//@ts-check

const fs = require("fs");
const Path = require("path");
const child_process = require("child_process");
const OS = require("os");

const os_totalmem = OS.totalmem();

const { argv_parse, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");

const { BlastnCoord, execAsync, exec_blastn, parseBlastnResults, blastn_coord } = require("./blastn_util.js");
const { readFasta, saveFasta, parseFasta, joinFastaSeq } = require("./fasta_util.js");
const { loadFragIdList } = require("./load_frag_list.js");
const { Dataset, MafftOptions } = require("./dataset.js");
const { loadSetting } = require("./setting.js");
const setting = loadSetting();

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

let {
	algorithm,
	default_algorithm,
	maxIterate: maxiterate,
	thread: num_thread,
} = mafftOptionFromDataset(dataset) || mafftOptionFromArgv(argv);

dataset.mafft = {};
dataset.mafft.algorithm = algorithm;
dataset.mafft.default_algorithm = default_algorithm;
dataset.mafft.maxIterate = maxiterate;
dataset.mafft.thread = num_thread;
fs.writeFileSync(argv_dataset_path, JSON.stringify(dataset, null, "\t"));

/**
 * @param {Dataset} dataset
 */
function mafftOptionFromDataset(dataset) {
	if (dataset.mafft) {
		let algorithm = dataset.mafft.algorithm || "";
		let default_algorithm = dataset.mafft.default_algorithm || "";
		let maxiterate = dataset.mafft.maxIterate;
		let num_thread = dataset.mafft.thread;
		return { algorithm, default_algorithm, maxiterate, num_thread };
	}
	else {
		return null;
	}
}
/**
 * @param {{[x:string]:string|boolean}} argv
 */
function mafftOptionFromArgv(argv) {
	//let mafft_options = new MafftOptions();
	let algorithm = String(argv["--mafft-algorithm"] || "");
	let default_algorithm = String(argv["--default-algorithm"] || "");
	let maxiterate = Number(argv["--mafft-maxiterate"]);
	let num_thread = Number(argv["--mafft-thread"]);
	if (!Number.isSafeInteger(maxiterate)) {
		maxiterate = null;
	}
	if (!Number.isSafeInteger(num_thread)) {
		num_thread = null;
	}
	return { algorithm, default_algorithm, maxiterate, num_thread };
}

const genome_info_list = dataset.loadGenomeInfoList();
const chr_info_list = genome_info_list.map(gInfo => gInfo.chr_list);

if (!fs.existsSync(`${dataset.tmp_path}/mafft_seq_frag`)) {
	fs.mkdirSync(`${dataset.tmp_path}/mafft_seq_frag`);
}
if (!fs.existsSync(`${dataset.tmp_path}/merged_fa`)) {
	fs.mkdirSync(`${dataset.tmp_path}/merged_fa`);
}
if (!fs.existsSync(`${dataset.tmp_path}/mem_usage`)) {
	fs.mkdirSync(`${dataset.tmp_path}/mem_usage`);
}

if (process.argv[1] == __filename) {
	main();
}

async function main() {
	const tetrad_analysis = dataset.mode == "tetrad";
	
	const all_chr_frag_list = loadFragIdList(dataset);

	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		let chr_frag_list = all_chr_frag_list[nChr];

		if (tetrad_analysis) {
			await loop_mafft(chr_frag_list.filter(coord => !coord.centromere), nChr);
			await loop_centromere_AT_mafft(all_chr_frag_list, nChr);
		}
		else {
			await loop_mafft(chr_frag_list, nChr);
		}
	}

	for (let i = 1; i <= genome_info_list[0].chr_list.length; ++i) {
		let input_path = `${dataset.tmp_path}/mafft_ch${i}.fa`;
		let output_path = `${dataset.output_path}/mafft_ch${i}.fa`;
		fs.createReadStream(input_path).pipe(fs.createWriteStream(output_path));//copy file
	}
}

async function loop_mafft(chr_frag_list, nChr) {
	for (let coord of chr_frag_list) {
		if (coord.list && coord.list.length) {
			console.log("ch", nChr, "is centromere ??", coord);
		}
		let r = await run_multialign(nChr, coord.id, algorithm);
		if (!r) {
			r = await run_multialign(nChr, coord.id, default_algorithm);
			if (!r) {
				console.log("error", algorithm, "nChr", nChr, "fileId", coord.id);
			}
		}
	}
}

/**
 * @param {number} nChr
 * @param {number} fragId
 * @param {string} [_algorithm]
 * @returns {Promise<boolean>}
 */
async function run_multialign(nChr, fragId, _algorithm) {
	const fasta_filename = `ch${nChr}_${fragId}.fa`;
	const input_path = `${dataset.tmp_path}/seq_frag/${fasta_filename}`;

	if (fs.existsSync(input_path)) {
		const output_file = `${dataset.tmp_path}/mafft_seq_frag/mafft_${fasta_filename}`;

		return await run_mafft(input_path, output_file, _algorithm, nChr, fragId);
	}
	else {
		const ref1_filename = `ch${nChr}_${fragId}_ref1.fa`;
		const ref2_filename = `ch${nChr}_${fragId}_ref2.fa`;
		const input_ref1_path = `${dataset.tmp_path}/seq_frag/${ref1_filename}`;
		const input_ref2_path = `${dataset.tmp_path}/seq_frag/${ref2_filename}`;
		const output_ref1_file = `${dataset.tmp_path}/mafft_seq_frag/mafft_${ref1_filename}`;
		const output_ref2_file = `${dataset.tmp_path}/mafft_seq_frag/mafft_${ref2_filename}`;

		if (fs.existsSync(input_ref1_path) && fs.existsSync(input_ref2_path)) {
			let result_1 = await run_mafft(input_ref1_path, output_ref1_file, _algorithm, nChr, fragId);
			let result_2 = await run_mafft(input_ref2_path, output_ref2_file, _algorithm, nChr, fragId);
			return result_1 && result_2;
		}
		else {
			// TODO: write log

			console.error({
				fragId: fragId,
				"ref1 file": fs.existsSync(input_ref1_path),
				"ref2 file": fs.existsSync(input_ref2_path),
			});
		}
	}

	return false;
}

/**
 * @param {string} input_path
 * @param {string} output_file
 * @param {string} _algorithm
 * @param {number} nChr
 * @param {number|string} fragId
 * @param {boolean} [reAlign]
 */
async function run_mafft(input_path, output_file, _algorithm, nChr, fragId, reAlign = false) {
	if (fs.existsSync(output_file)) {
		let stat = fs.statSync(output_file);
		if (stat.size == 0) {
			fs.unlinkSync(output_file);
			console.log("remove empty file", output_file);
		}
		else if (!reAlign) {
			console.log("skip exist", output_file);
			return true;
		}
	}
	
	const arg_algorithm = _algorithm ? `--${_algorithm}` : "";
	const arg_thread = num_thread >= 1 ? `--thread ${num_thread}` : "";
	const arg_maxiterate = maxiterate >= 0 ? `--maxiterate ${maxiterate}` : "";

	const mafft_cmd = `${setting.mafft_bin} --quiet ${arg_thread} ${arg_algorithm} ${arg_maxiterate} ${input_path} > ${output_file}`;

	console.log(mafft_cmd);

	let interval_id;

	let time1 = new Date();
	try {
		let proc = child_process.exec(mafft_cmd);

		interval_id = setInterval(function () {
			try {
				let time2 = new Date();
				let time_elapsed = time2.getTime() - time1.getTime();
				let text = `${time_elapsed}\t${os_totalmem - OS.freemem()}\n`;
				fs.writeFile(`${dataset.tmp_path}/mem_usage/${nChr}_${fragId}.txt`, text, { flag: "a" }, function (err) {
					void(err);
				});
			}
			catch (ex) {
			}
		}, 1000);

		let promise = new Promise(function (resolve, reject) {
			proc.on("exit", function (code, signal) {
				if (code || signal) {
					reject({
						code, signal,
					});
				}
				else {
					resolve();
				}
			});
		});
		await promise;

		let time2 = new Date();
		let time_elapsed = time2.getTime() - time1.getTime();

		let log_text = JSON.stringify({
			start_time: time1,
			elapsed: time_elapsed,
			cmd: mafft_cmd,
			input: input_path,
		}) + "\n";
		try {
			fs.writeFileSync(`${dataset.tmp_path}/loop-ma-log.txt`, log_text, { flag: "a" });
		}
		catch (ex) {
			console.error("error", ex, log_text);
		}
		console.log(time_elapsed);

		return true;
	}
	catch (ex) {
		let time2 = new Date();
		let time_elapsed = time2.getTime() - time1.getTime();
		let log_text = JSON.stringify({
			start_time: time1,
			elapsed: time_elapsed,
			cmd: mafft_cmd,
			error: ex.stack,
			input: input_path,
		}) + "\n";
		try {
			fs.writeFileSync(`${dataset.tmp_path}/loop-ma-log.txt`, log_text, { flag: "a" });//append to log file
		}
		catch (ex) {
			console.error("ex error", ex);
		}
		console.log("err", time_elapsed);

		throw ex;
	}
	finally {
		clearInterval(interval_id);
	}
}

/**
 * tetrad analysis
 * @param {*} all_chr_frag_list 
 * @param {number} nChr
 */
async function loop_centromere_AT_mafft(all_chr_frag_list, nChr) {
	let list = all_chr_frag_list[nChr];
	let chr_tasks = list.filter(coord => coord.centromere).map(async function (cen_coord) {
		let cen_fragId_list = cen_coord.list.map(a => a.id);

		cen_coord.list.forEach(a => console.log(": cen", "ch", nChr, "frag", a.id));

		/** @type {{[seqName:string]:string}} */
		const merged_fasta = joinFastaSeq(cen_fragId_list.map(a => `${dataset.tmp_path}/seq_frag/ch${nChr}_${a}.fa`).map(a => readFasta(a)));

		const src_filename = `${dataset.tmp_path}/merged_fa/ch${nChr}_${cen_fragId_list[0]}_${cen_fragId_list[cen_fragId_list.length - 1]}.fa`;
		if (!fs.existsSync(src_filename)) {
			saveFasta(src_filename, merged_fasta);
		}
		
		const qs_file_name = `ch${nChr}_${cen_fragId_list[0]}_${cen_fragId_list[cen_fragId_list.length - 1]}_ref1.fa`;
		const cs_file_name = `ch${nChr}_${cen_fragId_list[0]}_${cen_fragId_list[cen_fragId_list.length - 1]}_ref2.fa`;
		const qs_file_path = `${dataset.tmp_path}/merged_fa/${qs_file_name}`;
		const cs_file_path = `${dataset.tmp_path}/merged_fa/${cs_file_name}`;
		const output_qs_filename = `${dataset.tmp_path}/mafft_seq_frag/mafft_${qs_file_name}`;
		const output_cs_filename = `${dataset.tmp_path}/mafft_seq_frag/mafft_${cs_file_name}`;

		if (![qs_file_path, cs_file_path, output_qs_filename, output_cs_filename].every(a => fs.existsSync(a))) {
			const result_text = await exec_blastn(src_filename, src_filename);
			const rows = parseBlastnResults(result_text);
			
			const chr_name_list = chr_info_list.map(cInfo => cInfo[nChr - 1].chr);
			// const ref1_name = ref1_chr_list[nChr - 1].chr;
			// const ref2_name = ref2_chr_list[nChr - 1].chr;
			// const s1_name = s1_chr_list[nChr - 1].chr;
			// const s2_name = s2_chr_list[nChr - 1].chr;
			// const s3_name = s3_chr_list[nChr - 1].chr;
			// const s4_name = s4_chr_list[nChr - 1].chr;
			const [ref1_name, ref2_name, s1_name, s2_name, s3_name, s4_name] = chr_name_list;

			const q1 = rows.filter(a => a.query == ref1_name && a.subject == s1_name).sort((a, b) => b.score - a.score)[0];
			const q2 = rows.filter(a => a.query == ref1_name && a.subject == s2_name).sort((a, b) => b.score - a.score)[0];
			const q3 = rows.filter(a => a.query == ref1_name && a.subject == s3_name).sort((a, b) => b.score - a.score)[0];
			const q4 = rows.filter(a => a.query == ref1_name && a.subject == s4_name).sort((a, b) => b.score - a.score)[0];
			
			const c1 = rows.filter(a => a.query == ref2_name && a.subject == s1_name).sort((a, b) => b.score - a.score)[0];
			const c2 = rows.filter(a => a.query == ref2_name && a.subject == s2_name).sort((a, b) => b.score - a.score)[0];
			const c3 = rows.filter(a => a.query == ref2_name && a.subject == s3_name).sort((a, b) => b.score - a.score)[0];
			const c4 = rows.filter(a => a.query == ref2_name && a.subject == s4_name).sort((a, b) => b.score - a.score)[0];
			
			const s1 = q1 != null && c1 != null ? (q1.score > c1.score ? 0 : (c1.score > q1.score ? 1 : 2)) : (q1 != null ? 0 : (c1 != null ? 1 : 2));
			const s2 = q1 != null && c1 != null ? (q2.score > c2.score ? 0 : (c2.score > q2.score ? 1 : 2)) : (q2 != null ? 0 : (c2 != null ? 1 : 2));
			const s3 = q1 != null && c1 != null ? (q3.score > c3.score ? 0 : (c3.score > q3.score ? 1 : 2)) : (q3 != null ? 0 : (c3 != null ? 1 : 2));
			const s4 = q1 != null && c1 != null ? (q4.score > c4.score ? 0 : (c4.score > q4.score ? 1 : 2)) : (q4 != null ? 0 : (c4 != null ? 1 : 2));

			/** @type {string[][]} */
			const ss = [[], [], []];

			ss[s1].push(s1_name);
			ss[s2].push(s2_name);
			ss[s3].push(s3_name);
			ss[s4].push(s4_name);

			if (ss[2].length) {//if not 2:2
				console.log("not 2:2: ", src_filename, ss);
			}
			else {
				if (!fs.existsSync(output_qs_filename)) {
					const qs = {
						[ref1_name]: merged_fasta[ref1_name],
						[ss[0][0]]: merged_fasta[ss[0][0]],
						[ss[0][1]]: merged_fasta[ss[0][1]],
					};
					saveFasta(qs_file_path, qs);

					let result = await run_mafft(qs_file_path, output_qs_filename, algorithm, nChr, "centromere_ref1").catch(async function (err) {
						return await run_mafft(qs_file_path, output_qs_filename, default_algorithm, nChr, "centromere_ref1");
					});
					if (!result) {
						// TODO: write log
					}
				}
				if (!fs.existsSync(output_cs_filename)) {
					const cs = {
						[ref2_name]: merged_fasta[ref2_name],
						[ss[1][0]]: merged_fasta[ss[1][0]],
						[ss[1][1]]: merged_fasta[ss[1][1]],
					};
					saveFasta(cs_file_path, cs);
					
					let result = await run_mafft(cs_file_path, output_cs_filename, algorithm, nChr, "centromere_ref2").catch(async function (err) {
						return await run_mafft(cs_file_path, output_cs_filename, default_algorithm, nChr, "centromere_ref2");
					});
					if (!result) {
						// TODO: write log
					}
				}
				console.log("end cen ref1 ref2");
			}
		}//qs_filename, cs_filename
		else {
			console.log("done", nChr, {
				qs_filename: qs_file_path, cs_filename: cs_file_path, output_qs_filename, output_cs_filename
			});
		}

		//next step
		// ...
	});

	await Promise.all(chr_tasks);
}


module.exports.run_mafft = run_mafft;

