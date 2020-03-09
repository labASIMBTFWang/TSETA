//@ts-check

const fs = require("fs");
const Path = require("path");
const child_process = require("child_process");
const OS = require("os");

const { argv_parse, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");

const { BlastnCoord, execAsync, exec_blastn, parseBlastnResults, blastn_coord } = require("./blastn_util.js");
const { readFasta, saveFasta, parseFasta, joinFastaSeq } = require("./fasta_util.js");
const { loadFragIdList } = require("./load_frag_list.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const algorithm = dataset.mafft.algorithm || "";
const default_algorithm = dataset.mafft.default_algorithm || "";
const maxiterate = dataset.mafft.maxIterate;
const num_thread = dataset.mafft.thread;
const mafft_path = dataset.mafft.path || "mafft";

const genome_info_list = dataset.loadGenomeInfoList();
const chr_info_list = genome_info_list.map(gInfo => gInfo.chr_list);

if (!fs.existsSync("tmp/mafft_seq_frag")) {
	fs.mkdirSync("tmp/mafft_seq_frag");
}
if (!fs.existsSync("tmp/merged_fa")) {
	fs.mkdirSync("tmp/merged_fa");
}

if (process.argv[1] == __filename) {
	main();
}

async function main() {
	const all_chr_frag_list = loadFragIdList(dataset, true);
	
	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		let chr_frag_list = all_chr_frag_list[nChr];
		if (dataset.mode == "tetrad") {
			await loop_mafft(chr_frag_list.filter(coord => !coord.centromere), nChr);
			await loop_cen_AT_mafft(all_chr_frag_list, nChr);
		}
		else {
			await loop_mafft(chr_frag_list, nChr);
		}
	}
}

async function loop_mafft(chr_frag_list, nChr) {
	let tasks = chr_frag_list.map(async function (coord) {
		if (coord.list.length) {
			console.log("ch", nChr, "is centromere ??", coord);
		}
		let r = await run_multialign(nChr, coord.id, algorithm);
		if (!r) {
			r = await run_multialign(nChr, coord.id, default_algorithm);
			if (!r) {
				console.log("error", algorithm, "nChr", nChr, "fileId", coord.id);
			}
		}
	});
	await Promise.all(tasks);
}

/**
 * @param {number} nChr
 * @param {number} fragId
 * @param {string} [_algorithm]
 * @returns {Promise<boolean>}
 */
async function run_multialign(nChr, fragId, _algorithm) {
	const fasta_filename = `ch${nChr}_${fragId}.fa`;
	const input_path = `tmp/seq_frag/${fasta_filename}`;

	if (fs.existsSync(input_path)) {
		const output_file = `tmp/mafft_seq_frag/mafft_${fasta_filename}`;

		if (fs.existsSync(output_file)) {
			console.error("skip exist", output_file);
			return true;
		}
		else {
			return await run_mafft_(input_path, output_file, _algorithm, nChr, fragId);
		}
	}
	else {
		const ref1_filename = `ch${nChr}_${fragId}_ref1.fa`;
		const ref2_filename = `ch${nChr}_${fragId}_ref2.fa`;
		const input_ref1_path = `tmp/seq_frag/${ref1_filename}`;
		const input_ref2_path = `tmp/seq_frag/${ref2_filename}`;
		const output_ref1_file = `tmp/mafft_seq_frag/mafft_${ref1_filename}`;
		const output_ref2_file = `tmp/mafft_seq_frag/mafft_${ref2_filename}`;

		if (fs.existsSync(input_ref1_path) && fs.existsSync(input_ref2_path)) {
			let task1 = run_mafft_(input_ref1_path, output_ref1_file, _algorithm, nChr, fragId);
			let task2 = run_mafft_(input_ref2_path, output_ref2_file, _algorithm, nChr, fragId);
			await Promise.all([task1, task2]);
		}
	}

	return false;
}

/**
 * @param {string} input_path
 * @param {string} output_file
 * @param {string} _algorithm
 * @param {number} nChr
 * @param {number} fragId
 */
async function run_mafft_(input_path, output_file, _algorithm, nChr, fragId) {
	const arg_thread = num_thread >= 1 ? `--thread ${num_thread}` : "";
	const arg_maxiterate = maxiterate == null || maxiterate >= 0 ? `--maxiterate ${maxiterate}` : "";

	const mafft_cmd = `${mafft_path} --quiet ${arg_thread} ${_algorithm} ${arg_maxiterate} ${input_path} > ${output_file}`;

	console.log(mafft_cmd);

	let time1 = new Date();
	try {
		let proc = child_process.exec(mafft_cmd);

		let interval_id = setInterval(function () {
			try {
				let time2 = new Date();
				let time_elapsed = time2.getTime() - time1.getTime();
				let text = `${time_elapsed}\t${os_totalmem - OS.freemem()}`;
				fs.writeFileSync(`tmp/mem_usage/${nChr}_${fragId}.txt`, text, { flag: "a" });
			}
			catch (Ex) {
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

		clearInterval(interval_id);

		let time2 = new Date();
		let time_elapsed = time2.getTime() - time1.getTime();

		let log_text = JSON.stringify({
			start_time: time1,
			elapsed: time_elapsed,
			cmd: mafft_cmd,
			input: input_path,
		}) + "\n";
		try {
			fs.writeFileSync("tmp/loop-ma-log.txt", log_text, { flag: "a" });
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
			fs.writeFileSync("tmp/loop-ma-log.txt", log_text, { flag: "a" });//append to log file
		}
		catch (ex) {
			console.error("ex error", ex);
		}
		console.log("err", time_elapsed);
	}
}


async function loop_cen_AT_mafft(all_chr_frag_list, nChr) {
	let list = all_chr_frag_list[nChr];
	let chr_tasks = list.filter(coord => coord.centromere).map(async function (cen_coord) {
		let cen_fragId_list = cen_coord.list.map(a => a.id);

		cen_coord.list.forEach(a => console.log(": cen", "ch", nChr, "frag", a.id));

		const merged_fasta = joinFastaSeq(cen_fragId_list.map(a => `tmp/seq_frag/ch${nChr}_${a}.fa`).map(a => readFasta(a)));

		const src_filename = `tmp/merged_fa/ch${nChr}_${cen_fragId_list[0]}_${cen_fragId_list[cen_fragId_list.length - 1]}.fa`;
		if (!fs.existsSync(src_filename)) {
			saveFasta(src_filename, merged_fasta);
		}
		
		const qs_file_name = `ch${nChr}_${cen_fragId_list[0]}_${cen_fragId_list[cen_fragId_list.length - 1]}_ref1.fa`;
		const cs_file_name = `ch${nChr}_${cen_fragId_list[0]}_${cen_fragId_list[cen_fragId_list.length - 1]}_ref2.fa`;
		const qs_file_path = `tmp/merged_fa/${qs_file_name}`;
		const cs_file_path = `tmp/merged_fa/${cs_file_name}`;
		const output_qs_filename = `tmp/mafft_seq_frag/mafft_${qs_file_name}`;
		const output_cs_filename = `tmp/mafft_seq_frag/mafft_${cs_file_name}`;

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
				let tasks = [];
				if (!fs.existsSync(output_qs_filename)) {
					const qs = {
						[ref1_name]: merged_fasta[ref1_name],
						[ss[0][0]]: merged_fasta[ss[0][0]],
						[ss[0][1]]: merged_fasta[ss[0][1]],
					};
					saveFasta(qs_file_path, qs);

					tasks[0] = run_mafft_(qs_file_path, output_qs_filename, _algorithm, nChr, fragId);
				}
				if (!fs.existsSync(output_cs_filename)) {
					const cs = {
						[ref2_name]: merged_fasta[ref2_name],
						[ss[1][0]]: merged_fasta[ss[1][0]],
						[ss[1][1]]: merged_fasta[ss[1][1]],
					};
					saveFasta(cs_file_path, cs);
					
					tasks[1] = run_mafft_(qs_file_path, output_qs_filename, _algorithm, nChr, fragId);
				}
				await Promise.all(tasks);
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

