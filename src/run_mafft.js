//@ts-check

const fs = require("fs");
const child_process = require("child_process");
const Path = require("path");

const { argv_parse, loadChrLength, array_groupBy } = require("./util.js");
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
const maxiterate = dataset.mafft.maxIterate || 1000;
const num_thread = dataset.mafft.thread || 1;
const mafft_path = dataset.mafft.path || "mafft";

const ref1_name = dataset.ref;
const ref2_name = dataset.parental_list[1];
const spore_1_name = dataset.progeny_list[0];
const spore_2_name = dataset.progeny_list[1];
const spore_3_name = dataset.progeny_list[2];
const spore_4_name = dataset.progeny_list[3];

const ref1_chr_list = loadChrLength(`output/${ref1_name}.length.txt`).list;
const ref2_chr_list = loadChrLength(`output/${ref2_name}.length.txt`).list;
const s1_chr_list = loadChrLength(`output/${spore_1_name}.length.txt`).list;
const s2_chr_list = loadChrLength(`output/${spore_2_name}.length.txt`).list;
const s3_chr_list = loadChrLength(`output/${spore_3_name}.length.txt`).list;
const s4_chr_list = loadChrLength(`output/${spore_4_name}.length.txt`).list;

if (!fs.existsSync("tmp/mafft_seq_frag")) {
	fs.mkdirSync("tmp/mafft_seq_frag");
}
if (!fs.existsSync("tmp/merged_fa")) {
	fs.mkdirSync("tmp/merged_fa");
}

if (process.argv[1] == __filename) {
	main();
}

function main() {
	const all_chr_frag_list = loadFragIdList(dataset, true);
	
	for (let nChr = 1; nChr <= ref1_chr_list.length; ++nChr) {
		loop_mafft(all_chr_frag_list, nChr);

		loop_cen_AT_mafft(all_chr_frag_list, nChr);
	}
}

function loop_mafft(all_chr_frag_list, nChr) {
	let list = all_chr_frag_list[nChr];
	list.filter(coord => !coord.centromere).forEach(coord => {
		if (coord.list.length) {
			console.log("ch", nChr, "is centromere ??", coord);
		}
		let r = run_multialign(nChr, coord.id, algorithm);
		if (!r) {
			r = run_multialign(nChr, coord.id, default_algorithm);
			if (!r) {
				console.log("error", algorithm, "nChr", nChr, "fileId", coord.id);
			}
		}
	});
}

/**
 * @param {number} nChr
 * @param {number} fragId
 * @param {string} [_algorithm]
 * @returns {boolean}
 */
function run_multialign(nChr, fragId, _algorithm = "") {
	const fasta_filename = `ch${nChr}_${fragId}.fa`;
	const input_path = `tmp/seq_frag/${fasta_filename}`;

	if (fs.existsSync(input_path)) {
		const output_file = `tmp/mafft_seq_frag/mafft_${fasta_filename}`;

		if (fs.existsSync(output_file)) {
			console.error("skip exist", output_file);
			return true;
		}
		else {
			run_mafft_(input_path, output_file, _algorithm);
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
			let task1 = run_mafft_async(input_ref1_path, output_ref1_file, _algorithm);
			let task2 = run_mafft_async(input_ref2_path, output_ref2_file, _algorithm);
			//await Promise.all([task1, task2]);
		}
	}

	return false;
}

/**
 * @param {string} input_path
 * @param {string} output_file
 * @param {string} [_algorithm]
 */
function run_mafft_(input_path, output_file, _algorithm = "") {
	let arg_entries = [
		["--thread", num_thread],
		["--maxiterate", maxiterate],
	];

	let mafft_cmd = `${mafft_path} --quiet ${_algorithm} ${input_path} > ${output_file}`;

	console.log(mafft_cmd);

	let time1 = new Date();
	try {
		child_process.execSync(mafft_cmd);

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

/**
 * @param {string} input_path
 * @param {string} output_file
 * @param {string} [_algorithm]
 */
async function run_mafft_async(input_path, output_file, _algorithm = "") {
	let mafft_cmd = `${mafft_path} --quiet --thread ${num_thread} ${_algorithm} --maxiterate ${maxiterate} ${input_path} > ${output_file}`;

	console.log(mafft_cmd);

	return child_process.exec(mafft_cmd);
}

/**
 * @param {string} filename
 * @returns {{[nChr:number]:{start:number,end:number,length:number}[]}}
 */
function load_AT_island(filename, filter) {
	const text_tab = fs.readFileSync(filename).toString();
	let table = table_to_object_list(tsv_parse(text_tab), ["chr", "start", "end", "length"]);
	
	/** @type {{[nChr:number]:{start:number,end:number,length:number}[]}} */
	let group = {};

	//group by chr
	table.forEach(row => {
		if (!group[row.chr]) {
			group[row.chr] = [];
		}
		let data = {
			start: Number(row.start),
			end: Number(row.end),
			length: Number(row.length),
		};
		if (!filter || filter(data)) {
			group[Number(row.chr)].push(data);
		}
	});

	return group;
}

function loop_cen_AT_mafft(all_chr_frag_list, nChr) {
	let list = all_chr_frag_list[nChr];
	list.filter(coord => coord.centromere).forEach(cen_coord => {
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
			exec_blastn(src_filename, src_filename).then(function (result_text) {
				const rows = parseBlastnResults(result_text);
				
				const ref1_name = ref1_chr_list[nChr - 1].chr;
				const ref2_name = ref2_chr_list[nChr - 1].chr;
				const s1_name = s1_chr_list[nChr - 1].chr;
				const s2_name = s2_chr_list[nChr - 1].chr;
				const s3_name = s3_chr_list[nChr - 1].chr;
				const s4_name = s4_chr_list[nChr - 1].chr;

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

						let args = [
							"--quiet",
							"--thread", String(num_thread),
							"--maxiterate", String(maxiterate),
							qs_file_path,
						];
						let mafft_cmd = `${mafft_path} ${args.join(" ")} > ${output_qs_filename}`;
						console.log(mafft_cmd);

						child_process.execSync(mafft_cmd);
					}
					if (!fs.existsSync(output_cs_filename)) {
						const cs = {
							[ref2_name]: merged_fasta[ref2_name],
							[ss[1][0]]: merged_fasta[ss[1][0]],
							[ss[1][1]]: merged_fasta[ss[1][1]],
						};
						saveFasta(cs_file_path, cs);
						
						let args = [
							"--quiet",
							"--thread", String(num_thread),
							"--maxiterate", String(maxiterate),
							cs_file_path,
						];
						let mafft_cmd = `${mafft_path} ${args.join(" ")} > ${output_cs_filename}`;
						console.log(mafft_cmd);

						child_process.execSync(mafft_cmd);
					}
				}
			});
		}//qs_filename, cs_filename
		else {
			console.log("done", nChr, {
				qs_filename: qs_file_path, cs_filename: cs_file_path, output_qs_filename, output_cs_filename
			});
		}

		//next step
		// ...
	});
}

