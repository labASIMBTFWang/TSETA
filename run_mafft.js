//@ts-check

const fs = require("fs");
const child_process = require("child_process");
const Path = require("path");

const { argv_parse, loadChrLength, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");

const { BlastnCoord, execAsync, exec_blastn, parseBlastnResults, blastn_coord } = require("./blastn_util.js");
const { readFasta, saveFasta, parseFasta, joinFastaSeq } = require("./fasta_util.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const start_nChr = 1;
const end_nChr = 7;

const input_directory = String(argv["-i"] || "./tmp/seq_frag");
const mafft_output_directory = String(argv["-o"] || "./tmp/mafft_seq_frag");
const tmp_merged_fa = "./tmp/merged_fa";

const filePath_multi_coord_prefix = "tmp/";

const algorithm = dataset["variant calling"].args.algorithm || "--localpair";
const default_algorithm = dataset["variant calling"].args.default_algorithm || "";
const maxiterate = dataset["variant calling"].args.maxIterate || 1000;
const num_thread = dataset["variant calling"].args.thread || 20;

const ref1_name = dataset.ref;
const ref2_name = dataset.refs[1];
const spore_1_name = dataset.progeny_list[0];
const spore_2_name = dataset.progeny_list[1];
const spore_3_name = dataset.progeny_list[2];
const spore_4_name = dataset.progeny_list[3];

if (!fs.existsSync(mafft_output_directory)) {
	fs.mkdirSync(mafft_output_directory);
}
if (!fs.existsSync(tmp_merged_fa)) {
	fs.mkdirSync(tmp_merged_fa);
}

const ref1_chr_list = loadChrLength(`./${ref1_name}.length.txt`).list;
const ref2_chr_list = loadChrLength(`./${ref2_name}.length.txt`).list;
const s1_chr_list = loadChrLength(`./${spore_1_name}.length.txt`).list;
const s2_chr_list = loadChrLength(`./${spore_2_name}.length.txt`).list;
const s3_chr_list = loadChrLength(`./${spore_3_name}.length.txt`).list;
const s4_chr_list = loadChrLength(`./${spore_4_name}.length.txt`).list;

class MyCoord {
	constructor() {
		this.id = "";

		this.start = {
			search: 0,
			r1: 0,
			r2: 0,
			s1: 0,
			s2: 0,
			s3: 0,
			s4: 0
		};

		this.end = {
			search: 0,
			r1: 0,
			r2: 0,
			s1: 0,
			s2: 0,
			s3: 0,
			s4: 0
		};
		
		this.removed = false;

		/** @type {MyCoord[]} */
		this.list = [];

		this.centromere = false;
	}

	get length() {
		return {
			search: this.end.search - this.start.search,
			r1: this.end.r1 - this.start.r1,
			r2: this.end.r2 - this.start.r2,
			s1: this.end.s1 - this.start.s1,
			s2: this.end.s2 - this.start.s2,
			s3: this.end.s3 - this.start.s3,
			s4: this.end.s4 - this.start.s4
		};
	}
}

if (process.argv[1] == __filename) {
	main();
}

function main() {
	const all_chr_frag_list = load_frag_id_list();

	const stop_nChr = Math.min(end_nChr, ref1_chr_list.length);
	
	for (let nChr = start_nChr; nChr <= stop_nChr; ++nChr) {
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
	const input_path = Path.join(input_directory, fasta_filename);

	if (fs.existsSync(input_path)) {
		const output_file = Path.join(mafft_output_directory, `mafft_${fasta_filename}`);

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
		const input_ref1_path = Path.join(input_directory, ref1_filename);
		const input_ref2_path = Path.join(input_directory, ref2_filename);
		const output_ref1_file = Path.join(mafft_output_directory, `mafft_${ref1_filename}`);
		const output_ref2_file = Path.join(mafft_output_directory, `mafft_${ref2_filename}`);

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
	let mafft_cmd = `mafft --quiet --thread ${num_thread} ${_algorithm} --maxiterate ${maxiterate} ${input_path} > ${output_file}`;

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
	let mafft_cmd = `mafft --quiet --thread ${num_thread} ${_algorithm} --maxiterate ${maxiterate} ${input_path} > ${output_file}`;

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

/** @returns {{[nChr:number]:MyCoord[]}} */
function load_frag_id_list() {
	const stop_nChr = Math.min(end_nChr, ref1_chr_list.length);
	const AT_island = load_AT_island(dataset["AT-island"], data => data.length >= 3000);

	/** @type {{[nChr:number]:MyCoord[]}} */
	const all_chr_frag_list = {};

	for (let nChr = start_nChr; nChr <= stop_nChr; ++nChr) {
		try {
			const coords = load_ma_coord(`${filePath_multi_coord_prefix}multi_coord_ch${nChr}.txt`);

			const AT_desc_list = AT_island[nChr].sort((a, b) => b.length - a.length);
			const cen_range = AT_desc_list[0];
			{
				let [at1, at2] = [AT_desc_list[0], AT_desc_list[1]].sort((a, b) => a.start - b.start);
				console.log("ch", nChr, "at1, 500bp, at2", at1.start, at1.end, at2.start, at2.end);
				// 合併 2 個鄰近的 AT island，2 個 AT island 間最多能有 1 個 window (QM6a ChIV: 1482500-1559500,1560000-1659000)
				if ((at2.start - at1.end) <= Math.abs(at1.end - at1.start)) {
					cen_range.start = at1.start;
					cen_range.end = at2.end;
					cen_range.length = at2.end - at1.start;
				}
			}
			console.log("cen", cen_range.start, cen_range.end, cen_range.length);

			let cen_fragId_list = Object.keys(coords).filter(id => {
				let coord = coords[id];
				if (cen_range.start <= coord.start.r1 && coord.end.r1 <= cen_range.end) {
					return true;
				}
				// else if (coord.end.r1 > cen_range.end && (coord.end.r1 - cen_range.end) <= 5000) {
				// 	return true;
				// }
				else {
					return false;
				}
			});
			if (cen_fragId_list.length <= 0) {
				cen_fragId_list = Object.keys(coords).filter(id => {
					let coord = coords[id];
					if (coord.start.r1 <= cen_range.start && cen_range.end <= coord.end.r1) {
						return true;
					}
					else {
						return false;
					}
				});
				if (cen_fragId_list.length <= 0) {
					cen_fragId_list = Object.keys(coords).filter(id => {
						let coord = coords[id];
						if (coord.start.r1 <= cen_range.end && coord.end.r1 >= cen_range.start) {
							return true;
						}
						else {
							return false;
						}
					});
				}
			}

			if (cen_fragId_list.length > 1) {
				let { r1_len, r2_len, s1_len, s2_len, s3_len, s4_len } = cen_fragId_list.reduce((prev, id) => {
					let rs_length = coords[id].length;
					let [r1_len, r2_len, s1_len, s2_len, s3_len, s4_len] = [rs_length.r1, rs_length.r2, rs_length.s1, rs_length.s2, rs_length.s3, rs_length.s4];
					return {
						r1_len: prev.r1_len + r1_len,
						r2_len: prev.r2_len + r2_len,
						s1_len: prev.s1_len + s1_len,
						s2_len: prev.s2_len + s2_len,
						s3_len: prev.s3_len + s3_len,
						s4_len: prev.s4_len + s4_len
					};
				}, { r1_len: 0, r2_len: 0, s1_len: 0, s2_len: 0, s3_len: 0, s4_len: 0 });

				let s_len = [s1_len, s2_len, s3_len, s4_len];
				let qs_len = s_len.map(a => Number((a / r1_len).toFixed(1)));
				let cs_len = s_len.map(a => Number((a / r2_len).toFixed(1)));
				
				let qs_c = qs_len.reduce((prev, curr) => prev + (curr == 1 ? 1 : 0), 0);
				let cs_c = cs_len.reduce((prev, curr) => prev + (curr == 1 ? 1 : 0), 0);

				if (qs_c >= 2 && cs_c >= 2) {
					console.log("cen", "ch", nChr, "qs_c", qs_c);
					console.log("cen", "ch", nChr, "cs_c", cs_c);
				}
				else {
					let l_id = cen_fragId_list[cen_fragId_list.length - 1];
					let keys = Object.keys(coords);
					let next_id = keys[keys.indexOf(l_id) + 1];

					cen_fragId_list.push(next_id);
					console.log("cen", "ch", nChr, "next_id", next_id);
				}
			}

			cen_fragId_list.forEach(id => {
				coords[id].removed = true;
				console.log("cen", "ch", nChr, "frag", id);
			});

			coords[cen_fragId_list[0]].list = cen_fragId_list.map(id => coords[id]);
			coords[cen_fragId_list[0]] = Object.assign({}, coords[cen_fragId_list[0]]);//clone

			coords[cen_fragId_list[0]].removed = false;
			coords[cen_fragId_list[0]].id = `${cen_fragId_list[0]}_${cen_fragId_list[cen_fragId_list.length - 1]}`;
			coords[cen_fragId_list[0]].centromere = true;

			all_chr_frag_list[nChr] = Object.keys(coords).map(id => coords[id]).filter(coord => {
				if (coord.id && !coord.removed) {
					return true;
				}
				else {
					return false;
				}
			});
		}
		catch (ex) {
			console.error(ex);
		}
	}

	return all_chr_frag_list;
}

function loop_cen_AT_mafft(all_chr_frag_list, nChr) {
	let list = all_chr_frag_list[nChr];
	list.filter(coord => coord.centromere).forEach(cen_coord => {
		let cen_fragId_list = cen_coord.list.map(a => a.id);

		cen_coord.list.forEach(a => console.log(": cen", "ch", nChr, "frag", a.id));

		const merged_fasta = joinFastaSeq(cen_fragId_list.map(a => `${input_directory}/ch${nChr}_${a}.fa`).map(a => readFasta(a)));

		const src_filename = Path.join(tmp_merged_fa, `ch${nChr}_${cen_fragId_list[0]}_${cen_fragId_list[cen_fragId_list.length - 1]}.fa`);
		if (!fs.existsSync(src_filename)) {
			saveFasta(src_filename, merged_fasta);
		}
		
		const qs_file_name = `ch${nChr}_${cen_fragId_list[0]}_${cen_fragId_list[cen_fragId_list.length - 1]}_ref1.fa`;
		const cs_file_name = `ch${nChr}_${cen_fragId_list[0]}_${cen_fragId_list[cen_fragId_list.length - 1]}_ref2.fa`;
		const qs_file_path = Path.join(tmp_merged_fa, qs_file_name);
		const cs_file_path = Path.join(tmp_merged_fa, cs_file_name);
		const output_qs_filename = Path.join(mafft_output_directory, `mafft_${qs_file_name}`);
		const output_cs_filename = Path.join(mafft_output_directory, `mafft_${cs_file_name}`);

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
						let mafft_cmd = `mafft ${args.join(" ")} > ${output_qs_filename}`;
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
						let mafft_cmd = `mafft ${args.join(" ")} > ${output_cs_filename}`;
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

function load_ma_coord(filename) {
	const text = fs.readFileSync(filename).toString();
	let rows = text.trim().split("\n").map(a => a.split(/\t/).map(b => b.trim()).filter(b => b != "|"));
	/** @type {{[id:string]:MyCoord}} */
	let map = {};
	rows.forEach(row => {
		let [id, type, search, r1, r2, s1, s2, s3, s4] = row;
		if (!map[id]) {
			map[id] = Object.assign(new MyCoord(), {
				id: id,
			});
		}
		map[id][type] = {
			search, r1, r2, s1, s2, s3, s4,
		};
	});
	
	return map;
}

