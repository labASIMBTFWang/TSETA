//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, loadChrLength, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { Dataset } = require("./dataset.js");
const { BlastnCoord, execAsync, exec_blastn, parseBlastnResults, blastn_coord, isCollide, groupByOverlap } = require("./blastn_util.js");
const { readFasta, saveFasta } = require("./fasta_util.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

let mafft_output_directory = String(argv["-i"] || "./tmp/mafft_seq_frag");

const temp_direction_path = "./tmp";


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
else {
	debugger;
}

function main() {
	merge_chr_all_fa();
}

function merge_chr_all_fa() {
	const all_chr_frag_list = load_frag_id_list();

	for (let nChr = 1; nChr <= dataset.chrNames.length; ++nChr) {
		if (all_chr_frag_list[nChr]) {
			let id_list = all_chr_frag_list[nChr].map(a => a.id);

			let genome_name_list = [
				dataset.ref,
				dataset.refs[1],
				dataset.progeny_list[0],
				dataset.progeny_list[1],
				dataset.progeny_list[2],
				dataset.progeny_list[3],
			];
			let raw_seq = genome_name_list.map((genome_name, i) => {
				let chr_seq = readFasta(Path.resolve("./", `${genome_name}.genome.fa`));
				return chr_seq[id_list[i]];
			});
			
			let output_fa = {};
			id_list.forEach(i => {
				let in_filename = `${"./tmp/seq_frag"}/ch${nChr}_${i}.fa`;
		
				let in_fa;

				if (fs.existsSync(in_filename)) {
					in_fa = readFasta(in_filename);
		
					//console.log("ch", nChr, "frag", i);
				}
				else {
					let in_ref1_filename = `${"./tmp/seq_frag"}/ch${nChr}_${i}_ref1.fa`;
					let in_ref2_filename = `${"./tmp/seq_frag"}/ch${nChr}_${i}_ref2.fa`;
		
					if (fs.existsSync(in_ref1_filename) && fs.existsSync(in_ref2_filename)) {
						let in_ref1_fa = readFasta(in_ref1_filename);
						let in_ref2_fa = readFasta(in_ref2_filename);
		
						//console.log("ch", nChr, "frag", i, "ref1 ref2");

						in_fa = Object.assign({}, in_ref1_fa, in_ref2_fa);
					}
					else {
						console.log({
							in_filename, in_ref1_filename, in_ref2_filename,
						});
						throw new Error("ref1 ?? ref2 ??");
					}
				}

				Object.keys(in_fa).forEach(seq_name => {
					let seq = in_fa[seq_name];
					if (seq) {
						if (!output_fa[seq_name]) {
							output_fa[seq_name] = "";
						}
						output_fa[seq_name] += seq;
						//console.log("seq.len", seq_name, output_fa[seq_name].length, "+", seq.length);
					}
					else {
						//console.log("skip seq:", i, seq_name);
					}
				});
			});
			
			let a = genome_name_list.every((name, idx) => {
				if (output_fa[name] == raw_seq[idx]) {
					return true;
				}
				else {
					console.error("error", nChr, name);
				}
			});
			if (a) {
				console.log(nChr, "ok");
			}
		}
		else {
			console.log("skip", "ch", nChr);
		}
	}
}

/** @returns {{[nChr:number]:MyCoord[]}} */
function load_frag_id_list() {
	const AT_island = load_AT_island(dataset["AT-island"], data => data.length >= 3000);

	/** @type {{[nChr:number]:MyCoord[]}} */
	const all_chr_frag_list = {};

	for (let nChr = 1; nChr <= dataset.chrNames.length; ++nChr) {
		try {
			const coords = load_ma_coord(Path.join(temp_direction_path, `multi_coord_ch${nChr}.txt`));
			all_chr_frag_list[nChr] = Object.keys(coords).map(id => coords[id]);
		}
		catch (ex) {
			console.error(ex);
		}
	}

	return all_chr_frag_list;
}

//main();

/**
 * @param {number} nChr 
 * @param {string[]} file_id_list 
 * @param {string} output_name 
 */
function old_join_chr_frag(nChr, file_id_list, output_name) {
	let output_fa = {
	};
	
	file_id_list.forEach(i => {
		let in_filename = `${mafft_output_directory}/mafft_ch${nChr}_${i}.fa`;

		if (fs.existsSync(in_filename)) {
			let in_fa = readFasta(in_filename);

			console.log("ch", nChr, "frag", i);

			Object.keys(in_fa).forEach(seq_name => {
				let seq = in_fa[seq_name];
				if (seq) {
					if (!output_fa[seq_name]) {
						output_fa[seq_name] = "";
					}
					output_fa[seq_name] += seq;
					console.log("seq.len", seq_name, output_fa[seq_name].length, "+", seq.length);
				}
				else {
					console.log("skip seq:", i, seq_name);
				}
			});
		}
		else {
			let in_ref1_filename = `${mafft_output_directory}/mafft_ch${nChr}_${i}_ref1.fa`;
			let in_ref2_filename = `${mafft_output_directory}/mafft_ch${nChr}_${i}_ref2.fa`;

			if (fs.existsSync(in_ref1_filename) && fs.existsSync(in_ref2_filename)) {
				let in_ref1_fa = readFasta(in_ref1_filename);
				let in_ref2_fa = readFasta(in_ref2_filename);

				console.log("ch", nChr, "frag", i, "ref1 ref2");
				
				let ref1_max_length = Math.max(...Object.values(in_ref1_fa).map(seq => seq.length));
				let ref2_max_length = Math.max(...Object.values(in_ref2_fa).map(seq => seq.length));
				let multi_length = ref1_max_length + ref2_max_length;
				
				Object.keys(in_ref1_fa).forEach(seq_name => {
					in_ref1_fa[seq_name] = in_ref1_fa[seq_name].padEnd(ref1_max_length, "-");
				});
				Object.keys(in_ref2_fa).forEach(seq_name => {
					in_ref2_fa[seq_name] = in_ref2_fa[seq_name].padEnd(ref2_max_length, "-");
				});

				if (ref1_max_length >= ref2_max_length) {
					Object.keys(in_ref1_fa).forEach(seq_name => {
						in_ref1_fa[seq_name] = in_ref1_fa[seq_name].padStart(multi_length, "-");
					});
					Object.keys(in_ref2_fa).forEach(seq_name => {
						in_ref2_fa[seq_name] = in_ref2_fa[seq_name].padEnd(multi_length, "-");
					});
				}
				else {
					Object.keys(in_ref1_fa).forEach(seq_name => {
						in_ref1_fa[seq_name] = in_ref1_fa[seq_name].padEnd(multi_length, "-");
					});
					Object.keys(in_ref2_fa).forEach(seq_name => {
						in_ref2_fa[seq_name] = in_ref2_fa[seq_name].padStart(multi_length, "-");
					});
				}
				
				let in_fa = Object.assign({}, in_ref1_fa, in_ref2_fa);

				Object.keys(in_fa).forEach(seq_name => {
					let seq = in_fa[seq_name];
					if (seq) {
						if (!output_fa[seq_name]) {
							output_fa[seq_name] = "";
						}
						output_fa[seq_name] += seq;
						console.log("indel seq.len", seq_name, output_fa[seq_name].length, "+", seq.length);
					}
					else {
						console.log("skip seq:", i, seq_name);
					}
				});
			}
			else {
				console.log("no file:", i);
				//break;
			}
		}
	});

	saveFasta(output_name, output_fa);
}

// //let time_to_live = 5;
// for (let i = 1; i <= num_seq; ++i) {
// }

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

