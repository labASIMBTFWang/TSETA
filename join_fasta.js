//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, loadChrLength, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { Dataset } = require("./dataset.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { loadFragIdList } = require("./load_frag_list.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

let mafft_output_directory = String(argv["-i"] || "./tmp/mafft_seq_frag");

const temp_direction_path = "./tmp";

const ref1_chr_list = loadChrLength(`./${dataset.ref}.length.txt`).list;

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
	const all_chr_frag_list = loadFragIdList(dataset, true);//load_frag_id_list();

	for (let nChr = 1; nChr <= ref1_chr_list.length; ++nChr) {
		if (all_chr_frag_list[nChr]) {
			//let id_list = all_chr_frag_list[nChr].map(a => a.id);
	
			join_chr_frag(nChr, all_chr_frag_list[nChr], Path.resolve(get_make_temp_output_chr_path(nChr)));
		}
		else {
			console.log("skip", "ch", nChr);
		}
	}
}

/**
 * @param {number} nChr
 */
function get_make_temp_output_chr_path(nChr) {
	return Path.join(temp_direction_path, get_make_output_chr_file_name(nChr));
}

/**
 * @param {number} nChr
 */
function get_make_output_chr_file_name(nChr) {
	return `mafft_ch${nChr}.fa`;
}

/** @returns {{[nChr:number]:MyCoord[]}} */
function load_frag_id_list() {
	const AT_island = load_AT_island(dataset["AT-island"], data => data.length >= 3000);

	/** @type {{[nChr:number]:MyCoord[]}} */
	const all_chr_frag_list = {};

	for (let nChr = 1; nChr <= ref1_chr_list.length; ++nChr) {
		try {
			const coords = load_ma_coord(Path.join(temp_direction_path, `multi_coord_ch${nChr}.txt`));

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
					if (coord.start.r1 <= cen_range.start && cen_range.end <= coord.end.r1) {//if frag.length >= cen // cen in frag
						return true;
					}
					else {
						return false;
					}
				});
				if (cen_fragId_list.length <= 0) {
					cen_fragId_list = Object.keys(coords).filter(id => {
						let coord = coords[id];
						if (coord.start.r1 <= cen_range.end && coord.end.r1 >= cen_range.start) {//border
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

/**
 * @param {number} nChr 
 * @param {MyCoord[]} coord_list 
 * @param {string} output_name 
 */
function join_chr_frag(nChr, coord_list, output_name) {
	let output_fa = {
	};
	
	coord_list.forEach(coord => {
		let in_filename = `${mafft_output_directory}/mafft_ch${nChr}_${coord.id}.fa`;

		if (!coord.centromere && fs.existsSync(in_filename)) {
			let in_fa = readFasta(in_filename);

			console.log("ch", nChr, "frag", coord.id);

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
					console.log("skip seq:", coord.id, seq_name);
				}
			});
		}
		else {
			let in_ref1_filename = `${mafft_output_directory}/mafft_ch${nChr}_${coord.id}_ref1.fa`;
			let in_ref2_filename = `${mafft_output_directory}/mafft_ch${nChr}_${coord.id}_ref2.fa`;

			if (fs.existsSync(in_ref1_filename) && fs.existsSync(in_ref2_filename)) {
				let in_ref1_fa = readFasta(in_ref1_filename);
				let in_ref2_fa = readFasta(in_ref2_filename);

				console.log("ch", nChr, "frag", coord.id, "ref1 ref2");
				
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
						console.log("skip seq:", coord.id, seq_name);
					}
				});
			}
			else {
				console.log("no file:", coord.id);
				//break;
			}
		}
	});

	saveFasta(output_name, output_fa);
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

