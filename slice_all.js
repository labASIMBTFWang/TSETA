// @ts-check

const os = require("os");
const fs = require("fs");
const child_process = require("child_process");
const Path = require("path");

const { argv_parse, loadChrLength, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { readFasta, saveFasta, parseFasta } = require("./fasta_util.js");
const { BlastnCoord, execAsync, exec_blastn, parseBlastnResults, blastn_coord, isCollide } = require("./blastn_util.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const start_nChr = 1;

let argv_param_filename = (f => fs.existsSync(f) ? f : null)(process.argv[2]);
let _argv_param = argv_param_filename ? JSON.parse(fs.readFileSync(argv_param_filename).toString()) : {
	source_genome_seq_directory: Path.resolve("./"),
	output_dir: Path.resolve("./tmp"),
	source_chr_seq_directory: Path.resolve("./tmp/fasta"),
	output_frag_directory: Path.resolve("./tmp/seq_frag"),
	
	slice_size: 10000,
	init_max_delta: 10000 * 3,
	
	ref1_name: dataset.ref,
	ref2_name: dataset.refs[1],
	spore_1_name: dataset.progeny_list[0],
	spore_2_name: dataset.progeny_list[1],
	spore_3_name: dataset.progeny_list[2],
	spore_4_name: dataset.progeny_list[3],
};
argv_param_filename = "default.json";

let $save_argv_param_filename = false;
let argv_param = new Proxy(_argv_param, {
	set: function (target, propertyKey, value, receiver) {
		$save_argv_param_filename = true;
		//target[propertyKey] = value;
		return Reflect.set(target, propertyKey, value, receiver);
	},
});

{
	const header = ["Index#", "Chr#", "Length"].join("\t") + "\n";
	const genome_name_list = [
		argv_param.ref1_name,
		argv_param.ref2_name,
		argv_param.spore_1_name,
		argv_param.spore_2_name,
		argv_param.spore_3_name,
		argv_param.spore_4_name
	];
	
	genome_name_list.forEach(genome_name => {
		let fname = Path.join("./", `${genome_name}.length.txt`);
		if (!fs.existsSync(fname)) {
			let fa = readFasta(Path.join(argv_param.source_genome_seq_directory, `${genome_name}.genome.fa`));

			let out_text = "";
			out_text += header;
			out_text += Object.keys(fa).map((seq_name, idx) => [idx + 1, seq_name, fa[seq_name].length].join("\t")).join("\n");

			fs.writeFileSync(fname, out_text);

			console.log("output *.length.txt", fname);
		}
	})
}

const ref1_chr_list = loadChrLength(Path.join("./", `${argv_param.ref1_name}.length.txt`)).list;
const ref2_chr_list = loadChrLength(Path.join("./", `${argv_param.ref2_name}.length.txt`)).list;
const s1_chr_list = loadChrLength(Path.join("./", `${argv_param.spore_1_name}.length.txt`)).list;
const s2_chr_list = loadChrLength(Path.join("./", `${argv_param.spore_2_name}.length.txt`)).list;
const s3_chr_list = loadChrLength(Path.join("./", `${argv_param.spore_3_name}.length.txt`)).list;
const s4_chr_list = loadChrLength(Path.join("./", `${argv_param.spore_4_name}.length.txt`)).list;

{
	const genome_list = [
		[argv_param.ref1_name, ref1_chr_list, {}],
		[argv_param.ref2_name, ref2_chr_list, {}],
		[argv_param.spore_1_name, s1_chr_list, {}],
		[argv_param.spore_2_name, s2_chr_list, {}],
		[argv_param.spore_3_name, s3_chr_list, {}],
		[argv_param.spore_4_name, s4_chr_list, {}]
	];

	if (!fs.existsSync(argv_param.output_dir)) {
		fs.mkdirSync(argv_param.output_dir);
	}
	if (!fs.existsSync(Path.join(argv_param.output_dir, "fasta"))) {
		fs.mkdirSync(Path.join(argv_param.output_dir, "fasta"));
	}

	genome_list.forEach(([genome_name, chr_list, fa]) => {
		Object.assign(fa, readFasta(Path.join(argv_param.source_genome_seq_directory, `${genome_name}.genome.fa`)));
	});

	for (let nChr = start_nChr; nChr <= ref1_chr_list.length; ++nChr) {
		genome_list.forEach(([genome_name, chr_list, fa]) => {
			const chr_idx = nChr - 1;
			
			let chr_info = chr_list[chr_idx];
			
			let chr_name = chr_info.chr;
			
			let fname_fasta = Path.join(argv_param.source_chr_seq_directory, `${genome_name}_${chr_name}.fa`);
			
			if (!fs.existsSync(fname_fasta)) {
				console.log("output", genome_name, nChr, fname_fasta);
				saveFasta(fname_fasta, {
					[chr_name]: fa[chr_name],
				});
			}
		});
	}
}

if (!fs.existsSync(Path.join(argv_param.output_dir, "ma_util_blastn"))) {
	fs.mkdirSync(Path.join(argv_param.output_dir, "ma_util_blastn"));
}

if (process.argv[1] == __filename) {
	main();
}

async function main() {
	if (!fs.existsSync(argv_param.output_frag_directory)) {
		fs.mkdirSync(argv_param.output_frag_directory);
	}
	if (!fs.existsSync(Path.join(argv_param.output_dir, "./meta"))) {
		fs.mkdirSync(Path.join(argv_param.output_dir, "./meta"));
	}

	if ($save_argv_param_filename) {
		fs.writeFileSync(argv_param_filename, JSON.stringify(argv_param));
	}

	for (let nChr = start_nChr; nChr <= ref1_chr_list.length; ++nChr) {
		console.log("start ch", nChr);
		try {
			const count = await multi_coord(nChr);

			console.log(count);
		}
		catch (ex) {
			fs.writeFileSync(Path.join(argv_param.output_dir, "ma_util.error.txt"), "ch" + nChr + "\n" + ex.stack + "\n", { flag: "a" });
			console.error(ex);
		}
		console.log("end ch", nChr);
	}
}

/**
 * @param {number} nChr
 */
async function multi_coord(nChr) {
	const chr_idx = nChr - 1;
	let ref1_chr_info = ref1_chr_list[chr_idx];
	let ref2_chr_info = ref2_chr_list[chr_idx];
	let s1_chr_info = s1_chr_list[chr_idx];
	let s2_chr_info = s2_chr_list[chr_idx];
	let s3_chr_info = s3_chr_list[chr_idx];
	let s4_chr_info = s4_chr_list[chr_idx];

	let ref1_chr_name = ref1_chr_info.chr;
	let ref2_chr_name = ref2_chr_info.chr;
	let s1_chr_name = s1_chr_info.chr;
	let s2_chr_name = s2_chr_info.chr;
	let s3_chr_name = s3_chr_info.chr;
	let s4_chr_name = s4_chr_info.chr;

	let ref1_chr_fasta = Path.join(argv_param.source_chr_seq_directory, `${argv_param.ref1_name}_${ref1_chr_name}.fa`);
	let ref2_chr_fasta = Path.join(argv_param.source_chr_seq_directory, `${argv_param.ref2_name}_${ref2_chr_name}.fa`);
	let s1_chr_fasta = Path.join(argv_param.source_chr_seq_directory, `${argv_param.spore_1_name}_${s1_chr_name}.fa`);
	let s2_chr_fasta = Path.join(argv_param.source_chr_seq_directory, `${argv_param.spore_2_name}_${s2_chr_name}.fa`);
	let s3_chr_fasta = Path.join(argv_param.source_chr_seq_directory, `${argv_param.spore_3_name}_${s3_chr_name}.fa`);
	let s4_chr_fasta = Path.join(argv_param.source_chr_seq_directory, `${argv_param.spore_4_name}_${s4_chr_name}.fa`);

	let ref1_chr_seq = readFasta(ref1_chr_fasta)[ref1_chr_name];
	let ref2_chr_seq = readFasta(ref2_chr_fasta)[ref2_chr_name];
	let s1_chr_seq = readFasta(s1_chr_fasta)[s1_chr_name];
	let s2_chr_seq = readFasta(s2_chr_fasta)[s2_chr_name];
	let s3_chr_seq = readFasta(s3_chr_fasta)[s3_chr_name];
	let s4_chr_seq = readFasta(s4_chr_fasta)[s4_chr_name];

	let ref2_start = 1;
	let s1_start = 1;
	let s2_start = 1;
	let s3_start = 1;
	let s4_start = 1;

	let ref2_find_next_start = 0;
	let s1_find_next_start = 0;
	let s2_find_next_start = 0;
	let s3_find_next_start = 0;
	let s4_find_next_start = 0;

	let max_delta = argv_param.init_max_delta;

	let search_start = 1;
	let search_end = search_start + argv_param.slice_size;
	
	let ref1_start = search_start;
	let ref1_end = search_end;

	let $count = 1;

	let seq_from_ref1 = [true, null, null, null, null, null];

	let _global_search_align, _local_overlap_align;
	let ref2_search_start, s1_search_start, s2_search_start, s3_search_start, s4_search_start;

	for (; search_start <= ref1_chr_info.length && search_end <= ref1_chr_info.length; ++$count) {
		try {
			/**
			 * @param {BlastnCoord} coord_1
			 * @param {BlastnCoord} coord_2
			 */
			function check_q_start_end(coord_1, coord_2) {
				return coord_1.qstart <= coord_2.qend && coord_1.qend >= coord_2.qstart;
			}
			
			/**
			 * @param {BlastnCoord} ref2_coord
			 * @param {BlastnCoord} ref1_coord
			 */
			function check_q_s_start_end(ref2_coord, ref1_coord) {
				return (ref2_coord.qstart <= ref1_coord.qend && ref2_coord.qend >= ref1_coord.qstart);
			}
			
			/**
			 * @param {BlastnCoord} ref1_coord
			 */
			function check_q_no_repeats(ref1_coord) {
				return (ref1_coord.qstart == ref1_coord.sstart && ref1_coord.qend == ref1_coord.send);
			}

			/**
			 * @param {BlastnCoord} row
			 * @param {any[]} params
			 */
			function next_start_filter(row, params) {
				return (
					row.qstart < row.qend &&
					row.sstart < row.send
					// //row.sstart >= next_start
					// row.send >= params[0]
				);
			}

			if (search_start >=  21533 || search_end >=  21533) {
				debugger
			}

			ref2_search_start = Math.max(ref2_start, ref2_find_next_start);
			s1_search_start = Math.max(s1_start, s1_find_next_start);
			s2_search_start = Math.max(s2_start, s2_find_next_start);
			s3_search_start = Math.max(s3_start, s3_find_next_start);
			s4_search_start = Math.max(s4_start, s4_find_next_start);

			if (ref2_search_start > ref2_chr_info.length) {
				break;
			}
			if (s1_search_start > s1_chr_info.length) {
				break;
			}
			if (s2_search_start > s2_chr_info.length) {
				break;
			}
			if (s3_search_start > s3_chr_info.length) {
				break;
			}
			if (s4_search_start > s4_chr_info.length) {
				break;
			}

			let ref2_search_end = Math.min(ref2_search_start + max_delta, ref2_chr_info.length);
			let s1_search_end = Math.min(s1_search_start + max_delta, s1_chr_info.length);
			let s2_search_end = Math.min(s2_search_start + max_delta, s2_chr_info.length);
			let s3_search_end = Math.min(s3_search_start + max_delta, s3_chr_info.length);
			let s4_search_end = Math.min(s4_search_start + max_delta, s4_chr_info.length);
			
			// remove ref1-ref1 repeat and duplicate
			let ref1_coords_task = blastn_coord(ref1_chr_fasta, ref1_chr_fasta, search_start, search_end, search_start, search_end, next_start_filter, [], nChr + "_" + $count);
			
			let ref2_coords_task = blastn_coord(ref1_chr_fasta, ref2_chr_fasta, search_start, search_end, ref2_search_start, ref2_search_end, next_start_filter, [], nChr + "_" + $count);
			let s1_coords_task = blastn_coord(ref1_chr_fasta, s1_chr_fasta, search_start, search_end, s1_search_start, s1_search_end, next_start_filter, [], nChr + "_" + $count);
			let s2_coords_task = blastn_coord(ref1_chr_fasta, s2_chr_fasta, search_start, search_end, s2_search_start, s2_search_end, next_start_filter, [], nChr + "_" + $count);
			let s3_coords_task = blastn_coord(ref1_chr_fasta, s3_chr_fasta, search_start, search_end, s3_search_start, s3_search_end, next_start_filter, [], nChr + "_" + $count);
			let s4_coords_task = blastn_coord(ref1_chr_fasta, s4_chr_fasta, search_start, search_end, s4_search_start, s4_search_end, next_start_filter, [], nChr + "_" + $count);
			let [ref1_coords, ref2_coords, s1_coords, s2_coords, s3_coords, s4_coords] = _global_search_align = await Promise.all([ref1_coords_task, ref2_coords_task, s1_coords_task, s2_coords_task, s3_coords_task, s4_coords_task]);
			
			/**
			 * @param {BlastnCoord[]} ref1_coords
			 * @param {BlastnCoord[]} ref2_coords
			 * @param {BlastnCoord[]} s1_coords
			 * @param {BlastnCoord[]} s2_coords
			 * @param {BlastnCoord[]} s3_coords
			 * @param {BlastnCoord[]} s4_coords
			 */
			function check_overlap_ref2(ref1_coords, ref2_coords, s1_coords, s2_coords, s3_coords, s4_coords) {
				let _ref1_coords = ref1_coords.filter(r2c => r2c.align >= 1000);
				let _ref2_coords = ref2_coords.filter(r2c => r2c.align >= 1000);
				let _s1_coords = s1_coords.filter(r2c => r2c.align >= 1000);
				let _s2_coords = s2_coords.filter(r2c => r2c.align >= 1000);
				let _s3_coords = s3_coords.filter(r2c => r2c.align >= 1000);
				let _s4_coords = s4_coords.filter(r2c => r2c.align >= 1000);

				/**
				 * ch1,ch2: order by subject.end
				 * ch3: order by query.end
				 * new method: order 
				 */
				if (_ref2_coords && _ref2_coords.length) {
					let r2c = _ref2_coords.map(ref2_coord => {
						return {
							//191225//ref1_coord: ref1_coords.filter(coord2 => check_q_start_end(coord1, coord2)).sort((a, b) => b.send - a.send)[0],
							ref1_coord: _ref1_coords.filter(ref1_coord => check_q_no_repeats(ref1_coord)).filter(ref1_coord => check_q_s_start_end(ref2_coord, ref1_coord)).sort((a, b) => b.send - a.send)[0],//191225//
							ref2_coord: ref2_coord,
							s1_coord: _s1_coords.filter(coord2 => check_q_start_end(ref2_coord, coord2)).sort((a, b) => b.send - a.send)[0],
							s2_coord: _s2_coords.filter(coord2 => check_q_start_end(ref2_coord, coord2)).sort((a, b) => b.send - a.send)[0],
							s3_coord: _s3_coords.filter(coord2 => check_q_start_end(ref2_coord, coord2)).sort((a, b) => b.send - a.send)[0],
							s4_coord: _s4_coords.filter(coord2 => check_q_start_end(ref2_coord, coord2)).sort((a, b) => b.send - a.send)[0],
						};
					}).filter(a => a.ref1_coord && a.s1_coord && a.s2_coord && a.s3_coord && a.s4_coord).sort((a, b) => b.ref2_coord.send - a.ref2_coord.send);

					let all_match_group = _ref2_coords.map(ref2_coord => {
						return {
							ref1_coord: _ref1_coords.sort((a, b) => (Math.abs(a.qstart - ref2_coord.qstart) + Math.abs(a.qend - ref2_coord.qend)) - (Math.abs(b.qstart - ref2_coord.qstart) + Math.abs(b.qend - ref2_coord.qend)))[0],
							ref2_coord: ref2_coord,
							s1_coord: _s1_coords.sort((a, b) => (Math.abs(a.qstart - ref2_coord.qstart) + Math.abs(a.qend - ref2_coord.qend)) - (Math.abs(b.qstart - ref2_coord.qstart) + Math.abs(b.qend - ref2_coord.qend)))[0],
							s2_coord: _s2_coords.sort((a, b) => (Math.abs(a.qstart - ref2_coord.qstart) + Math.abs(a.qend - ref2_coord.qend)) - (Math.abs(b.qstart - ref2_coord.qstart) + Math.abs(b.qend - ref2_coord.qend)))[0],
							s3_coord: _s3_coords.sort((a, b) => (Math.abs(a.qstart - ref2_coord.qstart) + Math.abs(a.qend - ref2_coord.qend)) - (Math.abs(b.qstart - ref2_coord.qstart) + Math.abs(b.qend - ref2_coord.qend)))[0],
							s4_coord: _s4_coords.sort((a, b) => (Math.abs(a.qstart - ref2_coord.qstart) + Math.abs(a.qend - ref2_coord.qend)) - (Math.abs(b.qstart - ref2_coord.qstart) + Math.abs(b.qend - ref2_coord.qend)))[0],
						};
					}).sort((a, b) => b.ref2_coord.send - a.ref2_coord.send);

					return {
						best_match: all_match_group[0],
						all_match_groups: all_match_group,
					};
				}
				return {
					best_match: {
						ref1_coord: null,
						ref2_coord: null,
						s1_coord: null,
						s2_coord: null,
						s3_coord: null,
						s4_coord: null,
					},
					all_match_groups: [],
				};
			}

			let match_results = check_overlap_ref2(ref1_coords, ref2_coords, s1_coords, s2_coords, s3_coords, s4_coords);
			let { ref1_coord, ref2_coord, s1_coord, s2_coord, s3_coord, s4_coord } = _local_overlap_align = match_results.best_match;
			
			//check large In/Del
			{
				let _seq_from_ref1_p = [true, ref2_coord != null];
				let _seq_from_ref1_s = [s1_coord != null, s2_coord != null, s3_coord != null, s4_coord != null];
				
				let n_sp_p = _seq_from_ref1_p.reduce((prev, current) => prev + (current ? 1 : 0), 0);
				let n_sp_s = _seq_from_ref1_s.reduce((prev, current) => prev + (current ? 1 : 0), 0);

				if (n_sp_s != 2 && !(n_sp_p == 2 && n_sp_s == 4)) {
					//large 4:0 or 3:1
					let text = (JSON.stringify({
						$count,
						"Q:C": ["0:4", "1:3", "2:2", "3:1", "4:0"][n_sp_s],
						ref1_specific: n_sp_p == 1,
						seq_from_ref1: [true, ref2_coord != null, s1_coord != null, s2_coord != null, s3_coord != null, s4_coord != null],
						//
						search_start, search_end,
						ref2_start, ref2_find_next_start,
						s1_start, s1_find_next_start,
						s2_start, s2_find_next_start,
						s3_start, s3_find_next_start,
						s4_start, s4_find_next_start
					}));
					fs.writeFileSync(Path.join(argv_param.output_dir, "ref1_ch" + nChr + "_specific_seq.txt"), text + "\n", { flag: "a" });
				}
			}

			if (ref2_coord && s1_coord && s2_coord && s3_coord && s4_coord) {
				if ((ref2_coord.sstart - ref2_start) >= max_delta ||
					(s1_coord.sstart - s1_start) >= max_delta ||
					(s2_coord.sstart - s2_start) >= max_delta ||
					(s3_coord.sstart - s3_start) >= max_delta ||
					(s4_coord.sstart - s4_start) >= max_delta
				) {
					let oor = {
						r2D: Math.max(0, (ref2_coord.sstart - ref2_start) - max_delta),
						s1D: Math.max(0, (s1_coord.sstart - s1_start) - max_delta),
						s2D: Math.max(0, (s2_coord.sstart - s2_start) - max_delta),
						s3D: Math.max(0, (s3_coord.sstart - s3_start) - max_delta),
						s4D: Math.max(0, (s4_coord.sstart - s4_start) - max_delta)
					};
					
					console.log($count, "out of range:", "max_delta", max_delta);
					console.table(oor);

					search_end = Math.min(search_end + argv_param.slice_size, ref1_chr_info.length);//inc search range
					max_delta += argv_param.slice_size;
					continue;
				}
				ref2_find_next_start = ref2_coord.send + 1;//force ignore repeat, and any mistake
				s1_find_next_start = s1_coord.send + 1;
				s2_find_next_start = s2_coord.send + 1;
				s3_find_next_start = s3_coord.send + 1;
				s4_find_next_start = s4_coord.send + 1;

				//check qend
				const _qend_list = [ref2_coord.qend, s1_coord.qend, s2_coord.qend, s3_coord.qend, s4_coord.qend];
				let _min_qend = Math.min(..._qend_list);

				if (!(check_q_start_end(ref2_coord, s1_coord) &&
					check_q_start_end(ref2_coord, s2_coord) &&
					check_q_start_end(ref2_coord, s3_coord) &&
					check_q_start_end(ref2_coord, s4_coord)) || (
						ref2_coord.qstart == _min_qend ||
						s1_coord.qstart == _min_qend ||
						s2_coord.qstart == _min_qend ||
						s3_coord.qstart == _min_qend ||
						s4_coord.qstart == _min_qend
					)
				) {
					console.log($count, "out of q start end:", {
						search_start, search_end,
						ref2_coord, s1_coord, s2_coord, s3_coord, s4_coord
					});
					search_end = Math.min(search_end + argv_param.slice_size, ref1_chr_info.length);//inc search range
					max_delta += argv_param.slice_size;
					console.log("max_delta", max_delta);

					debugger;

					continue;
				}

				let max_iterate = 10;
				while ((--max_iterate) > 0) {
					/**
					 * @param {BlastnCoord} row
					 * @param {any[]} params
					 */
					function q_end_filter(row, params) {
						return (
							row.qstart < row.qend &&
							row.sstart < row.send &&
							//row.sstart >= next_start
							row.qend == params[0]
						);
					}

					let _ref2_coords_task = ref2_coord.qend != _min_qend ? blastn_coord(ref1_chr_fasta, ref2_chr_fasta, search_start, _min_qend, ref2_search_start, ref2_coord.send, q_end_filter, [_min_qend], nChr + "_" + $count) : null;
					let _s1_coords_task = s1_coord.qend != _min_qend ? blastn_coord(ref1_chr_fasta, s1_chr_fasta, search_start, _min_qend, s1_search_start, s1_coord.send, q_end_filter, [_min_qend], nChr + "_" + $count) : null;
					let _s2_coords_task = s2_coord.qend != _min_qend ? blastn_coord(ref1_chr_fasta, s2_chr_fasta, search_start, _min_qend, s2_search_start, s2_coord.send, q_end_filter, [_min_qend], nChr + "_" + $count) : null;
					let _s3_coords_task = s3_coord.qend != _min_qend ? blastn_coord(ref1_chr_fasta, s3_chr_fasta, search_start, _min_qend, s3_search_start, s3_coord.send, q_end_filter, [_min_qend], nChr + "_" + $count) : null;
					let _s4_coords_task = s4_coord.qend != _min_qend ? blastn_coord(ref1_chr_fasta, s4_chr_fasta, search_start, _min_qend, s4_search_start, s4_coord.send, q_end_filter, [_min_qend], nChr + "_" + $count) : null;
					let [_ref2_coords, _s1_coords, _s2_coords, _s3_coords, _s4_coords] = await Promise.all([_ref2_coords_task, _s1_coords_task, _s2_coords_task, _s3_coords_task, _s4_coords_task]);
					
					if (
						(_ref2_coords != null && _ref2_coords.length <= 0) ||
						(_s1_coords != null && _s1_coords.length <= 0) ||
						(_s2_coords != null && _s2_coords.length <= 0) ||
						(_s3_coords != null && _s3_coords.length <= 0) ||
						(_s4_coords != null && _s4_coords.length <= 0)
					) {
						console.log({ _min_qend });
						--_min_qend;
						continue;
					}

					let _ref2_coord = _ref2_coords ? _ref2_coords.sort((a, b) => b.send - a.send)[0] : null;
					let _s1_coord = _s1_coords ? _s1_coords.sort((a, b) => b.send - a.send)[0] : null;
					let _s2_coord = _s2_coords ? _s2_coords.sort((a, b) => b.send - a.send)[0] : null;
					let _s3_coord = _s3_coords ? _s3_coords.sort((a, b) => b.send - a.send)[0] : null;
					let _s4_coord = _s4_coords ? _s4_coords.sort((a, b) => b.send - a.send)[0] : null;

					ref2_coord = _ref2_coord || ref2_coord;
					s1_coord = _s1_coord || s1_coord;
					s2_coord = _s2_coord || s2_coord;
					s3_coord = _s3_coord || s3_coord;
					s4_coord = _s4_coord || s4_coord
					
					seq_from_ref1 = [true, ref2_coord != null, s1_coord != null, s2_coord != null, s3_coord != null, s4_coord != null];//last identical loaction

					//check qend
					const qend_list = [ref2_coord.qend, s1_coord.qend, s2_coord.qend, s3_coord.qend, s4_coord.qend];
					const min_qend = Math.min(...qend_list);
					if (!qend_list.every(a => min_qend == a)) {
						console.log($count, "gap or mis");
						continue;
					}

					if (ref2_coord.strand <= 0 || (ref2_start - 1) >= ref2_coord.send ||
						s1_coord.strand <= 0 || (s1_start - 1) >= s1_coord.send ||
						s2_coord.strand <= 0 || (s2_start - 1) >= s2_coord.send ||
						s3_coord.strand <= 0 || (s3_start - 1) >= s3_coord.send ||
						s4_coord.strand <= 0 || (s4_start - 1) >= s4_coord.send
					) {
						console.log($count, "inv:", {
							search_start, search_end,
							ref2_coord, s1_coord, s2_coord, s3_coord, s4_coord,
							ref2_start, s1_start, s2_start, s3_start, s4_start
						});
						continue;
					}

					let r2a = ref2_chr_seq[ref2_coord.send - 1];
					let s1a = s1_chr_seq[s1_coord.send - 1];
					let s2a = s2_chr_seq[s2_coord.send - 1];
					let s3a = s3_chr_seq[s3_coord.send - 1];
					let s4a = s4_chr_seq[s4_coord.send - 1];

					if (r2a == s1a && r2a == s2a && r2a == s3a && r2a == s4a) {
						ref1_start = search_start;
						ref1_end = min_qend;
						break;
					}
					else {
						console.log({ _min_qend, r2a, s1a, s2a, s3a, s4a });
						--_min_qend;
					}
				}
				if (max_iterate <= 0) {
					search_end = Math.min(search_end + argv_param.slice_size, ref1_chr_info.length);
					max_delta += argv_param.slice_size;
					console.log("max_delta", max_delta);
					continue;
				}

				const coord_text_1 = $count + "\tstart\t" + [search_start,	ref1_start,	ref2_start,			s1_start,		s2_start,		s3_start,		s4_start].join("\t|\t");
				const coord_text_2 = $count + "\t  end\t" + [search_end,	ref1_end,	ref2_coord.send,	s1_coord.send,	s2_coord.send,	s3_coord.send,	s4_coord.send].join("\t|\t");

				console.log(coord_text_1);
				console.log(coord_text_2);

				fs.writeFileSync(Path.join(argv_param.output_dir, "multi_coord_ch" + nChr + ".txt"), coord_text_1 + "\n" + coord_text_2 + "\n", { flag: "a" });
				
				let seq1 = ref1_chr_seq.slice(ref1_start - 1, ref1_end);
				let seq2 = ref2_chr_seq.slice(ref2_start - 1, ref2_coord.send);
				let ss1 = s1_chr_seq.slice(s1_start - 1, s1_coord.send);
				let ss2 = s2_chr_seq.slice(s2_start - 1, s2_coord.send);
				let ss3 = s3_chr_seq.slice(s3_start - 1, s3_coord.send);
				let ss4 = s4_chr_seq.slice(s4_start - 1, s4_coord.send);

				let min_len = Math.min(seq1.length, seq2.length, ss1.length, ss2.length, ss3.length, ss4.length);
				let max_len = Math.max(seq1.length, seq2.length, ss1.length, ss2.length, ss3.length, ss4.length);
				fs.writeFileSync(Path.join(argv_param.output_dir, "frag_length_ch" + nChr + ".txt"), [
					$count,
					seq1.length, seq2.length, ss1.length, ss2.length, ss3.length, ss4.length,
					(seq1.length / max_len).toFixed(2),
					(seq2.length / max_len).toFixed(2),
					(ss1.length / max_len).toFixed(2),
					(ss2.length / max_len).toFixed(2),
					(ss3.length / max_len).toFixed(2),
					(ss4.length / max_len).toFixed(2),
					min_len, max_len,
					(min_len / max_len).toFixed(2),
				].join("\t") + "\n", { flag: "a" });

				{
					let rp1 = (seq1.length * 10 / max_len);
					let rp2 = (seq2.length * 10 / max_len);
					let sp1 = (ss1.length * 10 / max_len);
					let sp2 = (ss2.length * 10 / max_len);
					let sp3 = (ss3.length * 10 / max_len);
					let sp4 = (ss4.length * 10 / max_len);

					fs.writeFileSync(Path.join(argv_param.output_dir, "meta", "meta_ch" + nChr + "_" + $count + ".json"), JSON.stringify({
						id: $count,
						match_group: match_results.all_match_groups,
						last_match: {
							//ref1_coord: null,
							ref2_coord: ref2_coord.send,
							s1_coord: s1_coord.send,
							s2_coord: s2_coord.send,
							s3_coord: s3_coord.send,
							s4_coord: s4_coord.send,
						},
						coord: {
							ref1_coord: [ref1_start - 1, ref1_end],
							ref2_coord: [ref2_start - 1, ref2_coord.send],
							s1_coord: [s1_start - 1, s1_coord.send],
							s2_coord: [s2_start - 1, s2_coord.send],
							s3_coord: [s3_start - 1, s3_coord.send],
							s4_coord: [s4_start - 1, s4_coord.send],
						},
						length: [
							seq1.length, seq2.length, ss1.length, ss2.length, ss3.length, ss4.length
						],
						"%": [
							rp1, rp2, sp1, sp2, sp3, sp4
						],
						"=": [
							true,
							rp2 == rp1,
							sp1 == rp1 || sp1 != rp2,
							sp2 == rp1 || sp2 != rp2,
							sp3 == rp1 || sp3 != rp2, 
							sp4 == rp1 || sp4 != rp2
						],
						"==": seqFrom(ref1_coords, ref2_coords, s1_coords, s2_coords, s3_coords, s4_coords),
					}, null, "\t"));
				}
				
				let seq1_fa, seq2_fa, seq3_fa, seq4_fa, seq5_fa, seq6_fa;

				seq1_fa = `>${ref1_chr_name} ${ref1_start}-${ref1_end} ${seq1.length}\n${seq1}\n`;
				seq2_fa = `>${ref2_chr_name} ${ref2_start}-${ref2_coord.send} ${seq2.length}\n${seq2}\n`;
				seq3_fa = `>${s1_chr_name} ${s1_start}-${s1_coord.send} ${ss1.length}\n${ss1}\n`;
				seq4_fa = `>${s2_chr_name} ${s2_start}-${s2_coord.send} ${ss2.length}\n${ss2}\n`;
				seq5_fa = `>${s3_chr_name} ${s3_start}-${s3_coord.send} ${ss3.length}\n${ss3}\n`;
				seq6_fa = `>${s4_chr_name} ${s4_start}-${s4_coord.send} ${ss4.length}\n${ss4}\n`;

				{
					let c1 = seq1_fa.slice(-2);
					let c2 = seq2_fa.slice(-2);
					let c3 = seq3_fa.slice(-2);
					let c4 = seq4_fa.slice(-2);
					let c5 = seq5_fa.slice(-2);
					let c6 = seq6_fa.slice(-2);
					if (c1 != c2 ||
						c1 != c3 ||
						c1 != c4 ||
						c1 != c5 ||
						c1 != c6
					) {
						console
						console.log("diff align len $count=", $count);
						debugger;
					}
				}

				const output_fasta_file_name = Path.join(argv_param.output_frag_directory, `ch${nChr}_${$count}.fa`);
				
				console.log(output_fasta_file_name);

				let text_fasta = [seq1_fa, seq2_fa, seq3_fa, seq4_fa, seq5_fa, seq6_fa].join("");
				fs.writeFileSync(output_fasta_file_name, text_fasta);
				
				//search next range
				if (search_end >= ref1_chr_info.length) {
					break;
				}
				search_start = Math.min(ref1_end + 1, ref1_chr_info.length);
				search_end = Math.min(search_start + argv_param.slice_size, ref1_chr_info.length);//
				if (search_start == search_end) {
					break;
				}
				//
				ref2_start = ref2_coord.send + 1;
				s1_start = s1_coord.send + 1;
				s2_start = s2_coord.send + 1;
				s3_start = s3_coord.send + 1;
				s4_start = s4_coord.send + 1;

				max_delta = argv_param.slice_size * 3;
			}
			else {
				//search next range
				if (search_end >= ref1_chr_info.length) {
					break;
				}

				let _seq_from_ref1 = [true, ref2_coords.length > 0, s1_coords.length > 0, s2_coords.length > 0, s3_coords.length > 0, s4_coords.length > 0];
				if (_seq_from_ref1.every(a => a)) {
					// every has align results, but all query loc no overlap
				}
				else {//save 1:3, 2:2, 3:1
					seq_from_ref1 = _seq_from_ref1;
				}
				search_end = Math.min(search_end + argv_param.slice_size, ref1_chr_info.length);//inc search range
				max_delta += argv_param.slice_size;
				
				console.log("not found");

				console.table({
					search_start, search_end,
					ref2_start, s1_start, s2_start, s3_start, s4_start,
					ref2_find_next_start, s1_find_next_start, s2_find_next_start, s3_find_next_start, s4_find_next_start,
					max_delta,
				});
				console.error("out meta.json");

				console.table({
					//r1: true,
					r2: ref2_coord != null,
					s1: s1_coord != null,
					s2: s2_coord != null,
					s3: s3_coord != null,
					s4: s4_coord != null
				});
			}
		}
		catch (ex) {
			console.error($count, ex.stack);
			throw ex;
		}
	}//for $count
	
	console.log({
		max_delta,
		search_start, search_end,
	});

	{
		let seq_len_list = [];
		let fa_text_1, fa_text_2, fa_text_3, fa_text_4, fa_text_5, fa_text_6;
		{
			ref1_start = ref1_end + 1;

			let seq1 = ref1_chr_seq.slice(ref1_start - 1, ref1_chr_info.length);
			fa_text_1 = `>${ref1_chr_name} ${ref1_start}-${ref1_chr_info.length} ${seq1.length}\n${seq1}\n`;
			seq_len_list.push(seq1.length);
		}
		{
			let seq2 = ref2_chr_seq.slice(ref2_start - 1, ref2_chr_info.length);
			fa_text_2 = `>${ref2_chr_name} ${ref2_start}-${ref2_chr_info.length} ${seq2.length}\n${seq2}\n`;
			seq_len_list.push(seq2.length);
		}
		{
			let ss1 = s1_chr_seq.slice(s1_start - 1, s1_chr_info.length);
			fa_text_3 = `>${s1_chr_name} ${s1_start}-${s1_chr_info.length} ${ss1.length}\n${ss1}\n`;
			seq_len_list.push(ss1.length);
		}
		{
			let ss2 = s2_chr_seq.slice(s2_start - 1, s2_chr_info.length);
			fa_text_4 = `>${s2_chr_name} ${s2_start}-${s2_chr_info.length} ${ss2.length}\n${ss2}\n`;
			seq_len_list.push(ss2.length);
		}
		{
			let ss3 = s3_chr_seq.slice(s3_start - 1, s3_chr_info.length);
			fa_text_5 = `>${s3_chr_name} ${s3_start}-${s3_chr_info.length} ${ss3.length}\n${ss3}\n`;
			seq_len_list.push(ss3.length);
		}
		{
			let ss4 = s4_chr_seq.slice(s4_start - 1, s4_chr_info.length);
			fa_text_6 = `>${s4_chr_name} ${s4_start}-${s4_chr_info.length} ${ss4.length}\n${ss4}\n`;
			seq_len_list.push(ss4.length);
		}

		{
			const coord_text_1 = $count + "\tstart\t" + [search_start,	ref1_start,				ref2_start,				s1_start,			s2_start,			s3_start,			s4_start].join("\t|\t");
			const coord_text_2 = $count + "\t  end\t" + [search_end,	ref1_chr_info.length,	ref2_chr_info.length,	s1_chr_info.length,	s2_chr_info.length,	s3_chr_info.length,	s4_chr_info.length].join("\t|\t");

			console.log(coord_text_1);
			console.log(coord_text_2);
			
			fs.writeFileSync(Path.join(argv_param.output_dir, "multi_coord_ch" + nChr + ".txt"), coord_text_1 + "\n" + coord_text_2 + "\n", { flag: "a" });
		}

		console.log("seq_from_ref1", seq_from_ref1);

		let fa_text_list = [fa_text_1, fa_text_2, fa_text_3, fa_text_4, fa_text_5, fa_text_6];

		console.log({
			_global_search_align: _global_search_align.map(a => a.length),
			_local_overlap_align: Object.keys(_local_overlap_align).map(key => _local_overlap_align[key] != null),
		});

		if (seq_len_list.every(len => len <= argv_param.init_max_delta) ||
			seq_from_ref1.every(a => a) ||
			seq_from_ref1.slice(1).every(a => !a)
		) {//seq_from_ref1[0] -> true
			const output_fasta_file_name = Path.join(argv_param.output_frag_directory, `ch${nChr}_${$count}.fa`);

			console.log(output_fasta_file_name);

			let text_fasta = fa_text_list.join("");
			fs.writeFileSync(output_fasta_file_name, text_fasta);
		}
		else {//check translocation
			let like_ref1_seq = [...fa_text_list.filter((seq, idx) => {
				return seq_from_ref1[idx];
			})].join("");

			let like_ref2_seq = [...fa_text_list.filter((seq, idx) => {
				return !seq_from_ref1[idx];
			})].join("");

			const output_ref1_fasta_file_name = Path.join(argv_param.output_frag_directory, `ch${nChr}_${$count}_ref1.fa`);
			const output_ref2_fasta_file_name = Path.join(argv_param.output_frag_directory, `ch${nChr}_${$count}_ref2.fa`);

			console.log(output_ref1_fasta_file_name);
			console.log(output_ref2_fasta_file_name);

			fs.writeFileSync(output_ref1_fasta_file_name, like_ref1_seq);
			fs.writeFileSync(output_ref2_fasta_file_name, like_ref2_seq);
		}
	}
	console.log("all fin");

	return $count;
}

/**
 * @param {BlastnCoord[]} ref1_list
 * @param {BlastnCoord[]} ref2_list
 * @param {BlastnCoord[]} a_list
 * @param {BlastnCoord[]} b_list
 * @param {BlastnCoord[]} c_list
 * @param {BlastnCoord[]} d_list
 */
function seqFrom(ref1_list, ref2_list, a_list, b_list, c_list, d_list) {
	let s1 = ref1_list.sort((a, b) => b.score - a.score)[0];
	let s2 = ref2_list.sort((a, b) => b.score - a.score)[0];
	let sa = a_list.sort((a, b) => b.score - a.score)[0];
	let sb = b_list.sort((a, b) => b.score - a.score)[0];
	let sc = c_list.sort((a, b) => b.score - a.score)[0];
	let sd = d_list.sort((a, b) => b.score - a.score)[0];
	return [
		true,
		s2.qlen == s2.slen && !s2.gap && !s2.mismatch,
		sa.identity == 100 || sa.identity != s2.identity,
		sb.identity == 100 || sb.identity != s2.identity,
		sc.identity == 100 || sc.identity != s2.identity,
		sd.identity == 100 || sd.identity != s2.identity
	];
}
