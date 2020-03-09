// @ts-check

const os = require("os");
const fs = require("fs");
const child_process = require("child_process");
const Path = require("path");

const { argv_parse, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { readFasta, saveFasta, parseFasta } = require("./fasta_util.js");
const { BlastnCoord, execAsync, exec_blastn, parseBlastnResults, blastn_coord, isCollide } = require("./blastn_util.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const argv_slice_size = 10000;
const argv_init_max_delta = argv_slice_size * 3;
const argv_minimum_align_length = 1000;

const genome_info_list = dataset.loadGenomeInfoList();

if (!fs.existsSync(`tmp/ma_util_blastn`)) {
	fs.mkdirSync(`tmp/ma_util_blastn`);
}

if (process.argv[1] == __filename) {
	main();
}

async function main() {
	if (!fs.existsSync(`tmp/seq_frag`)) {
		fs.mkdirSync(`tmp/seq_frag`);
	}
	if (!fs.existsSync(`tmp/meta`)) {
		fs.mkdirSync(`tmp/meta`);
	}

	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		console.log("start ch", nChr);
		try {
			const count = await multi_coord(nChr);

			console.log(count);
		}
		catch (ex) {
			fs.writeFileSync(`tmp/ma_util.error.txt`, "ch" + nChr + "\n" + ex.stack + "\n", { flag: "a" });
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
	
	const chr_info_list = genome_info_list.map(genome_info => genome_info.chr_list[chr_idx]);
	const chr_name_list = chr_info_list.map(chr_info => chr_info.chr);
	const chr_fastaPath_list = dataset.genomeNameList.map((gName, i) => `tmp/fasta/${gName}_${chr_name_list[i]}.fa`);

	const chr_seq_list = chr_fastaPath_list.map((filepath, i) => readFasta(filepath)[chr_name_list[i]]);

	const pos_start_list = dataset.genomeNameList.map(a => 1);
	const pos_end_list = dataset.genomeNameList.map(a => -1);
	const find_next_start_list = dataset.genomeNameList.map(a => 0);

	let max_delta = argv_init_max_delta;

	let search_start = 1;
	let search_end = search_start + argv_slice_size;
	
	pos_start_list[0] = search_start;
	pos_end_list[0] = search_end;

	let fragId = 1;

	const seq_from_ref1 = dataset.genomeNameList.map(gName => gName == dataset.ref);

	let _global_search_align, _local_overlap_align;
	let pos_search_start_list = dataset.genomeNameList.map(a => 1);

	for (; search_start <= chr_info_list[0].length && search_end <= chr_info_list[0].length; ++fragId) {
		try {
			/**
			 * @param {BlastnCoord} coord_1
			 * @param {BlastnCoord} coord_2
			 */
			function check_q_start_end(coord_1, coord_2) {
				return coord_1.qstart <= coord_2.qend && coord_1.qend >= coord_2.qstart;
			}
			
			/**
			 * @param {BlastnCoord} r2_coord
			 * @param {BlastnCoord} r1_coord
			 */
			function check_q_s_start_end(r2_coord, r1_coord) {
				return (r2_coord.qstart <= r1_coord.qend && r2_coord.qend >= r1_coord.qstart);
			}
			
			/**
			 * @param {BlastnCoord} r1_coord
			 */
			function check_q_no_repeats(r1_coord) {
				return (r1_coord.qstart == r1_coord.sstart && r1_coord.qend == r1_coord.send);
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

			for (let i = 1; i < pos_search_start_list.length; ++i) {
				pos_search_start_list[i] = Math.max(pos_start_list[i], find_next_start_list[i]);
			}

			//skip ref, pos_search_start[0] => 1
			if (pos_search_start_list.some((pos_search_start, i) => pos_search_start > chr_info_list[i].length)) {
				break;
			}

			let pos_search_end_list = pos_search_start_list.map((pos_search_start, i) => Math.min(pos_search_start + max_delta, chr_info_list[i].length))
			
			// remove ref1-ref1 repeat and duplicate
			let r1_coords_task = blastn_coord(chr_fastaPath_list[0], chr_fastaPath_list[0], search_start, search_end, search_start, search_end, next_start_filter, [], nChr + "_" + fragId);

			let coords_task_list = [];
			for (let i = 1; i < pos_search_start_list.length; ++i) {
				let task = blastn_coord(chr_fastaPath_list[0], chr_fastaPath_list[i], search_start, search_end, pos_search_start_list[i], pos_search_end_list[i], next_start_filter, [], nChr + "_" + fragId);
				coords_task_list.push(task);
			}
			let result_coords_list = await Promise.all([r1_coords_task, ...coords_task_list]);
			_global_search_align = result_coords_list;
			
			/**
			 * @param {BlastnCoord[]} ref_coords - r2[]
			 * @param {BlastnCoord[][]} coords_list - [r1[], r2[], s1[], s2[], s3[], s4[]]
			 */
			function check_overlap_ref(ref_coords, coords_list) {
				let _ref_coords = ref_coords.filter(r2c => r2c.align >= argv_minimum_align_length);
				let _coords_list = coords_list.map(coords => coords.filter(r2c => r2c.align >= argv_minimum_align_length));
				
				if (_coords_list[1] && _coords_list[1].length) {
					let all_match_group = _coords_list[1].map(r2_coord => {
						let _ret = [
							_ref_coords.sort((a, b) => (Math.abs(a.qstart - r2_coord.qstart) + Math.abs(a.qend - r2_coord.qend)) - (Math.abs(b.qstart - r2_coord.qstart) + Math.abs(b.qend - r2_coord.qend)))[0],
							r2_coord,
							..._coords_list.slice(2).map(_coords => _coords.sort((a, b) => (Math.abs(a.qstart - r2_coord.qstart) + Math.abs(a.qend - r2_coord.qend)) - (Math.abs(b.qstart - r2_coord.qstart) + Math.abs(b.qend - r2_coord.qend)))[0]),
						];
						return _ret;
					}).sort((a, b) => b[1].send - a[1].send);

					return {
						/** returns [r1, r2, s1, s2, s3, s4] */
						best_match: all_match_group[0],
						all_match_groups: all_match_group,
					};
				}
				return {
					best_match: [],
					all_match_groups: [],
				};
			}

			let match_results = check_overlap_ref(result_coords_list[1], result_coords_list);
			_local_overlap_align = match_results.best_match;

			if (match_results.best_match.every((match, i) => match)) {
				if (match_results.best_match.some((match, i) => (match.sstart - pos_start_list[i]) >= max_delta)) {

					console.log(fragId, "out of range:", "max_delta", max_delta);
					// let oor = {
					// 	r2D: Math.max(0, (r2_coord.sstart - r2_start) - max_delta),
					// 	s1D: Math.max(0, (s1_coord.sstart - s1_start) - max_delta),
					// 	s2D: Math.max(0, (s2_coord.sstart - s2_start) - max_delta),
					// 	s3D: Math.max(0, (s3_coord.sstart - s3_start) - max_delta),
					// 	s4D: Math.max(0, (s4_coord.sstart - s4_start) - max_delta)
					// };
					// console.table(oor);

					search_end = Math.min(search_end + argv_slice_size, chr_info_list[0].length);//inc search range
					max_delta += argv_slice_size;
					continue;
				}
				for (let i = 1; i < chr_name_list.length; ++i) {
					find_next_start_list[i] = match_results.best_match[i].send + 1;
				}

				//check qend
				let _min_qend = Math.min(...match_results.best_match.map(match => match.qend));

				if (!match_results.best_match.every(match => check_q_start_end(match_results.best_match[1], match)) ||
					match_results.best_match.some(match => match_results.best_match[1].qstart == _min_qend)
				) {
					// console.log(fragId, "out of q start end:", {
					// 	search_start, search_end,
					// 	r2_coord, s1_coord, s2_coord, s3_coord, s4_coord
					// });
					search_end = Math.min(search_end + argv_slice_size, chr_info_list[0].length);//inc search range
					max_delta += argv_slice_size;
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

					const _coords_tasks = match_results.best_match.map((match, i) => {
						if (match.qend != _min_qend) {
							return blastn_coord(chr_fastaPath_list[0], chr_fastaPath_list[i], search_start, _min_qend, pos_search_start_list[i], match.send, q_end_filter, [_min_qend], nChr + "_" + fragId);
						}
						else {
							return null;
						}
					});
					let _coords_list = await Promise.all(_coords_tasks);
					
 					if (_coords_list.some(_coords => _coords != null && _coords.length <= 0)) {
						console.log({ _min_qend });
						--_min_qend;
						continue;
					}

					_coords_list.forEach((_coords, i) => {
						if (_coords) {
							const n = _coords.sort((a, b) => b.send - a.send)[0];
							const p = match_results.best_match[i];
							match_results.best_match[i] = n || p;
						}
					});
					
					match_results.best_match.forEach((coord, i) => {
						seq_from_ref1[i] = coord != null;
					});//last identical loaction

					//check qend
					const qend_list = match_results.best_match.slice(1).map(match => match.qend);
					const min_qend = Math.min(...qend_list);
					if (!qend_list.every(a => min_qend == a)) {
						console.log(fragId, "gap or mis");
						continue;
					}
					
					if (match_results.best_match.every((match, i) => {
						return match.strand <= 0 || (pos_start_list[i] - 1) >= match.send;
					})) {
						console.log(fragId, "inv:", {
							search_start, search_end,
							best_match: match_results.best_match,
							pos_start_list,
						});
						continue;
					}

					const a_list = match_results.best_match.map((coord, i) => chr_seq_list[i][coord.send - 1]);
					if (a_list.slice(2).every(a => a_list[1] == a)) {
						pos_start_list[0] = search_start;
						pos_end_list[0] = min_qend;
						break;
					}
					else {
						console.log({ _min_qend, ...a_list });
						--_min_qend;
					}
				}
				if (max_iterate <= 0) {
					search_end = Math.min(search_end + argv_slice_size, chr_info_list[0].length);
					max_delta += argv_slice_size;
					console.log("max_delta", max_delta);
					continue;
				}
				
				match_results.best_match[0].send = pos_end_list[0];
				const output_coord_list = match_results.best_match.map((match, i) => [pos_start_list[i] - 1, match.send]);
				const extract_seq_list = output_coord_list.map(([start, end], i) => chr_seq_list[i].slice(start, end));
				const extract_length_list = extract_seq_list.map(seq => seq.length);

				{
					const coord_text_1 = fragId + "\tstart\t" + [search_start, ...match_results.best_match.map((match, i) => pos_start_list[i])].join("\t|\t");
					const coord_text_2 = fragId + "\t  end\t" + [search_end,   ...match_results.best_match.map((match, i) => match.send)].join("\t|\t");
					console.log(coord_text_1);
					console.log(coord_text_2);
					fs.writeFileSync(`tmp/multi_coord_ch${nChr}.txt`, coord_text_1 + "\n" + coord_text_2 + "\n", { flag: "a" });
				}

				const min_len = Math.min(...extract_length_list);
				const max_len = Math.max(...extract_length_list);
				fs.writeFileSync(`tmp/frag_length_ch${nChr}.txt`, [
					fragId,
					...extract_length_list,
					...extract_length_list.map(len => (len / max_len).toFixed(2)),
					min_len, max_len,
					(min_len / max_len).toFixed(2),
				].join("\t") + "\n", { flag: "a" });

				fs.writeFileSync(`tmp/meta/meta_ch${nChr}_${fragId}.json`, JSON.stringify({
					id: fragId,
					best_match: match_results.best_match,
					match_group: match_results.all_match_groups,
					coord: output_coord_list,
					length: extract_length_list,
				}, null, "\t"));
				
				const fa_seq_name_list = chr_name_list.map((chrName, i) => {
					const [start, end] = output_coord_list[i];
					return `${chrName} ${start}-${end} ${extract_seq_list[i].length}`;
				});
				
				{
					let cc = extract_seq_list.map(a => a.slice(-2));
					if (cc.some(a => cc[0] != a)) {
						console.log("diff align len fragId=", fragId);
						debugger;
					}
				}


				const output_fasta_file_name = `tmp/seq_frag/ch${nChr}_${fragId}.fa`;
				
				console.log(output_fasta_file_name);

				saveFasta(output_fasta_file_name, fa_seq_name_list.reduce((fasta, seqName, i) => {
					fasta[seqName] = extract_seq_list[i];
					return fasta;
				}, {}));
				
				//search next range
				if (search_end >= chr_info_list[0].length) {
					break;
				}
				search_start = Math.min(pos_end_list[0] + 1, chr_info_list[0].length);
				search_end = Math.min(search_start + argv_slice_size, chr_info_list[0].length);//
				if (search_start == search_end) {
					break;
				}
				//
				match_results.best_match.forEach((match, i) => {
					pos_start_list[i] = match.send + 1;
				});

				max_delta = argv_init_max_delta;
			}
			else {
				//search next range
				if (search_end >= chr_info_list[0].length) {
					break;
				}

				let _seq_from_ref1 = result_coords_list.map(coords => coords.length > 0);
				if (_seq_from_ref1.every(a => a)) {
					// every has align results, but all query loc no overlap
				}
				else {//save 1:3, 2:2, 3:1
					_seq_from_ref1.forEach((value, i) => {
						seq_from_ref1[i] = value;
					});//copy value
				}
				search_end = Math.min(search_end + argv_slice_size, chr_info_list[0].length);//inc search range
				max_delta += argv_slice_size;
				
				console.log("not found");

				console.table({
					search_start, search_end,
					pos_start_list,
					find_next_start_list,
					max_delta,
				});
				console.error("out meta.json");

				console.table(match_results.best_match.map(match => match != null));
			}
		}
		catch (ex) {
			console.error(fragId, ex.stack);
			throw ex;
		}
	}//for fragId
	
	console.log({
		max_delta,
		search_start, search_end,
	});

	{
		pos_start_list[0] = pos_end_list[0] + 1;
		let output_coord_list = pos_start_list.map((start, i) => {
			return [start - 1, chr_info_list[i].length];
		});
		const extract_seq_list = output_coord_list.map(([start, end], i) => chr_seq_list[i].slice(start, end));
		const extract_length_list = extract_seq_list.map(seq => seq.length);
		const fa_seq_name_list = chr_name_list.map((chrName, i) => {
			const [start, end] = output_coord_list[i];
			return `${chrName} ${start}-${end} ${extract_seq_list[i].length}`;
		});

		{
			const coord_text_1 = fragId + "\tstart\t" + [search_start, ...extract_seq_list.map((tuple) => tuple[0])].join("\t|\t");
			const coord_text_2 = fragId + "\t  end\t" + [search_end,   ...extract_seq_list.map((tuple) => tuple[1])].join("\t|\t");

			console.log(coord_text_1);
			console.log(coord_text_2);
			
			fs.writeFileSync(`tmp/multi_coord_ch${nChr}.txt`, coord_text_1 + "\n" + coord_text_2 + "\n", { flag: "a" });
		}

		console.log("seq_from_ref1", seq_from_ref1);

		console.log({
			_global_search_align: _global_search_align.map(a => a.length),
			local_overlap_align: Object.keys(_local_overlap_align).map(key => _local_overlap_align[key] != null),
		});

		if (extract_length_list.every(len => len <= argv_init_max_delta) ||
			seq_from_ref1.every(a => a) ||
			seq_from_ref1.slice(1).every(a => !a)
		) {//seq_from_ref1[0] -> true
			const output_fasta_file_name = `tmp/seq_frag/ch${nChr}_${fragId}.fa`;

			console.log(output_fasta_file_name);
			
			saveFasta(output_fasta_file_name, fa_seq_name_list.reduce((fasta, seqName, i) => {
				fasta[seqName] = extract_seq_list[i];
				return fasta;
			}, {}));
		}
		else {//check translocation
			let like_r1_seq = fa_seq_name_list.filter((seq, idx) => {
				return seq_from_ref1[idx];
			});
			let like_r2_seq = fa_seq_name_list.filter((seq, idx) => {
				return !seq_from_ref1[idx];
			});

			const output_r1_fasta_file_name = `tmp/seq_frag/ch${nChr}_${fragId}_ref1.fa`;
			const output_r2_fasta_file_name = `tmp/seq_frag/ch${nChr}_${fragId}_ref2.fa`;

			console.log(output_r1_fasta_file_name);
			console.log(output_r2_fasta_file_name);
			
			saveFasta(output_r1_fasta_file_name, like_r1_seq.reduce((fasta, seqName, i) => {
				fasta[seqName] = extract_seq_list[i];
				return fasta;
			}, {}));
			saveFasta(output_r2_fasta_file_name, like_r2_seq.reduce((fasta, seqName, i) => {
				fasta[seqName] = extract_seq_list[i];
				return fasta;
			}, {}));
		}
	}
	console.log("all fin");

	return fragId;
}

/**
 * @param {BlastnCoord[]} r1_list
 * @param {BlastnCoord[]} r2_list
 * @param {BlastnCoord[]} a_list
 * @param {BlastnCoord[]} b_list
 * @param {BlastnCoord[]} c_list
 * @param {BlastnCoord[]} d_list
 */
function seqFrom(r1_list, r2_list, a_list, b_list, c_list, d_list) {
	let s1 = r1_list.sort((a, b) => b.score - a.score)[0];
	let s2 = r2_list.sort((a, b) => b.score - a.score)[0];
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
