//usage: node --max-old-space-size=4096 crossover.js -dataset dataset.json
//node --max-old-space-size=4096 200115_crossover.js -dataset QM6a_WT_pacbio_mafft.json --segfile seg.txt --output-prefix test_ch6

//const NCO_MODE = "all";
const NCO_MODE = "all";
const IGNORE_INDEL_DUPLICATION_POSITION = true;

console.log({
	NCO_MODE,
	IGNORE_INDEL_DUPLICATION_POSITION,
});

const fs = require("fs");


const { argv_parse } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { Dataset } = require("./dataset.js");
const { loadSegFile, SegRow } = require("./SegFile.js");

const argv = argv_parse(process.argv);

const dataset_path = String(argv["-dataset"] || "");
const argv_segfile = String(argv["--segfile"] || "");//mafft_QM6a_WT_1_Ch1_segfile.txt
const argv_output_prefix = String(argv["--output-prefix"] || "");

const dataset = Dataset.loadFromFile(dataset_path);


class CrossoverAnalyser {
	constructor() {
		this.chrMinLength = Infinity;
		this.chrMaxLength = 0;
		
		//meta
		this.chrLength = {};
		
		/** @type {{[nChr:number]:SegRow[]}} */
		this.segFile = {};
		this.featureCounts = {};
		this.has_spo11_oligo_reads = false;
		
		//output
		this.co_list = [];
		this.nco_list = [];
	}
	clear() {
		this.chrMinLength = Infinity;
		this.chrMaxLength = 0;
		
		//meta
		this.chrLength = {};
		
		//data
		this.segFile = {};
		this.featureCounts = {};
	}
	clearOutput() {
		//output
		this.myBook = null;
		this.co_list = [];
		this.nco_list = [];
	}
	saveFile() {
		console.log("write to file");

		if (this.co_list.length || this.nco_list.length) {
			{
				const file_co = argv_output_prefix + "_co.txt";
				const file_nco = argv_output_prefix + "_nco.txt";
				
				//no tetrad

				fs.writeFileSync(file_co, "chr#	chr_len	pos(bp)	co_type	GCasso_tract_len	GCasso_marker#	snp_start_out	snp_start_in	snp_end_in	snp_end_out	type	why_remove	from\r\n");
				fs.writeFileSync(file_nco, "chr#	chr_len	pos(bp)	GCasso_tract_len	GCasso_marker#	snp_start_out	snp_start_in	snp_end_in	snp_end_out	type	why_remove\r\n");
				
				this.co_list.forEach(row => {
					fs.writeFileSync(file_co, row.join("\t") + "\r\n", { flag: "a" });
				});
				this.nco_list.forEach(row => {
					fs.writeFileSync(file_nco, row.join("\t") + "\r\n", { flag: "a" });
				});

				console.log("file_co", file_co);
				console.log("file_nco", file_nco);
			}
			{
				const file_co = argv_output_prefix + "_co.json";
				
				let props = [
					"chr", "chr_len", "pos", "co_type", "GCasso_tract_len", "GCasso_marker", "snp_start_out", "snp_start_in", "snp_end_in", "snp_end_out", "type", "why_remove", "before", "after"
				];

				let list = this.co_list.map(row => props.reduce((obj, key, idx) => {
					obj[key] = row[idx]
					return obj;
				}, {})).filter(co => !co.why_remove);

				fs.writeFileSync(file_co, JSON.stringify(list, null, "\t"));
				
				console.log("file_co", file_co);
			}
			// {
			// 	const file_nco = argv_output_prefix + "_nco.json";
			// 	fs.writeFileSync(file_nco, JSON.stringify(this.nco_list, null, "\t"));
			// 	console.log("file_nco", file_nco);
			// }
		}
	}
	load() {
		{	
			const group = loadSegFile(argv_segfile);

			//group
			//make index
			Object.keys(group).forEach(chr => {
				const nChr = Number(chr);
				const _list = group[nChr].sort((a, b) => a.pos - b.pos);
				
				/** @type {SegRow[]} */
				let list;
				if (IGNORE_INDEL_DUPLICATION_POSITION) {
					list = [];

					/** @type {Map<number, SegRow[]>} */
					let pos_map = new Map();
					_list.forEach(data => {
						let ll = pos_map.get(data.pos);
						if (!ll) {
							ll = [];
							pos_map.set(data.pos, ll);
						}
						ll.push(data);
					});

					pos_map.forEach(ll => {
						if (ll.length == 1) {
							list.push(ll[0]);
						}
						else {
							let nco_list = ll.filter(data => data.type == "type=3" || data.type == "type=4");
							let snv_list = ll.filter(data => !data.type || (!data.isInDel() && data.type == "type=2"));

							if (nco_list.length) {
								list.push(nco_list.pop());//nco
							}
							else if (snv_list.length) {
								list.push(snv_list.pop());//match or 2:2
							}
						}
					});
					list.sort((a, b) => a.pos - b.pos);
				}
				else {
					list = _list;
				}
				
				list.forEach((data, index) => {
					data.$index = index;
				});

				this.segFile[nChr] = list;
			});
		}

		{
			const text_tab = fs.readFileSync(dataset.ref + ".length.txt").toString();
			const header = [
				"nChr", "name", "length",
			];
			const rows = table_to_object_list(tsv_parse(text_tab), header, { start_row: 0 });
			
			this.chrMinLength = Infinity;
			this.chrMaxLength = 0;
			
			const group = {};
			rows.forEach(row => {
				row.nChr = Number(row.nChr);
				
				let list = group[row.nChr];
				if (!list) {
					list = group[row.nChr] = [];
				}
				
				row.length = Number(row.length);
				
				this.chrMinLength = Math.min(this.chrMinLength, row.length);
				this.chrMaxLength = Math.max(this.chrMaxLength, row.length);
				
				list.push(row);
			});
			
			//assign all
			Object.keys(group).forEach(chr => {
				this.chrLength[chr] = group[chr];
			});
		}
	}
	
	find_co_nco_simple() {
		console.log("step1: find all co / nco");

		function is_snp_equal(snp_1, snp_2) {
			return (
				(snp_1[1] == 2 || snp_2[1] == 2 || snp_1[1] == snp_2[1]) &&
				(snp_1[2] == 2 || snp_2[2] == 2 || snp_1[2] == snp_2[2]) &&
				(snp_1[3] == 2 || snp_2[3] == 2 || snp_1[3] == snp_2[3]) &&
				(snp_1[4] == 2 || snp_2[4] == 2 || snp_1[4] == snp_2[4])
			);
		}
		
		let results = Object.keys(this.chrLength).map(nChr => {
			nChr = Number(nChr);

			const seg = this.segFile[nChr];
			
			if (!seg || seg.length <= 0) {
				console.error("chr:", nChr, "if (!seg || seg.length <= 0) {", this.chrLength, Object.keys(this.segFile));
				return;
			}
			
			seg.forEach((data, index) => {
				let nco = [0, 0];//["red", "blue"]
				for (let row = 1; row <= 4; ++row) {
					let is_red = data[row];
					++nco[is_red];
				}
				
				data.$type = `${nco[1]}:${4 - nco[1]}`;
				//data.$snp = (data[4] << 3 | data[3] << 2 | data[2] << 1 | data[1] << 0);

				data.n_type = (nco[0] == 3 || nco[1] == 3) ? 3 : ((nco[0] == 4 || nco[1] == 4) ? 4 : 2);
			});

			function snp_distance(snp_1, snp_2) {
				return Math.abs(snp_2.pos - snp_1.pos);
			}
			
			const closeCOsMinDistance = dataset.crossover.closeCOsMinDistance;

			/** @type {{[start:number]:any}} */
			let snp_block_map = {};
			{
				let prev_snp;
				for (let i = 0; i < seg.length; ++i) {
					let snp_1 = seg[i];
					if (snp_1.n_type == 2) {
						prev_snp = snp_1;
						break;
					}
				}
				for (let i = 0; i < seg.length; ++i) {
					let snp_1 = seg[i];
					if (snp_1.n_type == 2) {
						if (!is_snp_equal(snp_1, prev_snp)) {
							snp_block_map[prev_snp.pos] = {
								start: prev_snp,
								end: snp_1,
								length: snp_1.pos - prev_snp.pos,
							};
							prev_snp = snp_1;
						}
					}
				}
			}
			
			//find simple CO

			let paired_co_list = [];
			{
				for (let i = 0; i < seg.length - 1; ++i) {
					let snp_1 = seg[i];
					let snp_2 = seg[i + 1];
					
					if (snp_1.n_type == 2 && snp_2.n_type == 2) {
						if (is_snp_equal(snp_1, snp_2)) {
						}
						else {
							let co_pair = [snp_1, snp_2];
							paired_co_list.push(co_pair);
						}
					}
				}
				paired_co_list.forEach(([snp_1, snp_2]) => {
					snp_1.$type_name = "CO";
					snp_2.$type_name = "CO";
				});
			}

			//find CO with GC

			let str = seg.map(raw => raw[4] << 3 | raw[3] << 2 | raw[2] << 1 | raw[1] << 0).map(a => a.toString(16)).join("");//make string
			let all_232 = [...str.matchAll(/((?:3|5|6|9|a|c)(?:0|1|2|4|7|8|b|d|e|f)+(?:3|5|6|9|a|c))/g)];//find 2:2 1:3|3:1|0:4|4:0 2:2
			let all_conco = all_232.filter(a => a[1][0] != a[1][a[1].length - 1]).map(a => seg.slice(a.index + 1, a.index + a[1].length - 1));//.map(a => [a.index, a.index + a[1].length - 1]);//co(nco) -> remove 2:2, save GC

			let all_nco = [...str.matchAll(/((?:0|1|2|4|7|8|b|d|e|f)+)/g)].map(a => seg.slice(a.index, a.index + a[1].length));
			
			all_conco.forEach(snp_list => {
				snp_list.forEach(snp => {
					snp.$type_name = "CO(NCO)";
				});
			});
			let paired_conco_list = all_conco;
			let merge_co_list = [
				...paired_conco_list,
				...paired_co_list
			].filter(snps => snps && snps.length).sort((a_snps, b_snps) => a_snps[0].pos - b_snps[0].pos);//sort co and co(nco)
			let remove_snps = {};
			{//remove 2:2 blocks if too short
				merge_co_list.forEach((list, snps_idx) => {
					if (list[0].$type_name == "CO") {
						let head = list[0];
						let foot = list[list.length - 1];
						let start_out = head;
						let end_out = foot;
						
						if (!snp_block_map[end_out.pos]) {
							//console.log("end_out.pos", end_out.pos);
						}
						else if (snp_block_map[end_out.pos].length < closeCOsMinDistance) {
							list.forEach(snp => {
								delete snp.$type_name;
							});
							remove_snps[snps_idx] = list;
							//console.log("remove co 2:2 blocks if too short:", head);
						}
					}
					else {
						let head = list[0];
						let foot = list[list.length - 1];
						let start_out = seg[head.$index - 1];
						let end_out = seg[foot.$index + 1];
						
						if (end_out) {
							if (!snp_block_map[end_out.pos]) {
								//console.log("end_out.pos", end_out.pos);
							}
							else if (snp_block_map[end_out.pos].length < closeCOsMinDistance) {
								list.forEach(snp => {
									delete snp.$type_name;
								});
								remove_snps[snps_idx] = list;
								//console.log("remove co(nco) 2:2 blocks if too short:", head);
							}
						}
					}
				});
			}
			console.log("merge_co_list.length", merge_co_list.length);
			console.log("remove_snps.length", Object.keys(remove_snps).length);
			{//remove co if too close
				let first_snaps = merge_co_list[0];
				for (let snps_idx = 1; snps_idx < merge_co_list.length; ++snps_idx) {
					const current_snps = merge_co_list[snps_idx];
					let foot;

					if (current_snps[0].pos == 1196305) {
						debugger;
					}
					
					if (current_snps[0].$type_name == "CO") {
						foot = current_snps[current_snps.length - 1];
					}
					else {
						foot = seg[current_snps[current_snps.length - 1].$index + 1];
					}

					let first_head;
					
					if (first_snaps[0].$type_name == "CO") {
						first_head = first_snaps[0];
					}
					else {
						first_head = seg[first_snaps[0].$index - 1];
					}

					if (!first_head || !foot) {
						continue;
					}

					let prev_co_dist = snp_distance(first_head, foot);
					if (prev_co_dist < closeCOsMinDistance) {
						const prev_snps_idx = snps_idx - 1;
						const prev_snaps = merge_co_list[prev_snps_idx];

						remove_snps[prev_snps_idx] = prev_snaps;
						prev_snaps.forEach(snp => {
							delete snp.$type_name;
						});
						//console.log("remove L:", prev_snaps[0], [first_head, foot]);

						if (is_snp_equal(first_head, foot)) {
							remove_snps[snps_idx] = current_snps;
							current_snps.forEach(snp => {
								delete snp.$type_name;
							});
							//console.log("remove R:", current_snps[0]);
						}
					}
					else {
						first_snaps = current_snps;
					}
				}
			}
			merge_co_list = merge_co_list.filter(list => list[0].$type_name);

			// check Q, C, S
			if (merge_co_list.length) {
				let co_prev_snp = seg[0];
				for (let index = 0; index < merge_co_list.length - 1; ++index) {
					const snp = merge_co_list[index][0];
					const next = merge_co_list[index + 1][0];
					before_after_crossover(snp, co_prev_snp, next);
					co_prev_snp = snp;
				}
				before_after_crossover(merge_co_list[merge_co_list.length - 1][0], co_prev_snp, seg[seg.length - 1]);
			}
			function before_after_crossover(snp, prev, next) {//QC_before_crossover
				let before = snp_qc(prev, snp);
				let after = snp_qc(snp, next);

				//snp.co_snp = [snp.pos, next.pos, aa[1], aa[2], aa[3], aa[4]].join(",");
				snp.co_x1 = [before[1], before[2], before[3], before[4]].join(",");
				snp.co_x2 = [after[1], after[2], after[3], after[4]].join(",");
			}

			function snp_qc(from_snp, to_snp) {
				let snps = seg.filter(a => a.pos >= from_snp.pos && a.pos <= to_snp.pos);
				
				let aa = {
					1: 0,
					2: 0,
					3: 0,
					4: 0
				};
				snps.forEach(_snp => {
					let [a, b, c, d] = [_snp[1] | 0, _snp[2] | 0, _snp[3] | 0, _snp[4] | 0];
					let n = a + b + c + d;
					if (n == 2) {
						aa[1] += a;
						aa[2] += b;
						aa[3] += c;
						aa[4] += d;
					}
				});

				let hl = snps.length / 2;
				let bb = {};
				Object.keys(aa).forEach(k => bb[k] = aa[k] >= hl ? 1 : 0);

				return bb;
			}

			// find NCO, 2NCO
			
			let paired_nco_list = [];
			let paired_2nco_list = [];
			all_nco.forEach(snp_list => {
				if (snp_list.every(a => a.rip)) {
					return;
				}

				if (snp_list.some(a => a.$type_name == "CO(NCO)")) {
					if (snp_list.every(a => a.$type_name == "CO(NCO)")) {
						//skip CO(NCO)
					}
					else {
						console.error("?? CO(NCO)", snp_list);
					}
				}
				else if (snp_list.every(a => a.n_type == 4)) {
					if (NCO_MODE == "all") {
						const list = snp_list;
						const _n_list = snp_list.filter(a => a.isInDel());
						list.forEach(snp => {
							if (snp.pos == 671437) {
								console.log(snp);
								console.log("_n_list.length", _n_list.length);
							}
							snp.$type_name = "2NCO" + (_n_list.length ? ("=" + _n_list.length + "/" + list.length) : "");
						});
						
						paired_2nco_list.push(list);
					}
					else {
						const snv_list = snp_list.filter(a => !a.isInDel()); //200203 // remove every indel
						if (snv_list.length) {
							snv_list.forEach(snp => {
								snp.$type_name = "2NCO";
							});
							paired_2nco_list.push(snv_list);
						}
					}
				}
				else {
					if (NCO_MODE == "all") {
						const list = snp_list;
						const _n_list = snp_list.filter(a => a.isInDel());
						list.forEach(snp => {
							if (snp.pos == 671437) {
								console.log(snp);
								console.log("_n_list.length", _n_list.length);
							}
							snp.$type_name = "NCO" + (_n_list.length ? ("=" + _n_list.length + "/" + list.length) : "");
						});
						
						paired_nco_list.push(list);
					}
					else {
						const snv_list = snp_list.filter(a => !a.isInDel()); //200203 // remove every indel
						if (snv_list.length) {
							snv_list.forEach(snp => {
								snp.$type_name = "NCO";
							});
							paired_nco_list.push(snv_list);
						}
					}
				}
			});

			let merge_list = [
				...paired_nco_list,
				...paired_2nco_list,
				...merge_co_list
			].filter(snps => snps && snps.length);

			let g_region_id = 0;
			let ret_list = [];
			function make_header(list) {
				++g_region_id;
				let head = list[0];
				let foot = list[list.length - 1];
				if (head.$type_name == "CO") {
					_make_header(head, head, foot, foot);
				}
				else {
					_make_header(seg[head.$index - 1], head, foot, seg[foot.$index + 1]);
				}
			}
			function _make_header(start_out, start_in, end_in, end_out) {
				start_out = start_out || start_in;
				end_out = end_out || end_in;
				
				head_comment.is_comment = true;
				
				head_comment.g_region_id = start_in.region_id;
				head_comment.start = head.pos;
				head_comment.end = foot.pos;
				head_comment.target = head;
				
				head_comment.idx_start_out = idx_start_out;
				head_comment.idx_start_in = idx_start_in;
				head_comment.idx_end_in = idx_end_in;
				head_comment.idx_end_out = idx_end_out;
				
				head_comment["snp_start_out"] = start_out.pos;
				head_comment["snp_start_in"] = head.pos;
				head_comment["snp_end_in"] = foot.pos;
				head_comment["snp_end_out"] = end_out.pos;
				head_comment["pos(bp)"] = Math.round((mid2 + mid1) * 0.5);//191024
				head_comment["tract_len"] = Math.round(Math.abs(mid2 - mid1));//191024
				head_comment["tract_min_len"] = Math.abs(end_in.pos - start_in.pos);
				head_comment["tract_max_len"] = Math.abs(end_out.pos - start_out.pos);
			}
			merge_list.forEach((list, list_idx) => {
				let region_id = list_idx + 1;
				let head = list[0];
				let foot = list[list.length - 1];
				
				let idx_start_out, idx_start_in, idx_end_in, idx_end_out;
				if (head.$type_name == "CO") {
					idx_start_out = head.$index;
					idx_start_in = head.$index;
					idx_end_in = foot.$index;
					idx_end_out = foot.$index;
				}
				else {
					idx_start_out = head.$index - 1;
					idx_start_in = head.$index;
					idx_end_in = foot.$index;
					idx_end_out = foot.$index + 1;
				}
				
				let start_in = seg[idx_start_in];
				let end_in = seg[idx_end_in];
				let start_out = seg[idx_start_out] || start_in;
				let end_out = seg[idx_end_out] || end_in;
				
				let mid1 = (start_in.pos + start_out.pos) * 0.5;
				let mid2 = (end_in.pos + end_out.pos) * 0.5;
				
				list.forEach(a => {
					a.region_id = region_id;
					a.$type_name = head.$type_name;
				});
				
				
				let head_comment;
				if (head.$type_name == "CO(NCO)") {
					head_comment = {
						...start_out,
						$type_name: head.$type_name,
						co_type: 1,
						"GCasso_marker#": list.length,//191021
					};//clone
				}
				else {
					head_comment = {
						...head,
						co_type: 0,
						"GCasso_marker#": head.$type_name == "CO" ? 0 : list.length,//191021
					};//clone
				}

				{
					if (head.$type_name == "CO") {
						//head_comment.co_snp = start_out.co_snp || start_in.co_snp;
						head_comment.co_x1 = start_out.co_x1 || start_in.co_x1;
						head_comment.co_x2 = start_out.co_x2 || start_in.co_x2;
					}
					if (head.$type_name == "CO(NCO)") {
						//head_comment.co_snp = start_out.co_snp || start_in.co_snp;
						head_comment.co_x1 = start_out.co_x1 || start_in.co_x1;
						head_comment.co_x2 = start_out.co_x2 || start_in.co_x2;
					}
				}
				
				head_comment.is_comment = true;
				delete head_comment.is_start;
				delete head_comment.is_end;
				delete head_comment.is_single;
				
				head_comment.region_id = region_id;
				head_comment.start = head.pos;
				head_comment.end = foot.pos;
				head_comment.target = head;
				
				head_comment.idx_start_out = idx_start_out;
				head_comment.idx_start_in = idx_start_in;
				head_comment.idx_end_in = idx_end_in;
				head_comment.idx_end_out = idx_end_out;
				
				head_comment["snp_start_out"] = start_out.pos;
				head_comment["snp_start_in"] = head.pos;
				head_comment["snp_end_in"] = foot.pos;
				head_comment["snp_end_out"] = end_out.pos;
				head_comment["pos(bp)"] = Math.round((mid2 + mid1) * 0.5);//191024
				head_comment["tract_len"] = Math.round(Math.abs(mid2 - mid1));//191024
				head_comment["tract_min_len"] = Math.abs(end_in.pos - start_in.pos);
				head_comment["tract_max_len"] = Math.abs(end_out.pos - start_out.pos);
				
				let remove_region_list = [];
				if (dataset.centromere) {
					remove_region_list.push({
						range: dataset.centromere[head_comment.chr],
						type: "centromere",
					});
				}
				try {
					if (dataset.telomere) {
						let [l_tele, r_tele] = dataset.telomere[head_comment.chr];
						remove_region_list.push({
							range: l_tele,
							type: "telomere",
						}, {
							range: r_tele,
							type: "telomere",
						});
					}
				}
				catch (ex) {
					console.error(ex);
					console.log(dataset.telomere);
					console.log(head_comment.chr);
				}
				if (!remove_region_list.length) {
					console.log(dataset.telomere);
					console.log(dataset.telomere[head_comment.chr]);
					throw 12;
				}
				remove_region_list.forEach(region => {
					let start = head_comment["snp_start_out"];
					let end = head_comment["snp_end_out"];
					let [r_start, r_end] = region.range;
					if (start <= r_end && end >= r_start) {
						head_comment.why_remove = region.type;
						//console.log(head_comment, region.type);
					}
				});
				
				ret_list.push(head_comment);
				
				if (head == foot && head.$type_name != "CO(NCO)" && head.$type_name != "CO") {//only nco 2nco//191024
					delete head.is_start;
					head_comment.is_start = true;
					head.is_single = true;
					foot.is_end = true;
					//head.region_id = region_id;
					
					ret_list.push(head);
				}
				else {
					if (head.$type_name == "CO(NCO)") {
						delete head.is_start;
						delete foot.is_end;
						//
						start_out.$type_name = head.$type_name;
						start_out.region_id = region_id;
						start_out.is_start = true;
						//
						end_out.$type_name = head.$type_name;
						end_out.region_id = region_id;
						end_out.is_end = true;
						//
						ret_list.push(start_out);
						//ret_list.push(head);
						//ret_list.push(foot);
						ret_list.push(end_out);
					}
					else {
						delete head.is_start;
						//delete foot.is_end;
						//
						head_comment.is_start = true;
						//
						//head.region_id = region_id;
						//head.is_start = true;
						//foot.region_id = region_id;
						foot.is_end = true;
						//
						ret_list.push(head);
						ret_list.push(foot);
					}
				}
			});
			
			return ret_list;
		});
		
		results = results.filter(a => a);
		if (results.length <= 0) {
			throw new Error("no result: results.length -> 0");
		}
		
		results.forEach(byChr => {
			byChr.filter(data => data.is_comment).sort((a, b) => a.snp_start_out - b.snp_start_out).forEach(data => {

				if (data.$type_name == "CO" || data.$type_name == "CO(NCO)") {
					let row = [
						//tetrad,//200217
						data.chr,
						this.chrLength[data.chr].length,
						data["pos(bp)"],
						data.co_type,
						data.tract_len,//GCasso_tract_len,
						data["GCasso_marker#"],
						data.snp_start_out, data.snp_start_in, data.snp_end_in, data.snp_end_out,
						data.$type_name,
						data.why_remove,
						//data.co_snp,
						data.co_x1,
						data.co_x2,
					];
					this.co_list.push(row);
				}
				else if (
					(NCO_MODE == "all" && (data.$type_name.indexOf("NCO") >= 0)) ||
					(NCO_MODE != "all" && (data.$type_name == "NCO" || data.$type_name == "2NCO"))
				) {
					let row = [
						//tetrad,
						data.chr,
						this.chrLength[data.chr].length,
						data["pos(bp)"],
						data.tract_len,//GCasso_tract_len,
						data["GCasso_marker#"],
						data.snp_start_out, data.snp_start_in, data.snp_end_in, data.snp_end_out,
						data.$type_name,
						data.why_remove,
					];
					this.nco_list.push(row);
				}
				else {
					//console.log("not co/nco", data);
					throw new Error("ls CO/NCO");
				}
			});
		});
	}
}

let conco = new CrossoverAnalyser();

conco.load();
conco.find_co_nco_simple();
conco.saveFile();

module.exports = CrossoverAnalyser;
