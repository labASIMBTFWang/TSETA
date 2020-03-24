// @ts-check

class RIP_data {
	constructor() {
		/** @type {number} pos index */
		this.pos = 0;
		
		this.ref1_pos = 0;
		
		/** @type {number} 0: no rip, 1: ref1, 2: ref2 */
		this.rip_ref = 0;

		this.ref1 = "";
		this.ref2 = "";
		this.a = "";
		this.b = "";
		this.c = "";
		this.d = "";
	}
}

let align_start_index;
let align_end_index;

let seg_snp = [];
/** @type {RIP_data[]} */
let rip_list = [];

///** @type {Uint8Array} */
//let parental_cmp_uint8array = null;
let region_rect = [];
/** @type {Uint32Array} ref1 pos map to multialign */
let ref1_pos_uint32array = null;
/** @type {Uint32Array} multialign pos map to ref1 */
let pos_ref1_uint32array = null;
/** @type {Uint32Array} ref2 pos map to multialign */
let ref2_pos_uint32array = null;
/** @type {Uint32Array} multialign pos map to ref2 */
let pos_ref2_uint32array = null;
// /** @type {Uint8Array} */
// let rip_map_uint8array = null;
/** @type {Uint32Array} ref1 ref2 score */
let ref1_ref2_score_uint32array = null;

class MarkerValue {
	constructor() {
		/** @type {number} */
		this.pos = 0;
		/** @type {number} */
		this.value = 0;
		/** @type {number} */
		this.ref1_pos = 0;
		/** @type {string} */
		this.type = "";
		/** @type {string} */
		this.name = "";
	}
}

class MarkerData {
	/**
	 * @param {string} name - display name
	 * @param {string} property - property name
	 * @param {number} order - display order
	 */
	constructor(name, property, order) {
		/** @type {string} - display name */
		this.name = name;

		/** @type {string} - property name */
		this.property = property;

		/** @type {number} - display order */
		this.order = order;

		/** @type {MarkerValue[]} */
		this.values = [];
	}
}

class MarkerList {
	constructor() {
		/** @type {MarkerData[]} */
		this.list = [];

		/** @type {{[propertyName:string]:MarkerData}} */
		this.map = {};
	}

	/**
	 * @param {string} name - display name
	 * @param {string} property - property name
	 * @param {any} [_]
	 */
	add_marker_type(name, property, _) {
		const marker_data = new MarkerData(name, property, this.list.length)
		this.list.push(marker_data);
		this.map[property] = marker_data;
	}

	/**
	 * @param {string} markerName
	 * @param {MarkerValue} value
	 */
	add_marker_value(markerName, value) {
		this.map[markerName].values.push(value);
	}
}

/** @type {MarkerList} */
let allMarker = null;

const ColorID = get_colorID_Data();
function get_colorID_Data() {
	let colorID = {
	};
	colorID.mask = 0b1111;
	colorID.dad = 0;
	colorID.mom = 1;
	colorID.identical = 2;
	colorID.dad_rip = 3;
	colorID.mom_rip = 4;
	colorID.diff = 5;
	colorID.none = 8;
	colorID.indel_bit = 0b1000_0000;
	colorID.indel_mask = 0b1000_0000;
	return colorID;
}
class AnalysisOptions {
	constructor() {
		/** @type {"tetrad"|"SNP"} */
		this.mode = "tetrad";

		/** @param {number} nChr */
		this.nChr = 0;

		/** @type {rDNA_info} */
		this.rDNA_info = null;

		/** @type {crossover_list[0]} - this chr co list*/
		this.co_list = null;

		this.fill_prev_color = true;

		/** @type {function(string):void} */
		this._onProgressReport = null;

		/** @type {boolean} */
		this.show_rDNA_snp = false;
	}
}

/**
 * @param {Partial<AnalysisOptions>} options
 */
function initAnalyser(options) {
	allMarker = new MarkerList();

	if (options.mode == "tetrad") {
		allMarker.add_marker_type("(2:2 makers)", "22", "MARKER_22");
		allMarker.add_marker_type("Noncrossors (3:1 makers)", "31", "MARKER_31");
		allMarker.add_marker_type("Noncrossors (4:0 makers)",  "40","MARKER_40");
		allMarker.add_marker_type("1n:3 makers", "1n3", "MARKER_1n3");
		allMarker.add_marker_type("2n:2 makers", "2n2", "MARKER_2n2");
		allMarker.add_marker_type("3n:1 makers", "3n1", "MARKER_3n1");
		allMarker.add_marker_type("4n:0 makers", "4n0", "MARKER_4n0");
		allMarker.add_marker_type("RIP", "rip", "MARKER_RIP");
		allMarker.add_marker_type("illegitimate mutation makers", "illegitimate_mutation", "MARKER_31");
	}
	else {
		allMarker.add_marker_type("(dummy makers)", "dummy", "MARKER_dummy");
		allMarker.add_marker_type("snp", "snp");
		allMarker.add_marker_type("snv", "snv");
	}
}

/**
 * @param {seq_list} seq_list
 * @param {Partial<AnalysisOptions>} options
 */
function loadCmpData(seq_list, options) {
	calc_seg_reg(seq_list, options);
	return {
		// argv_fill_prev_color: options.fill_prev_color,
		// argv_nChr: options.nChr,

		align_start_index,
		align_end_index,

		rip_list,
		seg_snp,

		region_rect,
		
		/** @type {Uint32Array} ref1 pos map to multialign */
		ref1_pos_uint32array: ref1_pos_uint32array,
		/** @type {Uint32Array} multialign pos map to ref1 */
		pos_ref1_uint32array: pos_ref1_uint32array,
		/** @type {Uint32Array} ref2 pos map to multialign */
		ref2_pos_uint32array: ref2_pos_uint32array,
		/** @type {Uint32Array} multialign pos map to ref2 */
		pos_ref2_uint32array: pos_ref2_uint32array,
		
		spore_cmp: {
			s22: allMarker.map["22"],
			s31: allMarker.map["31"],
			s40: allMarker.map["40"],
			s1n3: allMarker.map["1n3"],
			s2n2: allMarker.map["2n2"],
			s3n1: allMarker.map["3n1"],
			s4n0: allMarker.map["4n0"],
			illegitimate_mutation_list: allMarker.map["illegitimate_mutation"],
		},
		
		ColorID: ColorID,
	}
}

/**
 * @param {seq_list} seq_list
 * @param {Partial<AnalysisOptions>} options
 */
function calc_seg_reg(seq_list, options) {
	if (options.mode == null) {
		throw new TypeError("options.mode");
	}

	init_viewModel(seq_list);

	/** @type {{ start: number, end: number, state: number[], co_idx: number }[]} */
	let co_detail_list = options.co_list ? [] : null;
	if (options.co_list) {
		let prev_snp_r1_pos = 0;
		for (let i = 0; i < options.co_list.length; ++i) {
			let co = options.co_list[i];
			
			co_detail_list.push({
				start: prev_snp_r1_pos,
				end: co.snp_start_in,
				state: co.before,
				co_idx: i,
			});

			prev_snp_r1_pos = co.snp_end_in;
		}
		if (options.co_list.length > 1) {
			let last_co = options.co_list[options.co_list.length - 1];
			co_detail_list.push({
				start: prev_snp_r1_pos,
				end: Infinity,//pos_ref1_uint32array[pos_ref1_uint32array.length - 1],//ref1 last pos // pos_ref1_uint32array was not define yet
				state: last_co.after,
				co_idx: options.co_list.length,
			});
		}
	}
	console.log({
		co_detail_list,
	});

	/** @type {co_detail_list} */
	let work_co_detail_queue = co_detail_list ? [...co_detail_list] : null;
	/** @type {co_detail_list[0]} */
	let current_co_detail = co_detail_list ? work_co_detail_queue.shift() : null;

	seg_snp = new Array(seq_list[0].length);
	
	// check flanking => remove flanking => telomere
	for (let i = 0; i < seq_list[0].length; ++i) {
		if (seq_list.every(aa => aa[i] != "-")) {
			align_start_index = i;
			break;
		}
	}
	for (let i = seq_list[0].length - 1; i >= 0 ; --i) {
		if (seq_list.every(aa => aa[i] != "-")) {
			align_end_index = i;
			break;
		}
	}

	console.log({
		align_start_index,
		align_end_index,
	});
	if (align_start_index == -1 || align_end_index == -1 || align_start_index == align_end_index) {
		throw new Error("if (align_start_index == -1 || align_end_index == -1 || align_start_index == align_end_index) {");
	}

	// begin make_seg

	let ref1_pos = 1;
	let ref2_pos = 1;

	//args_colors
	let prev_22_color = [
		0, 1,
		0, 0, 0, 0
	];
	/** @type {[number, number, number, number, number, number]} */
	let all_prev_color = [
		0, 1,
		0, 0, 0, 0
	];
	let prev_has_rip = 0;

	//compare seq
	let ref1_seq = seq_list[0];
	let ref2_seq = seq_list[1];
	for (let pos = 0; pos < ref1_seq.length; ++pos) {
		// if (pos == 1645832) {
		// 	debugger;
		// }
		let ref1 = ref1_seq[pos];
		let ref2 = ref2_seq[pos];
		
		//AnalysisOptions
		if (options.mode == "tetrad") {
			push_tetrad_data(
				Math.max(1, ref1 != "-" ? ref1_pos : (ref1_pos - 1)),
				Math.max(1, ref2 != "-" ? ref2_pos : (ref2_pos - 1))
			);
		}
		else {
			push_SNP_data(Math.max(1, ref1 != "-" ? ref1_pos : (ref1_pos - 1)));
		}

		/**
		 * @param {number} ref1_pos
		 */
		function push_SNP_data(ref1_pos) {
			let row = {
				chr: options.nChr,
				pos: pos,
			};

			let has_snp = false;
			let has_indel = false;

			seq_list.slice(1).forEach((seq, idx) => {
				const sid = Number(idx) + 1;
				const va = seq[pos];
				
				if (va == ref1) {
					row[sid] = ColorID.dad;
				}
				else {
					row[sid] = ColorID.mom;
					has_snp = true;
				}

				let is_indel = va == "-";
				has_indel = is_indel || has_indel;

				row[sid] = row[sid] | (is_indel ? ColorID.indel_bit : 0);
			});
			let ref1_is_indel = ref1 == "-";
			row[0] = ColorID.dad | (ref1_is_indel ? ColorID.indel_bit : 0);

			seg_snp[pos] = row;
			
			if (has_snp) {
				if (ref1_is_indel || has_indel) {
					let marker = {
						type: "SNP",
						name: allMarker.map["snp"].name,
						pos: pos,
						value: 0,
						ref1_pos: ref1_pos,
					};
					allMarker.map["snp"].values.push(marker);
				}
				else {
					let marker = {
						type: "SNV",
						name: allMarker.map["snv"].name,
						pos: pos,
						value: 1,
						ref1_pos: ref1_pos,
					};
					allMarker.map["snv"].values.push(marker);
				}
			}

			if (!ref1_pos_uint32array[ref1_pos]) {
				ref1_pos_uint32array[ref1_pos] = pos + 1;
			}
			pos_ref1_uint32array[pos] = ref1_pos;
		}

		/**
		 * @param {number} ref1_pos
		 * @param {number} ref2_pos
		 */
		function push_tetrad_data(ref1_pos, ref2_pos) {
			if (co_detail_list) {
				if (ref1_pos < current_co_detail.end) {
					if (ref1_pos > current_co_detail.start) {
					}
					else {
						// console.log({
						// 	what: "CO(NCO) inner ??",
						// 	ref1_pos,
						// 	current_co_detail,
						// });
						// debugger;
					}
				}
				else {
					current_co_detail = work_co_detail_queue.shift();
				}
			}
			
			let a = seq_list[2][pos];
			let b = seq_list[3][pos];
			let c = seq_list[4][pos];
			let d = seq_list[5][pos];
			/** @type {[string, string, string, string]} */
			let spores = [a, b, c, d];

			const cmp_ref1_ref2 = ref1 == ref2;

			if (ref1 == "-" || ref2 == "-") {
				let prev_score = ref1_ref2_score_uint32array[pos - 1] | 0;
				ref1_ref2_score_uint32array[pos] = Math.max(0, prev_score - 1);
			}
			else {
				let prev_score = ref1_ref2_score_uint32array[pos - 1] | 0;
				ref1_ref2_score_uint32array[pos] = Math.max(0, prev_score + (cmp_ref1_ref2 ? 1 : -2));
			}
			
			let row = {
				chr: options.nChr,
				pos: pos,
			};

			let is_rDNA = (() => {
				if (options.show_rDNA_snp) {//ignore_rDNA
					return false;
				}
				else if (options.rDNA_info && options.nChr == options.rDNA_info.chr) {
					if (pos >= options.rDNA_info.alignment_start && pos <= options.rDNA_info.alignment_end) {
						return true;
					}
				}
			})();
			if (is_rDNA) {
				row[0] = 0 | (ref1 == "-" ? ColorID.indel_bit : 0);
				row[1] = 1 | (ref2 == "-" ? ColorID.indel_bit : 0);
				row[2] = 1 | (a == "-" ? ColorID.indel_bit : 0);
				row[3] = 1 | (b == "-" ? ColorID.indel_bit : 0);
				row[4] = 0 | (c == "-" ? ColorID.indel_bit : 0);
				row[5] = 0 | (d == "-" ? ColorID.indel_bit : 0);
			}
			else {
				push_not_rDNA();
			}//if (!is_rDNA) {

			seg_snp[pos] = row;

			if (!ref1_pos_uint32array[ref1_pos]) {
				ref1_pos_uint32array[ref1_pos] = pos + 1;
			}
			if (!ref2_pos_uint32array[ref2_pos]) {
				ref2_pos_uint32array[ref2_pos] = pos + 1;
			}
			pos_ref1_uint32array[pos] = ref1_pos;
			pos_ref2_uint32array[pos] = ref2_pos;

			function push_not_rDNA() {
				// if (pos == 2019326) {
				// 	debugger;
				// }
				let rip_ref = co_detail_list && isRIP(current_co_detail, ref1_pos, a, b, c, d, spores, ref1, ref2, cmp_ref1_ref2, row, all_prev_color, ColorID);
				if (rip_ref) {
					push_seg_RIP();//4:0 2:2 RIP
					
					let a_indel = spores.filter(a => a == "-");
					let spores_is_indel = a_indel.length > 0;
					let col_idx = allMarker.map["rip"].order;
					let marker = Object.assign({
						type: "RIP",
						name: allMarker.list[col_idx].name,
					}, make_marker(spores_is_indel));
					allMarker.list[col_idx].values.push(marker);
				}
				else {
					if (!cmp_ref1_ref2) {
						push_seg_SNP();
					}
					else {
						push_seg_not_SNP();
					}
					
					const label = push_spore_maker();
					if (label) {
						let { col_idx, marker } = label;
						if (!allMarker.list[col_idx]) {
							console.log({col_idx});
						}
						allMarker.list[col_idx].values.push(marker);
					}
				}
				prev_has_rip = rip_ref;

				function push_spore_maker() {
					let a_indel = spores.filter(a => a == "-");
					let spores_is_indel = a_indel.length > 0;
					let is_snp = ref1 != ref2;
					let is_snv = !spores_is_indel && is_snp;
					let a_1s = spores.filter(a => a == ref1);
					let a_2s = spores.filter(a => a == ref2);
					let n_1s = a_1s.length;
					let n_2s = a_2s.length;

					const snp_marker_order_map = [ [], [], [], [], [] ];

					snp_marker_order_map[0][4] = allMarker.map["40"].order;
					snp_marker_order_map[1][3] = allMarker.map["31"].order;
					snp_marker_order_map[2][2] = allMarker.map["22"].order;
					snp_marker_order_map[3][1] = allMarker.map["31"].order;
					snp_marker_order_map[4][0] = allMarker.map["40"].order;

					/** @type {number} */
					let col_idx = snp_marker_order_map[n_1s][n_2s];

					let marker;

					if (col_idx != null) {
						if (is_snv) {
							// 2:2, 3:1, 4:0
							marker = Object.assign({
								type: "SNV",
								name: allMarker.list[col_idx].name,
							}, make_marker(spores_is_indel));
						}
						else {
							// 1n3 2n2 3n1 4n0
							if (col_idx == allMarker.map["31"].order) {
								if (a_indel.length == 1) {
									col_idx = allMarker.map["1n3"].order;
								}
								else if (a_indel.length == 3) {
									col_idx = allMarker.map["3n1"].order;
								}
							}
							else {
								if (col_idx == allMarker.map["40"].order) {
									col_idx = allMarker.map["4n0"].order;
								}
								else if (col_idx == allMarker.map["22"].order) {
									col_idx = allMarker.map["2n2"].order;
								}
							}
							marker = Object.assign({
								type: "SNP",
								name: allMarker.list[col_idx].name,
							}, make_marker(spores_is_indel));
						}
					}
					else {
						let a_n1s = spores.filter(a => a != ref1);
						let a_n2s = spores.filter(a => a != ref2);

						let mut_n1s = a_n1s.filter(a => a != ref2 && ref2 != "-" && a != "-");	//mutation
						let mut_n2s = a_n2s.filter(a => a != ref1 && ref1 != "-" && a != "-");	//mutation
						let mut_arr = [
							mut_n1s, mut_n2s,
						];

						// let n1s_ins = a_n1s.filter(a => a != ref2 && ref2 == "-" && a != "-");	//progeny ins != r1, r2 del
						// let n2s_ins = a_n2s.filter(a => a != ref1 && ref1 == "-" && a != "-");	//progeny ins != r2, r1 del
						// let p_ins_arr = [
						// 	n1s_ins, n2s_ins,
						// ];

						if (mut_arr.some(aa => aa.length > 0)) {
							//is rip
							col_idx = allMarker.map["illegitimate_mutation"].order;

							// illegitimate mutation
							marker = Object.assign({
								type: "illegitimate_mutation",
								name: "illegitimate_mutation",
							}, make_marker(spores_is_indel));
						}
						// else if (p_ins_arr.some(aa => aa.length > 0)) {
						// 	//progeny ins 1n3 2n2 3n1 
						// 	col_idx = markerTypes.map["PROGENY_INS"];
						// 	// progeny ins
						// 	marker = Object.assign({
						// 		type: "progeny ins",
						// 		name: "progeny ins",
						// 	}, make_marker(is_indel));
						// }
						else {
							// progeny indel
							// 1n3 2n2 3n1 4n0
							let n_del = a_indel.length;
							if (n_del == 1) {
								col_idx = allMarker.map["1n3"].order;
							}
							else if (n_del == 2) {
								col_idx = allMarker.map["2n2"].order;
							}
							else if (n_del == 3) {
								col_idx = allMarker.map["3n1"].order;
							}
							else if (n_del == 4) {
								col_idx = allMarker.map["4n0"].order;
							}
							if (col_idx != null) {
								marker = Object.assign({
									type: "progeny indel",
									name: allMarker.list[col_idx]. name,
								}, make_marker(spores_is_indel));
							}
						}
					}

					if (col_idx >= 0 && marker) {
						return {
							col_idx, marker,
						};
					}

					// if ([0, 1, 2, 3, 4, 5].some(i => row[i] == ColorID.diff)) {
					// 	illegitimate_mutation.push({
					// 		pos: pos,//index
					// 		ref1_pos,
					// 	});
					// }
				}

				/**
				 * @param {boolean} is_indel
				 */
				function make_marker(is_indel) {
					return {
						pos: pos,
						ref1_pos,
						value: (row[5] << 3 | row[4] << 2 | row[3] << 1 | row[2] << 0) | (is_indel ? ColorID.indel_bit : 0),
					};
				}

				function push_seg_not_SNP() {
					if (options.fill_prev_color) {
						if (ref1 == "-" && ref2 == "-") {
							row[0] = ColorID.none;
							row[1] = ColorID.none;
							row[2] = a != "-" ? ColorID.diff : ColorID.none;
							row[3] = b != "-" ? ColorID.diff : ColorID.none;
							row[4] = c != "-" ? ColorID.diff : ColorID.none;
							row[5] = d != "-" ? ColorID.diff : ColorID.none;
							all_prev_color[0] = ColorID.none;
							all_prev_color[1] = ColorID.none;
							all_prev_color[2] = ColorID.none;
							all_prev_color[3] = ColorID.none;
							all_prev_color[4] = ColorID.none;
							all_prev_color[5] = ColorID.none;
						}
						else {
							// if (pos == 2328595) {
							// 	debugger;
							// }
							row[0] = ColorID.identical;
							row[1] = ColorID.identical;
							row[2] = prev_has_rip ? prev_22_color[2] : (a == ref1 ? ColorID.identical : ColorID.diff); //current_colorset["identical"]
							row[3] = prev_has_rip ? prev_22_color[3] : (b == ref1 ? ColorID.identical : ColorID.diff);
							row[4] = prev_has_rip ? prev_22_color[4] : (c == ref1 ? ColorID.identical : ColorID.diff);
							row[5] = prev_has_rip ? prev_22_color[5] : (d == ref1 ? ColorID.identical : ColorID.diff);
							all_prev_color[0] = ColorID.identical;
							all_prev_color[1] = ColorID.identical;
							all_prev_color[2] = ColorID.identical;
							all_prev_color[3] = ColorID.identical;
							all_prev_color[4] = ColorID.identical;
							all_prev_color[5] = ColorID.identical;
						}
					}
					else {
						if (ref1 == "-" && ref2 == "-") {
							row[0] = ColorID.none;
							row[1] = ColorID.none;
							row[2] = a != "-" ? ColorID.diff : ColorID.none;
							row[3] = b != "-" ? ColorID.diff : ColorID.none;
							row[4] = c != "-" ? ColorID.diff : ColorID.none;
							row[5] = d != "-" ? ColorID.diff : ColorID.none;
						}
						else {
							row[0] = ColorID.identical;
							row[1] = ColorID.identical;
							row[2] = ColorID.identical; //current_colorset["identical"]
							row[3] = ColorID.identical;
							row[4] = ColorID.identical;
							row[5] = ColorID.identical;
						}
					}
				}

				function push_seg_RIP() {
					let tmp = [];
					tmp[2] = (a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff)) | (a == "-" ? ColorID.indel_bit : 0);
					tmp[3] = (b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff)) | (b == "-" ? ColorID.indel_bit : 0);
					tmp[4] = (c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff)) | (c == "-" ? ColorID.indel_bit : 0);
					tmp[5] = (d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff)) | (d == "-" ? ColorID.indel_bit : 0);
					for (let i = 2; i <= 5; ++i) {
						if (tmp[i] == prev_22_color[i]) {
							row[i] = prev_22_color[i];
						}
					}
					row.is_rip = rip_ref;
					rip_list.push({
						pos,//index
						ref1_pos,
						rip_ref: rip_ref,
						ref1,
						ref2,
						a,
						b,
						c,
						d,
					});
				}

				function push_seg_SNP() {
					row[0] = ColorID.dad | (ref1 == "-" ? ColorID.indel_bit : 0);
					row[1] = ColorID.mom | (ref2 == "-" ? ColorID.indel_bit : 0);
					row[2] = (a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff)) | (a == "-" ? ColorID.indel_bit : 0);
					row[3] = (b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff)) | (b == "-" ? ColorID.indel_bit : 0);
					row[4] = (c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff)) | (c == "-" ? ColorID.indel_bit : 0);
					row[5] = (d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff)) | (d == "-" ? ColorID.indel_bit : 0);
					prev_22_color[0] = row[0];
					prev_22_color[1] = row[1];
					prev_22_color[2] = row[2];
					prev_22_color[3] = row[3];
					prev_22_color[4] = row[4];
					prev_22_color[5] = row[5];
					all_prev_color[0] = row[0];
					all_prev_color[1] = row[1];
					all_prev_color[2] = row[2];
					all_prev_color[3] = row[3];
					all_prev_color[4] = row[4];
					all_prev_color[5] = row[5];
					//seg.push(row);

					//2:2, 3:1, 4:0
					//push_SNP_maker();//snp 1n:3 2n:2 3n:1 4n:0 maker

					function push_SNP_maker() {
						if ((row[2] & ColorID.mask) <= ColorID.mom &&
							(row[3] & ColorID.mask) <= ColorID.mom &&
							(row[4] & ColorID.mask) <= ColorID.mom &&
							(row[5] & ColorID.mask) <= ColorID.mom) {
							//let va = (row[2] & ColorID.mask) + (row[3] & ColorID.mask) + (row[4] & ColorID.mask) + (row[5] & ColorID.mask);
							let indel = ref1 == "-" || ref2 == "-" || a == "-" || b == "-" || c == "-" || d == "-";
							//if (indel || va != 2) {
							//const cmp_column_idx = [4, 2, 0, 2, 4];//2:2,3:1,4:0
							let sq = spores.filter(v => v == ref1);
							let sc = spores.filter(v => v == ref2);
							let nsq = sq.length;
							let nsc = sc.length;
							let col_idx = -1;
							if (ref1 != "-" && ref2 != "-") {
								if (nsq == 2 && nsc == 2) {
									//2:2
									col_idx = 0;
								}
								else if (nsq == 1 || nsc == 1) {
									//3:1
									col_idx = 1;
								}
								else if (nsq == 4 || nsc == 4) {
									//4:1
									col_idx = 2;
								}
							}
							else if ((nsq == 1 && ref1 == "-") || nsc == 1 && ref2 == "-") {
								//1n:3
								col_idx = 3;
							}
							else if ((nsq == 2 || nsc == 2) && indel) {
								//2n:2
								col_idx = 4;
							}
							else if ((nsq == 3 && ref1 == "-") || nsc == 3 && ref2 == "-") {
								//3n:1
								col_idx = 5;
							}
							else if ((nsq == 4 && ref1 == "-") || nsc == 4 && ref2 == "-") {
								//4n:1
								col_idx = 6;
							}
							// spore_cmp_array[cmp_column_idx[va] + (indel ? 1 : 0)].push({
							// 	pos: pos,
							// 	value: (row[5] << 3 | row[4] << 2 | row[3] << 1 | row[2] << 0) | (indel ? ColorID.indel_bit : 0),
							// });
							if (col_idx >= 0 /* && pos >= align_start_index && pos <= align_end_index*/) {
								allMarker.list[col_idx].values.push({
									pos: pos,
									ref1_pos,
									value: (row[5] << 3 | row[4] << 2 | row[3] << 1 | row[2] << 0) | (indel ? ColorID.indel_bit : 0),
									type: "nco",
									name: "" + col_idx,
								});
							}
							//}
						}//is maker
					}//push_maker
				}//push_SNP
			}//notRDNA()
		}

		if (ref1 != "-") {
			++ref1_pos;
		}
		if (ref2 != "-") {
			++ref2_pos;
		}
	}

	// begin merge_seg

	_merge_seg(seq_list, seg_snp);
}

/**
 * @param {{ start: number, end: number, state: number[], co_idx: number }} co_detail
 * @param {number} ref1_pos
 * @param {string} a s1
 * @param {string} b s2
 * @param {string} c s3
 * @param {string} d s4
 * @param {[string, string, string, string]} spores [a, b, c, d]
 * @param {string} ref1 ref1
 * @param {string} ref2 ref2
 * @param {boolean} cmp_ref1_ref2 ref1 == ref2
 * @param {*} row seg row
 * @param {[number, number, number, number, number, number]} all_prev_color 
 * @param {{[id:string]:number}} ColorID
 * @returns {number} 1: ref1, 2: ref2
 */
function isRIP(co_detail, ref1_pos, a, b, c, d, spores, ref1, ref2, cmp_ref1_ref2, row, all_prev_color, ColorID) {
	if (co_detail.state) {
		const co_state = co_detail.state;

		let seg_ref1 = 0;
		let seg_ref2 = 0;
		let rip_ref1 = 0;
		let rip_ref2 = 0;
		
		let spores_rip = co_state.map((spore_val, spore_idx) => {
			let spore_s = spores[spore_idx];
			if (spore_val == 0) {
				if (ref1 == spore_s) {
					++seg_ref1;

					return {
						ref: -1,
						type: null
					};
				}
				else if (ref1 == "G" && spore_s == "A") {
					++rip_ref1;
					return {
						ref: 1,
						type: "A"
					};
				}
				else if (ref1 == "C" && spore_s == "T") {
					++rip_ref1;
					return {
						ref: 1,
						type: "T"
					};
				}
				return {
					//illeg
					ref: -1,
					type: null
				};
			}
			else {
				if (ref2 == spore_s) {
					++seg_ref2;
					return {
						ref: -2,
						type: null
					};
				}
				else if (ref2 == "G" && spore_s == "A") {
					++rip_ref2;
					return {
						ref: 2,
						type: "A"
					};
				}
				else if (ref2 == "C" && spore_s == "T") {
					++rip_ref2;
					return {
						ref: 2,
						type: "T"
					};
				}
				return {
					//illeg
					ref: -2,
					type: null
				};
			}
		});
		if (rip_ref1 >= 2 || rip_ref2 >= 2) {
			// if (ref1_pos == 1071253) {
			// 	debugger;
			// }

			// row[0] = (cmp_ref1_ref2 ? ColorID.identical : ColorID.dad) | (ref1 == "-" ? ColorID.indel_bit : 0);
			// row[1] = (cmp_ref1_ref2 ? ColorID.identical : ColorID.mom) | (ref2 == "-" ? ColorID.indel_bit : 0);

			row[0] = cmp_ref1_ref2 ? ColorID.identical : (ref1 != "-" ? ColorID.dad : all_prev_color[0]);
			row[1] = cmp_ref1_ref2 ? ColorID.identical : (ref2 != "-" ? ColorID.mom : all_prev_color[1]);

			spores_rip.forEach((rip, spore_idx) => {
				let bp = spores[spore_idx];
				if (rip.ref == -1) {
					row[2 + spore_idx] = (cmp_ref1_ref2 ? ColorID.identical : (bp == ref1 ? ColorID.dad : (bp == ref2 ? ColorID.mom : ColorID.diff))) | (a == "-" ? ColorID.indel_bit : 0);
				}
				else if (rip.ref == -2) {
					row[2 + spore_idx] = (cmp_ref1_ref2 ? ColorID.identical : (bp == ref1 ? ColorID.dad : (bp == ref2 ? ColorID.mom : ColorID.diff))) | (a == "-" ? ColorID.indel_bit : 0);
				}
				else if (rip.ref == 1) {
					row[2 + spore_idx] = ColorID.dad_rip;
				}
				else if (rip.ref == 2) {
					row[2 + spore_idx] = ColorID.mom_rip;
				}
			});

			return rip_ref1 >= 2 ? 1 : (rip_ref2 >= 2 ? 2 : 0);
		}
		else {
			return 0;
		}
		// let num_spores_rip = spores_rip.filter(a => a).length;
		// if (num_spores_rip == 2) {
		// 	// 2:2 rip
		// 	//return true;
		// }
		// else if (num_spores_rip > 2) {
		// 	spores.filter(v => ref1 == v).length == num_spores_rip;
		// 	spores.filter(v => ref2 == v).length == num_spores_rip;
		// 	//return true;
		// }
		// else if (num_spores_rip == 1) {
		// 	console.log({
		// 		type: "illeg", ref1_pos, ref1, ref2, a, b, c, d, co_state,
		// 	});
		// 	return false;
		// }
		// else {
		// 	return false;//not rip
		// }
	}
	else {
		let is_rip = 0;
		let rip;
		rip = check_rip(ref1, ref2, ColorID.dad_rip, ColorID.diff);
		if (rip == ColorID.dad_rip) {
			row[0] = cmp_ref1_ref2 ? ColorID.dad : ColorID.identical;
			row[1] = ref2 != "-" ? ColorID.mom : all_prev_color[1];
			is_rip = 1;
		}
		else if (rip == ColorID.mom_rip) {
			row[0] = ref1 != "-" ? ColorID.dad : all_prev_color[0];
			row[1] = cmp_ref1_ref2 ? ColorID.mom : ColorID.identical;
			is_rip = 2;
		}
		else {
			rip = check_rip(ref2, ref1, ColorID.mom_rip, ColorID.diff);
			if (rip == ColorID.dad_rip) {
				row[0] = cmp_ref1_ref2 ? ColorID.dad : ColorID.identical;
				row[1] = ref2 != "-" ? ColorID.mom : all_prev_color[1];
				is_rip = 1;
			}
			else if (rip == ColorID.mom_rip) {
				row[0] = ref1 != "-" ? ColorID.dad : all_prev_color[0];
				row[1] = cmp_ref1_ref2 ? ColorID.mom : ColorID.identical;
				is_rip = 2;
			}
		}

		return is_rip;
	}

	/**
	 * @param {string} ref_1
	 * @param {string} ref_2
	 * @param {number} col
	 * @param {number} col2
	 */
	function check_rip(ref_1, ref_2, col, col2) {
		let rip = 0;
		if (ref_1 == "G") { //GA AAAA
			if (ref_2 == "A") {
				if (a == "A" && b == "A" && c == "A" && d == "A") {
					//G->A
					rip = col;
					//parental_cmp_uint8array[pos] = col;
					row[2] = col;
					row[3] = col;
					row[4] = col;
					row[5] = col;
				}
			}
			else { //G? AA??
				let count_a = spores.filter(s => s == "A").length;
				let count_x = spores.filter(s => s == ref_2).length;
				if (count_a >= 2) {//2:2, 3:1
					if (count_x == (4 - count_a)) {//2:2, 3:1
						//G->A
						rip = col;
						//parental_cmp_uint8array[pos] = col;
						row[2] = a == "A" ? col : (a == ref_1 ? ColorID.dad : (a == ref_2 ? ColorID.mom : col2)) | (a == "-" ? ColorID.indel_bit : 0);
						row[3] = b == "A" ? col : (b == ref_1 ? ColorID.dad : (b == ref_2 ? ColorID.mom : col2)) | (b == "-" ? ColorID.indel_bit : 0);
						row[4] = c == "A" ? col : (c == ref_1 ? ColorID.dad : (c == ref_2 ? ColorID.mom : col2)) | (c == "-" ? ColorID.indel_bit : 0);
						row[5] = d == "A" ? col : (d == ref_1 ? ColorID.dad : (d == ref_2 ? ColorID.mom : col2)) | (d == "-" ? ColorID.indel_bit : 0);
					}
					else {
						//G??A
						rip = 0; //fill default color
					}
				}
			}
		}
		else if (ref_1 == "C") { //CT TTTT
			if (ref_2 == "T") {
				if (a == "T" && b == "T" && c == "T" && d == "T") {
					//C->T
					rip = col;
					//parental_cmp_uint8array[pos] = col;
					row[2] = col;
					row[3] = col;
					row[4] = col;
					row[5] = col;
				}
			}
			else { //C? TT??
				let count_t = spores.filter(s => s == "T").length;
				let count_x = spores.filter(s => s == ref_2).length;
				if (count_t >= 2) {//2:2, 3:1
					if (count_x == (4 - count_t)) {//2:2, 3:1
						//C->T
						rip = col;
						//parental_cmp_uint8array[pos] = col;
						row[2] = a == "T" ? col : (a == ref_1 ? ColorID.dad : (a == ref_2 ? ColorID.mom : col2)) | (a == "-" ? ColorID.indel_bit : 0);
						row[3] = b == "T" ? col : (b == ref_1 ? ColorID.dad : (b == ref_2 ? ColorID.mom : col2)) | (b == "-" ? ColorID.indel_bit : 0);
						row[4] = c == "T" ? col : (c == ref_1 ? ColorID.dad : (c == ref_2 ? ColorID.mom : col2)) | (c == "-" ? ColorID.indel_bit : 0);
						row[5] = d == "T" ? col : (d == ref_1 ? ColorID.dad : (d == ref_2 ? ColorID.mom : col2)) | (d == "-" ? ColorID.indel_bit : 0);
					}
					else {
						//C??T
						rip = 0; //fill default color
					}
				}
			}
		}
		return rip;
	}
}

function init_viewModel(seq_list) {
	region_rect = seq_list.map(_ => []);
	//parental_cmp_uint8array = new Uint8Array(seq_list[0].length);
	ref1_pos_uint32array = new Uint32Array(seq_list[0].length);
	pos_ref1_uint32array = new Uint32Array(seq_list[0].length);
	ref2_pos_uint32array = new Uint32Array(seq_list[0].length);
	pos_ref2_uint32array = new Uint32Array(seq_list[0].length);
	//rip_map_uint8array = new Uint8Array(seq_list[0].length);
	ref1_ref2_score_uint32array = new Uint32Array(seq_list[0].length);
	allMarker.list.forEach(a => a.values.splice(0));//clear all
	align_start_index = -1;
	align_end_index = -1;
}

function _merge_seg(seq_list, seg_snp) {
	for (let sid = 0; sid < seq_list.length; ++sid) {
		let start = 0;
		let prev = seg_snp[start][sid];
		let i = 1;
		for (; i < seq_list[0].length; ++i) {
			let col = seg_snp[i][sid];
			if (prev != col) {
				region_rect[sid].push({
					start: start,
					end: i,
					col: prev,
				});
				start = i;
				prev = col;
			}
		}
		if (region_rect[sid].length && i != region_rect[sid][region_rect[sid].length - 1].end) {
			region_rect[sid].push({
				start: start,
				end: i,
				col: prev,
			});
		}
	}
}

try {
	if (!globalThis.window) {
		module.exports.loadCmpData = loadCmpData;
		module.exports.initAnalyser = initAnalyser;
	}
}
catch (ex) {
}
