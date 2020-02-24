// @ts-check

const fs = require("fs");

const { LoadedDataset } = require("./dataset.js");


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


/**
 * @param {LoadedDataset} dataset
 * @param {boolean} [merge_centromere]
 * @returns {{[nChr:number]:MyCoord[]}}
 */
function loadFragIdList(dataset, merge_centromere = true) {
	/** @type {{[nChr:number]:MyCoord[]}} */
	const all_chr_frag_list = {};

	for (let nChr = 1; nChr <= dataset.chrNames.length; ++nChr) {
		try {
			const coords = load_my_coord(`tmp/multi_coord_ch${nChr}.txt`);

			if (merge_centromere) {
				const cen_range = (function get_centromere() {
					let { 0: start, 1: end } = dataset.centromere[nChr];
					return {
						start, end,
					};
				})();

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
			}//if (merge_centromere)

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

function load_my_coord(filename) {
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

module.exports.loadFragIdList = loadFragIdList;
module.exports.load_my_coord = load_my_coord;
module.exports.MyCoord = MyCoord;
