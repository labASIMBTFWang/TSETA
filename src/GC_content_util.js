
// @ts-check

class GC_Content_Data {
	constructor() {
		this.chr = 0;
		this.start = 0;
		this.end = 0;
		this.gc = 0;
	}

	/**
	 * @param {Partial<GC_Content_Data>} obj
	 * @returns {GC_Content_Data}
	 */
	static fromObject(obj) {
		let gc = new GC_Content_Data();
		gc.chr = Number(obj.chr);
		gc.gc = Number(obj.gc);
		gc.start = Number(obj.start);
		gc.end = Number(obj.end);
		return gc;
	}
}


/**
 * @param {{[key:string]:string}[]} table
 * @returns {{[parentalName:string]:{[nChr:number]:GC_Content_Data[]}}}
 */
function parse_GC_Content_table(table) {
	/** @type {{[parentalName:string]:GC_Content_Data[]}} */
	let name_group = {};

	/** @type {{[parentalName:string]:{[nChr:number]:GC_Content_Data[]}}} */
	let name_map = {};

	// group by name
	table.forEach(row => {
		if (!name_group[row.name]) {
			name_group[row.name] = [];
		}
		// @ts-ignore
		name_group[row.name].push(row);
	});
	//group by chr
	Object.keys(name_group).forEach(name => {
		/** @type {{[nChr:number]:GC_Content_Data[]}} */
		let chr_map = {};
		name_group[name].forEach(row => {
			if (!chr_map[row.chr]) {
				chr_map[row.chr] = [];
			}
			chr_map[row.chr].push(GC_Content_Data.fromObject(row));
		});
		name_map[name] = chr_map;
	});

	return name_map;
}

module.exports.parse_GC_Content_table = parse_GC_Content_table;
module.exports.GC_Content_Data = GC_Content_Data;
