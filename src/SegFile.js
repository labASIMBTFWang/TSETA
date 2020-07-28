
// @ts-check

const fs = require("fs");

const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");


class SegRow {
	constructor() {
		/** @type {number} */
		this.chr = 0;//1~N
		/** @type {number} ref position */
		this.pos = -1;
		/** @type {number} */
		this.a = -1;
		/** @type {number} */
		this.b = -1;
		/** @type {number} */
		this.c = -1;
		/** @type {number} */
		this.d = -1;
		{
			/** @type {number} */
			this[1] = -1;//a
			/** @type {number} */
			this[2] = -1;//b
			/** @type {number} */
			this[3] = -1;//c
			/** @type {number} */
			this[4] = -1;//d
		}
		/** @type {""|"type=2"|"type=3"|"type=4"} -  */
		this.type = "";
		/** @type {"align"|"indel=1"|"indel=2"|"indel=3"|"indel=4"} */
		this.indel = "align";
		
		/** @type {""|"-"|"A"|"C"|"G"|"T"} */
		this.r1 = "";
		/** @type {""|"-"|"A"|"C"|"G"|"T"} */
		this.r2 = "";
		/** @type {""|"-"|"A"|"C"|"G"|"T"} */
		this.s1 = "";
		/** @type {""|"-"|"A"|"C"|"G"|"T"} */
		this.s2 = "";
		/** @type {""|"-"|"A"|"C"|"G"|"T"} */
		this.s3 = "";
		/** @type {""|"-"|"A"|"C"|"G"|"T"} */
		this.s4 = "";

		/** @type {""|"G"|"C"|"?G"|"?C"} - is rip */
		this.rip = "";

		/** @type {number} - GC content % */
		this.gc = -1;//dataset

		/** @type {""|"centromere"|"telomere"} - is centromere or telomere */
		this.ct = "";//dataset
	}

	is_4n0_markers() {
		return (
			this.s1 == "-" &&
			this.s2 == "-" &&
			this.s3 == "-" &&
			this.s4 == "-"
		);
	}

	isInDel() {
		return this.indel.startsWith("indel");
	}
	
	/** @type {number} */
	get nChr() { return this.chr; }
	set nChr(val) { this.chr = val; }
	
	/** @type {number} */
	get [1]() { return this.a; }
	set [1](val) { this.a = val; }
	
	/** @type {number} */
	get [2]() { return this.b; }
	set [2](val) { this.b = val; }
	
	/** @type {number} */
	get [3]() { return this.c; }
	set [3](val) { this.c = val; }
	
	/** @type {number} */
	get [4]() { return this.d; }
	set [4](val) { this.d = val; }

	static keys() {
		const header = [
			//7	6
			"nChr", "pos",
			//0	0	1	0
			"a", "b", "c", "d",
			//type=3	indel
			"type", "indel",
			//-	C	C	C	-	C
			"r1", "r2", "s1", "s2", "s3", "s4",
			//	36.16	telomere
			"rip", "gc", "ct"
		];
		return header;
	}

	/**
	 * @param {Partial<SegRow>} obj
	 * @returns {SegRow}
	 */
	static formObject(obj) {
		let row = Object.assign(new SegRow(), obj);
		row.chr = row.nChr;
		row[1] = row.a;
		row[2] = row.b;
		row[3] = row.c;
		row[4] = row.d;
		return row;
	}
}

/**
 * @param {string} path_to_segfile
 * @returns {{[nChr:number]:SegRow[]}}
 */
function loadSegFile(path_to_segfile) {
	const text_tab = fs.readFileSync(path_to_segfile).toString();
	const _rows = table_to_object_list(tsv_parse(text_tab), SegRow.keys(), { start_row: 1 });
	console.log({ "row.length": _rows.length });
	
	const rows = _rows.map(data => SegRow.formObject(data));
	
	//console.log("loadSegFile", rows[0]);
	
	/** @type {{[nChr:number]:SegRow[]}} */
	const group = {};
	rows.forEach(row => {
		let list = group[row.nChr];
		if (!list) {
			list = group[row.nChr] = [];
		}
		row.nChr = Number(row.nChr);
		row.pos = Number(row.pos);
		list.push(row);
	});
	console.log("load segfile chr:", Object.keys(group));

	return group;
}

module.exports.loadSegFile = loadSegFile;
module.exports.SegRow = SegRow;

