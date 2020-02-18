//@ts-check

const fs = require("fs");

const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { readFasta } = require("./fasta_util.js");


class Crossover {
	constructor() {
		this.tetrad = "";
		this.chr = 0;
		this.chr_len = 0;
		this.pos = 0;
		this.co_type = 0;
		this.GCasso_tract_len = 0;
		this.GCasso_marker = 0;
		this.snp_start_out = 0;
		this.snp_start_in = 0;
		this.snp_end_in = 0;
		this.snp_end_out = 0;
		this.type = "";
		this.why_remove = "";
		this.from = ["", "", "", ""];
	}

	/**
	 * @param {Partial<Crossover>} obj
	 */
	static fromObject(obj) {
		let co = new Crossover();
		co.tetrad = obj.tetrad;
		co.chr = Number(obj.chr);
		co.chr_len = Number(obj.chr_len);
		co.pos = Number(obj.pos);
		co.co_type = Number(obj.co_type);
		co.GCasso_tract_len = Number(obj.GCasso_tract_len);
		co.GCasso_marker = Number(obj.GCasso_marker);
		co.snp_start_out = Number(obj.snp_start_out);
		co.snp_start_in = Number(obj.snp_start_in);
		co.snp_end_in = Number(obj.snp_end_in);
		co.snp_end_out = Number(obj.snp_end_out);
		co.type = obj.type;
		co.why_remove = obj.why_remove;
		co.from = String(obj.from).split(",");
		return co;
	}

	static get table_header_list() {
		return [
			"tetrad",
			"chr", "chr_len",
			"pos",
			"co_type", "GCasso_tract_len", "GCasso_marker",
			"snp_start_out", "snp_start_in", "snp_end_in", "snp_end_out",
			"type", "why_remove",
			"from"
		];
	}
	
	/**
	 * @param {string} dataName
	 * @param {string} tetrad
	 * @param {"co"|"nco"} type
	 */
	static getFilePath(dataName, tetrad, type) {
		return `${dataName}_${tetrad}_${type}.txt`.replace(/\//g, "_").replace(/ /g, "_");
	}
	/**
	 * @param {string} dataName
	 * @param {string} tetrad
	 * @param {string|number} chr
	 * @param {"co"|"nco"} type
	 */
	static getFilePathChr(dataName, tetrad, chr, type) {
		return `${dataName}_${tetrad}_${chr}_${type}.txt`.replace(/\//g, "_").replace(/ /g, "_");
	}
	/**
	 * @param {string} dataName
	 * @param {string} tetrad
	 * @param {string|number} chr
	 */
	static getFilePrefix(dataName, tetrad, chr) {
		return `${dataName}_${tetrad}_${chr}`.replace(/\//g, "_").replace(/ /g, "_");
	}

	/**
	 * @param {string} file_path
	 * @returns {Crossover[]}
	 */
	static loadTable(file_path) {
		let text = fs.readFileSync(file_path).toString();
		let rows = table_to_object_list(tsv_parse(text), Crossover.table_header_list, { start_row: 1 });
		return rows.map(a => Crossover.fromObject(a));
	}
}


module.exports.Crossover = Crossover;

