// @ts-check

const fs = require("fs");
const Path = require("path");

const { tsv_parse, _table_to_object_list, table_to_object_list } = require("./tsv_parser.js");
const { parse_GC_COntent_table, GC_Content_Data } = require("./GC_content_util.js");


/**
 * @param {string} filename
 * @returns {{[parentalName:string]:{[nChr:number]:GC_Content_Data[]}}}
 */
function loadGCContent(filename) {
	const text = fs.readFileSync(filename).toString();
	// @ts-ignore
	let table = table_to_object_list(tsv_parse(text), ["name", "chr", "start", "end", "gc"]);
	return parse_GC_COntent_table(table);
}

class RibosomalDNA_Data {
	constructor() {
		this.range = [2323505, 2387639];
		this.nChr = 6;
		this.sequence = "rDNA.fa";
	}
}

class RibosomalDNA_Position_info {
	constructor() {
		this.chr = 6;
		this.alignment_start = 2600866;
		this.alignment_end = 0;
	}
}

class Dataset {
	constructor() {
		//Chr#    start   end     GC%
		this.gc = {};

		this.ref = "";
		
		/** @type {string[]} */
		this.refs = [];

		/** @type {string[]} */
		this.progeny_list = [];

		/** @type {{[refName:string]:string}} - parental fasta file apth */
		this.parental = {};
		
		/** @type {{[refName:string]:string}} - progeny fasta file apth */
		this.progeny = {};
		
		/** @type {string[]} */
		this.chrNames = [];

		this.crossover = {
			closeCOsMinDistance: 5000,
		};

		/** file path */
		this["GC content"] = "";

		/** file path */
		this["AT-island"] = "";

		/** @type {{[nChr:number]:number[]}} */
		this.centromere = {};

		/** @type {{[nChr:number]:number[][]}} */
		this.telomere = {};

		/** @type {RibosomalDNA_Data} */
		this.rDNA = new RibosomalDNA_Data();

		/** @type {RibosomalDNA_Position_info} */
		this.rDNA_info = new RibosomalDNA_Position_info();
	}
	
	/**
	 * @param {number} nChr
	 * @param {number} pos
	 */
	isIn_rDNA(nChr, pos) {
		if (!this.rDNA_info) {
			throw new Error("this.rDNA_info");
		}
		const rDNA_nChr = Number(this.rDNA_info.chr);
		const rDNA_ma_start = Number(this.rDNA_info.alignment_start);
		const rDNA_ma_end = Number(this.rDNA_info.alignment_end);
		if (nChr == rDNA_nChr) {
			if (pos >= rDNA_ma_start && pos <= rDNA_ma_end) {
				return true;
			}
		}
		return false;
	}

	/**
	 * @param {string} ref_name
	 * @param {number} nChr
	 * @param {number} ref1_pos pos in bp
	 */
	getGCByPos(ref_name, nChr, ref1_pos) {
		return this["GC content"][ref_name][nChr].find(a => ref1_pos > a.start && ref1_pos <= a.end).gc;
	}
	
	/**
	 * @param {number} nChr
	 * @param {number} ref1_pos pos in bp
	 */
	isInTelomere(nChr, ref1_pos) {
		let [[start1, end1], [start2, end2]] = this.telomere[nChr];
		return (ref1_pos >= start1 && ref1_pos <= end1) || (ref1_pos >= start2 && ref1_pos <= end2);
	}	
	/**
	 * @param {number} nChr
	 * @param {number} ref1_pos pos in bp
	 */
	isInCentromere(nChr, ref1_pos) {
		let [start, end] = this.centromere[nChr];
		return ref1_pos >= start && ref1_pos <= end;
	}

	/**
	 * @param {string} dataset_path
	 * @returns {LoadedDataset}
	 */
	static loadFromFile(dataset_path) {
		let obj = JSON.parse(fs.readFileSync(dataset_path).toString());
		let dataset = Dataset.__fromObject(obj);

		/** @type {LoadedDataset} */
		// @ts-ignore
		let loaded = dataset;
		
		loaded["GC content"] = loadGCContent(dataset["GC content"]);
		
		loaded.rDNA_info = Dataset.__load_rDNA_info_fromDataset(dataset_path);

		loaded.refs = Object.keys(dataset.parental);
		loaded.progeny_list = Object.keys(dataset.progeny);

		return loaded;
	}

	/**
	 * @param {string} dataset_path
	 * @returns {RibosomalDNA_Position_info}
	 */
	static __load_rDNA_info_fromDataset(dataset_path) {
		let info_path = Path.join(Path.dirname(dataset_path), "rDNA_info.json");
		if (fs.existsSync(info_path)) {
			let data = JSON.parse(fs.readFileSync(info_path).toString());
			return data;
		}
		else {
			return null;
		}
	}

	/**
	 * @param {any} obj
	 * @returns {Dataset}
	 */
	static __fromObject(obj) {
		let ds = new Dataset();
		Object.assign(ds, obj);
		return ds;
	}
}

class LoadedDataset extends Dataset {
	constructor() {
		super();
		
		/** @type {{[parentalName:string]:{[nChr:number]:GC_Content_Data[]}}} */
		this["GC content"] = {};
	}
}


module.exports.Dataset = Dataset;
module.exports.loadGCContent = loadGCContent;
