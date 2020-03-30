// @ts-check

const fs = require("fs");
const Path = require("path");

const { tsv_parse, _table_to_object_list, table_to_object_list } = require("./tsv_parser.js");
const { parse_GC_Content_table, GC_Content_Data } = require("./GC_content_util.js");
const { readFasta } = require("./fasta_util.js");

const VERBOSE = process.argv.indexOf("--verbose") >= 0;

class ChromosomeData {
	constructor() {
		/** @type {number} - number of chromosome */
		this.index = null;
		/** @type {string} - chromosome name */
		this.chr = null;
		/** @type {number} - chromosome length*/
		this.length = null;
		/** @type {string} - chromosome fasta file path */
		this.path = null;
	}

	loadSeq() {
		return readFasta(this.path)[this.chr];
	}
}
class ChromosomeInfo {
	/**
	 * @param {ChromosomeData[]} list
	 * @param {{ [chrName:string]: ChromosomeData }} map
	 */
	constructor(list, map) {
		/** @type {ChromosomeData[]} */
		this.list = list;
		
		/** @type {{ [chrName:string]: ChromosomeData }} */
		this.map = map;
	}
}

class GenomeInfo {
	/**
	 * @param {Dataset} dataset
	 * @param {string} genome_name
	 */
	constructor(dataset, genome_name) {
		let info = GenomeInfo.loadChrLength(dataset, genome_name);
		
		this.name = genome_name;
		this.chr_list = info.list;
		this.chr_map = info.map;

		/** @type {{ [sChr:string]: string }} */
		this.fasta = {};
	}
	
	loadFasta() {
		return this.chr_list.reduce((obj, chrInfo) => Object.assign(obj, readFasta(chrInfo.path)), {});
	}
	
	/**
	 * @param {Dataset} dataset
	 * @param {string} genome_name
	 * @returns {ChromosomeInfo}
	 */
	static loadChrLength(dataset, genome_name) {
		let text_tab = fs.readFileSync(`${dataset.output_path}/${genome_name}.length.txt`).toString();
		let tab = tsv_parse(text_tab);
		let rows = table_to_object_list(tab, ["index", "chr", "length", "raw_chr_name"], { start_row: 1 });
		
		const _chrList = rows.map(row => {
			let data = new ChromosomeData();
			Object.assign(data, {
				index: Number(row.index),
				chr: String(row.chr),//seq name in fasta
				length: Number(row.length),
				path: `${dataset.tmp_path}/fasta/${String(row.chr)}.fa`,
			});
			return data;
		});
		const chrList = [..._chrList].sort((a, b) => a.index - b.index);

		/** @type {{ [chrName:string]: ChromosomeData }} */
		const chrMap = {};
		chrList.forEach(chrInfo => {
			chrMap[chrInfo.chr] = chrInfo;
		});
		
		return {
			map: chrMap,
			list: chrList,
		};
	}

	/**
	 * @param {string} genome_name
	 * @param {string} raw_chr_name
	 * @returns {string}
	 */
	static transformChrName(genome_name, raw_chr_name) {
		return `${genome_name}_${encodeURIComponent(raw_chr_name)}`;
	}

	/**
	 * @param {GenomeDataSet} dataset
	 * @param {string} genome_name
	 * @param {string} raw_chr_name
	 * @returns {string}
	 */
	static makeChrFilePath(dataset, genome_name, raw_chr_name) {
		return `${dataset.tmp_path}/fasta/${GenomeInfo.transformChrName(genome_name, raw_chr_name)}.fa`;
	}
}

class RibosomalDNA_Data {
	constructor() {
		/** @type {number} */
		this.nChr = null;
		/** @type {string} - rDNA fasta file path */
		this.sequence = null;
	}
}

class RibosomalDNA_Position_info {
	constructor() {
		/** @type {number} */
		this.chr = null;
		/** @type {number} */
		this.alignment_start = null;
		/** @type {number} */
		this.alignment_end = null;
	}
}

class GenomeDataSet {
	constructor() {
		this.name = "";

		/** @type {"tetrad"|"SNP"} */
		this.mode = "tetrad";

		this.ref = "";

		/** @type {string[]} */
		this.parental_list = [];

		/** @type {{[genomeName:string]:string}} - parental fasta file apth */
		this.parental = {};

		/** @type {string[]} */
		this.progeny_list = [];

		/** @type {{[genomeName:string]:string}} - progeny fasta file path*/
		this.progeny = {};
	}
	
	/** @type {string} */
	get output_path() {
		return encodeURIComponent(this.name);
	}
	
	/** @type {string} */
	get tmp_path() {
		return `${this.output_path}/tmp`;
	}
}

class MafftOptions {
	constructor() {
		this.algorithm = "localpair";
		this.default_algorithm = "";
		this.maxIterate = 1000;
		this.thread = 20;
	}
}

class Dataset extends GenomeDataSet {
	constructor() {
		super();

		//user input data
		
		/** @type {MafftOptions} */
		this.mafft = null;

		/** @type {{ closeCOsMinDistance: number }} */
		this.crossover = null;

		/** @type {{[parentalName:string]:{[nChr:number]:GC_Content_Data[]}}} */
		this.gc_content = null;
		// Object.defineProperty(this, "gc_content", {
		// 	get: function () {
		// 		throw new TypeError("gc");
		// 	},
		// });

		/** @type {number} */
		this.GC_Content_window = null;

		/** @type {string} file path */
		this.GC_Content_filePath = null;

		/** @type {{[nChr:number]:number[]}} */
		this.centromere = {};

		/** @type {{[nChr:number]:number[][]}} */
		this.telomere = {};

		/** @type {RibosomalDNA_Data} */
		this.rDNA = new RibosomalDNA_Data();

		// auto generate

		/** @type {RibosomalDNA_Position_info} */
		this.rDNA_info = null;

		// internal property

		/** @type {string[]} */
		this.genomeNameList = [];
	}

	/**
	 * @returns {"tetrad"|"SNP"}
	 */
	auto_detect_mode() {
		if (this.parental_list.length == 2) {
			if (this.progeny_list.length % 4 == 0) {
				//Tetrad analysis
				return "tetrad";
			}
		}
		return "SNP";
		// if (this.parental_list.length == 2) {
		// 	if (this.progeny_list.length == 4) {
		// 		return "SNP CO InDel";
		// 	}
		// 	else if (this.progeny_list.length > 0) {
		// 		return "2 parental +progeny SNP";
		// 	}
		// 	else {
		// 		return "2 parental SNP";
		// 	}
		// }
		// else if (this.parental_list.length == 1 && this.progeny_list.length > 0) {
		// 	return "1 parental +progeny SNP";
		// }
		// else {
		// 	throw new Error("");
		// }
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
		let row = this.gc_content[ref_name][nChr].find(a => ref1_pos >= a.start && ref1_pos <= a.end);
		this.gc_content[ref_name][nChr][ref1_pos / this.GC_Content_window];
		if (row) {
			return row.gc;
		}
		else {
			console.error({
				GC_Content_filePath: this.GC_Content_filePath,
				ref_name, nChr, ref1_pos,
				"s names": Object.keys(this.gc_content),
				"chr names": Object.keys(this.gc_content[ref_name]),
				"rows.length": Object.keys(this.gc_content[ref_name][nChr].length),
				"chr.rows.length": Object.keys(this.gc_content[ref_name]).map(chrName => this.gc_content[ref_name][chrName].length).join(","),
				"max": this.gc_content[ref_name][nChr].sort((a, b) => b.end - a.end)[0],
			});
			//throw new Error("getGCByPos(ref_name, nChr, ref1_pos) {");
		}
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

	loadGenomeInfoList() {
		return this.genomeNameList.map(gName => new GenomeInfo(this, gName));
	}
	
	loadGenomeInfoMap() {
		let map = {};
		this.genomeNameList.map(gName => map[gName] = new GenomeInfo(this, gName));
		return map;
	}

	/**
	 * @param {string} dataset_path
	 * @returns {Dataset}
	 */
	static loadFromFile(dataset_path) {
		if (!fs.existsSync(dataset_path)) {
			console.error({
				error: "No such file or directory",
				path: dataset_path,
				absolute: Path.resolve(dataset_path),
			});
			throw new Error("No such file or directory");
		}
		let obj = JSON.parse(fs.readFileSync(dataset_path).toString());
		let dataset = Dataset.__fromObject(obj);

		/** @type {Dataset} */
		// @ts-ignore
		let loaded = VERBOSE ? (new Proxy(dataset, {
			get: function (target, propertyKey, receiver) {
				let value = Reflect.get(target, propertyKey, receiver);
				if (value === null || value === undefined) {
					console.log(`dataset["${propertyKey.toString()}"] =>`, value);
					debugger;
				}
				return value;
			},
		})) : dataset;
		
		//loaded.gc_content = Dataset.__load_GC_content(dataset.GC_Content_filePath);

		//loaded.rDNA_info = Dataset.__load_rDNA_info_fromDataset(dataset_path);

		loaded.genomeNameList = [].concat(loaded.parental_list, loaded.progeny_list);

		Object.defineProperty(loaded, "$path", {
			enumerable: false,
			writable: false,
			configurable: false,
			value: dataset_path,
		});

		if (VERBOSE) {
			// @ts-ignore
			if (globalThis.dataset instanceof Dataset) {
				// @ts-ignore
				console.warn("dataset loaded:", globalThis.dataset.$path);
			}
		}
		Object.defineProperty(globalThis, "dataset", {
			value: loaded,
			configurable: true,
		});

		return loaded;
	}
	
	load_GC_content() {
		if (fs.existsSync(this.GC_Content_filePath)) {
			const text = fs.readFileSync(this.GC_Content_filePath).toString();
			// @ts-ignore
			let table = table_to_object_list(tsv_parse(text), ["name", "chr", "start", "end", "gc"]);
			this.gc_content = parse_GC_Content_table(table);
		}
	}

	load_rDNA_info() {
		if (fs.existsSync(`${dataset.output_path}/rDNA_info.json`)) {
			this.rDNA_info = JSON.parse(fs.readFileSync(`${dataset.output_path}/rDNA_info.json`).toString());
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

module.exports.MafftOptions = MafftOptions;

module.exports.GenomeDataSet = GenomeDataSet;
module.exports.Dataset = Dataset;

module.exports.GenomeInfo = GenomeInfo;

module.exports.RibosomalDNA_Data = RibosomalDNA_Data;
module.exports.RibosomalDNA_Position_info = RibosomalDNA_Position_info;

