// @ts-check

const fs = require("fs");
const Path = require("path");

const { tsv_parse, _table_to_object_list, table_to_object_list } = require("./tsv_parser.js");
const { parse_GC_Content_table, GC_Content_Data } = require("./GC_content_util.js");


class GenomeInfo {
	/**
	 * @param {string} name
	 */
	constructor(name) {
		this.name = name;
		
		this.chr_list = GenomeInfo.loadChrLength(`output/${name}.length.txt`).list;
		this.chr_map = GenomeInfo.loadChrLength(`output/${name}.length.txt`).map;

		/** @type {{ [sChr:string]: string }} */
		this.fasta = {};
	}
	
	/**
	 * @returns {{ list: { index: number, chr: string, length: number }[], map: {[chrName:string]:{ index: number, chr: string, length: number }} }}
	 */
	static loadChrLength(file_path) {
		let text_tab = fs.readFileSync(file_path).toString();
		let tab = tsv_parse(text_tab);
		let rows = table_to_object_list(tab, ["index", "chr", "length"], { start_row: 1 });
		
		/** @type {{ index: number, chr: string, length: number }[]} */
		let chrList = rows.map(row => {
			return {
				index: Number(row.index),
				chr: String(row.chr),
				length: Number(row.length),
			};
		});
		chrList = chrList.sort((a, b) => a.index - b.index);

		/** @type {{[chrName:string]:{ index: number, chr: string, length: number }}} */
		let chrMap = {};
		rows.forEach(row => {
			chrMap[row.chr] = {
				index: Number(row.index),
				chr: String(row.chr),
				length: Number(row.length),
			};
		});
		
		return {
			map: chrMap,
			list: chrList,
		};
	}
}

class RibosomalDNA_Data {
	constructor() {
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

class LoadingFlags {
	constructor() {
		this.ribosomal_dna = true;
		this.gc_content = true;
	}
}

class Dataset {
	constructor() {
		//user input data

		//Chr#    start   end     GC%
		this.gc = {};

		this.ref = "";
		
		/** @type {string[]} */
		this.parental_list = [];

		/** @type {string[]} */
		this.progeny_list = [];

		/** @type {{[refName:string]:string}} - parental fasta file apth */
		this.parental = {};
		
		/** @type {{[refName:string]:string}} - progeny fasta file apth */
		this.progeny = {};
		
		/** @type {string[]} ref chr / parental[0] chr */
		this.chrNames = [];

		this.mafft = {};
		this.mafft.algorithm = "--localpair";
		this.mafft.default_algorithm = "";
		this.mafft.maxIterate = 1000;
		this.mafft.thread = 20;

		this.crossover = {
			closeCOsMinDistance: 5000,
		};

		/** file path */
		this["GC content"] = "";

		/** file path */
		this["AT-island"] = "";

		/** file path */
		this["GC content file"] = "";

		/** file path */
		this["AT-island file"] = "";

		/** @type {{[nChr:number]:number[]}} */
		this.centromere = {};

		/** @type {{[nChr:number]:number[][]}} */
		this.telomere = {};

		/** @type {RibosomalDNA_Data} */
		this.rDNA = new RibosomalDNA_Data();

		/** @type {RibosomalDNA_Position_info} */
		this.rDNA_info = new RibosomalDNA_Position_info();

		// internal property

		/** @type {string[]} */
		this.genomeNameList = [];

		/** @type {{[name:string]:string}} */
		this.genomeFileMap = {};
	}

	/**
	 * @type {"tetrad"|"SNP"}
	 */
	get mode() {
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
		let row = this["GC content"][ref_name][nChr].find(a => ref1_pos >= a.start && ref1_pos <= a.end);
		if (row) {
			return row.gc;
		}
		else {
			console.error({
				"GC content file": this["GC content file"],
				ref_name, nChr, ref1_pos,
				"s names": Object.keys(this["GC content"]),
				"chr names": Object.keys(this["GC content"][ref_name]),
				"rows.length": Object.keys(this["GC content"][ref_name][nChr].length),
				"chr.rows.length": Object.keys(this["GC content"][ref_name]).map(chrName => this["GC content"][ref_name][chrName].length).join(","),
				"max": this["GC content"][ref_name][nChr].sort((a, b) => b.end - a.end)[0],
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
		return this.genomeNameList.map(name => new GenomeInfo(name));
	}

	/**
	 * @param {string} dataset_path
	 * @param {LoadingFlags} flags
	 * @returns {LoadedDataset}
	 */
	static loadFromFile(dataset_path, flags = null) {
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

		/** @type {LoadedDataset} */
		// @ts-ignore
		let loaded = dataset;
		
		if (flags && flags.gc_content) {
			loaded["GC content"] = Dataset.__load_GC_content(dataset["GC content"]);
		}
		else {
			console.warn(process.argv[1], "no loading:", "GC content");
		}

		if (flags && flags.ribosomal_dna) {
			loaded.rDNA_info = Dataset.__load_rDNA_info_fromDataset(dataset_path);
		}
		else {
			console.warn(process.argv[1], "no loading:", "rDNA_info");
		}

		loaded.genomeNameList = [].concat(loaded.parental_list, loaded.progeny_list);
		
		loaded.genomeFileMap = {};

		loaded.parental_list.forEach(name => {
			loaded.genomeFileMap[name] = loaded.parental[name];
		});
		loaded.progeny_list.forEach(name => {
			loaded.genomeFileMap[name] = loaded.progeny[name];
		});

		return loaded;
	}
	
	/**
	* @param {string} filename
	* @returns {{[parentalName:string]:{[nChr:number]:GC_Content_Data[]}}}
	*/
	static __load_GC_content(filename) {
		if (!fs.existsSync(filename)) {
			return null;
			// console.error({
			// 	error: "No such file or directory",
			// 	path: filename,
			// 	absolute: Path.resolve(filename),
			// });
			// throw new Error("No such file or directory");
		}
		const text = fs.readFileSync(filename).toString();
		// @ts-ignore
		let table = table_to_object_list(tsv_parse(text), ["name", "chr", "start", "end", "gc"]);
		return parse_GC_Content_table(table);
	}

	/**
	 * @param {string} dataset_path
	 * @returns {RibosomalDNA_Position_info}
	 */
	static __load_rDNA_info_fromDataset(dataset_path) {
		if (!fs.existsSync(dataset_path)) {
			return null;
			// console.error({
			// 	error: "No such file or directory",
			// 	path: dataset_path,
			// 	absolute: Path.resolve(dataset_path),
			// });
			// throw new Error("No such file or directory");
		}
		let info_path = Path.join(Path.dirname(dataset_path), "output/rDNA_info.json");
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

// @ts-ignore
class LoadedDataset extends Dataset {
	constructor() {
		super();
		
		/** @type {{[parentalName:string]:{[nChr:number]:GC_Content_Data[]}}} */
		this["GC content"] = {};
	}
}


module.exports.Dataset = Dataset;
module.exports.LoadedDataset = LoadedDataset;

module.exports.GenomeInfo = GenomeInfo;

module.exports.RibosomalDNA_Data = RibosomalDNA_Data;
module.exports.RibosomalDNA_Position_info = RibosomalDNA_Position_info;

