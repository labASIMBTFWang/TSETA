//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, loadChrLength, array_groupBy } = require("./util.js");
const { BlastnCoord, execAsync, exec_blastn, exec_blastn_Ex, parseBlastnResults, blastn_coord, isCollide, groupByOverlap } = require("./blastn_util.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);

const DEBUG = !!argv["--debug"];
const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const ref1_name = dataset.ref;
const ref2_name = dataset.refs[1];
const spore_1_name = dataset.progeny_list[0];
const spore_2_name = dataset.progeny_list[1];
const spore_3_name = dataset.progeny_list[2];
const spore_4_name = dataset.progeny_list[3];

let genome_name_list = [...dataset.refs, ...dataset.progeny_list];

const ref1_chr_list = loadChrLength(`./${ref1_name}.length.txt`).list;
const ref2_chr_list = loadChrLength(`./${ref2_name}.length.txt`).list;
const s1_chr_list = loadChrLength(`./${spore_1_name}.length.txt`).list;
const s2_chr_list = loadChrLength(`./${spore_2_name}.length.txt`).list;
const s3_chr_list = loadChrLength(`./${spore_3_name}.length.txt`).list;
const s4_chr_list = loadChrLength(`./${spore_4_name}.length.txt`).list;


if (process.argv[1] == __filename) {
	main();
}

function main() {
	const rDNA_nChr = Number(dataset.rDNA.nChr);

	for (let i = 1; i <= ref1_chr_list.length; ++i) {
		if (i != rDNA_nChr) {
			let input_path = `./tmp/mafft_ch${i}.fa`;
			let output_path = `./mafft_ch${i}.fa`;
			fs.createReadStream(input_path).pipe(fs.createWriteStream(output_path));//copy file
		}
	}

	restore_rdna();
}

async function restore_rdna() {
	const nChr = Number(dataset.rDNA.nChr);
	const chrIdx = nChr - 1;
	let rDNA_filePath = dataset.rDNA.sequence;

	let mafft_fa = readFasta(`./tmp/mafft_ch${nChr}.fa`);

	/** chrPos_to_multialign_posMap */
	let mafft_fa_posmap = {};
	/** multialign_to_chrPos_posMap */
	let chr_posmap = {};

	/** @type {{ [seqName: string]: rDNA_Data }} */
	let rdna_data = {};

	let chr_list = [
		ref1_chr_list[chrIdx].chr,
		ref2_chr_list[chrIdx].chr,
		s1_chr_list[chrIdx].chr,
		s2_chr_list[chrIdx].chr,
		s3_chr_list[chrIdx].chr,
		s4_chr_list[chrIdx].chr,
	];
	let promise = chr_list.map(async function (seqName, seqIndex) {
		const mafft_seq = mafft_fa[seqName];
		mafft_fa_posmap[seqName] = _chrPos_to_multialign_posMap(mafft_seq);
		chr_posmap[seqName] = _multialign_to_chrPos_posMap(mafft_seq);
		
		rdna_data[seqName] = await find_rDNA_use_blastn(rDNA_filePath, genome_name_list[seqIndex], nChr);
	});
	await Promise.all(promise);
	
	/** @type {{ [seqName: string]: string }} */
	let final_mafft_seq = {};

	function pos_to_ma(seqName) {
		const rd = rdna_data[seqName];

		let posmap = mafft_fa_posmap[seqName];

		let ma_start = posmap[rd.start - 1];
		let ma_end = posmap[rd.end - 1];

		const ma_seq = mafft_fa[seqName];

		while (ma_seq[ma_start - 1] == "-") {
			ma_start = ma_start - 1;
		}

		while (ma_seq[ma_end + 1] == "-") {
			ma_end = ma_end + 1;
		}

		return {
			ma_start,
			ma_end,
		};
	}

	let ma_rdna_range_list = chr_list.map(function (seqName, seqIndex) {
		let { ma_start, ma_end } = pos_to_ma(seqName);
		return {
			ma_start, ma_end
		};
	});
	let min_ma_start = Math.min(...ma_rdna_range_list.map(a => a.ma_start));
	let max_ma_end = Math.max(...ma_rdna_range_list.map(a => a.ma_end));//???

	let ma_start_delta = ma_rdna_range_list.map(a => a.ma_start - min_ma_start);
	
	// maybe rDNA
	let ex_range = [];
	let rdna_seq_list = chr_list.map(function (seqName, seqIndex) {//raw rDNA seq
		let posmap = chr_posmap[seqName];
		const rd = rdna_data[seqName];
		
		let start = posmap[min_ma_start - 1];
		let end = posmap[max_ma_end - 1];

		ex_range[seqIndex] = [start + 1, end];
		return rd.chrSeq.slice(start + 1, end);
	});
	let max_length = Math.max(...rdna_seq_list.map(a => a.length));
	
	chr_list.forEach(function (seqName, seqIndex) {
		const ma_seq = mafft_fa[seqName];
		const rd = rdna_data[seqName];

		console.log("rdna_seq_list[" + seqIndex + "].length", rdna_seq_list[seqIndex].length);

		let new_seq = rdna_seq_list[seqIndex].padEnd(max_length, "-");
		final_mafft_seq[seqName] = replace_seq_range(ma_seq, min_ma_start, max_ma_end, new_seq);

		//TODO: re-align IGS, rDNA repeat, IGS
	});
	
	console.log("rDNA multi align range", {
		start: min_ma_start,
		end: min_ma_start + max_length,
		before_end: max_ma_end,
	});

	let rDNA_info = {
		chr: nChr,
		
		//alignment index to position
		alignment_start: min_ma_start + 1,
		alignment_end: min_ma_start + max_length + 1,
		alignment_range_list: ma_rdna_range_list.map(a => [a.ma_start + 1, a.ma_end + 1]),

		alignment_delta_list: ma_start_delta,

		data: chr_list.map(function (seqName, seqIndex) {
			let info = rdna_data[seqName].info;

			//rebuild pos map
			{
				let posmap = __chrPos_to_multialign_posMap(final_mafft_seq[seqName]);

				info.alignment_repeats = [];
				for (let i = 0; i < info.repeats.length; ++i) {
					let start = info.repeats[i][0] - 1;
					let end = info.repeats[i][1] - 1;
					let ma_start = posmap[start - 1] + 1;//pos(bp)
					let ma_end = posmap[end - 1] + 1;//pos(bp)

					info.alignment_repeats[i] = [];
					info.alignment_repeats[i][0] = ma_start;
					info.alignment_repeats[i][1] = ma_end;
				}

				let ex_s = posmap[ex_range[seqIndex][0]];
				let ex_e = posmap[ex_range[seqIndex][1]];

				info.alignment_ex = [ex_s, ex_e];
			}

			return info;
		}),
	};
	fs.writeFileSync("rDNA_info.json", JSON.stringify(rDNA_info));

	if (DEBUG) {
		console.log("debug no output file");
	}
	else {
		let output_path = `./mafft_ch${nChr}.fa`;

		console.log("output fasta", Path.resolve(output_path));

		saveFasta(output_path, final_mafft_seq);
	}
}

class rDNA_Data {
	constructor() {
		//this.seq = "";
		this.start = 0;
		this.end = 0;
		this.chrSeq = "";
		
		/** @type {{ range: number[], repeats: number[][], alignment_repeats: number[][] }} */
		this.info = null;
	}
}

/**
 * @param {string} rDNA_filePath
 * @param {string} subject_genome_name
 * @param {number} nChr
 * @returns {Promise<rDNA_Data>}
 */
async function find_rDNA_use_blastn(rDNA_filePath, subject_genome_name, nChr) {
	const subject_chrInfo = loadChrLength(`./${subject_genome_name}.length.txt`).list;
	const subject_chr_name = subject_chrInfo[nChr - 1].chr;

	const subject_fa_filename = Path.join("./tmp/fasta", `${subject_genome_name}_${subject_chr_name}.fa`);

	let result_text = await exec_blastn_Ex(rDNA_filePath, subject_fa_filename, undefined, undefined, undefined, undefined, "-evalue 1e-5");
	let _table = parseBlastnResults(result_text);
	
	let max_len = Math.max(..._table.map(a => a.slen));
	/**
	 * TODO: Check Ribosomal DNA structure, IGS 18S ITS 5.8S ITS 26S IGS
	 * CBS1-1 min:7364, max:7835, 7835 / 7364 = 0.94
	 * range include IGS
	 */
	let group = _table.filter(a => (a.slen / max_len) >= 0.9);

	group = group.sort((a, b) => a.s_min - b.s_max);

	if (!DEBUG) {
		fs.writeFileSync(
			`blastn_rDNA_${subject_chr_name}.txt`,
			group.map(a => a.toArray().join("\t")).join("\n")
		);
	}

	let min_sstart = Math.min(...group.map(a => a.s_min));
	let max_send = Math.max(...group.map(a => a.s_max));

	console.log({
		subject: subject_genome_name,
		min_sstart, max_send,
		len: max_send - min_sstart,
		"_table.length": _table.length,
		"group.length": group.length,
	});

	if (!(max_send > min_sstart)) {
		throw new Error("if (!(max_send > min_sstart)) {");
	}
	else {		
		let raw_fa = readFasta(Path.join("./tmp/fasta", `${subject_genome_name}_${subject_chr_name}.fa`))[subject_chr_name];	

		// TDOD: save rDNA repeat position

		let info = {
			range: [min_sstart, max_send],
			repeats: group.map(a => [a.sstart, a.send]),
			alignment_repeats: null,
		};
		return {
			chrSeq: raw_fa,
			//seq: getSeq(raw_fa, min_sstart, max_send),
			start: min_sstart,
			end: max_send,
			info,
		};
	}
}

/**
 * index: 0 ~ length - 1
 * @param {string} ma_seq
 * @returns {number[]}
 */
function _chrPos_to_multialign_posMap(ma_seq) {
	let nPos = 0;
	let posmap = [];
	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		if (element != "-") {
			posmap[nPos] = index;
			++nPos;
		}
	}
	return posmap;
}

/**
 * index: 0 ~ length - 1
 * @param {string} ma_seq
 * @returns {number[]}
 */
function _multialign_to_chrPos_posMap(ma_seq) {
	let nPos = 0;
	let posmap = [];
	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		
		posmap[index] = nPos;

		if (element != "-") {
			++nPos;
		}
	}
	return posmap;
}



/**
 * pos: 1 ~ length
 * @param {string} ma_seq
 * @returns {number[]}
 */
function __chrPos_to_multialign_posMap(ma_seq) {
	let nPos = 1;
	let posmap = [];
	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		if (element != "-") {
			posmap[nPos] = index;
			++nPos;
		}
	}
	return posmap;
}

/**
 * pos: 1 ~ length
 * @param {string} ma_seq
 * @returns {number[]}
 */
function __multialign_to_chrPos_posMap(ma_seq) {
	let nPos = 1;
	let posmap = [];
	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		
		if (element == "-") {
			posmap[index] = Math.max(1, nPos - 1);
		}
		else {
			posmap[index] = nPos;
			++nPos;
		}
	}
	return posmap;
}

/**
 * @param {string} ma_seq
 * @param {number} ma_start
 * @param {number} ma_end
 * @param {string} replace_seq
 */
function replace_seq_range(ma_seq, ma_start, ma_end, replace_seq) {
	return getSeq(ma_seq, 1, ma_start) + replace_seq + getSeq(ma_seq, ma_end, ma_seq.length);
}

/**
 * @param {string} seq
 * @param {number} start - bp
 * @param {number} end - bp
 */
function getSeq(seq, start, end) {
	return seq.slice(start - 1, end);
}



