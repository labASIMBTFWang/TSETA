//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, array_groupBy, program_log } = require("./util.js");
const { BlastnCoord, execAsync, exec_blastn, exec_blastn_Ex, parseBlastnResults, blastn_coord, isCollide, groupByOverlap } = require("./blastn_util.js");
const { run_mafft } = require("./run_mafft.js");
const { readFasta, saveFasta, chrPos_to_multialign_posMap } = require("./fasta_util.js");
const { validation_chr } = require("./validation_seq.js");
const { Dataset, RibosomalDNA_Data } = require("./dataset.js");
const { loadFragIdList, MyCoord } = require("./load_frag_list.js");
const { join_chr_frag } = require("./join_chr_frag.js");

const argv = argv_parse(process.argv);

const DEBUG = !!argv["--debug"];
const VERBOSE = process.argv.indexOf("--verbose") >= 0;
const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const genome_info_list = dataset.loadGenomeInfoList();

main();

async function main() {
	program_log(`${dataset.name}.log.txt`, "start");

	if (argv["-chr"] != null && argv["-rDNA"] != null) {
		let nChr = Number(argv["-chr"]);
		let fasta_filepath = String(argv["-rDNA"] || "");

		if (!Number.isSafeInteger(nChr)) {
			console.error("error", "rDNA", "number of chromosome:", nChr);
		}

		if (!fs.existsSync(fasta_filepath)) {
			console.error("error", "rDNA", "fasta file path:", fasta_filepath);
			return;
		}

		dataset.rDNA = new RibosomalDNA_Data();
		dataset.rDNA.nChr = nChr;
		dataset.rDNA.sequence = fasta_filepath;
		fs.writeFileSync(argv_dataset_path, JSON.stringify(dataset, null, "\t"));
	}

	try {
		await align_rdna();
		
		// for (let i = 1; i <= genome_info_list[0].chr_list.length; ++i) {
		// 	validation_chr(i, dataset.output_path, true);
		// }
	}
	catch (ex) {
		throw ex;
	}
	finally {
		fs.writeFileSync(argv_dataset_path, JSON.stringify(dataset, null, "\t"));
		if (VERBOSE) {
			console.log("save:", argv_dataset_path);
		}
	}

	console.log("next step:", "analysis");
	console.log("snp calling command:", `node ./src/snp_summary.js -dataset ${argv_dataset_path}`);
	console.log("tetrad analysis command:", `node ./src/tetrad_summary.js -dataset ${argv_dataset_path} -min-co 5000`);

	program_log(`${dataset.name}.log.txt`, "exit");
}

async function align_rdna() {
	const all_chr_frag_list = loadFragIdList(dataset);

	const nChr = Number(dataset.rDNA.nChr);
	const chrIdx = nChr - 1;
	let rDNA_filePath = dataset.rDNA.sequence;

	/** @type {{ [seqName: string]: rDNA_Data }} */
	let rdna_data = {};

	let chr_list = genome_info_list.map(genome_info => genome_info.chr_list[chrIdx].chr);

	let promise = chr_list.map(async function (genomeChrSeqName, genomeIndex) {
		rdna_data[genomeChrSeqName] = await find_rDNA_use_blastn(rDNA_filePath, genomeIndex, nChr);
	});
	await Promise.all(promise);

	// get rDNA region frags

	const rDNA_range = chr_list.map(function (genomeChrSeqName, genomeIndex) {
		return {
			start: rdna_data[genomeChrSeqName].region[0],
			end: rdna_data[genomeChrSeqName].region[1],
		};
	});
	/** @type {number[]} */
	const rDNA_raw_frag_idx_list = [];

	const rDNA_raw_frags = all_chr_frag_list[nChr].filter((coord, coord_idx) => {
		const ks = Object.keys(coord.start);
		ks.shift();//remove search
		return ks.some((k, i) => {
			if (
				coord.start[k] <= rDNA_range[i].end &&
				coord.end[k] >= rDNA_range[i].start
			) {
				all_chr_frag_list[nChr][coord_idx].removed = true;
				rDNA_raw_frag_idx_list.push(coord_idx);

				if (VERBOSE) {
					console.log(coord.id);
				}

				return true;
			}
		});
	});

	// load frags

	/** @type {{ [chr: string]: string }} */
	const concat_raw_seq = {};
	//
	rDNA_raw_frags.forEach((coord, frag_idx) => {
		const fragId = coord.id;
		const fasta_filename = `ch${nChr}_${fragId}.fa`;
		const input_path = `${dataset.tmp_path}/seq_frag/${fasta_filename}`;
		const fa = readFasta(input_path);
		
		chr_list.forEach(function (genomeChrSeqName, genomeIndex) {
			if (!concat_raw_seq[genomeChrSeqName]) {
				concat_raw_seq[genomeChrSeqName] = "";
				// concat_raw_seq_aa[genomeChrSeqName] = [];
			}
			concat_raw_seq[genomeChrSeqName] += fa[genomeChrSeqName];
			// concat_raw_seq_aa[genomeChrSeqName] = concat_raw_seq_aa[genomeChrSeqName].concat(...)
		});
	});

	// seq.splice

	/** @type {{ [chr: string]: string[] }} */
	const extract_raw_seq = {};
	/** @type {{ [chr: string]: string }} */
	const extract_raw_last_seq = {};
	//
	chr_list.forEach(function (genomeChrSeqName, genomeIndex) {
		const ks = Object.keys(rDNA_raw_frags[0].start);
		ks.shift();//remove search
		const frag_start = rDNA_raw_frags[0].start[ks[genomeIndex]];

		// get repeat start end (min, max)
		const repeats = rdna_data[genomeChrSeqName].repeats.map(range => {
			const [start, end] = [
				range[0] - frag_start,
				range[1] - frag_start,
			].sort((a, b) => a - b);
			return {
				start,
				end,
			};
		});
		if (VERBOSE) {
			console.log({
				frag_start,
			});
		}

		/*
		 step 1: *
		 step 2: @
		 step 3: %
		 | <- <- <- < |
		 |* *  *  @  %
		   <-  *  @  %
			  <-  @  %
				 <-  %
					< |
		 */
		
		const concat_raw_seq_aa = [...concat_raw_seq[genomeChrSeqName]];
		
		// step 1
		// splice 0 to rDNA.start - 1
		repeats.reduce((prev_end, range, range_idx) => {
			const splice_end = range.start - 1;
			const length = splice_end - prev_end + 1;

			if (VERBOSE) {
				console.log({
					range_idx,
					prev_end,
					splice_end,
				});
			}

			const frag_a = concat_raw_seq_aa.splice(0, length);
			const frag_seq = frag_a.join("");
			
			if (length <= 0) {
				if (VERBOSE) {
					console.log({
						"splice": "0 to rDNA.start - 1",
						"length = splice_end - prev_end + 1": length,
						"frag_a.length": frag_a.length,
						"frag_seq.length": frag_seq.length,
						genomeChrSeqName,
						genomeIndex,
						prev_end,
						range,
					});
				}
			}

			if (!extract_raw_seq[genomeChrSeqName]) {
				extract_raw_seq[genomeChrSeqName] = [];
			}
			extract_raw_seq[genomeChrSeqName].push(frag_seq);

			return range.start;
		}, 0);

		// step 2
		// splice last repeats:
		// rDNA.start to rDNA.end
		const last_rDNA = repeats.slice(-1)[0];
		if (last_rDNA) {
			const range = last_rDNA;
			const prev_end = last_rDNA.start;
			const splice_end = range.end;
			const length = splice_end - prev_end + 1;

			if (VERBOSE) {
				console.log({
					range_idx: repeats.length,
					prev_end,
					splice_end,
				});
			}

			const frag_a = concat_raw_seq_aa.splice(0, length);
			const frag_seq = frag_a.join("");
		
			if (length <= 0) {
				if (VERBOSE) {
					console.log({
						"splice": "rDNA.start to rDNA.end",
						"length = splice_end - prev_end + 1": length,
						"frag_a.length": frag_a.length,
						"frag_seq.length": frag_seq.length,
						genomeChrSeqName,
						genomeIndex,
						prev_end,
						range,
					});
				}
			}
		
			if (!extract_raw_seq[genomeChrSeqName]) {
				extract_raw_seq[genomeChrSeqName] = [];
			}
			extract_raw_seq[genomeChrSeqName].push(frag_seq);
		}

		// step 3
		// splice last seq
		// rDNA.end to region.length
		const last_frag_a = concat_raw_seq_aa.splice(0);
		const last_frag_seq = last_frag_a.join("");
		if (VERBOSE) {
			console.log(genomeChrSeqName, "last_frag_seq.length", last_frag_seq.length);
		}
		extract_raw_last_seq[genomeChrSeqName] = last_frag_seq;

		if (VERBOSE) {
			const n_frag = extract_raw_seq[genomeChrSeqName].length;
			console.log(genomeChrSeqName, "n_frag", n_frag);
		}
	});

	// chr_list.forEach(function (genomeChrSeqName, genomeIndex) {
	// 	/** @type {{ [chr: string]: string }} */
	// 	const fag = {};
	// 	fa_frag_group.push(fag);
	// });

	const rDNA_repeatSeqList_genomeList = chr_list.map(seq_name => extract_raw_seq[seq_name]);

	const max_num_repeats = Math.max(...rDNA_repeatSeqList_genomeList.map(repeats => repeats.length));
	
	// /** @type {{ [chr: string]: string }[]} */
	// const fa_frag_group = [];

	// make group

	const fa_groups = [...Array(max_num_repeats)].map((_, repIdx) => {
		return rDNA_repeatSeqList_genomeList.map((repeats_list, seqIdx) => {
			return repeats_list[repIdx];
		})
	});
	fa_groups.push(chr_list.map(seq_name => extract_raw_last_seq[seq_name]));

	if (VERBOSE) {
		console.log(fa_groups.map(s => s.map(ss => ss ? ss.length : null)));
	}

	// save frag

	/** @type {string[]} */
	const frag_id_list = [];
	/** @type {string[]} */
	const fa_name_list = [];

	fa_groups.forEach(function save_frag(group, group_idx) {
		/** @type {{ [chr: string]: string }} */
		const fa = {};
		
		chr_list.forEach((seq_name, seq_idx) => {
			if (group[seq_idx]) {
				fa[seq_name] = group[seq_idx];
			}
		});

		const frag_id = `rDNA_${group_idx}`;
		const fa_name = `ch${nChr}_${frag_id}.fa`;
		const fa_path = `${dataset.tmp_path}/seq_frag/${fa_name}`;

		saveFasta(fa_path, fa);
		
		frag_id_list.push(frag_id);
		fa_name_list.push(fa_name);
	});

	// run mafft

	const reAlign = true;
	await fa_name_list.filter(a => a).reduce(async function (prev_promise, file_name, grpIdx) {
		await prev_promise;

		const input_path = `${dataset.tmp_path}/seq_frag/${file_name}`;
		const output_path = `${dataset.tmp_path}/mafft_seq_frag/mafft_${file_name}`;

		if (VERBOSE) {
			console.log({ grpIdx });
			console.log(fa_groups[grpIdx].length);
			console.log("fa_groups->length", fa_groups[grpIdx].map(s => s ? s.length : null));
			{
				chr_list.forEach((seq_name, seq_idx) => {
					if (fa_groups[grpIdx][seq_idx]) {
						console.log(seq_name, fa_groups[grpIdx][seq_idx] ? fa_groups[grpIdx][seq_idx].length : null);
					}
				});
			}
		}

		const n_fa_seq = fa_groups[grpIdx].filter(a => a).length;
		if (VERBOSE) {
			console.log("n_fa_seq", n_fa_seq);
		}
		
		if (n_fa_seq >= 2) {
			try {
				return await run_mafft(input_path, output_path, dataset.mafft.algorithm, nChr, `rDNA_${(grpIdx + 1)}`, reAlign);
			}
			catch (ex) {
				console.error(ex);
				
				try {
					return await run_mafft(input_path, output_path, dataset.mafft.default_algorithm, nChr, `rDNA_${(grpIdx + 1)}`, reAlign);
				}
				catch (ex) {
					console.error(ex);
				}
			}
		}
		else {
			// if (fs.existsSync(output_path)) {
			// 	if (VERBOSE) {
			// 		console.log("rDNA > skip exist:", output_path);
			// 	}
			// }
			// else {
				if (VERBOSE) {
					console.log("copy file:", output_path);
				}
				await fs.promises.copyFile(input_path, output_path);
			// }
		}
	}, Promise.resolve(true));

	// c f

	// modify list

	/**
	 * @type {MyCoord[]}
	 */
	const new_rDNA_frags = [];

	{
		const raw_first = rDNA_raw_frags[0];
		const raw_last = rDNA_raw_frags[rDNA_raw_frags.length - 1];

		if (VERBOSE) {
			console.log("raw_first.constructor.name", raw_first.constructor.name);
		}

		const template = MyCoord.prototype.clone.call(raw_first);

		delete template.removed;
		template.start.search = raw_first.start.search;
		template.end.search = raw_last.end.search;

		const ks = Object.keys(template.start);
		ks.shift();//remove search
		
		fa_groups.forEach((fa, grpIdx) => {
			const copy = template.clone();

			copy.id = frag_id_list[grpIdx];//`${raw_first.id}_${raw_last.id}_${grpIdx}`;
			
			ks.forEach((k, i) => {
				const genomeChrSeqName = chr_list[i];
				const repeat = rdna_data[genomeChrSeqName].repeats[grpIdx];

				if (repeat) {
					copy.start[k] = repeat[0];
					copy.end[k] = repeat[1];
				}
				else {
					copy.start[k] = null;
					copy.end[k] = null;
				}
			});

			new_rDNA_frags.push(copy);
		});
	}

	const new_frag_list = all_chr_frag_list[nChr].slice(0);// clone

	if (VERBOSE) {
		console.log("o frags", all_chr_frag_list[nChr].length);
	}

	// modify
	const splice_start = rDNA_raw_frag_idx_list[0];
	new_frag_list.splice(splice_start, rDNA_raw_frag_idx_list.length, ...new_rDNA_frags);

	if (VERBOSE) {
		console.log("n frags", new_frag_list.length);

		console.log("r frags", new_rDNA_frags.length);
	}
	
	const output_tmp_path = `${dataset.tmp_path}/re-align_ch${nChr}.fa`;
	const output_path = `${dataset.output_path}/mafft_ch${nChr}.fa`;
	const final_mafft_seq = join_chr_frag(nChr, new_frag_list, output_tmp_path, { padding: true, chr_list: chr_list });
	await fs.promises.copyFile(output_tmp_path, output_path);

	const rDNA_data_list = chr_list.map(function (seqName, seqIndex) {
		const data = rdna_data[seqName];

		//rebuild pos map
		{
			if (VERBOSE) {
				console.log(seqName, Object.keys(final_mafft_seq));
			}

			const posmap = chrPos_to_multialign_posMap(final_mafft_seq[seqName]);

			// arrow
			// info.alignment_repeats = [];

			data.alignment_region = data.region.map(a => posmap[a]);
			
			for (let i = 0; i < data.repeats.length; ++i) {
				const start = data.repeats[i][0];
				const end = data.repeats[i][1];
				const ma_start = posmap[start] + 1; //pos(bp)
				const ma_end = posmap[end] + 1; //pos(bp)

				data.alignment_repeats[i] = [];
				data.alignment_repeats[i][0] = ma_start;
				data.alignment_repeats[i][1] = ma_end;
			}
		}

		return data;
	});
	const rDNA_info = {
		chr: nChr,//??
		nChr: nChr,

		//alignment index to position

		// border
		region_start: Math.min(...rDNA_data_list.map(a => a.alignment_region[0])),
		region_end: Math.min(...rDNA_data_list.map(a => a.alignment_region[1])),
		
		// !important
		data: rDNA_data_list,
	};
	fs.writeFileSync(`${dataset.output_path}/rDNA_info.json`, JSON.stringify(rDNA_info));

	validation_chr(nChr, dataset.output_path, true);
}

// 20200727
class rDNA_Data {
	constructor() {
		/**
		 * nChr -> [1,n]
		 */
		this.nChr = 0;

		this.strand = 0;

		/**
		 * @typedef seq_range
		 * @type {[number, number] | number[]}
		 */

		/**
		 * region
		 * @type {seq_range}
		 */
		this.region = [0, 0];

		/**
		 * alignment region
		 * @type {seq_range}
		 */
		this.alignment_region = [0, 0];

		/**
		 * repeat
		 * @type {seq_range[]}
		 */
		this.repeats = [];
		
		/**
		 * alignment repeat
		 * @type {seq_range[]}
		 */
		this.alignment_repeats = [];
	}
}

/**
 * @param {string} rDNA_filePath
 * @param {number} genomeIndex
 * @param {number} nChr
 * @returns {Promise<rDNA_Data>}
 */
async function find_rDNA_use_blastn(rDNA_filePath, genomeIndex, nChr) {
	const subject_genome_name = dataset.genomeNameList[genomeIndex];
	const subject_chrInfo = genome_info_list[genomeIndex].chr_list;
	const subject_chr_name = subject_chrInfo[nChr - 1].chr;

	const subject_fa_filename = subject_chrInfo[nChr - 1].path;

	const result_text = await exec_blastn_Ex(rDNA_filePath, subject_fa_filename, undefined, undefined, undefined, undefined, "-evalue 1e-5");
	const _table = parseBlastnResults(result_text);
	
	const max_len = Math.max(..._table.map(a => a.slen));
	
	/**
	 * TODO: Check Ribosomal DNA structure, IGS 18S ITS 5.8S ITS 26S IGS
	 * CBS1-1 min:7364, max:7835, 7835 / 7364 = 0.94
	 * range include IGS
	 */
	const group = _table.filter(a => (a.slen / max_len) >= 0.9).sort((a, b) => a.s_min - b.s_max);
	if (!group.every(a => a.strand == group[0].strand)) {
		console.warn("found reverse rDNA");
	}
	
	const near_bp = 1;
	const region_list = _table.filter(aln => {
		const [ts, te] = [aln.sstart, aln.send].sort((a, b) => a - b);
		return group.some(rep => {
			const [rs, re] = [rep.sstart, rep.send].sort((a, b) => a - b);
			if (ts <= (re + near_bp) &&
				te >= (rs - near_bp)
			) {
				return true;
			}
		});
	})
	
	if (!DEBUG) {
		fs.writeFileSync(
			`${dataset.tmp_path}/blastn_rDNA_${subject_chr_name}.txt`,
			group.map(a => a.toArray().join("\t")).join("\n")
		);
	}

	let min_sstart = Math.min(...region_list.map(a => a.s_min));
	let max_send = Math.max(...region_list.map(a => a.s_max)) + 1;

	if (VERBOSE) {
		console.log({
			subject: subject_genome_name,
			min_sstart, max_send,
			len: max_send - min_sstart,
			"_table.length": _table.length,
			"group.length": group.length,
		});
	}

	if (!(max_send > min_sstart)) {
		throw new Error("if (!(max_send > min_sstart)) {");
	}
	else {
		// let raw_fa = subject_chrInfo[nChr - 1].loadSeq();
		// let info = {// save rDNA repeat position
		// 	region: [min_sstart, max_send],
		// 	repeats: group.map(a => [a.sstart, a.send]),
		// 	alignment_repeats: null,
		// };
		// return {
		// 	strand: group[0].strand,
		// 	// chrSeq: raw_fa,
		// 	//seq: getSeq(raw_fa, min_sstart, max_send),
		// 	// start: min_sstart, // rDNA region start
		// 	// end: max_send,     // rDNA region end
		// 	info,
		// };

		const data = new rDNA_Data();
		data.region = [min_sstart, max_send];
		data.repeats = group.map(a => [a.sstart, a.send]);
		data.strand = group[0].strand;
		data.nChr = nChr;
		return data;
	}
}

