//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, loadChrLength, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { Dataset } = require("./dataset.js");
const { BlastnCoord, execAsync, exec_blastn, parseBlastnResults, blastn_coord, isCollide, groupByOverlap } = require("./blastn_util.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { loadFragIdList } = require("./load_frag_list.js");

const argv = argv_parse(process.argv);

const argv_flags = String(argv["-flags"] || "") == "raw" ? "raw" : "";

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const r1_name = dataset.ref;
const r2_name = dataset.parental_list[1];
const s1_name = dataset.progeny_list[0];
const s2_name = dataset.progeny_list[1];
const s3_name = dataset.progeny_list[2];
const s4_name = dataset.progeny_list[3];

const r1_chr_list = loadChrLength(`output/${r1_name}.length.txt`).list;
const r2_chr_list = loadChrLength(`output/${r2_name}.length.txt`).list;
const s1_chr_list = loadChrLength(`output/${s1_name}.length.txt`).list;
const s2_chr_list = loadChrLength(`output/${s2_name}.length.txt`).list;
const s3_chr_list = loadChrLength(`output/${s3_name}.length.txt`).list;
const s4_chr_list = loadChrLength(`output/${s4_name}.length.txt`).list;

const genome_name_list = [
	r1_name,
	r2_name,
	s1_name,
	s2_name,
	s3_name,
	s4_name,
];
const genome_length_list = [
	r1_chr_list,
	r2_chr_list,
	s1_chr_list,
	s2_chr_list,
	s3_chr_list,
	s4_chr_list,
];

main();

function main() {
	validation_frag({ raw: argv_flags == "raw" });
}

function validation_frag(flags = { raw: true }) {
	const merge_centromere = !flags.raw;
	const all_chr_frag_list = loadFragIdList(dataset, merge_centromere);

	const input_path_dir = flags.raw ? "tmp/seq_frag" : "tmp/mafft_seq_frag";
	
	/** @type {{ [chrName: string]: string }[]} */
	let srcRaw_genome_list = [];
	genome_name_list.forEach((genome_name, i) => {
		let genome_seq = readFasta(dataset.genomeFileMap[genome_name]);
		srcRaw_genome_list[i] = genome_seq;
	});

	for (let nChr = 1; nChr <= dataset.chrNames.length; ++nChr) {
		/** @type {{ [chrName: string]: string }} */
		const output_src_fa = {};

		const chr_name_list = [];
		/** @type {{ [chrName: string]: string }} */
		const src_raw_seq = {};
		srcRaw_genome_list.forEach((src_raw_genome, i) => {
			const idxChr = nChr - 1;
			const chr_name = genome_length_list[i][idxChr].chr;
			src_raw_seq[chr_name] = src_raw_genome[chr_name];
			chr_name_list.push(chr_name);
			output_src_fa[chr_name] = "";//init
		});
		
		if (all_chr_frag_list[nChr]) {
			try {
				const final_list = [];

				for (let coord of all_chr_frag_list[nChr]) {
					const fragId = coord.id;
					const in_filename = `${input_path_dir}/${flags.raw ? "" : "mafft_"}ch${nChr}_${fragId}.fa`;
			
					/** @type {{ [chrName: string]: string }} */
					let in_fa;

					if (!coord.centromere && fs.existsSync(in_filename)) {
						if (fs.statSync(in_filename).size) {
							in_fa = readFasta(in_filename);
							final_list.push(in_filename);
						}
						else {
							console.log("invalid fasta", "ch", nChr, "frag", fragId, Path.resolve(in_filename));
							throw new Error("invalid fasta");
						}
						//console.log("ch", nChr, "frag", i);
					}
					else {
						const in_ref1_filename = `${input_path_dir}/${flags.raw ? "" : "mafft_"}ch${nChr}_${fragId}_ref1.fa`;
						const in_ref2_filename = `${input_path_dir}/${flags.raw ? "" : "mafft_"}ch${nChr}_${fragId}_ref2.fa`;
			
						if (fs.existsSync(in_ref1_filename) && fs.existsSync(in_ref2_filename)) {
							if (!fs.statSync(in_ref1_filename).size) {
								console.log("invalid fasta", "ch", nChr, "frag", fragId, "ref1", Path.resolve(in_ref1_filename));
								throw new Error("invalid fasta ref1");
							}
							if (!fs.statSync(in_ref2_filename).size) {
								console.log("invalid fasta", "ch", nChr, "frag", fragId, "ref2", Path.resolve(in_ref2_filename));
								throw new Error("invalid fasta ref2");
							}
							const in_ref1_fa = readFasta(in_ref1_filename);
							const in_ref2_fa = readFasta(in_ref2_filename);
							
							final_list.push([in_ref1_filename, in_ref2_filename]);
							
							//console.log("ch", nChr, "frag", i, "ref1 ref2");

							in_fa = Object.assign({}, in_ref1_fa, in_ref2_fa);
						}
						else {
							console.log({
								in_filename, in_ref1_filename, in_ref2_filename,
							});
							throw new Error("ref1 ?? ref2 ??");
						}
					}

					chr_name_list.forEach(chr_seq_name => {
						/** @type {string} */
						const _seq = in_fa[chr_seq_name];
						if (_seq) {
							const seq = flags.raw ? _seq : _seq.replace(/-/g, "");

							const raw_start = output_src_fa[chr_seq_name].length;
							const raw_end = raw_start + seq.length;
							const raw_frag = src_raw_seq[chr_seq_name].slice(raw_start, raw_end);
							for (let ii = 0; ii < seq.length; ++ii) {
								if (seq[ii] != raw_frag[ii]) {
									console.error("error at", in_filename, ii);
									console.error("src", raw_frag.slice(ii - 10, ii + 10));
									console.error("res", seq.slice(ii - 10, ii + 10));
									console.error(">>>", "----------^");
									throw new Error("invalid frag");
								}
							}
							
							output_src_fa[chr_seq_name] += seq;
						}
					});
				}

				let va_result = chr_name_list.every(chr_seq_name => {
					for (let ii = 0; ii < src_raw_seq[chr_seq_name].length; ++ii) {
						if (output_src_fa[chr_seq_name][ii] != src_raw_seq[chr_seq_name][ii]) {
							console.error("error at", ii);
							console.error("src", src_raw_seq[chr_seq_name].slice(ii - 10, ii + 10));
							console.error("res", output_src_fa[chr_seq_name].slice(ii - 10, ii + 10));
							console.error(">>>", "----------^");
							throw new Error("invalid full length");
						}
					}
				});
				
				fs.writeFileSync(`tmp/mafft_ch${nChr}.json`, JSON.stringify(final_list, null, "\t"));

				console.log("ch", nChr, "ok");
			}
			catch (ex) {
				console.error(ex);
			}
		}//if (all_chr_frag_list[nChr]) {
		else {
			console.log("skip", "ch", nChr);
		}
	}
}

