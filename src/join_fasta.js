//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { Dataset } = require("./dataset.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { loadFragIdList, MyCoord } = require("./load_frag_list.js");
const { validation_chr } = require("./validation_seq.js");

const argv = argv_parse(process.argv);

const VERBOSE = !!argv["--verbose"];

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

let mafft_output_directory = String(argv["-i"] || `${dataset.tmp_path}/mafft_seq_frag`);

const genome_info_list = dataset.loadGenomeInfoList();

if (process.argv[1] == __filename) {
	main();
}
else {
	debugger;
}

function main() {
	let genome_has_error = false;
	merge_chr_all_fa();
	
	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		let chr_has_error = validation_chr(nChr, dataset.tmp_path, false);
		genome_has_error = genome_has_error || chr_has_error;
	}

	if (genome_has_error) {
		console.log("has error, no output");
		return;
	}
	else {
		//clone all results to output path
		for (let i = 1; i <= genome_info_list[0].chr_list.length; ++i) {
			let input_path = `${dataset.tmp_path}/mafft_ch${i}.fa`;
			let output_path = `${dataset.output_path}/mafft_ch${i}.fa`;
			fs.createReadStream(input_path).pipe(fs.createWriteStream(output_path));//copy file
			console.log("output:", output_path);
		}
	}
	
	console.log("next step:", "define and align the 50S rDNA loci");
	console.log("command:", `node ./src/re_align_rDNA.js -dataset ${argv_dataset_path} -rDNA rDNA.fa -chr 6`);
}

function merge_chr_all_fa() {
	const all_chr_frag_list = loadFragIdList(dataset);//load_frag_id_list();

	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		if (all_chr_frag_list[nChr]) {
			//let id_list = all_chr_frag_list[nChr].map(a => a.id);
	
			join_chr_frag(nChr, all_chr_frag_list[nChr], `${dataset.tmp_path}/mafft_ch${nChr}.fa`);
		}
		else {
			console.log("skip", "ch", nChr);
		}
	}
}

/**
 * @param {number} nChr
 * @param {MyCoord[]} coord_list
 * @param {string} output_name
 */
function join_chr_frag(nChr, coord_list, output_name) {
	/** @type {{ [seqName:string]:string }} */
	let output_fa = {
	};
	
	coord_list.forEach(coord => {
		let in_filename = `${mafft_output_directory}/mafft_ch${nChr}_${coord.id}.fa`;

		if (!coord.centromere && fs.existsSync(in_filename)) {
			let in_fa = readFasta(in_filename);

			if (VERBOSE) {
				console.log("ch", nChr, "frag", coord.id);
			}

			Object.keys(in_fa).forEach(seq_name => {
				let seq = in_fa[seq_name];
				if (seq) {
					if (!output_fa[seq_name]) {
						output_fa[seq_name] = "";
					}
					output_fa[seq_name] += seq;
					
					if (VERBOSE) {
						console.log("seq.len", seq_name, output_fa[seq_name].length, "+", seq.length);
					}
				}
				else {
					console.log("skip seq:", coord.id, seq_name);
				}
			});
		}
		else {
			let in_ref1_filename = `${mafft_output_directory}/mafft_ch${nChr}_${coord.id}_ref1.fa`;
			let in_ref2_filename = `${mafft_output_directory}/mafft_ch${nChr}_${coord.id}_ref2.fa`;

			if (fs.existsSync(in_ref1_filename) && fs.existsSync(in_ref2_filename)) {
				let in_ref1_fa = readFasta(in_ref1_filename);
				let in_ref2_fa = readFasta(in_ref2_filename);

				if (VERBOSE) {
					console.log("ch", nChr, "frag", coord.id, "ref1 ref2");
				}

				let ref1_max_length = Math.max(...Object.values(in_ref1_fa).map(seq => seq.length));
				let ref2_max_length = Math.max(...Object.values(in_ref2_fa).map(seq => seq.length));
				let multi_length = ref1_max_length + ref2_max_length;
				
				Object.keys(in_ref1_fa).forEach(seq_name => {
					in_ref1_fa[seq_name] = in_ref1_fa[seq_name].padEnd(ref1_max_length, "-");
				});
				Object.keys(in_ref2_fa).forEach(seq_name => {
					in_ref2_fa[seq_name] = in_ref2_fa[seq_name].padEnd(ref2_max_length, "-");
				});

				if (ref1_max_length >= ref2_max_length) {
					Object.keys(in_ref1_fa).forEach(seq_name => {
						in_ref1_fa[seq_name] = in_ref1_fa[seq_name].padStart(multi_length, "-");
					});
					Object.keys(in_ref2_fa).forEach(seq_name => {
						in_ref2_fa[seq_name] = in_ref2_fa[seq_name].padEnd(multi_length, "-");
					});
				}
				else {
					Object.keys(in_ref1_fa).forEach(seq_name => {
						in_ref1_fa[seq_name] = in_ref1_fa[seq_name].padEnd(multi_length, "-");
					});
					Object.keys(in_ref2_fa).forEach(seq_name => {
						in_ref2_fa[seq_name] = in_ref2_fa[seq_name].padStart(multi_length, "-");
					});
				}
				
				let in_fa = Object.assign({}, in_ref1_fa, in_ref2_fa);
				
				Object.keys(in_fa).forEach(seq_name => {
					let seq = in_fa[seq_name];
					if (seq) {
						if (!output_fa[seq_name]) {
							output_fa[seq_name] = "";
						}
						output_fa[seq_name] += seq;
						if (VERBOSE) {
							console.log("indel seq.len", seq_name, output_fa[seq_name].length, "+", seq.length);//translocation, ...
						}
					}
					else {
						console.log("skip seq:", coord.id, seq_name);
					}
				});
			}
			else {
				console.log("no file:", coord.id);
				//break;
			}
		}
	});

	saveFasta(output_name, output_fa);
}
