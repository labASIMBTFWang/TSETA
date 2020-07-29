// @ts-check

const fs = require("fs");
const Path = require("path");

const { Dataset } = require("./dataset.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { loadFragIdList, MyCoord } = require("./load_frag_list.js");

const VERBOSE = process.argv.some(a => a == "--verbose");

const mafft_output_directory = `${dataset.tmp_path}/mafft_seq_frag`;

/**
 * @param {number} nChr
 * @param {MyCoord[]} coord_list
 * @param {string} output_name
 */
function join_chr_frag(nChr, coord_list, output_name, options = { padding: false, chr_list: null }) {
	/** @type {{ [seqName:string]:string }} */
	let output_fa = {
	};
	
	coord_list.forEach(coord => {
		let in_filename = `${mafft_output_directory}/mafft_ch${nChr}_${coord.id}.fa`;

		if (!coord.centromere && fs.existsSync(in_filename)) {
			let in_fa = readFasta(in_filename);

			if (VERBOSE) {
				console.log("ch", nChr, "frag", coord.id);
				console.log("in", in_filename);
			}

			const max_len = Math.max(...Object.values(in_fa).map(a => a ? a.length : 0));

			(options.chr_list || Object.keys(in_fa)).forEach(seq_name => {
				let seq = in_fa[seq_name];

				if (options.padding) {
					if (!seq) {
						seq = "";
					}

					const ori_len = seq.length;
					seq = seq.padEnd(max_len, "-");
					if (seq.length != ori_len) {
						if (VERBOSE) {
							console.log(in_filename, seq_name, "padEnd", ori_len, "->", seq.length);
						}
					}
				}

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
					
					console.log("in ref1", in_ref1_filename);
					console.log("in ref2", in_ref2_filename);
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
				console.log("no file:", coord);
				//break;
			}
		}
	});

	saveFasta(output_name, output_fa);

	return output_fa;
}

module.exports.join_chr_frag = join_chr_frag;

