
// @ts-check

const fs = require("fs");
const Path = require("path");
const child_process = require("child_process");

const { argv_parse, loadChrLength } = require("./util.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);

const argv_input_path = String(argv["-i"] || "");
const argv_dataset_patha = String(argv["-dataset"] || "");

const argv_find_border = String(argv["--find-border"] || "");

const dataset = Dataset.loadFromFile(argv_dataset_patha);


class ErrInfo {
	constructor() {
		this.start = -1;
		this.end = -1;

		this.min_raw = null;
		this.min_new = null;

		this.full_raw = null;
		this.full_new = null;
		this.old_fa = null;
	}
}

main();

function main() {
	for (let nChr = 1; nChr <= dataset.chrNames.length; ++nChr) {
		const input_path = `mafft_ch${nChr}.fa`;
		const input_fasta = readFasta(Path.join(argv_input_path, input_path));

		console.log("input_path", input_path);
		
		const seq_id_list = Object.keys(input_fasta);

		validation_seq(seq_id_list, input_fasta, nChr);

		if (argv_find_border) {
			find_border(seq_id_list, input_fasta, nChr);
		}
	}
}

function validation_seq(seq_id_list, fa, nChr) {
	let genome_name_list = [
		dataset.ref,
		dataset.refs[1],
		dataset.progeny_list[0],
		dataset.progeny_list[1],
		dataset.progeny_list[2],
		dataset.progeny_list[3],
	];
	
	let raw_seq = genome_name_list.map((genome_name, i) => {
		let chr_seq = readFasta(Path.resolve("./", `${genome_name}.genome.fa`));
		return chr_seq[seq_id_list[i]];
	});
	
	seq_id_list.forEach((name, idx) => {
		let s_name = seq_id_list[idx];
		let has_error = false;
		let to_raw_pos = 0;
		for (let i = 0; i < fa[name].length; ++i) {
			if (fa[name][i] != "-") {
				if (fa[name][i] == raw_seq[idx][to_raw_pos]) {
				}
				else {
					const ext = 10;
					console.log(s_name, "error at", i);
					console.log("in ", fa[name].slice(i - ext, i + 20));
					console.log("raw", raw_seq[idx].slice(to_raw_pos - ext, to_raw_pos + 20));
					let mark = [...fa[name].slice(i - ext, i)].fill(".").join("") + "^";
					console.log("mrk", mark);

					has_error = true;
					break;
				}
				++to_raw_pos;
			}
		}
		if (!has_error) {
			console.log(s_name, "ok");
		}
	});
}

function find_border(seq_id_list, fa, nChr) {
	console.log("validation seq");

	let seq_list = [
		"", "",
		"", "", "", "",
	];
	seq_id_list.forEach((id, i) => seq_list[i] = fa[id]);

	if (!fs.existsSync("err-ma")) {
		fs.mkdirSync("err-ma");
	}

	/** @type {ErrInfo[]} */
	let err_info = [];

	let count = 0;

	let ref1_pos = 1;
	for (let i = 0; i < seq_list[0].length - 1; ++i) {
		if (!dataset.isIn_rDNA(nChr, i)) {
			if (seq_list[0][i] != "-") {
				check(ref1_pos);
				++ref1_pos;
			}
			else if (ref1_pos - 1) {//indel
				check(ref1_pos - 1);
			}
		}
		function check(ref1_pos) {
			let valid = seq_list.length;
			let valid_2 = seq_list.length;
			if (nChr == 6 && i == 2086459) {
				console.log({
					chr: nChr,
					ma_pos: i,
					ref1_pos,
				});
				debugger;
			}
			if (dataset.isInCentromere(nChr, ref1_pos)) {
				//skip centromere
			}
			else {
				for (let j = 0; j < seq_list.length; ++j) {
					let current = seq_list[j][i];
					let next = seq_list[j][i + 1];
					if (current == "-" && next != "-") {
						--valid;
					}
					else if (current != "-" && next == "-") {
						--valid;
					}
					if (current == "-" && next == "-") {
						--valid_2;
					}
				}
				if (valid == 0 || (valid < seq_list.length && (seq_list.length - valid) == valid_2)) {
					console.log("error at:", i);
					let max_lower = 0;
					let max_upper = 0;
					for (let j = 0; j < seq_list.length; ++j) {
						for (let k = i - 1; k >= 0; --k) {
							let current = seq_list[j][k];
							if (current != "-") {
								max_lower = Math.max(max_lower, i - k);
								break;
							}
						}
						for (let k = i + 1; k < seq_list[j].length; ++k) {
							let current = seq_list[j][k];
							if (current != "-") {
								max_upper = Math.max(max_upper, k - i);
								break;
							}
						}
					}
					const show_more = 0;
					const end_200207 = 1;
					let min = Math.min(max_lower, max_upper);
					let max = Math.max(max_lower, max_upper);
					let err_start = i - max_lower - show_more + 1;
					let err_end = i + max_upper + show_more + end_200207;
					let err_min_start = i - min - show_more + 1;
					let err_min_end = i + min + show_more + end_200207;
					let err_min_fa = {};
					let err_fa = {};
					let old_fa = {};
					for (let j = 0; j < seq_list.length; ++j) {
						if (Math.abs(max_upper - max_lower) < 100) {
							let err_seq = seq_list[j].slice(err_start, err_end);
							console.log("seq" + j, max_lower, max_lower < max_upper ? "<" : ">", max_upper, min, err_seq);
							old_fa[seq_id_list[j]] = err_seq;
						}
						else {
							let err_seq;
							if (min < 100) {
								err_seq = seq_list[j].slice(i - min - show_more + 1, i + min + show_more);
								console.log("seq", max_lower, max_lower < max_upper ? "<" : ">", max_upper, min, err_seq, "...more", Math.abs(max_upper - max_lower));
								old_fa[seq_id_list[j]] = err_seq;
							}
							else {
								console.log("seq", max_lower, max_lower < max_upper ? "<" : ">", max_upper, min, "...more", Math.abs(max_upper - max_lower));
							}
						}
						err_min_fa[seq_id_list[j]] = seq_list[j].slice(err_min_start, err_min_end).replace(/-/g, "");
						err_fa[seq_id_list[j]] = seq_list[j].slice(err_start, err_end).replace(/-/g, "");
					}
					
					let ei = new ErrInfo();
					ei.start = err_min_start;
					ei.end = err_min_end;
					
					// let all_match = seq_id_list.slice(1).every(sname => err_min_fa[seq_id_list[0]] == err_min_fa[sname]);
					// if (all_match) {
					// }
					// else {
						const min_ma_name = `err-ch${nChr}-min-${err_min_start}-${err_min_end}.fa`;
						const min_ma_path = `err-ma/${min_ma_name}`;
						const min_output_path = `err-ma/mafft-${min_ma_name}`;
						
						const full_ma_name = `err-ch${nChr}-full-${err_start}-${err_end}.fa`;
						const full_ma_path = `err-ma/${full_ma_name}`;
						const full_output_path = `err-ma/mafft-${full_ma_name}`;

						ei.min_raw = err_min_fa;
						ei.full_raw = err_fa;
						ei.old_fa = old_fa;

						// saveFasta(min_ma_path, err_min_fa);
						// saveFasta(full_ma_path, err_fa);

						// let min_mafft_cmd = `mafft --quiet --thread ${num_thread} ${algorithm} --maxiterate ${maxiterate} ${min_ma_path} > ${min_output_path}`;
						// let full_mafft_cmd = `mafft --quiet --thread ${num_thread} ${algorithm} --maxiterate ${maxiterate} ${full_ma_path} > ${full_output_path}`;
						// try {
						// 	if (min < 100) {
						// 		child_process.execSync(min_mafft_cmd);
						// 		let ma_min_new = readFasta(min_output_path);
						// 		ei.min_new = ma_min_new;
						// 	}
						// 	else {
						// 		ei.min_new = "...more";
						// 	}
						// 	if (max < 100) {
						// 		child_process.execSync(full_mafft_cmd);
						// 		let ma_full_new = readFasta(full_output_path);
						// 		ei.full_new = ma_full_new;
						// 	}
						// 	else {
						// 		ei.full_new = "...more";
						// 	}
						// }
						// catch (ex) {
						// 	console.error(ex);
						// }
					// }

					err_info.push(ei);

					++count;
				}
			}
		}
	}

	//fs.writeFileSync(`err-ch${nChr}.json`, JSON.stringify(err_info, null, "\t"));

	console.log("count", count);
}
