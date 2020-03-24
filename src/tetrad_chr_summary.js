
// @ts-check

const fs = require("fs");

const { argv_parse } = require("./util.js");
const { Dataset } = require("./dataset.js");
const { Crossover } = require("./crossover_util.js");
const { readFasta, multialign_to_chrPos_posMap } = require("./fasta_util.js");
const { SegRow } = require("./SegFile.js");

const { initAnalyser, loadCmpData } = require("./analyser.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const argv_output_chr = Number(argv["-chr"]) | 0;
const argv_seq = String(argv["--seq"] || "");
const argv_co_list = String(argv["--co-list"] || "");
const argv_output_prefix = String(argv["--output-prefix"] || "");


if (!argv_dataset_path || !argv_seq || !argv_output_chr || !argv_output_prefix || !argv_co_list) {
	console.log("Usage: ", "node multi_align_to_segfile.js -dataset <dataset.json> -chr <nChr> --seq <multi_align.fa> --co-list <chrN_co_list.json> --output-prefix <output file prefix name>");
	console.log({
		argv_dataset_path,
		argv_seq,
		argv_output_chr,
		argv_output_prefix,
	})
	throw new Error("agv");
}

const co_list = ((co_list_filename) => {
	const text = fs.readFileSync(co_list_filename).toString();
	const list = JSON.parse(text);
	list.forEach(co => {
		co.before = co.before.split(",").map(n => Number(n));
		co.after = co.after.split(",").map(n => Number(n));
	});
	return list;
})(argv_co_list);

const input_fasta = readFasta(argv_seq);
const dataset = Dataset.loadFromFile(argv_dataset_path);
dataset.load_GC_content();
dataset.load_rDNA_info();

main();

function main() {
	console.log({
		ref: dataset.ref,
		chr: argv_output_chr
	});

	/** @type {string[]} */
	let seq_list = [];
	let seq_id_list = Object.keys(input_fasta);
	seq_id_list.forEach((id, i) => seq_list[i] = input_fasta[id]);

	let analysis_options = {
		get mode() { return dataset.mode; },
		get nChr() { return argv_output_chr },
		get rDNA_info() { return dataset.rDNA_info; },
		get show_rDNA_snp() { return false; },
		get co_list() { return co_list; },//in this stage
		get fill_prev_color() { return true; },
		get mode() { return dataset.mode; },
	};
	initAnalyser(analysis_options);
	let data = loadCmpData(seq_list, analysis_options);

	let ref_snp = [...seq_list[0]].map((q, i) => [data.pos_ref1_uint32array[i], q, seq_list[1][i], data.pos_ref2_uint32array[i]]).filter(([ref1_pos, q, c, ref2_pos]) => q != c);
	let ref_snv = ref_snp.filter(([ref1_pos, q, c, ref2_pos]) => q != c && q != "-" && c != "-");
	let ref_snp_indel = ref_snp.filter(([ref1_pos, q, c, ref2_pos]) => q != c && ((q == "-" && c != "-") || (q != "-" && c == "-")));//q is del or c is del

	fs.writeFileSync(`${dataset.output_path}/${argv_output_prefix}_snp.txt`, JSON.stringify(ref_snp));
	fs.writeFileSync(`${dataset.output_path}/${argv_output_prefix}_snv.txt`, JSON.stringify(ref_snv));
	fs.writeFileSync(`${dataset.output_path}/${argv_output_prefix}_snp_indel.txt`, JSON.stringify(ref_snp_indel));

	const output_head = [
		"Chromosome",
		"simple CO", "CO + NCO",
		"Q/C SNV", "Q/C SNP", "Q/C InDel",
		"RIP Q", "RIP C",
		"illegitimate mutation",
		"SNP 2:2", "NCO 3:1", "NCO 4:0",
		"1n:3", "2n:2", "3n:1", "4n:0"
	];

	const chr = argv_output_chr;

	const simpleCO_list = co_list.filter(co => co.type == "CO");
	const CO_NCO_list = co_list.filter(co => co.type == "CO(NCO)");

	const snv = ref_snv;
	const snp = ref_snp;
	const snp_indel = ref_snp_indel;
	
	const rip_ref1 = data.rip_list.filter(rip => rip.rip_ref == 1);
	const rip_ref2 = data.rip_list.filter(rip => rip.rip_ref == 2);

	const illeg = data.spore_cmp.illegitimate_mutation_list;

	const s22 = data.spore_cmp.s22;
	const s31 = data.spore_cmp.s31;
	const s40 = data.spore_cmp.s40;

	const s1n3 = data.spore_cmp.s1n3;
	const s2n2 = data.spore_cmp.s2n2;
	const s3n1 = data.spore_cmp.s3n1;
	const s4n0 = data.spore_cmp.s4n0;

	let output_text = "";

	output_text += output_head.join("\t") + "\n";
	output_text += [
		chr,
		simpleCO_list.length, CO_NCO_list.length,
		snv.length, snp.length, snp_indel.length,
		rip_ref1.length, rip_ref2.length,
		illeg.length,
		s22.length, s31.length, s40.length,
		s1n3.length, s2n2.length, s3n1.length, s4n0.length
	].join("\t") + "\n";

	fs.writeFileSync(`${dataset.output_path}/${argv_output_prefix}_summary.txt`, output_text);
}

