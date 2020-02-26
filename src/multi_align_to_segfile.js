// @ts-check

const fs = require("fs");

const { argv_parse } = require("./util.js");
const { Dataset } = require("./dataset.js");
const { Crossover } = require("./crossover_util.js");
const { readFasta } = require("./fasta_util.js");
const { SegRow } = require("./SegFile.js");

const { loadCmpData } = require("./analyser.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const argv_input = String(argv["-i"] || "");
const argv_output_chr = Number(argv["-chr"]) | 0;
const argv_output_prefix = String(argv["--output-prefix"] || "");


if (!argv_dataset_path || !argv_input || !argv_output_chr || !argv_output_prefix) {
	console.log("Usage: ", "node multi_align_to_segfile.js -dataset <dataset.json> -i <multi_align.fa> -chr <nChr> --output-prefix <output file prefix name>");
	console.log({
		argv_dataset_path,
		argv_input,
		argv_output_chr,
		argv_output_prefix,
	})
	throw new Error("agv");
}

const dataset = Dataset.loadFromFile(argv_dataset_path);
const input_fasta = readFasta(argv_input);


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

	let analyser_options = {
		get nChr() { return argv_output_chr },
		get rDNA_info() { return dataset.rDNA_info; },
		get co_list() { return null; },//not in this stage
		get fill_prev_color() { return true; },
	};
	let data = loadCmpData(seq_list, analyser_options);

	function clear_bits(input) {
		let bits = input & data.ColorID.mask;
		switch (bits) {
			case data.ColorID.dad:
			case data.ColorID.dad_rip:
				return 0;
			case data.ColorID.mom:
			case data.ColorID.mom_rip:
				return 1;
			default:
				return -1;//undefine result
		}
	}

	let segfile = data.seg_snp.map(cmp => {
		const r1 = seq_list[0][cmp.pos];
		const r2 = seq_list[1][cmp.pos];

		if (!cmp.is_rip && r1 == r2) {
			return null;
		}

		const s1 = seq_list[2][cmp.pos];
		const s2 = seq_list[3][cmp.pos];
		const s3 = seq_list[4][cmp.pos];
		const s4 = seq_list[5][cmp.pos];
		
		//Ch	pos(bp)	A	B	C	D	type	indel	Q	C	109	111	115	119	rip	gc	
		let row = new SegRow();
		row.chr = argv_output_chr;
		row.pos = data.pos_ref1_uint32array[cmp.pos];//cmp.ref1_pos;//

		row.a = clear_bits(cmp[2]);
		row.b = clear_bits(cmp[3]);
		row.c = clear_bits(cmp[4]);
		row.d = clear_bits(cmp[5]);
		if (row.a < 0 || row.b < 0 || row.c < 0 || row.d < 0) {
			return null;
		}

		const spore_indel_cnt = [s1, s2, s3, s4].filter(a => a == "-").length;

		/** @type {["indel=1", "indel=2", "indel=3", "indel=4"]} */
		const indel_type = ["indel=1", "indel=2", "indel=3", "indel=4"];

		if (row.a != 2 && row.b != 2 && row.c != 2 && row.d != 2) {
			let tv = Number(row.a) + Number(row.b) + Number(row.c) + Number(row.d);
			row.type = (tv == 3 || tv == 1) ? "type=3" : ((tv == 0 || tv == 4) ? "type=4" : "type=2");
		}
		else {
			row.type = "";
		}

		if (r1 == "-" || r2 == "-" || spore_indel_cnt >= 1) {
			if (r1 == "-") {
				row.indel =  "indel=1|" + spore_indel_cnt;
			}
			else if (r2 == "-") {
				row.indel =  "indel=2|" + spore_indel_cnt;
			}
			else {
				row.indel =  "indel=" + spore_indel_cnt;
			}
		}
		else {
			row.indel =  "align";
		}
		row.r1 = r1;
		row.r2 = r2;
		row.s1 = s1;
		row.s2 = s2;
		row.s3 = s3;
		row.s4 = s4;
		
		row.gc = dataset.getGCByPos(dataset.ref, argv_output_chr, row.pos);
		if (dataset.isInCentromere(argv_output_chr, row.pos)) {
			row.ct = "centromere";
		}
		else if (dataset.isInTelomere(argv_output_chr, row.pos)) {
			row.ct = "telomere";
		}
		
		return row;
	}).filter(a => a);

	save_seg(`output/${argv_output_prefix}.txt`, segfile);
	//save_seg(`output/rip_${argv_output_prefix}.txt`, rip_list);
}

function save_seg(output_name, segfile) {
	const output_header = SegRow.keys();

	let final_text = output_header.join("\t") + "\n" + segfile.map(row => {
		return output_header.map(key => row[key]).join("\t");
	}).join("\r\n");

	fs.writeFileSync(output_name, final_text);
}
