// @ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, loadChrLength, array_groupBy } = require("./util.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const r1_name = dataset.ref;
const r2_name = dataset.parental_list[1];
const s1_name = dataset.progeny_list[0];
const s2_name = dataset.progeny_list[1];
const s3_name = dataset.progeny_list[2];
const s4_name = dataset.progeny_list[3];

const genome_name_list = [
	r1_name,
	r2_name,
	s1_name,
	s2_name,
	s3_name,
	s4_name,
];

if (!fs.existsSync("output/")) {
	fs.mkdirSync("output/");
}

{
	const header = ["Index#", "Chr#", "Length"].join("\t") + "\n";

	genome_name_list.forEach(genome_name => {
		let fname = `output/${genome_name}.length.txt`;
		if (!fs.existsSync(fname)) {
			let fa = readFasta(dataset.genomeFileMap[genome_name]);

			let out_text = "";
			out_text += header;
			out_text += Object.keys(fa).map((seq_name, idx) => [idx + 1, seq_name, fa[seq_name].length].join("\t")).join("\n");

			fs.writeFileSync(fname, out_text);

			console.log("output *.length.txt", fname);
		}
	});
}

{
	const r1_chr_list = loadChrLength(`output/${r1_name}.length.txt`).list;
	const r2_chr_list = loadChrLength(`output/${r2_name}.length.txt`).list;
	const s1_chr_list = loadChrLength(`output/${s1_name}.length.txt`).list;
	const s2_chr_list = loadChrLength(`output/${s2_name}.length.txt`).list;
	const s3_chr_list = loadChrLength(`output/${s3_name}.length.txt`).list;
	const s4_chr_list = loadChrLength(`output/${s4_name}.length.txt`).list;

	const genome_list = [
		[r1_name, r1_chr_list, {}],
		[r2_name, r2_chr_list, {}],
		[s1_name, s1_chr_list, {}],
		[s2_name, s2_chr_list, {}],
		[s3_name, s3_chr_list, {}],
		[s4_name, s4_chr_list, {}]
	];

	if (!fs.existsSync(`tmp/`)) {
		fs.mkdirSync(`tmp/`);
	}
	if (!fs.existsSync(`tmp/fasta`)) {
		fs.mkdirSync(`tmp/fasta`);
	}

	genome_list.forEach(([genome_name, chr_list, fa]) => {
		Object.assign(fa, readFasta(dataset.genomeFileMap[genome_name]));
	});

	for (let nChr = 1; nChr <= r1_chr_list.length; ++nChr) {
		genome_list.forEach(([genome_name, chr_list, fa]) => {
			const chr_idx = nChr - 1;
			
			let chr_info = chr_list[chr_idx];
			
			let chr_name = chr_info.chr;
			
			let fname_fasta = `tmp/fasta/${genome_name}_${chr_name}.fa`;
			
			if (!fs.existsSync(fname_fasta)) {
				console.log("output", genome_name, nChr, fname_fasta);
				saveFasta(fname_fasta, {
					[chr_name]: fa[chr_name],
				});
			}
		});
	}
}

