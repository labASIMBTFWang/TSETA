// @ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, array_groupBy } = require("./util.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

if (!fs.existsSync("output/")) {
	fs.mkdirSync("output/");
}

{
	const header = ["Index", "Chromosome", "Length"].join("\t") + "\n";

	dataset.genomeNameList.forEach(genome_name => {
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
	const genome_info_list = dataset.loadGenomeInfoList();

	if (!fs.existsSync(`tmp/`)) {
		fs.mkdirSync(`tmp/`);
	}
	if (!fs.existsSync(`tmp/fasta`)) {
		fs.mkdirSync(`tmp/fasta`);
	}
	
	genome_info_list.forEach(genome_info => {
		Object.assign(genome_info.fasta, readFasta(dataset.genomeFileMap[genome_info.name]));
	});

	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		genome_info_list.forEach(genome_info => {
			const chr_idx = nChr - 1;
			
			let chr_info = genome_info.chr_list[chr_idx];
			
			let chr_name = chr_info.chr;
			
			let fname_fasta = `tmp/fasta/${genome_info.name}_${chr_name}.fa`;
			
			if (!fs.existsSync(fname_fasta)) {
				console.log("output", genome_info.name, nChr, fname_fasta);
				saveFasta(fname_fasta, {
					[chr_name]: genome_info.fasta[chr_name],
				});
			}
		});
	}
}

