//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { Dataset } = require("./dataset.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { loadFragIdList, MyCoord } = require("./load_frag_list.js");
const { validation_chr } = require("./validation_seq.js");
const { join_chr_frag } = require("./join_chr_frag.js");

const argv = argv_parse(process.argv);

const VERBOSE = !!argv["--verbose"];

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

let mafft_output_directory = String(argv["-i"] || `${dataset.tmp_path}/mafft_seq_frag`);

const genome_info_list = dataset.loadGenomeInfoList();

if (fs.realpathSync(process.argv[1]) == __filename) {
	main();
}
else {
	debugger;
}

function main() {
	let genome_has_error = false;
	merge_chr_all_fa();
	
	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		try {
			let chr_has_error = validation_chr(nChr, dataset.tmp_path, false);
			genome_has_error = genome_has_error || chr_has_error;
		}
		catch (ex) {
			console.error(ex);
		}
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

