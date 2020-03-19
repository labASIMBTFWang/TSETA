// @ts-check

const fs = require("fs");
const Path = require("path");

const { inputFile, inputDirectory, inputNumber, inputText, inputSelect } = require("./interact-util.js").userInput;
const { argv_parse, array_groupBy } = require("./util.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { GenomeDataSet, Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);


main();

async function main() {
	const { genome_dataset, output_filename } = await load_or_input();

	const loaded_data = make_genome_info(genome_dataset);
	explode_genome(genome_dataset, loaded_data);
	
	fs.writeFileSync(output_filename, JSON.stringify(genome_dataset, null, "\t"));

	console.log("next step:", "detect GC content %, AT-rich blocks, centromere, telomeres");
	console.log("command:", `node ./src/make_GC_AT_table.js -dataset ${output_filename} -w <window size> -m <???>`);
}

async function load_or_input() {
	const argv_dataset_path = String(argv["-dataset"] || "");

	let genome_dataset;
	try {
		if (fs.existsSync(argv_dataset_path)) {
			genome_dataset = Dataset.loadFromFile(argv_dataset_path);
		}
	}
	catch (ex) {
		genome_dataset = null;
		console.error(ex);
	}

	if (!genome_dataset) {
		genome_dataset = await userInput_genomeDataset();
	}
	
	if (!fs.existsSync(`${genome_dataset.output_path}/`)) {//make output directory
		fs.mkdirSync(`${genome_dataset.output_path}/`);
	}
	let output_filename = `${encodeURIComponent(genome_dataset.name)}.json`;

	return { genome_dataset, output_filename };
}

/**
 * @returns {Promise<GenomeDataSet>}
 */
async function userInput_genomeDataset() {
	let set_genomeName = new Set();

	let taskName = await inputText("Task name");

	/** @type {"tetrad"|"SNP"} */
	let mode = await inputSelect("mode", ["tetrad", "SNP"], "tetrad");
	
	// current version
	// genome \ mode | tetrad |  SNP
	// n parental    |    = 2 |  = 1
	// n progeny     |    = 4 | >= 1
	
	// future version
	// genome \ mode | tetrad |  SNP
	// n parental    |    = 2 | >= 1
	// n progeny     |   >= 4 | >= 1
	
	//let numParental = mode == "tetrad" ? 2 : await inputNumber("number of parental", { min: 1, max: 2, default: 1 });
	let numParental = mode == "tetrad" ? 2 : 1;
	/** @type {string[]} */
	let parental_list = [];
	/** @type {{[genomeName:string]:string}} */
	let parental_map = {};
	let refName = await inputText(`${mode == "tetrad" ? "ref (parental 1)" : "ref"} name`);
	let refFasta = await inputFile(`${mode == "tetrad" ? "ref (parental 1)" : "ref"} genome fasta file`, `${refName}.genome.fa`);
	parental_list.push(refName);
	set_genomeName.add(refName);
	parental_map[refName] = refFasta;

	for (let  parentalIdx = 1; parentalIdx < numParental; ++parentalIdx) {
		let pName = await inputText(`parental ${(parentalIdx + 1)} name`);
		if (!set_genomeName.has(pName)) {
			parental_list.push(pName);
			set_genomeName.add(pName);
			
			let filepath = await inputFile(`parental ${(parentalIdx + 1)} genome fasta file`, `${pName}.genome.fa`);
			parental_map[pName] = filepath;
		}
		else {
			console.warn(`exist ${pName}`);
		}
	}

	//let numProgeny = await inputNumber("number of progeny", { min: mode == "tetrad" ? 4 : (numParental > 1 ? 0 : 1), max: 16 });
	let numProgeny = mode == "tetrad" ? 4 : await inputNumber("number of subject", { min: 1, max: null });
	/** @type {string[]} */
	let progeny_list = [];
	/** @type {{[genomeName:string]:string}} */
	let progeny_map = {};

	for (let  progenyIdx = 0; progenyIdx < numProgeny; ++progenyIdx) {
		let pName = await inputText(`${mode == "tetrad" ? "progeny" : "subject"} ${(progenyIdx + 1)} name`);
		if (!set_genomeName.has(pName)) {
			progeny_list.push(pName);
			set_genomeName.add(pName);
			
			let filepath = await inputFile(`${mode == "tetrad" ? "progeny" : "subject"} ${(progenyIdx + 1)} fasta file`, `${pName}.genome.fa`);
			progeny_map[pName] = filepath;
		}
		else {
			console.warn(`exist ${pName}`);
		}
	}

	let genome_set = new GenomeDataSet();
	genome_set.name = taskName;
	genome_set.mode = mode;
	genome_set.ref = refName;
	genome_set.parental_list = parental_list;
	genome_set.progeny_list = progeny_list;
	genome_set.parental = parental_map;
	genome_set.progeny = progeny_map;

	console.log(genome_set);
	
	process.stdin.pause();//stop input

	return genome_set;
}

/**
 * @param {GenomeDataSet} genome_dataset
 * @returns {{[genomeName:string]:{genome_name:string,chr_name_list:string[],fasta:{[chrName:string]:string}}}}
 */
function make_genome_info(genome_dataset) {
	const header = ["Index", "Chromosome", "Length"].join("\t") + "\n";

	const genomeNameList = [...genome_dataset.parental_list, ...genome_dataset.progeny_list];
	const genomeFileMap = {
		...genome_dataset.parental,
		...genome_dataset.progeny,
	};

	/** @type {{[genomeName:string]:{genome_name:string,chr_name_list:string[],fasta:{[chrName:string]:string}}}} */
	const loaded_data = {};

	genomeNameList.forEach(genome_name => {
		let fname = `${genome_dataset.output_path}/${genome_name}.length.txt`;
		let fa = readFasta(genomeFileMap[genome_name]);

		const chr_name_list = Object.keys(fa);

		loaded_data[genome_name] = {
			genome_name: genome_name,
			fasta: fa,
			chr_name_list: chr_name_list,
		};

		if (!fs.existsSync(fname)) {
			let out_text = "";
			out_text += header;
			out_text += chr_name_list.map((seq_name, idx) => [idx + 1, seq_name, fa[seq_name].length].join("\t")).join("\n");

			fs.writeFileSync(fname, out_text);

			console.log("output *.length.txt", fname);
		}
	});

	return loaded_data;
}

/**
 * @param {GenomeDataSet} genome_dataset
 * @param {{[genomeName:string]:{genome_name:string,chr_name_list:string[],fasta:{[chrName:string]:string}}}} loaded_data
 */
function explode_genome(genome_dataset, loaded_data) {
	if (!fs.existsSync(`${genome_dataset.tmp_path}/`)) {
		fs.mkdirSync(`${genome_dataset.tmp_path}/`);
	}
	if (!fs.existsSync(`${genome_dataset.tmp_path}/fasta`)) {
		fs.mkdirSync(`${genome_dataset.tmp_path}/fasta`);
	}

	Object.keys(loaded_data).forEach(gName => {
		loaded_data[gName].chr_name_list.forEach((chrName, idx) => {
			const nChr = idx + 1;

			let fname_fasta = `${genome_dataset.tmp_path}/fasta/${gName}_${chrName}.fa`;
			
			if (!fs.existsSync(fname_fasta)) {
				console.log("output", gName, nChr, fname_fasta);
				saveFasta(fname_fasta, {
					[chrName]: loaded_data[gName].fasta[chrName],
				});
			}
		});
	});
}

