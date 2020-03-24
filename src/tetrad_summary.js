// @ts-check

const fs = require("fs");
const child_process = require("child_process");

const { tsv_parse, _table_to_object_list, table_to_object_list } = require("./tsv_parser.js");
const { argv_parse } = require("./util.js");
const { Dataset } = require("./dataset.js");
const { readFasta, saveFasta } = require("./fasta_util.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const argv_closeCOsMinDistance = Number(argv["-min-co"]);
const argv_output_summary_only = !!argv["--output-summary-only"];

const VERBOSE = !!argv["--verbose"];

const dataset = Dataset.loadFromFile(argv_dataset_path);
dataset.load_GC_content();
dataset.load_rDNA_info();

if (Number.isSafeInteger(argv_closeCOsMinDistance)) {
	dataset.crossover.closeCOsMinDistance = argv_closeCOsMinDistance;
	fs.writeFileSync(argv_dataset_path, JSON.stringify(dataset, null, "\t"));
}
else {
	console.error("Required: -min-co (number >= 0)");
	console.error("-min-co:", "The closest distance between two adjacent crossovers");
}

const genome_info_list = dataset.loadGenomeInfoList();

const output_prefix = encodeURIComponent(dataset.name);

let cmds = [];

main();

async function main() {
	if (!argv_output_summary_only) {
		let fin_list = [];
		
		//for each chromosome
		let tasks = genome_info_list[0].chr_list.map(async function (chrName, chrIdx) {
			const nChr = chrIdx + 1;
			const seq = `${dataset.output_path}/mafft_ch${nChr}.fa`;

			fin_list[chrIdx] = ["init"];

			let multi_align_to_segfile_args =  [
				"--max-old-space-size=4096",
				`${__dirname}/multi_align_to_segfile.js`,
				"-dataset", argv_dataset_path,
				"-i", seq,
				"-chr", nChr,
				"--output-prefix", `${output_prefix}_seg_ch${nChr}`,
			];
			try {
				await spawnNodeAsync(multi_align_to_segfile_args);

				if (VERBOSE) {
					console.log("nChr", nChr, "segfile");
				}

				fin_list[chrIdx].push("segfile");
			}
			catch (ex) {
				console.error(ex);
				console.error("error", nChr, ["node", ...multi_align_to_segfile_args].join(" "));
				return;
			}
			
			let crossover_ars = [
				"--max-old-space-size=4096",
				`${__dirname}/crossover.js`,
				"-dataset", argv_dataset_path,
				"--segfile", `${dataset.output_path}/${output_prefix}_seg_ch${nChr}.txt`,
				"--output-prefix", `${output_prefix}_co_ch${nChr}`
			];
			try {
				await spawnNodeAsync(crossover_ars);

				if (VERBOSE) {
					console.log("nChr", nChr, "crossover");
				}

				fin_list[chrIdx].push("crossover");
			}
			catch (ex) {
				console.error(ex);
				console.error("error", nChr, ["node", ...crossover_ars].join(" "));
				return;
			}

			let tetrad_chr_summary_args = [
				"--max-old-space-size=4096",
				`${__dirname}/tetrad_chr_summary.js`,
				"-dataset", argv_dataset_path,
				"-chr", nChr,
				"--seq", seq,
				"--co-list", `${dataset.output_path}/${output_prefix}_co_ch${nChr}_co.json`, 
				"--output-prefix", `${output_prefix}_ch${nChr}`
			];
			try {
				await spawnNodeAsync(tetrad_chr_summary_args);

				if (VERBOSE) {
					console.log("nChr", nChr, "summary");
				}

				fin_list[chrIdx].push("summary");
			}
			catch (ex) {
				console.error(ex);
				console.error("error", nChr, ["node", ...tetrad_chr_summary_args].join(" "));
				return;
			}
		});
		await Promise.all(tasks);

		if (VERBOSE) {
			console.log(fin_list);

			console.log(cmds.join("\n"));
		}
	}

	{
		let template_html = fs.readFileSync(`${__dirname}/template_viewer.html`).toString();
		let analyser_js = fs.readFileSync(`${__dirname}/analyser.js`).toString();
		let web_ui_js = fs.readFileSync(`${__dirname}/web_ui.js`).toString();
		
		let all_co = load_all_co();

		dataset.crossover_list = all_co;
		
		{
			dataset.results = genome_info_list[0].chr_list.map(function (chrName, chrIdx) {
				const nChr = chrIdx + 1;
				return `mafft_ch${nChr}.fa`;
			});
			const output_path = `${dataset.output_path}/debug_${output_prefix}.html`;

			output_html(template_html, output_path, false);
		
			if (VERBOSE) {
				console.log("output debug viewer:", output_path);
			}
		}

		//single html
		{
			let all_seq = load_all_seq();
			dataset.results = genome_info_list[0].chr_list.map(function (chrName, chrIdx) {
				return all_seq[chrIdx];
			});
			const output_path = `${dataset.output_path}/${output_prefix}.html`;

			output_html(template_html, output_path, true);
		
			console.log("output viewer:", output_path);
		}
		
		function output_html(input_html, output_path, inline_script) {
			let output_html = input_html.replace(`<script id="dataset.json" type="application/json"></script>`, `<script id="dataset.json" type="application/json">${JSON.stringify(dataset)}</script>`);

			if (inline_script) {
				output_html = output_html.replace(`<script src="analyser.js"></script>`, `<script>${analyser_js}</script>`);
				output_html = output_html.replace(`<script src="web_ui.js"></script>`, `<script>${web_ui_js}</script>`);
			}
			else {
				output_html = output_html.replace(`<script src="analyser.js"></script>`, `<script src="../src/analyser.js"></script>`);
				output_html = output_html.replace(`<script src="web_ui.js"></script>`, `<script src="../src/web_ui.js"></script>`);
			}

			fs.writeFileSync(output_path, output_html);
		}
	}

	output_final_table();
}

function load_all_seq() {
	let all_seq = [];
	genome_info_list[0].chr_list.forEach(function (chrName, chrIdx) {
		const nChr = chrIdx + 1;
		const seq = `${dataset.output_path}/mafft_ch${nChr}.fa`;
		
		const input_fasta = readFasta(seq);

		all_seq[chrIdx] = input_fasta;
	});
	return all_seq;
}

function load_all_co() {
	return genome_info_list[0].chr_list.map(function (chrName, chrIdx) {
		const nChr = chrIdx + 1;
		const co_list = load_co_list(`${dataset.output_path}/${output_prefix}_co_ch${nChr}_co.json`);
		return co_list;
	});
}

function load_co_list(co_list_filename) {
	const text = fs.readFileSync(co_list_filename).toString();
	const list = JSON.parse(text);
	list.forEach(co => {
		co.before = co.before.split(",").map(n => Number(n));
		co.after = co.after.split(",").map(n => Number(n));
	});
	return list;
}

function output_final_table() {
	const output_head = [
		"Chromosome",
		"simple CO", "CO + NCO",
		"Q/C SNV", "Q/C SNP", "Q/C InDel",
		"RIP Q", "RIP C",
		"illegitimate mutation",
		"SNP 2:2", "NCO 3:1", "NCO 4:0",
		"1n:3", "2n:2", "3n:1", "4n:0"
	];

	let final_table = [];

	genome_info_list[0].chr_list.map(function (chrName, chrIdx) {
		const nChr = chrIdx + 1;
		const tetrad_chr_summary = `${dataset.output_path}/${output_prefix}_ch${nChr}_summary.txt`;
		
		const text = fs.readFileSync(tetrad_chr_summary).toString();

		// @ts-ignore
		let table = table_to_object_list(tsv_parse(text), output_head, { start_row: 1 });

		final_table.push(...table);
	});

	let final_text = "";
	
	final_text += output_head.join("\t") + "\n";
	
	final_text += final_table.map(row => output_head.map(key => row[key]).join("\t")).join("\n");

	let output_path = `${dataset.output_path}/${output_prefix}_final_table.txt`;

	fs.writeFileSync(output_path, final_text);
	
	console.log("output table 1:", output_path);
}

function spawnNodeAsync(args) {
	cmds.push(["node", ...args].join(" "));

	return new Promise(function (resolve, reject) {
		let proc = child_process.spawn("node", args);
		proc.on("exit", function (code, signal) {
			if (code || signal) {
				reject({
					code, signal,
				});
			}
			else {
				resolve();
			}
		});
	});
}
