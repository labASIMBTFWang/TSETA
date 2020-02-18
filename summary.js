
// @ts-check

const fs = require("fs");
const child_process = require("child_process");

const { tsv_parse, _table_to_object_list, table_to_object_list } = require("./tsv_parser.js");
const { argv_parse } = require("./util.js");
const { Dataset } = require("./dataset.js");
const { readFasta, saveFasta } = require("./fasta_util.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const argv_output_prefix = String(argv["--output-prefix"] || "");

const dataset = Dataset.loadFromFile(argv_dataset_path);

let cmds = [];

main();

async function main() {
	let fin_list = [];
	
	//for each chromosome
	let tasks = dataset.chrNames.map(async function (chrName, chrIdx) {
		const nChr = chrIdx + 1;
		const seq = `mafft_ch${nChr}.fa`;

		return;

		fin_list[chrIdx] = ["init"];

		let multi_align_to_segfile_args =  [
			"--max-old-space-size=4096",
			"multi_align_to_segfile.js",
			"-dataset", argv_dataset_path,
			"-i", seq,
			"-chr", nChr,
			"--output-prefix", `${argv_output_prefix}_seg_ch${nChr}`,
		];
		try {
			await spawnNodeAsync(multi_align_to_segfile_args);

			console.log("nChr", nChr, "segfile");

			fin_list[chrIdx].push("segfile");
		}
		catch (ex) {
			console.error(ex);
			console.error("error", nChr, ["node", ...multi_align_to_segfile_args].join(" "));
			return;
		}
		
		let crossover_ars = [
			"--max-old-space-size=4096",
			"crossover.js",
			"-dataset", argv_dataset_path,
			"--segfile", `${argv_output_prefix}_seg_ch${nChr}.txt`,
			"--output-prefix", `${argv_output_prefix}_co_ch${nChr}`
		];
		try {
			await spawnNodeAsync(crossover_ars);

			console.log("nChr", nChr, "crossover");

			fin_list[chrIdx].push("crossover");
		}
		catch (ex) {
			console.error(ex);
			console.error("error", nChr, ["node", ...crossover_ars].join(" "));
			return;
		}

		let chr_summary_args = [
			"--max-old-space-size=4096",
			"chr_summary.js",
			"-dataset", argv_dataset_path,
			"-chr", nChr,
			"--seq", seq,
			"--co-list", `${argv_output_prefix}_co_ch${nChr}_co.json`, 
			"--output-prefix", `${argv_output_prefix}_ch${nChr}`
		];
		try {
			await spawnNodeAsync(chr_summary_args);

			console.log("nChr", nChr, "summary");

			fin_list[chrIdx].push("summary");
		}
		catch (ex) {
			console.error(ex);
			console.error("error", nChr, ["node", ...chr_summary_args].join(" "));
			return;
		}
	});

	await Promise.all(tasks);

	console.log(fin_list);

	console.log(cmds.join("\n"));

	output_final_table();

	{
		let template_html = fs.readFileSync("template.html").toString();
		let analyser_js = fs.readFileSync("analyser.js").toString();
		let web_ui_js = fs.readFileSync("web_ui.js").toString();
		
		let all_seq = load_all_seq();
		let all_co = load_all_co();

		dataset.crossover_list = all_co;
		dataset.results = dataset.chrNames.map(function (chrName, chrIdx) {
			const nChr = chrIdx + 1;
			const seq = `mafft_ch${nChr}.fa`;
			return seq;
		});

		template_html = template_html.replace(`<script id="dataset.json" type="application/json"></script>`, `<script id="dataset.json" type="application/json">${JSON.stringify(dataset)}</script>`);
		// template_html = template_html.replace(`<script id="all_seq.json" type="application/json"></script>`, `<script id="all_seq.json" type="application/json">${JSON.stringify(all_seq)}</script>`);

		// template_html = template_html.replace(`<script src="analyser.js"></script>`, `<script>${analyser_js}</script>`);
		// template_html = template_html.replace(`<script src="web_ui.js"></script>`, `<script>${web_ui_js}</script>`);

		console.log({
			argv_output_prefix
		})
		let output_html = `${argv_output_prefix}.html`;

		fs.writeFileSync(output_html, template_html);

		console.log("output html", output_html);
	}
}

function load_all_seq() {
	let all_seq = [];
	dataset.chrNames.forEach(function (chrName, chrIdx) {
		const nChr = chrIdx + 1;
		const seq = `mafft_ch${nChr}.fa`;
		
		const input_fasta = readFasta(seq);

		all_seq[chrIdx] = input_fasta;
	});
	return all_seq;
}

function load_all_co() {
	return dataset.chrNames.map(function (chrName, chrIdx) {
		const nChr = chrIdx + 1;
		const co_list = load_co_list(`${argv_output_prefix}_co_ch${nChr}_co.json`);
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

	dataset.chrNames.map(function (chrName, chrIdx) {
		const nChr = chrIdx + 1;
		const chr_summary = `${argv_output_prefix}_ch${nChr}_summary.txt`;
		
		const text = fs.readFileSync(chr_summary).toString();

		// @ts-ignore
		let table = table_to_object_list(tsv_parse(text), output_head, { start_row: 1 });

		final_table.push(...table);
	});

	let final_text = "";
	
	final_text += output_head.join("\t") + "\n";
	
	final_text += final_table.map(row => output_head.map(key => row[key]).join("\t")).join("\n");

	fs.writeFileSync(`${argv_output_prefix}_final_table.txt`, final_text);
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
