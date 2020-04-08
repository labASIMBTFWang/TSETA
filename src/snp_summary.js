//@ts-check

const fs = require("fs");
const Path = require("path");
const child_process = require("child_process");

const { argv_parse, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");

const { readFasta, saveFasta, joinFastaSeq } = require("./fasta_util.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);
const argv_output_summary_only = !!argv["--output-summary-only"];

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const VERBOSE = !!argv["--verbose"];

const genome_info_list = dataset.loadGenomeInfoList();
dataset.load_GC_content();
dataset.load_rDNA_info();

const output_prefix = encodeURIComponent(dataset.name);


main();

function main() {
	if (!argv_output_summary_only) {
		output_table();
	}
	output_viewer();
}

function output_table() {
	genome_info_list[0].chr_list.forEach((chrInfo, chrIdx) => {
		const nChr = chrIdx + 1;
		const input_fa = readFasta(`${dataset.output_path}/mafft_ch${nChr}.fa`);
		const chr_name_list = genome_info_list.map(genomeInfo => genomeInfo.chr_list[chrIdx].chr);
		const seq_list = chr_name_list.map(chrName => input_fa[chrName]);

		let snv_rows = [];
		let indel_rows = [];
		let _snp_list = [];

		let pos_list = chr_name_list.map(_ => 1);
		for (let i = 0; i < seq_list[0].length; ++i) {
			const ref1 = seq_list[0][i];
			const columns = [];
			let has_snp, is_indel;

			columns.push(pos_list[0], ref1);
			
			seq_list.forEach((target_seq, target_seqIdx) => {
				let pos;
				if (seq_list[target_seqIdx][i] == "-") {
					is_indel = true;
					pos = pos_list[target_seqIdx] - 1;
				}
				else {
					pos = pos_list[target_seqIdx];
					++pos_list[target_seqIdx];
				}

				if (target_seqIdx >= 1) {
					let target = seq_list[target_seqIdx][i];
					if (target != ref1) {
						columns.push(pos, target);
						has_snp = true;
					}
					else {
						columns.push(null, null);
					}
				}
			});

			if (has_snp) {
				if (is_indel) {
					indel_rows.push(columns);
				}
				else {//is snv
					snv_rows.push(columns);
				}
				_snp_list.push(columns);
			}
		}//for each bp
		
		try {
			output_table("snv", snv_rows);
		}
		catch (ex) {
			console.error(ex);
		}
		
		try {
			output_table("indel", indel_rows);
		}
		catch (ex) {
			console.error(ex);
		}
		
		try {
			output_table("snp", _snp_list);
		}
		catch (ex) {
			console.error(ex);
		}

		/**
		 * 
		 * @param {"snv"|"indel"|"snp"} type
		 * @param {(string|number)[][]} rows
		 */
		function output_table(type, rows) {
			if (VERBOSE) {
				console.log(nChr, type + ":", rows.length);
			}
	
			const output_path = `${dataset.output_path}/${type}_ch${nChr}.txt`;
			
			const output_header = ["ref_pos", "ref"];
			seq_list.slice(1).map((_, targetIdx) => output_header.push("target" + (1 + targetIdx) + "_pos", "target" + (1 + targetIdx)));
			
			const ws = fs.createWriteStream(output_path);
			ws.write(output_header.join("\t") + "\n");
			rows.forEach(cols => {
				let row = cols.map(v => v != null ? v : "").join("\t") + "\n";
				ws.write(row);
			});
			ws.end();
			console.log("output table:", output_path);

			//fs.writeFileSync(output_path, output_header.join("\t") + "\n");
			
			// rows.forEach(cols => {
			// 	let row = cols.map(v => v != null ? v : "").join("\t") + "\n";
			// 	fs.writeFile(output_path, row, { flag: "a" }, function (err) {
			// 		if (err) {
			// 			console.error(err);
			// 		}
			// 	});
			// });
			
			// const output_table = rows.map(cols => cols.map(v => v != null ? v : "").join("\t")).join("\n");
			// //console.log(output_table.length);
			// fs.writeFile(output_path, output_table, { flag: "a" }, function (err) {
			// 	if (err) {
			// 		console.error(err);
			// 	}
			// 	console.log("output:", output_path);
			// });
		}
	});
}

function output_viewer() {
	let template_html = fs.readFileSync(`${__dirname}/template_viewer.html`).toString();
	let analyser_js = fs.readFileSync(`${__dirname}/analyser.js`).toString();
	let web_ui_js = fs.readFileSync(`${__dirname}/web_ui.js`).toString();
	
	{
		dataset.results = genome_info_list[0].chr_list.map(function (chrName, chrIdx) {
			const nChr = chrIdx + 1;
			return `mafft_ch${nChr}.fa`;
		});
		const output_path = `${dataset.output_path}/debug_${output_prefix}_snp.html`;
		
		output_html(template_html, output_path, false);
		if (VERBOSE) {
			console.log("output debug viewer", output_path);
		}
	}

	//single html
	{
		let all_seq = load_all_seq();
		dataset.results = genome_info_list[0].chr_list.map(function (chrName, chrIdx) {
			return all_seq[chrIdx];
		});
		const output_path = `${dataset.output_path}/${output_prefix}_snp.html`;

		output_html(template_html, output_path, true);
	
		console.log("output viewer:", output_path);
	}
	
	function output_html(input_html, output_path, inline_script) {
		let output_html = input_html.replace(`<script id="dataset.json" type="application/json"></script>`, `<script id="dataset.json" type="application/json">${JSON.stringify(dataset)}</script>`);
		//output_html = output_html.replace(`<script id="all_seq.json" type="application/json"></script>`, `<script id="all_seq.json" type="application/json">${JSON.stringify(all_seq)}</script>`);

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

