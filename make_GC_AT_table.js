
// @ts-check

const fs = require("fs");

const { argv_parse, loadChrLength, array_groupBy } = require("./util.js");
const { Dataset } = require("./dataset.js");
const  { GC_Content_Data } = require("./GC_content_util.js");
const  { parseFasta } = require("./fasta_util.js");

const argv = argv_parse(process.argv);

const argv_window_size = Number(argv["--window"] || argv["-w"]) | 0;
const argv_AT_island_GC_minus = (Number(argv["--minus"] || argv["-m"]) | 0);//6 // AT_island => gc(window) <= (gc(all) - minus)

/**
A = [...seq_list[0]].filter(a => a != "-").slice(0, 5000).filter(a => a == "A").length;
T = [...seq_list[0]].filter(a => a != "-").slice(0, 5000).filter(a => a == "T").length;
G = [...seq_list[0]].filter(a => a != "-").slice(0, 5000).filter(a => a == "G").length;
C = [...seq_list[0]].filter(a => a != "-").slice(0, 5000).filter(a => a == "C").length;
(G + C) / (A + T + G + C) * 100
 */

// const argv_input_file = String(argv["--input"] || argv["-i"]);
// const argv_output_prifix = String(argv["--output"] || argv["-o"]);

console.log({
	argv_AT_island_GC_minus, argv_window_size,
});

{
	const argv_dataset_path = String(argv["-dataset"] || "");
	const dataset = Dataset.loadFromFile(argv_dataset_path);

	const ref1_name = dataset.ref;
	const ref2_name = dataset.refs[1];

	{
		let argv_param = {
			ref1_name: dataset.ref,
			ref2_name: dataset.refs[1],
			spore_1_name: dataset.progeny_list[0],
			spore_2_name: dataset.progeny_list[1],
			spore_3_name: dataset.progeny_list[2],
			spore_4_name: dataset.progeny_list[3],
		};
		
		const header = ["Index#", "Chr#", "Length"].join("\t") + "\n";
		const genome_name_list = [
			argv_param.ref1_name,
			argv_param.ref2_name,
		];
		
		genome_name_list.forEach(genome_name => {
			let fname = `${genome_name}.length.txt`;
			if (!fs.existsSync(fname)) {
				console.log(genome_name, dataset.parental[genome_name]);

				let fa = parseFasta(fs.readFileSync(dataset.parental[genome_name]).toString());

				let out_text = "";
				out_text += header;
				out_text += Object.keys(fa).map((seq_name, idx) => [idx + 1, seq_name, fa[seq_name].length].join("\t")).join("\n");

				fs.writeFileSync(fname, out_text);

				console.log("output *.length.txt", fname);
			}
		});
	}

	const ref1_chr_list = loadChrLength(`./${ref1_name}.length.txt`);
	const ref2_chr_list = loadChrLength(`./${ref2_name}.length.txt`);

	let { all_gc_table: ref1_all_gc_table, all_island_table: ref1_all_island_table } = make_table(dataset.parental[ref1_name], ref1_name);
	let { all_gc_table: ref2_all_gc_table, all_island_table: ref2_all_island_table } = make_table(dataset.parental[ref2_name], ref2_name);

	{
		ref1_all_gc_table.forEach(group => {
			group.forEach(gc => {
				const ref1_chr_info = ref1_chr_list.map[gc.chr];
				gc.name = ref1_name;
				gc.chr = ref1_chr_info.index;
				gc.start = gc.start + 1;
				gc.end = Math.min(gc.end + 1, ref1_chr_info.length);
			});
		});
		ref2_all_gc_table.forEach(group => {
			group.forEach(gc => {
				const ref2_chr_info = ref2_chr_list.map[gc.chr];
				gc.name = ref2_name;
				gc.chr = ref2_chr_info.index;
				gc.start = gc.start + 1;
				gc.end = Math.min(gc.end + 1, ref2_chr_info.length);
			});
		});
		
		const output_gc_table_header = [
			"name", "chr", "start", "end", "gc"
		];
		let all_gc_list = [].concat(...ref1_all_gc_table, ...ref2_all_gc_table);
		let text_all_gc_list = all_gc_list.map(row => output_gc_table_header.map(key => row[key]).join("\t")).join("\n");

		fs.writeFileSync(`${ref1_name}_${ref2_name}_GC_content.txt`, text_all_gc_list);
	}

	{
		ref1_all_island_table.forEach(group => {
			group.forEach(island => {
				island.chr = ref1_chr_list.map[island.chr].index;
				island.start = island.start + 1;
				island.end = island.end + 1;
			});
		});

		const output_at_table_header = [
			"chr", "start", "end", "length"
		];
		let all_island_list = [].concat(...ref1_all_island_table);
		let text_all_island_list = all_island_list.map(row => output_at_table_header.map(key => row[key]).join("\t")).join("\n");
		
		fs.writeFileSync(`${ref1_name}_AT_island.txt`, text_all_island_list);
	}

	{
		let ref1_cen = [];
		let ref1_tel = [];
		
		for (let nChr = 1; nChr <= dataset.chrNames.length; ++nChr) {
			const chrIdx = nChr - 1;
			const AT_desc_list = [...ref1_all_island_table[chrIdx]].sort((a, b) => b.length - a.length);
			const cen_range = AT_desc_list[0];
			ref1_cen.push(cen_range);
			
			let tel1 = Object.assign({}, ref1_all_island_table[chrIdx][0], {
				start: 1,
			});
			
			let tel2 = Object.assign({}, ref1_all_island_table[chrIdx][ref1_all_island_table[chrIdx].length - 1], {
				end: ref1_chr_list.list[ref1_all_island_table[chrIdx][ref1_all_island_table[chrIdx].length - 1].chr - 1].length,
			});

			ref1_tel.push([tel1, tel2]);
		}

		console.log("cen");
		console.log(ref1_cen.map(a => [a.chr, a.start, a.end].join("\t")).join("\n"));
		
		console.log("tel");
		console.log(ref1_tel.map(([a, b]) => [a.chr, a.start, a.end, b.start, b.end].join("\t")).join("\n"));
	}
}

class AT_island_Data {
	constructor() {
		this.chr = 0;
		this.start = 0;
		this.end = 0;
		this.length = 0;
	}
}

function make_table(input_filename, output_prifix) {
	console.log("input fa", input_filename);

	const input_file = fs.readFileSync(input_filename).toString();
	const input_seq = parseFasta(input_file);
	
	//console.log("input fa", Object.keys(input_seq));
	
	const all_gc_content = get_all_gc_content(input_seq);
	const at_island_maxmum_gc = all_gc_content - argv_AT_island_GC_minus;
	
	//console.log("AT island maxmum GC%", at_island_maxmum_gc.toFixed(2));

	/** @type {GC_Content_Data[][]} */
	let all_gc_table = [];
	/** @type {AT_island_Data[][]} */
	let all_island_table = [];

	Object.keys(input_seq).forEach(function (name) {
		/** @type {GC_Content_Data[]} */
		let gc_table = [];
		/** @type {AT_island_Data[]} */
		let island_table = [];

		let seq = input_seq[name];
		let counter = {
			A: 0,
			T: 0,
			G: 0,
			C: 0,
		};
		
		let i = 0, length = 0, start = 0;
		for (; i < seq.length; ++i) {
			const v = seq[i].toUpperCase();
			
			++counter[v];

			++length;
			if (length >= argv_window_size) {
				let gc = counter.G + counter.C;
				let gc_content = 100 * gc / (counter.A + counter.T + gc);

				gc_table.push({
					chr: name,
					start: start,
					end: i,
					gc: gc_content,
				});

				if (gc_content <= at_island_maxmum_gc) {
					island_table.push({
						chr: name,
						start: start,
						end: i,
						length: length,
					});
				}

				if ((i - start + 1) > argv_window_size) {
					console.error(start, i, argv_window_size);
					throw new Error("if ((start + 1) > argv_window_size) {");
				}

				length = 0;
				start = i + 1;//next
				counter = {
					A: 0,
					T: 0,
					G: 0,
					C: 0,
				};//clear
			}
		}
		if (length < argv_window_size) {
			let gc = counter.G + counter.C;
			let gc_content = 100 * gc / (counter.A + counter.T + gc);

			gc_table.push({
				chr: name,
				start: start,
				end: i,
				gc: gc_content,
			});
		}
		all_gc_table.push(gc_table);
		all_island_table.push(island_table);

		// console.log(name, "gc_table", gc_table.length);
		// console.log(name, "island_table", island_table.length);
	});

	{
		const output_gc_table_header = [
			"chr", "start", "end", "gc"
		];
		let all_gc_list = [].concat(...all_gc_table);
		let text_all_gc_list = all_gc_list.map(row => output_gc_table_header.map(key => row[key]).join("\t")).join("\n");

		fs.writeFileSync(`${output_prifix}_GC_content.txt`, text_all_gc_list);
	}

	function merge_island() {
		return all_island_table.map(island_list => {
			let megered_list = [];
			let list = island_list.slice(1);

			let current_island = island_list[0];
			list.forEach(next_island => {
				if ((current_island.end + 1) == next_island.start) {
					current_island.end = next_island.end;
					current_island.length = current_island.length + next_island.length;
				}
				else {
					megered_list.push(current_island);

					current_island = next_island;
				}
			});

			///.........

			return megered_list;
		});
	}

	{
		//console.log(all_island_table[0].length);
		all_island_table = merge_island();
		//console.log(all_island_table[0].length);

		const output_at_table_header = [
			"chr", "start", "end", "length"
		];
		let all_island_list = [].concat(...all_island_table);
		let text_all_island_list = all_island_list.map(row => output_at_table_header.map(key => row[key]).join("\t")).join("\n");
		
		fs.writeFileSync(`${output_prifix}_AT_island.txt`, text_all_island_list);
	}

	return {
		all_gc_table, all_island_table
	};
}

function get_all_gc_content(input_seq) {
	let counter = {
		A: 0,
		T: 0,
		G: 0,
		C: 0,
	};

	Object.keys(input_seq).forEach(function (name) {
		let seq = input_seq[name];
		for (let i = 0; i < seq.length; ++i) {
			const v = seq[i].toUpperCase();
			++counter[v];
		}
	});

	//console.log(counter);

	let all_gc = counter.G + counter.C;
	let all_atgc = counter.A + counter.T + all_gc;
	let all_gc_content = 100 * all_gc / all_atgc;

	// console.log("GC content", {
	// 	all_gc,
	// 	all_atgc,
	// 	all_gc_content,
	// 	"%": all_gc_content.toFixed(2),
	// });

	return all_gc_content;
}


