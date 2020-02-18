
const fs = require("fs");
const child_process = require("child_process");

const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");

/**
 * @param {{}[]} array
 * @param {string} groupBy
 */
function array_groupBy(array, groupBy) {
	let group = {};
	array.forEach(row => {
		if (!group[row[groupBy]]) {
			group[row[groupBy]] = [];
		}
		group[row[groupBy]].push(row);
	});

	return group;
}

/**
 * @param {string[]} argv
 * @return {{[key:string]:string|boolean}}
 */
function argv_parse(argv) {
	let a = argv.slice(2);
	/**
	 * @type {{[key:string]:string|boolean}}
	 */
	let paramMap = {};

	a.forEach((name, idx, arr) => {
		if (name.startsWith("-")) {
			let val = arr[idx + 1];
			paramMap[name] = !val || val.startsWith("-") ? true : val;
		}
	});

	return paramMap;
}

/**
 * very old
 * @param {string} cmd
 * @param {string} [stdout_fname]
 * @param {string} [stderr_fname]
 * @returns {Promise<number>} errno
 */
function execAsync(cmd, stdout_fname, stderr_fname) {
	const debug_print_cmd = process.env.debug_print_cmd != null;
	if (debug_print_cmd) {
		console.log(cmd);
		return Promise.resolve(0);
	}

	return new Promise(function (resolve, reject) {
		let start_time = new Date();

		if (stdout_fname) {
			fs.writeFileSync(stdout_fname, start_time + "\r\n" + cmd + "\r\n", { flag: "a" });//append file end
		}
		if (stderr_fname) {
			fs.writeFileSync(stderr_fname, start_time + "\r\n" + cmd + "\r\n", { flag: "a" });//append file end
		}
		let proc = child_process.exec(cmd);

		let stdout_text = "";
		let stderr_text = "";

		proc.stdout.on("data", function (chunk) {
			stdout_text += chunk;
		});
		proc.stderr.on("data", function (chunk) {
			stderr_text += chunk;
		});

		proc.on("exit", function (nCode, sSignal) {
			stdout_text += "exit code: " + nCode + "\r\n";

			if (stdout_fname) {
				fs.writeFileSync(stdout_fname, stdout_text, { flag: "a" });//append file end
			}
			if (stderr_fname) {
				fs.writeFileSync(stderr_fname, stderr_text, { flag: "a" });//append file end
			}

			resolve(nCode);
		});

		proc.on("error", function (err) {
			console.error(err.stack);
			if (stderr_fname) {
				fs.writeFileSync(stderr_fname, err.stack, { flag: "a" });//append file end
			}
		});
	});
}

/**
 * old
 * @param {string} cmd
 * @param {string} [stdout_fname]
 * @param {string} [stderr_fname]
 * @returns {Promise<number>} errno
 */
function _execAsync(cmd, stdout_fname = null, stderr_fname = null) {
	const debug_print_cmd = process.env.debug_print_cmd != null;
	if (debug_print_cmd) {
		console.log(cmd);
		return Promise.resolve(0);
	}
	return new Promise(function (resolve, reject) {
		try {
			let start_time = new Date();

			if (stdout_fname) {
				fs.writeFileSync(stdout_fname, start_time + "\r\n" + cmd + "\r\n", { flag: "a" });//append file end
			}
			if (stderr_fname) {
				fs.writeFileSync(stderr_fname, start_time + "\r\n" + cmd + "\r\n", { flag: "a" });//append file end
			}
	
			let proc = child_process.exec(cmd);
	
			proc.stdout.on("data", function (chunk) {
				if (stdout_fname) {
					fs.writeFileSync(stdout_fname, chunk.toString(), { flag: "a" });//append file end
				}
			});
			proc.stderr.on("data", function (chunk) {
				console.error(chunk);
				if (stderr_fname) {
					fs.writeFileSync(stderr_fname, chunk.toString(), { flag: "a" });//append file end
				}
			});
	
			proc.on("exit", function (nCode, sSignal) {
				if (stdout_fname) {
					fs.writeFileSync(stdout_fname, "exit code: " + nCode + "\r\n", { flag: "a" });//append file end
				}
				
				if (sSignal) {
					console.error({
						cmd, sSignal
					});
					reject({
						cmd, sSignal,
					});
				}
				else {
					resolve(nCode);
				}
			});
	
			proc.on("error", function (err) {
				console.error("child_process.exec", err);
				if (stderr_fname) {
					fs.writeFileSync(stderr_fname, JSON.stringify(err, null, "\t"), { flag: "a" });//append file end
				}
	
				reject("cmd: " + cmd + "\n\terr:" + err);
			});
		}
		catch (ex) {
			console.error("execAsync", ex);
			//fs.writeFileSync("test.error.txt", JSON.stringify(ex, null, "\t"), { flag: "a" });
		}
	});
}

async function enum_chr_name(genomeId) {
	await execAsync(`samtools faidx ${genomeId}.genome.fa`);
	
	const text_tab = fs.readFileSync(`${genomeId}.genome.fa.fai`).toString();
	const header = [
		"chr", "length",// "byte_start", "_1", "_2"
	];
	let rows = table_to_object_list(tsv_parse(text_tab), header, { start_row: 0 });
	
	let results = rows.map((data, index) => {
		return `${(index + 1)}\t${data.chr}\t${data.length}`;
	}).join("\r\n");
	
	results = `Index#\tChr#\tLength\r\n` + results;
	
	fs.writeFileSync(`${genomeId}.length.txt`, results);
	
	rows = rows.sort((a, b) => Number(a.index) - Number(b.index));
	
	if (rows.length <= 0) {
		throw new Error("ls_chr_name(" + genomeId + ")");
	}
	
	return rows;
}

/**
 * @returns {{ list: { index: number, chr: string, length: number }[], map: {[chrName:string]:{ index: number, chr: string, length: number }} }}
 */
function loadChrLength(file_path) {
	let text_tab = fs.readFileSync(file_path).toString();
	let tab = tsv_parse(text_tab);
	let rows = table_to_object_list(tab, ["index", "chr", "length"], { start_row: 1 });
	
	/** @type {{ index: number, chr: string, length: number }[]} */
	let chrList = rows.map(row => {
		return {
			index: Number(row.index),
			chr: String(row.chr),
			length: Number(row.length),
		};
	});
	chrList = chrList.sort((a, b) => a.index - b.index);

	/** @type {{[chrName:string]:{ index: number, chr: string, length: number }}} */
	let chrMap = {};
	rows.forEach(row => {
		chrMap[row.chr] = {
			index: Number(row.index),
			chr: String(row.chr),
			length: Number(row.length),
		};
	});
	
	return {
		map: chrMap,
		list: chrList,
		/**
		 * @param {number} num
		 */
		getByNum: function (num) {
			return this.list[n - 1];
		},
	};
}

// /**
//  * @returns {{[chrName:string]:{ index: number, chr: string, length: number }}}
//  */
// function loadLengthMapFile(file_path) {
// 	let text_tab = fs.readFileSync(file_path).toString();
// 	let tab = tsv_parse(text_tab);
// 	let rows = table_to_object_list(tab, ["index", "chr", "length"], { start_row: 1 });
//
// 	let chrMap = {};
//
// 	rows.forEach(row => {
// 		chrMap[row.chr] = {
// 			index: Number(row.index),
// 			chr: String(row.chr),
// 			length: Number(row.length),
// 		};
// 	});
//
// 	return chrMap;
// }

function loadDataSet(dataset_path) {
	try {
		return JSON.parse(fs.readFileSync(dataset_path).toString());
	}
	catch (ex) {
		throw ex;
	}
}

module.exports.array_groupBy = array_groupBy;
module.exports.execAsync = execAsync;
module.exports.argv_parse = argv_parse;

module.exports.loadDataSet = loadDataSet;
module.exports.enum_chr_name = enum_chr_name;
module.exports.loadChrLength = loadChrLength;
