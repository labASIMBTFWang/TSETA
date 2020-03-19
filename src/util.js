// @ts-check

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
			//fs.writeFileSync(`${dataset.tmp_path}/test.error.txt`, JSON.stringify(ex, null, "\t"), { flag: "a" });
		}
	});
}



module.exports.array_groupBy = array_groupBy;
module.exports.execAsync = execAsync;
module.exports.argv_parse = argv_parse;

