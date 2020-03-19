// @ts-check

const child_process = require("child_process");
const fs = require("fs");
const Path = require("path");
const { loadSetting } = require("./setting.js");

const setting = loadSetting();


class BlastnCoord {
	constructor() {
		this.query = "";
		this.subject = "";
		
		this.identity = 0;
		this.align = 0;
		this.mismatch = 0;
		this.gap = 0;
		
		this.qstart = 0;
		this.qend = 0;
		this.sstart = 0;
		this.send = 0;

		this.evalue = 0;
		this.score = 0;
	}

	get q_min() {
		return Math.min(this.qstart, this.qend);
	}
	get q_max() {
		return Math.max(this.qstart, this.qend);
	}
	
	get s_min() {
		return Math.min(this.sstart, this.send);
	}
	get s_max() {
		return Math.max(this.sstart, this.send);
	}

	get strand() {
		if (this.sstart > this.send || this.qstart > this.qend) {
			return -1;
		}
		else if (this.sstart < this.send && this.qstart < this.qend) {
			return 1;
		}
		else {
			//debugger;
			return 0;
		}
	}

	get slen() {
		return Math.abs(this.send - this.sstart);
	}
	
	get qlen() {
		return Math.abs(this.qend - this.qstart);
	}

	// /**
	//  * @param {number} start
	//  */
	// dist_sstart(start) {
	// 	return Math.abs(this.sstart - start);
	// }

	// /**
	//  * @param {number} start
	//  */
	// dist_send(start) {
	// 	return Math.abs(this.send - start);
	// }

	toArray() {
		return BlastnCoord.tableHeader.map(k => this[k]);
	}

	toString() {
		return JSON.stringify(this, null, "\t");
	}

	static get tableHeader() {
		return [
			"query", "subject",
			"identity", "align",
			"mismatch", "gap",
			"qstart", "qend", "sstart", "send",
			"evalue", "score"
		];
	}

	/**
	 * @param {Partial<BlastnCoord>} obj
	 * @returns {BlastnCoord}
	 */
	static fromObject(obj) {
		// let coord = new BlastnCoord();
		// coord.qstart = obj.qstart;
		// coord.qend = obj.qend;
		// coord.sstart = obj.sstart;
		// coord.send = obj.send;
		// return coord;
		return Object.assign(new BlastnCoord(), obj);
	}
}

/**
 * @param {string} cmd
 * @param {boolean} output_stdout
 * @param {boolean} output_stderr
 */
function execAsync(cmd, output_stdout, output_stderr) {
	// return new Promise(function (resolve, reject) {
	// 	let proc = child_process.spawn(cmd);
		
	// 	proc.on("error", function(err) {
	// 		reject(err);
	// 	});

	// 	let stdout_buf_list = [];
	// 	proc.stdout.on("data", function(buffer) {
	// 		if (output_stdout) {
	// 			stdout_buf_list.push(buffer);
	// 		}
	// 	});
		
	// 	let stderr_buf_list = [];
	// 	proc.stderr.on("data", function(buffer) {
	// 		if (output_stderr) {
	// 			stderr_buf_list.push(buffer);
	// 		}
	// 	});
		
	// 	proc.on("exit", function(code, signal) {
	// 		let stdout, stderr;
	// 		if (output_stdout) {
	// 			stdout = Buffer.concat(stdout_buf_list);
	// 		}
	// 		if (output_stderr) {
	// 			stderr = Buffer.concat(stderr_buf_list);
	// 		}
	// 		resolve({
	// 			code, signal, stdout, stderr
	// 		});
	// 	});
	// });

	return new Promise(function (resolve, reject) {
		child_process.exec(cmd, function (err, stdout, stderr) {
			if (output_stderr && stderr) {
				console.error(stderr.toString());
			}
		//
			if (err) {
				reject(err);
			}
			else if (output_stdout) {
				resolve(stdout.toString());
			}
			else {
				resolve();
			}
		});
	});
}

/**
 * @param {string} query_file
 * @param {string} subject_file
 * @param {number} [qstart]
 * @param {number} [qend]
 * @param {number} [sstart]
 * @param {number} [send]
 * @param {string} [_task_name]
 * @returns {Promise<string>}
 */
function exec_blastn(query_file, subject_file, qstart, qend, sstart, send, _task_name) {
	let query_loc = Number.isSafeInteger(qstart) && Number.isSafeInteger(qend) ? `-query_loc ${qstart}-${qend}` : "";
	let subject_loc = Number.isSafeInteger(sstart) && Number.isSafeInteger(send) ? `-subject_loc ${sstart}-${send}` : "";
	let cmd = `${setting.blastn_bin} -query ${query_file} ${query_loc} -subject ${subject_file} ${subject_loc} -outfmt 6`;

	_task_name = _task_name || "default";

	console.log(cmd);

	return new Promise(function (resolve, reject) {
		child_process.exec(cmd, function (err, stdout, stderr) {
			if (stderr) {
				console.error(cmd);
				console.error(stderr.toString());
			}

			// @ts-ignore
			const tmp_path = globalThis.dataset.tmp_path;

			if (err) {
				console.error(err, { query_file, subject_file, qstart, qend, sstart, send, _task_name });
				
				fs.writeFileSync(`${tmp_path}/ma_util_blastn/ma_util_blastn_${_task_name}.txt`, JSON.stringify(err) + "\n", { flag: "a" });
				reject(err);
			}
			else {
				// resolve({
				// 	stdout: stdout.toString(),
				// 	//stderr: stderr.toString(),
				// });
				let text = stdout.toString();

				fs.writeFileSync(`${tmp_path}/ma_util_blastn/ma_util_blastn_${_task_name}.txt`, cmd + "\n" + text + "\n", { flag: "a" });

				resolve(text);
			}
		});
	});
}
/**
 * @param {string} query_file
 * @param {string} subject_file
 * @param {number} [qstart]
 * @param {number} [qend]
 * @param {number} [sstart]
 * @param {number} [send]
 * @param {string} [args]
 * @returns {Promise<string>}
 */
function exec_blastn_Ex(query_file, subject_file, qstart, qend, sstart, send, args) {
	let query_loc = Number.isSafeInteger(qstart) && Number.isSafeInteger(qend) ? `-query_loc ${qstart}-${qend}` : "";
	let subject_loc = Number.isSafeInteger(sstart) && Number.isSafeInteger(send) ? `-subject_loc ${sstart}-${send}` : "";
	let cmd = `${setting.blastn_bin} -query ${query_file} ${query_loc} -subject ${subject_file} ${subject_loc} ${args} -outfmt 6`;

	console.log(cmd);

	return new Promise(function (resolve, reject) {
		child_process.exec(cmd, function (err, stdout, stderr) {
			if (stderr) {
				console.error(cmd);
				console.error(stderr.toString());
			}

			if (err) {
				console.error(err, { query_file, subject_file, qstart, qend, sstart, send });
				reject(err);
			}
			else {
				let text = stdout.toString();

				resolve(text);
			}
		});
	});
}

/**
 * @param {string} text
 * @returns {BlastnCoord[]}
 */
function parseBlastnResults(text) {
	let _table = text.trim().split("\n").filter(a => a && a.length).map(line => {
		let columns = line.split("\t");
		return BlastnCoord.fromObject({
			query: columns[0],
			subject: columns[1],
			identity: Number(columns[2]),
			align: Number(columns[3]),
			mismatch: Number(columns[4]),
			gap: Number(columns[5]),
			qstart: Number(columns[6]),
			qend: Number(columns[7]),
			sstart: Number(columns[8]),
			send: Number(columns[9]),
			evalue: Number(columns[10]),
			score: Number(columns[11]),
		});
	});
	return _table;
}

/**
 * return duplicate group
 * @param {BlastnCoord[]} rows
 * @returns {BlastnCoord[][]}
 */
function groupByOverlap(rows) {
	/** @type {BlastnCoord[][]} */
	let groups = [];

	const tempRows = rows.slice(0);

	while (tempRows.length) {
		const left = tempRows.splice(0, 1)[0];
		if (!left) {
			throw new Error("tempRows.splice(0, 1)[0] -> undefined");
		}
		const overlap_group = [left];

		for (let j = 0; j < tempRows.length; ++j) {
			const right = tempRows[j];
			if (isCollide(left.sstart, left.send, right.sstart, right.send) ||
				isCollide(left.qstart, left.qend, right.qstart, right.qend)
			) {
				tempRows.slice(j, 0);
				overlap_group.push(right);
			}
		}

		groups.push(overlap_group.sort((a, b) => b.score - a.score));
	}

	return groups.sort((a, b) => b[0].score - a[0].score);
}

function isCollide(x11, x12, x21, x22) {
	return x11 <= x22 && x12 >= x21;
}

/**
 * @param {BlastnCoord[]} rows
 */
function removeOverlap(rows) {
	/** @type {Set<BlastnCoord>} */
	let removes = new Set();

	for (let i = 0; i < rows.length; ++i) {
		const left = rows[i];
		if (!removes.has(left)) {
			if (left.strand <= 0) {
				removes.add(left);
			}
			else {
				for (let j = i + 1; j < rows.length; ++j) {
					const right = rows[j];
					if (isCollide(left.sstart, left.send, right.sstart, right.send) ||
						isCollide(left.qstart, left.qend, right.qstart, right.qend)
					) {
						if (left.score < right.score) {
							removes.add(left);
						}
						else if (right.score < left.score) {
							removes.add(right);
						}
					}
				}
			}
		}
	}

	let filtered = rows.filter(left => !removes.has(left));
	
	// let prev_sstart = 0;
	// filtered = rows = rows.filter(left => {
	// 	let delta = left.sstart - prev_sstart;
	// 	if (delta >= 0 && delta <= init_max_delta) {
	// 		prev_sstart = left.sstart;
	// 		return true;
	// 	}
	// });

	return filtered;
}

/**
 * @param {string} query_file
 * @param {string} subject_file
 * @param {number} qstart
 * @param {number} qend
 * @param {number} sstart
 * @param {number} send
 * @param {function(Partial<BlastnCoord>,any[]):boolean} filter
 * @param {any[]} filter_params
 * @param {string} _task_name
 * @returns {Promise<BlastnCoord[]>}
 */
async function blastn_coord(query_file, subject_file, qstart, qend, sstart, send, filter, filter_params, _task_name) {
	if (qstart >= qend || sstart >= send) {
		debugger;
		return null;
	}
	let text = await exec_blastn(query_file, subject_file, qstart, qend, sstart, send, _task_name);

	let _table = parseBlastnResults(text);

	// _table.forEach(row => {
	// 	if (row.send == 2641826) {//2004187
	// 		debugger;
	// 	}
	// });
	
	//_table = _table.sort((a, b) => b.score - a.score);
	_table = _table.sort((a, b) => a.send - b.send);

	// let table = _table.filter(row => {
	// 	return (
	// 		row.qstart < row.qend &&
	// 		row.sstart < row.send &&
	// 		//row.sstart >= next_start
	// 		row.send >= next_start
	// 	);
	// });

	let table = removeOverlap(_table).filter(row => filter(row, filter_params));
	//let table = _table.filter(row => filter(row, filter_params));

	return table;

	// //console.log(table);
	// if (table.length) {
	// 	return table[0];
	// }
	// //console.log(next_start, _table, query_file, subject_file, qstart, qend, sstart, send, next_start);
	// return null;
}


module.exports.BlastnCoord = BlastnCoord;
module.exports.execAsync = execAsync;
module.exports.exec_blastn = exec_blastn;
module.exports.exec_blastn_Ex = exec_blastn_Ex;
module.exports.parseBlastnResults = parseBlastnResults;
module.exports.blastn_coord = blastn_coord;
module.exports.isCollide = isCollide;
module.exports.removeOverlap = removeOverlap;
module.exports.groupByOverlap = groupByOverlap;

