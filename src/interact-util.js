//@ts-check

const fs = require("fs");
const { EOL } = require("os");


/**
 * @param {string} message
 * @param {string} [default_value]
 */
async function inputFile(message, default_value) {
	let default_msg = default_value != null ? ` (default ${default_value})` : "";
	for (;;) {
		try {
			process.stdout.write(`${message}${default_msg}: (file path)${EOL}> `);
			return await getUserInput(String, v => fs.existsSync(v) && !fs.statSync(v).isDirectory(), default_value);
		}
		catch (ex) {
			console.warn("invalid value, try again.");
		}
	}
}
/**
 * @param {string} message
 * @param {string} [default_value]
 */
async function inputDirectory(message, default_value) {
	let default_msg = default_value != null ? ` (default ${default_value})` : "";
	for (;;) {
		try {
			process.stdout.write(`${message}${default_msg}: (directory path)${EOL}> `);
			return await getUserInput(String, v => fs.existsSync(v) && fs.statSync(v).isDirectory(), default_value);
		}
		catch (ex) {
			console.warn("invalid value, try again.");
		}
	}
}

/**
 * @param {string} message
 * @param {string} [default_value]
 */
async function inputText(message, default_value) {
	let default_msg = default_value != null ? ` (default ${default_value})` : "";
	for (;;) {
		try {
			process.stdout.write(`${message}${default_msg}: (text)${EOL}> `);
			return await getUserInput(String, null, default_value);
		}
		catch (ex) {
			console.warn("invalid value, try again.");
		}
	}
}

/**
 * @param {string} message
 * @param {{min?:number, max?:number, default?:number}} [options]
 */
async function inputNumber(message, options) {
	let default_msg = options && options.default != null ? ` (default ${options.default})` : "";
	for (;;) {
		try {
			let type_msg = options && options.min < options.max ? ` (${options.min} ~ ${options.max})` : " (number)";

			process.stdout.write(`${message}${default_msg}:${type_msg}${EOL}> `);
			if (options) {
				let { min, max } = options;
				if (max != null) {
					if (min != null && min < max) {
						return await getUserInput(Number, v => {
							let n = Number(v);
							return min <= n && n <= max;
						}, options.default);
					}
					else {
						throw new TypeError("min != null && min < max");
					}
				}
				else {
					return await getUserInput(Number, v => {
						let n = Number(v);
						return min <= n;
					}, options.default);
				}
			}
			else {
				return await getUserInput(Number, v => /\d+.\d+/.test(v), options.default);
			}
		}
		catch (ex) {
			console.warn("invalid value, try again.");
		}
	}
}

/**
 * @template T
 * @param {string} message
 * @param {T extends string[]} option_list
 * @param {string} [default_value]
 * @returns {Promise<T>}
 */
async function inputSelect(message, option_list, default_value) {
	let default_msg = default_value != null ? ` (default ${default_value})` : "";
	for (;;) {
		try {
			process.stdout.write(`${message}${default_msg}: (${option_list.join(", ")}) ${EOL}> `);
			/** @type {option_list} */
			// @ts-ignore
			let value = await getUserInput(String, v => option_list.indexOf(v) >= 0, default_value);
			return value;
		}
		catch (ex) {
			console.warn("invalid value, try again.");
		}
	}
}

/**
 * @template T
 * @param {function(any):T} type
 * @param {function(string):boolean} [validator]
 * @param {T} [default_value]
 * @returns {Promise<T>}
 */
function getUserInput(type, validator, default_value) {
	return new Promise(function (resolve, reject) {
		process.stdin.once("data", function (chunk) {
			let raw_str = chunk.toString();
			if ((raw_str == EOL || raw_str == "\n") && default_value != null) {
				if (!validator || validator(default_value.toString())) {
					process.stdout.write("Use default value: " + default_value.toString() + EOL);
					resolve(default_value);
				}
				else {
					reject(default_value);
				}
			}
			else {
				let input_str = raw_str.trim();
				if (input_str && (!validator || validator(input_str))) {
					let value = type(input_str);
					resolve(value);
				}
				else {
					reject(input_str);
				}
			}
		});
	});
}


const userInput = {
	inputFile,
	inputDirectory,
	inputNumber,
	inputText,
	inputSelect,
};
module.exports.userInput = userInput;

