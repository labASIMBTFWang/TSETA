//@ts-check

const fs = require("fs");

/**
 * 
 * @param {string} filename
 * @param {{ [seqName:string]:string }} fasta
 */
function saveFasta(filename, fasta) {
	let str = "";
	Object.keys(fasta).forEach(key => {
		if (fasta[key] && fasta[key].length) {
			str += "> ";
			str += key;
			str += "\n";

			str += fasta[key];
			str += "\n";
		}
		else {
			console.warn("saveFasta:", key, "=>", fasta[key]);
		}
	});
	try {
		fs.writeFileSync(filename, str);
	}
	catch (ex) {
		console.error("write file error:", filename);
		throw ex;
	}
}

/**
 * @param {string} filename
 * @returns {{[chr:string]:string}}
 */
function readFasta(filename) {
	try {
		let text = fs.readFileSync(filename).toString();
		return parseFasta(text);
	}
	catch (ex) {
		console.error("read file error:", filename);
		throw ex;
	}
}

/**
 * @param {string} in_seq
 * @returns {{[chr:string]:string}}
 */
function parseFasta(in_seq) {
	let all = in_seq.split(">");

	/** @type {{[chr:string]:string}} */
	let results = {};

	all.filter(a => a.length).map(async (sub_fa) => {
		let li = sub_fa.indexOf("\n");
		let out_name = sub_fa.slice(0, li).trim().replace(/:/g, "");
		let out_file_name = out_name.match(/^([^ ]+)/)[1];
		let out_seq = sub_fa.slice(li).trim();

		out_seq = out_seq.replace(/\n/g, "");
		out_seq = out_seq.replace(/\r/g, "");
		out_seq = out_seq.replace(/ /g, "");

		results[out_file_name] = out_seq.toUpperCase();
	});

	return results;
}

function save_fasta_async(path, fasta) {
	return new Promise(function (resolve, reject) {
		try {
			let fp = fs.createWriteStream(path);
			Object.keys(fasta).forEach(key => {
				fp.write("> ");
				fp.write(key);
				fp.write("\n");

				fp.write(fasta[key]);
				fp.write("\n");
			});
			fp.end();
			fp.on("close", resolve);
		}
		catch (ex) {
			reject(ex);
		}
	});
}


/**
 * @param {{[chr:string]:string}[]} fasta_list
 * @returns {{[chr:string]:string}}
 */
function joinFastaSeq(fasta_list) {
	/** @type {{[chr:string]:string}} */
	let out = {};
	Object.keys(fasta_list[0]).forEach(k => {
		if (!out[k]) {
			out[k] = "";
		}
		fasta_list.forEach(seq => {
			out[k] += seq[k];
		});
	});
	return out;
}


/**
 * index: 0 ~ length - 1
 * @param {string} ma_seq
 * @returns {number[]}
 */
function chrPos_to_multialign_posMap(ma_seq) {
	let nPos = 1;
	let posmap = [];
	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		if (element != "-") {
			posmap[nPos] = index + 1;
			++nPos;
		}
	}
	return posmap;
}

/**
 * index: 0 ~ length - 1
 * @param {string} ma_seq
 * @returns {number[]}
 */
function multialign_to_chrPos_posMap(ma_seq) {
	let nPos = 0;
	let posmap = [];

	posmap[0] = 0;

	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		
		if (element != "-") {
			++nPos;
		}
		
		posmap[index + 1] = nPos;
	}

	return posmap;
}


module.exports.readFasta = readFasta;
module.exports.saveFasta = saveFasta;
module.exports.parseFasta = parseFasta;
module.exports.joinFastaSeq = joinFastaSeq;

module.exports.chrPos_to_multialign_posMap = chrPos_to_multialign_posMap;
module.exports.multialign_to_chrPos_posMap = multialign_to_chrPos_posMap;

