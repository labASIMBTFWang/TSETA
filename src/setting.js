// @ts-check

const fs = require("fs");

let default_value = {
	blastn_bin: "blastn",
	mafft_bin: "mafft",
}

/**
 * @returns {default_value}
 */
function loadSetting() {
	try {
		let text = fs.readFileSync("./setting.json").toString();
		
		let setting = JSON.parse(text);

		setting.blastn_bin = setting.blastn_bin || "blastn";
		setting.mafft_bin = setting.mafft_bin || "mafft";

		return setting;
	}
	catch (ex) {
		return default_value;
	}
}

module.exports.loadSetting = loadSetting;

