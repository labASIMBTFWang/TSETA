// @ts-check

const fs = require("fs");
const child_process = require("child_process");

const { inputFile, inputDirectory, inputNumber, inputText, inputSelect } = require(`${__dirname}/src/interact-util.js`).userInput;


main();

async function main() {
	let setting = {
		blastn_bin: "blastn",
		mafft_bin: "mafft",
	};
	
	if (fs.existsSync("./setting.json")) {
		let old_setting = JSON.parse(fs.readFileSync("./setting.json").toString());
		console.log("load setting", old_setting);

		setting.blastn_bin = old_setting.blastn_bin || "blastn";
		setting.mafft_bin = old_setting.mafft_bin || "mafft";
	}

	for (;;) {
		setting.blastn_bin = await inputText("path to blastn or blastn command", setting.blastn_bin);
		try {
			let cmd = `${setting.blastn_bin} -version`;
			let std_out = child_process.execSync(cmd);
			console.log(cmd);
			console.log(std_out.toString());
			break;
		}
		catch (ex) {
			console.error(ex);
		}
	}

	for (;;) {
		setting.mafft_bin = await inputText("path to mafft or mafft command", setting.mafft_bin);
		try {
			let cmd = `${setting.mafft_bin} --version`;
			let std_out = child_process.execSync(cmd);
			console.log(cmd);
			console.log(std_out.toString());
			break;
		}
		catch (ex) {
			console.error(ex);
		}
	}

	console.log("save setting", setting);

	fs.writeFileSync("setting.json", JSON.stringify(setting, null, "\t"));

	process.stdin.pause();
}

