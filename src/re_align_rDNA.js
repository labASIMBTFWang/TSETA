// @ts-check

const { argv_parse } = require("./util.js");

const argv = argv_parse(process.argv);

if (argv["--old-version"]) {
	require("./align_rDNA_v1.js");
}
else {
	require("./align_rDNA.js");
}
