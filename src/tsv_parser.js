

/**
 * @param {string} text
 * @param {string} splitter
 * @returns {string[][]}
 */
function tsv_parse(text, splitter = "\t") {
	let tab = [];
	let tr = [];
	let i = 0;

	text = text.trim().replace(/\r\n/g, "\n");

	while (i < text.length) {
		let c = text[i];
		
		if (c == '"') {
			let val = "";
			
			c = text[++i];//skip "
			while (c != '"' && c != null) {
				val += c;
				c = text[++i];
			}
			c = text[++i];//skip "
			
			tr.push(val);
			
			continue;
		}
		else if (c == splitter) {
			if (text[i - 1] == splitter) {//if this cell is empty
				tr.push("");
			}
			i++;//skip this splitter
			continue;
		}
		else if (c == "\n") {
			tab.push(tr);
			tr = [];
			i++;
			continue;
		}
		else {
			let val = "";
			
			while (c != splitter && c != '\n' && c != null) {
				val += c;
				c = text[++i];
			}
			
			tr.push(val);
			
			continue;
		}
	}
	if (tr.length) {
		tab.push(tr);
	}

	return tab;
}

/**
 * @param {any[][]} table
 * @param {number|string[]} header
 * @param {{start_row:number}} option
 * @returns {{[key:string]:string}[]}
 */
function _table_to_object_list(table, header, option = { start_row: 0, prototype: null }) {
	let start_row = option.start_row | 0;
	
	let _table = table.slice(start_row);

	if (header != null) {
		let head;
		if (Number.isSafeInteger(Number(header)) && header < _table.length) {
			for (let i = 0; i < header; ++i) {
				head = _table.shift();
			}
			head = _table.shift();
		}
		else if (Array.isArray(header)) {
			head = header;
		}
		
		//auto rename duplicate property
		head.reduce((prev, curr, i, arr) => {
			if (prev[arr[i]]) {
				prev[arr[i]]++;
				arr[i] += "_" + prev[arr[i]];
			}
			else {
				prev[arr[i]] = 1;
			}
			return prev;
		}, {});

		//to Object
		return _table.map((row) => {
			return row.reduce((prev, current, i) => {
				prev[head[i]] = current;
				return prev;
			}, option.prototype ? Object.create(option.prototype) : {});
		});
	}
	else {
		//@ts-ignore
		return _table;
	}
}

/**
 * @param {any[][]} table
 * @param {number|string[]} header
 * @param {{start_row:number}} option
 * @returns {{[key:string]:string}[]}
 */
function table_to_object_list(table, header, option = { start_row: 0 }) {
	let ret = _table_to_object_list(table, header, option);
	return ret;
}

if (typeof module == "object" && typeof module.exports == "object") {
	module.exports.tsv_parse = tsv_parse;
	module.exports.table_to_object_list = table_to_object_list;

	module.exports._table_to_object_list = _table_to_object_list;
}
