// @ts-check

let dataset = JSON.parse(document.getElementById("dataset.json").innerText);

let seq_list = [
	"", "",
	"", "", "", "",
];

class GC_Content {
	constructor() {
		this.name = "";
		this.chr = "";
		this.start = 0;
		this.end = 0;
		this.gc = 0;
	}
}
/** @type {{[name:string]:{[chr:string]:GC_Content[]}}} */
let gc_content = {};
/** @type {{[name:string]:{[chr:string]:number}}} */
let gc_content_average = {};

// @ts-ignore
window.get_r1a = function(start, end) {
	let ref1_pos = 1;
	let ra = seq_list.map(_ => "");

	[...seq_list[0]].forEach((ref1, pos) => {
		if (start <= ref1_pos && ref1_pos <= end) {
			for (let i = 0; i < seq_list.length; ++i) {
				ra[i] += seq_list[i][pos];
			}
		}
		if (ref1 != "-") {
			++ref1_pos;
		}
	});
	console.log(ra);
};

// @ts-ignore
window.get_ma = function(start, end) {
	for (let i = 0; i < seq_list.length; ++i) {
		console.log(seq_list[i].slice(start, end));
	}
};

{
	let el_markers_table = document.getElementById("markers_table");
	//let el_display_buttons_group = document.getElementById("display_buttons_group");

	MAKER_NAME.forEach(name => {
		appendMakerstable(name);
	});
	function appendMakerstable(name) {
		let elem = document.createElement("div");
		elem.classList.add("makers");
		elem.innerHTML = `<span>${name}</span>`;
		el_markers_table.append(elem);
	}
	// function appendDisplayButton() {
	// 	let elem = document.createElement("label");
	// 	elem.innerHTML = `<input type="checkbox" id="el_display_31" checked /> display ${name}`;
	// 	el_display_buttons_group.append(elem);
	// }
}

let el_display_bp_pos = document.getElementById("el_display_bp_pos");
let el_display_ref1_bp_pos = document.getElementById("el_display_ref1_bp_pos");
let el_display_ref2_bp_pos = document.getElementById("el_display_ref2_bp_pos");
let el_display_ref1_ref2_score = document.getElementById("el_display_ref1_ref2_score");

/** @type {HTMLInputElement} */
// @ts-ignore
let el_input_start = document.getElementById("el_input_start");
/** @type {HTMLInputElement} */
// @ts-ignore
let el_input_end = document.getElementById("el_input_end");
/** @type {HTMLInputElement} */
// @ts-ignore
let el_input_ref1_start = document.getElementById("el_input_ref1_start");
/** @type {HTMLInputElement} */
// @ts-ignore
let el_input_ref1_end = document.getElementById("el_input_ref1_end");

/** @type {HTMLCanvasElement} */
// @ts-ignore
let canvas = document.getElementById("canvas");
/** @type {CanvasRenderingContext2D} */
let ctx = canvas.getContext("2d");

{
	canvas.width = canvas.parentElement.clientWidth;
	canvas.style.width = canvas.width + "px";
	
	canvas.height = canvas.parentElement.clientHeight;
	canvas.style.height = canvas.height + "px";
}

let mouse_bp_pos = 1;
let bp_start = 0;
let bp_end = 0;

let getViewLength = () => (bp_end - bp_start + 1);
let getBPperPixel = () => getViewLength() / canvas.width;
let getPixelPerBP = () => canvas.width / getViewLength();
let g_maxPixelPerBP = 1;

class ViewerState {
	// /** @type {boolean} */
	// _display_13 = true;
	/** @type {boolean} */
	_display_31 = true;
	// /** @type {boolean} */
	// _display_04 = true;
	/** @type {boolean} */
	_display_40 = true;

	/** @type {boolean} */
	_display_22indel = true;
	/** @type {boolean} */
	_display_13indel = true;
	/** @type {boolean} */
	_display_31indel = true;
	/** @type {boolean} */
	_display_40indel = true;

	/** @type {boolean} */
	_display_illegitimate_mutation = true;

	/** @type {number} */
	_nChr = 1;

	_disable_max_length = true;
	_max_length = 7148324;

	_crossover_only = true;

	_padding_right = 16 * 4;

	/** @type {[number,number][]} */
	_range_makers = [];

	/** @type {number[]} */
	_position_makers = [];
	/** @type {number[]} */
	_position_ref1_makers = [];
	/** @type {number[]} */
	_position_ref2_makers = [];

	_rip_display_weight = 10000000;//100 is small

	_display_parent_rip = false;
	
	animationFrameId = 0;

	_gc_content_clip_indel = true;

	_display_rdna_border = false;

	constructor() {
		this.seg_row_height = 32;
		this.seg_row_separate = 5;

		this.fill_prev_color = true;
	}

	get_plot_height() {
		let { seg_row_height, seg_row_separate } = this;
		let height = 0;

		height += (seg_row_height + seg_row_separate) * seq_list.length;//seq

		height += (seg_row_height + seg_row_separate) * 1;//CO maker

		height += (seg_row_height + seg_row_separate) * spore_cmp_array.length;//nco makers

		height += (seg_row_height * 2 + seg_row_separate) * 2;//GC%

		return height;
	}

	get ref1() {
		return dataset.parental_list[0];
	}
	get ref2() {
		return dataset.parental_list[1];
	}
	get nChr() {
		return this._nChr;
	}
	set nChr(value) {
		if (this._nChr != value) {
			this._nChr = value;
			drawFrame();
		}
	}

	get disable_max_length() {
		return this._disable_max_length;
	}
	set disable_max_length(value) {
		this._disable_max_length = value;
		drawFrame();
	}

	get max_length() {
		if (this._disable_max_length) {
			return seq_list[0].length;
		}
		else {
			return this._max_length;
		}
	}
	set max_length(value) {
		//if (this._max_length != value) {
			this._max_length = value;
			drawFrame();
		//}
	}

	get max_view_width() {
		let all_scale = seq_list[0].length / viewerState.max_length;
		let max_view_width = canvas.clientWidth * all_scale;
		return max_view_width - this._padding_right;
	}

	get display_rdna_border() {
		return this._display_rdna_border;
	}
	set display_rdna_border(value) {
		if (this._display_rdna_border != value) {
			this._display_rdna_border = value;
			drawFrame();
		}
	}
	
	get display_22indel() {
		return this._display_22indel;
	}
	set display_22indel(value) {
		if (this._display_22indel != value) {
			this._display_22indel = value;
			drawFrame();
		}
	}
	
	// get display_13() {
	// 	return this._display_13;
	// }
	// set display_13(value) {
	// 	if (this._display_13 != value) {
	// 		this._display_13 = value;
	// 		drawFrame();
	// 	}
	// }
	get display_31() {
		return this._display_31;
	}
	set display_31(value) {
		if (this._display_31 != value) {
			this._display_31 = value;
			drawFrame();
		}
	}
	
	get display_13indel() {
		return this._display_13indel;
	}
	set display_13indel(value) {
		if (this._display_13indel != value) {
			this._display_13indel = value;
			drawFrame();
		}
	}
	
	get display_31indel() {
		return this._display_31indel;
	}
	set display_31indel(value) {
		if (this._display_31indel != value) {
			this._display_31indel = value;
			drawFrame();
		}
	}
	
	// get display_04() {
	// 	return this._display_04;
	// }
	// set display_04(value) {
	// 	if (this._display_04 != value) {
	// 		this._display_04 = value;
	// 		drawFrame();
	// 	}
	// }

	get display_40() {
		return this._display_40;
	}
	set display_40(value) {
		if (this._display_40 != value) {
			this._display_40 = value;
			drawFrame();
		}
	}

	get display_40indel() {
		return this._display_40indel;
	}
	set display_40indel(value) {
		if (this._display_40indel != value) {
			this._display_40indel = value;
			drawFrame();
		}
	}

	get display_illegitimate_mutation() {
		return this._display_illegitimate_mutation;
	}
	set display_illegitimate_mutation(value) {
		if (this._display_illegitimate_mutation != value) {
			this._display_illegitimate_mutation = value;
			drawFrame();
		}
	}

	get crossover_only() {
		return this._crossover_only;
	}
	set crossover_only(value) {
		if (this._crossover_only != value) {
			this._crossover_only = value;
			drawFrame();
		}
	}

	get padding_right() {
		return this._padding_right;
	}
	set padding_right(value) {
		if (this._padding_right != value) {
			this._padding_right = value;
			drawFrame();
		}
	}

	get display_parent_rip() {
		return this._display_parent_rip;
	}
	set display_parent_rip(value) {
		this._display_parent_rip = value;
		drawFrame();
	}

	get rip_display_weight() {
		return this._rip_display_weight;
	}
	set rip_display_weight(value) {
		this._rip_display_weight = value;
		drawFrame();
	}

	get gc_content_clip_indel() {
		return this._gc_content_clip_indel;
	}
	set gc_content_clip_indel(value) {
		this._gc_content_clip_indel = value;
		drawFrame();
	}

	/**
	 * @param {number} start
	 * @param {number} end
	 */
	pushRangeMaker(start, end) {
		this._range_makers.push([start, end]);
		drawFrame();
	}
	popRangeMaker() {
		let range = this._range_makers.pop();
		drawFrame();
		return range;
	}

	/**
	 * @param {number} position
	 */
	pushPositionMaker(position) {
		this._position_makers.push(position);
		
		if (pos_ref1_uint32array) {
			let ref1_pos = pos_ref1_uint32array[position - 1];
			if (ref1_pos != null) {
				this._position_ref1_makers.push(ref1_pos);
			}
		}
		
		if (pos_ref2_uint32array) {
			let ref2_pos = pos_ref2_uint32array[position - 1];
			if (ref2_pos != null) {
				this._position_ref2_makers.push(ref2_pos);
			}
		}
		
		drawFrame();
	}
	popPositionMaker() {
		let pos = this._position_makers.pop();
		let ref1_pos = this._position_ref1_makers.pop();
		let ref2_pos = this._position_ref2_makers.pop();
		drawFrame();
		return {
			pos, ref1_pos, ref2_pos,
		}
	}
}

const viewerState = new ViewerState();

// @ts-ignore
window.viewerState = viewerState;

const analyser_options = {
	get nChr() { return viewerState.nChr },
	get rDNA_info() { return dataset.rDNA_info; },
	get co_list() { return dataset.crossover_list[viewerState.nChr - 1]; },
	get fill_prev_color() { return true; },
};

/** @type {HTMLSelectElement} */
// @ts-ignore
let el_select_colorset = document.getElementById("el_select_colorset");
function onChangeColorSet(evt, colorset_name) {
	if (colorset_name) {
		el_select_colorset.value = colorset_name;
	}
	if (el_select_colorset.value == "view") {
		current_colorset = $color_set_view;
	}
	else {
		current_colorset = $color_set_print;
	}
	Object.keys(current_colorset).forEach(name => {
		let el = document.getElementById("el_" + name);
		if (el) {
			el.style.background = current_colorset[name];
		}
	});
	drawFrame();
};

{
	/** @type {HTMLSelectElement} */
	
	el_select_colorset.onchange = (evt) => onChangeColorSet(evt);
	el_select_colorset.onchange(null);

	{
		marker_order.forEach(label => {
			/** @type {HTMLInputElement} */
			// @ts-ignore
			let el_display_ = document.getElementById("el_display_" + label);
			if (el_display_) {
				// @ts-ignore
				el_display_.onchange = (evt) => viewerState["display_" + label] = evt.target.checked;
				el_display_.checked = viewerState["display_" + label];
			}
			else {
				console.log("el_display_" + label);
			}
		});
	}
	
	let el_input_chr = document.getElementById("el_input_chr");
	// @ts-ignore
	el_input_chr.value = viewerState._nChr;
	// @ts-ignore
	el_input_chr.oninput = (evt) => {
		viewerState.nChr = Number(evt.target.value);
	
		loadData();
	};

	let el_input_max_length = document.getElementById("el_input_max_length");
	// @ts-ignore
	el_input_max_length.value = viewerState._max_length;
	// @ts-ignore
	el_input_max_length.oninput = (evt) => viewerState.max_length = Number(evt.target.value);
	viewerState.max_length = viewerState._max_length;

	let el_input_disable_max_length = document.getElementById("el_input_disable_max_length");
	// @ts-ignore
	el_input_disable_max_length.value = viewerState._disable_max_length;
	// @ts-ignore
	el_input_disable_max_length.oninput = (evt) => viewerState.disable_max_length = !evt.target.checked;
	viewerState.disable_max_length = viewerState._disable_max_length;


	let el_input_rip_display_weight = document.getElementById("el_input_rip_display_weight");
	// @ts-ignore
	el_input_rip_display_weight.value = viewerState._rip_display_weight;
	// @ts-ignore
	el_input_rip_display_weight.oninput = (evt) => viewerState.rip_display_weight = Number(evt.target.value);
	viewerState.rip_display_weight = viewerState._rip_display_weight;
}

function set_ui_seq_pos(pageX) {
	mouse_bp_pos = screen_to_bp(pageX);
	
	el_display_bp_pos.innerText = mouse_bp_pos.toString();

	if (pos_ref1_uint32array) {
		let ref1_pos = pos_ref1_uint32array[mouse_bp_pos - 1];
		if (ref1_pos != null) {
			el_display_ref1_bp_pos.innerText = ref1_pos.toString();
		}
	}
	
	if (pos_ref2_uint32array) {
		let ref2_pos = pos_ref2_uint32array[mouse_bp_pos - 1];
		if (ref2_pos != null) {
			el_display_ref2_bp_pos.innerText = ref2_pos.toString();
		}
	}

	if (ref1_ref2_score_uint32array) {
		let score = ref1_ref2_score_uint32array[mouse_bp_pos - 1];
		if (score != null) {
			el_display_ref1_ref2_score.innerText = score.toString();
		}
	}
}

window.onmousemove = function (evt) {
	let rect = canvas.getBoundingClientRect();
	if (!(evt.pageY >= rect.top && evt.pageY <= rect.bottom)) {
		return;
	}

	set_ui_seq_pos(evt.pageX);
};

canvas.onclick = function (evt) {
	let rect = canvas.getBoundingClientRect();
	if (!(evt.pageY >= rect.top && evt.pageY <= rect.bottom)) {
		return;
	}
	
	set_ui_seq_pos(evt.pageX);
};

let mouse_bp_pos_0 = 0;
canvas.onmousedown = function (evt) {
	set_ui_seq_pos(evt.pageX);

	if (evt.which == 1) {
		if(evt.ctrlKey) {
			viewerState.pushPositionMaker(mouse_bp_pos);
			evt.preventDefault();
		}
		else if(evt.shiftKey) {
			viewerState.pushRangeMaker(bp_start, bp_end);
			evt.preventDefault();
		}
	}
	else if (evt.which == 3) {
		if(evt.ctrlKey) {
			viewerState.popPositionMaker();
			evt.preventDefault();
		}
		else if(evt.shiftKey) {
			viewerState.popRangeMaker();
			evt.preventDefault();
		}
	}

	mouse_bp_pos_0 = mouse_bp_pos;
};
canvas.onmousemove = function (evt) {
	if (evt.buttons) {
		set_ui_seq_pos(evt.pageX);

		let move_bp = mouse_bp_pos_0 - mouse_bp_pos;

		set_view_range(bp_start + move_bp, bp_end + move_bp);
	}
};
canvas.onmouseup = function (evt) {
	set_ui_seq_pos(evt.pageX);

	let move_bp = mouse_bp_pos_0 - mouse_bp_pos;

	set_view_range(bp_start + move_bp, bp_end + move_bp);
};


// @ts-ignore
canvas.onmousewheel = function (evt) {
	evt.preventDefault();
};

if (window.onwheel === null) {
	window.onwheel = onMouseWheel;
}
else {
	window.onmousewheel = onMouseWheel;
}
/**
 * @param {MouseWheelEvent} evt
 */
function onMouseWheel(evt) {
	let rect = canvas.getBoundingClientRect();
	if (!(evt.pageY >= rect.top && evt.pageY <= rect.bottom)) {
		return;
	}

	set_ui_seq_pos(evt.pageX);

	let start, end;

	let deltaY = 0;

	// @ts-ignore
	deltaY = (evt.deltaY || evt.wheelDelta || evt.wheelDeltaY);
	deltaY = deltaY > 0 ? 1 : (deltaY < 0 ? -1 : 0);
	
	const zoom_scale = 0.25;
	let d_left = (mouse_bp_pos - bp_start);
	let d_right = (bp_end - mouse_bp_pos);
	if (deltaY > 0) {
		start = bp_start + Math.ceil(Math.max(1, d_left * zoom_scale) * (-deltaY));
		end = bp_end - Math.ceil(Math.max(1, d_right * zoom_scale) * (-deltaY));
	}
	else if (deltaY < 0) {
		start = bp_start + Math.ceil(d_left * zoom_scale * (-deltaY));
		end = bp_end - Math.ceil(d_right * zoom_scale * (-deltaY));
	}

	set_view_range(start, end);
};

/**
 * @param {number} start
 * @param {number} end
 */
function set_view_range(start, end) {
	if (start < end) {
		if (start < 0) {
			end += -start;
			start = 1;
			end = Math.min(end, seq_list[0].length);
		}
		if (end > seq_list[0].length) {
			start -= end - seq_list[0].length;
			end = seq_list[0].length;
			start = Math.max(1, start);
		}
	
		start = Math.max(1, start);
		end = Math.min(end, seq_list[0].length);
	
		if (start != bp_start || bp_end != end) {
			bp_start = start;
			bp_end = end;
			
			el_input_start.value = bp_start.toString();
			el_input_end.value = bp_end.toString();
		
			drawFrame();
		}
	}
}

/**
 * @param {number} pos_x
 */
function screen_to_bp(pos_x) {
	const length = bp_end - bp_start;
	const scale = 1 / (length + 1);
	const bp_size = viewerState.max_view_width * scale;
	const rect = canvas.getBoundingClientRect();

	if (pos_x >= rect.left) {
		let center_bp_pos;

		if (pos_x <= rect.right) {
			center_bp_pos = bp_start + Math.trunc((pos_x - rect.left) / bp_size);
			center_bp_pos = Math.max(0, Math.min(center_bp_pos, seq_list[0].length));
		}
		else {
			center_bp_pos = seq_list[0].length;
		}
		
		return center_bp_pos;
	}
	else {
		return 0;
	}
}

async function loadData() {
	let fa = await (async function () {
		let obj = dataset.results[viewerState.nChr - 1];
		if (typeof obj == "string") {
			let raw_fa = await fetchData(obj, "text");
			let fa = await parseFasta(raw_fa);

			return fa;
		}
		else {
			return obj;
		}
	})();

	{
		let name_map = dataset["GC content"];
		
		let average = {};
		Object.keys(name_map).forEach(name => {
			average[name] = {};
			Object.keys(name_map[name]).forEach(chr =>{
				let first = name_map[name][chr][0];
				let w = Math.abs(first.end - first.start);
				average[name][chr] = name_map[name][chr].reduce((tt, curr) => {
					return tt + curr.gc / Math.abs(curr.end - curr.start);
				}, 0) / name_map[name][chr].length * w;
			});
		});
		// @ts-ignore
		gc_content = name_map;
		// @ts-ignore
		window.$gc_content = gc_content;
		// @ts-ignore
		gc_content_average = average;
		// @ts-ignore
		window.$gc_average = average;
	}

	//init
	bp_start = 1;
	bp_end = 20;
	el_input_start.value = bp_start.toString();
	el_input_end.value = bp_end.toString();

	el_input_ref1_start.oninput = (evt) => {
		if (ref1_pos_uint32array) {
			// @ts-ignore
			let newVal = Number(evt.target.value);
			bp_start = ref1_pos_uint32array[newVal] | 0;
			bp_start = Math.max(1, bp_start);
			el_input_start.value = bp_start.toString();
			if (bp_end > bp_start) {
				drawFrame();
			}
		}
	};
	el_input_ref1_end.oninput = (evt) => {
		if (ref1_pos_uint32array) {
			// @ts-ignore
			let newVal = Number(evt.target.value);
			bp_end = ref1_pos_uint32array[newVal] || seq_list[0].length;
			bp_end = Math.min(bp_end, seq_list[0].length);
			el_input_end.value = bp_end.toString();
			if (bp_end > bp_start) {
				drawFrame();
			}
		}
	};

	el_input_start.oninput = (evt) => {
		// @ts-ignore
		let newVal = Number(evt.target.value);
		bp_start = Math.max(1, newVal);
		el_input_start.value = bp_start.toString();
		if (bp_end > bp_start) {
			drawFrame();
		}
	};
	el_input_end.oninput = (evt) => {
		// @ts-ignore
		let newVal = Number(evt.target.value);
		bp_end = Math.min(newVal, seq_list[0].length);
		el_input_end.value = bp_end.toString();
		if (bp_end > bp_start) {
			drawFrame();
		}
	};

	//seq_list clear
	seq_id_list = Object.keys(fa);
	seq_id_list.forEach((id, i) => seq_list[i] = fa[id]);

	canvas.height = viewerState.get_plot_height();
	canvas.style.height = viewerState.get_plot_height() + "px";

	calc_seg_reg(seq_list, analyser_options);
	
	// // @ts-ignore
	// window.$seg = seg_snp;//no save
	// @ts-ignore
	//window.parental_cmp = parental_cmp_uint8array;
	// @ts-ignore
	window.spore_cmp_array = spore_cmp_array;
	// @ts-ignore
	window.ref1_pos_uint32array = ref1_pos_uint32array;
	// @ts-ignore
	window.pos_ref1_uint32array = pos_ref1_uint32array;
	// @ts-ignore
	window.ref2_pos_uint32array = ref2_pos_uint32array;
	// @ts-ignore
	window.pos_ref2_uint32array = pos_ref2_uint32array;
	// @ts-ignore
	window.ref1_ref2_score_uint32array = ref1_ref2_score_uint32array;

	//max length
	bp_end = Math.min(bp_end, seq_list[0].length);

	el_input_start.max = (seq_list[0].length - 1).toString();
	el_input_end.max = seq_list[0].length.toString();

	bp_start = 1;
	bp_end = seq_list[0].length;
	el_input_start.value = bp_start.toString();
	el_input_end.value = bp_end.toString();
	g_maxPixelPerBP = getPixelPerBP();

	drawFrame();
}

document.getElementById("el_show_all").onclick = function show_all() {
	bp_start = 1;
	bp_end = seq_list[0].length;
	el_input_start.value = bp_start.toString();
	el_input_end.value = bp_end.toString();

	drawFrame();
}

function drawFrame() {
	if (!viewerState.animationFrameId) {
		viewerState.animationFrameId = requestAnimationFrame(async function (time) {
			try {
				render(analyser_options);
			}
			catch (ex) {
				console.error(ex);
			}
			viewerState.animationFrameId = 0;
		});
		return true;
	}
	else {
		return false;
	}
}

//[].map();
/**
 * @param {AnalyserOptions} options
 */
function render(options) {
	if (!seq_list[0].length) {
		return;
	}
	_render(options);

	let v_scale = getPixelPerBP() / g_maxPixelPerBP;

	v_scale = seq_list[0].length / (bp_end - bp_start + 1);

	let list = [
		"scale: " + (v_scale).toFixed(2).toString(),
	];

	document.getElementById("scale").innerText = list.join(",");
}

//ctx.canvas.width / (bp_end - bp_start + 1)
/**
 * @param {AnalyserOptions} options
 */
function _render(options) {
	let max_view_width = viewerState.max_view_width;

	const scale = 1 / (bp_end - bp_start + 1);
	const bp_size = max_view_width * scale;
	const bp_per_px = 1 / bp_size;
	
	const bp_per_px_max = seq_list[0].length / max_view_width;
	const bp_per_px_min1 = Math.max(1, Math.max(1, bp_per_px) / bp_per_px_max * viewerState.rip_display_weight);

	ctx.lineJoin = "round";

	ctx.setTransform(1, 0, 0, 1, 0, 0);
	//ctx.translate(0.5, 0.5);
	//ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
	ctx.fillStyle = "#FFFFFF"
	ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);

	let allow_draw_small_object = bp_size / 4 * 3 >= 7;
	let font_size = Math.max(9, Math.min(bp_size / 4 * 3, 16));
	if (font_size == 9) {
		ctx.lineWidth = 0.5;
	}

	if (allow_draw_small_object) {
		onChangeColorSet(null, "view");
	}
	else {
		onChangeColorSet(null, "print");
	}

	let seg_row_height = viewerState.seg_row_height;//32
	let seg_row_separate = viewerState.seg_row_separate;//5
	let offset_y = 0;//seg_row_height * 2;
	
	let seg_row_hheight = Math.trunc(seg_row_height * 0.5);

	ctx.textAlign = "center";
	ctx.textBaseline = "middle";

	let prev_22_colorId = [
		ColorID.dad,
		ColorID.mom
	];
	/**
	 * @param {number} pos_start bp_start
	 */
	function calc_prev_parentalSNP_colorId(pos_start) {
		if (prev_22_colorId.length != 2) {
			return;
		}
		for (let i = pos_start - 1; i >= 0; --i) {
			let ref1 = seq_list[0][i];
			let ref2 = seq_list[1][i];
			let a = seq_list[2][i];
			let b = seq_list[3][i];
			let c = seq_list[4][i];
			let d = seq_list[5][i];
			let ss = [a, b, c, d];

			let rc1 = ss.filter(s => s == ref1).length;
			let rc2 = ss.filter(s => s == ref2).length;
			
			if (ref1 != ref2 && (rc1 + rc2) == 4) {//SNP 0:4 or 1:3 or 2:2 or 3:1 or 4:0
				prev_22_colorId[2] = (a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff));
				prev_22_colorId[3] = (b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff));
				prev_22_colorId[4] = (c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff));
				prev_22_colorId[5] = (d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff));
				return;
			}
		}
	}
	
	//let ct = 0;
	for (let seg_id = 0; seg_id < region_rect.length; ++seg_id) {
		let ss = region_rect[seg_id];

		let is_first_identical_seg = true;
		
		let prev_col_id = null;
		for (let si = 0; si < ss.length; ++si) {
			let start = ((1 - bp_start) + ss[si].start) * bp_size;
			let end = ((1 - bp_start) + ss[si].end) * bp_size;
			let length = end - start;
			let color_id = ss[si].col;

			if (color_id == null) {//no data
				continue;//transparent
			}

			if (!(start <= max_view_width)) {
				continue;
			}

			if (allow_draw_small_object) {
				if (end >= 0 && length >= 1) {
					let x1 = Math.max(0, start);
					let xlen = Math.min(end, max_view_width) - x1;

					let colId = color_id & ColorID.mask;

					ctx.filter = "brightness(0.8) saturate(2)";
					
					if (is_first_identical_seg && colId == ColorID.identical) {
						calc_prev_parentalSNP_colorId(bp_start);
						colId = prev_22_colorId[seg_id];
						is_first_identical_seg = false;
						prev_col_id = colId;
						ctx.filter = "saturate(0.5) brightness(1.5)";
					}
					if ((color_id & ColorID.mask) == ColorID.identical) {
						if (prev_col_id != null) {
							colId = prev_col_id;
							ctx.filter = "saturate(0.5) brightness(1.5)";
						}
					}
					else {
						prev_col_id = colId;
					}

					let col = args_colors[colId];

					ctx.beginPath();
					ctx.rect(x1, seg_id * (seg_row_height + seg_row_separate) + offset_y, xlen, seg_row_height);
					ctx.fillStyle = col;
					ctx.fill();
					ctx.globalAlpha = 1;
					ctx.filter = "none";
					//++ct;
				}
			}
			else if (length >= 1) {
				if (end >= 0) {
					if ((color_id & ColorID.indel_mask) == ColorID.diff) {
						debugger;
					}
					if (seg_id >= 2 && (color_id & ColorID.indel_mask)) {
						// nothing // large gap
						//prev_col_id = null;//reset to empty
						//skip
					}
					else {
						let b_fill;

						let colId = color_id & ColorID.mask;
						if (is_first_identical_seg && colId == ColorID.identical) {
							calc_prev_parentalSNP_colorId(ss[si].start);
							colId = prev_22_colorId[seg_id];
							let col = args_colors[colId];
							ctx.fillStyle = col;
							b_fill = true;
							is_first_identical_seg = false;
							prev_col_id = colId;
						}
						let col = args_colors[colId];

						let x1 = Math.max(0, start);
						let xlen = Math.min(end + 1, max_view_width) - x1;
						
						if ((color_id & ColorID.mask) == ColorID.identical) {
							if (prev_col_id != null) {
								if (args_colors[prev_col_id] == current_colorset["dad"]) {
									ctx.fillStyle = current_colorset["dad_bk"];
									b_fill = true;
								}
								else if (args_colors[prev_col_id] == current_colorset["mom"]) {
									ctx.fillStyle = current_colorset["mom_bk"];
									b_fill = true;
								}
								if (seg_id == 0) {
									ctx.fillStyle = current_colorset["dad_bk"];
									b_fill = true;
								}
								else if (seg_id == 1) {
									ctx.fillStyle = current_colorset["mom_bk"];
									b_fill = true;
								}
							}
						}
						else {
							if (seg_id < 2) {
								let indel = color_id & ColorID.indel_mask;
								if (!indel) {
									ctx.fillStyle = col;
									b_fill = true;
								}
							}
							else {
								if (viewerState.crossover_only) {
									if (col == current_colorset["dad"]) {
										ctx.fillStyle = current_colorset["dad_bk"];
										b_fill = true;
									}
									else if (col == current_colorset["mom"]) {
										ctx.fillStyle = current_colorset["mom_bk"];
										b_fill = true;
									}
									else {
										ctx.fillStyle = col;
										b_fill = true;
									}
								}
								else {
									ctx.fillStyle = col;
									b_fill = true;
								}
							}
							if (colId == ColorID.dad || colId == ColorID.mom) {
								prev_col_id = color_id & ColorID.mask;
							}
						}

						if (b_fill) {
							ctx.beginPath();
							ctx.rect(x1, seg_id * (seg_row_height + seg_row_separate) + offset_y, xlen, seg_row_height);
							ctx.fill();
						}
					}
				}
			}
			else {
				let sj = si;
				let next_start = 0;
				let next_end = 0;
				let n_col = [
					0, 0, 0, 0,
					0, 0, 0, 0,
					0, 0, 0, 0,
					0, 0, 0, 0
				];
				let next_length;
				let merge_length;
				let has_rip = false;
				//
				//console.time("merge_seg");
				for (; sj < ss.length; ++sj) {
					next_start = ((1 - bp_start) + ss[sj].start) * bp_size;
					next_end = ((1 - bp_start) + ss[sj].end) * bp_size;
					
					if (next_length >= 2) {
						break;
					}

					next_length = next_end - next_start;
					
					let ml = next_end - start;
					if (next_length >= 2) {
						break;
					}
					merge_length = ml;

					if (next_start > max_view_width) {
						break;
					}
					
					n_col[ss[sj].col] = (n_col[ss[sj].col] || 0) + next_length * [1, 1, 1, bp_per_px_min1, bp_per_px_min1, 1][ss[sj].col & ColorID.mask];

					has_rip = ss[sj].col == ColorID.dad_rip || ss[sj].col == ColorID.mom_rip;
					
					if (merge_length >= 1) {//merge size in px
						break;
					}
				}

				if (next_end >= 0) {
					//skip no diff if rect_length <= 5
					if (merge_length <= 2) {
						n_col[2] = 0;
						n_col[4] = 0;
					}
					
					let max_col = [...n_col].sort((a, b) => b - a)[0];
					if (max_col) {
						let color_id = n_col.findIndex(a => a == max_col);
						if (seg_id >= 2 && (color_id & ColorID.indel_mask)) {
							// nothing // large gap
							//prev_col_id = null;//reset to empty
							//skip
						}
						else {
							let colId = color_id & ColorID.mask;
							if (is_first_identical_seg && colId == ColorID.identical) {
								calc_prev_parentalSNP_colorId(ss[si].start);
								colId = prev_22_colorId[seg_id];
								is_first_identical_seg = false;
								prev_col_id = colId;
							}
							let col = args_colors[colId];

							let x1 = Math.max(0, start);
							let xlen = Math.min(x1 + merge_length + 1, max_view_width) - x1;
							
							let b_fill;

							if ((color_id & ColorID.mask) == ColorID.identical) {
								if (prev_col_id != null) {
									if (args_colors[prev_col_id] == current_colorset["dad"]) {
										ctx.fillStyle = current_colorset["dad_bk"];
										b_fill = true;
									}
									else if (args_colors[prev_col_id] == current_colorset["mom"]) {
										ctx.fillStyle = current_colorset["mom_bk"];
										b_fill = true;
									}
									if (seg_id == 0) {
										ctx.fillStyle = current_colorset["dad_bk"];
										b_fill = true;
									}
									else if (seg_id == 1) {
										ctx.fillStyle = current_colorset["mom_bk"];
										b_fill = true;
									}
								}
							}
							else {
								if (viewerState.crossover_only) {
									if (col == current_colorset["dad"]) {
										ctx.fillStyle = current_colorset["dad_bk"];
										b_fill = true;
									}
									else if (col == current_colorset["mom"]) {
										ctx.fillStyle = current_colorset["mom_bk"];
										b_fill = true;
									}
									else {
										ctx.fillStyle = col;
										b_fill = true;
									}
								}
								else {
									ctx.fillStyle = col;
									b_fill = true;
								}
								if (colId == ColorID.dad || colId == ColorID.mom) {
									prev_col_id = color_id & ColorID.mask;
								}
							}
							if (b_fill) {
								ctx.fillRect(x1, seg_id * (seg_row_height + seg_row_separate) + offset_y, xlen, seg_row_height);
								// if (colId != has_rip)  {
								// 	debugger;
								// } 
							}
						}
					}
					else {
						//large gap
					}
				}
				si = sj - 1;
			}
		}
	}

	if (allow_draw_small_object) {
		ctx.strokeStyle = "black";
		ctx.fillStyle = "black";

		// ctx.shadowBlur = 3;
		// ctx.shadowColor = "white";

		// bp text
		for (let seg_id = 0; seg_id < region_rect.length; ++seg_id) {
			for (let si = bp_start - 1; si < bp_end; ++si) {
				let x1 = ((1 - bp_start) + si) * bp_size;
				let x2 = ((1 - bp_start) + si + 1) * bp_size;

				ctx.beginPath();
				ctx.rect(x1, seg_id * (seg_row_height + seg_row_separate) + offset_y, bp_size, seg_row_height);
				ctx.stroke();
				
				ctx.font = Math.trunc(font_size) + "px Arial";
				
				//center text
				ctx.fillText(seq_list[seg_id][si], x2 - (bp_size * 0.5), (seg_id + 1) * (seg_row_height + seg_row_separate) - seg_row_hheight + offset_y);
			}
		}
	}

	if (options.rDNA_info.chr == viewerState.nChr) {
		let min_len_rDNA = Math.min(...options.rDNA_info.data.map(d => {
			return Math.min(...d.repeats.map((a, b) => Math.abs(a[1] - a[0])));
		}));
		if ((min_len_rDNA * bp_size) >= 1) {
			ctx.save();

			//if (x2 >= 0 && x1 <= max_view_width) {
			ctx.beginPath();
			ctx.rect(0, 0, max_view_width, 6 * (seg_row_height + seg_row_separate));
			ctx.clip();

			for (let seg_id = 0; seg_id < seq_list.length; ++seg_id) {
				let y = seg_id * (seg_row_height + seg_row_separate);
				let cy = y + seg_row_hheight;
				
				let ali_rep = options.rDNA_info.data[seg_id].alignment_repeats;
				let min_raw_start = Math.min(...ali_rep[0]);
				ali_rep.forEach((ar_data, rep_idx) => {
					let raw_start = Math.min(...ar_data);
					let raw_len = Math.abs(ar_data[1] - ar_data[0]);

					let strand = ar_data[1] < ar_data[0] ? -1 : 0;

					let x1 = ((1 - bp_start) + raw_start) * bp_size;
					let x2 = ((1 - bp_start) + raw_start + raw_len - 1) * bp_size;
					let len = x2 - x1;
					let hlen = len * 0.5;
					let ahlen = len * 0.3;
					let cx = (x1 + x2) * 0.5;
					
					const qy = seg_row_hheight * 0.5;

					ctx.save();
					{
						ctx.translate(cx, cy);
						ctx.scale(strand, 1);

						//arrow
						ctx.beginPath();
						ctx.moveTo(0 - hlen, 0 - qy);
						ctx.lineTo(0 + ahlen, 0 - qy);
						ctx.lineTo(0 + ahlen, 0 - seg_row_hheight);
						ctx.lineTo(0 + hlen, 0);
						ctx.lineTo(0 + ahlen, 0 + seg_row_hheight);
						ctx.lineTo(0 + ahlen, 0 + qy);
						ctx.lineTo(0 - hlen, 0 + qy);
						ctx.closePath();
						
						ctx.globalAlpha = Math.min(bp_size >= 1 ? (1.5 - Math.min(bp_size * 2, 32) / 32) : 1, 1);
						ctx.fillStyle = current_colorset.rDNA;
						ctx.fill();
						ctx.globalAlpha = 1;
					}
					ctx.restore();
					
					if (viewerState.display_rdna_border) {
						ctx.strokeStyle = "#7F7F00";
						ctx.beginPath();
						ctx.rect(x1, y + seg_row_hheight * 0.5, len, seg_row_hheight);
						ctx.stroke();
					}
				});
			}

			ctx.restore();
		}
		
		if (viewerState.display_rdna_border) {
			let y2 = seq_list.length * (seg_row_height + seg_row_separate);

			let start = ((1 - bp_start) + options.rDNA_info.alignment_start) * bp_size;
			let end = ((1 - bp_start) + options.rDNA_info.alignment_end) * bp_size;

			let x1 = Math.min(start, end);
			let x2 = Math.max(start, end);

			ctx.strokeStyle = "#7F7F00";
			ctx.beginPath();
			ctx.rect(x1, 0, x2 - x1, y2);
			ctx.stroke();
		}
	}

	// fill line in pixel
	ctx.lineWidth = 1;
	ctx.translate(0.5, 0.5);

	ctx.translate(0, 6 * (seg_row_height + seg_row_separate));

	// user input maker
	{
		const top = ctx.getTransform().f;

		for (let i = 0; i < viewerState._range_makers.length; ++i) {
			let { 0: start, 1: end } = viewerState._range_makers[i];
			let x1 = ((1 - bp_start) + (start - 1)) * bp_size;
			let x2 = ((1 - bp_start) + (end - 1)) * bp_size;
			ctx.beginPath();
			ctx.rect(x1 + 0.5, -top + 0.5, Math.max(1, x2 - x1), top - seg_row_separate);
			ctx.strokeStyle = "#000000DD";
			ctx.stroke();
		}
		for (let i = 0; i < viewerState._position_makers.length; ++i) {
			let pos = viewerState._position_makers[i] - 1;
			let x1 = ((1 - bp_start) + pos) * bp_size;
			ctx.beginPath();
			ctx.rect(x1 + 0.5, -top + 0.5, Math.max(1, bp_size), top - seg_row_separate);
			ctx.strokeStyle = "#000000DD";
			ctx.stroke();
		}
	}

	// crossover

	dataset.crossover_list[viewerState.nChr - 1].forEach(co_data => {
		if (co_data.chr == viewerState.nChr) {
			let x1 = ((1 - bp_start) + ref1_pos_uint32array[co_data.snp_start_out]) * bp_size;
			let x2 = ((1 - bp_start) + ref1_pos_uint32array[co_data.snp_start_in]) * bp_size;
			let x3 = ((1 - bp_start) + ref1_pos_uint32array[co_data.snp_end_in]) * bp_size;
			let x4 = ((1 - bp_start) + ref1_pos_uint32array[co_data.snp_end_out]) * bp_size;

			let in_length = x3 - x2;
			let out_length = x4 - x1;

			if (out_length < 5) {
				let cx = (x4 + x1) / 2;
				let ll = 5;
				ctx.fillStyle = "#FF0000";
				ctx.fillRect(cx - ll / 2, 0, ll, seg_row_height);
			}
			else {
				ctx.fillStyle = "#FF0000CF";
				ctx.fillRect(x1, 0, out_length, seg_row_height);
				ctx.fillStyle = "#FF0000";
				ctx.fillRect(x2, 0, in_length, seg_row_height);
			}
		}
	});

	ctx.translate(0, 1 * (seg_row_height + seg_row_separate));

	// pre-defined maker

	let maker_display_flag = [
		viewerState.display_31, viewerState.display_40,
		viewerState.display_13indel, viewerState.display_22indel, viewerState.display_31indel, viewerState.display_40indel,
		true, viewerState.display_illegitimate_mutation,
	];

	// no display 2:2 SNV
	for (let marker_idx = 1; marker_idx < spore_cmp_array.length; ++marker_idx) {
		let order = marker_idx - 1;
		
		if (maker_display_flag[order]) {
			for (let i = 0; i < spore_cmp_array[marker_idx].length; ++i) {
				let { pos, value } = spore_cmp_array[marker_idx][i];
				let cx = ((1 - bp_start) + pos) * bp_size;
				//let cx = ((1 - bp_start) + pos) * bp_size;
				let x1 = cx - 1;
				let x2 = cx + 1;
				if (!(x2 >= 0 && x1 <= max_view_width)) {
					continue;
				}
				let width = Math.min(x1 + Math.max(bp_size, x2 - x1), max_view_width) - x1;

				let col = marker_order[order];
				
				if (col) {
					ctx.fillStyle = current_colorset[col];
					ctx.fillRect(x1, 0, width, seg_row_height);
				}
			}
		}
		
		ctx.translate(0, 1 * (seg_row_height + seg_row_separate));
	}

	//begin GC%

	if (viewerState.ref1 && viewerState.ref2 && viewerState.nChr) {
		const dataset = [
			{
				refName: viewerState.ref1,
				ref_mapto_ma: ref1_pos_uint32array,
			},
			{
				refName: viewerState.ref2,
				ref_mapto_ma: ref2_pos_uint32array,
			}
		];

		const gc_max_height = seg_row_height * 2;

		dataset.forEach(({ refName, ref_mapto_ma }, sid) => {
			const gc_list = gc_content[refName][viewerState.nChr];
			if (!gc_list || gc_list.length <= 0) {
				return;
			}

			const avg_top = 0;
			const avg_mid = gc_max_height - gc_max_height * 0.5;
			const avg_bottom = gc_max_height;

			ctx.beginPath();
			ctx.moveTo(0, avg_mid);
			ctx.lineTo(max_view_width, avg_mid);

			ctx.strokeStyle = "darkgray";
			//ctx.setLineDash([2, 1]);
			ctx.stroke();
			//ctx.setLineDash([0, 0]);//reset

			ctx.fillStyle = "black";
			ctx.textAlign = "left";
			ctx.font = Math.trunc(gc_max_height / 3 / 1.5) + "px Arial";

			ctx.textBaseline = "middle";
			ctx.fillText((50).toFixed(1), max_view_width + 5, avg_mid);//50%
			
			ctx.textBaseline = "top";
			ctx.fillText((100).toFixed(1), max_view_width + 5, avg_top);//100%
			
			ctx.textBaseline = "bottom";
			ctx.fillText((0).toFixed(1), max_view_width + 5, avg_bottom - 1);//0%

			if (viewerState.gc_content_clip_indel) {
				ctx.save();//begin clip
				//
				ctx.beginPath();
				{//remove indel
					ctx.fillStyle = "#FFFFFF";
					let ss = region_rect[sid];
					for (let si = 0; si < ss.length; ++si) {
						let color_id = ss[si].col;
						if (color_id & ColorID.indel_mask) {
						}
						else {
							let start = ((1 - bp_start) + ss[si].start) * bp_size;
							let end = ((1 - bp_start) + ss[si].end) * bp_size;
							if (end >= 0 && start <= max_view_width) {
								let row_y_bottom = 1 * (gc_max_height + seg_row_separate);
								let row_y_top = row_y_bottom - gc_max_height;
								ctx.rect(start, row_y_top, end - start, gc_max_height);
								//ctx.fill();
							}
						}
					}
				}//gc remove indel
				ctx.clip();
			}

			let last_x = max_view_width;
			ctx.beginPath();
			ctx.moveTo(0, gc_max_height);
			for (let i = 0; i < gc_list.length; ++i) {
				let data = gc_list[i];
				let start = ref_mapto_ma[data.start];
				let end = ref_mapto_ma[data.end];
				if (start == null || end == null ) {
					break;
				}
				let x1 = ((1 - bp_start) + start) * bp_size;
				let x2 = ((1 - bp_start) + end) * bp_size;
				if (x2 >= 0 && x1 <= max_view_width) {
					x1 = Math.max(0, x1);
					x2 = Math.min(x2, max_view_width);
					let w = x2 - x1;

					//let cx = (end + start) >> 1;//=> Math.trunc((end + start) / 2)
					let gc = (data.gc || 0) * gc_max_height / 100;
					let y = avg_bottom - gc;

					ctx.lineTo(x1, y);
					ctx.lineTo(x2, y);
					last_x = x2;
				}
			}
			ctx.lineTo(last_x, gc_max_height);
			ctx.strokeStyle = "black";
			ctx.stroke();
			ctx.fillStyle = viewerState.gc_gap_color || "#FFFF0033";//"#FFFF0033"
			ctx.fill();

			if (viewerState.gc_content_clip_indel) {
				ctx.restore();//end clip => clear clip
			}

			ctx.textAlign = "left";//reset
			ctx.textBaseline = "alphabetic";//reset

			if (!viewerState.gc_content_clip_indel) {
				ctx.fillStyle = "#FFFFFF22";
				let ss = region_rect[sid];
				for (let si = 0; si < ss.length; ++si) {
					let start = ((1 - bp_start) + ss[si].start) * bp_size;
					let end = ((1 - bp_start) + ss[si].end) * bp_size;
					if (end >= 0 && start <= max_view_width) {
						let color_id = ss[si].col;
						if (color_id & ColorID.indel_mask) {
							let row_y_bottom = offset_y + (sid + 1) * (gc_max_height + seg_row_separate);
							let row_y_top = row_y_bottom - gc_max_height;
							ctx.beginPath();
							ctx.rect(start, row_y_top, end - start, gc_max_height);
							ctx.fill();
						}
						else {
						}
					}
				}
			}
			
			//TDOD: GC% plot border rectangle
			ctx.strokeStyle = "black";
			ctx.beginPath();

			ctx.rect(0, offset_y, max_view_width, gc_max_height);

			//ctx.moveTo(0, offset_y + sid * height);
			//ctx.lineTo(max_view_width, offset_y + sid * height);

			ctx.stroke();

			ctx.translate(0, 1 * (gc_max_height + seg_row_separate));
		});
	}//gc
}

function parseFasta(in_seq) {
	let all = in_seq.split(">");
	
	let results = {};
	
	all.filter(a => a.length).forEach((sub_fa) => {
		let li = sub_fa.indexOf("\n");
		let out_name = sub_fa.slice(0, li).trim().replace(/:/g, "");
		let out_file_name = out_name.match(/^([^ ]+)/)[1];
		let out_seq = sub_fa.slice(li).trim();
		
		out_seq = out_seq.replace(/\n/g, "");
		out_seq = out_seq.replace(/\r/g, "");
		out_seq = out_seq.replace(/ /g, "");
		
		results[out_file_name] = out_seq.toUpperCase();

		// @ts-ignore
		console.timeLog("loading", "seq:", out_file_name);
	});
	
	return results;
}

/**
 * @param {string} url
 * @param {""|"arraybuffer"|"blob"|"document"|"json"|"text"} responseType
 */
function fetchData(url, responseType) {
	return new Promise(function (resolve, reject) {
		let xhr = new XMLHttpRequest();
		xhr.open("GET", url, true);

		if (responseType) {
			xhr.responseType = responseType;;
		}

		xhr.timeout = 10 * 60 * 1000;//20000;

		xhr.onload = function () {
			if (this.status == 404 || this.status == 500) {
				alert("file: " + url);
				debugger;
				//resolve(null);
				reject(this.status + ": " + url);
			}
			else if (this.status == 200) {
				resolve(this.response);
			}
			else if (this.status == 304) {
				debugger
			}
		};

		xhr.ontimeout = function (e) {
			debugger;
			//resolve(null);
			reject("timeout: " + url);
		};

		xhr.onabort = function (e) {
			reject("abort: " + url);
		};

		xhr.send();
	});
}

