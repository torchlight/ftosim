/* ftosim.js - A face-turning octahedron twisty puzzle simulator
version 0.1 (2021-03-12)

Copyright (c) 2021 torchlight

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*

Internally, we treat the octahedron as being axis-aligned.

The axis convention used for 3D coordinates are right-handed:
positive x goes to the right of the screen
positive y goes to the top of the screen (NOTE: this is not the same as screen coordinates)
positive z goes out of the screen


Looking at the screen:

     +y

    /|\
   / | \
  /  |  \
 / U | R \
/___ |___ \  +x
\    |    /
 \ L | F /
  \  |  /
   \ | /
    \|/

Looking "through" the screen:

     +y

    /|\
   / | \
  /  |  \
 / B |BR \
/___ |___ \  +x
\    |    /
 \ BL| D /
  \  |  /
   \ | /
    \|/

The facelets are stored as a triple of points, with each point having integer coordinates.
In particular, |x| + |y| + |z| = size, always.

U face:   - x + y + z = size
R face:     x + y + z = size
F face:     x - y + z = size
L face:   - x - y + z = size
B face:   - x + y - z = size
BR face:    x + y - z = size
D face:     x - y - z = size
BL face:  - x - y - z = size

For each facelet, other than the logical 3D coordinates (which are not animated), we also store:

>>  previous 3D coordinates
>>  list of 3D rotations to apply (possibly empty)
    >>  each rotation is specified together with start and end times;
        rotations are specified as axis-angle
        rotations that are completed are applied "for real" then removed from this list
>>  intermediate 3D coordinates after applying the 3D rotations _in order_
>>  an SVG <g> element containing one or more <polygon>s
    corresponding to what to draw on the screen
>>  the colour (as a face label)
*/


'use strict';

const SVGNS = 'http://www.w3.org/2000/svg';

// cribbed from Aedan's FTOSIM v1.1
let colour_scheme = {
U:  '#FFFFFF',
F:  '#FF0000',
R:  '#00ED04',
L:  '#B515DC',
B:  '#0E1CFF',
D:  '#F6FD01',
BR: '#777777',
BL: '#FFA802',
};

const AXIS_VECTORS = [[-1,1,1], [-1,-1,1], [1,-1,1], [1,1,1]];
Object.freeze(AXIS_VECTORS);
AXIS_VECTORS.forEach(x => Object.freeze(x));

function cross(u, v)
{
	// cross product of two vectors
	return [u[1]*v[2] - u[2]*v[1],
	        u[2]*v[0] - u[0]*v[2],
	        u[0]*v[1] - u[1]*v[0]];
}

function dot(u, v)
{
	// dot product of two vectors
	let s = 0;
	for (let i = 0; i < u.length; i++) {s += u[i] * v[i];}
	return s;
}

function rotate(v, axis, angle)
{
	// axis must be a unit vector (we don't check here)
	let c = Math.cos(angle), s = Math.sin(angle);
	let cp = cross(axis, v);
	let sc = dot(axis, v) * (1-c);
	let rotated = [];
	for (let i = 0; i < 3; i++)
	{
		rotated[i] = v[i]*c + cp[i]*s + axis[i]*sc;
	}
	return rotated;
}

function calc_poly_normal(vertices)
{
	let normal = [0, 0, 0];
	let n = vertices.length;
	for (var i = 0; i < n; i++)
	{
		let tnormal = cross(vertices[i], vertices[(i+1)%n]);
		normal[0] += tnormal[0];
		normal[1] += tnormal[1];
		normal[2] += tnormal[2];
	}
	return normal;
}

class Facelet
{
	coords_list; // always integer
	old_coords_list; // always integer
	inter_coords_list; // might not be integer
	animations;
	g_el;
	value;
	needs_redraw = true;
	__verbose = false;

	constructor(coords_list, value)
	{
		this.coords_list = Array.from(coords_list);
		this.old_coords_list = this.coords_list.slice();
		this.inter_coords_list = this.coords_list.slice();
		this.animations = [];
		this.g_el = document.createElementNS(SVGNS, 'g');
		this.value = value;
	}

	update_coords(time)
	{
		let animations = this.animations;
		if (animations.length === 0)
		{
			// no animation; nothing to update
			return;
		}
		this.needs_redraw = true;
		if (this.__verbose) {console.log('needs redraw');}
		let num_expired = 0;
		this.inter_coords_list = this.old_coords_list.slice();
		for (let i = 0; i < animations.length; i++)
		{
			let animation = animations[i];
			for (let j = 0; j < this.inter_coords_list.length; j++)
			{
				this.inter_coords_list[j] = animation.transform(this.inter_coords_list[j], time);
			}
			if (time >= animation.end && i === num_expired)
			{
				if (this.__verbose)
				{
					console.log(`time: ${time}\nend: ${animation.end}`);
				}
				// all animations up to and including this one have finished
				num_expired++;
				// round the coordinates to integers to avoid accumulating float error
				for (let j = 0; j < this.inter_coords_list.length; j++)
				{
					for (let k = 0; k < 3; k++)
					{
						this.inter_coords_list[j][k] = Math.round(this.inter_coords_list[j][k]);
					}
				}
				this.old_coords_list = this.inter_coords_list.slice();
			}
		}
		animations.splice(0, num_expired);
	}

	draw(force=false)
	{
		if (!this.needs_redraw && !force) {return;}
		const inset = 0.1;
		let shrunk_coords_list = shrink_triangle(this.inter_coords_list, inset);
		let polys = [shrunk_coords_list];
		for (let axis_index = 0; axis_index < 3; axis_index++)
		{
			polys = polys.map(poly => split_polygon_by_axial_planes(poly, axis_index)).flat();
		}
		let polys_needed = polys.length;
		let g_el = this.g_el;
		let poly_els = Array.from(g_el.querySelectorAll('polygon'));
		while (poly_els.length < polys_needed)
		{
			let poly_el = document.createElementNS(SVGNS, 'polygon');
			poly_el.setAttribute('fill', colour_scheme[this.value] ?? 'none');
			poly_el.setAttribute('stroke', 'none');
			g_el.appendChild(poly_el);
			poly_els.push(poly_el);
		}
		if (this.__verbose)
		{
			console.log(polys);
			console.log(polys_needed);
			console.log(poly_els);
		}
		for (let i = 0; i < polys_needed; i++)
		{
			let poly = polys[i];
			let points = [];
			for (let j = 0; j < polys[i].length; j++)
			{
				points[j] = oct_project(polys[i][j]);
			}
			poly_els[i].setAttribute('points', points.map(p => p.map(x => x.toFixed(5)).join(',')).join(' '));
		}
		for (let i = polys_needed; i < poly_els.length; i++)
		{
			poly_els[i].setAttribute('points', '0,0');
		}
		this.needs_redraw = false;
	}
}

function shrink_triangle([u0, u1, u2], inset)
{
	let centre = [0, 0, 0];
	for (let i = 0; i < 3; i++)
	{
		centre[i] = (u0[i] + u1[i] + u2[i]) / 3;
	}
	return [u0, u1, u2].map(u => {
		let v = [];
		for (let i = 0; i < 3; i++)
		{
			v[i] = u[i]*(1-inset) + centre[i]*inset;
		}
		return v;
	});
}

function oct_project(vector)
{
	// octahedral projection
	let [x, y, z] = vector;
	let l1 = Math.abs(x) + Math.abs(y) + Math.abs(z);
	x /= l1;
	y /= l1;
	z /= l1;
	let xx = x+y, yy = x-y;
	if (z >= 0)
	{
		// in the U-R-F-L hemisphere
		return [xx, yy];
		// this is just a reflection and a scale
	}
	else if (x >= 0)
	{
		// BR octant
		if (y >= 0) {return [2-xx, yy];}
		// D octant
		else {return [xx, 2-yy];}
	}
	else // x < 0
	{
		// B octant
		if (y >= 0) {return [xx, -2-yy];}
		// BL octant
		else {return [-2-xx, yy];}
	}
}

function split_polygon_by_axial_planes(coords_list, axis_index)
{
	// axis_index = 0 ~ split by y-z plane
	//              1 ~ split by z-x plane
	//              2 ~ split by x-y plane
	const EPS = 1e-8;
	let plus = [];
	let minus = [];
	let n = coords_list.length;
	for (let i = 0; i < n; i++)
	{
		let u = coords_list[i];
		let v = coords_list[(i+1)%n];
		if (u[axis_index] >= 0) {plus.push(u);}
		else {minus.push(u);}
		if ((u[axis_index] >= 0 && v[axis_index] < 0) || (u[axis_index] < 0 && v[axis_index] >= 0))
		{
			// solve for where the edge intersects the plane
			let t = v[axis_index]/(v[axis_index] - u[axis_index]);
			let cut = [t*u[0] + (1-t)*v[0], t*u[1] + (1-t)*v[1], t*u[2] + (1-t)*v[2]];
			// add a bit of fudge to ensure that it has the correct sign
			cut[axis_index] = EPS;
			plus.push(cut.slice());
			cut[axis_index] = -EPS;
			minus.push(cut.slice());
		}
	}
	let out = [];
	if (plus.length > 1) {out.push(plus);}
	if (minus.length > 1) {out.push(minus);}
	return out;
}

class Animation
{
	axis; // unit vector
	angle; // angle in radians (>2pi allowed; that means rotating multiple times)
	start;
	end;

	constructor(axis, angle, start, end)
	{
		let l = dot(axis, axis) ** 0.5;
		this.axis = axis.map(x => x/l);
		this.angle = angle;
		this.start = start;
		this.end = end;
	}

	transform(vector, time)
	{
		let t = (time - this.start) / (this.end - this.start);
		if (time < this.start) {t = 0;}
		else if (time > this.end) {t = 1;}
		t = ease(t);
		let axis = this.axis;
		let angle = this.angle * t;
		return rotate(vector, axis, angle);
	}
}

function ease(x)
{
	if (x < 0) {return 0;}
	if (x > 1) {return 1;}
	return x * x * (3 - 2 * x);
}

class FTOSim
{
	size; // number of layers; the standard FTO has three layers
	num_facelets;
	facelets;
	last_updated_time;
	time_per_move = 180; // how long the move animation lasts in milliseconds
	queued_animations;

	constructor(size=3)
	{
		this.size = size;
		this.num_facelets = 8 * size**2;
		let facelets = [];
		// populate U facelets first
		for (let i = 0; i < size; i++)
		{
			let first = new Facelet([[0, i, size-i], [0, i+1, size-i-1], [-1, i, size-i-1]], 'U');
			facelets.push(first);
			// this loop is nontrivial iff i >= 1
			for (let j = 0; j < i; j++)
			{
				let down = new Facelet([[-j-1, i-j-1, size-i], [-j, i-j, size-i], [-j-1, i-j, size-i-1]], 'U');
				let up = new Facelet([[-j-1, i-j-1, size-i], [-j-1, i-j, size-i-1], [-j-2, i-j-1, size-i-1]], 'U');
				facelets.push(down, up);
			}
		}
		// copy them into F facelets
		facelets = facelets.concat(facelets.map(facelet => {
			let coords_list = facelet.coords_list;
			return new Facelet(coords_list.map(([x, y, z]) => [-x, -y, z]), 'F');
		}));
		// copy U/F into BR/BL
		facelets = facelets.concat(facelets.map(facelet => {
			let coords_list = facelet.coords_list;
			let value = {'U': 'BR', 'F': 'BL'}[facelet.value];
			return new Facelet(coords_list.map(([x, y, z]) => [-x, y, -z]), value);
		}));
		// copy U/F/BR/BL into L/R/D/B
		facelets = facelets.concat(facelets.map(facelet => {
			let coords_list = facelet.coords_list;
			let value = {'U': 'L', 'F': 'R', 'BR': 'D', 'BL': 'B'}[facelet.value];
			return new Facelet(coords_list.map(([x, y, z]) => [x, -y, z]).reverse(), value);
			// we reverse here to keep the normals of the facelets pointing outwards
			// (per right-hand rule)
		}));
		this.facelets = facelets;
		this.queued_animations = [];
	}

	attach(svg_el)
	{
		// this can be any "container" element that's in an <svg> element, not necessarily an
		// <svg> element itself
		for (let facelet of this.facelets)
		{
			if (facelet.g_el.parentNode !== svg_el)
			{
				svg_el.appendChild(facelet.g_el);
			}
		}
	}
	
	detach()
	{
		for (let facelet of this.facelets)
		{
			if (facelet.g_el.parentNode !== null)
			{
				facelet.g_el.parentNode.removeChild(facelet.g_el);
			}
		}
	}

	update_and_draw(time=performance.now(), force=false)
	{
		if (this.queued_animations.length === 0 && !force) {return;}
		for (let facelet of this.facelets)
		{
			facelet.update_coords(time);
			facelet.draw();
		}
		// clear expired animations
		while (this.queued_animations.length > 0)
		{
			//console.log(this.queued_animations[0].end, time);
			if (this.queued_animations[0].end < time) {this.queued_animations.shift();}
			else {break;}
		}
		this.last_updated_time = time;
	}

	apply_block_move(
		face,
		start_depth,
		end_depth,
		amount, // clockwise
		start_time=performance.now(),
		end_time=start_time+this.time_per_move
		)
	{
		/*
		the depths here refer to cutting depth, so e.g.
		start=0, end=1 corresponds to a normal single-layer turn,
		start=0, end=2 corresponds to a wide turn,
		start=0, end=size corresponds to a whole rotation,
		start=1, end=2 corresponds to an inner-slice turn.
		*/
		const size = this.size;
		let face_index = ['U', 'L', 'F', 'R', 'D', 'BR', 'B', 'BL'].indexOf(face);
		if (face_index === -1) {throw 'invalid move face';}
		if (!(0 <= start_depth && start_depth <= end_depth && end_depth <= size))
		{
			throw 'invalid move depth(s)';
		}
		if (amount !== Math.floor(amount) || Math.abs(amount) > 100)
		{
			throw 'invalid move amount';
		}
		if (face_index >= 4)
		{
			// rewrite D/BR/B/BL moves as moves on the opposite face
			[start_depth, end_depth] = [size - end_depth, size - start_depth];
			face_index -= 4;
			amount = -amount;
		}
		if (start_depth === end_depth)
		{
			/*
			we need to special-case this because otherwise we could have
			start=0, end=0 moving only the facelets on that face, which is not a legal move
			*/
			return;
		}
		const axis_vector = AXIS_VECTORS[face_index];
		const unit_axis_vector = axis_vector.map(x => x * 0.5773502691896257);
		const start = size - 2*end_depth;
		const end = size - 2*start_depth;
		const angle = -amount*2*Math.PI/3;
		const animation = new Animation(unit_axis_vector, angle, start_time, end_time);
		this.queued_animations.push(animation);
		const num_facelets = this.num_facelets;
		for (let i = 0; i < num_facelets; i++)
		{
			let facelet = this.facelets[i];
			let coords_list = facelet.coords_list;
			if (coords_list.map(v => dot(v, axis_vector)).every(x => (start <= x && x <= end)))
			{
				facelet.animations.push(animation);
				facelet.coords_list = coords_list.map(v => rotate(v, unit_axis_vector, angle).map(x => Math.round(x)));
				if (facelet.g_el.parentNode !== null)
				{
					/* move the <g> element to the end so it renders on top
					(SVG does not let you adjust z-index directly, not even through CSS)
					((the actually correct layering requires computing occlusions or
					whatever, which I don't know how to do, but this is a good-enough
					approximation))
					(((this is actually also the "wrong" place to adjust the layering
					because we should be adjusting it during the animation, not at
					the start. TODO figure this out maybe)))
					*/
					facelet.g_el.parentNode.appendChild(facelet.g_el);
				}
			}
		}
	}

	apply_full_puzzle_rotation(axis_vector, angle, start_time=performance.now(), end_time=start_time+this.time_per_move)
	{
		const axis_vector_length = dot(axis_vector, axis_vector)**0.5;
		const unit_axis_vector = axis_vector.map(x => x / axis_vector_length);
		for (let a = 0; a < 3; a++)
		{
			let v = [0, 0, 0];
			v[a] = this.size;
			let vv = rotate(v, unit_axis_vector, angle).map(x => Math.abs(x));
			vv.sort((x, y) => x-y);
			if (vv[0] > 0.49 || vv[1] > 0.49 || Math.abs(vv[2]-this.size) > 0.49)
			{
				throw 'invalid rotation (does not preserve axis alignment)';
			}
		}
		const animation = new Animation(unit_axis_vector, angle, start_time, end_time);
		this.queued_animations.push(animation);
		const num_facelets = this.num_facelets;
		for (let facelet of this.facelets)
		{
			facelet.animations.push(animation);
			facelet.coords_list = facelet.coords_list.map(v => rotate(v, unit_axis_vector, angle).map(x => Math.round(x)));
		}
	}

	apply_move_sequence_from_string(move_sequence_string, start_time=performance.now())
	{
		/*
		basic SiGN + Ben's notation is supported:

		U, L, F, R, D, BR, B, BL (outer layer)

		u, l, f, r, d, br, b, bl (two layers)
		Uw, Lw, etc. (two layers)

		nU, nL, etc. (nth slice only; 1 <= n <= size)
		nUw, nLw, etc. (n layers; 1 <= n <= size-1)
		(XXX currently the bounds aren't checked and `apply_block_move` has laxer requirements)

		Us, Ls, etc. (second slice only)

		Uo, Lo, etc. (whole puzzle rotation)
		Uv, Lv, etc. (whole puzzle rotation)
		T (rotate around U-R-F-L vertex)

		. (no-op; no suffixes allowed)

		// (single-line comment, but see below)

		numerical prefixes/suffixes must all be written in decimal without any leading zeros and
		without a decimal point. (XXX not currently fully enforced)

		time markers in comments:
		>> time must be specified in seconds; must be nonnegative
		>> no extraneous leading zeros allowed; scientific notation not allowed
		>> the time marker applies to all the moves on that line
		    (XXX does this make sense as a default behaviour?)
		>> the first line with a time marker is the zero-point
		>> there must be at most one time marker per line
		>> if there are multiple time markers, a parse error will be thrown
		    (nb: this means that not all single-line comments are valid here!)
		>> there is currently no notation for the end time of a move
		>> times currently do not have to be monotonous; however, the end state of the puzzle
		   will always have the moves applied in written order, rather than chronological order.
		>> lines without time markers specified will be interpolated (the exact interpolation
		   used is subject to change)

		acceptable whitespace:
		U+0009 Character tabulation (tab)
		U+000A Line feed (Unix newline)
		U+0020 Space

		NOT SUPPORTED:
		multi-line comments (slash-star star-slash)
		brackets
		commutator/conjugate notation
		*/
		let lines = move_sequence_string.replace(/[\t ]+/g, ' ').split(/\n+/g);
		let move_list = [];
		let time_list = [];
		for (let line of lines)
		{
			if (line.includes('//'))
			{
				let [moves, ...comments] = line.split(/ *\/\/ */);
				let markers = comments.filter(s => /^(0|[1-9][0-9]*)(\.[0-9]*)?$/.test(s));
				if (markers.length > 1) {throw 'multiple time markers in a line';}
				else if (markers.length === 0)
				{
					if (moves !== '')
					{
						move_list = move_list.concat(moves.split(' '));
						time_list.length = move_list.length;
					}
				}
				else /* markers.length === 1 */
				{
					let time = parseFloat(markers[0]);
					if (moves === '')
					{
						move_list.push('.');
						time_list.push(time);
					}
					else
					{
						let moves_split = moves.split(' ');
						move_list = move_list.concat(moves_split);
						time_list = time_list.concat(Array(moves_split.length).fill(time));
						// XXX does this make sense as a default behaviour?
					}
				}
			}
			else
			{
				if (line !== '')
					{
						move_list = move_list.concat(line.split(' '));
						time_list.length = move_list.length;
					}
			}
		}
		//console.log(move_list, time_list);
		let zero_point = time_list.find(x => x !== undefined);
		let zero_point_index = time_list.findIndex(x => x !== undefined);
		if (zero_point === undefined)
		{
			// none of the lines have time markers
			for (let i = 0; i < time_list.length; i++)
			{
				time_list[i] = i * this.time_per_move;
			}
		}
		else
		{
			// set all times before the first given time marker to that
			for (let i = 0; i < zero_point_index; i++)
			{
				time_list[i] = zero_point;
			}
			// interpolate all missing times up to the last given time marker
			let prev_index = zero_point_index;
			while (prev_index < time_list.length-1)
			{
				let next_index = time_list.findIndex((x, i) => x !== undefined && i > prev_index);
				if (next_index === -1) {break;} // no more time markers found
				let n = next_index - prev_index;
				let delta = (time_list[next_index] - time_list[prev_index]) / n;
				for (let i = prev_index+1; i < next_index; i++)
				{
					time_list[i] = (i - prev_index) * delta + time_list[prev_index];
				}
				prev_index = next_index;
			}
			// extrapolate all times beyond the last given time marker
			for (let i = prev_index+1; i < time_list.length; i++)
			{
				time_list[i] = time_list[i-1] + this.time_per_move;
			}
			// zero the zero point
			for (let i = 0; i < time_list.length; i++)
			{
				time_list[i] -= zero_point;
			}
		}
		for (let i = 0; i < move_list.length; i++)
		{
			let move_string = move_list[i];
			if (move_string === '.')
			{
				// this is a no-op
				continue;
			}
			// for all other types of moves, we allow suffixes
			let match = move_string.match(/^(?<prefix>(\d+-\d+)|(\d+)|())(?<move>[^0-9']+)(?<suffix>\d*'?)$/);
			if (match === null) {throw `invalid move ${move_string} (unrecognised format)`;}
			let {prefix, move, suffix} = match.groups;
			let face;
			let start_depth, end_depth;
			let amount;
			if (suffix === '') {amount = 1;}
			else if (suffix === "'") {amount = -1;}
			else if (!suffix.includes("'")) {amount = parseInt(suffix, 10);}
			else {amount = -parseInt(suffix, 10);}
			// note that parseInt strips non-digits at the end
			switch (move)
			{
				// normal moves / SiGN single-slice moves
				case 'U':
				case 'L':
				case 'F':
				case 'R':
				case 'D':
				case 'BR':
				case 'B':
				case 'BL':
					face = move;
					if (prefix === '') {start_depth = 0; end_depth = 1;}
					else if (prefix.includes('-'))
					{
						throw `invalid move ${move_string} (range prefix not allowed for SiGN slice)`;
					}
					else
					{
						let a = parseInt(prefix, 10);
						start_depth = a-1;
						end_depth = a;
					}
					break;
				// Ben notation slice moves
				case 'Us':
				case 'Ls':
				case 'Fs':
				case 'Rs':
				case 'Ds':
				case 'BRs':
				case 'Bs':
				case 'BLs':
					face = move.replace(/s/, '');
					if (prefix === '') {start_depth = 1; end_depth = 2;}
					else if (prefix.includes('-'))
					{
						throw `invalid move ${move_string} (range prefix not allowed for Ben slice)`;
					}
					else
					{
						/*
						let a = parseInt(prefix, 10);
						start_depth = a-1;
						end_depth = a;
						// XXX should we allow e.g. `3Us` to turn the third slice?
						// this is redundant with SiGN's `3U`.
						*/
						throw `invalid move ${move_string} (index prefix not allowed for Ben slice)`;
					}
					break;
				// SiGN wide moves (this is supported by a.c.n v2)
				case 'u':
				case 'l':
				case 'f':
				case 'r':
				case 'd':
				case 'br':
				case 'b':
				case 'bl':
				// Ben notation / WCA-style wide moves (this is not supported by a.c.n v2)
				case 'Uw':
				case 'Lw':
				case 'Fw':
				case 'Rw':
				case 'Dw':
				case 'BRw':
				case 'Bw':
				case 'BLw':
					face = move.replace(/w/, '').toUpperCase();
					if (prefix === '') {start_depth = 0; end_depth = 2;}
					else if (prefix.includes('-'))
					{
						throw `invalid move ${move_string}`;
						let [a, b] = prefix.split('-').map(s => parseInt(s, 10));
						if (a > b) {throw `invalid move ${move_string} (range prefix in wrong order)`;}
						// a.c.n behaves weirdly on e.g. 5-3Uw; we don't want to
						// replicate that behaviour.
						start_depth = a-1;
						end_depth = b;
					}
					else
					{
						let a = parseInt(prefix, 10);
						start_depth = 0;
						end_depth = a;
					}
					break;
				// Ben notation face rotations
				case 'Uo':
				case 'Lo':
				case 'Fo':
				case 'Ro':
				case 'Do':
				case 'BRo':
				case 'Bo':
				case 'BLo':
				// Twizzle notation face rotations
				case 'Uv':
				case 'Lv':
				case 'Fv':
				case 'Rv':
				case 'Dv':
				case 'BRv':
				case 'Bv':
				case 'BLv':
					face = move.replace(/o|v/, '');
					if (prefix !== '') {throw `invalid move ${move_string} (prefix not allowed)`;}
					start_depth = 0;
					end_depth = this.size;
					break;
				// Ben notation vertex rotations
				case 'T':
					if (prefix !== '') {throw `invalid move ${move_string} (prefix not allowed)`;}
					break;
				// XXX support Twizzle notation vertex rotations?
				// and maybe edge rotations too?
				default:
					throw `invalid move ${move_string}`;
			}
			let move_start_time = start_time + time_list[i];
			let move_end_time = move_start_time + this.time_per_move;
			if (move === 'T')
			{
				this.apply_full_puzzle_rotation([0, 0, 1], -amount*Math.PI/2, move_start_time, move_end_time);
			}
			else
			{
				this.apply_block_move(face, start_depth, end_depth, amount, move_start_time, move_end_time);
			}
		}
	}

	check_solved()
	{
		/* puzzle solved iff all facelets on a given face have the same colour;
		we determine which face a facelet lies on by calculating its normal vector.
		the sign convention used allows us to differentiate opposite faces since
		the normal vector chosen always points outwards.
		*/
		let values = {};
		for (let facelet of this.facelets)
		{
			let normal = calc_poly_normal(facelet.coords_list).join(',');
			if (!(normal in values))
			{
				values[normal] = facelet.value;
			}
			if (values[normal] !== facelet.value) {return false;}
		}
		return true;
	}
}
