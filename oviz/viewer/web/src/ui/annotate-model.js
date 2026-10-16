// Annotations: the data model, curve and shell geometry, and the grouping
// rules, kept free of the DOM and GL so Node can test them.
//
// A document is { groups, items }. Every item belongs to one group, and the
// groups are what the key lists (one switch each). Positions are points:
// { pos: [x, y, z] } fixed in the scene's frame, or { trace, index, offset }
// riding on an object, so they follow its orbit through time.
//
//   text   { at, text, size (px), weight, color, opacity, bg, dx, dy }
//   curve  { points: [point…], smooth (0–1), closed, width (px), color, opacity, dash, arrow }
//   sphere { center, radii [pc ×3], rot [deg ×3], style: bubble | shell | wire, color, opacity }
//   box    { center, radii (half sizes, pc ×3), rot, style: box | wire, color, opacity }
//
// Any item may set `present: true` to show only around the present day,
// like the dust maps it often marks. Shapes and curves may set `select`
// ("highlight" or "isolate") to select the objects inside them (a curve:
// within `reach` pc of it), coloured like the shape or shown alone.

import { clamp, lerp, DEG } from "../core/math.js";
import { parseColor } from "../core/color.js";

export const SPHERE_STYLES = ["bubble", "shell", "wire"];
export const BOX_STYLES = ["box", "wire"];
export const SELECT_MODES = ["highlight", "isolate"];
/** The twelve edges of a box, as pairs of corner indices (see boxCorners). */
export const BOX_EDGES = [[0, 1], [2, 3], [4, 5], [6, 7], [0, 2], [1, 3], [4, 6], [5, 7], [0, 4], [1, 5], [2, 6], [3, 7]];
export const CURVE_DASHES = ["solid", "dash", "dot"];
export const CURVE_ARROWS = ["none", "end", "start", "both"];

/** Default look of what each tool draws. */
export const ANNOT_DEFAULTS = {
  text: { size: 15, weight: 600, color: "#ffffff", opacity: 1, bg: false, dx: 0, dy: 0 },
  curve: { smooth: 1, closed: false, width: 2, color: "#ffd27a", opacity: 1, dash: "solid", arrow: "none" },
  arrow: { smooth: 1, closed: false, width: 2.5, color: "#ffd27a", opacity: 1, dash: "solid", arrow: "end" },
  bubble: { style: "bubble", color: "#7cc4ff", opacity: 0.6, rot: [0, 0, 0] },
  shell: { style: "shell", color: "#ff9e7a", opacity: 0.9, rot: [0, 0, 0] },
  box: { style: "box", color: "#9b7bff", opacity: 0.85, rot: [0, 0, 0] },
};

export function emptyAnnotations() {
  return { groups: [], items: [] };
}

let counter = 0;
/** A short unique id for groups and items. */
export function annotId(prefix = "a") {
  counter = (counter + 1) % 1e6;
  return `${prefix}${Date.now().toString(36).slice(-5)}${Math.random().toString(36).slice(2, 6)}${counter.toString(36)}`;
}

const num = (v, d) => (Number.isFinite(Number(v)) ? Number(v) : d);
const str = (v, d) => (typeof v === "string" && v.trim() ? v : d);
const vec3 = (v) => (Array.isArray(v) && v.length >= 3 && v.slice(0, 3).every((x) => Number.isFinite(Number(x))) ? v.slice(0, 3).map(Number) : null);

/** A point record, or null if it names no usable position. */
export function normalizePoint(p) {
  if (!p) return null;
  if (Array.isArray(p)) { const pos = vec3(p); return pos ? { pos } : null; }
  if (typeof p !== "object") return null;
  if (p.trace != null && Number.isInteger(Number(p.index)) && Number(p.index) >= 0) {
    const out = { trace: String(p.trace), index: Number(p.index) };
    const off = vec3(p.offset);
    if (off && off.some((x) => x !== 0)) out.offset = off;
    return out;
  }
  const pos = vec3(p.pos ?? p.position ?? [p.x, p.y, p.z]);
  return pos ? { pos } : null;
}

/** Validate a document (from a file, a link, a State or Python) into a clean copy. */
export function normalizeAnnotations(raw) {
  const doc = emptyAnnotations();
  if (!raw || typeof raw !== "object") return doc;
  const seen = new Set();
  const fresh = (id, prefix) => {
    let s = typeof id === "string" && id ? id : annotId(prefix);
    while (seen.has(s)) s = annotId(prefix);
    seen.add(s);
    return s;
  };
  for (const g of Array.isArray(raw.groups) ? raw.groups : []) {
    if (!g || typeof g !== "object") continue;
    const group = { id: fresh(g.id, "g"), name: String(g.name ?? ""), auto: g.auto !== false, visible: g.visible !== false };
    if (typeof g.nameFrom === "string") group.nameFrom = g.nameFrom;
    doc.groups.push(group);
  }
  const groupIds = new Set(doc.groups.map((g) => g.id));
  for (const it of Array.isArray(raw.items) ? raw.items : []) {
    const item = normalizeItem(it);
    if (!item) continue;
    item.id = fresh(it.id, "i");
    if (!groupIds.has(item.group)) {
      // An item names a group that is missing: give it one of its own.
      const g = { id: typeof it.group === "string" && it.group && !seen.has(it.group) ? fresh(it.group, "g") : fresh(null, "g"), name: "", auto: true, visible: true };
      doc.groups.push(g);
      groupIds.add(g.id);
      item.group = g.id;
    }
    doc.items.push(item);
  }
  // Groups left empty are dropped, and a name taken from a removed label goes.
  const used = new Set(doc.items.map((i) => i.group));
  const ids = new Set(doc.items.map((i) => i.id));
  doc.groups = doc.groups.filter((g) => used.has(g.id));
  for (const g of doc.groups) if (g.nameFrom && !ids.has(g.nameFrom)) delete g.nameFrom;
  return doc;
}

function normalizeItem(it) {
  if (!it || typeof it !== "object") return null;
  const kind = String(it.kind || "");
  const base = {
    id: "", group: typeof it.group === "string" ? it.group : "",
    color: str(it.color, null), opacity: clamp(num(it.opacity, 1), 0, 1),
  };
  if (it.present) base.present = true;
  if (it.visible === false) base.visible = false;
  if (kind !== "text" && SELECT_MODES.includes(it.select)) base.select = it.select;
  if (kind === "text") {
    const at = normalizePoint(it.at ?? it.anchor ?? it.position);
    if (!at) return null;
    const d = ANNOT_DEFAULTS.text;
    return {
      ...base, kind, at, text: String(it.text ?? "Label").slice(0, 400),
      size: clamp(num(it.size, d.size), 6, 96), weight: num(it.weight, d.weight) >= 550 ? 600 : 400,
      color: base.color || d.color, bg: !!it.bg, dx: clamp(num(it.dx, 0), -400, 400), dy: clamp(num(it.dy, 0), -400, 400),
    };
  }
  if (kind === "curve") {
    const points = (Array.isArray(it.points) ? it.points : []).map(normalizePoint).filter(Boolean);
    if (points.length < 2) return null;
    const d = it.arrow && it.arrow !== "none" ? ANNOT_DEFAULTS.arrow : ANNOT_DEFAULTS.curve;
    return {
      ...base, kind, points, smooth: clamp(num(it.smooth, d.smooth), 0, 1), closed: !!it.closed && points.length > 2,
      width: clamp(num(it.width, d.width), 0.5, 24), color: base.color || d.color,
      dash: CURVE_DASHES.includes(it.dash) ? it.dash : "solid", arrow: CURVE_ARROWS.includes(it.arrow) ? it.arrow : "none",
      ...(num(it.reach, 0) > 0 ? { reach: num(it.reach, 0) } : {}),
    };
  }
  if (kind === "box") {
    const center = normalizePoint(it.center ?? it.at);
    const r = Array.isArray(it.radii) ? it.radii : Array.isArray(it.size) ? it.size.map((x) => num(x, NaN) / 2) : [it.radius, it.radius, it.radius];
    const radii = r.slice(0, 3).map((x) => Math.max(1e-3, Math.abs(num(x, NaN))));
    if (!center || radii.length !== 3 || radii.some((x) => !Number.isFinite(x))) return null;
    const style = BOX_STYLES.includes(it.style) ? it.style : "box";
    return { ...base, kind, center, radii, rot: vec3(it.rot) || [0, 0, 0], style, color: base.color || ANNOT_DEFAULTS.box.color };
  }
  if (kind === "sphere") {
    const center = normalizePoint(it.center ?? it.at);
    const r = Array.isArray(it.radii) ? it.radii : [it.radius, it.radius, it.radius];
    const radii = r.slice(0, 3).map((x) => Math.max(1e-3, Math.abs(num(x, NaN))));
    if (!center || radii.length !== 3 || radii.some((x) => !Number.isFinite(x))) return null;
    const style = SPHERE_STYLES.includes(it.style) ? it.style : "bubble";
    const rot = vec3(it.rot) || [0, 0, 0];
    return { ...base, kind, center, radii, rot, style, color: base.color || ANNOT_DEFAULTS[style === "wire" ? "shell" : style].color };
  }
  return null;
}

/** Notes from earlier figures (and classic manual labels) as text annotations. */
export function notesToAnnotations(notes) {
  const doc = emptyAnnotations();
  if (!Array.isArray(notes) || !notes.length) return doc;
  const g = { id: annotId("g"), name: "Notes", auto: false, visible: true };
  doc.groups.push(g);
  for (const n of notes) {
    if (!n || typeof n !== "object") continue;
    const at = n.anchor?.pos ? normalizePoint(n.anchor.pos) : normalizePoint(n.anchor);
    if (!at) continue;
    const item = normalizeItem({ kind: "text", at, text: n.text, color: n.color, size: Number.isFinite(n.size) && n.size > 0 ? n.size : 13, weight: 600, visible: n.visible });
    if (!item) continue;
    // A note on an object sits just above it, as notes did.
    if (at.trace != null) item.dy = -14;
    item.id = String(n.id || annotId("i"));
    item.group = g.id;
    doc.items.push(item);
  }
  return doc.items.length ? doc : emptyAnnotations();
}

// ------------------------------------------------------------------ points

/** World position of a point now: `objectPosition(trace, index)` for attached ones. */
export function resolveAnnotPoint(p, objectPosition) {
  if (!p) return null;
  if (p.pos) return p.pos;
  const at = objectPosition?.(p.trace, p.index);
  if (!at) return null;
  const o = p.offset;
  return o ? [at[0] + o[0], at[1] + o[1], at[2] + o[2]] : at;
}

/** Whether anything in the document moves with time (attached or present-day only). */
export function annotationsTimeDependent(doc) {
  const moving = (p) => p && p.trace != null;
  for (const it of doc.items) {
    if (it.present) return true;
    if (it.kind === "text" && moving(it.at)) return true;
    if ((it.kind === "sphere" || it.kind === "box") && moving(it.center)) return true;
    if (it.kind === "curve" && it.points.some(moving)) return true;
  }
  return false;
}

// ------------------------------------------------------------------ curves

/**
 * Sample a curve through `pts` (world positions). `smooth` blends straight
 * segments (0) into a centripetal Catmull–Rom spline (1), which never
 * overshoots or loops between close points. Returns positions and the
 * cumulative length at each sample.
 */
export function curveSamples(pts, { smooth = 1, closed = false, perSegment = 24 } = {}) {
  const n = pts.length;
  const out = [];
  if (n < 2) return { points: n ? [pts[0].slice()] : [], arc: n ? [0] : [] };
  const segs = closed ? n : n - 1;
  const at = (i) => pts[closed ? ((i % n) + n) % n : clamp(i, 0, n - 1)];
  const steps = smooth > 0 ? perSegment : 1;
  for (let s = 0; s < segs; s++) {
    const p1 = at(s), p2 = at(s + 1);
    // Ends of an open curve mirror their neighbour, so the tangent is natural.
    const p0 = !closed && s === 0 ? mirror(p2, p1) : at(s - 1);
    const p3 = !closed && s === segs - 1 ? mirror(p1, p2) : at(s + 2);
    for (let k = 0; k < steps; k++) {
      const t = k / steps;
      const c = smooth > 0 ? catmullRom(p0, p1, p2, p3, t) : p1;
      out.push(smooth >= 1 ? c : [lerp(p1[0] + (p2[0] - p1[0]) * t, c[0], smooth), lerp(p1[1] + (p2[1] - p1[1]) * t, c[1], smooth), lerp(p1[2] + (p2[2] - p1[2]) * t, c[2], smooth)]);
    }
  }
  out.push((closed ? at(0) : at(n - 1)).slice());
  const arc = [0];
  for (let i = 1; i < out.length; i++) arc.push(arc[i - 1] + Math.hypot(out[i][0] - out[i - 1][0], out[i][1] - out[i - 1][1], out[i][2] - out[i - 1][2]));
  return { points: out, arc };
}

function mirror(a, b) {
  // The reflection of a through b.
  return [2 * b[0] - a[0], 2 * b[1] - a[1], 2 * b[2] - a[2]];
}

function catmullRom(p0, p1, p2, p3, t) {
  // Centripetal parameterisation (alpha = 0.5), Barry–Goldman form.
  const tj = (ti, a, b) => ti + Math.max(Math.sqrt(Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2])), 1e-6);
  const t0 = 0, t1 = tj(t0, p0, p1), t2 = tj(t1, p1, p2), t3 = tj(t2, p2, p3);
  const u = lerp(t1, t2, t);
  const L = (a, b, ta, tb) => {
    const w = (u - ta) / (tb - ta);
    return [a[0] + (b[0] - a[0]) * w, a[1] + (b[1] - a[1]) * w, a[2] + (b[2] - a[2]) * w];
  };
  const a1 = L(p0, p1, t0, t1), a2 = L(p1, p2, t1, t2), a3 = L(p2, p3, t2, t3);
  const b1 = L(a1, a2, t0, t2), b2 = L(a2, a3, t1, t3);
  return L(b1, b2, t1, t2);
}

// ------------------------------------------------------------------ spheres

/** Rotation matrix (row-major 3×3) of Euler angles in degrees, applied x, then y, then z. */
export function sphereRotation(rot = [0, 0, 0]) {
  const [ax, ay, az] = rot.map((d) => (Number(d) || 0) * DEG);
  const cx = Math.cos(ax), sx = Math.sin(ax), cy = Math.cos(ay), sy = Math.sin(ay), cz = Math.cos(az), sz = Math.sin(az);
  // Rz · Ry · Rx
  return [
    cz * cy, cz * sy * sx - sz * cx, cz * sy * cx + sz * sx,
    sz * cy, sz * sy * sx + cz * cx, sz * sy * cx - cz * sx,
    -sy, cy * sx, cy * cx,
  ];
}

/** The three semi-axes of an ellipsoid in world space ([x, y, z] scaled by the radii). */
export function sphereAxes(item) {
  const m = sphereRotation(item.rot);
  const r = item.radii;
  return [0, 1, 2].map((k) => [m[k] * r[k], m[3 + k] * r[k], m[6 + k] * r[k]]);
}

/** A unit quaternion [x, y, z, w] for a rotation matrix (for the GPU). */
export function rotationQuat(m) {
  const tr = m[0] + m[4] + m[8];
  let x, y, z, w;
  if (tr > 0) {
    const s = Math.sqrt(tr + 1) * 2;
    w = 0.25 * s; x = (m[7] - m[5]) / s; y = (m[2] - m[6]) / s; z = (m[3] - m[1]) / s;
  } else if (m[0] > m[4] && m[0] > m[8]) {
    const s = Math.sqrt(1 + m[0] - m[4] - m[8]) * 2;
    w = (m[7] - m[5]) / s; x = 0.25 * s; y = (m[1] + m[3]) / s; z = (m[2] + m[6]) / s;
  } else if (m[4] > m[8]) {
    const s = Math.sqrt(1 + m[4] - m[0] - m[8]) * 2;
    w = (m[2] - m[6]) / s; x = (m[1] + m[3]) / s; y = 0.25 * s; z = (m[5] + m[7]) / s;
  } else {
    const s = Math.sqrt(1 + m[8] - m[0] - m[4]) * 2;
    w = (m[3] - m[1]) / s; x = (m[2] + m[6]) / s; y = (m[5] + m[7]) / s; z = 0.25 * s;
  }
  const l = Math.hypot(x, y, z, w) || 1;
  return [x / l, y / l, z / l, w / l];
}

/** Euler angles (degrees, as sphereRotation takes them) of a rotation matrix. */
export function rotationToEuler(m) {
  const sy = clamp(-m[6], -1, 1);
  const y = Math.asin(sy);
  let x, z;
  if (Math.abs(sy) < 0.99999) {
    x = Math.atan2(m[7], m[8]);
    z = Math.atan2(m[3], m[0]);
  } else {
    // Gimbal lock: put the whole turn about x.
    x = Math.atan2(-m[5], m[4]);
    z = 0;
  }
  return [x / DEG, y / DEG, z / DEG].map((d) => Math.round(d * 1e6) / 1e6);
}

/** Matrix of a turn by `angle` (radians) about the unit axis `u`. */
export function axisRotation(u, angle) {
  const c = Math.cos(angle), s = Math.sin(angle), t = 1 - c;
  const [x, y, z] = u;
  return [
    t * x * x + c, t * x * y - s * z, t * x * z + s * y,
    t * x * y + s * z, t * y * y + c, t * y * z - s * x,
    t * x * z - s * y, t * y * z + s * x, t * z * z + c,
  ];
}

export function mul3(a, b) {
  const o = new Array(9);
  for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) o[r * 3 + c] = a[r * 3] * b[c] + a[r * 3 + 1] * b[3 + c] + a[r * 3 + 2] * b[6 + c];
  return o;
}

/** Euler angles after turning `rot` by the matrix `turn` (applied in world space). */
export function turnEuler(rot, turn) {
  return rotationToEuler(mul3(turn, sphereRotation(rot)));
}

/** The eight corners of a box (bit 0: +x, bit 1: +y, bit 2: +z). */
export function boxCorners(item, center) {
  const axes = sphereAxes(item);
  const out = [];
  for (let k = 0; k < 8; k++) {
    const sx = k & 1 ? 1 : -1, sy = k & 2 ? 1 : -1, sz = k & 4 ? 1 : -1;
    out.push([0, 1, 2].map((j) => center[j] + sx * axes[0][j] + sy * axes[1][j] + sz * axes[2][j]));
  }
  return out;
}

/** Whether world point `p` lies inside the box (grown by `margin`). */
export function insideBox(item, center, p, margin = 1) {
  const m = sphereRotation(item.rot);
  const d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
  for (let k = 0; k < 3; k++) {
    const local = m[k] * d[0] + m[3 + k] * d[1] + m[6 + k] * d[2];
    if (Math.abs(local) > item.radii[k] * margin) return false;
  }
  return true;
}

/** Distance in pc from `p` to a polyline of world points. */
export function polylineDistance3(p, pts) {
  let best = Infinity;
  for (let i = 1; i < pts.length; i++) {
    const a = pts[i - 1], b = pts[i];
    const ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    const l2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    const u = l2 > 0 ? clamp(((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1] + (p[2] - a[2]) * ab[2]) / l2, 0, 1) : 0;
    best = Math.min(best, Math.hypot(p[0] - a[0] - ab[0] * u, p[1] - a[1] - ab[1] * u, p[2] - a[2] - ab[2] * u));
  }
  return best;
}

/**
 * A test for "is this world point selected by `item`": inside a sphere or
 * box, or within its reach of a curve. `center` / `samples` are the item's
 * resolved centre or sampled curve. Returns null for items that select nothing.
 */
export function selectionTest(item, { center, samples } = {}) {
  if (!item.select) return null;
  if (item.kind === "sphere" && center) return (p) => insideSphere(item, center, p);
  if (item.kind === "box" && center) return (p) => insideBox(item, center, p);
  if (item.kind === "curve" && samples?.length > 1) {
    const r = curveReach(item, samples);
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (const q of samples) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], q[k] - r); hi[k] = Math.max(hi[k], q[k] + r); }
    return (p) => p[0] >= lo[0] && p[0] <= hi[0] && p[1] >= lo[1] && p[1] <= hi[1] && p[2] >= lo[2] && p[2] <= hi[2] && polylineDistance3(p, samples) <= r;
  }
  return null;
}

/** A curve's selection reach in pc: its own, else a twentieth of its length. */
export function curveReach(item, samples) {
  if (item.reach > 0) return item.reach;
  let len = 0;
  for (let i = 1; i < (samples?.length || 0); i++) len += Math.hypot(samples[i][0] - samples[i - 1][0], samples[i][1] - samples[i - 1][1], samples[i][2] - samples[i - 1][2]);
  return Math.max(1, Math.round(len / 20));
}

/** Whether world point `p` lies within `margin` × the ellipsoid. */
export function insideSphere(item, center, p, margin = 1) {
  const m = sphereRotation(item.rot);
  const d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
  let s = 0;
  for (let k = 0; k < 3; k++) {
    const local = m[k] * d[0] + m[3 + k] * d[1] + m[6 + k] * d[2];
    s += (local / (item.radii[k] * margin)) ** 2;
  }
  return s <= 1;
}

// ------------------------------------------------------------------ screen helpers

/** Distance from (x, y) to the polyline `pts` ([[x, y], …]) in pixels. */
export function polylineDistance(x, y, pts) {
  let best = Infinity;
  for (let i = 1; i < pts.length; i++) {
    const a = pts[i - 1], b = pts[i];
    if (!a || !b) continue;
    const dx = b[0] - a[0], dy = b[1] - a[1];
    const l2 = dx * dx + dy * dy;
    const u = l2 > 0 ? clamp(((x - a[0]) * dx + (y - a[1]) * dy) / l2, 0, 1) : 0;
    best = Math.min(best, Math.hypot(x - a[0] - dx * u, y - a[1] - dy * u));
  }
  return best;
}

// ------------------------------------------------------------------ names and groups

const KIND_NAMES = {
  text: ["Label", "Labels"], curve: ["Curve", "Curves"], arrow: ["Arrow", "Arrows"],
  bubble: ["Bubble", "Bubbles"], shell: ["Shell", "Shells"], wire: ["Sphere", "Spheres"],
  box: ["Box", "Boxes"],
};

/** What an item looks like to a reader: label, curve, arrow, bubble, shell or wire sphere. */
export function annotKind(item) {
  if (item.kind === "curve") return item.arrow && item.arrow !== "none" ? "arrow" : "curve";
  if (item.kind === "sphere") return item.style;
  if (item.kind === "box") return "box";
  return "text";
}

export function annotKindName(item, plural = false) {
  const k = KIND_NAMES[annotKind(item)] || KIND_NAMES.text;
  return k[plural ? 1 : 0];
}

/** The name a group shows: what its author typed, the label it is named after, or its kinds. */
export function annotGroupName(doc, group) {
  if (!group) return "";
  if (!group.auto && group.name) return group.name;
  const items = doc.items.filter((i) => i.group === group.id);
  if (group.nameFrom) {
    const t = items.find((i) => i.id === group.nameFrom);
    const line = t?.text?.split("\n")[0].trim();
    if (line) return line.length > 40 ? `${line.slice(0, 39)}…` : line;
  }
  if (!items.length) return group.name || "Annotations";
  // A lone label is its own name.
  if (items.length === 1 && items[0].kind === "text") {
    const line = items[0].text?.split("\n")[0].trim();
    if (line) return line.length > 40 ? `${line.slice(0, 39)}…` : line;
  }
  const kinds = new Set(items.map(annotKind));
  if (kinds.size === 1) return annotKindName(items[0], items.length > 1);
  return "Annotations";
}

/** The kind and colour a group's key entry shows: its first shape, else its first label. */
export function annotGroupLook(doc, group) {
  const items = doc.items.filter((i) => i.group === group.id);
  const lead = items.find((i) => i.kind !== "text") || items[0];
  return lead ? { kind: annotKind(lead), color: lead.color } : { kind: "text", color: "#ffffff" };
}

/**
 * Choose a group for a new item, the way people sort what they draw:
 *   1. a label put on or right beside a shape joins that shape (and names it);
 *   2. a shape drawn around or onto a label joins the label's group;
 *   3. otherwise it joins the group of the same kind and colour drawn
 *      before (a run of red arrows is one key entry);
 *   4. otherwise it starts a group of its own.
 * `ctx.resolve(point)` gives world positions; `ctx.screenDistance(item, pos)`,
 * when given, the pixel distance from a world point to a curve on screen.
 * Mutates `doc` (adds a group when needed) and returns the group id.
 */
export function assignGroup(doc, item, ctx = {}) {
  const resolve = ctx.resolve || ((p) => p?.pos || null);
  const near = ctx.nearPx ?? 36;
  const others = doc.items.filter((i) => i.id !== item.id && i.visible !== false && doc.groups.some((g) => g.id === i.group));
  const groupOf = (id) => doc.groups.find((g) => g.id === id);
  const nameAfter = (group, textItem) => {
    if (group.auto && !group.nameFrom && textItem?.text?.trim()) group.nameFrom = textItem.id;
  };
  if (item.kind === "text") {
    const at = resolve(item.at);
    let best = null;
    for (const o of others) {
      if (o.kind === "sphere" || o.kind === "box") {
        const c = resolve(o.center);
        if (!at || !c || !(o.kind === "box" ? insideBox(o, c, at, 1.25) : insideSphere(o, c, at, 1.25))) continue;
        const d = Math.hypot(at[0] - c[0], at[1] - c[1], at[2] - c[2]) / Math.max(...o.radii);
        if (!best || d < best.d) best = { o, d };
      } else if (o.kind === "curve" && at && ctx.screenDistance) {
        const px = ctx.screenDistance(o, at);
        if (px <= near && (!best || px / near < best.d)) best = { o, d: px / near };
      }
    }
    if (best) {
      const g = groupOf(best.o.group);
      nameAfter(g, item);
      return g.id;
    }
  } else {
    let best = null;
    for (const o of others) {
      if (o.kind !== "text") continue;
      const at = resolve(o.at);
      if (!at) continue;
      if (item.kind === "sphere" || item.kind === "box") {
        const c = resolve(item.center);
        if (!c || !(item.kind === "box" ? insideBox(item, c, at, 1.25) : insideSphere(item, c, at, 1.25))) continue;
        const d = Math.hypot(at[0] - c[0], at[1] - c[1], at[2] - c[2]) / Math.max(...item.radii);
        if (!best || d < best.d) best = { o, d };
      } else if (item.kind === "curve" && ctx.screenDistance) {
        const px = ctx.screenDistance(item, at);
        if (px <= near && (!best || px / near < best.d)) best = { o, d: px / near };
      }
    }
    // A label already grouped with another shape keeps that group.
    if (best) {
      const g = groupOf(best.o.group);
      const shapes = doc.items.filter((i) => i.group === g.id && i.kind !== "text" && i.id !== item.id);
      if (!shapes.length || shapes.every((s) => annotKind(s) === annotKind(item))) {
        nameAfter(g, best.o);
        return g.id;
      }
    }
  }
  // Same kind and colour as an earlier entry (newest first).
  const kind = annotKind(item);
  const color = String(item.color || "").toLowerCase();
  for (let i = others.length - 1; i >= 0; i--) {
    const o = others[i];
    const g = groupOf(o.group);
    if (!g || !g.auto || g.nameFrom || annotKind(o) !== kind || String(o.color || "").toLowerCase() !== color) continue;
    const members = doc.items.filter((x) => x.group === g.id && x.id !== item.id);
    if (members.every((x) => annotKind(x) === kind && String(x.color || "").toLowerCase() === color)) return g.id;
  }
  const g = { id: annotId("g"), name: "", auto: true, visible: true };
  doc.groups.push(g);
  return g.id;
}

// ------------------------------------------------------------------ transitions

const mixColor = (a, b, t) => {
  const ca = parseColor(a), cb = parseColor(b);
  const c = [0, 1, 2].map((k) => Math.round(clamp(lerp(ca[k], cb[k], t), 0, 1) * 255));
  return `#${c.map((x) => x.toString(16).padStart(2, "0")).join("")}`;
};

/**
 * The document part-way (t, 0–1) from `a` to `b`, for State transitions:
 * items in both (same id and kind) move, resize and recolour; items in one
 * only fade out or in. Each returned item carries `fade` (0–1) to multiply
 * its opacity by. Points are resolved to positions (attached ones too), so
 * the result is for drawing, not for editing.
 */
export function lerpAnnotations(a, b, t, resolve = (p) => p?.pos || null) {
  const out = emptyAnnotations();
  const gVis = (doc, id) => doc.groups.find((g) => g.id === id)?.visible !== false;
  const inA = new Map(a.items.map((i) => [i.id, i]));
  const inB = new Map(b.items.map((i) => [i.id, i]));
  const groups = new Map();
  for (const g of [...a.groups, ...b.groups]) groups.set(g.id, { ...g, visible: true });
  out.groups = [...groups.values()];
  const shown = (doc, it) => (it.visible !== false && gVis(doc, it.group) ? 1 : 0);
  const fixed = (p) => { const pos = resolve(p); return pos ? { pos: pos.slice() } : null; };
  const mixPoint = (p, q, u) => {
    const x = resolve(p), y = resolve(q);
    if (!x || !y) return fixed(u < 0.5 ? p : q);
    return { pos: [lerp(x[0], y[0], u), lerp(x[1], y[1], u), lerp(x[2], y[2], u)] };
  };
  for (const id of new Set([...inA.keys(), ...inB.keys()])) {
    const x = inA.get(id), y = inB.get(id);
    if (x && y && x.kind === y.kind) {
      const fa = shown(a, x), fb = shown(b, y);
      const fade = lerp(fa, fb, t);
      if (fade <= 1e-4) continue;
      const base = { ...(t < 0.5 ? x : y), fade, opacity: lerp(x.opacity, y.opacity, t), color: mixColor(x.color, y.color, t) };
      if (x.kind === "text") {
        base.at = mixPoint(x.at, y.at, t);
        base.size = lerp(x.size, y.size, t);
        base.dx = lerp(x.dx || 0, y.dx || 0, t);
        base.dy = lerp(x.dy || 0, y.dy || 0, t);
      } else if (x.kind === "sphere" || x.kind === "box") {
        base.center = mixPoint(x.center, y.center, t);
        base.radii = x.radii.map((r, k) => Math.exp(lerp(Math.log(r), Math.log(y.radii[k]), t)));
        base.rot = x.rot.map((r, k) => lerp(r, y.rot[k], t));
      } else if (x.points.length === y.points.length) {
        base.points = x.points.map((p, k) => mixPoint(p, y.points[k], t));
        base.width = lerp(x.width, y.width, t);
        base.smooth = lerp(x.smooth, y.smooth, t);
      } else {
        base.points = (t < 0.5 ? x : y).points.map(fixed);
      }
      if (base.at === null || base.center === null || (base.points && base.points.some((p) => !p))) continue;
      out.items.push(base);
    } else {
      const it = x || y;
      const fade = x ? shown(a, x) * (1 - t) : shown(b, y) * t;
      if (fade <= 1e-4) continue;
      const copy = { ...it, fade };
      if (it.kind === "text") copy.at = fixed(it.at);
      else if (it.kind === "sphere" || it.kind === "box") copy.center = fixed(it.center);
      else copy.points = it.points.map(fixed);
      if (copy.at === null || copy.center === null || (copy.points && copy.points.some((p) => !p))) continue;
      out.items.push(copy);
    }
  }
  return out;
}
