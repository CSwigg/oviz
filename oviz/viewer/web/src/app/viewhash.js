// Shareable view links:
//   #t=<time>&v=<mode>&c=<target,distance,yaw,pitch,fov>&a=<anchor>&g=<group>
//   &x=<state>                 the whole viewer State (see encodeStatePart): every
//                              property of the view, so the link opens exactly as sent
//   &w=<views>                 the saved views (see encodeViewsPart), so the receiver
//                              gets the same story to explore or present
//   &p=<n>                     open presenting the views, from view n (1-based)
// Links from before x= existed (or written without one) may instead keep:
//   l=<1|0 per layer>          the layers shown (traces and volumes, in the layers panel's order)
//   s=<h|d>[o],<trace>:<set>…  a lasso selection: hide (h) or dim (d) the rest, o = its filter off;
//                              per trace (its index in the figure) the selected objects
//   k=<outline>                the lasso's outline and camera, which clip the volumes
//   o=<trace>:<index>          the selected (clicked) object
// Object sets are base64url, a bitmap ("b…") or gaps between sorted indices as
// varints ("i…"), whichever is shorter.

import { encodeAnchor, decodeAnchor } from "../engine/anchor.js";

/**
 * Compact, human-readable hash: #t=-12.5&v=3d&c=x,y,z,d,yaw,pitch,fov&a=lsr.
 * `extra` adds what the link should also keep: `layers` (booleans),
 * `lasso` ({selection: [[traceIndex, indices]], isolate, filterOff, mask})
 * and `object` ({trace: traceIndex, index}).
 */
export function encodeViewHash(viewer, extra = {}) {
  const p = viewer.pose;
  const r = (x, d = 2) => Number(x.toFixed(d));
  const parts = [
    `t=${r(viewer.timeline.time, 3)}`,
    `v=${viewer.state.view.mode}`,
    `c=${[...p.target.map((x) => r(x, 1)), r(p.distance, 1), r(p.yaw, 4), r(p.pitch, 4), r(p.fov, 2)].join(",")}`,
  ];
  const anchor = encodeAnchor(viewer.state.view.anchor);
  if (anchor) parts.push(`a=${encodeURIComponent(anchor)}`);
  if (viewer.state.group) parts.push(`g=${encodeURIComponent(viewer.state.group)}`);
  if (Array.isArray(extra.layers) && extra.layers.length) parts.push(`l=${extra.layers.map((x) => (x ? "1" : "0")).join("")}`);
  const lasso = extra.lasso;
  if (lasso && (lasso.selection || lasso.mask)) {
    const sets = [];
    for (const [trace, idx] of lasso.selection || []) {
      const set = encodeIndexSet(idx);
      if (Number.isInteger(trace) && trace >= 0 && set) sets.push(`${trace}:${set}`);
    }
    parts.push(`s=${[(lasso.isolate === false ? "d" : "h") + (lasso.filterOff ? "o" : ""), ...sets].join(",")}`);
    const outline = lasso.mask ? encodeMask(lasso.mask) : "";
    if (outline) parts.push(`k=${outline}`);
  }
  const o = extra.object;
  if (o && Number.isInteger(o.trace) && Number.isInteger(o.index) && o.trace >= 0 && o.index >= 0) parts.push(`o=${o.trace}:${o.index}`);
  if (typeof extra.state === "string" && /^1[zj][A-Za-z0-9_-]*$/.test(extra.state)) parts.push(`x=${extra.state}`);
  if (typeof extra.views === "string" && /^1[zj][A-Za-z0-9_-]+$/.test(extra.views)) {
    parts.push(`w=${extra.views}`);
    const start = Number(extra.present);
    if (Number.isInteger(start) && start >= 1) parts.push(`p=${start}`);
  }
  return parts.join("&");
}

export function decodeViewHash(hash) {
  const out = {};
  for (const part of String(hash || "").replace(/^#/, "").split("&")) {
    const [k, val] = part.split("=");
    if (!k || val == null) continue;
    // A malformed escape in a pasted link must not stop the figure booting.
    try { out[k] = decodeURIComponent(val); } catch (_) { /* ignore this part */ }
  }
  if (!out.c && !out.t) return null;
  const view = {};
  if (out.t != null && Number.isFinite(Number(out.t))) view.time = Number(out.t);
  if (out.v) view.mode = out.v === "sky" ? "sky" : "3d";
  if (out.g) view.group = out.g;
  const anchor = out.a ? decodeAnchor(out.a) : null;
  if (anchor) view.anchor = anchor;
  if (out.c) {
    const n = out.c.split(",").map(Number);
    if (n.length === 7 && n.every(Number.isFinite)) {
      view.pose = { target: n.slice(0, 3), distance: n[3], yaw: n[4], pitch: n[5], fov: n[6] };
    }
  }
  if (out.l && /^[01]+$/.test(out.l)) view.layers = [...out.l].map((c) => c === "1");
  if (out.s != null) {
    const [flags, ...sets] = out.s.split(",");
    if (/^[hd]o?$/.test(flags)) {
      view.isolate = flags[0] !== "d";
      view.filterOff = flags.includes("o");
      view.selection = [];
      for (const s of sets) {
        const m = /^(\d+):(.+)$/.exec(s);
        const idx = m ? decodeIndexSet(m[2]) : null;
        if (idx) view.selection.push([Number(m[1]), idx]);
      }
    }
  }
  if (out.k) {
    const mask = decodeMask(out.k);
    if (mask) view.mask = mask;
  }
  const o = out.o ? /^(\d+):(\d+)$/.exec(out.o) : null;
  if (o) view.object = { trace: Number(o[1]), index: Number(o[2]) };
  // The whole State: decoded (asynchronously) against the figure's opening State.
  if (out.x && /^1[zj][A-Za-z0-9_-]+$/.test(out.x)) view.state = out.x;
  // The saved views, and whether to open presenting them.
  if (out.w && /^1[zj][A-Za-z0-9_-]+$/.test(out.w)) {
    view.views = out.w;
    const start = /^\d+$/.test(out.p || "") ? Number(out.p) : 0;
    if (start >= 1) view.present = start;
  }
  return view;
}

// ------------------------------------------------------------------ whole State

/**
 * A view link's x= part: the whole viewer State (`captureState`) as what it
 * changes in `base`, the State the figure opens with. The receiving figure
 * opens with the same `base` and rebuilds the State exactly, so every
 * property (camera, time, layer styles, volumes, the Sky stack and lens,
 * filters, the lasso, notes, widgets, panels) comes back as it was sent.
 * Numbers keep 7 significant digits, lasso object sets use the compact form
 * above and the JSON is deflated: "1z<base64url>", or "1j<base64url>" where
 * the browser cannot compress. Resolves to "" when nothing differs.
 */
export async function encodeStatePart(state, base) {
  const delta = diffJson(base, state);
  if (delta === undefined) return "";
  const bytes = new TextEncoder().encode(JSON.stringify(packValue(delta, "")));
  const z = await deflateBytes(bytes);
  return z ? `1z${toBase64Url(z)}` : `1j${toBase64Url(bytes)}`;
}

/**
 * A link's w= part: the saved views (name, caption, camera behaviour,
 * transition and State; thumbnails are left out and redrawn as each view is
 * visited), each State as what it changes in `base` like x=, deflated
 * together with the story's default transition. "" when there are none.
 */
export async function encodeViewsPart(items, base, defaultTransition = null) {
  if (!Array.isArray(items) || !items.length || !base) return "";
  const views = items.map((it) => {
    const o = { n: String(it.name || "") };
    if (it.caption) o.c = String(it.caption);
    if (it.camera === "keep") o.k = 1;
    if (isPlainObject(it.transition)) o.t = it.transition;
    const d = diffJson(base, it.state);
    o.s = d === undefined ? {} : packValue(d, "");
    return o;
  });
  const payload = { v: views };
  if (isPlainObject(defaultTransition)) payload.d = defaultTransition;
  const bytes = new TextEncoder().encode(JSON.stringify(payload));
  const z = await deflateBytes(bytes);
  return z ? `1z${toBase64Url(z)}` : `1j${toBase64Url(bytes)}`;
}

/**
 * The views a w= part describes: {items: [{name, caption, camera, transition,
 * state}], defaultTransition}, each State rebuilt over a copy of `base`;
 * null if it cannot be read.
 */
export async function decodeViewsPart(text, base) {
  const m = /^1([zj])([A-Za-z0-9_-]+)$/.exec(String(text || ""));
  if (!m || !base) return null;
  let bytes = fromBase64Url(m[2]);
  if (bytes && m[1] === "z") bytes = await inflateBytes(bytes);
  if (!bytes) return null;
  let payload;
  try { payload = JSON.parse(new TextDecoder().decode(bytes)); } catch (_) { return null; }
  if (!isPlainObject(payload) || !Array.isArray(payload.v)) return null;
  const items = [];
  for (const o of payload.v) {
    if (!isPlainObject(o) || !isPlainObject(o.s)) continue;
    const state = mergeJson(JSON.parse(JSON.stringify(base)), unpackValue(o.s));
    if (!isPlainObject(state) || !state.view) continue;
    items.push({
      name: String(o.n || `View ${items.length + 1}`),
      caption: String(o.c || ""),
      camera: o.k ? "keep" : "follow",
      transition: isPlainObject(o.t) ? o.t : null,
      state,
    });
  }
  return { items, defaultTransition: isPlainObject(payload.d) ? payload.d : null };
}

/** The State an x= part describes, rebuilt over (a copy of) `base`; null if it cannot be read. */
export async function decodeStatePart(text, base) {
  const m = /^1([zj])([A-Za-z0-9_-]+)$/.exec(String(text || ""));
  if (!m || !base) return null;
  let bytes = fromBase64Url(m[2]);
  if (bytes && m[1] === "z") bytes = await inflateBytes(bytes);
  if (!bytes) return null;
  let delta;
  try { delta = JSON.parse(new TextDecoder().decode(bytes)); } catch (_) { return null; }
  if (!isPlainObject(delta)) return null;
  return mergeJson(JSON.parse(JSON.stringify(base)), unpackValue(delta));
}

const isPlainObject = (x) => !!x && typeof x === "object" && !Array.isArray(x);

/**
 * What `b` changes in `a` (JSON values), or undefined when they are equal.
 * Objects recurse; equal-length arrays of objects change element by element
 * ({$arr: {index: change}}); other values are replaced whole; a key `b`
 * no longer has is null.
 */
function diffJson(a, b) {
  if (isPlainObject(a) && isPlainObject(b)) {
    const out = {};
    let any = false;
    for (const k of Object.keys(b)) {
      if (b[k] === undefined) continue;
      const d = diffJson(a[k], b[k]);
      if (d !== undefined) { out[k] = d; any = true; }
    }
    for (const k of Object.keys(a)) {
      if (a[k] !== undefined && (!(k in b) || b[k] === undefined)) { out[k] = null; any = true; }
    }
    return any ? out : undefined;
  }
  if (Array.isArray(a) && Array.isArray(b) && a.length === b.length && b.length && b.every(isPlainObject) && a.every(isPlainObject)) {
    const changes = {};
    let n = 0;
    b.forEach((x, i) => {
      const d = diffJson(a[i], x);
      if (d !== undefined) { changes[i] = d; n++; }
    });
    if (!n) return undefined;
    return n < b.length ? { $arr: changes } : b;
  }
  return JSON.stringify(a) === JSON.stringify(b) ? undefined : b === undefined ? null : b;
}

/** Apply a `diffJson` change to `a` (which it may modify). */
function mergeJson(a, d) {
  if (!isPlainObject(d)) return d;
  if (isPlainObject(d.$arr) && Object.keys(d).length === 1) {
    if (!Array.isArray(a)) return a;
    for (const [i, v] of Object.entries(d.$arr)) {
      const k = Number(i);
      if (Number.isInteger(k) && k >= 0 && k < a.length) a[k] = mergeJson(a[k], v);
    }
    return a;
  }
  const out = isPlainObject(a) ? a : {};
  for (const [k, v] of Object.entries(d)) {
    if (v === null) delete out[k];
    else out[k] = mergeJson(out[k], v);
  }
  return out;
}

/** Shrink a change for the link: 7 significant digits, lasso sets as compact strings, outlines simplified. */
function packValue(v, path) {
  if (typeof v === "number") return Number.isInteger(v) || !Number.isFinite(v) ? v : +v.toPrecision(7);
  if (Array.isArray(v)) {
    if (/\.lasso\.selection\.[^.]+$/.test(path) && v.every((i) => Number.isInteger(i) && i >= 0)) return { $is: encodeIndexSet(v) };
    if (/\.lasso\.mask\.polygon$/.test(path) && v.length > 4) return simplifyPolygon(v, 0.003).map(([x, y]) => [+x.toFixed(4), +y.toFixed(4)]);
    return v.map((x, i) => packValue(x, `${path}.${i}`));
  }
  if (isPlainObject(v)) {
    const out = {};
    for (const [k, x] of Object.entries(v)) out[k] = packValue(x, `${path}.${k}`);
    return out;
  }
  return v;
}

function unpackValue(v) {
  if (Array.isArray(v)) return v.map(unpackValue);
  if (isPlainObject(v)) {
    if (typeof v.$is === "string" && Object.keys(v).length === 1) return decodeIndexSet(v.$is) || [];
    const out = {};
    for (const [k, x] of Object.entries(v)) out[k] = unpackValue(x);
    return out;
  }
  return v;
}

async function deflateBytes(bytes) {
  if (typeof CompressionStream !== "function") return null;
  try {
    const stream = new Blob([bytes]).stream().pipeThrough(new CompressionStream("deflate-raw"));
    return new Uint8Array(await new Response(stream).arrayBuffer());
  } catch (_) {
    return null;
  }
}

async function inflateBytes(bytes) {
  if (typeof DecompressionStream !== "function") return null;
  try {
    const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("deflate-raw"));
    return new Uint8Array(await new Response(stream).arrayBuffer());
  } catch (_) {
    return null;
  }
}

// ------------------------------------------------------------------ object sets

/** Sorted unique indices as "b<bitmap>" or "i<gap varints>", base64url, whichever is shorter. */
function encodeIndexSet(indices) {
  const idx = [...new Set(Array.from(indices || [], Number))].filter((i) => Number.isInteger(i) && i >= 0).sort((a, b) => a - b);
  if (!idx.length) return "";
  const bits = new Uint8Array((idx[idx.length - 1] >> 3) + 1);
  for (const i of idx) bits[i >> 3] |= 1 << (i & 7);
  const gaps = [];
  let prev = -1;
  for (const i of idx) {
    let d = i - prev - 1;
    prev = i;
    while (d >= 0x80) { gaps.push((d & 0x7f) | 0x80); d = Math.floor(d / 0x80); }
    gaps.push(d);
  }
  return gaps.length < bits.length ? `i${toBase64Url(Uint8Array.from(gaps))}` : `b${toBase64Url(bits)}`;
}

function decodeIndexSet(text) {
  const bytes = fromBase64Url(text.slice(1));
  if (!bytes) return null;
  const out = [];
  if (text[0] === "b") {
    for (let i = 0; i < bytes.length * 8; i++) if (bytes[i >> 3] & (1 << (i & 7))) out.push(i);
  } else if (text[0] === "i") {
    let prev = -1, d = 0, shift = 1;
    for (const b of bytes) {
      d += (b & 0x7f) * shift;
      if (b & 0x80) { shift *= 0x80; continue; }
      prev += d + 1;
      out.push(prev);
      d = 0;
      shift = 1;
    }
  } else return null;
  return out;
}

// ------------------------------------------------------------------ lasso outline

const NDC_SPAN = 4; // outlines are stored for x, y in [−2, 2] (strokes can leave the canvas)

/**
 * The lasso outline in normalised device coordinates, simplified to about a
 * pixel and stored as 16-bit pairs, then the camera it was drawn with as 16
 * float32 (the same view-projection the volumes are clipped with).
 */
function encodeMask(mask) {
  const vp = Array.from(mask?.vp || [], Number);
  if (!Array.isArray(mask?.polygon) || vp.length !== 16 || !vp.every(Number.isFinite)) return "";
  const poly = simplifyPolygon(mask.polygon.filter((q) => Number.isFinite(q?.[0]) && Number.isFinite(q?.[1])), 0.003);
  if (poly.length < 3 || poly.length > 0xffff) return "";
  const dv = new DataView(new ArrayBuffer(2 + poly.length * 4 + 64));
  const q = (x) => Math.round(((Math.min(Math.max(x, -NDC_SPAN / 2), NDC_SPAN / 2) + NDC_SPAN / 2) / NDC_SPAN) * 0xffff);
  dv.setUint16(0, poly.length, true);
  poly.forEach(([x, y], i) => { dv.setUint16(2 + i * 4, q(x), true); dv.setUint16(4 + i * 4, q(y), true); });
  vp.forEach((x, i) => dv.setFloat32(2 + poly.length * 4 + i * 4, x, true));
  return toBase64Url(new Uint8Array(dv.buffer));
}

function decodeMask(text) {
  const bytes = fromBase64Url(text);
  if (!bytes || bytes.length < 2) return null;
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const n = dv.getUint16(0, true);
  if (n < 3 || bytes.length !== 2 + n * 4 + 64) return null;
  const dq = (u) => (u / 0xffff) * NDC_SPAN - NDC_SPAN / 2;
  const polygon = Array.from({ length: n }, (_, i) => [dq(dv.getUint16(2 + i * 4, true)), dq(dv.getUint16(4 + i * 4, true))]);
  const vp = Array.from({ length: 16 }, (_, i) => dv.getFloat32(2 + n * 4 + i * 4, true));
  return vp.every(Number.isFinite) ? { polygon, vp } : null;
}

/** Douglas–Peucker on a closed outline: the points that matter at `tol`. */
function simplifyPolygon(pts, tol) {
  if (pts.length <= 4) return pts.map((q) => [q[0], q[1]]);
  const keep = new Uint8Array(pts.length);
  // Split the loop at its first point and the point farthest from it.
  let far = 0, best = -1;
  for (let i = 1; i < pts.length; i++) {
    const d = Math.hypot(pts[i][0] - pts[0][0], pts[i][1] - pts[0][1]);
    if (d > best) { best = d; far = i; }
  }
  keep[0] = keep[far] = 1;
  const stack = [[0, far], [far, pts.length]];
  while (stack.length) {
    const [a, b] = stack.pop();
    const pa = pts[a], pb = pts[b % pts.length];
    let idx = -1, dmax = tol;
    for (let i = a + 1; i < b; i++) {
      const d = segmentDistance(pts[i], pa, pb);
      if (d > dmax) { dmax = d; idx = i; }
    }
    if (idx >= 0) { keep[idx] = 1; stack.push([a, idx], [idx, b]); }
  }
  return pts.filter((_, i) => keep[i]).map((q) => [q[0], q[1]]);
}

function segmentDistance(p, a, b) {
  const dx = b[0] - a[0], dy = b[1] - a[1];
  const len2 = dx * dx + dy * dy;
  const t = len2 > 0 ? Math.min(1, Math.max(0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2)) : 0;
  return Math.hypot(p[0] - (a[0] + t * dx), p[1] - (a[1] + t * dy));
}

// ------------------------------------------------------------------ base64url

function toBase64Url(bytes) {
  let s = "";
  for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  return btoa(s).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
}

function fromBase64Url(text) {
  if (!/^[A-Za-z0-9_-]*$/.test(text)) return null;
  try {
    const s = atob(text.replace(/-/g, "+").replace(/_/g, "/") + "=".repeat((4 - (text.length % 4)) % 4));
    return Uint8Array.from(s, (c) => c.charCodeAt(0));
  } catch (_) {
    return null;
  }
}
