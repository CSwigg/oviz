// AR model of the figure, for Apple's AR Quick Look (iPhone, iPad).
//
// A tabletop time-lapse of what is on screen:
// - every visible object as a small glowing sphere in its displayed colour,
//   moving along its track through the figure's timeline and growing in at
//   its birth, like the viewer's birth fade;
// - each visible dust volume as three stacks of see-through image slices
//   computed from the voxels with the renderer's own opacity model, which
//   fades in at the present day (dust maps are present-day);
// - a dark base plate with distance rings, the Galactic-centre direction and
//   a time readout that follows the animation.
// Quick Look plays the timeline; pinch to resize.
// Galactic x, y, z (pc) map to AR metres as (x, z, −y): +Y is Galactic north.

import { parseColor } from "../core/color.js";
import { clamp, niceFloor } from "../core/math.js";
import { birthFadeAt } from "../app/timeline.js";
import { usdzPackage } from "./usdz.js";

const TARGET_RADIUS_M = 0.4;   // most of the model fits within this radius
const FIT_QUANTILE = 0.9;      // …measured on this share of present-day positions
const VIEW_RADII = 1.45;       // objects farther than this (× fit radius) are out of view
const LIFT_M = 0.03;           // gap between the plate and the lowest object
const MIN_RADIUS_M = 0.0016;
const MAX_RADIUS_M = 0.01;
const SUN_RADIUS_M = 0.0045;
const TIME_CODES_PER_S = 12;   // animation speed: 12 timeline samples a second
const MAX_SAMPLES = 150;       // cap on time samples per track
const HOLD_START_S = 1;
const HOLD_END_S = 5;
const DUST_SLICES = 24;        // per axis (three orthogonal stacks)
const DUST_RES = 320;          // texels across the widest side of a slice
const DUST_STACK_WEIGHT = 0.6; // each stack's share of the optical depth
const PLATE_RES = 2048;

// ------------------------------------------------------------ geometry

/** Unit icosphere (level 0 = icosahedron), outward CCW faces. */
export function icosphere(level = 0) {
  const t = (1 + Math.sqrt(5)) / 2;
  const verts = [
    [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
    [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
    [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
  ].map(unit);
  let faces = [
    [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
    [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
    [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
    [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
  ];
  for (let l = 0; l < level; l++) {
    const mid = new Map();
    const midpoint = (a, b) => {
      const key = a < b ? `${a},${b}` : `${b},${a}`;
      if (!mid.has(key)) {
        mid.set(key, verts.length);
        verts.push(unit([(verts[a][0] + verts[b][0]) / 2, (verts[a][1] + verts[b][1]) / 2, (verts[a][2] + verts[b][2]) / 2]));
      }
      return mid.get(key);
    };
    const next = [];
    for (const [a, b, c] of faces) {
      const ab = midpoint(a, b), bc = midpoint(b, c), ca = midpoint(c, a);
      next.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]);
    }
    faces = next;
  }
  return { verts: Float32Array.from(verts.flat()), faces: Uint32Array.from(faces.flat()) };
}

function unit(v) {
  const l = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / l, v[1] / l, v[2] / l];
}

function sphereMesh(level) {
  const g = icosphere(level);
  return { points: g.verts, normals: g.verts, indices: g.faces };
}

/** A flat quad from four corners (counter-clockwise from the front), with st. */
function quadMesh(corners, normal, { doubleSided = true } = {}) {
  return {
    points: Float32Array.from(corners.flat()),
    normals: Float32Array.from([normal, normal, normal, normal].flat()),
    uvs: Float32Array.from([0, 0, 1, 0, 1, 1, 0, 1]),
    indices: Uint32Array.from([0, 1, 2, 0, 2, 3]),
    doubleSided,
  };
}

/** A flat disc (triangle fan) in the XZ plane at height y, st over its square. */
function discMesh(radius, y, segments = 128) {
  const p = [0, y, 0], n = [0, 1, 0], uv = [0.5, 0.5], idx = [];
  for (let i = 0; i <= segments; i++) {
    const a = (i / segments) * Math.PI * 2;
    const x = Math.cos(a) * radius, z = Math.sin(a) * radius;
    p.push(x, y, z);
    n.push(0, 1, 0);
    uv.push((x / radius + 1) / 2, 1 - (z / radius + 1) / 2);
    if (i > 0) idx.push(0, i + 1, i); // counter-clockwise seen from above (+Y)
  }
  return { points: Float32Array.from(p), normals: Float32Array.from(n), uvs: Float32Array.from(uv), indices: Uint32Array.from(idx), doubleSided: true };
}

// ------------------------------------------------------------ helpers

function lutColor(lut, u) {
  const w = lut.length / 4;
  const x = clamp(u, 0, 1) * (w - 1);
  const i = Math.floor(x), j = Math.min(i + 1, w - 1), t = x - i;
  const ch = (k) => (lut[i * 4 + k] * (1 - t) + lut[j * 4 + k] * t) / 255;
  return [ch(0), ch(1), ch(2), ch(3)];
}

function stretchValue(v, name) {
  if (name === "log" || name === "log10") return Math.log(1 + 999 * v) / Math.log(1000);
  if (name === "asinh") return Math.asinh(10 * v) / Math.asinh(10);
  return v;
}

/** A per-frame (frame-major) or static per-object value at a fractional frame. */
function frameValue(values, frames, count, i, frame, fallback) {
  if (!values) return fallback;
  if (!(frames > 1)) return values[i] === values[i] ? values[i] : fallback;
  const f = clamp(frame, 0, frames - 1);
  const f0 = Math.floor(f), f1 = Math.min(f0 + 1, frames - 1), t = f - f0;
  const a = values[f0 * count + i], b = values[f1 * count + i];
  const v = a * (1 - t) + b * t;
  return v === v ? v : a === a ? a : fallback;
}

async function pngFromRgba(rgba, width, height) {
  const c = document.createElement("canvas");
  c.width = width;
  c.height = height;
  c.getContext("2d").putImageData(new ImageData(rgba, width, height), 0, 0);
  return canvasPng(c);
}

async function canvasPng(canvas) {
  const blob = await new Promise((resolve) => canvas.toBlob(resolve, "image/png"));
  return new Uint8Array(await blob.arrayBuffer());
}

// ------------------------------------------------------------ time

/**
 * The timeline samples the animation steps through, in time order, at most
 * about `max` of them, always including the present-day frame.
 */
export function timeSamples(tl, max = MAX_SAMPLES) {
  const n = tl.count;
  if (n <= 1) return [{ frame: 0, time: tl.times[0] ?? 0 }];
  const step = Math.max(1, Math.ceil(n / max));
  const frames = new Set();
  for (let f = 0; f < n; f += step) frames.add(f);
  frames.add(n - 1);
  const zero = tl.zeroFrame?.();
  if (zero >= 0) frames.add(zero);
  return [...frames].map((f) => ({ frame: f, time: tl.times[f] })).sort((a, b) => a.time - b.time);
}

// ------------------------------------------------------------ objects

/**
 * Everything drawn as a point: its colour and brightness now, and for each
 * time sample its position (pc) and radius (pc; 0 before birth or when
 * absent). Hidden objects (filter, isolated lasso, hidden layers) stay out.
 */
function collectTracks(viewer, ui, samples, presentFrame) {
  const v = viewer;
  const params = v.points.params;
  const g = v.state.global;
  const out = [];
  for (const trace of v.traces) {
    const pts = trace.points;
    if (!pts || !trace.showInLegend) continue;
    const style = params?.styles.get(trace.key);
    if (!style || !style.visible) continue;
    const d = v.data.get(trace.key);
    const batch = v.points.byKey.get(trace.key);
    if (!d) continue;
    const N = pts.count;
    const byValue = style.colorMode === "by_value" && trace.colorBy;
    const lut = byValue ? ui.luts.get(style.colormap) : null;
    const override = style.colorOverride ? parseColor(style.colorOverride) : null;
    const base = parseColor(trace.color);
    const stars = trace.hasStars ? style.starsExp || 0 : 0;
    const sizeFrames = d.size ? pts.size.frames || 1 : 1;
    const opFrames = d.opacity ? pts.opacity.frames || 1 : 1;
    const scFrames = d.scalar ? pts.colorScalar?.frames || 1 : 1;
    const defSize = pts.size?.value ?? trace.pointSize ?? 1;
    const defOpacity = pts.opacity?.value ?? trace.opacity ?? 1;
    const sun = trace.key === v.world.sunTrace;
    for (let i = 0; i < N; i++) {
      const bits = batch ? batch.state[i] : 0;
      if (bits & 18) continue; // hidden by the filter or an isolated lasso
      const age = d.ageNow ? d.ageNow[i] : NaN;
      const starsK = stars && d.starsFactor ? Math.pow(Math.max(d.starsFactor[i], 1e-6), stars) : 1;
      const pos = new Array(samples.length);
      const size = new Float32Array(samples.length);
      let any = false;
      for (let k = 0; k < samples.length; k++) {
        const { frame, time } = samples[k];
        const p = v.objectPosition(trace.key, i, frame);
        if (!p) continue;
        pos[k] = p;
        const fade = birthFadeAt(time, age, g.fadeTime, g.fadeInOut);
        const sz = frameValue(d.size, sizeFrames, N, i, frame, defSize) * (style.sizeScale ?? 1) * g.pointSize * starsK;
        size[k] = 0.25 * sz * (v.world.pointScale || 1) * fade; // radius in pc
        if (size[k] > 0) any = true;
      }
      if (!any) continue;
      // Fill gaps so an absent object holds its last place while hidden.
      let last = pos.find(Boolean);
      for (let k = 0; k < samples.length; k++) { if (pos[k]) last = pos[k]; else pos[k] = last; }
      let rgb;
      if (lut) {
        const sc = frameValue(d.scalar, scFrames, N, i, presentFrame, pts.colorScalar?.value ?? 0);
        rgb = lutColor(lut, (sc - style.cmin) / ((style.cmax - style.cmin) || 1));
      } else if (override) rgb = override;
      else if (d.rgb) rgb = [d.rgb[i * 3] / 255, d.rgb[i * 3 + 1] / 255, d.rgb[i * 3 + 2] / 255];
      else rgb = base;
      let bright = frameValue(d.opacity, opFrames, N, i, presentFrame, defOpacity) * style.opacityScale * g.pointOpacity;
      if (bits & 9) bright *= style.dimOpacity ?? 0.16;
      out.push({ key: `${trace.key}_${i}`, pos, size, rgb, bright: clamp(bright, 0, 1), sun });
    }
  }
  return out;
}

// ------------------------------------------------------------ dust

/** The dust volumes drawn at the present day, as the renderer sets them up. */
function presentDustDraws(v) {
  const zero = v.timeline.zeroFrame?.();
  const frame = zero >= 0 ? zero : v.timeline.count - 1;
  return v._volumeParams(frame, v.timeline.times[frame] ?? 0).volumeDraws || [];
}

/**
 * Three orthogonal stacks of slices through each dust volume. Every slice
 * holds the optical depth of its slab, integrated with the renderer's model
 * (per-step alpha 1 − (1 − a·opacity)^(alphaCoef·step), so depth adds up as
 * −ln(1 − a·opacity)·alphaCoef per unit of texture length), and the
 * depth-weighted colour.
 */
async function dustSlices(v, ui, draws, toAr) {
  const nodes = [], materials = [], textures = {};
  let volumeIndex = 0;
  for (const d of draws) {
    const spec = v.volumeSpecs.find((x) => x.key === d.key);
    if (!spec) continue;
    const data = await v.store.get(spec.data.blob);
    const lut = ui.luts.get(d.state.colormap) || ui.luts.get(spec.colormaps?.[0]);
    if (!data || !lut) continue;
    const [nx, ny, nz] = spec.dims;
    const [lo, hi] = spec.bounds;
    // Per raw byte: optical depth per unit texture length, and colour value.
    const tauLut = new Float32Array(256), valLut = new Float32Array(256);
    const span = Math.max(d.high - d.low, 1e-6);
    const opacity = clamp(d.state.opacity ?? 1, 0, 1);
    for (let b = 0; b < 256; b++) {
      const w = (b / 255 - d.low) / span;
      if (w <= 0) continue;
      const val = stretchValue(Math.min(w, 0.999), d.state.stretch);
      const a = clamp(lutColor(lut, val)[3] * opacity, 0, 0.999);
      tauLut[b] = -Math.log(1 - a) * (d.state.alphaCoef ?? 1) * (d.fade ?? 1);
      valLut[b] = val;
    }
    const ext = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
    const widest = Math.max(...ext);
    const res = (e) => Math.max(8, Math.round((DUST_RES * e) / widest));
    const n = [nx, ny, nz];
    // Stacks: normal axis and the two in-plane axes (u, v), with their texels.
    const stacks = [
      { axis: 2, u: 0, v: 1 }, { axis: 0, u: 1, v: 2 }, { axis: 1, u: 0, v: 2 },
    ].map((st) => {
      const W = res(ext[st.u]), H = res(ext[st.v]), S = DUST_SLICES;
      const map = (count, bins) => Uint16Array.from({ length: count }, (_, i) => Math.min(bins - 1, Math.floor((i * bins) / count)));
      return {
        ...st, W, H, S,
        tau: new Float32Array(S * W * H), val: new Float32Array(S * W * H),
        sMap: map(n[st.axis], S), uMap: map(n[st.u], W), vMap: map(n[st.v], H),
        // Each texel averages the columns in its footprint; each voxel adds
        // its share of the slab's depth (one voxel = 1/n of texture length).
        norm: (W * H) / (n[st.u] * n[st.v]) / n[st.axis],
      };
    });
    const [sz, sx, sy] = stacks;
    for (let z = 0; z < nz; z++) {
      for (let y = 0; y < ny; y++) {
        const row = (z * ny + y) * nx;
        for (let x = 0; x < nx; x++) {
          const t = tauLut[data[row + x]];
          if (!t) continue;
          const val = valLut[data[row + x]] * t;
          let k = (sz.sMap[z] * sz.H + sz.vMap[y]) * sz.W + sz.uMap[x];
          sz.tau[k] += t; sz.val[k] += val;
          k = (sx.sMap[x] * sx.H + sx.vMap[z]) * sx.W + sx.uMap[y];
          sx.tau[k] += t; sx.val[k] += val;
          k = (sy.sMap[y] * sy.H + sy.vMap[z]) * sy.W + sy.uMap[x];
          sy.tau[k] += t; sy.val[k] += val;
        }
      }
    }
    const axisName = ["x", "y", "z"];
    const glow = lutColor(lut, 0.5);
    for (const st of stacks) {
      for (let slice = 0; slice < st.S; slice++) {
        const rgba = new Uint8ClampedArray(st.W * st.H * 4);
        let maxA = 0;
        for (let r = 0; r < st.H; r++) {
          // Image row 0 is the top: the high end of the in-plane v axis.
          const row = st.H - 1 - r;
          for (let c = 0; c < st.W; c++) {
            const k = (slice * st.H + row) * st.W + c;
            const tau = st.tau[k];
            if (!tau) continue;
            const alpha = 1 - Math.exp(-tau * st.norm * DUST_STACK_WEIGHT);
            const col = lutColor(lut, st.val[k] / tau);
            const o = (r * st.W + c) * 4;
            rgba[o] = col[0] * 255; rgba[o + 1] = col[1] * 255; rgba[o + 2] = col[2] * 255;
            rgba[o + 3] = alpha * 255;
            if (alpha > maxA) maxA = alpha;
          }
        }
        if (maxA < 2 / 255) continue; // an empty slab needs no slice
        const name = `Dust${volumeIndex}_${axisName[st.axis]}${String(slice).padStart(2, "0")}`;
        const file = `textures/${name}.png`;
        textures[file] = await pngFromRgba(rgba, st.W, st.H);
        // The slab's centre plane; corners in (u, v) order (lo,lo) (hi,lo) (hi,hi) (lo,hi)
        // to match st (0,0) (1,0) (1,1) (0,1), st's origin being the image's bottom-left.
        const at = lo[st.axis] + ((slice + 0.5) / st.S) * ext[st.axis];
        const corner = (uu, vv) => {
          const p = [0, 0, 0];
          p[st.axis] = at;
          p[st.u] = uu ? hi[st.u] : lo[st.u];
          p[st.v] = vv ? hi[st.v] : lo[st.v];
          return toAr(p);
        };
        const axisDir = [0, 0, 0];
        axisDir[st.axis] = 1;
        const mesh = quadMesh([corner(0, 0), corner(1, 0), corner(1, 1), corner(0, 1)], toAr(axisDir, true));
        materials.push({ name, texture: file, diffuse: [1, 1, 1], emissive: [glow[0] * 0.18, glow[1] * 0.18, glow[2] * 0.18], roughness: 0.95 });
        nodes.push({ name, material: name, mesh });
      }
    }
    volumeIndex++;
  }
  return { nodes, materials, textures };
}

// ------------------------------------------------------------ plate & labels

const FONT = "-apple-system, 'Helvetica Neue', Arial, sans-serif";

/**
 * The base plate texture: a dark translucent disc with distance rings
 * around the Sun, ring labels toward the front, the Galactic-centre
 * direction and the title. Texture x is +X; texture down is +Z (the front).
 */
async function plateTexture({ radiusPc, sunAt, title }) {
  const N = PLATE_RES, c = document.createElement("canvas");
  c.width = c.height = N;
  const ctx = c.getContext("2d");
  const R = N / 2 - 4, cx = N / 2, cy = N / 2;
  const px = (pc) => (pc / radiusPc) * R;
  const sx = cx + px(sunAt[0]), sy = cy + px(sunAt[1]);
  ctx.fillStyle = "rgba(7, 10, 17, 0.8)";
  ctx.beginPath(); ctx.arc(cx, cy, R, 0, Math.PI * 2); ctx.fill();
  ctx.strokeStyle = "rgba(255, 255, 255, 0.35)";
  ctx.lineWidth = 4;
  ctx.stroke();
  ctx.save();
  ctx.beginPath(); ctx.arc(cx, cy, R - 2, 0, Math.PI * 2); ctx.clip();
  const step = niceFloor(radiusPc / 3) || radiusPc / 3;
  ctx.font = `600 44px ${FONT}`;
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  for (let r = step; r < radiusPc * 1.4; r += step) {
    ctx.strokeStyle = "rgba(255, 255, 255, 0.2)";
    ctx.lineWidth = 3;
    ctx.beginPath(); ctx.arc(sx, sy, px(r), 0, Math.PI * 2); ctx.stroke();
    const a = Math.PI / 4; // front-right, where a viewer stands
    const lx = sx + Math.cos(a) * px(r), ly = sy + Math.sin(a) * px(r);
    if (Math.hypot(lx - cx, ly - cy) > R - 40) continue;
    const label = r >= 1000 ? `${+(r / 1000).toFixed(2)} kpc` : `${Math.round(r)} pc`;
    const w = ctx.measureText(label).width + 28;
    ctx.fillStyle = "rgba(7, 10, 17, 0.9)";
    ctx.fillRect(lx - w / 2, ly - 30, w, 60);
    ctx.fillStyle = "rgba(255, 255, 255, 0.78)";
    ctx.fillText(label, lx, ly);
  }
  ctx.restore();
  // The Sun, and the direction to the Galactic centre (+x).
  ctx.strokeStyle = "rgba(255, 224, 120, 0.9)";
  ctx.lineWidth = 5;
  ctx.beginPath(); ctx.arc(sx, sy, 26, 0, Math.PI * 2); ctx.stroke();
  ctx.strokeStyle = "rgba(255, 255, 255, 0.7)";
  ctx.fillStyle = "rgba(255, 255, 255, 0.8)";
  ctx.lineWidth = 6;
  ctx.beginPath(); ctx.moveTo(cx + R * 0.6, cy); ctx.lineTo(cx + R * 0.86, cy); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(cx + R * 0.91, cy); ctx.lineTo(cx + R * 0.85, cy - 20); ctx.lineTo(cx + R * 0.85, cy + 20); ctx.closePath(); ctx.fill();
  ctx.font = `600 38px ${FONT}`;
  ctx.fillText("Galactic centre", cx + R * 0.72, cy - 44);
  ctx.fillStyle = "rgba(255, 224, 120, 0.9)";
  ctx.fillText("Sun", sx, sy + 62);
  if (title) {
    ctx.font = `600 46px ${FONT}`;
    ctx.fillStyle = "rgba(255, 255, 255, 0.6)";
    ctx.fillText(title, cx, cy - R * 0.9);
  }
  return canvasPng(c);
}

async function labelTexture(text) {
  const c = document.createElement("canvas");
  c.width = 512;
  c.height = 128;
  const ctx = c.getContext("2d");
  ctx.fillStyle = "rgba(7, 10, 17, 0.88)";
  const r = 56;
  ctx.beginPath();
  ctx.moveTo(8 + r, 8);
  ctx.arcTo(504, 8, 504, 120, r); ctx.arcTo(504, 120, 8, 120, r); ctx.arcTo(8, 120, 8, 8, r); ctx.arcTo(8, 8, 504, 8, r);
  ctx.closePath();
  ctx.fill();
  ctx.fillStyle = "rgba(255, 255, 255, 0.96)";
  ctx.font = `700 64px ${FONT}`;
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.fillText(text, 256, 68);
  return canvasPng(c);
}

export function timeText(t) {
  if (Math.abs(t) < 1e-6) return "Today";
  const a = Math.abs(t);
  const s = a >= 100 ? a.toFixed(0) : String(+a.toFixed(1));
  return `${t < 0 ? "−" : "+"}${s} Myr`;
}

/**
 * Scale keys (time code, scale) that show a label from `from` to `to` and
 * hide it otherwise, with near-instant switches.
 */
export function windowKeys(from, to, end) {
  const keys = [];
  if (from > 0) keys.push([0, 0], [from - 0.01, 0]);
  keys.push([from, 1], [Math.max(from, to - 0.01), 1]);
  if (to < end) keys.push([to, 0], [end, 0]);
  return keys;
}

// ------------------------------------------------------------ model

/**
 * Build the USDZ for the current view: a time-lapse through the whole
 * timeline, or with `moment: true` a still of the current time. Returns
 * {blob, summary}: counts of objects, dust slices and time samples, the
 * animation length and the pc → m scale.
 */
export async function buildArModel(viewer, ui, { targetRadius = TARGET_RADIUS_M, moment = false } = {}) {
  const v = viewer;
  v._resolveParams?.(); // styles for exactly this state
  const tl = v.timeline;
  const samples = moment ? [{ frame: tl.frame, time: tl.time }] : timeSamples(tl);
  const zero = tl.zeroFrame?.();
  const presentFrame = moment ? tl.frame : zero >= 0 ? zero : tl.count - 1;
  let presentK = moment ? 0 : samples.findIndex((x) => x.frame === presentFrame);
  if (presentK < 0) presentK = samples.length - 1;
  const tracks = collectTracks(v, ui, samples, presentFrame);
  // A still shows the dust the viewer draws at that time; the time-lapse
  // reveals the present-day dust as it reaches today.
  const draws = moment ? v._volumeParams(tl.frame, tl.time).volumeDraws || [] : presentDustDraws(v);
  if (!tracks.length && !draws.length) throw new Error("Nothing is visible to place in AR.");

  // Fit like the figure frames itself: centred on the figure's own centre
  // (else the median object), with most present-day objects and the dust
  // inside the tabletop radius. Earlier in time objects may roam farther:
  // outside the view radius they are hidden until they come into view.
  const wc = v.world?.center;
  const nowPos = tracks.map((o) => o.pos[presentK]);
  const median = (k) => { const a = nowPos.map((p) => p[k]).sort((x, y) => x - y); return a[Math.floor(a.length / 2)] ?? 0; };
  const center = Array.isArray(wc) && wc.length === 3 && wc.every(Number.isFinite) ? wc.slice() : [median(0), median(1), median(2)];
  const dist = (p) => Math.hypot(p[0] - center[0], p[1] - center[1], p[2] - center[2]);
  const ds = nowPos.map(dist).sort((a, b) => a - b);
  let radius = ds.length ? ds[Math.floor(FIT_QUANTILE * (ds.length - 1))] : 0;
  let zmin = Infinity;
  for (const p of nowPos) zmin = Math.min(zmin, p[2]);
  for (const d of draws) {
    const spec = v.volumeSpecs.find((x) => x.key === d.key);
    if (!spec) continue;
    const [lo, hi] = spec.bounds;
    radius = Math.max(radius, Math.hypot(Math.max(Math.abs(lo[0] - center[0]), Math.abs(hi[0] - center[0])), Math.max(Math.abs(lo[1] - center[1]), Math.abs(hi[1] - center[1]))));
    zmin = Math.min(zmin, lo[2]);
  }
  radius = Math.max(radius, 1e-6);
  const viewRadius = radius * VIEW_RADII;
  for (const o of tracks) for (let k = 0; k < o.pos.length; k++) if (dist(o.pos[k]) > viewRadius) o.size[k] = 0;
  const s = targetRadius / radius;
  const lift = LIFT_M + (center[2] - (Number.isFinite(zmin) ? zmin : center[2])) * s;
  // Galactic (x, y, z) → USD (x, z, −y): right-handed, +Y up. Directions skip the offset.
  const toAr = (p, direction = false) => direction
    ? [p[0], p[2], -p[1]]
    : [(p[0] - center[0]) * s, (p[2] - center[2]) * s + lift, -(p[1] - center[1]) * s];

  // ---- time codes: a short hold at the start, the timeline, a long hold at the end
  const animated = samples.length > 1;
  const hold0 = animated ? Math.round(HOLD_START_S * TIME_CODES_PER_S) : 0;
  const tcOf = (k) => hold0 + k;
  const endTc = animated ? tcOf(samples.length - 1) + Math.round(HOLD_END_S * TIME_CODES_PER_S) : 0;
  const time = animated ? { start: 0, end: endTc, perSecond: TIME_CODES_PER_S } : null;

  const materials = [];
  const nodes = [];
  const textures = {};
  const geometries = { sphere: sphereMesh(1), sun: sphereMesh(2) };

  // ---- objects: one material per (colour, brightness); shared sphere geometry
  const groups = new Map();
  const same3 = (a, b) => Math.abs(a[0] - b[0]) < 1e-6 && Math.abs(a[1] - b[1]) < 1e-6 && Math.abs(a[2] - b[2]) < 1e-6;
  for (const o of tracks) {
    if (!o.size.some((x) => x > 0)) continue; // never within view
    const q = o.rgb.slice(0, 3).map((x) => Math.round(clamp(x, 0, 1) * 31));
    const level = Math.min(3, Math.floor(o.bright * 4));
    const key = o.sun ? "sun" : `${q.join("_")}_${level}`;
    if (!groups.has(key)) {
      const e = o.sun ? 1 : 0.3 + 0.7 * ((level + 0.5) / 4);
      const rgb = q.map((x) => x / 31);
      const name = o.sun ? "Sun" : `Points_${groups.size}`;
      materials.push({ name, diffuse: rgb.map((x) => x * 0.35), emissive: rgb.map((x) => x * e), roughness: 0.45 });
      groups.set(key, name);
    }
    const rad = Array.from(o.size, (pc) => (pc > 0 ? (o.sun ? Math.max(SUN_RADIUS_M, pc * s) : clamp(pc * s, MIN_RADIUS_M, MAX_RADIUS_M)) : 0));
    const ar = o.pos.map((p) => toAr(p));
    const node = { name: `O_${o.key}`, geometry: o.sun ? "sun" : "sphere", material: groups.get(key) };
    if (!animated || ar.every((p) => same3(p, ar[0]))) node.translate = ar[presentK];
    else node.translateSamples = [[0, ar[0]], ...ar.map((p, k) => [tcOf(k), p])];
    if (!animated || rad.every((r) => Math.abs(r - rad[0]) < 1e-7)) node.scale = rad[presentK];
    else {
      // Keys only where the size changes (and at both ends), for compactness.
      const keys = [[0, rad[0]]];
      for (let k = 0; k < rad.length; k++) {
        const prev = rad[k - 1] ?? rad[k], next = rad[k + 1] ?? rad[k];
        if (k === 0 || k === rad.length - 1 || Math.abs(rad[k] - prev) > 1e-7 || Math.abs(rad[k] - next) > 1e-7) keys.push([tcOf(k), rad[k]]);
      }
      node.scaleSamples = keys;
    }
    nodes.push(node);
  }

  // ---- dust: slices, revealed at the present day
  const dust = await dustSlices(v, ui, draws, toAr);
  if (dust.nodes.length) {
    // Grow from the figure's centre (the group's origin), not from the table.
    const pivot = toAr(center);
    for (const c of dust.nodes) {
      const p = c.mesh.points;
      for (let i = 0; i < p.length; i += 3) { p[i] -= pivot[0]; p[i + 1] -= pivot[1]; p[i + 2] -= pivot[2]; }
    }
    const dustNode = { name: "Dust", translate: pivot, children: dust.nodes };
    if (animated) {
      // Dust maps are present-day: they grow in as the timeline reaches
      // today (and shrink away again if it runs on into the future).
      const tp = tcOf(presentK);
      const keys = [[0, 0], [Math.max(0.01, tp - 6), 0], [tp, 1]];
      keys.push(presentK < samples.length - 1 ? [tcOf(presentK + 1), 0] : [endTc, 1]);
      if (presentK < samples.length - 1) keys.push([endTc, 0]);
      dustNode.scaleSamples = keys;
    }
    nodes.push(dustNode);
    materials.push(...dust.materials);
    Object.assign(textures, dust.textures);
  }

  // ---- the plate, and the time readout lying on its front edge
  const plateR = targetRadius * 1.15;
  const sunTrack = tracks.find((o) => o.sun);
  const sunPc = sunTrack ? sunTrack.pos[presentK] : center;
  const sunAt = [sunPc[0] - center[0], -(sunPc[1] - center[1])];
  textures["textures/plate.png"] = await plateTexture({ radiusPc: plateR / s, sunAt, title: v.manifest.title || "" });
  materials.push({ name: "Plate", texture: "textures/plate.png", diffuse: [1, 1, 1], emissive: [0.02, 0.025, 0.035], roughness: 0.9 });
  nodes.push({ name: "Plate", material: "Plate", mesh: discMesh(plateR, 0.002) });

  let labels = 0;
  // The time readout lies flat on the plate's front edge, readable from the front.
  const labelAt = (name) => ({
    name: "Label", material: name,
    mesh: quadMesh([[-0.08, 0, 0.02], [0.08, 0, 0.02], [0.08, 0, -0.02], [-0.08, 0, -0.02]], [0, 1, 0]),
  });
  if (moment && tl.count > 1) {
    // A still says when it is.
    const file = "textures/time_0.png";
    textures[file] = await labelTexture(timeText(tl.time));
    materials.push({ name: "Time_0", texture: file, diffuse: [1, 1, 1], emissive: [0.3, 0.3, 0.3], roughness: 0.9 });
    nodes.push({ name: "Time_0", translate: [0, 0.004, plateR * 0.84], children: [labelAt("Time_0")] });
    labels = 1;
  }
  if (animated) {
    const times = samples.map((x) => x.time);
    const t0 = times[0], t1 = times[times.length - 1];
    const step = niceFloor((t1 - t0) / 12) || (t1 - t0) / 12;
    const marks = [];
    for (let t = Math.ceil(t0 / step - 1e-9) * step; t <= t1 + 1e-9; t += step) marks.push(Math.abs(t) < 1e-9 ? 0 : t);
    // Time code of a time: interpolate between the samples around it.
    const tcAt = (t) => {
      let k = 0;
      while (k < times.length - 2 && times[k + 1] < t) k++;
      const f = times[k + 1] > times[k] ? (t - times[k]) / (times[k + 1] - times[k]) : 0;
      return tcOf(k + clamp(f, 0, 1));
    };
    const z0 = plateR * 0.84, y = 0.004;
    for (let j = 0; j < marks.length; j++) {
      const file = `textures/time_${j}.png`;
      textures[file] = await labelTexture(timeText(marks[j]));
      const name = `Time_${j}`;
      materials.push({ name, texture: file, diffuse: [1, 1, 1], emissive: [0.3, 0.3, 0.3], roughness: 0.9 });
      // Shown from halfway after the previous mark to halfway before the next.
      const from = j === 0 ? 0 : (tcAt(marks[j - 1]) + tcAt(marks[j])) / 2;
      const to = j === marks.length - 1 ? endTc : (tcAt(marks[j]) + tcAt(marks[j + 1])) / 2;
      nodes.push({ name, translate: [0, y, z0], scaleSamples: windowKeys(from, to, endTc), children: [labelAt(name)] });
      labels++;
    }
  }

  const title = v.manifest.title || "Oviz figure";
  const bytes = usdzPackage({ name: title, time, nodes, materials, geometries }, textures);
  return {
    blob: new Blob([bytes], { type: "model/vnd.usdz+zip" }),
    summary: {
      objects: nodes.filter((n) => n.geometry).length, dustSlices: dust.nodes.length, timeSamples: samples.length, labels,
      seconds: animated ? endTc / TIME_CODES_PER_S : 0, metresPerPc: s, bytes: bytes.length,
    },
  };
}
