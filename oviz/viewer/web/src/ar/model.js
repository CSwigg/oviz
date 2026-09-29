// AR model of the current view, for Apple's AR Quick Look (iPhone, iPad).
//
// What is on screen becomes a tabletop model: every visible object at the
// current time as a small glowing sphere in its displayed colour, and each
// visible dust volume as soft, crossed "cloudlet" quads where it is dense.
// Galactic x, y, z (pc) map to AR metres with Galactic north up (+Y) and
// the model resting a little above the detected table.

import { parseColor } from "../core/color.js";
import { clamp } from "../core/math.js";
import { birthFadeAt } from "../app/timeline.js";
import { usdzPackage } from "./usdz.js";

const TARGET_RADIUS_M = 0.4;   // most of the model fits within this radius
const FIT_QUANTILE = 0.9;      // …measured on this share of the objects
const KEEP_RADII = 1.8;        // far outliers beyond this many radii are left out
const LIFT_M = 0.02;           // gap between the table and the lowest object
const MIN_RADIUS_M = 0.0015;
const MAX_RADIUS_M = 0.009;
const SUN_RADIUS_M = 0.004;
const MAX_CLOUDLETS = 1500;
const CLOUD_BINS = 8;

// ------------------------------------------------------------ geometry

/** Unit icosphere (level 0 = icosahedron), outward CCW faces. */
export function icosphere(level = 0) {
  const t = (1 + Math.sqrt(5)) / 2;
  let verts = [
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

/** Growable mesh buffers. */
class MeshBuilder {
  constructor(name, material, { uvs = false, doubleSided = false } = {}) {
    this.name = name;
    this.material = material;
    this.doubleSided = doubleSided;
    this.p = [];
    this.n = [];
    this.uv = uvs ? [] : null;
    this.i = [];
  }

  sphere(geo, c, r) {
    const base = this.p.length / 3;
    const v = geo.verts;
    for (let k = 0; k < v.length; k += 3) {
      this.p.push(c[0] + v[k] * r, c[1] + v[k + 1] * r, c[2] + v[k + 2] * r);
      this.n.push(v[k], v[k + 1], v[k + 2]);
    }
    for (const f of geo.faces) this.i.push(base + f);
  }

  /** Three perpendicular soft quads through `c` (a cloudlet reads as volume from any side). */
  cloudlet(c, hx, hy, hz) {
    const quad = (corners, normal) => {
      const base = this.p.length / 3;
      for (const q of corners) { this.p.push(q[0], q[1], q[2]); this.n.push(normal[0], normal[1], normal[2]); }
      this.uv.push(0, 0, 1, 0, 1, 1, 0, 1);
      this.i.push(base, base + 1, base + 2, base, base + 2, base + 3);
    };
    const [x, y, z] = c;
    quad([[x - hx, y - hy, z], [x + hx, y - hy, z], [x + hx, y + hy, z], [x - hx, y + hy, z]], [0, 0, 1]);
    quad([[x - hx, y, z + hz], [x + hx, y, z + hz], [x + hx, y, z - hz], [x - hx, y, z - hz]], [0, 1, 0]);
    quad([[x, y - hy, z + hz], [x, y - hy, z - hz], [x, y + hy, z - hz], [x, y + hy, z + hz]], [1, 0, 0]);
  }

  finish() {
    return {
      name: this.name,
      material: this.material,
      doubleSided: this.doubleSided,
      points: Float32Array.from(this.p),
      normals: Float32Array.from(this.n),
      uvs: this.uv ? Float32Array.from(this.uv) : null,
      indices: Uint32Array.from(this.i),
    };
  }
}

// ------------------------------------------------------------ sampling

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

/** Everything drawn as a point right now: position (pc), colour, brightness, radius (pc). */
function collectObjects(viewer, ui) {
  const v = viewer;
  const params = v.points.params;
  const g = v.state.global;
  const t = v.timeline.time;
  const frame = v.timeline.frame;
  const out = [];
  for (const trace of v.traces) {
    const pts = trace.points;
    if (!pts || !trace.showInLegend) continue;
    const style = params?.styles.get(trace.key);
    if (!style || !style.visible || (style.presence ?? 1) <= 0.001) continue;
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
    for (let i = 0; i < N; i++) {
      const bits = batch ? batch.state[i] : 0;
      if (bits & 18) continue; // hidden by the filter or an isolated lasso
      const p = v.objectPosition(trace.key, i);
      if (!p) continue;
      const age = d.ageNow ? d.ageNow[i] : NaN;
      const fade = birthFadeAt(t, age, g.fadeTime, g.fadeInOut);
      if (fade <= 0.01) continue;
      let rgb;
      if (lut) {
        const s = frameValue(d.scalar, scFrames, N, i, frame, pts.colorScalar?.value ?? 0);
        const u = (s - style.cmin) / ((style.cmax - style.cmin) || 1);
        rgb = lutColor(lut, u);
      } else if (override) rgb = override;
      else if (d.rgb) rgb = [d.rgb[i * 3] / 255, d.rgb[i * 3 + 1] / 255, d.rgb[i * 3 + 2] / 255];
      else rgb = base;
      let bright = frameValue(d.opacity, opFrames, N, i, frame, defOpacity) * style.opacityScale * g.pointOpacity * (style.presence ?? 1);
      if (bits & 9) bright *= style.dimOpacity ?? 0.16;
      let size = frameValue(d.size, sizeFrames, N, i, frame, defSize) * (style.sizeScale ?? 1) * g.pointSize;
      if (g.fadeByOpacity) bright *= fade;
      else size *= fade;
      if (stars && d.starsFactor) size *= Math.pow(Math.max(d.starsFactor[i], 1e-6), stars);
      out.push({ p, rgb, bright: clamp(bright, 0, 1), radiusPc: 0.25 * size * (v.world.pointScale || 1), sun: trace.key === v.world.sunTrace });
    }
  }
  return out;
}

/** Soft samples of every dust volume drawn right now. */
async function collectClouds(viewer, ui, maxClouds) {
  const v = viewer;
  const draws = v.volumes.params?.volumeDraws || [];
  const out = [];
  for (const d of draws) {
    const spec = v.volumeSpecs.find((s) => s.key === d.key);
    if (!spec) continue;
    const data = await v.store.get(spec.data.blob);
    const lut = ui.luts.get(d.state.colormap) || ui.luts.get(spec.colormaps?.[0]);
    if (!data || !lut) continue;
    const [nx, ny, nz] = spec.dims;
    // A coarse grid of blocks: each keeps the mean density of a sparse
    // sample of its voxels (dense sheets read as clouds, not specks).
    const block = Math.max(2, Math.ceil(Math.cbrt((nx * ny * nz) / 60000)));
    const stride = Math.max(1, Math.floor(block / 3));
    const [lo, hi] = spec.bounds;
    const cell = [(hi[0] - lo[0]) / nx, (hi[1] - lo[1]) / ny, (hi[2] - lo[2]) / nz];
    const span = Math.max(d.high - d.low, 1e-6);
    const cx = (lo[0] + hi[0]) / 2, cy = (lo[1] + hi[1]) / 2;
    const ca = Math.cos(d.angle || 0), sa = Math.sin(d.angle || 0);
    const opacity = (d.state.opacity ?? 1) * (d.fade ?? 1);
    for (let bz = 0; bz < nz; bz += block) {
      for (let by = 0; by < ny; by += block) {
        for (let bx = 0; bx < nx; bx += block) {
          let sum = 0, n = 0;
          for (let z = bz; z < Math.min(bz + block, nz); z += stride) {
            for (let y = by; y < Math.min(by + block, ny); y += stride) {
              const row = (z * ny + y) * nx;
              for (let x = bx; x < Math.min(bx + block, nx); x += stride) {
                const w = (data[row + x] / 255 - d.low) / span;
                if (w > 0) sum += stretchValue(Math.min(w, 0.999), d.state.stretch);
                n++;
              }
            }
          }
          const w = n ? sum / n : 0;
          if (!(w > 0)) continue;
          const c = lutColor(lut, w);
          // Block centre in pc (co-rotating volumes turn about their centre).
          let px = lo[0] + (bx + Math.min(block, nx - bx) / 2) * cell[0];
          let py = lo[1] + (by + Math.min(block, ny - by) / 2) * cell[1];
          const pz = lo[2] + (bz + Math.min(block, nz - bz) / 2) * cell[2];
          if (d.angle) {
            const rx = px - cx, ry = py - cy;
            px = cx + rx * ca - ry * sa;
            py = cy + rx * sa + ry * ca;
          }
          out.push({ p: [px, py, pz], rgb: c, density: w, opacity, size: [block * cell[0], block * cell[1], block * cell[2]] });
        }
      }
    }
  }
  // The renderer integrates faint dust along each ray, so absolute block
  // densities are small: keep the densest blocks and weigh each against the
  // densest ones (a thin shell must still read as a shell).
  out.sort((a, b) => b.density - a.density);
  const top = out.length ? out[Math.min(out.length - 1, Math.floor(Math.min(out.length, maxClouds) * 0.02))].density : 1;
  // Relative, not absolute: dust windows differ between figures.
  const kept = out.filter((c) => c.density >= top * 0.03).slice(0, maxClouds);
  for (const c of kept) c.weight = clamp(c.density / top, 0, 1) * c.opacity;
  return kept;
}

/** A soft round sprite in `rgb`, peak opacity `alpha`, fading to the edge (PNG bytes). */
async function cloudTexture(rgb, alpha) {
  const c = document.createElement("canvas");
  c.width = c.height = 64;
  const ctx = c.getContext("2d");
  const g = ctx.createRadialGradient(32, 32, 0, 32, 32, 31.5);
  const [r, gg, b] = rgb.map((x) => Math.round(clamp(x, 0, 1) * 255));
  const stop = (at, a) => g.addColorStop(at, `rgba(${r},${gg},${b},${(a * alpha).toFixed(4)})`);
  stop(0, 1);
  stop(0.34, 0.8);
  stop(0.68, 0.25);
  stop(1, 0);
  ctx.fillStyle = g;
  ctx.fillRect(0, 0, 64, 64);
  const blob = await new Promise((resolve) => c.toBlob(resolve, "image/png"));
  return new Uint8Array(await blob.arrayBuffer());
}

// ------------------------------------------------------------ model

/**
 * Build the USDZ for the current view. Returns {blob, summary} where the
 * summary counts objects and cloudlets and gives the pc → m scale.
 */
export async function buildArModel(viewer, ui, { targetRadius = TARGET_RADIUS_M, maxClouds = MAX_CLOUDLETS } = {}) {
  const v = viewer;
  v._resolveParams?.(); // styles for exactly this time and state
  const objects = collectObjects(v, ui);
  const clouds = await collectClouds(v, ui, maxClouds);
  if (!objects.length && !clouds.length) throw new Error("Nothing is visible to place in AR.");

  // Fit like the figure frames itself: centred on the figure's own centre
  // (else the median object), with most objects and all the dust inside the
  // tabletop radius. A few far outliers must not shrink the rest to a speck.
  const all = [...objects.map((o) => o.p), ...clouds.map((c) => c.p)];
  const median = (k) => { const a = all.map((p) => p[k]).sort((x, y) => x - y); return a[Math.floor(a.length / 2)]; };
  const wc = v.world?.center;
  const center = Array.isArray(wc) && wc.length === 3 && wc.every(Number.isFinite) ? wc.slice() : [median(0), median(1), median(2)];
  const dist = (p) => Math.hypot(p[0] - center[0], p[1] - center[1], p[2] - center[2]);
  const od = objects.map((o) => dist(o.p)).sort((a, b) => a - b);
  let radius = od.length ? od[Math.floor(FIT_QUANTILE * (od.length - 1))] : 0;
  for (const c of clouds) radius = Math.max(radius, dist(c.p));
  radius = Math.max(radius, 1e-6);
  const keep = (p) => dist(p) <= radius * KEEP_RADII;
  const s = targetRadius / radius;
  let zmin = Infinity;
  for (const p of all) if (keep(p)) zmin = Math.min(zmin, p[2]);
  const lift = LIFT_M + (center[2] - zmin) * s;
  // Galactic (x, y, z) → USD (x, z, −y): right-handed, +Y up.
  const toAr = (p) => [(p[0] - center[0]) * s, (p[2] - center[2]) * s + lift, -(p[1] - center[1]) * s];

  const materials = [];
  const meshes = [];
  // Objects: one mesh per (colour, brightness) so each gets one material.
  // Small spheres read as round even as icosahedra; big ones need more faces.
  const small = icosphere(0), medium = icosphere(1), round = icosphere(2);
  const groups = new Map();
  let placed = 0;
  for (const o of objects) {
    if (!keep(o.p)) continue;
    placed++;
    const q = o.rgb.slice(0, 3).map((x) => Math.round(clamp(x, 0, 1) * 31));
    const level = Math.min(3, Math.floor(o.bright * 4));
    const key = o.sun ? "sun" : `${q.join("_")}_${level}`;
    if (!groups.has(key)) {
      const e = o.sun ? 1 : 0.3 + 0.7 * ((level + 0.5) / 4);
      const rgb = q.map((x) => x / 31);
      const name = `Points_${groups.size}`;
      materials.push({ name, diffuse: rgb.map((x) => x * 0.35), emissive: rgb.map((x) => x * e), roughness: 0.45 });
      groups.set(key, new MeshBuilder(name, name));
    }
    const r = o.sun ? Math.max(SUN_RADIUS_M, o.radiusPc * s) : clamp(o.radiusPc * s, MIN_RADIUS_M, MAX_RADIUS_M);
    groups.get(key).sphere(o.sun ? round : r > 0.004 ? medium : small, toAr(o.p), r);
  }
  for (const b of groups.values()) meshes.push(b.finish());

  // Dust: cloudlets binned by colour, faint and double-sided.
  const textures = {};
  if (clouds.length) {
    const bins = new Map();
    for (const c of clouds) {
      const bin = Math.min(CLOUD_BINS - 1, Math.floor(clamp(c.weight, 0, 1) * CLOUD_BINS));
      if (!bins.has(bin)) bins.set(bin, []);
      bins.get(bin).push(c);
    }
    for (const [bin, list] of bins) {
      const mean = (k) => list.reduce((acc, c) => acc + c.rgb[k], 0) / list.length;
      const rgb = [mean(0), mean(1), mean(2)];
      const name = `Dust_${bin}`;
      // Faint and layered, like the renderer's haze: dense sheets build up
      // where many cloudlets overlap, thin dust stays see-through.
      const file = `textures/dust_${bin}.png`;
      textures[file] = await cloudTexture(rgb, 0.06 + 0.24 * ((bin + 0.5) / CLOUD_BINS));
      materials.push({ name, texture: file, diffuse: rgb, emissive: rgb.map((x) => x * 0.12), roughness: 0.9 });
      const mb = new MeshBuilder(name, name, { uvs: true, doubleSided: true });
      for (const c of list) {
        const h = c.size.map((x) => x * s * 0.75); // overlap neighbours a little
        mb.cloudlet(toAr(c.p), h[0], h[2], h[1]);
      }
      meshes.push(mb.finish());
    }
  }
  const title = v.manifest.title || "Oviz figure";
  const bytes = usdzPackage({ name: title, meshes, materials }, textures);
  return {
    blob: new Blob([bytes], { type: "model/vnd.usdz+zip" }),
    summary: { objects: placed, clouds: clouds.length, metresPerPc: s, bytes: bytes.length },
  };
}
