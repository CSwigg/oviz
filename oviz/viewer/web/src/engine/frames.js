// Per-frame data textures and the GLSL that samples them.
//
// Positions for every frame live in one RGBA32F texture per channel set:
// texel k = frame * count + index, wrapped into rows of FRAME_TEX_WIDTH. The
// vertex shader interpolates between frames, so moving the time slider only
// changes a uniform: no CPU loops, no buffer uploads.

import { createDataTexture } from "./gl.js";

const FRAME_TEX_WIDTH = 4096;
const ABSENT = 3.0e38; // sentinel for "object not present in frame"

/**
 * Pack (frames, count, 3) positions plus an optional per-frame scalar into
 * an RGBA32F texture. Missing frames of the scalar replicate its static row.
 */
export function packFrameTexture(gl, { positions, posFrames, count, scalar = null, scalarFrames = 1, scalarDefault = 0, frames }) {
  const total = frames * count;
  const height = Math.max(1, Math.ceil(total / FRAME_TEX_WIDTH));
  const data = new Float32Array(FRAME_TEX_WIDTH * height * 4);
  for (let f = 0; f < frames; f++) {
    const pf = posFrames > 1 ? f : 0;
    const sf = scalarFrames > 1 ? f : 0;
    const pBase = pf * count * 3;
    const sBase = sf * count;
    const oBase = f * count * 4;
    for (let i = 0; i < count; i++) {
      const x = positions[pBase + i * 3];
      const o = oBase + i * 4;
      if (x !== x) { // NaN → absent
        data[o] = ABSENT; data[o + 1] = 0; data[o + 2] = 0;
      } else {
        data[o] = x;
        data[o + 1] = positions[pBase + i * 3 + 1];
        data[o + 2] = positions[pBase + i * 3 + 2];
      }
      const s = scalar ? scalar[sBase + i] : scalarDefault;
      data[o + 3] = s === s ? s : scalarDefault;
    }
  }
  const texture = createDataTexture(gl, data, FRAME_TEX_WIDTH, height, 4);
  return { texture, frames, count, height };
}

/** Pack up to four per-frame scalar channels into one RGBA32F texture. */
export function packScalarTexture(gl, { channels, count, frames }) {
  const total = frames * count;
  const height = Math.max(1, Math.ceil(total / FRAME_TEX_WIDTH));
  const data = new Float32Array(FRAME_TEX_WIDTH * height * 4);
  for (let c = 0; c < 4; c++) {
    const ch = channels[c];
    if (!ch) continue;
    const { values, frames: cf, fallback = 0 } = ch;
    for (let f = 0; f < frames; f++) {
      const src = (cf > 1 ? f : 0) * count;
      const dst = f * count * 4 + c;
      for (let i = 0; i < count; i++) {
        const v = values ? values[src + i] : fallback;
        data[dst + i * 4] = v === v ? v : fallback;
      }
    }
  }
  const texture = createDataTexture(gl, data, FRAME_TEX_WIDTH, height, 4);
  return { texture, frames, count, height };
}

/**
 * GLSL helpers. `frameFetch` reads one texel; `framePosition` interpolates
 * (Catmull-Rom when all four neighbours exist, else linear) and reports a
 * presence weight so objects that appear/disappear fade over one interval.
 */
export const FRAME_GLSL = `
const int FRAME_TEX_WIDTH = ${FRAME_TEX_WIDTH};
const float ABSENT_LIMIT = 1.0e37;

vec4 frameFetch(highp sampler2D tex, int frame, int index, int count) {
  int k = frame * count + index;
  return texelFetch(tex, ivec2(k % FRAME_TEX_WIDTH, k / FRAME_TEX_WIDTH), 0);
}

vec3 catmullRom(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float t) {
  float t2 = t * t;
  float t3 = t2 * t;
  return 0.5 * ((2.0 * p1) + (-p0 + p2) * t + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2 + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3);
}

// Returns xyz in .xyz, presence (0..1) in .w, and the interpolated
// scalar channel through 'scalarOut'.
vec4 framePosition(highp sampler2D tex, int index, int count, int frames, float frame, out float scalarOut) {
  if (frames <= 1) {
    vec4 a = frameFetch(tex, 0, index, count);
    scalarOut = a.w;
    return vec4(a.xyz, a.x > ABSENT_LIMIT ? 0.0 : 1.0);
  }
  float fc = clamp(frame, 0.0, float(frames - 1));
  int f1 = int(floor(fc));
  int f2 = min(f1 + 1, frames - 1);
  float t = fc - float(f1);
  vec4 a = frameFetch(tex, f1, index, count);
  vec4 b = frameFetch(tex, f2, index, count);
  bool hasA = a.x < ABSENT_LIMIT;
  bool hasB = b.x < ABSENT_LIMIT;
  if (!hasA && !hasB) {
    scalarOut = 0.0;
    return vec4(0.0);
  }
  if (!hasA) { scalarOut = b.w; return vec4(b.xyz, t); }
  if (!hasB) { scalarOut = a.w; return vec4(a.xyz, 1.0 - t); }
  scalarOut = mix(a.w, b.w, t);
  if (t <= 0.0) return vec4(a.xyz, 1.0);
  int f0 = max(f1 - 1, 0);
  int f3 = min(f2 + 1, frames - 1);
  vec4 p0 = frameFetch(tex, f0, index, count);
  vec4 p3 = frameFetch(tex, f3, index, count);
  if (p0.x > ABSENT_LIMIT || p3.x > ABSENT_LIMIT) {
    return vec4(mix(a.xyz, b.xyz, t), 1.0);
  }
  return vec4(catmullRom(p0.xyz, a.xyz, b.xyz, p3.xyz, t), 1.0);
}

vec4 frameScalars(highp sampler2D tex, int index, int count, int frames, float frame) {
  if (frames <= 1) return frameFetch(tex, 0, index, count);
  float fc = clamp(frame, 0.0, float(frames - 1));
  int f1 = int(floor(fc));
  int f2 = min(f1 + 1, frames - 1);
  return mix(frameFetch(tex, f1, index, count), frameFetch(tex, f2, index, count), fc - float(f1));
}
`;

/** CPU twin of framePosition for picking, labels and fly-to. */
export function cpuFramePosition(positions, posFrames, count, index, frame, out = [0, 0, 0]) {
  if (posFrames <= 1) {
    const o = index * 3;
    out[0] = positions[o]; out[1] = positions[o + 1]; out[2] = positions[o + 2];
    return out[0] === out[0] ? out : null;
  }
  const fc = Math.min(Math.max(frame, 0), posFrames - 1);
  const f1 = Math.floor(fc);
  const f2 = Math.min(f1 + 1, posFrames - 1);
  const t = fc - f1;
  const a = (f1 * count + index) * 3;
  const b = (f2 * count + index) * 3;
  const ax = positions[a], bx = positions[b];
  const hasA = ax === ax, hasB = bx === bx;
  if (!hasA && !hasB) return null;
  // Entering or leaving the catalogue the GPU fades the point by presence
  // (t, or 1 − t) and draws nothing at zero: labels and picks must agree.
  if (!hasA) {
    if (t <= 1e-3) return null;
    out[0] = bx; out[1] = positions[b + 1]; out[2] = positions[b + 2]; return out;
  }
  if (!hasB) {
    if (t >= 1 - 1e-3) return null;
    out[0] = ax; out[1] = positions[a + 1]; out[2] = positions[a + 2]; return out;
  }
  if (t <= 0) { out[0] = ax; out[1] = positions[a + 1]; out[2] = positions[a + 2]; return out; }
  const f0 = Math.max(f1 - 1, 0), f3 = Math.min(f2 + 1, posFrames - 1);
  const p0 = (f0 * count + index) * 3, p3 = (f3 * count + index) * 3;
  const t2 = t * t, t3 = t2 * t;
  const cr = positions[p0] === positions[p0] && positions[p3] === positions[p3];
  for (let k = 0; k < 3; k++) {
    const P1 = positions[a + k], P2 = positions[b + k];
    if (!cr) { out[k] = P1 + (P2 - P1) * t; continue; }
    const P0 = positions[p0 + k], P3 = positions[p3 + k];
    out[k] = 0.5 * (2 * P1 + (-P0 + P2) * t + (2 * P0 - 5 * P1 + 4 * P2 - P3) * t2 + (-P0 + 3 * P1 - 3 * P2 + P3) * t3);
  }
  return out;
}

/** Interpolate a per-frame (frames, 3) offset table on the CPU. */
export function frameOffset(offsets, frame, out = [0, 0, 0]) {
  if (!offsets) { out[0] = out[1] = out[2] = 0; return out; }
  const frames = offsets.length / 3;
  const fc = Math.min(Math.max(frame, 0), frames - 1);
  const f1 = Math.floor(fc), f2 = Math.min(f1 + 1, frames - 1), t = fc - f1;
  for (let k = 0; k < 3; k++) out[k] = offsets[f1 * 3 + k] + (offsets[f2 * 3 + k] - offsets[f1 * 3 + k]) * t;
  return out;
}

/**
 * Motion-trail span for a trail `myr` long: `span` is in frames, signed so
 * the trail points back along the playback `direction`; `myrPerFrame` is
 * signed by the time axis. Returns null when trails are off or undefined.
 */
export function trailParams(times, direction, myr) {
  const len = Number(myr) || 0;
  const n = times?.length || 0;
  if (len <= 0 || n < 2) return null;
  const myrPerFrame = (times[n - 1] - times[0]) / (n - 1);
  if (!(Math.abs(myrPerFrame) > 0)) return null;
  return { span: (len / Math.abs(myrPerFrame)) * (direction < 0 ? -1 : 1), myrPerFrame };
}
