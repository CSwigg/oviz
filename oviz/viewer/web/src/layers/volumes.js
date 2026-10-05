// Ray-marched scalar volumes (dust, ionised gas, densities).
//
// Optical model = the legacy Oviz renderer, so authored alpha coefficients
// keep their look: texture-space steps, per-step alpha
//   a = 1 − (1 − min(lut.a · opacity, 0.999))^(step · alpha_coef)
// front-to-back accumulation, and the legacy "premultiplied value blended
// with SRC_ALPHA" composite. What changes is the plumbing:
//   * one analytic fullscreen pass per volume (no proxy box geometry),
//   * block-max occupancy skipping (precomputed by the compiler),
//   * half-resolution marching while the view moves, full after it settles,
//   * a cached composite so hovers/labels never re-march the volume,
//   * time-series layers crossfade between neighbouring times.

import { createProgram, createTexture2D, RenderTarget, FULLSCREEN_VS } from "../engine/gl.js";

const FRAG = `
in vec2 vUv;
uniform mat4 uInvViewProj;
uniform vec3 uEye;
uniform vec3 uBoxMin;
uniform vec3 uBoxMax;
uniform vec3 uCenter;
uniform vec2 uRot;            // cos, sin of the co-rotation angle
uniform highp sampler3D uVolume;
uniform highp sampler3D uOcc;
// Optional colour field (e.g. velocity), stored as colour × density so
// trilinear filtering weights it by density: the LUT index is the ratio.
uniform highp sampler3D uColorVol;
uniform int uUseColor;
uniform vec3 uOccDims;
uniform int uUseOcc;
uniform sampler2D uCmap;
uniform float uLow;
uniform float uHigh;
uniform float uOpacity;
uniform float uSamples;
uniform float uAlphaCoef;
uniform int uStretch;         // 0 linear, 1 log10, 2 asinh
uniform float uFrameSeed;
uniform float uFade;
// Galactic-centre illumination (optional legacy lighting model).
uniform int uLighting;
uniform vec3 uGalCenter;
uniform float uGalIntensity;
uniform float uGalAmbient;
uniform float uGalExtinction;
uniform float uGalScattering;
uniform float uGalAnisotropy;
uniform float uGalWarmth;
// Lasso clip (classic "Lasso volumetric data"): samples that project
// outside the lasso, drawn with the camera of the moment it was drawn, fade
// to uMaskOutside (0 hides them). A State change between selections blends
// from the clip drawn when it began (…From) to the one it brings: what
// leaves fades out early in the flight, what joins fades in as it lands
// (uMaskMix). A change that interrupts another starts from that blend as
// it stood: …From itself blends from …Prev, held at uMaskFromMix.
// uMaskOn bits: 1 a clip now, 2 blending from another, 4 that one clips,
// 8 that one is a blend from a third, 16 the third clips.
uniform int uMaskOn;
uniform mat4 uMaskVP;
uniform sampler2D uMask;
uniform float uMaskOutside;
uniform mat4 uMaskFromVP;
uniform sampler2D uMaskFrom;
uniform float uMaskFromOutside;
uniform mat4 uMaskPrevVP;
uniform sampler2D uMaskPrev;
uniform float uMaskPrevOutside;
uniform vec2 uMaskMix;
uniform vec2 uMaskFromMix;
out vec4 outColor;

vec3 rotZ(vec3 p, vec2 cs) {
  return vec3(cs.x * p.x - cs.y * p.y, cs.y * p.x + cs.x * p.y, p.z);
}

float asinhF(float v) { return log(v + sqrt(v * v + 1.0)); }

float clipWeight(vec3 world, mat4 vp, sampler2D mask, float outside) {
  vec4 q = vp * vec4(world, 1.0);
  if (q.w <= 0.0) return outside;
  vec2 uv = q.xy / q.w * 0.5 + 0.5;
  if (any(lessThan(uv, vec2(0.0))) || any(greaterThan(uv, vec2(1.0)))) return outside;
  return texture(mask, uv).r > 0.5 ? 1.0 : outside;
}

float maskWeight(vec3 tc, vec3 scale) {
  vec3 world = rotZ(uBoxMin + tc * scale - uCenter, uRot) + uCenter;
  float b = (uMaskOn & 1) != 0 ? clipWeight(world, uMaskVP, uMask, uMaskOutside) : 1.0;
  if ((uMaskOn & 2) == 0) return b;
  float a = (uMaskOn & 4) != 0 ? clipWeight(world, uMaskFromVP, uMaskFrom, uMaskFromOutside) : 1.0;
  if ((uMaskOn & 8) != 0) {
    float p = (uMaskOn & 16) != 0 ? clipWeight(world, uMaskPrevVP, uMaskPrev, uMaskPrevOutside) : 1.0;
    a = mix(p, a, a < p ? uMaskFromMix.x : uMaskFromMix.y);
  }
  return mix(a, b, b < a ? uMaskMix.x : uMaskMix.y);
}

float stretch(float v) {
  v = clamp(v, 0.0, 1.0);
  if (uStretch == 1) return log(1.0 + 999.0 * v) / log(1000.0);
  if (uStretch == 2) return asinhF(10.0 * v) / asinhF(10.0);
  return v;
}

float density(vec3 tc) {
  if (any(lessThan(tc, vec3(0.0))) || any(greaterThan(tc, vec3(1.0)))) return 0.0;
  return clamp((texture(uVolume, tc).r - uLow) / max(uHigh - uLow, 1e-6), 0.0, 1.0);
}

vec3 galacticLight(vec3 premul, float a, float raw, float local, vec3 tc, vec3 scale, vec3 eyeLocal, vec3 gcLocal) {
  vec3 samplePhysical = (tc - 0.5) * scale + uCenter;
  vec3 toCenter = gcLocal - samplePhysical;
  float cd = max(length(toCenter), 1.0);
  vec3 tcDir = toCenter / cd;
  vec3 outgoing = normalize(eyeLocal - samplePhysical);
  float soft = sqrt(cd * cd + 1500.0 * 1500.0);
  float solar = max(length(uGalCenter), 1000.0);
  float falloff = clamp((solar * solar) / (soft * soft), 0.18, 4.0);
  vec3 gcTex = (gcLocal - uCenter) / max(scale, vec3(1e-6)) + 0.5;
  vec3 cDelta = gcTex - tc;
  float cDist = length(cDelta);
  float trans = 1.0;
  if (cDist > 1e-6 && uGalExtinction > 1e-6) {
    vec3 d = cDelta / cDist;
    float column = 0.5 * local + 0.32 * density(tc + d * 0.055) + 0.18 * density(tc + d * 0.16);
    trans = exp(-uGalExtinction * column * clamp(cDist, 0.2, 1.25));
  }
  float gs = 0.01;
  vec3 g = vec3(raw - texture(uVolume, tc + vec3(gs, 0, 0)).r, raw - texture(uVolume, tc + vec3(0, gs, 0)).r, raw - texture(uVolume, tc + vec3(0, 0, gs)).r);
  vec3 n = normalize(g / max(scale, vec3(1e-6)) + 1e-9);
  float facing = 0.30 + 0.70 * abs(dot(n, tcDir));
  float direct = uGalIntensity * falloff * trans * facing;
  float gA = clamp(uGalAnisotropy, 0.0, 0.9);
  float den = max(1.0 + gA * gA - 2.0 * gA * dot(-tcDir, outgoing), 0.04);
  float phase = clamp((1.0 - gA * gA) / pow(den, 1.5), 0.15, 4.5);
  float scattered = uGalScattering * falloff * trans * phase * 0.28;
  float field = clamp(uGalAmbient + direct + scattered, 0.025, 3.0);
  vec3 orig = premul / max(a, 1e-6);
  float lum = dot(orig, vec3(0.2126, 0.7152, 0.0722));
  vec3 incident = mix(vec3(0.72, 0.82, 1.0), vec3(1.0, 0.63, 0.30), uGalWarmth);
  vec3 lit = mix(vec3(0.24, 0.105, 0.045), incident, clamp(trans * (0.35 + 0.65 * field), 0.0, 1.0));
  vec3 dust = mix(orig, lit, 0.78);
  return dust * lum * field * mix(0.42, 1.0, trans) * a;
}

void main() {
  vec2 ndc = vUv * 2.0 - 1.0;
  vec4 nearP = uInvViewProj * vec4(ndc, -1.0, 1.0);
  vec4 farP = uInvViewProj * vec4(ndc, 1.0, 1.0);
  vec3 dirW = normalize(farP.xyz / farP.w - nearP.xyz / nearP.w);
  // Into the volume's (un-rotated) frame.
  vec2 inv = vec2(uRot.x, -uRot.y);
  vec3 eye = rotZ(uEye - uCenter, inv) + uCenter;
  vec3 dir = rotZ(dirW, inv);
  vec3 invDir = 1.0 / dir;
  vec3 t0 = (uBoxMin - eye) * invDir;
  vec3 t1 = (uBoxMax - eye) * invDir;
  vec3 tlo = min(t0, t1), thi = max(t0, t1);
  float tmin = max(max(tlo.x, tlo.y), tlo.z);
  float tmax = min(min(thi.x, thi.y), thi.z);
  if (tmax <= max(tmin, 0.0)) discard;
  tmin = max(tmin, 0.0);
  vec3 scale = uBoxMax - uBoxMin;
  vec3 tcStart = (eye + dir * tmin - uBoxMin) / scale;
  vec3 tcEnd = (eye + dir * tmax - uBoxMin) / scale;
  vec3 delta = tcEnd - tcStart;
  int count = min(int(length(delta) * uSamples), int(uSamples * 1.8));
  if (count <= 0) discard;
  delta /= float(count);
  float jitter = fract(52.9829189 * fract(dot(gl_FragCoord.xy + uFrameSeed, vec2(0.06711056, 0.00583715))));
  vec3 tc = tcStart - delta * (0.01 + 0.98 * jitter) - delta;
  float stepLen = length(delta);
  float stepAlpha = stepLen * uAlphaCoef;
  float invRange = 1.0 / max(uHigh - uLow, 1e-6);
  vec3 gcLocal = rotZ(uGalCenter - uCenter, inv) + uCenter;
  vec4 acc = vec4(0.0);
  for (int i = 0; i < 4096; i++) {
    if (i >= count) break;
    tc += delta;
    if (uUseOcc == 1 && texture(uOcc, tc).r - uLow <= 0.0) {
      vec3 cell = floor(clamp(tc, vec3(0.0), vec3(0.999999)) * uOccDims);
      vec3 exitPlanes = (cell + vec3(greaterThanEqual(delta, vec3(0.0)))) / uOccDims;
      vec3 toExit = mix(vec3(65536.0), (exitPlanes - tc) / delta, greaterThan(abs(delta), vec3(1e-12)));
      int skip = int(clamp(min(min(toExit.x, toExit.y), toExit.z), 0.0, 4096.0));
      i += skip;
      tc += float(skip) * delta;
      continue;
    }
    float raw = texture(uVolume, tc).r;
    float v = (raw - uLow) * invRange;
    if (v <= 0.0) continue;
    float mw = 1.0;
    if (uMaskOn != 0) {
      mw = maskWeight(tc, scale);
      if (mw <= 0.0) continue;
    }
    v = stretch(min(v, 0.999));
    float cv = uUseColor == 1 ? clamp(texture(uColorVol, tc).r / max(raw, 1e-6), 0.0, 1.0) : v;
    vec4 c = texture(uCmap, vec2(cv, 0.5));
    c.a *= mw;
    float a = 1.0 - pow(1.0 - clamp(c.a * uOpacity, 0.0, 0.999), stepAlpha);
    a *= (1.0 - acc.a);
    if (a <= 0.0) continue;
    vec3 rgb = c.rgb * a;
    if (uLighting == 1) rgb = galacticLight(rgb, a, raw, clamp(v, 0.0, 1.0), tc, scale, eye, gcLocal);
    acc += vec4(rgb, a);
    if (acc.a >= 0.99) { acc.a = 1.0; break; }
  }
  if (acc.a <= 0.001) discard;
  acc.rgb += vec3((jitter - 0.5) / (255.0 * max(acc.a, 0.05)));
  acc.rgb = max(acc.rgb, vec3(0.0));
  acc *= uFade;
  // Store the legacy composite contribution: rgb·a over (1 − a).
  outColor = vec4(acc.rgb * acc.a, acc.a);
}
`;

const COMPOSITE = `
in vec2 vUv;
uniform sampler2D uTex;
out vec4 outColor;
void main() { outColor = texture(uTex, vUv); }
`;

// Rebuilt densities kept on the GPU per volume (a 160×160×32 half-float
// sample is 1.6 MB; the classic runtime kept 24 float32 ones).
const VOLUME_DATA_CACHE = 32;

function stretchIndex(name) {
  const s = String(name || "linear").toLowerCase();
  if (s === "log10" || s === "log") return 1;
  if (s === "asinh" || s === "arcsinh") return 2;
  return 0;
}

export class VolumesLayer {
  constructor(gl, { caps }) {
    this.gl = gl;
    this.caps = caps;
    this.order = 10;
    this.program = createProgram(gl, FULLSCREEN_VS, FRAG, { label: "volume" });
    this.composite = createProgram(gl, FULLSCREEN_VS, COMPOSITE, { label: "volume-composite" });
    this.target = new RenderTarget(gl, {
      internalFormat: caps.halfFloatRT || caps.floatRT ? gl.RGBA16F : gl.RGBA8,
      type: caps.halfFloatRT || caps.floatRT ? gl.HALF_FLOAT : gl.UNSIGNED_BYTE,
    });
    this.volumes = new Map(); // key → gpu volume
    this.cmaps = new Map();
    this.params = null;
    this._cacheKey = "";
    this.lastMarchMs = 0;
    this.quality = 1;
    this.mask = null; // lasso clip: {texture, vp, outside, key}
    // While a State change blends selections: the clip drawn when it began
    // ({} for none), itself part-way from `prev` (held at `prevMix`) when
    // the change interrupted another.
    this.maskFrom = null;
    this.maskMix = [1, 1];
    this.maskVersion = 0;
  }

  _blankMask() {
    if (this._blank) return this._blank;
    const gl = this.gl;
    const t = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, t);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, 1, 1, 0, gl.RED, gl.UNSIGNED_BYTE, new Uint8Array([255]));
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.bindTexture(gl.TEXTURE_2D, null);
    this._blank = t;
    return t;
  }

  /**
   * Clip the volumes to a lasso: `polygon` in normalised device coordinates
   * (x right, y up, −1…1), `vp` the view-projection matrix it was drawn
   * with, `outside` the weight left outside (0 hides). Null clears it. A
   * blend in progress ends; an unchanged clip keeps its texture.
   */
  setMask(mask) {
    const key = clipKey(mask);
    const keep = this.mask && this.mask.key === key ? this.mask : null;
    if (keep && !this.maskFrom) return;
    this._dropMasks(keep);
    this.mask = keep || this._clip(mask, key);
  }

  /**
   * Blend from the clip drawn now (a clip, none, or a blend part-way) to
   * `to` (null: no clip) while a State change plays; `setMaskMix` moves it
   * along and `setMask` ends it. False when nothing would change.
   */
  blendMask(to) {
    const key = clipKey(to);
    if (!this.maskFrom && (this.mask?.key ?? null) === key) return false;
    const now = this.mask || {}, f = this.maskFrom;
    let from;
    if (!f) from = { ...now, prev: null };
    else if (!f.prev) from = { ...now, prev: f, prevMix: [...this.maskMix] };
    else {
      // Two blends deep at most: fold away whichever part is drawn least.
      const mo = (this.maskMix[0] + this.maskMix[1]) / 2;
      const mi = (f.prevMix[0] + f.prevMix[1]) / 2;
      if (mo < (1 - mo) * Math.min(mi, 1 - mi)) {
        // The latest blend has barely begun: carry on from the one before.
        this._deleteClip(now);
        from = f;
      } else {
        // The earlier blend is near one of its ends: keep that end.
        const end = mi >= 0.5 ? f : f.prev;
        this._deleteClip(end === f ? f.prev : f);
        from = { ...now, prev: { texture: end.texture, vp: end.vp, outside: end.outside, key: end.key }, prevMix: [...this.maskMix] };
      }
    }
    this.maskFrom = from;
    this.mask = this._clip(to, key);
    this.maskMix = [0, 0];
    this.maskVersion++;
    this._cacheKey = "";
    return true;
  }

  /** How far a clip blend has got: fade-out and fade-in, 0–1. */
  setMaskMix(out, inn) {
    if (!this.maskFrom || (this.maskMix[0] === out && this.maskMix[1] === inn)) return;
    this.maskMix = [out, inn];
    this.maskVersion++;
  }

  _deleteClip(c) {
    if (c?.texture) this.gl.deleteTexture(c.texture);
  }

  _dropMasks(keep = null) {
    for (const c of [this.mask, this.maskFrom, this.maskFrom?.prev]) {
      if (c && c.texture !== keep?.texture) this._deleteClip(c);
    }
    this.mask = null;
    this.maskFrom = null;
    this.maskMix = [1, 1];
    this.maskVersion++;
    this._cacheKey = "";
  }

  /** A lasso (`polygon` in NDC, the `vp` it was drawn with, `outside` weight) as a 512² clip texture; null if none. */
  _clip(mask, key = clipKey(mask)) {
    const gl = this.gl;
    if (key == null) return null;
    const N = 512;
    const canvas = document.createElement("canvas");
    canvas.width = canvas.height = N;
    const ctx = canvas.getContext("2d");
    ctx.fillStyle = "#fff";
    ctx.beginPath();
    mask.polygon.forEach(([x, y], i) => {
      const px = (x * 0.5 + 0.5) * N, py = (0.5 - y * 0.5) * N;
      if (i) ctx.lineTo(px, py); else ctx.moveTo(px, py);
    });
    ctx.closePath();
    ctx.fill();
    const img = ctx.getImageData(0, 0, N, N).data;
    // Texture rows run bottom-up (v = 0 at ndc y = −1); the canvas top-down.
    const data = new Uint8Array(N * N);
    for (let row = 0; row < N; row++) {
      const src = row * N * 4, dst = (N - 1 - row) * N;
      for (let col = 0; col < N; col++) data[dst + col] = img[src + col * 4];
    }
    const texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, N, N, 0, gl.RED, gl.UNSIGNED_BYTE, data);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.bindTexture(gl.TEXTURE_2D, null);
    return { texture, vp: Float32Array.from(mask.vp), outside: clipOutside(mask), key };
  }

  get hasContent() {
    return this.volumes.size > 0;
  }

  addColormap(name, bytes, width) {
    if (this.cmaps.has(name)) return;
    const gl = this.gl;
    this.cmaps.set(name, createTexture2D(gl, {
      width, height: 1, internalFormat: gl.RGBA8, format: gl.RGBA, type: gl.UNSIGNED_BYTE, data: bytes,
    }));
  }

  /**
   * Upload a volume. `data` is Uint8Array (z, y, x) or Float32Array of
   * values in [0, 1] (stored half-float); `occ` optional block maxima;
   * `color` an optional Uint8Array colour field of the same shape;
   * `dataId` names these voxels in the cache `useData` swaps through.
   */
  addVolume(spec, data, occ, color = null, dataId = "") {
    const gl = this.gl;
    // Replacing a volume frees the old textures (never leak GPU memory).
    const old = this.volumes.get(spec.key);
    if (old) this._freeVolume(old);
    const gpu = this._upload(spec, data, occ);
    let colorTex = null;
    if (color) colorTex = this._texture3D(spec, color, spec.interpolation === false ? gl.NEAREST : gl.LINEAR);
    const [lo, hi] = spec.bounds;
    const vol = {
      key: spec.key, spec, ...gpu, colorTex,
      center: [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2],
      boxMin: lo, boxMax: hi,
      // Rebuilt densities (event KDEs) by id: one texture per time sample.
      cache: new Map(), dataId: "",
    };
    if (dataId) {
      vol.cache.set(dataId, gpu);
      vol.dataId = dataId;
    }
    this.volumes.set(spec.key, vol);
    this._cacheKey = "";
    return vol;
  }

  /**
   * Show the voxels cached as `id` for volume `key`; with `data`, upload
   * and cache them first. False when `id` is not cached and no data was
   * given. Cached textures are capped (least recently used go first).
   */
  useData(key, id, data = null) {
    const vol = this.volumes.get(key);
    if (!vol) return false;
    let gpu = vol.cache.get(id);
    if (!gpu) {
      if (!data) return false;
      gpu = this._upload(vol.spec, data, null);
      vol.cache.set(id, gpu);
      while (vol.cache.size > VOLUME_DATA_CACHE) {
        const [oldest, g] = vol.cache.entries().next().value;
        if (oldest === vol.dataId) { vol.cache.delete(oldest); vol.cache.set(oldest, g); continue; }
        this._freeGpu(g);
        vol.cache.delete(oldest);
      }
    } else {
      vol.cache.delete(id);
      vol.cache.set(id, gpu);
    }
    if (vol.dataId !== id) {
      if (!vol.dataId) this._freeGpu(vol);
      vol.tex = gpu.tex;
      vol.occTex = gpu.occTex;
      vol.occDims = gpu.occDims;
      vol.dataId = id;
      this._cacheKey = "";
    }
    return true;
  }

  hasData(key, id) {
    return !!this.volumes.get(key)?.cache.has(id);
  }

  _texture3D(spec, data, filter) {
    const gl = this.gl;
    const [nx, ny, nz] = spec.dims;
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_3D, tex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    if (data instanceof Float32Array) gl.texImage3D(gl.TEXTURE_3D, 0, gl.R16F, nx, ny, nz, 0, gl.RED, gl.FLOAT, data);
    else gl.texImage3D(gl.TEXTURE_3D, 0, gl.R8, nx, ny, nz, 0, gl.RED, gl.UNSIGNED_BYTE, data);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, filter);
    for (const w of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T, gl.TEXTURE_WRAP_R]) gl.texParameteri(gl.TEXTURE_3D, w, gl.CLAMP_TO_EDGE);
    gl.bindTexture(gl.TEXTURE_3D, null);
    return tex;
  }

  _upload(spec, data, occ) {
    const gl = this.gl;
    const [nx, ny, nz] = spec.dims;
    const tex = this._texture3D(spec, data, spec.interpolation === false ? gl.NEAREST : gl.LINEAR);
    let occTex = null, occDims = [1, 1, 1];
    const grid = occ || computeOccupancy(data, nx, ny, nz);
    if (grid) {
      occTex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_3D, occTex);
      gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
      gl.texImage3D(gl.TEXTURE_3D, 0, gl.R8, grid.gx, grid.gy, grid.gz, 0, gl.RED, gl.UNSIGNED_BYTE, grid.data);
      gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      for (const w of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T, gl.TEXTURE_WRAP_R]) gl.texParameteri(gl.TEXTURE_3D, w, gl.CLAMP_TO_EDGE);
      gl.bindTexture(gl.TEXTURE_3D, null);
      occDims = [grid.gx, grid.gy, grid.gz];
    }
    return { tex, occTex, occDims };
  }

  _freeGpu(g) {
    if (g.tex) this.gl.deleteTexture(g.tex);
    if (g.occTex) this.gl.deleteTexture(g.occTex);
  }

  _freeVolume(vol) {
    // The drawn texture is either the volume's own or one of its cached ones.
    if (!vol.dataId) this._freeGpu(vol);
    for (const g of vol.cache.values()) this._freeGpu(g);
    vol.cache.clear();
    if (vol.colorTex) this.gl.deleteTexture(vol.colorTex);
  }

  prepare(frame) {
    const p = this.params;
    if (!p || !this.volumes.size) return;
    const gl = this.gl;
    const draws = p.volumeDraws || [];
    if (!draws.length) { this._visible = false; return; }
    this._visible = true;
    // Quality: march at half resolution while the view moves. Large
    // captures cap the march target (smooth volumes upsample cleanly and a
    // full 8K float target would need hundreds of MB).
    let scale = frame.interacting ? Math.min(0.5, this.quality) : this.quality;
    if (frame.capture) scale = Math.min(scale, 4096 / Math.max(frame.width, frame.height));
    const w = Math.max(1, Math.round(frame.width * scale));
    const h = Math.max(1, Math.round(frame.height * scale));
    const key = [frame.camera.version, w, h, p.frame.toFixed(5), p.volumeSignature, frame.width, frame.height, this.maskVersion].join("|");
    this._resized = this.target.resize(w, h);
    if (key === this._cacheKey && !this._resized) return;
    this._cacheKey = key;
    const t0 = performance.now();
    this.target.bind();
    gl.clearColor(0, 0, 0, 0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.disable(gl.DEPTH_TEST);
    gl.depthMask(false);
    gl.enable(gl.BLEND);
    // Accumulate volumes back-to-front exactly like sequential legacy draws.
    gl.blendFuncSeparate(gl.ONE, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    const prog = this.program.use();
    const cam = frame.camera;
    prog.m4("uInvViewProj", cam.invViewProj).v3("uEye", cam.eye);
    prog.f("uFrameSeed", 0);
    const mask = this.mask, from = this.maskFrom, prev = from?.prev;
    prog.i("uMaskOn", (mask ? 1 : 0) | (from ? 2 : 0) | (from?.texture ? 4 : 0) | (prev ? 8 : 0) | (prev?.texture ? 16 : 0));
    prog.f("uMaskOutside", mask ? mask.outside : 1).f("uMaskFromOutside", from?.texture ? from.outside : 1).f("uMaskPrevOutside", prev?.texture ? prev.outside : 1);
    if (mask) prog.m4("uMaskVP", mask.vp);
    if (from?.texture) prog.m4("uMaskFromVP", from.vp);
    if (prev?.texture) prog.m4("uMaskPrevVP", prev.vp);
    prog.v2("uMaskMix", this.maskMix[0], this.maskMix[1]);
    prog.v2("uMaskFromMix", from?.prevMix?.[0] ?? 1, from?.prevMix?.[1] ?? 1);
    // Always give the 2D mask samplers their own units: left unbound they
    // would share unit 0 with the 3D volume, and WebGL rejects the draw.
    prog.tex("uMask", mask ? mask.texture : this._blankMask());
    prog.tex("uMaskFrom", from?.texture ? from.texture : this._blankMask());
    prog.tex("uMaskPrev", prev?.texture ? prev.texture : this._blankMask());
    const sorted = draws
      .map((d) => ({ d, vol: this.volumes.get(d.key) }))
      .filter((x) => x.vol)
      .sort((a, b) => dist2(b.vol.center, cam.eye) - dist2(a.vol.center, cam.eye));
    for (const { d, vol } of sorted) {
      const st = d.state;
      const cm = this.cmaps.get(st.colormap) || this.cmaps.get(vol.spec.colormaps[0]);
      if (!cm) continue;
      prog.v3("uBoxMin", vol.boxMin).v3("uBoxMax", vol.boxMax).v3("uCenter", vol.center);
      prog.v2("uRot", Math.cos(d.angle || 0), Math.sin(d.angle || 0));
      prog.tex("uVolume", vol.tex, gl.TEXTURE_3D);
      prog.tex("uOcc", vol.occTex || vol.tex, gl.TEXTURE_3D);
      prog.tex("uColorVol", vol.colorTex || vol.tex, gl.TEXTURE_3D).i("uUseColor", vol.colorTex ? 1 : 0);
      prog.v3("uOccDims", vol.occDims).i("uUseOcc", vol.occTex && d.canSkip ? 1 : 0);
      prog.tex("uCmap", cm);
      prog.f("uLow", d.low).f("uHigh", d.high).f("uOpacity", st.opacity).f("uSamples", st.steps);
      prog.f("uAlphaCoef", st.alphaCoef).i("uStretch", stretchIndex(st.stretch)).f("uFade", d.fade);
      prog.i("uLighting", st.lightingMode === "galactic" ? 1 : 0);
      const gc = st.galacticCenter || [8122, 0, 0];
      prog.v3("uGalCenter", gc[0], gc[1], gc[2]);
      prog.f("uGalIntensity", st.galacticLightIntensity ?? 1.35).f("uGalAmbient", st.galacticAmbient ?? 0.22);
      prog.f("uGalExtinction", st.galacticExtinction ?? 2.4).f("uGalScattering", st.galacticScattering ?? 0.55);
      prog.f("uGalAnisotropy", st.galacticAnisotropy ?? 0.45).f("uGalWarmth", st.galacticWarmth ?? 0.72);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    this.lastMarchMs = performance.now() - t0;
  }

  draw(frame) {
    if (!this._visible || !this.params) return;
    const gl = this.gl;
    const prog = this.composite.use();
    gl.disable(gl.DEPTH_TEST);
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    prog.tex("uTex", this.target.texture);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    gl.enable(gl.DEPTH_TEST);
  }

  invalidate() {
    this._cacheKey = "";
  }

  dispose() {
    const gl = this.gl;
    for (const v of this.volumes.values()) this._freeVolume(v);
    for (const t of this.cmaps.values()) gl.deleteTexture(t);
    this._dropMasks();
    if (this._blank) gl.deleteTexture(this._blank);
    this.target.dispose();
    this.program.dispose();
    this.composite.dispose();
  }
}

/** A lasso clip's identity (null: no clip), so an unchanged one keeps its texture. */
function clipKey(mask) {
  if (!mask || !Array.isArray(mask.polygon) || mask.polygon.length < 3 || !mask.vp || mask.vp.length !== 16) return null;
  return JSON.stringify([mask.polygon, Array.from(mask.vp), clipOutside(mask)]);
}

function clipOutside(mask) {
  return Math.min(Math.max(Number(mask.outside) || 0, 0), 1);
}

function dist2(a, b) {
  const x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
  return x * x + y * y + z * z;
}

/**
 * Conservative block maxima (≤64 blocks per axis), dilated by one voxel so
 * trilinear tails at block borders are always covered. Used only when the
 * bundle did not ship a precomputed grid.
 */
export function computeOccupancy(data, nx, ny, nz) {
  const total = nx * ny * nz;
  if (total < (1 << 18)) return null;
  const gx = Math.min(64, nx), gy = Math.min(64, ny), gz = Math.min(64, nz);
  const blocks = new Uint8Array(gx * gy * gz);
  const isFloat = data instanceof Float32Array;
  const bx = new Int32Array(nx), by = new Int32Array(ny), bz = new Int32Array(nz);
  for (let x = 0; x < nx; x++) bx[x] = Math.min(gx - 1, Math.floor(((x + 0.5) * gx) / nx));
  for (let y = 0; y < ny; y++) by[y] = Math.min(gy - 1, Math.floor(((y + 0.5) * gy) / ny));
  for (let z = 0; z < nz; z++) bz[z] = Math.min(gz - 1, Math.floor(((z + 0.5) * gz) / nz));
  let any = false;
  let i = 0;
  for (let z = 0; z < nz; z++) {
    const z0 = bz[Math.max(z - 1, 0)], z1 = bz[Math.min(z + 1, nz - 1)];
    for (let y = 0; y < ny; y++) {
      const y0 = by[Math.max(y - 1, 0)], y1 = by[Math.min(y + 1, ny - 1)];
      for (let x = 0; x < nx; x++, i++) {
        // Float voxels (0–1) round up to the uint8 scale: never under-cover.
        const v = isFloat ? (data[i] > 0 ? Math.min(255, Math.ceil(data[i] * 255)) : 0) : data[i];
        if (!v) continue;
        any = true;
        const x0 = bx[Math.max(x - 1, 0)], x1 = bx[Math.min(x + 1, nx - 1)];
        for (let zz = z0; zz <= z1; zz++) {
          for (let yy = y0; yy <= y1; yy++) {
            const row = (zz * gy + yy) * gx;
            for (let xx = x0; xx <= x1; xx++) if (blocks[row + xx] < v) blocks[row + xx] = v;
          }
        }
      }
    }
  }
  return any ? { data: blocks, gx, gy, gz } : null;
}
