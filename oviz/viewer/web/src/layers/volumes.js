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
out vec4 outColor;

vec3 rotZ(vec3 p, vec2 cs) {
  return vec3(cs.x * p.x - cs.y * p.y, cs.y * p.x + cs.x * p.y, p.z);
}

float asinhF(float v) { return log(v + sqrt(v * v + 1.0)); }

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
    v = stretch(min(v, 0.999));
    vec4 c = texture(uCmap, vec2(v, 0.5));
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

  /** Upload a volume. `data` is Uint8Array (z, y, x); `occ` optional block maxima. */
  addVolume(spec, data, occ) {
    const gl = this.gl;
    // Replacing a volume frees the old textures (never leak GPU memory).
    const old = this.volumes.get(spec.key);
    if (old) {
      gl.deleteTexture(old.tex);
      if (old.occTex) gl.deleteTexture(old.occTex);
    }
    const [nx, ny, nz] = spec.dims;
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_3D, tex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage3D(gl.TEXTURE_3D, 0, gl.R8, nx, ny, nz, 0, gl.RED, gl.UNSIGNED_BYTE, data);
    const filter = spec.interpolation === false ? gl.NEAREST : gl.LINEAR;
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, filter);
    for (const w of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T, gl.TEXTURE_WRAP_R]) gl.texParameteri(gl.TEXTURE_3D, w, gl.CLAMP_TO_EDGE);
    let occTex = null, occDims = [1, 1, 1];
    const grid = occ || computeOccupancy(data, nx, ny, nz);
    if (grid) {
      occTex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_3D, occTex);
      gl.texImage3D(gl.TEXTURE_3D, 0, gl.R8, grid.gx, grid.gy, grid.gz, 0, gl.RED, gl.UNSIGNED_BYTE, grid.data);
      gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      for (const w of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T, gl.TEXTURE_WRAP_R]) gl.texParameteri(gl.TEXTURE_3D, w, gl.CLAMP_TO_EDGE);
      occDims = [grid.gx, grid.gy, grid.gz];
    }
    gl.bindTexture(gl.TEXTURE_3D, null);
    const [lo, hi] = spec.bounds;
    const vol = {
      key: spec.key, spec, tex, occTex, occDims,
      center: [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2],
      boxMin: lo, boxMax: hi,
    };
    this.volumes.set(spec.key, vol);
    this._cacheKey = "";
    return vol;
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
    const key = [frame.camera.version, w, h, p.frame.toFixed(5), p.volumeSignature, frame.width, frame.height].join("|");
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
    for (const v of this.volumes.values()) {
      gl.deleteTexture(v.tex);
      if (v.occTex) gl.deleteTexture(v.occTex);
    }
    for (const t of this.cmaps.values()) gl.deleteTexture(t);
    this.target.dispose();
    this.program.dispose();
    this.composite.dispose();
  }
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
        const v = data[i];
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
