// Instanced point sprites for every point trace.
//
// One draw call per trace. Each instance fetches its interpolated position
// from the trace's frame texture, applies the legacy Oviz sizing model
// (world-space size, size-by-n_stars, birth fade) and renders the stellar
// halo + core PSF in a single premultiplied pass: halo·a_h + core·a_c +
// dst·(1 − a_h) — identical to the legacy "normal halo, additive core" pair.

import { createProgram, createBuffer, createVAO, createTexture2D } from "../engine/gl.js";
import { FRAME_GLSL, packFrameTexture, packScalarTexture, frameOffset } from "../engine/frames.js";
import { parseColor } from "../core/color.js";

const VERT = `
${FRAME_GLSL}
layout(location = 0) in vec2 aCorner;
layout(location = 1) in float aIndex;
layout(location = 2) in vec3 aColor;
layout(location = 3) in vec4 aStatic;   // opacity, scalar, ageNow, starsFactor
layout(location = 4) in float aSizeMax;
layout(location = 5) in float aSymbol;
layout(location = 6) in float aState;   // 0 normal, 1 dimmed, 2 selected, 3 hidden

uniform mat4 uViewProj;
uniform mat4 uView;
uniform float uProjY;
uniform vec2 uViewport;
uniform highp sampler2D uFrameTex;   // xyz + size
uniform highp sampler2D uAuxTex;     // opacity, scalar (when varying)
uniform int uCount;
uniform int uFrames;
uniform int uAuxFrames;
uniform float uFrame;
uniform float uTime;
uniform vec3 uOffset;
uniform float uPointScale;
uniform float uSizeScale;      // trace size × global size
uniform float uGlobalSize;
uniform float uOpacityScale;   // trace opacity ratio × visibility × global
uniform float uGlow;
uniform vec3 uFade;            // fadeTime, fadeInOut, fadeByOpacity
uniform float uStarsExp;
uniform int uHover;
uniform float uHoverBoost;
uniform float uMinPx;
uniform float uDimOpacity;
uniform int uHasAux;
uniform int uHasAge;
uniform float uPickMinPx;
uniform float uMemberFade;

out vec2 vUv;
out vec3 vColor;
out float vScalar;
out float vGlowA;
out float vCoreA;
out float vCoreRatio;
out float vSymbol;
flat out int vIndex;

float smooth01(float x) { x = clamp(x, 0.0, 1.0); return x * x * (3.0 - 2.0 * x); }

float birthFade(float t, float ageNow, float fade, bool inOut) {
  float birth = -ageNow;
  if (fade <= 1e-9) {
    if (t < birth) return 0.0;
    if (!inOut) return 1.0;
    return abs(t - birth) <= 1e-9 ? 1.0 : 0.0;
  }
  if (t < birth - fade) return 0.0;
  if (t <= birth) return smooth01((t - (birth - fade)) / fade);
  if (!inOut) return 1.0;
  if (t <= birth + fade) return 1.0 - smooth01((t - birth) / fade);
  return 0.0;
}

void main() {
  int index = int(aIndex + 0.5);
  vIndex = index;
  float scalarFrame;
  vec4 p = framePosition(uFrameTex, index, uCount, uFrames, uFrame, scalarFrame);
  float frameSize = scalarFrame;
  float opacity = aStatic.x;
  float scalar = aStatic.y;
  if (uHasAux == 1) {
    vec4 aux = frameScalars(uAuxTex, index, uCount, uAuxFrames, uFrame);
    opacity = aux.x;
    scalar = aux.y;
  }
  float fade = 1.0;
  float size = frameSize;
  bool hasAge = uHasAge == 1 && aStatic.z == aStatic.z && aStatic.z > -1e29;
  if (hasAge) {
    fade = birthFade(uTime, aStatic.z, uFade.x, uFade.y > 0.5);
    if (uFade.z > 0.5) {
      size = aSizeMax > 0.0 ? aSizeMax : frameSize;
      opacity *= fade;
    } else {
      size = frameSize * fade;
    }
  }
  float stars = pow(max(aStatic.w, 1e-6), uStarsExp);
  float eff = opacity * uOpacityScale * p.w;
  // State bits: 1 dimmed (filter), 2 hidden (filter), 4 replaced by member
  // stars in Sky view.
  int st = int(aState + 0.5);
  if ((st & 2) != 0) eff = 0.0;
  else if ((st & 1) != 0) eff *= uDimOpacity;
  if ((st & 4) != 0) eff *= 1.0 - uMemberFade;
  vec4 center = uViewProj * vec4(p.xyz + uOffset, 1.0);
  if (size <= 0.0 || eff <= 0.001 || p.w <= 0.0 || center.w <= 0.0) {
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
    return;
  }
  float scale = max(size * uSizeScale * stars * uPointScale, uPointScale * 0.5 * uGlobalSize);
  bool hovered = index == uHover;
  float quadWorld;
  if (uGlow > 0.02) {
    float glowSize = scale * (3.15 + 0.70 * uGlow) * (hovered ? 1.08 * uHoverBoost : 1.0);
    float coreSize = scale * (2.65 + 0.18 * uGlow) * (hovered ? 1.22 * uHoverBoost : 1.0);
    quadWorld = max(glowSize, coreSize);
    vCoreRatio = quadWorld / max(coreSize, 1e-6);
    vGlowA = clamp(eff * (0.34 + 0.18 * uGlow), 0.0, 0.78);
    vCoreA = clamp(eff * (1.0 + 0.24 * uGlow), 0.0, 1.0);
  } else {
    quadWorld = scale * 0.36 * (hovered ? 1.45 * uHoverBoost : 1.0);
    vCoreRatio = 1.0;
    vGlowA = clamp(eff, 0.0, 1.0);
    vCoreA = 0.0;
  }
  float viewZ = -(uView * vec4(p.xyz + uOffset, 1.0)).z;
  float px = quadWorld * 0.5 * uViewport.y * uProjY / max(viewZ, 1e-6);
  #ifdef PICK
  px = max(px * (uGlow > 0.02 ? 0.45 : 1.0), uPickMinPx);
  #else
  // Keep sub-pixel stars stable: clamp the footprint and conserve flux.
  if (px < uMinPx) {
    float k = (px / uMinPx);
    vGlowA *= k * k;
    vCoreA *= k * k;
    px = uMinPx;
  }
  #endif
  vec2 ndc = aCorner * px / uViewport;
  gl_Position = center + vec4(ndc * center.w, 0.0, 0.0);
  vUv = aCorner;
  vColor = aColor;
  vScalar = scalar;
  vSymbol = aSymbol;
}
`;

const FRAG = `
in vec2 vUv;
in vec3 vColor;
in float vScalar;
in float vGlowA;
in float vCoreA;
in float vCoreRatio;
in float vSymbol;
flat in int vIndex;
uniform sampler2D uHaloTex;
uniform sampler2D uCoreTex;
uniform sampler2D uCmap;
uniform int uColorMode;       // 0 per-object/trace colour, 1 colormap
uniform vec3 uTraceColor;
uniform int uUseObjectColor;
uniform vec2 uCRange;
uniform float uGlow;
uniform int uPickTrace;
out vec4 outColor;

float sdSymbol(vec2 p, float kind) {
  // Returns signed distance (negative inside) for unit-radius symbols.
  if (kind < 0.5) return length(p) - 0.82;                          // circle
  if (kind < 1.5) { vec2 d = abs(p) - vec2(0.7); return max(d.x, d.y); } // square
  if (kind < 2.5) return (abs(p.x) + abs(p.y)) - 0.9;                // diamond
  if (kind < 3.5) { vec2 a = abs(p); return min(max(a.x - 0.85, a.y - 0.24), max(a.y - 0.85, a.x - 0.24)); } // cross
  if (kind < 4.5) { vec2 q = vec2(p.x + p.y, p.x - p.y) * 0.7071; vec2 a = abs(q); return min(max(a.x - 0.9, a.y - 0.2), max(a.y - 0.9, a.x - 0.2)); } // x
  if (kind < 5.5) { vec2 q = vec2(abs(p.x), -p.y + 0.25); return max(q.x * 0.866 + q.y * 0.5 - 0.5, -q.y - 0.55); } // triangle
  if (kind < 6.5) { float a = atan(p.y, p.x); float r = length(p); float k = 0.55 + 0.35 * cos(5.0 * a); return r - k; } // star
  return abs(length(p) - 0.7) - 0.14;                                  // ring
}

void main() {
  float r = length(vUv);
  if (r > 1.0) discard;
  #ifdef PICK
  if (uGlow <= 0.02) {
    if (sdSymbol(vUv, vSymbol) > 0.1) discard;
  }
  int id = vIndex + 1;
  outColor = vec4(float(uPickTrace) / 255.0, float((id >> 16) & 255) / 255.0, float((id >> 8) & 255) / 255.0, float(id & 255) / 255.0);
  return;
  #else
  vec3 color = uUseObjectColor == 1 ? vColor : uTraceColor;
  if (uColorMode == 1) {
    float u = clamp((vScalar - uCRange.x) / max(uCRange.y - uCRange.x, 1e-9), 0.0, 1.0);
    color = texture(uCmap, vec2(u, 0.5)).rgb;
  }
  if (uGlow > 0.02) {
    vec2 uv = vUv * 0.5 + 0.5;
    float halo = texture(uHaloTex, uv).r;
    vec2 cuv = vUv * vCoreRatio * 0.5 + 0.5;
    float core = (vCoreRatio * r <= 1.0) ? texture(uCoreTex, cuv).r : 0.0;
    float ah = halo * vGlowA;
    float ac = core * vCoreA;
    if (ah + ac < 0.0015) discard;
    outColor = vec4(color * ah + vec3(ac), ah);
  } else {
    float d = sdSymbol(vUv, vSymbol);
    float w = fwidth(d) * 0.75 + 1e-4;
    float a = (1.0 - smoothstep(-w, w, d)) * vGlowA;
    if (a < 0.004) discard;
    outColor = vec4(color * a, a);
  }
  #endif
}
`;

/** Radial profiles of the legacy Oviz star textures (r = 1 at the edge). */
const PSF = {
  halo: (r) => 0.82 * Math.exp(-0.5 * (r / 0.024) ** 2)
    + 0.30 * (1 + (r / 0.060) ** 2) ** -2.25
    + 0.050 * (1 + (r / 0.20) ** 2) ** -2.5
    + 0.055 * Math.exp(-0.5 * ((r - 0.22) / 0.11) ** 2),
  core: (r) => Math.exp(-0.5 * (r / 0.020) ** 2) + 0.28 * (1 + (r / 0.052) ** 2) ** -2.65,
  "member-halo": (r) => 0.72 * Math.exp(-0.5 * (r / 0.075) ** 2)
    + 0.24 * (1 + (r / 0.14) ** 2) ** -2.4
    + 0.035 * (1 + (r / 0.32) ** 2) ** -3.0,
  "member-core": (r) => Math.exp(-0.5 * (r / 0.055) ** 2) + 0.16 * (1 + (r / 0.12) ** 2) ** -2.8,
};

function psfTexture(gl, kind, size = 512) {
  // Tabulate the radial profile once, then fill the image by lookup; the
  // texture matches a direct evaluation to < 1/255.
  const LUT = 4096;
  const table = new Float32Array(LUT + 2);
  const f = PSF[kind];
  for (let i = 0; i <= LUT + 1; i++) {
    const r = (i / LUT) * 1.5;
    const t = Math.min(Math.max((r - 0.78) / 0.22, 0), 1);
    const taper = 1 - t * t * (3 - 2 * t);
    table[i] = Math.min(Math.max(f(r) * taper, 0), 1) * 255;
  }
  const data = new Uint8Array(size * size);
  const c = (size - 1) * 0.5;
  const invR = 2 / size;
  const k = LUT / 1.5;
  for (let y = 0; y < size; y++) {
    const dy = (y - c) * invR;
    const dy2 = dy * dy;
    const row = y * size;
    for (let x = 0; x < size; x++) {
      const dx = (x - c) * invR;
      const fi = Math.sqrt(dx * dx + dy2) * k;
      const i = fi | 0;
      const v = i >= LUT ? 0 : table[i] + (table[i + 1] - table[i]) * (fi - i);
      data[row + x] = v + 0.5;
    }
  }
  return createTexture2D(gl, {
    width: size, height: size, internalFormat: gl.R8, format: gl.RED, type: gl.UNSIGNED_BYTE,
    data, mipmaps: true,
  });
}

let sharedTextures = null;
/** Forget cached textures (they die with a lost WebGL context). */
export function resetStarTextures() {
  sharedTextures = null;
}

export function starTextures(gl) {
  if (!sharedTextures || sharedTextures.gl !== gl) {
    sharedTextures = {
      gl,
      halo: psfTexture(gl, "halo"),
      core: psfTexture(gl, "core"),
      memberHalo: psfTexture(gl, "member-halo", 256),
      memberCore: psfTexture(gl, "member-core", 256),
    };
  }
  return sharedTextures;
}

const QUAD = new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]);

/** One trace's GPU resources. */
class PointBatch {
  constructor(gl, trace, data, frames) {
    this.gl = gl;
    this.trace = trace;
    this.key = trace.key;
    const pts = trace.points;
    const N = pts.count;
    this.count = N;
    this.frames = frames;
    this.positions = data.position;
    this.posFrames = pts.position.frames || 1;
    this.offsets = data.offset || null;
    const sizeFrames = data.size ? (pts.size.frames || 1) : 1;
    const texFrames = Math.max(this.posFrames, sizeFrames);
    this.texFrames = texFrames;
    this.frameTex = packFrameTexture(gl, {
      positions: data.position,
      posFrames: this.posFrames,
      count: N,
      scalar: data.size,
      scalarFrames: sizeFrames,
      scalarDefault: pts.size.value ?? trace.pointSize ?? 1,
      frames: texFrames,
    });
    const opFrames = data.opacity ? (pts.opacity.frames || 1) : 1;
    const scFrames = data.scalar ? (pts.colorScalar.frames || 1) : 1;
    this.hasAux = opFrames > 1 || scFrames > 1;
    if (this.hasAux) {
      const auxFrames = Math.max(opFrames, scFrames);
      this.auxFrames = auxFrames;
      this.auxTex = packScalarTexture(gl, {
        count: N,
        frames: auxFrames,
        channels: [
          { values: data.opacity, frames: opFrames, fallback: pts.opacity.value ?? trace.opacity ?? 1 },
          { values: data.scalar, frames: scFrames, fallback: pts.colorScalar?.value ?? 0 },
        ],
      });
    }
    // Instance attributes.
    const index = new Float32Array(N);
    for (let i = 0; i < N; i++) index[i] = i;
    const color = new Uint8Array(N * 3);
    if (data.rgb) color.set(data.rgb);
    else {
      const c = parseColor(trace.color);
      for (let i = 0; i < N; i++) {
        color[i * 3] = c[0] * 255; color[i * 3 + 1] = c[1] * 255; color[i * 3 + 2] = c[2] * 255;
      }
    }
    this.objectColors = color;
    const stat = new Float32Array(N * 4);
    const opStatic = !data.opacity || opFrames > 1 ? null : data.opacity;
    const scStatic = !data.scalar || scFrames > 1 ? null : data.scalar;
    for (let i = 0; i < N; i++) {
      stat[i * 4] = opStatic ? opStatic[i] : (pts.opacity.value ?? trace.opacity ?? 1);
      stat[i * 4 + 1] = scStatic ? scStatic[i] : (pts.colorScalar?.value ?? 0);
      const age = data.ageNow ? data.ageNow[i] : NaN;
      stat[i * 4 + 2] = age === age ? age : -1e30;
      stat[i * 4 + 3] = data.starsFactor ? data.starsFactor[i] : 1;
    }
    this.hasAge = !!data.ageNow;
    const sizeMax = new Float32Array(N);
    if (data.sizeMax) sizeMax.set(data.sizeMax);
    else if (data.size && sizeFrames === 1) sizeMax.set(data.size);
    else sizeMax.fill(pts.size.value ?? trace.pointSize ?? 1);
    const symbol = new Float32Array(N);
    if (data.symbol) for (let i = 0; i < N; i++) symbol[i] = data.symbol[i];
    else symbol.fill(Math.max(0, (pts.symbols || []).indexOf(trace.symbol)));
    this.state = new Uint8Array(N);
    this.buffers = {
      quad: createBuffer(gl, QUAD),
      index: createBuffer(gl, index),
      color: createBuffer(gl, color),
      stat: createBuffer(gl, stat),
      sizeMax: createBuffer(gl, sizeMax),
      symbol: createBuffer(gl, symbol),
      state: createBuffer(gl, this.state, gl.DYNAMIC_DRAW),
    };
    const b = this.buffers;
    this.vao = createVAO(gl, [
      { loc: 0, buffer: b.quad, size: 2 },
      { loc: 1, buffer: b.index, size: 1, divisor: 1 },
      { loc: 2, buffer: b.color, size: 3, type: gl.UNSIGNED_BYTE, normalized: true, divisor: 1 },
      { loc: 3, buffer: b.stat, size: 4, divisor: 1 },
      { loc: 4, buffer: b.sizeMax, size: 1, divisor: 1 },
      { loc: 5, buffer: b.symbol, size: 1, divisor: 1 },
      { loc: 6, buffer: b.state, size: 1, type: gl.UNSIGNED_BYTE, divisor: 1 },
    ]);
  }

  /** Replace the bits in `mask` with `values` (Uint8Array, same length). */
  setStateBits(mask, values) {
    const gl = this.gl;
    const s = this.state;
    for (let i = 0; i < s.length; i++) s[i] = (s[i] & ~mask) | (values ? values[i] & mask : 0);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.state);
    gl.bufferSubData(gl.ARRAY_BUFFER, 0, s);
  }

  dispose() {
    const gl = this.gl;
    gl.deleteVertexArray(this.vao);
    for (const buf of Object.values(this.buffers)) gl.deleteBuffer(buf);
    gl.deleteTexture(this.frameTex.texture);
    if (this.auxTex) gl.deleteTexture(this.auxTex.texture);
  }
}

export class PointsLayer {
  constructor(gl, { frames }) {
    this.gl = gl;
    this.order = 30;
    this.frames = frames;
    this.batches = [];
    this.byKey = new Map();
    this.program = createProgram(gl, VERT, FRAG, { label: "points" });
    this.pickProgram = createProgram(gl, VERT, FRAG, { label: "points-pick", defines: { PICK: 1 } });
    this.textures = starTextures(gl);
    this.cmapTextures = new Map();
    this.fallbackCmap = createTexture2D(gl, {
      width: 2, height: 1, internalFormat: gl.RGBA8, format: gl.RGBA, type: gl.UNSIGNED_BYTE,
      data: new Uint8Array([40, 40, 200, 255, 250, 220, 40, 255]),
    });
    // Resolved each frame by the viewer.
    this.params = null;
  }

  addTrace(trace, data, pickId) {
    const batch = new PointBatch(this.gl, trace, data, this.frames);
    batch.pickId = pickId;
    this.batches.push(batch);
    this.byKey.set(trace.key, batch);
    return batch;
  }

  setColormap(name, lutBytes, width) {
    const gl = this.gl;
    if (this.cmapTextures.has(name)) return;
    this.cmapTextures.set(name, createTexture2D(gl, {
      width, height: 1, internalFormat: gl.RGBA8, format: gl.RGBA, type: gl.UNSIGNED_BYTE, data: lutBytes,
    }));
  }

  _uniforms(prog, frame, batch, style, global) {
    const cam = frame.camera;
    prog.m4("uViewProj", cam.viewProj).m4("uView", cam.view).f("uProjY", cam.proj[5]);
    prog.v2("uViewport", frame.width, frame.height);
    prog.tex("uFrameTex", batch.frameTex.texture);
    prog.i("uCount", batch.count).i("uFrames", batch.texFrames).f("uFrame", global.frame).f("uTime", global.time);
    if (batch.hasAux) {
      prog.tex("uAuxTex", batch.auxTex.texture).i("uAuxFrames", batch.auxFrames).i("uHasAux", 1);
    } else {
      prog.tex("uAuxTex", batch.frameTex.texture).i("uAuxFrames", 1).i("uHasAux", 0);
    }
    const off = frameOffset(batch.offsets, global.frame);
    prog.v3("uOffset", off[0], off[1], off[2]);
    prog.f("uPointScale", global.pointScale);
    prog.f("uSizeScale", style.sizeScale * global.pointSize);
    prog.f("uGlobalSize", global.pointSize);
    prog.f("uOpacityScale", style.opacityScale * global.pointOpacity);
    prog.f("uGlow", global.glow);
    prog.v3("uFade", global.fadeTime, global.fadeInOut ? 1 : 0, global.fadeByOpacity ? 1 : 0);
    prog.f("uStarsExp", style.starsExp);
    prog.i("uHover", style.hover ?? -1);
    prog.f("uHoverBoost", 1.0);
    prog.f("uMinPx", 1.25 * frame.dpr);
    prog.f("uPickMinPx", 9 * frame.dpr);
    prog.f("uDimOpacity", style.dimOpacity ?? 0.16);
    prog.i("uHasAge", batch.hasAge ? 1 : 0);
    prog.f("uMemberFade", global.memberFade ?? 0);
  }

  draw(frame) {
    const p = this.params;
    if (!p) return;
    const gl = this.gl;
    const prog = this.program.use();
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.enable(gl.DEPTH_TEST);
    gl.depthMask(false);
    prog.tex("uHaloTex", this.textures.halo);
    prog.tex("uCoreTex", this.textures.core);
    for (const batch of this.batches) {
      const style = p.styles.get(batch.key);
      if (!style || !style.visible || style.presence <= 0.001) continue;
      this._uniforms(prog, frame, batch, style, p);
      const colorMode = style.colorMode === "by_value" && batch.trace.colorBy ? 1 : 0;
      prog.i("uColorMode", colorMode);
      if (colorMode) {
        const cm = this.cmapTextures.get(style.colormap) || this.fallbackCmap;
        prog.tex("uCmap", cm);
        prog.v2("uCRange", style.cmin, style.cmax);
      } else {
        prog.tex("uCmap", this.fallbackCmap);
      }
      const c = style.colorOverride ? parseColor(style.colorOverride) : parseColor(batch.trace.color);
      prog.v3("uTraceColor", c[0], c[1], c[2]);
      prog.i("uUseObjectColor", style.colorOverride ? 0 : 1);
      gl.bindVertexArray(batch.vao);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, batch.count);
    }
    gl.bindVertexArray(null);
  }

  drawPick(frame) {
    const p = this.params;
    if (!p) return;
    const gl = this.gl;
    const prog = this.pickProgram.use();
    gl.disable(gl.BLEND);
    for (const batch of this.batches) {
      const style = p.styles.get(batch.key);
      if (!style || !style.visible || style.presence <= 0.001 || style.pickable === false) continue;
      this._uniforms(prog, frame, batch, style, p);
      prog.i("uPickTrace", batch.pickId);
      gl.bindVertexArray(batch.vao);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, batch.count);
    }
    gl.bindVertexArray(null);
  }

  dispose() {
    for (const b of this.batches) b.dispose();
    this.batches = [];
    this.program.dispose();
    this.pickProgram.dispose();
  }
}
