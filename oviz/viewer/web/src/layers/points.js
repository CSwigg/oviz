// Instanced point sprites for every point trace.
//
// One draw call per trace. Each instance fetches its interpolated position
// from the trace's frame texture, applies the legacy Oviz sizing model
// (world-space size, size-by-n_stars, birth fade) and renders the stellar
// halo + core PSF in a single premultiplied pass: halo·a_h + core·a_c +
// dst·(1 − a_h) — identical to the legacy "normal halo, additive core" pair.
//
// Optional motion trails: a tapered screen-space ribbon per point, sampled
// back along its own orbit from the same frame texture, so trails follow the
// exact interpolated path at no CPU cost.

import { createProgram, createBuffer, createVAO, createTexture2D } from "../engine/gl.js";
import { FRAME_GLSL, packFrameTexture, packScalarTexture, frameOffset } from "../engine/frames.js";
import { parseColor } from "../core/color.js";

// A click picks the point itself, not its glow: the pick disc is the core's,
// out to where the core profile has fallen to ~3% of its peak (in the core
// texture's radius units), at least PICK_MIN_PX across. pick()'s search
// radius then forgives near misses on small points.
export const POINT_PICK_CORE = 0.06;
export const PICK_MIN_PX = 3;

const VERT = `
${FRAME_GLSL}
layout(location = 0) in vec2 aCorner;
layout(location = 1) in float aIndex;
layout(location = 2) in vec3 aColor;
layout(location = 3) in vec4 aStatic;   // opacity, scalar, ageNow, starsFactor
layout(location = 4) in float aSizeMax;
layout(location = 5) in float aSymbol;
layout(location = 6) in float aState;   // state bits (below)
layout(location = 7) in vec2 aLasso;    // lasso level before / after a State change

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
uniform vec2 uLassoMix;

// The lasso: y is the object's level in the selection (1 shown, ½ dimmed,
// 0 hidden), x the weight it was drawn with when a State change began. The
// change blends one into the other: what leaves the selection fades out
// early in the flight, what joins fades in as it lands (uLassoMix:
// fade-out, fade-in progress; 1, 1 at rest, where only the level counts).
// Starting from the weight drawn, not the last level, lets an interrupted
// change carry on from where it was.
float lassoWeight(vec2 levels) {
  float b = levels.y > 0.75 ? 1.0 : levels.y > 0.25 ? uDimOpacity : 0.0;
  return mix(levels.x, b, b < levels.x ? uLassoMix.x : uLassoMix.y);
}

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
  // State bits: 1 dimmed / 2 hidden (distribution filter), 4 replaced by
  // member stars in Sky view, 32 dimmed (birth tree highlight). The lasso
  // (bits 8 / 16, for the CPU's readers) draws through aLasso, which can
  // blend; dimmed by both counts once.
  int st = int(aState + 0.5);
  float bitW = (st & 2) != 0 ? 0.0 : (st & 33) != 0 ? uDimOpacity : 1.0;
  eff *= min(bitW, lassoWeight(aLasso));
  if ((st & 4) != 0) eff *= 1.0 - uMemberFade;
  vec4 center = uViewProj * vec4(p.xyz + uOffset, 1.0);
  if (size <= 0.0 || eff <= 0.001 || p.w <= 0.0 || center.w <= 0.0) {
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
    return;
  }
  float scale = max(size * uSizeScale * stars * uPointScale, uPointScale * 0.5 * uGlobalSize);
  bool hovered = index == uHover;
  float quadWorld;
  float coreRatio = 1.0;  // quad radius over the core's
  if (uGlow > 0.02) {
    float glowSize = scale * (3.15 + 0.70 * uGlow) * (hovered ? 1.08 * uHoverBoost : 1.0);
    float coreSize = scale * (2.65 + 0.18 * uGlow) * (hovered ? 1.22 * uHoverBoost : 1.0);
    quadWorld = max(glowSize, coreSize);
    coreRatio = quadWorld / max(coreSize, 1e-6);
    vCoreRatio = coreRatio;
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
  // The point, not its glow (see POINT_PICK_CORE); a symbol is its own shape.
  px = max(uGlow > 0.02 ? px * ${POINT_PICK_CORE} / coreRatio : px, uPickMinPx);
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

const TRAIL_SEGMENTS = 24;

const TRAIL_VERT = `
${FRAME_GLSL}
layout(location = 0) in vec2 aTrail;    // s (0 head → 1 tail), side (±1)
layout(location = 1) in float aIndex;
layout(location = 2) in vec3 aColor;
layout(location = 3) in vec4 aStatic;   // opacity, scalar, ageNow, starsFactor
layout(location = 6) in float aState;
layout(location = 7) in vec2 aLasso;

uniform mat4 uViewProj;
uniform vec2 uViewport;
uniform highp sampler2D uFrameTex;
uniform int uCount;
uniform int uFrames;
uniform float uFrame;
uniform float uSpan;          // trail length in frames, signed toward the past
uniform float uTime;
uniform float uMyrPerFrame;   // signed
uniform vec3 uOffset;
uniform float uOpacityScale;
uniform float uWidthPx;
uniform float uDimOpacity;
uniform int uHasAge;
uniform float uMemberFade;
uniform highp sampler2D uAuxTex;     // per-frame opacity, colour value
uniform int uAuxFrames;
uniform int uHasAux;
uniform vec2 uLassoMix;

// The lasso: y is the object's level in the selection (1 shown, ½ dimmed,
// 0 hidden), x the weight it was drawn with when a State change began. The
// change blends one into the other: what leaves the selection fades out
// early in the flight, what joins fades in as it lands (uLassoMix:
// fade-out, fade-in progress; 1, 1 at rest, where only the level counts).
// Starting from the weight drawn, not the last level, lets an interrupted
// change carry on from where it was.
float lassoWeight(vec2 levels) {
  float b = levels.y > 0.75 ? 1.0 : levels.y > 0.25 ? uDimOpacity : 0.0;
  return mix(levels.x, b, b < levels.x ? uLassoMix.x : uLassoMix.y);
}

out vec3 vColor;
out float vScalar;
out float vAlpha;
out float vSide;

vec4 clipAt(int index, float f, out float presence) {
  float sc;
  vec4 p = framePosition(uFrameTex, index, uCount, uFrames, f, sc);
  presence = p.w;
  return uViewProj * vec4(p.xyz + uOffset, 1.0);
}

void main() {
  int index = int(aIndex + 0.5);
  float s = aTrail.x;
  float ds = 1.0 / float(${TRAIL_SEGMENTS});
  float last = float(uFrames - 1);
  float pres, pa, pb;
  float fs = clamp(uFrame - uSpan * s, 0.0, last);
  vec4 c = clipAt(index, fs, pres);
  vec4 ca = clipAt(index, clamp(uFrame - uSpan * max(s - ds, 0.0), 0.0, last), pa);
  vec4 cb = clipAt(index, clamp(uFrame - uSpan * min(s + ds, 1.0), 0.0, last), pb);
  // Opacity and colour value where the object was at this point of the trail.
  float opacity = aStatic.x;
  float scalar = aStatic.y;
  if (uHasAux == 1) {
    vec4 aux = frameScalars(uAuxTex, index, uCount, uAuxFrames, fs);
    opacity = aux.x;
    scalar = aux.y;
  }
  float eff = opacity * uOpacityScale * pres;
  int st = int(aState + 0.5);
  float bitW = (st & 2) != 0 ? 0.0 : (st & 33) != 0 ? uDimOpacity : 1.0;
  eff *= min(bitW, lassoWeight(aLasso));
  if ((st & 4) != 0) eff *= 1.0 - uMemberFade;
  if (uHasAge == 1 && aStatic.z == aStatic.z && aStatic.z > -1e29) {
    // No trail before the object was born.
    float birth = -aStatic.z;
    eff *= smoothstep(birth - 0.25, birth + 0.25, uTime - uSpan * s * uMyrPerFrame);
  }
  if (eff <= 0.002 || c.w <= 0.0 || ca.w <= 0.0 || cb.w <= 0.0) {
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
    vAlpha = 0.0;
    return;
  }
  vec2 hv = uViewport * 0.5;
  vec2 dir = (ca.xy / ca.w - cb.xy / cb.w) * hv;
  float len = length(dir);
  dir = len > 1e-4 ? dir / len : vec2(1.0, 0.0);
  vec2 nrm = vec2(-dir.y, dir.x);
  float w = uWidthPx * mix(1.0, 0.3, s);
  gl_Position = c + vec4(nrm * aTrail.y * w * 0.5 / hv * c.w, 0.0, 0.0);
  float tail = 1.0 - s;
  vAlpha = eff * tail * tail;
  vColor = aColor;
  vScalar = scalar;
  vSide = aTrail.y;
}
`;

const TRAIL_FRAG = `
in vec3 vColor;
in float vScalar;
in float vAlpha;
in float vSide;
uniform sampler2D uCmap;
uniform int uColorMode;
uniform vec3 uTraceColor;
uniform int uUseObjectColor;
uniform vec2 uCRange;
uniform float uStrength;
out vec4 outColor;

void main() {
  vec3 color = uUseObjectColor == 1 ? vColor : uTraceColor;
  if (uColorMode == 1) {
    float u = clamp((vScalar - uCRange.x) / max(uCRange.y - uCRange.x, 1e-9), 0.0, 1.0);
    color = texture(uCmap, vec2(u, 0.5)).rgb;
  }
  float a = vAlpha * uStrength * (1.0 - smoothstep(0.35, 1.0, abs(vSide)));
  if (a < 0.003) discard;
  // Mostly additive: trails brighten the sky rather than covering it.
  outColor = vec4(color * a, a * 0.45);
}
`;

function trailTemplate() {
  const out = new Float32Array((TRAIL_SEGMENTS + 1) * 4);
  for (let i = 0; i <= TRAIL_SEGMENTS; i++) {
    const s = i / TRAIL_SEGMENTS;
    out.set([s, -1, s, 1], i * 4);
  }
  return out;
}

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

/** A star sprite's radial profile (halo, core, member-halo, member-core) at r in sprite radii. */
export function starProfile(kind, r) {
  return PSF[kind](r);
}

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
    // Per object, the weight a lasso fade starts from and the selection's level (see aLasso).
    this.lasso = new Uint8Array(N * 2).fill(255);
    this.buffers = {
      quad: createBuffer(gl, QUAD),
      index: createBuffer(gl, index),
      color: createBuffer(gl, color),
      stat: createBuffer(gl, stat),
      sizeMax: createBuffer(gl, sizeMax),
      symbol: createBuffer(gl, symbol),
      state: createBuffer(gl, this.state, gl.DYNAMIC_DRAW),
      lasso: createBuffer(gl, this.lasso, gl.DYNAMIC_DRAW),
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
      { loc: 7, buffer: b.lasso, size: 2, type: gl.UNSIGNED_BYTE, normalized: true, divisor: 1 },
    ]);
  }

  /**
   * The lasso per object: `to` its level in the selection (255 shown, 128
   * dimmed, 0 hidden; null shows all) and, while a fade plays, `from` the
   * weight it starts from (0–255, as `lassoDrawn` gives; null at rest).
   */
  setLasso(to, from = null) {
    const n = this.count, l = this.lasso;
    for (let i = 0; i < n; i++) {
      const level = to ? to[i] : 255;
      l[i * 2] = from ? from[i] : level;
      l[i * 2 + 1] = level;
    }
    this.lassoVersion = (this.lassoVersion || 0) + 1;
    const gl = this.gl;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.lasso);
    gl.bufferSubData(gl.ARRAY_BUFFER, 0, l);
  }

  /** The lasso's opacity weight of object `i` now, as the shader draws it (for CPU readers such as member stars). */
  lassoWeight(i, mix, dim) {
    return lassoBlend(this.lasso[i * 2], this.lasso[i * 2 + 1], mix, dim);
  }

  /** Every object's lasso weight as drawn now (0–255): where a new fade starts. */
  lassoDrawn(mix, dim) {
    const out = new Uint8Array(this.count);
    for (let i = 0; i < out.length; i++) out[i] = Math.round(255 * this.lassoWeight(i, mix, dim));
    return out;
  }

  /**
   * Replace the bits in `mask` with `values` (Uint8Array, same length; null
   * clears them). Only a change uploads and bumps `stateVersion`; returns
   * whether anything changed.
   */
  setStateBits(mask, values) {
    const s = this.state;
    let changed = false;
    for (let i = 0; i < s.length; i++) {
      const b = (s[i] & ~mask) | (values ? values[i] & mask : 0);
      if (b !== s[i]) { s[i] = b; changed = true; }
    }
    if (!changed) return false;
    this.stateVersion = (this.stateVersion || 0) + 1;
    const gl = this.gl;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.state);
    gl.bufferSubData(gl.ARRAY_BUFFER, 0, s);
    return true;
  }

  dispose() {
    const gl = this.gl;
    gl.deleteVertexArray(this.vao);
    if (this.trailVao) gl.deleteVertexArray(this.trailVao);
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
    this.trailProgram = null; // compiled on first use
    this.trailBuffer = null;
    this.lassoMix = [1, 1];
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
    prog.f("uOpacityScale", style.opacityScale * global.pointOpacity * Math.min(1, style.presence ?? 1));
    prog.f("uGlow", global.glow);
    prog.v3("uFade", global.fadeTime, global.fadeInOut ? 1 : 0, global.fadeByOpacity ? 1 : 0);
    prog.f("uStarsExp", style.starsExp);
    prog.i("uHover", style.hover ?? -1);
    prog.f("uHoverBoost", 1.0);
    prog.f("uMinPx", 1.25 * frame.dpr);
    prog.f("uPickMinPx", PICK_MIN_PX * frame.dpr);
    prog.f("uDimOpacity", style.dimOpacity ?? 0.16);
    prog.i("uHasAge", batch.hasAge ? 1 : 0);
    prog.f("uMemberFade", global.memberFade ?? 0);
    prog.v2("uLassoMix", this.lassoMix[0], this.lassoMix[1]);
  }

  /** How far a State change's lasso fade has got: fade-out and fade-in, 0–1 (1, 1 at rest). */
  setLassoMix(out, inn) {
    this.lassoMix[0] = out;
    this.lassoMix[1] = inn;
  }

  /** Motion trails, drawn under the points. */
  drawTrails(frame) {
    const p = this.params;
    const tr = p?.trails;
    if (!tr || Math.abs(tr.span) < 0.05) return;
    const gl = this.gl;
    if (!this.trailProgram) {
      this.trailProgram = createProgram(gl, TRAIL_VERT, TRAIL_FRAG, { label: "point-trails" });
      this.trailBuffer = createBuffer(gl, trailTemplate());
    }
    const prog = this.trailProgram.use();
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.enable(gl.DEPTH_TEST);
    gl.depthMask(false);
    const cam = frame.camera;
    prog.m4("uViewProj", cam.viewProj).v2("uViewport", frame.width, frame.height);
    prog.f("uFrame", p.frame).f("uTime", p.time).f("uSpan", tr.span).f("uMyrPerFrame", tr.myrPerFrame);
    prog.f("uWidthPx", (1.2 + 0.8 * Math.min(p.pointSize, 3)) * frame.dpr);
    prog.f("uMemberFade", p.memberFade ?? 0);
    prog.v2("uLassoMix", this.lassoMix[0], this.lassoMix[1]);
    prog.f("uStrength", Math.min(1.6, 0.55 + 0.35 * p.glow));
    for (const batch of this.batches) {
      if (batch.posFrames <= 1) continue;
      const style = p.styles.get(batch.key);
      if (!style || !style.visible || style.presence <= 0.001 || style.trails === false) continue;
      if (!batch.trailVao) {
        const b = batch.buffers;
        batch.trailVao = createVAO(gl, [
          { loc: 0, buffer: this.trailBuffer, size: 2 },
          { loc: 1, buffer: b.index, size: 1, divisor: 1 },
          { loc: 2, buffer: b.color, size: 3, type: gl.UNSIGNED_BYTE, normalized: true, divisor: 1 },
          { loc: 3, buffer: b.stat, size: 4, divisor: 1 },
          { loc: 6, buffer: b.state, size: 1, type: gl.UNSIGNED_BYTE, divisor: 1 },
          { loc: 7, buffer: b.lasso, size: 2, type: gl.UNSIGNED_BYTE, normalized: true, divisor: 1 },
        ]);
      }
      prog.tex("uFrameTex", batch.frameTex.texture);
      prog.i("uCount", batch.count).i("uFrames", batch.texFrames);
      if (batch.hasAux) {
        prog.tex("uAuxTex", batch.auxTex.texture).i("uAuxFrames", batch.auxFrames).i("uHasAux", 1);
      } else {
        prog.tex("uAuxTex", batch.frameTex.texture).i("uAuxFrames", 1).i("uHasAux", 0);
      }
      const off = frameOffset(batch.offsets, p.frame);
      prog.v3("uOffset", off[0], off[1], off[2]);
      prog.f("uOpacityScale", style.opacityScale * p.pointOpacity * Math.min(1, style.presence ?? 1));
      prog.f("uDimOpacity", style.dimOpacity ?? 0.16);
      prog.i("uHasAge", batch.hasAge ? 1 : 0);
      const colorMode = style.colorMode === "by_value" && batch.trace.colorBy ? 1 : 0;
      prog.i("uColorMode", colorMode);
      if (colorMode) {
        prog.tex("uCmap", this.cmapTextures.get(style.colormap) || this.fallbackCmap);
        prog.v2("uCRange", style.cmin, style.cmax);
      } else {
        prog.tex("uCmap", this.fallbackCmap);
      }
      const c = style.colorOverride ? parseColor(style.colorOverride) : parseColor(batch.trace.color);
      prog.v3("uTraceColor", c[0], c[1], c[2]);
      prog.i("uUseObjectColor", style.colorOverride ? 0 : 1);
      gl.bindVertexArray(batch.trailVao);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, (TRAIL_SEGMENTS + 1) * 2, batch.count);
    }
    gl.bindVertexArray(null);
  }

  draw(frame) {
    const p = this.params;
    if (!p) return;
    this.drawTrails(frame);
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
    this.trailProgram?.dispose();
    if (this.trailBuffer) this.gl.deleteBuffer(this.trailBuffer);
    for (const tex of this.cmapTextures.values()) this.gl.deleteTexture(tex);
    this.cmapTextures.clear();
    this.gl.deleteTexture(this.fallbackCmap);
  }
}

/**
 * The shader's lasso weight on the CPU: from `a`, the weight a fade starts
 * from (0–255), toward the selection's level `b` (255 shown, 128 dimmed,
 * 0 hidden), by the fade-out (`mix[0]`, while the weight falls) or fade-in
 * progress (`mix[1]`).
 */
export function lassoBlend(a, b, mix, dim) {
  const wa = a / 255, wb = b > 191 ? 1 : b > 63 ? dim : 0;
  return wa + (wb - wa) * (wb < wa ? mix[0] : mix[1]);
}
