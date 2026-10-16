// Animated flow lines: static streamlines with light pulses running along
// them in integration time (see oviz/viewer/flow.py).
//
// Every vertex carries (x, y, z, line) and (tau, value, trust, edge). Each
// instance is one segment, read straight from the vertex buffer twice (the
// second read one vertex later), so the GPU holds the geometry once and a
// frame of animation costs one uniform: the flow's clock. Segments that
// would join two different lines collapse in the vertex shader.
//
// What readers can change lives in the trace's state (see flowStyle): the
// pulse rate, spacing and trail, the rails' brightness, the colour (solid
// or by value through a colormap) and the colour range.

import { createProgram, createBuffer, createVAO } from "../engine/gl.js";
import { parseColor } from "../core/color.js";

const VERT = `
layout(location = 0) in vec2 aCorner;     // x: 0 start / 1 end, y: -1..1 side
layout(location = 1) in vec4 aPosA;       // xyz, line id
layout(location = 2) in vec4 aPosB;
layout(location = 3) in vec4 aAttrA;      // tau (Myr), value, trust, edge
layout(location = 4) in vec4 aAttrB;
uniform mat4 uViewProj;
uniform vec2 uViewport;
uniform float uWidth;
uniform float uFocus;
out float vSide;
out float vHalf;
out float vTau;
out float vValue;
out float vTrust;
out float vEdge;
out float vPhase;
out float vDepth;

float hash(float n) { return fract(sin(n * 12.9898 + 4.1414) * 43758.5453); }

void main() {
  if (aPosA.w != aPosB.w) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }
  vec4 ca = uViewProj * vec4(aPosA.xyz, 1.0);
  vec4 cb = uViewProj * vec4(aPosB.xyz, 1.0);
  float near = 1e-4;
  if (ca.w < near && cb.w < near) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }
  if (ca.w < near) ca = mix(ca, cb, (near - ca.w) / (cb.w - ca.w));
  if (cb.w < near) cb = mix(cb, ca, (near - cb.w) / (ca.w - cb.w));
  vec2 sa = ca.xy / ca.w * uViewport * 0.5;
  vec2 sb = cb.xy / cb.w * uViewport * 0.5;
  vec2 dir = sb - sa;
  float len = length(dir);
  dir = len > 1e-6 ? dir / len : vec2(1.0, 0.0);
  vec2 n = vec2(-dir.y, dir.x);
  bool atB = aCorner.x > 0.5;
  vec4 c = atB ? cb : ca;
  // Nearer lines draw a little wider: a quiet depth cue.
  float depthScale = clamp(uFocus / max(c.w, 1e-3), 0.7, 1.35);
  float half_ = uWidth * 0.5 * depthScale + 1.0;
  vec2 sp = (atB ? sb : sa) + n * aCorner.y * half_ + dir * (atB ? 0.5 : -0.5);
  gl_Position = vec4(sp / (uViewport * 0.5) * c.w, c.z, c.w);
  vec4 attr = atB ? aAttrB : aAttrA;
  vSide = aCorner.y * half_;
  vHalf = half_;
  vTau = attr.x;
  vValue = attr.y;
  vTrust = attr.z;
  vEdge = attr.w;
  vPhase = hash(aPosA.w);
  vDepth = c.w;
}
`;

const FRAG = `
in float vSide;
in float vHalf;
in float vTau;
in float vValue;
in float vTrust;
in float vEdge;
in float vPhase;
in float vDepth;
uniform sampler2D uLut;
uniform int uUseLut;
uniform vec3 uColor;
uniform vec2 uRange;
uniform float uOpacity;
uniform float uRail;
uniform float uClock;
uniform float uPeriod;
uniform float uTail;
uniform float uTrustFloor;
uniform float uFocus;
out vec4 outColor;

void main() {
  float d = abs(vSide);
  float cover = clamp(vHalf - d, 0.0, 1.0);
  // Distance behind the nearest pulse head, in Myr of integration time.
  float behind = fract((uClock - vTau) / uPeriod + vPhase) * uPeriod;
  float tail = exp(-behind / uTail);
  // Soft leading edge so heads never alias into a hard step.
  float lead = smoothstep(uPeriod, uPeriod - max(fwidth(behind) * 1.5, 0.05 * uTail), behind);
  float pulse = tail * lead;
  float s = clamp((vValue - uRange.x) / max(uRange.y - uRange.x, 1e-6), 0.0, 1.0);
  vec3 col = uUseLut == 1 ? texture(uLut, vec2(0.02 + 0.96 * s, 0.5)).rgb : uColor;
  // The head burns toward white; the tail keeps the line's colour.
  col = mix(col, vec3(1.0), 0.3 * pulse * pulse * pulse);
  float trust = mix(uTrustFloor, 1.0, smoothstep(0.0, 1.0, vTrust));
  float fog = mix(1.0, 0.45, smoothstep(uFocus * 0.7, uFocus * 2.2, vDepth));
  // The core of the line is brighter than its rim.
  float core = 1.0 - 0.45 * smoothstep(0.0, vHalf, d);
  float a = uOpacity * cover * core * vEdge * trust * fog * (uRail + (1.0 - uRail) * pulse);
  if (a < 0.003) discard;
  outColor = vec4(col * a, a * 0.85);
}
`;

const CORNERS = new Float32Array([0, -1, 1, -1, 0, 1, 1, 1]);

const finite = (v, fb) => (v !== null && v !== undefined && v !== "" && Number.isFinite(Number(v)) ? Number(v) : fb);

/**
 * A flow trace's display settings: its authored defaults (`trace.flow`)
 * under the reader's choices in its trace state `ts`:
 *   colorMode "by_value" | "fixed", colormap, cmin, cmax, color,
 *   flowRate (Myr of integration time per second), flowPeriod (Myr between
 *   pulses), flowTail (Myr of trail), flowRail (0–1 brightness of the lines
 *   between pulses).
 */
export function flowStyle(trace, ts = {}) {
  const F = trace.flow || {};
  const A = F.animation || {};
  const [lo, hi] = F.speedRange || [0, 1];
  const colormap = ts.colormap || F.colormap || null;
  const mode = ts.colorMode || (F.fixedColor || !F.colormap ? "fixed" : "by_value");
  return {
    byValue: mode === "by_value" && !!colormap,
    colormap,
    color: ts.color || trace.color,
    range: [finite(ts.cmin, lo), finite(ts.cmax, hi)],
    rate: Math.max(0, finite(ts.flowRate, A.rate ?? 1.2)),
    period: Math.max(finite(ts.flowPeriod, A.period ?? 14), 0.05),
    tail: Math.max(finite(ts.flowTail, A.tail ?? 3.5), 0.02),
    rail: Math.min(Math.max(finite(ts.flowRail, F.railOpacity ?? 0.16), 0), 1),
  };
}

export class FlowLayer {
  /** `cmaps`: the figure's colormap textures by name (shared with volumes). */
  constructor(gl, { cmaps } = {}) {
    this.gl = gl;
    this.order = 35;
    // Animated: the renderer redraws only this layer while its pulses move.
    this.animated = true;
    this.batches = [];
    this.byKey = new Map();
    this.cmaps = cmaps || new Map();
    // Each flow's clock (Myr of integration time), advanced at its own rate
    // so a rate change never makes the pulses jump.
    this.clocks = new Map();
    this.program = createProgram(gl, VERT, FRAG, { label: "flow" });
    this.cornerBuffer = createBuffer(gl, CORNERS);
    this.params = null;
  }

  /** `data`: {position: Float32Array(V*4), attr: Float32Array(V*4)}. */
  addTrace(trace, data) {
    const gl = this.gl;
    const F = trace.flow;
    const V = F.vertices;
    const posBuf = createBuffer(gl, data.position.subarray(0, V * 4));
    const attrBuf = createBuffer(gl, data.attr.subarray(0, V * 4));
    const stride = 16;
    const vao = createVAO(gl, [
      { loc: 0, buffer: this.cornerBuffer, size: 2 },
      { loc: 1, buffer: posBuf, size: 4, stride, offset: 0, divisor: 1 },
      { loc: 2, buffer: posBuf, size: 4, stride, offset: stride, divisor: 1 },
      { loc: 3, buffer: attrBuf, size: 4, stride, offset: 0, divisor: 1 },
      { loc: 4, buffer: attrBuf, size: 4, stride, offset: stride, divisor: 1 },
    ]);
    const batch = { key: trace.key, trace, segments: V - 1, vao, buffers: [posBuf, attrBuf] };
    this.batches.push(batch);
    this.byKey.set(trace.key, batch);
    return batch;
  }

  /** Move every shown flow's pulses on by `seconds` of playback. */
  advance(seconds, styles) {
    for (const b of this.batches) {
      const s = styles?.get(b.key);
      if (!s?.flow || !this._shown(s)) continue;
      // Pulses repeat every period: wrapping keeps float32 exact forever.
      this.clocks.set(b.key, ((this.clocks.get(b.key) || 0) + seconds * s.flow.rate) % s.flow.period);
    }
  }

  _shown(s) {
    return !!s && s.visible && s.presence > 0.001 && s.opacityScale > 0.002;
  }

  /** Whether a flow is on screen (the renderer caches the rest of the scene then). */
  active() {
    return this.anyVisible(this.params?.styles);
  }

  /** Whether any flow is on screen (so the viewer knows to keep animating). */
  anyVisible(styles) {
    return this.batches.some((b) => this._shown(styles?.get(b.key)));
  }

  draw(frame) {
    const p = this.params;
    if (!p || !this.batches.length) return;
    const gl = this.gl;
    const prog = this.program.use();
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.enable(gl.DEPTH_TEST);
    gl.depthMask(false);
    const cam = frame.camera;
    // The orbit distance sets the depth scale; in Sky view the eye sits on
    // the target, so fall back to the scene scale.
    const focus = Math.max(cam.pose?.distance || 0, p.focusFallback || 500);
    prog.m4("uViewProj", cam.viewProj).v2("uViewport", frame.width, frame.height).f("uFocus", focus);
    for (const b of this.batches) {
      const style = p.styles.get(b.key);
      if (!this._shown(style)) continue;
      const F = b.trace.flow;
      const fs = style.flow || flowStyle(b.trace);
      const opacity = (b.trace.opacity ?? 1) * style.opacityScale * style.presence;
      if (opacity <= 0.002) continue;
      const lut = fs.byValue ? this.cmaps.get(fs.colormap) || this.cmaps.get(F.colormap) : null;
      const c = parseColor(fs.color || b.trace.color);
      if (lut) prog.tex("uLut", lut);
      prog.i("uUseLut", lut ? 1 : 0).v3("uColor", c[0], c[1], c[2]);
      prog.v2("uRange", fs.range[0], fs.range[1]);
      prog.f("uOpacity", Math.min(opacity, 1.5)).f("uRail", fs.rail);
      prog.f("uClock", (this.clocks.get(b.key) || 0) % fs.period);
      prog.f("uPeriod", fs.period).f("uTail", fs.tail).f("uTrustFloor", F.trustFloor ?? 0.35);
      prog.f("uWidth", Math.max(0.6, (F.width || 1.6) * (style.sizeScale || 1)) * frame.dpr);
      gl.bindVertexArray(b.vao);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, b.segments);
    }
    gl.bindVertexArray(null);
  }

  dispose() {
    const gl = this.gl;
    for (const b of this.batches) {
      gl.deleteVertexArray(b.vao);
      b.buffers.forEach((x) => gl.deleteBuffer(x));
    }
    this.batches = [];
    this.byKey.clear();
    gl.deleteBuffer(this.cornerBuffer);
    this.program.dispose();
  }
}
