// Anti-aliased screen-space thick lines (instanced segments).
// Each instance reads both endpoints from the trace's frame texture, so
// time-varying lines interpolate on the GPU exactly like points.

import { createProgram, createBuffer, createVAO } from "../engine/gl.js";
import { FRAME_GLSL, packFrameTexture, frameOffset } from "../engine/frames.js";
import { parseColor } from "../core/color.js";

const VERT = `
${FRAME_GLSL}
layout(location = 0) in vec2 aCorner;     // x: 0 start / 1 end, y: -1..1 side
layout(location = 1) in float aSegment;
layout(location = 2) in vec2 aArc;
layout(location = 3) in vec3 aColor;      // per-segment colour (when uSegColor)
uniform mat4 uViewProj;
uniform vec2 uViewport;
uniform highp sampler2D uFrameTex;
uniform int uCount;
uniform int uFrames;
uniform float uFrame;
uniform vec3 uOffset;
uniform float uWidth;
uniform vec3 uColor;
uniform int uSegColor;
out float vSide;
out vec3 vColor;
out float vArc;

void main() {
  int s = int(aSegment + 0.5);
  float d;
  vec4 a = framePosition(uFrameTex, 2 * s, uCount, uFrames, uFrame, d);
  vec4 b = framePosition(uFrameTex, 2 * s + 1, uCount, uFrames, uFrame, d);
  if (a.w <= 0.0 || b.w <= 0.0) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }
  vec4 ca = uViewProj * vec4(a.xyz + uOffset, 1.0);
  vec4 cb = uViewProj * vec4(b.xyz + uOffset, 1.0);
  // Clip the segment against the near plane so lines passing behind the
  // camera never explode across the screen.
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
  float half_ = uWidth * 0.5 + 1.0;
  vec4 c = aCorner.x < 0.5 ? ca : cb;
  vec2 sp = (aCorner.x < 0.5 ? sa : sb) + n * aCorner.y * half_ + dir * (aCorner.x < 0.5 ? -1.0 : 1.0) * 0.5;
  gl_Position = vec4(sp / (uViewport * 0.5) * c.w, c.z, c.w);
  vSide = aCorner.y * half_;
  vArc = aCorner.x < 0.5 ? aArc.x : aArc.y;
  vColor = uSegColor == 1 ? aColor : uColor;
}
`;

const FRAG = `
in float vSide;
in float vArc;
in vec3 vColor;
uniform float uOpacity;
uniform float uWidth;
uniform vec2 uDash;   // dash length, gap length (world units); 0 = solid
out vec4 outColor;
void main() {
  float d = abs(vSide);
  float a = clamp(uWidth * 0.5 + 0.5 - d, 0.0, 1.0);
  if (uWidth < 1.0) a *= uWidth;
  if (uDash.x > 0.0) {
    float period = uDash.x + uDash.y;
    float m = mod(vArc, period);
    float edge = fwidth(vArc) + 1e-6;
    a *= 1.0 - smoothstep(uDash.x - edge, uDash.x + edge, m);
  }
  a *= uOpacity;
  if (a < 0.003) discard;
  outColor = vec4(vColor * a, a);
}
`;

const CORNERS = new Float32Array([0, -1, 1, -1, 0, 1, 1, 1]);

export class LinesLayer {
  constructor(gl, { frames }) {
    this.gl = gl;
    this.order = 20;
    this.frames = frames;
    this.batches = [];
    this.byKey = new Map();
    this.program = createProgram(gl, VERT, FRAG, { label: "lines" });
    this.cornerBuffer = createBuffer(gl, CORNERS);
    this._offset = [0, 0, 0]; // scratch for per-frame offsets
    this.params = null;
  }

  addTrace(trace, data) {
    const gl = this.gl;
    const L = trace.lines;
    const S = L.count;
    const posFrames = L.position.frames || 1;
    const frameTex = packFrameTexture(gl, {
      positions: data.position, posFrames, count: 2 * S, frames: posFrames,
    });
    const seg = new Float32Array(S);
    for (let i = 0; i < S; i++) seg[i] = i;
    const arc = new Float32Array(S * 2);
    if (data.arc) arc.set(data.arc.subarray(0, S * 2));
    const segBuf = createBuffer(gl, seg);
    const arcBuf = createBuffer(gl, arc);
    const attrs = [
      { loc: 0, buffer: this.cornerBuffer, size: 2 },
      { loc: 1, buffer: segBuf, size: 1, divisor: 1 },
      { loc: 2, buffer: arcBuf, size: 2, divisor: 1 },
    ];
    const buffers = [segBuf, arcBuf];
    if (data.segmentColors) {
      const rgb = new Uint8Array(S * 3);
      rgb.set(data.segmentColors.subarray(0, S * 3));
      const colorBuf = createBuffer(gl, rgb);
      attrs.push({ loc: 3, buffer: colorBuf, size: 3, divisor: 1, type: gl.UNSIGNED_BYTE, normalized: true });
      buffers.push(colorBuf);
    }
    const vao = createVAO(gl, attrs);
    // Dash lengths relative to the line's total extent.
    let total = 0;
    for (let i = 0; i < arc.length; i++) total = Math.max(total, arc[i]);
    const unit = total > 0 ? total / 160 : 1;
    const dash = { solid: [0, 0], dash: [unit * 2.2, unit * 1.4], dot: [unit * 0.45, unit * 0.9], dashdot: [unit * 2, unit], longdash: [unit * 4, unit * 1.4] }[L.dash] || [0, 0];
    const batch = {
      key: trace.key, trace, count: S, frameTex, posFrames, vao, buffers, segmentColors: !!data.segmentColors,
      offsets: data.offset || null, dash, color: parseColor(L.color || trace.color),
    };
    this.batches.push(batch);
    this.byKey.set(trace.key, batch);
    return batch;
  }

  removeTrace(key) {
    const gl = this.gl;
    const b = this.byKey.get(key);
    if (!b) return;
    gl.deleteVertexArray(b.vao);
    b.buffers.forEach((x) => gl.deleteBuffer(x));
    gl.deleteTexture(b.frameTex.texture);
    this.byKey.delete(key);
    this.batches = this.batches.filter((x) => x !== b);
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
    prog.m4("uViewProj", frame.camera.viewProj).v2("uViewport", frame.width, frame.height);
    for (const b of this.batches) {
      const style = p.styles.get(b.key);
      if (!style || !style.visible || style.presence <= 0.001) continue;
      const lineStyle = b.trace.lines;
      const opacity = (lineStyle.opacity ?? 1) * style.opacityScale * style.presence * (p.lineOpacity ?? 1);
      if (opacity <= 0.002) continue;
      const c = style.colorOverride ? parseColor(style.colorOverride) : b.color;
      prog.tex("uFrameTex", b.frameTex.texture);
      prog.i("uCount", 2 * b.count).i("uFrames", b.posFrames).f("uFrame", p.frame);
      const off = frameOffset(b.offsets, p.frame, this._offset);
      prog.v3("uOffset", off[0], off[1], off[2]);
      prog.f("uWidth", Math.max(0.6, (lineStyle.width || 1) * (style.sizeScale || 1)) * frame.dpr);
      // A colour picked for the layer replaces its per-segment colours.
      prog.v3("uColor", c[0], c[1], c[2]).f("uOpacity", opacity).i("uSegColor", b.segmentColors && !style.colorOverride ? 1 : 0);
      prog.v2("uDash", b.dash[0], b.dash[1]);
      gl.bindVertexArray(b.vao);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, b.count);
    }
    gl.bindVertexArray(null);
  }

  dispose() {
    const gl = this.gl;
    for (const b of this.batches) {
      gl.deleteVertexArray(b.vao);
      b.buffers.forEach((x) => gl.deleteBuffer(x));
      gl.deleteTexture(b.frameTex.texture);
    }
    gl.deleteBuffer(this.cornerBuffer);
    this.program.dispose();
  }
}
