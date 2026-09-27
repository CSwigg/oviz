// Cluster member stars, drawn only in Sky view.
//
// Each star keeps its catalogue position (heliocentric pc at t = 0). At time
// t it is displaced rigidly with its parent cluster's orbit, so members
// inherit the parent's motion, birth fade, visibility and colour. Per-cluster
// state lives in two small textures updated on the CPU (thousands of
// clusters), while the hundreds of thousands of stars never leave the GPU.

import { createProgram, createBuffer, createVAO, createDataTexture, createTexture2D } from "../engine/gl.js";
import { starTextures } from "../layers/points.js";

const TEX_W = 1024;

const VERT = `
layout(location = 0) in vec2 aCorner;
layout(location = 1) in vec3 aPos;
layout(location = 2) in vec3 aVel;
layout(location = 3) in float aCluster;
uniform mat4 uViewProj;
uniform mat4 uView;
uniform float uProjY;
uniform vec2 uViewport;
uniform highp sampler2D uClusterState;   // dx, dy, dz, alpha
uniform sampler2D uClusterColor;         // rgb, size
uniform float uTime;
uniform float uInternal;
uniform float uStarSize;
uniform float uMinPx;
uniform float uReveal;
uniform float uGlow;
out vec2 vUv;
out vec3 vColor;
out float vA;
flat out int vId;
void main() {
  if (aCluster < -0.5) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }
  int c = int(aCluster + 0.5);
  ivec2 tc = ivec2(c % ${TEX_W}, c / ${TEX_W});
  vec4 st = texelFetch(uClusterState, tc, 0);
  vec4 col = texelFetch(uClusterColor, tc, 0);
  float a = st.w * uReveal;
  vec3 p = aPos + st.xyz + aVel * uTime * uInternal;
  vec4 clip = uViewProj * vec4(p, 1.0);
  vId = gl_InstanceID;
  if (a <= 0.002 || clip.w <= 0.0) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }
  float viewZ = -(uView * vec4(p, 1.0)).z;
  float world = uStarSize * (0.35 + col.a * 1.3);
  float px = world * 0.5 * uViewport.y * uProjY / max(viewZ, 1e-3);
  float k = 1.0;
  #ifdef PICK
  px = max(px, uMinPx * 3.0);
  #else
  if (px < uMinPx) { k = (px / uMinPx); k *= k; px = uMinPx; }
  #endif
  gl_Position = clip + vec4(aCorner * px / uViewport * clip.w, 0.0, 0.0);
  vUv = aCorner;
  vColor = col.rgb;
  vA = a * k;
}
`;

const FRAG = `
in vec2 vUv;
in vec3 vColor;
in float vA;
flat in int vId;
uniform sampler2D uHalo;
uniform sampler2D uCore;
uniform float uGlow;
out vec4 outColor;
void main() {
  float r = length(vUv);
  if (r > 1.0) discard;
  #ifdef PICK
  int id = vId + 1;
  outColor = vec4(200.0 / 255.0, float((id >> 16) & 255) / 255.0, float((id >> 8) & 255) / 255.0, float(id & 255) / 255.0);
  #else
  vec2 uv = vUv * 0.5 + 0.5;
  float halo = texture(uHalo, uv).r * vA * clamp(0.45 + 0.2 * uGlow, 0.0, 0.85);
  float core = texture(uCore, uv).r * vA;
  if (halo + core < 0.002) discard;
  outColor = vec4(vColor * halo + mix(vColor, vec3(1.0), 0.65) * core, halo);
  #endif
}
`;

const QUAD = new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]);

export class MembersLayer {
  constructor(gl, { count, clusters, position, velocity, valid, offsets }) {
    this.gl = gl;
    this.order = 32;
    this.count = count;
    this.clusterCount = clusters.length;
    this.clusterNames = clusters;
    this.program = createProgram(gl, VERT, FRAG, { label: "members" });
    this.pickProgram = createProgram(gl, VERT, FRAG, { label: "members-pick", defines: { PICK: 1 } });
    this.textures = starTextures(gl);
    this.positions = position;
    this.offsets = offsets;
    // Cluster index per star, and velocities relative to the cluster mean
    // (internal kinematics, as in the legacy motion arrows).
    const clusterOf = new Float32Array(count);
    const rel = new Float32Array(count * 3);
    this.clusterOfStar = new Uint32Array(count);
    for (let c = 0; c < this.clusterCount; c++) {
      const a = offsets[c], b = offsets[c + 1];
      let mx = 0, my = 0, mz = 0, n = 0;
      for (let i = a; i < b; i++) {
        if (!valid[i]) continue;
        mx += velocity[i * 3]; my += velocity[i * 3 + 1]; mz += velocity[i * 3 + 2]; n++;
      }
      if (n) { mx /= n; my /= n; mz /= n; }
      for (let i = a; i < b; i++) {
        clusterOf[i] = valid[i] ? c : -1;
        this.clusterOfStar[i] = c;
        rel[i * 3] = velocity[i * 3] - mx;
        rel[i * 3 + 1] = velocity[i * 3 + 1] - my;
        rel[i * 3 + 2] = velocity[i * 3 + 2] - mz;
      }
    }
    this.relVelocity = rel;
    const texH = Math.max(1, Math.ceil(this.clusterCount / TEX_W));
    this.stateData = new Float32Array(TEX_W * texH * 4);
    this.colorData = new Uint8Array(TEX_W * texH * 4);
    this.stateTex = createDataTexture(gl, this.stateData, TEX_W, texH, 4);
    this.colorTex = createTexture2D(gl, { width: TEX_W, height: texH, internalFormat: gl.RGBA8, format: gl.RGBA, type: gl.UNSIGNED_BYTE, data: this.colorData, filter: gl.NEAREST });
    this.texH = texH;
    const buf = {
      quad: createBuffer(gl, QUAD),
      pos: createBuffer(gl, position),
      vel: createBuffer(gl, rel),
      cluster: createBuffer(gl, clusterOf),
    };
    this.buffers = buf;
    this.vao = createVAO(gl, [
      { loc: 0, buffer: buf.quad, size: 2 },
      { loc: 1, buffer: buf.pos, size: 3, divisor: 1 },
      { loc: 2, buffer: buf.vel, size: 3, divisor: 1 },
      { loc: 3, buffer: buf.cluster, size: 1, divisor: 1 },
    ]);
    this.params = null;
    this.visible = false;
  }

  /** Upload per-cluster displacement/alpha (Float32 ×4) and colour (u8 ×4). */
  updateClusters(fill) {
    const gl = this.gl;
    fill(this.stateData, this.colorData);
    gl.bindTexture(gl.TEXTURE_2D, this.stateTex);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, TEX_W, this.texH, gl.RGBA, gl.FLOAT, this.stateData);
    gl.bindTexture(gl.TEXTURE_2D, this.colorTex);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, TEX_W, this.texH, gl.RGBA, gl.UNSIGNED_BYTE, this.colorData);
  }

  _uniforms(prog, frame) {
    const p = this.params;
    const cam = frame.camera;
    prog.m4("uViewProj", cam.viewProj).m4("uView", cam.view).f("uProjY", cam.proj[5]);
    prog.v2("uViewport", frame.width, frame.height);
    prog.tex("uClusterState", this.stateTex).tex("uClusterColor", this.colorTex);
    prog.f("uTime", p.time).f("uInternal", p.internalMotion ? 1 : 0);
    prog.f("uStarSize", p.starSize).f("uMinPx", p.minPx * frame.dpr).f("uReveal", p.reveal).f("uGlow", p.glow);
  }

  draw(frame) {
    const p = this.params;
    if (!p || p.reveal <= 0.002) return;
    const gl = this.gl;
    const prog = this.program.use();
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.enable(gl.DEPTH_TEST);
    gl.depthMask(false);
    this._uniforms(prog, frame);
    prog.tex("uHalo", this.textures.memberHalo).tex("uCore", this.textures.memberCore);
    gl.bindVertexArray(this.vao);
    gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, this.count);
    gl.bindVertexArray(null);
  }

  drawPick(frame) {
    const p = this.params;
    if (!p || p.reveal <= 0.5) return;
    const gl = this.gl;
    const prog = this.pickProgram.use();
    gl.disable(gl.BLEND);
    this._uniforms(prog, frame);
    gl.bindVertexArray(this.vao);
    gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, this.count);
    gl.bindVertexArray(null);
  }

  /** World position of star i at time t (matches the shader). */
  starPosition(i, time, internal) {
    const c = this.clusterOfStar[i];
    const s = this.stateData;
    const P = this.positions, V = this.relVelocity;
    const k = internal ? time : 0;
    return [
      P[i * 3] + s[c * 4] + V[i * 3] * k,
      P[i * 3 + 1] + s[c * 4 + 1] + V[i * 3 + 1] * k,
      P[i * 3 + 2] + s[c * 4 + 2] + V[i * 3 + 2] * k,
    ];
  }

  dispose() {
    const gl = this.gl;
    gl.deleteVertexArray(this.vao);
    for (const b of Object.values(this.buffers)) gl.deleteBuffer(b);
    gl.deleteTexture(this.stateTex);
    gl.deleteTexture(this.colorTex);
    this.program.dispose();
    this.pickProgram.dispose();
  }
}
