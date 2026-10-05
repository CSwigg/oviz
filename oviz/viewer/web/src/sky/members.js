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
// A click picks the star, not its glow: the member core's disc, out to where
// its profile has fallen to a few percent (in sprite radius units).
export const MEMBER_PICK_CORE = 0.14;

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
  px = max(px * ${MEMBER_PICK_CORE}, uMinPx * 3.0);
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

// Motion arrows (classic Oviz): each member star's sky-plane motion in its
// cluster's rest frame, drawn as 3 Myr of travel (at most 0.6 × its
// distance), half as wide as 4.5e-4 × distance, with a head 2.6× wider.
// Per cluster: length and width scales (0 length = no arrows).
const ARROW_VERT = `
layout(location = 1) in vec3 aPos;
layout(location = 2) in vec3 aVel;
layout(location = 3) in float aCluster;
uniform mat4 uViewProj;
uniform highp sampler2D uClusterState;   // dx, dy, dz, alpha
uniform sampler2D uClusterColor;         // rgb, size
uniform highp sampler2D uArrowStyle;     // length scale, width scale
uniform float uTime;
uniform float uInternal;
uniform float uReveal;
out vec4 vColor;
void main() {
  gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
  vColor = vec4(0.0);
  if (aCluster < -0.5) return;
  int c = int(aCluster + 0.5);
  ivec2 tc = ivec2(c % ${TEX_W}, c / ${TEX_W});
  vec4 style = texelFetch(uArrowStyle, tc, 0);
  vec4 st = texelFetch(uClusterState, tc, 0);
  float a = st.w * uReveal * 0.9;
  if (style.x <= 0.0 || a <= 0.002) return;
  vec3 p = aPos + st.xyz + aVel * uTime * uInternal;
  float dist = length(p);
  if (dist <= 1e-6) return;
  vec3 radial = p / dist;
  vec3 tang = aVel - radial * dot(aVel, radial);
  float speed = length(tang);
  if (speed <= 1e-9) return;
  vec3 dir = tang / speed;
  vec3 perp = cross(radial, dir);
  if (dot(perp, perp) <= 1e-12) return;
  perp = normalize(perp);
  float len = min(speed * 3.0 * style.x, dist * 0.6);
  float halfW = dist * 0.00045 * style.y;
  float head = min(max(0.30 * len, halfW * 2.0), len * 0.6);
  float shaft = max(len - head, 0.0);
  vec3 s = p + dir * shaft;
  vec3 tip = p + dir * len;
  vec3 q;
  int k = gl_VertexID;
  // Shaft: (p+w, p−w, s−w), (p+w, s−w, s+w); head: (s+2.6w, s−2.6w, tip).
  if (k == 0 || k == 3) q = p + perp * halfW;
  else if (k == 1) q = p - perp * halfW;
  else if (k == 2 || k == 4) q = s - perp * halfW;
  else if (k == 5) q = s + perp * halfW;
  else if (k == 6) q = s + perp * halfW * 2.6;
  else if (k == 7) q = s - perp * halfW * 2.6;
  else q = tip;
  gl_Position = uViewProj * vec4(q, 1.0);
  vColor = vec4(texelFetch(uClusterColor, tc, 0).rgb, a);
}
`;

const ARROW_FRAG = `
in vec4 vColor;
out vec4 outColor;
void main() {
  if (vColor.a <= 0.002) discard;
  outColor = vec4(vColor.rgb * vColor.a, vColor.a);
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
    this.arrowProgram = createProgram(gl, ARROW_VERT, ARROW_FRAG, { label: "member-arrows" });
    this.textures = starTextures(gl);
    this.positions = position;
    this.offsets = offsets;
    // Cluster index per star, and velocities relative to the cluster mean
    // (internal kinematics, as in the legacy motion arrows).
    const clusterOf = new Float32Array(count);
    const rel = new Float32Array(count * 3);
    this.clusterOfStar = new Uint32Array(count);
    // Clusters with no placeable star keep their own marker in Sky view.
    this.clusterHasStars = new Uint8Array(this.clusterCount);
    for (let c = 0; c < this.clusterCount; c++) {
      const a = offsets[c], b = offsets[c + 1];
      let mx = 0, my = 0, mz = 0, n = 0;
      for (let i = a; i < b; i++) {
        if (!valid[i]) continue;
        mx += velocity[i * 3]; my += velocity[i * 3 + 1]; mz += velocity[i * 3 + 2]; n++;
      }
      this.clusterHasStars[c] = n > 0 ? 1 : 0;
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
    this.arrowData = new Float32Array(TEX_W * texH * 4);
    this.arrowTex = createDataTexture(gl, this.arrowData, TEX_W, texH, 4);
    this.hasArrows = false;
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
    this.arrowVao = createVAO(gl, [
      { loc: 1, buffer: buf.pos, size: 3, divisor: 1 },
      { loc: 2, buffer: buf.vel, size: 3, divisor: 1 },
      { loc: 3, buffer: buf.cluster, size: 1, divisor: 1 },
    ]);
    this.params = null;
    this.visible = false;
  }

  /**
   * Motion-arrow scales per cluster: `fill(data)` writes length and width
   * scales at [c * 4], [c * 4 + 1] (0 length: no arrows for that cluster).
   */
  updateArrows(fill) {
    const gl = this.gl;
    this.arrowData.fill(0);
    fill(this.arrowData);
    let any = false;
    for (let c = 0; c < this.clusterCount && !any; c++) any = this.arrowData[c * 4] > 0;
    this.hasArrows = any;
    gl.bindTexture(gl.TEXTURE_2D, this.arrowTex);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, TEX_W, this.texH, gl.RGBA, gl.FLOAT, this.arrowData);
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
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.enable(gl.DEPTH_TEST);
    gl.depthMask(false);
    if (this.hasArrows) {
      // Arrows under the stars: 9 vertices (a shaft quad and a head) each.
      const ap = this.arrowProgram.use();
      ap.m4("uViewProj", frame.camera.viewProj);
      ap.tex("uClusterState", this.stateTex).tex("uClusterColor", this.colorTex).tex("uArrowStyle", this.arrowTex);
      ap.f("uTime", p.time).f("uInternal", p.internalMotion ? 1 : 0).f("uReveal", p.reveal);
      gl.bindVertexArray(this.arrowVao);
      gl.drawArraysInstanced(gl.TRIANGLES, 0, 9, this.count);
      gl.bindVertexArray(null);
    }
    const prog = this.program.use();
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
    gl.deleteVertexArray(this.arrowVao);
    for (const b of Object.values(this.buffers)) gl.deleteBuffer(b);
    gl.deleteTexture(this.stateTex);
    gl.deleteTexture(this.colorTex);
    gl.deleteTexture(this.arrowTex);
    this.program.dispose();
    this.pickProgram.dispose();
    this.arrowProgram.dispose();
  }
}
