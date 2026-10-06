// Annotations drawn into the figure: curves and arrows as anti-aliased
// screen-space ribbons with round joins, bubbles and shells as two-sided
// ellipsoids shaded by how obliquely they are seen, and text as textures
// drawn on top. Everything is in the WebGL frame, so screenshots and videos
// carry it. The annotate plugin hands over a resolved scene (world
// positions at the current time); this layer only uploads and draws it.

import { createProgram, createBuffer, createVAO, createTexture2D } from "../engine/gl.js";

// Per segment: a.xyz, b.xyz, arc0, arc1, rgba, width (CSS px), dash, gap, 0.
export const ANNOT_SEG_FLOATS = 16;
// Per arrowhead: tip.xyz, from.xyz, rgba, size (CSS px), width (CSS px).
export const ANNOT_ARROW_FLOATS = 12;
// Per ellipsoid: centre.xyz, radii.xyz, quaternion, rgba, style, highlight.
export const ANNOT_SPHERE_FLOATS = 16;
// Per mesh vertex (box faces, a curve's selection tube): xyz, normal, rgba, style, 0.
export const ANNOT_MESH_FLOATS = 12;

const SEG_VERT = `
layout(location = 0) in vec2 aCorner;   // x: 0 start / 1 end, y: -1..1 side
layout(location = 1) in vec3 aA;
layout(location = 2) in vec3 aB;
layout(location = 3) in vec2 aArc;
layout(location = 4) in vec4 aColor;
layout(location = 5) in vec4 aParams;   // width (CSS px), dash, gap, unused
uniform mat4 uViewProj;
uniform vec2 uViewport;                 // device pixels
uniform float uDpr;
out vec2 vP;
out vec2 vSa;
out vec2 vSb;
out float vArc;
out float vHalf;
out vec4 vColor;
out vec2 vDash;
void main() {
  vec4 ca = uViewProj * vec4(aA, 1.0);
  vec4 cb = uViewProj * vec4(aB, 1.0);
  float near = 1e-4;
  if (ca.w < near && cb.w < near) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }
  float arcA = aArc.x, arcB = aArc.y;
  if (ca.w < near) { float k = (near - ca.w) / (cb.w - ca.w); ca = mix(ca, cb, k); arcA = mix(arcA, arcB, k); }
  if (cb.w < near) { float k = (near - cb.w) / (ca.w - cb.w); cb = mix(cb, ca, k); arcB = mix(arcB, arcA, k); }
  vec2 sa = ca.xy / ca.w * uViewport * 0.5;
  vec2 sb = cb.xy / cb.w * uViewport * 0.5;
  vec2 dir = sb - sa;
  float len = length(dir);
  dir = len > 1e-6 ? dir / len : vec2(1.0, 0.0);
  vec2 n = vec2(-dir.y, dir.x);
  float half_ = aParams.x * uDpr * 0.5;
  float pad = half_ + 1.5;
  bool atEnd = aCorner.x > 0.5;
  vec2 sp = (atEnd ? sb + dir * pad : sa - dir * pad) + n * aCorner.y * pad;
  float z = atEnd ? cb.z / cb.w : ca.z / ca.w;
  // Screen-space output (w = 1) keeps the varyings linear on screen, so the
  // fragment's distance to the segment is exact.
  gl_Position = vec4(sp / (uViewport * 0.5), clamp(z, -1.0, 1.0), 1.0);
  vP = sp; vSa = sa; vSb = sb;
  vArc = atEnd ? arcB : arcA;
  vHalf = half_;
  vColor = aColor;
  vDash = aParams.yz;
}
`;

const SEG_FRAG = `
in vec2 vP;
in vec2 vSa;
in vec2 vSb;
in float vArc;
in float vHalf;
in vec4 vColor;
in vec2 vDash;
out vec4 outColor;
void main() {
  vec2 ab = vSb - vSa;
  float l2 = dot(ab, ab);
  float u = l2 > 1e-9 ? clamp(dot(vP - vSa, ab) / l2, 0.0, 1.0) : 0.0;
  float d = length(vP - vSa - ab * u);
  float a = clamp(vHalf + 0.5 - d, 0.0, 1.0);
  if (vHalf < 0.5) a *= vHalf * 2.0;
  if (vDash.x > 0.0) {
    float period = vDash.x + vDash.y;
    float m = mod(vArc, period);
    float e = fwidth(vArc) + 1e-6;
    a *= 1.0 - smoothstep(vDash.x - e, vDash.x + e, m) * (1.0 - smoothstep(period - e, period, m));
  }
  a *= vColor.a;
  if (a < 0.003) discard;
  outColor = vec4(vColor.rgb * a, a);
}
`;

const ARROW_VERT = `
layout(location = 0) in float aCorner;  // 0 tip, 1 left, 2 right
layout(location = 1) in vec3 aTip;
layout(location = 2) in vec3 aFrom;
layout(location = 3) in vec4 aColor;
layout(location = 4) in vec2 aSize;     // head length, line width (CSS px)
uniform mat4 uViewProj;
uniform vec2 uViewport;
uniform float uDpr;
out vec4 vColor;
void main() {
  vec4 ct = uViewProj * vec4(aTip, 1.0);
  vec4 cf = uViewProj * vec4(aFrom, 1.0);
  if (ct.w < 1e-4 || cf.w < 1e-4) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }
  vec2 st = ct.xy / ct.w * uViewport * 0.5;
  vec2 sf = cf.xy / cf.w * uViewport * 0.5;
  vec2 dir = st - sf;
  float len = length(dir);
  dir = len > 1e-6 ? dir / len : vec2(1.0, 0.0);
  vec2 n = vec2(-dir.y, dir.x);
  float L = aSize.x * uDpr;
  float W = max(L * 0.78, aSize.y * uDpr * 3.2);
  // The head starts a little behind the tip so a wide line never shows past it.
  vec2 tip = st + dir * aSize.y * uDpr * 0.5;
  vec2 base = tip - dir * L;
  vec2 p = aCorner < 0.5 ? tip : (aCorner < 1.5 ? base + n * W * 0.5 : base - n * W * 0.5);
  gl_Position = vec4(p / (uViewport * 0.5), clamp(ct.z / ct.w, -1.0, 1.0), 1.0);
  vColor = aColor;
}
`;

const ARROW_FRAG = `
in vec4 vColor;
out vec4 outColor;
void main() {
  outColor = vec4(vColor.rgb * vColor.a, vColor.a);
}
`;

const SPHERE_VERT = `
layout(location = 0) in vec3 aPos;      // unit sphere
layout(location = 1) in vec3 aCenter;
layout(location = 2) in vec3 aRadii;
layout(location = 3) in vec4 aQuat;
layout(location = 4) in vec4 aColor;
layout(location = 5) in vec2 aStyle;    // style (0 bubble, 1 shell, 2 wire), highlight
uniform mat4 uViewProj;
out vec3 vWorld;
out vec3 vNormal;
out vec3 vLocal;
out vec4 vColor;
flat out vec2 vStyle;
vec3 rotate(vec4 q, vec3 v) {
  vec3 t = 2.0 * cross(q.xyz, v);
  return v + q.w * t + cross(q.xyz, t);
}
void main() {
  vec3 world = aCenter + rotate(aQuat, aPos * aRadii);
  vWorld = world;
  vNormal = rotate(aQuat, aPos / max(aRadii, vec3(1e-6)));
  vLocal = aPos;
  vColor = aColor;
  vStyle = aStyle;
  gl_Position = uViewProj * vec4(world, 1.0);
}
`;

const SPHERE_FRAG = `
in vec3 vWorld;
in vec3 vNormal;
in vec3 vLocal;
in vec4 vColor;
flat in vec2 vStyle;
uniform vec3 uEye;
uniform float uBack;   // 1 while drawing the far side
out vec4 outColor;
float gridLine(float v, float step) {
  float f = abs(fract(v / step + 0.5) - 0.5) * step;
  float w = fwidth(v) * 1.1;
  return 1.0 - smoothstep(w * 0.5, w * 1.5, f);
}
void main() {
  vec3 N = normalize(vNormal);
  vec3 V = normalize(uEye - vWorld);
  float facing = abs(dot(N, V));
  float rim = 1.0 - facing;
  float a;
  int style = int(vStyle.x + 0.5);
  if (style == 0) {
    // A bubble: a soft body that brightens towards its edge.
    a = 0.16 + 0.84 * pow(rim, 2.2);
  } else if (style == 1) {
    // A thin shell: bright only where the line of sight runs along it.
    a = 0.035 + 0.965 * pow(rim, 3.2);
  } else {
    float lat = asin(clamp(vLocal.z, -1.0, 1.0));
    float lon = atan(vLocal.y, vLocal.x);
    float step = radians(22.5);
    float g = max(gridLine(lat, step), gridLine(lon, step) * smoothstep(0.0, 0.25, 1.0 - abs(vLocal.z)));
    a = g * (0.55 + 0.45 * rim) + 0.04 * pow(rim, 2.0);
  }
  a *= vColor.a * (uBack > 0.5 ? 0.5 : 1.0);
  if (vStyle.y > 0.5) a = max(a, 0.12 + 0.5 * pow(rim, 6.0));
  if (a < 0.002) discard;
  outColor = vec4(vColor.rgb * a, a);
}
`;

const MESH_VERT = `
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec4 aColor;
layout(location = 3) in vec2 aStyle;    // 0 box face, 1 selection tube
uniform mat4 uViewProj;
out vec3 vWorld;
out vec3 vNormal;
out vec4 vColor;
flat out float vStyle;
void main() {
  vWorld = aPos;
  vNormal = aNormal;
  vColor = aColor;
  vStyle = aStyle.x;
  gl_Position = uViewProj * vec4(aPos, 1.0);
}
`;

const MESH_FRAG = `
in vec3 vWorld;
in vec3 vNormal;
in vec4 vColor;
flat in float vStyle;
uniform vec3 uEye;
out vec4 outColor;
void main() {
  vec3 N = normalize(vNormal);
  vec3 V = normalize(uEye - vWorld);
  float rim = 1.0 - abs(dot(N, V));
  // A box face is a faint pane; a selection tube a glassy sleeve.
  float a = vStyle < 0.5 ? 0.07 + 0.11 * rim : 0.025 + 0.5 * pow(rim, 3.0);
  a *= vColor.a;
  if (a < 0.002) discard;
  outColor = vec4(vColor.rgb * a, a);
}
`;

const TEXT_VERT = `
layout(location = 0) in vec2 aCorner;
uniform vec2 uViewport;     // device pixels
uniform vec2 uOrigin;       // top-left in device pixels (from the top-left of the canvas)
uniform vec2 uSize;         // device pixels
out vec2 vUv;
void main() {
  vec2 p = uOrigin + aCorner * uSize;
  vec2 ndc = vec2(p.x / uViewport.x * 2.0 - 1.0, 1.0 - p.y / uViewport.y * 2.0);
  gl_Position = vec4(ndc, 0.0, 1.0);
  vUv = aCorner;
}
`;

const TEXT_FRAG = `
in vec2 vUv;
uniform sampler2D uTex;
uniform float uOpacity;
out vec4 outColor;
void main() {
  vec4 c = texture(uTex, vUv);
  outColor = c * uOpacity;
}
`;

const SEG_CORNERS = new Float32Array([0, -1, 1, -1, 0, 1, 1, 1]);
const QUAD = new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]);

/** Vertices and triangles of an icosphere (unit radius), `level` subdivisions. */
function unitIcosphere(level = 3) {
  const t = (1 + Math.sqrt(5)) / 2;
  let verts = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]];
  let faces = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
    [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]];
  const norm = (v) => { const l = Math.hypot(v[0], v[1], v[2]); return [v[0] / l, v[1] / l, v[2] / l]; };
  verts = verts.map(norm);
  for (let s = 0; s < level; s++) {
    const cache = new Map();
    const mid = (a, b) => {
      const key = a < b ? `${a},${b}` : `${b},${a}`;
      let i = cache.get(key);
      if (i == null) {
        const p = verts[a], q = verts[b];
        i = verts.length;
        verts.push(norm([(p[0] + q[0]) / 2, (p[1] + q[1]) / 2, (p[2] + q[2]) / 2]));
        cache.set(key, i);
      }
      return i;
    };
    const next = [];
    for (const [a, b, c] of faces) {
      const ab = mid(a, b), bc = mid(b, c), ca = mid(c, a);
      next.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]);
    }
    faces = next;
  }
  return { positions: new Float32Array(verts.flat()), indices: new Uint16Array(faces.flat()) };
}

export class AnnotationsLayer {
  constructor(gl) {
    this.gl = gl;
    // After the dust and lines, before the stars: points draw over shells.
    this.order = 25;
    this.segProgram = createProgram(gl, SEG_VERT, SEG_FRAG, { label: "annot-curves" });
    this.arrowProgram = createProgram(gl, ARROW_VERT, ARROW_FRAG, { label: "annot-arrows" });
    this.sphereProgram = createProgram(gl, SPHERE_VERT, SPHERE_FRAG, { label: "annot-shells" });
    this.textProgram = createProgram(gl, TEXT_VERT, TEXT_FRAG, { label: "annot-text" });
    this.meshProgram = createProgram(gl, MESH_VERT, MESH_FRAG, { label: "annot-mesh" });
    this.cornerBuffer = createBuffer(gl, SEG_CORNERS);
    this.quadBuffer = createBuffer(gl, QUAD);
    this.arrowCornerBuffer = createBuffer(gl, new Float32Array([0, 1, 2]));
    const ico = unitIcosphere(4);
    this.icoBuffer = createBuffer(gl, ico.positions);
    this.icoIndex = gl.createBuffer();
    // Index bindings belong to the bound VAO: make sure none is.
    gl.bindVertexArray(null);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.icoIndex);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, ico.indices, gl.STATIC_DRAW);
    this.icoCount = ico.indices.length;
    this.segBuffer = gl.createBuffer();
    this.arrowBuffer = gl.createBuffer();
    this.sphereBuffer = gl.createBuffer();
    this.meshBuffer = gl.createBuffer();
    const meshStride = ANNOT_MESH_FLOATS * 4;
    this.meshVao = this._vao([
      { loc: 0, buffer: this.meshBuffer, size: 3, stride: meshStride, offset: 0 },
      { loc: 1, buffer: this.meshBuffer, size: 3, stride: meshStride, offset: 12 },
      { loc: 2, buffer: this.meshBuffer, size: 4, stride: meshStride, offset: 24 },
      { loc: 3, buffer: this.meshBuffer, size: 2, stride: meshStride, offset: 40 },
    ]);
    this.segVao = this._vao([
      { loc: 0, buffer: this.cornerBuffer, size: 2 },
      ...this._instanced(this.segBuffer, ANNOT_SEG_FLOATS, [[1, 3], [2, 3], [3, 2], [4, 4], [5, 4]]),
    ]);
    this.arrowVao = this._vao([
      { loc: 0, buffer: this.arrowCornerBuffer, size: 1 },
      ...this._instanced(this.arrowBuffer, ANNOT_ARROW_FLOATS, [[1, 3], [2, 3], [3, 4], [4, 2]]),
    ]);
    this.sphereVao = createVAO(gl, [
      { loc: 0, buffer: this.icoBuffer, size: 3 },
      ...this._instanced(this.sphereBuffer, ANNOT_SPHERE_FLOATS, [[1, 3], [2, 3], [3, 4], [4, 4], [5, 2]]),
    ], this.icoIndex);
    this.textVao = this._vao([{ loc: 0, buffer: this.quadBuffer, size: 2 }]);
    this.counts = { seg: 0, arrow: 0, sphere: 0, mesh: 0 };
    this.texts = [];
    this.textures = new Map(); // style key → {tex, w, h, used}
    this._version = -1;
    this.scene = null;
  }

  _vao(attrs) {
    return createVAO(this.gl, attrs);
  }

  /** Interleaved per-instance attributes: [[loc, size], …] packed in order. */
  _instanced(buffer, floats, layout) {
    let offset = 0;
    return layout.map(([loc, size]) => {
      const a = { loc, buffer, size, stride: floats * 4, offset: offset * 4, divisor: 1 };
      offset += size;
      return a;
    });
  }

  /** A resolved scene from the plugin: typed arrays plus text records. */
  setScene(scene) {
    this.scene = scene;
  }

  _upload() {
    const s = this.scene;
    if (!s || s.version === this._version) return;
    this._version = s.version;
    const gl = this.gl;
    const put = (buf, data) => {
      gl.bindBuffer(gl.ARRAY_BUFFER, buf);
      gl.bufferData(gl.ARRAY_BUFFER, data && data.length ? data : new Float32Array(4), gl.DYNAMIC_DRAW);
    };
    put(this.segBuffer, s.segments);
    put(this.arrowBuffer, s.arrows);
    put(this.sphereBuffer, s.spheres);
    put(this.meshBuffer, s.mesh);
    this.counts = { seg: s.segmentCount || 0, arrow: s.arrowCount || 0, sphere: s.sphereCount || 0, mesh: s.meshCount || 0 };
  }

  /** A texture of `t`'s text at the current pixel ratio, cached by its look. */
  _textTexture(t, dpr) {
    const key = `${dpr}|${t.size}|${t.weight}|${t.color}|${t.bg ? 1 : 0}|${t.text}`;
    let entry = this.textures.get(key);
    if (!entry) {
      const r = rasterizeText(t, dpr);
      if (!r) return null;
      const gl = this.gl;
      gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, true);
      const tex = createTexture2D(gl, { width: r.canvas.width, height: r.canvas.height, internalFormat: gl.RGBA, format: gl.RGBA, type: gl.UNSIGNED_BYTE, data: null });
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, r.canvas);
      gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, false);
      entry = { tex, w: r.canvas.width, h: r.canvas.height, cssW: r.canvas.width / dpr, cssH: r.canvas.height / dpr };
      this.textures.set(key, entry);
    }
    entry.used = this._frameId;
    return entry;
  }

  /** On-screen box (CSS px) of each text item, for picking and handles. */
  textBoxes() {
    return this._boxes || new Map();
  }

  draw(frame) {
    const s = this.scene;
    if (!s) return;
    this._upload();
    const gl = this.gl;
    const cam = frame.camera;
    this._frameId = (this._frameId || 0) + 1;
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.enable(gl.DEPTH_TEST);
    gl.depthMask(false);
    // Ellipsoids: far side first, then the near side.
    if (this.counts.sphere) {
      const prog = this.sphereProgram.use();
      prog.m4("uViewProj", cam.viewProj).v3("uEye", cam.eye[0], cam.eye[1], cam.eye[2]);
      gl.bindVertexArray(this.sphereVao);
      gl.enable(gl.CULL_FACE);
      for (const back of [1, 0]) {
        gl.cullFace(back ? gl.FRONT : gl.BACK);
        prog.f("uBack", back);
        gl.drawElementsInstanced(gl.TRIANGLES, this.icoCount, gl.UNSIGNED_SHORT, 0, this.counts.sphere);
      }
      gl.disable(gl.CULL_FACE);
    }
    // Box faces and selection tubes: faint, so their order hardly shows.
    if (this.counts.mesh) {
      const prog = this.meshProgram.use();
      prog.m4("uViewProj", cam.viewProj).v3("uEye", cam.eye[0], cam.eye[1], cam.eye[2]);
      gl.bindVertexArray(this.meshVao);
      gl.drawArrays(gl.TRIANGLES, 0, this.counts.mesh);
    }
    if (this.counts.seg) {
      const prog = this.segProgram.use();
      prog.m4("uViewProj", cam.viewProj).v2("uViewport", frame.width, frame.height).f("uDpr", frame.dpr);
      gl.bindVertexArray(this.segVao);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, this.counts.seg);
    }
    if (this.counts.arrow) {
      const prog = this.arrowProgram.use();
      prog.m4("uViewProj", cam.viewProj).v2("uViewport", frame.width, frame.height).f("uDpr", frame.dpr);
      gl.bindVertexArray(this.arrowVao);
      gl.drawArraysInstanced(gl.TRIANGLES, 0, 3, this.counts.arrow);
    }
    // Text sits on top of everything, on whole device pixels.
    const boxes = new Map();
    if (s.texts?.length) {
      gl.disable(gl.DEPTH_TEST);
      const prog = this.textProgram.use();
      prog.v2("uViewport", frame.width, frame.height);
      gl.bindVertexArray(this.textVao);
      const scr = [0, 0, 0];
      for (const t of s.texts) {
        const p = cam.project(t.pos, scr);
        if (!p || t.opacity <= 0.003) continue;
        const entry = this._textTexture(t, frame.dpr);
        if (!entry) continue;
        const cx = p[0] + (t.dx || 0), cy = p[1] + (t.dy || 0);
        const x0 = Math.round((cx - entry.cssW / 2) * frame.dpr), y0 = Math.round((cy - entry.cssH / 2) * frame.dpr);
        boxes.set(t.id, { x: x0 / frame.dpr, y: y0 / frame.dpr, w: entry.cssW, h: entry.cssH, anchor: [p[0], p[1]] });
        prog.v2("uOrigin", x0, y0).v2("uSize", entry.w, entry.h).f("uOpacity", t.opacity).tex("uTex", entry.tex);
        gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
      }
      gl.enable(gl.DEPTH_TEST);
    }
    this._boxes = boxes;
    gl.bindVertexArray(null);
    // Let go of textures no longer drawn (edited text, another pixel ratio).
    if (this.textures.size > 64) {
      for (const [k, e] of this.textures) {
        if (e.used !== this._frameId) { gl.deleteTexture(e.tex); this.textures.delete(k); }
      }
    }
  }

  dispose() {
    const gl = this.gl;
    for (const e of this.textures.values()) gl.deleteTexture(e.tex);
    this.textures.clear();
    for (const b of [this.cornerBuffer, this.quadBuffer, this.arrowCornerBuffer, this.icoBuffer, this.icoIndex, this.segBuffer, this.arrowBuffer, this.sphereBuffer, this.meshBuffer]) gl.deleteBuffer(b);
    for (const v of [this.segVao, this.arrowVao, this.sphereVao, this.textVao, this.meshVao]) gl.deleteVertexArray(v);
    for (const p of [this.segProgram, this.arrowProgram, this.sphereProgram, this.textProgram, this.meshProgram]) p.dispose();
  }
}

/** Draw one label into a canvas at `dpr`: system type, a soft dark halo, an optional tag behind it. */
function rasterizeText(t, dpr) {
  if (typeof document === "undefined") return null;
  const lines = String(t.text || "").split("\n");
  const size = Math.max(4, t.size) * dpr;
  const font = `${t.weight >= 550 ? 600 : 400} ${size}px -apple-system, BlinkMacSystemFont, "SF Pro Text", "Inter", "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif`;
  const probe = (rasterizeText.probe ||= document.createElement("canvas").getContext("2d"));
  probe.font = font;
  const lineH = Math.ceil(size * 1.22);
  const widths = lines.map((l) => probe.measureText(l || " ").width);
  const padX = Math.ceil((t.bg ? 0.55 : 0.3) * size), padY = Math.ceil((t.bg ? 0.32 : 0.25) * size);
  const W = Math.ceil(Math.max(...widths) + 2 * padX), H = Math.ceil(lineH * lines.length + 2 * padY);
  const canvas = document.createElement("canvas");
  canvas.width = Math.max(2, W);
  canvas.height = Math.max(2, H);
  const ctx = canvas.getContext("2d");
  if (t.bg) {
    const r = Math.min(H / 2, 0.55 * size);
    ctx.fillStyle = "rgba(12, 15, 22, 0.62)";
    ctx.beginPath();
    ctx.roundRect ? ctx.roundRect(0.5, 0.5, W - 1, H - 1, r) : ctx.rect(0, 0, W, H);
    ctx.fill();
  }
  ctx.font = font;
  ctx.textBaseline = "middle";
  ctx.textAlign = "center";
  if (!t.bg) {
    // A dark halo keeps light text readable over bright dust.
    ctx.shadowColor = "rgba(0, 0, 0, 0.9)";
    ctx.shadowBlur = Math.max(2, size * 0.28);
  }
  ctx.fillStyle = t.color || "#ffffff";
  lines.forEach((l, i) => ctx.fillText(l, W / 2, padY + lineH * (i + 0.5)));
  return { canvas };
}
