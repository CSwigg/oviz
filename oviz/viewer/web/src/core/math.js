// Small, allocation-conscious linear algebra for the viewer.
// Matrices are column-major Float32Array(16), matching WebGL.

export const DEG = Math.PI / 180;

export function clamp(v, lo, hi) {
  return v < lo ? lo : v > hi ? hi : v;
}

export function lerp(a, b, t) {
  return a + (b - a) * t;
}

export function smoothstep(e0, e1, x) {
  const t = clamp((x - e0) / (e1 - e0 || 1e-9), 0, 1);
  return t * t * (3 - 2 * t);
}

export function easeInOutCubic(t) {
  return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
}

export function easeOutCubic(t) {
  return 1 - Math.pow(1 - t, 3);
}

/** Shortest signed angular difference b - a in radians. */
export function angleDelta(a, b) {
  let d = (b - a) % (2 * Math.PI);
  if (d > Math.PI) d -= 2 * Math.PI;
  if (d < -Math.PI) d += 2 * Math.PI;
  return d;
}

// ---- vec3 (plain arrays) -------------------------------------------------

export function v3(x = 0, y = 0, z = 0) {
  return [x, y, z];
}

export function v3add(a, b, out = [0, 0, 0]) {
  out[0] = a[0] + b[0]; out[1] = a[1] + b[1]; out[2] = a[2] + b[2];
  return out;
}

export function v3sub(a, b, out = [0, 0, 0]) {
  out[0] = a[0] - b[0]; out[1] = a[1] - b[1]; out[2] = a[2] - b[2];
  return out;
}

export function v3scale(a, s, out = [0, 0, 0]) {
  out[0] = a[0] * s; out[1] = a[1] * s; out[2] = a[2] * s;
  return out;
}

export function v3dot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

export function v3cross(a, b, out = [0, 0, 0]) {
  const x = a[1] * b[2] - a[2] * b[1];
  const y = a[2] * b[0] - a[0] * b[2];
  const z = a[0] * b[1] - a[1] * b[0];
  out[0] = x; out[1] = y; out[2] = z;
  return out;
}

export function v3len(a) {
  return Math.hypot(a[0], a[1], a[2]);
}

export function v3dist(a, b) {
  return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

export function v3norm(a, out = [0, 0, 0]) {
  const l = v3len(a) || 1;
  out[0] = a[0] / l; out[1] = a[1] / l; out[2] = a[2] / l;
  return out;
}

export function v3lerp(a, b, t, out = [0, 0, 0]) {
  out[0] = a[0] + (b[0] - a[0]) * t;
  out[1] = a[1] + (b[1] - a[1]) * t;
  out[2] = a[2] + (b[2] - a[2]) * t;
  return out;
}

// ---- mat4 ---------------------------------------------------------------

export function m4() {
  const m = new Float32Array(16);
  m[0] = m[5] = m[10] = m[15] = 1;
  return m;
}

export function m4mul(a, b, out = new Float32Array(16)) {
  const a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3];
  const a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7];
  const a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11];
  const a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15];
  for (let i = 0; i < 4; i++) {
    const b0 = b[i * 4], b1 = b[i * 4 + 1], b2 = b[i * 4 + 2], b3 = b[i * 4 + 3];
    out[i * 4] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    out[i * 4 + 1] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    out[i * 4 + 2] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    out[i * 4 + 3] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;
  }
  return out;
}

export function m4perspective(fovyRad, aspect, near, far, out = new Float32Array(16)) {
  const f = 1 / Math.tan(fovyRad / 2);
  const nf = 1 / (near - far);
  out.fill(0);
  out[0] = f / aspect;
  out[5] = f;
  out[10] = (far + near) * nf;
  out[11] = -1;
  out[14] = 2 * far * near * nf;
  return out;
}

export function m4lookAt(eye, center, up, out = new Float32Array(16)) {
  let zx = eye[0] - center[0], zy = eye[1] - center[1], zz = eye[2] - center[2];
  let l = Math.hypot(zx, zy, zz) || 1;
  zx /= l; zy /= l; zz /= l;
  let xx = up[1] * zz - up[2] * zy;
  let xy = up[2] * zx - up[0] * zz;
  let xz = up[0] * zy - up[1] * zx;
  l = Math.hypot(xx, xy, xz);
  if (l < 1e-9) {
    // Looking straight along `up`: pick any perpendicular axis.
    xx = 1; xy = 0; xz = 0;
    if (Math.abs(zx) > 0.9) { xx = 0; xy = 1; }
    const d = xx * zx + xy * zy + xz * zz;
    xx -= d * zx; xy -= d * zy; xz -= d * zz;
    l = Math.hypot(xx, xy, xz);
  }
  xx /= l; xy /= l; xz /= l;
  const yx = zy * xz - zz * xy;
  const yy = zz * xx - zx * xz;
  const yz = zx * xy - zy * xx;
  out[0] = xx; out[1] = yx; out[2] = zx; out[3] = 0;
  out[4] = xy; out[5] = yy; out[6] = zy; out[7] = 0;
  out[8] = xz; out[9] = yz; out[10] = zz; out[11] = 0;
  out[12] = -(xx * eye[0] + xy * eye[1] + xz * eye[2]);
  out[13] = -(yx * eye[0] + yy * eye[1] + yz * eye[2]);
  out[14] = -(zx * eye[0] + zy * eye[1] + zz * eye[2]);
  out[15] = 1;
  return out;
}

export function m4invert(a, out = new Float32Array(16)) {
  const a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3];
  const a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7];
  const a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11];
  const a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15];
  const b00 = a00 * a11 - a01 * a10, b01 = a00 * a12 - a02 * a10;
  const b02 = a00 * a13 - a03 * a10, b03 = a01 * a12 - a02 * a11;
  const b04 = a01 * a13 - a03 * a11, b05 = a02 * a13 - a03 * a12;
  const b06 = a20 * a31 - a21 * a30, b07 = a20 * a32 - a22 * a30;
  const b08 = a20 * a33 - a23 * a30, b09 = a21 * a32 - a22 * a31;
  const b10 = a21 * a33 - a23 * a31, b11 = a22 * a33 - a23 * a32;
  let det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
  if (!det) return null;
  det = 1 / det;
  out[0] = (a11 * b11 - a12 * b10 + a13 * b09) * det;
  out[1] = (a02 * b10 - a01 * b11 - a03 * b09) * det;
  out[2] = (a31 * b05 - a32 * b04 + a33 * b03) * det;
  out[3] = (a22 * b04 - a21 * b05 - a23 * b03) * det;
  out[4] = (a12 * b08 - a10 * b11 - a13 * b07) * det;
  out[5] = (a00 * b11 - a02 * b08 + a03 * b07) * det;
  out[6] = (a32 * b02 - a30 * b05 - a33 * b01) * det;
  out[7] = (a20 * b05 - a22 * b02 + a23 * b01) * det;
  out[8] = (a10 * b10 - a11 * b08 + a13 * b06) * det;
  out[9] = (a01 * b08 - a00 * b10 - a03 * b06) * det;
  out[10] = (a30 * b04 - a31 * b02 + a33 * b00) * det;
  out[11] = (a21 * b02 - a20 * b04 - a23 * b00) * det;
  out[12] = (a11 * b07 - a10 * b09 - a12 * b06) * det;
  out[13] = (a00 * b09 - a01 * b07 + a02 * b06) * det;
  out[14] = (a31 * b01 - a30 * b03 - a32 * b00) * det;
  out[15] = (a20 * b03 - a21 * b01 + a22 * b00) * det;
  return out;
}

/** Transform a point by a 4x4 matrix with perspective divide. */
export function m4project(m, p, out = [0, 0, 0, 0]) {
  const x = p[0], y = p[1], z = p[2];
  const w = m[3] * x + m[7] * y + m[11] * z + m[15];
  out[0] = (m[0] * x + m[4] * y + m[8] * z + m[12]) / w;
  out[1] = (m[1] * x + m[5] * y + m[9] * z + m[13]) / w;
  out[2] = (m[2] * x + m[6] * y + m[10] * z + m[14]) / w;
  out[3] = w;
  return out;
}

// ---- spherical / galactic helpers --------------------------------------

/** Unit direction for galactic longitude/latitude (degrees). x→GC, z→NGP. */
export function lbToDir(lDeg, bDeg, out = [0, 0, 0]) {
  const l = lDeg * DEG, b = bDeg * DEG;
  const cb = Math.cos(b);
  out[0] = cb * Math.cos(l);
  out[1] = cb * Math.sin(l);
  out[2] = Math.sin(b);
  return out;
}

export function dirToLb(d) {
  const r = Math.hypot(d[0], d[1], d[2]) || 1;
  let l = Math.atan2(d[1], d[0]) / DEG;
  if (l < 0) l += 360;
  const b = Math.asin(clamp(d[2] / r, -1, 1)) / DEG;
  return [l, b];
}

// ICRS → Galactic rotation T (row-major), Hipparcos/J2000 definition:
// gal = T · icrs and icrs = Tᵀ · gal.
const ICRS_TO_GAL = [
  -0.0548755604162154, -0.8734370902348850, -0.4838350155487132,
  0.4941094278755837, -0.4448296299600112, 0.7469822444972189,
  -0.8676661490190047, -0.1980763734312015, 0.4559837761750669,
];

export function galToIcrs(lDeg, bDeg) {
  const g = lbToDir(lDeg, bDeg);
  const T = ICRS_TO_GAL;
  const x = T[0] * g[0] + T[3] * g[1] + T[6] * g[2];
  const y = T[1] * g[0] + T[4] * g[1] + T[7] * g[2];
  const z = T[2] * g[0] + T[5] * g[1] + T[8] * g[2];
  let ra = Math.atan2(y, x) / DEG;
  if (ra < 0) ra += 360;
  return [ra, Math.asin(clamp(z, -1, 1)) / DEG];
}

export function icrsToGal(raDeg, decDeg) {
  const e = lbToDir(raDeg, decDeg);
  const T = ICRS_TO_GAL;
  return dirToLb([
    T[0] * e[0] + T[1] * e[1] + T[2] * e[2],
    T[3] * e[0] + T[4] * e[1] + T[5] * e[2],
    T[6] * e[0] + T[7] * e[1] + T[8] * e[2],
  ]);
}

// ---- formatting ---------------------------------------------------------

export function formatDistance(pc) {
  const a = Math.abs(pc);
  if (a >= 1000) return `${(pc / 1000).toFixed(a >= 10000 ? 1 : 2)} kpc`;
  if (a >= 100) return `${Math.round(pc)} pc`;
  if (a >= 10) return `${pc.toFixed(1)} pc`;
  return `${pc.toFixed(2)} pc`;
}

export function formatAngle(deg) {
  const a = Math.abs(deg);
  if (a >= 1) return `${deg.toFixed(a >= 10 ? 1 : 2)}°`;
  if (a * 60 >= 1) return `${(deg * 60).toFixed(1)}′`;
  return `${(deg * 3600).toFixed(1)}″`;
}

/** A "nice" number (1, 2, 5 × 10^n) not larger than `x`. */
export function niceFloor(x) {
  if (!(x > 0)) return 0;
  const p = Math.pow(10, Math.floor(Math.log10(x)));
  const m = x / p;
  return (m >= 5 ? 5 : m >= 2 ? 2 : 1) * p;
}
