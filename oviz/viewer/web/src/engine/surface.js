// Ray and screen geometry for double-click re-centring: where a click lands
// on a line, a volume or an image plane (classic Oviz re-centres its orbit
// on whatever is under the pointer, not only on points).

/**
 * Entry and exit distances of a ray through a box that is rotated by
 * `angle` (radians) about the z axis through `center`, or null on a miss.
 */
export function rayBox(origin, dir, boxMin, boxMax, center = null, angle = 0) {
  let o = origin, d = dir;
  if (angle && center) {
    const c = Math.cos(-angle), s = Math.sin(-angle);
    const ox = origin[0] - center[0], oy = origin[1] - center[1];
    o = [ox * c - oy * s + center[0], ox * s + oy * c + center[1], origin[2]];
    d = [dir[0] * c - dir[1] * s, dir[0] * s + dir[1] * c, dir[2]];
  }
  let tmin = -Infinity, tmax = Infinity;
  for (let k = 0; k < 3; k++) {
    if (Math.abs(d[k]) < 1e-12) {
      if (o[k] < boxMin[k] || o[k] > boxMax[k]) return null;
      continue;
    }
    const t0 = (boxMin[k] - o[k]) / d[k], t1 = (boxMax[k] - o[k]) / d[k];
    tmin = Math.max(tmin, Math.min(t0, t1));
    tmax = Math.min(tmax, Math.max(t0, t1));
  }
  if (!(tmax > Math.max(tmin, 0))) return null;
  return { near: Math.max(tmin, 0), far: tmax, local: { origin: o, dir: d } };
}

/**
 * Where along a ray through a volume half of its (windowed) density lies:
 * the middle of what the eye sees there. `voxels` is the grid (x fastest):
 * raw uint8, or a Float32Array already in 0–1 (an event KDE's density).
 * `low`/`high` is the display window in 0–1 of the raw range. Returns the
 * distance along the ray, or null when the ray sees no data.
 */
export function volumeMedianDepth(hit, voxels, dims, boxMin, boxMax, low = 0, high = 1, steps = 256) {
  if (!hit || !voxels) return null;
  const [nx, ny, nz] = dims;
  const { origin: o, dir: d } = hit.local;
  const span = hit.far - hit.near;
  const dt = span / steps;
  const inv = 1 / Math.max(high - low, 1e-6);
  const full = voxels instanceof Float32Array ? 1 : 255;
  const w = new Float32Array(steps);
  let total = 0;
  for (let i = 0; i < steps; i++) {
    const t = hit.near + (i + 0.5) * dt;
    const x = Math.floor(((o[0] + d[0] * t - boxMin[0]) / (boxMax[0] - boxMin[0])) * nx);
    const y = Math.floor(((o[1] + d[1] * t - boxMin[1]) / (boxMax[1] - boxMin[1])) * ny);
    const z = Math.floor(((o[2] + d[2] * t - boxMin[2]) / (boxMax[2] - boxMin[2])) * nz);
    if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) continue;
    const v = (voxels[x + nx * (y + ny * z)] / full - low) * inv;
    if (v <= 0) continue;
    w[i] = Math.min(v, 1);
    total += w[i];
  }
  // Too little along this ray to call it a hit (a few faint voxels).
  if (total < steps * 0.004) return null;
  let acc = 0;
  for (let i = 0; i < steps; i++) {
    acc += w[i];
    if (acc >= total / 2) return hit.near + (i + 0.5) * dt;
  }
  return null;
}

/** Distance along a ray to a horizontal rectangle (an image plane), or null. */
export function rayPlaneRect(origin, dir, center, width, height) {
  if (Math.abs(dir[2]) < 1e-9) return null;
  const t = (center[2] - origin[2]) / dir[2];
  if (!(t > 0)) return null;
  const x = origin[0] + dir[0] * t, y = origin[1] + dir[1] * t;
  if (Math.abs(x - center[0]) > width / 2 || Math.abs(y - center[1]) > height / 2) return null;
  return t;
}

/**
 * Closest point of screen segment a→b to (px, py): `{ dist, u }` with u in
 * 0–1 along the segment.
 */
export function screenSegment(px, py, a, b) {
  const dx = b[0] - a[0], dy = b[1] - a[1];
  const len2 = dx * dx + dy * dy;
  const u = len2 > 1e-9 ? Math.min(Math.max(((px - a[0]) * dx + (py - a[1]) * dy) / len2, 0), 1) : 0;
  return { dist: Math.hypot(a[0] + dx * u - px, a[1] + dy * u - py), u };
}
