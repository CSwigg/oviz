// Event-density volumes rebuilt in the viewer (classic "procedural KDE").
//
// A volume spec with `procedural.kind === "event-kde"` ships its events
// (float32 t, x, y, z rows sorted by t) instead of voxels. At each time
// sample the density is rebuilt exactly like the classic runtime did:
//   * the events in the trailing window (t − window, t] are deposited
//     cloud-in-cell on the grid (voxel centres at bounds + (i + ½)·d),
//   * smoothed by a core + halo Gaussian mixture (separable passes limited
//     to the deposited box plus the kernel margin),
//   * divided by `dataMax` and clamped to [0, 1].
// Time snaps DOWN to the sample grid, so a frame never shows events from
// its future. The result is float (the low tail lies below one uint8 step).

/** The KDE sample time for `time`: snapped down to `sampleMyr` (0 = exact). */
export function kdeSampleTime(proc, time) {
  const sample = Math.max(0, Number(proc.sampleMyr) || 0);
  const t = Number.isFinite(time) ? time : 0;
  if (sample <= 1e-6) return t;
  return Math.floor((t + 1e-6) / sample) * sample;
}

/** The trailing window a volume state asks for, else the spec's default. */
export function kdeWindow(proc, vs) {
  return Math.max(1e-6, Number(vs?.kdeTimeWindowMyr) || Number(proc.windowMyr) || 10);
}

/** Cache identity of one built density. */
export function kdeId(time, windowMyr) {
  return `kde:${time.toFixed(3)}:${windowMyr.toFixed(3)}`;
}

/** First row (×4) whose time is greater than `t` (events sorted by time). */
export function upperBound(events, t) {
  let lo = 0, hi = Math.floor(events.length / 4);
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (events[mid * 4] <= t) lo = mid + 1;
    else hi = mid;
  }
  return lo * 4;
}

function gaussianKernel(sigmaVox) {
  const sigma = Math.max(0, Number(sigmaVox) || 0);
  if (sigma <= 1e-6) return { values: new Float32Array([1]), radius: 0 };
  const radius = Math.max(1, Math.ceil(4 * sigma));
  const values = new Float32Array(2 * radius + 1);
  let total = 0;
  for (let o = -radius; o <= radius; o++) {
    const w = Math.exp((-0.5 * o * o) / (sigma * sigma));
    values[o + radius] = w;
    total += w;
  }
  for (let i = 0; i < values.length; i++) values[i] /= Math.max(total, 1e-12);
  return { values, radius };
}

let scratch64 = null;

/**
 * One separable pass along `axis` inside `box` (inclusive voxel bounds).
 * Accumulates in float64 offset by offset and stores float32 once, which
 * reproduces the classic kernel's rounding bit for bit.
 */
function blurAxis(input, nx, ny, nz, kernel, radius, axis, box) {
  if (radius <= 0) return input;
  const out = new Float32Array(input.length);
  const plane = nx * ny;
  const { x0, x1, y0, y1, z0, z1 } = box;
  if (!scratch64 || scratch64.length < input.length) scratch64 = new Float64Array(input.length);
  const acc = scratch64;
  for (let z = z0; z <= z1; z++) for (let y = y0; y <= y1; y++) {
    const row = z * plane + y * nx;
    acc.fill(0, row + x0, row + x1 + 1);
  }
  const stride = axis === 0 ? 1 : axis === 1 ? nx : plane;
  const limit = axis === 0 ? nx : axis === 1 ? ny : nz;
  const lo = axis === 0 ? x0 : axis === 1 ? y0 : z0;
  const hi = axis === 0 ? x1 : axis === 1 ? y1 : z1;
  for (let o = -radius; o <= radius; o++) {
    const w = kernel[o + radius];
    const shift = o * stride;
    const a = Math.max(lo, -o), b = Math.min(hi, limit - 1 - o);
    if (a > b) continue;
    const za = axis === 2 ? a : z0, zb = axis === 2 ? b : z1;
    const ya = axis === 1 ? a : y0, yb = axis === 1 ? b : y1;
    const xa = axis === 0 ? a : x0, xb = axis === 0 ? b : x1;
    for (let z = za; z <= zb; z++) {
      const zBase = z * plane;
      for (let y = ya; y <= yb; y++) {
        const row = zBase + y * nx;
        for (let x = xa; x <= xb; x++) acc[row + x] += w * input[row + x + shift];
      }
    }
  }
  for (let z = z0; z <= z1; z++) for (let y = y0; y <= y1; y++) {
    const row = z * plane + y * nx;
    for (let x = x0; x <= x1; x++) out[row + x] = acc[row + x];
  }
  return out;
}

function gaussianBlur(input, nx, ny, nz, sigmaVox, box) {
  const k = gaussianKernel(sigmaVox);
  let out = blurAxis(input, nx, ny, nz, k.values, k.radius, 0, box);
  out = blurAxis(out, nx, ny, nz, k.values, k.radius, 1, box);
  return blurAxis(out, nx, ny, nz, k.values, k.radius, 2, box);
}

/**
 * Build the normalised density (Float32Array, z-y-x) at `time` for a
 * trailing window of `windowMyr`. `dims` = [nx, ny, nz], `bounds` =
 * [[xmin, ymin, zmin], [xmax, ymax, zmax]]. Also returns the deposited box.
 */
export function buildEventKde(events, proc, dims, bounds, time, windowMyr) {
  const [nx, ny, nz] = dims;
  const [lo, hi] = bounds;
  const dx = (hi[0] - lo[0]) / nx, dy = (hi[1] - lo[1]) / ny, dz = (hi[2] - lo[2]) / nz;
  const plane = nx * ny;
  const scale = Math.max(1e-12, Number(proc.depositionScale) || 1);
  const density = new Float32Array(nx * ny * nz);
  let mx0 = nx, mx1 = -1, my0 = ny, my1 = -1, mz0 = nz, mz1 = -1;
  const start = upperBound(events, time - windowMyr);
  const end = upperBound(events, time + 1e-6);
  for (let e = start; e + 3 < end; e += 4) {
    const et = events[e];
    if (et <= time - windowMyr || et > time + 1e-6) continue;
    const fx = Math.min(nx - 1, Math.max(0, (events[e + 1] - lo[0]) / dx - 0.5));
    const fy = Math.min(ny - 1, Math.max(0, (events[e + 2] - lo[1]) / dy - 0.5));
    const fz = Math.min(nz - 1, Math.max(0, (events[e + 3] - lo[2]) / dz - 0.5));
    const x0 = Math.floor(fx), y0 = Math.floor(fy), z0 = Math.floor(fz);
    const tx = fx - x0, ty = fy - y0, tz = fz - z0;
    for (let oz = 0; oz <= 1; oz++) {
      const z = z0 + oz;
      if (z >= nz) continue;
      const wz = oz ? tz : 1 - tz;
      for (let oy = 0; oy <= 1; oy++) {
        const y = y0 + oy;
        if (y >= ny) continue;
        const wy = oy ? ty : 1 - ty;
        for (let ox = 0; ox <= 1; ox++) {
          const x = x0 + ox;
          if (x >= nx) continue;
          const wx = ox ? tx : 1 - tx;
          density[z * plane + y * nx + x] += wx * wy * wz * scale;
          if (x < mx0) mx0 = x;
          if (x > mx1) mx1 = x;
          if (y < my0) my0 = y;
          if (y > my1) my1 = y;
          if (z < mz0) mz0 = z;
          if (z > mz1) mz1 = z;
        }
      }
    }
  }
  const out = new Float32Array(density.length);
  if (mx1 < mx0) return { data: out, box: null };
  const core = Math.max(0, Number(proc.coreSigmaVox) || 0);
  const halo = Math.max(core, Number(proc.haloSigmaVox) || core);
  const coreWeight = Math.min(1, Math.max(0, Number(proc.coreWeight) || 1));
  const margin = Math.max(core <= 1e-6 ? 0 : Math.max(1, Math.ceil(4 * core)), halo <= 1e-6 ? 0 : Math.max(1, Math.ceil(4 * halo)));
  const box = {
    x0: Math.max(0, mx0 - margin), x1: Math.min(nx - 1, mx1 + margin),
    y0: Math.max(0, my0 - margin), y1: Math.min(ny - 1, my1 + margin),
    z0: Math.max(0, mz0 - margin), z1: Math.min(nz - 1, mz1 + margin),
  };
  const coreSmoothed = gaussianBlur(density, nx, ny, nz, core, box);
  let smoothed = coreSmoothed;
  if (halo > core + 1e-6 && coreWeight < 1 - 1e-6) {
    const haloSmoothed = gaussianBlur(density, nx, ny, nz, halo, box);
    smoothed = new Float32Array(coreSmoothed.length);
    const hw = 1 - coreWeight;
    for (let z = box.z0; z <= box.z1; z++) for (let y = box.y0; y <= box.y1; y++) {
      let i = z * plane + y * nx + box.x0;
      for (let x = box.x0; x <= box.x1; x++, i++) smoothed[i] = coreWeight * coreSmoothed[i] + hw * haloSmoothed[i];
    }
  }
  const inv = 1 / Math.max(1e-6, Number(proc.dataMax) || 1);
  for (let z = box.z0; z <= box.z1; z++) for (let y = box.y0; y <= box.y1; y++) {
    let i = z * plane + y * nx + box.x0;
    for (let x = box.x0; x <= box.x1; x++, i++) out[i] = Math.min(1, Math.max(0, smoothed[i] || 0) * inv);
  }
  return { data: out, box };
}

/** The sample times a timeline spans (snapped down, unique, ascending). */
export function kdeSampleTimes(proc, times) {
  const set = new Set();
  for (const t of times) set.add(kdeSampleTime(proc, t).toFixed(6));
  return [...set].map(Number).sort((a, b) => a - b);
}
