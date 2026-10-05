// Event-density volumes (classic "procedural KDE"), run with `node --test`.
import test from "node:test";
import assert from "node:assert/strict";

import { buildEventKde, kdeSampleTime, kdeSampleTimes, kdeWindow, upperBound } from "../../oviz/viewer/web/src/layers/kde.js";

const close = (a, b, eps = 1e-6) => assert.ok(Math.abs(a - b) <= eps, `${a} ≉ ${b}`);

// A 12³ grid of 1 pc voxels centred on the origin.
const DIMS = [12, 12, 12];
const BOUNDS = [[-6, -6, -6], [6, 6, 6]];
const PROC = { coreSigmaVox: 0.8, haloSigmaVox: 2.0, coreWeight: 0.7, dataMax: 1, depositionScale: 1 };

function events(rows) {
  return Float32Array.from(rows.flat());
}

/** Unbounded reference: deposit, then a dense separable blur of the whole grid. */
function reference(rows, proc, time, win) {
  const [nx, ny, nz] = DIMS;
  const grid = new Float64Array(nx * ny * nz);
  for (const [t, x, y, z] of rows) {
    if (t <= time - win || t > time + 1e-6) continue;
    const f = [x + 6 - 0.5, y + 6 - 0.5, z + 6 - 0.5].map((v, i) => Math.min(DIMS[i] - 1, Math.max(0, v)));
    const b = f.map(Math.floor);
    for (let oz = 0; oz < 2; oz++) for (let oy = 0; oy < 2; oy++) for (let ox = 0; ox < 2; ox++) {
      const X = b[0] + ox, Y = b[1] + oy, Z = b[2] + oz;
      if (X >= nx || Y >= ny || Z >= nz) continue;
      const w = (ox ? f[0] - b[0] : 1 - f[0] + b[0]) * (oy ? f[1] - b[1] : 1 - f[1] + b[1]) * (oz ? f[2] - b[2] : 1 - f[2] + b[2]);
      grid[(Z * ny + Y) * nx + X] += w;
    }
  }
  const blur = (src, sigma) => {
    const r = Math.max(1, Math.ceil(4 * sigma));
    const k = [];
    let tot = 0;
    for (let o = -r; o <= r; o++) { k.push(Math.exp(-0.5 * o * o / (sigma * sigma))); tot += k[k.length - 1]; }
    let cur = src;
    for (const axis of [0, 1, 2]) {
      const out = new Float64Array(cur.length);
      for (let z = 0; z < nz; z++) for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
        let acc = 0;
        for (let o = -r; o <= r; o++) {
          const p = [x, y, z];
          p[axis] += o;
          if (p[axis] < 0 || p[axis] >= DIMS[axis]) continue;
          acc += (k[o + r] / tot) * cur[(p[2] * ny + p[1]) * nx + p[0]];
        }
        out[(z * ny + y) * nx + x] = acc;
      }
      cur = out;
    }
    return cur;
  };
  const core = blur(grid, proc.coreSigmaVox), halo = blur(grid, proc.haloSigmaVox);
  return core.map((c, i) => Math.min(1, (proc.coreWeight * c + (1 - proc.coreWeight) * halo[i]) / proc.dataMax));
}

test("the event KDE matches a dense reference within float precision", () => {
  const rows = [[-3, 0.2, -1.4, 0.7], [-1, 2.5, 2.5, -2.1], [-0.5, -4, 3, 1], [-7, 1, 1, 1]];
  const ev = events(rows);
  const { data } = buildEventKde(ev, PROC, DIMS, BOUNDS, 0, 5);
  const ref = reference(rows, PROC, 0, 5);
  for (let i = 0; i < data.length; i++) close(data[i], ref[i], 2e-7);
});

test("the trailing window excludes old and future events", () => {
  const ev = events([[-10, 0, 0, 0], [-4, 0, 0, 0], [1, 0, 0, 0]]);
  const sum = (a) => a.reduce((s, x) => s + x, 0);
  // Only the −4 Myr event lies in (−5, 0].
  const one = sum(buildEventKde(ev, PROC, DIMS, BOUNDS, 0, 5).data);
  const only = sum(buildEventKde(events([[-4, 0, 0, 0]]), PROC, DIMS, BOUNDS, 0, 5).data);
  close(one, only, 1e-6);
  // An event exactly at t − window is outside; one exactly at t is inside.
  assert.equal(sum(buildEventKde(events([[-5, 0, 0, 0]]), PROC, DIMS, BOUNDS, 0, 5).data), 0);
  assert.ok(sum(buildEventKde(events([[0, 0, 0, 0]]), PROC, DIMS, BOUNDS, 0, 5).data) > 0);
  assert.equal(buildEventKde(events([]), PROC, DIMS, BOUNDS, 0, 5).box, null);
});

test("sample times snap down so a frame never shows its future", () => {
  const proc = { sampleMyr: 2 };
  assert.equal(kdeSampleTime(proc, 0), 0);
  assert.equal(kdeSampleTime(proc, -1), -2);
  assert.equal(kdeSampleTime(proc, -2), -2);
  assert.equal(kdeSampleTime(proc, -2.5), -4);
  assert.equal(kdeSampleTime({ sampleMyr: 0 }, -2.5), -2.5);
  assert.deepEqual(kdeSampleTimes(proc, [-5, -4, -3, -2, -1, 0]), [-6, -4, -2, 0]);
  assert.equal(kdeWindow({ windowMyr: 10 }, {}), 10);
  assert.equal(kdeWindow({ windowMyr: 10 }, { kdeTimeWindowMyr: 5 }), 5);
});

test("upperBound finds the first event after a time", () => {
  const ev = events([[-3, 0, 0, 0], [-2, 0, 0, 0], [-2, 0, 0, 0], [0, 0, 0, 0]]);
  assert.equal(upperBound(ev, -4), 0);
  assert.equal(upperBound(ev, -2), 12);
  assert.equal(upperBound(ev, 5), 16);
});
