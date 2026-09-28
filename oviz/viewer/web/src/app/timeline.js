// Continuous time over a discrete frame grid.
//
// The figure stores frames at fixed times; the viewer treats time as a real
// number and interpolates between frames on the GPU. `frame` is the
// fractional frame index, `time` the physical value (Myr).

import { clamp } from "../core/math.js";

export class Timeline {
  constructor(times, { initialIndex = 0, intervalMs = 240 } = {}) {
    this.times = Float64Array.from(times.length ? times : [0]);
    this.count = this.times.length;
    this.frame = clamp(initialIndex, 0, this.count - 1);
    this.playing = false;
    this.direction = 1;
    this.speed = 1; // multiplier on the authored pace
    this.loop = true;
    this.baseFps = 1000 / Math.max(intervalMs, 16); // frames per second at 1×
    this._last = 0;
    this.listeners = new Set();
  }

  get time() {
    return this.frameToTime(this.frame);
  }

  get min() {
    return this.times[0];
  }

  get max() {
    return this.times[this.count - 1];
  }

  frameToTime(f) {
    const n = this.count;
    if (n === 1) return this.times[0];
    const fc = clamp(f, 0, n - 1);
    const i = Math.floor(fc), j = Math.min(i + 1, n - 1);
    return this.times[i] + (this.times[j] - this.times[i]) * (fc - i);
  }

  timeToFrame(t) {
    const T = this.times, n = this.count;
    if (n === 1 || t <= T[0]) return 0;
    if (t >= T[n - 1]) return n - 1;
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (T[mid] <= t) lo = mid; else hi = mid;
    }
    const span = T[hi] - T[lo];
    return lo + (span > 0 ? (t - T[lo]) / span : 0);
  }

  /** Nearest frame index whose time is (approximately) zero, if any. */
  zeroFrame() {
    let best = -1, bestAbs = Infinity;
    for (let i = 0; i < this.count; i++) {
      const a = Math.abs(this.times[i]);
      if (a < bestAbs) { bestAbs = a; best = i; }
    }
    return best;
  }

  setFrame(f, { snap = false, silent = false } = {}) {
    let next = clamp(f, 0, this.count - 1);
    if (snap) next = Math.round(next);
    if (next === this.frame) return;
    this.frame = next;
    if (!silent) this._emit();
  }

  setTime(t, opts) {
    this.setFrame(this.timeToFrame(t), opts);
  }

  step(n = 1) {
    const target = Math.round(this.frame) + n;
    if (this.loop && this.count > 1) {
      this.setFrame(((target % this.count) + this.count) % this.count);
    } else {
      this.setFrame(target);
    }
  }

  play(direction = this.direction) {
    this.direction = direction >= 0 ? 1 : -1;
    if (this.count <= 1) return;
    // Restart from the far end when playing from a finished position.
    if (this.direction > 0 && this.frame >= this.count - 1) this.frame = 0;
    if (this.direction < 0 && this.frame <= 0) this.frame = this.count - 1;
    this.playing = true;
    this._last = 0;
    this._emit();
  }

  pause() {
    if (!this.playing) return;
    this.playing = false;
    this._emit();
  }

  toggle() {
    if (this.playing) this.pause(); else this.play();
  }

  /** Advance playback; returns true when the frame changed. */
  tick(now) {
    if (!this.playing) return false;
    if (!this._last) { this._last = now; return false; }
    const dt = Math.min(now - this._last, 100) / 1000;
    this._last = now;
    let f = this.frame + this.direction * this.baseFps * this.speed * dt;
    const last = this.count - 1;
    if (f > last || f < 0) {
      if (this.loop) {
        f = this.direction > 0 ? f - last : last + f;
        f = clamp(f, 0, last);
      } else {
        f = clamp(f, 0, last);
        this.playing = false;
      }
    }
    this.frame = f;
    this._emit();
    return true;
  }

  on(fn) {
    this.listeners.add(fn);
    return () => this.listeners.delete(fn);
  }

  _emit() {
    for (const fn of this.listeners) fn(this);
  }
}

export function formatTime(t, unit = "Myr") {
  const a = Math.abs(t);
  const digits = a >= 100 ? 0 : a >= 10 ? 1 : 2;
  let s = t.toFixed(digits);
  if (/^-0(\.0+)?$/.test(s)) s = s.slice(1);
  s = s.replace("-", "−");
  return `${s} ${unit}`;
}

/**
 * How many objects were born in each of `bins` equal slices of the
 * timeline (birth time = −age). Births outside the timeline are ignored.
 * Returns a Float64Array, lightly smoothed so the profile reads as a shape.
 */
export function birthProfile(timeline, ageArrays, bins = 96) {
  const counts = new Float64Array(bins);
  const tl = timeline;
  if (tl.count < 2) return counts;
  for (const ages of ageArrays) {
    for (let i = 0; i < ages.length; i++) {
      const birth = -ages[i];
      if (!(birth >= tl.min && birth <= tl.max)) continue;
      counts[Math.min(bins - 1, Math.floor((tl.timeToFrame(birth) / (tl.count - 1)) * bins))] += 1;
    }
  }
  const kernel = [0.06, 0.24, 0.4, 0.24, 0.06];
  const out = new Float64Array(bins);
  for (let i = 0; i < bins; i++) {
    let acc = 0;
    for (let j = -2; j <= 2; j++) acc += counts[Math.min(bins - 1, Math.max(0, i + j))] * kernel[j + 2];
    out[i] = acc;
  }
  return out;
}
