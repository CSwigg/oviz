// Timeline dock: transport controls, time readout and a continuous scrubber.

import { h, icon, iconButton } from "./dom.js";
import { niceFloor } from "../core/math.js";

const SPEEDS = [0.25, 0.5, 1, 2, 4, 8];

export class TimelineDock {
  constructor(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    const tl = (this.tl = v.timeline);
    this.el = h("div", { class: "ov-dock ov-glass" });
    if (tl.count <= 1 || v.manifest.time?.enabled === false) {
      this.el.hidden = true;
      return;
    }
    this.back = iconButton("stepBack", "Previous frame", () => { tl.pause(); tl.step(-1); }, { shortcut: "←" });
    this.playBtn = iconButton("play", "Play", () => this.togglePlay(), { cls: "ov-play", shortcut: "Space" });
    this.fwd = iconButton("stepFwd", "Next frame", () => { tl.pause(); tl.step(1); }, { shortcut: "→" });
    this.timeValue = h("div", { class: "ov-time-value" });
    this.timeCaption = h("div", { class: "ov-time-caption" }, "Time");
    const read = h("div", { class: "ov-time-read", "aria-live": "off" }, this.timeValue, this.timeCaption);
    this.scrub = this.buildScrubber();
    this.speed = h("button", { class: "ov-speed", type: "button", "data-tip": "Playback speed  < >", "aria-label": "Playback speed" }, "1×");
    this.speed.addEventListener("click", (e) => this.speedMenu(e.currentTarget));
    this.loop = iconButton("loop", "Loop playback", () => { tl.loop = !tl.loop; this.sync(); }, { pressed: tl.loop });
    this.el.append(this.back, this.playBtn, this.fwd, read, this.scrub.el, this.speed, this.loop);
    tl.on(() => this.sync());
    this.sync();
  }

  togglePlay() {
    const tl = this.tl;
    if (tl.playing) tl.pause();
    else {
      // Play toward the present when starting from the past, and vice versa.
      const dir = tl.frame >= tl.count - 1 ? -1 : tl.direction || 1;
      tl.play(dir);
    }
  }

  speedMenu(anchor) {
    this.ui.menu(anchor, SPEEDS.map((s) => ({
      label: `${s}×`,
      checked: this.tl.speed === s,
      run: () => { this.tl.speed = s; this.sync(); },
    })), { title: "Playback speed", align: "center", above: true });
  }

  cycleSpeed(dir) {
    const i = SPEEDS.indexOf(this.tl.speed);
    const j = Math.max(0, Math.min(SPEEDS.length - 1, (i < 0 ? 2 : i) + dir));
    this.tl.speed = SPEEDS[j];
    this.sync();
    this.ui.toast(`Playback ${SPEEDS[j]}×`, { ms: 900 });
  }

  buildScrubber() {
    const tl = this.tl;
    const el = h("div", { class: "ov-scrub", role: "slider", tabindex: "0", "aria-label": "Time", "aria-valuemin": String(tl.min), "aria-valuemax": String(tl.max) });
    const track = h("div", { class: "ov-scrub-track" });
    const fill = h("div", { class: "ov-scrub-fill" });
    track.append(fill);
    const ticks = h("div", { class: "ov-scrub-ticks" });
    const thumb = h("div", { class: "ov-scrub-thumb" });
    const hover = h("div", { class: "ov-scrub-hover" });
    const markers = h("div");
    el.append(track, ticks, markers, thumb, hover);
    // Ticks at nice time intervals, labels on majors.
    const span = tl.max - tl.min;
    if (span > 0) {
      const major = niceFloor(span / 4) || span;
      const minor = major / 5;
      const start = Math.ceil(tl.min / minor) * minor;
      for (let t = start; t <= tl.max + 1e-9; t += minor) {
        const f = (tl.timeToFrame(t)) / (tl.count - 1);
        const isMajor = Math.abs(t / major - Math.round(t / major)) < 1e-6;
        ticks.append(h("div", { class: `ov-scrub-tick${isMajor ? " ov-scrub-tick--major" : ""}`, style: { left: `${f * 100}%` } }));
        if (isMajor) ticks.append(h("div", { class: "ov-scrub-tick-label", style: { left: `${f * 100}%` } }, fmtTick(t)));
      }
      const z = tl.zeroFrame();
      if (z >= 0 && Math.abs(tl.times[z]) < 1e-9 && z > 0 && z < tl.count - 1) {
        el.append(h("div", { class: "ov-scrub-zero", style: { left: `${(z / (tl.count - 1)) * 100}%` }, title: "Present day" }));
      }
    }
    const frac = (clientX) => {
      const r = el.getBoundingClientRect();
      return Math.min(Math.max((clientX - r.left) / r.width, 0), 1);
    };
    let dragging = false, wasPlaying = false;
    el.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      dragging = true;
      wasPlaying = tl.playing;
      tl.pause();
      el.setPointerCapture(e.pointerId);
      el.dataset.dragging = "true";
      tl.setFrame(frac(e.clientX) * (tl.count - 1));
    });
    el.addEventListener("pointermove", (e) => {
      const f = frac(e.clientX);
      hover.style.left = `${f * 100}%`;
      hover.textContent = this.ui.formatTime(tl.frameToTime(f * (tl.count - 1)));
      if (dragging) tl.setFrame(f * (tl.count - 1));
    });
    const end = (e) => {
      if (!dragging) return;
      dragging = false;
      el.dataset.dragging = "false";
      // Snap gently to the nearest frame when released close to one.
      const f = tl.frame;
      if (Math.abs(f - Math.round(f)) < 0.12) tl.setFrame(Math.round(f));
      if (wasPlaying && e.type === "pointerup") { /* stay paused after a scrub */ }
    };
    el.addEventListener("pointerup", end);
    el.addEventListener("pointercancel", end);
    el.addEventListener("keydown", (e) => {
      if (e.key === "ArrowLeft" || e.key === "ArrowDown") { tl.step(e.shiftKey ? -5 : -1); e.preventDefault(); }
      else if (e.key === "ArrowRight" || e.key === "ArrowUp") { tl.step(e.shiftKey ? 5 : 1); e.preventDefault(); }
      else if (e.key === "Home") { tl.setFrame(0); e.preventDefault(); }
      else if (e.key === "End") { tl.setFrame(tl.count - 1); e.preventDefault(); }
      e.stopPropagation();
    });
    return { el, fill, thumb, markers };
  }

  setMarkers(times) {
    const s = this.scrub;
    if (!s) return;
    s.markers.replaceChildren(...times.map((t) => h("div", { class: "ov-scrub-marker", style: { left: `${(this.tl.timeToFrame(t) / (this.tl.count - 1)) * 100}%` } })));
  }

  sync() {
    const tl = this.tl;
    if (!this.scrub) return;
    const f = tl.count > 1 ? tl.frame / (tl.count - 1) : 1;
    this.scrub.fill.style.width = `${f * 100}%`;
    this.scrub.thumb.style.left = `${f * 100}%`;
    this.scrub.el.setAttribute("aria-valuenow", tl.time.toFixed(2));
    this.scrub.el.setAttribute("aria-valuetext", this.ui.formatTime(tl.time));
    const t = tl.time;
    const [num, unit] = splitTime(t);
    this.timeValue.replaceChildren(num, h("small", null, unit));
    this.timeCaption.textContent = Math.abs(t) < 1e-6 ? "Present day" : t < 0 ? "Before present" : "After present";
    this.playBtn.replaceChildren(icon(tl.playing ? "pause" : "play"));
    this.playBtn.setAttribute("aria-label", tl.playing ? "Pause" : "Play");
    this.playBtn.dataset.tip = `${tl.playing ? "Pause" : "Play"}  Space`;
    this.speed.textContent = `${tl.speed}×`;
    this.loop.setAttribute("aria-pressed", String(tl.loop));
  }
}

function splitTime(t) {
  const a = Math.abs(t);
  const digits = a >= 100 ? 0 : 1;
  let s = t.toFixed(digits);
  if (/^-0(\.0)?$/.test(s)) s = s.slice(1);
  return [s.replace("-", "−"), "Myr"];
}

function fmtTick(t) {
  const s = Math.abs(t) < 1e-9 ? "0" : String(+t.toFixed(2));
  return s.replace("-", "−");
}
