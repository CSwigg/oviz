// Timeline dock: transport controls, a time readout and the time slider.
//
// The slider shows the played part of the timeline, a faint silhouette of
// when the visible objects were born (brighter where time has already
// played), major ticks, a present-day mark and the saved views as dots
// you can click. A bubble reads out the time under the pointer, and
// follows the thumb while you drag.

import { h, icon, iconButton } from "./dom.js";
import { niceFloor } from "../core/math.js";
import { rafThrottle } from "../core/emitter.js";
import { birthProfile } from "../app/timeline.js";

const SPEEDS = [0.25, 0.5, 1, 2, 4, 8];
const HEAT_BINS = 96;
const HEAT_H = 14;

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
    // Click plays forward (Shift+click backward), like Space and ⇧ Space.
    this.playBtn = iconButton("play", "Play", (e) => this.togglePlay(e?.shiftKey ? -1 : 1), { cls: "ov-play", shortcut: "Space" });
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
    // The birth silhouette follows the visible layers and the slider width.
    const heat = rafThrottle(() => this.renderHeat());
    new ResizeObserver(heat).observe(this.scrub.el);
    v.on("style", heat);
    v.on("state-applied", heat);
    v.on("mode", heat);
    v.on("theme", heat);
  }

  /**
   * Space: play forward, or pause (classic Oviz). From the last frame
   * forward play starts over at the first; `direction` −1 (⇧ Space)
   * plays backward, from the first frame over at the last.
   */
  togglePlay(direction = 1) {
    const tl = this.tl;
    // Space pauses whatever is playing; ⇧ Space reverses forward playback.
    if (tl.playing && (direction >= 0 || tl.direction < 0)) tl.pause();
    else tl.play(direction < 0 ? -1 : 1);
  }

  speedMenu(anchor) {
    const tl = this.tl;
    this.ui.menu(anchor, [
      ...SPEEDS.map((s) => ({
        label: `${s}×`,
        checked: tl.speed === s,
        run: () => { tl.speed = s; this.sync(); },
      })),
      "-",
      { label: "Play forward", icon: "play", shortcut: "Space", run: () => tl.play(1) },
      { label: "Play backward", icon: "stepBack", shortcut: "⇧ Space", run: () => tl.play(-1) },
    ], { title: "Playback speed", align: "center", above: true });
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
    const heat = h("div", { class: "ov-scrub-heat", "aria-hidden": "true" });
    const heatLit = h("div", { class: "ov-scrub-heat-lit", "aria-hidden": "true" });
    const track = h("div", { class: "ov-scrub-track" });
    const fill = h("div", { class: "ov-scrub-fill" });
    track.append(fill);
    const ticks = h("div", { class: "ov-scrub-ticks" });
    const thumb = h("div", { class: "ov-scrub-thumb" });
    const bubble = h("div", { class: "ov-scrub-hover" });
    const markers = h("div");
    el.append(heat, heatLit, track, ticks, markers, thumb, bubble);
    // Major ticks at nice intervals (labels show in detailed mode).
    const span = tl.max - tl.min;
    if (span > 0) {
      const major = niceFloor(span / 4) || span;
      const start = Math.ceil(tl.min / major - 1e-9) * major;
      for (let t = start; t <= tl.max + 1e-9; t += major) {
        const left = `${(tl.timeToFrame(t) / (tl.count - 1)) * 100}%`;
        ticks.append(h("div", { class: "ov-scrub-tick", style: { left } }), h("div", { class: "ov-scrub-tick-label", style: { left } }, fmtTick(t)));
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
    let dragging = false;
    el.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      dragging = true;
      tl.pause();
      el.setPointerCapture(e.pointerId);
      el.dataset.dragging = "true";
      tl.setFrame(frac(e.clientX) * (tl.count - 1));
    });
    el.addEventListener("pointermove", (e) => {
      if (dragging) {
        tl.setFrame(frac(e.clientX) * (tl.count - 1));
        return;
      }
      const marker = e.target.closest?.(".ov-scrub-marker");
      if (marker) {
        this.setBubble(Number(marker.dataset.frac), marker.dataset.name || "Saved view", this.ui.formatTime(Number(marker.dataset.time)));
        return;
      }
      const f = frac(e.clientX);
      this.setBubble(f, this.ui.formatTime(tl.frameToTime(f * (tl.count - 1))));
    });
    const end = () => {
      if (!dragging) return;
      dragging = false;
      el.dataset.dragging = "false";
      // Snap gently to the nearest frame when released close to one.
      const f = tl.frame;
      if (Math.abs(f - Math.round(f)) < 0.12) tl.setFrame(Math.round(f));
    };
    el.addEventListener("pointerup", end);
    el.addEventListener("pointercancel", end);
    el.addEventListener("keydown", (e) => {
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      let handled = true;
      if (e.key === "ArrowLeft" || e.key === "ArrowDown") tl.step(e.shiftKey ? -5 : -1);
      else if (e.key === "ArrowRight" || e.key === "ArrowUp") tl.step(e.shiftKey ? 5 : 1);
      else if (e.key === "Home") tl.setFrame(0);
      else if (e.key === "End") tl.setFrame(tl.count - 1);
      else handled = false;
      // Every other key (Space, Esc, U, V…) still reaches the figure.
      if (handled) { e.preventDefault(); e.stopPropagation(); }
    });
    return { el, fill, thumb, markers, bubble, heat, heatLit };
  }

  setBubble(f, text, sub = "") {
    const b = this.scrub.bubble;
    b.style.left = `${f * 100}%`;
    b.replaceChildren(text, sub ? h("small", null, sub) : "");
  }

  /**
   * Saved views as clickable dots on the track. `items` are
   * {time, name, index} (plain times are accepted too).
   */
  setMarkers(items) {
    const s = this.scrub;
    if (!s) return;
    const tl = this.tl;
    s.markers.replaceChildren(...items.map((it) => {
      const m = typeof it === "number" ? { time: it } : it;
      const f = tl.timeToFrame(m.time) / (tl.count - 1);
      const b = h("button", {
        class: "ov-scrub-marker", type: "button",
        "aria-label": m.name ? `Go to view: ${m.name}` : "Saved view",
        style: { left: `${f * 100}%` },
        dataset: { frac: String(f), time: String(m.time), name: m.name || "" },
      });
      b.addEventListener("pointerdown", (e) => e.stopPropagation());
      b.addEventListener("click", (e) => {
        e.stopPropagation();
        if (m.index != null) this.ui.plugins.find((p) => p.name === "states")?.goTo(m.index);
        else tl.setTime(m.time);
      });
      return b;
    }));
  }

  /** Silhouette of when the visible objects were born along the timeline. */
  renderHeat() {
    const s = this.scrub;
    if (!s) return;
    const v = this.viewer;
    const tl = this.tl;
    const ageArrays = [];
    for (const t of v.traces) {
      if (!t.points || !t.showInLegend || v.state.traces[t.key]?.visible === false) continue;
      const ages = v.data.get(t.key)?.ageNow;
      if (ages) ageArrays.push(ages);
    }
    const smooth = birthProfile(tl, ageArrays, HEAT_BINS);
    let max = 0;
    for (const x of smooth) max = Math.max(max, x);
    const any = max > 0;
    s.heat.hidden = s.heatLit.hidden = !any;
    const w = s.el.clientWidth;
    if (!any || !w) return;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const color = getComputedStyle(this.ui.root).getPropertyValue("--ov-text").trim() || "#fff";
    for (const host of [s.heat, s.heatLit]) {
      const c = host.firstChild || host.appendChild(document.createElement("canvas"));
      c.width = Math.round(w * dpr);
      c.height = Math.round(HEAT_H * dpr);
      c.style.width = `${w}px`;
      const ctx = c.getContext("2d");
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, w, HEAT_H);
      ctx.beginPath();
      ctx.moveTo(0, HEAT_H);
      for (let i = 0; i < HEAT_BINS; i++) ctx.lineTo(((i + 0.5) / HEAT_BINS) * w, HEAT_H - (smooth[i] / max) * (HEAT_H - 1));
      ctx.lineTo(w, HEAT_H);
      ctx.closePath();
      ctx.fillStyle = color;
      ctx.fill();
    }
  }

  sync() {
    const tl = this.tl;
    if (!this.scrub) return;
    const f = tl.count > 1 ? tl.frame / (tl.count - 1) : 1;
    const pct = `${f * 100}%`;
    this.scrub.fill.style.width = pct;
    this.scrub.thumb.style.left = pct;
    this.scrub.heatLit.style.width = pct;
    this.scrub.el.setAttribute("aria-valuenow", tl.time.toFixed(2));
    this.scrub.el.setAttribute("aria-valuetext", this.ui.formatTime(tl.time));
    if (this.scrub.el.dataset.dragging === "true") this.setBubble(f, this.ui.formatTime(tl.time));
    // This runs every frame of playback: rebuild text and buttons (and
    // re-parse the play icon) only when what they show changes.
    const t = tl.time;
    const [num, unit] = splitTime(t);
    if (num !== this._shownTime) {
      this._shownTime = num;
      this.timeValue.replaceChildren(num, h("small", null, unit));
    }
    const caption = Math.abs(t) < 1e-6 ? "Present day" : t < 0 ? "Before present" : "After present";
    if (caption !== this.timeCaption.textContent) this.timeCaption.textContent = caption;
    const transport = `${tl.playing}|${tl.speed}|${tl.loop}`;
    if (transport === this._shownTransport) return;
    this._shownTransport = transport;
    this.ui.root.dataset.playing = String(tl.playing);
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
