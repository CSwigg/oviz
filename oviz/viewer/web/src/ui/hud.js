// Instrument skin HUD: a crosshair on the orbit target with a live readout
// of where the camera is aimed (target, range, azimuth, elevation in 3D;
// Galactic l, b and field in Sky). Hidden, and idle, in the other skins.

import { h } from "./dom.js";
import { rafThrottle } from "../core/emitter.js";

const DEG = 180 / Math.PI;

function num(v, digits = 0) {
  const s = v.toFixed(digits);
  return (/^-0(\.0+)?$/.test(s) ? s.slice(1) : s).replace("-", "−");
}

function dist(pc) {
  return pc >= 1000 ? `${num(pc / 1000, 2)} KPC` : `${num(pc, pc >= 100 ? 0 : 1)} PC`;
}

export class InstrumentHud {
  constructor(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.line1 = h("div");
    this.line2 = h("div");
    this.line3 = h("div");
    this.el = h("div", { class: "ov-hud", "aria-hidden": "true" },
      h("div", { class: "ov-hud-cross" }, h("i"), h("i"), h("b")),
      h("div", { class: "ov-hud-read" }, this.line1, this.line2, this.line3));
    ui.ui.prepend(this.el);
    const update = rafThrottle(() => this.update());
    v.on("camera", update);
    v.on("time", update);
    v.on("viewmode", update);
    v.renderer.afterRender.push(update);
  }

  get active() {
    return document.documentElement.dataset.ovizSkin === "instrument";
  }

  update() {
    if (!this.active) return;
    const v = this.viewer;
    const cam = v.renderer.camera;
    const canvas = v.canvas.getBoundingClientRect();
    const root = this.ui.root.getBoundingClientRect();
    this.el.style.setProperty("--hud-x", `${(canvas.left - root.left + canvas.width / 2).toFixed(1)}px`);
    this.el.style.setProperty("--hud-y", `${(canvas.top - root.top + canvas.height / 2).toFixed(1)}px`);
    const p = v.pose;
    const t = v.timeline.time;
    if (v.state.view.mode === "sky") {
      const f = cam.forward;
      const l = ((Math.atan2(f[1], f[0]) * DEG) + 360) % 360;
      const b = Math.asin(Math.max(-1, Math.min(1, f[2]))) * DEG;
      this.line1.textContent = `LOOK  L ${num(l, 2)}°  B ${b >= 0 ? "+" : ""}${num(b, 2)}°`;
      this.line2.textContent = `FOV   ${num(cam.hfov, 1)}°  FROM SUN`;
    } else {
      const [x, y, z] = p.target;
      const az = ((p.yaw * DEG) % 360 + 360) % 360;
      this.line1.textContent = `TGT   ${num(x)}  ${num(y)}  ${num(z)} PC`;
      this.line2.textContent = `RNG   ${dist(p.distance)}  AZ ${num(az, 1)}°  EL ${num(-p.pitch * DEG, 1)}°`;
    }
    this.line3.textContent = `T     ${Math.abs(t) < 1e-6 ? "0.0" : num(t, 1)} MYR${v.timeline.playing ? "  ▶ RUN" : ""}`;
  }
}
