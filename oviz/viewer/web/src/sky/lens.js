// Sky lens (the classic Oviz sky aperture): a circle pinned to a spot on the
// sky that shows it in another survey, a spyglass across wavelengths. A
// second Aladin view, kept on the main view's pose every frame, is clipped
// to the circle. Drag its ring to move it, scroll over it to resize, step
// through bands with its arrows, double-click to make that survey the
// background. Lives in state.sky.lens, so States keep it.

import { h, icon, iconButton } from "../ui/dom.js";
import { loadAladin, resolveSurvey } from "./aladin.js";
import { formatWavelength, surveyInfo } from "./catalog.js";
import { wavelengthStops } from "./stack.js";
import { lbToDir, dirToLb, clamp } from "../core/math.js";

const SVG = "http://www.w3.org/2000/svg";

export class SkyLens {
  constructor(plugin) {
    this.plugin = plugin;
    this.ui = plugin.ui;
    this.viewer = plugin.viewer;
    const ap = plugin.spec?.aperture || {};
    this.minDeg = Math.max(0.2, Number(ap.minSizeDeg) || 3);
    this.maxDeg = Math.max(this.minDeg, Number(ap.maxSizeDeg) || 48);
    this.defaultDeg = clamp(Number(ap.defaultSizeDeg) || 14, this.minDeg, this.maxDeg);
    this.aladin = null;
    this.survey = null;
  }

  get state() {
    return this.viewer.state.sky.lens || null;
  }

  /** Open a lens at the view centre (or toggle it off). */
  toggle() {
    const st = this.viewer.state.sky;
    if (st.lens) { this.close(); return; }
    const f = this.viewer.renderer.camera.forward;
    const [l, b] = dirToLb(f);
    const stops = wavelengthStops(st.layers || []);
    // Start on a band that differs from the background: dust if there is one.
    const base = (st.layers || []).filter((x) => x.visible !== false).pop();
    const pick = stops.find((s) => /planck|dust/i.test(`${s.id} ${s.label}`) && s.id !== base?.survey) || stops.find((s) => s.id !== base?.survey) || stops[0];
    st.lens = { l, b, sizeDeg: this.defaultDeg, survey: pick?.id || "P/PLANCK/R2/HFI/color", opacity: 1 };
    this.sync();
    this.ui.toast("Sky lens: drag its ring, scroll over it to resize, double-click to make it the background", { icon: icon("globe"), ms: 3800 });
  }

  close() {
    this.viewer.state.sky.lens = null;
    this.sync();
  }

  /** Match the DOM (and the second Aladin view) to state.sky.lens. */
  sync() {
    const lens = this.state;
    const inSky = this.viewer.state.view.mode === "sky";
    if (!lens || !inSky) {
      if (this.el) this.el.hidden = true;
      if (this.div) this.div.style.visibility = "hidden";
      this.plugin.changed?.();
      return;
    }
    this.ensureDom();
    this.el.hidden = false;
    this.ensureAladin().then(() => this.setSurvey(lens.survey));
    this.renderChip();
    this.plugin.changed?.();
  }

  ensureDom() {
    if (this.el) return;
    const ring = document.createElementNS(SVG, "svg");
    ring.setAttribute("class", "ov-lens-ring");
    this.circle = document.createElementNS(SVG, "circle");
    this.circleHit = document.createElementNS(SVG, "circle");
    this.circleHit.setAttribute("class", "ov-lens-hit");
    ring.append(this.circleHit, this.circle);
    this.chip = h("div", { class: "ov-lens-chip ov-glass" });
    this.el = h("div", { class: "ov-lens", hidden: true }, ring, this.chip);
    this.ring = ring;
    this.ui.ui.append(this.el);
    this.bind();
  }

  async ensureAladin() {
    if (this.aladin || this._starting) return this._starting;
    this._starting = (async () => {
      const A = await loadAladin();
      const div = document.createElement("div");
      div.className = "ov-aladin ov-aladin-lens";
      div.style.cssText = "position:absolute;inset:0;visibility:hidden;";
      this.plugin.host.append(div);
      this.div = div;
      const lens = this.state;
      this.aladin = A.aladin(div, {
        survey: resolveSurvey(lens?.survey || "P/DSS2/color"), fov: 60, target: "0 0", cooFrame: "galactic", projection: "TAN",
        showReticle: false, showZoomControl: false, showFullscreenControl: false, showLayersControl: false, showGotoControl: false,
        showShareControl: false, showSimbadPointerControl: false, showCooGridControl: false, showSettingsControl: false,
        showProjectionControl: false, showStatusBar: false, showCooLocation: false, showFrame: false, showFov: false,
        showContextMenu: false, showCooGrid: false, lockNorthUp: true, inertia: false, realFullscreen: false,
      });
      try { this.aladin.setProjection?.("TAN"); } catch (_) { /* older API */ }
      this.survey = resolveSurvey(lens?.survey || "P/DSS2/color");
    })();
    return this._starting;
  }

  setSurvey(id) {
    if (!this.aladin || !id) return;
    const url = resolveSurvey(id);
    if (url === this.survey) return;
    try { (this.aladin.setImageSurvey || this.aladin.setBaseImageLayer).call(this.aladin, url); this.survey = url; } catch (err) { console.warn("Oviz lens: survey failed", id, err); }
  }

  renderChip() {
    const lens = this.state;
    if (!lens || !this.chip) return;
    const info = surveyInfo(lens.survey);
    const stops = wavelengthStops(this.viewer.state.sky.layers || []);
    const step = (d) => {
      const i = stops.findIndex((s) => surveyInfo(s.id).id === info.id);
      const next = stops[((i < 0 ? 0 : i + d) % stops.length + stops.length) % stops.length];
      if (next) { lens.survey = next.id; this.setSurvey(next.id); this.renderChip(); this.plugin.changed?.(); }
    };
    this.chip.replaceChildren(
      iconButton("stepBack", "Shorter wavelength", () => step(-1), { cls: "ov-lens-btn" }),
      h("span", { class: "ov-lens-label", title: info.note || lens.survey },
        h("b", null, info.label), Number.isFinite(info.lambda) ? ` · ${formatWavelength(info.lambda)}` : info.band ? ` · ${info.band}` : ""),
      iconButton("stepFwd", "Longer wavelength", () => step(1), { cls: "ov-lens-btn" }),
      iconButton("close", "Close the lens", () => this.close(), { cls: "ov-lens-btn" }),
    );
  }

  bind() {
    const v = this.viewer;
    let drag = null;
    this.circleHit.addEventListener("pointerdown", (e) => {
      e.stopPropagation();
      drag = { id: e.pointerId };
      try { this.circleHit.setPointerCapture(e.pointerId); } catch (_) { /* synthetic */ }
    });
    this.circleHit.addEventListener("pointermove", (e) => {
      if (!drag) return;
      const r = v.canvas.getBoundingClientRect();
      const ray = v.renderer.camera.ray(e.clientX - r.left, e.clientY - r.top);
      const [l, b] = dirToLb(ray.dir);
      Object.assign(this.state, { l, b });
      this.place();
    });
    const end = () => { if (drag) { drag = null; this.plugin.changed?.(); } };
    this.circleHit.addEventListener("pointerup", end);
    this.circleHit.addEventListener("pointercancel", end);
    const wheel = (e) => {
      e.preventDefault();
      e.stopPropagation();
      const lens = this.state;
      if (!lens) return;
      lens.sizeDeg = clamp(lens.sizeDeg * Math.exp(clamp(e.deltaY, -120, 120) * 0.004), this.minDeg, this.maxDeg);
      this.place();
    };
    this.circleHit.addEventListener("wheel", wheel, { passive: false });
    this.chip.addEventListener("wheel", wheel, { passive: false });
    this.circleHit.addEventListener("dblclick", (e) => {
      e.stopPropagation();
      const lens = this.state;
      if (!lens) return;
      // Promote: the lens survey becomes the background (classic).
      const id = lens.survey;
      this.close();
      this.plugin.showOnly(id);
    });
  }

  /** Per frame (Sky view): keep the second view on the main pose and the circle on its spot. */
  tick(lDeg, bDeg, hfovDeg) {
    const lens = this.state;
    if (!lens || !this.aladin || this.viewer.state.view.mode !== "sky") return;
    const sig = `${lDeg.toFixed(5)}|${bDeg.toFixed(5)}|${hfovDeg.toFixed(4)}`;
    if (sig !== this._viewSig) {
      try {
        this.aladin.gotoPosition(((lDeg % 360) + 360) % 360, bDeg);
        this.aladin.setFoV(Math.min(Math.max(hfovDeg, 1e-4), 179));
        this._viewSig = sig;
      } catch (_) { /* warming up */ }
    }
    this.place();
  }

  /** Put the circle (clip and ring) where the lens's sky spot projects now. */
  place() {
    const lens = this.state;
    if (!lens || !this.el || !this.div) return;
    const cam = this.viewer.renderer.camera;
    const c = lbToDir(lens.l, lens.b);
    const far = 1e5;
    const s = cam.project([c[0] * far, c[1] * far, c[2] * far], [0, 0, 0]);
    // The radius: a point half the lens size away along the sky.
    const edge = lbToDir(lens.l, clamp(lens.b + lens.sizeDeg / 2, -89.9, 89.9));
    const s2 = cam.project([edge[0] * far, edge[1] * far, edge[2] * far], [0, 0, 0]);
    if (!s || !s2) { this.el.style.visibility = "hidden"; this.div.style.visibility = "hidden"; return; }
    const r = Math.max(8, Math.hypot(s2[0] - s[0], s2[1] - s[1]));
    this.el.style.visibility = "";
    this.div.style.visibility = "visible";
    this.div.style.clipPath = `circle(${r.toFixed(1)}px at ${s[0].toFixed(1)}px ${s[1].toFixed(1)}px)`;
    this.div.style.opacity = String(clamp(lens.opacity ?? 1, 0, 1));
    const d = 2 * r + 16;
    Object.assign(this.el.style, { left: `${(s[0] - r - 8).toFixed(1)}px`, top: `${(s[1] - r - 8).toFixed(1)}px`, width: `${d.toFixed(1)}px`, height: `${d.toFixed(1)}px` });
    this.ring.setAttribute("viewBox", `0 0 ${d} ${d}`);
    for (const el of [this.circle, this.circleHit]) {
      el.setAttribute("cx", String(d / 2));
      el.setAttribute("cy", String(d / 2));
      el.setAttribute("r", String(r));
    }
  }
}
