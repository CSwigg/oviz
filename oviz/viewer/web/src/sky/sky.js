// Sky view plugin: Aladin background, survey layer stack UI and member stars.

import { h, icon, iconButton, clear, MOD } from "../ui/dom.js";
import { slider, miniSeg, select, numberInput } from "../ui/controls.js";
import { AladinSky } from "./aladin.js";
import { MembersLayer } from "./members.js";
import { SkyPicker } from "./picker.js";
import { SkyIdentify, resolveName, lookAtPose } from "./identify.js";
import { SkyLens } from "./lens.js";
import { surveyInfo, surveyLabel, thumbnailUrl, formatWavelength, SKY_SURVEYS } from "./catalog.js";
import { showOnlySurvey, blendSurvey, hideSurvey, crossfadeToWavelength, layersInView } from "./stack.js";
import { parseColor } from "../core/color.js";
import { smoothstep, clamp, icrsToGal } from "../core/math.js";
import { cloneJson } from "../app/state.js";

// Aladin colormaps and stretches offered for a survey (classic Oviz list).
const SKY_COLORMAPS = [
  ["native", "Native colours"], ["grayscale", "Gray"], ["viridis", "Viridis"], ["plasma", "Plasma"], ["magma", "Magma"],
  ["inferno", "Inferno"], ["cividis", "Cividis"], ["rainbow", "Rainbow"], ["cubehelix", "Cubehelix"],
];

const DEG = 180 / Math.PI;

export class SkyPlugin {
  constructor() {
    this.name = "sky";
    this.expanded = null;
    this.group = "all";
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.spec = v.manifest.sky;
    if (!this.spec?.enabled) return;
    this.host = h("div", { class: "ov-sky-host", "aria-hidden": "true" });
    v.container.prepend(this.host);
    this.sky = new AladinSky(this.host);
    this.sky.onLayerStatus = (key, status) => this.onLayerStatus(key, status);
    this.sky.onNeedApply = () => this.applyLayers();
    this.picker = new SkyPicker(this);
    // Right-click (long-press) the sky: what SIMBAD lists there.
    this.identify = new SkyIdentify(this);
    // The sky lens (classic aperture): another survey inside a circle.
    this.lens = this.spec.aperture?.enabled === false ? null : new SkyLens(this);
    this.ensureBaseLayer();
    this.memberTraceKeys = new Set(v.traces.filter((t) => t.points?.members).map((t) => t.key));
    v.memberTraces = this.memberTraceKeys;
    v.memberReveal = 0;
    v.skyMembersActive = () => !!this.members && v.state.sky.members === "stars";
    v.skyPickDraw = (frame) => {
      if (this.members && v.state.view.mode === "sky") this.members.drawPick(frame);
    };
    v.skyPickResolve = (best) => this.resolveMemberPick(best);
    v.renderer.beforeRender.push((now) => this.tick(now));
    this.bgFade = 0;
    v.on("viewmode", ({ mode }) => this.onViewMode(mode));
    v.on("viewmode-settled", ({ mode }) => { if (mode === "sky") { this.revealMembers(true); this.preloadStateSurveys(); } });
    v.on("gpu-restored", () => {
      // The member layer's GL objects died with the context: rebuild lazily.
      this.members = null;
      this._arrowSig = null;
      this._membersLoading = null;
      this._clusterSig = "";
      if (v.state.view.mode === "sky") this.ensureMembers();
    });
    // Quick controls next to the 3D/Sky switch.
    this.memberSeg = miniSeg({
      label: "",
      value: v.state.sky.members,
      options: [{ value: "clusters", label: "Clusters" }, { value: "stars", label: "Stars" }],
      onChange: (m) => this.setMemberMode(m),
    });
    this.memberSeg.classList.add("ov-glass");
    this.memberSeg.style.cssText = "height:40px;border-radius:999px;padding:3px;align-items:center;";
    // The background button reads like a map-type switch: a thumbnail and
    // the name of the sky behind the data. It opens the picker (B still
    // shows or hides the background, Shift+B opens the picker anywhere).
    this.bgThumb = h("img", { class: "ov-skybg-thumb", alt: "" });
    this.bgThumb.addEventListener("error", () => { this.bgThumb.style.visibility = "hidden"; });
    this.bgLabel = h("span", { class: "ov-skybg-label" });
    this.bgBtn = h("button", { class: "ov-skybg ov-glass", type: "button", "aria-haspopup": "dialog", onclick: (e) => this.picker.toggle(e.currentTarget) },
      this.bgThumb, this.bgLabel, icon("chevron", "ov-skybg-chev"));
    ui.skyExtras.append(this.memberSeg, this.bgBtn);
    // Credit for the imagery, beside the scale while a survey is showing (classic).
    this.credit = h("a", { class: "ov-sky-credit", href: "https://aladin.cds.unistra.fr/AladinLite/", target: "_blank", rel: "noopener noreferrer", title: "Sky imagery from the CDS, Strasbourg, through Aladin Lite", hidden: true }, "Aladin Lite · CDS");
    ui.scaleRow?.append(this.credit);
    // Keep the picker's "today's sky" note current while time moves.
    v.on("time", () => { if (this.picker.open) this.picker.sync(); });
    this.syncExtras();
    ui.sky = this;
    // Warm the background up once the figure is idle so Sky opens instantly.
    const idle = window.requestIdleCallback || ((fn) => setTimeout(fn, 1500));
    idle(() => this.ensureSky().catch(() => {}), { timeout: 4000 });
  }

  afterBoot() {
    if (!this.spec?.enabled) return;
    this.syncExtras();
  }

  commands() {
    if (!this.spec?.enabled) return [];
    const v = this.viewer;
    const toSky = () => { if (v.state.view.mode !== "sky") this.ui.setViewMode("sky"); };
    return [
      { title: "Choose sky background…", icon: "stars", keywords: "survey hips aladin imagery wavelength map picker add", run: () => this.openPicker() },
      ...(this.lens ? [{ title: v.state.sky.lens ? "Close the sky lens" : "Sky lens: another wavelength in a circle", icon: "globe", keywords: "aperture spyglass compare band wavelength patch", run: () => { if (v.state.view.mode !== "sky") this.ui.setViewMode("sky"); this.lens.toggle(); } }] : []),
      { title: v.state.sky.backgroundVisible ? "Hide sky background" : "Show sky background", icon: "stars", shortcut: "B", run: () => this.setBackground(!v.state.sky.backgroundVisible) },
      ...(this.spec.members ? [{ title: v.state.sky.members === "stars" ? "Show cluster markers in Sky" : "Show member stars in Sky", icon: "sparkle", run: () => this.setMemberMode(v.state.sky.members === "stars" ? "clusters" : "stars") }] : []),
      ...(this.spec.members ? [this.arrowsCommand()] : []),
      ...this.pickerSurveys().map((t) => ({
        title: `Sky background: ${t.label}`,
        sub: [t.band, Number.isFinite(t.lambda) ? formatWavelength(t.lambda) : ""].filter(Boolean).join(" · "),
        icon: "globe",
        keywords: `survey hips aladin ${t.id} ${t.band || ""} ${t.note || ""}`,
        run: () => { toSky(); this.showOnly(t.id, { label: t.label }); },
      })),
    ];
  }

  /** Motion arrows for every visible trace with member stars (classic per-trace switch). */
  arrowsCommand() {
    const v = this.viewer;
    const keys = [...this.memberTraceKeys].filter((k) => v.state.traces[k]?.visible);
    const on = keys.some((k) => v.state.traces[k]?.memberArrows);
    return {
      title: on ? "Hide member motion arrows" : "Show member motion arrows",
      sub: "Each star's motion relative to its cluster, in Sky view",
      icon: "trail",
      keywords: "proper motion kinematics vectors arrows members stars sky",
      run: () => {
        for (const k of keys) v.setTraceStyle(k, { memberArrows: !on });
        if (!on && v.state.view.mode !== "sky") this.ui.toast("Motion arrows show on member stars in Sky view (V)", { icon: icon("trail") });
      },
    };
  }

  helpGroups() {
    if (!this.spec?.enabled) return [];
    return [["Sky view", [
      ["Choose or blend backgrounds", "✦"],
      ["Show / hide the background", "B"],
      ["What is here? (SIMBAD)", "Right-click"],
      ["Find a name on the sky", `${MOD} K`],
      ["Clusters ⇄ member stars", "Clusters Stars"],
    ]]];
  }

  /** Search-box entries for a query: find a named object on the sky (Sesame). */
  searchItems(q) {
    const text = String(q || "").trim();
    if (!this.spec?.enabled || text.length < 2) return [];
    return [{
      title: `Find “${text}” on the sky`,
      sub: "Look up the name in SIMBAD and turn the Sky view to it",
      icon: "globe",
      hint: "Sky",
      run: () => this.findOnSky(text),
    }];
  }

  /** Resolve a name (Sesame: SIMBAD, NED, VizieR) and turn the Sky view to it. */
  async findOnSky(name) {
    const v = this.viewer;
    this.ui.toast(`Looking up ${name}…`, { icon: icon("search"), ms: 1400 });
    let hit = null;
    try { hit = await resolveName(name); } catch (_) { hit = null; }
    if (!hit || !Number.isFinite(hit.ra) || !Number.isFinite(hit.dec)) {
      this.ui.toast(`SIMBAD does not know “${name}”`);
      return;
    }
    const [l, b] = icrsToGal(hit.ra, hit.dec);
    if (v.state.view.mode !== "sky") await this.ui.setViewMode("sky");
    await this.lookAt(l, b, { fov: Math.min(v.pose.fov, 20) });
    this.ui.toast(`${hit.name || name}${hit.type ? ` · ${hit.type}` : ""} · l ${l.toFixed(2)}°, b ${b.toFixed(2)}°`, { icon: icon("globe"), ms: 3600 });
  }

  /** Turn the Sky camera toward galactic (l, b). */
  lookAt(lDeg, bDeg, opts = {}) {
    const v = this.viewer;
    if (v.state.view.mode !== "sky") return Promise.resolve({ done: false });
    return v.animateTo(lookAtPose(v.pose, lDeg, bDeg, opts), { duration: 900 });
  }

  /** The surveys the picker lists, for the command palette. */
  pickerSurveys() {
    const seen = new Set();
    const out = [];
    for (const id of [...this.layers().map((l) => l.survey || l.key), ...SKY_SURVEYS.map((x) => x.id)]) {
      const info = surveyInfo(id);
      if (seen.has(info.id)) continue;
      seen.add(info.id);
      const own = this.layers().find((l) => (l.survey || l.key) === id);
      out.push({ ...info, label: !info.known && own?.label ? own.label : info.label });
    }
    return out;
  }

  /** Open the background picker (switching to Sky view first if needed). */
  openPicker() {
    const v = this.viewer;
    // The picker opens straight away; the camera finishes its flight to the Sun.
    if (v.state.view.mode !== "sky") this.ui.setViewMode("sky");
    if (this.ui.layersOpen && this.ui.layersFloating) this.ui.setLayersOpen(false);
    this.syncExtras();
    requestAnimationFrame(() => this.picker.show(this.bgBtn));
  }

  /** The more menu (phones and Focus mode): the background picker, from either view. */
  moreItems() {
    if (!this.spec?.enabled) return [];
    return [{ label: "Sky background…", icon: "stars", shortcut: "⇧B", run: () => this.openPicker() }];
  }

  onKey(e) {
    // Shift+B opens the background picker, from 3D view too.
    if (e.key === "B" && e.shiftKey && !e.metaKey && !e.ctrlKey && !e.altKey) {
      if (this.picker.open) this.ui.closeMenu();
      else this.openPicker();
      return true;
    }
    // Like the classic viewer, B belongs to Sky view.
    if ((e.key === "b" || e.key === "B") && !e.metaKey && !e.ctrlKey && this.viewer.state.view.mode === "sky") {
      this.setBackground(!this.viewer.state.sky.backgroundVisible);
      this.ui.toast(this.viewer.state.sky.backgroundVisible ? "Sky background on" : "Sky background off", { icon: icon("stars"), ms: 900 });
      return true;
    }
    return false;
  }

  // ------------------------------------------------------------ lifecycle

  async ensureSky() {
    if (this.sky.ready) return;
    await this.sky.init({ survey: this.baseSurvey(), fov: 60 });
    const r = this.viewer.renderer;
    // In Sky view the viewer renders inside Aladin's frame (see AladinSky).
    this.sky.onBeforeDraw = (t) => r.externalTick(t);
    r.externalDriver = () => this.sky.synchronized && this.viewer.state.view.mode === "sky";
    this.applyLayers();
  }

  /**
   * Once in Sky view, attach (hidden) the backgrounds the saved views use,
   * a few at a time, so going to a view never waits for its tiles.
   */
  preloadStateSurveys() {
    if (this._preloaded || !this.sky.ready) return;
    this._preloaded = true;
    const story = this.ui.plugins.find((p) => p.name === "states");
    const inView = new Set(layersInView(this.layers(), true).map((l) => l.key));
    const want = [];
    for (const it of story?.project?.items || []) {
      for (const l of it.state?.sky?.layers || []) {
        if (l.visible !== false && (l.opacity ?? 1) > 0.001 && !inView.has(l.key) && !want.some((x) => x.key === l.key)) want.push(l);
      }
    }
    want.slice(0, 4).forEach((l, i) => setTimeout(() => this.sky.preload(l), 1500 + i * 1200));
  }

  async ensureMembers() {
    if (this.members || this._membersLoading || !this.spec.members) return this._membersLoading;
    const v = this.viewer;
    const m = this.spec.members;
    this._membersLoading = (async () => {
      const [offsets, position, velocity, valid] = await Promise.all([
        v.store.get(m.offsets.blob), v.store.get(m.position.blob), v.store.get(m.velocity.blob), v.store.get(m.valid.blob),
      ]);
      // Links: trace object → member cluster.
      this.links = [];
      for (const t of v.traces) {
        if (!t.points?.members) continue;
        this.links.push({ key: t.key, link: await v.store.get(t.points.members.blob) });
      }
      this.members = new MembersLayer(v.gl, { count: m.count, clusters: m.clusters, position, velocity, valid, offsets });
      v.renderer.add(this.members);
      // Flag parent markers that member stars replace in Sky view.
      for (const { key, link } of this.links) {
        const batch = v.points.byKey.get(key);
        if (!batch) continue;
        const states = new Uint8Array(batch.count);
        const hasStars = this.members.clusterHasStars;
        for (let i = 0; i < link.length; i++) states[i] = link[i] >= 0 && hasStars[link[i]] ? 4 : 0;
        batch.setStateBits(4, states);
      }
      this.members.visible = true;
      v.renderer.invalidate();
    })();
    return this._membersLoading;
  }

  async onViewMode(mode) {
    const v = this.viewer;
    this.ui.root.dataset.view = mode;
    if (mode === "sky") {
      this.ensureSky().catch((err) => this.ui.toast(err.message || "Sky background unavailable"));
      this.ensureMembers();
    } else {
      this.revealMembers(false);
    }
    this.ui.layers.render();
    this.syncExtras();
    this.lens?.sync();
    v.invalidatePick();
  }

  revealMembers(on) {
    const v = this.viewer;
    const from = v.memberReveal ?? 0;
    const to = on && v.state.sky.members === "stars" ? 1 : 0;
    if (Math.abs(from - to) < 1e-3) return;
    const start = performance.now();
    const dur = 700;
    this._revealAnim?.();
    this._revealAnim = v.addAnimator((now) => {
      const t = clamp((now - start) / dur, 0, 1);
      v.memberReveal = from + (to - from) * smoothstep(0, 1, t);
      v.renderer.invalidate();
      v.invalidatePick();
      if (t >= 1) { this._revealAnim?.(); this._revealAnim = null; }
    });
    v.renderer.hold("reveal");
    setTimeout(() => v.renderer.release("reveal"), dur + 50);
  }

  setMemberMode(mode) {
    const v = this.viewer;
    v.state.sky.members = mode === "clusters" ? "clusters" : "stars";
    this.memberSeg.set(v.state.sky.members);
    if (v.state.view.mode === "sky") {
      this.ensureMembers().then(() => this.revealMembers(v.state.sky.members === "stars"));
    }
  }

  setBackground(on) {
    this.viewer.state.sky.backgroundVisible = !!on;
    this.applyLayers();
    this.changed();
  }

  syncExtras() {
    const v = this.viewer;
    const sky = v.state.view.mode === "sky";
    this.memberSeg.hidden = !sky || !this.spec.members;
    this.bgBtn.hidden = !sky;
    this.bgBtn.dataset.off = String(!v.state.sky.backgroundVisible);
    if (this.credit) this.credit.hidden = !sky || !v.state.sky.backgroundVisible;
    const on = v.state.sky.backgroundVisible;
    // The chip names the survey that dominates: the top one blended at half
    // strength or more, else the background at the bottom.
    const inView = layersInView(this.layers(), on);
    const bottom = [...inView].reverse().find((l) => (l.opacity ?? 1) >= 0.5) || inView[0];
    const name = bottom ? bottom.label || surveyLabel(bottom.survey || bottom.key) : "";
    const label = !on ? "Sky background (off)" : name ? `Sky background · ${name}` : "Sky background";
    this.bgBtn.setAttribute("aria-label", label);
    this.bgBtn.dataset.tip = `${label}\nChoose, blend or search surveys  ⇧B · B shows or hides`;
    this.bgLabel.textContent = !on ? "Sky off" : name || "Sky background";
    // The thumbnail of the background (the last one shown while it is off).
    const shown = bottom || layersInView(this.layers(), true)[0];
    const src = shown ? thumbnailUrl(shown.survey || shown.key, this.spec.thumbs) : "";
    if (src && this.bgThumb.getAttribute("src") !== src) {
      this.bgThumb.style.visibility = "";
      this.bgThumb.src = src;
    }
    this.ui.syncViewSeg?.();
  }

  // ------------------------------------------------------------ per frame

  tick(now = performance.now()) {
    const v = this.viewer;
    const r = v.renderer;
    const cam = r.camera;
    const pose = cam.pose;
    const inSky = v.state.view.mode === "sky";
    // The survey is registered to the data only when the eye is at the Sun
    // (any offset is parallax against the sky). Show it only after the
    // camera has arrived, and hide it as soon as the camera leaves.
    const atSun = pose.distance < 1e-3 && Math.hypot(pose.target[0], pose.target[1], pose.target[2]) < 1e-3;
    // The survey is today's sky: like the classic viewer it fades out over
    // the first Myr away from the present day.
    const today = 1 - smoothstep(0, 1, Math.min(Math.abs(v.timeline.time), 1));
    const target = inSky && atSun && this.sky.ready ? today : 0;
    const dt = this._lastTick ? Math.min(now - this._lastTick, 100) : 16;
    this._lastTick = now;
    if (Math.abs(this.bgFade - target) > 1e-4) {
      const rate = target > this.bgFade ? dt / 350 : dt / 120;
      this.bgFade = target > this.bgFade ? Math.min(target, this.bgFade + rate) : Math.max(target, this.bgFade - rate);
      r.hold("sky-fade");
    } else {
      this.bgFade = target;
      r.continuous.delete("sky-fade");
    }
    const op = smoothstep(0, 1, this.bgFade).toFixed(3);
    if (this.host.style.opacity !== op) this.host.style.opacity = op;
    if (inSky && this.sky.ready) {
      // A narrower window lowers the widest field Aladin can match.
      const cap = v.controls.maxFovFor("sky");
      if (pose.fov > cap) pose.fov = cap;
      // Push the pose every frame the camera is in Sky view; with the
      // Aladin loop hook this lands in the same frame as the WebGL render.
      cam.update(r.sceneRadius);
      this.sky.setView(pose.yaw * DEG, pose.pitch * DEG, cam.hfov);
      this.lens?.tick(pose.yaw * DEG, pose.pitch * DEG, cam.hfov);
    }
    if (this.members) this.updateMembers();
  }

  updateMembers() {
    const v = this.viewer;
    const m = this.members;
    const reveal = v.state.view.mode === "sky" ? v.memberReveal ?? 0 : 0;
    const t = v.timeline.time;
    const frame = v.timeline.frame;
    const zero = v.manifest.time?.zeroIndex ?? v.timeline.zeroFrame();
    const g = v.state.global;
    m.params = {
      time: t,
      reveal,
      internalMotion: true,
      starSize: (v.world.pointScale || 1) * 0.55 * g.pointSize,
      minPx: 1.1,
      glow: g.glow,
    };
    if (reveal <= 0.002) return;
    // Motion arrows: per-trace switches and scales, applied to their clusters.
    let arrowSig = "";
    for (const { key } of this.links) {
      const ts = v.state.traces[key];
      if (ts?.visible && ts.memberArrows) arrowSig += `${key}:${ts.memberArrowLength ?? 1}:${ts.memberArrowWidth ?? 1};`;
    }
    if (arrowSig !== this._arrowSig) {
      this._arrowSig = arrowSig;
      m.updateArrows((data) => {
        for (const { key, link } of this.links) {
          const ts = v.state.traces[key];
          if (!ts?.visible || !ts.memberArrows) continue;
          const L = clamp(Number(ts.memberArrowLength ?? 1) || 1, 0.1, 5), W = clamp(Number(ts.memberArrowWidth ?? 1) || 1, 0.2, 4);
          for (let i = 0; i < link.length; i++) {
            const c = link[i];
            if (c >= 0) { data[c * 4] = L; data[c * 4 + 1] = W; }
          }
        }
      });
    }
    // Parent state bits (filter / lasso) and lasso fades change without a new frame.
    let bitsVersion = 0;
    for (const { key } of this.links) {
      const b = v.points.byKey.get(key);
      bitsVersion += (b?.stateVersion || 0) + (b?.lassoVersion || 0);
    }
    const mix = v.points.lassoMix;
    const sig = `${frame.toFixed(4)}|${v._pickVersion}|${bitsVersion}|${mix[0].toFixed(3)},${mix[1].toFixed(3)}|${reveal > 0}`;
    if (sig === this._clusterSig) return;
    this._clusterSig = sig;
    const styles = v.points.params?.styles;
    m.updateClusters((state, color) => {
      state.fill(0);
      color.fill(0);
      for (const { key, link } of this.links) {
        const ts = v.state.traces[key];
        if (!ts?.visible) continue;
        const style = styles?.get(key);
        const trace = v.traceByKey.get(key);
        const d = v.data.get(key);
        const batch = v.points.byKey.get(key);
        const bits = batch?.state;
        const dim = style?.dimOpacity ?? 0.16;
        const base = parseColor(ts.color || trace.color);
        const opacityScale = ((ts.opacity ?? trace.opacity) / (trace.opacity || 1)) * g.pointOpacity * (style?.presence ?? 1);
        for (let i = 0; i < link.length; i++) {
          const c = link[i];
          if (c < 0 || state[c * 4 + 3] > 0) continue;
          // Members share their cluster's filter, lasso and birth-tree
          // state as the marker draws it: hidden or dimmed by the filter
          // (2, 1) or birth tree (32), and the lasso's level, which blends
          // while a State change plays.
          const b = bits ? bits[i] : 0;
          const w = Math.min(b & 2 ? 0 : b & 33 ? dim : 1, batch ? batch.lassoWeight(i, mix, dim) : 1);
          if (w <= 0.001) continue;
          const p = v.objectPosition(key, i, frame);
          const p0 = v.objectPosition(key, i, zero);
          if (!p || !p0) continue;
          let alpha = opacityScale * w;
          const age = d.ageNow ? d.ageNow[i] : NaN;
          if (age === age) alpha *= t >= -age ? 1 : 0;
          state[c * 4] = p[0] - p0[0];
          state[c * 4 + 1] = p[1] - p0[1];
          state[c * 4 + 2] = p[2] - p0[2];
          state[c * 4 + 3] = clamp(alpha, 0, 1);
          const rgb = d.rgb && !ts.color ? [d.rgb[i * 3] / 255, d.rgb[i * 3 + 1] / 255, d.rgb[i * 3 + 2] / 255] : base;
          color[c * 4] = rgb[0] * 255;
          color[c * 4 + 1] = rgb[1] * 255;
          color[c * 4 + 2] = rgb[2] * 255;
          const sf = d.starsFactor ? d.starsFactor[i] : 1;
          color[c * 4 + 3] = clamp(sf / 2.75, 0, 1) * 255;
        }
      }
    });
  }

  resolveMemberPick(best) {
    const m = this.members;
    if (!m) return null;
    const i = best.index;
    const v = this.viewer;
    const c = m.clusterOfStar[i];
    const pos = m.starPosition(i, v.timeline.time, true);
    const P = m.positions;
    const dist = Math.hypot(P[i * 3], P[i * 3 + 1], P[i * 3 + 2]);
    const V = m.relVelocity;
    const speed = Math.hypot(V[i * 3], V[i * 3 + 1], V[i * 3 + 2]) / 1.0227;
    const l = ((Math.atan2(P[i * 3 + 1], P[i * 3]) * DEG) + 360) % 360;
    const b = Math.asin(clamp(P[i * 3 + 2] / (dist || 1), -1, 1)) * DEG;
    const cr = m.colorData;
    return {
      kind: "member",
      trace: "__members",
      index: i,
      member: {
        cluster: String(m.clusterNames[c] || "").replace(/_/g, " "),
        position: pos, dist, l, b, speed,
        color: `rgb(${cr[c * 4]}, ${cr[c * 4 + 1]}, ${cr[c * 4 + 2]})`,
      },
    };
  }

  // ------------------------------------------------------------ layers

  layers() {
    return this.viewer.state.sky.layers || [];
  }

  /**
   * A figure (or State) with no Sky layer list still shows its survey, as
   * the classic viewer did; otherwise the background would stay blank with
   * no row in the panel to turn it on.
   */
  ensureBaseLayer() {
    const st = this.viewer.state.sky;
    if (Array.isArray(st.layers) && st.layers.length) return;
    const id = this.spec.survey || "P/DSS2/color";
    st.layers = [{ key: id, survey: id, label: surveyInfo(id).label || surveyLabel(id), opacity: 1, visible: true }];
  }

  /**
   * The survey Aladin should start with: the bottom visible layer of the
   * stack, so it does not first download (and then replace) a default.
   */
  baseSurvey() {
    const visible = this.layers().filter((l) => l.visible !== false && (l.opacity ?? 1) > 0.001);
    const bottom = visible[visible.length - 1];
    return bottom?.survey || bottom?.key || this.spec.survey;
  }

  throttledApply() {
    if (this._applyQueued) return;
    this._applyQueued = true;
    requestAnimationFrame(() => {
      this._applyQueued = false;
      this.applyLayers();
    });
  }

  applyLayers() {
    if (!this.sky.ready) return;
    this.sky.applyStack(this.layers(), { backgroundVisible: this.viewer.state.sky.backgroundVisible });
  }

  /** Something in the stack changed: refresh the panel row list, the picker and the button. */
  changed() {
    this.renderLayerList();
    this.picker?.sync();
    this.syncExtras();
  }

  setLayers(list, { quiet = false } = {}) {
    this.viewer.state.sky.layers = list;
    this.applyLayers();
    if (!quiet) this.changed();
  }

  layerByKey(key) {
    return this.layers().find((x) => x.key === key) || null;
  }

  patchLayer(key, patch) {
    const l = this.layerByKey(key);
    if (!l) return;
    Object.assign(l, patch);
    this.applyLayers();
    this.changed();
  }

  /** Show one survey (any HiPS ID or URL) alone as the background. */
  showOnly(id, meta = {}) {
    this.viewer.state.sky.backgroundVisible = true;
    this.setLayers(showOnlySurvey(this.layers(), id, meta));
  }

  /** Blend a survey over the background (half strength to start). */
  blendOnTop(id, meta = {}) {
    this.viewer.state.sky.backgroundVisible = true;
    this.setLayers(blendSurvey(this.layers(), id, meta));
  }

  hideLayer(key) {
    this.setLayers(hideSurvey(this.layers(), key));
  }

  removeLayer(key) {
    const list = this.layers().filter((l) => l.key !== key);
    if (!list.some((l) => l.visible !== false)) {
      this.ui.toast("The sky needs a background: choose another one first");
      return;
    }
    this.setLayers(list);
    this.sky.forget?.(key);
  }

  setLayerOpacity(key, x, { quiet = false } = {}) {
    const l = this.layerByKey(key);
    if (!l) return;
    l.opacity = clamp(Number(x) || 0, 0, 1);
    this.throttledApply();
    if (!quiet) this.changed();
  }

  /** Sweep the background across the spectrum (the picker's wavelength slider). */
  setWavelength(pos, stops) {
    this.viewer.state.sky.backgroundVisible = true;
    this.viewer.state.sky.layers = crossfadeToWavelength(this.layers(), stops, pos);
    this.throttledApply();
    this.picker?.sync();
    this.syncExtras();
    clearTimeout(this._listTimer);
    this._listTimer = setTimeout(() => this.renderLayerList(), 200);
  }

  /** Move a survey above (dir −1) or below (+1) the next survey in view. */
  moveLayer(key, dir) {
    const list = this.layers().slice();
    const i = list.findIndex((x) => x.key === key);
    if (i < 0) return;
    let j = i + dir;
    while (j >= 0 && j < list.length && !(list[j].visible !== false && (list[j].opacity ?? 1) > 0.001)) j += dir;
    if (j < 0 || j >= list.length) return;
    const [item] = list.splice(i, 1);
    list.splice(j, 0, item);
    this.setLayers(list);
  }

  /** A survey Aladin could not load: say so once and mark it in the picker. */
  onLayerStatus(key, status) {
    if (status === "error") {
      const l = this.layerByKey(key);
      (this._warned ||= new Set());
      if (!this._warned.has(key)) {
        this._warned.add(key);
        this.ui.toast(`Could not load ${l?.label || key}: the survey is unavailable or the ID is wrong`, { icon: icon("globe"), ms: 4200 });
      }
    }
    this.picker?.sync();
    this.renderLayerList();
  }

  /** The Layers panel's Sky section: what is in view, with the full editor. */
  layersSection() {
    if (!this.spec?.enabled) return null;
    const v = this.viewer;
    this.sectionEye = iconButton(v.state.sky.backgroundVisible ? "eye" : "eyeOff", v.state.sky.backgroundVisible ? "Hide the sky background" : "Show the sky background",
      () => this.setBackground(!v.state.sky.backgroundVisible), { shortcut: "B" });
    const add = iconButton("plus", "Choose, blend or search sky surveys", () => this.openPicker());
    const sec = h("div", { class: "ov-section ov-sky-section" },
      h("div", { class: "ov-section-head" }, h("span", { class: "ov-grow" }, "Sky background"), this.sectionEye, add));
    this.listHost = h("div");
    sec.append(this.listHost,
      h("div", { class: "ov-sky-choose" },
        h("button", { class: "ov-btn ov-btn--sm ov-btn--block", type: "button", onclick: () => this.openPicker() }, icon("stars"), "Choose background…")));
    this.renderLayerList();
    return sec;
  }

  renderLayerList() {
    if (!this.listHost || !this.listHost.isConnected) return;
    clear(this.listHost);
    const v = this.viewer;
    const st = v.state.sky;
    if (this.sectionEye) {
      this.sectionEye.replaceChildren(icon(st.backgroundVisible ? "eye" : "eyeOff"));
      const label = st.backgroundVisible ? "Hide the sky background" : "Show the sky background";
      this.sectionEye.setAttribute("aria-label", label);
      this.sectionEye.dataset.tip = `${label}  B`;
    }
    // What is in view, top first (the bottom one is the background).
    const view = layersInView(this.layers(), true).reverse();
    if (!view.length) {
      this.listHost.append(h("div", { class: "ov-note", style: { padding: "2px 14px 6px" } }, "No background: the sky is black."));
      return;
    }
    view.forEach((l, i) => {
      const isBase = i === view.length - 1;
      const key = l.key;
      const failed = this.sky.failed?.has(key);
      const badge = failed ? "Unavailable" : isBase ? "Background" : `Blend ${Math.round((l.opacity ?? 1) * 100)}%`;
      const eye = iconButton("eye", `Hide ${l.label || key}`, () => this.hideLayer(key), { cls: "ov-eye" });
      const chev = iconButton("chevron", `Adjust ${l.label || key}`, () => { this.expanded = this.expanded === key ? null : key; this.renderLayerList(); }, { cls: "ov-chev" });
      const name = h("button", { class: "ov-row-name", type: "button", title: `${l.survey || key}\nClick to adjust` }, l.label || surveyLabel(l.survey || key));
      name.addEventListener("click", () => { this.expanded = this.expanded === key ? null : key; this.renderLayerList(); });
      const thumb = h("img", { class: "ov-sky-thumb", alt: "", src: thumbnailUrl(l.survey || key, this.spec.thumbs), loading: "lazy" });
      thumb.addEventListener("error", () => thumb.remove(), { once: true });
      const row = h("div", { class: "ov-row", "data-visible": String(st.backgroundVisible), "data-expanded": String(this.expanded === key), "data-failed": String(!!failed) },
        h("span", { class: "ov-swatch-btn" }, thumb),
        name,
        h("span", { class: "ov-row-count" }, badge),
        h("div", { class: "ov-row-actions" }, eye, chev));
      this.listHost.append(row);
      if (this.expanded === key) this.listHost.append(this.layerEditor(key, i, view.length));
    });
  }

  /** Opacity, stretch, colormap, cuts and order of one survey. */
  layerEditor(key, i, n) {
    const l = this.layerByKey(key);
    const info = surveyInfo(l.survey || key);
    const stretch = l.stretch === "log10" ? "log" : l.stretch || "linear";
    const cut = (field) => numberInput({
      label: field === "cut_min" ? "Cut low" : "Cut high",
      value: Number.isFinite(l[field]) ? l[field] : "",
      onChange: (x) => { const cur = this.layerByKey(key); if (cur) { cur[field] = x; this.applyLayers(); } },
    });
    const clearCuts = h("button", { class: "ov-btn ov-btn--sm", type: "button", title: "Let Aladin choose the cuts", onclick: () => { this.patchLayer(key, { cut_min: undefined, cut_max: undefined }); } }, "Auto");
    return h("div", { class: "ov-editor" },
      slider({ label: "Opacity", min: 0, max: 1, value: l.opacity ?? 1, format: (x) => `${Math.round(x * 100)}%`, onInput: (x) => this.setLayerOpacity(key, x, { quiet: true }), onChange: () => this.changed() }),
      miniSeg({ label: "Stretch", value: stretch, options: [{ value: "linear", label: "Linear" }, { value: "log", label: "Log" }, { value: "asinh", label: "Asinh" }, { value: "sqrt", label: "Sqrt" }], onChange: (x) => this.patchLayer(key, { stretch: x }) }),
      select({ label: "Colours", value: l.colormap || "native", options: SKY_COLORMAPS.map(([value, label]) => ({ value, label })), onChange: (x) => this.patchLayer(key, { colormap: x === "native" ? undefined : x }) }),
      h("div", { class: "ov-field-row ov-sky-cuts" }, cut("cut_min"), cut("cut_max"), clearCuts),
      h("div", { class: "ov-field-row" },
        h("button", { class: "ov-btn ov-btn--sm", type: "button", disabled: i === 0, onclick: () => this.moveLayer(key, -1) }, "Bring forward"),
        h("button", { class: "ov-btn ov-btn--sm", type: "button", disabled: i === n - 1, onclick: () => this.moveLayer(key, 1) }, "Send back"),
        // Surveys added here (not the figure's own) can be taken out again.
        (this.spec.layers || []).some((x) => x.key === key) ? null
          : h("button", { class: "ov-btn ov-btn--sm ov-btn--danger", type: "button", disabled: n === 1, title: n === 1 ? "Choose another background first" : "Take this survey out of the figure", onclick: () => this.removeLayer(key) }, "Remove")),
      h("div", { class: "ov-note" }, [info.band, Number.isFinite(info.lambda) ? formatWavelength(info.lambda) : "", `HiPS ${l.survey || key}`].filter(Boolean).join(" · ")),
    );
  }

  // ------------------------------------------------------------ capture

  async captureBackground(ctx, w, h) {
    const v = this.viewer;
    ctx.fillStyle = "#000";
    ctx.fillRect(0, 0, w, h);
    if (v.state.view.mode === "sky" && this.sky.ready && v.state.sky.backgroundVisible) {
      await this.sky.drawInto(ctx, w, h);
    }
  }

  // ------------------------------------------------------------ states

  captureState() {
    return cloneJson(this.viewer.state.sky);
  }

  applyState(sky) {
    if (!sky) return;
    const v = this.viewer;
    const prevMembers = v.state.sky.members;
    v.state.sky = cloneJson(sky);
    this.ensureBaseLayer();
    this.applyLayers();
    this.changed();
    this.lens?.sync();
    this.memberSeg?.set(v.state.sky.members);
    if (prevMembers !== v.state.sky.members && v.state.view.mode === "sky") this.revealMembers(v.state.sky.members === "stars");
  }
}
