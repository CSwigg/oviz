// Sky view plugin: Aladin background, survey layer stack UI and member stars.

import { h, icon, iconButton, clear } from "../ui/dom.js";
import { slider, miniSeg } from "../ui/controls.js";
import { AladinSky } from "./aladin.js";
import { MembersLayer } from "./members.js";
import { parseColor } from "../core/color.js";
import { smoothstep, clamp } from "../core/math.js";
import { cloneJson } from "../app/state.js";

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
    this.memberTraceKeys = new Set(v.traces.filter((t) => t.points?.members).map((t) => t.key));
    v.memberTraces = this.memberTraceKeys;
    v.memberReveal = 0;
    v.skyMembersActive = () => !!this.members && v.state.sky.members === "stars";
    v.skyPickDraw = (frame) => {
      if (this.members && v.state.view.mode === "sky") this.members.drawPick(frame);
    };
    v.skyPickResolve = (best) => this.resolveMemberPick(best);
    v.renderer.beforeRender.push(() => this.tick());
    v.on("viewmode", ({ mode }) => this.onViewMode(mode));
    v.on("viewmode-settled", ({ mode }) => { if (mode === "sky") this.revealMembers(true); });
    v.on("gpu-restored", () => {
      // The member layer's GL objects died with the context: rebuild lazily.
      this.members = null;
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
    this.bgBtn = iconButton("stars", "Sky background", () => this.setBackground(!v.state.sky.backgroundVisible), { cls: "ov-glass", shortcut: "B", pressed: v.state.sky.backgroundVisible });
    ui.skyExtras.append(this.memberSeg, this.bgBtn);
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
    return [
      { title: v.state.sky.backgroundVisible ? "Hide sky background" : "Show sky background", icon: "stars", shortcut: "B", run: () => this.setBackground(!v.state.sky.backgroundVisible) },
      ...(this.spec.members ? [{ title: v.state.sky.members === "stars" ? "Show cluster markers in Sky" : "Show member stars in Sky", icon: "sparkle", run: () => this.setMemberMode(v.state.sky.members === "stars" ? "clusters" : "stars") }] : []),
      ...this.layers().map((l) => ({ title: `Sky: ${l.visible !== false ? "hide" : "show"} ${l.label || l.key}`, icon: "globe", keywords: "survey hips background", run: () => this.patchLayer(l.key, { visible: l.visible === false }) })),
    ];
  }

  onKey(e) {
    if ((e.key === "b" || e.key === "B") && !e.metaKey && !e.ctrlKey) {
      this.setBackground(!this.viewer.state.sky.backgroundVisible);
      return true;
    }
    return false;
  }

  // ------------------------------------------------------------ lifecycle

  async ensureSky() {
    if (this.sky.ready) return;
    await this.sky.init({ survey: this.spec.survey, fov: 60 });
    this.applyLayers();
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
        for (let i = 0; i < link.length; i++) states[i] = link[i] >= 0 ? 4 : 0;
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
    this.viewer.state.sky.backgroundVisible = on;
    this.applyLayers();
    this.syncExtras();
  }

  syncExtras() {
    const v = this.viewer;
    const sky = v.state.view.mode === "sky";
    this.memberSeg.hidden = !sky || !this.spec.members;
    this.bgBtn.hidden = !sky;
    this.bgBtn.setAttribute("aria-pressed", String(v.state.sky.backgroundVisible));
    this.ui.syncViewSeg?.();
  }

  // ------------------------------------------------------------ per frame

  tick() {
    const v = this.viewer;
    const cam = v.renderer.camera;
    const pose = cam.pose;
    const inSky = v.state.view.mode === "sky";
    // Registration only holds once the eye sits at the Sun; fade the
    // background in over the last few parsecs of the flight.
    const eye = cam.eye;
    const eyeDist = Math.hypot(eye[0], eye[1], eye[2]);
    const fade = inSky ? 1 - smoothstep(2, 40, eyeDist) : 0;
    const op = fade.toFixed(3);
    if (this.host.style.opacity !== op) this.host.style.opacity = op;
    if (fade > 0 && this.sky.ready) {
      this.sky.setView(pose.yaw * DEG, pose.pitch * DEG, cam.hfov);
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
    const sig = `${frame.toFixed(4)}|${v._pickVersion}|${reveal > 0}`;
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
        const base = parseColor(ts.color || trace.color);
        const opacityScale = ((ts.opacity ?? trace.opacity) / (trace.opacity || 1)) * g.pointOpacity * (style?.presence ?? 1);
        for (let i = 0; i < link.length; i++) {
          const c = link[i];
          if (c < 0 || state[c * 4 + 3] > 0) continue;
          const p = v.objectPosition(key, i, frame);
          const p0 = v.objectPosition(key, i, zero);
          if (!p || !p0) continue;
          let alpha = opacityScale;
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

  patchLayer(key, patch) {
    const l = this.layers().find((x) => x.key === key);
    if (!l) return;
    Object.assign(l, patch);
    this.applyLayers();
    this.renderLayerList();
  }

  /** Add a HiPS survey (CDS ID like "P/2MASS/color" or an https URL) on top. */
  addSurvey(raw) {
    const id = String(raw || "").trim();
    if (!id) return;
    if (!/^(https?:\/\/|[A-Za-z0-9_.-]+\/)/.test(id)) {
      this.ui.toast("Enter a HiPS ID such as P/2MASS/color, or an https:// URL");
      return;
    }
    const list = this.layers();
    const existing = list.find((l) => l.key === id || l.survey === id);
    if (existing) {
      existing.visible = true;
      existing.opacity = existing.opacity || 1;
    } else {
      const label = id.replace(/^https?:\/\/[^/]+\//, "").replace(/\//g, " ").trim();
      list.unshift({ key: id, label, survey: id, opacity: 0.7, visible: true });
    }
    this.applyLayers();
    this.renderLayerList();
    this.ui.toast(`Added ${id}`, { icon: icon("globe") });
  }

  moveLayer(key, dir) {
    const list = this.layers();
    const i = list.findIndex((x) => x.key === key);
    const j = i + dir;
    if (i < 0 || j < 0 || j >= list.length) return;
    [list[i], list[j]] = [list[j], list[i]];
    this.applyLayers();
    this.renderLayerList();
  }

  layersSection() {
    if (!this.spec?.enabled) return null;
    const sec = h("div", { class: "ov-section" });
    const groups = [{ key: "all", label: "All" }, ...(this.spec.layerGroups || []).filter((g) => g.key !== "all" && g.surveys?.length)];
    const head = h("div", { class: "ov-section-head" }, h("span", { class: "ov-grow" }, "Sky surveys"));
    sec.append(head);
    if (groups.length > 1) {
      const sel = h("select", { "aria-label": "Survey group" }, groups.map((g) => h("option", { value: g.key, selected: g.key === this.group }, g.label)));
      sel.addEventListener("change", () => { this.group = sel.value; this.renderLayerList(); });
      sec.append(h("div", { class: "ov-group", style: { marginTop: "0" } }, h("label", null, "Show"), sel));
    }
    this.listHost = h("div");
    sec.append(this.listHost);
    const input = h("input", { class: "ov-input", type: "text", placeholder: "Add survey: HiPS ID or URL", spellcheck: "false", "aria-label": "Add a HiPS survey" });
    input.addEventListener("keydown", (e) => {
      e.stopPropagation();
      if (e.key === "Enter") {
        this.addSurvey(input.value);
        input.value = "";
      }
    });
    sec.append(h("div", { style: { padding: "6px 12px 4px" } }, input));
    this.renderLayerList();
    return sec;
  }

  renderLayerList() {
    if (!this.listHost) return;
    clear(this.listHost);
    const all = this.layers();
    const g = (this.spec.layerGroups || []).find((x) => x.key === this.group);
    const shown = g?.surveys?.length ? all.filter((l) => g.surveys.includes(l.survey || l.key) || l.visible !== false) : all;
    const visibleKeys = all.filter((l) => l.visible !== false).map((l) => l.key);
    const baseKey = visibleKeys[visibleKeys.length - 1];
    for (const l of shown) {
      const vis = l.visible !== false;
      const badge = !vis ? "" : l.key === baseKey ? "Base" : "Overlay";
      const eye = iconButton(vis ? "eye" : "eyeOff", `${vis ? "Hide" : "Show"} ${l.label || l.key}`, () => this.patchLayer(l.key, { visible: !vis }), { cls: "ov-eye" });
      const chev = iconButton("chevron", `Adjust ${l.label || l.key}`, () => { this.expanded = this.expanded === l.key ? null : l.key; this.renderLayerList(); }, { cls: "ov-chev" });
      const name = h("button", { class: "ov-row-name", type: "button", title: l.survey || l.key }, l.label || l.key);
      name.addEventListener("click", () => this.patchLayer(l.key, { visible: !vis }));
      const row = h("div", { class: "ov-row", "data-visible": String(vis), "data-expanded": String(this.expanded === l.key) },
        h("span", { class: "ov-swatch-btn" }, icon("globe")),
        name,
        h("span", { class: "ov-row-count" }, badge),
        h("div", { class: "ov-row-actions" }, eye, chev));
      this.listHost.append(row);
      if (this.expanded === l.key) {
        const ed = h("div", { class: "ov-editor" },
          slider({ label: "Opacity", min: 0, max: 1, value: l.opacity ?? 1, format: (x) => `${Math.round(x * 100)}%`, onInput: (x) => { l.opacity = x; this.applyLayers(); } }),
          h("div", { class: "ov-field-row" },
            h("button", { class: "ov-btn ov-btn--sm", type: "button", onclick: () => this.moveLayer(l.key, -1) }, "Move up"),
            h("button", { class: "ov-btn ov-btn--sm", type: "button", onclick: () => this.moveLayer(l.key, 1) }, "Move down")),
          h("div", { class: "ov-note" }, `HiPS ${l.survey || l.key}`));
        this.listHost.append(ed);
      }
    }
  }

  // ------------------------------------------------------------ capture

  async captureBackground(ctx, w, h) {
    const v = this.viewer;
    ctx.fillStyle = "#04060a";
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
    this.applyLayers();
    this.renderLayerList();
    this.memberSeg?.set(v.state.sky.members);
    if (prevMembers !== v.state.sky.members && v.state.view.mode === "sky") this.revealMembers(v.state.sky.members === "stars");
    this.syncExtras();
  }
}
