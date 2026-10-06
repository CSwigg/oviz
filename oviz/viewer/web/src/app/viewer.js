// The Viewer: owns the renderer, layers, time and view state, and exposes
// the imperative API that the UI, States and the public JS API drive.

import { Emitter } from "../core/emitter.js";
import { clamp, v3dist, niceFloor } from "../core/math.js";
import { Renderer } from "../engine/renderer.js";
import { Controls } from "../engine/controls.js";
import { RenderTarget } from "../engine/gl.js";
import { makePose, clonePose, poseFromEyeTarget, poseTween, forwardFromAngles, anglesFromForward, sanitizePose } from "../engine/camera.js";
import { LSR_POINT, normalizeAnchor, sameAnchor, inferAnchor } from "../engine/anchor.js";
import { cpuFramePosition, frameOffset, trailParams } from "../engine/frames.js";
import { rayBox, volumeMedianDepth, rayPlaneRect, screenSegment } from "../engine/surface.js";
import { PointsLayer, resetStarTextures } from "../layers/points.js";
import { LinesLayer } from "../layers/lines.js";
import { FlowLayer, flowStyle } from "../layers/flow.js";
import { ImagesLayer } from "../layers/images.js";
import { VolumesLayer } from "../layers/volumes.js";
import { buildEventKde, kdeSampleTime, kdeWindow, kdeId, kdeSampleTimes } from "../layers/kde.js";
import { LabelsOverlay } from "../layers/labels.js";
import { Timeline, presentDayFade } from "./timeline.js";
import { initialViewerState, clampVolumeState, volumeColormap, volumeColorField } from "./state.js";

const SKY_FOV = 60;

export class Viewer extends Emitter {
  constructor(container, manifest, store) {
    super();
    this.manifest = manifest;
    this.store = store;
    this.container = container;
    this.canvas = document.createElement("canvas");
    this.canvas.className = "ov-canvas";
    this.canvas.setAttribute("aria-label", manifest.title ? `Interactive figure: ${manifest.title}` : "Interactive figure");
    this.canvas.tabIndex = -1;
    container.appendChild(this.canvas);
    this.renderer = new Renderer(this.canvas, { maxDpr: isMobile() ? 1.5 : 2 });
    this.gl = this.renderer.gl;
    const time = manifest.time || { values: [0] };
    this.timeline = new Timeline(time.values, { initialIndex: time.initialIndex, intervalMs: time.playbackIntervalMs });
    this.traces = manifest.traces || [];
    this.traceByKey = new Map(this.traces.map((t) => [t.key, t]));
    this.volumeSpecs = manifest.volumes || [];
    this.world = manifest.world || {};
    this.renderer.sceneRadius = Math.max(this.world.maxSpan || 1e4, 1000);
    this.state = initialViewerState(manifest);
    // The initial state may open on another frame (classic exports do).
    const openFrame = Number(this.state.time.frame);
    if (Number.isFinite(openFrame) && Math.abs(openFrame - this.timeline.frame) > 1e-9) this.timeline.setFrame(openFrame, { silent: true });
    this.hover = null;
    this.selection = [];
    this.tween = null;
    this.data = new Map(); // trace key → decoded arrays
    this.meta = new Map(); // trace key → metadata columns (lazy)
    this.pickIds = [];
    this.pickTarget = new RenderTarget(this.gl, { depth: true, filter: this.gl.NEAREST });
    this._pickVersion = 0;
    this.extraStyles = new Map(); // synthetic overlays (trails, measurements)
    this.stateExtensions = new Map(); // name → {capture, apply} for States
    this._buildLayers();
    this.controls = new Controls(this.canvas, this.renderer, {
      onChange: () => this._onCameraChange(),
      // Grabbing the camera ends any flight or State transition.
      onInteractStart: () => { this._cancelTween(); this.emit("camera-takeover", {}); this.emit("user-camera", {}); },
      anchored: () => this.anchored,
    });
    this.controls.maxDistance = Math.max(this.world.maxSpan || 1e4, 1e3) * 8;
    this.labels = new LabelsOverlay(container);
    this.renderer.overlays.push(this.labels);
    this.renderer.beforeRender.push((now) => this._beforeRender(now));
    this.renderer.onContextRestored = () => this._restoreGPU();
    this.timeline.on(() => {
      this._trackAnchor();
      this.renderer.invalidate();
      this.emit("time", this.timeline);
    });
    this._applyInitialCamera();
  }

  // -------------------------------------------------------------- setup

  _buildLayers() {
    const gl = this.gl;
    const frames = this.timeline.count;
    this.images = this.renderer.add(new ImagesLayer(gl));
    this.volumes = this.renderer.add(new VolumesLayer(gl, { caps: this.renderer.caps }));
    this.lines = this.renderer.add(new LinesLayer(gl, { frames }));
    this.flows = this.renderer.add(new FlowLayer(gl, { cmaps: this.volumes.cmaps }));
    this.points = this.renderer.add(new PointsLayer(gl, { frames }));
  }

  /** Upload every trace whose blobs are loaded (call after critical load). */
  async attachTraces({ labels = true } = {}) {
    const store = this.store;
    const cm = this.manifest.colormaps || {};
    for (const [name, entry] of Object.entries(cm)) {
      const lut = await store.get(entry.blob);
      this.points.setColormap(name, lut, entry.width);
      this.volumes.addColormap(name, lut, entry.width);
    }
    let pickId = 1;
    for (const trace of this.traces) {
      const d = {};
      if (trace.points) {
        const p = trace.points;
        d.position = await store.get(p.position.blob);
        if (p.position.offset) d.offset = await store.get(p.position.offset);
        if (p.size?.blob) d.size = await store.get(p.size.blob);
        if (p.sizeMax?.blob) d.sizeMax = await store.get(p.sizeMax.blob);
        if (p.opacity?.blob) d.opacity = await store.get(p.opacity.blob);
        if (p.colorScalar?.blob) d.scalar = await store.get(p.colorScalar.blob);
        if (p.rgb?.blob) d.rgb = await store.get(p.rgb.blob);
        if (p.symbol?.blob) d.symbol = await store.get(p.symbol.blob);
        if (p.ageNow?.blob) d.ageNow = await store.get(p.ageNow.blob);
        if (p.starsFactor?.blob) d.starsFactor = await store.get(p.starsFactor.blob);
        if (p.nStars?.blob) d.nStars = await store.get(p.nStars.blob);
        this.pickIds[pickId] = trace.key;
        this.points.addTrace(trace, d, pickId++);
        this.data.set(trace.key, d);
      }
      if (trace.lines) {
        const L = trace.lines;
        const ld = { position: await store.get(L.position.blob) };
        if (L.position.offset) ld.offset = await store.get(L.position.offset);
        if (L.arc?.blob) ld.arc = await store.get(L.arc.blob);
        if (L.segmentColors?.blob) ld.segmentColors = await store.get(L.segmentColors.blob);
        this.lines.addTrace(trace, ld);
        (this.lineData ||= new Map()).set(trace.key, ld);
      }
      if (trace.flow) {
        const F = trace.flow;
        this.flows.addTrace(trace, { position: await store.get(F.position.blob), attr: await store.get(F.attr.blob) });
      }
      if (trace.labels && labels) {
        const L = trace.labels;
        const ld = { position: await store.get(L.position.blob) };
        if (L.position.offset) ld.offset = await store.get(L.position.offset);
        this.labels.addTrace(trace, ld);
      }
    }
    this._settleInitialAnchor();
    this.renderer.invalidate();
  }

  /**
   * Rebuild every GPU resource after a lost WebGL context is restored.
   * Decoded arrays stay cached in the bundle store, so this re-uploads
   * without re-decoding; plugins rebuild their own layers on "gpu-restored".
   */
  async _restoreGPU() {
    // Loaders still streaming into the dead context's layers must stop.
    this._gpuGeneration = (this._gpuGeneration || 0) + 1;
    resetStarTextures();
    this.data.clear();
    this.lineData?.clear();
    this.pickIds = [];
    this.pickTarget = new RenderTarget(this.gl, { depth: true, filter: this.gl.NEAREST });
    this._pickVersion++;
    this.extraStyles.clear();
    this._buildLayers();
    this.renderer.overlays = this.renderer.overlays.filter((o) => o === this.labels);
    await this.attachTraces({ labels: false });
    this.emit("gpu-restored", {});
    this.attachImages();
    this.attachVolumes();
    this.renderer.invalidate();
  }

  async attachImages() {
    const gen = this._gpuGeneration || 0;
    const layer = this.images;
    for (const img of this.manifest.images || []) {
      try {
        const blob = await this.store.get(img.image.blob);
        const centers = img.center?.blob ? await this.store.get(img.center.blob) : null;
        const opac = img.opacity?.blob ? await this.store.get(img.opacity.blob) : null;
        if (gen !== (this._gpuGeneration || 0) || layer !== this.images) return;
        await layer.add(img, blob, centers, opac);
        this.renderer.invalidate();
      } catch (err) {
        console.warn("Oviz: image plane failed to load", img.key, err);
      }
    }
  }

  async attachVolumes(onEach) {
    // Smallest first so something appears quickly.
    const specs = [...this.volumeSpecs].sort((a, b) => a.dims[0] * a.dims[1] * a.dims[2] - b.dims[0] * b.dims[1] * b.dims[2]);
    const gen = this._gpuGeneration || 0;
    const layer = this.volumes;
    for (const spec of specs) {
      try {
        if (spec.procedural) {
          await this._attachProceduralVolume(spec, gen, layer);
          onEach?.(spec);
          this.emit("volume-ready", spec);
          continue;
        }
        const data = await this.store.get(spec.data.blob);
        const color = spec.colorData?.blob ? await this.store.get(spec.colorData.blob) : null;
        const fields = [];
        for (const f of spec.colorFields || []) fields.push(await this.store.get(f.data.blob));
        let occ = null;
        if (spec.occupancy?.blob) {
          const o = await this.store.get(spec.occupancy.blob);
          occ = { data: o, gx: spec.occupancy.dims[0], gy: spec.occupancy.dims[1], gz: spec.occupancy.dims[2] };
        }
        if (gen !== (this._gpuGeneration || 0) || layer !== this.volumes) return;
        layer.addVolume(spec, data, occ, color, "", fields);
        this.renderer.invalidate();
        onEach?.(spec);
        this.emit("volume-ready", spec);
        await new Promise((r) => setTimeout(r, 0));
      } catch (err) {
        console.warn("Oviz: volume failed to load", spec.key, err);
      }
    }
  }

  // ---------------------------------------------- rebuilt event densities

  /**
   * Load an event-KDE volume: its events, and the precomputed density for
   * its initial time when the bundle has one (the classic first frame).
   */
  async _attachProceduralVolume(spec, gen, layer) {
    const proc = spec.procedural;
    const events = await this.store.get(proc.events.blob);
    const initial = spec.data?.blob ? await this.store.get(spec.data.blob) : null;
    if (gen !== (this._gpuGeneration || 0) || layer !== this.volumes) return;
    (this._kde ||= new Map()).set(spec.key, { spec, events });
    const vs = this.state.volumes[spec.stateKey];
    const time = kdeSampleTime(proc, this.timeline.time);
    const win = kdeWindow(proc, vs);
    const initialId = initial && proc.initialTimeMyr != null ? kdeId(kdeSampleTime(proc, proc.initialTimeMyr), proc.initialWindowMyr ?? proc.windowMyr) : "";
    if (initialId) layer.addVolume(spec, initial, null, null, initialId);
    const id = kdeId(time, win);
    if (id !== initialId) {
      const built = this._buildKde(spec.key, time, win);
      if (initialId) layer.useData(spec.key, id, built);
      else layer.addVolume(spec, built, null, null, id);
    }
    this.renderer.invalidate();
    this._prefetchKde();
  }

  _buildKde(key, time, win) {
    const k = this._kde.get(key);
    const t0 = performance.now();
    const { data } = buildEventKde(k.events, k.spec.procedural, k.spec.dims, k.spec.bounds, time, win);
    this.kdeLastBuildMs = performance.now() - t0;
    this.canvas.dataset.kdeLastBuildMs = this.kdeLastBuildMs.toFixed(1);
    return data;
  }

  /** Swap an event-KDE volume to `time`'s sample, building it if needed; its cache id. */
  _ensureKde(spec, time, vs) {
    const k = this._kde?.get(spec.key);
    if (!k) return "";
    const proc = spec.procedural;
    const t = kdeSampleTime(proc, time), win = kdeWindow(proc, vs);
    const id = kdeId(t, win);
    if (!this.volumes.useData(spec.key, id)) {
      this.volumes.useData(spec.key, id, this._buildKde(spec.key, t, win));
      this._prefetchKde();
    }
    return id;
  }

  /**
   * While the page is idle, build every sample of the visible event-KDE
   * volumes nearest the current time first, so scrubbing only binds
   * textures. Phones skip it (memory), as the classic runtime did.
   */
  _prefetchKde() {
    if (this._kdePrefetching || isMobile() || !this._kde?.size) return;
    this._kdePrefetching = true;
    const idle = window.requestIdleCallback || ((fn) => setTimeout(fn, 200));
    const step = () => {
      this._kdePrefetching = false;
      for (const { spec } of this._kde.values()) {
        const vs = this.state.volumes[spec.stateKey];
        if (!vs?.visible || !this.volumes.volumes.has(spec.key)) continue;
        const proc = spec.procedural, win = kdeWindow(proc, vs);
        const now = this.timeline.time;
        const todo = kdeSampleTimes(proc, this.timeline.times)
          .filter((t) => !this.volumes.hasData(spec.key, kdeId(t, win)))
          .sort((a, b) => Math.abs(a - now) - Math.abs(b - now));
        if (!todo.length) continue;
        const keep = this.volumes.volumes.get(spec.key).dataId;
        this.volumes.useData(spec.key, kdeId(todo[0], win), this._buildKde(spec.key, todo[0], win));
        // Building caches without showing it: restore the drawn sample.
        this.volumes.useData(spec.key, keep);
        this._kdePrefetching = true;
        idle(step, { timeout: 4000 });
        return;
      }
    };
    idle(step, { timeout: 4000 });
  }

  _applyInitialCamera() {
    const cam = this.manifest.camera || {};
    const pos = cam.position || [0, -800, 1500];
    const tgt = cam.target || this.world.center || [0, 0, 0];
    const pose = poseFromEyeTarget(pos, tgt, cam.fov || 60);
    // The figure's camera anchor: the one it asks for, else the LSR when the
    // home view orbits it (Oviz scenes are LSR-centred), else none.
    // An explicit anchor, else the classic focus group (a layer's median).
    const focus = this.manifest.animation?.focusTrace;
    const asked = normalizeAnchor(this.manifest.viewer?.cameraAnchor) || (focus ? { kind: "layer", trace: String(focus) } : null);
    this.defaultAnchor = asked && this.anchorAvailable(asked) ? asked : this.inferAnchor(pose, { kind: "free" });
    if (this.defaultAnchor.kind === "lsr") pose.target = LSR_POINT.slice();
    this.state.view.anchor = { ...this.defaultAnchor };
    this.homePose = clonePose(pose);
    this.renderer.camera.pose = pose;
    this.state.view.pose = pose;
  }

  /** A Sun or object anchor needs its track: land the opening view on it. */
  _settleInitialAnchor() {
    const a = this.state.view.anchor;
    if (this._anchorSettled || !a || a.kind === "free" || a.kind === "lsr") return;
    const p = this.anchorPoint(a);
    if (!p) return;
    this._anchorSettled = true;
    this.homePose.target = p.slice();
    if (this.state.view.mode !== "sky" && !this.tween) {
      this.renderer.camera.pose.target = p.slice();
      this._onCameraChange();
    }
    this._anchorLast = p;
  }

  start() {
    this.renderer.start();
  }

  // -------------------------------------------------------------- per frame

  _beforeRender(now) {
    const tl = this.timeline;
    if (tl.playing) {
      tl.tick(now);
      this.renderer.hold("playback");
    } else {
      this.renderer.continuous.delete("playback");
    }
    if (this.tween) {
      const { pose, done } = this.tween.step(now);
      Object.assign(this.renderer.camera.pose, clonePose(pose));
      this.renderer.markInteraction();
      if (done) {
        const t = this.tween;
        this.tween = null;
        this.renderer.continuous.delete("tween");
        t.settle?.(true);
      }
      this._onCameraChange(true);
    }
    for (const fn of this._animators || []) fn(now);
    this._resolveParams();
    this._stepFlow(now);
  }

  // -------------------------------------------------------------- flow lines

  /** Whether the figure has animated flow layers. */
  get hasFlows() {
    return this.traces.some((t) => t.flow);
  }

  /**
   * Flow pulses run unless paused. The default follows the figure (and the
   * reader's reduced-motion preference); an explicit choice is kept in
   * `state.global.flowPaused`, so States restore it.
   */
  get flowPlaying() {
    const g = this.state.global;
    if (typeof g.flowPaused === "boolean") return !g.flowPaused;
    if (typeof window !== "undefined" && window.matchMedia?.("(prefers-reduced-motion: reduce)").matches) return false;
    return this.traces.some((t) => t.flow && t.flow.animation?.playing !== false);
  }

  setFlowPlaying(on) {
    this.state.global.flowPaused = !on;
    this.renderer.invalidate();
    this.emit("global", { flowPaused: !on });
  }

  /**
   * Advance the flow clocks, and keep frames coming only while a shown flow
   * is playing: a paused or hidden flow leaves the figure idle.
   */
  _stepFlow(now) {
    const layer = this.flows;
    if (!layer?.batches.length) return;
    const styles = layer.params?.styles;
    const running = this.flowPlaying && layer.anyVisible(styles);
    const last = this._flowLast;
    this._flowLast = running ? now : null;
    if (running && last != null) layer.advance(Math.min(Math.max(now - last, 0), 100) / 1000, styles);
    if (running) this.renderer.hold("flow");
    else this.renderer.continuous.delete("flow");
  }

  addAnimator(fn) {
    (this._animators ||= new Set()).add(fn);
    return () => this._animators.delete(fn);
  }

  _presenceWeight(trace, frame) {
    const pres = trace.presence;
    if (!pres) return 1;
    const set = (this._presenceSets ||= new Map());
    let s = set.get(trace.key);
    if (!s) set.set(trace.key, (s = new Set(pres)));
    const i = Math.floor(frame), j = Math.min(i + 1, this.timeline.count - 1), t = frame - i;
    return (s.has(i) ? 1 : 0) * (1 - t) + (s.has(j) ? 1 : 0) * t;
  }

  _resolveParams() {
    const st = this.state;
    const g = st.global;
    const frame = this.timeline.frame;
    const time = this.timeline.time;
    const styles = new Map();
    const skyMode = st.view.mode === "sky";
    const memberMode = skyMode && this.skyMembersActive?.();
    for (const trace of this.traces) {
      const ts = st.traces[trace.key] || {};
      // A legend-less trace keyed like a volume group (its outlines and
      // labels) shows with that volume, not with the grid.
      const companion = !trace.showInLegend ? st.volumes[trace.key] : null;
      const isReference = !trace.showInLegend && !companion;
      let visible = companion ? companion.visible === true && !skyMode : ts.visible !== false;
      if (isReference) visible = visible && g.grid !== false && !skyMode;
      // Sky view looks out from the Sun: it is not drawn there, nor is the
      // 3D-only context a figure marks (a shell around the Sun, say).
      if (skyMode && (trace.key === this.world.sunTrace || trace.skyHidden)) visible = false;
      const defaultOpacity = trace.opacity > 1e-6 ? trace.opacity : 1;
      let starsExp = 0;
      if (trace.hasStars) starsExp = g.sizeByStars ? (trace.sizeByStarsDefault ? 0 : 1) : (trace.sizeByStarsDefault ? -1 : 0);
      const hoverIdx = this.hover && this.hover.trace === trace.key ? this.hover.index : -1;
      let opacityScale = (ts.opacity ?? trace.opacity) / defaultOpacity;
      if (isReference) opacityScale *= g.gridOpacity ?? 1;
      // The longitude grid is drawn from today's Sun: it fades out within a
      // Myr of the present day (classic).
      if (trace.presentDayOnly) {
        opacityScale *= presentDayFade(time);
        visible = visible && opacityScale > 1e-3;
      }
      styles.set(trace.key, {
        visible,
        presence: this._presenceWeight(trace, frame),
        opacityScale,
        sizeScale: ts.sizeScale ?? 1,
        starsExp,
        colorMode: ts.colorMode || trace.colorBy?.defaultMode || "fixed",
        colormap: ts.colormap || trace.colorBy?.colormap,
        cmin: ts.cmin ?? trace.colorBy?.cmin ?? 0,
        cmax: ts.cmax ?? trace.colorBy?.cmax ?? 1,
        colorOverride: ts.color || null,
        hover: hoverIdx,
        dimOpacity: 0.16,
        // Samples of a surface or model are not picked (no hover, no click).
        pickable: trace.pickable !== false,
        flow: trace.flow ? flowStyle(trace, ts) : null,
      });
    }
    for (const [k, s] of this.extraStyles) styles.set(k, s);
    const pointParams = {
      styles,
      frame,
      time,
      pointScale: this.world.pointScale || 1,
      pointSize: g.pointSize,
      pointOpacity: g.pointOpacity,
      glow: g.glow,
      fadeTime: g.fadeTime,
      fadeInOut: g.fadeInOut,
      fadeByOpacity: g.fadeByOpacity,
      // Markers of clusters with member stars crossfade to the stars in Sky.
      memberFade: memberMode ? this.memberReveal ?? 0 : 0,
      trails: this._trailParams(g.trails),
    };
    this.points.params = pointParams;
    this.lines.params = { styles, frame, lineOpacity: 1 };
    this.flows.params = { styles, focusFallback: Math.max((this.world.maxSpan || 1000) * 0.5, 100) };
    this.labels.params = { styles, frame, labelSignature: this._labelSignature(), labelsVisible: g.labels !== false };
    // Images fade with the scale bar (zoom), like the legacy figure.
    const scaleBarPc = this.scaleBarPc();
    const imageStates = new Map();
    for (const img of this.manifest.images || []) {
      const s = st.images[img.key] || {};
      imageStates.set(img.key, { visible: s.visible !== false && !skyMode, opacity: s.opacity ?? 1 });
    }
    this.images.params = { frame, images: imageStates, scaleBarPc, galacticOpacity: 1 };
    this.volumes.params = this._volumeParams(frame, time);
  }

  _labelSignature() {
    let s = "";
    for (const trace of this.traces) {
      if (!trace.labels) continue;
      const ts = this.state.traces[trace.key] || {};
      s += `${trace.key}:${ts.visible !== false}:${this.state.volumes[trace.key]?.visible}:${ts.opacity}:${ts.sizeScale};`;
    }
    return s + this.state.global.grid + this.state.global.labels;
  }

  _volumeParams(frame, time) {
    const st = this.state;
    const draws = [];
    let sig = "";
    const T = this.timeline;
    const dtFrame = T.count > 1 ? Math.abs(T.frameToTime(Math.min(Math.floor(frame) + 1, T.count - 1)) - T.frameToTime(Math.floor(frame))) || 1 : 1;
    const groups = new Map();
    for (const spec of this.volumeSpecs) {
      if (!groups.has(spec.stateKey)) groups.set(spec.stateKey, []);
      groups.get(spec.stateKey).push(spec);
    }
    for (const [stateKey, specs] of groups) {
      const vs = st.volumes[stateKey];
      if (!vs || !vs.visible) continue;
      const timed = specs.filter((s) => s.timeMyr != null);
      const weights = [];
      if (timed.length) {
        // Crossfade the two layers bracketing the current time.
        const sorted = timed.slice().sort((a, b) => a.timeMyr - b.timeMyr);
        let lo = null, hi = null;
        for (const s of sorted) {
          if (s.timeMyr <= time + 1e-9) lo = s;
          if (s.timeMyr >= time - 1e-9 && !hi) hi = s;
        }
        if (lo && hi && lo !== hi) {
          const span = hi.timeMyr - lo.timeMyr || 1;
          const t = (time - lo.timeMyr) / span;
          if (span <= dtFrame * 1.5) {
            weights.push([lo, 1 - t], [hi, t]);
          } else {
            if (Math.abs(time - lo.timeMyr) < dtFrame) weights.push([lo, 1 - Math.abs(time - lo.timeMyr) / dtFrame]);
            if (Math.abs(time - hi.timeMyr) < dtFrame) weights.push([hi, 1 - Math.abs(time - hi.timeMyr) / dtFrame]);
          }
        } else if (lo || hi) {
          const s = lo || hi;
          const d = Math.abs(time - s.timeMyr);
          if (d < dtFrame) weights.push([s, 1 - d / dtFrame]);
        }
      }
      for (const spec of specs.filter((s) => s.timeMyr == null)) {
        let w = 1;
        const allTimes = vs.showAllTimes && spec.supportsShowAllTimes;
        if (!allTimes && spec.onlyAtT0 !== false) {
          w = clamp(1 - Math.abs(time) / dtFrame, 0, 1);
        }
        weights.push([spec, w]);
      }
      for (const [spec, w0] of weights) {
        // Sky view shows volumes only at the present day, like the classic
        // viewer: the sky is seen from where the Sun is today.
        const w = st.view.mode === "sky" ? w0 * clamp(1 - Math.abs(time) / dtFrame, 0, 1) : w0;
        if (w <= 0.002) continue;
        const vsClamped = clampVolumeState(spec, vs);
        // Members of one group keep their own colours (e.g. each cloud)
        // while the group shows one of its members' own colormaps.
        if (specs.length > 1 && spec.defaults?.colormap && specs.some((m) => m.defaults?.colormap === vs.colormap)) vsClamped.colormap = spec.defaults.colormap;
        const [dmin, dmax] = spec.dataRange;
        const span = dmax - dmin;
        let low = span > 0 ? clamp((vsClamped.vmin - dmin) / span, 0, 1) : 0;
        let high = span > 0 ? clamp((vsClamped.vmax - dmin) / span, 0, 1) : 1;
        if (!(high > low)) high = Math.min(1, low + 1e-6);
        const allTimes = vs.showAllTimes && spec.supportsShowAllTimes && spec.timeMyr == null;
        const angle = allTimes && spec.coRotate ? (spec.coRotationRate || 0) * (time - (spec.referenceTimeMyr || 0)) : 0;
        const dataId = spec.procedural ? this._ensureKde(spec, time, vs) : "";
        // A colour field (a KT map's velocity) colours it unless density is chosen.
        const field = volumeColorField(spec, vsClamped);
        const fieldIndex = field ? spec.colorFields.indexOf(field) : -1;
        const colormap = volumeColormap(spec, vsClamped);
        draws.push({ key: spec.key, state: vsClamped, low, high, fade: w, angle, canSkip: true, fieldIndex, colormap });
        sig += `${spec.key}:${dataId}:${w.toFixed(3)}:${low}:${high}:${vsClamped.opacity}:${vsClamped.steps}:${vsClamped.alphaCoef}:${vsClamped.stretch}:${colormap}:${fieldIndex}:${angle.toFixed(5)}:${vsClamped.lightingMode};`;
      }
    }
    return { frame, volumeDraws: draws, volumeSignature: sig };
  }

  // -------------------------------------------------------------- camera

  get pose() {
    return this.renderer.camera.pose;
  }

  _onCameraChange(fromTween = false) {
    this.state.view.pose = this.renderer.camera.pose;
    this.emit("camera", { fromTween });
  }

  /** Stop a running camera flight; its promise resolves with {done: false}. */
  _cancelTween() {
    const t = this.tween;
    if (!t) return;
    this.tween = null;
    this.renderer.continuous.delete("tween");
    this.renderer.invalidate();
    t.settle?.(false);
  }

  /**
   * Fly the camera to `pose`. Resolves {done: true} on arrival, or
   * {done: false} if the flight was interrupted (user input, a new flight).
   * The camera keeps its anchor: the orbit centre may sit off it, and the
   * camera keeps that offset as the anchor moves through time.
   */
  animateTo(pose, { duration = 1100, onDone, ease } = {}) {
    this.emit("camera-takeover", {});
    this._cancelTween();
    this.controls.stop();
    const d = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches ? 0 : duration;
    const tween = poseTween(this.renderer.camera.pose, sanitizePose(clonePose(pose)), d, ease);
    this.tween = tween;
    this.renderer.hold("tween");
    return new Promise((resolve) => {
      tween.settle = (done) => {
        tween.settle = null;
        if (done) onDone?.();
        resolve({ done });
      };
    });
  }

  setPose(pose) {
    this._cancelTween();
    Object.assign(this.renderer.camera.pose, sanitizePose(clonePose(pose)));
    this._anchorLast = this.anchorPoint();
    this._onCameraChange();
    this.renderer.invalidate();
  }

  /** Fly home: the opening view, anchored the way the figure opens. */
  resetView() {
    // Resetting is the reader taking the camera (it stops auto-orbit).
    this.emit("user-camera", { reset: true });
    if (this.state.view.mode === "sky") {
      const p = clonePose(this.pose);
      p.fov = SKY_FOV;
      return this.animateTo(p, { duration: 700 });
    }
    this.emit("camera-takeover", {});
    this._setAnchorState(this.defaultAnchor);
    return this.animateTo(this.homeView(), { duration: 1000 });
  }

  /** The home pose, centred on the default anchor where it is now. */
  homeView() {
    const p = clonePose(this.homePose);
    const at = this.anchorPoint(this.defaultAnchor);
    if (at) p.target = at;
    return p;
  }

  // -------------------------------------------------------------- camera anchor

  /** Whether the scene's origin is the LSR (figures centred on a cluster group are not). */
  get hasLsr() {
    return this.manifest.viewer?.lsrOrigin !== false;
  }

  /** How close the orbit centre must be to count as on an anchor (pc). */
  get anchorTolerance() {
    return Math.max((this.world.maxSpan || 1000) * 1e-4, 0.05);
  }

  get anchor() {
    return this.state.view.anchor || { kind: "free" };
  }

  /** The orbit centre is held on an anchor (never in Sky view, where the eye is at the Sun). */
  get anchored() {
    return this.state.view.mode !== "sky" && this.anchor.kind !== "free";
  }

  anchorAvailable(anchor) {
    const a = normalizeAnchor(anchor);
    if (!a) return false;
    if (a.kind === "lsr") return this.hasLsr;
    if (a.kind === "sun") return !!this.traceByKey.get(this.world.sunTrace)?.points;
    if (a.kind === "object") {
      const t = this.traceByKey.get(a.trace);
      return !!t?.points && a.index < t.points.count;
    }
    if (a.kind === "layer") return !!this.traceByKey.get(a.trace)?.points?.count;
    return true;
  }

  /** The anchor a pose implies (for States and links saved without one). */
  inferAnchor(pose, fallback = this.defaultAnchor || { kind: "free" }) {
    return inferAnchor(pose, this.hasLsr ? LSR_POINT : null, this.anchorTolerance, fallback);
  }

  /** Where an anchor is at a (fractional) frame; null for "free". */
  anchorPoint(anchor = this.anchor, frame = this.timeline.frame) {
    const a = normalizeAnchor(anchor);
    if (!a) return null;
    if (a.kind === "lsr") return this.hasLsr ? LSR_POINT.slice() : null;
    if (a.kind === "sun") return this.world.sunTrace ? this.objectPosition(this.world.sunTrace, 0, frame) : null;
    if (a.kind === "object") return this.objectPosition(a.trace, a.index, frame);
    if (a.kind === "layer") return this.layerMedian(a.trace, frame);
    return null;
  }

  /**
   * Anchor the camera (see engine/anchor.js). It flies its orbit centre onto
   * the anchor, keeping the viewing angle (and `distance`, if given), then
   * orbits, zooms and moves through time with it.
   */
  setAnchor(anchor, { fly = true, distance, duration = 900 } = {}) {
    const a = normalizeAnchor(anchor);
    if (!a || !this.anchorAvailable(a)) return Promise.resolve({ done: false });
    // A State transition in flight hands the camera over first.
    this.emit("camera-takeover", {});
    this._setAnchorState(a);
    const at = this._anchorLast;
    if (!fly || !at || this.state.view.mode === "sky") return Promise.resolve({ done: true });
    const p = clonePose(this.pose);
    p.target = at.slice();
    if (distance != null) p.distance = distance;
    return this.animateTo(p, { duration });
  }

  /** Let go of the anchor ("free"): the camera stays put and moves with nothing. */
  releaseAnchor() {
    if (this.anchor.kind === "free") return;
    this._setAnchorState({ kind: "free" }, { released: true });
  }

  _setAnchorState(anchor, { released = false } = {}) {
    const a = normalizeAnchor(anchor) || { kind: "free" };
    const changed = !sameAnchor(a, this.state.view.anchor);
    this.state.view.anchor = a;
    this._anchorLast = this.anchorPoint(a);
    if (changed) this.emit("anchor", { anchor: a, released });
  }

  /** After the camera was placed by other means: track the anchor from where it is now. */
  reconcileAnchor() {
    this._anchorLast = this.anchorPoint();
  }

  /**
   * Keep the camera with its anchor as time changes: move it by the anchor's
   * displacement (never snap), so a camera placed near the anchor keeps its
   * offset and nothing jumps. The LSR is fixed in the scene, so an
   * LSR-anchored camera stays put. In Sky view the pose to return to moves.
   */
  _trackAnchor() {
    const a = this.anchor;
    if (a.kind === "free" || a.kind === "lsr") return;
    const p = this.anchorPoint(a);
    if (!p) return; // not present at this time: catch up when it reappears
    const last = this._anchorLast;
    this._anchorLast = p;
    if (!last) return;
    const d = [p[0] - last[0], p[1] - last[1], p[2] - last[2]];
    if (!d[0] && !d[1] && !d[2]) return;
    if (this.state.view.mode === "sky") {
      if (this.returnPose) for (let i = 0; i < 3; i++) this.returnPose.target[i] += d[i];
      return;
    }
    if (this.tween) {
      this.tween.shiftTarget?.(d);
      return;
    }
    const t = this.renderer.camera.pose.target;
    t[0] += d[0]; t[1] += d[1]; t[2] += d[2];
    this._onCameraChange();
  }

  /** Current world position of an object (interpolated at the current time). */
  objectPosition(traceKey, index, frame = this.timeline.frame) {
    const trace = this.traceByKey.get(traceKey);
    const d = this.data.get(traceKey);
    if (!trace || !d || !trace.points) return null;
    const p = cpuFramePosition(d.position, trace.points.position.frames || 1, trace.points.count, index, frame);
    if (!p) return null;
    const off = frameOffset(d.offset, frame);
    return [p[0] + off[0], p[1] + off[1], p[2] + off[2]];
  }

  /**
   * The world point a double-click lands on when it misses every point, as
   * classic Oviz re-centred there: the nearest visible line within `radius`
   * pixels, the middle of the dust a volume shows along the ray, or an
   * image plane in view. Null on empty space.
   */
  pickSurface(x, y, { radius = 6 } = {}) {
    const cam = this.renderer.camera;
    const { origin: o, dir: d } = cam.ray(x, y);
    let best = null;
    const consider = (t) => {
      if (t > 1e-6 && (!best || t < best)) best = t;
    };
    const along = (p) => (p[0] - o[0]) * d[0] + (p[1] - o[1]) * d[1] + (p[2] - o[2]) * d[2];
    const frame = this.timeline.frame;
    // Lines: orbit tracks, the Galactic grid.
    const styles = this.lines.params?.styles;
    const sa = [0, 0, 0], sb = [0, 0, 0], wa = [0, 0, 0], wb = [0, 0, 0];
    for (const [key, ld] of this.lineData || []) {
      const st = styles?.get(key);
      if (!st?.visible || !(st.presence > 0.001) || !(st.opacityScale > 0.002)) continue;
      const L = this.traceByKey.get(key)?.lines;
      if (!L) continue;
      const S = L.count, F = L.position.frames || 1;
      const off = frameOffset(ld.offset, frame);
      for (let i = 0; i < S; i++) {
        if (!cpuFramePosition(ld.position, F, 2 * S, 2 * i, frame, wa) || !cpuFramePosition(ld.position, F, 2 * S, 2 * i + 1, frame, wb)) continue;
        for (let k = 0; k < 3; k++) { wa[k] += off[k]; wb[k] += off[k]; }
        if (!cam.project(wa, sa) || !cam.project(wb, sb)) continue;
        const hit = screenSegment(x, y, sa, sb);
        if (hit.dist > radius) continue;
        consider(along([wa[0] + (wb[0] - wa[0]) * hit.u, wa[1] + (wb[1] - wa[1]) * hit.u, wa[2] + (wb[2] - wa[2]) * hit.u]));
      }
    }
    // Volumes: where half of what the ray sees lies.
    for (const draw of this.volumes.params?.volumeDraws || []) {
      const vol = this.volumes.volumes.get(draw.key);
      if (!vol || !(draw.fade > 0.05)) continue;
      const hit = rayBox(o, d, vol.boxMin, vol.boxMax, vol.center, draw.angle || 0);
      const voxels = hit ? this._shownVoxels(vol.spec) : null;
      const t = volumeMedianDepth(hit, voxels, vol.spec.dims, vol.boxMin, vol.boxMax, draw.low, draw.high);
      if (t != null) consider(t);
    }
    // Image planes (the face-on Milky Way) while they show.
    for (const item of this.images.items || []) {
      const st = this.images.params?.images?.get(item.key);
      if (!st?.visible) continue;
      const c = item.centers ? frameOffset(item.centers, frame) : [0, 0, 0];
      const t = rayPlaneRect(o, d, c, item.spec.widthPc, item.spec.heightPc);
      if (t != null && this._imageShowing(item)) consider(t);
    }
    return best == null ? null : [o[0] + d[0] * best, o[1] + d[1] * best, o[2] + d[2] * best];
  }

  /**
   * The voxels a volume shows now, on the CPU: its uint8 grid, or an event
   * KDE's density at the current time, rebuilt (only the GPU keeps the
   * built ones). Null until loaded.
   */
  _shownVoxels(spec) {
    if (!spec.procedural) return (spec.data && this.store.peek(spec.data.blob)) || null;
    const k = this._kde?.get(spec.key);
    if (!k) return null;
    const proc = spec.procedural;
    const time = kdeSampleTime(proc, this.timeline.time);
    return buildEventKde(k.events, proc, spec.dims, spec.bounds, time, kdeWindow(proc, this.state.volumes[spec.stateKey])).data;
  }

  /** Whether an image plane is drawn at the current zoom (it fades as the view closes in). */
  _imageShowing(item) {
    const s = item.spec;
    const hide = s.hideBelowScalePc || 0, fade = s.fadeStartScalePc || 0;
    const bar = this.images.params?.scaleBarPc;
    return !(fade > hide && bar != null && bar <= hide);
  }

  /** Move the orbit centre to `point`, keeping the zoom and the viewing angle. */
  recenter(point, { duration = 550 } = {}) {
    if (!point || this.state.view.mode === "sky") return Promise.resolve({ done: false });
    const p = clonePose(this.pose);
    p.target = point.slice();
    return this.animateTo(p, { duration });
  }

  /**
   * Per-axis median of a layer's objects present at `frame` (the classic
   * focus-group centre). Cached for the last frame asked.
   */
  layerMedian(traceKey, frame = this.timeline.frame) {
    const key = `${traceKey}|${frame.toFixed(5)}`;
    if (this._medianKey === key) return this._median?.slice() ?? null;
    const trace = this.traceByKey.get(traceKey);
    const N = trace?.points?.count || 0;
    const xs = [], ys = [], zs = [];
    for (let i = 0; i < N; i++) {
      const p = this.objectPosition(traceKey, i, frame);
      if (!p) continue;
      xs.push(p[0]); ys.push(p[1]); zs.push(p[2]);
    }
    const med = (a) => {
      a.sort((x, y) => x - y);
      const m = a.length >> 1;
      return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2;
    };
    this._medianKey = key;
    this._median = xs.length ? [med(xs), med(ys), med(zs)] : null;
    return this._median?.slice() ?? null;
  }

  flyToObject(traceKey, index, { distance } = {}) {
    const pos = this.objectPosition(traceKey, index);
    if (!pos) return Promise.resolve();
    if (this.state.view.mode === "sky") {
      const eye = this.skyEye();
      const dir = [pos[0] - eye[0], pos[1] - eye[1], pos[2] - eye[2]];
      const a = anglesFromForward(dir);
      const p = clonePose(this.pose);
      p.yaw = a.yaw; p.pitch = a.pitch;
      p.fov = Math.min(p.fov, 20);
      return this.animateTo(p, { duration: 900 });
    }
    const p = clonePose(this.pose);
    p.target = pos;
    p.distance = distance ?? clamp(Math.min(p.distance, 900), 60, 4000);
    return this.animateTo(p, { duration: 1200 });
  }

  // -------------------------------------------------------------- view mode

  /**
   * Where the Sun is at `frame`. Positions are heliocentric at the present
   * day (seen from the origin, objects sit at their catalogue l, b and
   * distance), and the Sun moves away from there with its motion relative
   * to the LSR: the origin plus the Sun trace's displacement since t = 0
   * (the trace itself may carry a constant offset, e.g. its height above
   * the plane). Distances "from the Sun" are measured from here.
   */
  sunPosition(frame = this.timeline.frame) {
    const key = this.world.sunTrace;
    if (!key) return [0, 0, 0];
    const zero = this.manifest.time?.zeroIndex ?? this.timeline.zeroFrame();
    const p = this.objectPosition(key, 0, frame), p0 = this.objectPosition(key, 0, zero);
    return p && p0 ? [p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]] : [0, 0, 0];
  }

  skyEye() {
    return [0, 0, 0];
  }

  async setViewMode(mode, { animate = true } = {}) {
    // A State transition in flight hands over first (it may change the mode).
    this.emit("camera-takeover", {});
    if (mode === this.state.view.mode) return;
    const from = this.state.view.mode;
    this.state.view.mode = mode;
    this.emit("viewmode", { mode, from });
    if (mode === "sky") {
      this.returnPose = clonePose(this.pose);
      const f = forwardFromAngles(this.pose.yaw, this.pose.pitch);
      const a = anglesFromForward(f);
      const target = makePose({ target: this.skyEye(), distance: 0, yaw: a.yaw, pitch: a.pitch * 0.35, fov: SKY_FOV });
      this.controls.mode = "sky";
      const res = animate ? await this.animateTo(target, { duration: 1400 }) : null;
      if (!res || !res.done) {
        // Interrupted (or instant): land at the Sun, keeping the current
        // gaze, so Sky registration and member stars always come up.
        if (this.state.view.mode !== "sky") return;
        const p = clonePose(this.pose);
        p.target = this.skyEye();
        p.distance = 0;
        if (!animate) { p.yaw = target.yaw; p.pitch = target.pitch; p.fov = target.fov; }
        this.setPose(p);
      }
      if (this.state.view.mode !== "sky") return;
    } else {
      // The pose Sky was entered from rode along with the anchor meanwhile;
      // a figure that opened in Sky goes home, re-anchored.
      if (!this.returnPose) this._setAnchorState(this.defaultAnchor);
      const back = this.returnPose ? clonePose(this.returnPose) : this.homeView();
      this.controls.mode = "galactic";
      if (animate) await this.animateTo(back, { duration: 1300 });
      else this.setPose(back);
      if (this.state.view.mode !== "3d") return;
    }
    this.emit("viewmode-settled", { mode });
    this.renderer.invalidate();
  }

  // -------------------------------------------------------------- styles

  setTraceStyle(key, patch) {
    const cur = this.state.traces[key] || (this.state.traces[key] = {});
    Object.assign(cur, patch);
    this._pickVersion++;
    this.renderer.invalidate();
    this.emit("style", { key, patch });
  }

  setVolumeState(key, patch) {
    const cur = this.state.volumes[key];
    if (!cur) return;
    Object.assign(cur, patch);
    this.volumes.invalidate();
    this.renderer.invalidate();
    this.emit("volume", { key, patch });
  }

  /** Motion trails of `myr` Myr, pointing back along the playback direction. */
  _trailParams(myr) {
    return trailParams(this.timeline.times, this.timeline.direction, myr);
  }

  setGlobal(patch) {
    Object.assign(this.state.global, patch);
    this._pickVersion++;
    this.renderer.invalidate();
    this.emit("global", patch);
  }

  // -------------------------------------------------------------- picking

  /** GPU pick at CSS pixel (x, y); returns {trace, index} or null. */
  pick(x, y, { radius = 6 } = {}) {
    const gl = this.gl;
    const r = this.renderer;
    const w = r.width, h = r.height;
    const sig = `${r.camera.version}|${this.timeline.frame}|${w}x${h}|${this._pickVersion}|${this.state.view.mode}|${this.memberReveal ?? 0}`;
    if (sig !== this._pickSig) {
      this.pickTarget.resize(w, h);
      this.pickTarget.bind();
      // Layers leave the depth mask off, and a masked clear skips depth:
      // without this, old depth piles up and hides objects from picking.
      gl.enable(gl.DEPTH_TEST);
      gl.depthMask(true);
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      const frame = {
        gl, camera: r.camera, width: w, height: h, cssWidth: w, cssHeight: h, dpr: 1, interacting: false,
      };
      this.points.drawPick(frame);
      this.skyPickDraw?.(frame);
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      this._pickSig = sig;
    }
    const R = Math.max(1, Math.round(radius));
    const px = Math.round(x), py = Math.round(h - y);
    const x0 = clamp(px - R, 0, w - 1), y0 = clamp(py - R, 0, h - 1);
    const x1 = clamp(px + R, 0, w - 1), y1 = clamp(py + R, 0, h - 1);
    const bw = x1 - x0 + 1, bh = y1 - y0 + 1;
    if (bw <= 0 || bh <= 0) return null;
    const buf = new Uint8Array(bw * bh * 4);
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.pickTarget.fbo);
    gl.readPixels(x0, y0, bw, bh, gl.RGBA, gl.UNSIGNED_BYTE, buf);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    let best = null, bestD = Infinity;
    for (let j = 0; j < bh; j++) {
      for (let i = 0; i < bw; i++) {
        const o = (j * bw + i) * 4;
        const tr = buf[o];
        if (!tr) continue;
        const id = (buf[o + 1] << 16) | (buf[o + 2] << 8) | buf[o + 3];
        const dx = x0 + i - px, dy = y0 + j - py;
        const d = dx * dx + dy * dy;
        if (d < bestD) { bestD = d; best = { pickTrace: tr, index: id - 1 }; }
      }
    }
    if (!best) return null;
    if (best.pickTrace >= 200) return this.skyPickResolve?.(best) || null;
    const key = this.pickIds[best.pickTrace];
    return key ? { trace: key, index: best.index } : null;
  }

  invalidatePick() {
    this._pickVersion++;
  }

  setHover(hit) {
    const same = (hit && this.hover && hit.trace === this.hover.trace && hit.index === this.hover.index) || (!hit && !this.hover);
    if (same) return;
    this.hover = hit;
    this.renderer.invalidate();
    this.emit("hover", hit);
  }

  // -------------------------------------------------------------- metadata

  async objectMeta(traceKey) {
    if (this.meta.has(traceKey)) return this.meta.get(traceKey);
    const trace = this.traceByKey.get(traceKey);
    const blob = trace?.points?.meta?.blob;
    const cols = blob ? await this.store.get(blob) : { name: [] };
    this.meta.set(traceKey, cols);
    return cols;
  }

  objectMetaSync(traceKey) {
    return this.meta.get(traceKey) || null;
  }

  /** Structured info for the inspector/tooltip. */
  describeObject(traceKey, index) {
    const trace = this.traceByKey.get(traceKey);
    const d = this.data.get(traceKey);
    const meta = this.objectMetaSync(traceKey) || {};
    if (!trace) return null;
    const pos = this.objectPosition(traceKey, index);
    const name = (meta.name?.[index] || "").replace(/_/g, " ") || `${trace.name} #${index + 1}`;
    const ageNow = d?.ageNow ? d.ageNow[index] : NaN;
    const t = this.timeline.time;
    const nStars = d?.nStars ? d.nStars[index] : NaN;
    const color = d?.rgb ? `rgb(${d.rgb[index * 3]}, ${d.rgb[index * 3 + 1]}, ${d.rgb[index * 3 + 2]})` : trace.color;
    return {
      trace: traceKey,
      index,
      traceName: trace.name,
      name,
      color,
      aliases: (meta.aliases?.[index] || "").split(",").map((s) => s.trim().replace(/_/g, " ")).filter((s) => s && s !== name),
      ageNow: ageNow === ageNow ? ageNow : null,
      ageAt: ageNow === ageNow ? ageNow + t : null,
      nStars: nStars >= 0 ? nStars : null,
      l: meta.l?.[index] ?? null,
      b: meta.b?.[index] ?? null,
      dist: meta.dist?.[index] ?? null,
      ra: meta.ra?.[index] ?? null,
      dec: meta.dec?.[index] ?? null,
      hover: meta.hover?.[index] || "",
      position: pos,
      distanceFromSun: pos ? v3dist(pos, this.sunPosition()) : null,
    };
  }

  // -------------------------------------------------------------- scale

  /** Length in pc of a 120-px ruler at the orbit target (3D) or degrees (Sky). */
  scaleBarPc(px = 120) {
    const cam = this.renderer.camera;
    const p = cam.pose;
    if (p.distance <= 1e-6) return null;
    return cam.pixelScale(p.distance) * px;
  }

  scaleBar(px = 120) {
    const cam = this.renderer.camera;
    const p = cam.pose;
    if (this.state.view.mode === "sky" || p.distance <= 1e-6) {
      const deg = (p.fov * px) / Math.max(cam.height, 1);
      const nice = niceFloor(deg);
      return { value: nice, unit: "deg", px: (nice / deg) * px };
    }
    const pc = cam.pixelScale(p.distance) * px;
    const nice = niceFloor(pc);
    return { value: nice, unit: "pc", px: (nice / pc) * px };
  }

  dispose() {
    this.renderer.dispose();
  }
}

function isMobile() {
  const ua = navigator.userAgent || "";
  return /iPhone|iPad|iPod|Android/i.test(ua) || (navigator.maxTouchPoints > 1 && Math.min(screen.width, screen.height) < 820);
}

