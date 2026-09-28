// The Viewer: owns the renderer, layers, time and view state, and exposes
// the imperative API that the UI, States and the public JS API drive.

import { Emitter } from "../core/emitter.js";
import { clamp, DEG, v3dist, niceFloor } from "../core/math.js";
import { Renderer } from "../engine/renderer.js";
import { Controls } from "../engine/controls.js";
import { RenderTarget } from "../engine/gl.js";
import { makePose, clonePose, poseFromEyeTarget, poseTween, forwardFromAngles, anglesFromForward } from "../engine/camera.js";
import { cpuFramePosition, frameOffset } from "../engine/frames.js";
import { PointsLayer, resetStarTextures } from "../layers/points.js";
import { LinesLayer } from "../layers/lines.js";
import { ImagesLayer } from "../layers/images.js";
import { VolumesLayer } from "../layers/volumes.js";
import { LabelsOverlay } from "../layers/labels.js";
import { Timeline } from "./timeline.js";
import { initialViewerState, clampVolumeState } from "./state.js";

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
      onInteractStart: () => this._cancelTween(),
    });
    this.controls.maxDistance = Math.max(this.world.maxSpan || 1e4, 1e3) * 8;
    this.labels = new LabelsOverlay(container);
    this.renderer.overlays.push(this.labels);
    this.renderer.beforeRender.push((now) => this._beforeRender(now));
    this.renderer.onContextRestored = () => this._restoreGPU();
    this.timeline.on(() => {
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
        this.lines.addTrace(trace, ld);
      }
      if (trace.labels && labels) {
        const L = trace.labels;
        const ld = { position: await store.get(L.position.blob) };
        if (L.position.offset) ld.offset = await store.get(L.position.offset);
        this.labels.addTrace(trace, ld);
      }
    }
    this.renderer.invalidate();
  }

  /**
   * Rebuild every GPU resource after a lost WebGL context is restored.
   * Decoded arrays stay cached in the bundle store, so this re-uploads
   * without re-decoding; plugins rebuild their own layers on "gpu-restored".
   */
  async _restoreGPU() {
    resetStarTextures();
    this.data.clear();
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
    for (const img of this.manifest.images || []) {
      try {
        const blob = await this.store.get(img.image.blob);
        const centers = img.center?.blob ? await this.store.get(img.center.blob) : null;
        const opac = img.opacity?.blob ? await this.store.get(img.opacity.blob) : null;
        await this.images.add(img, blob, centers, opac);
        this.renderer.invalidate();
      } catch (err) {
        console.warn("Oviz: image plane failed to load", img.key, err);
      }
    }
  }

  async attachVolumes(onEach) {
    // Smallest first so something appears quickly.
    const specs = [...this.volumeSpecs].sort((a, b) => a.dims[0] * a.dims[1] * a.dims[2] - b.dims[0] * b.dims[1] * b.dims[2]);
    for (const spec of specs) {
      try {
        const data = await this.store.get(spec.data.blob);
        let occ = null;
        if (spec.occupancy?.blob) {
          const o = await this.store.get(spec.occupancy.blob);
          occ = { data: o, gx: spec.occupancy.dims[0], gy: spec.occupancy.dims[1], gz: spec.occupancy.dims[2] };
        }
        this.volumes.addVolume(spec, data, occ);
        this.renderer.invalidate();
        onEach?.(spec);
        this.emit("volume-ready", spec);
        await new Promise((r) => setTimeout(r, 0));
      } catch (err) {
        console.warn("Oviz: volume failed to load", spec.key, err);
      }
    }
  }

  _applyInitialCamera() {
    const cam = this.manifest.camera || {};
    const pos = cam.position || [0, -800, 1500];
    const tgt = cam.target || this.world.center || [0, 0, 0];
    const pose = poseFromEyeTarget(pos, tgt, cam.fov || 60);
    this.homePose = clonePose(pose);
    this.renderer.camera.pose = pose;
    this.state.view.pose = pose;
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
        const cb = this.tween.onDone;
        this.tween = null;
        this.renderer.continuous.delete("tween");
        cb?.();
      }
      this._onCameraChange(true);
    }
    for (const fn of this._animators || []) fn(now);
    this._resolveParams();
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
      const isReference = !trace.showInLegend;
      let visible = ts.visible !== false;
      if (isReference) visible = visible && g.grid !== false && !skyMode;
      const defaultOpacity = trace.opacity > 1e-6 ? trace.opacity : 1;
      let starsExp = 0;
      if (trace.hasStars) starsExp = g.sizeByStars ? (trace.sizeByStarsDefault ? 0 : 1) : (trace.sizeByStarsDefault ? -1 : 0);
      const hoverIdx = this.hover && this.hover.trace === trace.key ? this.hover.index : -1;
      let opacityScale = (ts.opacity ?? trace.opacity) / defaultOpacity;
      if (isReference) opacityScale *= g.gridOpacity ?? 1;
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
    };
    this.points.params = pointParams;
    this.lines.params = { styles, frame, lineOpacity: 1 };
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
      s += `${trace.key}:${ts.visible !== false}:${ts.opacity}:${ts.sizeScale};`;
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
      for (const [spec, w] of weights) {
        if (w <= 0.002) continue;
        const vsClamped = clampVolumeState(spec, vs);
        const [dmin, dmax] = spec.dataRange;
        const span = dmax - dmin;
        let low = span > 0 ? clamp((vsClamped.vmin - dmin) / span, 0, 1) : 0;
        let high = span > 0 ? clamp((vsClamped.vmax - dmin) / span, 0, 1) : 1;
        if (!(high > low)) high = Math.min(1, low + 1e-6);
        const allTimes = vs.showAllTimes && spec.supportsShowAllTimes && spec.timeMyr == null;
        const angle = allTimes && spec.coRotate ? (spec.coRotationRate || 0) * (time - (spec.referenceTimeMyr || 0)) : 0;
        draws.push({ key: spec.key, state: vsClamped, low, high, fade: w, angle, canSkip: true });
        sig += `${spec.key}:${w.toFixed(3)}:${low}:${high}:${vsClamped.opacity}:${vsClamped.steps}:${vsClamped.alphaCoef}:${vsClamped.stretch}:${vsClamped.colormap}:${angle.toFixed(5)}:${vsClamped.lightingMode};`;
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

  _cancelTween() {
    if (this.tween) {
      this.tween = null;
      this.renderer.continuous.delete("tween");
    }
  }

  animateTo(pose, { duration = 1100, onDone } = {}) {
    this.controls.stop();
    const d = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches ? 0 : duration;
    this.tween = poseTween(this.renderer.camera.pose, pose, d);
    this.tween.onDone = onDone;
    this.renderer.hold("tween");
    return new Promise((resolve) => {
      const prev = this.tween.onDone;
      this.tween.onDone = () => { prev?.(); resolve(); };
    });
  }

  setPose(pose) {
    this._cancelTween();
    Object.assign(this.renderer.camera.pose, clonePose(pose));
    this._onCameraChange();
    this.renderer.invalidate();
  }

  resetView() {
    if (this.state.view.mode === "sky") {
      const p = clonePose(this.pose);
      p.fov = SKY_FOV;
      return this.animateTo(p, { duration: 700 });
    }
    return this.animateTo(this.homePose, { duration: 1000 });
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

  skyEye() {
    return [0, 0, 0];
  }

  async setViewMode(mode, { animate = true } = {}) {
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
      if (animate) await this.animateTo(target, { duration: 1400 });
      else this.setPose(target);
    } else {
      const back = this.returnPose || this.homePose;
      this.controls.mode = "galactic";
      if (animate) await this.animateTo(back, { duration: 1300 });
      else this.setPose(back);
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
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      const frame = {
        gl, camera: r.camera, width: w, height: h, cssWidth: w, cssHeight: h, dpr: 1, interacting: false,
      };
      gl.enable(gl.DEPTH_TEST);
      gl.depthMask(true);
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
      distanceFromSun: pos ? v3dist(pos, this.skyEye()) : null,
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

  fovDeg() {
    return this.pose.fov;
  }

  dispose() {
    this.renderer.dispose();
  }
}

export function isMobile() {
  const ua = navigator.userAgent || "";
  return /iPhone|iPad|iPod|Android/i.test(ua) || (navigator.maxTouchPoints > 1 && Math.min(screen.width, screen.height) < 820);
}

