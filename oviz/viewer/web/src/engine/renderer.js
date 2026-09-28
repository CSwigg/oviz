// Frame scheduler and pass orchestration.
//
// Rendering is demand-driven: nothing is drawn unless something marked the
// frame dirty (camera motion, time change, style edit, async asset arrival)
// or an animation asked for continuous frames. Idle figures cost ~0 CPU/GPU.

import { createContext, enableExtensions } from "./gl.js";
import { Camera } from "./camera.js";

export class Renderer {
  constructor(canvas, { maxDpr = 2 } = {}) {
    this.canvas = canvas;
    const { gl, caps } = createContext(canvas);
    this.gl = gl;
    this.caps = caps;
    this.camera = new Camera();
    this.maxDpr = maxDpr;
    this.dpr = 1;
    this.width = 1;
    this.height = 1;
    this.layers = [];
    this.overlays = [];
    this.beforeRender = [];
    this.afterRender = [];
    this.dirty = true;
    this.continuous = new Set();
    this.frameCount = 0;
    this.lastFrameMs = 0;
    this.frameTimes = new Float32Array(90);
    this.frameTimeIndex = 0;
    this.sceneRadius = 1e4;
    this.clearColor = [0, 0, 0, 0];
    this.interacting = false;
    this._lastInteraction = 0;
    this._raf = 0;
    this._running = false;
    this._resizeObserver = new ResizeObserver(() => this.resize());
    this._resizeObserver.observe(canvas.parentElement || canvas);
    this.contextLost = false;
    this.onContextRestored = null;
    canvas.addEventListener("webglcontextlost", (e) => {
      // preventDefault asks the browser to restore the context later.
      e.preventDefault();
      this.contextLost = true;
    });
    canvas.addEventListener("webglcontextrestored", () => {
      this.contextLost = false;
      this.caps = enableExtensions(this.gl);
      // Every GL object died with the old context; the owner rebuilds them.
      this.layers = [];
      this.onContextRestored?.();
      this.invalidate();
    });
    this.resize();
  }

  add(layer) {
    this.layers.push(layer);
    this.layers.sort((a, b) => (a.order || 0) - (b.order || 0));
    layer.attach?.(this);
    this.invalidate();
    return layer;
  }

  remove(layer) {
    const i = this.layers.indexOf(layer);
    if (i >= 0) this.layers.splice(i, 1);
    layer.dispose?.();
    this.invalidate();
  }

  invalidate() {
    this.dirty = true;
    this._schedule();
  }

  /** Mark user interaction; volumes use this to drop to interactive LOD. */
  markInteraction() {
    this.interacting = true;
    this._lastInteraction = performance.now();
    this.invalidate();
  }

  /** Keep rendering every frame while `key` is held (animations, playback). */
  hold(key) {
    this.continuous.add(key);
    this._schedule();
  }

  release(key) {
    this.continuous.delete(key);
    this.invalidate();
  }

  resize() {
    const el = this.canvas.parentElement || this.canvas;
    // Layout size, not the painted box: a transform on the stage (the
    // entrance zoom) must not change the drawing buffer.
    const dpr = Math.min(window.devicePixelRatio || 1, this.maxDpr);
    const w = Math.max(1, el.clientWidth);
    const h = Math.max(1, el.clientHeight);
    if (w === this.width && h === this.height && dpr === this.dpr) return;
    this.width = w;
    this.height = h;
    this.dpr = dpr;
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    this.canvas.style.width = `${w}px`;
    this.canvas.style.height = `${h}px`;
    this.camera.setViewport(w, h);
    for (const layer of this.layers) layer.resize?.(this);
    this.invalidate();
  }

  start() {
    if (this._running) return;
    this._running = true;
    this._schedule();
  }

  stop() {
    this._running = false;
    cancelAnimationFrame(this._raf);
    this._raf = 0;
  }

  _schedule() {
    if (!this._running || this._raf) return;
    // While an external loop drives us (Aladin in Sky view) it calls
    // externalTick every frame; fall back to our own rAF if it stalls.
    if (this.externalDriver?.() && performance.now() - (this._lastExternal || 0) < 120) return;
    this._raf = requestAnimationFrame((t) => this._tick(t));
  }

  /**
   * Run one frame from an external animation loop. Aladin Lite calls this
   * right before drawing its survey, so the camera pose pushed to Aladin,
   * the WebGL data and the survey image are all produced in the same frame.
   */
  externalTick(now) {
    if (!this._running || !this.externalDriver?.()) return;
    this._lastExternal = performance.now();
    if (this._raf) {
      cancelAnimationFrame(this._raf);
      this._raf = 0;
    }
    this._tick(now);
  }

  _tick(now) {
    this._raf = 0;
    if (!this._running) return;
    if (this.interacting && now - this._lastInteraction > 160) {
      this.interacting = false;
      this.dirty = true; // settle: refine to full quality
    }
    for (const fn of this.beforeRender) fn(now);
    const needs = this.dirty || this.continuous.size > 0;
    if (needs && !this.contextLost) {
      const t0 = performance.now();
      this.dirty = false;
      this.render(now);
      const dt = performance.now() - t0;
      this.lastFrameMs = dt;
      this.frameTimes[this.frameTimeIndex++ % this.frameTimes.length] = dt;
      this.frameCount++;
      for (const fn of this.afterRender) fn(now);
    }
    if (this.continuous.size > 0 || this.dirty || this.interacting) this._schedule();
  }

  render(now = performance.now()) {
    const gl = this.gl;
    const cam = this.camera;
    cam.update(this.sceneRadius);
    const frame = {
      now,
      gl,
      camera: cam,
      renderer: this,
      width: this.canvas.width,
      height: this.canvas.height,
      cssWidth: this.width,
      cssHeight: this.height,
      dpr: this.dpr,
      interacting: this.interacting,
    };
    for (const layer of this.layers) if (layer.visible !== false) layer.prepare?.(frame);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, frame.width, frame.height);
    const c = this.clearColor;
    gl.clearColor(c[0], c[1], c[2], c[3]);
    gl.depthMask(true);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    for (const layer of this.layers) {
      if (layer.visible === false) continue;
      layer.draw(frame);
    }
    for (const overlay of this.overlays) overlay.update?.(frame);
  }

  /**
   * Render one frame into an image blob at `scale` × the current resolution.
   * An async background (e.g. the Sky survey) is prepared first; the WebGL
   * frame is then rendered and composited in the same task, before the
   * browser presents (and clears) the drawing buffer.
   */
  async capture({ scale = 1, background = null, type = "image/png", quality } = {}) {
    const w = this.width, h = this.height;
    const prevDpr = this.dpr;
    const dpr = Math.min(window.devicePixelRatio || 1, this.maxDpr) * scale;
    const W = Math.round(w * dpr), H = Math.round(h * dpr);
    let off = null, ctx = null;
    if (background) {
      off = document.createElement("canvas");
      off.width = W;
      off.height = H;
      ctx = off.getContext("2d");
      if (typeof background === "string") {
        ctx.fillStyle = background;
        ctx.fillRect(0, 0, W, H);
      } else if (typeof background === "function") {
        await background(ctx, W, H);
      }
    }
    this.canvas.width = W;
    this.canvas.height = H;
    this.dpr = dpr;
    for (const layer of this.layers) layer.resize?.(this);
    const prevInteracting = this.interacting;
    this.interacting = false;
    this.render();
    if (ctx) ctx.drawImage(this.canvas, 0, 0);
    const source = off || this.canvas;
    const blobPromise = new Promise((resolve) => source.toBlob(resolve, type, quality));
    this.dpr = prevDpr;
    this.interacting = prevInteracting;
    this.canvas.width = Math.round(w * prevDpr);
    this.canvas.height = Math.round(h * prevDpr);
    for (const layer of this.layers) layer.resize?.(this);
    this.invalidate();
    return blobPromise;
  }

  fps() {
    let sum = 0, n = 0;
    for (const t of this.frameTimes) if (t > 0) { sum += t; n++; }
    return n ? sum / n : 0;
  }

  dispose() {
    this.stop();
    this._resizeObserver.disconnect();
    for (const layer of this.layers) layer.dispose?.();
    this.layers = [];
  }
}
