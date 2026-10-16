// Frame scheduler and pass orchestration.
//
// Rendering is demand-driven: nothing is drawn unless something marked the
// frame dirty (camera motion, time change, style edit, async asset arrival)
// or an animation asked for continuous frames. Idle figures cost ~0 CPU/GPU.

import { createContext, enableExtensions, createProgram, RenderTarget, FULLSCREEN_VS } from "./gl.js";
import { Camera } from "./camera.js";

// Copies the cached static scene to the screen (premultiplied, one texel per pixel).
const COPY_FS = `
in vec2 vUv;
uniform sampler2D uTex;
out vec4 outColor;
void main() { outColor = texelFetch(uTex, ivec2(gl_FragCoord.xy), 0); }`;

// Largest capture we ask for (Chrome's drawing-buffer area limit).
const MAX_CAPTURE_PIXELS = 7680 * 4320;

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
    // Holds that only move "animated" layers (flow pulses): while nothing else
    // changes, those frames redraw the animated layers over a cached copy of
    // everything else instead of re-drawing the whole scene.
    this.animatedHolds = new Set(["flow"]);
    this._base = null;
    this._baseValid = false;
    this._copy = null;
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
    this._watchDpr();
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
      this._base = null;
      this._copy = null;
      this._baseValid = false;
      this.onContextRestored?.();
      this.invalidate();
    });
    this.resize();
  }

  /** Re-measure when the window moves to a display with another pixel ratio. */
  _watchDpr() {
    const mq = window.matchMedia?.(`(resolution: ${window.devicePixelRatio || 1}dppx)`);
    mq?.addEventListener?.("change", () => { this.resize(); this._watchDpr(); }, { once: true });
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
    // Resizing clears the canvas: redraw at once so animated layout changes
    // (a sidebar sliding open) never show a blank frame.
    if (this._running && !this.contextLost) this.render();
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
      // Nothing but animated layers moved since the last frame.
      let reuseBase = !this.dirty && !this.interacting && this.continuous.size > 0;
      if (reuseBase) for (const k of this.continuous) if (!this.animatedHolds.has(k)) { reuseBase = false; break; }
      this.dirty = false;
      this.render(now, { reuseBase });
      const dt = performance.now() - t0;
      this.lastFrameMs = dt;
      this.frameTimes[this.frameTimeIndex++ % this.frameTimes.length] = dt;
      this.frameCount++;
      for (const fn of this.afterRender) fn(now);
    }
    // While the context is lost nothing can draw: wait for the restore
    // handler (it invalidates) instead of spinning every frame.
    if (this.contextLost) return;
    if (this.continuous.size > 0 || this.dirty || this.interacting) this._schedule();
  }

  render(now = performance.now(), { reuseBase = false } = {}) {
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
      capture: !!this._capturing,
    };
    for (const layer of this.layers) if (layer.visible !== false) layer.prepare?.(frame);
    const animated = frame.capture ? [] : this.layers.filter((l) => l.animated && l.visible !== false && l.active?.() !== false);
    if (animated.length) {
      this._renderCached(frame, animated, reuseBase);
    } else {
      this._baseValid = false;
      this._clear(null, frame.width, frame.height);
      for (const layer of this.layers) {
        if (layer.visible === false) continue;
        layer.draw(frame);
      }
    }
    for (const overlay of this.overlays) overlay.update?.(frame);
  }

  _clear(fbo, width, height, color = this.clearColor) {
    const gl = this.gl;
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    gl.viewport(0, 0, width, height);
    gl.clearColor(color[0], color[1], color[2], color[3]);
    gl.depthMask(true);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  }

  /**
   * Draw the static layers into a cached target (only when something other
   * than the animated layers changed), copy it to the screen, then draw the
   * animated layers on top. No layer writes depth, so colour is all the
   * cache needs.
   */
  _renderCached(frame, animated, reuseBase) {
    const gl = this.gl;
    const base = (this._base ||= new RenderTarget(gl, { filter: gl.NEAREST }));
    const resized = base.resize(frame.width, frame.height);
    if (!reuseBase || resized || !this._baseValid || this.camera.version !== this._baseCamera) {
      this._clear(base.fbo, frame.width, frame.height);
      for (const layer of this.layers) {
        if (layer.visible === false || animated.includes(layer)) continue;
        layer.draw(frame);
      }
      this._baseValid = true;
      this._baseCamera = this.camera.version;
    }
    this._clear(null, frame.width, frame.height, [0, 0, 0, 0]);
    const copy = (this._copy ||= createProgram(gl, FULLSCREEN_VS, COPY_FS, { label: "base-copy" }));
    copy.use().tex("uTex", base.texture);
    gl.disable(gl.BLEND);
    gl.disable(gl.DEPTH_TEST);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    for (const layer of animated) layer.draw(frame);
  }

  /**
   * Render one frame into an image blob at `scale` × the current resolution.
   * An async background (e.g. the Sky survey) is prepared first; the WebGL
   * frame is then rendered and composited in the same task, before the
   * browser presents (and clears) the drawing buffer.
   */
  async capture({ scale = 1, background = null, type = "image/png", quality } = {}) {
    const gl = this.gl;
    const w = this.width, h = this.height;
    const prevDpr = this.dpr;
    let dpr = Math.min(window.devicePixelRatio || 1, this.maxDpr) * scale;
    // Browsers cap the drawing buffer (Chrome: 7680×4320 worth of pixels);
    // a larger canvas is silently cropped. Stay inside the caps.
    const maxDim = Math.min(gl.getParameter(gl.MAX_VIEWPORT_DIMS)[0], gl.getParameter(gl.MAX_VIEWPORT_DIMS)[1], gl.getParameter(gl.MAX_RENDERBUFFER_SIZE), 16384);
    const fit = Math.min(1, maxDim / (w * dpr), maxDim / (h * dpr), Math.sqrt(MAX_CAPTURE_PIXELS / (w * dpr * h * dpr)));
    dpr *= fit;
    let W = Math.floor(w * dpr), H = Math.floor(h * dpr);
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
    // Some browsers grant less than asked even inside the caps: render at
    // what was granted and scale it onto the background.
    if (gl.drawingBufferWidth < W || gl.drawingBufferHeight < H) {
      const k = Math.min(gl.drawingBufferWidth / W, gl.drawingBufferHeight / H);
      dpr *= k;
      W = Math.floor(w * dpr);
      H = Math.floor(h * dpr);
      this.canvas.width = W;
      this.canvas.height = H;
    }
    this.dpr = dpr;
    this.lastCapture = { width: W, height: H };
    for (const layer of this.layers) layer.resize?.(this);
    const prevInteracting = this.interacting;
    this.interacting = false;
    this._capturing = true;
    try { this.render(); } finally { this._capturing = false; }
    if (ctx) ctx.drawImage(this.canvas, 0, 0, off.width, off.height);
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
