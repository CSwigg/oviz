// Pointer, wheel, touch and keyboard navigation with inertia.
//
// Galactic mode: drag orbits, right/shift/two-finger drag pans, wheel and
// pinch zoom toward the cursor. Sky mode: drag grabs the sky (the point under
// the cursor stays under the cursor), wheel/pinch change the field of view.

import { DEG, clamp } from "../core/math.js";
import { forwardFromAngles } from "./camera.js";

const PITCH_LIMIT = 89.5 * DEG;

export class Controls {
  constructor(element, renderer, { onChange, onInteractStart, onInteractEnd } = {}) {
    this.el = element;
    this.renderer = renderer;
    this.camera = renderer.camera;
    this.onChange = onChange || (() => {});
    this.onInteractStart = onInteractStart || (() => {});
    this.onInteractEnd = onInteractEnd || (() => {});
    this.mode = "galactic"; // "galactic" | "sky"
    this.enabled = true;
    this.minDistance = 1;
    this.maxDistance = 2e5;
    this.minFov = 0.2;
    this.maxFov = 130;
    this.rotateSpeed = 1;
    this.damping = 0.14;
    this.velocity = { yaw: 0, pitch: 0, panX: 0, panY: 0, zoom: 0 };
    this.pointers = new Map();
    this.gesture = null;
    this.autoOrbit = 0; // rad/s
    this._lastMove = 0;
    this._bind();
  }

  _bind() {
    const el = this.el;
    el.style.touchAction = "none";
    el.addEventListener("pointerdown", (e) => this._down(e));
    el.addEventListener("pointermove", (e) => this._move(e));
    el.addEventListener("pointerup", (e) => this._up(e));
    el.addEventListener("pointercancel", (e) => this._up(e));
    el.addEventListener("lostpointercapture", (e) => this._up(e));
    el.addEventListener("wheel", (e) => this._wheel(e), { passive: false });
    el.addEventListener("contextmenu", (e) => e.preventDefault());
    this.renderer.beforeRender.push((now) => this._integrate(now));
  }

  get pose() {
    return this.camera.pose;
  }

  stop() {
    const v = this.velocity;
    v.yaw = v.pitch = v.panX = v.panY = v.zoom = 0;
  }

  _down(e) {
    if (!this.enabled) return;
    if (e.button === 1) return;
    this.el.setPointerCapture?.(e.pointerId);
    this.pointers.set(e.pointerId, { x: e.clientX, y: e.clientY, button: e.button, shift: e.shiftKey });
    this.stop();
    this._startGesture(e);
    this.onInteractStart(e);
  }

  _startGesture(e) {
    const pts = [...this.pointers.values()];
    if (pts.length >= 2) {
      const [a, b] = pts;
      this.gesture = {
        kind: "pinch",
        dist: Math.hypot(a.x - b.x, a.y - b.y),
        cx: (a.x + b.x) / 2,
        cy: (a.y + b.y) / 2,
      };
    } else if (pts.length === 1) {
      const p = pts[0];
      const pan = p.button === 2 || p.shift || (e && (e.ctrlKey || e.metaKey));
      this.gesture = { kind: pan && this.mode !== "sky" ? "pan" : "rotate", x: p.x, y: p.y };
    } else {
      this.gesture = null;
    }
  }

  _move(e) {
    const p = this.pointers.get(e.pointerId);
    if (!p || !this.enabled) return;
    const dt = Math.max(1, performance.now() - (this._lastMove || 0));
    this._lastMove = performance.now();
    const dx = e.clientX - p.x;
    const dy = e.clientY - p.y;
    p.x = e.clientX;
    p.y = e.clientY;
    const g = this.gesture;
    if (!g) return;
    this.renderer.markInteraction();
    if (g.kind === "pinch") {
      const pts = [...this.pointers.values()];
      if (pts.length < 2) return;
      const [a, b] = pts;
      const dist = Math.hypot(a.x - b.x, a.y - b.y);
      const cx = (a.x + b.x) / 2, cy = (a.y + b.y) / 2;
      const factor = g.dist > 0 ? g.dist / dist : 1;
      this.zoomAt(factor, cx, cy);
      if (this.mode !== "sky") this.pan(cx - g.cx, cy - g.cy);
      else this.lookDrag(cx - g.cx, cy - g.cy);
      g.dist = dist; g.cx = cx; g.cy = cy;
    } else if (g.kind === "pan") {
      this.pan(dx, dy);
      const k = 16 / dt;
      this.velocity.panX = dx * k * 0.5;
      this.velocity.panY = dy * k * 0.5;
    } else if (this.mode === "sky") {
      this.lookDrag(dx, dy);
      const k = 16 / dt;
      const perPx = this.pose.fov * DEG / this.camera.height;
      this.velocity.yaw = dx * perPx * k * 0.6;
      this.velocity.pitch = dy * perPx * k * 0.6;
    } else {
      const s = (2 * Math.PI) / Math.max(this.camera.height, 400) * this.rotateSpeed * 0.9;
      this.rotate(-dx * s, -dy * s);
      const k = 16 / dt;
      this.velocity.yaw = -dx * s * k * 0.55;
      this.velocity.pitch = -dy * s * k * 0.55;
    }
  }

  _up(e) {
    if (!this.pointers.has(e.pointerId)) return;
    this.pointers.delete(e.pointerId);
    this.el.releasePointerCapture?.(e.pointerId);
    // A pause before release means the user stopped; do not fling.
    if (performance.now() - this._lastMove > 60) this.stop();
    this._startGesture(null);
    if (this.pointers.size === 0) this.onInteractEnd(e);
  }

  _wheel(e) {
    if (!this.enabled) return;
    e.preventDefault();
    let dy = e.deltaY;
    if (e.deltaMode === 1) dy *= 16;
    else if (e.deltaMode === 2) dy *= 400;
    // Trackpad pinch arrives as ctrl+wheel with small deltas.
    const intensity = e.ctrlKey ? 0.012 : 0.0016;
    const factor = Math.exp(clamp(dy, -240, 240) * intensity);
    const rect = this.el.getBoundingClientRect();
    this.zoomAt(factor, e.clientX - rect.left, e.clientY - rect.top);
    this.renderer.markInteraction();
  }

  rotate(dYaw, dPitch) {
    const p = this.pose;
    p.yaw += dYaw;
    p.pitch = clamp(p.pitch + dPitch, -PITCH_LIMIT, PITCH_LIMIT);
    this._changed();
  }

  /** Sky-mode drag: rotate so the sky follows the cursor. */
  lookDrag(dx, dy) {
    const p = this.pose;
    const perPx = (p.fov * DEG) / this.camera.height;
    p.yaw += dx * perPx;
    p.pitch = clamp(p.pitch + dy * perPx, -PITCH_LIMIT, PITCH_LIMIT);
    this._changed();
  }

  pan(dx, dy) {
    const p = this.pose;
    const cam = this.camera;
    const scale = cam.pixelScale(Math.max(p.distance, 1));
    const r = cam.right, u = cam.camUp;
    for (let i = 0; i < 3; i++) p.target[i] += (-dx * r[i] + dy * u[i]) * scale;
    this._changed();
  }

  /** Zoom by `factor` (<1 zooms in) keeping the point under (x, y) fixed. */
  zoomAt(factor, x, y) {
    const p = this.pose;
    if (this.mode === "sky" || p.distance <= 1e-6) {
      const before = p.fov;
      p.fov = clamp(p.fov * factor, this.minFov, this.maxFov);
      // Keep the cursor direction steady while zooming.
      const cam = this.camera;
      if (x != null && cam.width > 0) {
        const k = 1 - p.fov / before;
        const perPx = (before * DEG) / cam.height;
        p.yaw += (x - cam.width / 2) * perPx * k * -1;
        p.pitch = clamp(p.pitch - (y - cam.height / 2) * perPx * k, -PITCH_LIMIT, PITCH_LIMIT);
      }
      this._changed();
      return;
    }
    const next = clamp(p.distance * factor, this.minDistance, this.maxDistance);
    const actual = next / p.distance;
    if (x != null) {
      // Move the target toward the cursor ray so zoom follows the pointer.
      const cam = this.camera;
      const ray = cam.ray(x, y);
      const f = forwardFromAngles(p.yaw, p.pitch);
      const denom = ray.dir[0] * f[0] + ray.dir[1] * f[1] + ray.dir[2] * f[2];
      if (denom > 1e-4) {
        const t = p.distance / denom;
        const hit = [ray.origin[0] + ray.dir[0] * t, ray.origin[1] + ray.dir[1] * t, ray.origin[2] + ray.dir[2] * t];
        for (let i = 0; i < 3; i++) p.target[i] = hit[i] + (p.target[i] - hit[i]) * actual;
      }
    }
    p.distance = next;
    this._changed();
  }

  _integrate(now) {
    const dtRaw = this._lastIntegrate ? now - this._lastIntegrate : 16;
    this._lastIntegrate = now;
    const dt = clamp(dtRaw, 0, 64);
    const v = this.velocity;
    let moved = false;
    if (this.pointers.size === 0) {
      const decay = Math.pow(1 - this.damping, dt / 16);
      if (Math.abs(v.yaw) + Math.abs(v.pitch) > 1e-5) {
        if (this.mode === "sky") {
          this.pose.yaw += v.yaw;
          this.pose.pitch = clamp(this.pose.pitch + v.pitch, -PITCH_LIMIT, PITCH_LIMIT);
        } else {
          this.pose.yaw += v.yaw;
          this.pose.pitch = clamp(this.pose.pitch + v.pitch, -PITCH_LIMIT, PITCH_LIMIT);
        }
        v.yaw *= decay; v.pitch *= decay;
        moved = true;
      } else {
        v.yaw = v.pitch = 0;
      }
      if (Math.abs(v.panX) + Math.abs(v.panY) > 0.05) {
        this.pan(v.panX, v.panY);
        v.panX *= decay; v.panY *= decay;
        moved = true;
      } else {
        v.panX = v.panY = 0;
      }
      if (this.autoOrbit && this.mode !== "sky") {
        this.pose.yaw += this.autoOrbit * dt / 1000;
        moved = true;
      }
    }
    if (moved) {
      this.renderer.markInteraction();
      this._changed();
    }
    if (this.autoOrbit && this.mode !== "sky") this.renderer.hold("auto-orbit");
    else this.renderer.continuous.delete("auto-orbit");
  }

  _changed() {
    this.renderer.invalidate();
    this.onChange();
  }
}
