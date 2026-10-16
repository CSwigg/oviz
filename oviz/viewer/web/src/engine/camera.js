// Camera model shared by the Galactic 3D view and the Sun-centred Sky view.
//
// The pose is (target, distance, yaw, pitch, fov). Galactic +z is always up,
// so the horizon never rolls. The Sky view is the same pose with distance 0:
// the eye sits on the target and looks along (yaw, pitch). Because a pinhole
// camera *is* a gnomonic (TAN) projection, this makes the Aladin TAN view and
// the WebGL scene register exactly when their centre, width-FOV and roll agree.

import {
  DEG, clamp, lerp, angleDelta, m4, m4mul, m4perspective, m4lookAt, m4invert,
  m4project, easeInOutCubic,
} from "../core/math.js";

const UP = [0, 0, 1];

/** Orbits stop just short of straight up/down, where yaw is undefined. */
export const PITCH_LIMIT = 89.5 * DEG;

/**
 * Keep a pose inside the orbit limits (in place). Exactly straight up or
 * down has no yaw; the look-at then put +x to the right, so keep that view.
 */
export function sanitizePose(p) {
  if (!Number.isFinite(p.pitch)) p.pitch = 0;
  if (Math.abs(p.pitch) >= Math.PI / 2 - 1e-9) p.yaw = Math.PI / 2;
  p.pitch = clamp(p.pitch, -PITCH_LIMIT, PITCH_LIMIT);
  return p;
}

export function makePose(p = {}) {
  return {
    target: (p.target || [0, 0, 0]).slice(),
    distance: p.distance ?? 1000,
    yaw: p.yaw ?? 0,
    pitch: p.pitch ?? 0.4,
    fov: p.fov ?? 45,
  };
}

export function clonePose(p) {
  return {
    target: p.target.slice(),
    distance: p.distance,
    yaw: p.yaw,
    pitch: p.pitch,
    fov: p.fov,
  };
}

/** Direction the camera looks along for a (yaw, pitch) pair. */
export function forwardFromAngles(yaw, pitch, out = [0, 0, 0]) {
  const cp = Math.cos(pitch);
  out[0] = cp * Math.cos(yaw);
  out[1] = cp * Math.sin(yaw);
  out[2] = Math.sin(pitch);
  return out;
}

export function anglesFromForward(d) {
  const r = Math.hypot(d[0], d[1], d[2]) || 1;
  return {
    yaw: Math.atan2(d[1], d[0]),
    pitch: Math.asin(clamp(d[2] / r, -1, 1)),
  };
}

/** Eye position for a pose (the camera looks from eye toward target). */
export function poseEye(p, out = [0, 0, 0]) {
  const f = forwardFromAngles(p.yaw, p.pitch);
  out[0] = p.target[0] - f[0] * p.distance;
  out[1] = p.target[1] - f[1] * p.distance;
  out[2] = p.target[2] - f[2] * p.distance;
  return out;
}

/** Build a pose from an eye/target pair (legacy camera specs use these). */
export function poseFromEyeTarget(eye, target, fov = 45) {
  const d = [target[0] - eye[0], target[1] - eye[1], target[2] - eye[2]];
  const distance = Math.hypot(d[0], d[1], d[2]);
  const a = anglesFromForward(distance > 0 ? d : [1, 0, 0]);
  // A top-down eye (no horizontal offset) has no yaw: treat it as ±90°.
  if (Math.hypot(d[0], d[1]) <= 1e-9 * Math.max(distance, 1)) a.pitch = Math.sign(a.pitch || -1) * Math.PI / 2;
  return sanitizePose(makePose({ target, distance, yaw: a.yaw, pitch: a.pitch, fov }));
}

export function lerpPose(a, b, t, out = makePose()) {
  out.target[0] = lerp(a.target[0], b.target[0], t);
  out.target[1] = lerp(a.target[1], b.target[1], t);
  out.target[2] = lerp(a.target[2], b.target[2], t);
  // Interpolate distance logarithmically so zooms feel uniform.
  const da = Math.max(a.distance, 1e-3), db = Math.max(b.distance, 1e-3);
  out.distance = (a.distance <= 1e-3 && b.distance <= 1e-3) ? 0 : Math.exp(lerp(Math.log(da), Math.log(db), t));
  if (b.distance <= 1e-3 && t >= 1) out.distance = 0;
  out.yaw = a.yaw + angleDelta(a.yaw, b.yaw) * t;
  out.pitch = lerp(a.pitch, b.pitch, t);
  out.fov = lerp(a.fov, b.fov, t);
  return out;
}

export class Camera {
  constructor() {
    this.pose = makePose();
    this.aspect = 1;
    this.near = 0.5;
    this.far = 1e6;
    this.view = m4();
    this.proj = m4();
    this.viewProj = m4();
    this.invViewProj = m4();
    this.invView = m4();
    this.invProj = m4();
    this.eye = [0, 0, 0];
    this.forward = [1, 0, 0];
    this.right = [0, 1, 0];
    this.camUp = [0, 0, 1];
    this.width = 1;
    this.height = 1;
    // CSS px along each edge that a host page covers (a sidebar over an
    // embedded figure). Like a map's padding, the view centres in what is
    // left: the projection shifts off-axis while the pose stays put, so
    // States and links mean the same view at any padding.
    this.padding = { top: 0, right: 0, bottom: 0, left: 0 };
    this.version = 0;
  }

  /** Horizontal field of view in degrees (Aladin's FoV is the view width). */
  get hfov() {
    const v = this.pose.fov * DEG;
    return 2 * Math.atan(Math.tan(v / 2) * this.aspect) / DEG;
  }

  setViewport(width, height) {
    this.width = Math.max(1, width);
    this.height = Math.max(1, height);
    this.aspect = this.width / this.height;
  }

  /** Set the padding in CSS px (edges left out keep their value). */
  setPadding(p = {}) {
    for (const k of ["top", "right", "bottom", "left"]) {
      const v = Number(p[k]);
      if (Number.isFinite(v)) this.padding[k] = Math.max(0, v);
    }
  }

  /** The padding in use: a fifth of each axis always stays uncovered. */
  usedPadding() {
    const p = this.padding;
    const fit = (a, b, size) => {
      const k = Math.min(1, (size * 0.8) / Math.max(a + b, 1e-9));
      return [a * k, b * k];
    };
    const [left, right] = fit(p.left, p.right, this.width);
    const [top, bottom] = fit(p.top, p.bottom, this.height);
    return { top, right, bottom, left };
  }

  /** The optical centre in CSS px: the middle of the uncovered view. */
  get center() {
    const p = this.usedPadding();
    return [(this.width + p.left - p.right) / 2, (this.height + p.top - p.bottom) / 2];
  }

  /**
   * Where a TAN background (Aladin in Sky view) must sit to register with
   * this view: centred on the optical centre, large enough to cover the
   * window, with the horizontal field its width spans (CSS px, degrees).
   * Without padding it is the window and `hfov`.
   */
  skyFrame() {
    const [cx, cy] = this.center;
    const width = 2 * Math.max(cx, this.width - cx);
    const height = 2 * Math.max(cy, this.height - cy);
    const f = this.height / 2 / Math.tan((this.pose.fov * DEG) / 2);
    return { left: cx - width / 2, top: cy - height / 2, width, height, hfov: (2 * Math.atan(width / 2 / f)) / DEG };
  }

  update(sceneRadius = 1e4) {
    const p = this.pose;
    // Only rebuild (and bump `version`, which downstream caches key on) when
    // the pose or viewport actually changed.
    const [cx, cy] = this.center;
    const sig = `${p.target[0]},${p.target[1]},${p.target[2]},${p.distance},${p.yaw},${p.pitch},${p.fov},${this.aspect},${this.height},${sceneRadius},${cx},${cy}`;
    if (sig === this._sig) return false;
    this._sig = sig;
    poseEye(p, this.eye);
    forwardFromAngles(p.yaw, p.pitch, this.forward);
    const tgt = p.distance > 1e-6
      ? p.target
      : [this.eye[0] + this.forward[0], this.eye[1] + this.forward[1], this.eye[2] + this.forward[2]];
    m4lookAt(this.eye, tgt, UP, this.view);
    // Row vectors of the view rotation are the camera basis.
    this.right[0] = this.view[0]; this.right[1] = this.view[4]; this.right[2] = this.view[8];
    this.camUp[0] = this.view[1]; this.camUp[1] = this.view[5]; this.camUp[2] = this.view[9];
    // Depth range adapts to the zoom level so near structure never clips and
    // precision is spent where the viewer is looking.
    const reach = Math.max(p.distance, 1) + sceneRadius;
    this.near = p.distance > 1e-6 ? clamp(p.distance * 0.002, 0.01, 50) : 0.05;
    this.far = Math.max(reach * 2.5, 1e5);
    m4perspective(p.fov * DEG, this.aspect, this.near, this.far, this.proj);
    // Padding: shift the image so the view axis lands on the optical centre.
    this.proj[8] = 1 - (2 * cx) / this.width;
    this.proj[9] = (2 * cy) / this.height - 1;
    m4mul(this.proj, this.view, this.viewProj);
    m4invert(this.viewProj, this.invViewProj);
    m4invert(this.view, this.invView);
    m4invert(this.proj, this.invProj);
    this.version++;
    return true;
  }

  /** World → CSS pixel coordinates. Returns null if behind the camera. */
  project(p, out = [0, 0, 0]) {
    const c = m4project(this.viewProj, p);
    if (c[3] <= 0) return null;
    out[0] = (c[0] * 0.5 + 0.5) * this.width;
    out[1] = (1 - (c[1] * 0.5 + 0.5)) * this.height;
    out[2] = c[3];
    return out;
  }

  /** World-space ray through a CSS pixel. */
  ray(x, y) {
    const nx = (x / this.width) * 2 - 1;
    const ny = 1 - (y / this.height) * 2;
    const a = m4project(this.invViewProj, [nx, ny, -1]);
    const b = m4project(this.invViewProj, [nx, ny, 1]);
    const d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    const l = Math.hypot(d[0], d[1], d[2]) || 1;
    return { origin: this.eye.slice(), dir: [d[0] / l, d[1] / l, d[2] / l] };
  }

  /** World units per CSS pixel at distance `depth` from the eye. */
  pixelScale(depth) {
    return (2 * depth * Math.tan((this.pose.fov * DEG) / 2)) / this.height;
  }
}

/** Animate between two poses; returns a stepper object. */
export function poseTween(from, to, durationMs, ease = easeInOutCubic) {
  const start = performance.now();
  const a = clonePose(from), b = clonePose(to);
  // Arc the path out a little when the target jumps far, so flights read
  // as a camera move rather than a zoom through the data.
  const jump = Math.hypot(b.target[0] - a.target[0], b.target[1] - a.target[1], b.target[2] - a.target[2]);
  const lift = Math.min(jump * 0.35, Math.max(a.distance, b.distance) * 1.5);
  const out = makePose();
  return {
    durationMs,
    /** Move the destination with a moving camera anchor, so the flight lands on it. */
    shiftTarget(d) {
      b.target[0] += d[0]; b.target[1] += d[1]; b.target[2] += d[2];
    },
    step(now = performance.now()) {
      const raw = durationMs > 0 ? clamp((now - start) / durationMs, 0, 1) : 1;
      const t = ease(raw);
      lerpPose(a, b, t, out);
      if (lift > 0 && a.distance > 1e-3 && b.distance > 1e-3) {
        out.distance += Math.sin(Math.PI * t) * lift;
      }
      return { pose: out, done: raw >= 1 };
    },
  };
}
