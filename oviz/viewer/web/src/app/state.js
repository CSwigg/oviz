// Viewer state: defaults, legacy conversion, volume clamping.
//
// The state object is plain JSON so it can be captured into a State,
// serialized into an export, diffed and interpolated. Shape:
//   { version, group, time: {frame, speed}, view: {mode, pose},
//     traces: {key: {visible, opacity, sizeScale, colorMode, colormap, color}},
//     volumes: {stateKey: {...}}, images: {key: {visible, opacity}},
//     global: {...}, sky: {...} }

import { clamp } from "../core/math.js";
import { makePose, poseFromEyeTarget, clonePose } from "../engine/camera.js";

export const STATE_VERSION = 2;

export function defaultVolumeState(spec) {
  const d = spec.defaults || {};
  const num = (v, fb) => (Number.isFinite(Number(v)) ? Number(v) : fb);
  return {
    visible: spec.visible !== false,
    opacity: num(d.opacity, 1),
    vmin: num(d.vmin, spec.dataRange[0]),
    vmax: num(d.vmax, spec.dataRange[1]),
    steps: num(d.steps, 192),
    alphaCoef: num(d.alpha_coef, 50),
    stretch: String(d.stretch || "linear"),
    colormap: String(d.colormap || spec.colormaps[0] || "inferno"),
    showAllTimes: Boolean(d.show_all_times),
    lightingMode: normalizeLighting(d.lighting_mode),
    galacticCenter: Array.isArray(d.galactic_center) ? d.galactic_center.map(Number) : [8122, 0, 0],
    galacticLightIntensity: num(d.galactic_light_intensity, 1.35),
    galacticAmbient: num(d.galactic_ambient, 0.22),
    galacticExtinction: num(d.galactic_extinction, 2.4),
    galacticScattering: num(d.galactic_scattering, 0.55),
    galacticAnisotropy: num(d.galactic_anisotropy, 0.45),
    galacticWarmth: num(d.galactic_warmth, 0.72),
  };
}

function normalizeLighting(mode) {
  const m = String(mode || "standard").toLowerCase().replace(/-/g, "_");
  return ["galactic", "galactic_center", "galactic_centre"].includes(m) ? "galactic" : "standard";
}

/** Port of the legacy clampVolumeStateForLayer (returns a clamped copy). */
export function clampVolumeState(spec, vs) {
  const out = { ...vs };
  const [dmin, dmax] = spec.dataRange || [0, 1];
  const d = spec.defaults || {};
  if (!Number.isFinite(out.vmin)) out.vmin = Number(d.vmin);
  if (!Number.isFinite(out.vmax)) out.vmax = Number(d.vmax);
  if (Number.isFinite(dmin) && Number.isFinite(dmax)) {
    out.vmin = clamp(out.vmin, dmin, dmax);
    out.vmax = clamp(out.vmax, dmin, dmax);
  }
  if (!(out.vmax > out.vmin)) out.vmax = Number.isFinite(dmax) && dmax > out.vmin ? dmax : out.vmin + 1e-6;
  out.opacity = clamp(Number(out.opacity), 0, 1);
  out.steps = Math.round(clamp(Number(out.steps), 24, 768));
  if (!Number.isFinite(out.steps)) out.steps = Number(d.steps || 100);
  out.alphaCoef = clamp(Number(out.alphaCoef), 1, 200);
  if (!Number.isFinite(out.alphaCoef)) out.alphaCoef = Number(d.alpha_coef || 50);
  return out;
}

function groupDefaults(manifest, group) {
  return (manifest.groups?.visibility || {})[group] || null;
}

/** Legend membership + visibility of a key in a group (legacy semantics). */
export function groupMode(manifest, group, key) {
  const defaults = groupDefaults(manifest, group);
  if (!defaults) return true;
  const mode = defaults[key];
  if (mode === undefined) return false;
  return mode;
}

export function initialViewerState(manifest) {
  const init = manifest.initialState || {};
  const gc = init.global_controls || {};
  const anim = manifest.animation || {};
  const groups = manifest.groups || {};
  const group = String(init.current_group || groups.default || groups.order?.[0] || "");
  const state = {
    version: STATE_VERSION,
    group,
    time: { frame: manifest.time?.initialIndex ?? 0, speed: 1 },
    view: { mode: "3d", pose: makePose() },
    traces: {},
    volumes: {},
    images: {},
    global: {
      pointSize: num(gc.point_size_scale, 1),
      pointOpacity: num(gc.point_opacity_scale, 1),
      glow: num(gc.point_glow_strength, 0.6),
      fadeTime: num(gc.fade_in_time_myr, anim.fadeTimeMyr ?? 0),
      fadeInOut: bool(gc.fade_in_and_out_enabled, anim.fadeInOut),
      fadeByOpacity: bool(gc.fade_opacity_by_birth_time_enabled, anim.fadeByOpacity),
      sizeByStars: bool(gc.size_points_by_stars_enabled, false),
      grid: bool(gc.galactic_reference_visible, init.galactic_coordinates_visible ?? true),
      gridOpacity: 1,
      labels: true,
      autoOrbit: bool(gc.camera_auto_orbit_enabled, false),
      trails: 0,
    },
    sky: {
      layers: cloneJson(manifest.sky?.layers || []),
      backgroundVisible: !gc.sky_background_hidden,
      members: gc.sky_member_display_mode === "clusters" ? "clusters" : "stars",
    },
  };
  applyGroupDefaults(state, manifest, group);
  const legend = init.legend_state || {};
  for (const [key, v] of Object.entries(legend)) {
    if (manifest.traces?.some((t) => t.key === key)) {
      (state.traces[key] ||= {}).visible = v === true;
    }
  }
  // Trace styles captured by the legacy runtime.
  for (const [key, ts] of Object.entries(init.trace_style_state || {})) {
    const trace = manifest.traces?.find((t) => t.key === key);
    if (!trace) continue;
    Object.assign((state.traces[key] ||= {}), legacyTraceStyle(ts, trace));
  }
  // Volumes: defaults, then legacy per-key state, then legend visibility.
  const seen = new Set();
  for (const spec of manifest.volumes || []) {
    if (seen.has(spec.stateKey)) continue;
    seen.add(spec.stateKey);
    const vs = defaultVolumeState(spec);
    const legacy = (init.volume_state_by_key || {})[spec.stateKey] || {};
    Object.assign(vs, pickVolumeFields(legacy));
    if (Object.prototype.hasOwnProperty.call(legend, spec.stateKey)) vs.visible = legend[spec.stateKey] === true;
    state.volumes[spec.stateKey] = vs;
  }
  for (const img of manifest.images || []) {
    state.images[img.key] = { visible: init.galaxy_image !== false, opacity: 1 };
  }
  if (Number.isFinite(Number(init.current_frame_value))) state.time.frame = Number(init.current_frame_value);
  else if (Number.isFinite(Number(init.current_frame_index))) state.time.frame = Number(init.current_frame_index);
  return state;
}

export function applyGroupDefaults(state, manifest, group) {
  state.group = group;
  for (const trace of manifest.traces || []) {
    const mode = groupMode(manifest, group, trace.key);
    const t = (state.traces[trace.key] ||= {});
    t.visible = mode === true;
    t.inGroup = mode !== false;
  }
}

function legacyTraceStyle(ts, trace) {
  const out = {};
  if (ts.opacity != null && Number.isFinite(Number(ts.opacity))) out.opacity = Number(ts.opacity);
  if (ts.sizeScale != null && Number.isFinite(Number(ts.sizeScale))) out.sizeScale = Number(ts.sizeScale);
  if (ts.colorMode) out.colorMode = ts.colorMode === "by_value" ? "by_value" : "fixed";
  if (ts.colormap) out.colormap = String(ts.colormap);
  if (ts.color && String(ts.color).toLowerCase() !== String(trace.color).toLowerCase()) {
    // Only an explicit user edit overrides per-object colours.
    if (ts.color !== (trace.legendColor || trace.color)) out.color = String(ts.color);
  }
  return out;
}

function pickVolumeFields(v) {
  const out = {};
  const keys = [
    "visible", "opacity", "vmin", "vmax", "steps", "alphaCoef", "stretch", "colormap", "showAllTimes",
    "lightingMode", "galacticCenter", "galacticLightIntensity", "galacticAmbient", "galacticExtinction",
    "galacticScattering", "galacticAnisotropy", "galacticWarmth",
  ];
  for (const k of keys) {
    if (v[k] === undefined || v[k] === null) continue;
    out[k] = k === "lightingMode" ? normalizeLighting(v[k]) : v[k];
  }
  return out;
}

/**
 * Convert a legacy State snapshot (captureRuntimeState()) into a v2 state,
 * layered over `base` so fields the snapshot lacks keep their values.
 */
export function legacySnapshotToState(snap, manifest, base) {
  const st = cloneJson(base);
  st.version = STATE_VERSION;
  if (snap.current_group) applyGroupDefaults(st, manifest, String(snap.current_group));
  for (const [key, v] of Object.entries(snap.legend_state || {})) {
    if (st.volumes[key]) st.volumes[key].visible = v === true;
    else if (manifest.traces?.some((t) => t.key === key)) (st.traces[key] ||= {}).visible = v === true;
  }
  for (const [key, ts] of Object.entries(snap.trace_style_state || {})) {
    const trace = manifest.traces?.find((t) => t.key === key);
    if (trace) Object.assign((st.traces[key] ||= {}), legacyTraceStyle(ts, trace));
  }
  for (const [key, vs] of Object.entries(snap.volume_state_by_key || {})) {
    if (st.volumes[key]) Object.assign(st.volumes[key], pickVolumeFields(vs));
  }
  const gc = snap.global_controls || {};
  const g = st.global;
  if (gc.point_size_scale != null) g.pointSize = num(gc.point_size_scale, g.pointSize);
  if (gc.point_opacity_scale != null) g.pointOpacity = num(gc.point_opacity_scale, g.pointOpacity);
  if (gc.point_glow_strength != null) g.glow = num(gc.point_glow_strength, g.glow);
  if (gc.fade_in_time_myr != null) g.fadeTime = num(gc.fade_in_time_myr, g.fadeTime);
  if (gc.fade_in_and_out_enabled != null) g.fadeInOut = !!gc.fade_in_and_out_enabled;
  if (gc.fade_opacity_by_birth_time_enabled != null) g.fadeByOpacity = !!gc.fade_opacity_by_birth_time_enabled;
  if (gc.size_points_by_stars_enabled != null) g.sizeByStars = !!gc.size_points_by_stars_enabled;
  if (gc.galactic_reference_visible != null) g.grid = !!gc.galactic_reference_visible;
  if (gc.camera_auto_orbit_enabled != null) g.autoOrbit = !!gc.camera_auto_orbit_enabled;
  if (gc.sky_member_display_mode) st.sky.members = gc.sky_member_display_mode === "clusters" ? "clusters" : "stars";
  if (gc.sky_background_hidden != null) st.sky.backgroundVisible = !gc.sky_background_hidden;
  if (Array.isArray(snap.sky_layers) && snap.sky_layers.length) st.sky.layers = cloneJson(snap.sky_layers);
  if (Array.isArray(snap.manual_labels)) {
    st.ext = { ...(st.ext || {}), notes: snap.manual_labels.map(legacyNote).filter(Boolean) };
  }
  const fv = Number(snap.current_frame_value ?? snap.current_frame_index);
  if (Number.isFinite(fv)) st.time.frame = fv;
  const mode = gc.camera_view_mode === "earth" ? "sky" : "3d";
  st.view.mode = mode;
  const cam = snap.camera;
  if (cam && cam.position && cam.target) {
    const p = [cam.position.x, cam.position.y, cam.position.z].map(Number);
    const t = [cam.target.x, cam.target.y, cam.target.z].map(Number);
    if (p.every(Number.isFinite) && t.every(Number.isFinite)) {
      const fov = num(gc.camera_fov, st.view.pose.fov || 60);
      if (mode === "sky") {
        const pose = poseFromEyeTarget(p, t, fov);
        pose.target = [0, 0, 0];
        pose.distance = 0;
        st.view.pose = pose;
      } else {
        st.view.pose = poseFromEyeTarget(p, t, fov);
      }
    }
  }
  return st;
}

export function cloneState(s) {
  const out = cloneJson(s);
  if (s.view?.pose) out.view.pose = clonePose(s.view.pose);
  return out;
}

export function cloneJson(v) {
  return v == null ? v : JSON.parse(JSON.stringify(v));
}

function num(v, fb) {
  const n = Number(v);
  return Number.isFinite(n) ? n : fb;
}

function bool(v, fb) {
  return v === undefined || v === null ? !!fb : !!v;
}

function legacyNote(l, i) {
  if (!l || typeof l !== "object") return null;
  const x = Number(l.x ?? l.position?.x), y = Number(l.y ?? l.position?.y), z = Number(l.z ?? l.position?.z);
  if (![x, y, z].every(Number.isFinite)) return null;
  return { id: String(l.id || `legacy-${i}`), text: String(l.text || "Note"), anchor: { pos: [x, y, z] }, color: String(l.color || "#ffffff") };
}
