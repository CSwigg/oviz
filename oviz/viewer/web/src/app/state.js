// Viewer state: defaults, legacy conversion, volume clamping.
//
// The state object is plain JSON so it can be captured into a State,
// serialized into an export, diffed and interpolated. Shape:
//   { version, group, time: {frame, speed}, view: {mode, pose, anchor},
//     traces: {key: {visible, opacity, sizeScale, colorMode, colormap, color}},
//     volumes: {stateKey: {...}}, images: {key: {visible, opacity}},
//     global: {...}, sky: {...} }

import { clamp } from "../core/math.js";
import { makePose, poseFromEyeTarget, clonePose } from "../engine/camera.js";

export const STATE_VERSION = 2;

function defaultVolumeState(spec) {
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
    // Event-density volumes: the trailing window their KDE counts (Myr).
    ...(spec.procedural ? { kdeTimeWindowMyr: num(spec.procedural.windowMyr, 10) } : {}),
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

/**
 * A volume's place in a group, the way a trace has one: true (shown),
 * "legendonly" (listed but hidden) or false (not in the group). Groups that
 * list no volume at all (figures grouped before they had volumes) leave
 * every volume listed and as it was: undefined.
 */
export function volumeGroupMode(manifest, group, spec) {
  const defaults = groupDefaults(manifest, group);
  if (!defaults || !spec) return undefined;
  const listed = (x) => defaults[x.stateKey] !== undefined || defaults[x.key] !== undefined;
  if (!(manifest.volumes || []).some(listed)) return undefined;
  const mode = defaults[spec.stateKey] ?? defaults[spec.key];
  return mode === undefined ? false : mode;
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
    view: { mode: "3d", pose: makePose(), anchor: { kind: "free" } }, // the viewer sets pose and anchor
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
      autoOrbitSpeed: Math.max(0, num(gc.camera_auto_orbit_speed, 1)),
      autoOrbitDirection: orbitDirection(gc.camera_auto_orbit_direction),
      navSpeed: clamp(num(gc.scroll_speed, 1), 0.2, 4),
      // The lasso also clips volumes (classic "Lasso volumetric data").
      lassoVolumes: bool(init.lasso_volume_selection_enabled, true),
      axes: bool(gc.axes_visible, manifest.world?.showAxes),
      trails: 0,
    },
    sky: {
      // Classic sky groups filter what is drawn: the layers outside the
      // active group start hidden, so the upgraded figure opens the same.
      layers: applyLegacySkyGroup(cloneJson(manifest.sky?.layers || []), manifest.sky?.layerGroups, init.active_sky_layer_group_key ?? gc.active_sky_layer_group_key ?? manifest.sky?.activeGroupKey),
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
    // Classic shows a volume only when both its state and its legend entry are on.
    if (Object.prototype.hasOwnProperty.call(legend, spec.stateKey)) vs.visible = vs.visible !== false && legend[spec.stateKey] === true;
    state.volumes[spec.stateKey] = vs;
  }
  for (const img of manifest.images || []) {
    state.images[img.key] = { visible: init.galaxy_image !== false, opacity: 1 };
  }
  if (Number.isFinite(Number(init.current_frame_value))) state.time.frame = Number(init.current_frame_value);
  else if (Number.isFinite(Number(init.current_frame_index))) state.time.frame = Number(init.current_frame_index);
  return state;
}

/**
 * Show or hide a set of layers at once, the way the classic legend's T key
 * and section eyes do. `items` are layer ids, `visible` the ids now shown.
 * If any is shown, all are hidden and the shown ones are remembered (the
 * returned `snapshot`). Otherwise the remembered ones come back, or every
 * layer when nothing (of these) was remembered. Returns
 * `{ visible: Map<id, boolean>, snapshot: Set<id> | null }`.
 */
export function toggleAllPlan(items, visible, snapshot) {
  const on = items.filter((id) => visible.has(id));
  if (on.length) return { visible: new Map(items.map((id) => [id, false])), snapshot: new Set(on) };
  const restore = snapshot && items.some((id) => snapshot.has(id)) ? snapshot : null;
  return { visible: new Map(items.map((id) => [id, restore ? restore.has(id) : true])), snapshot: null };
}

/**
 * Show one layer alone (double-click a name, ⇧ + its number), the way the
 * classic legend did: every other layer, volumes included, goes off and the
 * mix that was on is remembered. Doing it again on the soloed layer brings
 * that mix back (with the soloed layer kept on), or every layer when
 * nothing was remembered. Returns `{ visible, snapshot }` like
 * `toggleAllPlan`.
 */
export function soloPlan(items, visible, id, snapshot) {
  if (!items.includes(id)) return { visible: new Map(items.map((x) => [x, visible.has(x)])), snapshot: snapshot || null };
  const alone = items.every((x) => visible.has(x) === (x === id));
  if (alone) {
    const restore = snapshot && items.some((x) => x !== id && snapshot.has(x)) ? snapshot : null;
    return { visible: new Map(items.map((x) => [x, x === id || (restore ? restore.has(x) : true)])), snapshot: null };
  }
  return { visible: new Map(items.map((x) => [x, x === id])), snapshot: new Set(items.filter((x) => visible.has(x))) };
}

export function applyGroupDefaults(state, manifest, group) {
  state.group = group;
  for (const trace of manifest.traces || []) {
    const mode = groupMode(manifest, group, trace.key);
    const t = (state.traces[trace.key] ||= {});
    t.visible = mode === true;
    t.inGroup = mode !== false;
  }
  // Volumes are layers like the traces: a group that lists them shows or
  // hides them too, as the classic legend did.
  for (const spec of manifest.volumes || []) {
    const vs = state.volumes?.[spec.stateKey];
    const mode = volumeGroupMode(manifest, group, spec);
    if (vs && mode !== undefined) vs.visible = mode === true;
  }
}

/** The classic sky group a layer belongs to: its own group_key, else the first group listing its survey. */
function legacySkyGroupOf(layer, groups) {
  const slug = (x) => String(x || "").trim().toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-+|-+$/g, "");
  const own = slug(layer.group_key ?? layer.groupKey);
  if (own && own !== "all" && (own === "other" || groups.some((g) => g.key === own))) return own;
  const survey = String(layer.survey || layer.key || "").trim();
  const g = groups.find((d) => d.key !== "all" && (d.surveys || []).some((c) => String(c).trim() === survey));
  return g ? g.key : "other";
}

/**
 * Classic Sky layer groups drew only the active group's layers. The new
 * viewer has no groups in its stack, so a classic figure or State is read
 * with the other groups' layers hidden (the default group was "all-sky").
 */
export function applyLegacySkyGroup(layers, groups, activeKey) {
  const defs = Array.isArray(groups) ? groups.filter((g) => g && g.key) : [];
  if (!defs.length || !Array.isArray(layers)) return layers;
  const key = String(activeKey || (defs.some((g) => g.key === "all-sky") ? "all-sky" : "all"));
  if (key === "all" || !defs.some((g) => g.key === key)) return layers;
  return layers.map((l) => (legacySkyGroupOf(l, defs) === key ? l : { ...l, visible: false }));
}

function legacyTraceStyle(ts, trace) {
  const out = {};
  if (ts.opacity != null && Number.isFinite(Number(ts.opacity))) out.opacity = Number(ts.opacity);
  if (ts.sizeScale != null && Number.isFinite(Number(ts.sizeScale))) out.sizeScale = Number(ts.sizeScale);
  if (ts.colorMode) out.colorMode = ts.colorMode === "by_value" ? "by_value" : "fixed";
  if (ts.colormap) out.colormap = String(ts.colormap);
  // Classic Sky member motion arrows.
  if (ts.memberArrowsEnabled != null) out.memberArrows = !!ts.memberArrowsEnabled;
  if (Number.isFinite(Number(ts.memberArrowLength))) out.memberArrowLength = Number(ts.memberArrowLength);
  if (Number.isFinite(Number(ts.memberArrowWidth))) out.memberArrowWidth = Number(ts.memberArrowWidth);
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
    "galacticScattering", "galacticAnisotropy", "galacticWarmth", "kdeTimeWindowMyr",
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
  for (const [key, vs] of Object.entries(snap.volume_state_by_key || {})) {
    if (st.volumes[key]) Object.assign(st.volumes[key], pickVolumeFields(vs));
  }
  for (const [key, v] of Object.entries(snap.legend_state || {})) {
    // A volume shows only when its legend entry and its own state are both on.
    if (st.volumes[key]) st.volumes[key].visible = v === true && snap.volume_state_by_key?.[key]?.visible !== false;
    else if (manifest.traces?.some((t) => t.key === key)) (st.traces[key] ||= {}).visible = v === true;
  }
  for (const [key, ts] of Object.entries(snap.trace_style_state || {})) {
    const trace = manifest.traces?.find((t) => t.key === key);
    if (trace) Object.assign((st.traces[key] ||= {}), legacyTraceStyle(ts, trace));
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
  if (gc.camera_auto_orbit_speed != null) g.autoOrbitSpeed = Math.max(0, num(gc.camera_auto_orbit_speed, 1));
  if (gc.scroll_speed != null) g.navSpeed = clamp(num(gc.scroll_speed, 1), 0.2, 4);
  if (snap.lasso_volume_selection_enabled != null) g.lassoVolumes = !!snap.lasso_volume_selection_enabled;
  if (gc.axes_visible != null) g.axes = !!gc.axes_visible;
  if (gc.camera_auto_orbit_direction != null) g.autoOrbitDirection = orbitDirection(gc.camera_auto_orbit_direction);
  if (gc.sky_member_display_mode) st.sky.members = gc.sky_member_display_mode === "clusters" ? "clusters" : "stars";
  if (gc.sky_background_hidden != null) st.sky.backgroundVisible = !gc.sky_background_hidden;
  if (Array.isArray(snap.sky_layers) && snap.sky_layers.length) {
    const group = snap.active_sky_layer_group_key ?? gc.active_sky_layer_group_key;
    st.sky.layers = applyLegacySkyGroup(cloneJson(snap.sky_layers), manifest.sky?.layerGroups, group);
  }
  if (Array.isArray(snap.manual_labels)) {
    st.ext = { ...(st.ext || {}), notes: snap.manual_labels.map(legacyNote).filter(Boolean) };
  }
  // A classic click selection names its cluster; it is found by name on apply.
  const cs = snap.current_selection;
  const clicked = cs && typeof cs === "object" && cs.cluster_name ? { traceName: String(cs.trace_name || ""), name: String(cs.cluster_name) } : null;
  if (clicked || "current_selection" in snap || "zen_mode_enabled" in snap) {
    st.ext = { ...(st.ext || {}), ui: { selection: clicked, layers: null, details: null, zen: !!snap.zen_mode_enabled } };
  }
  // So does a lasso selection: its clusters' names (and its outline, which clips the dust).
  if ("selected_cluster_keys" in snap || "current_selection_mode" in snap) st.ext = { ...(st.ext || {}), lasso: legacyLasso(snap) };
  if ("cluster_filter_state" in snap) st.ext = { ...(st.ext || {}), filter: legacyFilter(snap.cluster_filter_state, manifest.widgets?.clusterFilter?.extents) };
  if ("widgets" in snap) st.ext = { ...(st.ext || {}), widgets: legacyWidgets(snap.widgets, snap.dendrogram_state) };
  if ("sky_aperture" in gc || "sky_apertures" in gc) {
    const lens = legacyLens(gc.sky_aperture ?? (Array.isArray(gc.sky_apertures) ? gc.sky_apertures[0] : null));
    if (lens) st.sky.lens = lens;
    else delete st.sky.lens;
  }
  const fv = Number(snap.current_frame_value ?? snap.current_frame_index);
  if (Number.isFinite(fv)) st.time.frame = fv;
  // Classic playback_state: direction 0 is paused, ±1 plays that way.
  const dir = Number(snap.playback_state?.direction);
  if (dir) { st.time.playing = true; st.time.direction = dir < 0 ? -1 : 1; }
  else { delete st.time.playing; delete st.time.direction; }
  const mode = gc.camera_view_mode === "earth" ? "sky" : "3d";
  st.view.mode = mode;
  // Classic States had no camera anchor: applying one infers it from the pose.
  delete st.view.anchor;
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
        // The 3D view V returns to (classic earth_view_return_camera_state).
        const back = gc.earth_view_return_camera_state;
        const bp = back?.position, bt = back?.target;
        if (bp && bt) {
          const e = [bp.x, bp.y, bp.z].map(Number), c = [bt.x, bt.y, bt.z].map(Number);
          if (e.every(Number.isFinite) && c.every(Number.isFinite)) st.view.returnPose = poseFromEyeTarget(e, c, num(back.fov, fov));
        }
      } else {
        st.view.pose = poseFromEyeTarget(p, t, fov);
      }
    }
  }
  return st;
}

/**
 * The Sky view a figure opens in, or null when it opens in 3D. A classic
 * "earth" initial view (a sky-start export) gives the gaze and field of
 * view, and the 3D view that V returns to.
 */
export function skyStartView(manifest, base) {
  const init = manifest.initialState || {};
  if (init.global_controls?.camera_view_mode !== "earth" || !manifest.sky?.enabled) return null;
  const st = legacySnapshotToState(init, manifest, base);
  if (st.view.mode !== "sky") return null;
  return { pose: st.view.pose, returnPose: st.view.returnPose || null };
}

export function cloneState(s) {
  const out = cloneJson(s);
  if (s.view?.pose) out.view.pose = clonePose(s.view.pose);
  return out;
}

export function cloneJson(v) {
  return v == null ? v : JSON.parse(JSON.stringify(v));
}

/**
 * Auto-orbit rate in rad/s added to the camera yaw. Classic Oviz orbited
 * like three.js OrbitControls at autoRotateSpeed 1.2 × speed: 2π/60 × 1.2
 * rad/s, clockwise seen from Galactic north (the yaw decreasing), with
 * direction −1 reversing it.
 */
export function autoOrbitRate(g) {
  if (!g || !g.autoOrbit) return 0;
  const speed = Number.isFinite(Number(g.autoOrbitSpeed)) ? Math.max(0, Number(g.autoOrbitSpeed)) : 1;
  return -((2 * Math.PI) / 60) * 1.2 * speed * (g.autoOrbitDirection < 0 ? -1 : 1);
}

function orbitDirection(v) {
  const s = String(v ?? "").trim().toLowerCase();
  if (["reverse", "clockwise", "cw", "-1", "negative"].includes(s)) return -1;
  return Number(v) < 0 ? -1 : 1;
}

function num(v, fb) {
  const n = Number(v);
  return Number.isFinite(n) ? n : fb;
}

function bool(v, fb) {
  return v === undefined || v === null ? !!fb : !!v;
}

/** The classic cluster-name key: lower case, underscores and spaces as one space. */
export function selectionKey(name) {
  return String(name ?? "").trim().toLowerCase().replace(/[_\s]+/g, " ");
}

/**
 * A classic lasso selection as ext.lasso: the selected clusters by name
 * (resolved to objects on apply) and the outline in normalised device
 * coordinates with the camera it was drawn with. Classic lassos hid the
 * rest while their filter was on; with it off the selection only stays.
 */
function legacyLasso(snap) {
  const names = [...new Set((Array.isArray(snap.selected_cluster_keys) ? snap.selected_cluster_keys : []).map(selectionKey).filter(Boolean))];
  const m = snap.lasso_selection_mask;
  const polygon = Array.isArray(m?.polygon_ndc) ? m.polygon_ndc.map((p) => [Number(p?.x), Number(p?.y)]).filter((p) => p.every(Number.isFinite)) : [];
  const vp = Array.isArray(m?.view_projection_matrix) ? m.view_projection_matrix.map(Number) : [];
  const mask = polygon.length >= 3 && vp.length === 16 && vp.every(Number.isFinite) ? { polygon, vp } : null;
  const mode = String(snap.current_selection_mode || (names.length || mask ? "lasso" : "none"));
  if (mode !== "lasso" || (!names.length && !mask)) return null;
  return { names, isolate: true, filterOff: snap.lasso_selection_filter_enabled === false, mask };
}

// Classic Cluster Filter parameters and the distribution filter's quantities.
const LEGACY_FILTER_PARAMS = { age_now_myr: { param: "ageNow" }, n_stars: { param: "nStars", log: true } };

/**
 * A classic Cluster Filter state as ext.filter. Only the active parameter
 * filtered, hiding what fell outside its range; a range spanning the whole
 * parameter (`extents`, from the classic widget) is no filter.
 */
function legacyFilter(cf, extents) {
  if (!cf || typeof cf !== "object") return null;
  const key = String(cf.parameter_key || "");
  const p = LEGACY_FILTER_PARAMS[key];
  const r = cf.ranges_by_key?.[key];
  if (!p || !r || typeof r !== "object") return null;
  const ext = Array.isArray(extents?.[key]) ? extents[key].map(Number) : null;
  let lo = r.min == null ? NaN : Number(r.min);
  let hi = r.max == null ? NaN : Number(r.max);
  if (ext && ext.every(Number.isFinite)) {
    lo = Number.isFinite(lo) ? clamp(lo, ext[0], ext[1]) : ext[0];
    hi = Number.isFinite(hi) ? clamp(hi, ext[0], ext[1]) : ext[1];
    if (lo <= ext[0] && hi >= ext[1]) return null;
  }
  if (!Number.isFinite(lo) || !Number.isFinite(hi)) return null;
  if (lo > hi) lo = hi = (lo + hi) / 2;
  // The filter keeps log-scaled quantities' ranges in log10.
  const tf = (x) => (p.log ? Math.log10(Math.max(x, 1e-9)) : x);
  return { param: p.param, range: [tf(lo), tf(hi)], mode: "hide" };
}

// Classic widget panels and the viewer widgets they became.
const LEGACY_WIDGETS = { age_kde: "relative-sfh", dendrogram: "birth-tree" };

/** Classic open widget panels (and the Dendrogram's settings) as ext.widgets. */
function legacyWidgets(w, dendrogram) {
  const open = {};
  for (const [classic, id] of Object.entries(LEGACY_WIDGETS)) {
    const p = w?.[classic];
    if (!p || typeof p !== "object" || !["normal", "fullscreen"].includes(p.mode)) continue;
    const r = p.rect || {};
    const box = [r.left, r.top, r.width, r.height].map((x) => (x == null ? NaN : Number(x)));
    open[id] = {
      ...(box.every(Number.isFinite) ? { rect: { x: box[0], y: box[1], w: box[2], h: box[3] } } : {}),
      state: id === "birth-tree" ? legacyBirthTree(dendrogram) : null,
    };
  }
  return Object.keys(open).length ? { open } : null;
}

function legacyBirthTree(d) {
  if (!d || typeof d !== "object") return null;
  const out = {};
  if (d.trace_key) out.trace = String(d.trace_key);
  if (d.threshold_mode) out.mode = d.threshold_mode === "birth_age_myr" ? "birth_age_myr" : "distance_pc";
  if (d.connection_mode) out.links = d.connection_mode === "birth_to_birth" ? "birth_to_birth" : "birth_to_older_track";
  if (d.threshold_pc != null && Number.isFinite(Number(d.threshold_pc))) out.thresholdPc = Number(d.threshold_pc);
  if (d.threshold_age_myr != null && Number.isFinite(Number(d.threshold_age_myr))) out.thresholdAgeMyr = Number(d.threshold_age_myr);
  return out;
}

/** A classic sky aperture as the Sky lens: its spot, size and the band it showed most. */
function legacyLens(ap) {
  if (!ap || typeof ap !== "object" || ap.open === false) return null;
  const l = Number(ap.sky_center?.l), b = Number(ap.sky_center?.b);
  if (ap.sky_center?.l == null || ap.sky_center?.b == null || !Number.isFinite(l) || !Number.isFinite(b)) return null;
  const blend = ap.blend && typeof ap.blend === "object" ? ap.blend : null;
  const band = String(ap.band_key || "");
  const survey = (blend && (Number(blend.fraction) >= 0.5 ? blend.upper_survey : blend.lower_survey)) || (band.includes("/") ? band : "") || "P/PLANCK/R2/HFI/color";
  const size = Number(ap.size_deg);
  return { l, b, sizeDeg: Number.isFinite(size) && size > 0 ? size : 14, survey: String(survey), opacity: clamp(num(ap.opacity, 1), 0, 1) };
}

function legacyNote(l, i) {
  if (!l || typeof l !== "object") return null;
  const x = Number(l.x ?? l.position?.x), y = Number(l.y ?? l.position?.y), z = Number(l.z ?? l.position?.z);
  if (![x, y, z].every(Number.isFinite)) return null;
  const size = Number(l.size ?? l.font_size);
  return {
    id: String(l.id || `legacy-${i}`), text: String(l.text || "Note"), anchor: { pos: [x, y, z] }, color: String(l.color || "#ffffff"),
    ...(Number.isFinite(size) && size > 0 ? { size } : {}),
    ...(l.visible === false ? { visible: false } : {}),
  };
}
