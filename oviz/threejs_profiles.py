from __future__ import annotations

import copy


def _deep_merge_dicts(base, overrides):
    merged = copy.deepcopy(base)
    for key, value in (overrides or {}).items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge_dicts(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def normalize_threejs_initial_state(initial_state=None):
    state = copy.deepcopy(initial_state or {})

    if "lite_mode_enabled" not in state and "minimal_mode_enabled" in state:
        state["lite_mode_enabled"] = bool(state.get("minimal_mode_enabled"))

    if "galactic_lite_mode_enabled" not in state and "galactic_simple_mode_enabled" in state:
        state["galactic_lite_mode_enabled"] = bool(state.get("galactic_simple_mode_enabled"))

    if "galactic_plane_image_path" in state:
        state["galaxy_image_path"] = state.get("galactic_plane_image_path")

    if "galactic_plane_size_pc" in state:
        state["galaxy_image_size_pc"] = state.get("galactic_plane_size_pc")

    if "galactic_plane_opacity" in state:
        state["galaxy_image_opacity"] = state.get("galactic_plane_opacity")

    if "galaxy_image" not in state:
        state["galaxy_image"] = bool(state.get("galaxy_image_path"))

    return state


def merge_threejs_profile(base_state=None, override_state=None):
    return normalize_threejs_initial_state(_deep_merge_dicts(base_state or {}, override_state or {}))


def threejs_profile(**overrides):
    return merge_threejs_profile({}, overrides)


def lite_profile(**overrides):
    return merge_threejs_profile({"lite_mode_enabled": True}, overrides)


def galactic_lite_profile(
    *,
    galaxy_image=False,
    galaxy_image_path=None,
    galaxy_image_size_pc=40000.0,
    galaxy_image_opacity=0.6,
    track_orbit_target_to_sun=False,
    initial_zoom_anchor=None,
    camera=None,
    global_controls=None,
    **overrides,
):
    profile = {
        "lite_mode_enabled": True,
        "galactic_lite_mode_enabled": True,
        "galaxy_image": bool(galaxy_image),
        "track_orbit_target_to_sun": bool(track_orbit_target_to_sun),
    }
    if galaxy_image_path is not None:
        profile["galaxy_image_path"] = str(galaxy_image_path)
    if galaxy_image_size_pc is not None:
        profile["galaxy_image_size_pc"] = float(galaxy_image_size_pc)
    if galaxy_image_opacity is not None:
        profile["galaxy_image_opacity"] = float(galaxy_image_opacity)
    if initial_zoom_anchor is not None:
        profile["initial_zoom_anchor"] = copy.deepcopy(initial_zoom_anchor)
    if camera is not None:
        profile["camera"] = copy.deepcopy(camera)
    if global_controls is not None:
        profile["global_controls"] = copy.deepcopy(global_controls)
    return merge_threejs_profile(profile, overrides)


def galactic_simple_profile(**overrides):
    return galactic_lite_profile(**overrides)


def website_background_profile(
    *,
    galaxy_image=False,
    galaxy_image_path=None,
    camera=None,
    global_controls=None,
    **overrides,
):
    base_camera = {
        "position": {"x": 500.0, "y": -4600.0, "z": 3300.0},
        "target": {"x": 0.0, "y": 0.0, "z": 0.0},
        "view_offset": {"x": 0.22, "y": 0.0},
        "up": {"x": 0.0, "y": 0.0, "z": 1.0},
    }
    base_controls = {
        "camera_fov": 60.0,
        "point_size_scale": 1.0,
        "point_glow_strength": 0.8,
        "fade_in_time_myr": 7.0,
        "camera_auto_orbit_enabled": True,
        "camera_auto_orbit_speed": 0.45,
        "camera_auto_orbit_direction": -1,
    }
    return galactic_lite_profile(
        galaxy_image=galaxy_image,
        galaxy_image_path=galaxy_image_path,
        track_orbit_target_to_sun=True,
        initial_zoom_anchor={"x": 0.0, "y": 0.0, "z": 0.0},
        click_selection_enabled=False,
        camera=_deep_merge_dicts(base_camera, camera or {}),
        global_controls=_deep_merge_dicts(base_controls, global_controls or {}),
        **overrides,
    )


_PROFILE_BUILDERS = {
    "full": threejs_profile,
    "standard": threejs_profile,
    "threejs": threejs_profile,
    "lite": lite_profile,
    "galactic_lite": galactic_lite_profile,
    "lite_galactic": galactic_lite_profile,
    "galactic_simple": galactic_simple_profile,
    "website_background": website_background_profile,
}


def build_threejs_profile(name=None, **overrides):
    if name is None:
        return threejs_profile(**overrides)

    normalized_name = str(name).strip().lower().replace("-", "_")
    builder = _PROFILE_BUILDERS.get(normalized_name)
    if builder is None:
        raise ValueError(
            "Unknown threejs profile. Expected one of: "
            + ", ".join(sorted(_PROFILE_BUILDERS))
        )
    return builder(**overrides)
