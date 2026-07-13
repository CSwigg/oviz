from __future__ import annotations

import hashlib

from oviz.threejs_figure import ThreeJSFigure


KNOWN_GOOD_STARTUP_HASHES = {
    "apply_initial": "d1e2211940073214e2409142a115f5522865530b04518097c695e7a1df7f8532",
    "marker": "e1ff90b3eb153b326af817e0eeb4eed8d79ca3cc5a3fa051fd16534bd7e890ff",
    "frame_scene": "11d826e0224e45691ccd6ee340ff4624e4d1b592a30279b28ff8402b1f0701d7",
    "render_frame": "3a6d2a6f01c07569d990dd7283f9e8fdf749de9c912acb35e4ae89c49860125a",
}


def _runtime_html() -> str:
    return ThreeJSFigure(
        {
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }
    ).to_html(compress_scene_spec=False)


def _digest(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def test_transition_optimizations_do_not_change_known_good_startup_renderer():
    html = _runtime_html()
    startup_sections = {
        "apply_initial": html.split(
            "function applyViewerStateSyncInternal(initialState, options = {})", 1
        )[1].split("function captureWidgetState", 1)[0],
        "marker": html.split("function addMarkerTrace(parent, trace)", 1)[1].split(
            "function addTextTrace", 1
        )[0],
        "frame_scene": html.split(
            "function renderFrameScene(frame, displayedTimeMyr, options = {})", 1
        )[1].split("let ovizRetainedTransitionScene", 1)[0],
        "render_frame": html.split("function renderFrame(index)", 1)[1].split(
            "function updateTopbarDensity", 1
        )[0],
    }

    assert {key: _digest(value) for key, value in startup_sections.items()} == (
        KNOWN_GOOD_STARTUP_HASHES
    )
