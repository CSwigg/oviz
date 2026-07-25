from __future__ import annotations

import hashlib

from oviz.threejs_figure import ThreeJSFigure


KNOWN_GOOD_STARTUP_HASHES = {
    "apply_initial": "03882b710aa2be5e479db7da4bfe59b96c134df27af35c903d4377242c996412",
    "marker": "c28a7d68b5b4dcf53c361811e495a8ee5c55b6b00ea8df922ca8742c6439579c",
    "frame_scene": "6844dc6f70c05bb8cff935ff33b0d2dedac9775b416d1bc856dd794ee3eae07e",
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


def test_startup_renderer_matches_visually_reviewed_optimized_runtime():
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
