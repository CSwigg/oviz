"""Focused contracts for the Paper and Mobile polish pass.

These tests intentionally exercise generated HTML rather than a stored main-
figure artifact.  They stay small, cover the integration seams that are easy
to regress, and do not require regenerating any of the large fixtures.
"""

from __future__ import annotations

import pytest

from oviz.paper import Paper
from oviz.threejs_figure import ThreeJSFigure


def _scene(*, mobile: bool = False, paper: dict | None = None) -> dict:
    scene = {
        "renderer": "threejs",
        "width": 640,
        "height": 480,
        "frames": [{"name": "0", "time": 0.0, "traces": []}],
        "initial_state": {"mobile_mode_enabled": mobile} if mobile else {},
        "states": {
            "schema_version": 1,
            "project_id": "paper-mobile-polish-tests",
            "items": [
                {
                    "id": "state-overview",
                    "name": "Overview",
                    "snapshot": {},
                }
            ],
        },
    }
    if paper is not None:
        scene["paper"] = paper
    return scene


def _paper_spec() -> dict:
    paper = Paper("Paper polish", enabled=True)
    paper.add_section("Introduction")
    paper.add_paragraph("Read along.", state="Overview", label="Overview")
    return paper.to_spec(
        states={"items": [{"id": "state-overview", "name": "Overview"}]}
    )


@pytest.fixture(scope="module")
def paper_html() -> str:
    return ThreeJSFigure(_scene(paper=_paper_spec())).to_html(
        compress_scene_spec=False
    )


@pytest.fixture(scope="module")
def mobile_html() -> str:
    return ThreeJSFigure(_scene(mobile=True)).to_html(compress_scene_spec=False)


@pytest.fixture(scope="module")
def desktop_html() -> str:
    return ThreeJSFigure(_scene()).to_html(compress_scene_spec=False)


def _function_region(source: str, function_name: str, next_function_name: str) -> str:
    start_marker = f"function {function_name}"
    end_marker = f"function {next_function_name}"
    start = source.index(start_marker)
    end = source.index(end_marker, start + len(start_marker))
    return source[start:end]


def test_paper_camera_offset_is_an_overlay_not_saved_camera_state(paper_html: str):
    """Opening Paper must not overwrite an authored State camera offset."""

    assert "let currentPaperCameraViewOffset" in paper_html
    assert "function combinedActionCameraViewOffset()" in paper_html
    assert "function applyActionCameraViewOffsetProjection()" in paper_html
    assert "function applyPaperCameraViewOffset(viewOffset)" in paper_html

    combined = _function_region(
        paper_html,
        "combinedActionCameraViewOffset()",
        "applyActionCameraViewOffsetProjection()",
    )
    assert "currentActionCameraViewOffset" in combined
    assert "currentPaperCameraViewOffset" in combined
    assert "logical.x + paper.x" in combined
    assert "logical.y + paper.y" in combined

    authored = _function_region(
        paper_html,
        "applyActionCameraViewOffset(viewOffset)",
        "applyPaperCameraViewOffset(viewOffset)",
    )
    overlay = _function_region(
        paper_html,
        "applyPaperCameraViewOffset(viewOffset)",
        "captureCurrentActionViewState()",
    )
    assert "currentActionCameraViewOffset =" in authored
    assert "currentPaperCameraViewOffset" not in authored
    assert "currentPaperCameraViewOffset =" in overlay
    assert "currentActionCameraViewOffset =" not in overlay
    assert "applyActionCameraViewOffsetProjection();" in authored
    assert "applyActionCameraViewOffsetProjection();" in overlay

    sky_geometry = _function_region(
        paper_html,
        "applySkyDomeFrameOffsetGeometry()",
        "updateSkyDomeBackgroundFrame(timestampMs = 0.0, options = {})",
    )
    assert 'typeof combinedActionCameraViewOffset === "function"' in sky_geometry
    assert "combinedActionCameraViewOffset()" in sky_geometry


def test_paper_scroll_retarget_is_last_wins_and_reports_status(paper_html: str):
    """Rapid reading-line changes coalesce and expose a useful status hook."""

    for contract in (
        "let ovizPaperJumpDebounceTimer",
        "let ovizPaperJumpGeneration",
        "function ovizPaperFlushStateJump(generation)",
        "function ovizPaperSetStatus(status, message = \"\")",
        "function ovizPaperSyncPanelAccessibility()",
        "function ovizPaperHandleTransitionLifecycle(event)",
        'class="oviz-three-paper-status"',
        'ovizPaperSetDataset("paperJumpPending"',
        'ovizPaperSetDataset("paperStatus"',
        'ovizPaperSetDataset("paperActiveAnchor"',
    ):
        assert contract in paper_html

    queue = _function_region(
        paper_html,
        "ovizPaperQueueStateJump(target)",
        "ovizPaperSampleTraces()",
    )
    assert "ovizPaperJumpGeneration += 1" in queue
    assert "window.clearTimeout(ovizPaperJumpDebounceTimer)" in queue
    assert "window.setTimeout" in queue
    assert "ovizPaperFlushStateJump(generation)" in queue

    flush = _function_region(
        paper_html,
        "ovizPaperFlushStateJump(generation)",
        "ovizPaperQueueStateJump(target)",
    )
    assert "generation !== ovizPaperJumpGeneration" in flush
    assert "ovizPaperPendingJumpTarget" in flush
    assert "await ovizGoToState" in flush

    for event_name in (
        "transition-complete",
        "transition-cancel",
        "transition-error",
    ):
        assert f'root.addEventListener("{event_name}"' in paper_html


def test_mobile_keeps_primary_authoring_actions_persistent(mobile_html: str):
    """Legend, Controls, and Lasso remain one tap away on mobile."""

    assert "oviz-three-earth-view-toggle oviz-three-view-segmented" in mobile_html
    assert 'class="oviz-three-search-shell"' in mobile_html
    assert 'class="oviz-three-mobile-legend"' in mobile_html
    assert 'class="oviz-three-mobile-controls"' in mobile_html
    assert 'class="oviz-three-mobile-lasso"' in mobile_html
    assert 'class="oviz-three-mobile-more"' in mobile_html
    assert 'data-mobile-action="search"' in mobile_html
    assert "min-height: 38px" in mobile_html or "height: 38px" in mobile_html

    assert '.oviz-three-widget-menu > * {' in mobile_html
    assert (
        '[data-search-open="true"] .oviz-three-widget-menu > .oviz-three-search-shell {'
        in mobile_html
    )
    assert (
        ".oviz-three-widget-menu > .oviz-three-mobile-legend," in mobile_html
    )
    assert (
        ".oviz-three-widget-menu > .oviz-three-mobile-controls," in mobile_html
    )
    assert (
        ".oviz-three-widget-menu > .oviz-three-mobile-lasso," in mobile_html
    )
    assert (
        ".oviz-three-widget-menu > .oviz-three-mobile-more {" in mobile_html
    )
    assert "grid-template-rows: 44px 44px" not in mobile_html
    assert "repeat(3, max-content)" not in mobile_html
    assert ".oviz-three-bottom-switches > * {" in mobile_html
    assert (
        ".oviz-three-bottom-switches > .oviz-three-earth-view-toggle {"
        in mobile_html
    )

    # The legacy Sky shortcut remains in the DOM for compatibility, while the
    # three requested authoring actions are the persistent mobile toolbar.
    assert (
        '[data-mobile="true"] .oviz-three-mobile-sky-view' in mobile_html
    )
    assert (
        '[data-mobile="true"] .oviz-three-mobile-lasso' in mobile_html
    )
    assert (
        '[data-mobile="true"] .oviz-three-mobile-legend' in mobile_html
    )


def test_mobile_panels_share_one_coordinated_sheet_and_sky_portal(mobile_html: str):
    for contract in (
        'class="oviz-three-mobile-sheet"',
        'class="oviz-three-mobile-sheet-backdrop"',
        'class="oviz-three-mobile-sheet-card"',
        'class="oviz-three-mobile-sheet-title"',
        'class="oviz-three-mobile-sheet-close"',
        'class="oviz-three-mobile-sheet-menu"',
        'class="oviz-three-mobile-sheet-content"',
        "function setMobileSheetPanel(panelName)",
        "function mountMobilePanelInSheet(panelName, panelEl)",
        "function restoreMobileSheetPanel()",
        "function syncMobileSheetAvailability()",
        "root.dataset.mobileSheetPanel",
    ):
        assert contract in mobile_html

    coordinator = _function_region(
        mobile_html,
        "setMobileSheetPanel(panelName)",
        "setManualLabelMenuOpen(isOpen)",
    )
    for panel_name in ("menu", "legend", "controls", "sky"):
        assert f'"{panel_name}"' in coordinator
    assert "mountMobilePanelInSheet" in coordinator
    assert "restoreMobileSheetPanel" in coordinator
    assert "skyControlsShellEl" in mobile_html
    assert "mountMobilePanelInSheet(nextName, skyControlsShellEl);" in coordinator


def test_mobile_volume_deferral_is_explicit_opt_in(mobile_html: str):
    assert (
        "const mobileDeferVolumes = mobileModeEnabled "
        "&& initialState.mobile_defer_volumes === true;"
    ) in mobile_html
    assert "initialState.mobile_defer_volumes !== false" not in mobile_html
    assert (
        'root.dataset.mobileVolumesDeferred = mobileDeferVolumes ? "true" : "false";'
        in mobile_html
    )


def test_desktop_keeps_existing_controls_and_ignores_mobile_sheet(desktop_html: str):
    assert 'data-mobile="false"' in desktop_html
    assert 'class="oviz-three-text-shell"' in desktop_html
    assert 'class="oviz-three-controls-shell"' in desktop_html
    assert 'class="oviz-three-reset-camera-view"' in desktop_html

    mount = _function_region(
        desktop_html,
        "mountMobilePanelInSheet(panelName, panelEl)",
        "syncMobileSheetAvailability()",
    )
    coordinator = _function_region(
        desktop_html,
        "setMobileSheetPanel(panelName)",
        "setManualLabelMenuOpen(isOpen)",
    )
    assert "!mobileModeEnabled" in mount
    assert "!mobileModeEnabled" in coordinator
    hidden_mobile_chrome = desktop_html[
        desktop_html.index(".oviz-three-mobile-sky-view,") :
        desktop_html.index(".oviz-three-fullscreen {")
    ]
    assert ".oviz-three-mobile-more," in hidden_mobile_chrome
    assert ".oviz-three-mobile-sheet," in hidden_mobile_chrome
    assert "display: none" in hidden_mobile_chrome
