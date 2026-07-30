"""Focused regression contracts for the coordinated mobile viewer shell."""

from __future__ import annotations

import re

import pytest

from oviz.threejs_figure import ThreeJSFigure


def _scene() -> dict:
    return {
        "renderer": "threejs",
        "width": 640,
        "height": 480,
        "frames": [{"name": "0", "time": 0.0, "traces": []}],
        "initial_state": {"mobile_mode_enabled": True},
        "states": {
            "schema_version": 1,
            "project_id": "mobile-shell-regressions",
            "items": [
                {
                    "id": "state-overview",
                    "name": "Overview",
                    "snapshot": {},
                }
            ],
        },
    }


@pytest.fixture(scope="module")
def mobile_html() -> str:
    return ThreeJSFigure(_scene()).to_html(compress_scene_spec=False)


def _function_region(source: str, name: str, next_name: str) -> str:
    start = source.index(f"function {name}")
    end = source.index(f"function {next_name}", start + len(name))
    return source[start:end]


def test_states_is_a_real_sheet_panel_and_resyncs_after_async_build(
    mobile_html: str,
):
    assert 'data-mobile-panel="states"' in mobile_html
    assert 'data-mobile-action="states"' not in mobile_html

    coordinator = _function_region(
        mobile_html,
        "setMobileSheetPanel(panelName)",
        "setManualLabelMenuOpen(isOpen)",
    )
    assert '"states"' in coordinator
    assert 'mountMobilePanelInSheet(nextName, statesShellEl);' in coordinator
    assert 'statesDrawerEl.removeAttribute("inert")' in coordinator
    assert "focusMobileSheetPanel(nextName);" in coordinator

    initialize_states = mobile_html.split(
        "async function initializeOvizStates()", 1
    )[1].split("initializeOvizStates().then", 1)[0]
    build_position = initialize_states.index("ovizBuildStatesDrawer();")
    sync_position = initialize_states.index("syncMobileSheetAvailability();")
    assert build_position < sync_position


def test_legacy_sky_shell_cannot_escape_the_mobile_sheet_or_zen(
    mobile_html: str,
):
    escaped_sky_rule = re.compile(
        r'\[data-mobile="true"\]\[data-camera-view-mode="earth"\]\s+'
        r'\.oviz-three-sky-controls-shell\[data-visible="true"\]\s*\{'
    )
    assert escaped_sky_rule.search(mobile_html) is None
    assert (
        ".oviz-three-mobile-sheet-content > .oviz-three-sky-controls-shell"
        in mobile_html
    )
    assert (
        '[data-mobile="true"][data-zen="true"] .oviz-three-mobile-sheet {'
        in mobile_html
    )


def test_mobile_legend_is_touch_scrollable_and_includes_sky_layers(
    mobile_html: str,
):
    assert (
        '.oviz-three-mobile-sheet[data-panel="legend"] .oviz-three-mobile-sheet-body {'
        in mobile_html
    )
    assert "overflow-y: auto !important;" in mobile_html
    assert "touch-action: pan-y !important;" in mobile_html
    assert (
        '.oviz-three-mobile-sheet[data-panel="legend"] '
        '.oviz-three-sky-controls-shell[data-visible="true"] {'
        in mobile_html
    )
    assert (
        '.oviz-three-mobile-sheet[data-panel="legend"] '
        '.oviz-three-sky-controls-shell[data-visible="false"] {'
        in mobile_html
    )

    coordinator = _function_region(
        mobile_html,
        "setMobileSheetPanel(panelName)",
        "setManualLabelMenuOpen(isOpen)",
    )
    assert 'cameraViewMode === "earth"' in coordinator
    assert 'setSkyControlsDrawerOpen(true, { keepLegendOpen: true });' in coordinator

    sky_drawer = _function_region(
        mobile_html,
        "setSkyControlsDrawerOpen(isOpen, options = {})",
        "syncSkyBackgroundDockVisibility()",
    )
    assert "options.keepLegendOpen !== true" in sky_drawer


def test_mobile_text_entry_does_not_trigger_ios_focus_zoom_or_page_refresh(
    mobile_html: str,
):
    # Do not solve input zoom by disabling user pinch zoom in the viewport.
    viewport = re.search(
        r'<meta name="viewport" content="([^"]+)"', mobile_html
    )
    assert viewport is not None
    assert "maximum-scale" not in viewport.group(1)
    assert "user-scalable=no" not in viewport.group(1)

    for selector in (
        'input[type="text"]',
        'input[type="search"]',
        'input[type="email"]',
        'input[type="url"]',
        'input[type="tel"]',
        'input[type="password"]',
        'input[type="number"]',
        "textarea",
        "select",
    ):
        assert f'#__ROOT_ID__[data-mobile="true"] {selector}' not in mobile_html
        assert f'[data-mobile="true"] {selector}' in mobile_html

    mobile_text_rule = mobile_html.split(
        '[data-mobile="true"] input[type="text"]', 1
    )[1].split("}", 1)[0]
    assert "font-size: 16px !important;" in mobile_text_rule
    assert "touch-action: manipulation;" in mobile_text_rule
    assert "scroll-margin-block: 18px;" in mobile_text_rule

    page_rule = mobile_html.split("html, body {", 1)[1].split("}", 1)[0]
    assert "overflow: hidden;" in page_rule
    assert "overscroll-behavior: none;" in page_rule


def test_mobile_refresh_guard_preserves_controls_and_scrollable_panels(
    mobile_html: str,
):
    guard = _function_region(
        mobile_html,
        "installMobileRefreshGuard()",
        "syncMobileViewportGeometry()",
    )
    assert 'root.addEventListener("touchmove"' in guard
    assert "{ passive: false }" in guard
    assert "mobileNativeTouchTarget(event.target)" in guard
    assert "mobileScrollableTouchAncestor(event.target)" in guard
    assert "deltaY > 0 && atTop" in guard
    assert "deltaY < 0 && atBottom" in guard
    assert "event.preventDefault();" in guard

    sync_mobile = _function_region(
        mobile_html,
        "syncMobileStaticUi()",
        "onWindowMessage(event)",
    )
    assert "installMobileRefreshGuard();" in sync_mobile


def test_mobile_zen_and_presentation_close_the_sheet_and_keep_touch_exit(
    mobile_html: str,
):
    zen = _function_region(
        mobile_html,
        "setZenMode(enabled)",
        "setPresentationMode(enabled)",
    )
    presentation = _function_region(
        mobile_html,
        "setPresentationMode(enabled)",
        "ovizNativeFullscreenElement()",
    )
    assert 'setMobileSheetPanel("closed");' in zen
    assert 'setMobileSheetPanel("closed");' in presentation
    assert 'presentation: ["Exit", "Exit presentation mode"]' in mobile_html
    assert (
        '[data-mobile="true"][data-presentation-mode="true"] '
        ".oviz-three-topbar {"
        in mobile_html
    )
    assert (
        '[data-mobile="true"][data-presentation-mode="true"] '
        ".oviz-three-mobile-more {"
        in mobile_html
    )

    more_handler = mobile_html.split(
        'mobileMoreButtonEl.addEventListener("click", (event) => {', 1
    )[1].split("if (mobileSheetBackdropEl)", 1)[0]
    assert more_handler.index("if (presentationModeEnabled)") < more_handler.index(
        "if (zenModeEnabled)"
    )
    assert "setPresentationMode(false);" in more_handler


def test_mobile_sheet_transfers_traps_and_restores_focus(mobile_html: str):
    assert (
        'class="oviz-three-mobile-sheet-card" role="dialog" aria-modal="true" '
        in mobile_html
    )
    assert 'tabindex="-1"' in mobile_html
    for helper in (
        "focusMobileSheetElement(element)",
        "mobileSheetFocusableElements(containerEl = mobileSheetEl)",
        "focusMobileSheetPanel(panelName)",
        "trapMobileSheetFocus(event)",
    ):
        assert f"function {helper}" in mobile_html

    coordinator = _function_region(
        mobile_html,
        "setMobileSheetPanel(panelName)",
        "setManualLabelMenuOpen(isOpen)",
    )
    assert "mobileSheetRestoreFocusEl = activeElement" in coordinator
    assert "focusMobileSheetElement(restoreFocusEl || mobileMoreButtonEl);" in coordinator
    assert "focusMobileSheetPanel(nextName);" in coordinator
    assert 'mobileSheetEl.addEventListener("keydown", trapMobileSheetFocus);' in mobile_html


def test_mobile_topbar_keeps_legend_controls_and_lasso_reachable(
    mobile_html: str,
):
    for class_name in (
        "oviz-three-mobile-legend",
        "oviz-three-mobile-controls",
        "oviz-three-mobile-lasso",
    ):
        assert f'class="{class_name}"' in mobile_html
        assert (
            f".oviz-three-widget-menu > .{class_name}" in mobile_html
        )

    controls_handler = mobile_html.split(
        'mobileControlsButtonEl.addEventListener("click", (event) => {', 1
    )[1].split("if (mobileArButtonEl)", 1)[0]
    assert 'mobileSheetPanelName === "controls"' in controls_handler
    assert ': "controls"' in controls_handler

    legend_handler = mobile_html.split(
        'mobileLegendButtonEl.addEventListener("click", (event) => {', 1
    )[1].split('legendPanelEl.addEventListener("click"', 1)[0]
    assert 'mobileSheetPanelName === "legend"' in legend_handler
    assert ': "legend"' in legend_handler

    lasso_handler = mobile_html.split(
        'mobileLassoButtonEl.addEventListener("click", () => {', 1
    )[1].split("if (mobileControlsButtonEl)", 1)[0]
    assert "lassoArmed = !lassoArmed;" in lasso_handler

    assert 'data-mobile-action="search">Search clusters</button>' in mobile_html
    mobile_menu_handler = mobile_html.split(
        'mobileSheetMenuEl.addEventListener("click", (event) => {', 1
    )[1].split('if (mobileModeEnabled) {', 1)[0]
    assert 'actionName === "search"' in mobile_menu_handler
    assert 'setOvizSearchOpen(true, { focus: true });' in mobile_menu_handler
