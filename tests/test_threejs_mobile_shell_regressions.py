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
