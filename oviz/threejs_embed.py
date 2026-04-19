from __future__ import annotations

import base64
import json
from typing import Any


def render_threejs_html(
    scene_spec: dict[str, Any],
    *,
    root_id: str,
    html_template: str,
    topbar_html: str,
    minimal_topbar_html: str,
    shell_html: str,
    legend_runtime_js: str,
    widget_runtime_js: str,
    widget_content_runtime_js: str,
    interaction_runtime_js: str,
    scene_runtime_js: str,
    sky_runtime_js: str,
    viewer_runtime_js: str,
) -> str:
    width = int(scene_spec.get("width") or 900)
    height = int(scene_spec.get("height") or 700)
    initial_state = scene_spec.get("initial_state") or {}
    minimal_mode = bool(
        initial_state.get("minimal_mode_enabled")
        or scene_spec.get("export_profile") == "minimal"
    )
    galactic_simple = bool((scene_spec.get("galactic_simple") or {}).get("enabled"))

    html = html_template.replace("__ROOT_ID__", root_id)
    html = html.replace("__WIDTH_PX__", str(width))
    html = html.replace("__HEIGHT_PX__", str(height))
    html = html.replace(
        "__ROOT_MINIMAL_ATTR__",
        'data-minimal="true"' if minimal_mode else 'data-minimal="false"',
    )
    html = html.replace(
        "__ROOT_GALACTIC_SIMPLE_ATTR__",
        'data-galactic-simple="true"' if galactic_simple else 'data-galactic-simple="false"',
    )
    html = html.replace(
        "__TOPBAR_HTML__",
        minimal_topbar_html if minimal_mode else topbar_html,
    )
    html = html.replace("__SHELL_HTML__", shell_html)
    html = html.replace("__LEGEND_RUNTIME_JS__", legend_runtime_js)
    html = html.replace("__WIDGET_RUNTIME_JS__", widget_runtime_js)
    html = html.replace("__WIDGET_CONTENT_RUNTIME_JS__", widget_content_runtime_js)
    html = html.replace("__INTERACTION_RUNTIME_JS__", interaction_runtime_js)
    html = html.replace("__SCENE_RUNTIME_JS__", scene_runtime_js)
    html = html.replace("__SKY_RUNTIME_JS__", sky_runtime_js)
    html = html.replace("__VIEWER_RUNTIME_JS__", viewer_runtime_js)
    html = html.replace("__SCENE_JSON__", json.dumps(scene_spec))
    return html


def threejs_data_url(html: str) -> str:
    encoded = base64.b64encode(html.encode("utf-8")).decode("ascii")
    return f"data:text/html;charset=utf-8;base64,{encoded}"


def threejs_iframe_html(data_url: str) -> str:
    return (
        f'<iframe src="{data_url}" '
        'style="width: 100vw; height: 100vh; border: 0; display: block; max-width: none; margin: 0 0 0 calc(50% - 50vw);" '
        'loading="eager" referrerpolicy="no-referrer"></iframe>'
    )
