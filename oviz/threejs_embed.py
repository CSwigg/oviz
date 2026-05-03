from __future__ import annotations

import base64
import gzip
import json
from typing import Any


DEFAULT_SCENE_SPEC_COMPRESSION_THRESHOLD_BYTES = 8 * 1024 * 1024
SCENE_SPEC_COMPRESSED_CHUNK_SIZE = 64 * 1024


def _normalize_scene_spec_compression_mode(compress_scene_spec: bool | str | None) -> bool | str:
    if compress_scene_spec is None:
        return "auto"
    if compress_scene_spec is True or compress_scene_spec is False:
        return bool(compress_scene_spec)
    if isinstance(compress_scene_spec, str):
        normalized = compress_scene_spec.strip().lower()
        if normalized == "auto":
            return "auto"
        if normalized == "true":
            return True
        if normalized == "false":
            return False
    raise ValueError('compress_scene_spec must be True, False, or "auto".')


def _normalize_scene_spec_compression_threshold(threshold: int | None) -> int:
    if threshold is None:
        return DEFAULT_SCENE_SPEC_COMPRESSION_THRESHOLD_BYTES
    normalized = int(threshold)
    if normalized < 0:
        raise ValueError("scene_spec_compression_threshold_bytes must be >= 0.")
    return normalized


def _scene_spec_json(scene_spec: dict[str, Any]) -> str:
    compact_payload = bool((scene_spec.get("initial_state") or {}).get("compact_payload_enabled"))
    return (
        json.dumps(scene_spec, separators=(",", ":"))
        if compact_payload
        else json.dumps(scene_spec)
    )


def _scene_spec_payload(
    scene_spec: dict[str, Any],
    *,
    root_id: str,
    compress_scene_spec: bool | str | None,
    scene_spec_compression_threshold_bytes: int | None,
) -> tuple[str, dict[str, Any], str]:
    scene_json = _scene_spec_json(scene_spec)
    scene_json_bytes = scene_json.encode("utf-8")
    raw_size = len(scene_json_bytes)
    compression_mode = _normalize_scene_spec_compression_mode(compress_scene_spec)
    compression_threshold = _normalize_scene_spec_compression_threshold(
        scene_spec_compression_threshold_bytes
    )
    should_compress = (
        compression_mode is True
        or (compression_mode == "auto" and raw_size > compression_threshold)
    )
    metadata: dict[str, Any] = {
        "compression_mode": compression_mode,
        "compression_method": "none",
        "compressed": False,
        "raw_scene_spec_size_bytes": raw_size,
        "compressed_size_bytes": None,
        "embedded_base64_size_bytes": None,
        "compression_threshold_bytes": compression_threshold,
    }

    if not should_compress:
        return scene_json, metadata, ""

    compressed_bytes = gzip.compress(scene_json_bytes, compresslevel=9, mtime=0)
    encoded_scene = base64.b64encode(compressed_bytes).decode("ascii")
    encoded_chunks = [
        encoded_scene[idx : idx + SCENE_SPEC_COMPRESSED_CHUNK_SIZE]
        for idx in range(0, len(encoded_scene), SCENE_SPEC_COMPRESSED_CHUNK_SIZE)
    ]
    payload_id = f"{root_id}-scene-spec-payload"
    payload_html = (
        f'<script type="application/octet-stream" id="{payload_id}">\n'
        + "\n".join(encoded_chunks)
        + "\n</script>"
    )
    metadata.update(
        {
            "compression_method": "gzip+base64",
            "compressed": True,
            "compressed_size_bytes": len(compressed_bytes),
            "embedded_base64_size_bytes": len(encoded_scene.encode("ascii")),
        }
    )
    return (
        f"await inflateOvizGzipBase64SceneSpec(readOvizSceneSpecPayload({json.dumps(payload_id)}))",
        metadata,
        payload_html,
    )


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
    action_runtime_js: str,
    compress_scene_spec: bool | str | None = "auto",
    scene_spec_compression_threshold_bytes: int | None = None,
) -> str:
    width = int(scene_spec.get("width") or 900)
    height = int(scene_spec.get("height") or 700)
    initial_state = scene_spec.get("initial_state") or {}
    minimal_mode = bool(
        initial_state.get("lite_mode_enabled")
        or
        initial_state.get("minimal_mode_enabled")
        or scene_spec.get("export_profile") in {"minimal", "lite"}
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
    html = html.replace("__ACTION_RUNTIME_JS__", action_runtime_js)
    scene_spec_expr, scene_spec_payload_metadata, scene_spec_payload_html = _scene_spec_payload(
        scene_spec,
        root_id=root_id,
        compress_scene_spec=compress_scene_spec,
        scene_spec_compression_threshold_bytes=scene_spec_compression_threshold_bytes,
    )
    html = html.replace("__SCENE_SPEC_PAYLOAD_HTML__", scene_spec_payload_html)
    html = html.replace(
        "__SCENE_SPEC_PAYLOAD_METADATA__",
        json.dumps(scene_spec_payload_metadata, separators=(",", ":")),
    )
    html = html.replace("__SCENE_SPEC_EXPR__", scene_spec_expr)
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
