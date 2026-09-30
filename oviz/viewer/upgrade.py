"""Upgrade existing Oviz HTML figures to the new viewer.

Every legacy figure embeds its complete scene spec, either as a gzip+base64
payload split across ``<script type="application/octet-stream">`` chunks or
as an inline ``const sceneSpec = {...};`` literal. This module extracts that
spec, compiles it into a bundle and writes a new standalone figure, so old
figures (and their saved States) move to the new viewer without re-running
any science pipeline.

Command line::

    python -m oviz.viewer.upgrade old_figure.html new_figure.html
"""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
import sys
import time
from pathlib import Path
from typing import Any

_PAYLOAD_ID_RE = re.compile(r"readOvizSceneSpecPayload\((.*?)\)")
_INLINE_RE = re.compile(r"/\*__SCENE_SPEC_START__\*/const sceneSpec = (.*?);/\*__SCENE_SPEC_END__\*/", re.S)


def read_legacy_scene_spec(path: str | Path) -> dict[str, Any]:
    """Return the scene spec embedded in a legacy Oviz HTML figure."""

    html = Path(path).read_text(encoding="utf-8")
    payload_id = None
    for arg in _PAYLOAD_ID_RE.findall(html):
        try:
            candidate = json.loads(arg)
        except json.JSONDecodeError:
            continue
        if isinstance(candidate, str) and candidate:
            payload_id = candidate
            break
    if payload_id:
        chunks = re.findall(
            r"<script\b(?=[^>]*type=[\"']application/octet-stream[\"'])"
            r"(?=[^>]*data-oviz-payload-id=[\"']" + re.escape(payload_id) + r"[\"'])"
            r"[^>]*data-oviz-payload-index=[\"'](\d+)[\"'][^>]*>(.*?)</script>",
            html,
            re.S,
        )
        if not chunks:
            # Older figures store the payload in one element keyed by id.
            single = re.search(
                r"<script\b(?=[^>]*type=[\"']application/octet-stream[\"'])"
                r"(?=[^>]*\bid=[\"']" + re.escape(payload_id) + r"[\"'])[^>]*>(.*?)</script>",
                html,
                re.S,
            )
            if not single:
                raise ValueError(f"{path}: compressed payload chunks for {payload_id!r} not found")
            chunks = [("0", single.group(1))]
        encoded = "".join(re.sub(r"\s+", "", c) for _, c in sorted(chunks, key=lambda x: int(x[0])))
        return json.loads(gzip.decompress(base64.b64decode(encoded)))
    m = _INLINE_RE.search(html)
    if m:
        return json.loads(m.group(1))
    raise ValueError(f"{path}: no embedded Oviz scene spec found")


def upgrade_html(src: str | Path, dst: str | Path, *, verbose: bool = False, mode: str = "focus",
                 camera_anchor: str | None = None, title: str | None = None) -> dict[str, Any]:
    """Convert the legacy figure ``src`` into a new-viewer figure at ``dst``.

    ``mode`` and ``camera_anchor`` are the :class:`~oviz.viewer.figure.OvizFigure`
    options (the anchor defaults to the LSR when the home view orbits it).
    """
    from .figure import OvizFigure
    from .compile import compile_scene_spec

    t0 = time.perf_counter()
    spec = read_legacy_scene_spec(src)
    t1 = time.perf_counter()
    bundle = compile_scene_spec(spec)
    t2 = time.perf_counter()
    fig = OvizFigure(bundle=bundle, mode=mode, camera_anchor=camera_anchor, title=title)
    out = fig.write_html(dst)
    t3 = time.perf_counter()
    report = {
        "source_bytes": Path(src).stat().st_size,
        "output_bytes": out.stat().st_size,
        "decode_s": round(t1 - t0, 2),
        "compile_s": round(t2 - t1, 2),
        "write_s": round(t3 - t2, 2),
        **bundle.size_report(),
    }
    if verbose:
        mb = 1 / (1024 * 1024)
        print(
            f"{src} ({report['source_bytes'] * mb:.1f} MiB) → {dst} ({report['output_bytes'] * mb:.1f} MiB) "
            f"in {t3 - t0:.1f}s"
        )
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Upgrade a legacy Oviz HTML figure to the new viewer.")
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=("focus", "detailed"), default="focus", help="mode the figure opens in")
    parser.add_argument("--camera-anchor", choices=("lsr", "sun", "free"), default=None,
                        help="what the camera orbits and moves with (default: the LSR when the home view orbits it)")
    args = parser.parse_args(argv)
    upgrade_html(args.source, args.output, verbose=True, mode=args.mode, camera_anchor=args.camera_anchor)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
