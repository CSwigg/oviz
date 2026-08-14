#!/usr/bin/env python3
"""Regenerate the July 25 non-paper Oviz figure with current cluster velocities.

The source scene is rebuilt with the July 4 Chronos ages and a selectable
cluster-velocity CSV, then converted to the clean July presentation while
keeping the editable Slides authoring option available. The embedded Three.js
runtime carries the volumetric quality and FPS upgrades developed for the
July 25 artifact.
"""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
from astropy.io import fits

from oviz.threejs_figure import ThreeJSFigure


TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from main_figure_chronos_july4 import (  # noqa: E402
    CLUSTER_MEMBERS_PATH,
    JULY4_CHRONOS_RESULTS_PATH,
)
from main_figure_new_chronos import (  # noqa: E402
    RATZENBOECK_SCOCEN_TRACE_NAME,
    run_main_figure as build_source_main_figure,
    run_main_figure_scene as build_source_main_figure_scene,
)


SOURCE_HTML = Path(__file__).with_name("main_figure_chronos_july4.html")
DEFAULT_OUTPUT_HTML = Path(__file__).with_suffix(".html")
DEFAULT_CLUSTER_VELOCITIES_PATH = Path(
    "/Users/swiggumc/Desktop/astro_research/cfa/velocity analysis/outputs/release/"
    "cluster_velocities_sdssv_covariance_audited.csv"
)
DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH = Path(
    "/Users/swiggumc/Desktop/astro_research/radcliffe/cluster_data/files_formatted/"
    "ratzi_2022.csv"
)
DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH = Path(
    "/Users/swiggumc/Downloads/"
    "ScoCen_SigMA2_groups_averages_34+KERR3_CenFar-to-Cham_37clsuters.fits"
)
DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH = Path(
    "/Users/swiggumc/Downloads/ScoCen_Sept2022.fits"
)
RATZENBOECK_SCOCEN_MEMBER_COUNT = 13_103
JULY25_LOOKBACK_MYR = 60
EDENHOFER_GALACTIC_LIGHTING = {
    "lighting_mode": "standard",
    "galactic_center": [8122.0, 0.0, 0.0],
    "galactic_light_intensity": 1.35,
    "galactic_ambient": 0.22,
    "galactic_extinction": 2.4,
    "galactic_scattering": 0.55,
    "galactic_anisotropy": 0.45,
    "galactic_warmth": 0.72,
}
SKY_BACKGROUND_GROUPS = [
    {
        "key": "all",
        "label": "All",
        "surveys": [],
    },
    {
        "key": "all-sky",
        "label": "All-Sky",
        "surveys": [
            "P/Mellinger/color",
            "P/DSS2/color",
            "P/PLANCK/R2/HFI/color",
            "P/Finkbeiner",
            "P/HI4PI/NHI",
            "P/Fermi/color",
        ],
    },
    {
        "key": "optical-surveys",
        "label": "Optical Surveys",
        "surveys": [
            "P/PanSTARRS/DR1/color-z-zg-g",
            "CDS/P/SDSS9/color",
            "P/DECaLS/DR5/color",
        ],
    },
    {
        "key": "infrared-ultraviolet",
        "label": "Infrared & UV",
        "surveys": [
            "P/2MASS/color",
            "P/allWISE/color",
            "CDS/P/IRIS/color",
            "CDS/P/AKARI/FIS/Color",
            "P/GALEXGR6/AIS/color",
        ],
    },
]
SKY_BACKGROUND_LAYERS = [
    {
        "key": "P/Mellinger/color",
        "label": "Mellinger Color",
        "survey": "P/Mellinger/color",
        "opacity": 1.0,
        "visible": True,
    },
    {
        "key": "P/DSS2/color",
        "label": "DSS2 Color",
        "survey": "P/DSS2/color",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "P/PanSTARRS/DR1/color-z-zg-g",
        "label": "Pan-STARRS DR1 Color",
        "survey": "P/PanSTARRS/DR1/color-z-zg-g",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "CDS/P/SDSS9/color",
        "label": "SDSS9 Color",
        "survey": "CDS/P/SDSS9/color",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "P/DECaLS/DR5/color",
        "label": "DECaLS DR5 Color",
        "survey": "P/DECaLS/DR5/color",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "P/2MASS/color",
        "label": "2MASS Color",
        "survey": "P/2MASS/color",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "P/allWISE/color",
        "label": "AllWISE Color",
        "survey": "P/allWISE/color",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "P/PLANCK/R2/HFI/color",
        "label": "Planck Dust Emission Color",
        "survey": "P/PLANCK/R2/HFI/color",
        "opacity": 1.0,
        "visible": True,
        "stretch": "asinh",
        "cut_max": 0.07,
    },
    {
        "key": "CDS/P/IRIS/color",
        "label": "IRIS Far IR",
        "survey": "CDS/P/IRIS/color",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "CDS/P/AKARI/FIS/Color",
        "label": "AKARI Far IR",
        "survey": "CDS/P/AKARI/FIS/Color",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "P/GALEXGR6/AIS/color",
        "label": "GALEX AIS Color",
        "survey": "P/GALEXGR6/AIS/color",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "P/Finkbeiner",
        "label": "Finkbeiner Hα",
        "survey": "P/Finkbeiner",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "P/HI4PI/NHI",
        "label": "HI4PI HI 21 cm",
        "survey": "P/HI4PI/NHI",
        "opacity": 1.0,
        "visible": False,
    },
    {
        "key": "P/Fermi/color",
        "label": "Fermi Gamma-ray",
        "survey": "P/Fermi/color",
        "opacity": 1.0,
        "visible": False,
    },
]
SKY_BACKGROUND_SURVEY_ALIASES = {
    "P/SDSS9/color": "CDS/P/SDSS9/color",
}
# Presentation defaults for the July 25 figure: Planck dust sits on top of
# the stack but starts hidden (opacity preset for when it is toggled on),
# Mellinger renders dimmed, and DSS2 is enabled as the second visible layer.
SKY_BACKGROUND_STACK_TOP_ORDER = [
    "P/PLANCK/R2/HFI/color",
    "P/Mellinger/color",
    "P/DSS2/color",
]
SKY_BACKGROUND_STACK_DEFAULTS = {
    "P/PLANCK/R2/HFI/color": {"visible": False, "opacity": 0.27},
    "P/Mellinger/color": {"visible": True, "opacity": 0.23},
    "P/DSS2/color": {"visible": True, "opacity": 1.0},
}
# The full catalog underlay is dropped from the July 25 presentation so the
# figure shows only the curated cluster samples.
FULL_CLUSTER_CATALOG_TRACE_NAME = "Full Cluster Catalog"
# Single-PNG volume atlases larger than this per side exceed iOS Safari's
# canvas limits (the decode silently fails, so the dust never appears on
# phones) and spike hundreds of MB of RGBA during read-back.
PNG_ATLAS_MAX_SIDE_PX = 2048
# Phones render this synchronous low-resolution cube instead of decoding
# the atlas, which avoids iOS Safari's canvas/async decode failure modes.
MOBILE_LOWRES_MAX_DIM = 128


def restore_raw_volume_payloads(scene: dict) -> int:
    """Put large volume payloads back on the synchronous raw-uint8 path.

    Reverts the PNG-slab encoding: the browser decodes ``data_b64`` with a
    plain base64 read at startup instead of pulling 26 atlas images through
    a canvas, which is what the desktop viewer used before the slab change.
    The small mobile low-res cube is kept, since phones decode that (also
    synchronously) rather than the full-resolution payload.
    """
    try:
        from PIL import Image
    except ImportError:
        return 0
    import io

    def _attach_mobile_lowres(layer, volume, nx, ny, nz):
        max_dim = max(nx, ny, nz)
        if max_dim <= MOBILE_LOWRES_MAX_DIM:
            return
        stride = int(np.ceil(max_dim / MOBILE_LOWRES_MAX_DIM))
        small = volume[::stride, ::stride, ::stride]
        layer["data_lowres_b64"] = base64.b64encode(
            np.ascontiguousarray(small, dtype=np.uint8).tobytes()
        ).decode("ascii")
        layer["lowres_shape"] = {
            "x": int(small.shape[2]),
            "y": int(small.shape[1]),
            "z": int(small.shape[0]),
        }

    restored = 0
    for layer in (scene.get("volumes") or {}).get("layers") or []:
        shape = layer.get("shape") or {}
        nx = max(1, int(shape.get("x") or 1))
        ny = max(1, int(shape.get("y") or 1))
        nz = max(1, int(shape.get("z") or 1))
        slabs = layer.get("data_b64_slabs") or None
        if not slabs:
            raw_payload = layer.get("data_b64") or ""
            if (
                raw_payload
                and str(layer.get("data_encoding") or "uint8") == "uint8"
                and not layer.get("data_lowres_b64")
            ):
                raw = np.frombuffer(base64.b64decode(raw_payload), dtype=np.uint8)
                expected = nx * ny * nz
                if raw.size == expected:
                    _attach_mobile_lowres(
                        layer,
                        raw.reshape((nz, ny, nx)),
                        nx,
                        ny,
                        nz,
                    )
            continue
        tiles = layer.get("data_atlas_tiles") or {}
        slab_cols = max(1, int(tiles.get("x") or 1))
        per_slab = max(
            1,
            int(tiles.get("slices_per_slab") or (slab_cols * max(1, int(tiles.get("y") or 1)))),
        )
        volume = np.zeros((nz, ny, nx), dtype=np.uint8)
        for slab_index, slab_b64 in enumerate(slabs):
            slab = np.asarray(
                Image.open(io.BytesIO(base64.b64decode(slab_b64))).convert("L"),
                dtype=np.uint8,
            )
            for offset in range(per_slab):
                z_index = (slab_index * per_slab) + offset
                if z_index >= nz:
                    break
                row, col = divmod(offset, slab_cols)
                volume[z_index] = slab[
                    row * ny:(row + 1) * ny,
                    col * nx:(col + 1) * nx,
                ]
        if not layer.get("data_lowres_b64"):
            _attach_mobile_lowres(layer, volume, nx, ny, nz)
        layer["data_b64"] = base64.b64encode(
            np.ascontiguousarray(volume, dtype=np.uint8).tobytes()
        ).decode("ascii")
        layer["data_b64_slabs"] = None
        layer["data_atlas_tiles"] = None
        layer["data_encoding"] = "uint8"
        restored += 1
    return restored

def read_embedded_scene_spec(path: Path) -> dict:
    html = Path(path).read_text(encoding="utf-8")
    payload_id = None
    for payload_arg in re.findall(r"readOvizSceneSpecPayload\((.*?)\)", html):
        try:
            candidate = json.loads(payload_arg)
        except json.JSONDecodeError:
            continue
        if isinstance(candidate, str) and candidate:
            payload_id = candidate
            break
    if not payload_id:
        raise ValueError(f"Could not locate the compressed scene payload in {path}.")
    chunks = re.findall(
        (
            r"<script\b(?=[^>]*type=[\"']application/octet-stream[\"'])"
            r"(?=[^>]*data-oviz-payload-id=[\"']"
            + re.escape(payload_id)
            + r"[\"'])[^>]*data-oviz-payload-index=[\"'](\d+)[\"'][^>]*>(.*?)</script>"
        ),
        html,
        re.S,
    )
    if not chunks:
        raise ValueError(f"Could not locate payload chunks for {payload_id}.")
    encoded = "".join(
        re.sub(r"\s+", "", chunk)
        for _, chunk in sorted(chunks, key=lambda item: int(item[0]))
    )
    return json.loads(gzip.decompress(base64.b64decode(encoded)))


def load_ratzenboeck_scocen_sky_members(
    members_path: Path,
    metadata_path: Path,
) -> dict[str, list[dict]]:
    """Load every Ratzenboeck/SigMA member star for the Sky renderer."""
    members_path = Path(members_path).expanduser().resolve()
    metadata_path = Path(metadata_path).expanduser().resolve()
    for label, path in [("member catalog", members_path), ("metadata", metadata_path)]:
        if not path.exists():
            raise FileNotFoundError(f"Missing Ratzenboeck/SigMA {label}: {path}")

    with fits.open(metadata_path, memmap=True) as hdus:
        metadata = hdus[1].data
        label_to_name = {
            int(row["cluster_label"]): str(row["cluster_name"]).strip()
            for row in metadata
        }
    if len(label_to_name) != 37:
        raise RuntimeError(
            "Expected 37 unique Ratzenboeck/SigMA cluster labels, "
            f"found {len(label_to_name)}."
        )

    required_columns = {
        "source_id",
        "labels",
        "ra",
        "dec",
        "l",
        "b",
        "parallax",
        "pmra",
        "pmdec",
    }
    with fits.open(members_path, memmap=True) as hdus:
        data = hdus[1].data
        available_columns = set(data.dtype.names or ())
        missing_columns = sorted(required_columns - available_columns)
        if missing_columns:
            raise RuntimeError(
                "Ratzenboeck/SigMA member catalog is missing "
                f"{missing_columns!r}."
            )
        if len(data) != RATZENBOECK_SCOCEN_MEMBER_COUNT:
            raise RuntimeError(
                f"Expected {RATZENBOECK_SCOCEN_MEMBER_COUNT:,} "
                f"Ratzenboeck/SigMA stars, found {len(data):,}."
            )
        unknown_labels = sorted(set(map(int, data["labels"])) - set(label_to_name))
        if unknown_labels:
            raise RuntimeError(
                "Ratzenboeck/SigMA stars contain unmapped cluster labels: "
                f"{unknown_labels!r}."
            )

        grouped = {cluster_name: [] for cluster_name in label_to_name.values()}
        for row in data:
            cluster_name = label_to_name[int(row["labels"])]
            parallax_mas = float(row["parallax"])
            member = {
                "ra": round(float(row["ra"]), 6),
                "dec": round(float(row["dec"]), 6),
                "l": round(float(row["l"]), 6),
                "b": round(float(row["b"]), 6),
                "label": cluster_name,
                "source_id": str(int(row["source_id"])),
                "is_cluster_member": True,
            }
            pmra = float(row["pmra"])
            pmdec = float(row["pmdec"])
            if np.isfinite(pmra):
                member["pmra_masyr"] = round(pmra, 6)
            if np.isfinite(pmdec):
                member["pmdec_masyr"] = round(pmdec, 6)
            if np.isfinite(parallax_mas) and parallax_mas > 0.0:
                member["distance_pc"] = round(1000.0 / parallax_mas, 3)
            grouped[cluster_name].append(member)

    loaded_count = sum(map(len, grouped.values()))
    if loaded_count != RATZENBOECK_SCOCEN_MEMBER_COUNT:
        raise RuntimeError(
            "Ratzenboeck/SigMA member grouping changed the row count: "
            f"{loaded_count:,}."
        )
    return grouped


def drop_trace_by_name(scene: dict, trace_name: str) -> set[str]:
    """Remove a trace from frames, legend, and group visibility maps.

    Trace keys are explicit (``trace-N``) rather than positional, so dropping
    one leaves the remaining keys and their visibility entries intact.
    """
    removed_keys: set[str] = set()
    for frame in scene.get("frames") or []:
        kept_traces = []
        for trace in frame.get("traces") or []:
            if str(trace.get("name") or "") == trace_name:
                key = str(trace.get("key") or "")
                if key:
                    removed_keys.add(key)
                continue
            kept_traces.append(trace)
        frame["traces"] = kept_traces
    legend = scene.get("legend")
    if isinstance(legend, dict) and isinstance(legend.get("items"), list):
        for item in legend["items"]:
            if str(item.get("name") or "") == trace_name:
                key = str(item.get("key") or "")
                if key:
                    removed_keys.add(key)
        legend["items"] = [
            item for item in legend["items"]
            if str(item.get("name") or "") != trace_name
        ]
    for visibility in (scene.get("group_visibility") or {}).values():
        if isinstance(visibility, dict):
            for key in removed_keys:
                visibility.pop(key, None)
    legend_state = (scene.get("initial_state") or {}).get("legend_state")
    if isinstance(legend_state, dict):
        for key in removed_keys:
            legend_state.pop(key, None)
    return removed_keys


def build_state_only_scene(
    scene_spec: dict,
    cluster_velocities_path: Path | None = None,
    ratzenboeck_scocen_catalog_path: Path | None = None,
    ratzenboeck_scocen_metadata_path: Path | None = None,
    ratzenboeck_scocen_members_path: Path | None = None,
    *,
    drop_full_cluster_catalog: bool = True,
) -> dict:
    scene = deepcopy(scene_spec)
    scene["title"] = ""
    if drop_full_cluster_catalog:
        drop_trace_by_name(scene, FULL_CLUSTER_CATALOG_TRACE_NAME)
    states = scene.get("states")
    if isinstance(states, dict):
        states["items"] = []
        states["assets"] = {}
        states["revision"] = 0
        states["synchronized_revision"] = 0
        # This figure is published read-only: never read back or write a
        # local IndexedDB authoring draft. A stale draft from an earlier
        # build was being auto-applied over the exported scene (hiding
        # volumes and spiking memory), which is not something a visitor to
        # a public URL should ever hit.
        states["autosave_drafts"] = False
    scene["deck"] = {
        "schema_version": 2,
        "available": True,
        "enabled": False,
        "embedded": False,
        "revision": 0,
        "aspect_ratio": "16:9",
        "guides": {"smart": True, "grid": False, "grid_size": 20},
        "slides": [],
    }
    initial_volume_state = (
        scene.setdefault("initial_state", {}).setdefault("volume_state_by_key", {})
    )
    scene["initial_state"].setdefault("global_controls", {}).update({
        "size_points_by_stars_enabled": True,
        "fade_opacity_by_birth_time_enabled": True,
    })
    # Volumes stay visible on phones; the runtime stride-downsamples large
    # cubes before the GPU upload (OVIZ_MOBILE_MAX_VOLUME_DIM), which keeps
    # the Edenhofer map on screen without tripping iOS Safari memory kills.
    scene["initial_state"].pop("mobile_defer_volumes", None)
    restore_raw_volume_payloads(scene)
    existing_sky_layers = scene["initial_state"].get("sky_layers") or []
    existing_sky_by_survey = {}
    for layer in existing_sky_layers:
        if not isinstance(layer, dict):
            continue
        survey = str(layer.get("survey") or layer.get("key") or "")
        canonical_survey = SKY_BACKGROUND_SURVEY_ALIASES.get(survey, survey)
        normalized_layer = deepcopy(layer)
        if canonical_survey != survey:
            normalized_layer["survey"] = canonical_survey
            normalized_layer["key"] = canonical_survey
        existing_sky_by_survey[canonical_survey] = normalized_layer
    configured_sky_layers = []
    configured_sky_surveys = set()
    for default_layer in SKY_BACKGROUND_LAYERS:
        survey = str(default_layer["survey"])
        merged_layer = deepcopy(default_layer)
        merged_layer.update(existing_sky_by_survey.get(survey, {}))
        configured_sky_layers.append(merged_layer)
        configured_sky_surveys.add(survey)
    configured_sky_layers.extend(
        deepcopy(layer)
        for layer in existing_sky_layers
        if isinstance(layer, dict)
        and SKY_BACKGROUND_SURVEY_ALIASES.get(
            str(layer.get("survey") or layer.get("key") or ""),
            str(layer.get("survey") or layer.get("key") or ""),
        ) not in configured_sky_surveys
    )
    def _layer_survey(layer: dict) -> str:
        return str(layer.get("survey") or layer.get("key") or "")

    for layer in configured_sky_layers:
        overrides = SKY_BACKGROUND_STACK_DEFAULTS.get(_layer_survey(layer))
        if overrides:
            layer.update(deepcopy(overrides))
    configured_sky_layers = [
        layer
        for survey in SKY_BACKGROUND_STACK_TOP_ORDER
        for layer in configured_sky_layers
        if _layer_survey(layer) == survey
    ] + [
        layer
        for layer in configured_sky_layers
        if _layer_survey(layer) not in SKY_BACKGROUND_STACK_TOP_ORDER
    ]
    scene["initial_state"]["sky_layers"] = configured_sky_layers
    scene["initial_state"].setdefault("active_sky_layer_key", "P/Mellinger/color")
    scene.setdefault("sky_dome", {})["layer_groups"] = deepcopy(SKY_BACKGROUND_GROUPS)
    for layer in (scene.get("volumes") or {}).get("layers") or []:
        layer_name = str(layer.get("name") or "").lower()
        state_key = str(layer.get("state_key") or layer.get("key") or "")
        state = initial_volume_state.setdefault(state_key, {}) if state_key else None
        if "edenhofer" in layer_name:
            defaults = layer.setdefault("default_controls", {})
            defaults.update(EDENHOFER_GALACTIC_LIGHTING)
            if state is None:
                continue
            state.update({
                "visible": True,
                "showAllTimes": False,
                "lightingMode": "standard",
                "galacticCenter": [8122.0, 0.0, 0.0],
                "galacticLightIntensity": 1.35,
                "galacticAmbient": 0.22,
                "galacticExtinction": 2.4,
                "galacticScattering": 0.55,
                "galacticAnisotropy": 0.45,
                "galacticWarmth": 0.72,
            })
        elif "mccallum" in layer_name or "vergely" in layer_name:
            layer["only_at_t0"] = True
            layer["supports_show_all_times"] = False
            if "vergely" in layer_name:
                layer.setdefault("default_controls", {})["colormap"] = "Greys"
            if state is not None:
                state.update({"visible": False, "showAllTimes": False})
                if "vergely" in layer_name:
                    state["colormap"] = "Greys"
    if cluster_velocities_path is not None:
        velocity_path = Path(cluster_velocities_path).expanduser().resolve()
        scene.setdefault("provenance", {})["cluster_velocities"] = {
            "filename": velocity_path.name,
            "path": str(velocity_path),
            "release": "SDSS-V covariance audited",
        }
    if ratzenboeck_scocen_catalog_path is not None:
        catalog_path = Path(
            ratzenboeck_scocen_catalog_path
        ).expanduser().resolve()
        metadata_path = (
            Path(ratzenboeck_scocen_metadata_path).expanduser().resolve()
            if ratzenboeck_scocen_metadata_path is not None
            else None
        )
        provenance = {
            "catalog_filename": catalog_path.name,
            "catalog_path": str(catalog_path),
            "metadata_filename": metadata_path.name if metadata_path else None,
            "metadata_path": str(metadata_path) if metadata_path else None,
            "n_clusters": 37,
            "membership_reference": "Ratzenboeck et al. 2023, A&A 677, A59",
            "age_reference": "Ratzenboeck et al. 2023, A&A 678, A71",
        }
        if ratzenboeck_scocen_members_path is not None:
            if metadata_path is None:
                raise ValueError(
                    "Ratzenboeck/SigMA members require the cluster metadata path."
                )
            members_path = Path(
                ratzenboeck_scocen_members_path
            ).expanduser().resolve()
            sigma_members = load_ratzenboeck_scocen_sky_members(
                members_path,
                metadata_path,
            )
            sky_panel = scene.setdefault("sky_panel", {})
            sky_panel["show_cluster_members_in_sky"] = True
            sky_panel.setdefault("members_by_cluster", {}).update(sigma_members)
            provenance.update({
                "member_catalog_filename": members_path.name,
                "member_catalog_path": str(members_path),
                "n_member_stars": RATZENBOECK_SCOCEN_MEMBER_COUNT,
                "member_identifier": "Gaia source_id",
            })
        scene.setdefault("provenance", {})["ratzenboeck_scocen"] = provenance

        present_frame = min(
            scene.get("frames") or [],
            key=lambda frame: abs(float(frame.get("time", 0.0))),
            default=None,
        )
        if present_frame is not None:
            sigma_trace = next(
                (
                    trace
                    for trace in present_frame.get("traces") or []
                    if trace.get("name") == RATZENBOECK_SCOCEN_TRACE_NAME
                ),
                None,
            )
            if sigma_trace is None:
                raise RuntimeError(
                    "The July 25 scene is missing the Ratzenboeck/SigMA trace."
                )
            sigma_trace_key = str(sigma_trace.get("key") or "")
            if not sigma_trace_key:
                raise RuntimeError("The Ratzenboeck/SigMA trace has no stable key.")

            sigma_track_name = f"{RATZENBOECK_SCOCEN_TRACE_NAME} Track"
            sigma_track_keys = {
                str(trace.get("key"))
                for frame in scene.get("frames") or []
                for trace in frame.get("traces") or []
                if trace.get("name") == sigma_track_name and trace.get("key")
            }
            for frame in scene.get("frames") or []:
                frame["traces"] = [
                    trace
                    for trace in frame.get("traces") or []
                    if trace.get("name") != sigma_track_name
                ]

            legend = scene.get("legend")
            if isinstance(legend, dict) and isinstance(legend.get("items"), list):
                legend["items"] = [
                    item
                    for item in legend["items"]
                    if item.get("name") != sigma_track_name
                    and str(item.get("key") or "") not in sigma_track_keys
                ]

            scene["group_order"] = [
                group_name
                for group_name in scene.get("group_order") or []
                if group_name != RATZENBOECK_SCOCEN_TRACE_NAME
            ]
            group_visibility = scene.setdefault("group_visibility", {})
            group_visibility.pop(RATZENBOECK_SCOCEN_TRACE_NAME, None)
            group_visibility.setdefault("All", {})
            for group_name, visibility in group_visibility.items():
                visibility[sigma_trace_key] = group_name == "All"
                for sigma_track_key in sigma_track_keys:
                    visibility.pop(sigma_track_key, None)
    return scene


def build_figure(
    source_html: Path,
    output_html: Path,
    cluster_velocities_path: Path | None = None,
    ratzenboeck_scocen_catalog_path: Path | None = None,
    ratzenboeck_scocen_metadata_path: Path | None = None,
    ratzenboeck_scocen_members_path: Path | None = None,
) -> Path:
    source_html = Path(source_html).expanduser().resolve()
    output_html = Path(output_html).expanduser().resolve()
    if source_html == output_html:
        raise ValueError("The July 25 figure must be a sibling copy, not the source artifact.")
    scene = build_state_only_scene(
        read_embedded_scene_spec(source_html),
        cluster_velocities_path=cluster_velocities_path,
        ratzenboeck_scocen_catalog_path=ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path=ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path=ratzenboeck_scocen_members_path,
    )
    html = ThreeJSFigure(scene, compress_scene_spec=True).to_html(compress_scene_spec=True)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(html, encoding="utf-8")
    return output_html


def resolve_source_input_paths(
    cluster_velocities_path: Path,
    ratzenboeck_scocen_catalog_path: Path,
    ratzenboeck_scocen_metadata_path: Path,
    ratzenboeck_scocen_members_path: Path,
) -> tuple[Path, Path, Path, Path]:
    cluster_velocities_path = Path(cluster_velocities_path).expanduser().resolve()
    if not cluster_velocities_path.exists():
        raise FileNotFoundError(
            f"Missing cluster velocity catalog: {cluster_velocities_path}"
        )
    ratzenboeck_scocen_catalog_path = Path(
        ratzenboeck_scocen_catalog_path
    ).expanduser().resolve()
    ratzenboeck_scocen_metadata_path = Path(
        ratzenboeck_scocen_metadata_path
    ).expanduser().resolve()
    ratzenboeck_scocen_members_path = Path(
        ratzenboeck_scocen_members_path
    ).expanduser().resolve()
    for label, path in [
        ("cluster catalog", ratzenboeck_scocen_catalog_path),
        ("cluster metadata", ratzenboeck_scocen_metadata_path),
        ("member catalog", ratzenboeck_scocen_members_path),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"Missing Ratzenboeck/SigMA {label}: {path}")
    return (
        cluster_velocities_path,
        ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path,
    )


def source_main_figure_kwargs(
    cluster_velocities_path: Path,
    ratzenboeck_scocen_catalog_path: Path,
    ratzenboeck_scocen_metadata_path: Path,
    *,
    lookback_myr: int = JULY25_LOOKBACK_MYR,
) -> dict:
    """Shared notebook-runner settings for the July source scene."""
    return {
        "mobile_mode": False,
        "compact_payload": True,
        "mobile_safe_mode": False,
        "lookback_myr": lookback_myr,
        "chronos_results_path": JULY4_CHRONOS_RESULTS_PATH,
        "chronos_model": "parsec",
        "include_spiral_arms": False,
        "jun6_catalog": True,
        "cluster_velocities_path": cluster_velocities_path,
        "include_background_cluster_trace": False,
        "cluster_members_file": CLUSTER_MEMBERS_PATH,
        "show_cluster_members_in_sky": True,
        "include_mccallum_halpha_volume": True,
        "include_vergely_dust_volume": True,
        "ratzenboeck_scocen_catalog_path": ratzenboeck_scocen_catalog_path,
        "ratzenboeck_scocen_metadata_path": ratzenboeck_scocen_metadata_path,
    }


def build_velocity_source_figure(
    cluster_velocities_path: Path,
    output_html: Path,
    ratzenboeck_scocen_catalog_path: Path = DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
    ratzenboeck_scocen_metadata_path: Path = DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    ratzenboeck_scocen_members_path: Path = DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
) -> Path:
    (
        cluster_velocities_path,
        ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path,
    ) = resolve_source_input_paths(
        cluster_velocities_path,
        ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path,
    )
    return build_source_main_figure(
        output_html=output_html,
        website_output_html=None,
        **source_main_figure_kwargs(
            cluster_velocities_path,
            ratzenboeck_scocen_catalog_path,
            ratzenboeck_scocen_metadata_path,
        ),
    )


def build_velocity_source_scene(
    cluster_velocities_path: Path,
    ratzenboeck_scocen_catalog_path: Path = DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
    ratzenboeck_scocen_metadata_path: Path = DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    ratzenboeck_scocen_members_path: Path = DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
    *,
    lookback_myr: int = JULY25_LOOKBACK_MYR,
) -> dict:
    """Build the shared source scene in-process, no HTML written or parsed."""
    (
        cluster_velocities_path,
        ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path,
    ) = resolve_source_input_paths(
        cluster_velocities_path,
        ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path,
    )
    return build_source_main_figure_scene(
        **source_main_figure_kwargs(
            cluster_velocities_path,
            ratzenboeck_scocen_catalog_path,
            ratzenboeck_scocen_metadata_path,
            lookback_myr=lookback_myr,
        ),
    )


def build_scene_from_velocity_catalog(
    cluster_velocities_path: Path,
    ratzenboeck_scocen_catalog_path: Path = DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
    ratzenboeck_scocen_metadata_path: Path = DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    ratzenboeck_scocen_members_path: Path = DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
) -> dict:
    """Return the finished July 25 presentation scene as a Python dict."""
    source_scene = build_velocity_source_scene(
        cluster_velocities_path,
        ratzenboeck_scocen_catalog_path=ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path=ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path=ratzenboeck_scocen_members_path,
    )
    return build_state_only_scene(
        source_scene,
        cluster_velocities_path=cluster_velocities_path,
        ratzenboeck_scocen_catalog_path=ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path=ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path=ratzenboeck_scocen_members_path,
    )


def build_figure_from_velocity_catalog(
    cluster_velocities_path: Path,
    output_html: Path,
    ratzenboeck_scocen_catalog_path: Path = DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
    ratzenboeck_scocen_metadata_path: Path = DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    ratzenboeck_scocen_members_path: Path = DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
) -> Path:
    output_html = Path(output_html).expanduser().resolve()
    scene = build_scene_from_velocity_catalog(
        cluster_velocities_path,
        ratzenboeck_scocen_catalog_path=ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path=ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path=ratzenboeck_scocen_members_path,
    )
    html = ThreeJSFigure(scene, compress_scene_spec=True).to_html(
        compress_scene_spec=True
    )
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(html, encoding="utf-8")
    return output_html


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-html",
        type=Path,
        default=None,
        help=(
            "Optional prebuilt source scene. If omitted, rebuild the source "
            "from the selected cluster velocity catalog."
        ),
    )
    parser.add_argument(
        "--cluster-velocities-path",
        type=Path,
        default=DEFAULT_CLUSTER_VELOCITIES_PATH,
        help="Cluster velocity CSV used when rebuilding the source scene.",
    )
    parser.add_argument(
        "--ratzenboeck-scocen-catalog-path",
        type=Path,
        default=DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
        help="Formatted 37-cluster Ratzenboeck/SigMA position/velocity catalog.",
    )
    parser.add_argument(
        "--ratzenboeck-scocen-metadata-path",
        type=Path,
        default=DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
        help="Ratzenboeck/SigMA cluster identifiers, regions, and member counts.",
    )
    parser.add_argument(
        "--ratzenboeck-scocen-members-path",
        type=Path,
        default=DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
        help="Ratzenboeck/SigMA 13,103-star membership catalog.",
    )
    parser.add_argument("--output-html", type=Path, default=DEFAULT_OUTPUT_HTML)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.source_html is not None:
        output = build_figure(
            args.source_html,
            args.output_html,
            cluster_velocities_path=args.cluster_velocities_path,
            ratzenboeck_scocen_catalog_path=args.ratzenboeck_scocen_catalog_path,
            ratzenboeck_scocen_metadata_path=args.ratzenboeck_scocen_metadata_path,
            ratzenboeck_scocen_members_path=args.ratzenboeck_scocen_members_path,
        )
    else:
        output = build_figure_from_velocity_catalog(
            args.cluster_velocities_path,
            args.output_html,
            ratzenboeck_scocen_catalog_path=args.ratzenboeck_scocen_catalog_path,
            ratzenboeck_scocen_metadata_path=args.ratzenboeck_scocen_metadata_path,
            ratzenboeck_scocen_members_path=args.ratzenboeck_scocen_members_path,
        )
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
