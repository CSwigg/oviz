"""make_plot(annotations=...) reaches the figure as the viewer's document."""

import numpy as np
import pandas as pd
import pytest

from oviz import Scene3D, Trace, TraceCollection
from oviz.annotations import normalize_annotations


def _scene():
    rng = np.random.default_rng(3)
    n = 12
    df = pd.DataFrame({
        "x": rng.normal(0, 100, n), "y": rng.normal(0, 100, n), "z": rng.normal(0, 20, n),
        "U": -10.0, "V": -15.0, "W": -7.0, "name": [f"c{i}" for i in range(n)], "age_myr": rng.uniform(5, 30, n),
    })
    return Scene3D(TraceCollection([Trace(df, data_name="Clusters", color="#ff5a5a")]), figure_theme="dark")


def test_annotations_compile_into_the_figure():
    fig = _scene().make_plot(time=np.arange(0, -6, -1.0), annotations=[
        {"kind": "shell", "center": [0, 0, 0], "radius": 150, "label": "Local Bubble", "present": True},
        {"kind": "arrow", "points": [[100, -300, 0], [50, -150, 20]], "color": "#ffd27a"},
        {"kind": "bubble", "center": [10, 0, 0], "radii": [30, 20, 10], "rot": [0, 0, 45]},
        {"kind": "text", "at": [300, 200, 0], "text": "Sco-Cen", "group": "Regions", "size": 16},
    ])
    doc = fig.bundle.manifest["annotations"]
    kinds = [i["kind"] for i in doc["items"]]
    assert kinds == ["sphere", "text", "curve", "sphere", "text"]
    shell, label = doc["items"][:2]
    assert shell["style"] == "shell" and shell["radii"] == [150.0] * 3 and shell["present"] is True
    assert label["text"] == "Local Bubble" and label["group"] == shell["group"]
    assert doc["items"][2]["arrow"] == "end"
    names = {g["id"]: g["name"] for g in doc["groups"]}
    assert names[shell["group"]] == "Local Bubble"
    assert names[doc["items"][4]["group"]] == "Regions"


def test_figures_without_annotations_are_unchanged():
    fig = _scene().make_plot(time=np.arange(0, -3, -1.0))
    assert "annotations" not in fig.bundle.manifest
    assert normalize_annotations(None) is None and normalize_annotations([]) is None


@pytest.mark.parametrize("bad", [
    {"kind": "blob"},
    {"kind": "text"},
    {"kind": "curve", "points": [[0, 0, 0]]},
    {"kind": "shell", "center": [0, 0, 0]},
    {"kind": "shell", "center": [0, 0], "radius": 3},
    {"kind": "bubble", "center": [0, 0, 0], "radii": [1, -1, 1]},
])
def test_bad_annotations_raise(bad):
    with pytest.raises(ValueError):
        normalize_annotations([bad])


def test_boxes_and_selections():
    doc = normalize_annotations([
        {"kind": "box", "center": [0, 0, 0], "size": [100, 50, 20], "rot": [0, 0, 30], "select": "isolate", "label": "Survey"},
        {"kind": "curve", "points": [[0, 0, 0], [100, 0, 0]], "select": True, "reach": 25},
        {"kind": "bubble", "center": [0, 0, 0], "radius": 10, "select": "highlight"},
    ])
    box, label, curve, bubble = doc["items"]
    assert box["kind"] == "box" and box["radii"] == [50.0, 25.0, 10.0] and box["select"] == "isolate"
    assert label["text"] == "Survey" and label["group"] == box["group"]
    assert curve["select"] == "highlight" and curve["reach"] == 25.0
    assert bubble["select"] == "highlight" and "reach" not in bubble


@pytest.mark.parametrize("bad", [
    {"kind": "box", "center": [0, 0, 0]},
    {"kind": "box", "center": [0, 0, 0], "size": [1, 0, 1]},
    {"kind": "bubble", "center": [0, 0, 0], "radius": 1, "select": "everything"},
    {"kind": "curve", "points": [[0, 0, 0], [1, 0, 0]], "select": "highlight", "reach": -1},
])
def test_bad_boxes_and_selections_raise(bad):
    with pytest.raises(ValueError):
        normalize_annotations([bad])
