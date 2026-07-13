from __future__ import annotations

import json
import subprocess

from oviz.threejs_figure import ThreeJSFigure
from oviz.threejs_runtime_states import THREEJS_STATE_RUNTIME_JS


def _run_mask_component_runtime(body: str) -> dict:
    start = THREEJS_STATE_RUNTIME_JS.index(
        "const OVIZ_MAX_VOLUME_MASK_SOURCE_COMPONENTS = 4;"
    )
    end = THREEJS_STATE_RUNTIME_JS.index("function ovizStatesUuid", start)
    helpers = THREEJS_STATE_RUNTIME_JS[start:end]
    script = f"""
function clampRange(value, low, high) {{
  return Math.min(high, Math.max(low, Number(value)));
}}
{helpers}
{body}
"""
    completed = subprocess.run(
        ["node", "-e", script],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(completed.stdout)


def test_repeated_retarget_freezes_four_weighted_runtime_masks() -> None:
    payload = _run_mask_component_runtime(
        """
const masks = ["A", "B", "C", "D", "E"].map((name) => ({ name }));
let source = ovizFreezeVolumeMaskSourceComponents(null, masks[0]);
source = ovizFreezeVolumeMaskSourceComponents({
  sourceMaskComponents: source,
  toMask: masks[1],
  progress: 0.25,
});
source = ovizFreezeVolumeMaskSourceComponents({
  sourceMaskComponents: source,
  toMask: masks[2],
  progress: 0.5,
});
source = ovizFreezeVolumeMaskSourceComponents({
  sourceMaskComponents: source,
  toMask: masks[3],
  progress: 0.2,
});
const nextTransition = {
  sourceMaskComponents: source,
  toMask: masks[4],
  progress: 0.0,
};
console.log(JSON.stringify({
  count: source.length,
  weights: Object.fromEntries(source.map((component) => [component.mask.name, component.weight])),
  nextSource: ovizVolumeMaskTransitionSourceComponents(nextTransition).map(
    (component) => component.mask.name
  ),
  sum: source.reduce((total, component) => total + component.weight, 0),
}));
"""
    )

    assert payload["count"] == 4
    assert payload["nextSource"] == ["A", "B", "C", "D"]
    expected = {"A": 0.3, "B": 0.1, "C": 0.4, "D": 0.2}
    assert payload["weights"].keys() == expected.keys()
    for key, value in expected.items():
        assert abs(payload["weights"][key] - value) < 1e-12
    assert abs(payload["sum"] - 1.0) < 1e-12


def test_unfiltered_endpoint_is_retained_as_a_weighted_component() -> None:
    payload = _run_mask_component_runtime(
        """
const mask = { name: "selected" };
const source = ovizFreezeVolumeMaskSourceComponents({
  sourceMaskComponents: [{ mask, weight: 1.0 }],
  toMask: null,
  progress: 0.4,
});
console.log(JSON.stringify(source.map((component) => ({
  name: component.mask ? component.mask.name : "all",
  weight: component.weight,
}))));
"""
    )

    assert payload == [
        {"name": "selected", "weight": 0.6},
        {"name": "all", "weight": 0.4},
    ]


def test_volume_shader_and_mask_lifetime_cover_all_retained_components() -> None:
    html = ThreeJSFigure(
        {"width": 320, "height": 240, "frames": [], "initial_state": {}}
    ).to_html(compress_scene_spec=False)

    assert "selectionSourceWeights" in html
    assert "selectionSourceTertiaryMaskTexture" in html
    assert "selectionSourceQuaternaryMaskTexture" in html
    assert "sourceSelectionWeight = dot(" in html
    assert "selectionTransition.sourceMaskComponents" in html
    assert "ovizStateSelectionTransition || ovizHeldSelectionTransition" in html
    assert "transition.sourceMaskComponents" in html
    assert "preservedSelectionMasks" in html
