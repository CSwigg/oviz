import shutil
import subprocess
import textwrap
import zipfile
from pathlib import Path

import pytest

from oviz.threejs_runtime_ar import THREEJS_AR_RUNTIME_JS


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not available")
def test_oviz_usdz_exporter_writes_transforms_and_aligned_store_zip(tmp_path):
    output_path = tmp_path / "oviz-test.usdz"
    script = f"""
    const fs = require("fs");
    {THREEJS_AR_RUNTIME_JS}

    function attribute(values, itemSize) {{
      return {{
        count: values.length / itemSize,
        getX(index) {{ return values[index * itemSize]; }},
        getY(index) {{ return values[(index * itemSize) + 1]; }},
        getZ(index) {{ return values[(index * itemSize) + 2]; }},
      }};
    }}

    const geometry = {{
      id: 3,
      userData: {{}},
      attributes: {{
        position: attribute([0, 0, 0, 1, 0, 0, 0, 1, 0], 3),
        normal: attribute([0, 0, 1, 0, 0, 1, 0, 0, 1], 3),
        uv: attribute([0, 0, 1, 0, 0, 1], 2),
      }},
      index: {{ count: 3, getX(index) {{ return index; }} }},
    }};
    const material = {{
      id: 4,
      color: {{ r: 0.2, g: 0.6, b: 0.9 }},
      emissive: {{ r: 0.02, g: 0.03, b: 0.04 }},
      emissiveIntensity: 1,
      opacity: 1,
      roughness: 0.5,
      metalness: 0,
      map: null,
    }};
    const mesh = {{
      id: 5,
      isMesh: true,
      geometry,
      material,
      matrixWorld: {{
        elements: [
          1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0.25, 0.5, -0.3, 1,
        ],
      }},
    }};
    const scene = {{
      updateMatrixWorld(force) {{ this.updated = force; }},
      traverseVisible(callback) {{ callback(mesh); }},
    }};

    (async () => {{
      const payload = await new OvizUSDZExporter().parse(scene);
      if (!scene.updated) throw new Error("world matrices were not updated");
      fs.writeFileSync({str(output_path)!r}, Buffer.from(new Uint8Array(payload)));
    }})().catch((error) => {{
      console.error(error && error.stack ? error.stack : error);
      process.exit(1);
    }});
    """
    subprocess.run(
        ["node"],
        input=textwrap.dedent(script),
        text=True,
        capture_output=True,
        check=True,
    )

    assert output_path.stat().st_size > 1024
    with zipfile.ZipFile(output_path) as archive:
        names = archive.namelist()
        assert names[0] == "model.usda"
        assert "geometries/Geometry_3.usda" in names
        assert all(info.compress_type == zipfile.ZIP_STORED for info in archive.infolist())
        model = archive.read("model.usda").decode("utf-8")
        assert 'defaultPrim = "Root"' in model
        assert 'prepend apiSchemas = ["MaterialBindingAPI"]' in model
        assert "0.2500000" in model
        assert "0.5000000" in model
        assert "-0.3000000" in model
        with output_path.open("rb") as handle:
            for info in archive.infolist():
                handle.seek(info.header_offset + 26)
                name_length = int.from_bytes(handle.read(2), "little")
                extra_length = int.from_bytes(handle.read(2), "little")
                data_offset = info.header_offset + 30 + name_length + extra_length
                assert data_offset % 64 == 0


@pytest.mark.skipif(shutil.which("/usr/bin/usdchecker") is None, reason="usdchecker is not available")
@pytest.mark.skipif(shutil.which("node") is None, reason="node is not available")
def test_oviz_usdz_exporter_passes_apple_arkit_checker(tmp_path):
    output_path = tmp_path / "oviz-test.usdz"
    script = f"""
    const fs = require("fs");
    {THREEJS_AR_RUNTIME_JS}
    const position = {{ count: 3, getX(i) {{ return [0, 1, 0][i]; }}, getY(i) {{ return [0, 0, 1][i]; }}, getZ() {{ return 0; }} }};
    const normal = {{ count: 3, getX() {{ return 0; }}, getY() {{ return 0; }}, getZ() {{ return 1; }} }};
    const mesh = {{
      id: 1,
      isMesh: true,
      geometry: {{ id: 1, userData: {{}}, attributes: {{ position, normal }}, index: {{ count: 3, getX(i) {{ return i; }} }} }},
      material: {{ id: 1, color: {{ r: 1, g: 0, b: 0 }}, emissive: {{ r: 0, g: 0, b: 0 }}, opacity: 1, roughness: 0.5, metalness: 0, map: null }},
      matrixWorld: {{ elements: [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0.1, 0.2, 0.3, 1] }},
    }};
    const scene = {{ updateMatrixWorld() {{}}, traverseVisible(callback) {{ callback(mesh); }} }};
    (async () => {{
      const payload = await new OvizUSDZExporter().parse(scene);
      fs.writeFileSync({str(output_path)!r}, Buffer.from(new Uint8Array(payload)));
    }})().catch((error) => {{ console.error(error); process.exit(1); }});
    """
    subprocess.run(["node"], input=textwrap.dedent(script), text=True, capture_output=True, check=True)

    result = subprocess.run(
        ["/usr/bin/usdchecker", "--arkit", str(output_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
