import importlib.util
import sys
import zipfile
from argparse import Namespace
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("main_figure_chronos_july4.py")


def _load_runner_module():
    tests_dir = str(SCRIPT_PATH.parent)
    if tests_dir not in sys.path:
        sys.path.insert(0, tests_dir)
    spec = importlib.util.spec_from_file_location("main_figure_chronos_july4_runner", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_write_zip_copy_contains_only_html_artifact(tmp_path):
    module = _load_runner_module()
    html_path = tmp_path / "figure.html"
    zip_path = tmp_path / "figure.zip"
    html_path.write_text("<!doctype html><html>oviz</html>", encoding="utf-8")

    result = module.write_zip_copy(html_path, zip_path)

    assert result == zip_path.resolve()
    with zipfile.ZipFile(zip_path) as archive:
        assert archive.namelist() == ["figure.html"]
        assert archive.read("figure.html") == html_path.read_bytes()


def test_main_writes_default_zip_for_rendered_html(tmp_path, monkeypatch):
    module = _load_runner_module()
    html_path = tmp_path / "chronos.html"

    def fake_run_main_figure(**_kwargs):
        html_path.write_text("<!doctype html><html>oviz</html>", encoding="utf-8")
        return html_path

    monkeypatch.setattr(module, "run_main_figure", fake_run_main_figure)
    monkeypatch.setattr(
        module,
        "parse_args",
        lambda: Namespace(
            output_html=html_path,
            include_spiral_arms=False,
            full_resolution=False,
            force_mobile_mode=False,
            full_payload=False,
            zip_output=None,
            no_zip=False,
        ),
    )

    module.main()

    with zipfile.ZipFile(html_path.with_suffix(".zip")) as archive:
        assert archive.namelist() == ["chronos.html"]
        assert archive.read("chronos.html") == html_path.read_bytes()
