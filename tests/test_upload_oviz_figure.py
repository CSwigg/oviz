import importlib.util
import sys
from pathlib import Path

import pytest


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "upload_oviz_figure.py"


def _load_upload_module():
    spec = importlib.util.spec_from_file_location("upload_oviz_figure", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_html(path: Path, body: str = "<!doctype html><html></html>") -> Path:
    path.write_text(body, encoding="utf-8")
    return path


def test_build_upload_plan_defaults_to_website_oviz_figures(tmp_path):
    module = _load_upload_module()
    source = _write_html(tmp_path / "figure.html")
    website_dir = tmp_path / "cam_website"

    plan = module.build_upload_plan(
        source,
        website_dir=website_dir,
        remote_url="https://github.com/CSwigg/cam_website.git",
    )

    assert plan.source == source.resolve()
    assert plan.destination == (website_dir / "oviz_figures" / "figure.html").resolve()
    assert plan.git_path == "oviz_figures/figure.html"
    assert plan.commit_message == "Upload Oviz figure figure.html"
    assert plan.push is True
    assert plan.public_url == "https://cswigg.github.io/cam_website/oviz_figures/figure.html"


def test_build_upload_plan_derives_project_pages_url_from_ssh_remote(tmp_path):
    module = _load_upload_module()
    source = _write_html(tmp_path / "figure with space.html")

    plan = module.build_upload_plan(
        source,
        website_dir=tmp_path / "site",
        remote_url="git@github.com:CSwigg/cam_website.git",
    )

    assert plan.public_url == "https://cswigg.github.io/cam_website/oviz_figures/figure%20with%20space.html"


def test_build_upload_plan_derives_user_pages_url(tmp_path):
    module = _load_upload_module()
    source = _write_html(tmp_path / "figure.html")

    plan = module.build_upload_plan(
        source,
        website_dir=tmp_path / "site",
        remote_url="https://github.com/CSwigg/cswigg.github.io.git",
    )

    assert plan.public_url == "https://cswigg.github.io/oviz_figures/figure.html"


def test_build_upload_plan_rejects_oversized_file(tmp_path):
    module = _load_upload_module()
    source = _write_html(tmp_path / "large.html", "x" * 2048)

    with pytest.raises(ValueError, match="above the configured"):
        module.build_upload_plan(source, website_dir=tmp_path / "site", max_size_mb=0.001)


def test_build_upload_plan_uses_cautious_default_size_limit(tmp_path):
    module = _load_upload_module()
    source = _write_html(tmp_path / "large.html")
    with source.open("ab") as handle:
        handle.truncate((module.DEFAULT_MAX_SIZE_MB * 1024 * 1024) + 1)

    with pytest.raises(ValueError, match="above the configured 25.0 MiB limit"):
        module.build_upload_plan(source, website_dir=tmp_path / "site")


def test_build_upload_plan_rejects_path_like_output_name(tmp_path):
    module = _load_upload_module()
    source = _write_html(tmp_path / "figure.html")

    with pytest.raises(ValueError, match="plain file name"):
        module.build_upload_plan(source, website_dir=tmp_path / "site", output_name="../figure.html")


def test_upload_oviz_figure_copies_and_runs_scoped_git_commands(tmp_path, monkeypatch):
    module = _load_upload_module()
    source = _write_html(tmp_path / "figure.html", "<html>oviz</html>")
    website_dir = tmp_path / "cam_website"
    plan = module.build_upload_plan(
        source,
        website_dir=website_dir,
        output_name="chronos.html",
        commit_message="Upload test figure",
        push=False,
    )
    commands = []

    def fake_run(command, *, cwd, check):
        commands.append((command, Path(cwd), check))

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    module.upload_oviz_figure(plan)

    assert plan.destination.read_text(encoding="utf-8") == "<html>oviz</html>"
    assert commands == [
        (["git", "add", "--", "oviz_figures/chronos.html"], website_dir.resolve(), True),
        (
            ["git", "commit", "-m", "Upload test figure", "--", "oviz_figures/chronos.html"],
            website_dir.resolve(),
            True,
        ),
    ]


def test_upload_oviz_figure_dry_run_does_not_copy_or_call_subprocess(tmp_path, monkeypatch):
    module = _load_upload_module()
    source = _write_html(tmp_path / "figure.html")
    website_dir = tmp_path / "cam_website"
    plan = module.build_upload_plan(source, website_dir=website_dir, dry_run=True)
    calls = []

    def fake_run(*args, **kwargs):
        calls.append((args, kwargs))

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    module.upload_oviz_figure(plan)

    assert not plan.destination.exists()
    assert calls == []
