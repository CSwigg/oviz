#!/usr/bin/env python3
"""Copy an Oviz HTML export into the website repo and publish it."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_HTML = REPO_ROOT / "tests" / "main_figure_chronos_july4.html"
DEFAULT_WEBSITE_DIR = Path.home() / "Desktop" / "astro_research" / "cam_website"
DEFAULT_TARGET_DIR = DEFAULT_WEBSITE_DIR / "oviz_figures"
DEFAULT_MAX_SIZE_MB = 50


@dataclass(frozen=True)
class UploadPlan:
    source: Path
    destination: Path
    website_dir: Path
    git_path: str
    size_bytes: int
    max_size_bytes: int
    commit_message: str
    push: bool
    dry_run: bool


def _max_size_bytes(max_size_mb: float) -> int:
    return int(float(max_size_mb) * 1024 * 1024)


def _validate_source(source: Path, max_size_bytes: int) -> int:
    if not source.exists():
        raise FileNotFoundError(f"Missing Oviz HTML export: {source}")
    if not source.is_file():
        raise ValueError(f"Oviz export is not a file: {source}")
    if source.suffix.lower() != ".html":
        raise ValueError(f"Oviz export must be an .html file: {source}")
    size_bytes = source.stat().st_size
    if size_bytes <= 0:
        raise ValueError(f"Oviz export is empty: {source}")
    if size_bytes > max_size_bytes:
        size_mb = size_bytes / (1024 * 1024)
        limit_mb = max_size_bytes / (1024 * 1024)
        raise ValueError(
            f"Oviz export is {size_mb:.1f} MiB, above the configured {limit_mb:.1f} MiB limit: {source}"
        )
    return size_bytes


def _destination_for(source: Path, target_dir: Path, output_name: str | None) -> Path:
    name = output_name or source.name
    if Path(name).name != name:
        raise ValueError("--name must be a plain file name, not a path")
    if Path(name).suffix.lower() != ".html":
        raise ValueError("--name must end in .html")
    return target_dir / name


def _git_relative_path(path: Path, website_dir: Path) -> str:
    try:
        return path.resolve().relative_to(website_dir.resolve()).as_posix()
    except ValueError as exc:
        raise ValueError(f"Destination must live inside the website repo: {path}") from exc


def build_upload_plan(
    source_html: Path,
    *,
    website_dir: Path = DEFAULT_WEBSITE_DIR,
    target_dir: Path | None = None,
    output_name: str | None = None,
    max_size_mb: float = DEFAULT_MAX_SIZE_MB,
    commit_message: str | None = None,
    push: bool = True,
    dry_run: bool = False,
) -> UploadPlan:
    source = source_html.expanduser().resolve()
    website_dir = website_dir.expanduser().resolve()
    target_dir = (target_dir or website_dir / "oviz_figures").expanduser().resolve()
    max_size_bytes = _max_size_bytes(max_size_mb)
    size_bytes = _validate_source(source, max_size_bytes)
    destination = _destination_for(source, target_dir, output_name)
    git_path = _git_relative_path(destination, website_dir)
    message = commit_message or f"Upload Oviz figure {destination.name}"
    return UploadPlan(
        source=source,
        destination=destination,
        website_dir=website_dir,
        git_path=git_path,
        size_bytes=size_bytes,
        max_size_bytes=max_size_bytes,
        commit_message=message,
        push=push,
        dry_run=dry_run,
    )


def _run_git(website_dir: Path, args: list[str], *, dry_run: bool) -> None:
    command = ["git", *args]
    if dry_run:
        print(f"DRY-RUN: cd {website_dir} && {' '.join(command)}")
        return
    subprocess.run(command, cwd=website_dir, check=True)


def upload_oviz_figure(plan: UploadPlan) -> None:
    size_mb = plan.size_bytes / (1024 * 1024)
    print(f"Source: {plan.source}")
    print(f"Destination: {plan.destination}")
    print(f"Size: {size_mb:.1f} MiB")
    if plan.dry_run:
        print("DRY-RUN: not copying or publishing")
    else:
        plan.destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(plan.source, plan.destination)

    _run_git(plan.website_dir, ["add", "--", plan.git_path], dry_run=plan.dry_run)
    _run_git(plan.website_dir, ["commit", "-m", plan.commit_message, "--", plan.git_path], dry_run=plan.dry_run)
    if plan.push:
        _run_git(plan.website_dir, ["push"], dry_run=plan.dry_run)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source_html",
        nargs="?",
        type=Path,
        default=DEFAULT_SOURCE_HTML,
        help="Oviz standalone HTML file to upload.",
    )
    parser.add_argument(
        "--name",
        help="Destination file name inside oviz_figures. Defaults to the source file name.",
    )
    parser.add_argument(
        "--website-dir",
        type=Path,
        default=DEFAULT_WEBSITE_DIR,
        help="Local cam_website checkout.",
    )
    parser.add_argument(
        "--target-dir",
        type=Path,
        help="Destination directory. Defaults to WEBSITE_DIR/oviz_figures.",
    )
    parser.add_argument(
        "--max-size-mb",
        type=float,
        default=DEFAULT_MAX_SIZE_MB,
        help="Refuse to upload files larger than this many MiB.",
    )
    parser.add_argument(
        "--message",
        help="Git commit message. Defaults to 'Upload Oviz figure <name>'.",
    )
    parser.add_argument(
        "--no-push",
        action="store_true",
        help="Copy, add, and commit, but do not push.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the copy/git actions without changing files or running git.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plan = build_upload_plan(
        args.source_html,
        website_dir=args.website_dir,
        target_dir=args.target_dir,
        output_name=args.name,
        max_size_mb=args.max_size_mb,
        commit_message=args.message,
        push=not bool(args.no_push),
        dry_run=bool(args.dry_run),
    )
    upload_oviz_figure(plan)


if __name__ == "__main__":
    main()
