#!/usr/bin/env python3
"""Copy an Oviz HTML export into the website repo and publish it."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from oviz.threejs_runtime_ar import THREEJS_AR_QUICKLOOK_SERVICE_WORKER_JS


DEFAULT_SOURCE_HTML = REPO_ROOT / "tests" / "main_figure_chronos_july4.html"
DEFAULT_WEBSITE_DIR = Path.home() / "Desktop" / "astro_research" / "cam_website"
DEFAULT_TARGET_DIR = DEFAULT_WEBSITE_DIR / "oviz_figures"
DEFAULT_MAX_SIZE_MB = 25
DEFAULT_VERIFY_ATTEMPTS = 6
DEFAULT_VERIFY_DELAY_SECONDS = 5.0
DEFAULT_VERIFY_TIMEOUT_SECONDS = 10.0
AR_QUICKLOOK_SERVICE_WORKER_NAME = "oviz-ar-quicklook-sw.js"


@dataclass(frozen=True)
class UploadPlan:
    source: Path
    destination: Path
    website_dir: Path
    git_path: str
    companion_destination: Path | None
    companion_git_path: str | None
    size_bytes: int
    max_size_bytes: int
    commit_message: str
    push: bool
    dry_run: bool
    public_url: str | None
    verify_url: bool
    verify_attempts: int
    verify_delay_seconds: float
    verify_timeout_seconds: float


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


def _origin_remote_url(website_dir: Path) -> str | None:
    if not website_dir.exists():
        return None
    result = subprocess.run(
        ["git", "remote", "get-url", "origin"],
        cwd=website_dir,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        return None
    remote_url = result.stdout.strip()
    return remote_url or None


def _github_pages_base_url(remote_url: str | None) -> str | None:
    if not remote_url:
        return None
    remote_url = remote_url.strip()
    owner_repo = None
    if remote_url.startswith("git@github.com:"):
        owner_repo = remote_url.removeprefix("git@github.com:")
    elif remote_url.startswith("https://github.com/"):
        owner_repo = remote_url.removeprefix("https://github.com/")
    elif remote_url.startswith("http://github.com/"):
        owner_repo = remote_url.removeprefix("http://github.com/")
    if owner_repo is None:
        return None
    owner_repo = owner_repo.removesuffix(".git").strip("/")
    parts = owner_repo.split("/")
    if len(parts) != 2 or not all(parts):
        return None
    owner, repo = parts
    host_owner = owner.lower()
    if repo.lower() == f"{host_owner}.github.io":
        return f"https://{host_owner}.github.io/"
    return f"https://{host_owner}.github.io/{quote(repo.strip('/'))}/"


def _public_url(base_url: str | None, git_path: str) -> str | None:
    if not base_url:
        return None
    normalized_base = base_url.rstrip("/") + "/"
    return normalized_base + "/".join(quote(part) for part in git_path.split("/"))


def build_upload_plan(
    source_html: Path,
    *,
    website_dir: Path = DEFAULT_WEBSITE_DIR,
    target_dir: Path | None = None,
    output_name: str | None = None,
    max_size_mb: float = DEFAULT_MAX_SIZE_MB,
    commit_message: str | None = None,
    public_base_url: str | None = None,
    remote_url: str | None = None,
    push: bool = True,
    dry_run: bool = False,
    verify_url: bool = False,
    verify_attempts: int = DEFAULT_VERIFY_ATTEMPTS,
    verify_delay_seconds: float = DEFAULT_VERIFY_DELAY_SECONDS,
    verify_timeout_seconds: float = DEFAULT_VERIFY_TIMEOUT_SECONDS,
) -> UploadPlan:
    source = source_html.expanduser().resolve()
    website_dir = website_dir.expanduser().resolve()
    target_dir = (target_dir or website_dir / "oviz_figures").expanduser().resolve()
    max_size_bytes = _max_size_bytes(max_size_mb)
    size_bytes = _validate_source(source, max_size_bytes)
    destination = _destination_for(source, target_dir, output_name)
    git_path = _git_relative_path(destination, website_dir)
    includes_ar_quicklook = b"oviz-ar-quicklook-sw.js" in source.read_bytes()
    companion_destination = target_dir / AR_QUICKLOOK_SERVICE_WORKER_NAME if includes_ar_quicklook else None
    companion_git_path = (
        _git_relative_path(companion_destination, website_dir)
        if companion_destination is not None
        else None
    )
    message = commit_message or f"Upload Oviz figure {destination.name}"
    base_url = public_base_url or _github_pages_base_url(remote_url or _origin_remote_url(website_dir))
    public_url = _public_url(base_url, git_path)
    if verify_url and not public_url:
        raise ValueError("Cannot verify the public URL because no GitHub Pages URL could be derived.")
    return UploadPlan(
        source=source,
        destination=destination,
        website_dir=website_dir,
        git_path=git_path,
        companion_destination=companion_destination,
        companion_git_path=companion_git_path,
        size_bytes=size_bytes,
        max_size_bytes=max_size_bytes,
        commit_message=message,
        push=push,
        dry_run=dry_run,
        public_url=public_url,
        verify_url=bool(verify_url),
        verify_attempts=max(int(verify_attempts), 1),
        verify_delay_seconds=max(float(verify_delay_seconds), 0.0),
        verify_timeout_seconds=max(float(verify_timeout_seconds), 0.1),
    )


def _run_git(website_dir: Path, args: list[str], *, dry_run: bool) -> None:
    command = ["git", *args]
    if dry_run:
        print(f"DRY-RUN: cd {website_dir} && {' '.join(command)}")
        return
    subprocess.run(command, cwd=website_dir, check=True)


def verify_public_url(
    url: str,
    *,
    attempts: int = DEFAULT_VERIFY_ATTEMPTS,
    delay_seconds: float = DEFAULT_VERIFY_DELAY_SECONDS,
    timeout_seconds: float = DEFAULT_VERIFY_TIMEOUT_SECONDS,
) -> None:
    last_error: Exception | None = None
    attempts = max(int(attempts), 1)
    for attempt in range(1, attempts + 1):
        request = urllib.request.Request(
            url,
            method="GET",
            headers={
                "Range": "bytes=0-2047",
                "User-Agent": "oviz-upload-helper/1.0",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=max(float(timeout_seconds), 0.1)) as response:
                status = int(getattr(response, "status", 200))
                content_type = str(response.headers.get("Content-Type", "")).lower()
                if status in {200, 206} and ("html" in content_type or not content_type):
                    print(f"Verified URL: {url}")
                    return
                last_error = RuntimeError(f"unexpected HTTP {status} with Content-Type {content_type!r}")
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc
        if attempt < attempts and delay_seconds > 0:
            time.sleep(float(delay_seconds))
    raise RuntimeError(f"Could not verify public URL after {attempts} attempts: {url}") from last_error


def upload_oviz_figure(plan: UploadPlan) -> None:
    size_mb = plan.size_bytes / (1024 * 1024)
    print(f"Source: {plan.source}")
    print(f"Destination: {plan.destination}")
    print(f"Size: {size_mb:.1f} MiB")
    if plan.public_url:
        print(f"Expected URL: {plan.public_url}")
    if plan.companion_destination:
        print(f"Quick Look worker: {plan.companion_destination}")
    if plan.dry_run:
        print("DRY-RUN: not copying or publishing")
    else:
        plan.destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(plan.source, plan.destination)
        if plan.companion_destination:
            plan.companion_destination.write_text(
                THREEJS_AR_QUICKLOOK_SERVICE_WORKER_JS + "\n",
                encoding="utf-8",
            )

    git_paths = [plan.git_path]
    if plan.companion_git_path:
        git_paths.append(plan.companion_git_path)
    _run_git(plan.website_dir, ["add", "--", *git_paths], dry_run=plan.dry_run)
    _run_git(plan.website_dir, ["commit", "-m", plan.commit_message, "--", *git_paths], dry_run=plan.dry_run)
    if plan.push:
        _run_git(plan.website_dir, ["push"], dry_run=plan.dry_run)
    if plan.verify_url and plan.public_url:
        if plan.dry_run:
            print(f"DRY-RUN: would verify {plan.public_url}")
        else:
            verify_public_url(
                plan.public_url,
                attempts=plan.verify_attempts,
                delay_seconds=plan.verify_delay_seconds,
                timeout_seconds=plan.verify_timeout_seconds,
            )


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
        help="Refuse to upload files larger than this many MiB. Defaults to a cautious mobile-safe limit.",
    )
    parser.add_argument(
        "--message",
        help="Git commit message. Defaults to 'Upload Oviz figure <name>'.",
    )
    parser.add_argument(
        "--public-base-url",
        help="Public base URL for the website. Defaults to one derived from the GitHub origin remote.",
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
    parser.add_argument(
        "--verify-url",
        action="store_true",
        help="After publishing, verify that the derived public URL serves an HTML response.",
    )
    parser.add_argument(
        "--verify-attempts",
        type=int,
        default=DEFAULT_VERIFY_ATTEMPTS,
        help="Number of public URL verification attempts.",
    )
    parser.add_argument(
        "--verify-delay-seconds",
        type=float,
        default=DEFAULT_VERIFY_DELAY_SECONDS,
        help="Delay between public URL verification attempts.",
    )
    parser.add_argument(
        "--verify-timeout-seconds",
        type=float,
        default=DEFAULT_VERIFY_TIMEOUT_SECONDS,
        help="Per-request public URL verification timeout.",
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
        public_base_url=args.public_base_url,
        push=not bool(args.no_push),
        dry_run=bool(args.dry_run),
        verify_url=bool(args.verify_url),
        verify_attempts=int(args.verify_attempts),
        verify_delay_seconds=float(args.verify_delay_seconds),
        verify_timeout_seconds=float(args.verify_timeout_seconds),
    )
    upload_oviz_figure(plan)


if __name__ == "__main__":
    main()
