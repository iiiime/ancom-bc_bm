#!/usr/bin/env python3
"""Download and verify benchmark datasets declared in data/manifest.yaml."""

from __future__ import annotations

import argparse
import ast
import hashlib
import os
import re
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = PROJECT_ROOT / "data" / "manifest.yaml"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data"
CHUNK_SIZE = 1024 * 1024
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class FileSpec:
    """One downloadable file from the dataset manifest."""

    dataset: str
    name: str
    sha256: str
    url: str


def _yaml_scalar(value: str) -> str:
    """Parse the quoted scalar style used by this repository's manifest."""
    value = value.strip()
    if value[:1] in {"'", '"'}:
        parsed = ast.literal_eval(value)
        if not isinstance(parsed, str):
            raise ValueError(f"Expected a string, got {value!r}")
        return parsed
    return value


def load_manifest(path: Path) -> dict[str, list[FileSpec]]:
    """Read the manifest without requiring a third-party YAML package."""
    path = Path(path)
    datasets: dict[str, list[FileSpec]] = {}
    in_datasets = False
    current_dataset: str | None = None
    current_file: dict[str, str] | None = None

    def finish_file() -> None:
        nonlocal current_file
        if current_file is None or current_dataset is None:
            return
        missing = {"name", "sha256", "url"} - current_file.keys()
        if missing:
            fields = ", ".join(sorted(missing))
            raise ValueError(f"Incomplete file entry in {current_dataset}: missing {fields}")
        name = current_file["name"]
        checksum = current_file["sha256"].lower()
        if Path(name).name != name:
            raise ValueError(f"Manifest filename must not contain a path: {name!r}")
        if not SHA256_PATTERN.fullmatch(checksum):
            raise ValueError(f"Invalid SHA-256 for {name}: {checksum!r}")
        datasets[current_dataset].append(
            FileSpec(current_dataset, name, checksum, current_file["url"])
        )
        current_file = None

    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        content = raw_line.split("#", 1)[0].rstrip()
        if not content.strip():
            continue
        indent = len(content) - len(content.lstrip())
        stripped = content.strip()

        if stripped == "datasets:":
            in_datasets = True
            continue
        if not in_datasets:
            continue

        if indent == 2 and stripped.endswith(":") and not stripped.startswith("-"):
            finish_file()
            current_dataset = stripped[:-1].strip()
            if not current_dataset or current_dataset in datasets:
                raise ValueError(f"Invalid or duplicate dataset at line {line_number}")
            datasets[current_dataset] = []
            continue

        if stripped.startswith("- name:"):
            if current_dataset is None:
                raise ValueError(f"File entry outside a dataset at line {line_number}")
            finish_file()
            current_file = {"name": _yaml_scalar(stripped.split(":", 1)[1])}
            continue

        if current_file is not None and ":" in stripped:
            key, value = stripped.split(":", 1)
            if key in {"sha256", "url"}:
                current_file[key] = _yaml_scalar(value)

    finish_file()
    if not datasets:
        raise ValueError(f"No datasets found in {path}")
    empty = [name for name, files in datasets.items() if not files]
    if empty:
        raise ValueError(f"Datasets contain no files: {', '.join(empty)}")
    return datasets


def sha256_file(path: Path) -> str:
    """Calculate a file's SHA-256 digest without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _format_size(size: int) -> str:
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB"):
        if value < 1024 or unit == "GiB":
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} GiB"


def download_file(
    spec: FileSpec,
    output_dir: Path,
    *,
    force: bool = False,
    retries: int = 2,
    timeout: float = 60,
    quiet: bool = False,
) -> str:
    """Download one file, verify it, and atomically move it into place."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / spec.name

    if destination.is_file() and not force:
        if sha256_file(destination) == spec.sha256:
            if not quiet:
                print(f"[cached] {spec.name}")
            return "cached"
        if not quiet:
            print(f"[replace] {spec.name}: existing checksum does not match")

    last_error: Exception | None = None
    for attempt in range(1, retries + 2):
        temporary: Path | None = None
        try:
            request = Request(
                spec.url,
                headers={"User-Agent": "ancom-bc-benchmark-data-fetcher/1.0"},
            )
            with tempfile.NamedTemporaryFile(
                mode="wb",
                prefix=f".{spec.name}.",
                suffix=".part",
                dir=output_dir,
                delete=False,
            ) as output:
                temporary = Path(output.name)
                digest = hashlib.sha256()
                downloaded = 0
                with urlopen(request, timeout=timeout) as response:
                    length = response.headers.get("Content-Length")
                    total = int(length) if length and length.isdigit() else None
                    while True:
                        chunk = response.read(CHUNK_SIZE)
                        if not chunk:
                            break
                        output.write(chunk)
                        digest.update(chunk)
                        downloaded += len(chunk)
                        if not quiet:
                            if total:
                                percent = downloaded / total * 100
                                progress = f"{percent:5.1f}% of {_format_size(total)}"
                            else:
                                progress = _format_size(downloaded)
                            print(f"\r[fetch]  {spec.name}: {progress}", end="", flush=True)
            if not quiet:
                print()
            actual = digest.hexdigest()
            if actual != spec.sha256:
                raise ValueError(
                    f"checksum mismatch for {spec.name}: expected {spec.sha256}, got {actual}"
                )
            os.replace(temporary, destination)
            if not quiet:
                print(f"[saved]  {destination}")
            return "downloaded"
        except (HTTPError, URLError, OSError, ValueError) as error:
            last_error = error
            if not quiet and attempt <= retries:
                print(f"[retry]  {spec.name}: {error} (attempt {attempt}/{retries + 1})")
                time.sleep(min(2**attempt, 5))
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    assert last_error is not None
    raise RuntimeError(f"Could not download {spec.name}: {last_error}") from last_error


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download benchmark datasets from Zenodo and verify SHA-256 checksums."
    )
    parser.add_argument(
        "datasets",
        nargs="*",
        metavar="DATASET",
        help="Dataset groups to fetch; omit to fetch every group.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="Manifest path (default: data/manifest.yaml beside this script).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Download directory (default: data/ beside this script).",
    )
    parser.add_argument("--force", action="store_true", help="Download even valid cached files.")
    parser.add_argument("--list", action="store_true", help="List dataset groups and exit.")
    parser.add_argument("--retries", type=int, default=2, help="Retries after the first attempt.")
    parser.add_argument("--timeout", type=float, default=60, help="Network timeout in seconds.")
    parser.add_argument("--quiet", action="store_true", help="Suppress progress output.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.retries < 0:
        parser.error("--retries must be zero or greater")
    if args.timeout <= 0:
        parser.error("--timeout must be greater than zero")

    try:
        manifest = load_manifest(args.manifest)
    except (OSError, SyntaxError, ValueError) as error:
        parser.error(str(error))

    if args.list:
        for dataset, files in manifest.items():
            print(f"{dataset} ({len(files)} files)")
            for spec in files:
                print(f"  {spec.name}")
        return 0

    selected = list(dict.fromkeys(args.datasets)) if args.datasets else list(manifest)
    unknown = [dataset for dataset in selected if dataset not in manifest]
    if unknown:
        parser.error(
            f"unknown dataset group(s): {', '.join(unknown)}; "
            f"choose from {', '.join(manifest)}"
        )

    failures: list[str] = []
    cached = 0
    downloaded = 0
    for dataset in selected:
        if not args.quiet:
            print(f"\nDataset: {dataset}")
        for spec in manifest[dataset]:
            try:
                status = download_file(
                    spec,
                    args.output_dir,
                    force=args.force,
                    retries=args.retries,
                    timeout=args.timeout,
                    quiet=args.quiet,
                )
                cached += status == "cached"
                downloaded += status == "downloaded"
            except RuntimeError as error:
                failures.append(str(error))
                print(f"[failed] {error}", file=sys.stderr)

    if failures:
        print(
            f"\nFinished with {len(failures)} failure(s); "
            f"downloaded {downloaded}, cached {cached}.",
            file=sys.stderr,
        )
        return 1
    if not args.quiet:
        print(f"\nDone: downloaded {downloaded}, cached {cached} in {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
