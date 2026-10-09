#!/usr/bin/env python3
"""Verify downloaded benchmark datasets against data/manifest.yaml."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from fetch_data import (
    DEFAULT_MANIFEST,
    DEFAULT_OUTPUT_DIR,
    FileSpec,
    load_manifest,
    sha256_file,
)


@dataclass(frozen=True)
class VerificationResult:
    """Checksum verification result for one manifest file."""

    spec: FileSpec
    path: Path
    status: str
    actual_sha256: str | None = None
    error: str | None = None


def verify_file(spec: FileSpec, data_dir: Path) -> VerificationResult:
    """Verify that one downloaded file exists and matches its manifest digest."""
    path = Path(data_dir) / spec.name
    if not path.exists():
        return VerificationResult(spec, path, "missing")
    if not path.is_file():
        return VerificationResult(spec, path, "error", error="path is not a regular file")

    try:
        actual = sha256_file(path)
    except OSError as error:
        return VerificationResult(spec, path, "error", error=str(error))

    if actual != spec.sha256:
        return VerificationResult(spec, path, "invalid", actual_sha256=actual)
    return VerificationResult(spec, path, "valid", actual_sha256=actual)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify downloaded benchmark data using manifest SHA-256 checksums."
    )
    parser.add_argument(
        "datasets",
        nargs="*",
        metavar="DATASET",
        help="Dataset groups to verify; omit to verify every group.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="Manifest path (default: data/manifest.yaml beside this script).",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory containing downloaded files (default: data/ beside this script).",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Show failures and the final summary, but hide valid files.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        manifest = load_manifest(args.manifest)
    except (OSError, SyntaxError, ValueError) as error:
        parser.error(str(error))

    selected = list(dict.fromkeys(args.datasets)) if args.datasets else list(manifest)
    unknown = [dataset for dataset in selected if dataset not in manifest]
    if unknown:
        parser.error(
            f"unknown dataset group(s): {', '.join(unknown)}; "
            f"choose from {', '.join(manifest)}"
        )

    results: list[VerificationResult] = []
    for dataset in selected:
        if not args.quiet:
            print(f"\nDataset: {dataset}")
        for spec in manifest[dataset]:
            result = verify_file(spec, args.data_dir)
            results.append(result)
            if result.status == "valid":
                if not args.quiet:
                    print(f"[valid]   {spec.name}")
            elif result.status == "missing":
                print(f"[missing] {result.path}")
            elif result.status == "invalid":
                print(f"[invalid] {result.path}")
                print(f"          expected: {spec.sha256}")
                print(f"          actual:   {result.actual_sha256}")
            else:
                print(f"[error]   {result.path}: {result.error}")

    counts = {
        status: sum(result.status == status for result in results)
        for status in ("valid", "missing", "invalid", "error")
    }
    print(
        "\nVerification summary: "
        f"{counts['valid']} valid, {counts['missing']} missing, "
        f"{counts['invalid']} invalid, {counts['error']} errors."
    )
    failures = counts["missing"] + counts["invalid"] + counts["error"]
    if failures:
        print(
            "Run fetch_data.py for the affected dataset groups to restore the files.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
