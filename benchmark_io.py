"""Shared, repository-relative I/O helpers for the benchmark scripts."""

from __future__ import annotations

import os
import pickle
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"


def input_path(filename: str, data_dir: Path = DATA_DIR) -> Path:
    """Return a dataset path, accepting the old repository-root layout too."""
    preferred = Path(data_dir) / filename
    legacy = PROJECT_ROOT / filename
    if preferred.is_file():
        return preferred
    if legacy.is_file():
        return legacy
    raise FileNotFoundError(
        f"Dataset not found: {preferred}. Download it using data/manifest.yaml."
    )


def output_path(filename: str, output_dir: Path = RESULTS_DIR) -> Path:
    """Return an output path and create its parent directory."""
    path = output_dir / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def atomic_write_csv(frame: pd.DataFrame, path: Path, **kwargs: Any) -> None:
    """Write a DataFrame without exposing a partially written checkpoint."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = "".join(path.suffixes) or ".csv"
    with tempfile.NamedTemporaryFile(
        mode="wb", prefix=f".{path.name}.", suffix=suffix, dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
    try:
        frame.to_csv(temporary, **kwargs)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_pickle_dump(value: Any, path: Path) -> None:
    """Atomically persist a trusted, locally generated Python cache."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{path.name}.",
            suffix=".pkl",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def load_pickle(path: Path) -> Any:
    """Load a trusted cache created by :func:`atomic_pickle_dump`."""
    with Path(path).open("rb") as handle:
        return pickle.load(handle)
