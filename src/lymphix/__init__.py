"""Lymphix — BCR/TCR clonality from custom-panel NGS."""

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

__all__ = ["__version__"]


def _read_version() -> str:
    """Version from installed metadata, falling back to the VERSION file.

    The VERSION file is the single source of truth (bin/bump_version.py
    propagates it to CITATION.cff and dxapp.json), so a source checkout that
    has never been installed still reports the right number.
    """
    try:
        return version("lymphix")
    except PackageNotFoundError:
        pass

    for parent in Path(__file__).resolve().parents:
        candidate = parent / "VERSION"
        if candidate.is_file():
            return candidate.read_text(encoding="utf-8").strip()

    return "0+unknown"


__version__ = _read_version()
