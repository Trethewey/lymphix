#!/usr/bin/env python3
"""Bump the project version and propagate to all declaration sites.

Sites kept in sync:
  - VERSION                 (source of truth)
  - CITATION.cff            (version + date-released)
  - dxapp.json              (version)

  nextflow.config reads VERSION at runtime (manifest { version = file('VERSION')... });
  no patching needed there since 2026-06-19.

Usage:
    python bin/bump_version.py {patch|minor|major} [--dry-run]
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VERSION_FILE = ROOT / "VERSION"
CITATION_FILE = ROOT / "CITATION.cff"
DXAPP_FILE = ROOT / "dxapp.json"


def read_version() -> tuple[int, int, int]:
    raw = VERSION_FILE.read_text(encoding="utf-8").strip()
    parts = raw.split(".")
    if len(parts) != 3:
        sys.exit(f"VERSION malformed (need MAJOR.MINOR.PATCH): {raw!r}")
    try:
        return int(parts[0]), int(parts[1]), int(parts[2])
    except ValueError:
        sys.exit(f"VERSION malformed (non-integer component): {raw!r}")


def bump(current: tuple[int, int, int], kind: str) -> str:
    major, minor, patch = current
    if kind == "major":
        return f"{major + 1}.0.0"
    if kind == "minor":
        return f"{major}.{minor + 1}.0"
    if kind == "patch":
        return f"{major}.{minor}.{patch + 1}"
    sys.exit(f"unknown bump kind: {kind}")


def write_version_file(new: str) -> None:
    VERSION_FILE.write_text(new + "\n", encoding="utf-8")


def update_citation(new: str, today: str) -> None:
    text = CITATION_FILE.read_text(encoding="utf-8")
    text, n_ver = re.subn(
        r'^version:\s*".*"$', f'version: "{new}"', text, count=1, flags=re.MULTILINE
    )
    text, n_date = re.subn(
        r'^date-released:\s*".*"$',
        f'date-released: "{today}"',
        text,
        count=1,
        flags=re.MULTILINE,
    )
    if not n_ver or not n_date:
        sys.exit("CITATION.cff: failed to locate version/date-released lines")
    CITATION_FILE.write_text(text, encoding="utf-8")


def update_dxapp(new: str) -> None:
    data = json.loads(DXAPP_FILE.read_text(encoding="utf-8"))
    if "version" not in data:
        sys.exit("dxapp.json: no 'version' key")
    data["version"] = new
    DXAPP_FILE.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


# nextflow.config no longer needs patching — it reads VERSION at runtime.


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("kind", choices=["patch", "minor", "major"])
    parser.add_argument("--dry-run", action="store_true", help="print the new version without writing")
    args = parser.parse_args()

    current = read_version()
    new = bump(current, args.kind)
    today = _dt.datetime.now().strftime("%Y-%m-%d")

    old_str = ".".join(str(p) for p in current)
    print(f"{old_str} -> {new}  ({today})")
    if args.dry_run:
        return

    write_version_file(new)
    update_citation(new, today)
    update_dxapp(new)

    print("Updated: VERSION, CITATION.cff, dxapp.json  (nextflow.config auto-reads VERSION)")
    print(f"Next steps:")
    print(f"  1. Review:  git -C app/lymphix diff   (optional — only if you locally git-track app/)")
    print(f"  2. Sync:    .\\sync_to_github.ps1")
    print(f"  3. In github/lymphix/:  git add -A && git commit -m 'Release v{new}' && git tag v{new}")
    print(f"  4. Push:    git push && git push --tags")


if __name__ == "__main__":
    main()
