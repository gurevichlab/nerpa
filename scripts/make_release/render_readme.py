#!/usr/bin/env python3

"""Render README.md from README_template.md."""

from __future__ import annotations

import argparse
import re
from datetime import date
from pathlib import Path


PLACEHOLDER_RE = re.compile(r"{{\s*([a-z][a-z0-9-]*)\s*}}")
VERSION_RE = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+$")


def parse_args(REPO_ROOT: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render the Nerpa README from README_template.md."
    )
    parser.add_argument(
        "--template",
        type=Path,
        default=REPO_ROOT / "docs" / "README_template.md",
        help="README template (default: README_template.md)",
    )
    parser.add_argument(
        "--version-file",
        type=Path,
        default=REPO_ROOT / "VERSION.txt",
        help="file containing the X.Y.Z release version (default: VERSION.txt)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "README.md",
        help="rendered README path (default: README.md)",
    )
    parser.add_argument(
        "--release-date",
        default=date.today().isoformat(),
        help="release date in ISO 8601 format (default: today's date)",
    )
    parser.add_argument(
        "--release-tag",
        help="GitHub release tag (default: nerpa-ms-X.Y.Z)",
    )
    return parser.parse_args()


def main() -> None:
    REPO_ROOT = Path(__file__).resolve().parents[2]
    assert (REPO_ROOT / "docs/README_template.md").exists(), "Expected to be run from the Nerpa repo root"
    args = parse_args(REPO_ROOT)

    version = args.version_file.read_text(encoding="utf-8").strip()
    if not VERSION_RE.fullmatch(version):
        raise SystemExit(
            f"Error: {args.version_file} does not contain an X.Y.Z version."
        )

    try:
        date.fromisoformat(args.release_date)
    except ValueError as error:
        raise SystemExit(
            f"Error: invalid ISO 8601 release date: {args.release_date!r}"
        ) from error

    stable_tag = f"nerpa-ms-{version}"
    release_tag = args.release_tag or stable_tag
    if release_tag != stable_tag and not release_tag.startswith(f"{stable_tag}-"):
        raise SystemExit(
            f"Error: release tag must be {stable_tag!r} or a prerelease of it."
        )

    major, minor, _patch = version.split(".")
    values = {
        "release-series": f"{major}.{minor}",
        "release-version": version,
        "release-date": args.release_date,
        "release-tag": release_tag,
    }

    template = args.template.read_text(encoding="utf-8")

    def replace_placeholder(match: re.Match[str]) -> str:
        name = match.group(1)
        try:
            return values[name]
        except KeyError:
            raise SystemExit(
                f"Error: unsupported placeholder in {args.template}: {name!r}"
            ) from None

    rendered = PLACEHOLDER_RE.sub(replace_placeholder, template)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = args.output.with_name(f".{args.output.name}.tmp")
    temporary_output.write_text(rendered, encoding="utf-8")
    temporary_output.replace(args.output)

    print(f"Created {args.output}")


if __name__ == "__main__":
    main()
