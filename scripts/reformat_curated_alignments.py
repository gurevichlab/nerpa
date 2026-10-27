#!/usr/bin/env python3

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO, TypeAlias, cast

import yaml


JSONScalar: TypeAlias = str | int | float | bool | None
JSONValue: TypeAlias = JSONScalar | list["JSONValue"] | dict[str, "JSONValue"]


@dataclass(frozen=True)
class Arguments:
    input: Path
    output: Path | None
    indent: int


def parse_args() -> Arguments:
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="Convert YAML/JSON to specially formatted JSON."
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Input .yaml, .yml, or .json file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output file (default: .json for YAML, .formatted.json for JSON)",
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="Indentation width (default: 2)",
    )

    namespace: argparse.Namespace = parser.parse_args()

    return Arguments(
        input=namespace.input,
        output=namespace.output,
        indent=namespace.indent,
    )


def load_data(path: Path) -> JSONValue:
    suffix: str = path.suffix.lower()

    with path.open(encoding="utf-8") as input_file:
        file: TextIO = input_file

        if suffix == ".json":
            data: Any = json.load(file)
        elif suffix in {".yaml", ".yml"}:
            data = yaml.safe_load(file)
        else:
            raise ValueError(f"Unsupported input format: {suffix}")

    return cast(JSONValue, data)


def format_alignment(
    alignment: list[JSONValue],
    indent: int,
    level: int,
) -> str:
    if not alignment:
        return "[]"

    prefix: str = " " * (indent * level)
    child_prefix: str = " " * (indent * (level + 1))

    rows: list[str] = [
        child_prefix
        + json.dumps(
            row,
            ensure_ascii=False,
            separators=(", ", ": "),
        )
        for row in alignment
    ]

    return "[\n" + ",\n".join(rows) + "\n" + prefix + "]"


def format_alternative_alignments(
    alignments: list[JSONValue],
    indent: int,
    level: int,
) -> str:
    if not alignments:
        return "[]"

    prefix: str = " " * (indent * level)
    child_prefix: str = " " * (indent * (level + 1))
    rendered_alignments: list[str] = []

    for alignment in alignments:
        if not isinstance(alignment, list):
            raise ValueError(
                "Each item in 'acceptable_alternative_alignments' "
                "must be an alignment list"
            )

        rendered: str = format_alignment(
            alignment,
            indent,
            level=level + 1,
        )
        rendered_alignments.append(child_prefix + rendered)

    return (
        "[\n"
        + ",\n".join(rendered_alignments)
        + "\n"
        + prefix
        + "]"
    )


def format_json(
    value: JSONValue,
    indent: int,
    level: int = 0,
    key: str | None = None,
) -> str:
    prefix: str = " " * (indent * level)
    child_prefix: str = " " * (indent * (level + 1))

    if key == "true_alignment":
        if not isinstance(value, list):
            raise ValueError("'true_alignment' must be a list")

        return format_alignment(value, indent, level)

    if key == "acceptable_alternative_alignments":
        if not isinstance(value, list):
            raise ValueError(
                "'acceptable_alternative_alignments' must be a list"
            )

        return format_alternative_alignments(value, indent, level)

    if isinstance(value, dict):
        if not value:
            return "{}"

        object_rows: list[str] = []

        for child_key, child_value in value.items():
            rendered: str = format_json(
                child_value,
                indent,
                level=level + 1,
                key=child_key,
            )

            object_rows.append(
                child_prefix
                + json.dumps(child_key, ensure_ascii=False)
                + ": "
                + rendered
            )

        return "{\n" + ",\n".join(object_rows) + "\n" + prefix + "}"

    if isinstance(value, list):
        if not value:
            return "[]"

        list_rows: list[str] = [
            child_prefix + format_json(item, indent, level=level + 1)
            for item in value
        ]

        return "[\n" + ",\n".join(list_rows) + "\n" + prefix + "]"

    return json.dumps(value, ensure_ascii=False)


def default_output_path(input_path: Path) -> Path:
    if input_path.suffix.lower() == ".json":
        return input_path.with_name(f"{input_path.stem}.formatted.json")

    return input_path.with_suffix(".json")


def main() -> None:
    args: Arguments = parse_args()

    if args.indent < 0:
        raise SystemExit("--indent must be non-negative")

    output: Path = args.output or default_output_path(args.input)

    if output.resolve() == args.input.resolve():
        raise SystemExit("Refusing to overwrite the input file")

    data: JSONValue = load_data(args.input)
    rendered: str = format_json(data, args.indent) + "\n"

    output.write_text(rendered, encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
