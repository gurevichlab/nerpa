from __future__ import annotations

import argparse
import ast
import csv
from pathlib import Path
from typing import cast


OLD_TAG = "alpha-hydroxy amino acid"
NEW_TAG = "alpha-hydroxy acid"

# These abbreviations describe the complete monomer, not its Nerpa residue.
ABBREVIATED_NAME_ENDINGS = (
    "hic",
    "hiv",
    "hmp",
    "hva",
    "lac",
)

# Names whose alpha-hydroxy character is not captured by an abbreviation.
DESCRIPTIVE_NAMES = {
    "2,3-dihydroxyhexadecanoic acid",
    "2,4-dihydroxypentanoic acid",
    "2-hydroxy-4-methylpentanoic acid",
    "2-hydroxypent-4-enoic acid",
    "2r-hydroxy-3-methylpentanoic acid",
    "2r-hydroxyisovaleric acid",
    "2s-hydroxyisocaproic acid",
    "2s-hydroxyisovaleric acid",
    "d-phenyllactic acid",
    "lactic acid",
    "malic acid",
    "p-hydroxymandelate",
    "phenyllactic acid",
}


def parse_args() -> argparse.Namespace:
    repository_root = Path(__file__).resolve().parents[1]
    default_table = (
        repository_root
        / "data"
        / "nerpa_internal"
        / "monomer_names_table.tsv"
    )

    parser = argparse.ArgumentParser(
        description="Derive the alpha-hydroxy acid tag from monomer names.",
    )
    parser.add_argument(
        "table",
        nargs="?",
        type=Path,
        default=default_table,
        help="Table to update in place",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Write to another path instead of updating the input table",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show changed rows without writing the table",
    )
    return parser.parse_args()


def normalize_name(name: str) -> str:
    return " ".join(name.lower().replace("_", " ").split())


def is_alpha_hydroxy_acid(name: str) -> bool:
    normalized_name = normalize_name(name)

    if normalized_name in DESCRIPTIVE_NAMES:
        return True

    compact_name = normalized_name.replace(" ", "")
    return compact_name.endswith(ABBREVIATED_NAME_ENDINGS)


def parse_tags(serialized_tags: str) -> list[str]:
    tags = ast.literal_eval(serialized_tags)
    return cast(list[str], tags)


def update_tags(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    for row in rows:
        old_tags = parse_tags(row["Tags"])
        new_tags = [
            tag
            for tag in old_tags
            if tag not in {OLD_TAG, NEW_TAG}
        ]

        if is_alpha_hydroxy_acid(row["MonomerName"]):
            new_tags.append(NEW_TAG)

        row["Tags"] = repr(new_tags)

    return rows


def main() -> None:
    args = parse_args()
    input_path: Path = args.table
    output_path: Path = args.output or input_path

    with input_path.open(encoding="utf-8", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)

    if fieldnames is None:
        raise ValueError(f"Missing header in {input_path}")

    original_tags = [row["Tags"] for row in rows]
    updated_rows = update_tags(rows)

    changed_rows = [
        row
        for row, old_tags in zip(updated_rows, original_tags, strict=True)
        if row["Tags"] != old_tags
    ]

    for row in changed_rows:
        print(
            f"{row['MonomerName']}\t"
            f"{row['NameFormat']}\t"
            f"{row['Tags']}"
        )

    print(f"\nChanged {len(changed_rows)} rows.")

    if args.dry_run:
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(f"{output_path.name}.tmp")

    with temporary_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(updated_rows)

    temporary_path.replace(output_path)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()