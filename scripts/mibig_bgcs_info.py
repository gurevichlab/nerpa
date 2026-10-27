#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable
from src.antismash_parsing.bgc_variant_types import BGC_Variant
import yaml

import polars as pl
from collections import defaultdict


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description='Build mibig_bgcs_info.tsv by combining BIG-SCAPE families, BGC variants, rBAN graphs, and approved matches.'
    )
    p.add_argument(
        '--bigscape-results-dir',
        type=Path,
        help='BIG-SCAPE clustering of MIBiG 4. Specifically, network_files results directory containing per-type subdirs',
    )
    p.add_argument(
        '--bigscape-cutoff',
        default='c0.30',
        help="BIG-SCAPE cutoff suffix used in '{type}_clustering_{cutoff}.tsv' (default: %(default)s)",
    )

    p.add_argument(
        '--preprocessed-bgcs',
        type=Path,
        help='YAML with Nerpa-preprocessed MIBiG BGCs',
    )

    p.add_argument(
        '--nrps-info',
        type=Path,
        help='TSV with NRP IDs',
    )

    p.add_argument(
        '--mibig3-jsons',
        type=Path,
        help=(
            'Path to a folder with MIBiG 3.1 JSONs. '
            'Used to assign the MIBiG version to each BGC. '
        ),
    )

    p.add_argument(
        '--out',
        type=Path,
        help='Output TSV path',
    )
    return p.parse_args()


def iter_bigscape_types(bigscape_results_dir: Path) -> list[str]:
    return sorted([
        p.name
        for p in bigscape_results_dir.iterdir()
        if p.is_dir()
    ])


BGC_ID = str
BIGSCAPE_FAMILY_FINGERPRINT = str


def load_bigscape_family_map(
    bigscape_results_dir: Path,
    bigscape_types: Iterable[str],
    cutoff: str,
) -> dict[BGC_ID, BIGSCAPE_FAMILY_FINGERPRINT]:
    """
    Build mapping: bgc_id -> 'TYPE_familyNum_TYPE_familyNum...'
    """
    families_by_bgc: dict[str, list[str]] = {}

    for bgc_type in bigscape_types:
        clustering_path = bigscape_results_dir / bgc_type / f'{bgc_type}_clustering_{cutoff}.tsv'
        df: pl.DataFrame = pl.read_csv(clustering_path, separator='\t')

        # Intent: replicate "if bgc_id in table then append '{type}_{family_number}'"
        # but in bulk by building a bgc->families list.
        for row in (
            df.select([pl.col('#BGC Name').alias('bgc_id'), pl.col('Family Number').alias('family_number')])
            .iter_rows(named=True)
        ):
            bgc_id = str(row['bgc_id'])
            family_number = row['family_number']
            fam = f'{bgc_type}_{family_number}'
            families_by_bgc.setdefault(bgc_id, []).append(fam)

    return {bgc_id: '_'.join(fams) for bgc_id, fams in families_by_bgc.items()}


def nrp_id_to_bgc_id(nrp_id: str) -> str:
    """
    Convert NRP ID to BGC ID by removing the last part after the last dot.
    Example: 'BGC000001.1' -> 'BGC000001'
    """
    return nrp_id.rsplit('.', 1)[0]


def main() -> None:
    args = parse_args()

    bigscape_types: list[str] = iter_bigscape_types(args.bigscape_results_dir)
    family_map: dict[str, str] = load_bigscape_family_map(
        bigscape_results_dir=args.bigscape_results_dir,
        bigscape_types=bigscape_types,
        cutoff=args.bigscape_cutoff,
    )
    bgc_variants = [
        BGC_Variant.from_dict(yaml_dict)
        for yaml_dict in yaml.safe_load(args.preprocessed_bgcs.open('r'))
    ]
    bgc_variant_ids = {
        bgc.bgc_variant_id._genome_id()
        for bgc in bgc_variants
    }
    missing_ids = bgc_variant_ids - set(family_map.keys())
    if missing_ids:
        print(
            'Warning: some BGCs from preprocessed MIBiG 4 set '
            f'are not present in the BIG-SCAPE results '
            f'({len(missing_ids)} in total):\n'
            f'{missing_ids}\n'
            'They will be assigned fictional unique families'
        )
    for i, missing_id in enumerate(missing_ids):
        family_map[missing_id] = f'fictional_{i}'

    bgc_ids = set(family_map.keys())

    mibig3_jsons = list(args.mibig3_jsons.glob('*.json'))
    mibig3_bgc_ids = {json_path.stem for json_path in mibig3_jsons}
    bgc_source = {
        bgc_id: (
            'MIBiG 3.1'
            if bgc_id in mibig3_bgc_ids
            else 'MIBiG 4.0'
        )
        for bgc_id in bgc_ids
    }
    num_a_domains_map: dict[str, int] = {
        bgc_variant.bgc_variant_id._genome_id(): len(bgc_variant.modules)
        for bgc_variant in bgc_variants
    }

    pnrpdb_df: pl.DataFrame = pl.read_csv(args.nrps_info, separator='\t')
    bgc_id_to_nrp_ids = defaultdict(list)
    for nrp_id in pnrpdb_df['ID'].to_list():
        bgc_id = nrp_id_to_bgc_id(nrp_id)
        bgc_id_to_nrp_ids[bgc_id].append(nrp_id)

    table = pl.DataFrame(
        {
            'bgc_id': list(bgc_ids),
            'source': [
                bgc_source[bgc_id]
                for bgc_id in bgc_ids
            ],
            'bigscape_families': [
                family_map[bgc_id]
                for bgc_id in bgc_ids
            ],
            'num_a_domains': [
                num_a_domains_map.get(bgc_id, None)
                for bgc_id in bgc_ids
            ],
            'nrp_ids': [
                ','.join(bgc_id_to_nrp_ids.get(bgc_id, []))
                for bgc_id in bgc_ids
            ],
        }
    )

    table.write_csv(args.out, separator='\t')


if __name__ == '__main__':
    main()
