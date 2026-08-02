import argparse
from pathlib import Path
from typing import (
    Dict,
    List,
    NamedTuple,
)

import polars as pl
import yaml
import networkx as nx

from src.config import load_monomer_names_helper
from src.monomer_names_helper import NOT_NRPS_RESIDUE, PKS_RESIDUE, UNKNOWN_RESIDUE
from src.pipeline.deduplication import (
    cluster_isomorphic_nrp_variants,
    reassign_representative_nrps,
)
from src.rban_parsing.nrp_variant_types import NRP_Variant, NRP_Variant_ID
from src.rban_parsing.rban_parser import (
    Parsed_rBAN_Record,
    NorineMonomerName
)
from src.rban_parsing.retrieve_nrp_variants import rban_records_to_nrp_variants
from collections import Counter



def cluster_isomorphic_rban_graphs(rban_records: List[Parsed_rBAN_Record]) -> Dict[str, str]:
    compound_id_to_representative: Dict[str, Parsed_rBAN_Record] = {}

    for i, rban_record in enumerate(rban_records, start=1):
        print(f'{i}/{len(rban_records)} Processing compound {rban_record.compound_id}...')

        representative = None
        for repr_candidate in compound_id_to_representative.values():
            if rban_record.is_isomorphic_to(repr_candidate):
                representative = repr_candidate
                break

        if representative is None:
            representative = rban_record  # if none found, use itself as representative

        compound_id_to_representative[rban_record.compound_id] = representative

    compound_id_to_repr_id: Dict[str, str] = {
        compound_id: repr_record.compound_id
        for compound_id, repr_record in compound_id_to_representative.items()
    }

    return reassign_representative_nrps(compound_id_to_repr_id)


def parse_args(nerpa_root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Compute various stats on PNRPDB2 compounds')
    default_pnrpdb2_preprocessed = (
        nerpa_root
        / 'data'
        / 'input'
        / 'preprocessed'
        / 'pnrpdb2_raw_parsed_rban_records.yaml'
    )
    parser.add_argument(
        '--pnrpdb2-preprocessed',
        type=Path,
        required=True,
        default=default_pnrpdb2_preprocessed,
        help=f'Path to the preprocessed PNRPDB2 data (default: {default_pnrpdb2_preprocessed})',
    )

    default_pnrpdb2_raw = (
        nerpa_root
        / 'data'
        / 'input'
        / 'pnrpdb2_raw.tsv'
    )
    parser.add_argument(
        '--pnrpdb2-raw',
        type=Path,
        default=default_pnrpdb2_raw,
        help=(
            f'Path to the raw PNRPDB2 data (default: {default_pnrpdb2_raw}). '
            'Used to extract IDs of compounds that rBAN was unable to process.'
        ),
    )
        
    default_output_table = (
        nerpa_root
        / 'data'
        / 'for_training_and_testing'
        / 'pnrpdb2_info.tsv'
    )
    parser.add_argument(
        '--out',
        type=Path,
        default=default_output_table,
        help=f'Path to the output table with PNRPDB2 info (default: {default_output_table})',
    )
    return parser.parse_args()


def main():
    nerpa_dir = Path(__file__).parent.parent
    assert (nerpa_dir / 'nerpa.py').exists(), f'Could not find nerpa.py in {nerpa_dir}'

    monomers_cfg_file = (
        nerpa_dir
        / 'configs'
        / 'monomers_config.yaml'
    )
    monomer_names_helper = load_monomer_names_helper(monomers_cfg_file, nerpa_dir)

    args = parse_args(nerpa_dir)

    # input
    preprocessed_dir = args.pnrpdb2_preprocessed.parent
    print(f'Loading preprocessed data from {preprocessed_dir}...')
    rban_records = [
        Parsed_rBAN_Record.from_dict(record_dict)
        for record_dict in yaml.safe_load(args.pnrpdb2_preprocessed.read_text())
    ]
    print(f'Loaded {len(rban_records)} rBAN records')
    nrp_variants = rban_records_to_nrp_variants(rban_records)
    print(f'Loaded {len(nrp_variants)} NRP variants')

    def is_rban_recognized_monomer(monomer_name: str) -> bool:
        return not any([monomer_name.startswith('X'),
                       ':' in monomer_name,])

    def is_nrps_monomer(monomer_name: str) -> bool:
        parsed_name = monomer_names_helper.parsed_name(
            monomer_name,
            name_format='rBAN/Norine'
        )
        return (
            is_rban_recognized_monomer(monomer_name) and
            parsed_name.residue not in (NOT_NRPS_RESIDUE, PKS_RESIDUE)
        )

    def is_nerpa_supported_monomer(monomer_name: str) -> bool:
        parsed_name = monomer_names_helper.parsed_name(
            monomer_name,
            name_format='rBAN/Norine'
        )
        return (
            is_rban_recognized_monomer(monomer_name) and
            parsed_name.residue not in (NOT_NRPS_RESIDUE, PKS_RESIDUE, UNKNOWN_RESIDUE)
        )

    print('Clustering isomorphic NRP variants...')
    nrp_id_to_repr_id = cluster_isomorphic_nrp_variants(nrp_variants)
    print(f'Found {len(set(nrp_id_to_repr_id.values()))} clusters of isomorphic NRP variants')
    print('Deduplicating rBAN graphs...')
    rban_id_to_repr_id = cluster_isomorphic_rban_graphs(rban_records)
    print(f'Found {len(set(rban_id_to_repr_id.values()))} clusters of isomorphic rBAN graphs')

    # prepare table
    pnrpdb_info_rows = []
    for record in rban_records:
        compound_id = record.compound_id
        monomers = [
            monomer_info.name for monomer_info in record.monomers.values()
        ]

        nrp_variant_id = NRP_Variant_ID(nrp_id=compound_id, variant_idx=0)
        if nrp_variant_id in nrp_id_to_repr_id:
            nrp_variant_representative = nrp_id_to_repr_id[nrp_variant_id].nrp_id
        else:
            nrp_variant_representative = None  # Nerpa was unable to build an NRP variant for this compound

        pnrpdb_info_rows.append({
            'compound_id': compound_id,
            'total_num_monomers':
                len(monomers),
            'num_rban_recognized_monomers':
                len(list(filter(is_rban_recognized_monomer, monomers))),
            'num_nrps_monomers':
                len(list(filter(is_nrps_monomer, monomers))),
            'num_nerpa_supported_monomers':
                len(list(filter(is_nerpa_supported_monomer, monomers))),
            'nrp_variant_iso_class_representative':
                nrp_variant_representative,
            'rban_graph_iso_class_representative':
                rban_id_to_repr_id[compound_id],
        })

    # Add compounds that Nerpa was unable to process
    pnrpdb_raw = pl.read_csv(args.pnrpdb2_raw, separator='\t')
    all_ids = set(pnrpdb_raw['ID'].to_list())
    missing_ids = all_ids - set(record.compound_id for record in rban_records)
    for compound_id in missing_ids:
        pnrpdb_info_rows.append({
            'compound_id': compound_id,
            'total_num_monomers': None,
            'num_rban_recognized_monomers': None,
            'num_nrps_monomers': None,
            'nrp_variant_iso_class_representative': None,
            'rban_graph_iso_class_representative': None,
        })

    pnrpdb_info = (
        pl.DataFrame(pnrpdb_info_rows)
        .sort('compound_id')
    )

    print(f'Writing table with info about PNRPDB2 NRP variants to {args.out}')
    args.out.write_text(pnrpdb_info.write_csv(separator='\t'))
    print('Done')


if __name__ == "__main__":
    main()
