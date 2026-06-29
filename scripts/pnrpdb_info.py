import argparse
from pathlib import Path
from typing import Dict, List

import polars as pl
import yaml
import networkx as nx

from src.config import load_monomer_names_helper
from src.monomer_names_helper import NOT_NRPS_RESIDUE, PKS_RESIDUE, UNKNOWN_RESIDUE
from src.pipeline.deduplication import cluster_isomorphic_nrp_variants
from src.rban_parsing.nrp_variant_types import NRP_Variant, NRP_Variant_ID
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.rban_parsing.retrieve_nrp_variants import rban_records_to_nrp_variants

def rban_record_to_nx_graph(rban_record: Parsed_rBAN_Record) -> nx.Graph:
    graph = nx.Graph()
    for monomer_idx, monomer_info in rban_record.monomers.items():
        graph.add_node(monomer_idx, name=monomer_info.name)
    for (mon1_idx, mon2_idx), edge_info in rban_record.monomer_bonds.items():
        bond_type = edge_info.bondType
        graph.add_edge(mon1_idx, mon2_idx, bond_type=bond_type)
    return graph

def cluster_isomorphic_rban_graphs(rban_records: List[Parsed_rBAN_Record]) -> Dict[str, str]:
    graphs = {
        record.compound_id: rban_record_to_nx_graph(record)
        for record in rban_records
    }
    compound_id_to_representative: Dict[str, str] = {}

    def node_match(node1_data: dict, node2_data: dict) -> bool:
        name1, name2 = node1_data['name'], node2_data['name']
        return any([
            name1 == name2,  # same amino acid
            name1.startswith('X') and name2.startswith('X'),  # both unknown
            ':' in name1 and ':' in name2   # both lipid tails
        ])

    def edge_match(edge1_data: dict, edge2_data: dict) -> bool:
        return edge1_data.get('bond_type') == edge2_data.get('bond_type')

    for i, (compound_id, graph) in enumerate(graphs.items()):
        print(f'{i}/{len(graphs)} Processing compound {compound_id}...')
        representative_id = next(
            (representative_id   # find first isomorphic graph from already found representatives
             for representative_id in compound_id_to_representative.values()
             if nx.is_isomorphic(
                     graph,
                     graphs[representative_id],
                     node_match=node_match,
                     edge_match=edge_match
             )),
            compound_id  # if none found, use itself as representative
        )

        compound_id_to_representative[compound_id] = representative_id

    return compound_id_to_representative


def parse_args(nerpa_root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Compute various stats on PNRPDB2 compounds')
    default_pnrpdb2_preprocessed = (
        nerpa_root
        / 'data'
        / 'input'
        / 'preprocessed'
        / 'pnrpdb2_parsed_rban_records.yaml'
    )
    parser.add_argument(
        '--pnrpdb2-preprocessed',
        type=Path,
        default=default_pnrpdb2_preprocessed,
        help=f'Path to the preprocessed PNRPDB2 data (default: {default_pnrpdb2_preprocessed})',
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
        parsed_name = monomer_names_helper.parsed_name(monomer_name,
                                                       name_format='rBAN/Norine')
        return (is_rban_recognized_monomer(monomer_name) and
                parsed_name.residue not in (NOT_NRPS_RESIDUE, PKS_RESIDUE))

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
            'nrp_variant_iso_class_representative':
                nrp_variant_representative,
            'rban_graph_iso_class_representative':
                rban_id_to_repr_id[compound_id],
        })

    pnrpdb_info = pl.DataFrame(pnrpdb_info_rows)

    print(f'Writing table with info about PNRPDB2 NRP variants to {args.out}')
    args.out.write_text(pnrpdb_info.write_csv(separator='\t'))
    print('Done')


if __name__ == "__main__":
    main()
