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

def rban_record_to_nx_graph(rban_record: Parsed_rBAN_Record) -> nx.Graph:
    graph = nx.Graph()
    for monomer_idx, monomer_info in rban_record.monomers.items():
        graph.add_node(monomer_idx, name=monomer_info.name)
    for (mon1_idx, mon2_idx), edge_info in rban_record.monomer_bonds.items():
        bond_type = edge_info.bondType
        graph.add_edge(mon1_idx, mon2_idx, bond_type=bond_type)
    return graph


def node_name_to_compare(node_data: dict) -> str:
    rban_name = node_data['name']
    if rban_name.startswith('X'):
        return 'X'       # unknown
    if ':' in rban_name:
        return 'LIPID'   # lipid tail
    return rban_name          # specific residue


def edge_name_to_compare(edge_data: dict) -> str | None:
    return edge_data.get('bond_type')

class rBAN_Graph_Signature(NamedTuple):
    num_nodes: int
    num_edges: int
    nodes_signature: Counter[str]
    edges_signature: Counter[str | None]

    
def rban_graph_signature(
    g: nx.Graph,
) -> rBAN_Graph_Signature:
    node_sig: Counter[str] = Counter(
        node_name_to_compare(data)
        for _, data in g.nodes(data=True)
    )
    edge_sig: Counter[str | None] = Counter(
        edge_name_to_compare(data)
        for _, _, data in g.edges(data=True)
    )
    return rBAN_Graph_Signature(
        g.number_of_nodes(),
        g.number_of_edges(),
        node_sig,
        edge_sig
    )

def rban_node_match(node1_data: dict, node2_data: dict) -> bool:
    return (
        node_name_to_compare(node1_data) ==
        node_name_to_compare(node2_data)
    )

def rban_edge_match(edge1_data: dict, edge2_data: dict) -> bool:
    return (
        edge_name_to_compare(edge1_data) ==
        edge_name_to_compare(edge2_data)
    )


def cluster_isomorphic_rban_graphs(rban_records: List[Parsed_rBAN_Record]) -> Dict[str, str]:
    graphs = {
        record.compound_id: rban_record_to_nx_graph(record)
        for record in rban_records
    }
    graph_signatures = {
        compound_id: rban_graph_signature(graph)
        for compound_id, graph in graphs.items()
    }
    compound_id_to_representative: Dict[str, str] = {}

    for i, (compound_id, graph) in enumerate(graphs.items()):
        print(f'{i}/{len(graphs)} Processing compound {compound_id}...')

        representative_id = None
        for repr_id_candidate in compound_id_to_representative.values():
            if graph_signatures[compound_id] != graph_signatures[repr_id_candidate]:
                continue  # signatures differ, cannot be isomorphic

            if nx.is_isomorphic(
                    graph,
                    graphs[repr_id_candidate],
                    node_match=rban_node_match,
                    edge_match=rban_edge_match
            ):
                representative_id = repr_id_candidate
                break

        if representative_id is None:
            representative_id = compound_id  # if none found, use itself as representative

        compound_id_to_representative[compound_id] = representative_id

    return reassign_representative_nrps(compound_id_to_representative)


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
