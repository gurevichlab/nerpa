from pathlib import Path
import polars as pl
from src.config import load_monomer_names_helper
from src.monomer_names_helper import (
    NOT_NRPS_RESIDUE,
    PKS_RESIDUE, UNKNOWN_RESIDUE)
import yaml
from src.pipeline.deduplication import cluster_isomorphic_nrp_variants
from src.rban_parsing.nrp_variant_types import NRP_Variant_ID, NRP_Variant
from typing import List, Dict
from networkx import DiGraph, is_isomorphic


def compact_dict_to_graph(compact: Dict) -> DiGraph:
    G = DiGraph()
    for mon_idx, name in compact['nodes'].items():
        G.add_node(mon_idx, name=name)
    for u, v, bond_type in compact['edges']:
        G.add_edge(u, v, bond_type=bond_type)
    return G

def node_match(n1, n2):
    name1, name2 = n1.get('name'), n2.get('name')
    return any([name1 == name2,  # same amino acid
                name1.startswith('X') and name2.startswith('X'),  # both unknown
                ':' in name1 and ':' in name2])  # both lipid tails

def edge_match(e1, e2):
    return e1.get('bond_type') == e2.get('bond_type')

def deduplicate_rban_graphs(compacts: List[Dict]) -> Dict[str, str]:
    graphs = {c['compound_id']: compact_dict_to_graph(c) for c in compacts}
    compound_id_to_representative: Dict[str, str] = {}

    for i, (compound_id, graph) in enumerate(graphs.items()):
        print(f'{i}/{len(graphs)} Processing compound {compound_id}...')
        representative_id = next(
            (representative_id   # find first isomorphic graph from already found representatives
             for representative_id in compound_id_to_representative.values()
             if is_isomorphic(graph, graphs[representative_id],
                              node_match=node_match,
                              edge_match=edge_match)
             ),
            compound_id  # if none found, use itself as representative
        )

        compound_id_to_representative[compound_id] = representative_id

    return compound_id_to_representative


def main():
    nerpa_dir = Path(__file__).parent.parent
    monomers_cfg_file = nerpa_dir / 'configs/monomers_config.yaml'
    monomer_names_helper = load_monomer_names_helper(monomers_cfg_file, nerpa_dir)

    # input
    preprocessed_dir = nerpa_dir / 'data/input/preprocessed'
    print(f'Loading preprocessed data from {preprocessed_dir}...')
    rban_graphs = yaml.safe_load((preprocessed_dir / 'pnrpdb2_preprocessed_rban_graphs.yaml').read_text())
    print(f'Loaded {len(rban_graphs)} rBAN graphs')
    nrp_variants_dicts = yaml.safe_load((preprocessed_dir / 'pnrpdb2_preprocessed.yaml').read_text())
    nrp_variants = [NRP_Variant.from_yaml_dict(nrp_variant_dict)
                    for nrp_variant_dict in nrp_variants_dicts]
    print(f'Loaded {len(nrp_variants)} NRP variants')
    #
    # output
    output_table_path = nerpa_dir / 'data/for_training_and_testing/pnrpdb2_info.tsv'
    #

    def is_rban_recognized_monomer(monomer_name: str) -> bool:
        return not any([monomer_name.startswith('X'),
                       ':' in monomer_name,])

    def is_nrps_monomer(monomer_name: str) -> bool:
        parsed_name = monomer_names_helper.parsed_name(monomer_name,
                                                       name_format='rBAN/Norine')
        return (is_rban_recognized_monomer(monomer_name) and
                parsed_name.residue not in (NOT_NRPS_RESIDUE, PKS_RESIDUE))

    print('Clustering isomorphic NRP variants...')
    nrp_representatives = cluster_isomorphic_nrp_variants(nrp_variants)
    print(f'Found {len(set(nrp_representatives.values()))} clusters of isomorphic NRP variants')
    print('Deduplicating rBAN graphs...')
    rban_id_to_repr = deduplicate_rban_graphs(rban_graphs)
    print(f'Found {len(set(rban_id_to_repr.values()))} clusters of isomorphic rBAN graphs')

    # prepare table
    pnrpdb_info_rows = []
    for record in rban_graphs:
        compound_id = record['compound_id']
        monomers = record['nodes'].values()

        nrp_variant_id = NRP_Variant_ID(nrp_id=compound_id, variant_idx=0)
        if nrp_variant_id in nrp_representatives:
            nrp_variant_representative = nrp_representatives[nrp_variant_id].nrp_id
        else:
            nrp_variant_representative = None

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
                rban_id_to_repr[compound_id],
        })

    pnrpdb_info = pl.DataFrame(pnrpdb_info_rows)

    print(f'Writing table with info about PNRPDB2 NRP variants to {output_table_path}')
    output_table_path.write_text(pnrpdb_info.write_csv(separator='\t'))
    print('Done')


if __name__ == "__main__":
    main()
