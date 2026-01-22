import csv
from itertools import combinations
from pathlib import Path
from typing import Callable

import networkx as nx
import yaml
from networkx.algorithms.isomorphism.isomorph import is_isomorphic

from src.config import load_config, load_monomer_names_helper
from src.generic.graphs import graphs_one_substitution_away
from src.monomer_names_helper import NRP_Monomer, Chirality, MonomerNamesHelper
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.rban_parsing.rban_parser import Parsed_rBAN_Record

def mon_cmp(mon1: NRP_Monomer, mon2: NRP_Monomer) -> bool:
    return all([mon1.residue == mon2.residue,
                mon1.methylated == mon2.methylated,
                mon1.chirality == mon2.chirality,
                mon1.is_pks_hybrid == mon2.is_pks_hybrid])

def mon_cmp_only_residues(mon1: NRP_Monomer, mon2: NRP_Monomer) -> bool:
    return mon1.residue == mon2.residue

def mon_cmp_unknown_chr_equal_known(mon1: NRP_Monomer, mon2: NRP_Monomer) -> bool:
    chr_equal = any([mon1.chirality == mon2.chirality,
                     mon1.chirality == Chirality.UNKNOWN,
                     mon2.chirality == Chirality.UNKNOWN])
    return all([mon1.residue == mon2.residue,
                mon1.methylated == mon2.methylated,
                chr_equal,
                mon1.is_pks_hybrid == mon2.is_pks_hybrid])

def parsed_record_to_graph(record: Parsed_rBAN_Record,
                           monomer_names_helper: MonomerNamesHelper) -> nx.DiGraph:
    G = nx.DiGraph()
    for mon_idx, mon_info in record.monomers.items():
        mon = monomer_names_helper.parsed_name(name=mon_info.name,
                                               name_format='rBAN/Norine')
        if mon_info.chirality != Chirality.UNKNOWN:
            mon = mon._replace(chirality=mon_info.chirality)
        if mon_info.is_pks_hybrid:
            mon = mon._replace(is_pks_hybrid=mon_info.is_pks_hybrid)

        G.add_node(mon.index, monomer=mon)

    for (u, v), edge_info in record.monomer_bonds.items():
        G.add_edge(u, v, bond_type=edge_info.bondType)

    return G


def main():
    nerpa_dir = Path(__file__).resolve().parent.parent
    nerpa_cfg = load_config()
    monomer_names_helper = load_monomer_names_helper(nerpa_cfg.monomers_config,
                                                     nerpa_dir)
    preprocessed_dir = nerpa_dir / 'data' / 'input' / 'preprocessed'
    pnrpdb_parsed_rban_records = preprocessed_dir / 'pnrpdb2_parsed_rban_records.yaml'
    pnrpdb_nrp_variants = preprocessed_dir / 'pnrpdb2_nrp_variants.yaml'

    print(f'Loading rBAN records and NRP variants from {preprocessed_dir}...')
    rban_records = [
        Parsed_rBAN_Record.from_dict(record)
        for record in yaml.safe_load(open(pnrpdb_parsed_rban_records, 'r'))
    ]
    nrp_variants = [
        NRP_Variant.from_yaml_dict(nrp_dict)
        for nrp_dict in yaml.safe_load(open(pnrpdb_nrp_variants, 'r'))
    ]

    rban_graphs_by_id = {
        record.compound_id: parsed_record_to_graph(record, monomer_names_helper)
        for record in rban_records
    }
    nrp_variants_by_id = {
        variant.nrp_variant_id.nrp_id: variant
        for variant in nrp_variants
    }


    num_pairs = len(rban_graphs_by_id) * (len(rban_graphs_by_id) - 1) // 2
    rows = []
    for i, (nrp1_id, nrp2_id) in enumerate(combinations(rban_graphs_by_id.keys(), 2)):
        if i % 1000 == 0:
            print(f'{i}/{num_pairs}. Comparing {nrp1_id} vs {nrp2_id}')

        nrp1_rban_graph = rban_graphs_by_id[nrp1_id]
        nrp2_rban_graph = rban_graphs_by_id[nrp2_id]
        nrp1_nerpa_graph = nrp_variants_by_id[nrp1_id]
        nrp2_nerpa_graph = nrp_variants_by_id[nrp2_id]

        row = {}
        for graph_type, (g1, g2) in [('rban', (nrp1_rban_graph, nrp2_rban_graph)),
                                     ('nerpa', (nrp1_nerpa_graph.to_nx_digraph(), nrp2_nerpa_graph.to_nx_digraph()))]:
            for (mon_cmp_name, _mon_cmp) in [('residue_only_cmp', mon_cmp_only_residues),
                                             ('unknown_chr_equal_known_cmp', mon_cmp_unknown_chr_equal_known)]:
                row[f'{graph_type}_isomorphic_{mon_cmp_name}'] = is_isomorphic(g1, g2,
                                                                               node_match=lambda n1, n2: _mon_cmp(n1['monomer'], n2['monomer']))
                row[f'{graph_type}_one_sub_away_{mon_cmp_name}'] = graphs_one_substitution_away(g1, g2,
                                                                                                nodes_comparator=_mon_cmp,
                                                                                                label_key='monomer')
        if any(row.values()):
            row['nrp1_id'] = nrp1_id
            row['nrp2_id'] = nrp2_id
            rows.append(row)

    output_path = (nerpa_dir
                   / 'data'
                   / 'for_training_and_testing'
                   / 'pnrpdb2_compound_similarity.tsv')

    with open(output_path, 'w') as f:
        csv_writer = csv.DictWriter(f, fieldnames=['nrp1_id', 'nrp2_id']
                                                  + [k for k in rows[0] if k not in ('nrp1_id', 'nrp2_id')])
        csv_writer.writeheader()
        csv_writer.writerows(rows)


if __name__ == '__main__':
    main()