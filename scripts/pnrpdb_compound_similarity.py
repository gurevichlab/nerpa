import csv
import os
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import networkx as nx
import yaml
from networkx.algorithms.isomorphism.isomorph import is_isomorphic

from src.config import load_config, load_monomer_names_helper
from src.generic.graphs import graphs_one_substitution_away
from src.monomer_names_helper import NRP_Monomer, Chirality, MonomerNamesHelper
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.rban_parsing.retrieve_nrp_variants import build_monomer
from functools import partial
from itertools import combinations, islice
from joblib import Parallel, delayed

def rban_mon_cmp_key(mon: rBAN_Monomer) -> tuple:
    return (mon.rban_name,
            mon.methylated,
            mon.chirality,
            mon.is_pks_hybrid)

def rban_mon_cmp(mon1: rBAN_Monomer, mon2: rBAN_Monomer) -> bool:
    return all([mon1.rban_name == mon2.rban_name,
                mon1.methylated == mon2.methylated,
                mon1.chirality == mon2.chirality,
                mon1.is_pks_hybrid == mon2.is_pks_hybrid])

def nerpa_mon_cmp(mon1: rBAN_Monomer, mon2: rBAN_Monomer) -> bool:
    return all([mon1.residue == mon2.residue,
                mon1.methylated == mon2.methylated,
                mon1.chirality == mon2.chirality,
                mon1.is_pks_hybrid == mon2.is_pks_hybrid])

def mon_cmp_only_residues(mon1: rBAN_Monomer, mon2: rBAN_Monomer) -> bool:
    return mon1.residue == mon2.residue

def mon_cmp_wo_chr_key(mon: rBAN_Monomer) -> tuple:
    return (mon.residue,
            mon.methylated,
            mon.is_pks_hybrid)

def unknown_chr_equal_known_cmp(mon1: rBAN_Monomer, mon2: rBAN_Monomer) -> bool:
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
        mon = build_monomer(mon_info, mon_idx, monomer_names_helper)
        G.add_node(mon_idx, monomer=mon)

    for (u, v), edge_info in record.monomer_bonds.items():
        G.add_edge(u, v, bond_type=edge_info.bondType)

    return G

class IgnoreTagsLoader(yaml.SafeLoader):
    pass

def unknown_tag(loader, tag_suffix, node):
    if isinstance(node, yaml.ScalarNode):
        return loader.construct_scalar(node)
    elif isinstance(node, yaml.SequenceNode):
        return loader.construct_sequence(node)
    elif isinstance(node, yaml.MappingNode):
        return loader.construct_mapping(node)

IgnoreTagsLoader.add_multi_constructor("!", unknown_tag)
IgnoreTagsLoader.add_constructor(
    "tag:yaml.org,2002:python/tuple",
    lambda loader, node: loader.construct_sequence(node),
)


def _node_match(n1, n2, *, cmp):
    return cmp(n1["monomer"], n2["monomer"])


def batched(iterable, batch_size: int):
    it = iter(iterable)
    while True:
        batch = list(islice(it, batch_size))
        if not batch:
            return
        yield batch


def add_sim_info(graphs_by_id: Dict[str, nx.DiGraph],
                 sim_dict: Dict[Tuple[str, str], Dict[str, bool]],
                 pref: str):
    print(f'Adding similarity info for {len(graphs_by_id)} {pref} graphs...')
    graphs_by_key = defaultdict(list)
    for nrp_id, g in graphs_by_id.items():
        key = tuple(sorted(mon_cmp_wo_chr_key(g.nodes[node_idx]['monomer'])
                           for node_idx in g.nodes))
        graphs_by_key[key].append(nrp_id)

    unique_keys = len(graphs_by_key)
    for i, (key, ids) in enumerate(graphs_by_key.items()):
        print(f'Processing cluster {i}/{unique_keys}')
        for nrp1_id, nrp2_id in combinations(ids, 2):
            for mon_cmp in [nerpa_mon_cmp, unknown_chr_equal_known_cmp]:
                g1 = graphs_by_id[nrp1_id]
                g2 = graphs_by_id[nrp2_id]

                isomorphic = is_isomorphic(g1, g2, node_match=lambda m1, m2: mon_cmp(m1['monomer'], m2['monomer']))
                one_sub_away = graphs_one_substitution_away(
                    g1, g2, nodes_comparator=mon_cmp, label_key="monomer"
                )

                if (nrp1_id, nrp2_id) not in sim_dict:
                    sim_dict[(nrp1_id, nrp2_id)] = {}
                sim_dict[(nrp1_id, nrp2_id)][f"{pref}_equal_{mon_cmp.__name__}"] = isomorphic
                sim_dict[(nrp1_id, nrp2_id)][f"{pref}_one_sub_away_{mon_cmp.__name__}"] = one_sub_away


def main():
    nerpa_dir = Path(__file__).resolve().parent.parent
    nerpa_cfg = load_config()
    monomer_names_helper = load_monomer_names_helper(nerpa_cfg.monomers_config,
                                                     nerpa_dir)
    preprocessed_dir = nerpa_dir / 'data' / 'input' / 'preprocessed'
    pnrpdb_parsed_rban_records = preprocessed_dir / 'pnrpdb2_parsed_rban_records.yaml'
    pnrpdb_nrp_variants = preprocessed_dir / 'pnrpdb2_nrp_variants.yaml'

    print(f'Loading rBAN records and NRP variants from {preprocessed_dir}...')
    # rban_records = [
    #     Parsed_rBAN_Record.from_dict(record)
    #     for record in yaml.load(open(pnrpdb_parsed_rban_records, 'r'),
    #                             Loader=IgnoreTagsLoader)
    # ]
    nrp_variants = [
        NRP_Variant.from_yaml_dict(nrp_dict)
        for nrp_dict in yaml.safe_load(open(pnrpdb_nrp_variants, 'r'))
    ]
    # Keep the original if you still need it elsewhere; otherwise this is enough:
    nerpa_graphs_by_id = {
        variant.nrp_variant_id.nrp_id: variant.to_nx_digraph()
        for variant in nrp_variants
    }

    # rban_graphs_by_id = {
    #     record.compound_id: parsed_record_to_graph(record, monomer_names_helper)
    #     for record in rban_records
    #     if record.compound_id in nerpa_graphs_by_id
    # }
    #
    print(f'Comparing {len(nerpa_graphs_by_id)} compounds...')
    sim_dict = {}
    # add_sim_info(rban_graphs_by_id, sim_dict, pref="rban")
    add_sim_info(nerpa_graphs_by_id, sim_dict, pref="nerpa")

    rows = [
        {"nrp1_id": nrp1_id,
         "nrp2_id": nrp2_id,
         **sim_info}
        for (nrp1_id, nrp2_id), sim_info in sim_dict.items()
    ]


    output_path = (
            nerpa_dir
            / "data"
            / "for_training_and_testing"
            / "pnrpdb2_compound_similarity.tsv"
    )
    # Stable header (don’t depend on rows[0])
    '''
    fieldnames = [
        "nrp1_id",
        "nrp2_id",
        "rban_isomorphic_rban_mon_cmp",
        "rban_one_sub_away_rban_mon_cmp",
        "rban_isomorphic_unknown_chr_equal_known_cmp",
        "rban_one_sub_away_unknown_chr_equal_known_cmp",
        "nerpa_isomorphic_rban_mon_cmp",
        "nerpa_one_sub_away_rban_mon_cmp",
        "nerpa_isomorphic_unknown_chr_equal_known_cmp",
        "nerpa_one_sub_away_unknown_chr_equal_known_cmp",
    ]
    '''

    fieldnames = ['nrp1_id', 'nrp2_id'] + [k for k in rows[0].keys()
                                           if k not in ('nrp1_id', 'nrp2_id')]
    with open(output_path, "w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


    print(f"Wrote comparison results to {output_path}")


if __name__ == '__main__':
    main()