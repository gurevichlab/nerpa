import csv
import os
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
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.rban_parsing.retrieve_nrp_variants import build_monomer
from functools import partial
from itertools import combinations, islice
from joblib import Parallel, delayed


def rban_mon_cmp(mon1: rBAN_Monomer, mon2: rBAN_Monomer) -> bool:
    return all([mon1.rban_name == mon2.rban_name,
                mon1.methylated == mon2.methylated,
                mon1.chirality == mon2.chirality,
                mon1.is_pks_hybrid == mon2.is_pks_hybrid])

def mon_cmp_only_residues(mon1: rBAN_Monomer, mon2: rBAN_Monomer) -> bool:
    return mon1.residue == mon2.residue

def mon_cmp_unknown_chr_equal_known(mon1: rBAN_Monomer, mon2: rBAN_Monomer) -> bool:
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


def compare_pair(
    nrp1_id: str,
    nrp2_id: str,
    rban_graphs_by_id,
    nerpa_graphs_by_id,
):
    g1_rban = rban_graphs_by_id[nrp1_id]
    g2_rban = rban_graphs_by_id[nrp2_id]
    g1_nerpa = nerpa_graphs_by_id[nrp1_id]
    g2_nerpa = nerpa_graphs_by_id[nrp2_id]

    row = {}

    for graph_type, (g1, g2) in [
        ("rban", (g1_rban, g2_rban)),
        ("nerpa", (g1_nerpa, g2_nerpa)),
    ]:
        for mon_cmp_name, mon_cmp in [
            ("rban_mon_cmp", rban_mon_cmp),
            ("unknown_chr_equal_known_cmp", mon_cmp_unknown_chr_equal_known),
        ]:
            nm = partial(_node_match, cmp=mon_cmp)

            row[f"{graph_type}_isomorphic_{mon_cmp_name}"] = is_isomorphic(
                g1, g2, node_match=nm
            )
            row[f"{graph_type}_one_sub_away_{mon_cmp_name}"] = graphs_one_substitution_away(
                g1, g2, nodes_comparator=mon_cmp, label_key="monomer"
            )

    if any(row.values()):
        row["nrp1_id"] = nrp1_id
        row["nrp2_id"] = nrp2_id
        return row

    return None


def batched(iterable, batch_size: int):
    it = iter(iterable)
    while True:
        batch = list(islice(it, batch_size))
        if not batch:
            return
        yield batch


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
        for record in yaml.load(open(pnrpdb_parsed_rban_records, 'r'),
                                Loader=IgnoreTagsLoader)
    ]
    nrp_variants = [
        NRP_Variant.from_yaml_dict(nrp_dict)
        for nrp_dict in yaml.safe_load(open(pnrpdb_nrp_variants, 'r'))
    ]
    # Keep the original if you still need it elsewhere; otherwise this is enough:
    nerpa_graphs_by_id = {
        variant.nrp_variant_id.nrp_id: variant.to_nx_digraph()
        for variant in nrp_variants
    }

    rban_graphs_by_id = {
        record.compound_id: parsed_record_to_graph(record, monomer_names_helper)
        for record in rban_records
    }

    print(f'Comparing {len(rban_graphs_by_id)} compounds...')
    ids = list(rban_graphs_by_id.keys())
    num_pairs = len(ids) * (len(ids) - 1) // 2

    output_path = (
            nerpa_dir
            / "data"
            / "for_training_and_testing"
            / "pnrpdb2_compound_similarity.tsv"
    )

    # Stable header (don’t depend on rows[0])
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

    pairs_iter = combinations(ids, 2)

    # Tune this: larger batches reduce joblib overhead; smaller batches update progress more often.
    batch_size = 10000

    done = 0
    with open(output_path, "w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for batch in batched(pairs_iter, batch_size):
            # Progress (batch-level)
            print(f"{done}/{num_pairs} pairs processed...")

            results = Parallel(
                n_jobs=os.cpu_count() // 2,
                backend="multiprocessing",  # fork on Linux -> cheap sharing of big dicts/graphs
                verbose=0,
            )(
                delayed(compare_pair)(a, b, rban_graphs_by_id, nerpa_graphs_by_id)
                for (a, b) in batch
            )

            # Write only positives
            for row in results:
                if row is not None:
                    writer.writerow(row)

            done += len(batch)

    print(f"Wrote comparison results to {output_path}")


if __name__ == '__main__':
    main()