import csv
from itertools import combinations
from pathlib import Path

import yaml

from src.generic.graphs import graphs_one_substitution_away
from src.monomer_names_helper import NRP_Monomer, Chirality
from src.rban_parsing.nrp_variant_types import NRP_Variant

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

def main():
    nerpa_dir = Path(__file__).resolve().parent.parent
    pnrpdb_preprocessed_path = (
            nerpa_dir
            / 'input'
            / 'preprocessed'
            / 'pnrpdb2_preprocessed_nrp_variants.yaml'
    )

    nrp_variants = [
        NRP_Variant.from_yaml_dict(nrp_variant_dict)
        for nrp_variant_dict in yaml.safe_load(open(pnrpdb_preprocessed_path, 'r'))
    ]

    nrp_variants = sorted(nrp_variants, key=lambda nrp: nrp.nrp_variant_id.nrp_id)
    nrp_variants_by_id = {nrp.nrp_variant_id.nrp_id: nrp for nrp in nrp_variants}
    nrp_nx_graphs_by_id = {nrp.nrp_variant_id.nrp_id: nrp.to_nx_digraph(node_label_key='monomer')
                           for nrp in nrp_variants}

    rows = []
    for nrp1_id, nrp2_id in combinations(nrp_variants_by_id.keys(), 2):
        nrp1_variant, nrp1_graph = nrp_variants_by_id[nrp1_id], nrp_nx_graphs_by_id[nrp1_id]
        nrp2_variant, nrp2_graph = nrp_variants_by_id[nrp2_id], nrp_nx_graphs_by_id[nrp2_id]

        row = {}
        for (mon_cmp_name, _mon_cmp) in [('residue_only_cmp', mon_cmp_only_residues),
                                         ('unknown_chr_equal_known_cmp', mon_cmp_unknown_chr_equal_known)]:
            row['isomorphic_' + mon_cmp_name] = nrp1_variant.is_isomorphic_to(nrp2_variant,
                                                                              monomers_comparator=_mon_cmp)
            row['one_sub_away_' + mon_cmp_name] = graphs_one_substitution_away(nrp1_graph,
                                                                               nrp2_graph,
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