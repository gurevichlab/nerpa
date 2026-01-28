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

from scripts.pnrpdb_compound_similarity import nerpa_mon_cmp

def main():
    nerpa_dir = Path(__file__).resolve().parent
    nerpa_cfg = load_config()
    # preprocessed_dir = Path('/home/ilianolhin/git/nerpa2/nerpa_results/parsed_compounds/nerpa_results/preprocessed_input')
    # pnrpdb_nrp_variants = preprocessed_dir / 'NRP_variants.yaml'
    preprocessed_dir = nerpa_dir / 'data' / 'input' / 'preprocessed'
    pnrpdb_nrp_variants = preprocessed_dir / 'pnrpdb2_nrp_variants.yaml'

    print(f'Loading rBAN records and NRP variants from {preprocessed_dir}...')
    nrp_variants = [
        NRP_Variant.from_yaml_dict(nrp_dict)
        for nrp_dict in yaml.safe_load(open(pnrpdb_nrp_variants, 'r'))
        if nrp_dict['nrp_variant_id']['nrp_id'] in ('NPA028802', 'NPA028880')
    ]

    nerpa_graphs_by_id = {
        variant.nrp_variant_id.nrp_id: variant.to_nx_digraph()
        for variant in nrp_variants
    }

    for (nrp1_id, graph1), (nrp2_id, graph2) in combinations(nerpa_graphs_by_id.items(), 2):
        are_isomorphic = is_isomorphic(graph1, graph2,
                                       node_match=lambda n1, n2: nerpa_mon_cmp(n1['monomer'], n2['monomer']))
        if are_isomorphic:
            print(f'Graphs for {nrp1_id} and {nrp2_id} are isomorphic!')

if __name__ == '__main__':
    main()
