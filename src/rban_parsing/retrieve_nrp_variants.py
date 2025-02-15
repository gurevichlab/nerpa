from functools import partial
from itertools import (
    chain,
    count
)
from typing import (
    List
)

import networkx as nx

from src.data_types import (
    NRP_Fragment,
    NRP_Variant,
)
from src.monomer_names_helper import (
    MonomerNamesHelper,
    NRP_Monomer,
    UNKNOWN_RESIDUE
)
from src.rban_parsing.rban_parser import MonomerInfo
from src.rban_parsing.monomer_features import get_monomer_features
from src.generic.graphs import (
    BackboneSequence,
    permuted_backbones,
    putative_backbones
)
from src.pipeline.logger import NerpaLogger
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.config import rBAN_Processing_Config
from src.rban_parsing.rban_monomer import rBAN_Monomer

NerpaMonomerGraph = nx.DiGraph


def build_monomer(mon_info: MonomerInfo,
                  idx: int,
                  backbone_sequence: BackboneSequence,
                  monomer_names_helper: MonomerNamesHelper) -> rBAN_Monomer:
    parsed_name = monomer_names_helper.parsed_name(mon_info.name, name_format='norine')
    monomer_features = get_monomer_features(monomer_idx=idx,
                                            monomer_info=mon_info,
                                            backbone_sequence=backbone_sequence,
                                            names_helper=monomer_names_helper)
    return rBAN_Monomer(residue=parsed_name.residue,
                        methylated=parsed_name.methylated,
                        chirality=mon_info.chirality,
                        is_pks_hybrid=mon_info.is_pks_hybrid,
                        rban_name=mon_info.name,
                        rban_idx=idx)


def build_nx_graph(rban_record: Parsed_rBAN_Record,
                   backbone_bonds: List[str]) -> NerpaMonomerGraph:
    lipid_monomers = {monomer_idx for monomer_idx, monomer_info in rban_record.monomers.items()
                      if ':' in monomer_info.name}

    graph = nx.DiGraph()
    graph.add_nodes_from((monomer_idx, {'data': monomer_info})
                         for monomer_idx, monomer_info in rban_record.monomers.items()
                         if monomer_idx not in lipid_monomers)

    graph.add_edges_from((end, start, {'data': edge_info})
                         for (start, end), edge_info in rban_record.monomer_bonds.items()
                         if all([start not in lipid_monomers,
                                 end not in lipid_monomers,
                                 edge_info.bondType in backbone_bonds]))

    return graph


def backbone_sequence_to_fragment(backbone_sequence: BackboneSequence,
                                  G: NerpaMonomerGraph,
                                  monomer_names_helper: MonomerNamesHelper) -> NRP_Fragment:
    return NRP_Fragment(is_cyclic=backbone_sequence.is_cyclic,
                        monomers=[build_monomer(mon_info=G.nodes[idx]['data'],
                                                idx=idx,
                                                backbone_sequence=backbone_sequence,
                                                monomer_names_helper=monomer_names_helper)
                                  for idx in backbone_sequence.node_idxs])


def sufficiently_covered(sequence: BackboneSequence,
                         G: NerpaMonomerGraph,
                         min_recognized_nodes=2) -> bool:
    return sum(G.nodes[node_idx]['data'].residue != UNKNOWN_RESIDUE
               for node_idx in sequence.node_idxs) >= min_recognized_nodes


def trim_unrecognized_monomers(backbone: BackboneSequence, G: NerpaMonomerGraph) -> BackboneSequence:
    if backbone.is_cyclic:
        return backbone

    leftmost_monomer = G.nodes[backbone.node_idxs[0]]['data']
    rightmost_monomer = G.nodes[backbone.node_idxs[-1]]['data']
    left_bound = 1 if leftmost_monomer.residue == UNKNOWN_RESIDUE else 0
    right_bound = len(backbone.node_idxs) - 1 if rightmost_monomer.residue == UNKNOWN_RESIDUE \
        else len(backbone.node_idxs)
    return BackboneSequence(node_idxs=backbone.node_idxs[left_bound:right_bound],
                            is_cyclic=False)


def process_single_record(rban_record: Parsed_rBAN_Record,
                          backbone_bond_types: List[str],
                          monomer_names_helper: MonomerNamesHelper,
                          min_recognized_nodes=2) -> NRP_Variant:
    graph = build_nx_graph(rban_record, backbone_bond_types)

    backbone_to_fragment = partial(backbone_sequence_to_fragment, G=graph, monomer_names_helper=monomer_names_helper)

    all_fragments = [backbone_to_fragment(backbone)
                     for backbone in putative_backbones(graph)]
    proper_fragments = [fragment
                        for fragment in all_fragments
                        if len(fragment.monomers) > 1 or fragment.monomers[0].residue != UNKNOWN_RESIDUE]
    isolated_unknown_monomers = [fragment.monomers[0]
                                 for fragment in all_fragments
                                 if len(fragment.monomers) == 1 and fragment.monomers[0].residue == UNKNOWN_RESIDUE]

    return NRP_Variant(variant_idx=0,
                       nrp_id=rban_record.compound_id,
                       fragments=proper_fragments,
                       isolated_unknown_monomers=isolated_unknown_monomers)


def retrieve_nrp_variants(rban_records: List[Parsed_rBAN_Record],
                          monomer_names_helper: MonomerNamesHelper,
                          config: rBAN_Processing_Config,
                          log: NerpaLogger) -> List[NRP_Variant]:

    nrp_variants = []
    for rban_record in rban_records:
        #log.info(f'Processing {rban_record.compound_id}')
        new_variant = process_single_record(rban_record, config.PNP_BONDS,
                                            monomer_names_helper=monomer_names_helper,
                                            min_recognized_nodes=config.MIN_RECOGNIZED_NODES)
        if not new_variant.fragments:
            log.warning(f'Structure "{rban_record.compound_id}": unable to determine backbone sequence. '
                        f'Skipping "{rban_record.compound_id}".')
        nrp_variants.append(new_variant)

    return nrp_variants