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
    NRP_Monomer,
    NRP_Variant,
)
from src.monomer_names_helper import MonomerNamesHelper, UNKNOWN_RESIDUE
from src.rban_parsing.rban_parser import MonomerInfo
from src.generic.graph_manipulations import (
    BackboneSequence,
    permuted_backbones,
    putative_backbones
)
from src.pipeline.logger import NerpaLogger
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.config import rBAN_Processing_Config

NerpaMonomerGraph = nx.DiGraph


def build_monomer(mon_info: MonomerInfo,
                  idx: int,
                  monomer_names_helper: MonomerNamesHelper) -> NRP_Monomer:
    parsed_name = monomer_names_helper.parsed_name(mon_info.name, name_format='norine')

    return NRP_Monomer(residue=parsed_name.residue,
                       modifications=parsed_name.modifications,
                       chirality=mon_info.chirality,
                       rban_name=mon_info.name,
                       rban_idx=idx,
                       is_hybrid=mon_info.is_hybrid)


def build_nx_graph(rban_record: Parsed_rBAN_Record,
                   backbone_bonds: List[str],
                   monomer_names_helper: MonomerNamesHelper) -> NerpaMonomerGraph:
    lipid_monomers = {monomer_idx for monomer_idx, monomer_info in rban_record.monomers.items()
                      if ':' in monomer_info.name}

    graph = nx.DiGraph()
    graph.add_nodes_from((monomer_idx, {'data': build_monomer(monomer_info, monomer_idx, monomer_names_helper)})
                         for monomer_idx, monomer_info in rban_record.monomers.items()
                         if monomer_idx not in lipid_monomers)

    graph.add_edges_from((end, start, {'data': edge_info})
                         for (start, end), edge_info in rban_record.monomer_bonds.items()
                         if all([start not in lipid_monomers,
                                 end not in lipid_monomers,
                                 edge_info.bondType in backbone_bonds]))

    return graph


def backbone_sequence_to_fragment(backbone_sequence: BackboneSequence, G: NerpaMonomerGraph) -> NRP_Fragment:
    return NRP_Fragment(is_cyclic=backbone_sequence.is_cyclic,
                        monomers=[G.nodes[idx]['data'] for idx in backbone_sequence.node_idxs])


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
                          min_recognized_nodes=2) -> List[NRP_Variant]:
    graph = build_nx_graph(rban_record, backbone_bond_types, monomer_names_helper)

    backbone_to_fragment = partial(backbone_sequence_to_fragment, G=graph, monomer_names_helper=monomer_names_helper)

    # all fragments individually
    return [NRP_Variant(variant_idx=0,
                        nrp_id=rban_record.compound_id,
                        fragments=[backbone_to_fragment(backbone)
                                   for backbone in putative_backbones(graph)])]


def retrieve_nrp_variants(rban_records: List[Parsed_rBAN_Record],
                          monomer_names_helper: MonomerNamesHelper,
                          config: rBAN_Processing_Config,
                          log: NerpaLogger) -> List[NRP_Variant]:

    nrp_variants = []
    for rban_record in rban_records:
        log.info(f'Processing {rban_record.compound_id}')
        new_variants = process_single_record(rban_record, config.PNP_BONDS,
                                             monomer_names_helper=monomer_names_helper,
                                             min_recognized_nodes=2)
        if not new_variants:
            log.warning(f'Structure "{rban_record.compound_id}": unable to determine backbone sequence. '
                        f'Skipping "{rban_record.compound_id}".')
        nrp_variants.extend(new_variants)

    return nrp_variants