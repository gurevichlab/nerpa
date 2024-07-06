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
    rBAN_Residue_Name
)
from src.generic.graph_manipulations import (
    BackboneSequence,
    permuted_backbones,
    putative_backbones
)
from src.pipeline.logger import NerpaLogger
from src.rban_parsing.rban_names_helper import rBAN_Names_Helper
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.config import rBAN_Processing_Config

NerpaMonomerGraph = nx.DiGraph


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


def build_monomer(G: NerpaMonomerGraph, mon_idx: int,
                  names_helper: rBAN_Names_Helper) -> NRP_Monomer:
    mon_info = G.nodes[mon_idx]['data']
    parsed_name = names_helper.parsed_rban_name[mon_info.name]

    return NRP_Monomer(residue=parsed_name.residue,
                       modifications=parsed_name.modifications,
                       chirality=mon_info.chirality,
                       rban_name=mon_info.name,
                       rban_idx=mon_idx)


def backbone_sequence_to_fragment(backbone_sequence: BackboneSequence, G: NerpaMonomerGraph,
                                  names_helper: rBAN_Names_Helper) -> NRP_Fragment:
    return NRP_Fragment(is_cyclic=backbone_sequence.is_cyclic,
                        monomers=[build_monomer(G, idx, names_helper=names_helper)
                                  for idx in backbone_sequence.node_idxs])


def sufficiently_covered(sequence: BackboneSequence,
                         G: NerpaMonomerGraph,
                         recognized_monomers: List[rBAN_Residue_Name],
                         min_recognized_nodes=2) -> bool:
    return sum(G.nodes[node_idx]['data'].name in recognized_monomers
               for node_idx in sequence.node_idxs) >= min_recognized_nodes


def trim_unrecognized_monomers(backbone: BackboneSequence, G: NerpaMonomerGraph) -> BackboneSequence:
    if backbone.is_cyclic:
        return backbone

    leftmost_monomer = G.nodes[backbone.node_idxs[0]]['data']
    rightmost_monomer = G.nodes[backbone.node_idxs[-1]]['data']
    left_bound = 1 if leftmost_monomer.name.startswith('X') else 0
    right_bound = len(backbone.node_idxs) - 1 if rightmost_monomer.name.startswith('X') \
        else len(backbone.node_idxs)
    return BackboneSequence(node_idxs=backbone.node_idxs[left_bound:right_bound],
                            is_cyclic=False)


def process_single_record(rban_record: Parsed_rBAN_Record,
                          recognized_monomers: List[rBAN_Residue_Name],
                          backbone_bond_types: List[str],
                          rban_names_helper: rBAN_Names_Helper,
                          min_recognized_nodes=2) -> List[NRP_Variant]:
    graph = build_nx_graph(rban_record, backbone_bond_types)

    # Split the graph into paths and simple cycles
    backbones = [trim_unrecognized_monomers(backbone, graph)
                 for backbone in putative_backbones(graph, min_nodes=2)
                 if sufficiently_covered(backbone, graph, recognized_monomers, min_recognized_nodes)]

    perm_backbones = permuted_backbones(backbones)

    variant_idx_counter = count()
    backbone_to_fragment = partial(backbone_sequence_to_fragment, G=graph, names_helper=rban_names_helper)

    # all fragments individually
    nrp_variants = [NRP_Variant(variant_idx=next(variant_idx_counter),
                                nrp_id=rban_record.compound_id,
                                fragments=[backbone_to_fragment(backbone)])
                    for backbone in chain(backbones, perm_backbones)]

    # all fragments together
    if len(backbones) > 1:
        nrp_variants.append(NRP_Variant(variant_idx=next(variant_idx_counter),
                                        nrp_id=rban_record.compound_id,
                                        fragments=list(map(backbone_to_fragment, backbones))))

    return nrp_variants


def retrieve_nrp_variants(rban_records: List[Parsed_rBAN_Record],
                          recognized_monomers: List[rBAN_Residue_Name],
                          names_helper: rBAN_Names_Helper,
                          config: rBAN_Processing_Config,
                          log: NerpaLogger) -> List[NRP_Variant]:

    nrp_variants = []
    for rban_record in rban_records:
        log.info(f'Processing {rban_record.compound_id}')
        new_variants = process_single_record(rban_record, recognized_monomers, config.PNP_BONDS,
                                             names_helper, min_recognized_nodes=2)
        if not new_variants:
            log.warning(f'Structure "{rban_record.compound_id}": unable to determine backbone sequence. '
                        f'Skipping "{rban_record.compound_id}".')
        nrp_variants.extend(new_variants)

    return nrp_variants