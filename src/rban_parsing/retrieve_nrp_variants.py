from functools import partial
from typing import (
    Dict,
    List,
    Optional
)

import networkx as nx

from src.rban_parsing.nrp_variant_types import (
    NRP_Fragment,
    NRP_Variant,
    NRP_Variant_ID
)
from src.monomer_names_helper import (
    MonomerNamesHelper,
    UNKNOWN_RESIDUE
)
from src.rban_parsing.rban_parser import (
    AtomId,
    MonomerIdx,
    MonomerInfo
)
from src.generic.graphs import (
    BackboneSequence,
    putative_backbones
)
from src.pipeline.logging.logger import NerpaLogger
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.config import rBAN_Processing_Config, load_config, load_monomer_names_helper
from src.rban_parsing.rban_monomer import rBAN_Monomer

NerpaMonomerGraph = nx.DiGraph


def build_monomer(mon_info: MonomerInfo,
                  idx: int,
                  monomer_names_helper: MonomerNamesHelper) -> rBAN_Monomer:
    parsed_name = monomer_names_helper.parsed_name(mon_info.name, name_format='rBAN/Norine')
    return rBAN_Monomer(residue=parsed_name.residue,
                        methylated=parsed_name.methylated,
                        chirality=mon_info.chirality,
                        is_pks_hybrid=mon_info.is_pks_hybrid,
                        rban_name=mon_info.name,
                        rban_idx=idx)

def _get_bond_type(edge_info: List[Dict[MonomerIdx, AtomId]],
                   rban_record: Parsed_rBAN_Record) -> str:
    atom1, atom2 = edge_info[0].values()  # get the first
    if (atom1, atom2) in rban_record.atomic_bonds:
        return rban_record.atomic_bonds[(atom1, atom2)].bond_type
    elif (atom2, atom1) in rban_record.atomic_bonds:
        return rban_record.atomic_bonds[(atom2, atom1)].bond_type
    else:
        raise ValueError(f'No atomic bond found for monomer bond between atoms {atom1} and {atom2}')

def build_nx_graph(rban_record: Parsed_rBAN_Record,
                   backbone_bonds: List[str]) -> NerpaMonomerGraph:
    lipid_monomers = {monomer_idx for monomer_idx, monomer_info in rban_record.monomers.items()
                      if ':' in monomer_info.name}

    graph = nx.DiGraph()
    graph.add_nodes_from((monomer_idx, {'data': monomer_info})
                         for monomer_idx, monomer_info in rban_record.monomers.items()
                         if monomer_idx not in lipid_monomers)

    for (start, end), edge_info in rban_record.monomer_bonds.items():
        bond_type = _get_bond_type(edge_info, rban_record)
        if all([start not in lipid_monomers,
                end not in lipid_monomers,
                bond_type in backbone_bonds]):

            if all([bond_type == 'AMINO',
                    rban_record.atoms[start].name == 'N',
                    rban_record.atoms[end].name == 'C']):
                start, end = end, start  # ensure the direction of amino bond to go from C to N

            graph.add_edge(start, end, bond_type=bond_type)
            if bond_type != 'AMINO':
                graph.add_edge(end, start, bond_type=bond_type)  # non-amino bonds are bidirectional

    return graph


def backbone_sequence_to_fragment(backbone_sequence: BackboneSequence,
                                  G: NerpaMonomerGraph,
                                  monomer_names_helper: MonomerNamesHelper) -> NRP_Fragment:
    return NRP_Fragment(is_cyclic=backbone_sequence.is_cyclic,
                        monomers=[build_monomer(mon_info=G.nodes[idx]['data'],
                                                idx=idx,
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
                     for backbone in putative_backbones(graph,
                                                        compound_id=rban_record.compound_id,)]
    proper_fragments = [fragment
                        for fragment in all_fragments
                        if len(fragment.monomers) > 1 or fragment.monomers[0].residue != UNKNOWN_RESIDUE]
    isolated_unknown_monomers = [fragment.monomers[0]
                                 for fragment in all_fragments
                                 if len(fragment.monomers) == 1
                                 and fragment.monomers[0].residue not in monomer_names_helper.supported_residues]

    return NRP_Variant(nrp_variant_id=NRP_Variant_ID(variant_idx=0,
                                                     nrp_id=rban_record.compound_id),
                       fragments=proper_fragments,
                       isolated_unknown_monomers=isolated_unknown_monomers,
                       metadata=rban_record.metadata)


def rban_records_to_nrp_variants(rban_records: List[Parsed_rBAN_Record],
                                 monomer_names_helper: Optional[MonomerNamesHelper] = None,
                                 config: Optional[rBAN_Processing_Config] = None,
                                 log: Optional[NerpaLogger] = None) -> List[NRP_Variant]:
    # For usage outside the main pipeline -- load default configs if not provided as arguments
    if config is None or monomer_names_helper is None:
        main_config = load_config()
        if config is None:
            config = main_config.rban_processing_config
        if monomer_names_helper is None:
            monomer_names_helper = load_monomer_names_helper(monomers_cfg_file=main_config.monomers_config,
                                                             nerpa_dir=main_config.nerpa_dir,)

    nrp_variants = []
    for rban_record in rban_records:
        #log.info(f'Processing {rban_record.compound_id}')
        new_variant = process_single_record(rban_record, config.PNP_BONDS,
                                            monomer_names_helper=monomer_names_helper,
                                            min_recognized_nodes=config.MIN_RECOGNIZED_NODES)
        if new_variant.fragments:
            nrp_variants.append(new_variant)
        else:
            if log is not None:
                log.warning(f'Structure "{rban_record.compound_id}": unable to determine backbone sequence. '
                            f'Skipping "{rban_record.compound_id}".')

    return nrp_variants
