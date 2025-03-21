from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM
from typing import Dict, List, Tuple
from src.matching.genes_fragments_intervals import get_genes_intervals, get_fragments_intervals
from src.matching.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType
)
from itertools import pairwise
from src.matching.alignment_type import Alignment
from src.matching.alignment_step_type import (
    AlignmentStep,
    AlignmentStep_BGC_Module_Info,
)
from src.rban_parsing.rban_monomer import rBAN_Monomer
from pathlib import Path
from itertools import count


# "long" edges are those that involve multiple modules
LONG_EDGE_TYPES = {
    DetailedHMMEdgeType.SKIP_GENE,
    DetailedHMMEdgeType.SKIP_FRAGMENT,
    DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START,
    DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END,
}

def get_long_edges_ranges(hmm: DetailedHMM, path: List[int]) -> Dict[Tuple[int, int], Tuple[int, int]]:
    # returns (edge_from, edge_to) -> (fst_module_idx, last_module_idx)
    genes_intervals = get_genes_intervals(hmm.bgc_variant.modules)
    fragments_intervals = get_fragments_intervals(hmm.bgc_variant.modules)

    edges_ranges = {}
    fst_module_idx, last_module_idx = -1, -1
    for edge_from, edge_to in pairwise(path):
        edge = hmm.transitions[edge_from][edge_to]
        if edge_from in hmm.state_idx_to_module_idx:
            last_module_idx = hmm.state_idx_to_module_idx[edge_from] - 1
        if edge.edge_type not in LONG_EDGE_TYPES:
            continue

        fst_module_idx = last_module_idx + 1
        gene_id = hmm.bgc_variant.modules[fst_module_idx].gene_id
        fragment_idx = hmm.bgc_variant.modules[fst_module_idx].fragment_idx
        if edge.edge_type == DetailedHMMEdgeType.SKIP_GENE:
            last_module_idx = genes_intervals[gene_id][1]
        else:
            last_module_idx = fragments_intervals[fragment_idx][1]

        edges_ranges[(edge_from, edge_to)] = (fst_module_idx, last_module_idx)

    return edges_ranges


# NOTE: I didn't bother to optimize it in any way but maybe that's ok
def hmm_path_to_alignment(hmm: DetailedHMM,
                          path: List[int],
                          nrp_monomers: List[rBAN_Monomer]) -> Alignment:
    #free_idx = next(i for i in count()
    #                if not Path(f'{hmm.bgc_variant.genome_id}_path_{i}.png').exists())
    # hmm.draw(Path(f'{hmm.bgc_variant.genome_id}_path_{free_idx}.png'), path)  # for debugging
    #hmm.draw(Path(f'hmm_path.png'), path)  # for debugging
    #print('bgc: ', hmm.bgc_variant.genome_id)
    #print('path: ', path)
    #print('nrp monomers: ', [mon.residue for mon in nrp_monomers])

    long_edges_ranges = get_long_edges_ranges(hmm, path)
    alignment = []
    mon_idx = 0
    for i, (edge_from, edge_to) in enumerate(pairwise(path)):
        edge = hmm.transitions[edge_from][edge_to]

        # these edge types are handled separately because they involve multiple modules
        # and hence yield multiple alignment steps
        if edge.edge_type in LONG_EDGE_TYPES:
            fst_module_idx, last_module_idx = long_edges_ranges[(edge_from, edge_to)]
            alignment.extend(AlignmentStep(
                bgc_module=AlignmentStep_BGC_Module_Info.from_bgc_module(hmm.bgc_variant.modules[module_idx]),
                nrp_monomer=None,
                score=edge.weight,
                match_detailed_score=None,
                step_type=edge.edge_type)
                             for module_idx in range(fst_module_idx, last_module_idx + 1))
            continue

        bgc_module = None
        nrp_monomer = None
        match_detailed_score = None
        transition_score = edge.weight
        emission_score = 0

        match edge.edge_type:
            # with emissions
            case DetailedHMMEdgeType.MATCH:
                bgc_module = hmm.bgc_variant.modules[hmm.state_idx_to_module_idx[edge_from]]
                nrp_monomer = nrp_monomers[mon_idx]
                match_detailed_score = hmm.hmm_helper.match_detailed_score(bgc_module,
                                                                           nrp_monomers[mon_idx].to_base_mon(),
                                                                           hmm.bgc_variant.has_pks_domains())
                emission_score = sum(match_detailed_score)
                mon_idx += 1
            case DetailedHMMEdgeType.START_INSERTING | DetailedHMMEdgeType.INSERT:
                # bgc_module = hmm.bgc_variant.modules[hmm.state_idx_to_module_idx[edge_from]]
                nrp_monomer = nrp_monomers[mon_idx]
                # in the current setup, emission_score is always 0 for insertions, as if the emission is always the same
                mon_idx += 1
            case DetailedHMMEdgeType.START_INSERTING_AT_START | DetailedHMMEdgeType.INSERT_AT_START:
                #next_module_state_idx = next(next_vertex
                #                             for next_vertex in hmm.adj_list[edge_to]
                #                             if hmm.states[next_vertex].state_type == DetailedHMMStateType.MODULE_START)
                #bgc_module = hmm.bgc_variant.modules[next_module_state_idx]
                # in the current setup, emission_score is always 0 for insertions, as if the emission is always the same
                nrp_monomer = nrp_monomers[mon_idx]
                mon_idx += 1
            # skips
            case DetailedHMMEdgeType.SKIP_MODULE:
                next_module_state_idx = next(next_vertex
                                             for next_vertex in hmm.transitions[edge_from]
                                             if hmm.states[next_vertex].state_type in (DetailedHMMStateType.MODULE_START,
                                                                                       DetailedHMMStateType.FINAL))
                module_idx = hmm.state_idx_to_module_idx[next_module_state_idx] - 1 \
                    if next_module_state_idx != edge_to else len(hmm.bgc_variant.modules) - 1
                bgc_module = hmm.bgc_variant.modules[module_idx]

        score = transition_score + emission_score
        bgc_module_info = AlignmentStep_BGC_Module_Info.from_bgc_module(bgc_module) \
            if bgc_module is not None else None
        alignment.append(AlignmentStep(bgc_module=bgc_module_info,
                                       nrp_monomer=nrp_monomer,
                                       score=score,
                                       match_detailed_score=match_detailed_score,
                                       step_type=edge.edge_type))

    return alignment