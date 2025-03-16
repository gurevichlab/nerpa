from typing import (
    Callable,
    Dict,
    List,
    TYPE_CHECKING,
    Tuple
)
if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM
from src.data_types import BGC_Module, NRP_Monomer
from src.matching.hmm_auxiliary_types import StateIdx
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.generic.combinatorics import longest_increasing_subsequence_without_collisions
from itertools import chain


def bgc_module_matches_nrp_monomer(bgc_module_emissions: Dict[NRP_Monomer, float],
                                   nrp_monomer: NRP_Monomer) -> bool:
    predicted_residues = []
    for monomer, score in sorted(bgc_module_emissions.items(), key=lambda x: x[1], reverse=True):
        if monomer.residue not in predicted_residues:
            predicted_residues.append(monomer.residue)

    return nrp_monomer.residue == predicted_residues[0]


def _get_checkpoints(len_bgc_modules: int,
                     len_nrp_monomers: int,
                     are_equal: Callable[[int, int], bool]) -> List[Tuple[int, int]]:  # (i, j) means i-th module corresponds to j-th monomer
    all_matches = []

    for i in range(len_bgc_modules - 1):
        for j in range(len_nrp_monomers - 1):
            if are_equal(i, j) and are_equal(i + 1, j + 1):
                all_matches.append((i, j))

    retained_matches = max_not_intersecting_edges(all_matches)

    return retained_matches


def max_not_intersecting_edges(edges: list[tuple[int, int]]) -> list[tuple[int, int]]:
    '''
    finds the maximum subset of not intersecting edges in a bipartile graph,
    edges (i1, j1), (i2, j2) are called not intersecting if i2 > i1 and j2 > j1 or vice versa
    '''
    lis_indices = longest_increasing_subsequence_without_collisions(edges)

    retained = [edges[i] for i in lis_indices]

    return retained


def get_checkpoints(hmm: 'DetailedHMM',
                    nrp_monomers: List[rBAN_Monomer]) -> List[Tuple[StateIdx, int]]:
    bgc_modules_predictions = [hmm.hmm_helper.get_emissions(module, hmm.bgc_variant.has_pks_domains())
                               for module in hmm.bgc_variant.modules]

    def are_equal(module_idx: int,
                  monomer_idx: int) -> bool:
        return bgc_module_matches_nrp_monomer(bgc_modules_predictions[module_idx],
                                              nrp_monomers[monomer_idx].to_base_mon())

    matched_pairs = _get_checkpoints(len(hmm.bgc_variant.modules),
                                     len(nrp_monomers),
                                     are_equal)

    inner_checkpoints = [(hmm._module_idx_to_match_state_idx[module_idx], monomer_idx)
                         for module_idx, monomer_idx in matched_pairs]
    return list(chain([(hmm.start_state_idx, 0)],
                      inner_checkpoints,
                      [(hmm.final_state_idx, len(nrp_monomers))]))