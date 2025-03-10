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
from itertools import chain
from bisect import bisect_left


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

    retained_matches = minimize_crossings_binary(all_matches)

    return retained_matches


def minimize_crossings_binary(alignment: list[tuple[int, int]]) -> list[tuple[int, int]]:
    lis_indices = longest_increasing_subsequence_binary(alignment)

    retained = [alignment[i] for i in lis_indices]

    return retained


def longest_increasing_subsequence_binary(pairs: list[tuple[int, int]]) -> list[int]:
    if not pairs:
        return []

    indexed_pairs = sorted(enumerate(pairs), key=lambda x: (x[1][0], x[1][1]))

    pairs_sorted = [pair for _, pair in indexed_pairs]
    original_indices = [idx for idx, _ in indexed_pairs]

    lis = []
    predecessors = [-1] * len(pairs_sorted)

    used_i = set()
    used_j = set()

    for current_index, (i, j) in enumerate(pairs_sorted):
        if i in used_i or j in used_j:
            continue

        pos = bisect_left([pairs_sorted[lis[k]][1] for k in range(len(lis))], j) if lis else 0

        if pos < len(lis):
            lis[pos] = current_index
        else:
            lis.append(current_index)

        if pos > 0:
            predecessors[current_index] = lis[pos - 1]

        used_i.add(i)
        used_j.add(j)

    lis_result = []
    last_index = lis[-1] if lis else -1
    while last_index != -1:
        lis_result.append(original_indices[last_index])
        last_index = predecessors[last_index]

    return lis_result[::-1]


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