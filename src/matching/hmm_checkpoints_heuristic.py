from typing import (
    Callable,
    Dict,
    List,
    TYPE_CHECKING,
    Tuple
)
if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM
from src.data_types import NRP_Monomer
from src.matching.hmm_auxiliary_types import StateIdx
from src.rban_parsing.rban_monomer import rBAN_Monomer
from itertools import chain


def nrp_monomer_matches_prediction(nrp_monomer: rBAN_Monomer,
                                   predicted_emissions: Dict[NRP_Monomer, float]) -> bool:
    predicted_residues = []
    for monomer, score in sorted(predicted_emissions.items(), key=lambda x: x[1], reverse=True):
        if monomer.residue not in predicted_residues:
            predicted_residues.append(monomer.residue)

    return nrp_monomer.residue == predicted_residues[0]


def _get_checkpoints(len_bgc_modules: int,
                     len_nrp_monomers: int,
                     are_equal: Callable[[int, int], bool]) -> List[Tuple[int, int]]:  # (i, j) means i-th module corresponds to j-th monomer
    return []


def get_checkpoints(hmm: 'DetailedHMM',
                    nrp_monomers: List[rBAN_Monomer]) -> List[Tuple[StateIdx, int]]:
    bgc_modules_predictions = [hmm.hmm_helper.get_emissions(module, hmm.bgc_variant.has_pks_domains())
                               for module in hmm.bgc_variant.modules]

    def are_equal(module_idx: int,
                    monomer_idx: int) -> bool:
        return nrp_monomer_matches_prediction(nrp_monomers[monomer_idx],
                                              bgc_modules_predictions[module_idx])

    matched_pairs = _get_checkpoints(len(hmm.bgc_variant.modules),
                                     len(nrp_monomers),
                                     are_equal)

    inner_checkpoints = [(hmm._module_idx_to_match_state_idx[module_idx], monomer_idx)
                         for module_idx, monomer_idx in matched_pairs]
    return list(chain([(hmm.start_state_idx, 0)],
                      inner_checkpoints,
                      [(hmm.final_state_idx, len(nrp_monomers))]))