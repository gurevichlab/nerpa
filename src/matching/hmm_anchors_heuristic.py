from typing import Dict, List, TYPE_CHECKING
if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM
from src.data_types import NRP_Monomer
from src.matching.match_type import Match
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.matching.hmm_to_alignment import hmm_path_to_alignment
from src.matching.alignment_type import show_alignment
from pathlib import Path


def nrp_monomer_matches_prediction(nrp_monomer: rBAN_Monomer,
                                   predicted_emissions: Dict[NRP_Monomer, float]) -> bool:
    predicted_residues = []
    for monomer, score in sorted(predicted_emissions.items(), key=lambda x: x[1], reverse=True):
        if monomer.residue not in predicted_residues:
            predicted_residues.append(monomer.residue)

    return nrp_monomer.residue == predicted_residues[0]


def heuristic_opt_path(hmm, # type: DetailedHMM, # can't import DetailedHMM here because of circular imports
                       bgc_predictions: List[Dict[NRP_Monomer, float]],
                       nrp_monomers: List[rBAN_Monomer]) -> List[int]:
    bgc_module_idx_to_match_state_idx = {module_idx: hmm._module_idx_to_state_idx[module_idx] + 2
                                         for module_idx in range(len(hmm.bgc_variant.modules))}

    path_with_emissions = hmm.get_opt_path_with_emissions(start=hmm.start_state_idx,
                                                          finish=hmm.final_state_idx,
                                                          emitted_monomers=nrp_monomers)
    opt_path = [state for state, _ in path_with_emissions]

    '''
    hmm.draw(Path('hmm.png'), opt_path)  # for debugging
    with open('alignment.txt', 'w') as f:  # for debugging
        f.write(show_alignment(hmm_path_to_alignment(hmm, opt_path, nrp_monomers)))
    '''

    return opt_path