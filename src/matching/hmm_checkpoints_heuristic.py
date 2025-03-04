from typing import (
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


def nrp_monomer_matches_prediction(nrp_monomer: rBAN_Monomer,
                                   predicted_emissions: Dict[NRP_Monomer, float]) -> bool:
    predicted_residues = []
    for monomer, score in sorted(predicted_emissions.items(), key=lambda x: x[1], reverse=True):
        if monomer.residue not in predicted_residues:
            predicted_residues.append(monomer.residue)

    return nrp_monomer.residue == predicted_residues[0]


def get_checkpoints(hmm, # type: DetailedHMM, # can't import DetailedHMM here because of circular imports
                    nrp_monomers: List[rBAN_Monomer]) -> List[Tuple[StateIdx, int]]:

    return [(hmm.start_state_idx, 0), (hmm.finish_state_idx, len(nrp_monomers))]