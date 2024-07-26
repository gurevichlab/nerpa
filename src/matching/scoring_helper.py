import math

from typing import (
    List,
    Tuple,
    Optional
)
from src.data_types import (
    BGC_Module,
    BGC_Module_Modification,
    BGC_Fragment,
    NRP_Monomer,
    NRP_Fragment,
    LogProb,
    NRP_Monomer_Modification,
    UNKNOWN_RESIDUE
)
from src.matching.dp_types import DP_State
from src.matching.scoring_config import ScoringConfig, ModMatch, ChiralityMatch
from src.matching.alignment_types import AlignmentStepType
from dataclasses import dataclass
from enum import Enum, auto

FUNCTION_NAME_TO_STEP_TYPE = {'bgc_module_skip': AlignmentStepType.BGC_MODULE_SKIP,
                              'nrp_mon_skip': AlignmentStepType.NRP_MONOMER_SKIP,
                              'match': AlignmentStepType.MATCH,
                              'iterate_module': AlignmentStepType.ITERATE_MODULE,
                              'iterate_gene': AlignmentStepType.ITERATE_GENE}

@dataclass
class ScoringHelper:
    scoring_config: ScoringConfig
    bgc_modules: List[BGC_Module] = None
    nrp_monomers: List[NRP_Monomer] = None

    def set_bgc_modules(self, bgc_modules: List[BGC_Module]):
        self.bgc_modules = bgc_modules

    def set_nrp_monomers(self, nrp_monomers: List[NRP_Monomer]):
        self.nrp_monomers = nrp_monomers

    def bgc_module_skip(self, bgc_module_idx: int, dp_state: Optional[DP_State] = None) -> LogProb:  # argument is kept for better alignment backtracking
        return self.scoring_config.bgc_module_skip_score

    def nrp_mon_skip(self, mon_idx: int, dp_state: Optional[DP_State] = None) -> LogProb:
        mon = self.nrp_monomers[mon_idx]
        if dp_state.module_pos in (0, len(self.bgc_modules)) and mon.residue == UNKNOWN_RESIDUE:  # skipping an unknown monomer at an end
            return -2  # TODO: load from config

        return self.scoring_config.nrp_monomer_skip_score[mon.residue]  # should we take into account methylation and chirality as well?

    def match(self, bgc_pred: BGC_Module,
              nrp_mon: NRP_Monomer,
              dp_state: Optional[DP_State] = None) -> LogProb:
        if dp_state.module_pos == 2 and dp_state.monomer_pos == 2:
            pass
        residue_score = bgc_pred.residue_score[nrp_mon.residue]
        mod_match = ModMatch(mod=NRP_Monomer_Modification.METHYLATION,
                             bgc_mod=BGC_Module_Modification.METHYLATION in bgc_pred.modifications,
                             nrp_mod=NRP_Monomer_Modification.METHYLATION in nrp_mon.modifications)
        chirality_match = ChiralityMatch(bgc_epim=BGC_Module_Modification.EPIMERIZATION in bgc_pred.modifications,
                                         nrp_chr=nrp_mon.chirality)
        result = sum([self.scoring_config.match_score,
                      residue_score,
                      self.scoring_config.mod_score[mod_match],
                      self.scoring_config.chirality_score[chirality_match]])
        result = max(result, -6)  # TODO: remove this line
        return result

    def match_detailed_score(self, bgc_pred: BGC_Module,
                             nrp_mon: NRP_Monomer) -> Tuple[LogProb, LogProb, LogProb, LogProb]:
        residue_score = bgc_pred.residue_score[nrp_mon.residue]
        mod_match = ModMatch(mod=NRP_Monomer_Modification.METHYLATION,
                             bgc_mod=BGC_Module_Modification.METHYLATION in bgc_pred.modifications,
                             nrp_mod=NRP_Monomer_Modification.METHYLATION in nrp_mon.modifications)
        chirality_match = ChiralityMatch(bgc_epim=BGC_Module_Modification.EPIMERIZATION in bgc_pred.modifications,
                                         nrp_chr=nrp_mon.chirality)
        result = sum([self.scoring_config.match_score,
                      residue_score,
                      self.scoring_config.mod_score[mod_match],
                      self.scoring_config.chirality_score[chirality_match]])
        if result < -6:
            residue_score += -6 - result
        return (residue_score,
                self.scoring_config.mod_score[mod_match],
                self.scoring_config.chirality_score[chirality_match],
                self.scoring_config.match_score)

    def iterate_module(self, dp_state: Optional[DP_State] = None) -> LogProb:
        return 0

    def iterate_gene(self, dp_state: Optional[DP_State] = None) -> LogProb:
        return 0

    def null_hypothesis_score(self, nrp_monomer: NRP_Monomer) -> LogProb:
        return self.scoring_config.null_hypothesis_residue_score[nrp_monomer.residue] \
            + self.scoring_config.null_hypothesis_chirality_score[nrp_monomer.chirality] \
            + self.scoring_config.null_hypothesis_mod_score[(NRP_Monomer_Modification.METHYLATION,
                                                             NRP_Monomer_Modification.METHYLATION in nrp_monomer.modifications)]

    def skip_bgc_fragment_score(self, bgc_fragment: BGC_Fragment) -> LogProb:  # TODO: load from scoring_config
        return -3 * len(bgc_fragment)

    def skip_nrp_fragment_score(self, nrp_fragment: NRP_Fragment) -> LogProb:  # TODO: load from scoring_config
        return -3 * len(nrp_fragment.monomers)