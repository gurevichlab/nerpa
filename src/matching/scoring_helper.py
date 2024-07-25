import math

from typing import Tuple
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

    def bgc_module_skip(self, bgc_pred: BGC_Module) -> LogProb:  # argument is kept for better alignment backtracking
        return self.scoring_config.bgc_module_skip_score

    def nrp_mon_skip(self, mon: NRP_Monomer) -> LogProb:
        return self.scoring_config.nrp_monomer_skip_score[mon.residue]  # should we take into account methylation and chirality as well?

    def match(self, bgc_pred: BGC_Module,
              nrp_mon: NRP_Monomer) -> LogProb:
        residue_score = bgc_pred.residue_score[nrp_mon.residue]
        if nrp_mon.residue == UNKNOWN_RESIDUE:
            residue_score -= math.log(self.scoring_config.num_unknown_residues)
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
        if nrp_mon.residue == UNKNOWN_RESIDUE:
            residue_score -= math.log(self.scoring_config.num_unknown_residues)
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

    def iterate_module(self) -> LogProb:
        return 0

    def iterate_gene(self) -> LogProb:
        return 0

    def null_hypothesis_score(self, nrp_monomer: NRP_Monomer) -> LogProb:
        return self.scoring_config.null_hypothesis_residue_score[nrp_monomer.residue] \
            + self.scoring_config.null_hypothesis_chirality_score[nrp_monomer.chirality] \
            + self.scoring_config.null_hypothesis_mod_score[(NRP_Monomer_Modification.METHYLATION,
                                                             NRP_Monomer_Modification.METHYLATION in nrp_monomer.modifications)]

    def skip_bgc_fragment_score(self, bgc_fragment: BGC_Fragment) -> LogProb:  # TODO: load from scoring_config
        return 0

    def skip_nrp_fragment_score(self, nrp_fragment: NRP_Fragment) -> LogProb:  # TODO: load from scoring_config
        return 0