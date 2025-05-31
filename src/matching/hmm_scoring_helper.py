import math

from typing import (
    Dict,
    Tuple,
)
from src.data_types import (
    BGC_Module,
    BGC_Module_Modification,
    Chirality,
    NRP_Monomer,
    LogProb,
)
from src.matching.hmm_config import (
    ChiralityMatch,
    MethylationMatch,
    HMMScoringConfig
)
from src.matching.alignment_step_type import MatchDetailedScore
from src.monomer_names_helper import MonomerNamesHelper, UNKNOWN_RESIDUE, PKS
from src.matching.hmm_edge_weights import get_edge_weights
from dataclasses import dataclass


@dataclass
class HMMHelper:
    scoring_config: HMMScoringConfig
    monomer_names_helper: MonomerNamesHelper

    def match_residue_score(self,
                            bgc_module: BGC_Module,
                            nrp_monomer: NRP_Monomer,
                            pks_domains_in_bgc: bool) -> LogProb:
        if nrp_monomer.residue == PKS:
            return -math.inf  # PKS residues are never matched to BGC modules

        # if there are no PKS domains in the BGC, PKS hybrids are treated as unknown residues
        if nrp_monomer.is_pks_hybrid and not pks_domains_in_bgc:
            nrp_residue = UNKNOWN_RESIDUE
        else:
            nrp_residue = nrp_monomer.residue
        return bgc_module.residue_score[nrp_residue]

    def match_methylation_score(self,
                                bgc_module: BGC_Module,
                                nrp_monomer: NRP_Monomer) -> LogProb:
        mod_match = MethylationMatch(BGC_Module_Modification.METHYLATION in bgc_module.modifications,
                                     nrp_monomer.methylated)
        return self.scoring_config.methylation_score[mod_match]

    def match_chirality_score(self,
                              bgc_module: BGC_Module,
                              nrp_monomer: NRP_Monomer) -> LogProb:
        chirality_match = ChiralityMatch(BGC_Module_Modification.EPIMERIZATION in bgc_module.modifications,
                                         nrp_monomer.chirality)
        return self.scoring_config.chirality_score[chirality_match]

    def match_detailed_score(self,
                             bgc_module: BGC_Module,
                             nrp_monomer: NRP_Monomer,
                             pks_domains_in_bgc: bool = False) -> MatchDetailedScore:
        residue_score = self.match_residue_score(bgc_module, nrp_monomer, pks_domains_in_bgc)
        methylation_score = self.match_methylation_score(bgc_module, nrp_monomer)
        chirality_score = self.match_chirality_score(bgc_module, nrp_monomer)

        return MatchDetailedScore(residue_score,
                                  methylation_score,
                                  chirality_score)

    def match(self,
              bgc_module: BGC_Module,
              nrp_monomer: NRP_Monomer,
              pks_domains_in_bgc: bool = False) -> LogProb:
        return sum(self.match_detailed_score(bgc_module, nrp_monomer, pks_domains_in_bgc))

    def monomer_detailed_default_score(self,
                                       nrp_monomer: NRP_Monomer) -> MatchDetailedScore:
        return self.scoring_config.monomer_detailed_default_score[nrp_monomer]

    def monomer_default_score(self,
                              nrp_monomer: NRP_Monomer) -> LogProb:
        return sum(self.scoring_config.monomer_detailed_default_score[nrp_monomer])

    def normalized_match_detailed_score(self,
                                        bgc_module: BGC_Module,
                                        nrp_monomer: NRP_Monomer,
                                        pks_domains_in_bgc: bool = False) -> MatchDetailedScore:
          """
          Returns the match detailed score normalized by the sum of the default scores.
          """
          match_score = self.match_detailed_score(bgc_module, nrp_monomer, pks_domains_in_bgc)
          default_score = self.monomer_detailed_default_score(nrp_monomer)
          return MatchDetailedScore(
                match_score.residue_score - default_score.residue_score,
                match_score.methylation_score - default_score.methylation_score,
                match_score.chirality_score - default_score.chirality_score
          )

    def get_emissions(self,
                      bgc_module: BGC_Module,
                      pks_domains_in_bgc: bool = False) -> Dict[NRP_Monomer, LogProb]:
        return {nrp_monomer: self.match(bgc_module, nrp_monomer, pks_domains_in_bgc)
                for nrp_monomer in self.monomer_names_helper.mon_to_int}

    def get_insert_emissions(self,
                             bgc_module: BGC_Module,
                             pks_domains_in_bgc: bool = False) -> Dict[NRP_Monomer, LogProb]:
        return {nrp_monomer: 0 for nrp_monomer in self.monomer_names_helper.mon_to_int}

    def get_insert_at_start_emissions(self,
                                      bgc_module: BGC_Module,
                                      pks_domains_in_bgc: bool = False) -> Dict[NRP_Monomer, LogProb]:
        return {nrp_monomer: 0 for nrp_monomer in self.monomer_names_helper.mon_to_int}

    def get_edge_weights(self, hmm) -> Dict[Tuple[int, int], LogProb]:
        return get_edge_weights(hmm, self.scoring_config)
