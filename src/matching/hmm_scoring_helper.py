import math

from typing import (
    Dict,
    Tuple, Literal,
)
from src.data_types import (
    BGC_Module,
    BGC_Module_Modification,
    Chirality,
    NRP_Monomer,
    LogProb,
)
from src.matching.hmm_scoring_config import (
    ChiralityMatch,
    MethylationMatch,
    HMMScoringConfig
)
from src.matching.alignment_step_type import MatchDetailedScore
from src.monomer_names_helper import MonomerNamesHelper, UNKNOWN_RESIDUE, PKS
from src.matching.hmm_edge_weights import get_edge_weights
from dataclasses import dataclass
from math import log


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
        if nrp_monomer.is_pks_hybrid:
            # PKS hybrids are treated as unknown residues for default scoring
            mon = NRP_Monomer(residue=UNKNOWN_RESIDUE,
                              chirality=nrp_monomer.chirality,
                              methylated=nrp_monomer.methylated,
                              is_pks_hybrid=False)
            return self.scoring_config.monomer_detailed_default_score[mon]
        else:
            return self.scoring_config.monomer_detailed_default_score[nrp_monomer]

    def monomer_default_score(self,
                              nrp_monomer: NRP_Monomer) -> LogProb:
        return sum(self.monomer_detailed_default_score(nrp_monomer))

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
                      scoring_type: Literal['LogProb', 'LogOdds'] = 'LogProb',
                      allow_unknown_chiralities: bool = True,
                      pks_domains_in_bgc: bool = False) -> Dict[NRP_Monomer, LogProb]:
        assert scoring_type in ['LogProb'], \
            f'Invalid scoring type: {scoring_type}. Expected "LogProb". "LogOdds" is deprecated.'
        emission_scores = {}
        for nrp_monomer in self.monomer_names_helper.mon_to_int:
            if (nrp_monomer.chirality == Chirality.UNKNOWN
                    and not allow_unknown_chiralities):
                emission_scores[nrp_monomer] = float('-inf')
                continue

            if scoring_type == 'LogProb':
                emission_scores[nrp_monomer] = self.match(bgc_module, nrp_monomer, pks_domains_in_bgc)
            elif scoring_type == 'LogOdds':
                detailed_score = self.normalized_match_detailed_score(bgc_module, nrp_monomer, pks_domains_in_bgc)
                emission_scores[nrp_monomer] = (detailed_score.residue_score
                                          + detailed_score.methylation_score
                                          + detailed_score.chirality_score)

        return emission_scores

    # TODO: precompute once and cache
    def get_insert_emissions(self,
                             bgc_module: BGC_Module,
                             scoring_type: Literal['LogProb', 'LogOdds'] = 'LogProb',
                             allow_unknown_chiralities: bool = True,
                             pks_domains_in_bgc: bool = False) -> Dict[NRP_Monomer, LogProb]:
        # adjust default monomer frequencies to force the probability of UNKNOWN_RESIDUE to be UNKNOWN_INSERT_FREQ
        UNKNOWN_INSERT_FREQ = 0.77  # TODO: make this configurable
        scale_factor = (1 - UNKNOWN_INSERT_FREQ) / (1 - self.monomer_names_helper.default_frequencies.residue[UNKNOWN_RESIDUE])
        default_insert_freqs = {res: freq * scale_factor
                                for res, freq in self.monomer_names_helper.default_frequencies.residue.items()}
        default_insert_freqs[UNKNOWN_RESIDUE] = UNKNOWN_INSERT_FREQ

        assert scoring_type in ['LogProb'], \
            (f'Invalid scoring type: {scoring_type}. Expected "LogProb" '
             f'("LogOdds" deprecated for now).')

        emission_scores = {}
        for nrp_monomer in self.monomer_names_helper.mon_to_int:
            if (not allow_unknown_chiralities
                    and (nrp_monomer.chirality == Chirality.UNKNOWN
                         or nrp_monomer.is_pks_hybrid)):
                emission_scores[nrp_monomer] = float('-inf')
                continue

            # If the monomer is a PKS hybrid, we treat it as an unknown residue
            res = nrp_monomer.residue if not nrp_monomer.is_pks_hybrid else UNKNOWN_RESIDUE

            res_freq = default_insert_freqs[res]
            res_score = log(res_freq)
            meth_freq = (self.monomer_names_helper.default_frequencies.methylation
                            if nrp_monomer.methylated else
                            1 - self.monomer_names_helper.default_frequencies.methylation)
            meth_score = log(meth_freq)

            d_chr_freq = self.monomer_names_helper.default_frequencies.d_chirality
            match nrp_monomer.chirality:
                case Chirality.D:
                    chr_score = log(d_chr_freq)
                case Chirality.L:
                    chr_score = log(1 - d_chr_freq)
                case Chirality.UNKNOWN:
                    chr_score = (
                            d_chr_freq * log(d_chr_freq)
                            + (1 - d_chr_freq) * log(1 - d_chr_freq)
                    )

            emission_scores[nrp_monomer] = res_score + meth_score + chr_score


        # assert that all scores sum to 1
        total_score = sum(math.e ** score
                          for mon, score in emission_scores.items()
                          if mon.chirality != Chirality.UNKNOWN
                          and not mon.is_pks_hybrid)
        assert math.isclose(total_score, 1.0), \
            f'Total score of insert emissions is {total_score}, expected 1.0'

        return emission_scores

    def get_insert_at_start_emissions(self,
                                      bgc_module: BGC_Module,
                                      scoring_type: Literal['LogProb', 'LogOdds'] = 'LogProb',
                                      allow_unknown_chiralities: bool = True,
                                      pks_domains_in_bgc: bool = False) -> Dict[NRP_Monomer, LogProb]:
        return self.get_insert_emissions(bgc_module,
                                         scoring_type,
                                         allow_unknown_chiralities,
                                         pks_domains_in_bgc)

    def get_edge_weights(self, hmm) -> Dict[Tuple[int, int], LogProb]:
        return get_edge_weights(hmm, self.scoring_config)
