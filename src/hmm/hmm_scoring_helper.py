import math

from typing import (
    Dict,
    Tuple, Literal, Optional,
)
from src.antismash_parsing.bgc_variant_types import (
    BGC_Module,
    BGC_Module_Modification,
)
from src.monomer_names_helper import Chirality, NRP_Monomer, PKS_RESIDUE, NOT_NRPS_RESIDUE
from src.general_type_aliases import (
    LogProb,
)
from src.hmm.hmm_scoring_config import (
    ChiralityMatch,
    MethylationMatch,
    HMMScoringConfig
)
from src.matching.alignment_step_type import MatchDetailedScore
from src.monomer_names_helper import MonomerNamesHelper, UNKNOWN_RESIDUE
from dataclasses import dataclass
from math import log, isnan
from enum import Enum, auto

class ScoringScheme(Enum):
    NERPA = auto()  # LogProb, unknown chiralities allowed
    NERPA_MS = auto()  # LogOdds, chirality score component removed


@dataclass
class HMMHelper:
    scoring_config: HMMScoringConfig
    monomer_names_helper: MonomerNamesHelper

    def match_residue_score(self,
                            bgc_module: BGC_Module,
                            nrp_monomer: NRP_Monomer,
                            pks_domains_in_bgc: bool) -> LogProb:
        if nrp_monomer.residue in (PKS_RESIDUE, NOT_NRPS_RESIDUE):
            return float('-inf')  # PKS residues are never matched to BGC modules

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
        if nrp_monomer.residue in (PKS_RESIDUE, NOT_NRPS_RESIDUE):
            return MatchDetailedScore(float('-inf'), 0, 0)
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
          res_score = (
              float('-inf')
              if match_score.residue_score == float('-inf')
              else match_score.residue_score - default_score.residue_score
          )
          return MatchDetailedScore(
                res_score,
                match_score.methylation_score - default_score.methylation_score,
                match_score.chirality_score - default_score.chirality_score
          )

    def get_match_emissions(
            self,
            bgc_module: BGC_Module,
            scoring_scheme: ScoringScheme = ScoringScheme.NERPA,
            pks_domains_in_bgc: bool = False
    ) -> Dict[NRP_Monomer, LogProb]:
        emission_scores = {}
        for nrp_monomer in self.monomer_names_helper.mon_to_int:
            match scoring_scheme:
                case ScoringScheme.NERPA:
                     emission_scores[nrp_monomer] = self.match(bgc_module, nrp_monomer, pks_domains_in_bgc)
                case ScoringScheme.NERPA_MS:
                    detailed_score = self.normalized_match_detailed_score(bgc_module, nrp_monomer, pks_domains_in_bgc)
                    emission_scores[nrp_monomer] = detailed_score.residue_score + detailed_score.methylation_score
                case _:
                    raise ValueError(f'Invalid scoring scheme: {scoring_scheme}')

            if isnan(emission_scores[nrp_monomer]):
                print(f'Warning: NaN score for module {(bgc_module.gene_id, bgc_module.a_domain_idx)} and monomer {nrp_monomer} with scoring scheme {scoring_scheme}')

        return emission_scores

    # TODO: precompute once and cache
    def get_insert_emissions(self,
                             bgc_module: Optional[BGC_Module],
                             scoring_scheme: ScoringScheme = ScoringScheme.NERPA,
                             pks_domains_in_bgc: bool = False) -> Dict[NRP_Monomer, LogProb]:
        # adjust default monomer frequencies to force the probability of UNKNOWN_RESIDUE to be UNKNOWN_INSERT_FREQ
        UNKNOWN_INSERT_FREQ = self.scoring_config.insert_unknown_freq
        scale_factor = (1 - UNKNOWN_INSERT_FREQ) / (1 - self.monomer_names_helper.default_frequencies.residue[UNKNOWN_RESIDUE])
        default_insert_freqs = {res: freq * scale_factor
                                for res, freq in self.monomer_names_helper.default_frequencies.residue.items()}
        default_insert_freqs[UNKNOWN_RESIDUE] = UNKNOWN_INSERT_FREQ

        assert math.isclose(sum(default_insert_freqs.values()), 1.0), \
            f'Default insert frequencies do not sum to 1: {sum(default_insert_freqs.values())}'

        mon_helper_residues = {mon.residue for mon in self.monomer_names_helper.mon_to_int}
        assert set(default_insert_freqs.keys()) == mon_helper_residues, \
            ('Default insert frequencies do not match monomer names helper residues.\n'
             'Probably some nerpa residues are missing in the Norine database.\n'
             f'insert_residues - nerpa_residues: {set(default_insert_freqs.keys()) - mon_helper_residues}\n'
             f'nerpa_residues - insert_residues: {mon_helper_residues - set(default_insert_freqs.keys())}')

        emission_scores = {}
        for nrp_monomer in self.monomer_names_helper.mon_to_int:
            res = (
                nrp_monomer.residue 
                if not nrp_monomer.is_pks_hybrid or pks_domains_in_bgc
                else UNKNOWN_RESIDUE
            )

            res_freq = default_insert_freqs[res]
            res_score = log(res_freq)
            meth_freq = (
                self.monomer_names_helper.default_frequencies.methylation
                if nrp_monomer.methylated else
                1 - self.monomer_names_helper.default_frequencies.methylation
            )
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
            match scoring_scheme:
                case ScoringScheme.NERPA:
                    emission_scores[nrp_monomer] = res_score + meth_score + chr_score
                case ScoringScheme.NERPA_MS:
                    emission_scores[nrp_monomer] = res_score + meth_score
                case _:
                    raise ValueError(f'Invalid scoring scheme: {scoring_scheme}')

        distinct_proper_monomers = []
        match scoring_scheme:
            case ScoringScheme.NERPA:
                distinct_proper_monomers = [
                    mon
                    for mon in emission_scores
                    if mon.chirality != Chirality.UNKNOWN and not mon.is_pks_hybrid
                ]
            case ScoringScheme.NERPA_MS:
                mon_by_res_meth = {
                    (mon.residue, mon.methylated): mon
                    for mon in emission_scores
                    if not mon.is_pks_hybrid
                }
                distinct_proper_monomers = list(mon_by_res_meth.values())               
            case _:
                raise ValueError(f'Invalid scoring scheme: {scoring_scheme}')

        # assert that all scores sum to 1
        total_score = sum(math.e ** emission_scores[mon]
                          for mon in distinct_proper_monomers)

        assert math.isclose(total_score, 1.0), \
            f'Total score of insert emissions is {total_score}, expected 1.0'

        return emission_scores

    def get_insert_at_start_emissions(self,
                                      bgc_module: BGC_Module,
                                      scoring_scheme: ScoringScheme = ScoringScheme.NERPA,
                                      pks_domains_in_bgc: bool = False) -> Dict[NRP_Monomer, LogProb]:
        return self.get_insert_emissions(bgc_module,
                                         scoring_scheme=scoring_scheme,
                                         pks_domains_in_bgc=pks_domains_in_bgc)

