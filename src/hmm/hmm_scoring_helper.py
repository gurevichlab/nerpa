import math

from typing import (
    Dict,
    Tuple, Literal, Optional,
)
from src.antismash_parsing.bgc_variant_types import (
    BGC_Module,
    BGC_Module_Modification,
)
from src.monomer_names_helper.monomer_names_helper import Chirality, NRP_Monomer, PKS_RESIDUE, NOT_NRPS_RESIDUE
from src.general_type_aliases import (
    LogProb,
)
from src.hmm.hmm_scoring_config import (
    ChiralityMatch,
    MethylationMatch,
    HMMScoringConfig
)
from src.hmm.hmm_scoring_types import (
    LogOddsScoreComponents,
    LogProbScoreComponents,
    ScoreComponentProbabilities,
    ScoringScheme,
    ScoringSchemeName,
)
from src.monomer_names_helper.monomer_names_helper import MonomerNamesHelper, UNKNOWN_RESIDUE
from src.general_type_aliases import Prob
from src.training.hmm_parameters.emission_probabilities.modifications_schema import MODIFICATION_TYPE as MT

from dataclasses import dataclass
from math import log, isnan


@dataclass
class HMMScoringHelper:
    scoring_config: HMMScoringConfig
    monomer_names_helper: MonomerNamesHelper

    def match_residue_score(
            self,
            bgc_module: BGC_Module,
            nrp_monomer: NRP_Monomer,
            pks_domains_in_bgc: bool
    ) -> LogProb:
        if nrp_monomer.residue in (PKS_RESIDUE, NOT_NRPS_RESIDUE):
            return float('-inf')  # PKS residues are never matched to BGC modules

        # if there are no PKS domains in the BGC, PKS hybrids are treated as unknown residues
        if nrp_monomer.is_pks_hybrid and not pks_domains_in_bgc:
            nrp_residue = UNKNOWN_RESIDUE
        else:
            nrp_residue = nrp_monomer.residue
        return bgc_module.residue_score[nrp_residue]

    def match_methylation_score(
            self,
            bgc_module: BGC_Module,
            nrp_monomer: NRP_Monomer
    ) -> LogProb:
        mod_match = MethylationMatch(BGC_Module_Modification.METHYLATION in bgc_module.modifications,
                                     nrp_monomer.methylated)
        return self.scoring_config.methylation_score[mod_match]

    def match_chirality_score(
            self,
            bgc_module: BGC_Module,
            nrp_monomer: NRP_Monomer
    ) -> LogProb:
        chirality_match = ChiralityMatch(BGC_Module_Modification.EPIMERIZATION in bgc_module.modifications,
                                         nrp_monomer.chirality)
        return self.scoring_config.chirality_score[chirality_match]

    def match_log_prob_components(
            self,
            bgc_module: BGC_Module,
            nrp_monomer: NRP_Monomer,
            pks_domains_in_bgc: bool = False
    ) -> LogProbScoreComponents:
        residue_score = self.match_residue_score(
            bgc_module,
            nrp_monomer,
            pks_domains_in_bgc,
        )
        methylation_score = self.match_methylation_score(
            bgc_module,
            nrp_monomer,
        )
        chirality_score = self.match_chirality_score(
            bgc_module,
            nrp_monomer,
        )

        return LogProbScoreComponents(
            residue_score=residue_score,
            methylation_score=methylation_score,
            chirality_score=chirality_score,
        )

    def match_log_prob(
            self,
            bgc_module: BGC_Module,
            nrp_monomer: NRP_Monomer,
            pks_domains_in_bgc: bool = False
    ) -> LogProb:
        components = self.match_log_prob_components(
            bgc_module,
            nrp_monomer,
            pks_domains_in_bgc,
        )
        return components.total_score()


    def match_log_odds_components(
            self,
            bgc_module: BGC_Module,
            nrp_monomer: NRP_Monomer,
            pks_domains_in_bgc: bool = False
    ) -> LogOddsScoreComponents:
        """Return match components under the ordinary Nerpa scheme."""
        log_prob_components = self.match_log_prob_components(
            bgc_module,
            nrp_monomer,
            pks_domains_in_bgc,
        )
        score_components = self.score_components(
            log_prob_components,
            nrp_monomer,
            ScoringSchemeName.NERPA,
        )
        if not isinstance(score_components, LogOddsScoreComponents):
            raise AssertionError(
                "The NERPA scoring scheme must produce log-odds"
            )
        return score_components

    def _background_probabilities_for_monomer(
            self,
            nrp_monomer: NRP_Monomer,
            scoring_scheme: ScoringScheme,
    ) -> Optional[ScoreComponentProbabilities]:
        background_probabilities = scoring_scheme.background_probabilities
        if background_probabilities is None:
            return None

        if nrp_monomer.residue in (PKS_RESIDUE, NOT_NRPS_RESIDUE):
            insert_components = self.insert_log_prob_components(nrp_monomer)
            return ScoreComponentProbabilities(
                residue=math.exp(insert_components.residue_score),
                methylation=math.exp(
                    insert_components.methylation_score
                ),
                chirality=math.exp(insert_components.chirality_score),
            )

        background_monomer = (
            NRP_Monomer(
                residue=UNKNOWN_RESIDUE,
                chirality=nrp_monomer.chirality,
                methylated=nrp_monomer.methylated,
                is_pks_hybrid=False,
            )
            if nrp_monomer.is_pks_hybrid
            else nrp_monomer
        )
        return background_probabilities[background_monomer]

    def background_log_prob(
            self,
            nrp_monomer: NRP_Monomer,
            scoring_scheme_name: ScoringSchemeName,
    ) -> LogProb:
        scoring_scheme = self.scoring_config.scoring_schemes[
            scoring_scheme_name
        ]
        background = self._background_probabilities_for_monomer(
            nrp_monomer,
            scoring_scheme,
        )
        if background is None:
            return 0.0

        chirality_log_prob = (
            log(background.chirality)
            if scoring_scheme.include_chirality
            else 0.0
        )
        return (
            log(background.residue)
            + log(background.methylation)
            + chirality_log_prob
        )

    def score_components(
            self,
            log_prob_components: LogProbScoreComponents,
            nrp_monomer: NRP_Monomer,
            scoring_scheme_name: ScoringSchemeName,
    ) -> LogProbScoreComponents | LogOddsScoreComponents:
        scoring_scheme = self.scoring_config.scoring_schemes[
            scoring_scheme_name
        ]
        background = self._background_probabilities_for_monomer(
            nrp_monomer,
            scoring_scheme,
        )

        if background is None:
            return LogProbScoreComponents(
                residue_score=log_prob_components.residue_score,
                methylation_score=log_prob_components.methylation_score,
                chirality_score=(
                    log_prob_components.chirality_score
                    if scoring_scheme.include_chirality
                    else 0.0
                ),
            )

        return LogOddsScoreComponents(
            residue_score=(
                float("-inf")
                if log_prob_components.residue_score == float("-inf")
                else (
                    log_prob_components.residue_score
                    - log(background.residue)
                )
            ),
            methylation_score=(
                log_prob_components.methylation_score
                - log(background.methylation)
            ),
            chirality_score=(
                log_prob_components.chirality_score
                - log(background.chirality)
                if scoring_scheme.include_chirality
                else 0.0
            ),
        )

    def get_match_emissions(
            self,
            bgc_module: BGC_Module,
            scoring_scheme_name: ScoringSchemeName = ScoringSchemeName.LOG_PROBS,
            pks_domains_in_bgc: bool = False
    ) -> Dict[NRP_Monomer, LogProb]:
        emission_scores: Dict[NRP_Monomer, LogProb] = {}

        for nrp_monomer in self.monomer_names_helper.mon_to_int:
            log_prob_components = self.match_log_prob_components(
                bgc_module,
                nrp_monomer,
                pks_domains_in_bgc,
            )
            score_components = self.score_components(
                log_prob_components,
                nrp_monomer,
                scoring_scheme_name,
            )
            emission_scores[nrp_monomer] = score_components.total_score()

            if isnan(emission_scores[nrp_monomer]):
                print(
                    "Warning: NaN score for module "
                    f"{(bgc_module.gene_id, bgc_module.a_domain_idx)} "
                    f"and monomer {nrp_monomer} with scoring scheme "
                    f"{scoring_scheme_name}"
                )

        return emission_scores

    def average_nrp_monomer_probability(
            self,
            nrp_monomer: NRP_Monomer,
    ) -> Prob:
        residue_probability = self.scoring_config.norine_residue_freqs[
            nrp_monomer.residue
        ]

        methylation_frequency = self.scoring_config.norine_mods_freqs[
            MT.METHYLATION
        ]
        methylation_probability = (
            methylation_frequency
            if nrp_monomer.methylated
            else 1 - methylation_frequency
        )

        d_chirality_frequency = self.scoring_config.norine_mods_freqs[
            MT.EPIMERIZATION
        ]
        match nrp_monomer.chirality:
            case Chirality.D:
                chirality_probability = d_chirality_frequency
            case Chirality.L:
                chirality_probability = 1 - d_chirality_frequency
            case Chirality.UNKNOWN:
                raise ValueError(
                    "Average-NRP probability requires known chirality"
                )

        return (
            residue_probability
            * methylation_probability
            * chirality_probability
        )

    def insert_log_prob_components(
            self,
            nrp_monomer: NRP_Monomer,
            pks_domains_in_bgc: bool = False,
    ) -> LogProbScoreComponents:
        residue = (
            nrp_monomer.residue
            if not nrp_monomer.is_pks_hybrid or pks_domains_in_bgc
            else UNKNOWN_RESIDUE
        )
        residue_score = log(
            self.scoring_config.insert_residue_probabilities[residue]
        )

        methylation_freq = self.scoring_config.norine_mods_freqs[
            MT.METHYLATION
        ]
        methylation_score = log(
            methylation_freq
            if nrp_monomer.methylated
            else 1 - methylation_freq
        )

        d_chirality_freq = self.scoring_config.norine_mods_freqs[
            MT.EPIMERIZATION
        ]
        match nrp_monomer.chirality:
            case Chirality.D:
                chirality_score = log(d_chirality_freq)
            case Chirality.L:
                chirality_score = log(1 - d_chirality_freq)
            case Chirality.UNKNOWN:
                chirality_score = (
                    d_chirality_freq * log(d_chirality_freq)
                    + (1 - d_chirality_freq)
                    * log(1 - d_chirality_freq)
                )

        return LogProbScoreComponents(
            residue_score=residue_score,
            methylation_score=methylation_score,
            chirality_score=chirality_score,
        )

    # TODO: precompute once and cache
    def get_insert_emissions(
            self,
            bgc_module: Optional[BGC_Module],
            scoring_scheme_name: ScoringSchemeName = ScoringSchemeName.LOG_PROBS,
            pks_domains_in_bgc: bool = False
    ) -> Dict[NRP_Monomer, LogProb]:
        scoring_scheme = self.scoring_config.scoring_schemes[
            scoring_scheme_name
        ]
        emission_scores: Dict[NRP_Monomer, LogProb] = {}

        for nrp_monomer in self.monomer_names_helper.mon_to_int:
            log_prob_components = self.insert_log_prob_components(
                nrp_monomer,
                pks_domains_in_bgc,
            )
            score_components = self.score_components(
                log_prob_components,
                nrp_monomer,
                scoring_scheme_name,
            )
            emission_scores[nrp_monomer] = score_components.total_score()

        # Only raw log-probability emissions form a probability distribution.
        if scoring_scheme.background_probabilities is None:
            if scoring_scheme.include_chirality:
                distinct_monomers = [
                    monomer
                    for monomer in emission_scores
                    if (
                        monomer.chirality != Chirality.UNKNOWN
                        and not monomer.is_pks_hybrid
                    )
                ]
            else:
                monomers_by_residue_and_methylation = {
                    (monomer.residue, monomer.methylated): monomer
                    for monomer in emission_scores
                    if not monomer.is_pks_hybrid
                }
                distinct_monomers = list(
                    monomers_by_residue_and_methylation.values()
                )

            total_probability = sum(
                math.exp(emission_scores[monomer])
                for monomer in distinct_monomers
            )
            assert math.isclose(total_probability, 1.0), (
                "Total probability of insert emissions is "
                f"{total_probability}, expected 1.0"
            )

        return emission_scores

    def get_insert_at_start_emissions(
            self,
            bgc_module: BGC_Module,
            scoring_scheme_name: ScoringSchemeName = ScoringSchemeName.LOG_PROBS,
            pks_domains_in_bgc: bool = False
    ) -> Dict[NRP_Monomer, LogProb]:
        return self.get_insert_emissions(
            bgc_module,
            scoring_scheme_name=scoring_scheme_name,
            pks_domains_in_bgc=pks_domains_in_bgc,
        )

