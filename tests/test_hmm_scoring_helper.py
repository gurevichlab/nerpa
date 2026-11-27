import math
import unittest
from types import SimpleNamespace
from typing import Any, cast

from src.antismash_parsing.antismash_parser_types import GeneId
from src.antismash_parsing.bgc_variant_types import (
    BGC_Module,
    BGC_Module_Modification,
)
from src.hmm.hmm_scoring_config import ChiralityMatch, MethylationMatch
from src.hmm.hmm_scoring_helper import HMMScoringHelper
from src.hmm.hmm_scoring_types import (
    ScoreComponentProbabilities,
    ScoringScheme,
    ScoringSchemeName,
)
from src.monomer_names_helper.monomer_names_helper import (
    Chirality,
    NRP_Monomer,
    NerpaResidue,
    PKS_RESIDUE,
)
from src.training.hmm_parameters.emission_probabilities.modifications_schema import (
    MODIFICATION_TYPE,
)


class HMMScoringHelperTest(unittest.TestCase):
    def setUp(self) -> None:
        self.alanine = NRP_Monomer(
            residue=NerpaResidue("Ala"),
            methylated=True,
            chirality=Chirality.D,
        )
        self.pks_monomer = NRP_Monomer(
            residue=PKS_RESIDUE,
            methylated=False,
            chirality=Chirality.L,
        )

        average_bgc_background = {
            self.alanine: ScoreComponentProbabilities(
                residue=0.08,
                methylation=0.25,
                chirality=0.30,
            ),
        }
        scoring_config = SimpleNamespace(
            methylation_score={
                MethylationMatch(False, False): math.log(0.9),
                MethylationMatch(False, True): math.log(0.1),
                MethylationMatch(True, False): math.log(0.2),
                MethylationMatch(True, True): math.log(0.8),
            },
            chirality_score={
                ChiralityMatch(False, Chirality.L): math.log(0.9),
                ChiralityMatch(False, Chirality.D): math.log(0.1),
                ChiralityMatch(False, Chirality.UNKNOWN): math.log(0.5),
                ChiralityMatch(True, Chirality.L): math.log(0.2),
                ChiralityMatch(True, Chirality.D): math.log(0.8),
                ChiralityMatch(True, Chirality.UNKNOWN): math.log(0.5),
            },
            insert_residue_probabilities={
                NerpaResidue("Ala"): 0.7,
                PKS_RESIDUE: 0.3,
            },
            norine_mods_freqs={
                MODIFICATION_TYPE.METHYLATION: 0.2,
                MODIFICATION_TYPE.EPIMERIZATION: 0.4,
            },
            scoring_schemes={
                ScoringSchemeName.LOG_PROBS: ScoringScheme(
                    include_chirality=True,
                    background_probabilities=None,
                ),
                ScoringSchemeName.NERPA: ScoringScheme(
                    include_chirality=True,
                    background_probabilities=average_bgc_background,
                ),
            },
        )
        monomer_names_helper = SimpleNamespace(
            mon_to_int={
                self.alanine: 0,
                self.pks_monomer: 1,
            },
        )
        self.helper = HMMScoringHelper(
            scoring_config=cast(Any, scoring_config),
            monomer_names_helper=cast(Any, monomer_names_helper),
        )
        self.module = BGC_Module(
            gene_id=GeneId("gene"),
            fragment_idx=0,
            a_domain_idx=0,
            genomic_context=frozenset(),
            residue_score={
                NerpaResidue("Ala"): math.log(0.6),
            },
            modifications=(
                BGC_Module_Modification.METHYLATION,
                BGC_Module_Modification.EPIMERIZATION,
            ),
            iterative_module=False,
            iterative_gene=False,
        )

    def assert_inside_outside_normalization_equal(
            self,
            raw_log_prob: float,
            normalized_score: float,
            monomer: NRP_Monomer,
    ) -> None:
        background_log_prob = self.helper.background_log_prob(
            monomer,
            ScoringSchemeName.NERPA,
        )
        self.assertAlmostEqual(
            raw_log_prob - background_log_prob,
            normalized_score,
        )

    def test_match_normalization_inside_equals_outside(self) -> None:
        raw_components = self.helper.match_log_prob_components(
            self.module,
            self.alanine,
        )
        normalized_components = self.helper.score_components(
            raw_components,
            self.alanine,
            ScoringSchemeName.NERPA,
        )

        self.assert_inside_outside_normalization_equal(
            raw_components.total_score(),
            normalized_components.total_score(),
            self.alanine,
        )

    def test_match_log_odds_uses_nerpa_scheme(self) -> None:
        raw_components = self.helper.match_log_prob_components(
            self.module,
            self.alanine,
        )
        expected_components = self.helper.score_components(
            raw_components,
            self.alanine,
            ScoringSchemeName.NERPA,
        )

        self.assertEqual(
            self.helper.match_log_odds_components(
                self.module,
                self.alanine,
            ),
            expected_components,
        )

    def test_insert_normalization_inside_equals_outside(self) -> None:
        for monomer in (self.alanine, self.pks_monomer):
            with self.subTest(monomer=monomer):
                raw_components = self.helper.insert_log_prob_components(
                    monomer
                )
                normalized_components = self.helper.score_components(
                    raw_components,
                    monomer,
                    ScoringSchemeName.NERPA,
                )

                self.assert_inside_outside_normalization_equal(
                    raw_components.total_score(),
                    normalized_components.total_score(),
                    monomer,
                )


if __name__ == "__main__":
    unittest.main()