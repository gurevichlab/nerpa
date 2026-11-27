from enum import StrEnum
from typing import TYPE_CHECKING, NamedTuple

from src.general_type_aliases import LogOdds, LogProb, Prob

if TYPE_CHECKING:
    from src.monomer_names_helper.monomer_names_helper import NRP_Monomer


class ScoringSchemeName(StrEnum):
    LOG_PROBS = "LOG_PROBS"
    NERPA = "NERPA"
    NERPA_MS = "NERPA_MS"


class ScoreComponentProbabilities(NamedTuple):
    residue: Prob
    methylation: Prob
    chirality: Prob


class ScoringScheme(NamedTuple):
    include_chirality: bool
    background_probabilities: (
        dict["NRP_Monomer", ScoreComponentProbabilities] | None
    )


class LogProbScoreComponents(NamedTuple):
    residue_score: LogProb
    methylation_score: LogProb
    chirality_score: LogProb

    def total_score(self) -> LogProb:
        return (
            self.residue_score
            + self.methylation_score
            + self.chirality_score
        )


class LogOddsScoreComponents(NamedTuple):
    residue_score: LogOdds
    methylation_score: LogOdds
    chirality_score: LogOdds

    def total_score(self) -> LogOdds:
        return (
            self.residue_score
            + self.methylation_score
            + self.chirality_score
        )