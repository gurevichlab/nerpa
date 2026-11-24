from collections import defaultdict
import math
from math import log
from typing import (
    Dict,
    List,
    NamedTuple
)

import polars as pl
from src.aa_specificity_prediction_model.specificity_prediction_helper import SpecificityPredictionHelper
from src.general_type_aliases import (
    LogProb,
     Prob,
)
from src.antismash_parsing.genomic_context import (
    ModuleGenomicContext,
    ModuleGenomicContextFeature
)
from src.generic.numeric import safe_log
from src.hmm.hmm_auxiliary_types import DetailedHMMEdgeType, DetailedHMMStateType
from src.matching.alignment_step_type import MatchDetailedScore

from dataclasses import dataclass
from pathlib import Path
import yaml

from src.monomer_names_helper.monomer_names_helper import (
    Chirality,
    NRP_Monomer,
    NerpaResidue,
    MonomerNamesHelper,
    MonomersDefaultFrequencies
)
from src.training.hmm_parameters.norine_monomers_info import NorineMonomersInfo
from src.training.hmm_parameters.emission_probabilities.modifications_schema import (
    ModificationsProbs_df,
    MODIFICATION_TYPE,
    ModificationsProbsSchema
)

ST = DetailedHMMStateType
ET = DetailedHMMEdgeType
MGF = ModuleGenomicContextFeature
GC = ModuleGenomicContext

class ChiralityMatch(NamedTuple):
    bgc_epim: bool
    nrp_chr: Chirality

class MethylationMatch(NamedTuple):
    bgc_meth: bool
    nrp_meth: bool


            
def load_methylation_score(
        mod_freqs_df: ModificationsProbs_df
) -> Dict[MethylationMatch, LogProb]:
    MT = MODIFICATION_TYPE
    MPS = ModificationsProbsSchema

    meth_freqs_df = mod_freqs_df.df.filter(
        pl.col(MPS.MODIFICATION) == MT.METHYLATION
    )
    return {
        MethylationMatch(
            bgc_meth=row[MPS.PRESENT_IN_BGC],
            nrp_meth=row[MPS.PRESENT_IN_NRP]
        ):
        log(row[MPS.PROBABILITY])

        for row in meth_freqs_df.iter_rows(named=True)
    }


def load_chirality_score(
        mod_freqs_df: ModificationsProbs_df,
        norine_monomers_info: NorineMonomersInfo,
) -> Dict[ChiralityMatch, LogProb]:
    MT = MODIFICATION_TYPE
    MPS = ModificationsProbsSchema

    chr_freqs_df = mod_freqs_df.df.filter(
        pl.col(MPS.MODIFICATION) == MT.EPIMERIZATION
    )

    chr_score: Dict[ChiralityMatch, LogProb] = {}
    # for known chirality, we use the estimated probabilities directly
    for row in chr_freqs_df.iter_rows(named=True):
        nrp_chr = (
            Chirality.D
            if row[MPS.PRESENT_IN_NRP]
            else Chirality.L
        )
        chr_match = ChiralityMatch(
            bgc_epim=row[MPS.PRESENT_IN_BGC],
            nrp_chr=nrp_chr
        )
        chr_score[chr_match] = log(row[MPS.PROBABILITY])

    # for unknown chirality, we compute the weighted average,
    # based on the adundance of D chirality in the Norine monomers database
    default_d_chr_freq = (
        norine_monomers_info
        .modification_frequencies()
        [MT.EPIMERIZATION]
    )

    for bgc_epim in (False, True):
        chr_match_unk = ChiralityMatch(
            bgc_epim=bgc_epim, nrp_chr=Chirality.UNKNOWN
        )
        chr_match_nrp_false = ChiralityMatch(
            bgc_epim=bgc_epim, nrp_chr=Chirality.L
        )
        chr_match_nrp_true = ChiralityMatch(
            bgc_epim=bgc_epim, nrp_chr=Chirality.D
        )
        chr_score[chr_match_unk] = (
            default_d_chr_freq * chr_score[chr_match_nrp_true]
            + (1 - default_d_chr_freq) * chr_score[chr_match_nrp_false]
        )

    return chr_score

import math
from src.generic.numeric import clamp

def load_edge_weight_params(cfg: dict) -> Dict[ST, Dict[ModuleGenomicContext, Dict[ET, LogProb]]]:
    edge_weights_cfg = pl.read_csv(cfg['edge_weight_parameters'], separator='\t')
    ET = DetailedHMMEdgeType
    MGF = ModuleGenomicContextFeature
    min_allowed_log_prob = -3  # a quick fix because I am too lazy to update training right now
    max_allowed_log_prob = math.log(1 - math.e ** min_allowed_log_prob)

    parsed_data: Dict[ST, Dict[GC, Dict[ET, LogProb]]] = (
        defaultdict(lambda: defaultdict(dict))
    )
    for row in edge_weights_cfg.iter_rows(named=True):
        state_type = ST[row["STATE"]]
        edge_type = ET[row["EDGE_TYPE"]]
        gc = frozenset(MGF[feature_str] for feature_str in eval(row["RELEVANT_GC"]))

        prob: float = row["PROB"]
        match prob:
            case 0.0:
                log_prob = -math.inf
            case 1.0:
                log_prob = 0.0
            case _:
                log_prob = clamp(
                    math.log(prob),
                    min_allowed_log_prob,
                    max_allowed_log_prob,
                )

        parsed_data[state_type][gc][edge_type] = log_prob

    return parsed_data


def compute_monomers_default_detailed_score(
        paras_default_residue_freqs: Dict[NerpaResidue, Prob],
        methylation_score: Dict[MethylationMatch, LogProb],
        chirality_score: Dict[ChiralityMatch, LogProb],
        methylation_freq: Prob,
        epimerization_freq: Prob
) -> Dict[NRP_Monomer, MatchDetailedScore]:
    monomer_default_score = {}
    for residue, residue_freq in paras_default_residue_freqs.items():
        for nrp_meth in (False, True):
            for nrp_chr in Chirality:
                residue_score = safe_log(residue_freq)
                default_methylation_score = (
                    methylation_freq * methylation_score[MethylationMatch(bgc_meth=True, nrp_meth=nrp_meth)]
                    + (1 - methylation_freq) * methylation_score[MethylationMatch(bgc_meth=False, nrp_meth=nrp_meth)]
                )
                default_chirality_score = (
                    epimerization_freq * chirality_score[ChiralityMatch(bgc_epim=True, nrp_chr=nrp_chr)]
                    + (1 - epimerization_freq) * chirality_score[ChiralityMatch(bgc_epim=False, nrp_chr=nrp_chr)]
                )

                monomer = NRP_Monomer(
                    residue=residue,
                    methylated=nrp_meth,
                    chirality=nrp_chr
                )
                monomer_default_score[monomer] = MatchDetailedScore(
                    residue_score,
                    default_methylation_score,
                    default_chirality_score
                )

    return monomer_default_score


def load_paras_avg_residue_scores(
        paras_results_tsv: Path,
        specificity_prediction_helper: SpecificityPredictionHelper,
) -> Dict[NerpaResidue, Prob]:
    assert specificity_prediction_helper.DEFAULT_MODEL == 'paras', \
        'PARAS residue frequencies can be loaded only if PARAS is used as the default specificity prediction model.'
    paras_results: pl.DataFrame = (
        pl.read_csv(paras_results_tsv, separator='\t')
    )
    predictions_by_id_df: pl.DataFrame = (
        paras_results
        .group_by("id", maintain_order=True)
        .agg(
            substrates=pl.col("substrate"),
            scores=pl.col("score"),
        )
    )
    paras_predictions: List[Dict[str, float]] = [
        dict(zip(row["substrates"], row["scores"]))
        for row in predictions_by_id_df.iter_rows(named=True)
    ]

    predictions = [
        specificity_prediction_helper.paras_predictions_to_nerpa_predictions(a_domain_prediction)
        for a_domain_prediction in paras_predictions
    ]

    calibrated_predictions = [
        specificity_prediction_helper.calibrate_scores(
            prediction,
            specificity_prediction_helper._calibration_step_function_paras
        )
        for prediction in predictions
    ]

    residues = calibrated_predictions[0].keys()
    averaged_predictions = {
        res: (
            sum(pred[res] for pred in calibrated_predictions)
            / len(calibrated_predictions)
        )
        for res in residues
    }

    assert math.isclose(sum(averaged_predictions.values()), 1.0, rel_tol=1e-4), \
        f'Total score of averaged PARAS predictions is {sum(averaged_predictions.values())}, expected 1.0'

    return averaged_predictions


class HMMScoringConfig:
    methylation_score: Dict[MethylationMatch, LogProb]
    chirality_score: Dict[ChiralityMatch, LogProb]

    edge_weight_parameters: Dict[ST, Dict[GC, Dict[ET, LogProb]]]

    # the score of a monomer averaged over all BGC modules from MIBiG
    # used for normalization
    monomer_avg_detailed_score: Dict[NRP_Monomer, MatchDetailedScore]

    # probability that an inserted monomer has unknown residue
    insert_unknown_freq: Prob

    # residue frequencies estimated from the Norine database
    norine_residue_freqs: Dict[NerpaResidue, Prob]
    # METHYLATION and EPIMERIZATION frequencies estimated from the Norine database
    norine_mods_freqs: Dict[MODIFICATION_TYPE, Prob]

    def __init__(
            self,
            nerpa_dir: Path,
            path_to_config: Path,
            specificity_prediction_helper: SpecificityPredictionHelper,
            norine_monomers_info: NorineMonomersInfo,
    ):
        cfg = yaml.safe_load(path_to_config.open('r'))
        mod_freqs_df = ModificationsProbs_df(
            pl.read_csv(
                cfg['emission_parameters']['modification_likelihoods'],
                separator='\t'
            )
        )

        self.methylation_score: Dict[MethylationMatch, LogProb] = (
            load_methylation_score(mod_freqs_df)
        )
        self.chirality_score: Dict[ChiralityMatch, LogProb] = (
            load_chirality_score(mod_freqs_df, norine_monomers_info)
        )

        paras_prediction_tsv: Path = nerpa_dir / cfg['paras_predictions_table']
        paras_avg_residue_scores: Dict[NerpaResidue, Prob] = (
            load_paras_avg_residue_scores(
                paras_prediction_tsv,
                specificity_prediction_helper,
            )
        )
        self.monomer_avg_detailed_score: Dict[NRP_Monomer, MatchDetailedScore] = (
            compute_monomers_default_detailed_score(
                paras_default_residue_freqs=paras_avg_residue_scores,
                methylation_freq=cfg['emission_parameters']['mibig_methylation_freq'],
                epimerization_freq=cfg['emission_parameters']['mibig_epimerization_freq'],
                methylation_score=self.methylation_score,
                chirality_score=self.chirality_score
            )
        )

        monomer_names_helper: MonomerNamesHelper = (
            specificity_prediction_helper.monomer_names_helper
        )
        self.norine_residue_freqs: Dict[NerpaResidue, Prob] = (
            norine_monomers_info.residue_frequencies(
                supported_residues=set(monomer_names_helper.all_residues()),
                min_allowed_freq=0.01, # TODO: put in config
            )
        )
        self.norine_mods_freqs: Dict[MODIFICATION_TYPE, Prob] = (
            norine_monomers_info.modification_frequencies()
        )
        self.insert_unknown_freq: Prob = cfg['emission_parameters']['insert_unknown_freq']
        self.edge_weight_parameters: Dict[ST, Dict[GC, Dict[ET, LogProb]]] = (
            load_edge_weight_params(cfg)
        )

