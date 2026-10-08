from typing import (
    Callable,
    Dict,
    List,
    Tuple,
    NamedTuple,
    Optional,
)
from collections import Counter, defaultdict

from src.aa_specificity_prediction_model.specificity_prediction_helper import create_step_function
from src.antismash_parsing.bgc_variant_types import BGC_Module, BGC_Module_Modification
from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import INSERT_STATE_TYPES, MATCHING_STATE_TYPES
from src.training.hmm_parameters.emission_probabilities.step_function import (
    create_bins,
    fit_step_function_to_bins,
    plot_step_function,
    plot_step_function_stacked,
)
from src.monomer_names_helper import NRP_Monomer, UNKNOWN_RESIDUE, MonomerNamesHelper, NerpaResidue
from src.general_type_aliases import LogProb, Prob
from src.training.hmm_parameters.norine_monomers_info import NorineMonomersInfo
from src.training.hmm_parameters.training_types import HMM_Chosen_Steps_df, HMM_Chosen_Steps_Schema
from src.training.hmm_parameters.short_aliases import *
from logging import Logger
from src.training.hmm_parameters.emission_probabilities.modifications_schema import ModificationsProbs_df
from src.training.hmm_parameters.emission_probabilities.modifications_schema import MODIFICATION_TYPE
from src.training.hmm_parameters.norine_monomers_info import NorineMonomersInfo
from src.training.hmm_parameters.statistics_helpers import get_profile_likelihood_ci

from pathlib import Path
from math import log
import numpy as np
import matplotlib.pyplot as plt
import polars as pl
from math import e as math_e
from scipy.optimize import minimize_scalar, brentq
from scipy.stats import chi2

class StepFunctionEstimationConfig(NamedTuple):
    num_bins: int = 20
    step_range: int = 1000

def get_match_emissions(hmm_chosen_steps_df: HMM_Chosen_Steps_df) -> List[Tuple[BGC_Module, NRP_Monomer]]:
    matching_state_types_names = [
        state_type.name for state_type in MATCHING_STATE_TYPES
    ]
    emissions_info_df = (
        hmm_chosen_steps_df.df
        .filter(
            pl.col(SC.STATE_TYPE).is_in(matching_state_types_names)
        )
        .select([SC.MODULE, SC.EMISSION])
    )
    match_emissions: List[Tuple[BGC_Module, NRP_Monomer]] = [
        (row[SC.MODULE], row[SC.EMISSION])
        for row in emissions_info_df.iter_rows(named=True)
    ]
    return match_emissions


def get_score_correctness(
        match_emissions: List[Tuple[BGC_Module, NRP_Monomer]]
) -> List[Tuple[LogProb, bool]]:
    score_correctness: List[Tuple[LogProb, bool]] = []
    for bgc_module, nrp_monomer in match_emissions:
        for predicted_residue, score in bgc_module.residue_score.items():
            score_correctness.append(
                (score, predicted_residue == nrp_monomer.residue)
            )
            # TODO: check that PKS-hybrids are handled correctly
    return score_correctness

def fit_step_function(
        match_emissions: List[Tuple[BGC_Module, NRP_Monomer]],
        cfg: StepFunctionEstimationConfig,
        output_dir: Optional[Path] = None,
) -> List[float]:
    score_correctness: List[Tuple[LogProb, bool]] = (
        get_score_correctness(match_emissions)
    )
    score_correctness_bins: List[Tuple[int, int]] = (
        create_bins(score_correctness, cfg.num_bins)
    )
    step_function: List[float] = (
        fit_step_function_to_bins(score_correctness_bins, cfg.step_range)
    )

    if output_dir is not None:
        plot_step_function(
            score_correctness_bins,
            step_function,
            out_file=output_dir / Path('step_function.png')
        )
        plot_step_function_stacked(
            score_correctness_bins,
            step_function,
            out_file=output_dir / Path('step_function_stacked.png')
        )

    return step_function


def get_insert_unknown_prob(hmm_chosen_steps: HMM_Chosen_Steps_df) -> Prob:
    insert_state_types_names = [
        state_type.name for state_type in INSERT_STATE_TYPES
    ]
    inserts_df = (
        hmm_chosen_steps.df
        .filter(
            pl.col(SC.STATE_TYPE).is_in(insert_state_types_names)
        )
        .with_columns(
            pl.col(SC.EMISSION)
            .map_elements(
                lambda emission: emission.to_base_mon().residue == UNKNOWN_RESIDUE,
                return_dtype=pl.Boolean
            )
            .alias("inserted_unknown")
        )
    )

    return Prob(
        inserts_df["inserted_unknown"].sum() / len(inserts_df)
    )


def get_unknown_because_pks_prob(
        match_emissions: List[Tuple[BGC_Module, NRP_Monomer]],
        calibration_function: Callable[[Prob], Prob],
        confidence_level: float = 0.95
) -> Tuple[Prob, Tuple[Prob, Prob]]:
    """
    Model:
    modules m1, m2, ..., mn attract UNKNOWN residue with known probabilities p1, p2, ..., pn
    PKS domains can independently turn a residue into UNKNOWN with an unknown probability x

    For a module with baseline UNKNOWN probability p:
        P(UNKNOWN | PKS_DOWNSTREAM) = p + (1 - p) * x

    Goal:
    estimate x given the observed correctness of the predictions for the UNKNOWN residue
    Returns:
        Tuple[Prob, Tuple[Prob, Prob]]: (MLE of x, (lower bound of CI, upper bound of CI))
    """

    # 1. Build a polars table with the relevant information for each emission
    rows = [
        {
        "PKS_DOWNSTREAM": (
            ModuleGenomicContextFeature.PKS_DOWNSTREAM in bgc_module.genomic_context
        ),
        "UNKNOWN_RESIDUE": (
            nrp_monomer.residue == UNKNOWN_RESIDUE
        ),
        "UNKNOWN_PREDICTED_PROB": (
            calibration_function(
                math_e ** bgc_module.residue_score[UNKNOWN_RESIDUE]
            )
        )
        }
        for bgc_module, nrp_monomer in match_emissions
    ]

    df = pl.DataFrame(rows)
    pks_df = df.filter(pl.col("PKS_DOWNSTREAM"))

    unknown_observed = pks_df["UNKNOWN_RESIDUE"].to_numpy()
    predicted_probs = pks_df["UNKNOWN_PREDICTED_PROB"].to_numpy()

    def log_likelihood(x: float) -> float:
        unknown_probs = (
            predicted_probs
            + (1.0 - predicted_probs) * x
        )

        return float(
            np.where(
                unknown_observed,
                np.log(unknown_probs),
                np.log1p(-unknown_probs)
            ).sum()
        )

    # Maximum-likelihood estimate.
    result = minimize_scalar(
        lambda x: -log_likelihood(x),
        bounds=(0.0, 1.0),
        method="bounded"
    )

    if not result.success:
        raise RuntimeError(
            f"Failed to estimate PKS UNKNOWN probability: {result.message}"
        )

    x_mle = float(result.x)
    ci_lower, ci_upper = get_profile_likelihood_ci(
        log_likelihood=log_likelihood,
        mle=x_mle,
        parameter_bounds=(0.0, 1.0),
        confidence_level=confidence_level
    )
    return (
        Prob(x_mle),
        (Prob(ci_lower), Prob(ci_upper))
    )



class EmissionParams(NamedTuple):
    step_function: List[float]
    modifications_frequences: List[dict]
    unknown_because_pks_prob: Prob
    insert_unknown_prob: Prob  # probability of inserting UNKNOWN residue
    default_modifications_frequencies: Dict[str, Prob]  # default frequencies for modifications, used for score normalization
    default_residue_frequencies: Dict[str, Prob]  # default probabilities for residues, used for score normalization


def infer_emission_params(
        hmm_chosen_steps_df: HMM_Chosen_Steps_df,
        norine_monomers_info: NorineMonomersInfo,
        step_function_cfg: StepFunctionEstimationConfig,
        logger: Logger,
        output_dir: Optional[Path] = None,
) -> EmissionParams:
    logger.info('Building step function...')
    match_emissions: List[Tuple[BGC_Module, NRP_Monomer]] = (
        get_match_emissions(hmm_chosen_steps_df)
    )
    step_function_steps = fit_step_function(
        match_emissions=match_emissions,
        cfg=step_function_cfg,
        output_dir=output_dir
    )

    logger.info('Calculating modifications frequencies...')
    modification_probs: ModificationsProbs_df = (
        ModificationsProbs_df.from_match_emissions(match_emissions)
    )
    default_mod_freqs: Dict[MODIFICATION_TYPE, Prob] = (
        norine_monomers_info.modification_frequencies()
    )
    default_residue_freqs: Dict[NerpaResidue, Prob] = (
        norine_monomers_info.residue_frequencies()
    )

    logger.info('Calculating PKS probability...')
    unknown_because_pks_prob, (ci_lb, ci_ub) = get_unknown_because_pks_prob(
        match_emissions,
        create_step_function(step_function_steps)
    )
    logger.info(
        f"Estimated PKS UNKNOWN probability: {unknown_because_pks_prob:.4f} "
        f"(95% CI: {ci_lb:.4f}, {ci_ub:.4f})"
    )
    insert_unknown_prob = get_insert_unknown_prob(hmm_chosen_steps_df)

    return EmissionParams(
        step_function=step_function_steps,
        modifications_frequences=(
            list(modification_probs.df.iter_rows(named=True))
        ),
        default_modifications_frequencies=default_mod_freqs,
        default_residue_frequencies=default_residue_freqs,
        unknown_because_pks_prob=unknown_because_pks_prob,
        insert_unknown_prob=insert_unknown_prob
    )
