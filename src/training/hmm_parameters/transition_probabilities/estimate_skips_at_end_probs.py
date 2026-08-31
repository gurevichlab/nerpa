from logging import Logger
from typing import (
    Dict,
)
from src.antismash_parsing.genomic_context import (
    ModuleGenomicContextFeature,
    ModuleGenomicContext
)
from src.general_type_aliases import LogProb, Prob
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    all_possible_genomic_contexts
)
from src.training.hmm_parameters.training_types import HMM_Chosen_Steps_df, HMM_Chosen_Steps_Schema
from src.training.hmm_parameters.training_types import SuccessCnts
from src.training.hmm_parameters.transition_probabilities.estimate_event_prob_given_context import (
    EstimateEventProbGivenContextConfig,
    estimate_event_prob_given_context_conditioned
)
from src.training.hmm_parameters.transition_probabilities.estimate_skips_at_start_probs import estimate_prob_non_marker_skipping
from src.training.hmm_parameters.transition_probabilities.estimate_transition_probs_conditioned import (
    estimate_transition_probs_for_state_conditioned,
)

from math import log
from scipy.stats import beta
from scipy.special import betaln, expit
import polars as pl

MCF = ModuleGenomicContextFeature
GC = ModuleGenomicContext
ST = DetailedHMMStateType
ET = DetailedHMMEdgeType
SC = HMM_Chosen_Steps_Schema


def estimate_probs_finish_matching(
        edge_choices_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger
) -> Dict[GC, Dict[ET, Prob]]:
    transition_probs: Dict[GC, Dict[ET, Prob]] = {}

    # 1. With the finish marker
    # standard estimation conditioned on the presence of the finish marker in the genomic context
    transition_probs.update(
        estimate_transition_probs_for_state_conditioned(
            state_type=ST.CHOOSE_IF_FINISH_MATCHING,
            conditioning_feature=MCF.ASSEMBLY_LINE_FINISH_MARKER,
            include_conditioning_feature_in_context=True,
            edge_choices_df=edge_choices_df,
            cfg=cfg,
            logger=logger
        )
    )

    # 1.2. Without the finish marker
    # finishing is disallowed; the only possible edge is CONTINUE_MATCHING
    transition_probs.update(
        {
            gc: {
                ET.FINISH_MATCHING: 0.0,
                ET.CONTINUE_MATCHING: 1.0
            }
            for gc in all_possible_genomic_contexts(ST.CHOOSE_IF_FINISH_MATCHING)
            if MCF.ASSEMBLY_LINE_FINISH_MARKER not in gc
        }
    ) 

    return transition_probs


def estimate_probs_skip_at_end_with_marker(
        edge_choices_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger
) -> Dict[GC, Dict[ET, Prob]]:
    # 1. Define auxiliary columns for filtering edge_choices_df
    aux_columns: Dict[str, pl.Expr] = {}
    aux_columns["finish_marker_present"] = (
            pl.col(SC.GENOMIC_CONTEXT)
            .map_elements(
                lambda gc: MCF.ASSEMBLY_LINE_FINISH_MARKER in gc,
                return_dtype=pl.Boolean
            )
    )
    aux_columns["matching_finished"] = (
            (pl.col(SC.STATE_TYPE) == ST.CHOOSE_IF_FINISH_MATCHING.name)
            & (pl.col(SC.EDGE_TYPE) == ET.FINISH_MATCHING.name)
    )
    aux_columns["skip_at_end"] = (
            (pl.col(SC.STATE_TYPE) == ST.SKIPPING_MODULES_AT_END.name)
            & (pl.col(SC.EDGE_TYPE) == ET.SKIP_MODULES_AT_END.name)
    )

    # 2. Gather relevant transitions:
    # (ST.SKIPPING_MODULES_AT_END --> ET.SKIP_MODULES_AT_END) as positive observations
    # (ST.CHOOSE_IF_FINISH_MATCHING -> ET.FINISH_MATCHING) as negative observations
    marker_decisions_df: pl.DataFrame = (
        edge_choices_df.df
        .with_columns(**aux_columns)
        .filter(
            pl.col("finish_marker_present")
            & (pl.col("matching_finished") | pl.col("skip_at_end"))
        )
    )

    marker_skip_counts_df: pl.DataFrame = (
        marker_decisions_df
        .group_by(SC.GENOMIC_CONTEXT)
        .agg(
            skips_cnt=pl.sum("skip_at_end"),
            finish_matching_cnt=pl.sum("matching_finished")
        )
    )
    
    marker_skip_cnts_by_gc: Dict[GC, SuccessCnts] = {}
    for row in marker_skip_counts_df.iter_rows(named=True):
        gc = row[SC.GENOMIC_CONTEXT]
        marker_skip_cnts_by_gc[gc] = SuccessCnts(
            successes=row["skips_cnt"],
            failures=row["finish_matching_cnt"]
        )

    # 3. Estimate the probability of skipping marker modules at the end,
    # conditioned on the presence of the finish marker in the genomic context
    marker_skip_probs_by_gc: Dict[GC, Prob] = (
        estimate_event_prob_given_context_conditioned(
            success_cnts=marker_skip_cnts_by_gc,
            conditioning_feature=MCF.ASSEMBLY_LINE_FINISH_MARKER,
            include_conditioning_feature_in_context=True,
            cfg=cfg,
            logger=logger,
        )
    )

    return {
        gc: {ET.SKIP_MODULES_AT_END: skip_prob}
        for gc, skip_prob in marker_skip_probs_by_gc.items()
    }


def estimate_skips_at_end_probs(
        edge_choices_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger
) -> Dict[ST, Dict[GC, Dict[ET, Prob]]]:
    """
    Estimate probabilities associated with skipping modules at the end
    of an assembly line.

    Marker and non-marker modules play different roles.

    At a module carrying ASSEMBLY_LINE_FINISH_MARKER, the model makes the
    actual biological decision: either finish matching at this module or
    continue matching (at least) until the next possible assembly-line finish.

    Say, the sequence of modules ends with "..., G1, g2, g3, G4, g5, G6",
    with G1, G4, G6 being marker modules and g2, g3, g5 being non-marker modules.
    The probability of finishing matching at G1 and thus skipping "g2, g3, G4, g5, G6",
    is split into three parts:
    1. The probability of finishing matching at G1 (the actual biological decision);
       it is estimated in the standard way, conditioned on the presence of the finish marker ST.CHOOSE_IF_FINISH_MATCHING.
    2-3. The length penalty for skipping the suffix "g2, g3, G4, g5, G6".
    To compute the penalty, we view the sequence in reverse: "G6, g5, G4, g3, g2, G1, ...".
    We pretend that the model makes a decision at every marker module, whether to skip it or to start (i.e. "finish" in forward view) matching.
    Here, G6 and G4 are skipped, while G1 is not.
    2. We estimate the probability of skipping marker modules at end using this reverse view, with
    skipped modules like G6 and G4 be positive observations, and the module where matching is started (G1) being a negative observation.
    3. For non-marker modules we calculate a length penalty based on the geometric distribution,
    treating each run (e.g. "g2, g3") as a sequence of Bernoulli trials with probability q of the next module being another non-marker module.
    """
    transition_probs: Dict[ST, Dict[GC, Dict[ET, Prob]]] = {}

    # 1. ST.CHOOSE_IF_FINISH_MATCHING
    # choose whether to finish matching; forward view
    transition_probs[ST.CHOOSE_IF_FINISH_MATCHING] = estimate_probs_finish_matching(
        edge_choices_df=edge_choices_df,
        cfg=cfg,
        logger=logger
    )

    # 2. ST.SKIPPING_MODULES_AT_END with the finish marker
    # choose whether to skip or "start matching"; reverse view
    transition_probs[ST.SKIPPING_MODULES_AT_END] = estimate_probs_skip_at_end_with_marker(
        edge_choices_df=edge_choices_df,
        cfg=cfg,
        logger=logger
    )

    # 3. ST.SKIPPING_MODULES_AT_END without the finish marker
    # a length penalty based on the geometric distribution
    skip_prob_non_marker = estimate_prob_non_marker_skipping(
        edge_choices_df=edge_choices_df,
        state_type=ST.SKIPPING_MODULES_AT_END,
        skip_edge_type=ET.SKIP_MODULES_AT_END,
        marker_feature=MCF.ASSEMBLY_LINE_FINISH_MARKER,
        cfg=cfg,
        logger=logger
    )

    transition_probs[ST.SKIPPING_MODULES_AT_END].update(
        {
            gc: {
                ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START:
                    skip_prob_non_marker,
            }
            for gc in all_possible_genomic_contexts(
                ST.SKIPPING_MODULES_AT_END
            )
            if MCF.ASSEMBLY_LINE_START_MARKER not in gc
        }
    )

    return transition_probs

        
