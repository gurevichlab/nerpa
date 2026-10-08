import math
from functools import partial
from logging import Logger
from typing import (
    Dict,
    NamedTuple,
    Tuple
)
from src.antismash_parsing.genomic_context import (
    ModuleGenomicContextFeature,
    ModuleGenomicContext,
)
from src.general_type_aliases import LogProb, Prob
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType
)
from src.generic.numeric import safe_log
from src.training.hmm_parameters.training_types import (
    HMM_Chosen_Steps_df,
    HMM_Chosen_Steps_Schema
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    STATE_TYPE_TO_EDGE_TYPES,
    RELEVANT_GENOMIC_CONTEXT,
    NORMAL_FLOW_EDGE_TYPES,
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import all_possible_genomic_contexts
from math import log as math_ln
from scipy.stats import beta
from scipy.special import betaln, expit
import polars as pl
from src.training.hmm_parameters.transition_probabilities.estimate_event_prob_given_context import (
    EstimateEventProbGivenContextConfig,
    estimate_event_prob_given_context
)
from src.training.hmm_parameters.transition_probabilities.estimate_skips_at_start_probs import estimate_skips_at_start_probs
from src.training.hmm_parameters.transition_probabilities.estimate_skips_at_end_probs import estimate_skips_at_end_probs
from src.training.hmm_parameters.transition_probabilities.validate import transition_weights_errors
from src.training.hmm_parameters.short_aliases import *
from src.training.hmm_parameters.training_types import (
    SuccessCnts,
    Estimated_Transition_Probs,
    Estimated_Transition_Probs_Schema,
    ProbabilityEstimate,
)

def estimate_transition_probs_for_state(
        state_type: DetailedHMMStateType,
        hmm_chosen_steps_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger
) -> Estimated_Transition_Probs:
    logger.info(f"Estimating transition probabilities for state type: {state_type.name}")
    edge_types = list(STATE_TYPE_TO_EDGE_TYPES[state_type])
    if len(edge_types) == 0:
        return Estimated_Transition_Probs()
    if len(edge_types) == 1:
        return Estimated_Transition_Probs.from_dict({
            state_type: {
                gc: { edge_types[0]: 1.0 }
                for gc in all_possible_genomic_contexts(state_type)
            }
        })

    assert len(edge_types) == 2, \
        "each state should have at most two outgoing edges"

    if (edge_types[0] in NORMAL_FLOW_EDGE_TYPES
        and edge_types[1] not in NORMAL_FLOW_EDGE_TYPES):
        normal_flow_edge_t = edge_types[0]
        abnormal_flow_edge_t = edge_types[1]
    elif (edge_types[1] in NORMAL_FLOW_EDGE_TYPES
          and edge_types[0] not in NORMAL_FLOW_EDGE_TYPES):
        normal_flow_edge_t = edge_types[1]
        abnormal_flow_edge_t = edge_types[0]
    else:
        raise ValueError(
            "each state should have one normal flow edge, "
            "and one abnormal flow edge"
        )

    edge_counts_df = (
        hmm_chosen_steps_df.df
        .filter(pl.col(SC.STATE_TYPE) == state_type.name)
        .with_columns(
            normal_flow_edge_chosen=(
                (pl.col(SC.EDGE_TYPE) == normal_flow_edge_t.name)
                .cast(pl.Int64)
            )
        )
        .group_by(SC.RELEVANT_GC)
        .agg(
            normal_flow_edge_chosen_count=pl.sum("normal_flow_edge_chosen"),
            total_count=pl.count()
        )
        .with_columns(
            abnormal_flow_edge_chosen_count=pl.col("total_count") - pl.col("normal_flow_edge_chosen_count")
        )
    )

    edge_counts_by_gc: Dict[GC, SuccessCnts] = {
        SC.RELEVANT_GC.from_polars_dtype(row[SC.RELEVANT_GC]):
        SuccessCnts(
            failures=row["normal_flow_edge_chosen_count"],
            successes=row["abnormal_flow_edge_chosen_count"]
        )
        for row in edge_counts_df.iter_rows(named=True)
    }

    abnormal_edge_probs_by_gc: Dict[GC, ProbabilityEstimate] = (
        estimate_event_prob_given_context(
            success_cnts=edge_counts_by_gc,
            cfg=cfg,
            logger=logger,
            event_label=state_type.name
        )
    )

    # Convert to a polars table for easier inspection and logging
    ETPS = Estimated_Transition_Probs_Schema
    rows = []
    for gc in all_possible_genomic_contexts(state_type):
        prob_est = abnormal_edge_probs_by_gc[gc]
        rows.append({
            ETPS.STATE_TYPE: state_type,
            ETPS.RELEVANT_GC: gc,
            ETPS.EDGE_TYPE: abnormal_flow_edge_t,
            ETPS.PROB: prob_est.prob,
            ETPS.COMMENTS: prob_est.comments,
        })
        rows.append({
            ETPS.STATE_TYPE: state_type,
            ETPS.RELEVANT_GC: gc,
            ETPS.EDGE_TYPE: normal_flow_edge_t,
            ETPS.PROB: 1 - prob_est.prob,
            ETPS.COMMENTS: prob_est.comments,
        })

    return Estimated_Transition_Probs.from_rows(rows)


def estimate_hmm_transition_probs(
        hmm_chosen_steps_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger,
) -> Estimated_Transition_Probs:
    """
    Infer log probabilities of HMM transitions from the chosen steps data.

    Normally the transition probabilities depend only on the state type and the genomic context.
    However, there are two special cases: skipping at the start and skipping at the end.

    I allow a BGC to start/end only at modules with the
    ASSEMBLY_LINE_START_MARKER/ASSEMBLY_LINE_FINISH_MARKER, respectively,
    so skipping at the start/end means that all modules until the next marker are skipped.

    I want the probability of skipping to depend on
    the number and genomic context of the modules skipped,
    so I have to treat this events specially.
    """

    transition_probs: Estimated_Transition_Probs = Estimated_Transition_Probs()

    # Special cases
    transition_probs.update(
        estimate_skips_at_start_probs(
            hmm_chosen_steps_df,
            cfg,
            logger
        )
    )

    transition_probs.update(
        estimate_skips_at_end_probs(
            hmm_chosen_steps_df,
            cfg,
            logger
        )
    )

    # General case
    for state_type in STATE_TYPE_TO_EDGE_TYPES:
        state_processed: bool = (
            transition_probs.df
            .get_column(SC.STATE_TYPE)
            .eq(state_type.name)
            .any()
        )
        if not state_processed:
            transition_probs.update(
                estimate_transition_probs_for_state(
                    state_type,
                    hmm_chosen_steps_df,
                    cfg,
                    logger
                )
            )


    if (errors := transition_weights_errors(transition_probs.to_dict())):
        for error_msg in errors:
            logger.warning(error_msg)
        raise ValueError("Couldn't estimate transition probabilities")

    return transition_probs
            



