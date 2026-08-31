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

MCF = ModuleGenomicContextFeature
GC = ModuleGenomicContext
ST = DetailedHMMStateType
ET = DetailedHMMEdgeType
SC = HMM_Chosen_Steps_Schema

def estimate_transition_probs_for_state(
        state_type: DetailedHMMStateType,
        hmm_chosen_steps_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger
) -> Dict[GC, Dict[ET, Prob]]:
    edge_types = list(STATE_TYPE_TO_EDGE_TYPES[state_type])
    if len(edge_types) == 0:
        return {}
    if len(edge_types) == 1:
        return {
            gc: { edge_types[0]: 1.0 }
            for gc in all_possible_genomic_contexts(state_type)
        }

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
        .filter(pl.col(SC.STATE_TYPE) == state_type)
        .with_columns(
            normal_flow_edge_chosen=(pl.col(SC.EDGE_TYPE) == normal_flow_edge_t).cast(pl.Int64)
        )
        .group_by(SC.GENOMIC_CONTEXT)
        .agg(
            normal_flow_edge_chosen_count=pl.sum("normal_flow_edge_chosen"),
            total_count=pl.count()
        )
        .with_columns(
            abnormal_flow_edge_chosen_count=pl.col("total_count") - pl.col("normal_flow_edge_chosen_count")
        )
    )

    edge_counts_by_gc: Dict[GenomicContext, SuccessCnts] = {
        row[SC.GENOMIC_CONTEXT]: SuccessCnts(
            failures=row["normal_flow_edge_chosen_count"],
            successes=row["abnormal_flow_edge_chosen_count"]
        )
        for row in edge_counts_df.iter_rows(named=True)
    }

    abnormal_edge_probs_by_gc: Dict[GC, Prob] = (
        estimate_event_prob_given_context(
            success_cnts=edge_counts_by_gc,
            cfg=cfg,
            logger=logger,
            state_type=state_type
        )
    )

    return {
        gc: {
            normal_flow_edge_t: 1 - abnormal_edge_probs_by_gc[gc],
            abnormal_flow_edge_t: abnormal_edge_probs_by_gc[gc]
        }
        for gc in all_possible_genomic_contexts(state_type)
    }



def estimate_hmm_transition_probs(
        hmm_chosen_steps_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger,
) -> Dict[ST, Dict[GC, Dict[ET, Prob]]]:
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

    transtition_probs: Dict[ST, Dict[GC, Dict[ET, Prob]]] = {}

    # Special cases
    transtition_probs.update(
        estimate_skips_at_start_probs(
            hmm_chosen_steps_df,
            cfg,
            logger
        )
    )

    transtition_probs.update(
        estimate_skips_at_end_probs(
            hmm_chosen_steps_df,
            cfg,
            logger
        )
    )

    # General case
    for state_type in STATE_TYPE_TO_EDGE_TYPES:
        if state_type not in transtition_probs:
            transtition_probs[state_type] = estimate_transition_probs_for_state(
                state_type,
                hmm_chosen_steps_df,
                cfg,
                logger
            )

    if (errors := transition_weights_errors(transtition_probs)):
        for error_msg in errors:
            logger.warning(error_msg)
        raise ValueError("Couldn't estimate transition probabilities")

    return transtition_probs
            

def estimate_hmm_transition_logprobs(
        hmm_chosen_steps_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger,
) -> Dict[ST, Dict[GC, Dict[ET, LogProb]]]:
    transition_probs = estimate_hmm_transition_probs(
        hmm_chosen_steps_df,
        cfg,
        logger
    )
    transition_logprobs: Dict[ST, Dict[GC, Dict[ET, LogProb]]] = {
        state_type: {
            gc: {
                edge_type: math_ln(prob)
                for edge_type, prob in edge_probs.items()
            }
            for gc, edge_probs in gc_edge_probs.items()
        }
        for state_type, gc_edge_probs in transition_probs.items()
    }
    return transition_logprobs


