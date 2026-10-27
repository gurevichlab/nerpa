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
from src.training.hmm_parameters.training_types import HMM_Chosen_Steps_df, HMM_Chosen_Steps_Schema
from src.training.hmm_parameters.training_types import (
    SuccessCnts,
    Estimated_Transition_Probs,
    Estimated_Transition_Probs_Schema,
)
from src.training.hmm_parameters.transition_probabilities.estimate_event_prob_given_context import (
    EstimateEventProbGivenContextConfig,
    estimate_event_prob_given_context_conditioned
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    STATE_TYPE_TO_EDGE_TYPES,
    RELEVANT_GENOMIC_CONTEXT,
    NORMAL_FLOW_EDGE_TYPES,
    get_normal_flow_edge,
    get_abnormal_flow_edge,
)
from src.training.hmm_parameters.short_aliases import *
from src.training.hmm_parameters.transition_probabilities.estimate_event_prob_given_context import (
    ProbabilityEstimate,
)

import polars as pl

def estimate_transition_probs_for_state_conditioned(
        state_type: DetailedHMMStateType,
        conditioning_feature: ModuleGenomicContextFeature,
        edge_choices_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger,
        include_conditioning_feature_in_context: bool = True,
) -> Estimated_Transition_Probs:
    logger.info(
        f"Estimating transition probabilities for state type {state_type} "
        f"conditioned on genomic context feature {conditioning_feature}"
    )
    edge_choices_state_only: pl.DataFrame = edge_choices_df.df.filter(
        pl.col(SC.STATE_TYPE) == state_type.name
    )

    edge_choices_with_cond_col: pl.DataFrame = edge_choices_state_only.with_columns(
        pl.lit(conditioning_feature.name)
        .is_in(pl.col(SC.GENOMIC_CONTEXT))
        .alias("_cond_feature_present")
    )

    edge_choices_with_cond_feature: pl.DataFrame = edge_choices_with_cond_col.filter(
        pl.col("_cond_feature_present")
    )

    normal_flow_edge: ET = get_normal_flow_edge(state_type)
    abnormal_flow_edge: ET = get_abnormal_flow_edge(state_type)

    # Count how many times the SKIP_MODULE_AT_START edge was (not) chosen for each genomic context
    edge_choices_cnts: pl.DataFrame = (
        edge_choices_with_cond_feature
        .with_columns(
            (pl.col(SC.EDGE_TYPE) == abnormal_flow_edge.name)
            .alias("_abnormal_flow_edge_chosen")
        )
        .group_by(SC.RELEVANT_GC)
        .agg(
            abnormal_edge_cnt=pl.sum("_abnormal_flow_edge_chosen"),
            total=pl.count(),
        )
        .with_columns(
            normal_edge_cnt=pl.col("total") - pl.col("abnormal_edge_cnt")
        )
    )

    abnormal_choices_by_gc: Dict[GC, SuccessCnts] = {}
    for row in edge_choices_cnts.iter_rows(named=True):
        gc = SC.RELEVANT_GC.from_polars_dtype(row[SC.RELEVANT_GC])
        abnormal_choices_by_gc[gc] = SuccessCnts(
            failures=int(row["normal_edge_cnt"]),
            successes=int(row["abnormal_edge_cnt"]),
        )

    abnormal_edge_probs_by_gc: Dict[GC, ProbabilityEstimate] = (
        estimate_event_prob_given_context_conditioned(
            success_cnts=abnormal_choices_by_gc,
            conditioning_feature=conditioning_feature,
            include_conditioning_feature_in_context=include_conditioning_feature_in_context,
            event_label=state_type.name,
            cfg=cfg,
            logger=logger,
        )
    )

    # Convert to a polars table for easier inspection and logging
    ETPS = Estimated_Transition_Probs_Schema
    rows = []
    for gc, prob_est in abnormal_edge_probs_by_gc.items():
        rows.append({
            ETPS.STATE_TYPE: state_type,
            ETPS.RELEVANT_GC: gc,
            ETPS.EDGE_TYPE: abnormal_flow_edge,
            ETPS.PROB: prob_est.prob,
            ETPS.COMMENTS: prob_est.comments,
        })
        rows.append({
            ETPS.STATE_TYPE: state_type,
            ETPS.RELEVANT_GC: gc,
            ETPS.EDGE_TYPE: normal_flow_edge,
            ETPS.PROB: 1 - prob_est.prob,
            ETPS.COMMENTS: prob_est.comments,
        })

    return Estimated_Transition_Probs.from_rows(rows)
