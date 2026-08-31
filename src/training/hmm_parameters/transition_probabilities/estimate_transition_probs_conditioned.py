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
from src.training.hmm_parameters.training_types import SuccessCnts
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

import polars as pl

MCF = ModuleGenomicContextFeature
GC = ModuleGenomicContext
ST = DetailedHMMStateType
ET = DetailedHMMEdgeType
SC = HMM_Chosen_Steps_Schema

def estimate_transition_probs_for_state_conditioned(
        state_type: DetailedHMMStateType,
        conditioning_feature: ModuleGenomicContextFeature,
        edge_choices_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger,
        include_conditioning_feature_in_context: bool = True,
) -> Dict[GC, Dict[ET, Prob]]:
    logger.info(
        f"Estimating transition probabilities for state type {state_type} "
        f"conditioned on genomic context feature {conditioning_feature}"
    )
    # Filter the rows corresponding to the state type state_type
    # and with conditioning_feature in the genomic context
    logger.debug(f"edge_choices_df size before filtering: {edge_choices_df.df.shape}")
    # Debug: filter step-by-step so we can see which condition is dropping rows
    edge_choices_state_only: pl.DataFrame = edge_choices_df.df.filter(
        pl.col(SC.STATE_TYPE) == state_type.name
    )
    logger.debug(
        "edge_choices_df size after filtering by state type: "
        f"{edge_choices_state_only.shape}"
    )

    edge_choices_with_cond_col: pl.DataFrame = edge_choices_state_only.with_columns(
        pl.col(SC.GENOMIC_CONTEXT)
        .map_elements(
            lambda gc: conditioning_feature in gc,
            return_dtype=pl.Boolean,
        )
        .alias("_cond_feature_present")
    )

    edge_choices_with_cond_feature: pl.DataFrame = edge_choices_with_cond_col.filter(
        pl.col("_cond_feature_present")
    )
    logger.debug(
        "edge_choices_df size after filtering by conditioning feature: "
        f"{edge_choices_with_cond_feature.shape}"
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
            total=pl.len(),
        )
        .with_columns(
            normal_edge_cnt=pl.col("total") - pl.col("abnormal_edge_cnt")
        )
    )

    abnormal_choices_by_gc: Dict[GC, SuccessCnts] = {}
    for row in edge_choices_cnts.iter_rows(named=True):
        gc = row[SC.RELEVANT_GC]
        abnormal_choices_by_gc[gc] = SuccessCnts(
            failures=int(row["normal_edge_cnt"]),
            successes=int(row["abnormal_edge_cnt"]),
        )

    abnormal_edge_probs_by_gc: Dict[GC, Prob] = (
        estimate_event_prob_given_context_conditioned(
            success_cnts=abnormal_choices_by_gc,
            conditioning_feature=conditioning_feature,
            include_conditioning_feature_in_context=include_conditioning_feature_in_context,
            state_type=state_type,
            cfg=cfg,
            logger=logger,
        )
    )

    return {
        gc: {
            normal_flow_edge: 1 - prob,
            abnormal_flow_edge: prob
        }
        for gc, prob in abnormal_edge_probs_by_gc.items()
    }
