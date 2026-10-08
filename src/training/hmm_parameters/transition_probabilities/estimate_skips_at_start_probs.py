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
from src.training.hmm_parameters.training_types import (
    HMM_Chosen_Steps_df,
    HMM_Chosen_Steps_Schema,
    Estimated_Transition_Probs,
    Estimated_Transition_Probs_Schema,
)
from src.training.hmm_parameters.transition_probabilities.estimate_event_prob_given_context import (
    EstimateEventProbGivenContextConfig,
    estimate_event_prob_given_context
)
from src.training.hmm_parameters.transition_probabilities.estimate_transition_probs_conditioned import (
    estimate_transition_probs_for_state_conditioned,
) 
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    all_possible_genomic_contexts
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

def estimate_prob_non_marker_skipping(
        edge_choices_df: HMM_Chosen_Steps_df,
        state_type: ST,
        skip_edge_type: ET,
        marker_feature: MCF,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger
) -> Prob:
    """
    There's a chain of skipped modules: G1, g2, g3, G4, g5, G6, ...,
    with G1, G4, G6, ... being marker modules and g2, g3, g5, ... being
    Here, we estimate the probability that the next module in the chain is a non-marker module

    We model the number K of such consecutive non-marker modules with a
    geometric distribution:

        P(K = k) = (1 - q) * q**k,

    Consequently:

      * every skipped non-marker module contributes one Bernoulli
        observation with probability q;
      * every completed run contributes one observation with probability
        1 - q.

    Each skipped marker module starts exactly one such run. Therefore the
    number of skipped marker modules is the number of observed runs, while
    the total number of skipped non-marker modules is the total number of
    q-events across those runs.

    A Beta(1, 1) prior is used for q. The returned value is its posterior
    mean. Multiplying this value once for every skipped non-marker module
    gives q**k, i.e. the length-dependent part of the geometric penalty.

    Note that the complete probability of a run of exactly k non-marker
    modules is (1 - q) * q**k; the terminating (1 - q) factor must therefore
    be accounted for separately if the full geometric probability is needed.
    """

    modules_skipped_df: pl.DataFrame = (
        edge_choices_df.df
        .filter(
            (pl.col(SC.STATE_TYPE) == state_type.name)
            & (pl.col(SC.EDGE_TYPE) == skip_edge_type.name)
        )
        .with_columns(
            pl.lit(marker_feature.name)
            .is_in(pl.col(SC.GENOMIC_CONTEXT))
            .alias("marker_present")
        )
    )

    # Total number of non-marker modules appearing inside the skipped runs.
    #
    # Each such module contributes one factor q to the geometric likelihood.
    non_marker_skips: int = (
        modules_skipped_df
        .filter(~pl.col("marker_present"))
        .height
    )

    # Every skipped marker module starts one run of zero or more following
    # non-marker modules. Hence this is also the number of geometric runs
    # observed in the training data.
    num_skip_runs: int = (
        modules_skipped_df
        .filter(pl.col("marker_present"))
        .height
    )

    # Uniform Beta prior for
    #
    #     q = P(next module is a non-marker).
    #
    # The likelihood contributed by all runs is proportional to
    #
    #     q**non_marker_skips * (1 - q)**num_skip_runs,
    #
    # giving the posterior
    #
    #     q | data ~ Beta(
    #         alpha + non_marker_skips,
    #         beta  + num_skip_runs
    #     ).
    alpha_prior = 1.0
    beta_prior = 1.0

    posterior_mean = (
        (alpha_prior + non_marker_skips)
        / (
            alpha_prior
            + beta_prior
            + non_marker_skips
            + num_skip_runs
        )
    )

    return posterior_mean
    

def estimate_skips_at_start_probs(
        edge_choices_df: HMM_Chosen_Steps_df,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger
) -> Estimated_Transition_Probs:
    """
    Estimate probabilities associated with skipping modules at the beginning
    of an assembly line.

    Marker and non-marker modules play different roles.

    At a module carrying ASSEMBLY_LINE_START_MARKER, the model makes the
    actual biological decision: either start matching at this module or
    skip it and continue looking for a later possible assembly-line start.

    If a marker module is skipped, all immediately following non-marker
    modules must also be skipped. These are therefore not independent
    skip/start decisions. Instead, they contribute a geometric length
    penalty: longer stretches of skipped non-marker modules are progressively
    less likely.

    The marker-skip probabilities and the non-marker length penalty are
    estimated separately.
    """

    # 1) Marker modules.
    #
    # At a marker module there is a genuine choice between:
    #
    #   * starting matching here, and
    #   * skipping this candidate start and continuing to the next one.
    #
    # Estimate this probability conditioned on the genomic context of the
    # 1. marker module.
    logger.info(
        "Estimating probabilities of skipping assembly-line "
        "for modules with start markers."
    )
    transition_probs: Estimated_Transition_Probs = (
        estimate_transition_probs_for_state_conditioned(
            state_type=ST.CHOOSE_IF_START_MATCHING,
            conditioning_feature=MCF.ASSEMBLY_LINE_START_MARKER,
            edge_choices_df=edge_choices_df,
            cfg=cfg,
            logger=logger
        )
    )

    # 2) Non-marker modules.
    #
    # Once a marker module has been skipped, any following non-marker
    # modules must also be skipped. Their contribution is therefore a
    # length penalty rather than an independent skip/start decision.
    #
    # I treat all non-marker modules as equivalent, regardless of their genomic context, due to the lack of data.
    #
    # If q is the estimated probability of extending the run by one
    # additional non-marker module, a run containing k such modules
    # contributes q**k.
    logger.info(
        "Estimating probabilities of skipping assembly-line "
        "for modules without start markers."
    )
        
    skip_prob_non_marker: Prob = estimate_prob_non_marker_skipping(
        edge_choices_df=edge_choices_df,
        state_type=ST.CHOOSE_IF_START_MATCHING,
        skip_edge_type=ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START,
        marker_feature=MCF.ASSEMBLY_LINE_START_MARKER,
        cfg=cfg,
        logger=logger
    )

    nonmarker_probs: Dict[GC, Dict[ET, Prob]] = (
        {
            gc: {
                ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START:
                    skip_prob_non_marker,

                # START_MATCHING at a non-marker module is forbidden by the model.
                # When actually assigning weights, the code will use skip_prob_non_marker
                # to add a length penalty to the edge of the corresponding marker module.
                # the weights at non-marker modules will be set to 1.0 and 0.0, respectively.
                ET.START_MATCHING:
                    1.0 - skip_prob_non_marker,
            }
            for gc in all_possible_genomic_contexts(
                ST.CHOOSE_IF_START_MATCHING
            )
            if MCF.ASSEMBLY_LINE_START_MARKER not in gc
        }
    )

    transition_probs.update(
        Estimated_Transition_Probs.from_dict(
            {ST.CHOOSE_IF_START_MATCHING: nonmarker_probs}
        )
    )
    return transition_probs

        
