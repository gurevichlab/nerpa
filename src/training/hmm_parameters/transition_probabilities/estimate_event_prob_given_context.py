import math
from functools import partial
from logging import Logger
from typing import (
    Dict,
    NamedTuple,
    Tuple,
    Optional
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
    SuccessCnts,
    ProbabilityEstimate
)
from src.training.hmm_parameters.statistics_helpers import (
    beta_posterior_mean,
    calc_posterior_hyp_is_true
)
from collections import defaultdict
from more_itertools import powerset

MCF = ModuleGenomicContextFeature
GC = ModuleGenomicContext
ST = DetailedHMMStateType
ET = DetailedHMMEdgeType


class EstimateEventProbGivenContextConfig(NamedTuple):
    prior_gc_has_effect: Prob = 0.5
    prior_p_success_if_gc_has_effect: Tuple[float, float] = (1, 1)  # Beta(alpha, beta)
    min_num_cnts_for_empty_context: int = 50
    min_num_cnts_for_nonempty_context: int = 10
    unclear_effect_prob_range: Tuple[Prob, Prob] = (0.05, 0.95)
    min_allowed_prob: Prob = 0.01
    try_indirectly_estimate_uncertain_effects: bool = True


def estimate_event_prob_given_context_indirectly(
        gc: ModuleGenomicContext,
        success_cnts: Dict[ModuleGenomicContext, SuccessCnts],
        posterior_p_success_with_gc: Dict[ModuleGenomicContext, Prob],
        all_features: set[ModuleGenomicContextFeature],
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger,
) -> Optional[Prob]:
    """
    Ad-hoc crude solution:
    estimate fa by comparing (fa, fb) and (fb) for some other feature fb,
    assumes the model
    P(success | fa_present, fb_present) = 1 - (1 - p_success_fa)^fa_present * (1 - p_success_fb)^fb_present * (1 - p_success_base)
    """
    if not (len(gc) == 1 and len(all_features) == 2):
        return None

    this_f = next(iter(gc))
    other_f = next(iter(all_features - gc))

    # Estimating the effect of fa by comparing (fa, fb) and (fb)
    cnts_other = success_cnts.get(frozenset({other_f}), SuccessCnts(0, 0))
    cnts_both = success_cnts.get(frozenset({this_f, other_f}), SuccessCnts(0, 0))
    p_success_other = posterior_p_success_with_gc[frozenset({other_f})]
    p_success_both = posterior_p_success_with_gc[frozenset({this_f, other_f})]

    if (cnts_other.total() < cfg.min_num_cnts_for_nonempty_context
        or cnts_both.total() < cfg.min_num_cnts_for_nonempty_context):
        logger.warning(
            f'Too little data for indirect estimation of effect of context {this_f.name}:'
            f' counts for {(this_f.name, other_f.name)} = {cnts_both.total()},'
            f' counts for {other_f.name} = {cnts_other.total()}.'
        )
        return None

    if p_success_both < p_success_other:
        logger.warning(
            f'Unexpected probabilities for indirect estimation of effect of context {this_f.name}:'
            f' p_success({(this_f.name, other_f.name)}) = {p_success_both:.3f},'
            f' p_success({other_f.name}) = {p_success_other:.3f}.'
        )
        return None
    if p_success_other > 1 - 1e-3:
        logger.warning(
            f"Cannot estimate effect of context {this_f.name} "
            f"because p_success({other_f.name}) is too close to 1."
        )
        return None

    # Under independence assumption, p_success_fa_fb = 1 - (1 - p_success_fa) * (1 - p_success_fb)
    p_success_this_wo_base = (p_success_both - p_success_other) / (1 - p_success_other)
    p_success_base = posterior_p_success_with_gc[frozenset()]
    p_success_this = 1 - (1 - p_success_this_wo_base) * (1 - p_success_base)
    p_success_this = max(p_success_this, cfg.min_allowed_prob)

    return p_success_this
    

def estimate_event_prob_given_context(
        success_cnts: Dict[ModuleGenomicContext, SuccessCnts],
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger,
        event_label: str = "NA",  # for logging purposes
) -> Dict[ModuleGenomicContext, ProbabilityEstimate]:
    """
    Model:
    each nonempty context has a probability of having any effect at all 
    prior_gc_has_effect ~ Bernoulli(0.5)
    if a gc has effect, then the probability of success is unknown and has a uniform prior Beta(1, 1)

    The empty context is assumed to always have effect, and its probability of success is estimated from the counts for the empty context.
    """
    all_features = set().union(*success_cnts.keys())

    comments_by_gc: Dict[ModuleGenomicContext, list[str]] = defaultdict(list)

    # 1. Estimate the success probability for the empty context (base case)
    base_cnts: SuccessCnts = success_cnts.get(frozenset(), SuccessCnts(0, 0))
    if base_cnts.total() < cfg.min_num_cnts_for_empty_context:
        comments_by_gc[frozenset()].append("too little data")
        logger.warning(
            f'{event_label} '
            'Too little data for base case (empty context): '
            f'{base_cnts.total()} counts. '
            f'Minimum required: {cfg.min_num_cnts_for_empty_context}.'
        )
    if base_cnts.successes == 0:
        comments_by_gc[frozenset()].append("no successes")
        logger.warning(
            f'{event_label}: '
            'No successes for base case (empty context).'
        )

    p_success_base = max(base_cnts.estimate_prob(), cfg.min_allowed_prob)

    # 2. Estimate the probabilities that each context has an effect on the event, and if so, what is the probability of success given that context
    posterior_gc_has_effect = {frozenset(): 1.0}
    posterior_p_success_with_gc = {frozenset(): p_success_base}
    for gc in map(frozenset, powerset(all_features)):
        if not gc:
            continue  # already handled empty context
        
        gc_cnts = success_cnts.get(gc, SuccessCnts(0, 0))
        if gc_cnts.total() < cfg.min_num_cnts_for_nonempty_context:
            comments_by_gc[gc].append("too little data")
            logger.warning(
                f'{event_label}: '
                f'Too little data for context {[f.name for f in gc]}: '
                f'{gc_cnts.total()} counts.'
            )

        posterior_gc_has_effect[gc] = calc_posterior_hyp_is_true(   
            p_success_if_hyp_false=p_success_base,
            cnts=gc_cnts,
            prior_hyp_is_true=cfg.prior_gc_has_effect,  # uniform prior
            prior_p_success_if_hyp_true=cfg.prior_p_success_if_gc_has_effect
        )
        if (cfg.unclear_effect_prob_range[0]
            < posterior_gc_has_effect[gc]
            < cfg.unclear_effect_prob_range[1]):
            comments_by_gc[gc].append("unclear effect")
            logger.warning(
                f'{event_label}: '
                f'Unclear effect of context {[f.name for f in gc]}:'
                f'posterior probability of effect = {posterior_gc_has_effect[gc]:.3f}.'
            )

        p_success_if_gc_has_effect = beta_posterior_mean(
            prior_params=cfg.prior_p_success_if_gc_has_effect,
            cnts=gc_cnts
        )
        posterior_p_success_with_gc[gc] = (
            (1 - posterior_gc_has_effect[gc]) * p_success_base
            + posterior_gc_has_effect[gc] * p_success_if_gc_has_effect
        )
        posterior_p_success_with_gc[gc] = max(posterior_p_success_with_gc[gc], cfg.min_allowed_prob)

    # 3. Try to estimate the probabilities of contexts with uncertain effect using other contexts
    # This assumes that the effects of different contexts are independent, which may not be true in general
    if cfg.try_indirectly_estimate_uncertain_effects:
        for gc in posterior_p_success_with_gc:
            little_data_unclear_effect = (
                success_cnts.get(gc, SuccessCnts(0, 0)).total() < cfg.min_num_cnts_for_nonempty_context
                and cfg.unclear_effect_prob_range[0] < posterior_gc_has_effect[gc] < cfg.unclear_effect_prob_range[1]
            )
            if little_data_unclear_effect:
                logger.info(
                    f'{event_label}: '
                    f'Trying to estimate effect of context {[f.name for f in gc]} indirectly.'
                )
            else:
                continue

            indirect_p_success_estimate = estimate_event_prob_given_context_indirectly(
                gc=gc,
                success_cnts=success_cnts,
                posterior_p_success_with_gc=posterior_p_success_with_gc,
                all_features=all_features,
                cfg=cfg,
                logger=logger
            )
            if indirect_p_success_estimate is not None:
                comments_by_gc[gc].append("estimated indirectly")
                posterior_p_success_with_gc[gc] = indirect_p_success_estimate
            else:
                comments_by_gc[gc].append("could not estimate indirectly")
                logger.warning(
                    f'{event_label}: '
                    f'Could not estimate effect of context {[f.name for f in gc]} indirectly. '
                    f'Using direct estimate with little data and unclear effect.'
                )

    # 4. Merge the probabilities with the comments into a single result
    results: Dict[ModuleGenomicContext, ProbabilityEstimate] = {}
    for gc in posterior_p_success_with_gc:
        results[gc] = ProbabilityEstimate(
            prob=posterior_p_success_with_gc[gc],
            comments=comments_by_gc[gc]
        )

    return results

def gc_to_str_short(gc: GC) -> str:
    return str([f.name for f in gc])

def estimate_event_prob_given_context_conditioned(
        success_cnts: Dict[ModuleGenomicContext, SuccessCnts],
        conditioning_feature: ModuleGenomicContextFeature,
        cfg: EstimateEventProbGivenContextConfig,
        logger: Logger,
        event_label: str = "NA",  # for logging purposes
        include_conditioning_feature_in_context: bool = True
) -> Dict[ModuleGenomicContext, ProbabilityEstimate]:
    logger.debug(
        f"Estimating even prob given {conditioning_feature.name}.\n"
        f"state_type: {event_label}\n"
    )
    logger.debug(
        "successs_cnts:\n" +
        '\n'.join(f"{gc_to_str_short(k)}: {v}" for k, v in success_cnts.items())
    )
    success_cnts_conditioned = {}
    for gc, cnts in success_cnts.items():
        if conditioning_feature in gc:
            gc_wo_cond_feature = gc - {conditioning_feature}
            success_cnts_conditioned[gc_wo_cond_feature] = cnts
        else:
            logger.warning(
                f'{event_label}: '
                f'Context {[f.name for f in gc]} does not contain conditioning feature {conditioning_feature.name}.'
            )
    logger.debug(
        "successs_cnts_conditioned:\n" +
        '\n'.join(f"{gc_to_str_short(k)}: {v}" for k, v in success_cnts_conditioned.items())
    )

    posterior_p_success_with_gc_conditioned: Dict[GC, ProbabilityEstimate] = (
        estimate_event_prob_given_context(
            success_cnts=success_cnts_conditioned,
            cfg=cfg,
            logger=logger,
            event_label=event_label
        )
    )

    if include_conditioning_feature_in_context:
        posterior_p_success_with_gc_conditioned = {
            gc | {conditioning_feature}: p_success
            for gc, p_success in posterior_p_success_with_gc_conditioned.items()
        }

    
    return posterior_p_success_with_gc_conditioned
    
