from typing import (
    Callable,
    Tuple,
)
from math import log
from scipy.stats import beta
from scipy.special import betaln, expit
from scipy.optimize import brentq
from scipy.stats import chi2
from src.training.hmm_parameters.training_types import SuccessCnts
from src.general_type_aliases import Prob
import math

def beta_posterior_params(
        prior_params: Tuple[float, float],
        cnts: SuccessCnts
) -> Tuple[float, float]:
    """
    Returns the parameters of the posterior Beta distribution given the counts and a Beta prior.
    """
    alpha_post = prior_params[0] + cnts.successes
    beta_post = prior_params[1] + cnts.failures
    return alpha_post, beta_post

def beta_posterior_mean(
        prior_params: Tuple[float, float],
        cnts: SuccessCnts
) -> float:
    """
    Returns the mean of the posterior Beta distribution given the counts and a Beta prior.
    """
    alpha_post, beta_post = beta_posterior_params(prior_params, cnts)
    return alpha_post / (alpha_post + beta_post)


def calc_posterior_hyp_is_true(
        p_success_if_hyp_false: Prob,
        cnts: SuccessCnts,
        prior_hyp_is_true: Prob,
        prior_p_success_if_hyp_true: Tuple[float, float]  # Beta(alpha, beta)
) -> Prob:
    """
    Null hypothesis: event occurs with probability p_success_if_hyp_false
    Hypothesis: event occurs with an unknown probability p_success_if_hyp_true ~ Beta(alpha, beta)
    Hypothesis is true with prior probability prior_hyp_is_true
    Output: posterior probability that hypothesis is true given the counts of successes and failures
    """
    if p_success_if_hyp_false < 1e-3 or p_success_if_hyp_false > 1-1e-3:
        raise ValueError(f'p_success_if_hyp_false too close to 0 or 1: {p_success_if_hyp_false}')
    
    # Calculate the likelihood of the data given the feature is inactive
    log_likelihood_if_hyp_false = (
        cnts.successes * math.log(p_success_if_hyp_false)
        + cnts.failures * math.log(1 - p_success_if_hyp_false)
    )

    # Calculate the likelihood of the data given the feature is active
    log_likelihood_if_hyp_true = (
        betaln(cnts.successes + prior_p_success_if_hyp_true[0],
               cnts.failures + prior_p_success_if_hyp_true[1])
        - betaln(prior_p_success_if_hyp_true[0], prior_p_success_if_hyp_true[1])
    )

    log_prior_odds = (
        math.log(prior_hyp_is_true)
        - math.log1p(-prior_hyp_is_true)
    )

    log_posterior_odds = (
        log_prior_odds
        + log_likelihood_if_hyp_true
        - log_likelihood_if_hyp_false
    )

    return expit(log_posterior_odds)




def get_profile_likelihood_ci(
        log_likelihood: Callable[[float], float],
        mle: float,
        parameter_bounds: Tuple[float, float],
        confidence_level: float = 0.95,
) -> Tuple[float, float]:
    """
    Compute a profile-likelihood confidence interval for a scalar parameter.

    The confidence interval contains parameter values whose log-likelihood
    differs from the maximum by at most

        chi2.ppf(confidence_level, df=1) / 2.

    If the likelihood does not fall below this threshold before reaching a
    parameter boundary, that boundary is used as the corresponding CI limit.
    """
    lower_bound, upper_bound = parameter_bounds

    max_log_likelihood = log_likelihood(mle)
    max_log_likelihood_drop = chi2.ppf(confidence_level, df=1) / 2
    target_log_likelihood = max_log_likelihood - max_log_likelihood_drop

    def ci_equation(parameter: float) -> float:
        return log_likelihood(parameter) - target_log_likelihood

    if ci_equation(lower_bound) >= 0:
        ci_lower = lower_bound
    else:
        ci_lower = brentq(ci_equation, lower_bound, mle)

    if ci_equation(upper_bound) >= 0:
        ci_upper = upper_bound
    else:
        ci_upper = brentq(ci_equation, mle, upper_bound)

    return ci_lower, ci_upper
