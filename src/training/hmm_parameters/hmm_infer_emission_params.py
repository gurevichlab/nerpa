from typing import Dict, List, Tuple, NamedTuple
from collections import Counter, defaultdict

from src.aa_specificity_prediction_model.specificity_prediction_helper import create_step_function
from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.training.hmm_parameters.step_function import (
    create_bins,
    fit_step_function_to_bins,
    plot_step_function,
    plot_step_function_stacked,
)
from src.monomer_names_helper import NRP_Monomer, UNKNOWN_RESIDUE, MonomerNamesHelper, MonomerResidue
from src.data_types import LogProb, Prob
from src.data_types import (
    BGC_Module,
    BGC_Module_Modification
)
from src.training.hmm_parameters.training_types import MatchEmissionInfo
from src.training.hmm_parameters.norine_stats import NorineStats

from pathlib import Path
from math import log
import numpy as np
import matplotlib.pyplot as plt


def get_score_correctness(emissions: List[MatchEmissionInfo]) -> List[Tuple[LogProb, bool]]:
    score_correctness = []
    for bgc_id, bgc_module, nrp_monomer in emissions:
        for predicted_residue, score in bgc_module.residue_score.items():
            score_correctness.append((score, predicted_residue == nrp_monomer.residue))
            # TODO: check that PKS-hybrids are handled correctly
    return score_correctness

def fit_step_function(emissions: List[MatchEmissionInfo],
                      num_bins: int,
                      step_range: int,
                      output_dir: Path) -> List[float]:
    score_correctness = get_score_correctness(emissions)
    score_correctness_bins = create_bins(score_correctness, num_bins)
    step_function = fit_step_function_to_bins(score_correctness_bins, step_range)
    plot_step_function(score_correctness_bins, step_function,
                       out_file=output_dir / Path('step_function.png'))
    plot_step_function_stacked(score_correctness_bins, step_function,
                               out_file=output_dir / Path('step_function_stacked.png'))
    return step_function


def get_modifications_frequencies(emissions: List[Tuple[BGC_Module, NRP_Monomer]]) -> Dict[str, Dict[str, float]]:
    # (METHYLATION/EPIMERIZATION) -> (BGC_{True/False}_NRP_{True/False} -> frequency)
    mt_cnt = Counter()
    ep_cnt = Counter()
    for bgc_module, nrp_monomer in emissions:
        bgc_mods = bgc_module.modifications
        bgc_meth = BGC_Module_Modification.METHYLATION in bgc_mods
        bgc_chir = BGC_Module_Modification.EPIMERIZATION in bgc_mods

        nrp_meth = nrp_monomer.methylated
        nrp_chir = nrp_monomer.chirality

        mt_cnt[f'BGC_{bgc_meth}_NRP_{nrp_meth}'] += 1
        ep_cnt[f'BGC_{bgc_chir}_NRP_{nrp_chir.name}'] += 1

    result = {}
    # Laplace rule of succession (add 1 to numerator and 2 to denominator)
    result['METHYLATION'] = \
        {'BGC_True_NRP_True': (mt_cnt['BGC_True_NRP_True'] + 1) / (mt_cnt['BGC_True_NRP_True'] + mt_cnt['BGC_True_NRP_False'] + 2),
         'BGC_False_NRP_True': (mt_cnt['BGC_False_NRP_True'] + 1) / (mt_cnt['BGC_False_NRP_True'] + mt_cnt['BGC_False_NRP_False'] + 2)}

    result['METHYLATION']['BGC_True_NRP_False'] = 1 - result['METHYLATION']['BGC_True_NRP_True']
    result['METHYLATION']['BGC_False_NRP_False'] = 1 - result['METHYLATION']['BGC_False_NRP_True']

    # TODO: check if this is correct (some amino acids are stereosymmetric, etc)
    result['EPIMERIZATION'] = \
        {'BGC_True_NRP_True': (ep_cnt['BGC_True_NRP_D'] + 1) / (ep_cnt['BGC_True_NRP_D'] + ep_cnt['BGC_True_NRP_L'] + 2),
         'BGC_False_NRP_True': (ep_cnt['BGC_False_NRP_D'] + 1) / (ep_cnt['BGC_False_NRP_D'] + ep_cnt['BGC_False_NRP_L'] + 2)}

    result['EPIMERIZATION']['BGC_True_NRP_False'] = 1 - result['EPIMERIZATION']['BGC_True_NRP_True']
    result['EPIMERIZATION']['BGC_False_NRP_False'] = 1 - result['EPIMERIZATION']['BGC_False_NRP_True']

    return result


def get_modifications_scores(emissions: List[Tuple[BGC_Module, NRP_Monomer]],
                             default_frequencies: Dict[str, float]) -> Dict[str, float]:
    mod_freqs = get_modifications_frequencies(emissions)

    def log_odds(p1, p2, inverse=False):
        return log(1 - p1) - log(1 - p2) if inverse else log(p1) - log(p2)

    def nrp_label(mod, nrp):
        if mod == 'METHYLATION':
            return str(nrp)
        else:
            return 'D' if nrp else 'L'

    return {f'Mod_{mod}_BGC_{bgc}_NRP_{nrp_label(mod, nrp)}':
                log_odds(mod_freqs[mod][f'BGC_{bgc}'], default_frequencies[mod], inverse=not nrp)
            for mod in ('METHYLATION', 'EPIMERIZATION')
            for bgc in (True, False)
            for nrp in (True, False)}


def get_unknown_because_pks_prob(emissions: List[MatchEmissionInfo],
                                 step_function_steps: List[float]) -> Prob:
    emissions_with_pks = [emission for emission in emissions
                          if ModuleGenomicContextFeature.PKS_DOWNSTREAM in emission.bgc_module.genomic_context]
    calibration_function = create_step_function(step_function_steps)
    unknown_preds_correctness =[(calibration_function(emission.bgc_module.residue_score[UNKNOWN_RESIDUE]),
                                 emission.nrp_monomer.residue == UNKNOWN_RESIDUE)
                                for emission in emissions_with_pks]
    print(f'Unknown residue predictions correctness: {unknown_preds_correctness}')

    # in our model UNKNOWN_RESIDUE can appear due to two independent reasons:
    # 1. the unknown residue is attracted by A domain with probability p
    # 2. the residue modified by PKS domains nearby with probability x
    # here we estimate x given p-s

    probs = np.array([pred for pred, correct in unknown_preds_correctness])
    successes = np.array([correct for pred, correct in unknown_preds_correctness], dtype=int)

    # 2. Prior hyper-parameters  (here: uniform Beta(1,1))  ─────────────────────
    alpha = 1.0  # shape parameter α
    beta = 1.0  # shape parameter β

    # 3. Log-posterior kernel as a Python function ─────────────────────────────-
    def logposterior(x, probs, successes, alpha, beta):
        """
        Returns log π(x | data) up to a normalising constant.
        ------------------------------------------------------
        x : scalar in (0,1)
        probs : array of p_i values
        successes: array of e_i (0/1) values
        alpha, beta : Beta(a,b) prior hyper-parameters
        """
        # Likelihood part: product over i of  Bernoulli(q_i)
        q = probs + (1.0 - probs) * x
        ll = (successes * np.log(q) + (1.0 - successes) * np.log1p(-q)).sum()

        # Prior part: Beta(a,b)  ∝ x^{a-1} (1-x)^{b-1}
        lp = (alpha - 1.0) * np.log(x) + (beta - 1.0) * np.log1p(-x)

        return ll + lp

    # 4. Build an x-grid and evaluate the log-posterior everywhere ──────────────
    #    A 2000-point grid is already plenty for 3–4 decimal places.
    x_grid = np.linspace(0.0005, 0.9995, 2000)
    log_post = np.array([logposterior(x, probs, successes, alpha, beta)
                         for x in x_grid])

    # 5. Convert log-values to normalised weights  w_j  on the grid ─────────────
    #    • subtract max(log_post)  → prevents numerical underflow
    #    • exponentiate            → back to probability scale
    #    • divide by sum           → now ∑ w_j = 1  (proper mass function)
    weights = np.exp(log_post - log_post.max())
    weights /= weights.sum()

    # 6. Posterior summaries  (mean & 95 % central credible interval) ──────────
    mean_x = np.sum(weights * x_grid)

    # cumulative weights for the empirical CDF on the grid
    cdf = np.cumsum(weights)
    lo = x_grid[np.searchsorted(cdf, 0.025)]  # lower 2.5 %
    hi = x_grid[np.searchsorted(cdf, 0.975)]  # upper 97.5 %

    print(f"Posterior mean      = {mean_x:0.4f}")
    print(f"95 % credible band  = ({lo:0.4f}, {hi:0.4f})")

    # 7. (Nice-to-have)  plot the posterior density  ───────────────────────────
    plt.figure()
    # convert the discrete masses to a “density” by dividing by grid spacing
    plt.plot(x_grid, weights / (x_grid[1] - x_grid[0]))
    plt.axvline(mean_x, linestyle="--")  # mark the posterior mean
    plt.xlabel("x")
    plt.ylabel("posterior density")
    plt.title("Posterior of x via grid integration")
    plt.tight_layout()
    #plt.show()

    return mean_x


class EmissionParams(NamedTuple):
    step_function: List[float]
    modifications_frequences: Dict[str, Dict[str, Prob]]  # (modification_type -> (BGC_{True/False}_NRP_{True/False} -> frequency))
    unknown_because_pks_prob: Prob
    default_modifications_frequencies: Dict[str, Prob]  # default frequencies for modifications, used for score normalization
    default_residue_frequencies: Dict[str, Prob]  # default probabilities for residues, used for score normalization


def infer_emission_params(emissions: List[MatchEmissionInfo],
                          norine_stats: NorineStats,
                          output_dir: Path,
                          monomer_names_helper: MonomerNamesHelper) -> EmissionParams:
    print('Building step function...')
    step_function = fit_step_function(emissions, 20, 1000,
                                      output_dir)  # TODO: put in config

    print('Calculating modifications frequencies...')
    emissions_ = [(bgc_module, nrp_monomer.to_base_mon())
                  for bgc_info, bgc_module, nrp_monomer in emissions]
    #modifications_scores = get_modifications_scores(emissions_, default_freqs)
    modification_frequencies = get_modifications_frequencies(emissions_)

    default_mod_freqs = {'METHYLATION': norine_stats.methylated / norine_stats.total_monomers,
                         'EPIMERIZATION': norine_stats.d_chirality / norine_stats.total_monomers}
    default_residue_freqs: Dict[MonomerResidue, Prob] = defaultdict(float)
    for norine_residue, freq in norine_stats.residue_frequencies.items():
        nerpa_residue = monomer_names_helper.parsed_name(norine_residue, 'norine').residue
        default_residue_freqs[nerpa_residue] += freq

    default_residue_freqs = dict(default_residue_freqs)  # to dump as YAML

    print('Calculating PKS probability...')
    unknown_because_pks_prob = get_unknown_because_pks_prob(emissions, step_function)

    return EmissionParams(step_function=step_function,
                          modifications_frequences=modification_frequencies,
                          default_modifications_frequencies=default_mod_freqs,
                          default_residue_frequencies=default_residue_freqs,
                          unknown_because_pks_prob=unknown_because_pks_prob)
