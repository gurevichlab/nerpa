from typing import Dict, List, NamedTuple, Tuple, Optional
from collections import Counter, defaultdict
from src.training.step_function import (
    create_bins,
    fit_step_function_to_bins,
    plot_step_function,
    plot_step_function_stacked
)
from src.monomer_names_helper import MonomerResidue, NRP_Monomer
from src.data_types import LogProb, ModuleGenomicContextFeature, ModuleGenomicContext
from src.data_types import (
    BGC_Module,
    BGC_Module_Modification,
    NRP_Monomer_Modification,
    Chirality
)
from src.training.training_types import MatchWithBGCNRP, MatchEmissionInfo
from src.training.norine_stats import NorineStats

from dataclasses import asdict, dataclass
from pathlib import Path
from math import log, e



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
        {'BGC_True': (mt_cnt['BGC_True_NRP_True'] + 1) / (mt_cnt['BGC_True_NRP_True'] + mt_cnt['BGC_True_NRP_False'] + 2),
         'BGC_False': (mt_cnt['BGC_False_NRP_True'] + 1) / (mt_cnt['BGC_False_NRP_True'] + mt_cnt['BGC_False_NRP_False'] + 2)}
    # TODO: check if this is correct (some amino acids are stereosymmetric, etc)
    result['EPIMERIZATION'] = \
        {'BGC_True': (ep_cnt['BGC_True_NRP_D'] + 1) / (ep_cnt['BGC_True_NRP_D'] + ep_cnt['BGC_True_NRP_L'] + 2),
         'BGC_False': (ep_cnt['BGC_False_NRP_D'] + 1) / (ep_cnt['BGC_False_NRP_D'] + ep_cnt['BGC_False_NRP_L'] + 2)}
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


def infer_emission_params(emissions: List[MatchEmissionInfo],
                          norine_stats: NorineStats,
                          output_dir: Path) -> Tuple[List[float], Dict[str, float]]:
    results = {}
    print('Building step function...')
    results['step_function'] = fit_step_function(emissions, 20, 1000,
                                                 output_dir)  # TODO: put in config
    print('Calculating modifications frequencies...')
    default_freqs = {'METHYLATION': norine_stats.methylated / norine_stats.total_monomers,
                     'EPIMERIZATION': norine_stats.d_chirality / norine_stats.total_monomers}

    emissions_ = [(bgc_module, nrp_monomer) for bgc_info, bgc_module, nrp_monomer in emissions]
    results['modifications_scores'] = get_modifications_scores(emissions_, default_freqs)

    return results['step_function'], results['modifications_scores']
