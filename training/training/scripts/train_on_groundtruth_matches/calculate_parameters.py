from typing import Dict, List, NamedTuple, Tuple, Optional
from collections import Counter, defaultdict
from step_function import create_bins, fit_step_function_to_bins, plot_step_function
from src.monomer_names_helper import MonomerResidue
from src.data_types import LogProb, ModuleLocFeature, ModuleLocFeatures
from src.matching.matching_types_alignment_step import AlignmentStepType
from src.matching.matching_types_match import Match
from src.data_types import (
    NRP_Variant,
    BGC_Variant,
    BGC_Module_Modification,
    NRP_Monomer_Modification,
   Chirality
)
from auxilary_types import MatchWithBGCNRP
from norine_stats import NorineStats

from src.matching.heuristic_matching import HeuristicMatchingConfig
from heuristics_thresholds import calculate_heuristic_parameters
from indel_params import (
    get_insert_probabilities_all,
    get_skip_probs,
    SkipsProbs,
    BGC_Fragment_Loc_Features,
    GeneLocFeatures
)
from extract_steps import get_steps_info
from extract_matches_skips_steps import ModuleMatchStepInfo
from dataclasses import asdict, dataclass
from pathlib import Path
from math import log
from itertools import groupby


'''
def get_mismatched_pairs_with_perfect_prediction(alignment_steps_info: List[AlignmentStepInfo]) \
        -> Dict[Tuple[MonomerResidue, MonomerResidue], int]:  # (predicted, actual) -> count
    mismatched_pairs = Counter()
    for specificity_predictions, step_info, step_location, _ in alignment_steps_info:
        if step_info['Alignment_step'] != 'MATCH':
            continue
        if specificity_predictions[step_info['NRP_residue']] == 0.0:
            continue
        for predicted_residue, score in specificity_predictions.items():
            if score == 0.0:
                mismatched_pairs[(predicted_residue, step_info['NRP_residue'])] += 1
    return mismatched_pairs
'''

def get_score_correctness(match_steps: List[ModuleMatchStepInfo]) -> List[Tuple[LogProb, bool]]:
    score_correcntess = []
    for match_step in match_steps:
        for predicted_residue, score in match_step.residue_scores.items():
            score_correcntess.append((score, predicted_residue == match_step.step.nrp_monomer_info.residue))
    return score_correcntess


def fit_step_function(match_steps: List[ModuleMatchStepInfo],
                      num_bins: int,
                      step_range: int,
                      output_dir: Path) -> List[float]:
    score_correctness = get_score_correctness(match_steps)
    score_correctness_bins = create_bins(score_correctness, num_bins)
    step_function = fit_step_function_to_bins(score_correctness_bins, step_range)
    plot_step_function(score_correctness_bins, step_function, output_dir)
    return step_function


def get_modifications_frequencies(match_steps_info: List[ModuleMatchStepInfo]) -> Dict[str, Dict[str, float]]:
    # (METHYLATION/EPIMERIZATION) -> (BGC_{True/False}_NRP_{True/False} -> frequency)
    mt_cnt = Counter()
    ep_cnt = Counter()
    for step_info in match_steps_info:
        bgc_mods = step_info.step.bgc_module_info.modifying_domains
        bgc_meth = BGC_Module_Modification.METHYLATION in bgc_mods
        bgc_chir = BGC_Module_Modification.EPIMERIZATION in bgc_mods

        nrp_mods = step_info.step.nrp_monomer_info.modifications
        nrp_meth = NRP_Monomer_Modification.METHYLATION in nrp_mods
        nrp_chir = step_info.step.nrp_monomer_info.chirality

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


def get_modifications_scores(alignment_steps_info: List[ModuleMatchStepInfo],
                             default_frequencies: Dict[str, float]) -> Dict[str, float]:
    mod_freqs = get_modifications_frequencies(alignment_steps_info)

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


@dataclass
class TrainedParameters:
    mismatched_pairs_with_perfect_prediction: List[Tuple[Tuple[MonomerResidue, MonomerResidue], int]]  # (predicted, actual) -> count. Dict is displated incorrectly, so list
    step_function: List[float]
    insert_after_prob: Dict[ModuleLocFeatures, float]
    insert_at_start_prob: Dict[ModuleLocFeatures, float]
    skip_probs: dict  # SkipsProbs.to_dict() for dumping
    modifications_scores: Dict[str, Dict[str, float]]
    heuristic_matching_cfg: HeuristicMatchingConfig


def calculate_training_parameters(matches_with_bgcs_nrps: List[MatchWithBGCNRP],
                                  norine_stats: NorineStats,
                                  output_dir: Path) -> TrainedParameters:
    print('Extracting alignment steps...')
    steps_info = get_steps_info(matches_with_bgcs_nrps)
    match_steps = steps_info.matches_and_skips.module_matches

    results = {}
    '''
    print('Writing mismatched pairs with perfect prediction...')
    mismatched_pairs_cnts = sorted(get_mismatched_pairs_with_perfect_prediction(matching_steps_with_modules).items(),
                                   key=lambda x: x[1], reverse=True)
    results['mismatched_pairs_with_perfect_prediction'] = mismatched_pairs_cnts
    '''
    print('Building step function...')
    results['step_function'] = fit_step_function(match_steps, 100, 1000,
                                                 output_dir)  # TODO: put in config
    print('Calculating indel frequencies...')
    results['insert_after_prob'], results['insert_at_start_prob'] = get_insert_probabilities_all(steps_info.insert_runs)
    results['skip_probs'] = get_skip_probs(steps_info.matches_and_skips).to_dict()

    print('Calculating modifications frequencies...')
    default_freqs = {'METHYLATION': norine_stats.methylated / norine_stats.total_monomers,
                     'EPIMERIZATION': norine_stats.d_chirality / norine_stats.total_monomers}
    results['modifications_scores'] = get_modifications_scores(match_steps, default_freqs)

    '''
    print('Creating HeuristicMatchingConfig...')
    results['heuristic_matching_cfg'] = calculate_heuristic_parameters(matches_with_bgcs, output_dir)
    '''

    print('Done!')

    return TrainedParameters(**results)
