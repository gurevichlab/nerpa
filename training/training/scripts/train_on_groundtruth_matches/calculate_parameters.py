from typing import Dict, List, NamedTuple, Tuple, Optional
from collections import Counter, defaultdict
from step_function import create_bins, fit_step_function_to_bins, plot_step_function
from src.monomer_names_helper import MonomerResidue
from src.data_types import LogProb, ModuleLocFeature, ModuleLocFeatures
from src.matching.alignment_types import AlignmentStepType
from src.matching.heuristic_matching import HeuristicMatchingConfig
from heuristics_thresholds import calculate_heuristic_parameters
from alignment_steps import AlignmentStepInfo, get_steps_info, MatchDict, BGC_VariantDict, StepType
from indel_params import (
    get_indel_frequencies,
    get_indel_frequencies_cnt,
    get_large_indel_frequencies,
    BGC_Fragment_Loc_Features,
    GeneLocFeatures
)
from dataclasses import asdict, dataclass
from pathlib import Path
from math import log

NRP_VariantDict = dict

def get_mismatched_pairs_with_perfect_prediction(alignment_steps_info: List[AlignmentStepInfo]) \
        -> Dict[Tuple[MonomerResidue, MonomerResidue], int]:  # (predicted, actual) -> count
    mismatched_pairs = Counter()
    for specificity_predictions, step_info, step_location in alignment_steps_info:
        if step_info['Alignment_step'] != 'MATCH':
            continue
        if specificity_predictions[step_info['NRP_residue']] == 0.0:
            continue
        for predicted_residue, score in specificity_predictions.items():
            if score == 0.0:
                mismatched_pairs[(predicted_residue, step_info['NRP_residue'])] += 1
    return mismatched_pairs


def get_score_correctness(alignment_steps_info: List[AlignmentStepInfo]) -> List[Tuple[LogProb, bool]]:
    score_correcntess = []
    for specificity_predicion, step_info, step_location in alignment_steps_info:
        if step_info['Alignment_step'] != 'MATCH':
            continue
        for predicted_residue, score in specificity_predicion.items():
            score_correcntess.append((score, predicted_residue == step_info['NRP_residue']))
    return score_correcntess


def fit_step_function(matching_steps_info: List[AlignmentStepInfo],
                      num_bins: int,
                      step_range: int,
                      output_dir: Path) -> List[float]:
    score_correctness = get_score_correctness(matching_steps_info)
    score_correctness_bins = create_bins(score_correctness, num_bins)
    step_function = fit_step_function_to_bins(score_correctness_bins, step_range)
    plot_step_function(score_correctness_bins, step_function, output_dir)
    return step_function


def get_modifications_frequencies(alignment_steps_info: List[AlignmentStepInfo]) -> Dict[str, Dict[str, float]]:
    # (METHYLATION/EPIMERIZATION) -> (BGC_{True/False}_NRP_{True/False} -> frequency)
    mt_cnt = Counter()
    ep_cnt = Counter()
    for specificity_predicion, step_info, step_location in alignment_steps_info:
        if step_info['Alignment_step'] != 'MATCH':
            continue
        bgc_mods = step_info['Modifying_domains'].split(',')
        bgc_meth = 'METHYLATION' in bgc_mods
        bgc_chir = 'EPIMERIZATION' in bgc_mods

        nrp_mods = step_info['NRP_modifications'].split(',')
        nrp_meth = 'METHYLATION' in nrp_mods
        nrp_chir = step_info['NRP_chirality']

        mt_cnt[f'BGC_{bgc_meth}_NRP_{nrp_meth}'] += 1
        ep_cnt[f'BGC_{bgc_chir}_NRP_{nrp_chir}'] += 1

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


def get_modifications_scores(alignment_steps_info: List[AlignmentStepInfo],
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
    indel_frequencies: Dict[ModuleLocFeatures, Dict[StepType, float]]
    indel_frequencies_cnt: Dict[ModuleLocFeatures, Dict[StepType, int]]
    indel_frequencies_bgc_fragments: Dict[BGC_Fragment_Loc_Features, Tuple[int, int]]
    indel_frequencies_genes: Dict[GeneLocFeatures, Tuple[int, int]]
    indel_frequencies_nrp_fragments: Tuple[int, int]
    modifications_scores: Dict[str, Dict[str, float]]
    heuristic_matching_cfg: HeuristicMatchingConfig

    def to_dict(self):
        # Custom method to handle the conversion of ModuleLocation keys
        return {
            'mismatched_pairs_with_perfect_prediction': self.mismatched_pairs_with_perfect_prediction,
            'step_function': self.step_function,
            'insert_probabilities': {
                tuple(loc): freq_dict for loc, freq_dict in self.insert_probabilities.items()
            },
            'insert_cnts': {
                tuple(loc): dict(counter) for loc, counter in self.indel_cnts.items()
            },
            'skip_probabilities': {
                tuple(f.name for f in loc): {'MATCH': cnt[0], 'SKIP': cnt[1]}
                for loc, cnt in self.skip_probabilities.items()
            },
            'indel_frequencies_genes': {
                tuple(f.name for f in loc): {'MATCH': cnt[0], 'SKIP': cnt[1]}
                for loc, cnt in self.indel_frequencies_genes.items()
            },
            'indel_frequencies_nrp_fragments': {
                'MATCH': self.indel_frequencies_nrp_fragments[0],
                'SKIP': self.indel_frequencies_nrp_fragments[1]
            },
            'modifications_scores': self.modifications_scores,
            'heuristic_matching_cfg': asdict(self.heuristic_matching_cfg)
        }


def calculate_training_parameters(matches_with_bgcs_nrps: List[Tuple[MatchDict, BGC_VariantDict, NRP_VariantDict]],
                                  norine_stats: dict,
                                  output_dir: Path) -> TrainedParameters:
    print('Extracting unique alignment steps...')
    matches_with_bgcs = [(match, bgc_variant) for match, bgc_variant, _ in matches_with_bgcs_nrps]
    alignment_steps_info = get_steps_info(matches_with_bgcs)

    results = {}
    print('Writing mismatched pairs with perfect prediction...')
    mismatched_pairs_cnts = sorted(get_mismatched_pairs_with_perfect_prediction(alignment_steps_info).items(),
                                   key=lambda x: x[1], reverse=True)
    results['mismatched_pairs_with_perfect_prediction'] = mismatched_pairs_cnts
    print('Building step function...')
    results['step_function'] = fit_step_function(alignment_steps_info, 100, 1000,
                                                 output_dir)  # TODO: put in config
    print('Calculating indel frequencies...')
    results['insert_probabilities'] = get_insert_probabilities(alignment_steps_info)
    results['skip_probabilities'] = get_skip_probabilities(matches_with_bgcs_nrps)

    print('Calculating modifications frequencies...')
    default_freqs = {'METHYLATION': norine_stats['methylated'] / norine_stats['total_monomers'],
                     'EPIMERIZATION': norine_stats['epimerized'] / norine_stats['total_monomers']}
    results['modifications_scores'] = get_modifications_scores(alignment_steps_info, default_freqs)

    print('Creating HeuristicMatchingConfig...')
    results['heuristic_matching_cfg'] = calculate_heuristic_parameters(matches_with_bgcs, output_dir)

    print('Done!')

    return TrainedParameters(**results)
