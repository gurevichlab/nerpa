from typing import Dict, List, NamedTuple, Tuple, Optional
from collections import Counter, defaultdict
from step_function import create_bins, fit_step_function_to_bins, plot_step_function
from src.monomer_names_helper import MonomerResidue
from src.data_types import LogProb
from src.matching.alignment_types import AlignmentStepType
from alignment_steps import AlignmentStepInfo, get_steps_info, StepLocation, MatchDict, BGC_VariantDict
from dataclasses import dataclass
from pathlib import Path


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


def get_indel_frequencies(alignment_steps_info: List[AlignmentStepInfo]) -> Dict[AlignmentStepType, Dict[StepLocation, float]]:
    indel_cnt: Dict[Tuple[AlignmentStepType, StepLocation], int] = Counter()
    for specificity_predicion, step_info, step_location in alignment_steps_info:
        if step_info['Alignment_step'] in ('NRP_MONOMER_SKIP', 'BGC_MODULE_SKIP'):
            indel_cnt[AlignmentStepType[step_info['Alignment_step']], step_location] += 1

    steps_total = len(alignment_steps_info)
    return {step_type: {location: indel_cnt[step_type, location] / steps_total
                        for location in StepLocation}
            for step_type in (AlignmentStepType.BGC_MODULE_SKIP, AlignmentStepType.NRP_MONOMER_SKIP)}


def get_modifications_frequencies(alignment_steps_info: List[AlignmentStepInfo]) -> Dict[str, Dict[str, float]]:
    # (METHYLATION/EPIMERIZATION) -> (BGC_{True/False}_NRP_{True/False} -> frequency)
    modifications_cnt: Dict[Tuple[str, str], int] = Counter()
    for specificity_predicion, step_info, step_location in alignment_steps_info:
        if step_info['Alignment_step'] != 'MATCH':
            continue
        bgc_mods = step_info['Modifying_domains'].split(',')
        bgc_meth = 'METHYLATION' in bgc_mods
        bgc_chir = 'EPIMERIZATION' in bgc_mods

        nrp_mods = step_info['NRP_modifications'].split(',')
        nrp_meth = 'METHYLATION' in nrp_mods
        nrp_chir = step_info['NRP_chirality']

        modifications_cnt['METHYLATION', f'BGC_{bgc_meth}_NRP_{nrp_meth}'] += 1
        modifications_cnt['EPIMERIZATION', f'BGC_{bgc_chir}_NRP_{nrp_chir}'] += 1

    steps_total = len(alignment_steps_info)
    result = {'METHYLATION': {}, 'EPIMERIZATION': {}}
    for (mod_type, bgc_nrp_presense), cnt in modifications_cnt.items():
        result[mod_type][bgc_nrp_presense] = cnt / steps_total
    return result


@dataclass
class TrainedParameters:
    mismatched_pairs_with_perfect_prediction: List[Tuple[Tuple[MonomerResidue, MonomerResidue], int]]  # (predicted, actual) -> count. Dict is displated incorrectly, so list
    step_function: List[float]
    indel_frequencies: Dict[AlignmentStepType, Dict[StepLocation, float]]
    modifications_frequencies: Dict[str, Dict[str, float]]


def calculate_training_parameters(matches_with_variants_for_bgc: List[Tuple[MatchDict, BGC_VariantDict]],
                                  output_dir: Path) -> TrainedParameters:
    print('Extracting unique alignment steps...')
    alignment_steps_info = get_steps_info(matches_with_variants_for_bgc)

    results = {}
    print('Writing mismatched pairs with perfect prediction...')
    mismatched_pairs_cnts = sorted(get_mismatched_pairs_with_perfect_prediction(alignment_steps_info).items(),
                                   key=lambda x: x[1], reverse=True)
    results['mismatched_pairs_with_perfect_prediction'] = mismatched_pairs_cnts
    print('Building step function...')
    results['step_function'] = fit_step_function(alignment_steps_info, 100, 1000,
                                                 output_dir)  # TODO: put in config
    print('Calculating indel frequencies...')
    results['indel_frequencies'] = get_indel_frequencies(alignment_steps_info)
    print('Calculating modifications frequencies...')
    results['modifications_frequencies'] = get_modifications_frequencies(alignment_steps_info)
    print('Done!')

    return TrainedParameters(**results)
