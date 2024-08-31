from typing import Dict, List, NamedTuple, Tuple, Optional
from collections import Counter, defaultdict
from enum import Enum, auto
from pathlib import Path
import yaml
from math import e
from training.training.scripts.train_on_groundtruth_matches.fit_increasing_step_function import fit_step_function_to_bins


def specificity_predictions_for_step(step, bgc_variant) -> Optional[Dict[str, float]]:
    if step['Alignment_step'] != 'MATCH':
        return None
    return next(module['residue_score']
                for module in bgc_variant['tentative_assembly_line']
                if
                module['gene_id'] == step['Gene'] and module['module_idx'] == step['A-domain_idx'])


class StepLocation(Enum):
    AT_END = auto()
    BETWEEN_GENES = auto()
    INSIDE_GENE = auto()


class AlignmentStepType(Enum):  # TODO: import from src/matching/alignment_types.py
    MATCH = auto()
    NRP_MONOMER_SKIP = auto()
    BGC_MODULE_SKIP = auto()
    ITERATE_MODULE = auto()
    ITERATE_GENE = auto()

def enum_representer(dumper, e: Enum):
    return dumper.represent_scalar(f'!{e.__class__.__name__}', e.name)


yaml.add_representer(StepLocation, enum_representer)
yaml.add_representer(AlignmentStepType, enum_representer)

def get_step_location(step_idx: int, alignment: List[dict]) -> StepLocation:
    previous_match_step = next((alignment[i]
                                for i in range(step_idx - 1, -1, -1)
                                if alignment[i]['Alignment_step'] == 'MATCH'),
                               None)
    next_match_step = next((alignment[i]
                            for i in range(step_idx + 1, len(alignment))
                            if alignment[i]['Alignment_step'] == 'MATCH'),
                           None)
    if previous_match_step is None or next_match_step is None:
        return StepLocation.AT_END
    if any([previous_match_step['Gene'] == next_match_step['Gene'],  # for nrp_monomer_skip
            previous_match_step['Gene'] == alignment[step_idx]['Gene'],  # for bgc_module_skip
            next_match_step['Gene'] == alignment[step_idx]['Gene']]):  # for bgc_module_skip
        return StepLocation.INSIDE_GENE
    return StepLocation.BETWEEN_GENES


def get_step_key(step_idx: int, alignment: List[dict]) -> tuple:
    '''
    steps matching the same A-domain to the same NRP monomer are considered the same and not counted twice
    to avoid overfitting to large BGC families
    '''
    previous_match_step_idx = next((i for i in range(step_idx - 1, -1, -1)
                                    if alignment[i]['Alignment_step'] == 'MATCH'),
                                   None)
    next_match_step_idx = next((i for i in range(step_idx + 1, len(alignment))
                                if alignment[i]['Alignment_step'] == 'MATCH'),
                           None)
    step = alignment[step_idx]
    if step['Alignment_step'] == 'MATCH':
        return (AlignmentStepType.MATCH,
                step['Gene'],
                step['A-domain_idx'],
                step['Modifying_domains'],
                step['NRP_chirality'],
                step['NRP_modifications'],
                step['NRP_residue'])
    elif step['Alignment_step'] == 'NRP_MONOMER_SKIP':
        if previous_match_step_idx is None:
            return (AlignmentStepType.NRP_MONOMER_SKIP, None, next_match_step_idx - step_idx)
        else:
            previous_match_step = alignment[previous_match_step_idx]
            return (previous_match_step['Gene'],
                    previous_match_step['A-domain_idx'],
                    step_idx - previous_match_step_idx)
    elif step['Alignment_step'] == 'BGC_MODULE_SKIP':
        return (AlignmentStepType.BGC_MODULE_SKIP, step['Gene'], step['A-domain_idx'])


class MathingStepInfo(NamedTuple):
    specificity_predictions: Optional[Dict[str, float]]
    step_info: dict
    step_location: StepLocation


def get_steps_info_bgc(matches_with_variants_for_bgc: List[Tuple[dict, dict]]) -> List[MathingStepInfo]:
    matching_steps_info: List[MathingStepInfo] = []
    processed_steps = set()
    for match, bgc_variant in matches_with_variants_for_bgc:
        for alignment in match['Alignments']:
            for i, step in enumerate(alignment):
                step_key = get_step_key(i, alignment)
                if step_key in processed_steps:
                    continue
                processed_steps.add(step_key)

                step_location = get_step_location(i, alignment)
                specificity_predictions = specificity_predictions_for_step(step, bgc_variant)
                matching_steps_info.append(MathingStepInfo(specificity_predictions, step, step_location))
    return matching_steps_info


def get_mismatched_pairs_with_perfect_prediction(matching_steps_info: List[MathingStepInfo]) -> Dict[Tuple[str, str], int]:
    mismatched_pairs = Counter()
    for specificity_predictions, step_info, step_location in matching_steps_info:
        if step_info['Alignment_step'] != 'MATCH':
            continue
        if specificity_predictions[step_info['NRP_residue']] == 0.0:
            continue
        for predicted_residue, score in specificity_predictions.items():
            if score == 0.0:
                mismatched_pairs[(predicted_residue, step_info['NRP_residue'])] += 1
    return mismatched_pairs


def get_score_correctness(matching_steps_info: List[MathingStepInfo]) -> List[Tuple[float, bool]]:
    score_correcntess = []
    for specificity_predicion, step_info, step_location in matching_steps_info:
        if step_info['Alignment_step'] != 'MATCH':
            continue
        for predicted_residue, score in specificity_predicion.items():
            score_correcntess.append((score, predicted_residue == step_info['NRP_residue']))
    return score_correcntess


def create_bins(score_correctness_log: List[Tuple[float, bool]],
                num_bins: int) -> List[Tuple[int, int]]:
    score_correctness = [(e ** score, correctness) for score, correctness in score_correctness_log]
    bins: List[List[bool]] = [[] for _ in range(num_bins)]
    bin_size = 1 / num_bins
    for score, correctness in score_correctness:
        bin_idx = min(int(score / bin_size), num_bins - 1)
        bins[bin_idx].append(correctness)
    return [(bin.count(False), bin.count(True)) for bin in bins]


def plot_step_function(score_correctness_bins: List[Tuple[int, int]],
                       step_function: List[float]) -> None:
    assert len(score_correctness_bins) == len(step_function)
    import matplotlib.pyplot as plt
    bin_size = 1 / len(score_correctness_bins)
    bin_centers = [i * bin_size + bin_size / 2 for i in range(len(score_correctness_bins))]
    bin_fractions = [0 if not num_true else num_true / (num_true + num_false)
                     for num_false, num_true in score_correctness_bins]
    plt.bar(bin_centers, bin_fractions, align='center', width=bin_size)
    plt.plot(bin_centers, step_function, color='red')
    plt.show()


def fit_step_function(matching_steps_info: List[MathingStepInfo],
                      num_bins: int,
                      step_range: int) -> List[float]:
    score_correctness = get_score_correctness(matching_steps_info)
    score_correctness_bins = create_bins(score_correctness, num_bins)
    step_function = fit_step_function_to_bins(score_correctness_bins, step_range)
    print('Plotting...')
    plot_step_function(score_correctness_bins, step_function)
    return step_function


def get_indel_frequencies(matching_steps_info: List[MathingStepInfo]) -> Dict[AlignmentStepType, Dict[StepLocation, float]]:
    indel_cnt: Dict[Tuple[AlignmentStepType, StepLocation], int] = Counter()
    for specificity_predicion, step_info, step_location in matching_steps_info:
        if step_info['Alignment_step'] in ('NRP_MONOMER_SKIP', 'BGC_MODULE_SKIP'):
            indel_cnt[AlignmentStepType[step_info['Alignment_step']], step_location] += 1

    steps_total = len(matching_steps_info)
    return {step_type: {location: indel_cnt[step_type, location] / steps_total
                        for location in StepLocation}
            for step_type in (AlignmentStepType.BGC_MODULE_SKIP, AlignmentStepType.NRP_MONOMER_SKIP)}


def get_modifications_frequencies(matching_steps_info: List[MathingStepInfo]) -> Dict[str, Dict[str, float]]:
    modifications_cnt: Dict[Tuple[str, str], int] = Counter()
    for specificity_predicion, step_info, step_location in matching_steps_info:
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

    steps_total = len(matching_steps_info)
    result = {'METHYLATION': {}, 'EPIMERIZATION': {}}
    for (mod_type, bgc_nrp_presense), cnt in modifications_cnt.items():
        result[mod_type][bgc_nrp_presense] = cnt / steps_total
    return result


def load_matches_with_bgc_variants(matches_file: Path, nerpa_results_dir: Path) -> List[Tuple[dict, dict]]:
    matches = yaml.safe_load(matches_file.open('r'))
    bgc_variants = defaultdict(list)
    for bgc_dir in nerpa_results_dir.iterdir():
        bgc_id = bgc_dir.name
        try:
            yaml_files = [f for f in (bgc_dir / 'BGC_variants').iterdir() if f.name.endswith('.yaml')]
        except FileNotFoundError:
            print(f'No BGC variants for {bgc_id}')
            continue
        for yaml_file in yaml_files:
            bgc_variants[bgc_id].extend(yaml.safe_load(yaml_file.read_text()))

    matches_with_variants = []
    for match in matches:
        bgc_id = match['NRP'].split('.')[0]
        bgc_variant = next(variant for variant in bgc_variants[bgc_id]
                           if variant['variant_idx'] == match['BGC_variant_idx'])
        matches_with_variants.append((match, bgc_variant))
    return matches_with_variants


def write_yaml(data, out_file: Path):  # TODO: import from src/write_results.py
    # dirty hack to erase information about types and make output less verbose
    # https://github.com/yaml/pyyaml/issues/408
    yaml.emitter.Emitter.prepare_tag = lambda self, tag: ''

    # another hack (albeit less dirty) to forbid yaml creating references
    # https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml
    yaml.Dumper.ignore_aliases = lambda *args: True

    with open(out_file, 'w') as out:
        yaml.dump(data, out,
                  default_flow_style=None, sort_keys=False)

if __name__ == '__main__':
    nerpa_results_dir = Path('/home/ilianolhin/git/nerpa2/training/training/nerpa_results/currently_inspecting/individual')
    matches_yaml = Path('/home/ilianolhin/git/nerpa2/training/training/matches_inspection_results/approved_matches.yaml')

    print('Loading matches with BGC variants...')
    matches_with_variants = load_matches_with_bgc_variants(matches_yaml, nerpa_results_dir)

    print('Extracting unique alignment steps...')
    matching_steps_info = get_steps_info_bgc(matches_with_variants)

    results = {}
    print('Writing mismatched pairs with perfect prediction...')
    results['mismatched_pairs_with_perfect_prediction'] = get_mismatched_pairs_with_perfect_prediction(matching_steps_info)
    results['mismatched_pairs_with_perfect_prediction'] = sorted(results['mismatched_pairs_with_perfect_prediction'].items(), key=lambda x: x[1], reverse=True)
    print('Building step function...')
    results['step_function'] = fit_step_function(matching_steps_info, 20, 1000)
    print('Calculating indel frequencies...')
    results['indel_frequencies'] = get_indel_frequencies(matching_steps_info)
    print('Calculating modifications frequencies...')
    results['modifications_frequencies'] = get_modifications_frequencies(matching_steps_info)
    print('Done!')

    write_yaml(results, Path('/home/ilianolhin/git/nerpa2/training/training/matches_inspection_results/calculated_statistics.yaml'))
