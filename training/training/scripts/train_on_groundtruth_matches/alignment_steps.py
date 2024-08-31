from typing import Dict, List, Optional, Tuple, NamedTuple
from enum import Enum, auto
from src.data_types import LogProb
from src.monomer_names_helper import MonomerResidue
from src.matching.alignment_types import AlignmentStepType
import yaml


AlignmentStepDict = dict
BGC_VariantDict = dict
MatchDict = dict


def specificity_predictions_for_step(step: AlignmentStepDict,
                                     bgc_variant: BGC_VariantDict) -> Optional[Dict[MonomerResidue, LogProb]]:
    if step['Alignment_step'] != 'MATCH':
        return None
    return next(module['residue_score']
                for module in bgc_variant['tentative_assembly_line']
                if module['gene_id'] == step['Gene'] and module['module_idx'] == step['A-domain_idx'])


class StepLocation(Enum):
    AT_END = auto()
    BETWEEN_GENES = auto()
    INSIDE_GENE = auto()


def enum_representer(dumper, e: Enum):
    return dumper.represent_scalar(f'!{e.__class__.__name__}', e.name)


yaml.add_representer(StepLocation, enum_representer)
yaml.add_representer(AlignmentStepType, enum_representer)


def get_step_location(step_idx: int, alignment: List[AlignmentStepDict]) -> StepLocation:
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


# steps with the same StepKey are considered the same and not counted twice
class StepKey(NamedTuple):
    step_type: AlignmentStepType
    gene: Optional[str] = None
    a_domain_idx: Optional[int] = None
    modifying_domains: Optional[str] = None
    nrp_chirality: Optional[str] = None
    nrp_modifications: Optional[str] = None
    nrp_residue: Optional[str] = None
    dist_to_previous_match_step: Optional[int] = None
    dist_to_next_match_step: Optional[int] = None


def get_step_key(step_idx: int, alignment: List[dict]) -> StepKey:
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

    match step['Alignment_step']:
        case 'MATCH':
            return StepKey(step_type=AlignmentStepType.MATCH,
                           gene=step['Gene'],
                           a_domain_idx=step['A-domain_idx'],
                           modifying_domains=step['Modifying_domains'],
                           nrp_chirality=step['NRP_chirality'],
                           nrp_modifications=step['NRP_modifications'],
                           nrp_residue=step['NRP_residue'])
        case 'NRP_MONOMER_SKIP':
            if previous_match_step_idx is None:
                return StepKey(step_type=AlignmentStepType.NRP_MONOMER_SKIP,
                               dist_to_previous_match_step=None,
                               dist_to_next_match_step=next_match_step_idx - step_idx)
            else:
                previous_match_step = alignment[previous_match_step_idx]
                return StepKey(step_type=AlignmentStepType.NRP_MONOMER_SKIP,
                               gene=previous_match_step['Gene'],
                               a_domain_idx=previous_match_step['A-domain_idx'],
                               dist_to_previous_match_step=step_idx - previous_match_step_idx)
        case 'BGC_MODULE_SKIP':
            return StepKey(step_type=AlignmentStepType.BGC_MODULE_SKIP,
                           gene=step['Gene'],
                           a_domain_idx=step['A-domain_idx'])


class AlignmentStepInfo(NamedTuple):
    specificity_predictions: Optional[Dict[MonomerResidue, LogProb]]
    step_info: AlignmentStepDict
    step_location: StepLocation


def get_steps_info(matches_with_variants_for_bgc: List[Tuple[MatchDict, BGC_VariantDict]]) -> List[AlignmentStepInfo]:
    matching_steps_info: List[AlignmentStepInfo] = []
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
                matching_steps_info.append(AlignmentStepInfo(specificity_predictions, step, step_location))
    return matching_steps_info
