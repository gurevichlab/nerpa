from typing import Dict, List, Optional, Tuple, NamedTuple
from enum import Enum, auto

import dacite

from src.data_types import LogProb, ModuleLocFeature, ModuleLocFeatures
from src.monomer_names_helper import MonomerResidue
from src.matching.matching_types_alignment_step import AlignmentStepType
import yaml


AlignmentStepDict = dict
BGC_VariantDict = dict
BGC_ModuleDict = dict
MatchDict = dict


def module_for_step(step: AlignmentStepDict,
                    bgc_variant: BGC_VariantDict) -> Optional[BGC_ModuleDict]:
    if step['Alignment_step'] == 'NRP_MONOMER_SKIP':
        return None
    return next(module
                for module in bgc_variant['modules']
                if module['gene_id'] == step['Gene'] and module['a_domain_idx'] == step['A-domain_idx'])

def specificity_predictions_for_step(step: AlignmentStepDict,
                                     bgc_variant: BGC_VariantDict) -> Optional[Dict[MonomerResidue, LogProb]]:
    if step['Alignment_step'] != 'MATCH':
        return None
    else:
        return module_for_step(step, bgc_variant)['residue_score']


def enum_representer(dumper, e: Enum):
    return dumper.represent_scalar(f'!{e.__class__.__name__}', e.name)


class StepType(Enum):
    MATCH = auto()
    INSERTION_BEFORE = auto()
    INSERTION_AFTER = auto()
    BGC_MODULE_SKIP = auto()

yaml.add_representer(AlignmentStepType, enum_representer)
yaml.add_representer(StepType, enum_representer)


def get_mod_loc_to_step_type(step_idx: int,
                             alignment: List[AlignmentStepDict],
                             bgc_variant: BGC_VariantDict) -> List[Tuple[ModuleLocFeatures, StepType]]:
    step_tp = alignment[step_idx]['Alignment_step']
    if step_tp in ('ITERATE_MODULE', 'ITERATE_GENE'):
        return []
    if step_tp in ('MATCH', 'BGC_MODULE_SKIP'):
        module = module_for_step(alignment[step_idx], bgc_variant)
        module_loc = tuple(ModuleLocFeature[loc_feature]
                           for loc_feature in module['module_loc'])
        return [(module_loc, StepType[step_tp])]

    assert step_tp == 'NRP_MONOMER_SKIP'
    previous_match_step = next((alignment[i]
                                for i in range(step_idx - 1, -1, -1)
                                if alignment[i]['Alignment_step'] == 'MATCH'),
                               None)
    next_match_step = next((alignment[i]
                            for i in range(step_idx + 1, len(alignment))
                            if alignment[i]['Alignment_step'] == 'MATCH'),
                           None)
    results = []
    if previous_match_step is not None:
        module = module_for_step(previous_match_step, bgc_variant)
        mod_loc = tuple(ModuleLocFeature[loc_feature]
                        for loc_feature in module['module_loc'])
        results.append((mod_loc, StepType.INSERTION_AFTER))

    if next_match_step is not None:
        module = module_for_step(next_match_step, bgc_variant)
        mod_loc = tuple(ModuleLocFeature[loc_feature]
                        for loc_feature in module['module_loc'])
        results.append((mod_loc, StepType.INSERTION_BEFORE))
    return results


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
                return StepKey(step_type=AlignmentStepType.NRP_MONOMER_INSERT,
                               dist_to_previous_match_step=None,
                               dist_to_next_match_step=next_match_step_idx - step_idx)
            else:
                previous_match_step = alignment[previous_match_step_idx]
                return StepKey(step_type=AlignmentStepType.NRP_MONOMER_INSERT,
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
    mod_loc_to_step_type: List[Tuple[ModuleLocFeatures, StepType]]
    nrp_id: str


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

                mod_loc_to_step_type = get_mod_loc_to_step_type(i, alignment, bgc_variant)
                specificity_predictions = specificity_predictions_for_step(step, bgc_variant)
                matching_steps_info.append(AlignmentStepInfo(specificity_predictions,
                                                             step,
                                                             mod_loc_to_step_type,
                                                             match['NRP']))
    return matching_steps_info
