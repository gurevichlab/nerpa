from typing import (
    List,
    Dict,
    Tuple,
    TypeVar,
    NamedTuple,
    Union
)
from collections import Counter, defaultdict
from src.write_results import write_yaml
from pathlib import Path

import yaml

from src.data_types import ModuleLocFeatures
from src.antismash_parsing.location_features import (
    ModuleLocFeature, ModuleLocFeatures,
    GeneLocFeature, GeneLocFeatures,
    BGC_Fragment_Loc_Feature, BGC_Fragment_Loc_Features
)
from training.training.scripts.train_on_groundtruth_matches.alignment_steps import AlignmentStepInfo, StepType
from alignment_steps import AlignmentStepInfo, StepType
from itertools import groupby
from enum import Enum, auto
from extract_insert_steps import (
    InsertRunsInfo,
    MonomerInsertRunInfo,
    MonomerInsertAtStartRunInfo
)
from extract_matches_skips_steps import (
    MatchesSkipsSteps
)


MatchDict = dict
NRP_VariantDict = dict
BGC_VariantDict = dict

T = TypeVar('T')


def group_by_module_context(insert_runs: Union[List[MonomerInsertRunInfo], List[MonomerInsertAtStartRunInfo]]) \
        -> Dict[ModuleLocFeatures, List[Tuple[int, str]]]:  # module_context -> [(num_inserts, nrp_id)]; nrp_ids are for debug only
    insert_runs_grouped = defaultdict(list)
    for insert_run in insert_runs:
        insert_runs_grouped[insert_run.module_context].append((len(insert_run.inserted_monomers), insert_run.nrp_id))
    return insert_runs_grouped


def get_insert_prob(insert_run_lens: List[int],
                    pseudocounts: bool = False,
                    min_prob: float = 0,
                    max_prob: float = 1) -> Tuple[float, float]:  # start_prob, p_continue
    num_nonempty_runs = sum(run_len > 0 for run_len in insert_run_lens)
    start_prob = num_nonempty_runs / len(insert_run_lens)
    mean_run_len = sum(insert_run_lens) / num_nonempty_runs if num_nonempty_runs > 0 else 0
    p_finish = 1 / (mean_run_len + 1)  # geometric distribution (at each step we have p_finish chance to finish)
    p_continue = 1 - p_finish  # probability to continue inserting
    return start_prob, p_continue


def get_insert_probabilities(inserts_by_module_context: Dict[ModuleLocFeatures, List[Tuple[int, str]]]) \
        -> Dict[ModuleLocFeatures, Tuple[float, float]]:  # module_context -> (start_run, prolong_run)
    '''
    total_runs_for_module_context = {module_context: len(insert_runs)
                                     for module_context, insert_runs in inserts_by_module_context.items()}

    frequent_module_contexts = {module_context for module_context, total_runs in total_runs_for_module_context.items()
                                if total_runs >= 50 or
                                total_runs >= 7 and any(num_inserts > 0 for num_inserts, _ in
                                                        inserts_by_module_context[module_context])}
    '''
    return {module_context: get_insert_prob([num_inserts for num_inserts, _ in insert_runs])
            for module_context, insert_runs in inserts_by_module_context.items()}



def get_insert_probabilities_all(insert_runs: InsertRunsInfo) \
        -> Tuple[Dict[ModuleLocFeatures, Tuple[float, float]], Dict[ModuleLocFeatures, Tuple[float, float]]]:  # insert_after, insert_before
    insert_by_module_context = group_by_module_context(insert_runs.inserts)
    insert_by_module_context_at_start = group_by_module_context(insert_runs.inserts_at_start)


    # dump insert_by_module_context and insert_by_module_context_at_start to yaml
    insert_by_module_context_nonempty = {module_context:
                                             [
                                                 insert_run
                                                 for insert_run in insert_runs
                                                 if insert_run[0] > 0
                                              ]
                                         for module_context, insert_runs in insert_by_module_context.items()}
    insert_by_module_context_at_start_nonempty = {module_context:
                                                        [
                                                            insert_run
                                                            for insert_run in insert_runs
                                                            if insert_run[0] > 0
                                                        ]
                                                    for module_context, insert_runs in insert_by_module_context_at_start.items()}
    write_yaml(insert_by_module_context_nonempty, Path('insert_by_module_context.yaml'))
    write_yaml(insert_by_module_context_at_start_nonempty, Path('insert_by_module_context_at_start.yaml'))

    return get_insert_probabilities(insert_by_module_context), get_insert_probabilities(insert_by_module_context_at_start)


class SkipsProbs(NamedTuple):
    bgc_fragment: Dict[BGC_Fragment_Loc_Features, float]
    gene: Dict[GeneLocFeatures, float]
    module: Dict[ModuleLocFeatures, float]

    def to_dict(self):
        return {
            'bgc_fragment': {tuple(f.name for f in loc_features): prob for loc_features, prob in self.bgc_fragment.items()},
            'gene': {tuple(f.name for f in loc_features): prob for loc_features, prob in self.gene.items()},
            'module': {tuple(f.name for f in loc_features): prob for loc_features, prob in self.module.items()},
        }


def calculate_skip_probability(matches, skips):
    contexts = {step.context for step in matches} | {step.context for step in skips}
    return {
        context: len([step for step in skips if step.context == context]) /
                 (len([step for step in matches if step.context == context]) +
                  len([step for step in skips if step.context == context]))
        for context in contexts
    }


def get_skip_probs(matches_and_skips_steps: MatchesSkipsSteps) -> SkipsProbs:
    bgc_fragment_skip_prob = calculate_skip_probability(
        matches_and_skips_steps.bgc_fragment_matches, matches_and_skips_steps.bgc_fragment_skips
    )
    gene_skip_prob = calculate_skip_probability(
        matches_and_skips_steps.gene_matches, matches_and_skips_steps.gene_skips
    )
    module_skip_prob = calculate_skip_probability(
        matches_and_skips_steps.module_matches, matches_and_skips_steps.module_skips
    )

    return SkipsProbs(bgc_fragment=bgc_fragment_skip_prob, gene=gene_skip_prob, module=module_skip_prob)