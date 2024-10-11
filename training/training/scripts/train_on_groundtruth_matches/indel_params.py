from typing import List, Dict, Tuple, TypeVar, NamedTuple
from collections import Counter, defaultdict
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

MatchDict = dict
NRP_VariantDict = dict
BGC_VariantDict = dict

T = TypeVar('T')


def get_insert_cnts(alignment_steps_info: List[AlignmentStepInfo]) -> Dict[ModuleLocFeatures, Dict[StepType, int]]:
    counter = Counter((mod_loc, step_type)
                      for step_info in alignment_steps_info
                      for mod_loc, step_type in step_info.mod_loc_to_step_type)
    cnt_per_mod_loc = defaultdict(lambda: defaultdict(lambda: 0))
    for (mod_loc, step_type), cnt in counter.items():
        cnt_per_mod_loc[mod_loc][step_type] = cnt
    return cnt_per_mod_loc


def get_insert_probabilities(alignment_steps_info: List[AlignmentStepInfo]) \
        -> Tuple[Dict[ModuleLocFeatures, float], Dict[ModuleLocFeatures, float]]:  # insert_after, insert_before
    insert_cnts = get_insert_cnts(alignment_steps_info)
    frequent_locs = {mod_loc for mod_loc, cnts in insert_cnts.items()
                     if sum(cnts.values()) >= 10 or all(v > 0 for v in cnts.values())}
    rare_locs = {mod_loc for mod_loc in insert_cnts.keys() if mod_loc not in frequent_locs}
    ins_prob_before, ins_prob_after = dict(), dict()
    for loc in frequent_locs:
        bgc_mod_skips = insert_cnts[loc].get(StepType.BGC_MODULE_SKIP, 0)
        insert_after = insert_cnts[loc].get(StepType.INSERTION_AFTER, 0)
        insert_before = insert_cnts[loc].get(StepType.INSERTION_BEFORE, 0)
        matches = insert_cnts[loc].get(StepType.MATCH, 0)
        ins_prob_after[loc] = (insert_after + 1) / (matches + 2)
        if ModuleLocFeature.START_OF_BGC in loc:
            ins_prob_before[loc] = (insert_before + 1) / (matches + 2)

    for loc in rare_locs:
        ins_prob_after = max(ins_prob_after[subloc]
                             for subloc in frequent_locs
                             if all(subloc_feature in loc for subloc_feature in subloc))
        if ModuleLocFeature.START_OF_BGC in loc:
            ins_prob_before[loc] = max(ins_prob_before[subloc]
                                       for subloc in frequent_locs
                                       if all(subloc_feature in loc for subloc_feature in subloc))

    return ins_prob_after, ins_prob_before


def get_gene_to_bgc_fragment(bgc_variant: Dict) -> Dict[str, int]:  # gene_id -> fragment_idx
    gene_to_fragment = dict()
    for fragment_idx, fragment in enumerate(bgc_variant['fragments']):
        for module in fragment:
            gene_to_fragment[module['gene_id']] = fragment_idx
    return gene_to_fragment


def get_rban_idx_to_nrp_fragment(nrp_variant: Dict) -> Dict[int, int]:  # rban_idx -> fragment_idx
    rban_idx_to_fragment = dict()
    for fragment_idx, fragment in enumerate(nrp_variant['fragments']):
        for mon in fragment['monomers']:
            rban_idx_to_fragment[mon['rban_idx']] = fragment_idx
    return rban_idx_to_fragment


def get_bgc_fragment_features(fragment_idx: int, bgc_variant: BGC_VariantDict) -> BGC_Fragment_Loc_Features:
    fst_module = bgc_variant['fragments'][fragment_idx][0]
    lst_module = bgc_variant['fragments'][fragment_idx][-1]
    pairs = [(BGC_Fragment_Loc_Feature.START_OF_BGC, fragment_idx == 0),
             (BGC_Fragment_Loc_Feature.END_OF_BGC, fragment_idx == len(bgc_variant['fragments']) - 1),
             (BGC_Fragment_Loc_Feature.PKS_UPSTREAM, 'PKS_UPSTREAM' in fst_module['module_loc']),
             (BGC_Fragment_Loc_Feature.PKS_DOWNSTREAM, 'PKS_DOWNSTREAM' in lst_module['module_loc'])]
    return tuple(feature for feature, is_present in pairs if is_present)


def get_gene_features(gene_id: str, bgc_variant: BGC_VariantDict) -> GeneLocFeatures:
    fst_module = next(module
                      for fragment in bgc_variant['fragments']
                      for module in fragment
                      if module['gene_id'] == gene_id)
    lst_module = next(module
                        for fragment in reversed(bgc_variant['fragments'])
                        for module in reversed(fragment)
                        if module['gene_id'] == gene_id)
    pairs = [(GeneLocFeature.START_OF_BGC, 'START_OF_BGC' in fst_module['module_loc']),
             (GeneLocFeature.END_OF_BGC, 'END_OF_BGC' in lst_module['module_loc']),
             (GeneLocFeature.START_OF_FRAGMENT, 'START_OF_FRAGMENT' in fst_module['module_loc']),
             (GeneLocFeature.END_OF_FRAGMENT, 'END_OF_FRAGMENT' in lst_module['module_loc']),
             (GeneLocFeature.PKS_UPSTREAM, 'PKS_UPSTREAM' in fst_module['module_loc']),
             (GeneLocFeature.PKS_DOWNSTREAM, 'PKS_DOWNSTREAM' in lst_module['module_loc'])]
    return tuple(feature for feature, is_present in pairs if is_present)


def get_module_features(gene_id: str,
                        a_domain_idx: int,
                        bgc_variant: BGC_VariantDict) -> ModuleLocFeatures:
    module = next(module
                  for fragment in bgc_variant['fragments']
                  for module in fragment
                  if module['gene_id'] == gene_id and module['a_domain_idx'] == a_domain_idx)
    return tuple(ModuleLocFeature[loc_feature] for loc_feature in module['module_loc'])


class MatchesSkipsCnt(NamedTuple):
    matches: int
    skips: int


def sum_pairs_dicts(_d1: Dict[T, MatchesSkipsCnt],
                    _d2: Dict[T, MatchesSkipsCnt]) -> Dict[T, MatchesSkipsCnt]:
    d1 = defaultdict(lambda: MatchesSkipsCnt(0, 0), _d1)
    d2 = defaultdict(lambda: MatchesSkipsCnt(0, 0), _d2)
    joined_keys = set(d1.keys()) | set(d2.keys())
    return {k: MatchesSkipsCnt(d1[k][0] + d2[k][0], d1[k][1] + d2[k][1]) for k in joined_keys}


def get_modules_skips_and_matches(gene_steps: List[Dict],
                                  bgc_variant: BGC_VariantDict) \
        -> Dict[ModuleLocFeatures, MatchesSkipsCnt]:  # module features -> (matches, skips)
    modules_params = defaultdict(lambda: MatchesSkipsCnt(0, 0))
    for step in gene_steps:
        if step['Gene'] == '---':
            continue
        module_features = get_module_features(step['Gene'], step['A-domain_idx'], bgc_variant)
        if step['Alignment_step'] == 'BGC_MODULE_SKIP':
            modules_params[module_features] += MatchesSkipsCnt(matches=0, skips=1)
        else:
            modules_params[module_features] += MatchesSkipsCnt(matches=1, skips=0)

    return modules_params


def get_genes_and_modules_skips_and_matches(fragment_steps: List[Dict],
                                            bgc_variant: BGC_VariantDict) \
        -> Tuple[Dict[GeneLocFeatures, MatchesSkipsCnt], Dict[ModuleLocFeatures, MatchesSkipsCnt]]:
    # gene features -> (matches, skips), module features -> (matches, skips)
    steps_grouped_by_gene = groupby(filter(lambda step: step['Gene'] != '---', fragment_steps),
                                    key=lambda step: step['Gene'])
    genes_params = defaultdict(lambda: MatchesSkipsCnt(0, 0))
    modules_params = defaultdict(lambda: MatchesSkipsCnt(0, 0))
    for gene_id, _gene_steps in steps_grouped_by_gene:
        gene_steps = list(_gene_steps)
        gene_features = get_gene_features(gene_id, bgc_variant)
        if all(step['Alignment_step'] == 'BGC_MODULE_SKIP' for step in gene_steps):
            genes_params[gene_features] += MatchesSkipsCnt(matches=0, skips=1)
        else:
            genes_params[gene_features] += MatchesSkipsCnt(matches=1, skips=0)
            modules_params = sum_pairs_dicts(modules_params, get_modules_skips_and_matches(gene_steps, bgc_variant))
    return genes_params, modules_params


def get_nrp_skips_and_matches(steps: List[Dict],
                                        rban_idx_to_nrp_fragment: Dict[int, int])\
        -> Tuple[MatchesSkipsCnt, MatchesSkipsCnt]:
    steps_grouped_by_fragment = groupby(filter(lambda step: step['rBAN_idx'] != '---', steps),
                                        key=lambda step: rban_idx_to_nrp_fragment[step['rBAN_idx']])
    fragment_matches, fragment_skips = 0, 0
    monomer_matches, monomer_skips = 0, 0
    for fragment_idx, _fragment_steps in steps_grouped_by_fragment:
        fragment_steps = list(_fragment_steps)
        if all(step['Alignment_step'] == 'NRP_MONOMER_SKIP' for step in fragment_steps):
            fragment_skips += 1
        else:
            fragment_matches += 1
            monomer_matches += len([step for step in fragment_steps if step['Alignment_step'] == 'MATCH'])
            monomer_skips += len([step for step in fragment_steps if step['Alignment_step'] == 'NRP_MONOMER_SKIP'])
    return MatchesSkipsCnt(fragment_matches, fragment_skips), MatchesSkipsCnt(monomer_matches, monomer_skips)


class SkipsAndMatchesCnts(NamedTuple):
    bgc_fragments_cnts: Dict[BGC_Fragment_Loc_Features, MatchesSkipsCnt]
    genes_cnts: Dict[GeneLocFeatures, MatchesSkipsCnt]
    modules_cnts: Dict[ModuleLocFeatures, MatchesSkipsCnt]
    nrp_fragments_cnts: MatchesSkipsCnt
    monomers_cnts: MatchesSkipsCnt

    def __add__(self, other):
        return SkipsAndMatchesCnts(
            sum_pairs_dicts(self.bgc_fragments_cnts, other.bgc_fragments_cnts),
            sum_pairs_dicts(self.genes_cnts, other.genes_cnts),
            sum_pairs_dicts(self.modules_cnts, other.modules_cnts),
            MatchesSkipsCnt(self.nrp_fragments_cnts.matches + other.nrp_fragments_cnts.matches,
                            self.nrp_fragments_cnts.skips + other.nrp_fragments_cnts.skips),
            MatchesSkipsCnt(self.monomers_cnts.matches + other.monomers_cnts.matches,
                            self.monomers_cnts.skips + other.monomers_cnts.skips)
        )


def get_skips_and_matches_cnts(steps: List[Dict],
                               bgc_variant: BGC_VariantDict,
                               gene_to_bgc_fragment: Dict[str, int],
                               rban_idx_to_nrp_fragment: Dict[int, int]) -> SkipsAndMatchesCnts:
    steps_grouped_by_fragment = groupby(filter(lambda step: step['Gene'] != '---', steps),
                                        key=lambda step: gene_to_bgc_fragment[step['Gene']])

    fragments_params = defaultdict(lambda: MatchesSkipsCnt(0, 0))
    genes_params = defaultdict(lambda: MatchesSkipsCnt(0, 0))
    modules_params = defaultdict(lambda: MatchesSkipsCnt(0, 0))
    for fragment_idx, _fragment_steps in steps_grouped_by_fragment:
        fragment_steps = list(_fragment_steps)
        fragment_features = get_bgc_fragment_features(fragment_idx, bgc_variant)
        if all(step['Alignment_step'] == 'BGC_MODULE_SKIP' for step in fragment_steps):
            fragments_params[fragment_features] += MatchesSkipsCnt(matches=0, skips=1)
        else:
            fragments_params[fragment_features] += MatchesSkipsCnt(matches=1, skips=0)
            fragment_genes_params, fragment_modules_params = get_genes_and_modules_skips_and_matches(fragment_steps, bgc_variant)
            genes_params = sum_pairs_dicts(genes_params, fragment_genes_params)
            modules_params = sum_pairs_dicts(modules_params, fragment_modules_params)
    return SkipsAndMatchesCnts(fragments_params, genes_params, modules_params,
                               *get_nrp_skips_and_matches(steps, rban_idx_to_nrp_fragment))


def cnts_to_freqs(cnts: Dict[Tuple[T, ...], Tuple[int, int]]) \
    -> Dict[Tuple[T, ...], Tuple[float, float]]:  # (matches, skips) -> (match_freq, skip_freq) for each list of features T
    frequent_features = set(features for features, cnts in cnts.items()
                            if sum(cnts) >= 7)
    result = dict()
    for features in frequent_features:
        matches, skips = cnts[features]
        result[features] = (skips + 1) / (matches + skips + 2)  # Laplace rule of succession

    for features in set(cnts.keys()) - frequent_features:
        result[features] = max(result[subfeatures][1]
                               for subfeatures in frequent_features
                               if all(subfeature in features for subfeature in subfeatures))
    return result


class SkipsProbs(NamedTuple):
    bgc_fragment: Dict[BGC_Fragment_Loc_Features, float]
    gene: Dict[GeneLocFeatures, float]
    module: Dict[ModuleLocFeatures, float]
    nrp_fragment: float
    monomer: float


def get_skip_probs(matches_with_bgcs_nrps: List[Tuple[MatchDict, BGC_VariantDict, NRP_VariantDict]]) -> SkipsProbs:

    skips_matches_counts = SkipsAndMatchesCnts(defaultdict(lambda: MatchesSkipsCnt(0, 0)),
                                               defaultdict(lambda: MatchesSkipsCnt(0, 0)),
                                               defaultdict(lambda: MatchesSkipsCnt(0, 0)),
                                               MatchesSkipsCnt(0, 0),
                                               MatchesSkipsCnt(0, 0))
    for match, bgc_variant, nrp_variant in matches_with_bgcs_nrps:
        gene_to_bgc_fragment = get_gene_to_bgc_fragment(bgc_variant)
        rban_idx_to_nrp_fragment = get_rban_idx_to_nrp_fragment(nrp_variant)
        match_steps = [step
                       for alignment in match['Alignments']
                       for step in alignment
                       if step['Alignment_step'] not in ('ITERATE_GENE', 'ITERATE_MODULE')]
        skips_matches_counts += get_skips_and_matches_cnts(match_steps, bgc_variant, gene_to_bgc_fragment, rban_idx_to_nrp_fragment)
    return SkipsProbs(
        cnts_to_freqs(skips_matches_counts.bgc_fragments_cnts),
        cnts_to_freqs(skips_matches_counts.genes_cnts),
        cnts_to_freqs(skips_matches_counts.modules_cnts),
        (skips_matches_counts.nrp_fragments_cnts.skips + 1) / (skips_matches_counts.nrp_fragments_cnts.matches + skips_matches_counts.nrp_fragments_cnts.skips + 2),
        (skips_matches_counts.monomers_cnts.skips + 1) / (skips_matches_counts.monomers_cnts.matches + skips_matches_counts.monomers_cnts.skips + 2)
    )