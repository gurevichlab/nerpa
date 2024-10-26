from __future__ import annotations
from typing import (
    Dict,
    Iterable,
    List,
    NamedTuple,
    Tuple,
    Set
)
from dataclasses import dataclass
from src.antismash_parsing.location_features import (
    ModuleLocFeatures,
    ModuleLocFeature,
    GeneLocFeatures,
    GeneLocFeature,
    BGC_Fragment_Loc_Features,
    BGC_Fragment_Loc_Feature,
    module_features_to_fragment_features,
    module_features_to_gene_features
)
from src.matching.matching_types_alignment_step import (
    AlignmentStep,
    AlignmentStepType,
    AlignmentStep_BGC_Module_Info
)
from src.matching.matching_types_match import Match
from src.monomer_names_helper import MonomerResidue
from src.data_types import (
    BGC_Variant,
    GeneId,
    LogProb,
    NRP_Variant
)
from dataclasses import dataclass
from itertools import groupby

@dataclass
class ModuleMatchStepInfo:
    step: AlignmentStep
    nrp_id: str
    bgc_module_context: ModuleLocFeatures
    residue_scores: Dict[MonomerResidue, LogProb]

    def __eq__(self, other):
        return all([
            self.step.bgc_module_info.gene_id == other.step.bgc_module_info.gene_id,
            self.step.bgc_module_info.a_domain_idx == other.step.bgc_module_info.a_domain_idx,
            self.step.bgc_module_info.modifying_domains == other.step.bgc_module_info.modifying_domains,
            self.step.nrp_monomer_info.rban_name == other.step.nrp_monomer_info.rban_name
        ])


@dataclass
class ModuleSkipStepInfo:
    step: AlignmentStep
    nrp_id: str
    bgc_module_context: ModuleLocFeatures

    def __eq__(self, other):
        return all([
            self.step.bgc_module_info.gene_id == other.step.bgc_module_info.gene_id,
            self.step.bgc_module_info.a_domain_idx == other.step.bgc_module_info.a_domain_idx
        ])


@dataclass
class GeneSkipInfo:
    nrp_id: str
    gene_id: GeneId
    bgc_gene_context: GeneLocFeatures

    def __eq__(self, other):
        return self.gene_id == other.gene_id


@dataclass
class GeneMatchInfo:
    nrp_id: str
    gene_id: GeneId
    bgc_gene_context: GeneLocFeatures

    def __eq__(self, other):
        return self.gene_id == other.gene_id

@dataclass
class BGCFragmentSkipInfo:
    nrp_id: str
    fragment_genes: List[GeneId]
    bgc_fragment_context: BGC_Fragment_Loc_Features

    def __eq__(self, other):
        return self.fragment_genes == other.fragment_genes  # account for the order of genes?


@dataclass
class BGCFragmentMatchInfo:
    nrp_id: str
    fragment_genes: List[GeneId]
    bgc_fragment_context: BGC_Fragment_Loc_Features

    def __eq__(self, other):
        return self.fragment_genes == other.fragment_genes  # account for the order of genes?


@dataclass
class MatchesSkipsSteps:
    module_matches: List[ModuleMatchStepInfo]
    module_skips: List[ModuleSkipStepInfo]
    gene_matches: List[GeneMatchInfo]
    gene_skips: List[GeneSkipInfo]
    bgc_fragment_matches: List[BGCFragmentMatchInfo]
    bgc_fragment_skips: List[BGCFragmentSkipInfo]


    def join(self, other: MatchesSkipsSteps,
             unique: bool = False) -> MatchesSkipsSteps:
        joined_data = {}
        for attr in self.__annotations__.keys():
            current_list = getattr(self, attr)[:]
            other_list = getattr(other, attr)
            if unique:
                for item in other_list:
                    if item not in current_list:
                        current_list.append(item)
            else:
                current_list.extend(other_list)
            joined_data[attr] = current_list
        return MatchesSkipsSteps(**joined_data)


def fragments_matches_skips_info(steps_wo_ins: List[AlignmentStep],
                                 bgc_variant: BGC_Variant,
                                 nrp_id: str,
                                 gene_to_fragment: Dict[GeneId, int]) \
    -> Tuple[List[BGCFragmentMatchInfo], List[BGCFragmentSkipInfo]]:
    steps_by_fragment = groupby(steps_wo_ins,
                                key=lambda step: gene_to_fragment[step.bgc_module_info.gene_id])
    fragments_matches, fragments_skips = [], []
    for fragment_idx, fragment_steps in steps_by_fragment:
        fragment_steps = list(fragment_steps)
        fragment_modules = [module for module in bgc_variant.modules
                            if module.fragment_idx == fragment_idx]
        fragment_context = module_features_to_fragment_features(fragment_modules[0].module_loc,
                                                                fragment_modules[-1].module_loc)
        fragment_genes = []
        for module in fragment_modules:
            if module.gene_id not in fragment_genes:
                fragment_genes.append(module.gene_id)

        if all(step.nrp_monomer_info is None for step in fragment_steps):
            fragments_skips.append(BGCFragmentSkipInfo(nrp_id=nrp_id,
                                                       fragment_genes=fragment_genes,
                                                       bgc_fragment_context=fragment_context))
        else:
            fragments_matches.append(BGCFragmentMatchInfo(nrp_id=nrp_id,
                                                          fragment_genes=fragment_genes,
                                                          bgc_fragment_context=fragment_context))
    return fragments_matches, fragments_skips


def genes_matches_skips_info(steps_wo_ins: List[AlignmentStep],
                             bgc_variant: BGC_Variant,
                             nrp_id: str,
                             gene_to_fragment: Dict[GeneId, int],
                             skipped_fragments: Set[int]) -> Tuple[List[GeneMatchInfo], List[GeneSkipInfo]]:
        steps_by_gene = groupby(steps_wo_ins,
                                key=lambda step: step.bgc_module_info.gene_id)
        genes_matches, genes_skips = [], []
        for gene_id, gene_steps in steps_by_gene:
            if gene_to_fragment[gene_id] in skipped_fragments:
                continue
            gene_steps = list(gene_steps)
            gene_context = module_features_to_gene_features(gene_steps[0].bgc_module_info.module_loc,
                                                            gene_steps[-1].bgc_module_info.module_loc)
            if all(step.nrp_monomer_info is None for step in gene_steps):
                genes_skips.append(GeneSkipInfo(nrp_id=nrp_id,
                                                gene_id=gene_id,
                                                bgc_gene_context=gene_context))
            else:
                genes_matches.append(GeneMatchInfo(nrp_id=nrp_id,
                                                   gene_id=gene_id,
                                                   bgc_gene_context=gene_context))
        return genes_matches, genes_skips



def extract_matches_skips_info(alignment: List[AlignmentStep],
                               bgc_variant: BGC_Variant,
                               nrp_id: str) -> MatchesSkipsSteps:
    gene_to_fragment = {module.gene_id: module.fragment_idx
                        for module in bgc_variant.modules}
    steps_wo_ins = [step for step in alignment if step.step_type != AlignmentStepType.NRP_MONOMER_INSERT]
    fragments_matches, fragments_skips = fragments_matches_skips_info(steps_wo_ins, bgc_variant, nrp_id, gene_to_fragment)
    skipped_fragments = {gene_to_fragment[gene_id]
                         for fragment in fragments_skips
                         for gene_id in fragment.fragment_genes}
    genes_matches, genes_skips = genes_matches_skips_info(steps_wo_ins, bgc_variant, nrp_id, gene_to_fragment, skipped_fragments)
    skipped_genes = {skipped_gene.gene_id for skipped_gene in genes_skips}

    module_matches, module_skips = [], []
    for step in steps_wo_ins:
        gene_id = step.bgc_module_info.gene_id
        a_domain_idx = step.bgc_module_info.a_domain_idx
        if gene_id in skipped_genes or gene_to_fragment[gene_id] in skipped_fragments:
            continue

        module = next(module for module in bgc_variant.modules
                      if module.gene_id == gene_id and module.a_domain_idx == a_domain_idx)

        if step.nrp_monomer_info is None:
            module_skips.append(ModuleSkipStepInfo(step=step,
                                                   nrp_id=nrp_id,
                                                   bgc_module_context=module.module_loc))
        else:
            module_matches.append(ModuleMatchStepInfo(step=step,
                                                      nrp_id=nrp_id,
                                                      bgc_module_context=module.module_loc,
                                                      residue_scores=module.residue_score))
    return MatchesSkipsSteps(module_matches=module_matches,
                            module_skips=module_skips,
                            gene_matches=genes_matches,
                            gene_skips=genes_skips,
                            bgc_fragment_matches=fragments_matches,
                            bgc_fragment_skips=fragments_skips)
