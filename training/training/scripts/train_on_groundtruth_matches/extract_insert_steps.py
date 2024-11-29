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
from itertools import groupby, pairwise
from more_itertools import split_at

@dataclass
class MonomerInsertRunInfo:
    nrp_id: str
    inserted_monomers: List[MonomerResidue]
    module_context: ModuleLocFeatures  # the module in the alignment matched to a monomer right before the inserted monomers
    module_info: AlignmentStep_BGC_Module_Info

    def __eq__(self, other):
        return all([
            self.module_info.gene_id == other.module_info.gene_id,
            self.module_info.a_domain_idx == other.module_info.a_domain_idx,
            self.inserted_monomers == other.inserted_monomers
            ])


@dataclass
class MonomerInsertAtStartRunInfo:
    nrp_id: str
    inserted_monomers: List[MonomerResidue]
    module_context: ModuleLocFeatures  # the first module in the alignment matched to a monomer
    module_info: AlignmentStep_BGC_Module_Info

    def __eq__(self, other):
        return all([
            self.module_info.gene_id == other.module_info.gene_id,
            self.module_info.a_domain_idx == other.module_info.a_domain_idx,
            self.inserted_monomers == other.inserted_monomers
        ])


@dataclass
class InsertRunsInfo:
    inserts_at_start: List[MonomerInsertAtStartRunInfo]
    inserts: List[MonomerInsertRunInfo]

    def join(self, other: InsertRunsInfo,
             unique: bool = False) -> InsertRunsInfo:
        inserts_at_start = self.inserts_at_start[:]
        inserts = self.inserts[:]
        if unique:
            for insert in other.inserts_at_start:
                if insert not in inserts_at_start:
                    inserts_at_start.append(insert)
            for insert in other.inserts:
                if insert not in inserts:
                    inserts.append(insert)
        else:
            inserts_at_start.extend(other.inserts_at_start)
            inserts.extend(other.inserts)

        return InsertRunsInfo(
            inserts_at_start=inserts_at_start,
            inserts=inserts
        )


def extract_inserts_info(alignment: List[AlignmentStep],
                         bgc_variant: BGC_Variant,
                         nrp_id: str) -> InsertRunsInfo:
    steps_wo_skips = [step for step in alignment if step.nrp_monomer_info is not None]
    insert_runs = [steps_list
                   for steps in split_at(steps_wo_skips, lambda step: step.step_type == AlignmentStepType.MATCH,
                                         keep_separator=True)
                   if (steps_list := list(steps))]  # to remove empty list in the beginning if the first step is a match
    insert_at_start_steps_info, insert_steps_info = [], []
    if insert_runs[0][0].step_type == AlignmentStepType.NRP_MONOMER_INSERT:
        gene_id = insert_runs[1][0].bgc_module_info.gene_id
        a_domain_idx = insert_runs[1][0].bgc_module_info.a_domain_idx
        module = next(module for module in bgc_variant.modules
                      if module.gene_id == gene_id and module.a_domain_idx == a_domain_idx)
        insert_at_start_steps_info.append(MonomerInsertAtStartRunInfo(
            nrp_id=nrp_id,
            inserted_monomers=[step.nrp_monomer_info.residue
                               for step in insert_runs[0]],
            module_context=module.module_loc,
            module_info=AlignmentStep_BGC_Module_Info.from_bgc_module(module)
        ))
    else:
        gene_id = insert_runs[0][0].bgc_module_info.gene_id
        a_domain_idx = insert_runs[0][0].bgc_module_info.a_domain_idx
        module = next(module for module in bgc_variant.modules
                      if module.gene_id == gene_id and module.a_domain_idx == a_domain_idx)
        insert_at_start_steps_info.append(MonomerInsertAtStartRunInfo(
            nrp_id=nrp_id,
            inserted_monomers=[],
            module_context=module.module_loc,
            module_info=AlignmentStep_BGC_Module_Info.from_bgc_module(module)
        ))

    for steps1, steps2 in pairwise(insert_runs):
        if steps1[0].step_type != AlignmentStepType.MATCH:
            continue
        prev_match_step = steps1[0]
        module = next(module for module in bgc_variant.modules
                      if module.gene_id == prev_match_step.bgc_module_info.gene_id and
                      module.a_domain_idx == prev_match_step.bgc_module_info.a_domain_idx)
        if steps2[0].step_type == AlignmentStepType.NRP_MONOMER_INSERT:
            inserted_monomers = [step.nrp_monomer_info.residue
                                 for step in steps2]
        else:
            inserted_monomers = []
        insert_steps_info.append(MonomerInsertRunInfo(
            nrp_id=nrp_id,
            inserted_monomers=inserted_monomers,
            module_context=module.module_loc,
            module_info=AlignmentStep_BGC_Module_Info.from_bgc_module(module)
        ))

    return InsertRunsInfo(
        inserts_at_start=insert_at_start_steps_info,
        inserts=insert_steps_info
    )

