from __future__ import annotations
from typing import (
    Dict,
    NamedTuple,
    List,
    Optional,
    Tuple
)
from src.data_types import (
    BGC_Module,
    BGC_Module_Modification,
    LogProb,
    NRP_Monomer,
    Chirality
)
from src.antismash_parsing.antismash_parser_types import GeneId
from src.monomer_names_helper import MonomerResidue, NRP_Monomer_Modification
from enum import Enum, auto
from collections import OrderedDict
from dataclasses import dataclass


class AlignmentStep_BGC_Module_Info(NamedTuple):
    gene_id: GeneId
    a_domain_idx: int
    top_scoring_residues: List[MonomerResidue]
    modifying_domains: List[BGC_Module_Modification]
    aa10_code: str
    aa34_code: str

    @classmethod
    def from_bgc_module(cls, bgc_module: BGC_Module) -> AlignmentStep_BGC_Module_Info:
        top_score = max(bgc_module.residue_score.values())
        top_residues = [res for res, score in bgc_module.residue_score.items()
                        if score == top_score]
        return cls(gene_id=bgc_module.gene_id,
                   a_domain_idx=bgc_module.a_domain_idx,
                   top_scoring_residues=top_residues,
                   modifying_domains=list(bgc_module.modifications),
                   aa10_code=bgc_module.aa10_code,
                   aa34_code=bgc_module.aa34_code)


class AlignmentStep_NRP_Monomer_Info(NamedTuple):
    residue: MonomerResidue
    chirality: Chirality
    modifications: List[NRP_Monomer_Modification]
    rban_name: str
    rban_idx: int

    @classmethod
    def from_nrp_monomer(cls, nrp_monomer: NRP_Monomer) -> AlignmentStep_NRP_Monomer_Info:
        return cls(residue=nrp_monomer.residue,
                   chirality=nrp_monomer.chirality,
                   modifications=list(nrp_monomer.modifications),
                   rban_name=nrp_monomer.rban_name,
                   rban_idx=nrp_monomer.rban_idx)


class AlignmentStepType(Enum):
    MATCH = auto()
    NRP_MONOMER_INSERT = auto()
    BGC_MODULE_SKIP = auto()
    GENE_SKIP = auto()
    BGC_FRAGMENT_SKIP = auto()
    ITERATE_MODULE = auto()
    ITERATE_GENE = auto()


class MatchDetailedScore(NamedTuple):
    residue_score: LogProb
    methylation_score: LogProb
    chirality_score: LogProb


@dataclass
class AlignmentStep:
    bgc_module_info: Optional[AlignmentStep_BGC_Module_Info]
    nrp_monomer_info: Optional[AlignmentStep_NRP_Monomer_Info]
    score: LogProb
    match_detailed_score: Optional[MatchDetailedScore]  # only for MATCH, sum(match_detailed_score) == score
    step_type: AlignmentStepType

    NA = '---'

    def get_score(self) -> LogProb:  # TODO: this method is kept for backward compatibility, remove it in the future
        return self.score


    def to_dict(self) -> Dict[str, str]:
        NA = AlignmentStep.NA
        return OrderedDict({'Gene': self.bgc_module_info.gene_id if self.bgc_module_info else NA,
                            'A-domain_idx': self.bgc_module_info.a_domain_idx if self.bgc_module_info else NA,
                            'Top_scoring_residues': ','.join(self.bgc_module_info.top_scoring_residues)
                            if self.bgc_module_info else NA,
                            'Modifying_domains': ','.join(mod.name for mod in self.bgc_module_info.modifying_domains)
                            if self.bgc_module_info and self.bgc_module_info.modifying_domains else NA,
                            'NRP_residue': self.nrp_monomer_info.residue if self.nrp_monomer_info else NA,
                            'NRP_chirality': self.nrp_monomer_info.chirality.name if self.nrp_monomer_info else NA,
                            'NRP_modifications': ','.join(mod.name for mod in self.nrp_monomer_info.modifications)
                            if self.nrp_monomer_info and self.nrp_monomer_info.modifications else NA,
                            'rBAN_name': self.nrp_monomer_info.rban_name if self.nrp_monomer_info else NA,
                            'rBAN_idx': self.nrp_monomer_info.rban_idx if self.nrp_monomer_info else NA,
                            'Alignment_step': self.step_type.name,
                            'Score': round(self.score, 3),
                            'ResidueScore': round(self.match_detailed_score.residue_score, 3)
                            if self.match_detailed_score else NA,
                            'MethylationScore': round(self.match_detailed_score.methylation_score, 3)
                            if self.match_detailed_score else NA,
                            'ChiralityScore': round(self.match_detailed_score.chirality_score, 3)
                            if self.match_detailed_score else NA,
                            'aa10_code': self.bgc_module_info.aa10_code if self.bgc_module_info else NA,
                            'aa34_code': self.bgc_module_info.aa34_code if self.bgc_module_info else NA})

    @classmethod
    def from_dict(cls, data: dict) -> AlignmentStep:
        # get bgc module info
        if data['Gene'] == AlignmentStep.NA:
            bgc_module_info = None
        else:
            modifying_domains = [BGC_Module_Modification[mod]
                                 for mod in data['Modifying_domains'].split(',')] \
                if data['Modifying_domains'] != AlignmentStep.NA else []
            bgc_module_info = AlignmentStep_BGC_Module_Info(gene_id=data['Gene'],
                                                            a_domain_idx=data['A-domain_idx'],
                                                            top_scoring_residues=data['Top_scoring_residues'].split(','),
                                                            modifying_domains=modifying_domains,
                                                            aa10_code=data['aa10_code'],
                                                            aa34_code=data['aa34_code'])

        # get nrp monomer info
        if data['NRP_residue'] == AlignmentStep.NA:
            nrp_monomer_info = None
        else:
            modifications = [NRP_Monomer_Modification[mod]
                             for mod in data['NRP_modifications'].split(',')] \
                if data['NRP_modifications'] != AlignmentStep.NA else []
            nrp_monomer_info = AlignmentStep_NRP_Monomer_Info(residue=data['NRP_residue'],
                                                              chirality=Chirality[data['NRP_chirality']],
                                                              modifications=modifications,
                                                              rban_name=data['rBAN_name'],
                                                              rban_idx=data['rBAN_idx'])

        # get match detailed score
        step_type = AlignmentStepType[data['Alignment_step']]
        if step_type == AlignmentStepType.MATCH:
            match_detailed_score = MatchDetailedScore(residue_score=data['ResidueScore'],
                                                      methylation_score=data['MethylationScore'],
                                                      chirality_score=data['ChiralityScore'])
        else:
            match_detailed_score = None

        return cls(bgc_module_info=bgc_module_info,
                   nrp_monomer_info=nrp_monomer_info,
                   score=data['Score'],
                   match_detailed_score=match_detailed_score,
                   step_type=step_type)
