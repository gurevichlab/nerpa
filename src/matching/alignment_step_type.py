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
)
from src.antismash_parsing.antismash_parser_types import GeneId
from src.monomer_names_helper import (
    MonomerResidue,
    Chirality
)
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType
from src.generic.functional import make_optional
from enum import Enum, auto
from collections import OrderedDict
from dataclasses import dataclass


# whole BGC_Module is too big
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


class MatchDetailedScore(NamedTuple):
    residue_score: LogProb
    methylation_score: LogProb
    chirality_score: LogProb


@dataclass
class AlignmentStep:
    bgc_module: Optional[AlignmentStep_BGC_Module_Info]
    nrp_monomer: Optional[rBAN_Monomer]
    score: LogProb
    match_detailed_score: Optional[MatchDetailedScore]  # only for MATCH, sum(match_detailed_score) == score
    step_type: Optional[DetailedHMMEdgeType]  # optional for backward compatibility -- when reading old alignments with different step types

    NA = '---'

    def get_score(self) -> LogProb:  # TODO: this method is kept for backward compatibility, remove it in the future
        return self.score


    def to_dict(self) -> Dict[str, str]:
        NA = AlignmentStep.NA
        return OrderedDict({'Gene': self.bgc_module.gene_id if self.bgc_module else NA,
                            'A-domain_idx': self.bgc_module.a_domain_idx if self.bgc_module else NA,
                            'Top_scoring_residues': ','.join(self.bgc_module.top_scoring_residues)
                            if self.bgc_module else NA,
                            'Modifying_domains': ','.join(mod.name for mod in self.bgc_module.modifying_domains)
                            if self.bgc_module and self.bgc_module.modifying_domains else NA,
                            'NRP_residue': self.nrp_monomer.residue if self.nrp_monomer else NA,
                            'NRP_chirality': self.nrp_monomer.chirality.name if self.nrp_monomer else NA,
                            'NRP_methylated': str(self.nrp_monomer.methylated)
                            if self.nrp_monomer else NA,
                            'rBAN_name': self.nrp_monomer.rban_name if self.nrp_monomer else NA,
                            'rBAN_idx': self.nrp_monomer.rban_idx if self.nrp_monomer else NA,
                            'Alignment_step': self.step_type.name if self.step_type else NA,
                            'Score': round(self.score, 3),
                            'ResidueScore': round(self.match_detailed_score.residue_score, 3)
                            if self.match_detailed_score else NA,
                            'MethylationScore': round(self.match_detailed_score.methylation_score, 3)
                            if self.match_detailed_score else NA,
                            'ChiralityScore': round(self.match_detailed_score.chirality_score, 3)
                            if self.match_detailed_score else NA,
                            'aa10_code': self.bgc_module.aa10_code if self.bgc_module else NA,
                            'aa34_code': self.bgc_module.aa34_code if self.bgc_module else NA})

    @classmethod
    def from_dict(cls, data: dict) -> AlignmentStep:
        # get bgc module info
        if data['Gene'] == AlignmentStep.NA:
            bgc_module = None
        else:
            modifying_domains = [BGC_Module_Modification[mod]
                                 for mod in data['Modifying_domains'].split(',')] \
                if data['Modifying_domains'] != AlignmentStep.NA else []
            bgc_module = AlignmentStep_BGC_Module_Info(gene_id=data['Gene'],
                                                       a_domain_idx=int(data['A-domain_idx']),
                                                       top_scoring_residues=data['Top_scoring_residues'].split(','),
                                                       modifying_domains=modifying_domains,
                                                       aa10_code=data.get('aa10_code', AlignmentStep.NA),
                                                       aa34_code=data.get('aa34_code', AlignmentStep.NA))

        # get nrp monomer info
        if data['NRP_residue'] == AlignmentStep.NA:
            nrp_monomer = None
        else:
            if 'NRP_modifications' in data:  # backward compatibility
                methylated = 'METHYLATION' in data['NRP_modifications'].split(',')
            else:
                methylated = data['NRP_methylated'] == 'True'
            if 'is_pks_hybrid' in data:
                is_pks_hybrid = data['is_pks_hybrid']
            else:
                is_pks_hybrid = False
            nrp_monomer = rBAN_Monomer(residue=data['NRP_residue'],
                                       chirality=Chirality[data['NRP_chirality']],
                                       methylated=methylated,
                                       rban_name=data['rBAN_name'],
                                       rban_idx=data['rBAN_idx'],
                                       is_pks_hybrid=is_pks_hybrid)

        # get match detailed score
        if data['Alignment_step'] == 'NRP_MONOMER_SKIP':
            data['Alignment_step'] = 'NRP_MONOMER_INSERT'  # for backward compatibility, to be removed in the future
        step_type = DetailedHMMEdgeType[data['Alignment_step']] \
            if data['Alignment_step'] in DetailedHMMEdgeType.__members__ else None  # for backward compatibility
        if step_type == DetailedHMMEdgeType.MATCH:
            try:
                match_detailed_score = MatchDetailedScore(residue_score=float(data['ResidueScore']),
                                                          methylation_score=float(data['MethylationScore']),
                                                          chirality_score=float(data['ChiralityScore']))
            except:
                match_detailed_score = None  # for backward compatibility
        else:
            match_detailed_score = None

        return cls(bgc_module=bgc_module,
                   nrp_monomer=nrp_monomer,
                   score=float(data['Score']),
                   match_detailed_score=match_detailed_score,
                   step_type=step_type)
