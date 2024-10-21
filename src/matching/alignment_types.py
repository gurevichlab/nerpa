from __future__ import annotations
from typing import (
    Dict,
    NamedTuple,
    List,
    Optional,
    Tuple,
    Union
)
from src.data_types import (
    BGC_Module,
    NRP_Monomer,
    BGC_Variant,
    LogProb,
    NRP_Variant,
    NRP_Fragment,
    BGC_Module_Modification,
    Chirality
)
from src.antismash_parsing.antismash_parser_types import GeneId
from src.monomer_names_helper import MonomerResidue, NRP_Monomer_Modification
from dataclasses import dataclass
from prettytable import PrettyTable
from io import StringIO
from enum import Enum, auto
from collections import OrderedDict
import more_itertools


class AlignmentStep_BGC_Module_Info(NamedTuple):
    gene_id: GeneId
    a_domain_idx: int
    top_scoring_residues: List[MonomerResidue]
    modifying_domains: List[BGC_Module_Modification]
    aa10_code: str
    aa34_code: str

    @classmethod
    def from_bgc_module(cls, BGC_Module) -> AlignmentStep_BGC_Module_Info:
        top_score = max(BGC_Module.residue_score.values())
        top_residues = [res for res, score in BGC_Module.residue_score.items()
                        if score == top_score]
        return cls(gene_id=BGC_Module.gene_id,
                   a_domain_idx=BGC_Module.a_domain_idx,
                   top_scoring_residues=top_residues,
                   modifying_domains=BGC_Module.modifications,
                   aa10_code=BGC_Module.aa10_code,
                   aa34_code=BGC_Module.aa34_code)


class AlignmentStep_NRP_Monomer_Info(NamedTuple):
    residue: MonomerResidue
    chirality: Chirality
    modifications: List[NRP_Monomer_Modification]
    rban_name: str
    rban_idx: int

    @classmethod
    def from_nrp_monomer(cls, NRP_Monomer) -> AlignmentStep_NRP_Monomer_Info:
        return cls(residue=NRP_Monomer.residue,
                   chirality=NRP_Monomer.chirality,
                   modifications=NRP_Monomer.modifications,
                   rban_name=NRP_Monomer.rban_name,
                   rban_idx=NRP_Monomer.rban_idx)


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


class AlignmentStep(NamedTuple):
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
            bgc_module_info = AlignmentStep_BGC_Module_Info(gene_id=data['Gene'],
                                                            a_domain_idx=data['A-domain_idx'],
                                                            top_scoring_residues=data['Top_scoring_residues'].split(','),
                                                            modifying_domains=[BGC_Module_Modification[mod]
                                                                               for mod in data['Modifying_domains'].split(',')],
                                                            aa10_code=data['aa10_code'],
                                                            aa34_code=data['aa34_code'])

        # get nrp monomer info
        if data['NRP_residue'] == AlignmentStep.NA:
            nrp_monomer_info = None
        else:
            nrp_monomer_info = AlignmentStep_NRP_Monomer_Info(residue=data['NRP_residue'],
                                                              chirality=Chirality[data['NRP_chirality']],
                                                              modifications=[NRP_Monomer_Modification[mod]
                                                                             for mod in data['NRP_modifications'].split(',')],
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


Alignment = List[AlignmentStep]

def alignment_score(alignment: Alignment) -> LogProb:
    return sum(alignment_step.get_score() for alignment_step in alignment)


def combined_alignments_score(alignments: List[Alignment]) -> LogProb:
    return sum(alignment_score(alignment) for alignment in alignments)


def show_alignment(alignment: Alignment) -> str:
    rows = [alignment_step.to_dict()
            for alignment_step in alignment]

    t = PrettyTable(rows[0].keys(), align='l', border=False)
    t.add_rows(row.values() for row in rows)
    return str(t)


def alignment_from_str(alignment_str: str) -> Alignment:
    lines = alignment_str.strip().split('\n')
    step_lines = more_itertools.split_at(lines, lambda line: line.startswith('Gene'))
    return [AlignmentStep(bgc_module=BGC_Module.from_str(step_lines[0][0].split('=')[1]),
                          nrp_monomer=NRP_Monomer.from_str(step_lines[0][1].split('=')[1]),
                          score=eval(step_lines[0][2].split('=')[1]),
                          action=AlignmentStepType[step_lines[0][3].split('=')[1])]
            for step_lines in step_lines]


def alignment_to_str(alignment_str: str) -> Alignment:

@dataclass
class Match:
    bgc_variant: BGC_Variant
    nrp_variant: NRP_Variant
    alignments: List[Alignment]  # alignments of each fragment
    normalized_score: float

    def raw_score(self) -> LogProb:
        return sum(map(alignment_score, self.alignments))

    def __str__(self):
        out = StringIO()
        out.write('\n'.join([f'Genome={self.bgc_variant.genome_id}',
                             f'BGC={self.bgc_variant.bgc_idx}',
                             f'BGC_variant={self.bgc_variant.variant_idx}',
                             f'NRP={self.nrp_variant.nrp_id}',
                             f'NRP_variant={self.nrp_variant.variant_idx}',
                             f'NormalisedScore={self.normalized_score}',
                             f'Score={self.raw_score()}']))
        out.write('\n')

        for i, alignment in enumerate(self.alignments):
            if len(self.alignments) > 1:
                out.write(f'Fragment_#{i}\n')
            out.write(show_alignment(alignment) + '\n')

        return out.getvalue()

    def to_dict_light(self) -> dict:  # because full Match is too big to write
        return {'Genome': self.bgc_variant.genome_id,
                'BGC': self.bgc_variant.bgc_idx,
                'BGC_variant_idx': self.bgc_variant.variant_idx,
                'NRP': self.nrp_variant.nrp_id,
                'NRP_variant_idx': self.nrp_variant.variant_idx,
                'NormalisedScore': self.normalized_score,
                'Score': self.raw_score(),
                'Alignments': [[dict(alignment_step.to_dict())  # for some reason yaml.dump treats OrderedDict as list of pairs
                                for alignment_step in alignment]
                               for alignment in self.alignments]}

    @classmethod
    def from_yaml_dict(cls, data: dict) -> Match:
        return cls(bgc_variant=BGC_Variant.from_yaml_dict(data['bgc_variant']),
                   nrp_variant=NRP_Variant.from_yaml_dict(data['nrp_variant']),
                   alignments=[[AlignmentStep.from_yaml_dict(alignment_step_data)
                                for alignment_step_data in alignment_data]
                               for alignment_data in data['alignments']],
                   normalized_score=data['normalised_score'])


@dataclass
class MatchLight:
    genome_id: str
    bgc_idx: int
    bgc_variant_idx: int
    nrp_id: str
    nrp_variant_idx: int
    normalized_score: LogProb
    alignments: List[Alignment]

    def raw_score(self) -> LogProb:
        return sum(map(alignment_score, self.alignments))

    def __str__(self):
        out = StringIO()
        out.write('\n'.join([f'Genome={self.genome_id}',
                             f'BGC={self.bgc_idx}',
                             f'BGC_variant={self.bgc_variant_idx}',
                             f'NRP={self.nrp_id}',
                             f'NRP_variant={self.nrp_variant_idx}',
                             f'NormalizedScore={self.normalized_score}',
                             f'Score={self.raw_score()}']))
        out.write('\n')

        for i, alignment in enumerate(self.alignments):
            if len(self.alignments) > 1:
                out.write(f'Fragment_#{i}\n')
            out.write(show_alignment(alignment) + '\n')

        return out.getvalue()

    @classmethod
    def from_str(cls, string: str) -> MatchLight:
        lines = string.strip().split('\n')
        genome_id = lines[0].split('=')[1]
        bgc_idx = int(lines[1].split('=')[1])
        bgc_variant_idx = int(lines[2].split('=')[1])
        nrp_id = lines[3].split('=')[1]
        nrp_variant_idx = int(lines[4].split('=')[1])
        normalized_score = float(lines[5].split('=')[1])

        alignments = []
        if lines[6].startswith('Fragment'):
            fragments_strs = more_itertools.split_at(lines[6:], lambda line: line.startswith('Fragment'))
        else:
            fragments_strs = [lines[6:]]
        alignments = [alignment_from_str(fragment_str)
                        for fragment_str in fragments_strs]
        return cls(genome_id=genome_id,
                   bgc_idx=bgc_idx,
                   bgc_variant_idx=bgc_variant_idx,
                   nrp_id=nrp_id,
                   nrp_variant_idx=nrp_variant_idx,
                   normalized_score=normalized_score,
                   alignments=alignments)
