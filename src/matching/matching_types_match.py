from __future__ import annotations
from dataclasses import dataclass
from io import StringIO
from typing import List, NamedTuple
from src.data_types import (
    BGC_Variant,
    LogProb,
    NRP_Variant,
    GeneId
)
from src.matching.matching_types_alignment_step import AlignmentStep
from src.matching.matching_types_alignment import Alignment, alignment_score, show_alignment, alignment_from_str
from more_itertools import split_at
from collections import defaultdict
from itertools import takewhile


class Match_BGC_Variant_Info(NamedTuple):
    genome_id: str
    contig_idx: int
    bgc_idx: int
    variant_idx: int

    @classmethod
    def from_bgc_variant(cls, bgc_variant: BGC_Variant) -> Match_BGC_Variant_Info:
        return cls(genome_id=bgc_variant.genome_id,
                   contig_idx=bgc_variant.contig_idx,
                   bgc_idx=bgc_variant.bgc_idx,
                   variant_idx=bgc_variant.variant_idx)

    def get_antismash_id(self) -> str:
        return f'r{self.contig_idx}c{self.bgc_idx}'


class Match_NRP_Variant_Info(NamedTuple):
    nrp_id: str
    variant_idx: int

    @classmethod
    def from_nrp_variant(cls, nrp_variant: NRP_Variant) -> Match_NRP_Variant_Info:
        return cls(nrp_id=nrp_variant.nrp_id,
                   variant_idx=nrp_variant.variant_idx)

@dataclass
class Match:
    bgc_variant_info: Match_BGC_Variant_Info
    nrp_variant_info: Match_NRP_Variant_Info
    normalized_score: float
    alignments: List[Alignment]  # alignments of each fragment

    def raw_score(self) -> LogProb:
        return sum(map(alignment_score, self.alignments))

    def to_dict(self) -> dict:  # because full Match is too big to write
        return {'Genome': self.bgc_variant_info.genome_id,
                'Contig_idx': self.bgc_variant_info.contig_idx,
                'BGC_idx': self.bgc_variant_info.bgc_idx,
                'BGC_variant_idx': self.bgc_variant_info.variant_idx,
                'NRP': self.nrp_variant_info.nrp_id,
                'NRP_variant_idx': self.nrp_variant_info.variant_idx,
                'NormalisedScore': self.normalized_score,
                'Score': self.raw_score(),
                'Alignments': [[dict(alignment_step.to_dict())  # for some reason yaml.dump treats OrderedDict as list of pairs
                                for alignment_step in alignment]
                               for alignment in self.alignments]}

    @classmethod
    def from_dict(cls, data: dict) -> Match:
        return cls(bgc_variant_info=Match_BGC_Variant_Info(genome_id=GeneId(data['Genome']),
                                                           contig_idx=data.get('Contig_idx', 0),
                                                           bgc_idx=data.get('BGC_idx', 0),
                                                           variant_idx=data['BGC_variant_idx']),
                   nrp_variant_info=Match_NRP_Variant_Info(nrp_id=data['NRP'],
                                                           variant_idx=data['NRP_variant_idx']),
                   normalized_score=data['NormalisedScore'],
                   alignments=[[AlignmentStep.from_dict(alignment_step_data)
                                for alignment_step_data in alignment_data]
                               for alignment_data in data['Alignments']])

    def __str__(self):
        out = StringIO()
        out.write('\n'.join([f'Genome: {self.bgc_variant_info.genome_id}',
                             f'Contig_idx: {self.bgc_variant_info.contig_idx}',
                             f'BGC_idx: {self.bgc_variant_info.bgc_idx}',
                             f'BGC_variant: {self.bgc_variant_info.variant_idx}',
                             f'NRP: {self.nrp_variant_info.nrp_id}',
                             f'NRP_variant: {self.nrp_variant_info.variant_idx}',
                             f'NormalisedScore: {self.normalized_score}',
                             f'Score: {self.raw_score()}',
                             'Alignment:']))
        out.write('\n')

        for i, alignment in enumerate(self.alignments):
            if len(self.alignments) > 1:
                out.write(f'Fragment_#{i}\n')
            out.write(show_alignment(alignment) + '\n')

        return out.getvalue()

    @classmethod
    def from_str(cls, match_str: str,
                 print_warnings: bool = False) -> Match:
        lines = match_str.splitlines()
        # q: remove empty lines at the beginning and end
        fst_non_empty_line = next(i for i, line in enumerate(lines) if line.strip())
        last_non_empty_line = next(i for i, line in enumerate(reversed(lines)) if line.strip())
        lines_iter = iter(lines[fst_non_empty_line:len(lines) - last_non_empty_line])

        # parsing is a bit weird because of backward compatibility: some fields are optional
        data = defaultdict(lambda: None)
        field_type = {
            'Genome': str,
            'Contig_idx': int,
            'BGC_idx': int,
            'BGC_variant_idx': int,
            'NRP': str,
            'NRP_variant_idx': int,
            'NormalisedScore': float,
            'Score': float
        }
        field_lines = takewhile(lambda x: not x.startswith('Alignment:'), lines_iter)
        for line in field_lines:
            field_name, value = line.split(': ')
            if field_name in field_type:
                try:
                    data[field_name] = field_type[field_name](value)
                except:
                    if print_warnings:
                        print('Error parsing field:', field_name, value)
                        print('Assigning None')
                    data[field_name] = None

        bgc_variant_info = Match_BGC_Variant_Info(genome_id=data['Genome'],
                                                  contig_idx=data['Contig_idx'],
                                                  bgc_idx=data['BGC_idx'],
                                                  variant_idx=data['BGC_variant_idx'])
        nrp_variant_info = Match_NRP_Variant_Info(nrp_id=data['NRP'],
                                                  variant_idx=data['NRP_variant_idx'])

        fst_alignment_line = next(lines_iter)
        if fst_alignment_line.startswith('Fragment'):
            fragments_blocks = split_at(lines_iter, lambda x: x.startswith('Fragment'))
        else:
            fragments_blocks = [[fst_alignment_line] + list(lines_iter)]
        alignments = [alignment_from_str('\n'.join(fragment_block))
                      for fragment_block in fragments_blocks]
        return cls(bgc_variant_info=bgc_variant_info,
                   nrp_variant_info=nrp_variant_info,
                   alignments=alignments,
                   normalized_score=data['NormalisedScore'])