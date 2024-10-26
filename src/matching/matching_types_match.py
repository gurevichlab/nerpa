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


class Match_BGC_Variant_Info(NamedTuple):
    genome_id: str
    bgc_idx: int
    variant_idx: int


class Match_NRP_Variant_Info(NamedTuple):
    nrp_id: str
    variant_idx: int

@dataclass
class Match:
    bgc_variant_info: Match_BGC_Variant_Info
    nrp_variant_info: Match_NRP_Variant_Info
    normalized_score: float
    alignments: List[Alignment]  # alignments of each fragment

    def __init__(self,
                 bgc_variant: BGC_Variant,
                 nrp_variant: NRP_Variant,
                 normalized_score: float,
                 alignments: List[Alignment]):
        self.bgc_variant_info = Match_BGC_Variant_Info(genome_id=bgc_variant.genome_id,
                                                       bgc_idx=bgc_variant.bgc_idx,
                                                       variant_idx=bgc_variant.variant_idx)
        self.nrp_variant_info = Match_NRP_Variant_Info(nrp_id=nrp_variant.nrp_id,
                                                         variant_idx=nrp_variant.variant_idx)
        self.normalized_score = normalized_score
        self.alignments = alignments


    def raw_score(self) -> LogProb:
        return sum(map(alignment_score, self.alignments))

    def to_dict(self) -> dict:  # because full Match is too big to write
        return {'Genome': self.bgc_variant_info.genome_id,
                'BGC': self.bgc_variant_info.bgc_idx,
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
                                                           bgc_idx=data['BGC'],
                                                           variant_idx=data['BGC_variant_idx']),
                   nrp_variant_info=Match_NRP_Variant_Info(nrp_id=data['NRP'],
                                                           variant_idx=data['NRP_variant_idx']),
                   normalized_score=data['Normalised_score'],
                   alignments=[[AlignmentStep.from_dict(alignment_step_data)
                                for alignment_step_data in alignment_data]
                               for alignment_data in data['Alignments']])

    def __str__(self):
        out = StringIO()
        out.write('\n'.join([f'Genome={self.bgc_variant_info.genome_id}',
                             f'BGC={self.bgc_variant_info.bgc_idx}',
                             f'BGC_variant={self.bgc_variant_info.variant_idx}',
                             f'NRP={self.nrp_variant_info.nrp_id}',
                             f'NRP_variant={self.nrp_variant_info.variant_idx}',
                             f'NormalisedScore={self.normalized_score}',
                             f'Score={self.raw_score()}']))
        out.write('\n')

        for i, alignment in enumerate(self.alignments):
            if len(self.alignments) > 1:
                out.write(f'Fragment_#{i}\n')
            out.write(show_alignment(alignment) + '\n')

        return out.getvalue()

    @classmethod
    def from_str(cls, match_str: str) -> Match:
        lines = match_str.splitlines()
        # q: remove empty lines at the beginning and end
        fst_non_empty_line = next(i for i, line in enumerate(lines) if line.strip())
        last_non_empty_line = next(i for i, line in enumerate(reversed(lines)) if line.strip())
        lines = lines[fst_non_empty_line:len(lines) - last_non_empty_line]

        genome_id = lines[0].split('=')[1]
        bgc_idx = int(lines[1].split('=')[1])
        bgc_variant_idx = int(lines[2].split('=')[1])
        nrp_id = lines[3].split('=')[1]
        nrp_variant_idx = int(lines[4].split('=')[1])
        normalized_score = float(lines[5].split('=')[1])
        score = float(lines[6].split('=')[1])

        bgc_variant_info = Match_BGC_Variant_Info(genome_id=genome_id,
                                                  bgc_idx=bgc_idx,
                                                  variant_idx=bgc_variant_idx)
        nrp_variant_info = Match_NRP_Variant_Info(nrp_id=nrp_id,
                                                  variant_idx=nrp_variant_idx)

        if lines[7].startswith('Fragment'):
            fragments_blocks = split_at(lines[7:], lambda x: x.startswith('Fragment'))
        else:
            fragments_blocks = [lines[7:]]
        alignments = [alignment_from_str('\n'.join(fragment_block))
                      for fragment_block in fragments_blocks]
        return cls(bgc_variant_info=bgc_variant_info,
                   nrp_variant_info=nrp_variant_info,
                   alignments=alignments,
                   normalized_score=normalized_score)