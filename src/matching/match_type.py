from __future__ import annotations
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import List, NamedTuple, Optional
from src.antismash_parsing.bgc_variant_types import (
    BGC_ID,
    BGC_Variant_ID,
    BGC_Variant,
    GeneId,
)
from src.rban_parsing.nrp_variant_types import NRP_Variant, NRP_Variant_ID
from src.general_type_aliases import (
    LogProb,
    Prob,
)
from src.matching.alignment_step_type import AlignmentStep
from src.matching.alignment_type import Alignment, alignment_score, show_alignment, alignment_from_str
from more_itertools import split_at
from collections import defaultdict
from itertools import takewhile


@dataclass
class Match:
    genome_id: str
    bgc_variant_id: BGC_Variant_ID
    nrp_variant_id: NRP_Variant_ID
    raw_score: LogProb  # log probability of the match in HMM
    log_odds_vs_avg_nrp: LogProb  # the score of matching the BGC to an "average" NRP of the same length
    log_odds_vs_avg_bgc: LogProb  # the score of matching the NRP to an "average" BGC of the same length
    p_value: Prob
    alignments: List[Alignment]  # alignments of each fragment

    def to_dict(self) -> dict:
        return {'genome_id': self.genome_id,
                'bgc_variant_id': self.bgc_variant_id.to_dict(),
                'nrp_variant_id': self.nrp_variant_id._asdict(),
                'raw_score': self.raw_score,
                'log_odds_vs_avg_nrp': self.log_odds_vs_avg_nrp,
                'log_odds_vs_avg_bgc': self.log_odds_vs_avg_bgc,
                'p_value': self.p_value,
                'alignments': [[dict(alignment_step.to_dict())  # for some reason yaml.dump treats OrderedDict as list of pairs
                                for alignment_step in alignment]
                               for alignment in self.alignments]}

    @classmethod
    def from_dict(cls, data: dict) -> Match:
        return cls(genome_id=data['genome_id'],
                   bgc_variant_id=BGC_Variant_ID.from_dict(data['bgc_variant_id']),
                   nrp_variant_id=NRP_Variant_ID(**data['nrp_variant_id']),
                   raw_score=data['raw_score'],
                   log_odds_vs_avg_nrp=data['log_odds_vs_avg_nrp'],
                   log_odds_vs_avg_bgc=data['log_odds_vs_avg_bgc'],
                   p_value=data['p_value'],
                   alignments=[[AlignmentStep.from_dict(alignment_step_data)
                                for alignment_step_data in alignment_data]
                               for alignment_data in data['alignments']])

    def __str__(self):
        out = StringIO()
        out.write('\n'.join([f'Genome: {self.genome_id}',
                             f'antiSMASH_file: {self.bgc_variant_id.bgc_id.antiSMASH_file}',
                             f'Contig_idx: {self.bgc_variant_id.bgc_id.contig_idx}',
                             f'BGC_idx: {self.bgc_variant_id.bgc_id.bgc_idx}',
                             f'BGC_variant: {self.bgc_variant_id.variant_idx}',
                             f'NRP: {self.nrp_variant_id.nrp_id}',
                             f'NRP_variant: {self.nrp_variant_id.variant_idx}',
                             f'RawScore: {self.raw_score}',
                             f'LogOdds_vs_avg_BGC: {self.log_odds_vs_avg_bgc}',
                             f'LogOdds_vs_avg_NRP: {self.log_odds_vs_avg_nrp}',
                             f'P-value: {self.p_value}',
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
            'antiSMASH_file': Path,
            'Contig_idx': int,
            'BGC_idx': int,
            'BGC_variant_idx': int,
            'NRP': str,
            'NRP_variant_idx': int,
            'RawScore': float,
            'LogOdds_vs_avg_BGC': float,
            'LogOdds_vs_avg_NRP': float,
            'P-value': float,
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

        bgc_variant_id = BGC_Variant_ID(bgc_id=BGC_ID(antiSMASH_file=data['antiSMASH_file'],
                                                      contig_idx=data['Contig_idx'],
                                                      bgc_idx=data['BGC_idx']),
                                        variant_idx=data['BGC_variant_idx'])
        nrp_variant_id = NRP_Variant_ID(nrp_id=data['NRP'],
                                        variant_idx=data['NRP_variant_idx'])

        fst_alignment_line = next(lines_iter)
        if fst_alignment_line.startswith('Fragment'):
            fragments_blocks = split_at(lines_iter, lambda x: x.startswith('Fragment'))
        else:
            fragments_blocks = [[fst_alignment_line] + list(lines_iter)]
        alignments = [alignment_from_str('\n'.join(fragment_block))
                      for fragment_block in fragments_blocks]

        return cls(genome_id=data['Genome'],
                   bgc_variant_id=bgc_variant_id,
                   nrp_variant_id=nrp_variant_id,
                   alignments=alignments,
                   raw_score=data['RawScore'],
                   log_odds_vs_avg_bgc=data['LogOdds_vs_avg_bgc'],
                   log_odds_vs_avg_nrp=data['LogOdds_vs_avg_nrp'],
                   p_value=data['P-value'])