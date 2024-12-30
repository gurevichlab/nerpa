from typing import List
from prettytable import PrettyTable
import more_itertools
from src.data_types import LogProb, BGC_Module
from src.matching.matching_types_alignment_step import AlignmentStep
from src.generic.other import parse_pretty_table

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
    table = parse_pretty_table(alignment_str)
    return [AlignmentStep.from_dict(row)
            for row in table]
