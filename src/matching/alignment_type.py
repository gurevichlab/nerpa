from typing import List
from prettytable import PrettyTable
import more_itertools
from src.general_type_aliases import LogProb
from src.antismash_parsing.bgc_variant_types import (
    BGC_ID,
    BGC_Variant_ID,
    BGC_Variant,
    BGC_Module,
    GeneId,
)
from src.rban_parsing.nrp_variant_types import NRP_Variant, NRP_Variant_ID
from src.matching.alignment_step_type import AlignmentStep
from src.generic.other import parse_pretty_table

Alignment = List[AlignmentStep]

def alignment_score(alignment: Alignment) -> LogProb:
    return sum(alignment_step.get_score() for alignment_step in alignment)


def combined_alignments_score(alignments: List[Alignment]) -> LogProb:
    return sum(alignment_score(alignment) for alignment in alignments)


def show_alignment(alignment: Alignment,
                   skip_auxiliary_steps: bool = True) -> str:
    #skip_auxiliary_steps = False
    rows = [alignment_step.to_dict()
            for alignment_step in alignment
            if any([alignment_step.bgc_module is not None,
                    alignment_step.nrp_monomer is not None,
                    not skip_auxiliary_steps])]

    t = PrettyTable(rows[0].keys(), align='l', border=False)
    t.add_rows([row.values() for row in rows])
    return str(t)


def alignment_from_str(alignment_str: str) -> Alignment:
    table = parse_pretty_table(alignment_str)
    return [AlignmentStep.from_dict(row)
            for row in table]
