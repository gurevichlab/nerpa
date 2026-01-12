from typing import List, Optional
from prettytable import PrettyTable
import more_itertools
from src.general_type_aliases import LogProb
from src.antismash_parsing.bgc_variant_types import (
    BGC_ID,
    BGC_Variant_ID,
    BGC_Variant,
    BGC_Module,
    GeneId, A_Domain_ID,
)
from src.rban_parsing.nrp_variant_types import NRP_Variant, NRP_Variant_ID
from src.matching.alignment_step_type import AlignmentStep, AlignmentStepLight, AlignmentStep_BGC_Module_Info
from src.generic.other import parse_pretty_table
from src.rban_parsing.rban_monomer import rBAN_Monomer

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


AlignmentLight = List[AlignmentStepLight]

def alignment_to_light_alignment(alignment: Alignment) -> List[AlignmentStepLight]:
    def bgc_module_info_to_a_domain_id(bgc_module_info: Optional[AlignmentStep_BGC_Module_Info]) -> Optional[A_Domain_ID]:
        if bgc_module_info is None:
            return None
        return A_Domain_ID(gene_id=bgc_module_info.gene_id,
                           a_domain_idx=bgc_module_info.a_domain_idx)

    return [AlignmentStepLight(bgc_module_info_to_a_domain_id(step.bgc_module),
                               step.nrp_monomer)
            for step in alignment]

def light_alignment_to_str(alignment_light: AlignmentLight) -> str:
    def show_a_domain_id(a_domain_id: Optional[A_Domain_ID]) -> str:
        if a_domain_id is None:
            return 'None'
        return f'{a_domain_id.gene_id}:{a_domain_id.a_domain_idx}'
    def show_nrp_monomer(nrp_monomer: Optional[rBAN_Monomer]) -> str:
        if nrp_monomer is None:
            return 'None'
        return nrp_monomer.rban_name

    return '\n'.join(f'{show_a_domain_id(step.a_domain_id)}\t{show_nrp_monomer(step.nrp_monomer)}'
                     for step in alignment_light)