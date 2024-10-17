from typing import (
    Any,
    Dict,
    List,
    NamedTuple,
    Tuple,
    Union)

from src.data_types import (
    BGC_Module,
    BGC_MODULE_DUMMY,
    Chirality,
    GeneId,
    LogProb,
    NRP_Monomer,
    NRP_MONOMER_DUMMY
)
from src.matching.scoring_helper import ScoringHelper, FUNCTION_NAME_TO_STEP_TYPE
from src.matching.alignment_types import Alignment, AlignmentStep, AlignmentStepType
from src.matching.dp_types import (
    DP_State,
    DP_Value,
    DP_Table,
    START_STATE,
    valid_state
)
import numpy as np
from itertools import groupby
from functools import partial


def dp_recalc(dp_table: DP_Table,
              transitions: List[Tuple[DP_State, LogProb, AlignmentStepType]]) -> Union[DP_Value, None]:
    return max((DP_Value(dp_table[state].score + score, state, transition_type)
                for state, score, transition_type in transitions),
               default=None)


def get_genes_intervals(assembly_line: List[BGC_Module]) -> Dict[GeneId, Tuple[int, int]]:  # gene_id -> (start, end)
    return {gene_id: (next(i
                           for i, module in enumerate(assembly_line)
                           if module.gene_id == gene_id),
                      next(len(assembly_line) - 1 - i
                           for i, module in enumerate(reversed(assembly_line))
                           if module.gene_id == gene_id))
            for gene_id in {module.gene_id for module in assembly_line}}


def get_fragments_intervals(assembly_line: List[BGC_Module]) -> Dict[int, Tuple[int, int]]:  # fragment_idx -> (start, end)
    return {fragment_idx: (next(i
                                for i, module in enumerate(assembly_line)
                                if module.fragment_idx == fragment_idx),
                            next(len(assembly_line) - 1 - i
                                for i, module in enumerate(reversed(assembly_line))
                                if module.fragment_idx == fragment_idx))
            for fragment_idx in {module.fragment_idx for module in assembly_line}}


def calculate_dp(assembly_line: List[BGC_Module],
                 nrp_monomers: List[NRP_Monomer],
                 dp_helper: ScoringHelper) -> DP_Table:  # functions computing scores and other parameters
    dp_helper.set_bgc_modules(assembly_line)
    dp_helper.set_nrp_monomers(nrp_monomers)

    genes_intervals = get_genes_intervals(assembly_line)
    fragments_intervals = get_fragments_intervals(assembly_line)
    max_gene_reps = dp_helper.scoring_config.max_gene_reps if any(module.iterative_gene for module in assembly_line) \
        else 0
    max_module_reps = dp_helper.scoring_config.max_module_reps if any(module.iterative_module for module in assembly_line) \
        else 0

    dp_table = np.empty((len(assembly_line)+1, len(nrp_monomers)+1, max_gene_reps+1, max_module_reps+1),
                        dtype=DP_Value)
    dp_table[START_STATE] = DP_Value(0, START_STATE, AlignmentStepType.MATCH)  # step type is not important

    recalc = partial(dp_recalc, dp_table)
    for gene_reps in range(max_gene_reps + 1):
        for module_reps in range(max_module_reps + 1):
            for i, bgc_module in enumerate([BGC_MODULE_DUMMY] + assembly_line):
                for j, nrp_mon in enumerate([NRP_MONOMER_DUMMY] + nrp_monomers):
                    if i == 0 and j == 0:
                        continue

                    # regular match/bgc_module_skip/nrp_mon_skip
                    transitions = [(DP_State(i - 1, j, gene_reps, module_reps), dp_helper.bgc_module_skip, (i - 1,)),
                                   (DP_State(i, j - 1, gene_reps, module_reps), dp_helper.nrp_mon_insert, (j - 1,)),
                                   (DP_State(i - 1, j - 1, gene_reps, module_reps), dp_helper.match, (i - 1, j - 1))]

                    prev_bgc_module = assembly_line[i - 2] if i > 1 else None
                    next_bgc_module = assembly_line[i] if i < len(assembly_line) else None

                    # gene skip
                    if prev_bgc_module is not None and bgc_module.gene_id != prev_bgc_module.gene_id:
                        prev_gene_start, prev_gene_end = genes_intervals[prev_bgc_module.gene_id]
                        transitions.append((DP_State(prev_gene_start, j, gene_reps - 1, module_reps),
                                            dp_helper.gene_skip, (prev_gene_start, prev_gene_end)))

                    # fragment skip
                    if prev_bgc_module is not None and bgc_module.fragment_idx != prev_bgc_module.fragment_idx:
                        prev_fragment_start, prev_fragment_end = fragments_intervals[prev_bgc_module.fragment_idx]
                        transitions.append((DP_State(prev_fragment_start, j, gene_reps, module_reps - 1),
                                            dp_helper.fragment_skip, (prev_fragment_start, prev_fragment_end)))

                    # iterate module
                    if next_bgc_module is not None and next_bgc_module.iterative_module:
                        transitions.append((DP_State(i + 1, j, gene_reps, module_reps - 1),
                                            dp_helper.iterate_module, ()))

                    # iterate gene
                    if next_bgc_module is not None \
                            and next_bgc_module.gene_id != bgc_module.gene_id \
                            and next_bgc_module.iterative_gene:
                        next_gene_end = genes_intervals[next_bgc_module.gene_id][1]
                        transitions.append((DP_State(next_gene_end + 1, j, gene_reps - 1, module_reps),
                                            dp_helper.iterate_gene, ()))
                    dp_table[DP_State(i, j, gene_reps, module_reps)] = recalc((dp_state, transition(*(args + (dp_state,))),
                                                                               FUNCTION_NAME_TO_STEP_TYPE[transition.__name__])
                                                                              for dp_state, transition, args in transitions
                                                                              if valid_state(dp_table, dp_state))

    return dp_table


def retrieve_alignment(dp_table: DP_Table,
                       state: DP_State,
                       assembly_line: List[BGC_Module],
                       nrp_monomers: List[NRP_Monomer],
                       scoring_helper: ScoringHelper) -> Alignment:
    if state == START_STATE:
        return []

    parent = dp_table[state].parent
    bgc_module = assembly_line[parent.module_pos] if parent.module_pos < state.module_pos else None
    nrp_monomer = nrp_monomers[parent.monomer_pos] if parent.monomer_pos < state.monomer_pos else None
    score = scoring_helper.match_detailed_score(bgc_module, nrp_monomer) if dp_table[state].action == AlignmentStepType.MATCH \
            else dp_table[state].score - dp_table[parent].score

    return (retrieve_alignment(dp_table, parent, assembly_line, nrp_monomers, scoring_helper)
            + [AlignmentStep(bgc_module=bgc_module,
                             nrp_monomer=nrp_monomer,
                             score=score,
                             action=dp_table[state].action)])


def get_alignment(assembly_line: List[BGC_Module],
                  nrp_monomers: List[NRP_Monomer],
                  dp_helper: ScoringHelper) -> Alignment:
    dp_table = calculate_dp(assembly_line, nrp_monomers, dp_helper)

    def last_state(gene_reps: int, module_reps: int) -> DP_State:
        return DP_State(len(assembly_line), len(nrp_monomers), gene_reps, module_reps)

    final_state = max([last_state(gene_reps, module_reps)
                       for gene_reps in range(dp_table.shape[2])
                       for module_reps in range(dp_table.shape[3])
                       if dp_table[last_state(gene_reps, module_reps)] is not None],
                      key=lambda dp_state: dp_table[dp_state])
    return retrieve_alignment(dp_table, final_state,
                              assembly_line, nrp_monomers,
                              dp_helper)
