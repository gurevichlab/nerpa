from __future__ import annotations

from itertools import pairwise

from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    StateIdx,
)

ET = DetailedHMMEdgeType
ST = DetailedHMMStateType


def add_skipping_at_start_finish_subgraph(constructor,  # constructor: HMM_Constructor
                                          module_idx: int) -> StateIdx:
    # returns the index of the SKIPPING_MODULES_AT_START_FINISHED state
    module = constructor.bgc_variant.modules[module_idx]
    constructor.states.extend([
        DetailedHMMState(state_type=ST.SKIPPING_MODULES_AT_START_FINISHED,
                         related_module_idx=module_idx),
        DetailedHMMState(state_type=ST.INSERT_AT_START,
                         emissions=constructor.get_insert_emissions(module),
                         related_module_idx=module_idx)
    ])
    skipping_finished = len(constructor.states) - 2
    insert_at_start = len(constructor.states) - 1
    module_subgraph_root = constructor.module_idx_to_subgraph_root[module_idx]

    constructor.add_edge(skipping_finished, insert_at_start, ET.START_INSERTING_AT_START)
    constructor.add_edge(skipping_finished, module_subgraph_root, ET.NO_INSERTIONS)

    constructor.add_edge(insert_at_start, insert_at_start, ET.INSERT_AT_START)
    constructor.add_edge(insert_at_start, module_subgraph_root, ET.FINISH_INSERTING)

    return skipping_finished


def add_skips_at_start(constructor, cls):  # constructor: HMM_Constructor
    if constructor.construction_cfg.ASSEMBLY_LINE_BORDERS_ONLY_AT_MARKERS:
        tentative_assembly_starts = [
           module_idx
           for module_idx, module in enumerate(constructor.bgc_variant.modules)
           if (ModuleGenomicContextFeature.ASSEMBLY_LINE_START_MARKER in module.genomic_context
               or module_idx == 0)
        ]
    else:
        tentative_assembly_starts = list(range(len(constructor.bgc_variant.modules)))

    # only one tentative assembly line start -> skipping at start impossible
    if len(tentative_assembly_starts) == 1:
        skipping_finished = add_skipping_at_start_finish_subgraph(constructor, 0)
        constructor.add_edge(constructor.initial_state, skipping_finished, ET.INITIAL)
        return

    # multiple tentative assembly line starts -> skipping at start possible
    constructor.states.append(
        DetailedHMMState(state_type=ST.CHOOSE_IF_START_MATCHING,
                         related_module_idx=0)
    )
    choose_if_start_matching = len(constructor.states) - 1
    skipping_finished = add_skipping_at_start_finish_subgraph(constructor, 0)

    constructor.add_edge(constructor.initial_state, choose_if_start_matching, ET.INITIAL)
    constructor.add_edge(choose_if_start_matching, skipping_finished, ET.START_MATCHING)

    # invariant: choose_if_start_matching points to the CHOOSE_IF_START_MATCHING state
    # which decides whether to start matching at tent_start or skip to next_tent_start
    # the "SKIPPING_MODULES_AT_START_FINISHED" path for choose_if_start_matching is already added
    # but the "SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START" path is not yet added

    for tent_start, next_tent_start in pairwise(tentative_assembly_starts):  # all but last
        skipping_finished = add_skipping_at_start_finish_subgraph(constructor, next_tent_start)
        constructor.states.append(
            DetailedHMMState(state_type=ST.CHOOSE_IF_START_MATCHING,
                             related_module_idx=next_tent_start)
        )
        new_choose_if_start_matching = len(constructor.states) - 1
        constructor.add_edge(choose_if_start_matching,
                             new_choose_if_start_matching,
                             ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START)
        constructor.add_edge(new_choose_if_start_matching,
                             skipping_finished,
                             ET.START_MATCHING)
        choose_if_start_matching = new_choose_if_start_matching

    #constructor.convert_to_detailed_hmm(cls).draw(Path("adding_skips_at_start.png"))
    constructor.add_edge(choose_if_start_matching,
                         constructor.final_state,
                         ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START)




