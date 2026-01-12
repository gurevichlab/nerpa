from __future__ import annotations

from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    StateIdx,
)

ET = DetailedHMMEdgeType
ST = DetailedHMMStateType


def add_ordinary_module_subgraph(constructor,
                                  module_idx: StateIdx):  # constructor: HMM_Constructor
    module = constructor.bgc_variant.modules[module_idx]

    constructor.states.extend([
        DetailedHMMState(state_type=ST.MATCH,
                         emissions=constructor.get_emissions(module),
                         related_module_idx=module_idx),
        DetailedHMMState(state_type=ST.INSERT,
                         emissions=constructor.get_insert_emissions(module),
                         related_module_idx=module_idx)
    ])
    match_state_idx = len(constructor.states) - 2
    insert_state_idx = len(constructor.states) - 1
    root = constructor.module_idx_to_subgraph_root[module_idx]

    constructor.add_edge(root, match_state_idx, ET.MATCH)
    constructor.add_edge(match_state_idx, insert_state_idx, ET.START_INSERTING)
    constructor.add_edge(insert_state_idx, insert_state_idx, ET.INSERT)

    next_module_root = constructor.module_idx_to_subgraph_root[module_idx + 1]
    iteration_subgraph_root = constructor.add_iteration_subgraph(module_idx=module_idx,
                                                                 node_after_subgraph=next_module_root)

    constructor.add_edge(root, iteration_subgraph_root, ET.SKIP_MODULE)
    constructor.add_edge(match_state_idx, iteration_subgraph_root, ET.NO_INSERTIONS)
    constructor.add_edge(insert_state_idx, iteration_subgraph_root, ET.FINISH_INSERTING)


def add_tentative_assembly_end_module_subgraph(constructor,
                                                module_idx: StateIdx):  # constructor: HMM_Constructor
    module = constructor.bgc_variant.modules[module_idx]
    if constructor.construction_cfg.ASSEMBLY_LINE_BORDERS_ONLY_AT_MARKERS:
        assert ModuleGenomicContextFeature.ASSEMBLY_LINE_FINISH_MARKER in module.genomic_context

    constructor.states.extend([
        DetailedHMMState(state_type=ST.MATCH_POSSIBLE_ASSEMBLY_FINISH,
                         emissions=constructor.get_emissions(module),
                         related_module_idx=module_idx),
        DetailedHMMState(state_type=ST.CHOOSE_IF_FINISH_MATCHING,
                         related_module_idx=module_idx),
        DetailedHMMState(state_type=ST.MATCHING_FINISHED,
                         related_module_idx=module_idx),
        DetailedHMMState(state_type=ST.MATCHING_CONTINUED,
                         related_module_idx=module_idx),
        DetailedHMMState(state_type=ST.INSERT,
                         emissions=constructor.get_insert_emissions(module),
                         related_module_idx=module_idx)
    ])
    match_state_idx = len(constructor.states) - 5
    choose_if_finish_matching = len(constructor.states) - 4
    matching_finished = len(constructor.states) - 3
    matching_continued = len(constructor.states) - 2
    insert_state_idx = len(constructor.states) - 1

    root = constructor.module_idx_to_subgraph_root[module_idx]

    constructor.add_edge(root, match_state_idx, ET.MATCH)
    constructor.add_edge(root, choose_if_finish_matching, ET.SKIP_MODULE)
    constructor.add_edge(match_state_idx, choose_if_finish_matching, ET.CHOOSE_IF_FINISH_MATCHING)

    constructor.add_edge(choose_if_finish_matching, matching_finished, ET.FINISH_MATCHING)
    constructor.add_edge(choose_if_finish_matching, matching_continued, ET.CONTINUE_MATCHING)
    constructor.add_edge(matching_continued, insert_state_idx, ET.START_INSERTING)
    constructor.add_edge(insert_state_idx, insert_state_idx, ET.INSERT)

    next_module_root = constructor.module_idx_to_subgraph_root[module_idx + 1]
    iteration_subgraph_root = constructor.add_iteration_subgraph(module_idx=module_idx,
                                                                 node_after_subgraph=next_module_root)

    constructor.add_edge(matching_continued, iteration_subgraph_root, ET.NO_INSERTIONS)
    constructor.add_edge(insert_state_idx, iteration_subgraph_root, ET.FINISH_INSERTING)


def add_last_module_subgraph(constructor):  # constructor: HMM_Constructor
    last_module_idx = len(constructor.bgc_variant.modules) - 1
    last_module = constructor.bgc_variant.modules[-1]
    constructor.states.extend([
        DetailedHMMState(state_type=ST.MATCH_LAST_MODULE,
                         emissions=constructor.get_emissions(last_module),
                         related_module_idx=last_module_idx),
        DetailedHMMState(state_type=ST.MATCHING_FINISHED,
                         related_module_idx=last_module_idx),
        DetailedHMMState(state_type=ST.INSERT_AT_END,
                         emissions=constructor.get_insert_emissions(last_module),
                         related_module_idx=last_module_idx)
    ])
    match_state_idx = len(constructor.states) - 3
    matching_finished = len(constructor.states) - 2
    insert_state_idx = len(constructor.states) - 1

    root = constructor.module_idx_to_subgraph_root[last_module_idx]

    constructor.add_edge(root, match_state_idx, ET.MATCH)
    constructor.add_edge(root, matching_finished, ET.SKIP_MODULE)

    constructor.add_edge(match_state_idx, matching_finished, ET.FINISH_MATCHING)

    constructor.add_edge(matching_finished, insert_state_idx, ET.START_INSERTING)
    constructor.add_edge(insert_state_idx, insert_state_idx, ET.INSERT_AT_END)

    iteration_subgraph_root = constructor.add_iteration_subgraph(module_idx=last_module_idx,
                                                                 node_after_subgraph=constructor.final_state)

    constructor.add_edge(matching_finished, iteration_subgraph_root, ET.NO_INSERTIONS)
    constructor.add_edge(insert_state_idx, iteration_subgraph_root, ET.FINISH_INSERTING)

MCF = ModuleGenomicContextFeature

def add_module_subgraphs(constructor):  # constructor: HMM_Constructor
    for module_idx, module in enumerate(constructor.bgc_variant.modules[:-1]):  # the order does not matter
        if (constructor.construction_cfg.ASSEMBLY_LINE_BORDERS_ONLY_AT_MARKERS
                and MCF.ASSEMBLY_LINE_FINISH_MARKER not in module.genomic_context):
                add_ordinary_module_subgraph(constructor, module_idx)
        else:
            add_tentative_assembly_end_module_subgraph(constructor, module_idx)

    add_last_module_subgraph(constructor)
