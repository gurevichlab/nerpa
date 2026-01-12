from __future__ import annotations

import math
from typing import Tuple, Literal

from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.general_type_aliases import LogProb
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    DetailedHMMEdge,
)

from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import RELEVANT_GENOMIC_CONTEXT

GenomicContext = Tuple[ModuleGenomicContextFeature, ...]
GC = GenomicContext
ET = DetailedHMMEdgeType
ST = DetailedHMMStateType
MCF = ModuleGenomicContextFeature

def make_edge(helper,  # type: HMM_Constructor,
              state: DetailedHMMState,
              edge_type: DetailedHMMEdgeType,
              ) -> DetailedHMMEdge:
    # 0. Aliases for brevity
    bgc_variant = helper.bgc_variant
    hmm_helper = helper.hmm_helper

    WEIGHT_PARAMS = hmm_helper.scoring_config.edge_weight_parameters  # for brevity
    # stub for debugging
    # WEIGHT_PARAMS: Dict[ST, Dict[GC, Dict[ET, LogProb]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: -1)))

    current_module = (
        bgc_variant.modules[state.related_module_idx]
        if state.related_module_idx is not None
        else None
    )

    full_genomic_context = current_module.genomic_context if current_module is not None else ()
    relevant_genomic_context = tuple(f
                                     for f in full_genomic_context
                                     if f in RELEVANT_GENOMIC_CONTEXT[state.state_type])

    # 1. Sanity checks
    if edge_type == ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START:
        start_marker_present = (MCF.ASSEMBLY_LINE_START_MARKER in full_genomic_context
                                or state.related_module_idx == 0)
        if (not start_marker_present
                and helper.construction_cfg.ASSEMBLY_LINE_BORDERS_ONLY_AT_MARKERS):
            raise ValueError('Skipping until the end of assembly line '
             'should start after a tentative assembly line finish')

    if edge_type == ET.SKIP_MODULES_AT_END:
        finish_marker_present = MCF.ASSEMBLY_LINE_FINISH_MARKER in full_genomic_context
        if (not finish_marker_present
                and helper.construction_cfg.ASSEMBLY_LINE_BORDERS_ONLY_AT_MARKERS):
            raise ValueError('Skipping until the end of assembly line '
                             'should start after a tentative assembly line finish')

    # 2. Helper function
    def skip_at_start_cost(module_idx: int) -> LogProb:
        module_ctxt = tuple(f
                            for f in bgc_variant.modules[module_idx].genomic_context
                            if f in RELEVANT_GENOMIC_CONTEXT[ST.CHOOSE_IF_START_MATCHING])
        if (MCF.START_OF_BGC not in module_ctxt
                and MCF.ASSEMBLY_LINE_START_MARKER not in module_ctxt):
            module_ctxt = tuple(sorted(module_ctxt + (MCF.ASSEMBLY_LINE_START_MARKER,)))
        return WEIGHT_PARAMS[ST.CHOOSE_IF_START_MATCHING][module_ctxt][ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START]

    def skip_at_end_cost(module_idx: int) -> LogProb:
        module_ctxt = tuple(f
                            for f in bgc_variant.modules[module_idx].genomic_context
                            if f in RELEVANT_GENOMIC_CONTEXT[ST.SKIPPING_MODULES_AT_END])
        return WEIGHT_PARAMS[ST.SKIPPING_MODULES_AT_END][module_ctxt][ET.SKIP_MODULES_AT_END]


    def get_full_skip_cost(fst_skipped_module_idx: int,
                           last_skipped_module_idx: int,
                           skip_state_type: Literal[ST.CHOOSE_IF_START_MATCHING, ST.CHOOSE_IF_FINISH_MATCHING]) -> LogProb:
        # That's a little hack: all skips of individual modules are compressed into
        # a single skip edge, but I pretend that there's a chain of skip_module_type states (they don't exist in the HMM)
        # to compute the total skip cost while maintaining WEIGHT_PARAMS structure.
        full_skip_cost = 0
        for skipped_module_idx in range(fst_skipped_module_idx, last_skipped_module_idx + 1):
            full_skip_cost += (
                skip_at_start_cost(skipped_module_idx)
                if skip_state_type == ST.CHOOSE_IF_START_MATCHING
                else skip_at_end_cost(skipped_module_idx)
            )

        return full_skip_cost

    # 3. Compute weight

    # 3.1 special case 1 -- if we choose to skip until next tentative assembly line start
    # we need to aggregate the penalty over all the modules skipped
    if (state.state_type == ST.CHOOSE_IF_START_MATCHING
            and MCF.ASSEMBLY_LINE_START_MARKER in full_genomic_context):
        fst_skipped_module_idx = state.related_module_idx
        lst_skipped_module_idx = next((module_idx - 1
                                       for module_idx, module in enumerate(bgc_variant.modules[state.related_module_idx + 1:],
                                                                           start=state.related_module_idx + 1)
                                       if MCF.ASSEMBLY_LINE_START_MARKER in module.genomic_context),
                                      len(bgc_variant.modules) - 1)

        skip_cost = get_full_skip_cost(fst_skipped_module_idx,
                                       lst_skipped_module_idx,
                                       ST.CHOOSE_IF_START_MATCHING)
        weight = (
            skip_cost
            if edge_type == ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START
            else math.log(1 - math.exp(skip_cost))
        )
    elif (state.state_type == ST.CHOOSE_IF_FINISH_MATCHING):
        # 3.2 special case 2 -- if we choose to skip until the end
        #  we need to aggregate the penalty over all the modules skipped
        fst_skipped_module_idx = state.related_module_idx + 1
        lst_skipped_module_idx = len(bgc_variant.modules) - 1

        # skip cost at end consists of two parts:
        # first, initiation cost
        skip_cost = WEIGHT_PARAMS[ST.CHOOSE_IF_FINISH_MATCHING][relevant_genomic_context][ET.FINISH_MATCHING]
        # second, elongation
        skip_cost += get_full_skip_cost(fst_skipped_module_idx,
                                       lst_skipped_module_idx,
                                       ST.CHOOSE_IF_FINISH_MATCHING)
        weight = (
            skip_cost
            if edge_type == ET.FINISH_MATCHING
            else math.log(1 - math.exp(skip_cost))
        )
    elif (state.state_type == ST.SKIPPING_MODULES_AT_END):
        # in weight params for ST.SKIPPING_MODULES_AT_END
        # I store the skip cost for a module.
        # It is accounted for in the ET.FINISH_MATCHING edge
        # but the weight of the edge ET.SKIP_MODULES_AT_END in the HMM is zero (because it's the only edge)
        weight = 0
    else:
        # general case -- just lookup the parameters
        weight = WEIGHT_PARAMS[state.state_type][relevant_genomic_context][edge_type]


    return DetailedHMMEdge(edge_type=edge_type, weight=weight)
