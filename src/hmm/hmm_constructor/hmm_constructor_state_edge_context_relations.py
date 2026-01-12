from __future__ import annotations
from typing import Tuple, Dict, Set, Iterable

from more_itertools import powerset

from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
)
from collections import defaultdict

GenomicContext = Tuple[ModuleGenomicContextFeature, ...]
GC = GenomicContext
ET = DetailedHMMEdgeType
ST = DetailedHMMStateType
MCF = ModuleGenomicContextFeature

STATE_TYPE_TO_EDGE_TYPES: Dict[DetailedHMMStateType, Set[DetailedHMMEdgeType]] = {
    # 0. Before matching
    ST.INITIAL: {ET.INITIAL},
    ST.CHOOSE_IF_START_MATCHING: {
        ET.START_MATCHING,
        ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START
    },
    ST.SKIPPING_MODULES_AT_START_FINISHED: {
        ET.START_INSERTING_AT_START,
        ET.NO_INSERTIONS
    },
    ST.INSERT_AT_START: {
        ET.INSERT_AT_START,
        ET.FINISH_INSERTING
    },

    # 1. During matching
    # 1.0 Ordinary flow
    ST.MODULE_SUBGRAPH_ROOT: {
        ET.MATCH,
        ET.SKIP_MODULE
    },
#    ST.MATCH: {
#        ET.START_INSERTING,
#        ET.NO_INSERTIONS
#    },
    ST.INSERT: {
        ET.INSERT,
        ET.FINISH_INSERTING
    },
    # 1.1 Iterations
    ST.CHOOSE_IF_ITERATE_GENE: {
        ET.ITERATE_GENE,
        ET.NO_ITERATION,
    },
    ST.CHOOSE_IF_ITERATE_MODULE: {
        ET.ITERATE_MODULE,
        ET.NO_ITERATION,
    },
    # 1.2 Tentative assembly line ends
    ST.MATCH_POSSIBLE_ASSEMBLY_FINISH: {
        ET.CHOOSE_IF_FINISH_MATCHING,
    },
    ST.CHOOSE_IF_FINISH_MATCHING: {
        ET.FINISH_MATCHING,
        ET.CONTINUE_MATCHING,
    },
    ST.MATCHING_CONTINUED: {
        ET.START_INSERTING,
        ET.NO_INSERTIONS
    },
    ST.MATCH_LAST_MODULE: {
        ET.FINISH_MATCHING,
    },
    # 2. After matching
    ST.MATCHING_FINISHED: {
        ET.START_INSERTING,
        ET.NO_INSERTIONS,
    },
    ST.INSERT_AT_END: {
        ET.INSERT_AT_END,
        ET.FINISH_INSERTING,
    },
    ST.SKIPPING_MODULES_AT_END: {
        ET.SKIP_MODULES_AT_END,
    },
    ST.FINAL: set(),

    # 3. Aunxiliary
    #ST.SKIP_MODULE_AT_START: {ET.SKIP_MODULE},
    #ST.SKIP_MODULE_AT_END: {ET.SKIP_MODULE},
}


NORMAL_FLOW_EDGE_TYPES: Set[DetailedHMMEdgeType] = {
    ET.INITIAL,
    ET.START_MATCHING,
    ET.NO_INSERTIONS,
    ET.FINISH_INSERTING,
    ET.MATCH,
    ET.NO_ITERATION,
    ET.CONTINUE_MATCHING,
}

RELEVANT_GENOMIC_CONTEXT: Dict[DetailedHMMStateType, Set[MCF]] = \
    defaultdict(lambda: set(),
                {
                    ST.CHOOSE_IF_START_MATCHING: {MCF.START_OF_BGC, MCF.ONLY_A_DOMAIN, MCF.ASSEMBLY_LINE_START_MARKER},  # choosing if start matching
                    #ST.SKIP_MODULE_AT_START: {MCF.ONLY_A_DOMAIN},
                    ST.SKIPPING_MODULES_AT_START_FINISHED: {MCF.PKS_UPSTREAM},  # choosing if insert at start
                    ST.INSERT_AT_START: {MCF.PKS_UPSTREAM},  # choosing if insert at start

                    ST.MODULE_SUBGRAPH_ROOT: {MCF.ONLY_A_DOMAIN},  # choosing if skip module
                    #ST.MATCH: {MCF.PKS_DOWNSTREAM, MCF.END_OF_GENE},  # choosing if insert
                    ST.MATCHING_CONTINUED: {MCF.PKS_DOWNSTREAM, MCF.END_OF_GENE},  # choosing if insert
                    ST.MATCHING_FINISHED: {MCF.PKS_DOWNSTREAM},  # choosing if insert
                    ST.INSERT: {MCF.PKS_DOWNSTREAM, MCF.END_OF_GENE},  # choosing if insert
                    ST.INSERT_AT_END: {MCF.PKS_DOWNSTREAM},  # choosing if insert at end

                    #ST.SKIP_MODULE_AT_END: {MCF.ONLY_A_DOMAIN},  # choosing if skip module
                    ST.CHOOSE_IF_FINISH_MATCHING: {MCF.ASSEMBLY_LINE_FINISH_MARKER},  # choosing if finish matching
                    ST.SKIPPING_MODULES_AT_END: {MCF.ONLY_A_DOMAIN},  # context not needed for construction but is used in training
                }
                )

def genomic_context_is_possible(state_type: DetailedHMMStateType,
                                gc: GenomicContext) -> bool:
    if any(f not in RELEVANT_GENOMIC_CONTEXT[state_type]
           for f in gc):
        return False

    match state_type:
        case ST.CHOOSE_IF_START_MATCHING:
            if MCF.START_OF_BGC in gc and MCF.ASSEMBLY_LINE_START_MARKER not in gc:
                return False
        case _:
            pass

    return True


def all_possible_genomic_contexts(state_type: DetailedHMMStateType) -> Iterable[GenomicContext]:
    mgfs_sets = powerset(RELEVANT_GENOMIC_CONTEXT[state_type])
    gcs = [tuple(sorted(mgfs)) for mgfs in mgfs_sets]
    return filter(lambda gc: genomic_context_is_possible(state_type, gc),
                  gcs)


MATCHING_STATE_TYPES: Set[DetailedHMMStateType] = {
    DetailedHMMStateType.MATCH,
    DetailedHMMStateType.MATCH_POSSIBLE_ASSEMBLY_FINISH,
    DetailedHMMStateType.MATCH_LAST_MODULE,
}

INSERT_STATE_TYPES: Set[DetailedHMMStateType] = {
    DetailedHMMStateType.INSERT,
    DetailedHMMStateType.INSERT_AT_START,
    DetailedHMMStateType.INSERT_AT_END,
}