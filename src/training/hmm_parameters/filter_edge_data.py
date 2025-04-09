from pathlib import Path
from typing import (
    List,
    Tuple,
    Optional
)
from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.matching.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    GenomicContext,
)
from src.training.hmm_parameters.training_types import EdgeInfo
from src.write_results import write_yaml

# which genomic context features are relevant for the edge probability for each edge type
EDGE_TYPE_DEPENDENCIES = {
    # Insertions
    DetailedHMMEdgeType.START_INSERTING_AT_START: {
        ModuleGenomicContextFeature.PKS_UPSTREAM,
    },
    DetailedHMMEdgeType.INSERT_AT_START: set(),
    DetailedHMMEdgeType.START_INSERTING: {
        ModuleGenomicContextFeature.END_OF_GENE,
        ModuleGenomicContextFeature.PKS_DOWNSTREAM,
    },
    DetailedHMMEdgeType.INSERT: set(),
    DetailedHMMEdgeType.END_INSERTING: set(),

    # Skips
    DetailedHMMEdgeType.START_SKIP_MODULES_AT_START: {ModuleGenomicContextFeature.ONLY_A_DOMAIN},
    DetailedHMMEdgeType.SKIP_MODULE_AT_START: {ModuleGenomicContextFeature.ONLY_A_DOMAIN},
    DetailedHMMEdgeType.SKIP_MODULE: {ModuleGenomicContextFeature.ONLY_A_DOMAIN},
    DetailedHMMEdgeType.SKIP_GENE: set(),
    DetailedHMMEdgeType.START_SKIP_MODULES_AT_END: {ModuleGenomicContextFeature.ONLY_A_DOMAIN},
    DetailedHMMEdgeType.SKIP_MODULE_AT_END: {ModuleGenomicContextFeature.ONLY_A_DOMAIN},

    # Iterations
    DetailedHMMEdgeType.ITERATE_MODULE: set(),
    DetailedHMMEdgeType.ITERATE_GENE: set(),

    # Other (not needed for training)
    #DetailedHMMEdgeType.START_MATCHING: set(),
    #DetailedHMMEdgeType.MATCH: set(),
    #DetailedHMMEdgeType.NO_INSERTIONS: set(),
}


# leave only relevant genomic context features for each edge type
def get_filtered_edge_choices(edge_choices: List[Tuple[EdgeInfo, bool]]) \
        -> List[Tuple[EdgeInfo, bool]]:
    filtered_edge_choices: List[Tuple[EdgeInfo, bool]] = []
    for (edge_type, genomic_context), chosen in edge_choices:
        # weights of edges of these types are automatically computed based on other edges,
        # so they are not relevant
        if edge_type not in EDGE_TYPE_DEPENDENCIES:
            continue
            #raise KeyError(f'Edge type {edge_type} is not in the list of edge types')

        if genomic_context is not None:
            genomic_context = tuple(feature
                                    for feature in genomic_context
                                    if feature in EDGE_TYPE_DEPENDENCIES[edge_type])

        new_edge_info = EdgeInfo(edge_type=edge_type, genomic_context=genomic_context)
        filtered_edge_choices.append((new_edge_info, chosen))

    #write_yaml(filtered_edge_choices, Path('filtered_data.yaml'))
    return filtered_edge_choices
