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


# TODO: create GenomicContext already filtered?
def remove_redundancy(genomic_context: Optional[GenomicContext]) -> Optional[GenomicContext]:
    # if ModuleLocFeature.START_OF_BGC is present, then MODULE_LOC_FEATURE.START_OF_FRAGMENT is automatically present
    # I want to keep only ModuleLocFeature.START_OF_BGC.
    # The same for other pairs of features
    if genomic_context is None:
        return None

    if ModuleGenomicContextFeature.START_OF_FRAGMENT in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleGenomicContextFeature.START_OF_GENE)
    if ModuleGenomicContextFeature.END_OF_FRAGMENT in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleGenomicContextFeature.END_OF_GENE)
    if ModuleGenomicContextFeature.START_OF_BGC in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleGenomicContextFeature.START_OF_FRAGMENT)
    if ModuleGenomicContextFeature.END_OF_BGC in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleGenomicContextFeature.END_OF_FRAGMENT)
    return genomic_context


# which genomic context features are relevant for the edge probability for each edge type
EDGE_TYPE_DEPENDENCIES = {
    # Insertion edge dependencies
    DetailedHMMEdgeType.START_INSERTING_AT_START: {
        ModuleGenomicContextFeature.START_OF_BGC,
        ModuleGenomicContextFeature.START_OF_FRAGMENT,
        ModuleGenomicContextFeature.START_OF_GENE,
        ModuleGenomicContextFeature.PKS_UPSTREAM_PREV_GENE,
        ModuleGenomicContextFeature.PKS_UPSTREAM_SAME_GENE
    },
    DetailedHMMEdgeType.INSERT_AT_START: set(),
    DetailedHMMEdgeType.START_INSERTING: {
        ModuleGenomicContextFeature.END_OF_BGC,
        ModuleGenomicContextFeature.END_OF_FRAGMENT,
        ModuleGenomicContextFeature.END_OF_GENE,
        ModuleGenomicContextFeature.PKS_DOWNSTREAM_NEXT_GENE,
        ModuleGenomicContextFeature.PKS_DOWNSTREAM_SAME_GENE
    },
    DetailedHMMEdgeType.INSERT: set(),
    DetailedHMMEdgeType.END_INSERTING: set(),

    # Skip edge dependencies (empty for now)
    DetailedHMMEdgeType.START_SKIP_MODULES_AT_START: set(),
    DetailedHMMEdgeType.START_SKIP_GENES_AT_START: set(),
    DetailedHMMEdgeType.START_SKIP_FRAGMENTS_AT_START: set(),
    DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START: set(),
    DetailedHMMEdgeType.SKIP_MODULE: set(),
    DetailedHMMEdgeType.SKIP_GENE: set(),
    DetailedHMMEdgeType.SKIP_FRAGMENT: set(),
    DetailedHMMEdgeType.START_SKIPPING_AT_END: set(),
    DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END: set(),

    # Iteration edge dependencies (empty for now)
    DetailedHMMEdgeType.ITERATE_MODULE: set(),
    DetailedHMMEdgeType.ITERATE_GENE: set()
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

        genomic_context = remove_redundancy(genomic_context)
        if genomic_context is not None:
            genomic_context = tuple(feature
                                    for feature in genomic_context
                                    if feature in EDGE_TYPE_DEPENDENCIES[edge_type])

        new_edge_info = EdgeInfo(edge_type=edge_type, genomic_context=genomic_context)
        filtered_edge_choices.append((new_edge_info, chosen))

    # write_yaml(filtered_data, Path('filtered_data.yaml'))
    return filtered_edge_choices
