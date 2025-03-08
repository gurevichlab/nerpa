from typing import (
    Dict,
    List,
    Tuple,
    Optional,
    Union
)
from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature, GeneGenomicContextFeature, FragmentGenomicContextFeature
from src.matching.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    GenomicContext,
)


def filter_context(genomic_context: GenomicContext) -> GenomicContext:
    # if ModuleLocFeature.START_OF_BGC is present, then MODULE_LOC_FEATURE.START_OF_FRAGMENT is automatically present
    # I want to keep only ModuleLocFeature.START_OF_BGC.
    # The same for other pairs of features
    if ModuleGenomicContextFeature.START_OF_FRAGMENT in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleGenomicContextFeature.START_OF_GENE)
    if ModuleGenomicContextFeature.END_OF_FRAGMENT in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleGenomicContextFeature.END_OF_GENE)
    if ModuleGenomicContextFeature.START_OF_BGC in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleGenomicContextFeature.START_OF_FRAGMENT)
    if ModuleGenomicContextFeature.END_OF_BGC in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleGenomicContextFeature.END_OF_FRAGMENT)
    return genomic_context


# parameters are estimated for single features only
# for multiple features, parameters are computed based on the single features
SingleFeatureContext = Union[ModuleGenomicContextFeature, GeneGenomicContextFeature, FragmentGenomicContextFeature, None]


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

def get_filtered_edge_data(edge_choices: List[Tuple[DetailedHMMEdgeType, GenomicContext, bool]]) \
        -> Dict[DetailedHMMEdgeType, Dict[Tuple[SingleFeatureContext, ...], Tuple[int, int]]]:
    filtered_data = {}
    for edge_type, context, chosen in edge_choices:
        if edge_type not in EDGE_TYPE_DEPENDENCIES:
            continue
        if context is None:
            context = ()
        context = filter_context(context)
        context = tuple(feature for feature in context if feature in EDGE_TYPE_DEPENDENCIES[edge_type])
        if edge_type not in filtered_data:
            filtered_data[edge_type] = {}
        if context not in filtered_data[edge_type]:
            filtered_data[edge_type][context] = (0, 0)
        filtered_data[edge_type][context] = (filtered_data[edge_type][context][0] + 1 - int(chosen),
                                             filtered_data[edge_type][context][1] + int(chosen))

    # write_yaml(filtered_data, Path('filtered_data.yaml'))
    return filtered_data
