from typing import (
    Dict,
    List,
    Tuple,
    Optional,
    Union
)
from src.antismash_parsing.location_features import ModuleLocFeature, GeneLocFeature, BGC_Fragment_Loc_Feature
from src.matching.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    GenomicContext,
)


def filter_context(genomic_context: GenomicContext) -> GenomicContext:
    # if ModuleLocFeature.START_OF_BGC is present, then MODULE_LOC_FEATURE.START_OF_FRAGMENT is automatically present
    # I want to keep only ModuleLocFeature.START_OF_BGC.
    # The same for other pairs of features
    if ModuleLocFeature.START_OF_FRAGMENT in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleLocFeature.START_OF_GENE)
    if ModuleLocFeature.END_OF_FRAGMENT in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleLocFeature.END_OF_GENE)
    if ModuleLocFeature.START_OF_BGC in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleLocFeature.START_OF_FRAGMENT)
    if ModuleLocFeature.END_OF_BGC in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleLocFeature.END_OF_FRAGMENT)
    return genomic_context


# parameters are estimated for single features only
# for multiple features, parameters are computed based on the single features
SingleFeatureContext = Union[ModuleLocFeature, GeneLocFeature, BGC_Fragment_Loc_Feature, None]


EDGE_TYPE_DEPENDENCIES = {
    # Insertion edge dependencies
    DetailedHMMEdgeType.START_INSERTING_AT_START: {
        ModuleLocFeature.START_OF_BGC,
        ModuleLocFeature.START_OF_FRAGMENT,
        ModuleLocFeature.START_OF_GENE,
        ModuleLocFeature.PKS_UPSTREAM_PREV_GENE,
        ModuleLocFeature.PKS_UPSTREAM_SAME_GENE
    },
    DetailedHMMEdgeType.INSERT_AT_START: set(),
    DetailedHMMEdgeType.START_INSERTING: {
        ModuleLocFeature.END_OF_BGC,
        ModuleLocFeature.END_OF_FRAGMENT,
        ModuleLocFeature.END_OF_GENE,
        ModuleLocFeature.PKS_DOWNSTREAM_NEXT_GENE,
        ModuleLocFeature.PKS_DOWNSTREAM_SAME_GENE
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
