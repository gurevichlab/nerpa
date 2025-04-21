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
from src.matching.hmm_config import EDGE_TYPE_DEPENDENCIES
from src.training.hmm_parameters.training_types import EdgeInfo
from src.write_results import write_yaml


# leave only relevant genomic context features for each edge type
def get_filtered_edge_choices(edge_choices: List[Tuple[EdgeInfo, bool]]) \
        -> List[Tuple[EdgeInfo, bool]]:
    filtered_edge_choices: List[Tuple[EdgeInfo, bool]] = []
    for (edge_type, genomic_context), chosen in edge_choices:
        # weights of edges of these types are automatically computed based on other edges,
        # so they are not relevant
        if edge_type not in EDGE_TYPE_DEPENDENCIES:
            continue

        genomic_context = tuple(feature
                                for feature in genomic_context
                                if feature in EDGE_TYPE_DEPENDENCIES[edge_type])

        new_edge_info = EdgeInfo(edge_type=edge_type, genomic_context=genomic_context)
        filtered_edge_choices.append((new_edge_info, chosen))

    #write_yaml(filtered_edge_choices, Path('filtered_data.yaml'))
    return filtered_edge_choices
