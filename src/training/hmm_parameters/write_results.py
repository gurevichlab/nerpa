from typing import Dict

import yaml

from src.data_types import LogProb
from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType, GenomicContext
from src.training.hmm_parameters.training_types import DataForTraining
from src.write_results import write_yaml


def edge_params_to_yaml(edge_params: Dict[DetailedHMMEdgeType,
                                     Dict[GenomicContext, LogProb]]) -> dict:
    edge_params_yaml: Dict[str, Dict[str, float]] = {}
    for edge_type, weights_dict in edge_params.items():
        edge_params_yaml[edge_type.name] = {}
        for gc, weight in weights_dict.items():
            gc_str = f'("{gc[0].name}",)' if len(gc) == 1 else \
                '(' + ', '.join(f'"{feature.name}"' for feature in gc) + ')'
            edge_params_yaml[edge_type.name][gc_str] = weight

    return edge_params_yaml


def write_params(edge_params,
                 emission_params,
                 data_for_training: DataForTraining,
                 output_dir):
    output_dir.mkdir(exist_ok=True, parents=True)
    write_yaml(data_for_training, output_dir / 'data_for_training.yaml')

    if edge_params is not None:
        yaml.dump(edge_params_to_yaml(edge_params),
                  (output_dir / 'edge_params.yaml').open('w'))
    if emission_params is not None:
        write_yaml(emission_params, output_dir / 'emission_params.yaml')
    print(f'Edge and emission params written to {output_dir}')
