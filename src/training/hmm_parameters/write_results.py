from pathlib import Path
from typing import Dict, Optional

import yaml

from src.data_types import LogProb
from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType, GenomicContext
from src.training.hmm_parameters.hmm_infer_emission_params import EmissionParams
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


def create_new_configs(old_configs_dir: Path,
                          config_updates_dir: Path,
                          new_configs_dir: Path):
     for old_config_path in old_configs_dir.glob('*.yaml'):
        old_config = yaml.safe_load(old_config_path.read_text())
        new_config = old_config.copy()

        update_path = config_updates_dir / old_config_path.name
        updates = yaml.safe_load(update_path.read_text()) \
            if update_path.exists() else {}

        new_config.update(updates)

        new_config_path = new_configs_dir / old_config_path.name
        write_yaml(new_config, new_config_path)


def write_params(edge_params,
                 emission_params: EmissionParams,
                 data_for_training: DataForTraining,
                 old_configs_dir: Path,
                 output_dir: Path):
    output_dir.mkdir(exist_ok=True, parents=True)
    write_yaml(data_for_training, output_dir / 'data_for_training.yaml')

    config_updates_dir = output_dir / 'config_updates'
    config_updates_dir.mkdir(exist_ok=True, parents=True)

    write_yaml({'CALIBRATION_STEP_FUNCTION_STEPS_PARAS': emission_params.step_function},
               config_updates_dir / 'specificity_prediction_config.yaml')
    hmm_scoring_config = {
        'emission_parameters': {
            'modification_frequencies': emission_params.modifications_frequences
        },
        'edge_weight_parameters': edge_params_to_yaml(edge_params)
    }
    write_yaml(hmm_scoring_config, config_updates_dir / 'hmm_scoring_config.yaml')

    new_configs_dir = output_dir / 'new_configs'
    new_configs_dir.mkdir(exist_ok=True, parents=True)

    create_new_configs(old_configs_dir,
                       config_updates_dir,
                       new_configs_dir)