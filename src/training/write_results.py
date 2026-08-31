from pathlib import Path
from typing import Dict

import yaml

from src.general_type_aliases import LogProb
from src.hmm.hmm_auxiliary_types import DetailedHMMEdgeType
from src.build_output.write_results import write_yaml
from src.training.hmm_parameters.short_aliases import *
from src.training.hmm_parameters.emission_probabilities.hmm_infer_emission_params import EmissionParams
from src.training.hmm_parameters.training_types import HMM_Chosen_Steps_df


def edge_params_to_yaml(
        edge_params: Dict[ST, Dict[GC, Dict[ET, LogProb]]]
) -> dict:
    edge_params_yaml: Dict[str, Dict[str, Dict[str, float]]] = {}

    for state_type, edges_by_gc in edge_params.items():
        edge_params_yaml[state_type.name] = {}
        for gc, edge_weights in edges_by_gc.items():
            features = sorted(f'"{f.name}"' for f in gc)
            gc_str = (
                f'({features[0]},)'
                if len(features) == 1 else
                '(' + ', '.join(features) + ')'
            )

            edge_params_yaml[state_type.name][gc_str] = {
                edge_type.name: weight
                for edge_type, weight in edge_weights.items()
            }


    return edge_params_yaml


def create_new_configs(old_configs_dir: Path,
                       config_updates_dir: Path,
                       new_configs_dir: Path):
     for old_config_path in old_configs_dir.glob('*.yaml'):
        old_config = yaml.safe_load(old_config_path.read_text())
        new_config = old_config.copy()

        update_path = config_updates_dir / old_config_path.name
        updates = (
            yaml.safe_load(update_path.read_text())
            if update_path.exists()
            else {}
        )

        new_config.update(updates)

        new_config_path = new_configs_dir / old_config_path.name
        write_yaml(new_config, new_config_path)


def write_params(
        edge_params: Dict[ST, Dict[GC, Dict[ET, LogProb]]],
        emission_params: EmissionParams,
        data_for_training: HMM_Chosen_Steps_df,
        old_configs_dir: Path,
        output_dir: Path
):
    output_dir.mkdir(exist_ok=True, parents=True)
    data_for_training.write_tsv(output_dir / 'data_for_training.tsv')

    write_yaml(emission_params._asdict(), output_dir / 'emission_params.yaml')

    config_updates_dir = output_dir / 'config_updates'
    config_updates_dir.mkdir(exist_ok=True, parents=True)

    write_yaml(
        {
        'CALIBRATION_STEP_FUNCTION_STEPS_PARAS': emission_params.step_function
        },
        config_updates_dir / 'specificity_prediction_config.yaml'
    )
    hmm_scoring_config = {
        'emission_parameters': {
            'modification_frequencies': emission_params.modifications_frequences
        },
        'edge_weight_parameters': edge_params_to_yaml(edge_params)
    }
    write_yaml(hmm_scoring_config, config_updates_dir / 'hmm_scoring_config.yaml')

    new_configs_dir = output_dir / 'new_configs'
    new_configs_dir.mkdir(exist_ok=True, parents=True)

    create_new_configs(
        old_configs_dir,
        config_updates_dir,
        new_configs_dir
    )
