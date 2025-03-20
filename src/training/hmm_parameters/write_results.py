from src.training.hmm_parameters.training_types import DataForTraining
from src.write_results import write_yaml


def write_params(edge_params,
                 emission_params,
                 data_for_training: DataForTraining,
                 output_dir):
    output_dir.mkdir(exist_ok=True, parents=True)
    write_yaml(data_for_training, output_dir / 'data_for_training.yaml')

    if edge_params is not None:
        write_yaml(edge_params, output_dir / 'edge_params.yaml')
    if emission_params is not None:
        write_yaml(emission_params, output_dir / 'emission_params.yaml')
    print(f'Edge and emission params written to {output_dir}')
