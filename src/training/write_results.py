from src.write_results import write_yaml
from pathlib import Path
import yaml


def write_params(edge_params, emission_params, output_dir):
    output_dir.mkdir(exist_ok=True, parents=True)
    write_yaml(edge_params, output_dir / 'edge_params.yaml')
    write_yaml(emission_params, output_dir / 'emission_params_old.yaml')
    print(f'Edge and emission params written to {output_dir}')
