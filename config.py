from typing import List
from dataclasses import dataclass
from pathlib import Path
import nerpa_init
from command_line_args_helper import CommandLineArgs
from src.nerpa_pipeline.rban_parser import rBAN_Residue_Name
import yaml
import dacite


@dataclass
class rBAN_Config:
    rban_dir: Path
    rban_jar: Path
    rban_out_dir: Path

    default_monomers_file: Path
    default_input_file: Path
    default_output_file_name: str  # this is how rBAN is invoked for some reason

    putative_hybrids_input_file: Path
    putative_hybrids_output_file_name: str


@dataclass
class rBAN_Processing_Config:
    SUPPORTED_RESIDUES: List[rBAN_Residue_Name]
    MIN_RECOGNIZED_NODES: int
    CUT_LIPIDS: bool
    PNP_BONDS: List[str]


@dataclass
class ConfigPaths:
    main_out_dir: Path
    antismash_out_dir: Path
    nerpa_monomers: Path  # TODO: rename
    nerpa_monomers_info: Path
    configs_input: Path
    configs_output: Path
    scoring_config: Path
    report: Path
    matches_details: Path


@dataclass
class Config:
    paths: ConfigPaths
    rban_config: rBAN_Config
    rban_processing_config: rBAN_Processing_Config



def load_config(args: CommandLineArgs) -> Config:
    rban_dir = Path(nerpa_init.external_tools_dir) / Path('rBAN')
    rban_jar = rban_dir / Path('rBAN-1.0.jar')

    main_out_dir = args.output_dir.resolve()

    nerpa_dir = Path(__file__).parent.resolve()
    cfg = yaml.safe_load((nerpa_dir / Path('config.yaml')).open('r'))

    paths_cfg_dict = cfg['paths']
    paths = ConfigPaths(main_out_dir=main_out_dir,
                        antismash_out_dir=main_out_dir / Path(paths_cfg_dict['antismash_out_dir']),
                        nerpa_monomers=nerpa_dir / Path(paths_cfg_dict['nerpa_monomers']),
                        nerpa_monomers_info=nerpa_dir / Path(paths_cfg_dict['nerpa_monomers_info']),
                        configs_input=args.configs_dir if args.configs_dir else \
                            nerpa_dir / Path(paths_cfg_dict['configs_input']),
                        configs_output = main_out_dir / Path('configs_output'),
                        scoring_config=nerpa_dir / Path(paths_cfg_dict['scoring_config']),
                        report=main_out_dir / Path(paths_cfg_dict['report']),
                        matches_details=main_out_dir / Path(paths_cfg_dict['matches_details']))

    rban_cfg_dict = cfg['rban_config']
    rban_config = rBAN_Config(rban_dir=rban_dir,
                              rban_jar=rban_jar,
                              rban_out_dir=main_out_dir,
                              default_monomers_file=main_out_dir / Path(rban_cfg_dict['default_monomers_file']),
                              default_input_file=main_out_dir / Path(rban_cfg_dict['default_input_file']),
                              default_output_file_name=rban_cfg_dict['default_output_file_name'],
                              putative_hybrids_input_file=main_out_dir / Path(rban_cfg_dict['putative_hybrids_input_file']),
                              putative_hybrids_output_file_name= rban_cfg_dict['putative_hybrids_output_file_name'])

    rban_processing_cfg = dacite.from_dict(rBAN_Processing_Config,
                                           cfg['rban_processing_config'])

    return Config(paths=paths,
                  rban_config=rban_config,
                  rban_processing_config=rban_processing_cfg)
