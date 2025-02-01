from typing import Dict, List, Literal, Optional
from dataclasses import dataclass
from pathlib import Path
from src.pipeline.command_line_args_helper import CommandLineArgs
import yaml
import dacite


@dataclass
class antiSMASH_Processing_Config:
    ANTISMASH_DOMAINS_NAMES_MAPPING: Dict[str, str]
    MAX_DISTANCE_BETWEEN_GENES: int
    MAX_VARIANTS_PER_BGC: int
    MAX_PERMUTATIONS_PER_BGC: int
    MAX_BGC_SPLITS_INTO_FRAGMENTS: int


@dataclass
class SpecificityPredictionConfig:
    specificity_prediction_model: Path
    a_domains_signatures: Path
    aa_codes: Path
    SVM_SUBSTRATES: List[str]
    SVM_NOT_SUPPORTED_SCORE: float
    SVM_NO_PREDICTION_SCORE: float
    SCORING_TABLE_INDEX: str
    SCORING_TABLE_COLUMNS: List[str]
    apriori_residue_prob: Dict[str, float]
    calibration_step_function_steps: List[float]
    pseudo_count_fraction: float

    def __init__(self,
                 nerpa_dir: Path,
                 cfg_dict: dict):
        for k, v in cfg_dict.items():
            setattr(self, k, v)
        self.specificity_prediction_model = nerpa_dir / Path(cfg_dict['specificity_prediction_model'])
        self.a_domains_signatures = nerpa_dir / Path(cfg_dict['a_domains_signatures'])
        self.aa_codes = nerpa_dir / Path(cfg_dict['aa_codes'])


@dataclass
class rBAN_Config:
    rban_jar: Path
    nerpa_monomers: Path

    def __init__(self,
                 rban_cfg_dict: dict,
                 nerpa_dir: Path):
        self.rban_jar = nerpa_dir / Path(rban_cfg_dict['rban_jar'])
        self.nerpa_monomers = nerpa_dir / Path(rban_cfg_dict['nerpa_monomers'])


@dataclass
class rBAN_Output_Config:
    default_monomers_file: Path
    default_input_file: Path
    default_output_file_name: Path
    putative_hybrids_input_file: Path
    putative_hybrids_output_file_name: Path

    def __init__(self,
                 rban_output_cfg_dict: dict,
                 main_out_dir: Path):
        for k, v in rban_output_cfg_dict.items():
            setattr(self, k, main_out_dir / Path(v))


@dataclass
class rBAN_Processing_Config:
    MIN_RECOGNIZED_NODES: int
    PNP_BONDS: List[str]


@dataclass
class HMM_Scoring_Config:
    edge_weight_parameters: Dict[str, Dict[Optional[str], float]]
    emission_parameters: Dict[str, Dict[str, float]]


@dataclass
class OutputConfig:
    main_out_dir: Path
    default_results_root_dirname: Path
    default_results_dirname_prefix: Path
    configs_output: Path
    antismash_out_dir: Path
    bgc_variants_dir: Path
    nrp_variants_dir: Path
    nrp_images_dir: Path
    rban_graphs: Path
    matches_details: Path
    report: Path
    html_report: Path
    logo: Path
    rban_output_config: rBAN_Output_Config

    def __init__(self,
                 output_cfg_dict: dict,
                 main_out_dir: Path):
        for k, v in output_cfg_dict.items():
            if k != 'rban_output_config':
                print(k, v)
                setattr(self, k, main_out_dir / Path(v))
        self.rban_output_config = rBAN_Output_Config(output_cfg_dict['rban_output_config'],
                                                     main_out_dir)
        self.main_out_dir = main_out_dir


@dataclass
class Config:
    configs_dir: Path
    monomer_names_table: Path
    antismash_processing_config: antiSMASH_Processing_Config
    specificity_prediction_config: SpecificityPredictionConfig
    rban_config: rBAN_Config
    rban_processing_config: rBAN_Processing_Config
    hmm_scoring_config: Path
    output_config: OutputConfig


def load_config(args: CommandLineArgs) -> Config:
    nerpa_dir = Path(__file__).parent.parent.resolve()
    main_out_dir = args.output_dir.resolve()
    configs_dir = nerpa_dir / Path('configs')
    cfg = yaml.safe_load((configs_dir / 'config.yaml').open('r'))

    antismash_processing_cfg_dict = yaml.safe_load((nerpa_dir / cfg['antismash_processing_config']).open('r'))
    antismash_processing_cfg = dacite.from_dict(antiSMASH_Processing_Config, antismash_processing_cfg_dict)

    specificity_prediction_cfg_dict = yaml.safe_load((nerpa_dir / cfg['specificity_prediction_config']).open('r'))
    specificity_prediction_cfg = SpecificityPredictionConfig(nerpa_dir, specificity_prediction_cfg_dict)

    rban_cfg_dict = yaml.safe_load((nerpa_dir / cfg['rban_config']).open('r'))
    rban_cfg = rBAN_Config(rban_cfg_dict, nerpa_dir)

    rban_processing_cfg_dict = yaml.safe_load((nerpa_dir / cfg['rban_processing_config']).open('r'))
    rban_processing_cfg = dacite.from_dict(rBAN_Processing_Config, rban_processing_cfg_dict)

    output_cfg_dict = yaml.safe_load((nerpa_dir / cfg['output_config']).open('r'))
    output_cfg = OutputConfig(output_cfg_dict, main_out_dir)

    return Config(antismash_processing_config=antismash_processing_cfg,
                  configs_dir=configs_dir,
                  specificity_prediction_config=specificity_prediction_cfg,
                  rban_config=rban_cfg,
                  rban_processing_config=rban_processing_cfg,
                  hmm_scoring_config=nerpa_dir / cfg['hmm_scoring_config'],
                  output_config=output_cfg,
                  monomer_names_table=nerpa_dir / cfg['monomer_names_table'])
