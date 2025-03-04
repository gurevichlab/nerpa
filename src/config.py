import time
from typing import Dict, List, Literal, Optional
from dataclasses import dataclass
from pathlib import Path
from src.pipeline.command_line_args_helper import CommandLineArgs
import yaml
import dacite
from datetime import datetime


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
    KNOWN_AA10_CODES: Dict[str, List[str]]
    KNOWN_AA34_CODES: Dict[str, List[str]]
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
        aa_codes_dict = yaml.safe_load((nerpa_dir / Path(cfg_dict['aa_codes'])).open('r'))
        self.KNOWN_AA10_CODES = aa_codes_dict['aa10']
        self.KNOWN_AA34_CODES = aa_codes_dict['aa34']


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
    default_output_file_name: str  # I have name instead of file due to rBAN quirks
    putative_hybrids_input_file: Path
    putative_hybrids_output_file_name: str
    rban_output_dir: Path

    def __init__(self,
                 rban_output_cfg_dict: dict,
                 main_out_dir: Path):
        for k, v in rban_output_cfg_dict.items():
            if not k.endswith('name'):
                setattr(self, k, main_out_dir / Path(v))
        self.default_output_file_name = rban_output_cfg_dict['default_output_file_name']
        self.putative_hybrids_output_file_name = rban_output_cfg_dict['putative_hybrids_output_file_name']


@dataclass
class rBAN_Processing_Config:
    MIN_RECOGNIZED_NODES: int
    PNP_BONDS: List[str]


@dataclass
class HMM_Scoring_Config:
    edge_weight_parameters: Dict[str, Dict[Optional[str], float]]
    emission_parameters: Dict[str, Dict[str, float]]


@dataclass
class CppIOConfig:
    hmms_json: Path
    nrp_linearizations_json: Path
    cpp_output_json: Path

    def __init__(self,
                 cpp_cfg_dict: dict,
                 main_out_dir: Path):
            for k, v in cpp_cfg_dict.items():
                setattr(self, k, main_out_dir / Path(v))

    def to_dict(self):
        return {k: str(v) for k, v in self.__dict__.items()}


@dataclass
class MatchingConfig:
    max_num_matches_per_bgc: int
    max_num_matches_per_nrp: int  # 0 means no limit
    max_num_matches: int
    checkpoints_heuristic: bool

    def __init__(self,
                 cfg_dict: dict,
                 args: Optional[CommandLineArgs]):
        for k, v in cfg_dict.items():
            setattr(self, k, v)
        if args is None:
            return

        if args.max_num_matches_per_bgc is not None:
            self.max_num_matches_per_bgc = args.max_num_matches_per_bgc
        if args.max_num_matches_per_nrp is not None:
            self.max_num_matches_per_nrp = args.max_num_matches_per_nrp
        if args.max_num_matches is not None:
            self.max_num_matches = args.max_num_matches
        if args.fast_matching is not None:
            self.checkpoints_heuristic = args.fast_matching


@dataclass
class OutputConfig:
    main_out_dir: Path
    symlink_to_latest: Path
    configs_output: Path
    antismash_out_dir: Path
    bgc_variants_dir: Path
    nrp_variants_dir: Path
    nrp_images_dir: Path
    rban_graphs: Path
    draw_molecules: bool
    matches_details: Path
    report: Path
    html_report: Path
    logo: Path  # seems a bit out of place here
    rban_output_config: rBAN_Output_Config
    cpp_io_config: CppIOConfig

    def __init__(self,
                 output_cfg_dict: dict,
                 nerpa_dir: Path,
                 main_out_dir: Path,
                 args: Optional[CommandLineArgs]):
        for k, v in output_cfg_dict.items():
            if k not in ('rban_output_config',
                         'draw_molecules',
                         'cpp_io_config'):
                setattr(self, k, main_out_dir / Path(v))
        self.rban_output_config = rBAN_Output_Config(output_cfg_dict['rban_output_config'],
                                                     main_out_dir)
        self.cpp_io_config = CppIOConfig(output_cfg_dict['cpp_io_config'],
                                         main_out_dir)
        if args is not None and args.skip_molecule_drawing is not None:
            self.draw_molecules = not args.skip_molecule_drawing
        else:
            self.draw_molecules = True

        self.main_out_dir = main_out_dir
        self.logo = nerpa_dir / Path(output_cfg_dict['logo'])
        self.symlink_to_latest = nerpa_dir / Path(output_cfg_dict['symlink_to_latest'])
        if self.symlink_to_latest == self.main_out_dir:
            raise ValueError(f'Invalid output directory: the path {self.symlink_to_latest} is reserved')


@dataclass
class Config:
    nerpa_dir: Path
    configs_dir: Path
    monomers_config: Path
    antismash_processing_config: antiSMASH_Processing_Config
    specificity_prediction_config: SpecificityPredictionConfig
    rban_config: rBAN_Config
    rban_processing_config: rBAN_Processing_Config
    hmm_scoring_config: Path
    matching_config: MatchingConfig
    cpp_matcher_exec: Path
    output_config: OutputConfig


def get_default_output_dir(nerpa_dir: Path,
                           cfg: dict) -> Path:
    while (out_dir := (nerpa_dir / Path(cfg['default_results_root_dirname']) /
              Path(datetime.now().strftime('%Y-%m-%d_%H-%M-%S')))).exists():
        time.sleep(1)
    return out_dir



def load_config(args: Optional[CommandLineArgs] = None) -> Config:
    nerpa_dir = Path(__file__).parent.parent.resolve()
    configs_dir = nerpa_dir / Path('configs')
    cfg = yaml.safe_load((configs_dir / 'config.yaml').open('r'))
    main_out_dir = args.output_dir.resolve() \
        if args is not None and args.output_dir is not None else get_default_output_dir(nerpa_dir, cfg)

    antismash_processing_cfg_dict = yaml.safe_load((nerpa_dir / cfg['antismash_processing_config']).open('r'))
    antismash_processing_cfg = dacite.from_dict(antiSMASH_Processing_Config, antismash_processing_cfg_dict)

    specificity_prediction_cfg_dict = yaml.safe_load((nerpa_dir / cfg['specificity_prediction_config']).open('r'))
    specificity_prediction_cfg = SpecificityPredictionConfig(nerpa_dir, specificity_prediction_cfg_dict)

    rban_cfg_dict = yaml.safe_load((nerpa_dir / cfg['rban_config']).open('r'))
    rban_cfg = rBAN_Config(rban_cfg_dict, nerpa_dir)

    rban_processing_cfg_dict = yaml.safe_load((nerpa_dir / cfg['rban_processing_config']).open('r'))
    rban_processing_cfg = dacite.from_dict(rBAN_Processing_Config, rban_processing_cfg_dict)

    output_cfg_dict = yaml.safe_load((nerpa_dir / cfg['output_config']).open('r'))
    output_cfg = OutputConfig(output_cfg_dict, nerpa_dir, main_out_dir, args)

    matching_cfg_dict = yaml.safe_load((nerpa_dir / cfg['matching_config']).open('r'))
    matching_cfg = MatchingConfig(matching_cfg_dict, args)

    return Config(antismash_processing_config=antismash_processing_cfg,
                  configs_dir=configs_dir,
                  nerpa_dir=nerpa_dir,
                  specificity_prediction_config=specificity_prediction_cfg,
                  rban_config=rban_cfg,
                  rban_processing_config=rban_processing_cfg,
                  hmm_scoring_config=nerpa_dir / cfg['hmm_scoring_config'],
                  output_config=output_cfg,
                  monomers_config=nerpa_dir / cfg['monomers_config'],
                  matching_config=matching_cfg,
                  cpp_matcher_exec=nerpa_dir / cfg['cpp_matcher_exec'])