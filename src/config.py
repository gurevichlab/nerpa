from typing import Dict, List, Literal
from dataclasses import dataclass
from pathlib import Path
from src.pipeline.command_line_args_helper import CommandLineArgs
from src.monomer_names_helper import NorineMonomerName, antiSMASH_MonomerName, MonomerResidue
from src.matching.scoring_config import load_scoring_config, ScoringConfig
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

    def __init__(self,
                 rban_cfg_dict: dict,
                 rban_jar: Path,
                 main_out_dir: Path):
        self.rban_dir = rban_jar.parent
        self.rban_jar = rban_jar
        self.rban_out_dir = main_out_dir
        self.default_monomers_file = main_out_dir / Path(rban_cfg_dict['default_monomers_file'])
        self.default_input_file = main_out_dir / Path(rban_cfg_dict['default_input_file'])
        self.default_output_file_name = rban_cfg_dict['default_output_file_name']
        self.putative_hybrids_input_file = main_out_dir / Path(rban_cfg_dict['putative_hybrids_input_file'])
        self.putative_hybrids_output_file_name = rban_cfg_dict['putative_hybrids_output_file_name']


@dataclass
class rBAN_Processing_Config:
    SUPPORTED_RESIDUES: List[NorineMonomerName]
    MIN_RECOGNIZED_NODES: int
    CUT_LIPIDS: bool
    PNP_BONDS: List[str]

    def __init__(self,
                 rban_processing_cfg_dict: dict):
        for k, v in rban_processing_cfg_dict.items():
            setattr(self, k, v)


@dataclass
class antiSMASH_Parsing_Config:
    ANTISMASH_DOMAINS_NAMES: Dict[str, str]
    MAX_DISTANCE_BETWEEN_GENES: int
    MAX_VARIANTS_PER_BGC: int
    MAX_PERMUTATIONS_PER_BGC: int
    MAX_BGC_SPLITS_INTO_FRAGMENTS: int
    SCORING_TABLE_COLUMNS: List[str]
    SCORING_TABLE_INDEX: str
    SVM_SUBSTRATES: List[MonomerResidue]
    SVM_NOT_SUPPORTED_SCORE: float
    SVM_NO_PREDICTION_SCORE: float
    KNOWN_AA10_CODES: Dict[MonomerResidue, List[str]] = None
    KNOWN_AA34_CODES: Dict[MonomerResidue, List[str]] = None


    def __init__(self,
                 antismash_parsing_cfg_dict: dict,
                 aa_codes: dict):
        for k, v in antismash_parsing_cfg_dict.items():
            setattr(self, k, v)
        self.KNOWN_AA10_CODES = aa_codes['aa10']
        self.KNOWN_AA34_CODES = aa_codes['aa34']


@dataclass
class ConfigPaths:
    main_out_dir: Path
    antismash_out_dir: Path
    aa_codes: Path
    specificity_prediction_model: Path
    nerpa_monomers: Path  # TODO: rename
    nerpa_monomers_info: Path
    configs_input: Path
    configs_output: Path
    scoring_config: Path
    report: Path
    rban_dir: Path
    matches_details: Path
    default_results_root_dirname: str
    default_results_dirname_prefix: str

    def __init__(self,
                 paths_cfg_dict: dict,
                 nerpa_dir: Path,
                 main_out_dir: Path,
                 configs_dir: Path = None):
        self.main_out_dir = main_out_dir
        self.antismash_out_dir = main_out_dir / Path(paths_cfg_dict['antismash_out_dir'])
        self.aa_codes = nerpa_dir / Path(paths_cfg_dict['aa_codes'])
        self.specificity_prediction_model = nerpa_dir / Path(paths_cfg_dict['specificity_prediction_model'])
        self.nerpa_monomers = nerpa_dir / Path(paths_cfg_dict['nerpa_monomers'])
        self.nerpa_monomers_info = nerpa_dir / Path(paths_cfg_dict['nerpa_monomers_info'])
        self.configs_input = configs_dir if configs_dir else nerpa_dir / Path(paths_cfg_dict['configs_input'])
        self.configs_output = main_out_dir / Path('configs_output')
        self.scoring_config = nerpa_dir / Path(paths_cfg_dict['scoring_config'])
        self.report = main_out_dir / Path(paths_cfg_dict['report'])
        self.html_report = main_out_dir / Path(paths_cfg_dict['html_report'])
        self.matches_details = main_out_dir / Path(paths_cfg_dict['matches_details'])
        self.rban_dir = nerpa_dir / Path(paths_cfg_dict['rban_dir'])
        self.default_results_root_dirname = paths_cfg_dict['default_results_root_dirname']
        self.default_results_dirname_prefix = paths_cfg_dict['default_results_dirname_prefix']


@dataclass
class SpecificityPredictionConfig:
    calibration: bool
    apply_step_function: bool
    calibration_step_function_steps: List[float]
    normalize_scores: bool
    apriori_residue_prob: Dict[str, float]
    pseudo_counts: bool
    pseudo_count_fraction: float
    compute_evidence: bool


@dataclass
class MatchingConfig:
    scoring_config: ScoringConfig
    heuristic_discard_on: bool


@dataclass
class Config:
    paths: ConfigPaths
    antismash_parsing_config: antiSMASH_Parsing_Config
    specificity_prediction_config: SpecificityPredictionConfig
    rban_config: rBAN_Config
    rban_processing_config: rBAN_Processing_Config
    matching_config: MatchingConfig



def load_config(args: CommandLineArgs) -> Config:
    nerpa_dir = Path(__file__).parent.parent.resolve()
    main_out_dir = args.output_dir.resolve()
    cfg = yaml.safe_load((nerpa_dir / Path('configs/config.yaml')).open('r'))

    paths_config = ConfigPaths(paths_cfg_dict=cfg['paths'],
                               nerpa_dir=nerpa_dir,
                               main_out_dir=main_out_dir,
                               configs_dir=args.configs_dir)
    rban_config = rBAN_Config(rban_cfg_dict=cfg['rban_config'],
                              rban_jar=paths_config.rban_dir / Path('rBAN-1.0.jar'),
                              main_out_dir=main_out_dir)

    rban_processing_config = rBAN_Processing_Config(cfg['rban_processing_config'])

    antismash_parsing_config = antiSMASH_Parsing_Config(cfg['antismash_parsing_config'],
                                                        yaml.safe_load(paths_config.aa_codes.open('r')))
    specificity_prediction_config = dacite.from_dict(SpecificityPredictionConfig,
                                                     cfg['specificity_prediction_config'])
    matching_config = MatchingConfig(heuristic_discard_on=args.heuristic_discard,
                                     scoring_config=load_scoring_config(paths_config.scoring_config))
    return Config(paths=paths_config,
                  antismash_parsing_config=antismash_parsing_config,
                  rban_config=rban_config,
                  rban_processing_config=rban_processing_config,
                  specificity_prediction_config=specificity_prediction_config,
                  matching_config=matching_config)
