from typing import (
    Dict,
    List,
    Tuple
)

import yaml

from src.pipeline.command_line_args_helper import (
    CommandLineArgs,
    get_command_line_args,
    ValidationError
)
from src.pipeline.logger import NerpaLogger
from src.config import Config, load_config
from src.data_types import (
    BGC_Variant,
    NRP_Variant
)
import src.pipeline.nerpa_utils as nerpa_utils
from src.rban_parsing.rban_parser import (
    Parsed_rBAN_Record
)
from src.monomer_names_helper import MonomerNamesHelper

from src.matching.match_type import Match, Match_NRP_Variant_Info
from src.matching.matcher import get_matches
from src.matching.detailed_hmm import DetailedHMM
from src.pipeline.pipeline_helper_rban import PipelineHelper_rBAN
from src.matching.hmm_config import load_hmm_scoring_config
from src.matching.hmm_scoring_helper import HMMHelper
import src.write_results as report

import shutil
import pandas as pd
from src.pipeline.pipeline_helper_antismash import PipelineHelper_antiSMASH
from src.rban_parsing.get_linearizations import (
    NRP_Linearizations,
    nrp_linearizations_to_str
)
from src.matching.hmm_match import HMM_Match, convert_to_detailed_matches
import json
from pathlib import Path
from dataclasses import asdict
import subprocess
from dataclasses import dataclass


@dataclass
class PipelineHelperCpp:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    monomer_names_helper: MonomerNamesHelper

    def dump_hmms(self, detailed_hmms: List[DetailedHMM]) -> Path:
        out_file = self.config.output_config.cpp_output_config.hmms_json
        with open(out_file, 'w') as f:
            json.dump([detailed_hmm.to_hmm().to_json() for detailed_hmm in detailed_hmms], f, indent=4)
        return out_file


    def dump_nrp_linearizations(self, nrp_linearizations: List[NRP_Linearizations]) -> Path:
        out_file = self.config.output_config.cpp_output_config.nrp_linearizations_json
        with open(out_file, 'w') as f:
            json.dump([nrp_linearization.to_mon_codes_json(self.monomer_names_helper)
                       for nrp_linearization in nrp_linearizations], f, indent=4)
        return out_file

    def run_cpp_matcher(self,
                        hmms_json: Path,
                        nrp_linearizations_json: Path) -> List[HMM_Match]:
        # 1. Dump cpp config
        cpp_config_dict = asdict(self.config.output_config.cpp_output_config)
        cpp_config_json = self.config.output_config.cpp_output_config.config_json
        with open(cpp_config_json, 'w') as f:
            json.dump(cpp_config_dict, f, indent=4)

        # 2. Run cpp matcher
        # Construct the command
        cmd = [
            self.config.cpp_matcher_exec,  # Path to compiled C++ executable
            "--hmms_json", hmms_json,
            "--nrps_json", nrp_linearizations_json,
            "--config_json", cpp_config_json,
            "--threads", str(self.args.num_threads),
            "--output", self.config.output_config.cpp_output_config.cpp_output_json
        ]
        print(f"Running C++ matcher with command: {' '.join(cmd)}")

        # Run the C++ executable
        result = subprocess.run(cmd, check=True)

        # 3. Collect the results
        with open(self.config.output_config.cpp_output_config.cpp_output_json, 'r') as f:
            hmm_matches = [HMM_Match.from_json(json_dict)
                           for json_dict in json.load(f)]
        return hmm_matches

    def get_hmm_matches(self,
                        detailed_hmms: List[DetailedHMM],
                        nrp_linearizations: List[NRP_Linearizations]) -> List[HMM_Match]:
        hmms_json = self.dump_hmms(detailed_hmms)
        nrp_linearizations_json = self.dump_nrp_linearizations(nrp_linearizations)
        return self.run_cpp_matcher(hmms_json, nrp_linearizations_json)

