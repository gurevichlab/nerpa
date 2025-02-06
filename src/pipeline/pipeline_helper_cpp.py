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



class PipelineHelperCpp:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    monomer_names_helper: MonomerNamesHelper

    def dump_hmms(self, hmms: List[DetailedHMM]) -> Path:
        with open(self.config.output_config.hmms_file, 'w') as f:
            f.write('\n\n'.join(str(hmm) for hmm in hmms))

    def dump_nrp_linearizations(self, nrp_linearizations: Dict[str, NRP_Linearizations]) -> Path:
        with open(self.config.output_config.nrp_linearizations_file, 'w') as f:
            for nrp_id, linearizations in nrp_linearizations.items():
                f.write(f'{nrp_id}\n')
                f.write(nrp_linearizations_to_str(linearizations, self.monomer_names_helper))
                f.write('\n\n')

    def run_cpp_matcher(self,
                        hmms_file: Path,
                        nrp_linearizations_file: Path):
        pass
