from collections import defaultdict
from typing import (
    List,
    Tuple,
    Union
)
from command_line_args_helper import (
    CommandLineArgs,
    get_command_line_args,
    ValidationError
)
from src.nerpa_pipeline.logger import NerpaLogger
import nerpa_init
from config import Config, load_config
from src.data_types import (
    BGC_Variant,
    SMILES,
    NRP_Variant
)

import src.nerpa_pipeline.nerpa_utils as nerpa_utils
from src.nerpa_pipeline import handle_rban
from src.nerpa_pipeline.rban_parser import (
    Parsed_rBAN_Record,
    parsed_chiralities,
    Raw_rBAN_Record
)
from src.nerpa_pipeline.rban_names_helper import rBAN_Names_Helper
from src.nerpa_pipeline.handle_monomers import get_monomers_chirality

from src.NewMatcher.scoring_config import load_scoring_config
from src.NewMatcher.scoring_helper import ScoringHelper
from src.NewMatcher.alignment_types import Match
from src.NewMatcher.matcher import get_matches
from pipeline_helper_rban import PipelineHelper_rBAN
import src.write_results as report

from rban_helper import rBAN_Helper
from pathlib import Path
import csv
import json
import shutil
import yaml
from dataclasses import dataclass
from itertools import chain
from src.pipeline_helpers.pipeline_helper_antismash import PipelineHelper_antiSMASH




class PipelineHelper:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    pipeline_helper_rban: PipelineHelper_rBAN
    pipeline_helper_antismash: PipelineHelper_antiSMASH
    scoring_helper: ScoringHelper

    def __init__(self, log: NerpaLogger):
        nerpa_init.init()
        self.log = log

        try:
            self.args = get_command_line_args()
        except ValidationError as e:
            self.log.error(str(e), to_stderr=True)
            raise e

        self.config = load_config(self.args)
        nerpa_utils.set_up_output_dir(output_dirpath=self.config.paths.main_out_dir,
                                      crash_if_exists=not self.args.output_dir_reuse,
                                      log=self.log)

        self.log.set_up_file_handler(self.config.paths.main_out_dir)
        self.log.start()

        shutil.copytree(self.config.paths.configs_input, self.config.paths.configs_output, copy_function=shutil.copy)

        self.pipeline_helper_rban = PipelineHelper_rBAN(self.config, self.args, self.log)
        self.pipeline_helper_antismash = PipelineHelper_antiSMASH(self.config, self.args, self.log)
        scoring_config = load_scoring_config(self.config.paths.scoring_config)  # TODO: this is ugly
        self.scoring_helper = ScoringHelper(scoring_config)

    def get_bgc_variants(self) -> List[BGC_Variant]:
        return self.pipeline_helper_antismash.get_bgc_variants()

    def get_nrp_variants_and_rban_records(self) -> Tuple[List[NRP_Variant], List[Parsed_rBAN_Record]]:
        if self.pipeline_helper_rban.preprocessed_nrp_variants():
            rban_records = []
            nrp_variants = self.pipeline_helper_rban.load_nrp_variants()
        else:
            rban_records = self.pipeline_helper_rban.get_rban_results()
            nrp_variants = self.pipeline_helper_rban.get_nrp_variants(rban_records)
        return nrp_variants, rban_records

    def get_matches(self,
                    bgc_variants: List[BGC_Variant],
                    nrp_variants: List[NRP_Variant]) -> List[Match]:
        self.log.info("\n======= Nerpa matching")
        return get_matches(bgc_variants, nrp_variants, self.scoring_helper,
                           min_score=self.args.min_score,
                           max_num_matches=self.args.num_matches,
                           num_threads=self.args.threads,
                           log=self.log)

    def write_results(self,
                      matches: List[Match],
                      bgc_variants: List[BGC_Variant],
                      nrp_variants: List[NRP_Variant],
                      rban_records: List[Parsed_rBAN_Record],
                      matches_details:  bool = True):
        self.log.info("RESULTS:")
        self.log.info("Main report is saved to " + str(self.config.paths.report), indent=1)
        self.log.info("Detailed reports are saved to " + str(self.config.paths.matches_details), indent=1)
        report.write_results(matches, self.config.paths.main_out_dir,
                             bgc_variants, nrp_variants,
                             rban_records,
                             matches_details)
        self.log.finish()


