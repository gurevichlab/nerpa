from typing import (
    List,
    Tuple
)
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

from src.matching.scoring_config import load_scoring_config
from src.matching.scoring_helper import ScoringHelper
from src.matching.alignment_types import Match
from src.matching.matcher import get_matches
from src.pipeline.pipeline_helper_rban import PipelineHelper_rBAN
import src.write_results as report

import shutil
import pandas as pd
from src.pipeline.pipeline_helper_antismash import PipelineHelper_antiSMASH




class PipelineHelper:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    monomer_names_helper: MonomerNamesHelper
    pipeline_helper_rban: PipelineHelper_rBAN
    pipeline_helper_antismash: PipelineHelper_antiSMASH
    scoring_helper: ScoringHelper

    def __init__(self, log: NerpaLogger):
        self.log = log

        try:
            self.args = get_command_line_args()
        except ValidationError as e:
            self.log.error(str(e), to_stderr=True)
            raise e

        self.config = load_config(self.args)
        nerpa_utils.set_up_output_dir(output_dirpath=self.config.paths.main_out_dir,
                                      config=self.config,
                                      crash_if_exists=not self.args.output_dir_reuse,
                                      log=self.log)

        self.log.set_up_file_handler(self.config.paths.main_out_dir)
        self.log.start()

        shutil.copytree(self.config.paths.configs_input, self.config.paths.configs_output, copy_function=shutil.copy)

        monomer_names_helper = MonomerNamesHelper(pd.read_csv(self.config.paths.nerpa_monomers_info, sep='\t'))
        self.pipeline_helper_rban = PipelineHelper_rBAN(self.config, self.args, self.log, monomer_names_helper)
        self.pipeline_helper_antismash = PipelineHelper_antiSMASH(self.config, self.args, monomer_names_helper, self.log)
        self.scoring_helper = ScoringHelper(scoring_config=self.config.matching_config.scoring_config,
                                            heuristic_discard_on=self.config.matching_config.heuristic_discard_on)

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
        self.log.info("HTML report is saved to " + str(self.config.paths.html_report), indent=1)
        self.log.info("Detailed reports are saved to " + str(self.config.paths.matches_details), indent=1)
        report.write_results(matches, self.config.paths.main_out_dir,
                             bgc_variants, nrp_variants,
                             rban_records,
                             matches_details)
        self.log.finish()


