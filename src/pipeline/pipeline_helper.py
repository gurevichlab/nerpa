from typing import (
    Dict,
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

from src.matching.matching_types_match import Match, Match_NRP_Variant_Info
from src.matching.hmm_matcher import get_matches
from src.matching.matcher_viterbi_detailed_hmm import DetailedHMM
from src.pipeline.pipeline_helper_rban import PipelineHelper_rBAN
from src.matching.hmm_config import load_hmm_scoring_config
from src.matching.hmm_scoring_helper import HMMHelper
import src.write_results as report

import shutil
import pandas as pd
from src.pipeline.pipeline_helper_antismash import PipelineHelper_antiSMASH
from src.rban_parsing.get_linearizations import get_all_nrp_linearizations, NRP_Linearizations



class PipelineHelper:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    monomer_names_helper: MonomerNamesHelper
    pipeline_helper_rban: PipelineHelper_rBAN
    pipeline_helper_antismash: PipelineHelper_antiSMASH

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
        self.monomer_names_helper = monomer_names_helper
        self.pipeline_helper_rban = PipelineHelper_rBAN(self.config, self.args, self.log, monomer_names_helper)
        self.pipeline_helper_antismash = PipelineHelper_antiSMASH(self.config, self.args, monomer_names_helper, self.log)
        hmm_scoring_config = load_hmm_scoring_config(self.config.paths.hmm_scoring_config)
        DetailedHMM.hmm_helper = HMMHelper(hmm_scoring_config, monomer_names_helper)

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

    def construct_hmms(self, bgc_variants: List[BGC_Variant]) -> List[DetailedHMM]:
        self.log.info("\n======= Constructing HMMs")
        return [DetailedHMM.from_bgc_variant(bgc_variant) for bgc_variant in bgc_variants]

    def get_nrp_linearizations(self, nrp_variants: List[NRP_Variant]) -> Dict[Match_NRP_Variant_Info, NRP_Linearizations]:
        self.log.info("\n======= Generating NRP linearizations")
        return get_all_nrp_linearizations(nrp_variants)

    def get_matches(self,
                    hmms: List[DetailedHMM],
                    nrp_linearizations: Dict[str, NRP_Linearizations]) -> List[Match]:
        self.log.info("\n======= Nerpa matching")
        return get_matches(hmms, nrp_linearizations,
                           max_num_matches_per_bgc_variant=self.args.num_matches,
                           num_threads=self.args.threads,
                           log=self.log)

    def write_results(self,
                      matches: List[Match],
                      bgc_variants: List[BGC_Variant],
                      nrp_variants: List[NRP_Variant],
                      rban_records: List[Parsed_rBAN_Record],
                      matches_details: bool = True):
        self.log.info("RESULTS:")
        self.log.info("Main report is saved to " + str(self.config.paths.report), indent=1)
        self.log.info("HTML report is saved to " + str(self.config.paths.html_report), indent=1)
        self.log.info("Detailed reports are saved to " + str(self.config.paths.matches_details), indent=1)
        report.write_results(matches, self.config.paths,
                             bgc_variants, nrp_variants,
                             rban_records,
                             matches_details,
                             log=self.log,
                             draw_molecules=self.args.draw_molecules)
        self.log.finish()  # TODO: create a separate method for this and "cleaning up"


