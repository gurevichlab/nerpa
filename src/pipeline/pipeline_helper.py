from typing import (
    Dict,
    List,
    Tuple
)

import yaml

from src.generic.functional import timing_decorator
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
from src.pipeline.pipeline_helper_cpp import PipelineHelperCpp
from src.matching.hmm_config import load_hmm_scoring_config
from src.matching.hmm_scoring_helper import HMMHelper
import src.write_results as report

import shutil
import pandas as pd
from src.pipeline.pipeline_helper_antismash import PipelineHelper_antiSMASH
from src.rban_parsing.get_linearizations import get_all_nrp_linearizations, NRP_Linearizations
from src.matching.hmm_match import HMM_Match, convert_to_detailed_matches
from pathlib import Path



class PipelineHelper:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    monomer_names_helper: MonomerNamesHelper
    pipeline_helper_rban: PipelineHelper_rBAN
    pipeline_helper_antismash: PipelineHelper_antiSMASH
    pipeline_helper_cpp: PipelineHelperCpp

    def __init__(self, log: NerpaLogger):
        self.log = log

        try:
            self.args = get_command_line_args()
        except ValidationError as e:
            self.log.error(str(e), to_stderr=True)
            raise e

        try:
            self.config = load_config(self.args)
        except ValueError as e:
            self.log.error(str(e), to_stderr=True)
            raise e
        nerpa_utils.set_up_output_dir(output_cfg=self.config.output_config,
                                      crash_if_exists=not self.args.output_dir_reuse,
                                      log=self.log)

        self.log.set_up_file_handler(self.config.output_config.main_out_dir)
        self.log.start()

        shutil.copytree(self.config.configs_dir, self.config.output_config.configs_output, copy_function=shutil.copy)

        monomers_cfg = yaml.safe_load(self.config.monomers_config.open('r'))
        monomers_table_tsv = self.config.nerpa_dir / monomers_cfg['monomer_names_table']
        monomer_names_helper = MonomerNamesHelper(pd.read_csv(monomers_table_tsv, sep='\t'),
                                                  monomers_cfg['supported_residues'],
                                                  monomers_cfg['pks_names'])
        self.monomer_names_helper = monomer_names_helper
        self.pipeline_helper_rban = PipelineHelper_rBAN(self.config, self.args, self.log, monomer_names_helper)
        self.pipeline_helper_antismash = PipelineHelper_antiSMASH(self.config, self.args, monomer_names_helper, self.log)
        self.pipeline_helper_cpp = PipelineHelperCpp(self.config, self.args, self.log, monomer_names_helper)
        hmm_scoring_config = load_hmm_scoring_config(self.config.hmm_scoring_config)
        DetailedHMM.hmm_helper = HMMHelper(hmm_scoring_config, monomer_names_helper)

    @timing_decorator
    def get_bgc_variants(self) -> List[BGC_Variant]:
        return self.pipeline_helper_antismash.get_bgc_variants()

    @timing_decorator
    def get_nrp_variants_and_rban_records(self) -> Tuple[List[NRP_Variant], List[Parsed_rBAN_Record]]:
        if self.pipeline_helper_rban.preprocessed_nrp_variants():
            rban_records = []
            nrp_variants = self.pipeline_helper_rban.load_nrp_variants()
        else:
            rban_records = self.pipeline_helper_rban.get_rban_results()
            nrp_variants = self.pipeline_helper_rban.get_nrp_variants(rban_records)
        return nrp_variants, rban_records

    @timing_decorator
    def construct_hmms(self, bgc_variants: List[BGC_Variant]) -> List[DetailedHMM]:
        self.log.info("\n======= Constructing HMMs")
        return [DetailedHMM.from_bgc_variant(bgc_variant) for bgc_variant in bgc_variants]

    @timing_decorator
    def get_nrp_linearizations(self, nrp_variants: List[NRP_Variant]) \
            -> List[NRP_Linearizations]:
        self.log.info("\n======= Generating NRP linearizations")
        return get_all_nrp_linearizations(nrp_variants)

    @timing_decorator
    def get_matches(self,
                    hmms: List[DetailedHMM],
                    nrp_linearizations: List[NRP_Linearizations],
                    nrp_variants: List[NRP_Variant]) -> List[Match]:
        #for i, hmm in enumerate(hmms):
        #    hmm.draw(Path(f'{hmm.bgc_variant.genome_id}.png'))
        self.log.info("\n======= Nerpa matching")
        hmm_matches = self.pipeline_helper_cpp.get_hmm_matches(hmms, nrp_linearizations)
        return convert_to_detailed_matches(hmms, nrp_variants, hmm_matches)

    def write_results(self,
                      matches: List[Match],
                      bgc_variants: List[BGC_Variant],
                      nrp_variants: List[NRP_Variant],
                      rban_records: List[Parsed_rBAN_Record],
                      matches_details: bool = True):
        self.log.info("RESULTS:")
        self.log.info("Main report is saved to " + str(self.config.output_config.report), indent=1)
        self.log.info("HTML report is saved to " + str(self.config.output_config.html_report), indent=1)
        self.log.info("Detailed reports are saved to " + str(self.config.output_config.matches_details), indent=1)
        report.write_results(matches, self.config.output_config,
                             bgc_variants, nrp_variants,
                             rban_records,
                             matches_details,
                             log=self.log)
        # q: create symlink to results in self.config.default_results_root_dirname

        self.log.finish()  # TODO: create a separate method for this and "cleaning up"


