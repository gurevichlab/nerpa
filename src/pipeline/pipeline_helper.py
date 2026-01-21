from pathlib import Path
from typing import (
    List
)

import yaml

from src.aa_specificity_prediction_model.specificity_prediction_helper import SpecificityPredictionHelper
from src.generic.functional import timing_decorator
from src.pipeline.command_line_args_helper import (
    CommandLineArgs,
    get_command_line_args,
    ValidationError
)
from src.pipeline.logging.logger import NerpaLogger, PreliminaryLogger
from src.config import (
    Config,
    load_config,
    load_monomer_names_helper
)
from src.antismash_parsing.bgc_variant_types import BGC_Variant
from src.rban_parsing.nrp_variant_types import NRP_Variant
import src.pipeline.nerpa_utils as nerpa_utils
from src.pipeline.deduplication import cluster_isomorphic_nrp_variants, cluster_isomorphic_bgc_variants
from src.monomer_names_helper import MonomerNamesHelper

from src.matching.match_type import Match
from src.matching.matcher import get_hmm_matches
from src.hmm.detailed_hmm import DetailedHMM
from src.pipeline.pipeline_helper_rban import PipelineHelper_rBAN
from src.rban_parsing.nrp_variant_types import NRP_Variants_Info
from src.pipeline.pipeline_helper_cpp import PipelineHelperCpp
from src.hmm.hmm_scoring_config import load_hmm_scoring_config
from src.hmm.hmm_scoring_helper import HMMHelper
import src.build_output.write_results as write_results

import shutil
from src.pipeline.pipeline_helper_antismash import PipelineHelper_antiSMASH
from src.antismash_parsing.bgc_variant_types import BGC_Variants_Info
from src.pipeline.paras_parsing import get_paras_results_all
from src.rban_parsing.get_linearizations import get_all_nrp_linearizations, NRP_Linearizations
from src.matching.hmm_match import HMM_Match, convert_to_detailed_matches
from src.rban_parsing.rban_parser import Parsed_rBAN_Record


def load_parsed_rban_records(parsed_rban_records_path: Path) -> List[Parsed_rBAN_Record]:
    with parsed_rban_records_path.open('r') as f:
        return [Parsed_rBAN_Record.from_dict(record)
                for record in yaml.safe_load(f)]

class PipelineHelper:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    monomer_names_helper: MonomerNamesHelper
    hmm_helper: HMMHelper
    pipeline_helper_rban: PipelineHelper_rBAN
    pipeline_helper_antismash: PipelineHelper_antiSMASH
    pipeline_helper_cpp: PipelineHelperCpp

    def __init__(self, pre_logger: PreliminaryLogger):

        default_cfg = load_config()
        try:
            self.args = get_command_line_args(default_cfg)
        except ValidationError as e:
            raise

        try:
            self.config = load_config(self.args)
        except ValueError as e:
            raise e
        nerpa_utils.set_up_output_dir(output_cfg=self.config.output_config,
                                      crash_if_exists=not self.args.output_dir_reuse,
                                      log=pre_logger)

        self.log = NerpaLogger(self.config.logging_config)
        self.log.start()

        shutil.copytree(self.config.configs_dir, self.config.output_config.configs_output, copy_function=shutil.copy)

        self.monomer_names_helper = load_monomer_names_helper(self.config.monomers_config,
                                                              self.config.nerpa_dir)

        external_specificity_predictions = get_paras_results_all(self.args.paras_results,
                                                                 self.monomer_names_helper,
                                                                 self.log) \
            if self.args.paras_results is not None else None
        specificity_prediction_helper = SpecificityPredictionHelper(self.config.specificity_prediction_config,
                                                                    self.monomer_names_helper,
                                                                    external_specificity_predictions)

        hmm_scoring_config = load_hmm_scoring_config(self.config.nerpa_dir,
                                                     self.config.hmm_scoring_config,
                                                     specificity_prediction_helper,
                                                     self.monomer_names_helper)
        self.hmm_helper = HMMHelper(hmm_scoring_config, self.monomer_names_helper)

        self.pipeline_helper_rban = PipelineHelper_rBAN(self.config, self.args, self.log, self.monomer_names_helper)
        self.pipeline_helper_antismash = PipelineHelper_antiSMASH(self.config, self.args,
                                                                  self.monomer_names_helper,
                                                                  specificity_prediction_helper,
                                                                  self.log)
        self.pipeline_helper_cpp = PipelineHelperCpp(self.config, self.args, self.log, self.monomer_names_helper)

    @timing_decorator('Getting BGC variants')
    def get_bgc_variants(self) -> BGC_Variants_Info:
        bgc_variants = self.pipeline_helper_antismash.get_bgc_variants()
        if not bgc_variants:
            self.log.info("No BGC variants found. Exiting.")
            self.finish()
            exit(0)

        if not self.args.disable_bgc_deduplication:
            self.log.info("Deduplicating extracted BGC variants...")
            bgc_id_to_repr_id = cluster_isomorphic_bgc_variants(bgc_variants)
        else:
            bgc_id_to_repr_id = {bgc_variant.bgc_variant_id: bgc_variant.bgc_variant_id
                                 for bgc_variant in bgc_variants}
        return BGC_Variants_Info(bgc_variants=bgc_variants,
                                 bgc_id_to_repr_id=bgc_id_to_repr_id)

    @timing_decorator('Getting NRP variants')
    def get_nrp_variants_and_rban_records(self) -> NRP_Variants_Info:
        self.log.info("\n======= Getting rBAN output")
        if self.args.parsed_rban_records is not None:
            self.log.info(f"Loading preprocessed rBAN records from {self.args.parsed_rban_records}...")
            rban_records = load_parsed_rban_records(self.args.parsed_rban_records)
        else:
            rban_records = self.pipeline_helper_rban.get_rban_results()
        nrp_variants = self.pipeline_helper_rban.get_nrp_variants_info(rban_records)

        if not nrp_variants:
            self.log.info("No NRP variants found. Exiting.")
            self.finish()
            exit(0)

        if not self.args.disable_nrp_deduplication:
            self.log.info("Deduplicating extracted NRP variants...")
            nrp_id_to_repr_id = cluster_isomorphic_nrp_variants(nrp_variants)
        else:
            nrp_id_to_repr_id = {nrp_variant.nrp_variant_id: nrp_variant.nrp_variant_id
                                 for nrp_variant in nrp_variants}
        return NRP_Variants_Info(nrp_variants=nrp_variants,
                                 rban_records=rban_records,
                                 nrp_id_to_repr_id=nrp_id_to_repr_id)

    @timing_decorator('Constructing HMMs')
    def construct_hmms(self, bgc_variants: List[BGC_Variant]) -> List[DetailedHMM]:
        self.log.info("\n======= Constructing HMMs")
        hmms = []
        for bgc_variant in bgc_variants:
            try:
                hmms.append(DetailedHMM.from_bgc_variant(bgc_variant, self.hmm_helper))
            except Exception as e:
                if self.args.let_it_crash:
                    raise
                self.log.warning(f"Error constructing HMM for BGC variant "
                                 f"{bgc_variant.bgc_variant_id.bgc_id.antiSMASH_file}")

        return hmms

    @timing_decorator('Generating linearizations')
    def get_nrp_linearizations(self, nrp_variants: List[NRP_Variant]) \
            -> List[NRP_Linearizations]:
        self.log.info("\n======= Generating NRP linearizations")
        return get_all_nrp_linearizations(nrp_variants)

    @timing_decorator('Matching')
    def get_hmm_matches(self,
                        hmms: List[DetailedHMM],
                        nrp_linearizations: List[NRP_Linearizations]) -> List[HMM_Match]:
        #for i, hmm in enumerate(hmms):
        #    hmm.draw(Path(f'{hmm.bgc_variant.genome_id}.png'))
        self.log.info("\n======= Nerpa matching")
        if self.args.fast_matching:
            cpp_output = self.pipeline_helper_cpp.get_hmm_matches_and_p_values(hmms, nrp_linearizations)

            # Set p-value estimators for hmms (TODO: that's ugly, it breaks encapsulation)
            hmm_by_id = {hmm.bgc_variant.bgc_variant_id: hmm for hmm in hmms}
            assert hmm_by_id.keys() == cpp_output.p_values_by_bgc_variant.keys(), \
                "Mismatch in BGC variant IDs between HMMs and precomputed p-values"
            '''
            for bgc_variant_id, p_values in cpp_output.p_values_by_bgc_variant.items():
                hmm_by_id[bgc_variant_id]._p_value_estimator = PValueEstimator._from_precomputed_p_values(p_values)
            '''
            return cpp_output.matches
        else:
            return get_hmm_matches(hmms,
                                   nrp_linearizations,
                                   self.config.matching_config,
                                   self.args.threads,
                                   self.log)

    def get_matches(self,
                    hmms: List[DetailedHMM],
                    nrp_linearizations: List[NRP_Linearizations],
                    nrp_variants: List[NRP_Variant]) -> List[Match]:
        hmm_matches = self.get_hmm_matches(hmms, nrp_linearizations)
        if self.args.draw_hmms:
            self.log.info("======== Drawing HMMs with optimal paths")
            write_results.draw_hmms_with_optimal_paths(hmms, hmm_matches,
                                                       self.config.output_config.main_out_dir)
        self.log.info("\n======= Reconstructing alignments for matches")
        return convert_to_detailed_matches(hmms, nrp_variants, hmm_matches)

    @timing_decorator('Writing results')
    def write_results(self,
                      matches: List[Match],
                      bgc_variants_info: BGC_Variants_Info,
                      nrp_variants_info: NRP_Variants_Info,
                      write_only_what_is_matched: bool = True,
                      matches_details: bool = True):
        self.log.info("\n======= Writing results")
        write_results.write_results(matches,
                                    bgc_variants_info,
                                    nrp_variants_info,
                                    self.config.output_config,
                                    write_only_what_is_matched=write_only_what_is_matched,
                                    matches_details=matches_details,
                                    log=self.log)
        self.log.info("RESULTS:")
        self.log.info("Main report is saved to " + str(self.config.output_config.report), indent=1)
        self.log.info("HTML report is saved to " + str(self.config.output_config.html_report), indent=1)
        self.log.info("Detailed reports are saved to " + str(self.config.output_config.matches_details), indent=1)

    def finish(self):
        if not self.args.keep_intermediate_files:
            self.log.info("Removing intermediate files...")
            shutil.rmtree(self.config.output_config.rban_output_config.rban_output_dir,
                          ignore_errors=True)
            self.config.output_config.cpp_io_config.cpp_output_json.unlink(missing_ok=True)
            self.config.output_config.cpp_io_config.hmms_json.unlink(missing_ok=True)
            self.config.output_config.cpp_io_config.nrp_linearizations_json.unlink(missing_ok=True)

        self.log.finish()


