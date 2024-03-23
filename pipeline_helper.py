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
import src.write_results as report

from rban_helper import rBAN_Helper
from pathlib import Path
import csv
import json
import shutil
import yaml
from dataclasses import dataclass
from itertools import chain


class PipelineHelper:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    rban_helper: rBAN_Helper
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

        self.set_rban_helper()
        scoring_config = load_scoring_config(self.config.paths.scoring_config)  # TODO: this is ugly
        self.scoring_helper = ScoringHelper(scoring_config)

    def set_rban_helper(self):
        custom_monomers = list(chain(*(json.load(path_to_monomers.open('r'))
                                       for path_to_monomers in (self.config.paths.nerpa_monomers,
                                                                self.args.rban_monomers)
                                       if path_to_monomers is not None)))
        self.rban_helper = rBAN_Helper(self.config.rban_config, custom_monomers)

    def get_input_for_rban(self) -> List[dict]:
        def default_id(i: int) -> str:
            return f'compound_{i:06d}'

        if self.args.smiles_tsv:
            reader = csv.DictReader(self.args.smiles_tsv.open('r'),
                                    delimiter=self.args.sep, quoting=csv.QUOTE_NONE)
            return [{'id': row[self.args.col_id] if self.args.col_id else default_id(i),
                     'smiles': row[self.args.col_smiles]}
                    for i, row in enumerate(reader)]
        elif self.args.smiles:
            return [{'id': default_id(i), 'smiles': smiles}
                    for i, smiles in enumerate(self.args.smiles)]

    def get_rban_results(self) -> List[Parsed_rBAN_Record]:
        smiles_with_ids = self.get_input_for_rban()
        self.log.info('\n======= Structures preprocessing with rBAN')
        rban_records = self.rban_helper.run_rban(smiles_with_ids,
                                                 self.log, report_not_processed=True)
        hybrid_monomers_dict = self.rban_helper.get_hybrid_monomers(rban_records, self.log)
        chiralities_dict = {}
        for i, rban_record in enumerate(rban_records):
            try:
                chiralities_dict[i] = parsed_chiralities(get_monomers_chirality(rban_record))
            except:
                self.log.warning(
                    f'Structure "{rban_record["id"]}": unexpected error while determining stereoisomeric configuration '
                    f'of monomers. Stereoisomeric configuration will be omitted.')

        self.log.info("\n======= Done with Structures preprocessing with rBAN")
        return [Parsed_rBAN_Record(rban_record,
                                   hybrid_monomers_dict[i],
                                   chiralities_dict[i])
                for i, rban_record in enumerate(rban_records)]

    def get_nrp_variants(self, parsed_rban_records: List[Parsed_rBAN_Record]) -> List[NRP_Variant]:
        self.log.info('\n======= Processing rBAN output')
        self.log.info(f'results will be in {self.config.paths.main_out_dir / Path("NRP_Variants")}')
        recognized_monomers = [row['Code']
                               for row in csv.DictReader(self.config.paths.nerpa_monomers_info.open('r'),
                                                         delimiter='\t')]
        rban_names_helper = rBAN_Names_Helper(self.config.paths.nerpa_monomers_info)
        nrp_variants = handle_rban.retrieve_nrp_variants(parsed_rban_records,
                                                         recognized_monomers, rban_names_helper,
                                                         self.log)
        self.log.info('\n======= Done with Processing rBAN output')
        return nrp_variants

    def preprocessed_nrp_variants(self) -> bool:
        return self.args.structures is not None

    def load_nrp_variants(self) -> List[NRP_Variant]:
        self.log.info('Loading preprocessed NRP variants')
        nrp_variants = []
        for file_with_nrp_variants in filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                             Path(self.args.structures).iterdir()):
            nrp_variants.extend(NRP_Variant.from_yaml_dict(yaml_record)
                                for yaml_record in yaml.safe_load(file_with_nrp_variants.read_text()))
        return nrp_variants

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
