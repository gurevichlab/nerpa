from typing import List
from command_line_args_helper import CommandLineArgs
from src.nerpa_pipeline.logger import NerpaLogger
from config import Config
from src.data_types import NRP_Variant

from src.nerpa_pipeline import handle_rban
from src.nerpa_pipeline.rban_parser import (
    Parsed_rBAN_Record,
    parsed_chiralities,
)
from src.nerpa_pipeline.handle_monomers import get_monomers_chirality
from src.nerpa_pipeline.nerpa_utils import sys_call, get_path_to_program
from src.nerpa_pipeline.rban_names_helper import rBAN_Names_Helper

from rban_helper import rBAN_Helper
from pathlib import Path
import csv
import json
import yaml
from dataclasses import dataclass
from itertools import chain


@dataclass
class PipelineHelper_antiSMASH:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    antismash_helper: antiSMASH_Helper = None

    def __post_init__(self):
        self.set_antismash_helper()

    def set_antismash_helper(self):
        if self.args.antismash_path:
            antismash_exec = get_path_to_program('run_antismash.py', dirpath=self.args.antismash_path,
                                                 min_version='5.0')
        else:
            antismash_exec = get_path_to_program('antismash', min_version='5.0')
        self.antismash_helper = antiSMASH_Helper(antismash_exec, self.config.paths.antismash_out_dir)

    def preprocessed_bgc_variants(self) -> bool:
        return self.args.predictions is not None

    def load_bgc_variants(self) -> List[BGC_Variant]:
        self.log.info('Loading preprocessed BGC variants')
        bgc_variants = []
        # q: iterate over yaml files in self.args.predictions directory
        for path_to_predictions in self.args.predictions:
            for file_with_bgc_variants in filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                                 Path(path_to_predictions).iterdir()):
                bgc_variants.extend(BGC_Variant.from_yaml_dict(yaml_record)
                                    for yaml_record in yaml.safe_load(file_with_bgc_variants.read_text()))

        return bgc_variants

    def get_bgc_variants(self) -> List[BGC_Variant]:
        if self.preprocessed_bgc_variants():
            return self.load_bgc_variants()

        antismash_results = self.args.antismash if self.args.antismash is not None else []
        if self.args.seqs:
            new_results = self.antismash_helper.run_antismash(self.args.seqs,
                                                              self.args.threads,
                                                              self.config.paths.antismash_out_dir)
            if new_results is not None:
                antismash_results.append(new_results)


        self.log.info('\n======= Predicting BGC variants')
        bgc_variants = self.antismash_helper.get_bgc_variants()
        self.log.info('\n======= Done with Predicting BGC variants')
        return bgc_variants

    def run_antismash(self):
        antismash_results_dirs = self.args.antismash if self.args.antismash is not None else []
        if self.args.seqs:
            cur_antismash_out = os.path.join(output_dir, 'antismash_output')
            if args.antismash_path:
                antismash_exe = nerpa_utils.get_path_to_program('run_antismash.py', dirpath=args.antismash_path, min_version='5.0')
            else:
                antismash_exe = nerpa_utils.get_path_to_program('antismash', min_version='5.0')
            if antismash_exe is None:
                log.error("Can't find antismash 5.x executable. Please make sure that you have antismash 5.x installed "
                          "in your system or provide path to antismash source directory via --antismash-path option.")
            command = [antismash_exe,
                       '--genefinding-tool', 'prodigal',
                       '--output-dir', cur_antismash_out,
                       '--minimal', '--skip-zip-file', '--enable-nrps-pks',
                       '--cpus', str(args.threads), args.seqs]
            nerpa_utils.sys_call(command, log, cwd=output_dir)
            antismash_out_dirs.append(cur_antismash_out)

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
