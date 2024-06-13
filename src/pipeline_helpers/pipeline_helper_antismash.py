from typing import (
    Iterable,
    List,
    Union
)
from command_line_args_helper import CommandLineArgs
from src.nerpa_pipeline.logger import NerpaLogger
from src.nerpa_pipeline.nerpa_utils import sys_call
from config import Config
from src.data_types import BGC_Variant

from src.nerpa_pipeline import handle_rban
from src.nerpa_pipeline.rban_parser import (
    Parsed_rBAN_Record,
    parsed_chiralities,
)
from src.nerpa_pipeline.handle_monomers import get_monomers_chirality
from src.nerpa_pipeline.nerpa_utils import sys_call, get_path_to_program
from src.nerpa_pipeline.rban_names_helper import rBAN_Names_Helper
from src.antismash_parsing.antismash_parser import parse_antismash_json
from src.antismash_parsing.bgcs_split_and_reorder import split_and_reorder
from src.antismash_parsing.build_bgc_variants import build_bgc_variants
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
    antismash_exec: Union[Path, None] = None

    def __post_init__(self):
        if self.args.antismash_path:
            self.antismash_exec = get_path_to_program('run_antismash.py', dirpath=self.args.antismash_path,
                                                 min_version='5.0')
        else:
            self.antismash_exec = get_path_to_program('antismash', min_version='5.0')

    def preprocessed_bgc_variants(self) -> bool:
        return self.args.predictions is not None

    def load_bgc_variants(self) -> List[BGC_Variant]:
        self.log.info('Loading preprocessed BGC variants')
        bgc_variants = []
        for path_to_predictions in self.args.predictions:
            for file_with_bgc_variants in filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                                 Path(path_to_predictions).iterdir()):
                bgc_variants.extend(BGC_Variant.from_yaml_dict(yaml_record)
                                    for yaml_record in yaml.safe_load(file_with_bgc_variants.read_text()))

        return bgc_variants

    def get_antismash_results(self) -> Iterable[dict]:
        antismash_results = self.args.antismash if self.args.antismash is not None else []

        if self.args.seqs:
            try:
                new_results = self.run_antismash(self.args.seqs,
                                                 self.args.threads,
                                                 self.config.paths.antismash_out_dir,
                                                 self.log)
            except Exception as e:
                self.log.error(f'Error while running antismash on {self.args.seqs}, aborting')
                raise e
            antismash_results.append(new_results)

        return (json.loads(antismash_json_file.read_text())
                for antismash_dir in antismash_results
                for antismash_json_file in antismash_dir.glob('**/*.json'))

    def get_bgc_variants(self) -> List[BGC_Variant]:
        if self.preprocessed_bgc_variants():
            return self.load_bgc_variants()

        antismash_results = self.get_antismash_results()

        self.log.info('\n======= Predicting BGC variants')

        bgc_variants = []
        for antismash_record in antismash_results:
            try:
                bgc_variants.extend(self.extract_bgc_variants_from_antismash(antismash_record))
            except:
                self.log.error(f'Error while parsing antismash record: {antismash_record}, skipping')

        self.log.info('\n======= Done with Predicting BGC variants')
        return bgc_variants

    def run_antismash(self, dna_sequences: Path,
                      threads: int, output_dir: Path,
                      log: NerpaLogger) -> Path:
        command = [str(self.antismash_exec),
                   '--genefinding-tool', 'prodigal',
                   '--output-dir', str(output_dir),
                   '--minimal', '--skip-zip-file', '--enable-nrps-pks',
                   '--cpus', str(threads),
                   str(dna_sequences)]
        sys_call(command, log, cwd=self.config.paths.main_out_dir)
        return output_dir

    def extract_bgc_variants_from_antismash(self, antismash_json: dict) -> List[BGC_Variant]:
        antismash_bgcs = parse_antismash_json(antismash_json, self.config.antismash_parsing_config)
        bgcs = split_and_reorder(antismash_bgcs)
        return chain.from_iterable(build_bgc_variants(bgc, self.config.antismash_parsing_config) for bgc in bgcs)


