from typing import (
    Callable,
    Dict,
    Iterable,
    List,
    Optional,
    Union
)

from src.aa_specificity_prediction_model.specificity_prediction_helper import SpecificityPredictionHelper
from src.pipeline.command_line_args_helper import CommandLineArgs
from src.pipeline.logger import NerpaLogger
from src.pipeline.download_antismash_results import download_antismash_results
from src.config import Config
from src.data_types import AA34, BGC_Variant, LogProb

from src.pipeline.nerpa_utils import sys_call, get_path_to_program
from src.monomer_names_helper import MonomerNamesHelper, MonomerResidue
from src.antismash_parsing.antismash_parser import parse_antismash_json
from src.antismash_parsing.antismash_parser_types import antiSMASH_record
from src.antismash_parsing.build_bgc_variants import build_bgc_variants
from src.aa_specificity_prediction_model.model_wrapper import ModelWrapper
from src.write_results import write_bgc_variants
from pathlib import Path
import json
import yaml
from dataclasses import dataclass
from itertools import chain
from copy import deepcopy


@dataclass
class PipelineHelper_antiSMASH:
    config: Config
    args: CommandLineArgs
    monomer_names_helper: MonomerNamesHelper
    specificity_prediction_helper: SpecificityPredictionHelper
    log: NerpaLogger
    antismash_exec: Union[Path, None] = None

    def __post_init__(self):
        if self.args.antismash_path:
            self.antismash_exec = get_path_to_program('run_antismash.py', dirpath=self.args.antismash_path,
                                                      min_version='5.0')
        else:
            self.antismash_exec = get_path_to_program('antismash', min_version='5.0')

    def preprocessed_bgc_variants(self) -> bool:
        return self.args.bgc_variants is not None

    def load_bgc_variants(self) -> List[BGC_Variant]:
        self.log.info('Loading preprocessed BGC variants')
        bgc_variants: List[BGC_Variant] = []
        for file_with_bgc_variants in filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                             self.args.bgc_variants.iterdir()):
            bgc_variants.extend(BGC_Variant.from_yaml_dict(yaml_record)
                                for yaml_record in yaml.safe_load(file_with_bgc_variants.read_text()))

        return bgc_variants

    # TODO: watch out for name collisions
    def create_symlinks_to_antismash_results(self, antismash_dirs: Iterable[Path]):
        for antismash_dir in antismash_dirs:
            if antismash_dir.parent == self.config.output_config.antismash_out_dir:
                continue
            symlink = self.config.output_config.antismash_out_dir / antismash_dir.name
            if not symlink.exists():  # in case of a name collision keep only the first path (TODO: handle properly)
                symlink.symlink_to(antismash_dir.resolve())

    # TODO: refactor: this function is doing too much
    def get_antismash_results(self) -> Iterable[antiSMASH_record]:
        self.config.output_config.antismash_out_dir.mkdir(exist_ok=True)
        antismash_results: List[Path] = []
        if self.args.antismash is not None:
            antismash_results.extend(self.args.antismash)
        if self.args.antismash_outpaths_file is not None:
            antismash_results.extend(map(Path, self.args.antismash_outpaths_file.read_text().strip().splitlines()))

        if self.args.seqs:
            for seq in self.args.seqs:
                try:
                    new_results = self.run_antismash(seq,
                                                     self.args.threads,
                                                     self.config.output_config.antismash_out_dir,
                                                     self.log)
                except Exception as e:
                    self.log.error(f'Error while running antismash on {seq}, aborting')
                    raise e
                antismash_results.append(new_results)

        if self.args.antismash_job_ids is not None:
            antismash_results.extend(download_antismash_results(job_id,
                                                                self.config.output_config.antismash_out_dir,
                                                                self.log)
                                     for job_id in self.args.antismash_job_ids)
        antismash_jsons = [antismash_json_file
                           for antismash_dir in antismash_results
                           for antismash_json_file in antismash_dir.glob('**/*.json')]
        self.create_symlinks_to_antismash_results(antismash_json.parent
                                                  for antismash_json in antismash_jsons)
        return (antiSMASH_record(json.loads(antismash_json_file.read_text()))
                for antismash_dir in antismash_results
                for antismash_json_file in antismash_dir.glob('**/*.json'))

    def get_bgc_variants(self) -> List[BGC_Variant]:
        if self.preprocessed_bgc_variants():
            self.log.info('Loading preprocessed BGC variants. All other inputs will be ignored')
            return self.load_bgc_variants()

        antismash_results = self.get_antismash_results()

        self.log.info('\n======= Predicting BGC variants')

        bgc_variants = []
        for antismash_record in antismash_results:
            try:
                bgc_variants.extend(self.extract_bgc_variants_from_antismash(antismash_record))
            except Exception as e:
                self.log.info(f'Error while parsing antismash record: {antismash_record["input_file"]}, skipping')
                raise e  # TODO: remove

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
        sys_call(command, log, cwd=self.config.output_config.main_out_dir)
        return output_dir

    def extract_bgc_variants_from_antismash(self,
                                            antismash_json: antiSMASH_record) -> List[BGC_Variant]:
        antismash_bgcs = parse_antismash_json(antismash_json,
                                              self.config.antismash_processing_config,
                                              self.log)
        bgc_variants = list(chain.from_iterable(build_bgc_variants(bgc,
                                                                   self.specificity_prediction_helper,
                                                                   self.config.antismash_processing_config,
                                                                   self.log)
                                                for bgc in antismash_bgcs))

        if self.args.debug:
            debug_specificity_prediction_cfg = deepcopy(self.config.specificity_prediction_config)
            debug_specificity_prediction_cfg.ENABLE_CALIBRATION = False
            debug_specificity_prediction_cfg.ENABLE_DICTIONARY_LOOKUP = False

            debug_specificity_prediction_cfg.APRIORI_RESIDUE_PROB = \
                {res: 1.0 for res in self.monomer_names_helper.supported_residues}  # this effectively forces computing logprobs instead of log-odds
            debug_spec_helper = SpecificityPredictionHelper(debug_specificity_prediction_cfg,
                                                            self.monomer_names_helper,
                                                            self.log)

            bgc_variants_no_cal = list(chain.from_iterable(build_bgc_variants(bgc,
                                                                              debug_spec_helper,
                                                                              self.config.antismash_processing_config,
                                                                              self.log)
                                                    for bgc in antismash_bgcs))
            write_bgc_variants(bgc_variants_no_cal,
                               self.config.output_config.bgc_variants_no_calibration_dir)


        return bgc_variants


