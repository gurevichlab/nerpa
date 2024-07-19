from typing import (
    Callable,
    Dict,
    Iterable,
    List,
    Union
)
from src.pipeline.command_line_args_helper import CommandLineArgs
from src.pipeline.logger import NerpaLogger
from src.config import Config, SpecificityPredictionConfig
from src.aa_specificity_prediction_model.model_wrapper import ModelWrapper
from src.data_types import BGC_Variant

from src.pipeline.nerpa_utils import sys_call, get_path_to_program
from src.rban_parsing.rban_names_helper import rBAN_Names_Helper
from src.antismash_parsing.antismash_parser import parse_antismash_json
from src.antismash_parsing.build_bgc_variants import build_bgc_variants
from src.write_results import write_bgc_variants
from pathlib import Path
import json
import yaml
from dataclasses import dataclass
from itertools import chain
import math


def create_step_function(steps: List[float]) -> Callable[[float], float]:
    step_len = 1 / len(steps)
    def step_function(x: float) -> float:
        for i, step in enumerate(steps):
            if x < (i + 1) * step_len:
                return steps[i]
        return steps[-1]
    return step_function


def calibrate_scores(predictions: Dict[str, float],
                     config: SpecificityPredictionConfig):
    predictions = {res: math.e ** score for res, score in predictions.items()}
    if config.apply_step_function:
        step_function = create_step_function(config.calibration_step_function_steps)
        predictions = {res: step_function(score)
                       for res, score in predictions.items()}
    if config.normalize_scores:
        total_prob = sum(predictions.values())
        if total_prob < 1:
            predictions = {res: score + (1 - total_prob) * config.apriori_residue_prob[res]
                           for res, score in predictions.items()}
        else:
            predictions = {res: score / total_prob
                           for res, score in predictions.items()}

    if config.pseudo_counts:
        predictions = {res: score * (1 - config.pseudo_count_fraction) +
                            config.pseudo_count_fraction * config.apriori_residue_prob[res]
                       for res, score in predictions.items()}

    if config.compute_evidence:
        predictions = {res: score / config.apriori_residue_prob[res]
                       for res, score in predictions.items()}

    predictions = {res: math.log(score) for res, score in predictions.items()}
    return predictions


@dataclass
class PipelineHelper_antiSMASH:
    config: Config
    args: CommandLineArgs
    monomer_names_helper: rBAN_Names_Helper
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
        sys_call(command, log, cwd=self.config.paths.main_out_dir)
        return output_dir

    def extract_bgc_variants_from_antismash(self,
                                            antismash_json: dict) -> List[BGC_Variant]:
        antismash_bgcs = parse_antismash_json(antismash_json,
                                              self.monomer_names_helper,
                                              self.config.antismash_parsing_config)
        specificity_prediction_model = ModelWrapper(self.config.paths.specificity_prediction_model)
        bgc_variants = list(chain.from_iterable(build_bgc_variants(bgc,
                                                                   self.log,
                                                                   specificity_prediction_model,
                                                                   self.config.antismash_parsing_config)
                                                for bgc in antismash_bgcs))

        if self.config.specificity_prediction_config.calibration:
            (self.config.paths.main_out_dir / 'BGC_variants_before_calibration').mkdir()
            write_bgc_variants(bgc_variants,
                               self.config.paths.main_out_dir / 'BGC_variants_before_calibration')

            for bgc_variant in bgc_variants:
                for module in bgc_variant.tentative_assembly_line:
                    module.residue_score = calibrate_scores(module.residue_score,
                                                            self.config.specificity_prediction_config)
        return bgc_variants


