from __future__ import annotations
from typing import (
    Iterable,
    List,
    Union, Tuple, NamedTuple, Dict
)

from more_itertools import chunked

from src.aa_specificity_prediction_model.specificity_prediction_helper import SpecificityPredictionHelper
from src.pipeline.logging.buffered_logger import BufferedLogger
from src.pipeline.command_line_args_helper import CommandLineArgs
from src.pipeline.logging.logger import NerpaLogger
from src.pipeline.download_antismash_results import download_antismash_results
from src.config import Config, antiSMASH_Processing_Config
from src.data_types import BGC_Variant

from src.pipeline.nerpa_utils import sys_call, get_path_to_program
from src.monomer_names_helper import MonomerNamesHelper
from src.antismash_parsing.antismash_parser import parse_antismash_json
from src.antismash_parsing.antismash_parser_types import antiSMASH_record, BGC_ID
from src.antismash_parsing.build_bgc_variants import build_bgc_variants
from src.pipeline.parse_antismash_parallel import extract_bgc_variants_from_antismash_batch
from pathlib import Path
import json
import yaml
from dataclasses import dataclass
from itertools import chain
from joblib import Parallel, delayed


class BGC_Variant_ID(NamedTuple):
    bgc_id: BGC_ID
    variant_idx: int

    @classmethod
    def from_dict(cls, data: dict) -> BGC_Variant_ID:
        return cls(
            bgc_id=BGC_ID.from_dict(data["bgc_id"]),
            variant_idx=data["variant_idx"]
        )

    def to_dict(self) -> dict:
        return {'bgc_id': self.bgc_id.to_dict(),
                'variant_idx': self.variant_idx}

class BGC_Variants_Info(NamedTuple):
    bgc_variants: List[BGC_Variant]
    bgc_id_to_repr_id: Dict[BGC_Variant_ID, BGC_Variant_ID]

    def get_representative_bgcs(self) -> List[BGC_Variant]:
        repr_ids = set(self.bgc_id_to_repr_id.values())
        return [bgc for bgc in self.bgc_variants if bgc.bgc_variant_id in repr_ids]


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

        if self.args.bgc_variants.is_dir():
            files_with_bgc_variants = filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                             self.args.bgc_variants.iterdir())
        else:
            files_with_bgc_variants = [self.args.bgc_variants]

        return [BGC_Variant.from_yaml_dict(yaml_record)
                for file_with_bgc_variants in files_with_bgc_variants
                for yaml_record in yaml.safe_load(file_with_bgc_variants.read_text())]


    # TODO: watch out for name collisions
    def create_symlinks_to_antismash_results(self, antismash_dirs: Iterable[Path]):
        for antismash_dir in antismash_dirs:
            if antismash_dir.parent == self.config.output_config.antismash_out_dir:
                continue
            symlink = self.config.output_config.antismash_out_dir / antismash_dir.name
            if not symlink.exists():  # in case of a name collision keep only the first path (TODO: handle properly)
                symlink.symlink_to(antismash_dir.resolve())

    # TODO: refactor: this function is doing too much
    # I am returning an iterator of antiSMASH_record to avoid loading all JSONs into memory at once
    def get_antismash_results(self) -> List[Path]:
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
        antismash_jsons: List[Path] = []
        for antismash_dir in antismash_results:
            if antismash_dir.is_file() and antismash_dir.suffix == '.json':
                antismash_jsons.append(antismash_dir)
                continue

            antismash_jsons_in_dir = list(antismash_dir.glob('**/*.json'))
            if not antismash_jsons_in_dir:
                self.log.warning(f'No antiSMASH json files found in {antismash_dir}, skipping')
            antismash_jsons.extend(antismash_jsons_in_dir)

        self.create_symlinks_to_antismash_results(antismash_json.parent
                                                  for antismash_json in antismash_jsons)
        return antismash_jsons

    def get_bgc_variants(self) -> List[BGC_Variant]:
        if self.preprocessed_bgc_variants():
            self.log.info('Loading preprocessed BGC variants. All other inputs will be ignored')
            return self.load_bgc_variants()

        as_json_paths = self.get_antismash_results()
        if not as_json_paths:
            self.log.error('No antiSMASH results provided, aborting')
            raise ValueError('No antiSMASH results provided')

        self.log.info(f'\n======= Predicting BGC variants from antiSMASH results'
                     f' using {self.args.threads} threads')

        bgc_variants_in_chunks = Parallel(n_jobs=self.args.threads,
                                          prefer='processes')(
            delayed(extract_bgc_variants_from_antismash_batch)(
                antismash_paths_chunk,
                self.args)
            for antismash_paths_chunk in chunked(as_json_paths,
                                                 len(as_json_paths) // self.args.threads + 1)
        )
        bgc_variants_batches, loggers = zip(*bgc_variants_in_chunks)
        BufferedLogger.replay(self.log, loggers)


        self.log.info('\n======= Done with Predicting BGC variants')
        return list(chain(*bgc_variants_batches))

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





