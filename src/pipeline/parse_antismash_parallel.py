from typing import (
    Iterable,
    List,
    Tuple
)

from src.aa_specificity_prediction_model.specificity_prediction_helper import SpecificityPredictionHelper
from src.pipeline.logging.buffered_logger import BufferedLogger
from src.pipeline.command_line_args_helper import CommandLineArgs
from src.config import antiSMASH_Processing_Config, load_monomer_names_helper
from src.data_types import BGC_Variant

from src.antismash_parsing.antismash_parser import parse_antismash_json
from src.antismash_parsing.build_bgc_variants import build_bgc_variants
from src.config import load_config
from src.pipeline.logging.logging_types import AnyLogger
from src.pipeline.paras_parsing import get_paras_results_all
from pathlib import Path
import json
from itertools import chain

# I load configs in each worker to avoid issues with multiprocessing and large objects pickling
def _load_configs(args: CommandLineArgs) -> Tuple[antiSMASH_Processing_Config, SpecificityPredictionHelper]:
    config = load_config(args)

    monomer_names_helper = load_monomer_names_helper(config.monomers_config,
                                                     config.nerpa_dir)

    external_specificity_predictions = get_paras_results_all(args.paras_results,
                                                             monomer_names_helper,
                                                             log=AnyLogger()) \
        if args.paras_results is not None else None
    specificity_prediction_helper = SpecificityPredictionHelper(config.specificity_prediction_config,
                                                                monomer_names_helper,
                                                                external_specificity_predictions)

    return config.antismash_processing_config, specificity_prediction_helper


def extract_bgc_variants_from_antismash_batch(antismash_paths: Iterable[Path],
                                              args: CommandLineArgs)\
        -> Tuple[List[BGC_Variant], BufferedLogger]:
    antismash_processing_config, specificity_prediction_helper = _load_configs(args)
    log = BufferedLogger()
    bgc_variants: List[BGC_Variant] = []
    for antismash_json_file in antismash_paths:
        try:
            antismash_bgcs = parse_antismash_json(antismash_json_file,
                                                  antismash_processing_config,
                                                  log)
            new_bgc_variants = chain.from_iterable(build_bgc_variants(bgc,
                                                                      specificity_prediction_helper,
                                                                      antismash_processing_config,
                                                                      log)
                                                   for bgc in antismash_bgcs)
            bgc_variants.extend(new_bgc_variants)
        except Exception as e:
            log.error(f'Unexpected error while parsing antiSMASH JSON {antismash_json_file}: {e}'
                      f'\nSkipping this file.')

    return bgc_variants, log
