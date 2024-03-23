#!/usr/bin/env python3

import sys
import os
import argparse
import csv
import site
import shutil
import json
import yaml

import nerpa_init
from pipeline_helper import PipelineHelper
from run_antismash import get_antismash_v3_compatible_input_paths

site.addsitedir(nerpa_init.python_modules_dir)

import src.nerpa_pipeline.predictions_preprocessor as predictions_preprocessor
import src.nerpa_pipeline.nerpa_utils as nerpa_utils
import src.nerpa_pipeline.handle_rban
from src.nerpa_pipeline.logger import NerpaLogger

from src.nerpa_pipeline.rban_names_helper import rBAN_Names_Helper
from pathlib import Path

# for detecting and processing antiSMASH v.5 output
site.addsitedir(os.path.join(nerpa_init.python_modules_dir, 'NRPSPredictor_utils'))
from src.nerpa_pipeline.NRPSPredictor_utils.json_handler import get_main_json_fpath
from src.nerpa_pipeline.NRPSPredictor_utils.main import main as convert_antiSMASH_v5
from src.NewMatcher.scoring_helper import ScoringHelper
from src.NewMatcher.scoring_config import load_scoring_config as load_scoring_config
from src.NewMatcher.matcher import get_matches

from src.data_types import (
    BGC_Variant,
    NRP_Variant,
    UNKNOWN_RESIDUE
)
from command_line_args_helper import CommandLineArgs

from src.write_results import write_results, write_nrp_variants, write_bgc_variants







def main(log: NerpaLogger):
    pipeline_helper = PipelineHelper(log)
    args = pipeline_helper.args
    output_dir = pipeline_helper.config.paths.main_out_dir

    if args.predictions is not None:
            bgc_variants = []
            for path_to_predictions in args.predictions:
                for file_with_bgc_variants in filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                                     Path(path_to_predictions).iterdir()):
                    bgc_variants.extend(BGC_Variant.from_yaml_dict(yaml_record)
                                        for yaml_record in yaml.safe_load(file_with_bgc_variants.read_text()))
    else:
        antismash_out_dirs = args.antismash if args.antismash is not None else []
        if args.seqs:
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

        bgc_variants = predictions_preprocessor.parse_antismash_output(get_antismash_v3_compatible_input_paths(
                listing_fpath=args.antismash_out, list_of_paths=antismash_out_dirs,
                output_dir=output_dir, log=log), output_dir, args.debug, log)


    if pipeline_helper.preprocessed_nrp_variants():
        rban_records = []
        nrp_variants = pipeline_helper.load_nrp_variants()
    else:
        rban_records = pipeline_helper.get_rban_results()
        nrp_variants = pipeline_helper.get_nrp_variants(rban_records)

    matches = pipeline_helper.get_matches(bgc_variants, nrp_variants)
    pipeline_helper.write_results(matches, bgc_variants, nrp_variants, rban_records)


if __name__ == "__main__":
    log = NerpaLogger()
    try:
        main(log)
    except Exception:
        _, exc_value, _ = sys.exc_info()
        log.exception(exc_value)
    finally:
        # TODO: clean up: remove all intermediate files
        pass
