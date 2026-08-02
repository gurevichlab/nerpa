#!/usr/bin/env python3

import sys
from typing import Optional

from src.pipeline.pipeline_helper import PipelineHelper
from src.pipeline.logging.logger import PreliminaryLogger
from src.build_output.write_results import write_nrp_variants
import traceback
import argparse


def nerpa(
        pre_logger: PreliminaryLogger,
        args: Optional[argparse.Namespace] = None
):
    pipeline_helper = PipelineHelper(pre_logger=pre_logger, args=args)

    bgc_variants_info = pipeline_helper.get_bgc_variants()
    representative_bgcs = bgc_variants_info.get_representative_bgc_variants()
    hmms = pipeline_helper.construct_hmms(representative_bgcs)

    nrp_variants_info = pipeline_helper.get_nrp_variants_and_rban_records()
    representative_nrps = nrp_variants_info.get_representative_nrp_variants()
    nrp_linearizations = pipeline_helper.get_nrp_linearizations(representative_nrps)

    matches = pipeline_helper.get_matches(
        hmms,
        nrp_linearizations,
        representative_nrps,
        nrp_variants_info.rban_records
    )
    pipeline_helper.write_results(
        matches,
        bgc_variants_info,
        nrp_variants_info,
        write_only_what_is_matched=not pipeline_helper.args.dump_all_preprocessed
    )

    pipeline_helper.finish()

def main(args: Optional[argparse.Namespace] = None):  # log is passed as an argument to make it easier to write log in case of exception
    pre_logger = PreliminaryLogger()
    try:
        nerpa(pre_logger=pre_logger, args=args)
    except Exception as e:
        pre_logger.exception(traceback.format_exc())
        exit(1)


if __name__ == "__main__":
    main()
