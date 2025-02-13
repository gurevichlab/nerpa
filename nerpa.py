#!/usr/bin/env python3

import sys
from src.pipeline.pipeline_helper import PipelineHelper
from src.pipeline.logger import NerpaLogger
from src.write_results import write_bgc_variants
from pathlib import Path


def main(log: NerpaLogger):  # log is passed as an argument to make it easier to write log in case of exception
    pipeline_helper = PipelineHelper(log)

    bgc_variants = pipeline_helper.get_bgc_variants()
    #write_bgc_variants(bgc_variants, Path('/home/ilianolhin/git/nerpa2/input'))
    hmms = pipeline_helper.construct_hmms(bgc_variants)

    nrp_variants, rban_records = pipeline_helper.get_nrp_variants_and_rban_records()
    nrp_linearizations = pipeline_helper.get_nrp_linearizations(nrp_variants)

    matches = pipeline_helper.get_matches(hmms, nrp_linearizations, nrp_variants)
    pipeline_helper.write_results(matches, bgc_variants, nrp_variants, rban_records)


if __name__ == "__main__":
    log = NerpaLogger()
    try:
        main(log)
    except Exception as e:
        _, exc_value, _ = sys.exc_info()
        log.exception(exc_value)
    finally:
        # TODO: clean up: remove all intermediate files
        pass
