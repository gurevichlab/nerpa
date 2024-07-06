import sys
from src.pipeline.pipeline_helper import PipelineHelper
from src.pipeline.logger import NerpaLogger


def main(log: NerpaLogger):  # log is passed as an argument to make it easier to write log in case of exception
    pipeline_helper = PipelineHelper(log)

    bgc_variants = pipeline_helper.pipeline_helper_antismash.get_bgc_variants()
    nrp_variants, rban_records = pipeline_helper.get_nrp_variants_and_rban_records()

    matches = pipeline_helper.get_matches(bgc_variants, nrp_variants)
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
