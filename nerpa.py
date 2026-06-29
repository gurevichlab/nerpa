#!/usr/bin/env python3

import sys

from src.pipeline.pipeline_helper import PipelineHelper
from src.pipeline.logging.logger import PreliminaryLogger
from src.build_output.write_results import write_nrp_variants
import traceback


def build_data_for_lars(pipeline_helper: PipelineHelper):
    import polars as pl
    import json
    from pathlib import Path
    df: pl.DataFrame = pl.read_csv(
        pipeline_helper.config.specificity_prediction_config.a_domains_signatures,
        separator='\t',
        has_header=True,
    )

    # Add PARAS_predictions by running foo on the aa34 sequence for each row.
    df = df.with_columns(
        pl.col('aa34')
        .map_elements(lambda aa34: pipeline_helper.specificity_prediction_helper.paras_model_wrapper.predict(aa34),
                      return_dtype=pl.Object)
        .alias('PARAS_predictions')
    )

    rows: list[dict[str, Any]] = df.to_dicts()

    # Write as a JSON list of dicts.
    output_json = Path('paras_predictions_known_specificities.json')
    output_json.write_text(json.dumps(rows, indent=2), encoding='utf-8')


def main(pre_logger: PreliminaryLogger):  # log is passed as an argument to make it easier to write log in case of exception
    pipeline_helper = PipelineHelper(pre_logger)

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
    pipeline_helper.write_results(matches, bgc_variants_info, nrp_variants_info,
                                  write_only_what_is_matched=not pipeline_helper.args.dump_all_preprocessed)

    pipeline_helper.finish()


if __name__ == "__main__":
    log = PreliminaryLogger()
    try:
        main(log)
    except Exception as e:
        log.exception(traceback.format_exc())
        exit(1)
