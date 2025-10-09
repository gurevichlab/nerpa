#!/usr/bin/env python3

import sys

from src.pipeline.pipeline_helper import PipelineHelper
from src.pipeline.logging.logger import PreliminaryLogger
from src.build_output.write_results import write_nrp_variants
import traceback


def main(pre_logger: PreliminaryLogger):  # log is passed as an argument to make it easier to write log in case of exception
    pipeline_helper = PipelineHelper(pre_logger)

    nrp_variants_info = pipeline_helper.get_nrp_variants_and_rban_records()
    #representative_nrps = nrp_variants_info.get_representative_nrp_variants()
    representative_nrps = nrp_variants_info.nrp_variants
    nrp_linearizations = pipeline_helper.get_nrp_linearizations(representative_nrps)

    bgc_variants_info = pipeline_helper.get_bgc_variants()
    #representative_bgcs = bgc_variants_info.get_representative_bgc_variants()
    representative_bgcs = bgc_variants_info.bgc_variants
    #compute_modification_freqs(bgc_variants)
    hmms = pipeline_helper.construct_hmms(representative_bgcs)

    matches = pipeline_helper.get_matches(hmms, nrp_linearizations, representative_nrps)
    pipeline_helper.write_results(matches, bgc_variants_info, nrp_variants_info)

    pipeline_helper.finish()


if __name__ == "__main__":
    log = PreliminaryLogger()
    try:
        main(log)
    except Exception as e:
        log.exception(traceback.format_exc())
        exit(1)

'''
def compute_modification_freqs(bgc_variants: List[BGC_Variant]):
    processed_bgc_ids = set()

    total_modules = 0
    total_methylated = 0
    total_epimerized = 0
    for bgc_variant in bgc_variants:
        if bgc_variant.bgc_variant_id.bgc_id in processed_bgc_ids:
            continue
        processed_bgc_ids.add(bgc_variant.bgc_variant_id.bgc_id)

        for module in bgc_variant.modules:
            total_modules += 1
            if BGC_Module_Modification.METHYLATION in module.modifications:
                total_methylated += 1
            if BGC_Module_Modification.EPIMERIZATION in module.modifications:
                total_epimerized += 1

    methylation_freq = total_methylated / total_modules
    epimerization_freq = total_epimerized / total_modules
    print(f"Total modules: {total_modules}")
    print(f"Methylation frequency: {methylation_freq:.4f}")
    print(f"Epimerization frequency: {epimerization_freq:.4f}")
'''