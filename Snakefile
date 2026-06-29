from __future__ import annotations

from pathlib import Path
from typing import Any
import yaml

from snakemake.exceptions import WorkflowError

NERPA_DIR = Path(workflow.basedir)
PATHS_TO_EXTERNAL_DATA = "paths_to_external_data.yaml"


def ext_path(key: str) -> Path:
    config = yaml.safe_load(open(PATHS_TO_EXTERNAL_DATA))
    val: Any = config.get(key)
    if val in (None, ''):
        raise WorkflowError(
            f'Missing external data path: config["{key}"] is not set.\n'
            f'Edit {PATHS_TO_EXTERNAL_DATA} and set "{key}: /path/to/data".'
        )
    return Path(str(val)).expanduser()

# ===== Source data files -- these are not produced by any rule
APPROVED_MATCHES = NERPA_DIR / 'data/for_training_and_testing/approved_matches.yaml'
PNRPDB2_INITIAL = NERPA_DIR / 'data' / 'input' / 'pnrpdb2.tsv'
ANTISMASH_RESULTS = ext_path('antiSMASH results on MIBiG 4')

# ===== Important intermediate files (not to rewrite paths in the rules)
# PNRPDB2 without the compounds Nerpa is unable to process
PNRPDB2 = NERPA_DIR / 'data' / 'input' / 'pnrpdb2_filtered.tsv'
# Deduplicated compounds from MIBiG and Norine, aka "confirmed NRPs"
PNRPDB2_MIBIG_NORINE = NERPA_DIR / 'data' / 'input' / 'pnrpdb2_mibig_norine_deduplicated.tsv'
MIBIG_NORINE_PREPROCESSED = NERPA_DIR / 'data' / 'input' / 'preprocessed' / 'pnrpdb2_mibig_norine_deduplicated_parsed_rban_records.yaml'
# Clustering and stats on PNRPDB2 compounds
PNRPDB2_INFO = NERPA_DIR / 'data' / 'for_training_and_testing' / 'pnrpdb2_info.tsv'

# ==== Important directories
TEST_RESULTS_DIR = NERPA_DIR / 'test_results'
SMILES_DIR = NERPA_DIR / 'data' / 'input'
NRP_PREPROCESSED_DIR = NERPA_DIR / 'data' / 'input' / 'preprocessed'

rule all:
    input:
        TEST_RESULTS_DIR / 'test_nerpa.log',


rule test_nerpa:
    params:
        test_results_dir=TEST_RESULTS_DIR,
    input:
        matches=APPROVED_MATCHES,
        mibig_norine_preprocessed=MIBIG_NORINE_PREPROCESSED,
        antismash=ANTISMASH_RESULTS,
    output:
        log=NERPA_DIR / 'test_results' / 'test_nerpa.log',
        wrong_matches=TEST_RESULTS_DIR / 'wrong_matches.yaml',
    shell:
        r"""
        mkdir -p test_results
        python test_nerpa.py \
          --test-matches {input.matches} \
          --parsed-rban-records {input.mibig_norine_preprocessed} \
          --antismash-results {input.antismash} \
          --output-dir {params.test_results_dir} \
          > {output.log} 2>&1
        """

rule preprocess_pnrpdb2:
    params:
        tables_dir=NERPA_DIR / 'data' / 'input',
        preprocessed_dir=NERPA_DIR / 'data' / 'input' / 'preprocessed',
        preprocess_script=NERPA_DIR / 'scripts' / 'preprocess_pnrpdb2.py',
        build_pnrpdb2_info_script=NERPA_DIR / 'scripts' / 'pnrpdb_info.py',
        nerpa_output_dir=NERPA_DIR / 'nerpa_results' / 'preprocessed_pnrpdb2',
    input:
        pnrpdb2=PNRPDB2_INITIAL,
    output:
        # "filtered" means ID for which Nerpa is able build NRP variants
        pnrpdb2_filtered=PNRPDB2,
        pnrpdb2_preprocessed=NRP_PREPROCESSED_DIR / 'pnrpdb2_parsed_rban_records.yaml',
        pnrpdb2_deduplicated=SMILES_DIR / 'pnrpdb2_deduplicated.tsv',
        pnrpdb2_deduplicated_preprocessed=NRP_PREPROCESSED_DIR / 'pnrpdb2_deduplicated_parsed_rban_records.yaml',
        pnrpdb2_mibig_norine=SMILES_DIR / 'pnrpdb2_mibig_norine.tsv',
        pnrpdb2_mibig_norine_preprocessed=NRP_PREPROCESSED_DIR / 'pnrpdb2_mibig_norine_parsed_rban_records.yaml',
        pnrpdb2_mibig_norine_deduplicated=PNRPDB2_MIBIG_NORINE,
        pnrpdb2_mibig_norine_deduplicated_preprocessed=MIBIG_NORINE_PREPROCESSED,
        pnrpdb2_info=PNRPDB2_INFO,
    shell:
        r"""
        python {params.preprocess_script}  \
          --nrp-database {input.pnrpdb2} \
          --filtered-tables-dir {params.tables_dir} \
          --preprocessed-dir {params.preprocessed_dir} \
          run-nerpa \
          --output-for-nerpa-run {params.nerpa_output_dir} \
          --deduplicate

        python {params.build_pnrpdb2_info_script} \
            --pnrpdb2-preprocessed {output.pnrpdb2_preprocessed} \
            --out {output.pnrpdb2_info}
        """

# used for training
NERPA_RESULTS_ON_MIBIG_NO_CALIBRATION = NERPA_DIR / 'nerpa_results' / 'mibig_no_calibration'

rule preprocess_mibig_no_calibration:
    params:
        nerpa_exec=NERPA_DIR / 'nerpa.py',
    input:
        ANTISMASH_RESULTS
    output:
        NERPA_RESULTS_ON_MIBIG_NO_CALIBRATION / 'preprocessed_input' / 'BGC_variants.yaml',
    shell:
        r"""
        python {params.nerpa_exec} \
            --antismash-results {input} \
            # a random SMILES. It doesn't matter -- we just want to get prepocessed BGCs
            --smiles CC(=O)OC1=CC=CC=C1C(=O)O \
            --output-dir {NERPA_RESULTS_ON_MIBIG_NO_CALIBRATION} \
            --force-output-dir \
            --min-num-matches-per-bgc 1 \
            # 0 means unlimited
            --max-num-matches 0 \
            --disable-calibration \
            --dump-all-preprocessed \
            --let-it-crash
        """
        
        
    
