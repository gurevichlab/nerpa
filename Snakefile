from __future__ import annotations

from pathlib import Path
from typing import Any

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
ANTISMASH_RESULTS = ext_path('as_results_mibig4_nrps')

# ===== Important intermediate files (not to rewrite paths in the rules)
PNRPDB2 = NERPA_DIR / 'data' / 'input' / 'pnrpdb2_filtered.tsv'
PNRPDB2_MIBIG_NORINE = NERPA_DIR / 'data' / 'input' / 'pnrpdb2_mibig_norine_deduplicated.tsv'
MIBIG_NORINE_PREPROCESSED = NERPA_DIR / 'data' / 'input' / 'preprocessed' / 'pnrpdb2_mibig_norine_parsed_rban_records.yaml'

# ==== Important directories
TEST_RESULTS_DIR = NERPA_DIR / 'test_results'

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
        wrong_matches=params.test_results_dir / 'wrong_matches.yaml',
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
        script=NERPA_DIR / 'scripts' / 'preprocess_pnrpdb2.py',
        nerpa_output_dir=NERPA_DIR / 'nerpa_results' / 'preprocessed_pnrpdb2',
    input:
        pnrpdb2=PNRPDB2_INITIAL,
    output:
        # "filtered" means ID for which Nerpa is able build NRP variants
        pnrpdb2_filtered=params.tables_dir / 'pnrpdb2_filtered.tsv',
        pnrpdb2_preprocessed=params.preprocessed_dir / 'pnrpdb2_parsed_rban_records.yaml',
        pnrpdb2_deduplicated=PNRPDB2,
        pnrpdb2_deduplicated_preprocessed=params.preprocessed_dir / 'pnrpdb2_deduplicated_parsed_rban_records.yaml',
        pnrpdb2_mibig_norine=PNRPDB2_MIBIG_NORINE,
        pnrpdb2_mibig_norine_preprocessed=params.preprocessed_dir / 'pnrpdb2_mibig_norine_parsed_rban_records.yaml',
        pnrpdb2_mibig_norine_deduplicated=params.tables_dir / 'pnrpdb2_mibig_norine_deduplicated.tsv',
        pnrpdb2_mibig_norine_deduplicated_preprocessed=PNRPDB2_MIBIG_NORINE_PREPROCESSED,
    shell:
        r"""
        python {params.script}  \
          --nrp-database {input.pnrpdb2} \
          --filtered-tables-dir {params.tables_dir} \
          --preprocessed-dir {params.preprocessed_dir} \
          run-nerpa \
          --output-for-nerpa-run {params.nerpa_output_dir} \
          --deduplicate
        """
