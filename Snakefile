from __future__ import annotations

from pathlib import Path
from typing import Any
import yaml

from snakemake.exceptions import WorkflowError

NERPA_ROOT = Path(workflow.basedir)
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
APPROVED_MATCHES = NERPA_ROOT / 'data/for_training_and_testing/approved_matches.yaml'
PNRPDB2_INITIAL = NERPA_ROOT / 'data' / 'input' / 'pnrpdb2_raw.tsv'
ANTISMASH_RESULTS = ext_path('antiSMASH results on MIBiG 4')
NERPA1_RESULTS_ON_MIBIG = ext_path('NERPA1 results on MIBiG 4 vs MIBIG+NORINE NRPs')
BIOCAT_RESULTS_ON_MIBIG = ext_path('BioCAT results on MIBiG 4 vs MIBIG+NORINE NRPs')

# ===== Important intermediate files (not to rewrite paths in the rules)
# PNRPDB2 without the compounds Nerpa is unable to process
PNRPDB2 = NERPA_ROOT / 'data' / 'input' / 'pnrpdb2.tsv'
# Deduplicated compounds from MIBiG and Norine, aka "confirmed NRPs"
PNRPDB2_MIBIG_NORINE = NERPA_ROOT / 'data' / 'input' / 'pnrpdb2_mibig_norine_deduplicated.tsv'
MIBIG_NORINE_PREPROCESSED = NERPA_ROOT / 'data' / 'input' / 'preprocessed' / 'pnrpdb2_mibig_norine_deduplicated_parsed_rban_records.yaml'
# Clustering and stats on PNRPDB2 compounds
PNRPDB2_INFO = NERPA_ROOT / 'data' / 'for_training_and_testing' / 'pnrpdb2_info.tsv'
NERPA2_RESULTS_ON_MIBIG_NO_CALIBRATION = NERPA_ROOT / 'nerpa_results' / 'mibig_no_calibration'
NERPA2_RESULTS_ON_MIBIG_VS_MIBIG_NORINE = NERPA_ROOT / 'nerpa_results' / 'mibig_vs_mibig_norine' 
# ===== Important directories
TEST_RESULTS_DIR = NERPA_ROOT / 'test_results'
TRAINING_RESULTS_DIR = NERPA_ROOT / 'training_results'
BENCHMARKING_PLOTS_DIR = NERPA_ROOT / 'benchmarking_results'
SMILES_DIR = NERPA_ROOT / 'data' / 'input'
NRP_PREPROCESSED_DIR = NERPA_ROOT / 'data' / 'input' / 'preprocessed'

rule all:
    input:
        TEST_RESULTS_DIR / '.done',
        TRAINING_RESULTS_DIR / '.done',
        BENCHMARKING_PLOTS_DIR / 'mibig4_wo_training_bgcs' / '.done',
        BENCHMARKING_PLOTS_DIR / 'training_bgcs' / '.done',


rule test_nerpa:
    params:
        test_results_dir=TEST_RESULTS_DIR,
    input:
        matches=APPROVED_MATCHES,
        mibig_norine_preprocessed=MIBIG_NORINE_PREPROCESSED,
        pnrpdb2_info=PNRPDB2_INFO,
        antismash=ANTISMASH_RESULTS,
    output:
        test_nerpa_done=TEST_RESULTS_DIR / '.done',
    log:
        TEST_RESULTS_DIR / 'test_nerpa.log',
    shell:
        r"""
        mkdir -p test_results
        python test_nerpa.py \
          --test-matches {input.matches} \
          --preprocessed-nrps {input.mibig_norine_preprocessed} \
          --nrps-info {input.pnrpdb2_info} \
          --antismash-results {input.antismash} \
          --output-dir {params.test_results_dir} \
          > {log} 2>&1
        touch {output.test_nerpa_done}
        """

rule preprocess_pnrpdb2:
    params:
        tables_dir=NERPA_ROOT / 'data' / 'input',
        preprocessed_dir=NERPA_ROOT / 'data' / 'input' / 'preprocessed',
        preprocess_script=NERPA_ROOT / 'scripts' / 'preprocess_pnrpdb2.py',
        nerpa_output_dir=NERPA_ROOT / 'nerpa_results' / 'preprocessed_pnrpdb2',
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
    shell:
        r"""
        PYTHONPATH={NERPA_ROOT} python {params.preprocess_script}  \
          --nrp-database {input.pnrpdb2} \
          --filtered-tables-dir {params.tables_dir} \
          --preprocessed-dir {params.preprocessed_dir} \
          run-nerpa \
          --output-for-nerpa-run {params.nerpa_output_dir} \
          --deduplicate
        """

rule compute_pnrpdb2_info:
    params:
        build_pnrpdb2_info_script=NERPA_ROOT / 'scripts' / 'pnrpdb_info.py',
    input:
        pnrpdb2_preprocessed=NRP_PREPROCESSED_DIR / 'pnrpdb2_parsed_rban_records.yaml',
    output:
        pnrpdb2_info=PNRPDB2_INFO,
    shell:
        r"""
        PYTHONPATH={NERPA_ROOT} python {params.build_pnrpdb2_info_script} \
            --pnrpdb2-preprocessed {input.pnrpdb2_preprocessed} \
            --out {output.pnrpdb2_info}
        """

rule preprocess_mibig_no_calibration:
    params:
        nerpa_exec=NERPA_ROOT / 'nerpa.py',
        output_dir=NERPA2_RESULTS_ON_MIBIG_NO_CALIBRATION,
    input:
        ANTISMASH_RESULTS
    output:
        NERPA2_RESULTS_ON_MIBIG_NO_CALIBRATION / 'preprocessed_input' / 'BGC_variants.yaml',
    shell:
        r"""
        python {params.nerpa_exec} \
            --antismash {input} \
            # a random SMILES. It doesn't matter -- we just want to get prepocessed BGCs
            --smiles CC(=O)OC1=CC=CC=C1C(=O)O \
            --output-dir {params.output_dir} \
            --force-output-dir \
            --min-num-matches-per-bgc 1 \
            # 0 means unlimited
            --max-num-matches 0 \
            --disable-calibration \
            --dump-all-preprocessed \
            --let-it-crash
        """
        
    
rule train_nerpa:
    params:
        train_script=NERPA_ROOT / 'scripts' / 'train_nerpa.py',
        nerpa_results_on_mibig_no_calibration=NERPA2_RESULTS_ON_MIBIG_NO_CALIBRATION,
        output_dir=TRAINING_RESULTS_DIR,
    input:
        approved_matches=APPROVED_MATCHES,
        # ensure we have preprocessed BGC variants available for training
        mibig_no_calib_variants=NERPA2_RESULTS_ON_MIBIG_NO_CALIBRATION / 'preprocessed_input' / 'BGC_variants.yaml',
        pnrpdb2_info=PNRPDB2_INFO,
        mibig_norine_preprocessed=MIBIG_NORINE_PREPROCESSED,
    output:
        train_nerpa_done=NERPA_ROOT / 'training_results' / '.done',
    log:
        NERPA_ROOT / 'training_results' / 'train_nerpa.log',
    shell:
        r"""
        mkdir -p {output.outdir}
        PYTHONPATH={NERPA_ROOT} python {params.train_script} \
          --approved-matches {input.approved_matches} \
          --nerpa-results-on-mibig-no-calibration {params.nerpa_results_on_mibig_no_calibration} \
          --pnrpdb2-info {input.pnrpdb2_info} \
          --mibig-norine-preprocessed {input.mibig_norine_preprocessed} \
          --output-dir {output.outdir} \
          > {output.log} 2>&1
        touch {output.train_nerpa_done}
        """
   
rule nerpa2_on_mibig_vs_mibig_norine:
    params:
        nerpa_exec=NERPA_ROOT / 'nerpa.py',
        outdir=NERPA2_RESULTS_ON_MIBIG_VS_MIBIG_NORINE,
    input:
        antismash=ANTISMASH_RESULTS,
        # confirmed NRPs (MIBiG+Norine), already parsed by rBAN
        parsed_rban_records=MIBIG_NORINE_PREPROCESSED,
    output:
        NERPA2_RESULTS_ON_MIBIG_VS_MIBIG_NORINE / 'report.tsv',
    shell:
        r"""
        python {params.nerpa_exec} \
          --antismash {input.antismash} \
          --parsed-rban-records {input.parsed_rban_records} \
          --output-dir {params.outdir} \
          --force-output-dir \
          --min-num-matches-per-bgc 10 \
          --max-num-matches-per-bgc 10 \
          --min-num-matches-per-nrp 0 \
          --max-num-matches-per-nrp 0 \
          --max-num-matches 0 \
          --disable-bgc-deduplication \
          --process-hybrids \
          --threads 8 \
          --skip-molecule-drawing \
          --dump-all-preprocessed
        """

rule benchmarking_plots:
    params:
        benchmarking_script=NERPA_ROOT / 'benchmarking_plots.py',
        outdir=BENCHMARKING_PLOTS_DIR / '{bgc_test_set}',
    input:
        nerpa1_report=NERPA1_RESULTS_ON_MIBIG,
        nerpa2_report=NERPA2_RESULTS_ON_MIBIG_VS_MIBIG_NORINE / 'report.tsv',
        biocat_report=BIOCAT_RESULTS_ON_MIBIG,
    output:
        benchmarking_done=BENCHMARKING_PLOTS_DIR / '{bgc_test_set}' / '.done',
    log:
        BENCHMARKING_PLOTS_DIR / '{bgc_test_set}' / 'benchmarking_plots.log',
    shell:
        r"""
         mkdir -p {params.outdir}
         PYTHONPATH={NERPA_ROOT} python {params.benchmarking_script} \
           --nerpa1-report {input.nerpa1_report} \
           --nerpa2-report {input.nerpa2_report} \
           --biocat-report {input.biocat_report} \
           --output-dir {params.outdir} \
           --bgc-test-set {wildcards.bgc_test_set} \
           > {log} 2>&1
         touch {output.benchmarking_done}
        """
