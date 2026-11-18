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
            f'Edit/Create {PATHS_TO_EXTERNAL_DATA} and set "{key}: /path/to/data".'
        )
    return Path(str(val)).expanduser()

# ===== Source data files -- these are not produced by any rule
MONOMER_NAMES_TABLE = NERPA_ROOT / 'data' / 'nerpa_internal' / 'monomer_names_table.tsv'
CURATED_ALIGNMENTS = NERPA_ROOT / 'data/for_training_and_testing/curated_alignments.json'
PNRPDB2_RAW = NERPA_ROOT / 'data' / 'input' / 'pnrpdb2_raw.tsv'
NORINE_GRAPHS = NERPA_ROOT / 'data/for_training_and_testing/norine_graphs.csv'
ANTISMASH_RESULTS = ext_path('antiSMASH results on MIBiG 4')
BIGSCAPE_RESULTS_ON_MIBIG = ext_path('BiG-SCAPE clustering of MIBiG 4')
MIBIG3_JSONS = ext_path('MIBiG 3 JSONs')

NERPA1_RESULTS_ON_MIBIG = NERPA_ROOT / 'data' / 'for_training_and_testing' / 'nerpa1_report_mibig4_vs_mibig_norine.csv'
BIOCAT_RESULTS_ON_MIBIG = ext_path('BioCAT results on MIBiG 4 vs MIBiG+NORINE NRPs')

# ===== Important intermediate files (not to rewrite paths in the rules)
# information about the monomers in the Norine database
NORINE_INFO = NERPA_ROOT / 'data/nerpa_internal/norine_monomers_info.tsv'
MIBIG_INFO = NERPA_ROOT / 'data' / 'for_training_and_testing' / 'mibig_bgcs_info.tsv'
# PNRPDB2 without the compounds Nerpa is unable to process
PNRPDB2 = NERPA_ROOT / 'data' / 'input' / 'pnrpdb2.tsv'
# PNRPDB2 with alternative IDs expanded, so that each ID has its own row, and each compound has a unique ID
PNRPDB2_EXPANDED = NERPA_ROOT / 'data' / 'input' / 'pnrpdb2_expanded.tsv'
# Deduplicated compounds from MIBiG and Norine, aka "confirmed NRPs"
PNRPDB2_MIBIG_NORINE = NERPA_ROOT / 'data' / 'input' / 'pnrpdb2_mibig_norine_deduplicated.tsv'
MIBIG_NORINE_PREPROCESSED = NERPA_ROOT / 'data' / 'input' / 'preprocessed' / 'pnrpdb2_mibig_norine_deduplicated_parsed_rban_records.yaml'
MIBIG_NORINE_PREPROCESSED_NO_DEDUP = NERPA_ROOT / 'data' / 'input' / 'preprocessed' / 'pnrpdb2_mibig_norine_parsed_rban_records.yaml'

# Clustering and stats on PNRPDB2 compounds
PNRPDB2_INFO = NERPA_ROOT / 'data' / 'for_training_and_testing' / 'pnrpdb2_info.tsv'
PNRPDB2_COMPOUND_SIMILARITY = (
    NERPA_ROOT / 'data' / 'for_training_and_testing' / 'pnrpdb2_compound_similarity.tsv'
)
NERPA2_RESULTS_ON_MIBIG_NO_CALIBRATION = NERPA_ROOT / 'nerpa_results' / 'mibig_no_calibration'
NERPA2_RESULTS_ON_MIBIG_VS_MIBIG_NORINE = NERPA_ROOT / 'nerpa_results' / 'mibig_vs_mibig_norine' 

# ===== BGC and NRP sets for training and testing
MIN_NUM_A_DOMAINS = 3
MIN_NUM_NRPS_MONOMERS = 3

TRAINING_BGCs_AS_RESULTS_PATHS = NERPA_ROOT / 'data' / 'for_training_and_testing' / 'training_bgcs_as_results_paths.txt'
TRAINING_NRPs = NERPA_ROOT / 'data' / 'for_training_and_testing' / 'training_nrps.tsv'
TEST_BGCs_AS_RESULTS_PATHS = NERPA_ROOT / 'data' / 'for_training_and_testing' / 'test_bgcs_as_results_paths.txt'
TEST_NRPs = NERPA_ROOT / 'data' / 'for_training_and_testing' / 'test_nrps.tsv'

NERPA2_RESULTS_ON_TRAINING_BGCs_vs_their_NRPs = (
    NERPA_ROOT
    / 'nerpa_results'
    / f'training_bgcs_vs_their_nrps'
)
NERPA2_RESULTS_ON_TRAINING_BGCs_vs_MIBIG_NORINE = (
    NERPA_ROOT
    / 'nerpa_results'
    / f'training_bgcs_vs_mibig_norine'
)

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
        PNRPDB2_COMPOUND_SIMILARITY,
        NERPA2_RESULTS_ON_TRAINING_BGCs_vs_their_NRPs / '.done',

rule release:
    input:
        version=NERPA_ROOT / "VERSION.txt",
        make_release_script=(
            NERPA_ROOT
            / "scripts"
            / "make_release"
            / "make_release.sh"
        ),
        readme_template=NERPA_ROOT / "docs" / "README_template.md",
        render_readme_script=(
            NERPA_ROOT
            / "scripts"
            / "make_release"
            / "render_readme.py"
        ),
    params:
        generated_readme=NERPA_ROOT / "tmp" / "release" / "README.md",
    shell:
        r"""
        mapfile -t commit_tags < <(git tag --points-at HEAD)
        release_tag_args=()

        case "${{#commit_tags[@]}}" in
          0)
            ;;
          1)
            release_tag_args=(--release-tag "${{commit_tags[0]}}")
            ;;
          *)
            echo "Error: multiple tags are attached to HEAD:" >&2
            printf '  %s\n' "${{commit_tags[@]}}" >&2
            exit 1
            ;;
        esac

        python {input.render_readme_script:q} \
          --template {input.readme_template:q} \
          --version-file {input.version:q} \
          --output {params.generated_readme:q} \
          "${{release_tag_args[@]}}"

        bash {input.make_release_script:q} \
          {input.version:q} \
          {params.generated_readme:q}
        """

rule create_training_bgcs_nrps:
    params:
        script=(
            NERPA_ROOT
            / "scripts"
            / "build_training_and_testing_sets"
            / "get_training_bgcs_nrps.py"
        ),
    input:
        curated_alignments=CURATED_ALIGNMENTS,
        pnrpdb_expanded=PNRPDB2_EXPANDED,
        antismash_results=ANTISMASH_RESULTS,
    output:
        bgc_paths=TRAINING_BGCs_AS_RESULTS_PATHS,
        nrps=TRAINING_NRPs
    shell:
        r"""
        PYTHONPATH={NERPA_ROOT:q} python {params.script:q} \
          --curated-alignments {input.curated_alignments:q} \
          --antismash-results {input.antismash_results:q} \
          --pnrpdb-expanded {input.pnrpdb_expanded:q} \
          --output-bgcs {output.bgc_paths:q} \
          --output-nrps {output.nrps:q}
        """

rule create_test_bgcs_as_results_paths:
    params:
        script=(
            NERPA_ROOT
            / "scripts"
            / "build_training_and_testing_sets"
            / "get_testing_bgcs.py"
        ),
    input:
        curated_alignments=CURATED_ALIGNMENTS,
        antismash_results=ANTISMASH_RESULTS,
        pnrpdb_info=PNRPDB2_INFO,
        mibig_info=MIBIG_INFO,
    output:
        paths=TEST_BGCs_AS_RESULTS_PATHS,
    shell:
        r"""
        PYTHONPATH={NERPA_ROOT:q} python {params.script:q} \
          --curated-alignments {input.curated_alignments:q} \
          --antismash-results {input.antismash_results:q} \
          --pnrpdb-info {input.pnrpdb_info:q} \
          --mibig-bgcs-info {input.mibig_info:q} \
          --min-a-domains {MIN_NUM_A_DOMAINS} \
          --min-nrps-monomers {MIN_NUM_NRPS_MONOMERS} \
          --output {output.paths:q}
        """

rule create_test_nrps:
    params:
        script=(
            NERPA_ROOT
            / "scripts"
            / "build_training_and_testing_sets"
            / "get_testing_nrps.py"
        ),
    input:
        pnrpdb_expanded=PNRPDB2_RAW,
        pnrpdb_info=PNRPDB2_INFO,
    output:
        nrps=TEST_NRPs,
    shell:
        r"""
        PYTHONPATH={NERPA_ROOT:q} python {params.script:q} \
          --pnrpdb-expanded {input.pnrpdb_expanded:q} \
          --pnrpdb-info {input.pnrpdb_info:q} \
          --min-nrps-monomers {MIN_NUM_NRPS_MONOMERS} \
          --output {output.nrps:q}
        """

rule run_nerpa_on_training_bgcs_vs_their_nrps:
    params:
        script=NERPA_ROOT / "scripts/nerpa_runs/training_bgcs_vs_their_nrps.py",
        output_dir=NERPA2_RESULTS_ON_TRAINING_BGCs_vs_their_NRPs,
    input:
        antismash=TRAINING_BGCs_AS_RESULTS_PATHS,
        nrps=TRAINING_NRPs,
        mibig_info=MIBIG_INFO,
    output:
        done=NERPA2_RESULTS_ON_TRAINING_BGCs_vs_their_NRPs / '.done',
    threads: 8
    log:
        NERPA2_RESULTS_ON_TRAINING_BGCs_vs_their_NRPs / 'run.log',
    shell:
        r"""
        mkdir -p {params.output_dir:q}
        PYTHONPATH={NERPA_ROOT:q} python {params.script:q} \
          --output-dir {params.output_dir:q} \
          --num-threads {threads} \
          --antismash-results {input.antismash:q} \
          --nrps {input.nrps:q} \
          --mibig-bgcs-info {input.mibig_info:q} \
          > {log:q} 2>&1
        touch {output.done:q}
        """

rule nerpa2_on_training_bgs_vs_mibig_norine:
    params:
        nerpa_exec=NERPA_ROOT / 'nerpa.py',
        outdir=NERPA2_RESULTS_ON_TRAINING_BGCs_vs_MIBIG_NORINE,
    input:
        antismash=TRAINING_BGCs_AS_RESULTS_PATHS,
        # confirmed NRPs (MIBiG+Norine), already parsed by rBAN
        parsed_rban_records=MIBIG_NORINE_PREPROCESSED,
    output:
         NERPA2_RESULTS_ON_TRAINING_BGCs_vs_MIBIG_NORINE / 'report.tsv',
    shell:
        r"""
        python {params.nerpa_exec} \
          --antismash-paths-file {input.antismash} \
          --parsed-rban-records {input.parsed_rban_records} \
          --output-dir {params.outdir} \
          --force-output-dir \
          --min-num-matches-per-bgc 50 \
          --max-num-matches-per-bgc 50 \
          --min-num-matches-per-nrp 0 \
          --max-num-matches-per-nrp 0 \
          --max-num-matches 0 \
          --disable-bgc-deduplication \
          --process-hybrids \
          --threads 8 \
          --dump-all-preprocessed
        """

rule test_nerpa:
    params:
        output_dir=TEST_RESULTS_DIR,
        nerpa_results_individual=NERPA2_RESULTS_ON_TRAINING_BGCs_vs_their_NRPs,
        nerpa_results_combined=NERPA2_RESULTS_ON_TRAINING_BGCs_vs_MIBIG_NORINE,
    input:
        curated_alignments=CURATED_ALIGNMENTS,
        nerpa_results_individual_done=NERPA2_RESULTS_ON_TRAINING_BGCs_vs_their_NRPs / ".done",
        nerpa_results_combined_report=NERPA2_RESULTS_ON_TRAINING_BGCs_vs_MIBIG_NORINE / "report.tsv",
        pnrpdb_info=PNRPDB2_INFO,
    output:
        test_nerpa_done=TEST_RESULTS_DIR / ".done",
    log:
        TEST_RESULTS_DIR / "test_nerpa.log",
    shell:
        r"""
        mkdir -p {params.output_dir:q}
        python test_nerpa.py \
          --curated-alignments {input.curated_alignments:q} \
          --nerpa-results-individual {params.nerpa_results_individual:q} \
          --nerpa-results-combined {params.nerpa_results_combined:q} \
          --pnrpdb-info {input.pnrpdb_info:q} \
          --output-dir {params.output_dir:q} \
          > {log:q} 2>&1
        touch {output.test_nerpa_done:q}
        """

# creates several tables by filtering PNRPDB2_RAW,
# and corresponding preprocessed rBAN records for each of the tables
rule preprocess_pnrpdb2:
    params:
        tables_dir=NERPA_ROOT / 'data' / 'input',
        preprocessed_dir=NERPA_ROOT / 'data' / 'input' / 'preprocessed',
        preprocess_script=(
            NERPA_ROOT
            / 'scripts'
            / 'data_preprocessing'
            / 'preprocess_pnrpdb2.py'
        ),
        nerpa_output_dir=NERPA_ROOT / 'nerpa_results' / 'preprocessed_pnrpdb2',

        # some antiSMASH results just to give Nerpa something to work with, so that it can produce preprocessed rBAN records
        # Any would work
        antismash_results_example=NERPA_ROOT / 'test_data' / 'nerpa' / 'antismash',
    input:
        pnrpdb2_raw=PNRPDB2_RAW,
    output:
        # no filtering
        NRP_PREPROCESSED_DIR / f'{PNRPDB2_RAW.stem}_parsed_rban_records.yaml',

        # expanded alternative IDs, so that each ID has its own row, and each compound has a unique ID
        PNRPDB2_EXPANDED,
        NRP_PREPROCESSED_DIR / f'{PNRPDB2_EXPANDED.stem}_parsed_rban_records.yaml',

        # filter out compounds for which Nerpa is unable build NRP variants
        PNRPDB2,
        NRP_PREPROCESSED_DIR / f'{PNRPDB2.stem}_parsed_rban_records.yaml',

        # only one compound per isomorphic NRP_Variant group (deduplicated)
        SMILES_DIR / 'pnrpdb2_deduplicated.tsv',
        NRP_PREPROCESSED_DIR / 'pnrpdb2_deduplicated_parsed_rban_records.yaml',

        # only compounds that are present in MIBiG and/or Norine (confirmed NRPs)
        mibig_norine=SMILES_DIR / 'pnrpdb2_mibig_norine.tsv',
        mibig_norine_readme=SMILES_DIR / "pnrpdb2_mibig_norine.tsv.readme",
        MIBIG_NORINE_PREPROCESSED_NO_DEDUP,

        # copies of non-deduplicated MIBiG+Norine table and its README with user-friendly names, for convenience
        friendly_mibig_norine=SMILES_DIR / "mibig_norine.tsv",
        friendly_mibig_norine_readme=SMILES_DIR / "mibig_norine.tsv.readme",

        # only one compound per isomorphic NRP_Variant group (deduplicated) among the confirmed NRPs
        PNRPDB2_MIBIG_NORINE,
        MIBIG_NORINE_PREPROCESSED,
    shell:
        r"""
        PYTHONPATH={NERPA_ROOT} python {params.preprocess_script}  \
          --pnrpdb2-raw {input.pnrpdb2_raw} \
          --antismash-results {params.antismash_results_example} \
          --filtered-tables-dir {params.tables_dir} \
          --preprocessed-dir {params.preprocessed_dir} \
          --output-for-nerpa-run {params.nerpa_output_dir}

       cp {output.mibig_norine:q} {output.friendly_mibig_norine:q}
       cp {output.mibig_norine_readme:q} {output.friendly_mibig_norine_readme:q}
        """

rule compute_pnrpdb2_info:
    params:
        build_pnrpdb2_info_script=(
            NERPA_ROOT
            / 'scripts'
            / 'data_preprocessing'
            / 'pnrpdb_info.py'
        ),
    input:
        pnrpdb2_expanded=PNRPDB2_EXPANDED,
        pnrpdb2_preprocessed=NRP_PREPROCESSED_DIR / 'pnrpdb2_expanded_parsed_rban_records.yaml',
    output:
        pnrpdb2_info=PNRPDB2_INFO,
    shell:
        r"""
        PYTHONPATH={NERPA_ROOT} python {params.build_pnrpdb2_info_script} \
            --pnrpdb2-preprocessed {input.pnrpdb2_preprocessed} \
            --pnrpdb2-expanded {input.pnrpdb2_expanded} \
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
            --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" \
            --output-dir {params.output_dir} \
            --force-output-dir \
            --min-num-matches-per-bgc 1 \
            --max-num-matches 0 \
            --disable-calibration \
            --dump-all-preprocessed \
            --let-it-crash
        """
        
    
rule train_nerpa:
    params:
        train_script=NERPA_ROOT / 'train_nerpa.py',
        nerpa_results_on_mibig_no_calibration=NERPA2_RESULTS_ON_MIBIG_NO_CALIBRATION,
        output_dir=TRAINING_RESULTS_DIR,
    input:
        curated_alignments=CURATED_ALIGNMENTS,
        # ensure we have preprocessed BGC variants available for training
        mibig_no_calib_variants=NERPA2_RESULTS_ON_MIBIG_NO_CALIBRATION / 'preprocessed_input' / 'BGC_variants.yaml',
        pnrpdb2_info=PNRPDB2_INFO,
        mibig_norine_preprocessed=MIBIG_NORINE_PREPROCESSED,
        norine_monomers_info=NORINE_INFO,
    output:
        train_nerpa_done=NERPA_ROOT / 'training_results' / '.done',
    shell:
        r"""
        PYTHONPATH={NERPA_ROOT} python {params.train_script} \
          --curated-alignments {input.curated_alignments} \
          --nerpa-results-on-mibig-no-calibration {params.nerpa_results_on_mibig_no_calibration} \
          --pnrpdb2-info {input.pnrpdb2_info} \
          --mibig-norine-preprocessed {input.mibig_norine_preprocessed} \
          --norine-monomers-info {input.norine_monomers_info} \
          --output-dir {params.output_dir}

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
        outdir=lambda wildcards: BENCHMARKING_PLOTS_DIR / wildcards.bgc_test_set,
    input:
        curated_alignments=CURATED_ALIGNMENTS,
        nerpa1_report=NERPA1_RESULTS_ON_MIBIG,
        nerpa2_report=NERPA2_RESULTS_ON_MIBIG_VS_MIBIG_NORINE / 'report.tsv',
        biocat_report=BIOCAT_RESULTS_ON_MIBIG,
        pnrpdb2_info=PNRPDB2_INFO,
        pnrpdb2_compound_similarity=PNRPDB2_COMPOUND_SIMILARITY,
        mibig_info=MIBIG_INFO,
    output:
        benchmarking_done=BENCHMARKING_PLOTS_DIR / '{bgc_test_set}' / '.done',
    log:
        BENCHMARKING_PLOTS_DIR / '{bgc_test_set}' / 'benchmarking_plots.log',
    shell:
        r"""
         mkdir -p {params.outdir}
         MPLBACKEND=Agg PYTHONPATH={NERPA_ROOT} python {params.benchmarking_script} \
           --nerpa1-report {input.nerpa1_report} \
           --nerpa2-report {input.nerpa2_report} \
           --biocat-report {input.biocat_report} \
           --mibig-bgcs-info {input.mibig_info} \
           --pnrpdb2-info {input.pnrpdb2_info} \
           --pnrpdb2-compound-similarity {input.pnrpdb2_compound_similarity} \
           --curated-alignments {input.curated_alignments} \
           --output-dir {params.outdir} \
           --bgc-test-set {wildcards.bgc_test_set} \
           > {log} 2>&1
         touch {output.benchmarking_done}
        """
        
rule compute_compound_similarity:
    params:
        compound_similarity_script=(
            NERPA_ROOT
            / 'scripts'
            / 'data_preprocessing'
            / 'pnrpdb_compound_similarity.py'
        ),
        nerpa_root=NERPA_ROOT,
    input:
        # Use the preprocessed PNRPDB2 rBAN records
        preprocessed_nrps=MIBIG_NORINE_PREPROCESSED_NO_DEDUP,
    output:
        out=PNRPDB2_COMPOUND_SIMILARITY,
    shell:
        r"""
        PYTHONPATH={params.nerpa_root} python {params.compound_similarity_script} \
            --preprocessed-nrps {input.preprocessed_nrps} \
            --out {output.out}
        """

rule compute_mibig_info:
    params:
        mibig_info_script=(
            NERPA_ROOT
            / 'scripts'
            / 'data_preprocessing'
            / 'mibig_bgcs_info.py'
        ),
        nerpa_root=NERPA_ROOT,
    input:
        mibig3_jsons=MIBIG3_JSONS,
        bigscape_results_on_mibig4=BIGSCAPE_RESULTS_ON_MIBIG,
        pnrpdb2_info=PNRPDB2_EXPANDED,
        preprocessed_bgcs=NERPA2_RESULTS_ON_MIBIG_NO_CALIBRATION / 'preprocessed_input' / 'BGC_variants.yaml',
    output:
        mibig_info=MIBIG_INFO,
    shell:
        r"""
        PYTHONPATH={params.nerpa_root} python {params.mibig_info_script} \
            --bigscape-results-dir {input.bigscape_results_on_mibig4} \
            --preprocessed-bgcs {input.preprocessed_bgcs} \
            --nrps-info {input.pnrpdb2_info} \
            --mibig3-jsons {input.mibig3_jsons} \
            --out {output.mibig_info}
        """
        
rule compute_norine_monomers_info:
    params:
        script=NERPA_ROOT / 'scripts' / 'data_preprocessing' / 'compute_norine_monomers_info.py',
    input:
        names_table=MONOMER_NAMES_TABLE,
        norine_graphs=NORINE_GRAPHS,
    output:
        norine_info=NORINE_INFO,
    shell:
        r"""
        PYTHONPATH={NERPA_ROOT} python {params.script} \
            --norine-graphs {input.norine_graphs} \
            --monomer-names-table {input.names_table} \
            --out {output.norine_info}
        """
    
