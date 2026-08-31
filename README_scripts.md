# Nerpa 2.1 — Training / Testing / Benchmarking README

This directory contains a Snakemake workflow for **training**, **testing**, and **benchmarking** of Nerpa.

## Requirements

- **Nerpa 2.1** installed (the main Nerpa repository)
- **Conda environment for Nerpa** created and activated (e.g., `nerpa-env`)
- **Snakemake** available in the active environment (or otherwise installed and runnable as `snakemake`)
- Linux recommended (tested on typical Linux setups)

## Installation / Setup

1. **Unpack the archive into the Nerpa root directory** (top-level directory of the Nerpa repository):

   ```bash
   tar -xzf nerpa-2.1-scripts.tar.gz
   ```

   Make sure the extracted files/folders end up inside the Nerpa repository root.

2. **Activate the Nerpa conda environment**:

   ```bash
   conda activate nerpa-env
   ```

## Running the workflow

From the **Nerpa root directory**, run:

```bash
snakemake --cores 1
```

That’s it — Snakemake will execute the training/testing/benchmarking pipeline with a single core.

## Output

The workflow produces several output directories. The most useful ones are:

- `benchmarking_results/{bgcs_test_set}/nerpa1_vs_nerpa2_vs_biocat_plots/`  
  Various benchmarking plots comparing **Nerpa 1**, **Nerpa 2**, and **BioCAT**.  
  Here `{bgcs_test_set}` is one of:
  - `training_bgcs` — BGCs from the manually curated alignments used for training
  - `mibig4_wo_training_bgcs` — MIBiG 4 BGCs dissimilar from those used for training (**main benchmarking set**)

- `training_results/`  
  Training output: configuration files/parameters calculated from the manually curated alignments.

- `test_results/`  
  Testing output: Nerpa errors on manually curated alignments.
