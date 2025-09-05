# Nerpa 2.0 Manual

1. [About Nerpa](#sec_about) </br>
   1.1 [Nerpa pipeline](#sec_about_pipeline)</br>
   1.2 [Supported data types](#sec_about_data)</br>
2. [Installation](#sec_install)</br>
    2.1. [Prerequisites](#sec_install_prereq)</br>
    2.2. [Installation from tarball](#sec_install_source)</br>
    2.3. [Verifying your installation](#sec_install_verify)</br>
3. [Running Nerpa](#sec_run)</br>
    3.1. [Quick start](#sec_run_quick)</br>
    3.2. [Command-line options](#sec_run_options)</br>
    3.3. [Output files](#sec_run_results)</br>
4. [Citation](#sec_cite)</br>
5. [Feedback and bug reports](#sec_feedback)</br>

<a name="sec_about"></a>
# About Nerpa

Nerpa is a tool for linking biosynthetic gene clusters (BGCs) to known nonribosomal peptide (NRP) structures.
You can read more about the Nerpa algorithm and the practical applications of the tool in [our papers](#sec_cite).
Nerpa is currently developed and maintained by [Gurevich Lab](https://helmholtz-hips.de/en/hmsb) 
at the [Helmholtz Institute for Pharmaceutical Research Saarland (HIPS)](https://helmholtz-hips.de/en/) 
and the [Center for Bioinformatics Saar (CBI)](https://zbi-www.bioinf.uni-sb.de/en/).

This manual will help you to install and run the tool. Nerpa version 2.0.0 was released on 19.03.2025. 
The tool is dual-licensed and is available under GPLv3 or Creative Commons BY-NC-SA 4.0, see [LICENSE.txt](LICENSE.txt).

<a name="sec_about_pipeline"></a>
## Nerpa pipeline
The simplified Nerpa pipeline is depicted in the figure below.

![Nerpa pipeline](docs/img/pipeline.png "Nerpa 2 pipeline")

Nerpa takes as input an NRP structure database and genome sequences.
The pipeline goes as follows:  
1. Construct tentative NRP synthetase assembly lines 
along with respective sequences of genome-predicted residues
(using [antiSMASH](https://academic.oup.com/nar/article/51/W1/W46/7151336)).  
2. Construct representations of the database structures as monomer graphs
(using [rBAN](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0335-x)).  
3. Build HMMs for genome-predicted NRP synthetase assembly lines
as described in [the Nerpa 2 paper](https://doi.org/10.1101/2024.11.19.624380).
4. Extract NRP linearizations from the monomer graphs.  
5. Score the NRP linearizations against the HMMs all-vs-all manner 
(using the [Viterbi algorithm](https://en.wikipedia.org/wiki/Viterbi_algorithm)).  
6. Create an interactive report with the best matches and detailed alignments.  

 
<a name="sec_about_data"></a>
## Supported data types
For **genome sequences**:  

* *Recommended*: complete **antiSMASH output** after processing your raw genome sequence
(e.g., downloaded from the [antiSMASH web server](https://antismash.secondarymetabolites.org/));
or antiSMASH job IDs (in this case, Nerpa will download it automatically).
* *Also accepted*: raw genome sequences in the **FASTA and GenBank** formats; 
in this case, Nerpa will predict NRP BGCs in them with antiSMASH 
(should be installed separately and present in `PATH` or provided to Nerpa via `--antismash-installation-dir`).

For **NRP structures**:

* *Recommended*: [**isomeric SMILES format**](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system);
Nerpa distinguishes between L- and D-configurations of amino acids, so the use of the isomeric format leads to more accurate results.
* *Also accepted*: any other [**SMILES**](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system), 
i.e., without stereochemistry information.

*Note:* you can use free online converters to get (isomeric) SMILES from other popular chemical formats such as MDL MOL or InChI,
e.g., [this one from UNM](https://datascience.unm.edu/tomcat/biocomp/convert). Alternatively, there are many command-line convertors, e.g. [molconvert](https://docs.chemaxon.com/display/docs/molconvert.md), or programming libraries, e.g. [RDKit](https://www.rdkit.org/).


<a name="sec_install"></a>
# Installation
<a name="sec_install_prereq"></a>
## Prerequisites

* **(Required)** Nerpa relies on **Java** (to run the embedded rBAN), **Python v3.10 or higher**, and a number of **Python dependencies specified in the [environment.yml](environment.yml)** file.  
  We highly recommend installing [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to easily set up all dependencies, as demonstrated below.

* **(Optional)** If you plan to use Nerpa with raw genome sequences (FASTA or GenBank) rather than antiSMASH-processed files, 
you will also need to [install antiSMASH locally](https://docs.antismash.secondarymetabolites.org/install/).  
Alternatively, you can use the [antiSMASH web server](https://antismash.secondarymetabolites.org/).

* **(Optional)** Nerpa is quite fast by default, but we provide an even faster C++ implementation. 
To use it, you will need a **[C++20 compiler](https://isocpp.org/get-started)** and **[CMake v3.10](https://cmake.org/)** or higher.  

<a name="sec_install_source"></a>
## Installation from tarball

First, download and unpack the release tarball:

```commandline
wget https://github.com/gurevichlab/nerpa/releases/download/nerpa_2.0.0/nerpa-2.0.0.tar.gz
tar -xzf nerpa-2.0.0.tar.gz
cd nerpa-2.0.0
```
Next, install all required dependencies. We recommend creating and activating a Conda environment:  

```commandline
conda env create -f environment.yml
conda activate nerpa-env
```
Finally, if you want to use the fast C++ version (optional), compile it by running:

```commandline
bash install.sh
```

<a name="sec_install_verify"></a>
## Verifying your installation

We recommend adding the `nerpa` directory to `PATH`. In this case, you can run Nerpa simply as `nerpa.py` from anywhere;
otherwise, you would need to specify path from the current directory to `./nerpa.py`.
All running examples below assume that Nerpa is in `PATH`.

To test your installation, first, try to get the list of the Nerpa command-line options:  

```commandline
nerpa.py -h
```

Then, try any example from the [Quick start](#sec_run_quick) section and ensure the log contains no error messages.

If you have any problems, please do not hesitate to [contact us](#sec_feedback).

<a name="sec_run"></a>
# Running Nerpa

<a name="sec_run_quick"></a>
## Quick start

Sample test data with three antiSMASH-processed BGCs and three NRP structures in the SMILES format is included in the release tarball.  
Alternatively, you can download it from [here](https://github.com/gurevichlab/nerpa/releases/download/nerpa_2.0.0/test_data.tar.gz) and unpack it in your current working directory.

To run Nerpa on the test data, execute:

```commandline
nerpa.py -a test_data/antismash --smiles-tsv test_data/smiles.tsv --col-id ID --output-dir nerpa_results/test_run
```

The output will be saved in the `nerpa_results/{CURRENT_TIME}` directory and symlinked to `nerpa_results/latest` for your convenience.  
For details on the output directory contents and their interpretation refer to the [corresponding section](#sec_run_results).

<a name="sec_run_options"></a>
## Command-line options

To see the full list of available options, type

```commandline
nerpa.py -h
```
	
All options are divided into four categories. The most important options in each category are listed below.

### Genomic input (genome sequences or BGCs)

The most convenient way to obtain antiSMASH predictions of BGCs in your genomic data is to upload your  
FASTA or GenBank file to their [web server](https://antismash.secondarymetabolites.org/).  
Once the server job is completed, download the results (`Download -> Download all results`), unpack the archive,  
and provide the path to the unpacked directory using the `-a` option.

Alternatively, you can provide Nerpa with the server job ID (e.g., `bacteria-2a9bb79e-e804-42c9-bb62-516cac47eca2`)  
via the `--antismash-job-ids` option, and Nerpa will download everything automatically.  
For multiple jobs, specify them as a space-separated list of IDs.

You can also use the [command-line version of antiSMASH](https://docs.antismash.secondarymetabolites.org/install/).  
Nerpa has been tested with outputs from antiSMASH version 7 (7.0.0 and 7.1.0).  
If antiSMASH is installed on your system, you can provide raw genome sequences in FASTA or GenBank format via the `--genome` option,  
and Nerpa will run antiSMASH automatically.  
To enable this, antiSMASH should be available in your system’s `PATH` variable, or the path to its installation directory should be specified via the `--antismash-installation-path` option.

Note that you can specify an unlimited number of antiSMASH output files by either:  

- Using the `-a` option multiple times.  
- Specifying a root directory containing many inputs.  
- Writing paths to all antiSMASH outputs in a single text file and providing it via the `--antismash-paths-file` option.  

### Chemical input (compounds)

NRP molecules should be specified in the [SMILES format](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system).  
You can provide them in one of the following ways:

- As a space-separated list of SMILES strings using the `--smiles` option.
- In a multi-column file specified via the `--smiles-tsv` option.

In the latter case, the default column separator (`\t`), names of the SMILES column (`SMILES`) and the column with molecule IDs 
(*row index*) could be adjusted using the `--sep`, `--col-smiles`, and `--col-id` options, respectively.

The Nerpa release package comes with a set of NRP databases in the SMILES format:  

- Compounds from [MIBiG 4.0](https://mibig.secondarymetabolites.org/) and [Norine](https://bioinfo.cristal.univ-lille.fr/norine/index.jsp), available in [data/mibig_norine.tsv](data/mibig_norine.tsv).  
- Our own database of putative NRP structures, pNRPdb, available in [data/pnrpdb2rc1_summary.tsv](data/input/pnrpdb2rc1_summary.tsv).  

### Advanced Input

You can reuse preprocessed BGCs and/or chemical structures from a previous Nerpa run.  
This is useful, for example, if you want to screen the same BGCs against different NRPs, or vice versa.

The preprocessed files are stored in the Nerpa output directory in:  

- `BGC_variants/` (for BGCs).  
- `NRP_variants/` (for NRPs).  

To reuse them, provide the corresponding paths via the `--bgc-variants` and `--nrp-variants` options.

### Pipeline Options

- `--output_dir <DIR>, -o <DIR>`  
  Path to the output directory.  
  If the directory already exists, Nerpa will exit with an error unless `--force-output-dir` is specified.  
  If not set, Nerpa will create the directory `nerpa_results/{CURRENT_TIME}` and symlink it to `nerpa_results/latest`.

- `--process-hybrids`  
  Process NRP-polyketide hybrid monomers (requires rBAN to be used). Recommended.

- `--threads`  
  Number of threads for running Nerpa. Default: `1`.

- `--skip-molecule-drawing`  
  Disable drawing of NRP compounds (they will not appear in the HTML report). Enabling this option speeds up the run and reduces output size.

- `--fast-matching`  
  Enable the fast C++-based matching (requires pre-compilation; see the [Installation](#sec_install) section).


<a name="sec_run_results"></a>
## Output Files

The key files and directories inside the Nerpa output directory (`--output-dir`) are:

- `report.html`  
  An interactive HTML report showing the best Nerpa matches, along with the corresponding annotated BGCs and NRPs.

- `report.tsv`  
  A tab-separated file containing matched NRP-BGC pairs with their corresponding scores.

- `BGC_variants/`  
  Directory containing preprocessed antiSMASH outputs. These can be reused for another run via the `--bgc-variants` option.

- `NRP_variants/`  
  Directory containing preprocessed compounds.  These can be reused for another run via the `--nrp-variants` option.

<a name="sec_cite"></a>
## Citation

If you use Nerpa in your research, please cite our papers:  
Nerpa v.2 is described in [Olkhovskii et al, bioRxiv 2024](https://doi.org/10.1101/2024.11.19.624380).  
Nerpa v.1 is published in [Kunyavskaya, Tagirdzhanov et al., Metabolites 2021](https://doi.org/10.3390/metabo11100693).

<a name="sec_feedback"></a>
## Feedback and bug reports
You can leave your comments and bug reports at [https://github.com/gurevichlab/nerpa/issues](https://github.com/gurevichlab/nerpa/issues) (*recommended way*) 
or sent it via e-mail to [alexey.gurevich@helmholtz-hips.de](alexey.gurevich@helmholtz-hips.de).

Your comments, bug reports, and suggestions are **very welcomed**.
They will help us to improve Nerpa further.
In particular, we would love to hear your thought on desired features of the future Nerpa web service.

If you have any troubles running Nerpa, please attach `nerpa.log` from the output directory.

