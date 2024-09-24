## How to augment signatures.tsv with SVM predictions

This instruction will help you to run `scripts/signatures_extender.py` and create `data/extended_signatures.tsv` from `data/signatures.tsv`.
The extended file 
- adds the header line, 
- picks only one out of two columns with A domain specificities 
(by default: column #3, i.e., all experimentally validated substrates that can be attracted by the A domain;
column #4 contains substrates observed in the corresponding NRP final structures) and converts them to a proper Python list,
- substitutes the column with A domain IDs with just their total number,
- adds 28 extra columns corresponding to each of the SVM-supported substrates. The values are 4-tuples 
with SVM prediction scores from 'single aa' to 'small cluster' to 'large cluster' to 'physicochemical_class' 
(or 0.0 if the corresponding score was not reported, e.g., only one substrate has a non-zero 'single aa' score).

### Prerequisites

1. Source code of antiSMASH (e.g., from <https://github.com/antismash/antismash/>). You should add it to `sys.path` or mark it as a Source Root in PyCharm (see comments inside the script)
2. The signatures file (available in this folder as `signatures.tsv`, downloaded from <https://dl.secondarymetabolites.org/releases/stachelhaus/1.1/signatures.tsv.xz>)
3. SVM models (available in this folder as `models.tar.xz`, downloaded from <https://dl.secondarymetabolites.org/releases/nrps_svm/2.0/models.tar.xz>; should be unpacked before use!)
4. Python3 with `nrpys` package installed (<https://pypi.org/project/nrpys/>)

Note: the most actual versions of (2) & (3) are likely specified in <https://github.com/antismash/antismash/blob/master/antismash/download_databases.py> (as `STACHELHAUS_URL` and `NRPS_SVM_URL`).

### Running

Note: the script contains hard-coded paths to input and output file locations.
The running is simply   
`python scripts/signatures_extender.py`  

The resulting file (`extended_signatures.tsv`) will be saved near the original input file (`data/signatures.tsv`)

### Important note

It seems there is a bug in the antiSMASH wrapper of NRPys.
(<https://github.com/antismash/antismash/blob/7-1-stable/antismash/modules/nrps_pks/nrpys.py#L189>)
One chemical class misses two its amino acids ('Sal' and 'Trp'). 
As a result, their scores are not reported at the 'class' level, only at the 'large_cluster' level and below.

## How to create the big training table for the specificity prediction model

### Prerequisites
1. Several data files available in the `data` direct ory (see more details in the help message).
2. Python3 with basic packages (like `pandas` and `pyyaml`) installed. 

### Running

Note: the script requires many input/output paths that can be changed via options. 
Please use `-h` to show all of them. 
Two of these files, `monomers_unique.tsv` (aka Monomers Table) and `core_frequency.tsv` (frequencies of the core Nerpa monomers in the `signatures.tsv` file) are discussed (and sometimes updated) [here](https://github.com/gurevichlab/nerpa/issues/18).
If unsure, use default values: all required files are available in the `data` directory already. 

The running is simply   
`python scripts/training_table_generator.py --num-monomers 49`  (the default num monomers is `-1` that correspond to "all"). 

The resulting files (`monomer_signatures.yaml` and `nerpa_scoring_table.tsv`) will be saved in the `data` directory by default.