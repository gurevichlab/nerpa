Install dependencies (Python 3):
```
pip install dacite parse more_itertools pandas rdkit networkx  
```


Use case example (run from the nerpa main directory):

```
python nerpa.py -a ./test_data/MIBiG_subset/genome_predictions_v7 --smiles-tsv ./test_data/MIBiG_subset/structures_info.tsv --col-id ID -o ./nerpa_results/test_example --force-existing-outdir
```
