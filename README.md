##Dependencies

### Core

* Linux or macOS  
* Python 3 (tested with Python 3.10)  
* Java (tested with OpenJDK 21)  

### Required Python modules 
(the list is probably incomplete)

They can be installed via:

```
pip install dacite parse more_itertools pandas rdkit networkx matplotlib  
```




## Use case example 

(run from the nerpa root directory):

```
python nerpa.py -a ./test_data/MIBiG_subset/genome_predictions_v7 --smiles-tsv ./test_data/MIBiG_subset/structures_info.tsv --col-id ID -o ./nerpa_results/test_example --force-existing-outdir
```
