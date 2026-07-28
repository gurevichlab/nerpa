# ARGS=""
python nerpa_ms.py \
       -a test_data/NC_008095.1 \
       --smiles-tsv test_data/myxochromides_wo_A.tsv \
       --col-id ID \
       --spectra test_data/myxococcus_xanthux.mgf \
       --output-dir nerpa_results/nerpa_ms_test_run \
       --force-output-dir \
       --let-it-crash \

