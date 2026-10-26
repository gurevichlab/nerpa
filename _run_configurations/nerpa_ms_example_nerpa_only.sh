# ARGS=""
(
    cd -- "$(dirname -- "$0")/.." || exit
    python nerpa.py \
           -a test_data/NC_008095.1 \
           --smiles-tsv test_data/myxochromides_wo_A.tsv \
           --col-id ID \
           --output-dir nerpa_results/nerpa_ms_test_nerpa_only \
           --force-output-dir \
           --let-it-crash
)
