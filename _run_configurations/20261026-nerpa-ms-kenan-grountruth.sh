# ARGS=""
(
    cd -- "$(dirname -- "$0")/.." || exit
    python nerpa_ms.py \
           --nerpa-results nerpa_results/for_kenan_2026-10-08 \
           --spectra /home/ilianolhin/work/data/kenan-msms/msms \
           --output-dir nerpa_results/20261026-nerpa-ms-kenan-groundtruth \
           --max-edits 1 \
           --nerpa-match-ids-for-candidate-generation "4, 6" \
           --num-variants-per-num-edits 5
    # --max-nerpa-matches-for-candidate-generation 5 \
    #       --smiles-tsv test_data/myxochromides_wo_A.tsv \

)
