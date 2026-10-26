# ARGS=""
(
    cd -- "$(dirname -- "$0")/.." || exit
    python nerpa_ms.py \
           --nerpa-results nerpa_results/nerpa_ms_test_nerpa_only \
           --spectra test_data/myxococcus_xanthux.mgf \
           --output-dir nerpa_results/nerpa_ms_test_run_reuse_nerpa_output \
           --max-edits 1 \
           --nerpa-match-ids-for-candidate-generation "0, 2" \
           --num-variants-per-num-edits 5
    # --max-nerpa-matches-for-candidate-generation 5 \
    #       --smiles-tsv test_data/myxochromides_wo_A.tsv \

)
