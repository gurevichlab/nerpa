# ARGS=""
(
    cd -- "$(dirname -- "$0")/.." || exit
    OUTPUT_FOLDER="nerpa_results/nerpa_ms_test_run_reuse_nerpa_output"
    mkdir -p -- "$OUTPUT_FOLDER"

    python -u nerpa_ms.py \
           --nerpa-results nerpa_results/nerpa_ms_test_nerpa_only \
           --spectra test_data/myxococcus_xanthux.mgf \
           --output-dir "$OUTPUT_FOLDER" \
           --max-edits 1 \
           --nerpa-match-ids-for-candidate-generation "0, 2" \
           --num-variants-per-num-edits 5 \
           2>&1 | tee "${OUTPUT_FOLDER}.log"

    mv -- "${OUTPUT_FOLDER}.log" "$OUTPUT_FOLDER/nerpa_ms.log"

    # --max-nerpa-matches-for-candidate-generation 5 \
    #       --smiles-tsv test_data/myxochromides_wo_A.tsv \

)
