python nerpa.py \
-a /share/home/azat/nerpa/asdb-jsons \
--nrp-variants ./data/input/preprocessed/pnrpdb2_deduplicated_preprocessed.yaml \
-o ./nerpa_results/antismashdb/ \
--force-output-dir \
--fast-matching \
--process-hybrids \
--max-num-matches-per-bgc 10 \
--max-num-matches 10000 \
--threads 32