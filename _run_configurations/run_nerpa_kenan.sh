python nerpa.py \
-a /share/home/azat/nerpa/asdb-jsons \
--smiles-tsv \
/share/home/azat/nerpa/nerpa2-antismashdb-kenan/compounds_kenan.tsv \
--col-id ID \
-o ./nerpa_results/antismashdb_kenan/ \
--force-output-dir \
--fast-matching \
--process-hybrids \
--min-num-matches-per-nrp 10 \
--max-num-matches 1000 \
--threads 32 \
--skip-molecule-drawing
