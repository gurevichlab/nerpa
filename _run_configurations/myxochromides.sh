# ARGS=""
python nerpa.py \
       -a /home/ilianolhin/work/data/NC_008095.1 \
       --smiles-tsv /home/ilianolhin/work/data/myxochromides.tsv \
       --col-id ID \
       --output-dir nerpa_results/myxochromides \
       --force-output-dir \
       --max-num-matches 100 \
       --let-it-crash

