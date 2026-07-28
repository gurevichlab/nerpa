# ARGS=""
python nerpa.py \
       -a /home/ilianolhin/work/data/NC_008095.1 \
       --smiles-tsv /home/ilianolhin/work/data/myxoschromides_wo_A.tsv \
       --col-id ID \
       --output-dir nerpa_results/myxochromides_wo_A \
       --force-output-dir \
       --max-num-matches 5 \
       --let-it-crash

