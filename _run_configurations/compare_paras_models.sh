nerpa_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
(
    cd "$nerpa_dir" || exit 1
    PYTHONPATH=. python scripts/compare_paras_models.py \
        --table training_results/hmm_chosen_steps.tsv \
        --model-1 external_tools/paras/model.paras.bak \
        --python-1 envs/paras-sklearn-1.2/bin/python \
        --sklearn-version-1 1.2.0 \
        --model-2 external_tools/paras/model.paras \
        --python-2 envs/nerpa/bin/python \
        --sklearn-version-2 1.8.0 \
        --nerpa-root "$nerpa_dir" \
        --outdir tmp/paras_model_comparison
)
