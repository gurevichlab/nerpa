#!/usr/bin/env bash
# Install the PARAS model paired with scikit-learn 1.8.0 in environment.yml.
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
paras_model_path="${1:-$repo_root/external_tools/paras/model.paras}"
paras_model_sha256='1f97ea44a9bf7a59a9fc83b1e9e06cc76d0c942318fbcb1534fe94034064fd7f'
paras_model_url="${PARAS_MODEL_URL:-}"

matches_model() {
    [[ -f "$1" ]] && [[ "$(sha256sum -- "$1" | cut -d ' ' -f 1)" == "$paras_model_sha256" ]]
}

if matches_model "$paras_model_path"; then
    echo 'PARAS model checksum verified.'
    exit 0
fi

mkdir -p -- "$(dirname -- "$paras_model_path")"
paras_tmp_dir="$(mktemp -d "${paras_model_path}.download.XXXXXX")"
trap 'rm -rf -- "$paras_tmp_dir"' EXIT

# Reuse a previously downloaded archive when upgrading an offline checkout.
if [[ -f "${paras_model_path}.gz" ]] && \
        gzip -dc -- "${paras_model_path}.gz" > "$paras_tmp_dir/model.paras" && \
        matches_model "$paras_tmp_dir/model.paras"; then
    echo 'Using the verified local PARAS archive.'
else
    if [[ -z "$paras_model_url" ]]; then
        echo 'Set PARAS_MODEL_URL to the matching model.paras.gz download URL, or place the archive beside model.paras.' >&2
        exit 1
    fi
    curl --location --fail --retry 2 --connect-timeout 30 \
        --output "$paras_tmp_dir/model.paras.gz" "$paras_model_url"
    gzip -dc -- "$paras_tmp_dir/model.paras.gz" > "$paras_tmp_dir/model.paras"
    if ! matches_model "$paras_tmp_dir/model.paras"; then
        echo 'PARAS model checksum mismatch; the existing model has been preserved.' >&2
        exit 1
    fi
fi

if [[ -e "$paras_model_path" ]]; then
    paras_backup="$(mktemp "${paras_model_path}.previous.XXXXXX")"
    cp -p -- "$paras_model_path" "$paras_backup"
    echo "Previous PARAS model saved to $paras_backup"
fi
# The temporary file is on the same filesystem, so installation is atomic.
mv -- "$paras_tmp_dir/model.paras" "$paras_model_path"
echo 'PARAS model installed and checksum verified.'
