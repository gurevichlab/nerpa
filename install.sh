#!/usr/bin/env bash

# Exit immediately if a command fails
set -euo pipefail

# Set the repository root directory to the directory of this script
repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo_root"

# Exit immediately if a command fails
set -euo pipefail

# ===== 1. Clone kakapo repository if it doesn't exist
kakapo_repo_dir="$repo_root/external_tools/kakapo"

if [[ ! -d "$kakapo_repo_dir" ]]; then
    echo "Cloning kakapo release 0.1.0b1 to $kakapo_repo_dir..."
    git clone --branch 0.1.0b1 https://github.com/gurevichlab/kakapo.git "$kakapo_repo_dir"
    echo "Kakapo release 0.1.0b1 cloned"
fi

# ===== 2. Set up environments

nerpa_env_prefix="$repo_root/envs/nerpa"
kakapo_env_prefix="$repo_root/envs/kakapo"
bin_dir="$repo_root/bin"

mkdir -p "$bin_dir"

# 2.1 Install the pinned runtime, including when upgrading an existing checkout.
if [[ ! -d "$nerpa_env_prefix" ]]; then
  conda env create -p "$nerpa_env_prefix" -f "$repo_root/environment.yml"
else
  conda env update --prune -p "$nerpa_env_prefix" -f "$repo_root/environment.yml"
fi

# 2.2 Create the Kakapo environment if it doesn't exist
if [[ ! -d "$kakapo_env_prefix" ]]; then
  conda create -y -p "$kakapo_env_prefix" python=3.12 pip
fi

"$kakapo_env_prefix/bin/python" -m pip install -U pip
"$kakapo_env_prefix/bin/python" -m pip install -e "$repo_root/external_tools/kakapo"


## ===== 3. Build C++ code

# Remove the build directory if it exists
build_dir="$repo_root/build"

if [[ -d "$build_dir" ]]; then
    echo "Removing existing build directory..."
    rm -rf "$build_dir"
 fi

# Create and enter the build directory
mkdir build && cd build

# Run CMake and Make
echo "Running CMake..."
cmake -S "$repo_root" -B "$build_dir"

echo "Building project..."
cmake --build "$build_dir"

echo "Build completed successfully!"


# ===== 4. Install and validate the PARAS model paired with this runtime.
bash "$repo_root/scripts/installation/install_paras.sh"
PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}" "$nerpa_env_prefix/bin/python" - <<'PY'
from src.config import load_config
from src.paras.paras_wrapper import ParasWrapper

model = ParasWrapper(load_config().specificity_prediction_config.paras_model, sklearn_version='1.8.0')
model.predict('A' * 34)
print('PARAS model loading and prediction verified.')
PY

echo "Installation finished successfully!"
