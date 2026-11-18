#!/usr/bin/env bash

# Exit immediately if a command fails
set -euo pipefail

# Set the repository root directory to the directory of this script
repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo_root"

# Exit immediately if a command fails
set -euo pipefail

# ===== 1. Set up the Nerpa environment

nerpa_env_prefix="$repo_root/envs/nerpa"

if command -v mamba >/dev/null 2>&1; then
  environment_manager="mamba"
elif command -v conda >/dev/null 2>&1; then
  environment_manager="conda"
else
  echo "Error: neither mamba nor conda is available." >&2
  exit 1
fi

echo "Using $environment_manager to set up the Nerpa environment."

# Install the pinned runtime, including when upgrading an existing checkout.
if [[ ! -d "$nerpa_env_prefix" ]]; then
  "$environment_manager" env create -p "$nerpa_env_prefix" -f "$repo_root/environment.yml"
else
  "$environment_manager" env update --prune -p "$nerpa_env_prefix" -f "$repo_root/environment.yml"
fi


# ===== 2. Build C++ code

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


# ===== 3. Install and validate the PARAS model paired with this runtime.
bash "$repo_root/scripts/installation/install_paras.sh"
PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}" "$nerpa_env_prefix/bin/python" - <<'PY'
from src.config import load_config
from src.paras.paras_wrapper import ParasWrapper

model = ParasWrapper(load_config().specificity_prediction_config.paras_model, sklearn_version='1.8.0')
model.predict('A' * 34)
print('PARAS model loading and prediction verified.')
PY

echo "Installation finished successfully!"
