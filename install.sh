#!/usr/bin/env bash

# Exit immediately if a command fails
set -euo pipefail

# Set the repository root directory to the directory of this script
repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo_root"

# Exit immediately if a command fails
set -euo pipefail


# ===== 1. Set up Nerpa and Kakapo environments

nerpa_env_prefix="$repo_root/envs/nerpa"
kakapo_env_prefix="$repo_root/envs/kakapo"
bin_dir="$repo_root/bin"

mkdir -p "$bin_dir"

# 1.1 Create the Nerpa environment if it doesn't exist
if [[ ! -d "$nerpa_env_prefix" ]]; then
  conda env create -p "$nerpa_env_prefix" -f "$repo_root/environment.yml"
fi

# 1.2 Create the Kakapo environment if it doesn't exist
if [[ ! -d "$kakapo_env_prefix" ]]; then
  conda create -y -p "$kakapo_env_prefix" python=3.12 pip
fi

"$kakapo_env_prefix/bin/python" -m pip install -U pip
"$kakapo_env_prefix/bin/python" -m pip install -e "$repo_root/external_tools/kakapo"


## ===== 2. Build C++ code

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


# ===== 3. Download PARAS model if not present
paras_model_link='https://zenodo.org/records/13165500/files/model.paras?download=1'

paras_model_path="$repo_root/external_tools/paras/model.paras"

if [[ ! -f "$paras_model_path" ]]; then
    echo "Downloading model.paras from $paras_model_link..."
    mkdir -p "$(dirname -- "$paras_model_path")"
    curl -L --fail -o "$paras_model_path" "$paras_model_link"
    echo "PARAS model downloaded"
fi

echo "Installation finished successfully!"

