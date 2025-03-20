#!/bin/bash

# Exit immediately if a command fails
set -e

paras_model_link='https://zenodo.org/records/13165500/files/model.paras?download=1'
# q: if ./external_tools/paras/model.paras doesn't exist, download it
if [ ! -f "./external_tools/paras/model.paras" ]; then
    echo "Downloading model.paras from $paras_model_link..."
    wget -O ./external_tools/paras/model.paras $paras_model_link
fi

# Remove the build directory if it exists
if [ -d "./build" ]; then
    echo "Removing existing build directory..."
    rm -rf ./build
fi

# Create and enter the build directory
mkdir build && cd build

# Run CMake and Make
echo "Running CMake..."
cmake ..

echo "Building project..."
make

echo "Build completed successfully!"