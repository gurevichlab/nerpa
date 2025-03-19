#!/bin/bash

# Exit immediately if a command fails
set -e

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