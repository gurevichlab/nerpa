#!/usr/bin/env python3
"""
Script to compare Nerpa performance with and without anchors heuristic.
Single-threaded version to avoid multiprocessing issues.
"""
import os
import sys
import time
import subprocess
from pathlib import Path
import argparse


def run_nerpa(antismash_path, rban_path, output_dir, num_matches, use_anchors=False):
    """
    Run Nerpa with specified parameters and measure execution time.
    """
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build the command
    command = [
        "python", "nerpa.py",
        "--antismash", str(antismash_path),
        "--rban-json", str(rban_path),
        "--output_dir", str(output_dir),
        "--force-existing-outdir",
        "--threads", "1",  # Force single thread for reliability
        "--num-matches", str(num_matches)
    ]

    # Add anchors heuristic flag if specified
    if use_anchors:
        command.append("--use-anchors-heuristic")

    # Print command for reference
    print(f"Running command: {' '.join([str(c) for c in command])}")

    # Time the execution
    start_time = time.time()

    try:
        # Run the command
        subprocess.run(command, check=True)

        end_time = time.time()
        execution_time = end_time - start_time

        # Print success message
        print(f"Command completed successfully in {execution_time:.2f} seconds")

        return execution_time

    except subprocess.CalledProcessError as e:
        print(f"Error running Nerpa: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Compare Nerpa performance with and without anchors heuristic")
    parser.add_argument("--antismash", required=True, help="Path to antiSMASH input directory")
    parser.add_argument("--rban-json", required=True, help="Path to rBAN JSON file")
    parser.add_argument("--output-dir", default="./output", help="Base output directory")
    parser.add_argument("--num-matches", type=int, default=10, help="Maximum number of matches per BGC")

    args = parser.parse_args()

    # Convert paths to Path objects
    antismash_path = Path(args.antismash)
    rban_path = Path(args.rban_json)
    output_base_dir = Path(args.output_dir)

    # Create output directories
    standard_dir = output_base_dir / "standard"
    anchors_dir = output_base_dir / "with_anchors"

    # Run without anchors heuristic
    print("\n=== Running without anchors heuristic ===")
    standard_time = run_nerpa(
        antismash_path,
        rban_path,
        standard_dir,
        args.num_matches,
        use_anchors=False
    )

    # Run with anchors heuristic
    print("\n=== Running with anchors heuristic ===")
    anchors_time = run_nerpa(
        antismash_path,
        rban_path,
        anchors_dir,
        args.num_matches,
        use_anchors=True
    )

    # Compare results
    print("\n=== Comparison Results ===")
    print(f"Without anchors: {standard_time:.2f} seconds")
    print(f"With anchors: {anchors_time:.2f} seconds")

    if anchors_time < standard_time:
        speedup = (standard_time - anchors_time) / standard_time * 100
        print(f"Anchors heuristic is faster by {speedup:.2f}%")
    else:
        slowdown = (anchors_time - standard_time) / standard_time * 100
        print(f"Anchors heuristic is slower by {slowdown:.2f}%")


if __name__ == "__main__":
    main()