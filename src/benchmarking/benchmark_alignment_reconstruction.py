import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import yaml
from typing import List, Dict
from src.testing.testing_nerpa1 import SimplifiedMatch, load_nerpa1_matches
from src.testing.testing_types import TestMatch
from src.testing.simplified_alignment import get_mismatched_steps, simplified_alignment_from_match
from src.matching.match_type import Match
from src.testing.simplified_match import nerpa2_match_to_simplified_match

def calculate_error_counts(matches: List[SimplifiedMatch],
                           matches_label: str,
                           tests: List[TestMatch],
                           remove_outliers: bool = False) -> List[int]:
    """
    Calculate the number of mismatches between each test alignment and the corresponding match.

    Args:
        matches: List of SimplifiedMatch objects
        tests: List of TestMatch objects containing true alignments

    Returns:
        List of error counts for each match
    """
    # Create a dict of TestMatch by (bgc_id, nrp_id) for quick lookup
    matches_by_bgc_nrp = {(m.bgc_id, m.nrp_id): m for m in matches}

    error_counts = []
    for tm in tests:
        key = (tm.bgc_id, tm.nrp_id)
        m = matches_by_bgc_nrp.get(key)
        if m is None:
            print(f"Warning: No match for ({tm.bgc_id}, {tm.nrp_id}) in {matches_label}")
            continue

        mismatches = get_mismatched_steps(tm.true_alignment, m.alignment)
        error_counts.append(len(mismatches))

    if remove_outliers:
        q1 = np.percentile(error_counts, 25)
        q3 = np.percentile(error_counts, 75)
        iqr = q3 - q1
        #lower_bound = q1 - 5 * iqr
        #upper_bound = q3 + 5 * iqr
        lower_bound = 0  # stub
        upper_bound = 7
        error_counts = [ec for ec in error_counts if lower_bound <= ec <= upper_bound]

    return error_counts



def plot_error_histograms(error_counts_dict: Dict[str, List[int]],
                         filename: Path):
    """
    Plot histograms for multiple sets of error counts.

    Args:
        error_counts_dict: Dictionary mapping match set names to lists of error counts
        filename: Path to save the histogram image
    """
    plt.figure(figsize=(12, 6))

    max_errors = max(max(error_counts) for error_counts in error_counts_dict.values())
    bins = range(0, max_errors + 2)

    num_sets = len(error_counts_dict)
    bar_width = 0.8 / num_sets

    for idx, (name, error_counts) in enumerate(error_counts_dict.items()):
        counts, _ = np.histogram(error_counts, bins=bins)
        x_positions = [x + idx * bar_width for x in range(len(bins) - 1)]
        plt.bar(x_positions, counts, width=bar_width, label=name, alpha=0.8)

    plt.title('Histogram of Alignment Errors')
    plt.xlabel('Number of Errors')
    plt.ylabel('Frequency')
    plt.xticks(range(len(bins) - 1), range(max_errors + 1))
    plt.legend()
    plt.grid(True, alpha=0.3, axis='y')
    filename.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(filename, bbox_inches='tight')
    plt.close()


def benchmark_alignments(match_sets: Dict[str, List[SimplifiedMatch]],
                         tests: List[TestMatch],
                         plots_dir: Path) -> None:
    """
    Benchmark multiple sets of matches against true alignments and visualize the results.

    Args:
        match_sets: Dictionary mapping names to lists of SimplifiedMatch objects
        tests: List of TestMatch objects containing true alignments
        plots_dir: Path to save the plots
    """
    error_counts_dict = {}
    for name, matches in match_sets.items():
        error_counts = calculate_error_counts(matches, name, tests,
                                              remove_outliers=True)
        error_counts_dict[name] = error_counts

    histogram_file = plots_dir / 'alignment_reconstruction_histogram.png'
    boxplots_file = plots_dir / 'alignment_reconstruction_boxplots.png'

    plot_error_histograms(error_counts_dict, histogram_file)


def load_nerpa2_matches(results_dir: Path) -> List[SimplifiedMatch]:
    matches_file = results_dir / 'matches_details' / 'matches.yaml'
    match_dicts = yaml.safe_load(matches_file.read_text())
    simplified_matches = [
        nerpa2_match_to_simplified_match(Match.from_dict(match_dict))
        for match_dict in match_dicts
    ]

    return simplified_matches


def nerpa1_vs_nerpa2(nerpa_dir: Path,
                     args: argparse.Namespace) -> None:
    # Load nerpa1 matches
    nerpa1_matches = load_nerpa1_matches(args.nerpa1_results)

    # Load nerpa2 matches
    nerpa2_matches = load_nerpa2_matches(args.nerpa2_results)

    test_matches = [TestMatch.from_yaml_dict(item)
                    for item in yaml.safe_load(args.approved_matches_yaml.read_text())]

    # Create a dictionary of match sets
    match_sets = {
        'Nerpa1': nerpa1_matches,
        'Nerpa2': nerpa2_matches,
    }

    # Benchmark and plot
    benchmark_alignments(match_sets, test_matches, args.plots_dir)

def parse_args() -> argparse.Namespace:
    """
    Parse CLI arguments for nerpa paths.

    - You can provide a single --nerpa-dir to derive the original defaults, or
      override any path individually with the specific flags.
    - Returns an argparse.Namespace with attributes:
      nerpa1_results, nerpa2_results, approved_matches_yaml, plots_dir
    """
    parser = argparse.ArgumentParser(description="Override nerpa paths")
    parser.add_argument(
        "--nerpa1-results",
        type=Path,
        default=None,
        help="path to nerpa1 results",
    )
    parser.add_argument(
        "--nerpa2-results",
        type=Path,
        default=None,
        help="path to nerpa2 results",
    )
    parser.add_argument(
        "--approved-matches-yaml",
        type=Path,
        default=None,
        help="path to approved_matches.yaml",
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=None,
        help="path to output plots dir",
    )

    ns = parser.parse_args()
    return ns

if __name__ == '__main__':
    args = parse_args()
    nerpa_dir = Path(__file__).resolve().parent.parent.parent
    assert nerpa_dir.name.startswith('nerpa'), f'Wrong nerpa_dir: {nerpa_dir}'

    nerpa1_vs_nerpa2(nerpa_dir=nerpa_dir,
                     args=args)
