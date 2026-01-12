import polars as pl
import numpy as np
from matplotlib import pyplot as plt

from src.benchmarking.data_helper_external_methods import interval_str_to_tuple
from src.benchmarking.plots_data_helper import PlotsDataHelper, NerpaReport

def score_correctness_per_bgc_length(
    report: NerpaReport,
    data_helper: PlotsDataHelper,
    num_bins: int = 3,
) -> dict[str, pl.DataFrame]:
    """Compute score correctness grouped by BGC length."""
    grouped = data_helper.group_by_bgc_length(report, num_bins)

    # Compute for all data first to get consistent score bins
    all_correctness = data_helper.compute_score_correctness(report)

    # Get bin edges from the score_bin column
    score_bins = sorted(map(interval_str_to_tuple, all_correctness['score_bin']))

    # Extract numeric boundaries
    bin_edges = [score_bins[0][0]] + [score_bin[1] for score_bin in score_bins]

    # Compute per length bin using same score bins
    out: dict[str, pl.DataFrame] = {"all": all_correctness}
    for interval_key, group in grouped.items():
        interval_str = str(interval_key)
        out[interval_str] = data_helper.compute_score_correctness(group, bin_edges=bin_edges)

    return out


def plot_score_correctness_graph(
    df: pl.DataFrame,
    label: str | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot correctness proportion per score bin."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    x_labels = df['score_bin'].to_list()
    x = range(len(df))

    correct = df['Correct'].to_numpy()
    incorrect = df['Incorrect'].to_numpy()
    total = correct + incorrect
    correctness = np.where(total > 0, correct / total, np.nan)

    ax.plot(x, correctness, marker='o', label=label or "Score Correctness")
    ax.set_xticks(x)
    ax.set_xticklabels([str(lbl) for lbl in x_labels], rotation=60, ha="right")
    ax.set_xlabel("Score Bins")
    ax.set_ylabel("Proportion Correct")
    ax.set_ylim(0, 1)

    return ax


def plot_compare_score_correctness(datasets: dict[str, pl.DataFrame],
                                   title: str = 'Score Correctness Comparison',
                                   ax: plt.Axes = None) -> plt.Axes:
    """Plot multiple score correctness datasets on same axes."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 6))
    for name, df in datasets.items():
        plot_score_correctness_graph(df, label=name, ax=ax)
    ax.set_title(title)
    ax.legend()
    return ax


'''
def plot_score_correctness_hist(df: pd.DataFrame,
                                ax: plt.Axes) -> plt.Axes:
    """
    Plot a histogram of correct vs incorrect counts per score bin.

    Parameters
    ----------
    df : pd.DataFrame
        Output of `compute_score_correctness()`.
        Indexed by score_bin, with 'Correct' and 'Incorrect' columns.
    """
    # Convert score bins to string labels for plotting
    x_labels = df.index.astype(str)
    x = range(len(df))  # numeric positions for bins

    correct = df["Correct"]
    incorrect = df["Incorrect"]

    width = 0.4  # width of each bar

    #fig, ax = plt.subplots(figsize=(10, 6))

    # Bars for incorrect counts
    ax.bar([i - width/2 for i in x], incorrect, width=width, label="Incorrect", alpha=0.7)
    # Bars for correct counts
    ax.bar([i + width/2 for i in x], correct, width=width, label="Correct", alpha=0.7)

    # Labels and legend
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.set_xlabel("Score Bins")
    ax.set_ylabel("Count")
    ax.set_title("Correct vs Incorrect Counts per Score Bin")
    ax.legend()
    return ax
'''
