from typing import Dict, Tuple

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from scripts.benchmarking.plots_data_helper import PlotsDataHelper, NerpaReport

def score_correctness_per_bgc_length(
    report: NerpaReport,
    data_helper: PlotsDataHelper,
    num_bins: int = 3,
) -> dict[str, pd.DataFrame]:
    # Build map from BGC ID to length directly (allowing duplicates in table)
    bgc_len = {
        row['bgc_id']: row["num_a_domains"]
        for _, row in data_helper.mibig_bgcs_info.iterrows()
    }

    df = report.copy()
    df["_bgc_len"] = df[NerpaReport.BGC_ID_COL].map(bgc_len)

    # Bin by BGC length
    df["_len_bin"] = pd.cut(df["_bgc_len"], bins=num_bins, include_lowest=True, right=True)
    len_bins = df["_len_bin"].cat.categories

    # Compute 'all' once to fix score bins
    all_correctness = data_helper.compute_score_correctness(df)
    score_idx = all_correctness.index  # CategoricalIndex or IntervalIndex
    intervals = getattr(score_idx, "categories", score_idx)  # IntervalIndex either way
    bin_edges = intervals.left.tolist() + [intervals.right[-1]]

    # Compute per length bin using the same score bins
    out: dict[str, pd.DataFrame] = {"all": all_correctness}
    for interval in len_bins:
        sub = df.loc[df["_len_bin"] == interval]
        out[str(interval)] = data_helper.compute_score_correctness(sub, bin_edges=bin_edges)

    return out


def plot_score_correctness_graph(
    df: pd.DataFrame,
    label: str | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Plot correctness proportion per score bin.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame returned by `compute_score_correctness`.
    label : str, optional
        Label for this dataset (used in legend).
    ax : plt.Axes, optional
        Axes to plot into. If None, a new figure and axes are created.

    Returns
    -------
    plt.Axes
        The axes with the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    x_labels = df.index.astype(str)
    x = range(len(df))

    total = df["Correct"] + df["Incorrect"]
    correctness = np.where(total > 0, df["Correct"] / total, np.nan)

    ax.plot(x, correctness, marker='o', label=label or "Score Correctness")
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation=60, ha="right")
    ax.set_xlabel("Score Bins")
    ax.set_ylabel("Proportion Correct")
    ax.set_ylim(0, 1)

    return ax


def plot_compare_score_correctness(datasets: dict[str, pd.DataFrame],
                                   title: str = 'Score Correctness Comparision',
                                   ax: plt.Axes = None) -> plt.Axes:
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
