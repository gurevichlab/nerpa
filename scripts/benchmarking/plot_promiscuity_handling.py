from __future__ import annotations
from itertools import chain
from pathlib import Path
import os
from typing import Literal, Dict, Optional, Sequence, List, Tuple

import numpy as np
import pandas as pd
import polars as pl
from collections import defaultdict

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from sklearn.linear_model import LinearRegression

from scripts.benchmarking.plots_data_helper import PlotsDataHelper, NerpaReport
from scripts.benchmarking.score_correctness_per_bgc_length import (
    score_correctness_per_bgc_length,
    plot_compare_score_correctness,
)

def plot_promiscuity_handling(reports: List[NerpaReport],
                              data_helper,
                              ax: Axes) -> Axes:
    """
    Overlay multiple reports: grouped bars per num_nrp_classes.
    """
    # Collect data
    all_dfs = []
    for report in reports:
        df = data_helper.promiscuity_handling_stats(report)
        df['report'] = report.name
        all_dfs.append(df)

    merged = pd.concat(all_dfs, ignore_index=True)

    # Unique x values
    x_vals = sorted(merged["num_nrp_classes"].unique())
    n_reports = len(reports)
    bar_width = 0.8 / n_reports

    for i, report in enumerate(reports):
        subset = merged[merged["report"] == report.name]
        xs = subset["num_nrp_classes"]
        ys = subset["mean_fraction_identified"]
        # Shift bars for grouped style
        ax.bar(
            [x + i * bar_width - (0.8 - bar_width) / 2 for x in xs],
            ys,
            width=bar_width,
            alpha=0.8,
            label=report.name,
        )

    ax.set_xlabel("num_nrp_classes")
    ax.set_ylabel("average of (nrp_classes_identified - 1) / (num_nrp_classes - 1)")
    ax.set_xticks(x_vals)
    ax.set_xticklabels([str(x) for x in x_vals])
    ax.set_ylim(0, 1)
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    ax.legend()
    return ax


