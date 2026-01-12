from __future__ import annotations
from typing import Optional, List

import numpy as np

from matplotlib.axes import Axes

from src.benchmarking.plots_data_helper import PlotsDataHelper, NerpaReport
from src.benchmarking.plots_utils import highlight_point


def plot_precision_recall_curve(nerpa_reports: List[NerpaReport],
                                ax: Axes,
                                data_helper: PlotsDataHelper,
                                top_matches_per_bgc: Optional[int] = None):
    """
    Plot Precision-Recall curves for all reports.

    Starting from top-scoring pairs, gradually lower the score threshold
    and compute precision and recall at each threshold.
    """
    for report in nerpa_reports:
        pr_points = data_helper.compute_precision_recall_curve(
            report,
            top_matches_per_bgc=top_matches_per_bgc
        )
        precision_list, recall_list, score_list = zip(*[
            (point.precision, point.recall, point.score)
            for point in pr_points
        ])

        # Plot precision vs recall
        ax.plot(
            recall_list,
            precision_list,
            label=report.name,
            marker='o',
            markersize=3,
            linewidth=2
        )

        # Highlight points at recall intervals of 0.1
        recall_thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        for recall_threshold in recall_thresholds:
            # Only highlight if the recall threshold is within the curve's range
            if min(recall_list) <= recall_threshold <= max(recall_list):
                # Interpolate to find precision and score at this recall value
                precision_at_threshold = np.interp(
                    recall_threshold,
                    recall_list,
                    precision_list
                )
                score_at_threshold = np.interp(
                    recall_threshold,
                    recall_list,
                    score_list
                )
                
                highlight_point(
                    ax,
                    x=recall_threshold,
                    y=precision_at_threshold,
                    label=f'{score_at_threshold:.2f}',
                    color='red'
                )

    ax.set_xlabel('Recall (TP / (TP + FN))', fontsize=14)
    ax.set_ylabel('Precision (TP / (TP + FP))', fontsize=14)
    ax.set_title(f'Precision-Recall Curve (top {top_matches_per_bgc} matches per BGC)', fontsize=16)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12)
