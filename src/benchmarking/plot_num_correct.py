from __future__ import annotations
from typing import Literal, Optional, Tuple

import polars as pl

from matplotlib.axes import Axes
from src.benchmarking.plots_data_helper import NerpaReport
from src.benchmarking.plots_utils import highlight_point


def _with_row_index(df: pl.DataFrame,
                    name: str = 'index',
                    offset: int = 0) -> pl.DataFrame:
    if hasattr(df, "with_row_index"):
        return df.with_row_index(name, offset=offset)
    else:
        return df.with_row_count(name, offset=offset)


def _plot_num_correct_matches(ax: Axes,
                              _report: NerpaReport,
                              y_axis: Literal['Count', 'Percentage'],
                              title: str,
                              label: Optional[str] = None,
                              top_matches_per_bgc: Optional[int] = 10,
                              max_num_matches: Optional[int] = 500,
                              pct_thresholds: Tuple[int, ...] = (50,), ):  # label scores for given accuracy thresholds
    # 1. Remove matches beyond top_matches_per_bgc
    #print(f"top matches per bgc: {top_matches_per_bgc}")
    #print(f"Plotting num correct matches for report: {_report.name}")
    #print(f"Total matches before filtering: {_report.height}")
    #print(f"Total number of BGCs: {_report[NerpaReport.BGC_ID].n_unique()}")
    #print(f'Number of unique scores: {_report[NerpaReport.SCORE].n_unique()}')
    report = _report
    if top_matches_per_bgc is not None:
        report = report.filter(
            pl.col(NerpaReport.MATCH_RANK_FOR_BGC) <= top_matches_per_bgc
        )

    #print(f"Total matches after filtering: {report.height}")
    #print(f'Number of unique scores: {report[NerpaReport.SCORE].n_unique()}')

    report = _with_row_index(report, offset=1)

    # 2. Compute cumulative num correct matches
    report = report.with_columns(
        pl.col(NerpaReport.IS_CORRECT).cast(pl.Int64)
        .cum_sum()
        .alias("_num_correct_matches")
    )
    report = report.with_columns(
        (pl.col("_num_correct_matches") / pl.col('index') * 100.0)
        .alias("_pct_correct_matches")
    )

    # 3. Define xs and ys to plot
    if y_axis == 'Percentage':
        ys = report["_pct_correct_matches"].to_list()
    else:
        ys = report["_num_correct_matches"].to_list()
    if max_num_matches is not None:
        ys = ys[:max_num_matches]

    xs = range(1, len(ys) + 1)

    # 4. Plot the graph
    if label:
        ax.plot(xs, ys, label=label)
    else:
        ax.plot(xs, ys)
    ax.set_xlabel('Number of matches')
    ax.set_ylabel(f'{y_axis} of correct matches')
    ax.set_title(title)
    ax.grid()
    ax.set_aspect('equal', adjustable='box')

    # 5. Add score threshold annotations
    for threshold in pct_thresholds:
        # Find the first x where score >= threshold
        filtered = report.filter(pl.col('_pct_correct_matches') >= threshold)
        last_row = next(reversed(list(filtered.iter_rows(named=True))),
                         None)
        if last_row is None:
            continue  # threshold not reached

        x_val = last_row['index']
        y_val = last_row['_num_correct_matches'] \
            if y_axis == 'Count' else last_row['_pct_correct_matches']
        highlight_point(ax, x_val, y_val, label=f'score={last_row[NerpaReport.SCORE]:.2f}')
