from __future__ import annotations
from itertools import chain
from pathlib import Path
import os
from typing import Literal, Dict, Optional, Sequence, List, Tuple, TYPE_CHECKING

import numpy as np
import pandas as pd
import polars as pl
from collections import defaultdict

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from sklearn.linear_model import LinearRegression

from scripts.benchmarking.nerpa_report import (
    NerpaReport,
    NerpaReportRow,
    load_nerpa_report,
    is_bgc_or_norine_nrp, OutputSizeConfig, check_reports_discrepancy
)
if TYPE_CHECKING:
    from scripts.benchmarking.plots_data_helper import PlotsDataHelper


def compute_num_correct_matches(data_helper,
                                nerpa_report: NerpaReport) -> pd.Series[int]:
    """
    Returns an array where the i-th element is the number of correct matches
    among the first i matches in the nerpa_report, sorted by the specified column.
    """
    num_correct = (nerpa_report
                   .sort_values(by=NerpaReport.SCORE_COL, ascending=False)
                   .apply(data_helper.match_is_correct, axis=1)
                   .cumsum())
    # Prepend the “0 correct at N=0” point
    return pd.concat([pd.Series([0]), num_correct], ignore_index=True)


def compute_score_correctness(data_helper,
                              nerpa_report: NerpaReport,
                              num_bins: int = 20,
                              bin_edges: List[float] = None) -> pd.DataFrame:
    """
    Returns a summary table with counts of correct and incorrect matches per score bin.

    The result is a DataFrame indexed by score intervals, with columns:
    - 'Correct': number of correctly matched rows in each bin
    - 'Incorrect': number of incorrectly matched rows in each bin

    Example output:
                      Incorrect  Correct
    score_bin
    (0.0, 0.2]               3        1
    (0.2, 0.4]               2        4
    ...
    """
    df = nerpa_report.copy()

    # --- Bin scores ------------------------------------------------------------
    # Cut scores into equal-width bins; include lowest edge to capture min score
    if bin_edges is not None:
        df["score_bin"] = pd.cut(df[NerpaReport.SCORE_COL], bins=bin_edges, include_lowest=True)
    else:
        df["score_bin"] = pd.cut(df[NerpaReport.SCORE_COL], bins=num_bins, include_lowest=True)

    return (
        df.groupby(["score_bin", "is_correct"], observed=False)
        .size()
        .unstack(fill_value=0)  # columns: False / True
        .rename(columns={False: "Incorrect", True: "Correct"})
    )


def compute_num_identified(data_helper,
                           nerpa_report: NerpaReport,
                           id_column: str,
                           top_k: int = 1) -> pd.Series[int]:
    """
    1. Sorts all BGCs/NRPs by their best match score (uses sort_matches_by column).
    2. For each N, counts how many BGCs/NRPs have been identified among the first N.
       A BGC/NRP is considered identified if it has at least one correct match
       among the top_k matches for it.
    """
    assert id_column in (NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL), \
        f"id_column must be either '{NerpaReport.BGC_ID_COL}' or '{NerpaReport.NRP_ISO_CLASS_COL}'"

    # 1) Sort IDs by their maximum score (so that we count in descending‐score order)
    best_scores = nerpa_report.groupby(id_column)[NerpaReport.SCORE_COL].max()
    sorted_ids = best_scores.sort_values(ascending=False).index

    # 2) For each ID, check if it has at least one correct match among the top_k matches
    identified_flags = (
        nerpa_report
        .groupby(id_column)[
            nerpa_report.columns]  # explicitly use all columns so that stupid pandas does not drop id_column
        .apply(
            lambda g: (
                g.nlargest(top_k, NerpaReport.SCORE_COL)
                .apply(data_helper.match_is_correct, axis=1)
                .any()
            )
        )
    )

    # 3) cumulative sum
    assert set(sorted_ids) == set(identified_flags.index)
    cumulative = identified_flags.reindex(sorted_ids).cumsum()

    # 4) Prepend the “0 identified at N=0” point
    return pd.concat([pd.Series([0]), cumulative], ignore_index=True)


def compute_total_identified(data_helper: 'PlotsDataHelper',
                             nerpa_report: NerpaReport,
                             id_column: str,
                             max_top_k: int = 10,
                             y_axis: Literal['Count', 'Percentage'] = 'Count') -> pd.Series[int]:
    assert id_column in (NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL), \
        f"id_column must be either '{NerpaReport.BGC_ID_COL}' or '{NerpaReport.NRP_ISO_CLASS_COL}'"

    total_identified = []
    for top_k in range(1, max_top_k + 1):
        num_identified = data_helper.compute_num_identified(nerpa_report, id_column, top_k)
        total_ids = len(data_helper.test_bgcs) \
            if id_column == NerpaReport.BGC_ID_COL \
            else len(data_helper.nrp_classes_with_matches)

        value = 100 * num_identified.iloc[-1] / total_ids \
            if y_axis == 'Percentage' else num_identified.iloc[-1]
        total_identified.append(value)

    return pd.Series(total_identified, index=range(1, max_top_k + 1))


# TODO: now redundant. Use data_helper.test_bgcs and data_helper.nrp_classes_with_matches instead
def get_identified_ids(data_helper,
                       report: NerpaReport,
                       id_column: str) -> Sequence[str]:
    assert id_column in (NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL), \
        f"id_column must be either '{NerpaReport.BGC_ID_COL}' or '{NerpaReport.NRP_ISO_CLASS_COL}'"

    def id_identified(report: NerpaReport, id_value: str) -> bool:
        return (report[report[id_column] == id_value]
                .apply(data_helper.match_is_correct, axis=1)
                .any())

    return [id_value for id_value in set(report[id_column])
            if id_identified(report, id_value)]

def promiscuity_handling_stats(data_helper,
                               nerpa_report: NerpaReport,
                               identified_if: Literal['default', 'top_N_variants_matches'] = 'top_N_variants_matches') -> pd.DataFrame:
    report = pl.from_pandas(nerpa_report)

    # add a column 'num_nrp_classes' with the total number of NRP iso-classes per BGC
    report = report.with_columns(
        [
            pl.col(NerpaReport.BGC_ID_COL)
            .map_elements(lambda bgc_id: len(data_helper.bgc_to_nrp_iso_classes.get(bgc_id, [])))
            .alias('num_nrp_classes'),

            pl.col(NerpaReport.SCORE_COL)
            .rank(method="dense", descending=True)  # rank best scores as 1, ties share the same rank
            .over(NerpaReport.BGC_ID_COL)
            .alias("rank_in_bgc")
        ]
    )

    if identified_if == 'top_N_variants_matches':
        # only keep top N matches per BGC, where N is the number of NRP iso-classes for that BGC
        report = report.filter(
            pl.col('rank_in_bgc') <= pl.col('num_nrp_classes')
        )

    # Count identified iso-classes per BGC
    cnts = (
        report.group_by(NerpaReport.BGC_ID_COL)
        .agg([
            pl.col(NerpaReport.NRP_ISO_CLASS_COL)
            .filter(pl.col('is_correct') == True)
            .n_unique()
            .alias("num_identified_nrp_classes"),

            pl.col('num_nrp_classes').first().alias('num_nrp_classes')
        ])
    )

    fraction_variants_identified = (
        cnts.filter((pl.col('num_identified_nrp_classes') > 0) & (pl.col('num_nrp_classes') > 1))
        .with_columns(
            ((pl.col('num_identified_nrp_classes') - 1) / (pl.col('num_nrp_classes') - 1))
            .alias('fraction_variants_identified')
        )
    )

    result = (
        fraction_variants_identified
        .group_by('num_nrp_classes')
        .agg([
            pl.count().alias('count_bgcs'),
            pl.col('fraction_variants_identified').mean().alias('mean_fraction_identified'),
            pl.col('fraction_variants_identified').median().alias('median_fraction_identified'),
            pl.col('fraction_variants_identified').std().alias('std_fraction_identified'),
        ])
    )

    return result.sort('num_nrp_classes').to_pandas()


def group_by_bgc_length(data_helper: 'PlotsDataHelper',
                        report: NerpaReport,
                        num_bins: int = 5,) -> Dict[pd.IntervalIndex, pd.DataFrame]:
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

    # Split into groups
    out: dict[pd.IntervalIndex, pd.DataFrame] = {
        interval: df.loc[df["_len_bin"] == interval]
        for interval in len_bins
    }
    return out

