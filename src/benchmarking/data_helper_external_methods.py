from __future__ import annotations

from math import ceil, floor
from typing import Literal, Sequence, TYPE_CHECKING, NamedTuple, List, Optional

import polars as pl

from src.benchmarking.data_frames import MIBiG_BGCs_Info
from src.benchmarking.nerpa_report import NerpaReport

if TYPE_CHECKING:
    from src.benchmarking.plots_data_helper import PlotsDataHelper

def interval_str_to_tuple(interval_str: str) -> tuple[float, float]:
    """Convert interval string to tuple of floats."""
    parts = interval_str.strip('()[]{}"').split(',')
    left = float(parts[0])
    right = float(parts[1])
    return (left, right)


def compute_num_correct_matches(data_helper: 'PlotsDataHelper',
                                nerpa_report: NerpaReport) -> pl.Series:
    """Returns cumulative count of correct matches."""
    sorted_report = nerpa_report.sort(NerpaReport.SCORE, descending=True)
    cumsum = sorted_report[NerpaReport.IS_CORRECT].cast(pl.Int64).cum_sum()

    # Prepend 0 for N=0 point
    return pl.Series([0]).extend(cumsum)


def compute_score_correctness(data_helper: 'PlotsDataHelper',
                              nerpa_report: NerpaReport,
                              num_bins: int = 20,
                              bin_edges: list[float] = None) -> pl.DataFrame:
    """Returns counts of correct and incorrect matches per score bin."""
    if bin_edges is not None:
        df = nerpa_report.with_columns(
            pl.col(NerpaReport.SCORE)
            .cut(breaks=bin_edges, left_closed=True)
            .alias('score_bin')
        )
    else:
        # Create equal-width bins
        min_score = nerpa_report[NerpaReport.SCORE].min()
        max_score = nerpa_report[NerpaReport.SCORE].max()
        step = (max_score - min_score) / num_bins
        edges = [min_score + i * step for i in range(num_bins + 2)]  # extra empty bin at the end to avoid inf
        df = nerpa_report.with_columns(
            pl.col(NerpaReport.SCORE).cut(breaks=edges, left_closed=True).alias('score_bin')
        )

    # Group by score bin and correctness
    result = (
        df.group_by(['score_bin', NerpaReport.IS_CORRECT])
        .agg(pl.count().alias('count'))
        .pivot(index='score_bin', columns=NerpaReport.IS_CORRECT, values='count')
        .fill_null(0)
    )
    # Rename columns if they exist
    if 'false' in result.columns:
        result = result.rename({'false': 'Incorrect'})
    if 'true' in result.columns:
        result = result.rename({'true': 'Correct'})

    # Ensure both columns exist (add missing ones with 0s)
    if 'Correct' not in result.columns:
        result = result.with_columns(pl.lit(0).alias('Correct'))
    if 'Incorrect' not in result.columns:
        result = result.with_columns(pl.lit(0).alias('Incorrect'))

    return result


def compute_num_identified(data_helper: 'PlotsDataHelper',
                           nerpa_report: NerpaReport,
                           id_column: str,
                           top_k: int = 1) -> pl.Series:
    """Count cumulative identified BGCs/NRPs."""
    assert id_column in (NerpaReport.BGC_ID, NerpaReport.NRP_ISO_CLASS), \
        f"Wrong id_column: {id_column}. Must be either BGC_ID or NRP_ID"

    # Compute table (ID -> best score)
    best_scores = (
        nerpa_report
        .group_by(id_column)
        .agg(
            pl.col(NerpaReport.SCORE)
            .max()
            .alias('max_score')
        )
        .sort('max_score', descending=True)
    )

    # Compute table (ID -> identified) (if any of top_k matches is correct)
    identified = (
        nerpa_report
        .sort(NerpaReport.SCORE, descending=True)
        .group_by(id_column)
        .agg(
            pl.col(NerpaReport.IS_CORRECT)
            .head(top_k)
            .any()
            .alias('identified')
        )
    )

    # Join and compute cumsum
    result = (
        best_scores
        .join(identified, on=id_column, validate="1:1", maintain_order="left")  # ids are sorted by max_score
        .select(
            pl.col('identified').cast(pl.Int64)
            .cum_sum()
        )
        .to_series()
    )

    # Prepend 0
    return pl.Series([0]).extend(result)


def compute_total_identified(data_helper: 'PlotsDataHelper',
                             nerpa_report: NerpaReport,
                             id_column: str,
                             max_top_k: int = 10,
                             y_axis: Literal['Count', 'Percentage'] = 'Count') -> pl.Series:
    """Compute total identified for various top_k values."""
    total_ids = (
        len(data_helper.test_bgcs)
        if id_column == NerpaReport.BGC_ID
        else len(data_helper.nrp_classes_with_matches)
    )

    results = []
    for top_k in range(1, max_top_k + 1):
        num_identified = compute_num_identified(data_helper, nerpa_report, id_column, top_k)
        value = num_identified[-1]
        if y_axis == 'Percentage':
            value = 100 * value / total_ids
        results.append(value)

    return pl.Series(results)


def get_identified_ids(data_helper: 'PlotsDataHelper',
                       report: NerpaReport,
                       id_column: str) -> Sequence[str]:
    """Get list of identified IDs."""
    id_to_identified = (
        report.group_by(id_column)
        .agg(
            pl.col(NerpaReport.IS_CORRECT)
            .any()
            .alias('identified')
        )
    )
    return (id_to_identified
            .filter(pl.col('identified'))
            .select(id_column)
            .to_series())


def promiscuity_handling_stats(data_helper: 'PlotsDataHelper',
                               nerpa_report: NerpaReport) -> pl.DataFrame:
    """Compute promiscuity handling statistics."""
    # Add num_nrp_classes and rank columns
    report = nerpa_report.with_columns([
        pl.col(NerpaReport.BGC_ID)
        .map_elements(lambda bgc_id: len(data_helper.bgc_to_nrp_iso_classes.get(bgc_id, [])), return_dtype=pl.Int64)
        .alias('num_nrp_classes'),

        pl.col(NerpaReport.SCORE)
        .rank(method="dense", descending=True)
        .over(NerpaReport.BGC_ID)
        .alias("match_rank_in_bgc")
    ])

    # Filter top N matches per BGC
    report = report.filter(
        pl.col('match_rank_in_bgc') <= pl.col('num_nrp_classes')
    )

    # Count identified iso-classes per BGC
    cnts = (
        report.group_by(NerpaReport.BGC_ID)
        .agg([
            pl.col(NerpaReport.NRP_ISO_CLASS)
            .filter(pl.col(NerpaReport.IS_CORRECT))
            .count()
            .alias("num_identified_nrp_classes"),

            pl.col('num_nrp_classes').first()
        ])
    )

    # Compute fraction of variants identified
    result = (
        cnts
        .filter((pl.col('num_identified_nrp_classes') > 0) & (pl.col('num_nrp_classes') > 1))
        .with_columns(
            ((pl.col('num_identified_nrp_classes') - 1) / (pl.col('num_nrp_classes') - 1))
            .alias('fraction_variants_identified')
        )

        .group_by('num_nrp_classes')
        .agg([
            pl.count().alias('count_bgcs'),
            pl.col('fraction_variants_identified').mean().alias('mean_fraction_identified'),
            pl.col('fraction_variants_identified').median().alias('median_fraction_identified'),
            pl.col('fraction_variants_identified').std().alias('std_fraction_identified'),
        ])
        .sort('num_nrp_classes')
    )

    return result


def group_by_bgc_length(data_helper: 'PlotsDataHelper',
                        report: NerpaReport,
                        num_bins: int = 5) -> dict[tuple, pl.DataFrame]:
    """Group report by BGC length bins."""
    # Build BGC length mapping
    bgc_len_map = {
        row[MIBiG_BGCs_Info.BGC_ID]: row[MIBiG_BGCs_Info.NUM_A_DOMAINS]
        for row in data_helper.mibig_bgcs_info.iter_rows(named=True)
    }

    # Add BGC length column
    df = report.with_columns(
        pl.col(NerpaReport.BGC_ID)
        .map_elements(lambda bgc_id: bgc_len_map.get(bgc_id), return_dtype=pl.Int64)
        .alias('_bgc_len')
    )

    # Create bins
    min_len = df['_bgc_len'].min()
    max_len = df['_bgc_len'].max()
    step = (max_len - min_len) / num_bins

    # Bin and group
    df = df.with_columns(
        pl.col('_bgc_len')
        .cut(breaks=[min_len + i * step for i in range(num_bins + 2)],  # extra empty bin at the end to avoid inf
             left_closed=True,)
        .alias('_len_bin')
    )

    def parse_interval(interval_str) -> tuple[int, int]:
        parts = interval_str.strip('()[]{}"').split(',')
        left = ceil(float(parts[0]))
        right = floor(float(parts[1]))
        return (left, right)

    # Create dictionary of groups
    groups = {}
    for interval_str in df['_len_bin'].unique().sort().cast(str):
        group = (df
                 .filter(pl.col('_len_bin') == interval_str)
                 .drop(['_bgc_len', '_len_bin']))
        groups[parse_interval(interval_str)] = group

    return groups

class PrecisionRecallPoint(NamedTuple):
    precision: float
    recall: float
    score: float

def compute_precision_recall_curve(data_helper: 'PlotsDataHelper',
                                   nerpa_report: NerpaReport,
                                   top_matches_per_bgc: Optional[int] = None) -> List[PrecisionRecallPoint]:
    """
    Compute precision and recall at various score thresholds.
    Assuming each BGC-NRP pair is scored, for each score threshold T
    we say the prediction is positive if score >= T, negative otherwise.
    Then we compute precision and recall based on these predictions:

        precision = TP / (TP + FP)
        recall = TP / (TP + FN)

    For pairs missing in the report we assume they have score -inf (i.e., negative prediction).
    Returns a list of PrecisionRecallPoint named tuples.

    normalize_total_tp: if BGC doesn't appear in the report,
    we ignore it when computing total number of true pairs.
    """
    # Sort by score descending

    _report = nerpa_report.sort(NerpaReport.SCORE, descending=True)
    if top_matches_per_bgc is not None:
        _report = _report.filter(
            pl.col(NerpaReport.MATCH_RANK_FOR_BGC) <= top_matches_per_bgc
        )

    # Total number of true pairs (all correct matches that exist)
    total_true_pairs = sum(len(data_helper.bgc_to_nrp_iso_classes[bgc_id])
                           for bgc_id in data_helper.test_bgcs)

    precision_recall_points = []
    for score_threshold in _report[NerpaReport.SCORE].unique():
        _report_positive = _report.filter(pl.col(NerpaReport.SCORE) >= score_threshold)
        num_positive = _report_positive.height
        num_true_positives = _report_positive.filter(pl.col(NerpaReport.IS_CORRECT)).height
        num_false_negatives = total_true_pairs - num_true_positives

        precision = num_true_positives / num_positive
        recall = num_true_positives / (num_true_positives + num_false_negatives)

        precision_recall_points.append(PrecisionRecallPoint(precision=precision,
                                                            recall=recall,
                                                            score=score_threshold))

    return precision_recall_points


def extra_false_positives(data_helper: 'PlotsDataHelper',
                          report1: NerpaReport,
                          report2: NerpaReport) -> pl.DataFrame:
    """
    Assess which matches are handled better in report1 than in report2.
    For each prefix length N, includes the N-th match from report1 if:
    1. The match is wrong (IS_CORRECT == False)
    2. The match is absent in the N-th prefix of report2 (i.e., either not in report2 at all,
       or its index in report2 is >= N)

    Returns a DataFrame with extra false positives, their positions in both reports,
    and their prefix positions.
    """
    # Sort both reports by score descending and add indexes
    sorted_report1 = report1.sort(NerpaReport.SCORE, descending=True).with_row_index(f'index_in_{report1.name}')
    sorted_report2 = report2.sort(NerpaReport.SCORE, descending=True).with_row_index(f'index_in_{report2.name}')

    extra_fps = []

    for row_report1 in (sorted_report1
                        .filter(pl.col(~NerpaReport.IS_CORRECT))
                        .iter_rows(named=True)):
        # Get the N-th match from report1 (0-indexed, so n-1)

        rows_report2 = sorted_report2.filter(
            (pl.col(NerpaReport.BGC_ID) == row_report1[NerpaReport.BGC_ID]) &
            (pl.col(NerpaReport.NRP_ISO_CLASS) == row_report1[NerpaReport.NRP_ISO_CLASS])
        )
        row_report2 = rows_report2[0] if rows_report2.height > 0 else None
        if (row_report2 is None or
                row_report2[f'index_in_{report2.name}'] > row_report1[f'index_in_{report1.name}']):
            extra_fps.append({
                'bgc_id': row_report1[NerpaReport.BGC_ID],
                'nrp_iso_class': row_report1[NerpaReport.NRP_ISO_CLASS],
                f'index_in_{report1.name}': row_report1[f'index_in_{report1.name}'],
                f'index_in_{report2.name}': row_report2[f'index_in_{report2.name}'] if row_report2 is not None else None,
                'score_in_' + report1.name: row_report1[NerpaReport.SCORE],
                'score_in_' + report2.name: row_report2[NerpaReport.SCORE] if row_report2 is not None else None,
            })

    return pl.DataFrame(extra_fps)



