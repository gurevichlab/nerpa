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

from scripts.benchmarking.plot_promiscuity_handling import plot_promiscuity_handling
from scripts.benchmarking.plots_data_helper import PlotsDataHelper, NerpaReport
from scripts.benchmarking.score_correctness_per_bgc_length import (
    score_correctness_per_bgc_length,
    plot_compare_score_correctness,
)


def cnts_to_percentages(cnts: pd.Series[int]) -> pd.Series[float]:
    """
    Convert counts to percentages.
    """
    percentages = cnts[1:] / np.arange(1, len(cnts)) * 100.0  # skip the first value (0 identified at N=0)
    percentages.index = pd.RangeIndex(start=1,
                                      stop=len(percentages) + 1)  # adjust index to match the percentage calculation
    if not ((percentages > -1e-9) & (percentages < 100 + 1e-9)).all():
        raise ValueError("Percentages should be in the range [0, 100]")
    return percentages


class PlotsHelper:
    def __init__(self,
                 bgc_test_set: Literal['mibig4_wo_training_bgcs', 'training_bgcs'] = 'training_bgcs'):
        self.data_helper = PlotsDataHelper(bgc_test_set)

    def plot_num_correct_matches(self,
                                 nerpa_reports: List[NerpaReport],
                                 output_dir: Path,
                                 in_one_plot: bool = True,
                                 y_axis: Literal['Count', 'Percentage'] = 'Count') -> List[Path]:

        output_dir.mkdir(parents=True, exist_ok=True)
        graphs = {nerpa_report.name: self.data_helper.compute_num_correct_matches(nerpa_report)
                  for nerpa_report in nerpa_reports}

        def _plot_single_graph(ax: Axes,
                               _values: pd.Series,
                               title: str,
                               label: Optional[str] = None) -> None:
            if y_axis == 'Percentage':
                values = cnts_to_percentages(_values)
            else:
                values = _values
            if label:
                ax.plot(values.index, values, label=label)
            else:
                ax.plot(values.index, values)
            ax.set_xlabel('Number of matches')
            ax.set_ylabel(f'{y_axis} of correct matches')
            ax.set_title(title)
            ax.grid()
            ax.set_aspect('equal', adjustable='box')
            ax.set_ylim(bottom=0)

        if in_one_plot:
            out_file = output_dir / f'{y_axis}_correct_all.png'

            fig, ax = plt.subplots(figsize=(10, 6))
            for report_name, num_correct_array in graphs.items():
                _plot_single_graph(ax, num_correct_array, '', label=report_name)

            ax.set_title(f'{y_axis} of Correct Matches among the top N vs N')
            ax.legend()
            fig.tight_layout()
            fig.savefig(out_file)
            plt.close(fig)
            return [out_file]
        else:
            out_files = {report_name: output_dir / f'{y_axis}_correct_{report_name}.png'
                         for report_name in graphs.keys()}

            for report_name, num_correct_array in graphs.items():
                fig, ax = plt.subplots(figsize=(10, 6))
                _plot_single_graph(
                    ax,
                    num_correct_array,
                    f'{y_axis} of correct matches {report_name}'
                )
                fig.tight_layout()
                fig.savefig(out_files[report_name])
                plt.close(fig)

            return list(out_files.values())

    def plot_score_correctness(self,
                               nerpa_reports: List[NerpaReport],
                               output_dir: Path,
                               num_bins: int = 20) -> List[Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_files = {report.name: output_dir / f'score_correctness_{report.name}.png'
                     for report in nerpa_reports}

        for nerpa_report in nerpa_reports:
            correct_incorrect_counts = self.data_helper.compute_score_correctness(nerpa_report, num_bins)

            ax = correct_incorrect_counts[["Correct", "Incorrect"]].plot(
                kind="bar",
                stacked=True,
                figsize=(10, 4),
                width=0.9,
                rot=45,
                xlabel="Score range",
                ylabel="Number of matches",
                title=f"{nerpa_report.name}: Correct vs. incorrect matches per score bin (N={num_bins})"
            )
            ax.set_xticklabels(
                [f"{iv.left:.2f}–{iv.right:.2f}" for iv in correct_incorrect_counts.index],
                rotation=45, ha="right"
            )

            plt.tight_layout()
            plt.savefig(out_files[nerpa_report.name])
            plt.close()

        return list(out_files.values())

    def plot_score_correctness_per_bgc_len(self,
                                           nerpa_report: NerpaReport,
                                           output_dir: Path,
                                           num_len_bins: int = 5,
                                           ) -> Path:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / f'score_correctness_per_bgc_len_{nerpa_report.name}.png'
        data = score_correctness_per_bgc_length(nerpa_report,
                                                self.data_helper,
                                                num_len_bins)
        fig, ax = plt.subplots(figsize=(10, 6))
        plot_compare_score_correctness(data, ax=ax)
        ax.set_title(f'Score Correctness per BGC Length: {nerpa_report.name}')
        fig.tight_layout()
        fig.savefig(out_file)
        plt.close(fig)
        return out_file

    def box_plots_score_per_bgc_len(self,
                                    nerpa_report: NerpaReport,
                                    output_dir: Path,
                                    num_len_bins: int = 5):
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / f'score_boxplot_per_bgc_len_{nerpa_report.name}.png'
        grouped = self.data_helper.group_by_bgc_length(nerpa_report, num_len_bins)

        correct_scores_per_len = {
            len_interval: (
                pl.from_pandas(group)
                .filter(pl.col('is_correct'))
                .get_column(NerpaReport.SCORE_COL).to_list()
            )
            for len_interval, group in grouped.items()
        }

        wrong_scores_per_len = {
            len_interval: (
                pl.from_pandas(group)
                .filter(~pl.col('is_correct'))
                .get_column(NerpaReport.SCORE_COL).to_list()
            )
            for len_interval, group in grouped.items()
        }

        scores_to_plot = list(chain.from_iterable(
            (correct_scores_per_len[interval], wrong_scores_per_len[interval])
            for interval in grouped.keys()
        ))
        labels = list(chain.from_iterable(
            (f'{int(interval.left)}-{int(interval.right)}\nCorrect',
             f'{int(interval.left)}-{int(interval.right)}\nWrong')
            for interval in grouped.keys()
        ))

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.boxplot(scores_to_plot,
                   labels=labels,
                   showfliers=False)
        ax.set_title(f'Score Distribution per BGC Length: {nerpa_report.name}')
        ax.set_xlabel('BGC Length')
        ax.set_ylabel('Score')
        fig.tight_layout()
        fig.savefig(out_file)
        plt.close(fig)
        return out_file

    def bgc_len_histogram(self,
                          nerpa_report: NerpaReport,
                          output_dir: Path,
                          num_len_bins: int = 5) -> Path:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / f'bgc_length_histogram_{nerpa_report.name}.png'

        bgc_len = {
            row['bgc_id']: row["num_a_domains"]
            for _, row in self.data_helper.mibig_bgcs_info.iterrows()
        }

        df = pl.from_pandas(nerpa_report)
        df = df.with_columns(
            pl.col(NerpaReport.BGC_ID_COL)
            .map_elements(lambda bgc_id: bgc_len[bgc_id], return_dtype=pl.Int64)
            .alias("bgc_len")
        )

        print(df.head())

        bgc_lengths = list(
            df.group_by(NerpaReport.BGC_ID_COL)
            .agg(pl.col("bgc_len").first())["bgc_len"]
        )

        print(bgc_lengths[:10])
        fig, ax = plt.subplots(figsize=(10, 6))
        min_len, max_len = min(bgc_lengths), max(bgc_lengths)
        bins = np.arange(min_len, max_len + 2) - 0.5  # shift by 0.5 so bars center on integers

        ax.hist(bgc_lengths,
                bins=bins,
                edgecolor='black')
        ax.set_xticks(range(min_len, max_len + 1))

        ax.set_title(f'BGC Length Histogram: {nerpa_report.name}')
        ax.set_xlabel('BGC Length')
        ax.set_ylabel('Count')

        fig.tight_layout()
        fig.savefig(out_file)
        plt.close(fig)
        return out_file





    def plot_num_identified(self,
                            nerpa_reports: List[NerpaReport],
                            id_column: str,
                            output_dir: Path,
                            y_axis: Literal['Count', 'Percentage'] = 'Count',
                            top_ks: Tuple[int, ...] = (1,)) -> List[Path]:
        assert id_column in (NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL), \
            f'id_column should be one of {NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL}'

        output_dir.mkdir(parents=True, exist_ok=True)
        out_files_per_top_k = {top_k: output_dir / f'{y_axis}_identified_{id_column}_top{top_k}.png'
                               for top_k in top_ks}
        out_files_per_id = {report.name: output_dir / f'{y_axis}_identified_{id_column}_{report.name}.png'
                            for report in nerpa_reports}

        graphs = {
            report.name: {
                top_k: self.data_helper.compute_num_identified(
                    report,
                    id_column,
                    top_k
                )
                for top_k in top_ks
            }
            for report in nerpa_reports
        }

        # Plot per report (all top_ks)
        for report_name, topk_results in graphs.items():
            fig, ax = plt.subplots(figsize=(10, 6))
            for top_k, _values in topk_results.items():
                if y_axis == 'Percentage':
                    values = cnts_to_percentages(_values)
                else:
                    values = _values
                ax.plot(values.index, values.values, label=f"top-{top_k}")
            ax.set_title(f"Num Identified by {id_column}: {report_name}")
            ax.set_xlabel(f"Num top {id_column}")
            ax.set_ylabel(f"{y_axis} identified")
            ax.grid()
            ax.legend()
            fig.tight_layout()
            fig.savefig(out_files_per_id[report_name])
            plt.close(fig)

        # Plot per top_k (all reports)
        for top_k in top_ks:
            fig, ax = plt.subplots(figsize=(10, 6))
            for report_name, topk_results in graphs.items():
                _values = topk_results[top_k]
                if y_axis == 'Percentage':
                    values = cnts_to_percentages(_values)
                else:
                    values = _values
                ax.plot(values.index, values.values, label=f'{report_name}')
            ax.set_title(f"{y_axis} Identified {id_column}: top-{top_k}")
            ax.set_xlabel("Identifier")
            ax.set_ylabel(y_axis)
            ax.grid()
            ax.legend()
            fig.tight_layout()
            fig.savefig(out_files_per_top_k[top_k])
            plt.close(fig)

        return list(chain(out_files_per_top_k.values(), out_files_per_id.values()))

    def plot_total_identified(self,
                              nerpa_reports: List[NerpaReport],
                              id_column: str,
                              output_dir: Path,
                              max_top_k: int = 10,
                              y_axis: Literal['Count', 'Percentage'] = 'Count') -> List[Path]:
        """
        Plot the total number of identified BGCs/NRPs for each report
        """
        assert id_column in (NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL), \
            f'id_column should be one of {NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL}'

        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / f'total_{y_axis}_identified_{id_column}.png'
        graphs = {report.name: self.data_helper.compute_total_identified(
            report,
            id_column,
            max_top_k=max_top_k,
            y_axis=y_axis
        ) for report in nerpa_reports}

        fig, ax = plt.subplots(figsize=(10, 6))
        for report_name, total_identified in graphs.items():
            ax.step(
                total_identified.index,
                total_identified.values,
                where="post",  # "mid" makes the step flat across each x interval
                label=report_name
            )
            ax.set_xticks(total_identified.index)

        ax.set_title(f'{y_axis} of identified {id_column}',
                     fontsize=22)
        ax.set_xlabel(f'Num top matches considered for each {id_column}',
                      fontsize=18)
        ax.set_ylabel(f'{y_axis} Identified', fontsize=18)
        ax.grid()
        ax.legend(fontsize=18)
        fig.tight_layout()
        fig.savefig(out_file)
        plt.close(fig)

        return [out_file]

    def plot_promiscuity_handling(self,
                                  reports: List[NerpaReport],
                                  output_dir: Path) -> Path:
        """
        Create a grouped bar chart for all reports and save it to a PNG.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / f"promiscuity_hist.png"

        fig, ax = plt.subplots(figsize=(10, 6))
        plot_promiscuity_handling(reports, self.data_helper, ax)
        ax.set_title("Mean fraction of identified variants by num_nrp_classes")
        fig.tight_layout()
        fig.savefig(out_file, dpi=200)
        plt.close(fig)
        return out_file

    def plot_all(self,
                 nerpa_reports: List[NerpaReport],
                 output_dir: Path,
                 num_score_bins: int = 20,
                 top_ks: Tuple[int, ...] = (1, 3, 5, 10)) -> List[Path]:
        """
        Generate all plots and return paths to the generated files.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        plot_files = []

        plot_files.extend(self.plot_score_correctness(nerpa_reports,
                                                      output_dir,
                                                      num_score_bins))
        plot_files.extend(self.plot_score_correctness_per_bgc_len(nerpa_report, output_dir)
                          for nerpa_report in nerpa_reports)
        plot_files.append(self.plot_promiscuity_handling(nerpa_reports, output_dir))
        plot_files.extend(self.box_plots_score_per_bgc_len(nerpa_report, output_dir, num_len_bins=20)
                          for nerpa_report in nerpa_reports)
        plot_files.extend(self.bgc_len_histogram(nerpa_report, output_dir, num_len_bins=20)
                          for nerpa_report in nerpa_reports)

        #return plot_files
        for y_axis in ['Count', 'Percentage']:
            plot_files.extend(self.plot_num_correct_matches(nerpa_reports,
                                                            y_axis=y_axis,
                                                            output_dir=output_dir))
            for id_column in (NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL):
                plot_files.extend(self.plot_num_identified(nerpa_reports,
                                                           id_column=id_column,
                                                           y_axis=y_axis,
                                                           output_dir=output_dir,
                                                           top_ks=top_ks))
                plot_files.extend(self.plot_total_identified(nerpa_reports,
                                                             id_column=id_column,
                                                             y_axis=y_axis,
                                                             output_dir=output_dir))


        return plot_files
