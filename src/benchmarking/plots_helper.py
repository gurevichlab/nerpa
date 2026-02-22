from __future__ import annotations
from itertools import chain
from pathlib import Path
from typing import Literal, Optional, List, Tuple

import numpy as np
import pandas as pd
import polars as pl

from matplotlib import pyplot as plt

from src.benchmarking.data_frames import PNRPDB_Compound_Similarity
from src.benchmarking.data_helper_external_methods import interval_str_to_tuple
from src.benchmarking.plot_num_correct import _plot_num_correct_matches
from src.benchmarking.plot_precision_recall import plot_precision_recall_curve
from src.benchmarking.plot_promiscuity_handling import plot_promiscuity_handling
from src.benchmarking.plots_data_helper import PlotsDataHelper, NerpaReport
from src.benchmarking.score_correctness_per_bgc_length import (
    score_correctness_per_bgc_length,
    plot_compare_score_correctness,
)

def cnts_to_percentages(cnts: pl.Series) -> pl.Series:
        """
        Convert counts to percentages.
        """
        if cnts.len() <= 1:
            return pl.Series(name="percentages", values=[], dtype=pl.Float64)

        percentages = cnts[1:] / pl.Series(values=np.arange(1, cnts.len())) * 100.0  # skip the first value (0 identified at N=0)

        if not ((percentages > -1e-9) & (percentages < 100 + 1e-9)).all():
            raise ValueError("Percentages should be in the range [0, 100]")
        return percentages


class PlotsHelper:
    def __init__(self,
                 bgc_test_set: Literal['mibig4_wo_training_bgcs', 'training_bgcs'] = 'training_bgcs'):
        self.data_helper = PlotsDataHelper(bgc_test_set)
        self.dpi = 300
        self.height_px = 1000
        self.width_px = 1600
        self.axis_fontsize = 14
        self.title_fontsize = 15
        self.legend_fontsize = 14

    def plot_num_correct_matches(self,
                                 nerpa_reports: List[NerpaReport],
                                 output_dir: Path,
                                 y_axis: Literal['Count', 'Percentage'] = 'Count',
                                 top_matches_per_bgc_vals: Tuple[int, ...] = (1,3,5,10),
                                 pct_thresholds: Tuple[int, ...] = (50, ),
                                 max_num_matches: Optional[int] = 500) -> List[Path]:

        output_dir.mkdir(parents=True, exist_ok=True)
        out_files = []
        for top_matches_per_bgc in top_matches_per_bgc_vals:
            fig, ax = plt.subplots(figsize=(10, 6))

            for report in nerpa_reports:
                _plot_num_correct_matches(ax=ax,
                                          _report=report,
                                          title='',
                                          label=report.name,
                                          y_axis=y_axis,
                                          top_matches_per_bgc=top_matches_per_bgc,
                                          max_num_matches=max_num_matches,
                                          pct_thresholds=pct_thresholds)

            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)

            out_file = output_dir / f'{y_axis}_correct_top_{top_matches_per_bgc}.png'
            out_files.append(out_file)
            ax.set_title(f'{y_axis} of Correct Matches among the top N vs N')
            ax.legend()
            fig.tight_layout()
            fig.savefig(out_file)
            plt.close(fig)


        return out_files

    def plot_score_correctness(self,
                               nerpa_reports: List[NerpaReport],
                               output_dir: Path,
                               num_bins: int = 20) -> List[Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_files = {report.name: output_dir / f'score_correctness_{report.name}.png'
                     for report in nerpa_reports}

        for nerpa_report in nerpa_reports:
            correct_incorrect_counts = self.data_helper.compute_score_correctness(nerpa_report, num_bins)

            intervals_list = correct_incorrect_counts['score_bin'].to_list()

            # Extract lower and upper bounds from interval objects
            intervals = [
                interval_str_to_tuple(iv)
                for iv in intervals_list
            ]

            # Convert Polars to Pandas for plotting
            df_pandas: pd.DataFrame = correct_incorrect_counts.to_pandas()

            # Set score_bin as index and select only the count columns
            df_pandas = df_pandas.set_index('score_bin')[["Correct", "Incorrect"]]

            ax = df_pandas.plot(
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
                [f"{iv[0]:.2f}–{iv[1]:.2f}" for iv in intervals],
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
                group
                .filter(pl.col(NerpaReport.IS_CORRECT))
                .get_column(NerpaReport.SCORE).to_list()
            )
            for len_interval, group in grouped.items()
        }

        wrong_scores_per_len = {
            len_interval: (
                group
                .filter(~pl.col(NerpaReport.IS_CORRECT))
                .get_column(NerpaReport.SCORE).to_list()
            )
            for len_interval, group in grouped.items()
        }

        scores_to_plot = list(chain.from_iterable(
            (correct_scores_per_len[interval], wrong_scores_per_len[interval])
            for interval in grouped.keys()
        ))
        labels = list(chain.from_iterable(
            (f'{int(interval[0])}-{int(interval[1])}\nCorrect',
             f'{int(interval[0])}-{int(interval[1])}\nWrong')
            for interval in grouped.keys()
        ))

        fig, ax = plt.subplots(figsize=(10, 6))
        boxplots = ax.boxplot(scores_to_plot,
                              labels=labels,
                              showfliers=False,
                              patch_artist=True,  # needed to color the boxes
                              )

        # Alternate colors for odd/even
        colors = [(0, 1, 0, 0.2),  # green, 20% opaque
                  (1, 0, 0, 0.2)]  # red, 20% opaque

        for i, box in enumerate(boxplots['boxes']):
            box.set_facecolor(colors[i % 2])  # odd/even coloring

        # You may also want to color the medians, whiskers, etc.
        for i, median in enumerate(boxplots['medians']):
            median.set_color("black")

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
            for row in self.data_helper.mibig_bgcs_info.iter_rows(named=True)
        }

        # build mapping DF
        bgc_len_df = pl.DataFrame({
            NerpaReport.BGC_ID: list(bgc_len.keys()),
            "bgc_len": list(bgc_len.values()),
        }).with_columns(pl.col("bgc_len").cast(pl.Int64))

        # attach lengths
        df = nerpa_report.join(bgc_len_df, on=NerpaReport.BGC_ID, how="left")

        bgc_lengths = list(
            df.group_by(NerpaReport.BGC_ID)
            .agg(pl.col("bgc_len").first())["bgc_len"]
        )

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
        assert id_column in (NerpaReport.BGC_ID, NerpaReport.NRP_ISO_CLASS), \
            f'id_column should be one of {NerpaReport.BGC_ID, NerpaReport.NRP_ISO_CLASS}'

        output_dir.mkdir(parents=True, exist_ok=True)
        out_files_per_top_k = {top_k: output_dir / f'{y_axis}_identified_{id_column}_top{top_k}.png'
                               for top_k in top_ks}
        out_files_per_id = {report.name: output_dir / f'{y_axis}_identified_{id_column}_{report.name}.png'
                            for report in nerpa_reports}

        num_identified_graphs = {
            report.name: {
                top_k: {
                    cmp_method:
                        self.data_helper.compute_num_identified(
                            report,
                            id_column,
                            cmp_method,
                            top_k
                        )
                    for cmp_method in [PNRPDB_Compound_Similarity.NERPA_EQUAL_ALLOW_UNK_CHR,]
                                       #PNRPDB_Compound_Similarity.NERPA_NO_MORE_ONE_SUB_ALLOW_UNK_CHR,]
                }
                for top_k in top_ks
            }
            for report in nerpa_reports
        }

        colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray']
        # Plot per report (all top_ks)
        for report_name, topk_results in num_identified_graphs.items():
            fig, ax = plt.subplots()
            for (top_k, cmp_methods), color in zip(topk_results.items(), colors):
                for cmp_method, _values in cmp_methods.items():
                    values = (cnts_to_percentages(_values)
                              if y_axis == 'Percentage'
                              else _values)
                    linestyle = ('-'
                                  if cmp_method == PNRPDB_Compound_Similarity.NERPA_EQUAL_ALLOW_UNK_CHR
                                  else '--')


                    xs = range(len(values))
                    ax.plot(xs, values,
                            label=f"top-{top_k}({cmp_method})",
                            color=color,
                            linestyle=linestyle)
            ax.set_title(f"Num Identified by {id_column}: {report_name}")
            ax.set_xlabel(f"Num top {id_column}")
            ax.set_ylabel(f"{y_axis} identified")
            ax.grid()
            ax.legend(fontsize=self.legend_fontsize)
            fig.tight_layout()
            fig.savefig(out_files_per_id[report_name])
            plt.close(fig)

        # Plot per top_k (all reports)
        for top_k in top_ks:
            fig, ax = plt.subplots()
            for (report_name, topk_results), color in zip(num_identified_graphs.items(), colors):
                for cmp_method, linestyle in [(PNRPDB_Compound_Similarity.NERPA_EQUAL_ALLOW_UNK_CHR, '-'),]:
                                              #(PNRPDB_Compound_Similarity.NERPA_NO_MORE_ONE_SUB_ALLOW_UNK_CHR, '--')]:
                    _values = topk_results[top_k][cmp_method]
                    values = (cnts_to_percentages(_values)
                              if y_axis == 'Percentage'
                              else _values)
                    xs = range(len(values))
                    ax.plot(xs, values,
                            label=f'{report_name}',
                            #label=f'{report_name}({cmp_method})',
                            color=color,
                            linestyle=linestyle)

            name_identified = "BGCs" if id_column == NerpaReport.BGC_ID else "NRP iso classes"
            ax.set_title(f"{y_axis} of true hits present among top {top_k} matches",
                         fontsize=self.title_fontsize)
            ax.set_xlabel(f"Num top {name_identified}", fontsize=self.axis_fontsize)
            ax.set_ylabel(y_axis, fontsize=self.axis_fontsize)
            ax.grid()
            ax.legend(fontsize=self.legend_fontsize,
                      loc='lower right')
            fig.set_figheight(self.height_px / self.dpi)
            fig.tight_layout()
            fig.savefig(out_files_per_top_k[top_k], dpi=self.dpi)
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
        assert id_column in (NerpaReport.BGC_ID, NerpaReport.NRP_ISO_CLASS), \
            f'id_column should be one of {NerpaReport.BGC_ID, NerpaReport.NRP_ISO_CLASS}'

        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / f'total_{y_axis}_identified_{id_column}.png'
        total_identified_graphs = {
            report.name: {
                cmp_method: self.data_helper.compute_total_identified(report,
                                                                      id_column,
                                                                      max_top_k=max_top_k,
                                                                      y_axis=y_axis,
                                                                      cmp_method=cmp_method)
                for cmp_method in [PNRPDB_Compound_Similarity.NERPA_EQUAL_ALLOW_UNK_CHR,]
                                   #PNRPDB_Compound_Similarity.NERPA_NO_MORE_ONE_SUB_ALLOW_UNK_CHR,]
            }
            for report in nerpa_reports
        }

        fig, ax = plt.subplots()
        colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray']
        for report_name, color in zip(total_identified_graphs.keys(), colors):
            for cmp_method, total_identified in total_identified_graphs[report_name].items():
                linestyle = ('-'
                              if cmp_method == PNRPDB_Compound_Similarity.NERPA_EQUAL_ALLOW_UNK_CHR
                              else '--')
                ax.step(
                    range(len(total_identified)),
                    total_identified,
                    where="post",  # "mid" makes the step flat across each x interval
                    label=report_name,
                    linestyle=linestyle,
                    color=color,
                )
                ax.set_xticks(range(len(total_identified)))

        name_identified = "BGCs" if id_column == NerpaReport.BGC_ID else "NRP iso classes"
        ax.set_title(f'{y_axis} of true hits present',
                     fontsize=self.title_fontsize)
        ax.set_xlabel(f'Num top matches considered for each {name_identified}',
                      fontsize=self.axis_fontsize)
        ax.set_ylabel(f'{y_axis}', fontsize=self.axis_fontsize)
        ax.set_ylim(bottom=0)

        ax.grid()
        ax.legend(fontsize=self.legend_fontsize,
                  loc='lower right')
        fig.set_figheight(self.height_px / self.dpi)
        fig.tight_layout()
        fig.savefig(out_file, dpi=self.dpi)
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
        fig.set_figheight(self.height_px / self.dpi)
        fig.savefig(out_file, dpi=self.dpi)
        plt.close(fig)
        return out_file

    def plot_precision_recall_curve(self,
                                    nerpa_reports: List[NerpaReport],
                                    output_dir: Path,
                                    top_matches_per_bgc_vals: Tuple[Optional[int], ...] = (None, 1, 3, 5, 10)) -> List[Path]:
        """
        Plot Precision-Recall curves for all reports.

        Starting from top-scoring pairs, gradually lower the score threshold
        and compute precision and recall at each threshold.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        out_files = []

        for top_matches_per_bgc in top_matches_per_bgc_vals:
            fig, ax = plt.subplots()

            plot_precision_recall_curve(nerpa_reports,
                                        ax,
                                        self.data_helper,
                                        top_matches_per_bgc=top_matches_per_bgc,
                                        axis_fontsize=self.axis_fontsize,
                                        title_fontsize=self.title_fontsize,
                                        legend_fontsize=self.legend_fontsize)

            out_file = output_dir / f'precision_recall_curve_top_{top_matches_per_bgc}.png'
            fig.set_figheight(self.height_px / self.dpi)
            fig.tight_layout()
            fig.savefig(out_file, dpi=self.dpi)
            plt.close(fig)
            out_files.append(out_file)

        return out_files

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
        plot_files.extend(self.plot_precision_recall_curve(nerpa_reports, output_dir))
        plot_files.extend(self.box_plots_score_per_bgc_len(nerpa_report, output_dir, num_len_bins=20)
                          for nerpa_report in nerpa_reports)
        plot_files.extend(self.bgc_len_histogram(nerpa_report, output_dir, num_len_bins=20)
                          for nerpa_report in nerpa_reports)

        #return plot_files
        for y_axis in ['Count', 'Percentage']:
            plot_files.extend(self.plot_num_correct_matches(nerpa_reports,
                                                            y_axis=y_axis,
                                                            output_dir=output_dir))
            for id_column in (NerpaReport.BGC_ID, NerpaReport.NRP_ISO_CLASS):
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
