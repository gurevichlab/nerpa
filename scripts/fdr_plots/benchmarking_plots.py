from __future__ import annotations
from itertools import chain
from pathlib import Path
import os
from typing import Literal, Dict, Optional, Sequence, List, Tuple

import numpy as np
import pandas as pd
from collections import defaultdict

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from sklearn.linear_model import LinearRegression

NerpaReport = pd.DataFrame
NerpaReportRaw = pd.Series

class PlotsDataHelper:
    mibig_nrps_table: pd.DataFrame
    rban_graphs_table: pd.DataFrame
    bgc_to_nrp_iso_classes: defaultdict[str, set[int]]

    def __init__(self):
        nerpa_dir = Path(os.getcwd()).parent.parent
        self.mibig_nrps_table = pd.read_csv(nerpa_dir / 'scripts' / 'build_mibig_info_table' / 'mibig_bgcs_info.tsv',
                                       sep='\t')
        self.rban_graphs_table = pd.read_csv(nerpa_dir / 'scripts' / 'build_mibig_info_table' / 'rban_graphs_filtered.tsv',
                                        sep='\t')
        self.mibig_norine_nrps = set(pd.read_csv(nerpa_dir / 'data' / 'mibig_norine_compounds.tsv', sep='\t')['ID'])

        self.bgc_to_nrp_iso_classes = defaultdict(set)
        for _, row in self.mibig_nrps_table.iterrows():
            self.bgc_to_nrp_iso_classes[row['bgc_id']].add(row['iso_class_idx'])

    def nrp_id_to_iso_class(self, nrp_id: str) -> int:
        return self.rban_graphs_table.loc[self.rban_graphs_table['compound_id'] == nrp_id, 'iso_class_idx'].values[0]

    def match_is_correct(self, report_row: NerpaReportRaw):  # row in report.tsv
        try:
            nrp_id = report_row['NRP_ID']
            bgc_id = report_row['Genome_ID']
            nrp_class = self.rban_graphs_table.loc[self.rban_graphs_table['compound_id'] == nrp_id, 'iso_class_idx'].values[0]
            return nrp_class in self.bgc_to_nrp_iso_classes[bgc_id]
        except KeyError:
            print('Invalid row: ', report_row)
            raise KeyError("Row does not contain valid NRP_ID or Genome_ID")

    def compute_num_correct_matches(self,
                                    nerpa_report: NerpaReport,
                                    sort_by: Literal['LogOdds score', 'p_value'] = 'LogOdds score') -> pd.Series[int]:
        """
        Returns an array where the i-th element is the number of correct matches
        among the first i matches in the nerpa_report, sorted by the specified column.
        """
        num_correct = (nerpa_report
                       .sort_values(by=sort_by, ascending=(sort_by == 'p_value'))
                       .apply(self.match_is_correct, axis=1)
                       .cumsum())
        # Prepend the “0 correct at N=0” point
        return pd.concat([pd.Series([0]), num_correct], ignore_index=True)

    def compute_score_correctness(self,
                                  nerpa_report: NerpaReport,
                                  score_column: Literal['LogOdds score', 'p_value'] = 'LogOdds score',
                                  num_bins: int = 20) -> pd.DataFrame:
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
        df["is_correct"] = df.apply(self.match_is_correct, axis=1)

        # --- Bin scores ------------------------------------------------------------
        # Cut scores into equal-width bins; include lowest edge to capture min score
        df["score_bin"] = pd.cut(df[score_column], bins=num_bins, include_lowest=True)
        return (
            df.groupby(["score_bin", "is_correct"], observed=False)
            .size()
            .unstack(fill_value=0)               # columns: False / True
            .rename(columns={False: "Incorrect", True: "Correct"})
        )

    def compute_num_identified(self,
                               nerpa_report: NerpaReport,
                               id_column: Literal['Genome_ID', 'NRP_ID'],
                               score_column: Literal['LogOdds score', 'p_value'] = 'LogOdds score',
                               top_k: int = 1) -> pd.Series[int]:
        """
        1. Sorts all BGCs/NRPs by their best match score (uses sort_matches_by column).
        2. For each N, counts how many BGCs/NRPs have been identified among the first N.
           A BGC/NRP is considered identified if it has at least one correct match
           among the top_k matches for it.
        """
        if score_column == 'p_value':
            nerpa_report['1-p_value'] = 1.0 - nerpa_report['p_value']
            score_column = '1-p_value'

        if id_column == 'NRP_ID':  # TODO: refactor upstream to avoid this hack
            id_column = 'iso_class_idx'
            
        # 1) Sort IDs by their maximum score (so that we count in descending‐score order)
        best_scores = nerpa_report.groupby(id_column)[score_column].max()
        sorted_ids = best_scores.sort_values(ascending=False).index

        # 2) For each ID, check if it has at least one correct match among the top_k matches
        identified_flags = (
            nerpa_report
            .groupby(id_column)[nerpa_report.columns]  # explicitly use all columns so that stupid pandas does not drop id_column
            .apply(
                lambda g: (
                    g.nlargest(top_k, score_column)
                    .apply(self.match_is_correct, axis=1)
                    .any()
                )
            )
        )

        # 3) cumulative sum
        assert set(sorted_ids) == set(identified_flags.index)
        cumulative = identified_flags.reindex(sorted_ids).cumsum()

        # 4) Prepend the “0 identified at N=0” point
        return pd.concat([pd.Series([0]), cumulative], ignore_index=True)

    def compute_total_identified(self,
                                 nerpa_report: NerpaReport,
                                 id_column: Literal['Genome_ID', 'NRP_ID'],
                                 score_column: Literal['LogOdds score', 'p_value'] = 'LogOdds score',
                                 max_top_k: int = 10,
                                 y_axis: Literal['Count', 'Percentage'] = 'Count') -> pd.Series[int]:
        total_identified = []
        for top_k in range(1, max_top_k + 1):
            num_identified = self.compute_num_identified(nerpa_report, id_column, score_column, top_k)
            value = 100 * num_identified.iloc[-1] / len(num_identified) \
                if y_axis == 'Percentage' else num_identified.iloc[-1]
            total_identified.append(value)

        return pd.Series(total_identified, index=range(1, max_top_k + 1))

    def get_identified_ids(self,
                           report: NerpaReport,
                           id_column: Literal['Genome_ID', 'NRP_ID']) -> Sequence[str]:
        def id_identified(report: NerpaReport, id_value: str) -> bool:
            return (report[report[id_column] == id_value]
                    .apply(self.match_is_correct, axis=1)
                    .any())

        return [id_value for id_value in set(report[id_column])
                if id_identified(report, id_value)]

    def leave_only_data_for_benchmarking(self,
                                         _report: NerpaReport,
                                         report_name: str = 'REPORT_NAME_MISSING',
                                         score_column: str = 'LogOdds score',
                                         top_matches_per_bgc: int = 10) -> NerpaReport:
        report = _report.copy()
        mibig_nrps_table = self.mibig_nrps_table
        bgc_ids_appr = mibig_nrps_table[mibig_nrps_table['in_approved_matches'] == True]['bgc_id'].unique()

        # Filter report to only include BGCs present in the approved matches
        report = report[report['Genome_ID'].isin(bgc_ids_appr)]
        if set(report['Genome_ID']) != set(bgc_ids_appr):
            print(f"Warning! Some BGC IDs from the approved matches are missing in the report {report_name}")
            print(f"Missing BGC IDs: {set(bgc_ids_appr) - set(report['Genome_ID'])}")

        # Filter report to only include NRPs from MIBiG + Norine
        report = report[report['NRP_ID'].map(lambda nrp_id: nrp_id in self.mibig_norine_nrps)]

        # Leave only one of isomorphic NRPs per BGC
        report['iso_class_idx'] = report['NRP_ID'].map(self.nrp_id_to_iso_class)
        report.drop_duplicates(subset=['Genome_ID', 'iso_class_idx'], inplace=True)

        '''
        mibig_norine_iso_classes = set(map(self.nrp_id_to_iso_class, self.mibig_norine_nrps))
        if set(report['iso_class_idx']) != mibig_norine_iso_classes:
            print(f"Warning! Some NRP iso-classes from MIBiG + Norine are missing in the report {report_name}")
            print(f"Missing iso classes: {mibig_norine_iso_classes - set(report['iso_class_idx'])}")
        '''

        # Keep only top matches per BGC
        report = (report
                    .sort_values(by=score_column, ascending=False)
                    .groupby('Genome_ID')
                    .head(top_matches_per_bgc)
                    .reset_index(drop=True))

        return report

    def load_nerpa_report(self,
                          report_path: Path,
                          report_name: str = 'REPORT_NAME_MISSING',
                          score_column: str = 'LogOdds score',
                          nerpa_version: Literal['Nerpa 1', 'Nerpa 2'] = 'Nerpa 2') -> NerpaReport:
        report = pd.read_csv(report_path, sep='\t' if nerpa_version == 'Nerpa 2' else ',')
        if nerpa_version == 'Nerpa 1':
            report = nerpa1_report_to_nerpa2_compatible(report)

        report['LogOdds score'] = report[score_column]
        report = self.leave_only_data_for_benchmarking(report, report_name)
        return report

    def check_report_discrepancy(self,
                                 report1: NerpaReport,
                                 report2: NerpaReport,
                                 report1_name: str,
                                 report2_name: str,
                                 check_nrp_ids: bool = True) -> None:
        """
        Check for discrepancies between two reports.
        """
        print('Checking report discrepancies between '
              f'{report1_name} and {report2_name}...')

        if set(report1['Genome_ID']) != set(report2['Genome_ID']):
            print(f"Warning! Genome IDs in reports "
                  f"{report1_name} and {report2_name} do not match")
            print(f'{report1_name} \ {report2_name} = {set(report1["Genome_ID"]) - set(report2["Genome_ID"])}')
            print(f'{report2_name} \ {report1_name} = {set(report2["Genome_ID"]) - set(report1["Genome_ID"])}')
        else:
            print("Genome IDs sets match")

        if check_nrp_ids:
            if set(report1['NRP_ID']) != set(report2['NRP_ID']):
                print(f"Warning! NRP IDs in reports "
                      f"{report1_name} and {report2_name} do not match")
                print(f'{report1_name} \ {report2_name} = {set(report1["NRP_ID"]) - set(report2["NRP_ID"])}')
                print(f'{report2_name} \ {report1_name} = {set(report2["NRP_ID"]) - set(report1["NRP_ID"])}')
            else:
                print("NRP IDs sets match")

        id_columns = ['Genome_ID', 'NRP_ID'] if check_nrp_ids else ['Genome_ID']
        for id_column in id_columns:
            identified_ids_in_report1 = self.get_identified_ids(report1, id_column)
            identified_ids_in_report2 = self.get_identified_ids(report2, id_column)
            if set(identified_ids_in_report1) != set(identified_ids_in_report2):
                print(f"Warning! Identified {id_column}s in reports "
                      f"{report1_name} and {report2_name} do not match")
                print(f'Identified in {report1_name} but not in {report2_name}:'
                      f'\n{set(identified_ids_in_report1) - set(identified_ids_in_report2)}')
                print(f'Identified in {report2_name} but not in {report1_name}:'
                      f'\n{set(identified_ids_in_report2) - set(identified_ids_in_report1)}')
            else:
                print(f'Identified {id_column}s match')


def cnts_to_percentages(cnts: pd.Series[int]) -> pd.Series[float]:
    """
    Convert counts to percentages.
    """
    percentages = cnts[1:] / np.arange(1, len(cnts)) * 100.0  # skip the first value (0 identified at N=0)
    percentages.index = pd.RangeIndex(start=1, stop=len(percentages) + 1)  # adjust index to match the percentage calculation
    if not ((percentages > -1e-9) & (percentages < 100 + 1e-9)).all():
        raise ValueError("Percentages should be in the range [0, 100]")
    return percentages

class PlotsHelper:
    def __init__(self):
        self.data_helper = PlotsDataHelper()

    def plot_num_correct_matches(self,
                                 nerpa_reports: Dict[str, NerpaReport],
                                 output_dir: Path,
                                 in_one_plot: bool = True,
                                 score_columns: Optional[Dict[str, Literal['LogOdds score', 'p_value']]] = None,
                                 y_axis: Literal['Count', 'Percentage'] = 'Count') -> List[Path]:
        if score_columns is None:
            score_columns = {report_name: 'LogOdds score' for report_name in nerpa_reports.keys()}

        output_dir.mkdir(parents=True, exist_ok=True)
        graphs = {report_name: self.data_helper.compute_num_correct_matches(nerpa_report,
                                                                            sort_by=score_columns[report_name])
                  for report_name, nerpa_report in nerpa_reports.items()}

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
                               nerpa_reports: Dict[str, NerpaReport],
                               score_columns: Dict[str, Literal['LogOdds score', 'p_value']],
                               output_dir: Path,
                               num_bins: int = 20) -> List[Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_files = {report_name: output_dir / f'score_correctness_{report_name}.png'
                     for report_name in nerpa_reports.keys()}

        for report_name, nerpa_report in nerpa_reports.items():
            correct_incorrect_counts = self.data_helper.compute_score_correctness(nerpa_report, 
                                                                                  score_columns[report_name],
                                                                                  num_bins)

            ax = correct_incorrect_counts[["Correct", "Incorrect"]].plot(
                kind="bar",
                stacked=True,
                figsize=(10, 4),
                width=0.9,
                rot=45,
                xlabel="Score range",
                ylabel="Number of matches",
                title=f"{report_name}: Correct vs. incorrect matches per score ({score_columns[report_name]}) bin (N={num_bins})"
            )
            ax.set_xticklabels(
                [f"{iv.left:.2f}–{iv.right:.2f}" for iv in correct_incorrect_counts.index],
                rotation=45, ha="right"
            )

            plt.tight_layout()
            plt.savefig(out_files[report_name])
            plt.close()
    
        return list(out_files.values())

    def plot_num_identified(self,
                            nerpa_reports: Dict[str, 'NerpaReport'],
                            id_column: Literal['Genome_ID', 'NRP_ID'],
                            output_dir: Path,
                            score_columns: Dict[str, Literal['LogOdds score', 'p_value']],
                            y_axis: Literal['Count', 'Percentage'] = 'Count',
                            top_ks: Tuple[int, ...] = (1,)) -> List[Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_files_per_top_k = {top_k: output_dir / f'{y_axis}_identified_{id_column}_top{top_k}.png'
                               for top_k in top_ks}
        out_files_per_id = {report_name: output_dir / f'{y_axis}_identified_{id_column}_{report_name}.png'
                            for report_name in nerpa_reports.keys()}

        graphs = {
            report_name: {
                top_k: self.data_helper.compute_num_identified(
                    nerpa_report,
                    id_column,
                    score_columns[report_name],
                    top_k
                )
                for top_k in top_ks
            }
            for report_name, nerpa_report in nerpa_reports.items()
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
            ax.set_title(f"Num Identified by {id_column} using {score_columns[report_name]}: {report_name}")
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
                ax.plot(values.index, values.values, label=f'{report_name} ({score_columns[report_name]})')
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
                              nerpa_reports: Dict[str, NerpaReport],
                              id_column: Literal['Genome_ID', 'NRP_ID'],
                              output_dir: Path,
                              score_columns: Dict[str, Literal['LogOdds score', 'p_value']],
                              max_top_k: int = 10,
                              y_axis: Literal['Count', 'Percentage'] = 'Count') -> List[Path]:
        """
        Plot the total number of identified BGCs/NRPs for each report
        """

        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / f'total_{y_axis}_identified_{id_column}.png'
        graphs = {report_name: self.data_helper.compute_total_identified(
            nerpa_report,
            id_column,
            score_columns[report_name],
            max_top_k=max_top_k,
            y_axis=y_axis
        ) for report_name, nerpa_report in nerpa_reports.items()}

        fig, ax = plt.subplots(figsize=(10, 6))
        for report_name, total_identified in graphs.items():
            ax.plot(total_identified.index,
                    total_identified.values,
                    label=f"{report_name}")
            ax.set_xticks(total_identified.index)

        ax.set_title(f'{y_axis} of identified {"BGCs" if id_column == "Genome_ID" else "NRPs"}',
                     fontsize=22)
        ax.set_xlabel(f'Num top matches considered for each {"BGC" if id_column == "Genome_ID" else "NRP"}',
                      fontsize=18)
        ax.set_ylabel(f'{y_axis} Identified', fontsize=18)
        ax.grid()
        ax.legend(fontsize=18)
        fig.tight_layout()
        fig.savefig(out_file)
        plt.close(fig)

        return [out_file]

    def plot_all(self,
                    nerpa_reports: Dict[str, NerpaReport],
                    output_dir: Path,
                    num_score_bins: int = 20,
                    top_ks: Tuple[int, ...] = (1, 3, 5, 10)) -> List[Path]:
            """
            Generate all plots and return paths to the generated files.
            """
            output_dir.mkdir(parents=True, exist_ok=True)
            score_columns = {report_name: 'LogOdds score' for report_name in nerpa_reports.keys()}
            plot_files = []

            plot_files.extend(self.plot_score_correctness(nerpa_reports,
                                                          score_columns,
                                                          output_dir,
                                                          num_score_bins))
            for y_axis in ['Count', 'Percentage']:
                plot_files.extend(self.plot_num_correct_matches(nerpa_reports,
                                                                y_axis=y_axis,
                                                                output_dir=output_dir))
                for id_column in ('Genome_ID', 'NRP_ID'):
                    plot_files.extend(self.plot_num_identified(nerpa_reports,
                                                               id_column=id_column,
                                                               y_axis=y_axis,
                                                               output_dir=output_dir,
                                                               score_columns=score_columns,
                                                               top_ks=top_ks))
                    plot_files.extend(self.plot_total_identified(nerpa_reports,
                                                                 id_column=id_column,
                                                                 y_axis=y_axis,
                                                                 output_dir=output_dir,
                                                                 score_columns=score_columns))

            return plot_files
    


def nerpa1_report_to_nerpa2_compatible(nerpa1_report: pd.DataFrame,
                                       genome_ids_to_filter: Optional[List[str]] = None) -> NerpaReport:
    """
    Convert a Nerpa1 report to a Nerpa2 compatible format.
    Nerpa1 report should have columns: 'NRP_ID', 'Genome_ID', 'LogOdds score', 'p_value'.
    """
    def pred_file_name_to_genome_id(pred_file_name: str) -> str:
        """
        Extract genome ID from the prediction
        """
        name = pred_file_name.split('/')[-1]  # Assuming the first part is the genome ID
        return name.split('_')[0]  # Extract the genome ID from the file name

    nerpa2_report = pd.DataFrame()
    nerpa2_report['LogOdds score'] = nerpa1_report['Score']
    nerpa2_report['NRP_ID'] = nerpa1_report['StructureId'].apply(lambda s: s.split('_')[0])  # Extract NRP ID from 'Structure ID'
    nerpa2_report['Genome_ID'] = nerpa1_report['PredictionFileName'].apply(pred_file_name_to_genome_id)  # Extract Genome ID from 'BGC ID'

    '''
    # q: group by NRP_ID and Genome_ID, and take the first row with the highest LogOdds score
    nerpa2_report = (nerpa2_report
                     .sort_values(by='LogOdds score', ascending=False)
                     .drop_duplicates(subset=['NRP_ID', 'Genome_ID'], keep='first')
                     .reset_index(drop=True))

    if genome_ids_to_filter is not None:
        # q: filter by genome IDs
        nerpa2_report = nerpa2_report[nerpa2_report['Genome_ID'].isin(genome_ids_to_filter)]

    # q: leave only 10 best matches per Genome_ID
    nerpa2_report = (nerpa2_report
                     .groupby('Genome_ID')
                     .head(10)
                     .reset_index(drop=True))
    '''

    return nerpa2_report


def nerpa1_vs_nerpa2():
    nerpa2_report_tsv = Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine/report.tsv')
    nerpa2_report = pd.read_csv(nerpa2_report_tsv, sep='\t')
    nerpa2_new_report_tsv = Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine2/report.tsv')
    nerpa2_new_report = pd.read_csv(nerpa2_new_report_tsv, sep='\t')

    nerpa1_report_csv = Path('/data/for_training_and_testing/nerpa1_report_mibig4_vs_mibig_norine.csv')
    nerpa1_report = pd.read_csv(nerpa1_report_csv)

    helper = PlotsHelper()
    mibig_nrps_table = helper.data_helper.mibig_nrps_table
    bgc_ids = mibig_nrps_table[mibig_nrps_table['in_approved_matches'] == True]['bgc_id'].unique()

    nerpa1_report = nerpa1_report_to_nerpa2_compatible(nerpa1_report,
                                                       genome_ids_to_filter=bgc_ids.tolist())

    print(nerpa1_report.head())
    print(nerpa2_report.head())

    """
    assert set(nerpa1_report['Genome_ID']) == set(nerpa2_report['Genome_ID']), \
        ("Genome IDs in Nerpa1 and Nerpa2 reports do not match\n"
         f'nerpa1 \ nerpa2 = {set(nerpa1_report["Genome_ID"]) - set(nerpa2_report["Genome_ID"])}\n'
         f'nerpa2 \ nerpa1 = {set(nerpa2_report["Genome_ID"]) - set(nerpa1_report["Genome_ID"])}')
    """
    output_dir = Path('/home/ilianolhin/git/nerpa2/benchmarking/nerpa1_vs_nerpa2_plots')
    helper.plot_all({'Nerpa 1': nerpa1_report,
                     'Nerpa 2': nerpa2_report,
                     'Nerpa 2 new': nerpa2_new_report},
                    output_dir=output_dir)

    def bgc_id_identified(report: NerpaReport, genome_id: str) -> bool:
        """
        Check if a BGC with the given genome_id is identified in the report.
        """
        return (report[report['Genome_ID'] == genome_id]
                .apply(helper.data_helper.match_is_correct, axis=1)
                .any())

    # q: print all Genome_IDs in identified in nerpa1_report but not in nerpa2_report
    identified_in_nerpa1 = {bgc_id for bgc_id in set(nerpa1_report['Genome_ID'])
                            if bgc_id_identified(nerpa1_report, bgc_id)}
    identified_in_nerpa2 = {bgc_id for bgc_id in set(nerpa2_report['Genome_ID'])
                            if bgc_id_identified(nerpa2_report, bgc_id)}
    print('Identified in Nerpa1 but not in Nerpa2:'
          f'\n{identified_in_nerpa1 - identified_in_nerpa2}')



def nerpa1_vs_nerpa2_vs_nerpa2new():
    helper = PlotsHelper()

    # 1. Load reports
    nerpa1_report = helper.data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/data/for_training_and_testing/nerpa1_report_mibig4_vs_mibig_norine.csv'),
                                                         nerpa_version='Nerpa 1',
                                                         report_name='Nerpa 1')
    nerpa2_report = helper.data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine/report.tsv'),
                                                         report_name='Nerpa 2')
    nerpa2_new_vs_avg_bgc = helper.data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine_new/report.tsv'),
                                                                 score_column='LogOdds_vs_avg_BGC',
                                                                 report_name='Nerpa 2 new')
    nerpa2_new_vs_avg_nrp = helper.data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine_new/report.tsv'),
                                                                 score_column='LogOdds_vs_avg_NRP',
                                                                 report_name='Nerpa 2 new')

    # 2. Check for reports discrepancies
    helper.data_helper.check_report_discrepancy(nerpa1_report, nerpa2_report,
                                                'Nerpa 1', 'Nerpa 2',
                                                check_nrp_ids=False)
    helper.data_helper.check_report_discrepancy(nerpa2_report, nerpa2_new_vs_avg_bgc,
                                                'Nerpa 2', 'Nerpa 2 new vs avg BGC',
                                                check_nrp_ids=False)


    # 4. Generate plots
    print('Generating plots...')
    output_dir = Path('/home/ilianolhin/git/nerpa2/benchmarking/nerpa1_vs_nerpa2_plots')
    helper.plot_all({'Nerpa 1': nerpa1_report,
                     'Nerpa 2': nerpa2_report,
                     'Nerpa 2 score vs avg NRP': nerpa2_new_vs_avg_nrp,
                     'Nerpa 2 score vs avg BGC': nerpa2_new_vs_avg_bgc},
                    output_dir=output_dir)


def cross_validation():
    benchmarking_results_dir = Path('/home/ilianolhin/git/nerpa2/benchmarking/cross_validation_results')
    nerpa_reports = {
        f'subset_{i}': pd.read_csv(benchmarking_results_dir / f'subset_{i}' / 'report.tsv', sep='\t')
        for i in range(5)  # Assuming you have 5 samples
    }
    output_dir = Path(benchmarking_results_dir / 'plots')

    helper = PlotsHelper()
    helper.plot_all(nerpa_reports, output_dir)


if __name__ == "__main__":
    # Example usage
    #nerpa1_vs_nerpa2()
    nerpa1_vs_nerpa2_vs_nerpa2new()
    #cross_validation()
    #log_odds_vs_p_values()

'''
def log_odds_vs_p_values():
    nerpa_report_tsv = Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine/report.tsv')
    nerpa_report = pd.read_csv(nerpa_report_tsv, sep='\t')
    output_dir = Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine/plots')
    helper = PlotsHelper()
    helper.plot_log_scores_vs_p_values(nerpa_report, output_dir)
    
    def plot_log_scores_vs_p_values(self,
                                    nerpa_report: NerpaReport,
                                    output_dir: Path,
                                    num_bins: int = 20,
                                    top_ks: Tuple[int, ...] = (1, 3, 5, 10)) -> List[Path]:
        nerpa_reports = {'nerpa_results_log_score': nerpa_report,
                         'nerpa_results_p_value': nerpa_report}
        score_columns = {
            'nerpa_results_log_score': 'LogOdds score',
            'nerpa_results_p_value': 'p_value'
        }

        plots = []
        plots.extend(self.plot_score_correctness(nerpa_reports,
                                                 score_columns,
                                                 output_dir,
                                                 num_bins))
        for y_axis in ('Count', 'Percentage'):
            plots.extend(self.plot_num_identified(nerpa_reports,
                                                  id_column='Genome_ID',
                                                  output_dir=output_dir,
                                                  score_columns=score_columns,
                                                  top_ks=top_ks,
                                                  y_axis=y_axis))
            plots.extend(self.plot_num_identified(nerpa_reports,
                                                  id_column='NRP_ID',
                                                  output_dir=output_dir,
                                                  score_columns=score_columns,
                                                  top_ks=top_ks,
                                                  y_axis=y_axis))
            plots.extend(self.plot_num_correct_matches(nerpa_reports,
                                                       output_dir,
                                                       score_columns=score_columns,
                                                       y_axis=y_axis,
                                                       in_one_plot=True))
        return plots
'''