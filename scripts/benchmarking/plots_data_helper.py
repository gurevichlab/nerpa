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

from scripts.benchmarking.nerpa_report import (
    NerpaReport,
    NerpaReportRow,
    load_nerpa_report,
    is_bgc_or_norine_nrp
)


class PlotsDataHelper:
    mibig_bgcs_info: pd.DataFrame
    pnrpdb_info: pd.DataFrame

    bgc_to_nrp_iso_classes: defaultdict[str, set[int]]
    nrp_id_to_iso_class: Dict[str, int]

    def __init__(self):
        nerpa_dir = Path(os.getcwd()).parent.parent
        self.mibig_bgcs_info = pd.read_csv(nerpa_dir / 'data/for_training_and_testing/mibig_bgcs_info.tsv',
                                           sep='\t')
        self.pnrpdb_info = pd.read_csv(nerpa_dir / 'data/for_training_and_testing/pnrpdb2_additional_info.tsv',
                                       sep='\t')
        self.mibig_norine_nrps = set(
            self.pnrpdb_info[self.pnrpdb_info['compound_id'].map(is_bgc_or_norine_nrp)]['compound_id'])

        self.bgc_to_nrp_iso_classes = defaultdict(set)
        for _, row in self.mibig_bgcs_info.iterrows():
            self.bgc_to_nrp_iso_classes[row['bgc_id']].add(row['iso_class_idx'])

        # q: compute nrp_id_to_iso_class only once
        self.nrp_id_to_iso_class = {row['compound_id']: row['iso_class_idx']
                                    for _, row in self.pnrpdb_info.iterrows()}

    def match_is_correct(self, report_row: NerpaReportRow):  # row in report.tsv
        try:
            nrp_iso_class = report_row[NerpaReport.NRP_ISO_CLASS_COL]
            bgc_id = report_row[NerpaReport.BGC_ID_COL]
            return nrp_iso_class in self.bgc_to_nrp_iso_classes[bgc_id]
        except KeyError:
            print('Invalid row: ', report_row)
            raise KeyError("Row does not contain valid NRP_ID or Genome_ID")

    def compute_num_correct_matches(self,
                                    nerpa_report: NerpaReport) -> pd.Series[int]:
        """
        Returns an array where the i-th element is the number of correct matches
        among the first i matches in the nerpa_report, sorted by the specified column.
        """
        num_correct = (nerpa_report
                       .sort_values(by=NerpaReport.SCORE_COL, ascending=False)
                       .apply(self.match_is_correct, axis=1)
                       .cumsum())
        # Prepend the “0 correct at N=0” point
        return pd.concat([pd.Series([0]), num_correct], ignore_index=True)

    def compute_score_correctness(self,
                                  nerpa_report: NerpaReport,
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
        df["score_bin"] = pd.cut(df[NerpaReport.SCORE_COL], bins=num_bins, include_lowest=True)
        return (
            df.groupby(["score_bin", "is_correct"], observed=False)
            .size()
            .unstack(fill_value=0)  # columns: False / True
            .rename(columns={False: "Incorrect", True: "Correct"})
        )

    def compute_num_identified(self,
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
                                 id_column: str,
                                 max_top_k: int = 10,
                                 y_axis: Literal['Count', 'Percentage'] = 'Count') -> pd.Series[int]:
        assert id_column in (NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL), \
            f"id_column must be either '{NerpaReport.BGC_ID_COL}' or '{NerpaReport.NRP_ISO_CLASS_COL}'"

        total_identified = []
        for top_k in range(1, max_top_k + 1):
            num_identified = self.compute_num_identified(nerpa_report, id_column, top_k)
            value = 100 * num_identified.iloc[-1] / len(num_identified) \
                if y_axis == 'Percentage' else num_identified.iloc[-1]
            total_identified.append(value)

        return pd.Series(total_identified, index=range(1, max_top_k + 1))

    def get_identified_ids(self,
                           report: NerpaReport,
                           id_column: str) -> Sequence[str]:
        assert id_column in (NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL), \
            f"id_column must be either '{NerpaReport.BGC_ID_COL}' or '{NerpaReport.NRP_ISO_CLASS_COL}'"

        def id_identified(report: NerpaReport, id_value: str) -> bool:
            return (report[report[id_column] == id_value]
                    .apply(self.match_is_correct, axis=1)
                    .any())

        return [id_value for id_value in set(report[id_column])
                if id_identified(report, id_value)]

    def check_reports_discrepancy(self, report1: NerpaReport, report2: NerpaReport) -> None:
        """
        Check for discrepancies between two reports.
        """
        print('Checking report discrepancies between '
              f'{report1.name} and {report2.name}...')

        if set(report1[NerpaReport.BGC_ID_COL]) != set(report2[NerpaReport.BGC_ID_COL]):
            print(f"Warning! Genome IDs in reports "
                  f"{report1.name} and {report2.name} do not match")
            print(f'{report1.name} \ {report2.name} = '
                  f'{set(report1[NerpaReport.BGC_ID_COL]) - set(report2[NerpaReport.BGC_ID_COL])}')
            print(f'{report2.name} \ {report1.name} = '
                  f'{set(report2[NerpaReport.BGC_ID_COL]) - set(report1[NerpaReport.BGC_ID_COL])}')
        else:
            print("Genome IDs sets match")


        for id_column in (NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL):
            identified_ids_in_report1 = self.get_identified_ids(report1, id_column)
            identified_ids_in_report2 = self.get_identified_ids(report2, id_column)
            if set(identified_ids_in_report1) != set(identified_ids_in_report2):
                print(f"Warning! Identified {id_column}s in reports "
                      f"{report1.name} and {report2.name} do not match")
                print(f'Identified in {report1.name} but not in {report2.name}:'
                      f'\n{set(identified_ids_in_report1) - set(identified_ids_in_report2)}')
                print(f'Identified in {report2.name} but not in {report1.name}:'
                      f'\n{set(identified_ids_in_report2) - set(identified_ids_in_report1)}')
            else:
                print(f'Identified {id_column}s match')

    def load_nerpa_report(self,
                          report_path: Path,
                          report_name: str = 'REPORT_NAME_MISSING',
                          score_column: str = 'LogOdds score') -> NerpaReport:
        appr_matches_rows = self.mibig_bgcs_info[self.mibig_bgcs_info['is_approved_match']]
        bgc_ids_to_keep = appr_matches_rows['bgc_id'].dropduplicates()
        nrp_ids_to_keep = appr_matches_rows['compound_id'].dropduplicates()
        out_cfg = NerpaReport.OutputSizeConfig(
            max_num_matches=float('inf'),
            max_num_matches_per_bgc=10,
            max_num_matches_per_nrp=10,
            min_num_matches_per_bgc=10,
            min_num_matches_per_nrp=10,
        )

        return load_nerpa_report(report_path=report_path,
                                 bgc_ids_to_keep=bgc_ids_to_keep,
                                 nrp_ids_to_keep=nrp_ids_to_keep,
                                 nrp_id_to_iso_class=self.nrp_id_to_iso_class,
                                 out_size_cfg=out_cfg,
                                 report_name=report_name,
                                 score_column=score_column)
