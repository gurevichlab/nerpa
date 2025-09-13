from __future__ import annotations
from itertools import chain
from pathlib import Path
import os
from typing import Literal, Dict, Optional, Sequence, List, Tuple

import numpy as np
import pandas as pd
from collections import defaultdict
import polars as pl

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from sklearn.linear_model import LinearRegression

from scripts.benchmarking.data_helper_external_methods import compute_num_correct_matches, compute_score_correctness, \
    compute_num_identified, compute_total_identified, get_identified_ids, promiscuity_handling_stats
from scripts.benchmarking.nerpa_report import (
    NerpaReport,
    NerpaReportRow,
    load_nerpa_report,
    is_bgc_or_norine_nrp, OutputSizeConfig, check_reports_discrepancy
)


class PlotsDataHelper:
    mibig_bgcs_info: pd.DataFrame
    pnrpdb_info: pd.DataFrame

    mibig_norine_nrps: set[str]
    bgc_to_nrp_iso_classes: defaultdict[str, set[int]]
    nrp_id_to_iso_class: Dict[str, int]

    def __init__(self):
        nerpa_dir = Path(__file__).parent.parent.parent
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
        return compute_num_correct_matches(self, nerpa_report)

    def compute_score_correctness(self,
                                  nerpa_report: NerpaReport,
                                  num_bins: int = 20,
                                  bin_edges: List[float] = None) -> pd.DataFrame:
        return compute_score_correctness(self, nerpa_report, num_bins, bin_edges)

    def compute_num_identified(self,
                               nerpa_report: NerpaReport,
                               id_column: str,
                               top_k: int = 1) -> pd.Series[int]:
        return compute_num_identified(self, nerpa_report, id_column, top_k)

    def compute_total_identified(self,
                                 nerpa_report: NerpaReport,
                                 id_column: str,
                                 max_top_k: int = 10,
                                 y_axis: Literal['Count', 'Percentage'] = 'Count') -> pd.Series[int]:
        return compute_total_identified(self, nerpa_report, id_column, max_top_k, y_axis)

    def get_identified_ids(self,
                           report: NerpaReport,
                           id_column: str) -> Sequence[str]:
        return get_identified_ids(self, report, id_column)

    def promiscuity_handling_stats(self,
                                   nerpa_report: NerpaReport) -> pd.DataFrame:
        return promiscuity_handling_stats(self, nerpa_report)

    def check_reports_discrepancy(self,  # type : PlotsDataHelper. Can't import due to circular dependency
                                  report1: NerpaReport, report2: NerpaReport) -> None:
        check_reports_discrepancy(self, report1, report2)

    def load_nerpa_report(self,
                          report_path: Path,
                          report_name: str = 'REPORT_NAME_MISSING',
                          score_column: str = 'LogOdds score',
                          tool_version: Literal['Nerpa 1, Nerpa 2'] = 'Nerpa 2') -> NerpaReport:
        appr_matches_rows = self.mibig_bgcs_info[self.mibig_bgcs_info['in_approved_matches']]
        bgc_ids_to_keep = set(appr_matches_rows['bgc_id'].drop_duplicates())
        nrp_ids_to_keep = self.mibig_norine_nrps
        out_cfg = OutputSizeConfig(
            max_num_matches=float('inf'),  # no limit on total matches
            max_num_matches_per_bgc=10,
            max_num_matches_per_nrp=10,
            min_num_matches_per_bgc=10,
            min_num_matches_per_nrp=10,
        )

        report = load_nerpa_report(report_path=report_path,
                                   bgc_ids_to_keep=bgc_ids_to_keep,
                                   nrp_ids_to_keep=nrp_ids_to_keep,
                                   nrp_id_to_iso_class=self.nrp_id_to_iso_class,
                                   out_size_cfg=out_cfg,
                                   report_name=report_name,
                                   score_column=score_column,
                                   tool_version=tool_version)

        report['is_correct'] = report.apply(self.match_is_correct, axis=1)
        return report
