from __future__ import annotations
from itertools import chain
from pathlib import Path
import os
from typing import Literal, Dict, Optional, Sequence, List, Tuple, NamedTuple, Set

import numpy as np
import pandas as pd
from collections import defaultdict, Counter

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from sklearn.linear_model import LinearRegression


class NerpaReport(pd.DataFrame):
    NRP_ISO_CLASS_COL = 'nrp_iso_class'  # absent in the original report, but added here for convenience
    NRP_ID_COL = 'NRP_ID'
    BGC_ID_COL = 'Genome_ID'
    SCORE_COL = 'Score'

    name: str

    def __init__(self, *args, report_name: str = 'REPORT_NAME_MISSING', **kwargs):
        super().__init__(*args, **kwargs)
        self.name = report_name


class OutputSizeConfig(NamedTuple):
    min_num_matches_per_bgc: int
    max_num_matches_per_bgc: int
    min_num_matches_per_nrp: int
    max_num_matches_per_nrp: int
    max_num_matches: int

NerpaReportRow = pd.Series


def is_bgc_or_norine_nrp(nrp_id: str) -> bool:
    return nrp_id.startswith('BGC') or nrp_id.startswith('NOR')


def nerpa1_report_to_nerpa2_compatible(nerpa1_report: pd.DataFrame) -> pd.DataFrame:
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
    nerpa2_report[NerpaReport.SCORE_COL] = nerpa1_report['Score']
    nerpa2_report[NerpaReport.NRP_ID_COL] = nerpa1_report['StructureId'].apply(lambda s: s.split('_')[0])  # Extract NRP ID from 'Structure ID'
    nerpa2_report[NerpaReport.BGC_ID_COL] = nerpa1_report['PredictionFileName'].apply(pred_file_name_to_genome_id)  # Extract Genome ID from 'BGC ID'

    return nerpa2_report


def truncate_output(report: NerpaReport,
                    out_size_cfg: OutputSizeConfig) -> NerpaReport:
    """
    Truncate the report to satisfy the output size configuration.
    """
    num_matches_per_bgc = Counter()
    num_matches_per_nrp = Counter()

    truncated_report = []
    for _, row in report.iterrows():
        bgc_id = row[NerpaReport.BGC_ID_COL]
        nrp_id = row[NerpaReport.NRP_ID_COL]

        can_add_bgc = num_matches_per_bgc[bgc_id] < out_size_cfg.max_num_matches_per_bgc
        can_add_nrp = num_matches_per_nrp[nrp_id] < out_size_cfg.max_num_matches_per_nrp

        should_add_bgc = num_matches_per_bgc[bgc_id] < out_size_cfg.min_num_matches_per_bgc
        should_add_nrp = num_matches_per_nrp[nrp_id] < out_size_cfg.min_num_matches_per_nrp

        if (can_add_bgc and can_add_nrp) or should_add_bgc or should_add_nrp:
            truncated_report.append(row)
            num_matches_per_bgc[bgc_id] += 1
            num_matches_per_nrp[nrp_id] += 1

            if len(truncated_report) >= out_size_cfg.max_num_matches:
                break

    return NerpaReport(pd.DataFrame(truncated_report).reset_index(drop=True),
                       report_name=report.name)


def prepare_data_for_benchmarking(_report: pd.DataFrame,
                                  bgc_ids_to_keep: Set[str],
                                  nrp_ids_to_keep: Set[str],
                                  nrp_id_to_iso_class: Dict[str, int],
                                  out_size_cfg: OutputSizeConfig,
                                  report_name: str = 'REPORT_NAME_MISSING') -> NerpaReport:
    report = _report.copy()

    # Filter report to only include BGCs present in the approved matches
    report = report[report[NerpaReport.BGC_ID_COL].isin(bgc_ids_to_keep)]
    if set(report[NerpaReport.BGC_ID_COL]) != set(bgc_ids_to_keep):
        print(f"Warning! Some BGC IDs are missing in the report {report_name}")
        print(f"Missing BGC IDs: {set(bgc_ids_to_keep) - set(report['Genome_ID'])}")


    # Filter report to only include NRPs from MIBiG + Norine
    report = (
        report[report[NerpaReport.NRP_ID_COL]
        .map(lambda nrp_id: nrp_id in nrp_ids_to_keep)]
    )

    # Leave only one of isomorphic NRPs per BGC
    report[NerpaReport.NRP_ISO_CLASS_COL] = (
        report[NerpaReport.NRP_ID_COL]
        .map(nrp_id_to_iso_class)
    )
    report.drop_duplicates(subset=[NerpaReport.BGC_ID_COL, NerpaReport.NRP_ISO_CLASS_COL],
                           inplace=True)

    report = report.sort_values(by=NerpaReport.SCORE_COL, ascending=False)
    return truncate_output(NerpaReport(report, report_name=report_name),
                           out_size_cfg)


def load_nerpa_report(report_path: Path,
                      bgc_ids_to_keep: Set[str],
                      nrp_ids_to_keep: Set[str],
                      nrp_id_to_iso_class: Dict[str, int],
                      out_size_cfg: OutputSizeConfig,
                      report_name: str = 'REPORT_NAME_MISSING',
                      score_column: str = NerpaReport.SCORE_COL,
                      tool_version: Literal['Nerpa 1', 'Nerpa 2'] = 'Nerpa 2') -> NerpaReport:
    report = pd.read_csv(report_path, sep='\t' if tool_version == 'Nerpa 2' else ',')
    if tool_version == 'Nerpa 1':
        report = nerpa1_report_to_nerpa2_compatible(report)

    report[NerpaReport.SCORE_COL] = report[score_column]
    return prepare_data_for_benchmarking(_report=report,
                                         report_name=report_name,
                                         bgc_ids_to_keep=bgc_ids_to_keep,
                                         nrp_ids_to_keep=nrp_ids_to_keep,
                                         nrp_id_to_iso_class=nrp_id_to_iso_class,
                                         out_size_cfg=out_size_cfg)


def check_reports_discrepancy(data_helper,  # type : PlotsDataHelper. Can't import due to circular dependency
                              report1: NerpaReport, report2: NerpaReport) -> None:
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
        identified_ids_in_report1 = data_helper.get_identified_ids(report1, id_column)
        identified_ids_in_report2 = data_helper.get_identified_ids(report2, id_column)
        if set(identified_ids_in_report1) != set(identified_ids_in_report2):
            print(f"Warning! Identified {id_column}s in reports "
                  f"{report1.name} and {report2.name} do not match")
            print(f'Identified in {report1.name} but not in {report2.name}:'
                  f'\n{set(identified_ids_in_report1) - set(identified_ids_in_report2)}')
            print(f'Identified in {report2.name} but not in {report1.name}:'
                  f'\n{set(identified_ids_in_report2) - set(identified_ids_in_report1)}')
        else:
            print(f'Identified {id_column}s match')

'''
mibig_norine_iso_classes = set(map(self.nrp_id_to_iso_class, self.mibig_norine_nrps))
if set(report['iso_class_idx']) != mibig_norine_iso_classes:
    print(f"Warning! Some NRP iso-classes from MIBiG + Norine are missing in the report {report_name}")
    print(f"Missing iso classes: {mibig_norine_iso_classes - set(report['iso_class_idx'])}")
'''
