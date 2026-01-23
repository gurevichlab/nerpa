from __future__ import annotations
from itertools import chain
from pathlib import Path
from typing import Literal, Dict, Sequence, List, Optional

import polars as pl
from collections import defaultdict

from src.benchmarking.data_frames import MIBiG_BGCs_Info, PNRPDB_Info, PNRPDB_Compound_Similarity
from src.benchmarking.data_helper_external_methods import (
    compute_num_correct_matches, compute_score_correctness,
    compute_num_identified, compute_total_identified, get_identified_ids,
    promiscuity_handling_stats, group_by_bgc_length, compute_precision_recall_curve, PrecisionRecallPoint,
    extra_false_positives
)
from src.benchmarking.nerpa_report import (
    NerpaReport, load_nerpa_report, is_bgc_or_norine_nrp,
    OutputSizeConfig, check_reports_discrepancy
)

PCS = PNRPDB_Compound_Similarity

def sanity_check_similarity_table(pnrpdb_compound_similarity: PNRPDB_Compound_Similarity,
                                  nrp_id_to_iso_class: Dict[str, str]) -> None:
    for nrp_id, repr_id in nrp_id_to_iso_class.items():
        similarity_info = (
            pnrpdb_compound_similarity
            .filter(
                ((pl.col(PCS.FST_COMPOUND_ID) == nrp_id) & (pl.col(PCS.SND_COMPOUND_ID) == repr_id)) |
                ((pl.col(PCS.FST_COMPOUND_ID) == repr_id) & (pl.col(PCS.SND_COMPOUND_ID) == nrp_id))
            )
            .select(pl.col(PCS.rBAN_ISO_ALLOW_UNK_CHR))
        )
        assert not similarity_info.is_empty() and similarity_info[PCS.rBAN_ISO_ALLOW_UNK_CHR][0], \
            f'Similarity table inconsistency for NRP {nrp_id} and its iso-class representative {repr_id}.'


class PlotsDataHelper:
    mibig_bgcs_info: MIBiG_BGCs_Info
    pnrpdb_info: PNRPDB_Info
    pnrpdb_compound_similarity: PNRPDB_Compound_Similarity

    test_nrp_classes: set[str]
    bgc_to_nrp_iso_classes: defaultdict[str, set[str]]
    nrp_id_to_iso_class: Dict[str, str]

    test_bgcs: set[str]
    nrp_classes_with_matches: set[str]

    def __init__(self,
                 bgc_test_set: Literal['mibig4_wo_training_bgcs', 'training_bgcs'] = 'training_bgcs'):
        nerpa_dir = Path(__file__).parent.parent.parent
        self.mibig_bgcs_info = MIBiG_BGCs_Info.from_csv(
            nerpa_dir / 'data/for_training_and_testing/mibig_bgcs_info.tsv'
        )
        self.pnrpdb_info = PNRPDB_Info.from_csv(
            nerpa_dir / 'data/for_training_and_testing/pnrpdb2_info.tsv'
        )
        self.pnrpdb_compound_similarity = PNRPDB_Compound_Similarity.from_csv(
            nerpa_dir / 'data/for_training_and_testing/pnrpdb2_compound_similarity.tsv'
        )

        bgc_to_nrps = defaultdict(set)
        for row in self.mibig_bgcs_info.iter_rows(named=True):
            bgc_to_nrps[row[MIBiG_BGCs_Info.BGC_ID]].add(row[MIBiG_BGCs_Info.NRP_ID])

        self.nrp_id_to_iso_class = {
            row[PNRPDB_Info.COMPOUND_ID]: row[PNRPDB_Info.ISO_CLASS_ID]
            for row in self.pnrpdb_info.iter_rows(named=True)
        }
        sanity_check_similarity_table(self.pnrpdb_compound_similarity,
                                      self.nrp_id_to_iso_class)

        self.bgc_to_nrp_iso_classes = defaultdict(set)
        for row in self.mibig_bgcs_info.iter_rows(named=True):
            bgc_id = row[MIBiG_BGCs_Info.BGC_ID]
            nrp_id = row[MIBiG_BGCs_Info.NRP_ID]
            nrp_class = self.nrp_id_to_iso_class[nrp_id]
            self.bgc_to_nrp_iso_classes[bgc_id].add(nrp_class)

        training_bgcs_info = self.mibig_bgcs_info.filter(pl.col(MIBiG_BGCs_Info.IN_APPROVED_MATCHES))
        training_bgcs_fams = set(training_bgcs_info[MIBiG_BGCs_Info.BIGSCAPE_FAMILIES])

        dissimilar_to_training_bgcs_info = (
            self.mibig_bgcs_info
            .filter(~pl.col(MIBiG_BGCs_Info.BIGSCAPE_FAMILIES).is_in(training_bgcs_fams))
        )

        match bgc_test_set:
            case 'training_bgcs':
                test_bgcs_info = training_bgcs_info
            case 'mibig4_wo_training_bgcs':
                test_bgcs_info = dissimilar_to_training_bgcs_info
            case _:
                raise ValueError(f'Unknown bgc_test_set: {bgc_test_set}')

        # MIBiG or NORINE NRPs with at least 3 NRPS monomers
        test_nrps = set(
            self.pnrpdb_info
            .filter(is_bgc_or_norine_nrp() &
                    (pl.col(PNRPDB_Info.NUM_NRPS_MONOMERS) >= 3))
            [PNRPDB_Info.COMPOUND_ID]
        )

        self.test_nrp_classes = set(
            self.nrp_id_to_iso_class[nrp_id]
            for nrp_id in test_nrps
        )

        _test_bgcs = set(
            test_bgcs_info
            .filter(pl.col(MIBiG_BGCs_Info.NUM_A_DOMAINS) >= 2)
            [MIBiG_BGCs_Info.BGC_ID]
        )

        # Keep only BGCs that have at least one NRP iso-class in the test NRP classes
        self.test_bgcs = {bgc for bgc in _test_bgcs
                          if any(nrp_class in self.test_nrp_classes
                                 for nrp_class in self.bgc_to_nrp_iso_classes[bgc])}


        self.nrp_classes_with_matches = set(
            chain.from_iterable(
                self.bgc_to_nrp_iso_classes[bgc_id]
                for bgc_id in self.test_bgcs
            )
        )

    def match_is_correct(self, nrp_iso_class: str, bgc_id: str,
                         cmp_mode: str = PCS.rBAN_ISO_ALLOW_UNK_CHR) -> bool:
        """Check if a match between NRP iso-class and BGC is correct."""
        for bgc_iso_class in self.bgc_to_nrp_iso_classes[bgc_id]:
            if bgc_iso_class == nrp_iso_class:
                return True
            similarity_info = (
                self.pnrpdb_compound_similarity
                .filter(
                    ((pl.col(PCS.FST_COMPOUND_ID) == nrp_iso_class) & (pl.col(PCS.SND_COMPOUND_ID) == bgc_iso_class)) |
                    ((pl.col(PCS.FST_COMPOUND_ID) == bgc_iso_class) & (pl.col(PCS.SND_COMPOUND_ID) == nrp_iso_class))
                )
                .select(pl.col(cmp_mode))
            )
            if not similarity_info.is_empty() and similarity_info[cmp_mode][0]:
                return True
        return False

    def compute_num_correct_matches(self, nerpa_report: NerpaReport) -> pl.Series:
        return compute_num_correct_matches(self, nerpa_report)

    def compute_score_correctness(self, nerpa_report: NerpaReport,
                                  num_bins: int = 20,
                                  bin_edges: list[float] = None) -> pl.DataFrame:
        return compute_score_correctness(self, nerpa_report, num_bins, bin_edges)

    def compute_num_identified(self, nerpa_report: NerpaReport,
                               id_column: str,
                               cmp_mode: str,
                               top_k: int = 1) -> pl.Series:
        return compute_num_identified(self, nerpa_report, id_column, cmp_mode, top_k)

    def compute_total_identified(self, nerpa_report: NerpaReport,
                                 id_column: str,
                                 max_top_k: int = 10,
                                 y_axis: Literal['Count', 'Percentage'] = 'Count') -> pl.Series:
        return compute_total_identified(self, nerpa_report, id_column, max_top_k, y_axis)

    def get_identified_ids(self, report: NerpaReport, id_column: str) -> Sequence[str]:
        return get_identified_ids(self, report, id_column)

    def promiscuity_handling_stats(self, nerpa_report: NerpaReport) -> pl.DataFrame:
        return promiscuity_handling_stats(self, nerpa_report)

    def group_by_bgc_length(self, nerpa_report: NerpaReport,
                           num_bins: int = 5) -> dict[tuple, pl.DataFrame]:
        return group_by_bgc_length(self, nerpa_report, num_bins)

    def check_reports_discrepancy(self, report1: NerpaReport, report2: NerpaReport) -> None:
        check_reports_discrepancy(self, report1, report2)

    def load_nerpa_report(self, report_path: Path,
                          report_name: str = 'REPORT_NAME_MISSING',
                          score_column: str = NerpaReport.SCORE,
                          tool_version: Literal['Nerpa 1', 'Nerpa 2'] = 'Nerpa 2') -> NerpaReport:
        nrp_classes_to_keep = self.test_nrp_classes
        bgc_ids_to_keep = self.test_bgcs

        #
        '''
        out_cfg = OutputSizeConfig(
            max_num_matches=float('inf'),  # no limit on total matches
            max_num_matches_per_bgc=10,
            max_num_matches_per_nrp=10,
            min_num_matches_per_bgc=10,
            min_num_matches_per_nrp=10,
        )
        '''

        # No limits
        out_cfg = OutputSizeConfig(
            max_num_matches=float('inf'),  # no limit on total matches
            max_num_matches_per_bgc=float('inf'),#10,
            max_num_matches_per_nrp=float('inf'),#10,
            min_num_matches_per_bgc=0,#10,
            min_num_matches_per_nrp=0,#10,
        )

        report = load_nerpa_report(report_path=report_path,
                                   bgc_ids_to_keep=bgc_ids_to_keep,
                                   nrp_classes_to_keep=nrp_classes_to_keep,
                                   nrp_id_to_iso_class=self.nrp_id_to_iso_class,
                                   out_size_cfg=out_cfg,
                                   report_name=report_name,
                                   score_column=score_column,
                                   tool_version=tool_version)

        # Add is_correct column
        report = report.with_columns([
            pl.struct([NerpaReport.NRP_ISO_CLASS, NerpaReport.BGC_ID])
            .map_elements(lambda row: self.match_is_correct(row[NerpaReport.NRP_ISO_CLASS], row[NerpaReport.BGC_ID]),
                         return_dtype=pl.Boolean)
            .alias(NerpaReport.IS_CORRECT),

            pl.col(NerpaReport.SCORE)
            .rank(method="ordinal", descending=True)
            .over(NerpaReport.BGC_ID)
            .alias(NerpaReport.MATCH_RANK_FOR_BGC),

            pl.col(NerpaReport.SCORE)
            .rank(method="ordinal", descending=True)
            .over(NerpaReport.NRP_ISO_CLASS)
            .alias(NerpaReport.MATCH_RANK_FOR_NRP)
        ]
        )

        # add row index
        report = report.with_row_count(name=NerpaReport.MATCH_RANK)

        return NerpaReport(report, report_name=report_name)

    def compute_precision_recall_curve(self,
                                       nerpa_report: NerpaReport,
                                       top_matches_per_bgc: Optional[int],) -> List[PrecisionRecallPoint]:
        return compute_precision_recall_curve(self, nerpa_report, top_matches_per_bgc)

    def extra_false_positives(self,
                              report1: NerpaReport,
                              report2: NerpaReport) -> pl.DataFrame:
        """Count extra false positives in report1 compared to report2."""
        return extra_false_positives(self, report1, report2)
