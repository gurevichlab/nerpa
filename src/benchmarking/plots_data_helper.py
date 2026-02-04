from __future__ import annotations
from itertools import chain, combinations
from pathlib import Path
from typing import Literal, Dict, Sequence, List, Optional, Set, Tuple

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
from src.generic.graphs import DSU

PCS = PNRPDB_Compound_Similarity

def sanity_check_similarity_table(similarity_dict: Dict[COMPARISION_METHOD, Set[Tuple[NRP_ID, NRP_ID]]],
                                  pnrpdb_info: PNRPDB_Info) -> None:

    nrp_id_to_iso_class = {
        row[PNRPDB_Info.COMPOUND_ID]: row[PNRPDB_Info.ISO_CLASS_ID]
        for row in pnrpdb_info.iter_rows(named=True)
    }
    ids_in_similarity_table = set(
        chain.from_iterable(
            (nrp1_id, nrp2_id)
            for cmp_mode, similar_nrp_pairs in similarity_dict.items()
            for nrp1_id, nrp2_id in similar_nrp_pairs
        )
    )
    check_successful = True
    out_file = Path('./similarity_table_inconsistencies.txt')
    out_file.parent.mkdir(parents=True, exist_ok=True)

    with open(out_file, 'w') as f:
        for nrp_id, repr_id in nrp_id_to_iso_class.items():
            if nrp_id == repr_id or nrp_id not in ids_in_similarity_table:
                continue
            if (repr_id, nrp_id) not in similarity_dict[PCS.NERPA_EQUAL_ALLOW_UNK_CHR]:
                f.write(f'{(nrp_id, repr_id)} is absent in similarity_dict[PCS.NERPA_EQUAL_ALLOW_UNK_CHR].\n')
                check_successful = False

        for similar_nrp_pairs in similarity_dict[PCS.NERPA_EQUAL]:
            for nrp1_id, nrp2_id in similar_nrp_pairs:
                if nrp_id_to_iso_class[nrp1_id] != nrp_id_to_iso_class[nrp2_id]:
                    f.write(f'{(nrp1_id, nrp2_id)} from similarity_dict are in different iso-classes in pnrpdb_info.\n')
                    check_successful = False

    if check_successful:
        print('Sanity check of similarity table successful.')
    else:
        print(f'Sanity check of similarity table failed. See {out_file} for details.')



BGC_ID = str
NRP_ID = str
COMPARISION_METHOD = str

def get_similarity_dict(pnrpdb_compound_similarity: PNRPDB_Compound_Similarity)\
        -> Dict[COMPARISION_METHOD, Set[Tuple[NRP_ID, NRP_ID]]]:
    '''
    comparison_method -> set of (nrp1_id, nrp2_id) tuples that are similar according to this method
    '''
    similarity_dict = defaultdict(set)
    for row in pnrpdb_compound_similarity.iter_rows(named=True):
        for cmp_method in [
            PCS.NERPA_EQUAL_ALLOW_UNK_CHR,
            PCS.NERPA_ONE_SUB_ALLOW_UNK_CHR,
        ]:
            if row[cmp_method]:
                similarity_dict[cmp_method].add(
                    (row[PCS.FST_COMPOUND_ID], row[PCS.SND_COMPOUND_ID])
                )
                similarity_dict[cmp_method].add(
                    (row[PCS.SND_COMPOUND_ID], row[PCS.FST_COMPOUND_ID])
                )

    return similarity_dict


def get_match_correct_dict(bgc_to_nrps: Dict[BGC_ID, Set[NRP_ID]],
                           similarity_dict: Dict[COMPARISION_METHOD, Set[Tuple[NRP_ID, NRP_ID]]]) \
    -> Dict[COMPARISION_METHOD, Set[Tuple[BGC_ID, NRP_ID]]]:
    '''
    comparison_method -> set of (bgc_id, nrp_id) tuples that are correct matches according to this method
    '''
    match_correct_dict = defaultdict(set)
    nrp_id_to_bgc_id = {
        nrp_id: bgc_id
        for bgc_id, nrp_ids in bgc_to_nrps.items()
        for nrp_id in nrp_ids
    }

    for cmp_method, similar_nrp_pairs in similarity_dict.items():
        for nrp1_id, nrp2_id in similar_nrp_pairs:
            if nrp1_id in nrp_id_to_bgc_id:  # similarity_dict is symmetric so no need to check the reverse
                bgc_id = nrp_id_to_bgc_id[nrp1_id]
                match_correct_dict[cmp_method].add((bgc_id, nrp2_id))

    # Because reflexivity is not included in similarity_dict, add (bgc_id, nrp_id) for each known correct match
    for nrp_id, bgc_id in nrp_id_to_bgc_id.items():
        for cmp_method in (
            PCS.NERPA_EQUAL,
            PCS.NERPA_EQUAL_ALLOW_UNK_CHR,
            PCS.NERPA_NO_MORE_ONE_SUB,
            PCS.NERPA_NO_MORE_ONE_SUB_ALLOW_UNK_CHR,
        ):
            match_correct_dict[cmp_method].add((bgc_id, nrp_id))

    return match_correct_dict


def get_nrp_id_to_iso_class(similarity_dict: Dict[COMPARISION_METHOD, Set[Tuple[NRP_ID, NRP_ID]]],
                            pnrpdb_info: PNRPDB_Info,
                            cmp: str = PCS.NERPA_EQUAL) \
    -> Dict[NRP_ID, NRP_ID]:
    dsu = DSU()
    for nrp1_id, nrp2_id in similarity_dict[cmp]:
        dsu.union(nrp1_id, nrp2_id)

    nrp_id_to_iso_class = {nrp_id: repr_id for nrp_id, repr_id in dsu.items()}
    for nrp_id in pnrpdb_info[PNRPDB_Info.COMPOUND_ID]:
        if nrp_id not in nrp_id_to_iso_class:
            nrp_id_to_iso_class[nrp_id] = nrp_id  # singleton iso-class

    return nrp_id_to_iso_class


def add_similarities_from_pnrpdb_info(similarity_dict: Dict[COMPARISION_METHOD, Set[Tuple[NRP_ID, NRP_ID]]],
                                      pnrpdb_info: PNRPDB_Info) -> None:
    iso_class_to_nrp_ids = defaultdict(set)
    for row in pnrpdb_info.iter_rows(named=True):
        iso_class_to_nrp_ids[row[PNRPDB_Info.ISO_CLASS_ID]].add(row[PNRPDB_Info.COMPOUND_ID])

    for nrp_ids in iso_class_to_nrp_ids.values():
        for nrp1_id, nrp2_id in combinations(nrp_ids, 2):
                similarity_dict[PCS.NERPA_EQUAL_ALLOW_UNK_CHR].add((nrp1_id, nrp2_id))
                similarity_dict[PCS.NERPA_EQUAL_ALLOW_UNK_CHR].add((nrp2_id, nrp1_id))
                similarity_dict[PCS.NERPA_EQUAL].add((nrp1_id, nrp2_id))
                similarity_dict[PCS.NERPA_EQUAL].add((nrp2_id, nrp1_id))


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
        self.similarity_dict = get_similarity_dict(self.pnrpdb_compound_similarity)
        sanity_check_similarity_table(self.similarity_dict,
                                      self.pnrpdb_info)
        # right now similarity table is not complete, so I update it from pnrpdb_info
        # to ensure that all compounds in the same iso-class are similar
        add_similarities_from_pnrpdb_info(self.similarity_dict, self.pnrpdb_info)
        self.similarity_dict[PCS.NERPA_NO_MORE_ONE_SUB_ALLOW_UNK_CHR] = (
            self.similarity_dict[PCS.NERPA_EQUAL_ALLOW_UNK_CHR] |
            self.similarity_dict[PCS.NERPA_ONE_SUB_ALLOW_UNK_CHR]
        )
        self.nrp_id_to_iso_class = get_nrp_id_to_iso_class(self.similarity_dict,
                                                           self.pnrpdb_info)

        bgc_to_nrps = defaultdict(set)
        for row in self.mibig_bgcs_info.iter_rows(named=True):
            bgc_to_nrps[row[MIBiG_BGCs_Info.BGC_ID]].add(row[MIBiG_BGCs_Info.NRP_ID])


        self.bgc_to_nrp_iso_classes = defaultdict(set)
        for row in self.mibig_bgcs_info.iter_rows(named=True):
            bgc_id = row[MIBiG_BGCs_Info.BGC_ID]
            nrp_id = row[MIBiG_BGCs_Info.NRP_ID]
            nrp_class = self.nrp_id_to_iso_class[nrp_id]
            self.bgc_to_nrp_iso_classes[bgc_id].add(nrp_class)

        self.match_correct_dict = get_match_correct_dict(bgc_to_nrps, self.similarity_dict)
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
                         cmp_mode: str = PCS.NERPA_EQUAL_ALLOW_UNK_CHR) -> bool:
        """Check if a match between NRP iso-class and BGC is correct."""
        return (bgc_id, nrp_iso_class) in self.match_correct_dict[cmp_mode]

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
                                 y_axis: Literal['Count', 'Percentage'] = 'Count',
                                 cmp_method: str = PCS.NERPA_EQUAL_ALLOW_UNK_CHR) -> pl.Series:
        return compute_total_identified(self, nerpa_report, id_column, max_top_k, y_axis, cmp_method)

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
