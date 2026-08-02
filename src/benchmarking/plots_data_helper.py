from __future__ import annotations
from itertools import chain, combinations
from pathlib import Path
from typing import Literal, Dict, Sequence, List, Optional, Set, Tuple

import polars as pl
from collections import defaultdict
from src.testing.testing_types import CuratedAlignment

from src.benchmarking.data_frames import (
    MIBiG_BGCs_Info,
    PNRPDB_Info,
    PNRPDB_Compound_Similarity,
    NRP_CMP_METHOD,
    NRP_CMP_METHODS,
)
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
from src.rban_parsing.nrp_variant_types import (
    NRP_ID,
    is_mibig_nrp,
    is_norine_nrp,
)

from src.generic.graphs import DSU
import yaml

PCS = PNRPDB_Compound_Similarity
BGC_ID = str


def sanity_check_similarity_table(
        similarity_dict: Dict[NRP_CMP_METHOD, Set[Tuple[NRP_ID, NRP_ID]]],
        pnrpdb_info: PNRPDB_Info,
        out_file: Optional[Path] = None
) -> List[str]:
    """
    Validate consistency between:
      1) `pnrpdb_info` iso-class assignments (COMPOUND_ID -> REPR_ID),
      2) `similarity_dict` contents for specific comparison modes.

    Returns:
      - A list of error messages,

    Checks performed:
      A) For every compound such that nrp_id != repr_id,
         ensure the pair (repr_id, nrp_id) exists in the NERPA_EQUAL_ALLOW_UNK_CHR similarity set
      B) For every pair in NERPA_EQUAL, ensure both compounds are assigned to the same iso-class
         according to `pnrpdb_info`.
    """
    errors = []
    nrp_id_to_iso_class = {
        row[PNRPDB_Info.COMPOUND_ID]: row[PNRPDB_Info.REPR_ID]
        for row in pnrpdb_info.iter_rows(named=True)
    }

    for nrp_id, repr_id in nrp_id_to_iso_class.items():
        if any([
                # similarity dict is only computed for MIBiG and NORINE compounds
                not is_mibig_nrp(nrp_id) and not is_norine_nrp(nrp_id),
                # Nerpa was unable to build an NRP variant for this compound
                repr_id is None,
                # reflexive pairs are absent from similarity_dict
                nrp_id == repr_id,
        ]):
            continue

        if all([
                (repr_id, nrp_id) not in similarity_dict[PCS.NERPA_ISO_WEAK_CHR],
                (nrp_id, repr_id) not in similarity_dict[PCS.NERPA_ISO_WEAK_CHR]
        ]):
            errors.append(
                f'{(nrp_id, repr_id)} is absent in '
                'similarity_dict[PCS.NERPA_ISO_WEAK_CHR].'
            )

    for nrp1_id, nrp2_id in similarity_dict[PCS.NERPA_ISO]:
        if nrp_id_to_iso_class[nrp1_id] != nrp_id_to_iso_class[nrp2_id]:
            errors.append(
                f'{(nrp1_id, nrp2_id)} from similarity_dict (cmp mode {PCS.NERPA_ISO}) '
                'are in different iso-classes in pnrpdb_info:\n'
                f'  {nrp1_id} -> {nrp_id_to_iso_class[nrp1_id]}\n'
                f'  {nrp2_id} -> {nrp_id_to_iso_class[nrp2_id]}'
            )

    if not errors:
        print('No errors found in similarity table.')
    else:
        print(
            f'Errors were found in similarity table '
            f'({len(errors)} in total)'
        )
        if out_file is not None:
            out_file.parent.mkdir(parents=True, exist_ok=True)
            with open(out_file, 'w') as f:
                f.write('\n'.join(errors))
                print(f'Errors log written to {out_file.resolve()}')

    return errors



def get_similarity_dict(
        pnrpdb_compound_similarity: PNRPDB_Compound_Similarity,
        pnrpdb_info: Optional[PNRPDB_Info] = None,  # used only for sanity check
        log_file: Optional[Path] = None
) -> Dict[NRP_CMP_METHOD, Set[Tuple[NRP_ID, NRP_ID]]]:
    '''
    pnrpdb_compound_similarity has structure:
    nrp1_id, nrp2_id, cmp_method1_result, cmp_method2_result, ...

    we convert it to a dictionary of sets for easier access:
    comparison_method -> set of (nrp1_id, nrp2_id) tuples that are similar according to this method
    '''
    similarity_dict = defaultdict(set)
    for row in pnrpdb_compound_similarity.iter_rows(named=True):
        for cmp_method in NRP_CMP_METHODS:
            if row[cmp_method]:
                similarity_dict[cmp_method].add(
                    (row[PCS.FST_COMPOUND_ID], row[PCS.SND_COMPOUND_ID])
                )
                similarity_dict[cmp_method].add(
                    (row[PCS.SND_COMPOUND_ID], row[PCS.FST_COMPOUND_ID])
                )

    if pnrpdb_info is not None and log_file is not None:
        sanity_check_similarity_table(
            similarity_dict,
            pnrpdb_info,
            out_file=log_file
        )

    return similarity_dict


def get_match_correct_dict(
        bgc_to_nrps: Dict[BGC_ID, Set[NRP_ID]],
        similarity_dict: Dict[NRP_CMP_METHOD, Set[Tuple[NRP_ID, NRP_ID]]]
) -> Dict[NRP_CMP_METHOD, Set[Tuple[BGC_ID, NRP_ID]]]:
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
            # if nrp1 is a known correct match to a bgc, and nrp2 is close to nrp1 then (bgc_id, nrp2) is also a correct match
            # no need to check the other direction because similar_nrp_pairs is symmetric
            if nrp1_id in nrp_id_to_bgc_id:
                bgc_id = nrp_id_to_bgc_id[nrp1_id]
                match_correct_dict[cmp_method].add((bgc_id, nrp2_id))

    # Because reflexivity is not included in similarity_dict, add (bgc_id, nrp_id) for each known correct match
    for nrp_id, bgc_id in nrp_id_to_bgc_id.items():
        for cmp_method in NRP_CMP_METHODS:
            match_correct_dict[cmp_method].add((bgc_id, nrp_id))

    return match_correct_dict


class PlotsDataHelper:
    mibig_bgcs_info: MIBiG_BGCs_Info
    pnrpdb_info: PNRPDB_Info
    pnrpdb_compound_similarity: PNRPDB_Compound_Similarity

    test_nrp_classes: set[str]
    bgc_to_nrp_iso_classes: defaultdict[str, set[str]]
    nrp_id_to_iso_class: Dict[str, str]

    test_bgcs: set[str]
    nrp_classes_with_matches: set[str]

    def __init__(
            self,
            mibig_bgcs_info_path: Path,
            pnrpdb_info_path: Path,
            pnrpdb_compound_similarity_path: Path,
            curated_alignments_path: Path,
            out_dir: Path,
            bgc_test_set: Literal[
                "mibig4_wo_training_bgcs", "training_bgcs"
            ] = "training_bgcs",
    ):
        self.mibig_bgcs_info = MIBiG_BGCs_Info.from_csv(mibig_bgcs_info_path)
        self.pnrpdb_info = PNRPDB_Info.from_csv(pnrpdb_info_path)
        self.pnrpdb_compound_similarity = PNRPDB_Compound_Similarity.from_csv(pnrpdb_compound_similarity_path)
        curated_alignments = [
            CuratedAlignment.from_yaml_dict(test_match_dict)
            for test_match_dict in yaml.safe_load(curated_alignments_path.open('r'))
        ]
        self.training_bgcs = set(
            test_match.bgc_id
            for test_match in curated_alignments
        )

        self.similarity_dict = get_similarity_dict(
            self.pnrpdb_compound_similarity,
            self.pnrpdb_info,  # needed only for sanity check
            log_file=out_dir.parent / 'similarity_table_errors.txt'
        )

        self.nrp_id_to_iso_class = {
            row[PNRPDB_Info.COMPOUND_ID]: row[PNRPDB_Info.REPR_ID]
            for row in self.pnrpdb_info.iter_rows(named=True)
        }

        # with open(nerpa_dir / 'tmp/nrp_id_to_iso_class.txt', 'w') as f:
        #     for nrp_id, iso_class in sorted(self.nrp_id_to_iso_class.items()):
        #         f.write(f'{nrp_id}: {iso_class}' + '\n')

        bgc_to_nrps = defaultdict(set)
        for row in self.mibig_bgcs_info.iter_rows(named=True):
            bgc_id = row[MIBiG_BGCs_Info.BGC_ID]
            nrp_ids = (
                row[MIBiG_BGCs_Info.NRP_IDS].split(',')
                if row[MIBiG_BGCs_Info.NRP_IDS]
                else []
            )
            bgc_to_nrps[bgc_id] = set(nrp_ids)


        self.bgc_to_nrp_iso_classes = defaultdict(set)
        for bgc_id, nrp_ids in bgc_to_nrps.items():
            for nrp_id in nrp_ids:
                iso_class = self.nrp_id_to_iso_class.get(nrp_id)
                # missing iso_class means that Nerpa was unable
                # to build an NRP variant for this compound, so we skip it
                if iso_class is not None:
                    self.bgc_to_nrp_iso_classes[bgc_id].add(iso_class)

        print(f'{len(self.bgc_to_nrp_iso_classes)} BGCs in total (MIBiG 3+4)')

        self.match_correct_dict = get_match_correct_dict(bgc_to_nrps, self.similarity_dict)
        training_bgcs_info = self.mibig_bgcs_info.filter(
            pl.col(MIBiG_BGCs_Info.BGC_ID).is_in(self.training_bgcs)
        )

        print(f'{len(set(training_bgcs_info[MIBiG_BGCs_Info.BGC_ID].to_list()))} training BGCs (MIBiG 3+4)')
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
        # with open(nerpa_dir / 'tmp/test_nrps.txt', 'w') as f:
        #     for nrp in sorted(test_nrps):
        #         f.write(f'{nrp}: {self.nrp_id_to_iso_class[nrp]}' + '\n')

        self.test_nrp_classes = set(
            self.nrp_id_to_iso_class[nrp_id]
            for nrp_id in test_nrps
        )

        _test_bgcs = set(
            test_bgcs_info
            .filter(pl.col(MIBiG_BGCs_Info.NUM_A_DOMAINS) >= 3)
            [MIBiG_BGCs_Info.BGC_ID]
        )

        # Keep only BGCs that have at least one NRP iso-class in the test NRP classes
        self.test_bgcs = {bgc for bgc in _test_bgcs
                          if any(nrp_class in self.test_nrp_classes
                                 for nrp_class in self.bgc_to_nrp_iso_classes[bgc])}
        self.test_bgcs -= {'BGC0001561'}  # a mistake in MIBiG

        print(f'Number of test BGCs: {len(self.test_bgcs)}')

        # with open(nerpa_dir / 'tmp/test_bgcs.txt', 'w') as f:
        #     for bgc_id in sorted(self.test_bgcs):
        #         f.write(f'{bgc_id}: {",".join(self.bgc_to_nrp_iso_classes[bgc_id])}' + '\n')


        self.nrp_classes_with_matches = set(
            chain.from_iterable(
                self.bgc_to_nrp_iso_classes[bgc_id]
                for bgc_id in self.test_bgcs
            )
        )

    def match_is_correct(self, nrp_iso_class: str, bgc_id: str,
                         cmp_mode: str = PCS.NERPA_ISO_WEAK_CHR) -> bool:
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
                                 cmp_method: str = PCS.NERPA_ISO_WEAK_CHR) -> pl.Series:
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

    def load_nerpa_report(
            self,
            report_path: Path,
            report_name: str = 'REPORT_NAME_MISSING',
            score_column: str = NerpaReport.SCORE,
            tool_version: Literal['Nerpa 1', 'Nerpa 2'] = 'Nerpa 2'
    ) -> NerpaReport:
        nrp_classes_to_keep = self.test_nrp_classes
        bgc_ids_to_keep = self.test_bgcs


        out_cfg = OutputSizeConfig(
            max_num_matches=10000,  # no limit on total matches
            max_num_matches_per_bgc=10,
            max_num_matches_per_nrp=10000,  # no limit
            min_num_matches_per_bgc=10,
            min_num_matches_per_nrp=0,
        )

        report = load_nerpa_report(
            report_path=report_path,
            bgc_ids_to_keep=bgc_ids_to_keep,
            nrp_classes_to_keep=nrp_classes_to_keep,
            nrp_id_to_iso_class=self.nrp_id_to_iso_class,
            out_size_cfg=out_cfg,
            report_name=report_name,
            score_column=score_column,
            tool_version=tool_version
        )

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

        # q: take first five BGCs in the report and write all the corresponding rows to a tsv
        first_bgcs = (
            report.select(NerpaReport.BGC_ID)
            .unique(maintain_order=True)
            .head(5)[NerpaReport.BGC_ID]
            .to_list()
        )
        if first_bgcs:
            safe_report_name = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in report_name)
            out_path = Path(f"{safe_report_name}_first5_bgcs.tsv")
            # report.filter(pl.col(NerpaReport.BGC_ID).is_in(first_bgcs)).write_csv(out_path, separator="\t")
            # q: write the report to out_path. Sort the rows by BGC_ID
            report.filter(pl.col(NerpaReport.BGC_ID).is_in(first_bgcs)).sort(by=NerpaReport.BGC_ID).write_csv(out_path, separator="\t")

        for cmp_method in [
            PCS.NERPA_ISO_WEAK_CHR,
            #PCS.NERPA_NO_MORE_ONE_SUB_ALLOW_UNK_CHR,
        ]:
            report = report.with_columns(
                pl.struct([NerpaReport.NRP_ISO_CLASS, NerpaReport.BGC_ID])
                .map_elements(lambda row: self.match_is_correct(nrp_iso_class=row[NerpaReport.NRP_ISO_CLASS],
                                                                bgc_id=row[NerpaReport.BGC_ID],
                                                                cmp_mode=cmp_method),
                             return_dtype=pl.Boolean)
                .alias(NerpaReport.is_correct_col(cmp_method))
            )

        # # check that if NerpaReport.IS_CORRECT is true, then IS_CORRECT_COL is also true
        # incorrect_matches = report.filter(
        #     pl.col(NerpaReport.is_correct_col(PCS.NERPA_EQUAL_ALLOW_UNK_CHR))
        #     &
        #     ~pl.col(NerpaReport.is_correct_col(PCS.NERPA_NO_MORE_ONE_SUB_ALLOW_UNK_CHR))
        # )
        # if incorrect_matches.height > 0:
        #     row = incorrect_matches.row(0, named=True)
        #     raise ValueError(f"Inconsistent correctness columns\n"
        #                      f"Row: {row}\n")

        # add row index
        report = report.with_row_count(name=NerpaReport.MATCH_RANK)

        return NerpaReport(report, report_name=report_name)

    def compute_precision_recall_curve(self,
                                       nerpa_report: NerpaReport,
                                       top_matches_per_bgc: Optional[int],
                                       cmp_method: str = PCS.NERPA_ISO_WEAK_CHR) -> List[PrecisionRecallPoint]:
        return compute_precision_recall_curve(self, nerpa_report, top_matches_per_bgc, cmp_method)

    def extra_false_positives(self,
                              report1: NerpaReport,
                              report2: NerpaReport) -> pl.DataFrame:
        """Count extra false positives in report1 compared to report2."""
        return extra_false_positives(self, report1, report2)
