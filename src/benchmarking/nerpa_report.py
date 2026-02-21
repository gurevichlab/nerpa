from __future__ import annotations
from pathlib import Path
from typing import Literal, Dict, NamedTuple, Set

import polars as pl

from src.benchmarking.data_frames import PNRPDB_Info, PNRPDB_Compound_Similarity

PCS = PNRPDB_Compound_Similarity

class NerpaReport(pl.DataFrame):
    NRP_ISO_CLASS = 'nrp_iso_class_representative'
    NRP_ID = 'NRP_ID'
    BGC_ID = 'Genome_ID'
    SCORE = 'Score'
    IS_CORRECT = 'is_correct'
    MATCH_RANK_FOR_BGC = 'match_rank_for_bgc'
    MATCH_RANK_FOR_NRP = 'match_rank_for_nrp'
    MATCH_RANK = 'match_rank'

    name: str

    def __init__(self, *args, report_name: str = 'REPORT_NAME_MISSING', **kwargs):
        if len(args) == 1 and isinstance(args[0], pl.DataFrame) and not kwargs:
            # Convert to something older Polars constructors accept
            data = args[0].to_dict(as_series=False)  # dict[str, list]
            super().__init__(data)
        else:
            super().__init__(*args, **kwargs)
        self.name = report_name

    @classmethod
    def is_correct_col(cls, cmp_method: str) -> str:
        assert cmp_method in [PCS.NERPA_EQUAL_ALLOW_UNK_CHR,
                              PCS.NERPA_NO_MORE_ONE_SUB_ALLOW_UNK_CHR], \
            f"Unsupported comparison method: {cmp_method}"

        return f'{cls.IS_CORRECT}_{cmp_method}'


class OutputSizeConfig(NamedTuple):
    min_num_matches_per_bgc: int
    max_num_matches_per_bgc: int
    min_num_matches_per_nrp: int
    max_num_matches_per_nrp: int
    max_num_matches: int


def is_bgc_or_norine_nrp() -> pl.Expr:
    return (pl.col(PNRPDB_Info.COMPOUND_ID).str.starts_with('BGC') |
            pl.col(PNRPDB_Info.COMPOUND_ID).str.starts_with('NOR'))

def nerpa1_report_to_nerpa2_compatible(nerpa1_report: pl.DataFrame) -> pl.DataFrame:
    """Convert a Nerpa1 report to a Nerpa2 compatible format."""
    def pred_file_name_to_genome_id(pred_file_name: str) -> str:
        name = pred_file_name.split('/')[-1]
        return name.split('_')[0]

    return nerpa1_report.select([
        pl.col('Score')
        .alias(NerpaReport.SCORE),

        pl.col('StructureId').str
        .split('_').list
        .first()
        .alias(NerpaReport.NRP_ID),

        pl.col('PredictionFileName')
        .map_elements(pred_file_name_to_genome_id,
                      return_dtype=pl.Utf8)
        .alias(NerpaReport.BGC_ID)
    ])

def biocat_to_nerpa2_compatible(biocat_report: pl.DataFrame) -> pl.DataFrame:
        """Convert a Biocat report to a Nerpa2 compatible format."""
        # 1. Rename 'Relative score' -> NerpaReport.SCORE, 'Chromosome ID' -> NerpaReport.BGC_ID, 'Substance' -> NerpaReport.NRP_ID
        # 2. In case of multiple matches between the same BGC and NRP, keep only the one with the highest score
        report = biocat_report.rename({
            'Relative score': NerpaReport.SCORE,
            'Chromosome ID': NerpaReport.BGC_ID,
            'Substance': NerpaReport.NRP_ID
        })
        report = report.group_by([NerpaReport.BGC_ID, NerpaReport.NRP_ID]).agg(
            pl.col(NerpaReport.SCORE).max()
        )
        return report

def truncate_output(report: NerpaReport,
                    out_size_cfg: OutputSizeConfig) -> NerpaReport:
    """Truncate the report to satisfy the output size configuration."""
    name = report.name
    # Sort by score
    report = report.sort(NerpaReport.SCORE, descending=True)

    # Add counters for each BGC and NRP
    report = report.with_columns([
        pl.col(NerpaReport.BGC_ID)
        .cum_count()
        .over(NerpaReport.BGC_ID)
        .alias('_match_rank_for_bgc'),

        pl.col(NerpaReport.NRP_ID)
        .cum_count()
        .over(NerpaReport.NRP_ID)
        .alias('_match_rank_for_nrp')
    ])

    # For the case when out_size_cfg.max_num_matches_per_bgc = float('inf')
    max_num_matches = min(out_size_cfg.max_num_matches, len(report))

    # Filter based on configuration
    filtered = report.filter(
        (
            (pl.col('_match_rank_for_bgc') < out_size_cfg.max_num_matches_per_bgc)  # < because match_rank starts from 0
            &
            (pl.col('_match_rank_for_nrp') < out_size_cfg.max_num_matches_per_nrp)
        )
        |
        (pl.col('_match_rank_for_bgc') < out_size_cfg.min_num_matches_per_bgc)
        |
        (pl.col('_match_rank_for_nrp') < out_size_cfg.min_num_matches_per_nrp)
    ).head(max_num_matches)

    print(f"Truncated report {name} from {len(report)} matches to {len(filtered)} matches based on output size configuration")
    return NerpaReport(filtered.drop(['_match_rank_for_bgc', '_match_rank_for_nrp']),
                       report_name=name)


def prepare_data_for_benchmarking(_report: pl.DataFrame,
                                  bgc_ids_to_keep: Set[str],
                                  nrp_classes_to_keep: Set[str],
                                  nrp_id_to_iso_class: Dict[str, str],
                                  out_size_cfg: OutputSizeConfig,
                                  report_name: str = 'REPORT_NAME_MISSING') -> NerpaReport:
    # Filter BGCs and NRPs
    report = _report.filter(
        pl.col(NerpaReport.BGC_ID).is_in(list(bgc_ids_to_keep))
    )

    missing_bgcs = set(bgc_ids_to_keep) - set(report[NerpaReport.BGC_ID].unique())
    if missing_bgcs:
        print(f"Warning! Some BGC IDs are missing in the report {report_name}")
        print(f"Missing BGC IDs: {sorted(missing_bgcs)}")
    else:
        print(f'All {len(bgc_ids_to_keep)} BGC IDs are present in the report {report_name}')

    # Add NRP iso-class column
    report = report.with_columns(
        pl.col(NerpaReport.NRP_ID)
        .map_elements(lambda x: nrp_id_to_iso_class.get(x), return_dtype=pl.Utf8)
        .alias(NerpaReport.NRP_ISO_CLASS)
    )

    report = report.filter(
        pl.col(NerpaReport.NRP_ISO_CLASS).is_in(nrp_classes_to_keep)
    )

    # Sort by score
    report = report.sort(NerpaReport.SCORE, descending=True)

    # Remove duplicates
    report = report.unique(subset=[NerpaReport.BGC_ID, NerpaReport.NRP_ISO_CLASS],
                           keep='first')


    return truncate_output(NerpaReport(report, report_name=report_name), out_size_cfg)


def load_nerpa_report(report_path: Path,
                      bgc_ids_to_keep: Set[str],
                      nrp_classes_to_keep: Set[str],
                      nrp_id_to_iso_class: Dict[str, str],
                      out_size_cfg: OutputSizeConfig,
                      report_name: str = 'REPORT_NAME_MISSING',
                      score_column: str = NerpaReport.SCORE,
                      tool_version: Literal['Nerpa 1', 'Nerpa 2', 'BioCAT'] = 'Nerpa 2') -> NerpaReport:
    print(f"Loading report from {report_path} with tool version {tool_version}...")
    match tool_version:
        case 'Nerpa 2':
            report = pl.read_csv(report_path, separator='\t')
        case 'Nerpa 1':
            report = pl.read_csv(report_path)
            report = nerpa1_report_to_nerpa2_compatible(report)
        case 'BioCAT':
            report = pl.read_csv(report_path, separator='\t')
            report = biocat_to_nerpa2_compatible(report)
        case _:
            raise ValueError(f"Unsupported tool version: {tool_version}")

    report = report.with_columns(
        pl.col(score_column)
        .alias(NerpaReport.SCORE)
    )

    return prepare_data_for_benchmarking(
        _report=report,
        report_name=report_name,
        bgc_ids_to_keep=bgc_ids_to_keep,
        nrp_classes_to_keep=nrp_classes_to_keep,
        nrp_id_to_iso_class=nrp_id_to_iso_class,
        out_size_cfg=out_size_cfg
    )


def check_reports_discrepancy(data_helper,
                              report1: NerpaReport,
                              report2: NerpaReport) -> None:
    """Check for discrepancies between two reports."""
    print(f'Checking report discrepancies between {report1.name} and {report2.name}...')

    bgc_ids_1 = set(report1[NerpaReport.BGC_ID])
    bgc_ids_2 = set(report2[NerpaReport.BGC_ID])

    if bgc_ids_1 != bgc_ids_2:
        print(f"Warning! Genome IDs in reports {report1.name} and {report2.name} do not match")
        print(f'{report1.name} \\ {report2.name} = {bgc_ids_1 - bgc_ids_2}')
        print(f'{report2.name} \\ {report1.name} = {bgc_ids_2 - bgc_ids_1}')
    else:
        print("Genome IDs sets match")

    for id_column in (NerpaReport.BGC_ID, NerpaReport.NRP_ISO_CLASS):
        identified_ids_1 = data_helper.get_identified_ids(report1, id_column)
        identified_ids_2 = data_helper.get_identified_ids(report2, id_column)
        if set(identified_ids_1) != set(identified_ids_2):
            print(f"Warning! Identified {id_column}s in reports {report1.name} and {report2.name} do not match")
            print(f'Identified in {report1.name} but not in {report2.name}:\n{set(identified_ids_1) - set(identified_ids_2)}')
            print(f'Identified in {report2.name} but not in {report1.name}:\n{set(identified_ids_2) - set(identified_ids_1)}')
        else:
            print(f'Identified {id_column}s match')
