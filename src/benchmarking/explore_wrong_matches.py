from src.benchmarking.plots_data_helper import PlotsDataHelper
from src.benchmarking.nerpa_report import NerpaReport
from src.benchmarking.data_frames import PNRPDB_Compound_Similarity
from pathlib import Path
import polars as pl


def main():
    data_helper = PlotsDataHelper()

    # nerpa1_report = data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/data/for_training_and_testing/nerpa1_report_mibig4_vs_mibig_norine.csv'),
    #                                                      tool_version='Nerpa 1',
    #                                                      report_name='Nerpa 1')
    nerpa2_report = data_helper.load_nerpa_report(Path('/nerpa_results/approved_vs_mibig_norine_all/report.tsv'),
                                                  report_name='Nerpa 2.1',
                                                  score_column='LogOdds_vs_avg_BGC')
    print(nerpa2_report.head())

    nerpa2_report = nerpa2_report.with_columns(
        pl.struct([NerpaReport.BGC_ID, NerpaReport.NRP_ID])
        .map_elements(lambda s: data_helper.match_is_correct(
            bgc_id=s[NerpaReport.BGC_ID],
            nrp_iso_class=s[NerpaReport.NRP_ID],
            cmp_mode=PNRPDB_Compound_Similarity.NERPA_NO_MORE_ONE_SUB_ALLOW_UNK_CHR,
        ))
        .alias("is_correct_tmp")
    )

    incorrect_matches = (
        nerpa2_report
        .filter(pl.col("is_correct_tmp") == False)
        .select([
            NerpaReport.SCORE,
            NerpaReport.BGC_ID,
            NerpaReport.NRP_ID,
            NerpaReport.MATCH_RANK,
            NerpaReport.MATCH_RANK_FOR_BGC,
            NerpaReport.MATCH_RANK_FOR_NRP,
        ])
    )

    # q: write incorrect matches to tsv
    incorrect_matches.write_csv('incorrect_matches.tsv', separator='\t')


if __name__ == "__main__":
    main()
