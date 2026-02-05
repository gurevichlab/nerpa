from src.benchmarking.plots_data_helper import PlotsDataHelper
from src.benchmarking.nerpa_report import NerpaReport
from src.benchmarking.data_frames import PNRPDB_Compound_Similarity
from pathlib import Path
import polars as pl
import shutil


def main():
    data_helper = PlotsDataHelper()

    # nerpa1_report = data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/data/for_training_and_testing/nerpa1_report_mibig4_vs_mibig_norine.csv'),
    #                                                      tool_version='Nerpa 1',
    #                                                      report_name='Nerpa 1')
    nerpa2_report = data_helper.load_nerpa_report(Path('nerpa_results/mibig4_vs_mibig_norine/report.tsv'),
                                                  report_name='Nerpa 2.1',
                                                  score_column='LogOdds_vs_avg_BGC')
    # Compute table (ID -> best score)
    best_scores = (
        nerpa2_report
        .group_by(NerpaReport.BGC_ID)
        .agg(
            pl.col(NerpaReport.SCORE)
            .max()
            .alias('max_score')
        )
        .sort('max_score', descending=True)
    )

    # Compute table (ID -> identified) (if any of top_k matches is correct)
    identified = (
        nerpa2_report
        .sort(NerpaReport.SCORE, descending=True)
        .group_by(NerpaReport.BGC_ID)
        .agg(
            pl.col(NerpaReport.is_correct_col(PNRPDB_Compound_Similarity.NERPA_EQUAL_ALLOW_UNK_CHR))
            .head(10)
            .any()
            .alias('identified')
        )
    )

    # print number of not identified
    num_not_identified = (
        identified
        .filter(pl.col('identified') == False)
        .height
    )
    print(f"Number of not identified BGCs: {num_not_identified}")
    # print total number of bgcs
    total_bgcs = identified.height
    print(f"Total number of BGCs: {total_bgcs}")

    # get top 10 not identified bgc ids
    num_not_identified_to_inspect = 10
    not_identified_bgc_ids = (
        identified
        .filter(pl.col('identified') == False)
        .join(best_scores, on=NerpaReport.BGC_ID)
        .sort('max_score', descending=True)
        .select(NerpaReport.BGC_ID)
        .head(num_not_identified_to_inspect)
        .to_series()
    ).to_list()

    print(f"Top {num_not_identified_to_inspect} not identified BGC IDs:\n"
          '\n'.join(not_identified_bgc_ids))
    exit(0)

    nerpa2_report = nerpa2_report.with_columns(
        pl.struct([NerpaReport.BGC_ID, NerpaReport.NRP_ID])
        .map_elements(lambda s: data_helper.match_is_correct(
            bgc_id=s[NerpaReport.BGC_ID],
            nrp_iso_class=s[NerpaReport.NRP_ID],
            cmp_mode=PNRPDB_Compound_Similarity.NERPA_EQUAL_ALLOW_UNK_CHR,
        ))
        .alias("is_correct_tmp"),
        
        # iso classes for each bgc
        pl.col(NerpaReport.BGC_ID)
        .map_elements(lambda bgc_id: ','.join(data_helper.bgc_to_nrp_iso_classes.get(bgc_id)))
        .alias("true_nrp_ids"),
    )

    incorrect_matches = (
        nerpa2_report
        .filter(pl.col("is_correct_tmp") == False)
        .select([
            NerpaReport.SCORE,
            NerpaReport.BGC_ID,
            NerpaReport.NRP_ID,
            NerpaReport.MATCH_RANK,
            'true_nrp_ids'
        ])
    )

    # q: write incorrect matches to tsv
    out_file = Path('incorrect_matches.tsv')
    incorrect_matches.write_csv(out_file, separator='\t')
    print(f'The incorrect matches are written to in {out_file.resolve()}')

    # paths
    nerpa_dir = Path(__file__).parent.parent.parent
    debug_root = nerpa_dir / "debug" / "incorrect_matches"
    molecules_dir = nerpa_dir / "nerpa_results" / "pnrpdb2_vs_antismash_db5" / "NRP_images" / "molecules"
    graphs_dir = nerpa_dir / "nerpa_results" / "pnrpdb2_vs_antismash_db5" / "NRP_images" / "graphs"

    debug_root.mkdir(parents=True, exist_ok=True)

    # take first 10 incorrect rows
    num_rows_to_inspect = 10
    rows = incorrect_matches.head(num_rows_to_inspect).to_dicts()

    def copy_png(src: Path, dst: Path):
        if not src.exists():
            print(f"[WARN] Missing image: {src}")
            return
        shutil.copy2(src, dst)

    for i, r in enumerate(rows):
        wrong_nrp_id = r[NerpaReport.NRP_ID]
        true_nrp_ids = r["true_nrp_ids"].split(",")

        # folder per incorrect match
        folder = debug_root / f"{i}"
        folder.mkdir(parents=True, exist_ok=True)

        # ---- wrong ----
        wrong_mol_src = molecules_dir / f"{wrong_nrp_id}.png"
        wrong_graph_src = graphs_dir / f"{wrong_nrp_id}.png"
        wrong_mol_dst = folder / f"{wrong_nrp_id}_molecule_wrong.png"
        wrong_graph_dst = folder / f"{wrong_nrp_id}_graph_wrong.png"
        copy_png(wrong_mol_src, wrong_mol_dst)
        copy_png(wrong_graph_src, wrong_graph_dst)

        # ---- true (possibly multiple) ----
        for true_id in true_nrp_ids:
            true_id = true_id.strip()
            true_mol_src = molecules_dir / f"{true_id}.png"
            true_graph_src = graphs_dir / f"{true_id}.png"
            true_mol_dst = folder / f"{true_id}_molecule_true.png"
            true_graph_dst = folder / f"{true_id}_graph_true.png"
            copy_png(true_mol_src, true_mol_dst)
            copy_png(true_graph_src, true_graph_dst)

    print(f"True and wrong NRP molecules for top {num_rows_to_inspect} are saved to {debug_root}")


if __name__ == "__main__":
    main()
