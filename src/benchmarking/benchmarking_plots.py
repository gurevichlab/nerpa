from __future__ import annotations
from pathlib import Path

from src.benchmarking.plots_helper import PlotsHelper
from PIL import Image


def join_pngs_side_by_side(output_path, *image_paths):
    # Load all images
    images = [Image.open(p) for p in image_paths]

    # Find the smallest height
    min_height = min(img.height for img in images)

    # Rescale all images to the minimum height (preserving aspect ratio)
    resized = []
    for img in images:
        ratio = min_height / img.height
        new_width = int(img.width * ratio)
        resized.append(img.resize((new_width, min_height), Image.LANCZOS))

    # Compute total width after resizing
    total_width = sum(img.width for img in resized)

    # Create new blank image
    result = Image.new("RGBA", (total_width, min_height))

    # Paste images side by side
    x_offset = 0
    for img in resized:
        result.paste(img, (x_offset, 0))
        x_offset += img.width

    # Save
    result.save(output_path)
    print(f"Saved: {output_path}")


def nerpa1_vs_nerpa2():
    helper = PlotsHelper(bgc_test_set='mibig4_wo_training_bgcs')

    nerpa1_report = helper.data_helper.load_nerpa_report(Path(
        '/data/for_training_and_testing/nerpa1_report_mibig4_vs_mibig_norine.csv'),
                                                         tool_version='Nerpa 1',
                                                         report_name='Nerpa 1')
    nerpa2_report = helper.data_helper.load_nerpa_report(Path('/nerpa_results/mibig4_vs_mibig_norine/report.tsv'),
                                                         report_name='Nerpa 2',
                                                         score_column='LogOdds_vs_avg_BGC')

    extra_fp = helper.data_helper.extra_false_positives(nerpa2_report, nerpa1_report)
    with open('/benchmarking/nerpa1_vs_nerpa2_plots/extra_false_positives.tsv', 'w') as f:
        extra_fp.write_csv(f, separator='\t')


    # output_dir = Path('/home/ilianolhin/git/nerpa2/benchmarking/nerpa1_vs_nerpa2_plots')
    # helper.plot_all([nerpa1_report, nerpa2_report],
    #                 output_dir=output_dir)

def nerpa1_vs_nerpa2_vs_biocat():
    #helper = PlotsHelper(bgc_test_set='mibig4_wo_training_bgcs')
    helper = PlotsHelper(bgc_test_set='mibig4_wo_training_bgcs')

    nerpa1_report = helper.data_helper.load_nerpa_report(Path(
        '/data/for_training_and_testing/nerpa1_report_mibig4_vs_mibig_norine.csv'),
                                                         tool_version='Nerpa 1',
                                                         report_name='Nerpa 1')
    nerpa2_report = helper.data_helper.load_nerpa_report(Path('/nerpa_results/mibig4_vs_mibig_norine/report.tsv'),
                                                         report_name='Nerpa 2',
                                                         score_column='LogOdds_vs_avg_BGC')
    biocat_report = helper.data_helper.load_nerpa_report(Path('/data/for_training_and_testing/biocat_results.txt'),
                                                         report_name='BioCAT',
                                                         tool_version='BioCAT')


    output_dir = Path('/benchmarking/nerpa1_vs_nerpa2_vs_biocat_plots')
    helper.plot_all([nerpa1_report, nerpa2_report, biocat_report],
                    output_dir=output_dir)


def plots_for_paper():
    nerpa_dir = Path(__file__).parent.parent.parent
    
    # nerpa1_vs_nerpa2_vs_biocat()
    # nerpa1_vs_nerpa2_alignment()
    plots_for_paper_dir = nerpa_dir / 'benchmarking' / 'plots_for_paper'
    plots_for_paper_dir.mkdir(exist_ok=True)

    # copy the following files to plots_for_paper_dir:
    # nerpa_dir / 'benchmarking' / 'nerpa1_vs_nerpa2_vs_biocat_plots' / 'precision_recall_curve_top_10.png'
    # nerpa_dir / 'benchmarking' / 'nerpa1_vs_nerpa2_vs_biocat_plots' / 'total_Percentage_identified_Genome_ID.png'
    # nerpa_dir / 'benchmarking' / 'nerpa1_vs_nerpa2_plots' / 'alignment_reconstruction_histogram.png'
    files_to_copy = [
        nerpa_dir / 'benchmarking' / 'nerpa1_vs_nerpa2_vs_biocat_plots' / 'total_Percentage_identified_Genome_ID.png',
        nerpa_dir / 'benchmarking' / 'nerpa1_vs_nerpa2_vs_biocat_plots' / 'Percentage_identified_Genome_ID_top10.png',
        nerpa_dir / 'benchmarking' / 'nerpa1_vs_nerpa2_vs_biocat_plots' / 'precision_recall_curve_top_10.png',
        nerpa_dir / 'benchmarking' / 'nerpa1_vs_nerpa2_plots' / 'alignment_reconstruction_histogram.png',
    ]
    for file_path in files_to_copy:
        dest_path = plots_for_paper_dir / file_path.name
        dest_path.write_bytes(file_path.read_bytes())
        
    join_pngs_side_by_side(plots_for_paper_dir / 'combined_figure.png',
                           *files_to_copy)


if __name__ == "__main__":
    # Example usage
    plots_for_paper()
    #nerpa1_vs_nerpa2()
    #nerpa1_vs_nerpa2_vs_biocat()
    #nerpa1_vs_nerpa2_vs_nerpa2new()
    #cross_validation()
    #log_odds_vs_p_values()
    #nerpa1_vs_nerpa2_mibig4_wo_training_bgcs()

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

'''
def cross_validation():
    benchmarking_results_dir = Path('/home/ilianolhin/git/nerpa2/benchmarking/cross_validation_results')
    nerpa_reports = {
        f'subset_{i}': pd.read_csv(benchmarking_results_dir / f'subset_{i}' / 'report.tsv', sep='\t')
        for i in range(5)  # Assuming you have 5 samples
    }
    output_dir = Path(benchmarking_results_dir / 'plots')

    helper = PlotsHelper()
    helper.plot_all(nerpa_reports, output_dir)
'''
