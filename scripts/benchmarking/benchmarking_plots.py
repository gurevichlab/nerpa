from __future__ import annotations
from itertools import chain
from pathlib import Path
from typing import Literal, Dict, Optional, Sequence, List, Tuple

import numpy as np
import pandas as pd
from scripts.benchmarking.plots_helper import PlotsHelper



def nerpa1_vs_nerpa2():
    helper = PlotsHelper()

    nerpa1_report = helper.data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/data/for_training_and_testing/nerpa1_report_mibig4_vs_mibig_norine.csv'),
                                                         nerpa_version='Nerpa 1',
                                                         report_name='Nerpa 1')
    nerpa2_report = helper.data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine/report.tsv'),
                                                         report_name='Nerpa 2',
                                                         score_column='LogOdds_score')


    output_dir = Path('/home/ilianolhin/git/nerpa2/benchmarking/nerpa1_vs_nerpa2_plots')
    helper.plot_all([nerpa1_report, nerpa2_report],
                    output_dir=output_dir)



def nerpa1_vs_nerpa2_vs_nerpa2new():
    helper = PlotsHelper()

    # 1. Load reports
    nerpa1_report = helper.data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/data/for_training_and_testing/nerpa1_report_mibig4_vs_mibig_norine.csv'),
                                                         nerpa_version='Nerpa 1',
                                                         report_name='Nerpa 1')
    nerpa2_report = helper.data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine/report.tsv'),
                                                         report_name='Nerpa 2')
    nerpa2_new_vs_avg_bgc = helper.data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine_new/report.tsv'),
                                                                 score_column='LogOdds_vs_avg_BGC',
                                                                 report_name='Nerpa 2 new')

    # 2. Check for reports discrepancies
    helper.data_helper.check_report_discrepancy(nerpa1_report, nerpa2_report)
    helper.data_helper.check_report_discrepancy(nerpa2_report, nerpa2_new_vs_avg_bgc)


    # 4. Generate plots
    print('Generating plots...')
    output_dir = Path('/home/ilianolhin/git/nerpa2/benchmarking/nerpa1_vs_nerpa2_plots')
    helper.plot_all([nerpa1_report, nerpa2_report, nerpa2_new_vs_avg_bgc],
                    output_dir=output_dir)




if __name__ == "__main__":
    # Example usage
    #nerpa1_vs_nerpa2()
    nerpa1_vs_nerpa2_vs_nerpa2new()
    #cross_validation()
    #log_odds_vs_p_values()

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
