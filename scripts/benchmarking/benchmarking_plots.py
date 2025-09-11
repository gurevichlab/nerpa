

def cnts_to_percentages(cnts: pd.Series[int]) -> pd.Series[float]:
    """
    Convert counts to percentages.
    """
    percentages = cnts[1:] / np.arange(1, len(cnts)) * 100.0  # skip the first value (0 identified at N=0)
    percentages.index = pd.RangeIndex(start=1, stop=len(percentages) + 1)  # adjust index to match the percentage calculation
    if not ((percentages > -1e-9) & (percentages < 100 + 1e-9)).all():
        raise ValueError("Percentages should be in the range [0, 100]")
    return percentages

class PlotsHelper:
    def __init__(self):
        self.data_helper = PlotsDataHelper()

    def plot_num_correct_matches(self,
                                 nerpa_reports: Dict[str, NerpaReport],
                                 output_dir: Path,
                                 in_one_plot: bool = True,
                                 score_columns: Optional[Dict[str, Literal['LogOdds score', 'p_value']]] = None,
                                 y_axis: Literal['Count', 'Percentage'] = 'Count') -> List[Path]:
        if score_columns is None:
            score_columns = {report_name: 'LogOdds score' for report_name in nerpa_reports.keys()}

        output_dir.mkdir(parents=True, exist_ok=True)
        graphs = {report_name: self.data_helper.compute_num_correct_matches(nerpa_report,
                                                                            sort_by=score_columns[report_name])
                  for report_name, nerpa_report in nerpa_reports.items()}

        def _plot_single_graph(ax: Axes,
                               _values: pd.Series,
                               title: str,
                               label: Optional[str] = None) -> None:
            if y_axis == 'Percentage':
                values = cnts_to_percentages(_values)
            else:
                values = _values
            if label:
                ax.plot(values.index, values, label=label)
            else:
                ax.plot(values.index, values)
            ax.set_xlabel('Number of matches')
            ax.set_ylabel(f'{y_axis} of correct matches')
            ax.set_title(title)
            ax.grid()
            ax.set_aspect('equal', adjustable='box')
            ax.set_ylim(bottom=0)

        if in_one_plot:
            out_file = output_dir / f'{y_axis}_correct_all.png'

            fig, ax = plt.subplots(figsize=(10, 6))
            for report_name, num_correct_array in graphs.items():
                _plot_single_graph(ax, num_correct_array, '', label=report_name)

            ax.set_title(f'{y_axis} of Correct Matches among the top N vs N')
            ax.legend()
            fig.tight_layout()
            fig.savefig(out_file)
            plt.close(fig)
            return [out_file]
        else:
            out_files = {report_name: output_dir / f'{y_axis}_correct_{report_name}.png'
                         for report_name in graphs.keys()}

            for report_name, num_correct_array in graphs.items():
                fig, ax = plt.subplots(figsize=(10, 6))
                _plot_single_graph(
                    ax,
                    num_correct_array,
                    f'{y_axis} of correct matches {report_name}'
                )
                fig.tight_layout()
                fig.savefig(out_files[report_name])
                plt.close(fig)

            return list(out_files.values())
        
    def plot_score_correctness(self,
                               nerpa_reports: Dict[str, NerpaReport],
                               score_columns: Dict[str, Literal['LogOdds score', 'p_value']],
                               output_dir: Path,
                               num_bins: int = 20) -> List[Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_files = {report_name: output_dir / f'score_correctness_{report_name}.png'
                     for report_name in nerpa_reports.keys()}

        for report_name, nerpa_report in nerpa_reports.items():
            correct_incorrect_counts = self.data_helper.compute_score_correctness(nerpa_report, 
                                                                                  score_columns[report_name],
                                                                                  num_bins)

            ax = correct_incorrect_counts[["Correct", "Incorrect"]].plot(
                kind="bar",
                stacked=True,
                figsize=(10, 4),
                width=0.9,
                rot=45,
                xlabel="Score range",
                ylabel="Number of matches",
                title=f"{report_name}: Correct vs. incorrect matches per score ({score_columns[report_name]}) bin (N={num_bins})"
            )
            ax.set_xticklabels(
                [f"{iv.left:.2f}–{iv.right:.2f}" for iv in correct_incorrect_counts.index],
                rotation=45, ha="right"
            )

            plt.tight_layout()
            plt.savefig(out_files[report_name])
            plt.close()
    
        return list(out_files.values())

    def plot_num_identified(self,
                            nerpa_reports: Dict[str, 'NerpaReport'],
                            id_column: Literal['Genome_ID', 'NRP_ID'],
                            output_dir: Path,
                            score_columns: Dict[str, Literal['LogOdds score', 'p_value']],
                            y_axis: Literal['Count', 'Percentage'] = 'Count',
                            top_ks: Tuple[int, ...] = (1,)) -> List[Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_files_per_top_k = {top_k: output_dir / f'{y_axis}_identified_{id_column}_top{top_k}.png'
                               for top_k in top_ks}
        out_files_per_id = {report_name: output_dir / f'{y_axis}_identified_{id_column}_{report_name}.png'
                            for report_name in nerpa_reports.keys()}

        graphs = {
            report_name: {
                top_k: self.data_helper.compute_num_identified(
                    nerpa_report,
                    id_column,
                    score_columns[report_name],
                    top_k
                )
                for top_k in top_ks
            }
            for report_name, nerpa_report in nerpa_reports.items()
        }

        # Plot per report (all top_ks)
        for report_name, topk_results in graphs.items():
            fig, ax = plt.subplots(figsize=(10, 6))
            for top_k, _values in topk_results.items():
                if y_axis == 'Percentage':
                    values = cnts_to_percentages(_values)
                else:
                    values = _values
                ax.plot(values.index, values.values, label=f"top-{top_k}")
            ax.set_title(f"Num Identified by {id_column} using {score_columns[report_name]}: {report_name}")
            ax.set_xlabel(f"Num top {id_column}")
            ax.set_ylabel(f"{y_axis} identified")
            ax.grid()
            ax.legend()
            fig.tight_layout()
            fig.savefig(out_files_per_id[report_name])
            plt.close(fig)

        # Plot per top_k (all reports)
        for top_k in top_ks:
            fig, ax = plt.subplots(figsize=(10, 6))
            for report_name, topk_results in graphs.items():
                _values = topk_results[top_k]
                if y_axis == 'Percentage':
                    values = cnts_to_percentages(_values)
                else:
                    values = _values
                ax.plot(values.index, values.values, label=f'{report_name} ({score_columns[report_name]})')
            ax.set_title(f"{y_axis} Identified {id_column}: top-{top_k}")
            ax.set_xlabel("Identifier")
            ax.set_ylabel(y_axis)
            ax.grid()
            ax.legend()
            fig.tight_layout()
            fig.savefig(out_files_per_top_k[top_k])
            plt.close(fig)
        
        return list(chain(out_files_per_top_k.values(), out_files_per_id.values()))

    def plot_total_identified(self,
                              nerpa_reports: Dict[str, NerpaReport],
                              id_column: Literal['Genome_ID', 'NRP_ID'],
                              output_dir: Path,
                              score_columns: Dict[str, Literal['LogOdds score', 'p_value']],
                              max_top_k: int = 10,
                              y_axis: Literal['Count', 'Percentage'] = 'Count') -> List[Path]:
        """
        Plot the total number of identified BGCs/NRPs for each report
        """

        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / f'total_{y_axis}_identified_{id_column}.png'
        graphs = {report_name: self.data_helper.compute_total_identified(
            nerpa_report,
            id_column,
            score_columns[report_name],
            max_top_k=max_top_k,
            y_axis=y_axis
        ) for report_name, nerpa_report in nerpa_reports.items()}

        fig, ax = plt.subplots(figsize=(10, 6))
        for report_name, total_identified in graphs.items():
            ax.plot(total_identified.index,
                    total_identified.values,
                    label=f"{report_name}")
            ax.set_xticks(total_identified.index)

        ax.set_title(f'{y_axis} of identified {"BGCs" if id_column == "Genome_ID" else "NRPs"}',
                     fontsize=22)
        ax.set_xlabel(f'Num top matches considered for each {"BGC" if id_column == "Genome_ID" else "NRP"}',
                      fontsize=18)
        ax.set_ylabel(f'{y_axis} Identified', fontsize=18)
        ax.grid()
        ax.legend(fontsize=18)
        fig.tight_layout()
        fig.savefig(out_file)
        plt.close(fig)

        return [out_file]

    def plot_all(self,
                    nerpa_reports: Dict[str, NerpaReport],
                    output_dir: Path,
                    num_score_bins: int = 20,
                    top_ks: Tuple[int, ...] = (1, 3, 5, 10)) -> List[Path]:
            """
            Generate all plots and return paths to the generated files.
            """
            output_dir.mkdir(parents=True, exist_ok=True)
            score_columns = {report_name: 'LogOdds score' for report_name in nerpa_reports.keys()}
            plot_files = []

            plot_files.extend(self.plot_score_correctness(nerpa_reports,
                                                          score_columns,
                                                          output_dir,
                                                          num_score_bins))
            for y_axis in ['Count', 'Percentage']:
                plot_files.extend(self.plot_num_correct_matches(nerpa_reports,
                                                                y_axis=y_axis,
                                                                output_dir=output_dir))
                for id_column in ('Genome_ID', 'NRP_ID'):
                    plot_files.extend(self.plot_num_identified(nerpa_reports,
                                                               id_column=id_column,
                                                               y_axis=y_axis,
                                                               output_dir=output_dir,
                                                               score_columns=score_columns,
                                                               top_ks=top_ks))
                    plot_files.extend(self.plot_total_identified(nerpa_reports,
                                                                 id_column=id_column,
                                                                 y_axis=y_axis,
                                                                 output_dir=output_dir,
                                                                 score_columns=score_columns))

            return plot_files
    


def nerpa1_report_to_nerpa2_compatible(nerpa1_report: pd.DataFrame,
                                       genome_ids_to_filter: Optional[List[str]] = None) -> NerpaReport:
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
    nerpa2_report['LogOdds score'] = nerpa1_report['Score']
    nerpa2_report['NRP_ID'] = nerpa1_report['StructureId'].apply(lambda s: s.split('_')[0])  # Extract NRP ID from 'Structure ID'
    nerpa2_report['Genome_ID'] = nerpa1_report['PredictionFileName'].apply(pred_file_name_to_genome_id)  # Extract Genome ID from 'BGC ID'

    '''
    # q: group by NRP_ID and Genome_ID, and take the first row with the highest LogOdds score
    nerpa2_report = (nerpa2_report
                     .sort_values(by='LogOdds score', ascending=False)
                     .drop_duplicates(subset=['NRP_ID', 'Genome_ID'], keep='first')
                     .reset_index(drop=True))

    if genome_ids_to_filter is not None:
        # q: filter by genome IDs
        nerpa2_report = nerpa2_report[nerpa2_report['Genome_ID'].isin(genome_ids_to_filter)]

    # q: leave only 10 best matches per Genome_ID
    nerpa2_report = (nerpa2_report
                     .groupby('Genome_ID')
                     .head(10)
                     .reset_index(drop=True))
    '''

    return nerpa2_report


def nerpa1_vs_nerpa2():
    nerpa2_report_tsv = Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine/report.tsv')
    nerpa2_report = pd.read_csv(nerpa2_report_tsv, sep='\t')
    nerpa2_new_report_tsv = Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine2/report.tsv')
    nerpa2_new_report = pd.read_csv(nerpa2_new_report_tsv, sep='\t')

    nerpa1_report_csv = Path('/data/for_training_and_testing/nerpa1_report_mibig4_vs_mibig_norine.csv')
    nerpa1_report = pd.read_csv(nerpa1_report_csv)

    helper = PlotsHelper()
    mibig_nrps_table = helper.data_helper.mibig_bgcs_info
    bgc_ids = mibig_nrps_table[mibig_nrps_table['in_approved_matches'] == True]['bgc_id'].unique()

    nerpa1_report = nerpa1_report_to_nerpa2_compatible(nerpa1_report,
                                                       genome_ids_to_filter=bgc_ids.tolist())

    print(nerpa1_report.head())
    print(nerpa2_report.head())

    """
    assert set(nerpa1_report['Genome_ID']) == set(nerpa2_report['Genome_ID']), \
        ("Genome IDs in Nerpa1 and Nerpa2 reports do not match\n"
         f'nerpa1 \ nerpa2 = {set(nerpa1_report["Genome_ID"]) - set(nerpa2_report["Genome_ID"])}\n'
         f'nerpa2 \ nerpa1 = {set(nerpa2_report["Genome_ID"]) - set(nerpa1_report["Genome_ID"])}')
    """
    output_dir = Path('/home/ilianolhin/git/nerpa2/benchmarking/nerpa1_vs_nerpa2_plots')
    helper.plot_all({'Nerpa 1': nerpa1_report,
                     'Nerpa 2': nerpa2_report,
                     'Nerpa 2 new': nerpa2_new_report},
                    output_dir=output_dir)

    def bgc_id_identified(report: NerpaReport, genome_id: str) -> bool:
        """
        Check if a BGC with the given genome_id is identified in the report.
        """
        return (report[report['Genome_ID'] == genome_id]
                .apply(helper.data_helper.match_is_correct, axis=1)
                .any())

    # q: print all Genome_IDs in identified in nerpa1_report but not in nerpa2_report
    identified_in_nerpa1 = {bgc_id for bgc_id in set(nerpa1_report['Genome_ID'])
                            if bgc_id_identified(nerpa1_report, bgc_id)}
    identified_in_nerpa2 = {bgc_id for bgc_id in set(nerpa2_report['Genome_ID'])
                            if bgc_id_identified(nerpa2_report, bgc_id)}
    print('Identified in Nerpa1 but not in Nerpa2:'
          f'\n{identified_in_nerpa1 - identified_in_nerpa2}')



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
    nerpa2_new_vs_avg_nrp = helper.data_helper.load_nerpa_report(Path('/home/ilianolhin/git/nerpa2/nerpa_results/approved_vs_mibig_norine_new/report.tsv'),
                                                                 score_column='LogOdds_vs_avg_NRP',
                                                                 report_name='Nerpa 2 new')

    # 2. Check for reports discrepancies
    helper.data_helper.check_report_discrepancy(nerpa1_report, nerpa2_report,
                                                'Nerpa 1', 'Nerpa 2',
                                                check_nrp_ids=False)
    helper.data_helper.check_report_discrepancy(nerpa2_report, nerpa2_new_vs_avg_bgc,
                                                'Nerpa 2', 'Nerpa 2 new vs avg BGC',
                                                check_nrp_ids=False)


    # 4. Generate plots
    print('Generating plots...')
    output_dir = Path('/home/ilianolhin/git/nerpa2/benchmarking/nerpa1_vs_nerpa2_plots')
    '''
    helper.plot_all({'Nerpa 1': nerpa1_report,
                     'Nerpa 2': nerpa2_report,
                     'Nerpa 2 score vs avg NRP': nerpa2_new_vs_avg_nrp,
                     'Nerpa 2 score vs avg BGC': nerpa2_new_vs_avg_bgc},
                    output_dir=output_dir)
    '''
    helper.plot_all({'Nerpa 1': nerpa1_report,
                     'Nerpa 2 score vs avg BGC': nerpa2_new_vs_avg_bgc},
                    output_dir=output_dir)


def cross_validation():
    benchmarking_results_dir = Path('/home/ilianolhin/git/nerpa2/benchmarking/cross_validation_results')
    nerpa_reports = {
        f'subset_{i}': pd.read_csv(benchmarking_results_dir / f'subset_{i}' / 'report.tsv', sep='\t')
        for i in range(5)  # Assuming you have 5 samples
    }
    output_dir = Path(benchmarking_results_dir / 'plots')

    helper = PlotsHelper()
    helper.plot_all(nerpa_reports, output_dir)


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