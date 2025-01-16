import argparse
from pathlib import Path
from typing import List, Tuple, Union, Set, Iterable, Optional
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Input group
    input_group = parser.add_argument_group("input", "Input options")
    input_group.add_argument("--nerpa2", action='append', 
                             help="Nerpa2 output dir (could be used multiple times)")
    input_group.add_argument("--nerpa1", type=str, default=None,
                             help="Nerpa1 output dir")
    input_group.add_argument("--exclude", type=str, default='',
                             help="Comma-separated list of MiBIG accessions to exclude.")
    input_group.add_argument("--max-rank", type=int, default=10,
                             help="Maximum rank to still consider the match correct.")
    
    # Output group
    output_group = parser.add_argument_group("output", "Output options")
    output_group.add_argument("--fig", type=str, default="./fdrplot_nerpa1_vs_nerpa2.png",
                              help="Path to PNG.")
    output_group.add_argument("--tsv-pref", type=str, default="./fps",
                              help="Path with a prefix for a TSV file to store false positives.")

    return parser.parse_args()


def parse_nerpa1_results(nerpa_dir: Path) -> pd.DataFrame:
    report = pd.read_csv(nerpa_dir / 'report.csv')
    report['AccStruct'] = report['StructureId'].map(lambda x: x.rsplit('_', 1)[0].split('.')[0])
    report['AccGen'] = report['PredictionFileName'].map(lambda x: x.rsplit('/', 1)[-1].split('_', 1)[0])
    return report

def parse_nerpa2_results(nerpa_dir: Path) -> pd.DataFrame:
    report = pd.read_csv(nerpa_dir / 'report.tsv', sep='\t').rename(columns={'Genome_ID': 'AccGen'})
    report['AccStruct'] = report['NRP_ID'].map(lambda x: x.split('.')[0])
    return report

def calculate_fdr_data(report_: pd.DataFrame,
                       rank: int,
                       exclude: Optional[Set[str]] = None,
                       keep_only: Optional[Set[str]] = None) -> Tuple[Iterable, pd.DataFrame]:
    report = report_.copy()
    if exclude is None:
        exclude = set()
    if keep_only is not None:
        report = report[report['AccGen'].isin(keep_only - exclude)]

    print('Total BGCs: ', report['AccGen'].nunique())

    report['Correct'] = report['AccStruct'] == report['AccGen']
    report['Rank'] = report.groupby('AccGen')['Score'].rank('average', ascending=False)

    mask = report['Rank'] < rank + 1
    mask &= ~report['AccGen'].isin(exclude)
    top = report.loc[mask].sort_values(['Correct', 'Score'], ascending=[False, False])\
        .groupby('AccGen').head(1).sort_values('Score', ascending=False)
    print(top.shape)
    top = report.sort_values(['Correct', 'Score'], ascending=[False, False]) \
        .groupby('AccGen').head(1).sort_values('Score', ascending=False)
    print()
    fdr = np.zeros((3, top.shape[0]))
    for i in range(top.shape[0]):
        n = i + 1
        _df = top.head(n)
        fdr[0, i] = n
        n_correct = _df['Correct'].sum()
        fdr[1, i] = 1 - (n_correct / n)
        fdr[2, i] = _df['Score'].min()
    
    return fdr, top

def sample_nerpa_results(nerpa_results: pd.DataFrame, sample_size: int,
                         exclude: Optional[Set[str]] = None) -> pd.DataFrame:
    if exclude is None:
        exclude = set()
    bgc_ids = list(set(nerpa_results['AccGen'].unique()) - exclude)
    bgc_ids = np.random.choice(bgc_ids, sample_size, replace=False)
    return nerpa_results[nerpa_results['AccGen'].isin(bgc_ids)]


def run_default(args):
    exclude = set(x.strip() for x in args.exclude.split(','))
    mibig3_ids = {'BGC0000064', 'BGC0000081', 'BGC0000186', 'BGC0000287', 'BGC0000288', 'BGC0000292', 'BGC0000293', 'BGC0000298', 'BGC0000302', 'BGC0000303', 'BGC0000304', 'BGC0000308', 'BGC0000311', 'BGC0000320', 'BGC0000323', 'BGC0000326', 'BGC0000327', 'BGC0000329', 'BGC0000330', 'BGC0000332', 'BGC0000333', 'BGC0000337', 'BGC0000341', 'BGC0000342', 'BGC0000345', 'BGC0000348', 'BGC0000350', 'BGC0000352', 'BGC0000353', 'BGC0000355', 'BGC0000360', 'BGC0000362', 'BGC0000371', 'BGC0000372', 'BGC0000373', 'BGC0000376', 'BGC0000378', 'BGC0000380', 'BGC0000394', 'BGC0000395', 'BGC0000396', 'BGC0000404', 'BGC0000406', 'BGC0000408', 'BGC0000409', 'BGC0000410', 'BGC0000414', 'BGC0000416', 'BGC0000418', 'BGC0000419', 'BGC0000421', 'BGC0000422', 'BGC0000428', 'BGC0000430', 'BGC0000438', 'BGC0000439', 'BGC0000443', 'BGC0000444', 'BGC0000446', 'BGC0000448', 'BGC0000454', 'BGC0000458', 'BGC0000465', 'BGC0000816', 'BGC0000818', 'BGC0000819', 'BGC0000949', 'BGC0000954', 'BGC0000958', 'BGC0000959', 'BGC0000960', 'BGC0000961', 'BGC0000962', 'BGC0000964', 'BGC0000966', 'BGC0000967', 'BGC0000968', 'BGC0000970', 'BGC0000973', 'BGC0000974', 'BGC0000975', 'BGC0000976', 'BGC0000977', 'BGC0000982', 'BGC0000983', 'BGC0000988', 'BGC0000989', 'BGC0000990', 'BGC0000991', 'BGC0000993', 'BGC0000994', 'BGC0000995', 'BGC0000996', 'BGC0000997', 'BGC0000999', 'BGC0001001', 'BGC0001010', 'BGC0001011', 'BGC0001012', 'BGC0001013', 'BGC0001015', 'BGC0001016', 'BGC0001017', 'BGC0001020', 'BGC0001021', 'BGC0001022', 'BGC0001023', 'BGC0001024', 'BGC0001025', 'BGC0001026', 'BGC0001035', 'BGC0001036', 'BGC0001037', 'BGC0001039', 'BGC0001040', 'BGC0001042', 'BGC0001043', 'BGC0001044', 'BGC0001048', 'BGC0001049', 'BGC0001051', 'BGC0001052', 'BGC0001053', 'BGC0001054', 'BGC0001055', 'BGC0001057', 'BGC0001059', 'BGC0001069', 'BGC0001070', 'BGC0001071', 'BGC0001084', 'BGC0001088', 'BGC0001098', 'BGC0001099', 'BGC0001100', 'BGC0001101', 'BGC0001103', 'BGC0001106', 'BGC0001108', 'BGC0001109', 'BGC0001110', 'BGC0001111', 'BGC0001114', 'BGC0001115', 'BGC0001116', 'BGC0001122', 'BGC0001123', 'BGC0001125', 'BGC0001136', 'BGC0001142', 'BGC0001147', 'BGC0001149', 'BGC0001153', 'BGC0001165', 'BGC0001178', 'BGC0001182', 'BGC0001187', 'BGC0001192', 'BGC0001193', 'BGC0001216', 'BGC0001218', 'BGC0001220', 'BGC0001249', 'BGC0001250', 'BGC0001255', 'BGC0001268', 'BGC0001290', 'BGC0001327', 'BGC0001328', 'BGC0001329', 'BGC0001331', 'BGC0001332', 'BGC0001342', 'BGC0001344', 'BGC0001351', 'BGC0001367', 'BGC0001368', 'BGC0001370', 'BGC0001371', 'BGC0001394', 'BGC0001425', 'BGC0001426', 'BGC0001427', 'BGC0001428', 'BGC0001429', 'BGC0001430', 'BGC0001431', 'BGC0001434', 'BGC0001435', 'BGC0001437', 'BGC0001442', 'BGC0001445', 'BGC0001449', 'BGC0001455', 'BGC0001457', 'BGC0001459', 'BGC0001460', 'BGC0001467', 'BGC0001477', 'BGC0001502', 'BGC0001519', 'BGC0001545', 'BGC0001553', 'BGC0001560', 'BGC0001563', 'BGC0001609', 'BGC0001635', 'BGC0001652', 'BGC0001667', 'BGC0001671', 'BGC0001672', 'BGC0001692', 'BGC0001702', 'BGC0001705', 'BGC0001711', 'BGC0001717', 'BGC0001718', 'BGC0001737', 'BGC0001738', 'BGC0001740', 'BGC0001741', 'BGC0001750', 'BGC0001751', 'BGC0001768', 'BGC0001780', 'BGC0001800', 'BGC0001801', 'BGC0001816', 'BGC0001820', 'BGC0001845', 'BGC0001857', 'BGC0001873', 'BGC0001881', 'BGC0001902', 'BGC0001919', 'BGC0001942', 'BGC0001970', 'BGC0001971', 'BGC0002014', 'BGC0002035', 'BGC0002037', 'BGC0002046', 'BGC0002050', 'BGC0002055', 'BGC0002058', 'BGC0002059', 'BGC0002066', 'BGC0002067', 'BGC0002068', 'BGC0002075', 'BGC0002082', 'BGC0002091', 'BGC0002093', 'BGC0002095', 'BGC0002114', 'BGC0002118', 'BGC0002119', 'BGC0002131', 'BGC0002135', 'BGC0002157', 'BGC0002158', 'BGC0002160', 'BGC0002164', 'BGC0002166', 'BGC0002169', 'BGC0002197', 'BGC0002198', 'BGC0002204', 'BGC0002205', 'BGC0002210', 'BGC0002230', 'BGC0002239', 'BGC0002242', 'BGC0002248', 'BGC0002250', 'BGC0002256', 'BGC0002259', 'BGC0002264', 'BGC0002268', 'BGC0002272', 'BGC0002273', 'BGC0002287', 'BGC0002296', 'BGC0002300', 'BGC0002304', 'BGC0002314', 'BGC0002340', 'BGC0002344', 'BGC0002367', 'BGC0002373', 'BGC0002378', 'BGC0002408', 'BGC0002412', 'BGC0002434', 'BGC0002437', 'BGC0002456', 'BGC0002468', 'BGC0002491', 'BGC0002492', 'BGC0002495', 'BGC0002505', 'BGC0002509', 'BGC0002542', 'BGC0002543', 'BGC0002546', 'BGC0002555', 'BGC0002565', 'BGC0002611', 'BGC0002613', 'BGC0002649', 'BGC0002735'}
    mibig_nrps = set(Path('/home/ilianolhin/programs/antismash7/mibig_only_nrps.txt').read_text().splitlines())

    fdr_data = {}
    tops = {}

    if args.nerpa1:
        k = 'Nerpa 1'
        fdr_data[k], tops[k] = calculate_fdr_data(
            parse_nerpa1_results(Path(args.nerpa1)),
            args.max_rank, exclude=mibig3_ids, keep_only=mibig_nrps
        )
    if args.nerpa2:
        for i, p in enumerate(args.nerpa2):
            k = f'Nerpa 2 run {i}' if len(args.nerpa2) > 1 else 'Nerpa 2'
            nerpa2_results = parse_nerpa2_results(Path(p))

            fdr_data[k], tops[k] = calculate_fdr_data(
                nerpa2_results,
                args.max_rank, exclude=mibig3_ids, keep_only=mibig_nrps
            )

    f, ax = plt.subplots(figsize=(8, 6))
    for k, v in fdr_data.items():
        color = 'blue' if 'Nerpa 2' in k else 'orange'
        _ = ax.plot(v[0], v[1], label=k, color=color)

    _ = ax.legend()
    _ = ax.set_ylim((0, 1))
    ax.set_xlim(left=0)
    _ = ax.locator_params(axis='y', nbins=10)
    _ = ax.locator_params(axis='x', nbins=14)
    _ = ax.set_xlabel('Number of BGCs')
    _ = ax.set_ylabel('FDR')
    _ = ax.grid()

    f.savefig(args.fig)
    for k, df in tops.items():
        df.loc[~df['Correct']].to_csv(f'{args.tsv_pref}_{"_".join(k.split())}.tsv', index=False, sep='\t')


def run_compare_mibig3_and_mibig4(args, num_samples: int = 100):
    mibig3_results = parse_nerpa2_results(Path(args.nerpa2[0]))
    mibig4_results = parse_nerpa2_results(Path(args.nerpa2[1]))

    exclude = set(x.strip() for x in args.exclude.split(','))
    exclude = {'BGC0000404', 'BGC0000435', 'BGC0001570', 'BGC0000403', 'BGC0001461', 'BGC0001999', 'BGC0002203', 'BGC0002336', 'BGC0001297', 'BGC0002752', 'BGC0000408', 'BGC0002653'}
    mibig3_ids = {'BGC0000064', 'BGC0000081', 'BGC0000186', 'BGC0000287', 'BGC0000288', 'BGC0000292', 'BGC0000293', 'BGC0000298', 'BGC0000302', 'BGC0000303', 'BGC0000304', 'BGC0000308', 'BGC0000311', 'BGC0000320', 'BGC0000323', 'BGC0000326', 'BGC0000327', 'BGC0000329', 'BGC0000330', 'BGC0000332', 'BGC0000333', 'BGC0000337', 'BGC0000341', 'BGC0000342', 'BGC0000345', 'BGC0000348', 'BGC0000350', 'BGC0000352', 'BGC0000353', 'BGC0000355', 'BGC0000360', 'BGC0000362', 'BGC0000371', 'BGC0000372', 'BGC0000373', 'BGC0000376', 'BGC0000378', 'BGC0000380', 'BGC0000394', 'BGC0000395', 'BGC0000396', 'BGC0000404', 'BGC0000406', 'BGC0000408', 'BGC0000409', 'BGC0000410', 'BGC0000414', 'BGC0000416', 'BGC0000418', 'BGC0000419', 'BGC0000421', 'BGC0000422', 'BGC0000428', 'BGC0000430', 'BGC0000438', 'BGC0000439', 'BGC0000443', 'BGC0000444', 'BGC0000446', 'BGC0000448', 'BGC0000454', 'BGC0000458', 'BGC0000465', 'BGC0000816', 'BGC0000818', 'BGC0000819', 'BGC0000949', 'BGC0000954', 'BGC0000958', 'BGC0000959', 'BGC0000960', 'BGC0000961', 'BGC0000962', 'BGC0000964', 'BGC0000966', 'BGC0000967', 'BGC0000968', 'BGC0000970', 'BGC0000973', 'BGC0000974', 'BGC0000975', 'BGC0000976', 'BGC0000977', 'BGC0000982', 'BGC0000983', 'BGC0000988', 'BGC0000989', 'BGC0000990', 'BGC0000991', 'BGC0000993', 'BGC0000994', 'BGC0000995', 'BGC0000996', 'BGC0000997', 'BGC0000999', 'BGC0001001', 'BGC0001010', 'BGC0001011', 'BGC0001012', 'BGC0001013', 'BGC0001015', 'BGC0001016', 'BGC0001017', 'BGC0001020', 'BGC0001021', 'BGC0001022', 'BGC0001023', 'BGC0001024', 'BGC0001025', 'BGC0001026', 'BGC0001035', 'BGC0001036', 'BGC0001037', 'BGC0001039', 'BGC0001040', 'BGC0001042', 'BGC0001043', 'BGC0001044', 'BGC0001048', 'BGC0001049', 'BGC0001051', 'BGC0001052', 'BGC0001053', 'BGC0001054', 'BGC0001055', 'BGC0001057', 'BGC0001059', 'BGC0001069', 'BGC0001070', 'BGC0001071', 'BGC0001084', 'BGC0001088', 'BGC0001098', 'BGC0001099', 'BGC0001100', 'BGC0001101', 'BGC0001103', 'BGC0001106', 'BGC0001108', 'BGC0001109', 'BGC0001110', 'BGC0001111', 'BGC0001114', 'BGC0001115', 'BGC0001116', 'BGC0001122', 'BGC0001123', 'BGC0001125', 'BGC0001136', 'BGC0001142', 'BGC0001147', 'BGC0001149', 'BGC0001153', 'BGC0001165', 'BGC0001178', 'BGC0001182', 'BGC0001187', 'BGC0001192', 'BGC0001193', 'BGC0001216', 'BGC0001218', 'BGC0001220', 'BGC0001249', 'BGC0001250', 'BGC0001255', 'BGC0001268', 'BGC0001290', 'BGC0001327', 'BGC0001328', 'BGC0001329', 'BGC0001331', 'BGC0001332', 'BGC0001342', 'BGC0001344', 'BGC0001351', 'BGC0001367', 'BGC0001368', 'BGC0001370', 'BGC0001371', 'BGC0001394', 'BGC0001425', 'BGC0001426', 'BGC0001427', 'BGC0001428', 'BGC0001429', 'BGC0001430', 'BGC0001431', 'BGC0001434', 'BGC0001435', 'BGC0001437', 'BGC0001442', 'BGC0001445', 'BGC0001449', 'BGC0001455', 'BGC0001457', 'BGC0001459', 'BGC0001460', 'BGC0001467', 'BGC0001477', 'BGC0001502', 'BGC0001519', 'BGC0001545', 'BGC0001553', 'BGC0001560', 'BGC0001563', 'BGC0001609', 'BGC0001635', 'BGC0001652', 'BGC0001667', 'BGC0001671', 'BGC0001672', 'BGC0001692', 'BGC0001702', 'BGC0001705', 'BGC0001711', 'BGC0001717', 'BGC0001718', 'BGC0001737', 'BGC0001738', 'BGC0001740', 'BGC0001741', 'BGC0001750', 'BGC0001751', 'BGC0001768', 'BGC0001780', 'BGC0001800', 'BGC0001801', 'BGC0001816', 'BGC0001820', 'BGC0001845', 'BGC0001857', 'BGC0001873', 'BGC0001881', 'BGC0001902', 'BGC0001919', 'BGC0001942', 'BGC0001970', 'BGC0001971', 'BGC0002014', 'BGC0002035', 'BGC0002037', 'BGC0002046', 'BGC0002050', 'BGC0002055', 'BGC0002058', 'BGC0002059', 'BGC0002066', 'BGC0002067', 'BGC0002068', 'BGC0002075', 'BGC0002082', 'BGC0002091', 'BGC0002093', 'BGC0002095', 'BGC0002114', 'BGC0002118', 'BGC0002119', 'BGC0002131', 'BGC0002135', 'BGC0002157', 'BGC0002158', 'BGC0002160', 'BGC0002164', 'BGC0002166', 'BGC0002169', 'BGC0002197', 'BGC0002198', 'BGC0002204', 'BGC0002205', 'BGC0002210', 'BGC0002230', 'BGC0002239', 'BGC0002242', 'BGC0002248', 'BGC0002250', 'BGC0002256', 'BGC0002259', 'BGC0002264', 'BGC0002268', 'BGC0002272', 'BGC0002273', 'BGC0002287', 'BGC0002296', 'BGC0002300', 'BGC0002304', 'BGC0002314', 'BGC0002340', 'BGC0002344', 'BGC0002367', 'BGC0002373', 'BGC0002378', 'BGC0002408', 'BGC0002412', 'BGC0002434', 'BGC0002437', 'BGC0002456', 'BGC0002468', 'BGC0002491', 'BGC0002492', 'BGC0002495', 'BGC0002505', 'BGC0002509', 'BGC0002542', 'BGC0002543', 'BGC0002546', 'BGC0002555', 'BGC0002565', 'BGC0002611', 'BGC0002613', 'BGC0002649', 'BGC0002735'}

    fdr_data = {}
    tops = {}

    sample_size = len(set(mibig4_results['AccGen'].unique()) - exclude)

    nerpa2_results_samples = [sample_nerpa_results(mibig3_results, sample_size, exclude=mibig3_ids)
                              for _ in range(num_samples)]
    fdr_data_samples = [calculate_fdr_data(
        nerpa2_results_sample,
        args.max_rank
    ) for nerpa2_results_sample in nerpa2_results_samples]

    fdr_data_samples_values = [fdr[1] for fdr, _ in fdr_data_samples]
    fdr_data_mibig4, _ = calculate_fdr_data(
        mibig4_results,
        args.max_rank, exclude
    )
    fdr_data_mibig4 = fdr_data_mibig4[1]
    stacked = np.vstack(fdr_data_samples_values)

    # Compute mean and standard deviation for each position
    mean = np.mean(stacked, axis=0)
    std = np.std(stacked, axis=0)

    # Calculate the upper and lower bounds
    upper_bound = mean + std
    lower_bound =  mean - std
    lower_bound = np.where(lower_bound < 0, 0, lower_bound)

    # Plotting
    # q: set y range to 0-1
    plt.ylim(0, 1)
    plt.locator_params(axis='y', nbins=10)
    plt.locator_params(axis='x', nbins=14)
    plt.grid()
    plt.plot(mean, label='MIBiG 3: Mean', color='blue', linestyle='--')
    plt.plot(fdr_data_mibig4, label='FDR for MIBiG 4', color='red')
    plt.fill_between(range(len(mean)), lower_bound, upper_bound, color='lightblue', alpha=0.5, label='FDR for MIBiG 3: Mean ± Std Dev')
    plt.xlabel('Number of BGCs')
    plt.ylabel('FDR')
    #plt.title('Mean ± Standard Deviation')
    plt.legend()
    plt.savefig('fdrplot_nerpa2_mibig3_vs_nerpa2_mibig4.png')
    plt.show()
    # q: save the plot


def main():
    args = parse_args()

    # mode = 'compare_mibig3_and_mibig4'
    mode = 'default'
    if mode == 'compare_mibig3_and_mibig4':
        run_compare_mibig3_and_mibig4(args)

    # mode == 'default'
    run_default(args)


if __name__ == '__main__':
    main()