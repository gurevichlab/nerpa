import argparse
from pathlib import Path
from typing import List, Tuple, Union, Set, Iterable
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
    output_group.add_argument("--fig", type=str, default="./fdrplot.png",
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
    report['AccGen'] = report['AccGen'].map(lambda x: x.split('.')[0])
    report['AccStruct'] = report['NRP_ID'].map(lambda x: x.split('.')[0])
    return report


def calculate_fdr_data(report: pd.DataFrame,
                       rank: int,
                       exclude: Set[str] = []) -> Tuple[Iterable, pd.DataFrame]:
    report['Correct'] = report['AccStruct'] == report['AccGen']
    report['Rank'] = report.groupby('AccGen')['Score'].rank('average', ascending=False)

    mask = report['Rank'] < rank + 1
    mask &= ~report['AccGen'].isin(exclude)
    top = report.loc[mask].sort_values(['Correct', 'Score'], ascending=[False, False]) \
        .groupby('AccGen').head(1).sort_values('Score', ascending=False)

    fdr = np.zeros((3, top.shape[0]))
    for i in range(top.shape[0]):
        n = i + 1
        _df = top.head(n)
        fdr[0, i] = n
        n_correct = _df['Correct'].sum()
        fdr[1, i] = 1 - (n_correct / n)
        fdr[2, i] = _df['Score'].min()

    return fdr, top


def main():
    args = parse_args()
    exclude = set(x.strip() for x in args.exclude.split(','))

    fdr_data = {}
    tops = {}

    if args.nerpa1:
        k = 'Nerpa 1.1'
        fdr_data[k], tops[k] = calculate_fdr_data(
            parse_nerpa1_results(Path(args.nerpa1)),
            args.max_rank, exclude
        )
    if args.nerpa2:
        for i, p in enumerate(args.nerpa2):
            k = f'Nerpa 2 run {i}'
            fdr_data[k], tops[k] = calculate_fdr_data(
                parse_nerpa2_results(Path(p)),
                args.max_rank, exclude
            )

    f, ax = plt.subplots(figsize=(8, 6))
    for k, v in fdr_data.items():
        _ = ax.plot(v[0], v[1], label=k)

    _ = ax.legend()
    _ = ax.set_ylim((0, 1))
    _ = ax.locator_params(axis='y', nbins=10)
    _ = ax.locator_params(axis='x', nbins=14)
    _ = ax.set_xlabel('Num BGCs')
    _ = ax.set_ylabel('FDR')
    _ = ax.grid()

    f.savefig(args.fig)
    for k, df in tops.items():
        df.loc[~df['Correct']].to_csv(f'{args.tsv_pref}_{"_".join(k.split())}.tsv', index=False, sep='\t')


if __name__ == '__main__':
    main()