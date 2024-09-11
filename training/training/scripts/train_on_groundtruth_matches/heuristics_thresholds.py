from typing import Dict, List, Tuple, Optional, Set
from pathlib import Path
from dataclasses import dataclass
from src.data_types import MonomerResidue, LogProb
from src.generic.combinatorics import intersection_with_repeats
from src.generic.statistics import fit_between_lines
from src.matching.heuristic_matching import HeuristicMatchingConfig
from itertools import chain, pairwise, product
from collections import defaultdict


def nrp_residues(match: dict,
                 MISSING_VALUE: str = '---') -> List[MonomerResidue]:
    return [alignment_step['NRP_residue']
            for alignment in match['Alignments']
            for alignment_step in alignment
            if alignment_step['NRP_residue'] != MISSING_VALUE]


def top_n_predictions(residue_scores: Dict[MonomerResidue, LogProb],
                      num_top_predictions: int) -> List[MonomerResidue]:
    return sorted(residue_scores, key=residue_scores.get, reverse=True)[:num_top_predictions]


def bgc_residues(bgc_variant: dict,
                 num_top_predictions: int = 3) -> List[MonomerResidue]:
    return list(chain(*[top_n_predictions(bgc_module['residue_score'], num_top_predictions)
                        for bgc_fragment in bgc_variant['fragments']
                        for bgc_module in bgc_fragment]))


def num_common_residues(match: dict,
                        bgc_variant: dict,
                        num_top_predictions: int = 3,
                        MISSING_VALUE: str = '---') -> int:
    return len(intersection_with_repeats(nrp_residues(match, MISSING_VALUE),
                                         bgc_residues(bgc_variant, num_top_predictions)))


def get_common_cons_pairs_cnt(match: dict, bgc_variant: dict) -> int:
    nrp_pairs = pairwise(nrp_residues(match))

    bgc_set_pairs = pairwise(top_n_predictions(bgc_module['residue_score'], 3)
                             for bgc_fragment in bgc_variant['fragments']
                             for bgc_module in bgc_fragment)
    bgc_pairs = chain(*(product(s1, s2) for s1, s2 in bgc_set_pairs))
    return len(intersection_with_repeats(nrp_pairs, bgc_pairs))


# TODO: refactor
def plot_box_plots(crc_for_nrp_len: Dict[int, List[int]],
                   nrp_bgc_lengths: Dict[int, List[int]],
                   cons_pairs_for_nrp_len: Dict[int, List[int]],
                   params: HeuristicMatchingConfig,
                   output_dir: Path):
    from matplotlib import pyplot as plt
    plt.boxplot([crc_for_nrp_len[bgc_len] for bgc_len in sorted(crc_for_nrp_len.keys())],
                positions=sorted(crc_for_nrp_len.keys()))
    # q: plot the line from the params
    k, b, d = params.LINEAR_DISCARD_PARAMS_AA_CONTENTS
    plt.plot(sorted(crc_for_nrp_len.keys()),
             [k * x + b - d for x in sorted(crc_for_nrp_len.keys())],
             color='red')
    plt.plot(sorted(crc_for_nrp_len.keys()),
             [k * x + b + d for x in sorted(crc_for_nrp_len.keys())],
             color='red')
    plt.xticks(range(1, len(crc_for_nrp_len) + 1), sorted(crc_for_nrp_len.keys()))
    plt.xlabel('NRP length')
    plt.ylabel('Common residues count')
    plt.savefig(output_dir / 'common_residues_counts_boxplot.png')
    plt.close()

    # q: make box and whisker plot for bgc and nrp length
    plt.boxplot([nrp_bgc_lengths[nrp_len] for nrp_len in sorted(nrp_bgc_lengths.keys())],
                positions=sorted(nrp_bgc_lengths.keys()))
    k, b, d = params.LINEAR_DISCARD_PARAMS_LENGTHS
    plt.plot(sorted(nrp_bgc_lengths.keys()),
             [k * x + b - d for x in sorted(nrp_bgc_lengths.keys())],
             color='red')
    plt.plot(sorted(nrp_bgc_lengths.keys()),
             [k * x + b + d for x in sorted(nrp_bgc_lengths.keys())],
             color='red')
    plt.xlabel('NRP length')
    plt.ylabel('BGC length')
    plt.savefig(output_dir / 'bgc_nrp_lengths_boxplot.png')
    plt.close()

    plt.boxplot([cons_pairs_for_nrp_len[nrp_len] for nrp_len in sorted(cons_pairs_for_nrp_len.keys())],
                positions=sorted(cons_pairs_for_nrp_len.keys()))
    k, b, d = params.LINEAR_DISCARD_PARAMS_CONS_PAIRS
    plt.plot(sorted(nrp_bgc_lengths.keys()),
             [k * x + b - d for x in sorted(cons_pairs_for_nrp_len.keys())],
             color='red')
    plt.plot(sorted(cons_pairs_for_nrp_len.keys()),
             [k * x + b + d for x in sorted(cons_pairs_for_nrp_len.keys())],
             color='red')
    plt.xlabel('NRP length')
    plt.ylabel('# common consequtive residue pairs')
    plt.savefig(output_dir / 'common_consequtive_residues_boxplot.png')
    plt.close()


def calculate_heuristic_parameters(matches_with_bgc_variants: List[Tuple[dict, dict]],
                                   output_dir: Path) -> HeuristicMatchingConfig:
    def bgc_len(bgc_variant: dict) -> int:
        return sum(len(fragment) for fragment in bgc_variant['fragments'])

    def nrp_len(match: dict) -> int:
        return len(nrp_residues(match))

    crc_for_nrp_len = defaultdict(list)
    for match, bgc_variant in matches_with_bgc_variants:
        crc = num_common_residues(match, bgc_variant)
        crc_for_nrp_len[nrp_len(match)].append(crc)

    nrp_bgc_lengths = defaultdict(list)
    for match, bgc_variant in matches_with_bgc_variants:
        nrp_bgc_lengths[nrp_len(match)].append(bgc_len(bgc_variant))

    common_cons_pairs_cnt_for_nrp_len = defaultdict(list)
    for match, bgc_variant in matches_with_bgc_variants:
        common_cons_pairs_cnt_for_nrp_len[nrp_len(match)].append(
            get_common_cons_pairs_cnt(match, bgc_variant)
        )
    params = HeuristicMatchingConfig(
        LINEAR_DISCARD_PARAMS_LENGTHS=fit_between_lines(nrp_bgc_lengths, 0.05),  # TODO: put in config
        LINEAR_DISCARD_PARAMS_AA_CONTENTS=fit_between_lines(crc_for_nrp_len, 0.05),
        LINEAR_DISCARD_PARAMS_CONS_PAIRS=fit_between_lines(common_cons_pairs_cnt_for_nrp_len, 0.05),
        NUM_TOP_PREDICTIONS=3  # TODO: magic number, bad
    )
    plot_box_plots(crc_for_nrp_len, nrp_bgc_lengths, common_cons_pairs_cnt_for_nrp_len, params, output_dir)

    return params

