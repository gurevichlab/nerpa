from typing import Callable, Dict, List, Optional, Tuple
from math import log
import matplotlib.pyplot as plt
from pathlib import Path
import math
from src.data_types import LogProb

def step_cost(step_y_norm: float, bin: Tuple[int, int]) -> float:
    '''
    Each bin is treated as a realization of a binomial distribution with bin[0] false and bin[1] true outcomes.
    The cost of a step is the log-likelihood of the binomial distribution with the probability of success equal to
    step_y_norm.
    '''
    num_false, num_true = bin
    step_y_norm = min(1 - 1e-10, max(1e-10, step_y_norm))
    return num_true * log(step_y_norm) + num_false * log(1 - step_y_norm)


def fit_step_function_to_bins(bins: List[Tuple[int, int]],  # TODO: pseudocounts?
                              step_range: int) -> List[float]:
    '''
    Each bin is treated as a realization of a binomial distribution with bin[0] false and bin[1] true outcomes.
    The function finds the non-decreasing step function that maximizes the total log-likelihood of the binomial distributions.
    step_range determines the range of possible step values: each step value is of the form i / precision
    The function returns the list of step values.
    '''

    minus_infinity = float('-inf')
    # dp[(bin_idx, step_y)] is the maximum log-likelihood of the bins[:bin_idx + 1]
    # with the last step equal to step_y / step_range
    dp: Dict[Tuple[int, int], float] = {}
    # parent[(bin_idx, step_y)] is the last step value that maximizes dp[(bin_idx, step_y)]
    # used for backtracking
    parent: Dict[Tuple[int, int], int] = {}
    for step_y in range(step_range + 1):  # initialize the first column
        dp[(-1, step_y)] = 0

    #  calculate dp
    for bin_idx, bin in enumerate(bins):
        if bin[0] + bin[1] == 0:  # for empty bins just copy the previous column
            for step_y in range(step_range + 1):
                dp[(bin_idx, step_y)] = dp[(bin_idx - 1, step_y)]
                parent[(bin_idx, step_y)] = step_y
            continue

        for step_y in range(step_range + 1):
            dp[(bin_idx, step_y)] = minus_infinity
            cost = step_cost(step_y / step_range, bin)
            for prev_step_y in range(step_y + 1):
                if dp[(bin_idx - 1, prev_step_y)] + cost > dp[(bin_idx, step_y)]:
                    dp[(bin_idx, step_y)] = dp[(bin_idx - 1, prev_step_y)] + cost
                    parent[(bin_idx, step_y)] = prev_step_y

    # backtracking
    last_step_y = max(range(step_range + 1), key=lambda step_y: dp[(len(bins) - 1, step_y)])
    steps = []
    for bin_idx in range(len(bins) - 1, -1, -1):
        steps.append(last_step_y)
        last_step_y = parent[(bin_idx, last_step_y)]

    return [step_y / step_range for step_y in steps[::-1]]


def create_bins(score_correctness_log: List[Tuple[LogProb, bool]],
                num_bins: int) -> List[
    Tuple[int, int]]:  # (num_false, num_true) for each bin [i/numbins, (i+1)/numbins)
    score_correctness = [(math.e ** score, correctness) for score, correctness in score_correctness_log]
    bins: List[List[bool]] = [[] for _ in range(num_bins)]
    bin_size = 1 / num_bins
    for score, correctness in score_correctness:
        bin_idx = min(int(score / bin_size), num_bins - 1)
        bins[bin_idx].append(correctness)
    return [(bin.count(False), bin.count(True)) for bin in bins]


def plot_step_function(score_correctness_bins: List[Tuple[int, int]],
                       # (num_false, num_true) for each bin [i/numbins, (i+1)/numbins)
                       step_function: List[float],
                       out_file: Optional[Path] = None,
                       show: bool = False) -> None:
    assert len(score_correctness_bins) == len(step_function)
    bin_size = 1 / len(score_correctness_bins)
    bin_centers = [i * bin_size + bin_size / 2 for i in range(len(score_correctness_bins))]
    bin_fractions = [0 if not num_true else num_true / (num_true + num_false)
                     for num_false, num_true in score_correctness_bins]
    plt.bar(bin_centers, bin_fractions, align='center', width=bin_size)
    plt.plot(bin_centers, step_function, color='red')
    if show:
        plt.show()
    if out_file is not None:
        plt.savefig(out_file)
    plt.close()
    return out_file


def plot_step_function_stacked(score_correctness_bins: List[Tuple[int, int]],
                               step_function: List[float],
                               out_file: Optional[Path] = None,
                               show: bool = False) -> Optional[Path]:
    assert len(score_correctness_bins) == len(step_function)

    bin_size = 1 / len(score_correctness_bins)
    bin_centers = [i * bin_size + bin_size / 2 for i in range(len(score_correctness_bins))]

    num_false = [num_false for num_false, num_true in score_correctness_bins]
    num_true = [num_true for num_false, num_true in score_correctness_bins]

    plt.bar(bin_centers, num_false, align='center', width=bin_size, color='orange', label='Incorrect')
    plt.bar(bin_centers, num_true, align='center', width=bin_size, bottom=num_false, color='blue', label='Correct')

    plt.plot(bin_centers, step_function, color='red', label='Step Function')

    plt.legend()

    if show:
        plt.show()
    if out_file is not None:
        plt.savefig(out_file)

    plt.close()
    return out_file


def fit_step_function(score_correctness: List[Tuple[LogProb, bool]],
                      num_bins: int,
                      step_range: int) -> List[float]:
    score_correctness_bins = create_bins(score_correctness, num_bins)
    step_function = fit_step_function_to_bins(score_correctness_bins, step_range)
    print('Plotting...')
    plot_step_function(score_correctness_bins, step_function)
    return step_function


def create_step_function(steps: List[float], step_len: float) -> Callable[[float], float]:
    def step_function(x: float) -> float:
        for i, step in enumerate(steps):
            if x < (i + 1) * step_len:
                return steps[i]
        return steps[-1]
    return step_function
