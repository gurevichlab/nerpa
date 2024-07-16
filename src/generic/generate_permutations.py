from typing import Iterable, List, TypeVar
from itertools import combinations, chain


T = TypeVar('T')


def permutations_no_consecutive(n: int) -> Iterable[List[int]]:
    '''
    generates all permutations of {0,...,n-1} with no consecutive numbers
    also, 0 and n-1 should change their positions (they are consequent with imaginary -1 and n)
    '''
    def extend_permutation(perm: List[int]) -> Iterable[List[int]]:
        if len(perm) == n:
            if perm[-1] != n - 1:
                yield perm
            return
        for num in range(n):
            if num not in perm and num != (perm[-1] + 1 if perm else 0):
                yield from extend_permutation(perm + [num])

    yield from extend_permutation([])


def generate_permutations_idxs(n: int, max_blocks: int = None) -> Iterable[List[int]]:
    '''
    generates all permutations of {0,...,n-1} in the order of increasing number of contiguous blocks
    we also add imaginary -1 and n at the ends that never move under a permutation
    '''
    yield list(range(0, n))
    if max_blocks is None:
        max_blocks = n + 2
    for blocks_num in range(4, max_blocks + 1):  # numbers of contiguous blocks. 2, 3 are impossible to make non consecutive (given that blocks at the ends are stationary)
        permutations = list(permutations_no_consecutive(blocks_num - 2))
        for block_starts in combinations(range(n + 1), blocks_num - 1):  # -1 and n are imaginary blocks at the ends
            block_starts = [-1] + list(block_starts) + [n + 1]
            blocks = [list(range(block_starts[i], block_starts[i + 1])) for i in range(blocks_num)]
            for perm in permutations:
                joined_blocks = blocks[0] + list(chain(*(blocks[perm[i] + 1] for i in range(blocks_num - 2)))) + blocks[-1]
                yield joined_blocks[1:-1]  # remove -1 and n


def generate_permutations(xs_: Iterable[T], max_blocks: int = None) -> Iterable[List[T]]:
    xs = list(xs_)
    print(xs_)
    print(list(generate_permutations_idxs(len(xs), max_blocks)))
    return (list(xs[i] for i in perm)
            for perm in generate_permutations_idxs(len(xs), max_blocks))