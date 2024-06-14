from typing import Iterable, List, Callable, TypeVar
from itertools import chain, combinations

T = TypeVar('T')

def hamming_distance(s1: str, s2: str) -> int:
    if len(s1) != len(s2):
        raise ValueError("Strings must be of equal length")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def list_monad_compose(*functions: Callable) -> Callable[[List], List]:
    def composed(lst: List) -> List:
        for f in functions:
            lst = list(chain(*(map(f, lst))))
        return lst
    return composed


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
    max_blocks = n + 2 if max_blocks is not None else max_blocks + 2  # +2 for imaginary blocks at the ends
    for blocks_num in range(4, max_blocks + 1):  # numbers of contiguous blocks. 2, 3 are impossible to make non consecutive (given that blocks at the ends are stationary)
        permutations = list(permutations_no_consecutive(blocks_num - 2))
        for block_starts in combinations(range(n + 1), blocks_num - 1):  # -1 and n are imaginary blocks at the ends
            block_starts = [-1] + list(block_starts) + [n + 1]
            blocks = [list(range(block_starts[i], block_starts[i + 1])) for i in range(blocks_num)]
            for perm in permutations:
                joined_blocks = blocks[0] + list(chain(*(blocks[perm[i] + 1] for i in range(blocks_num - 2)))) + blocks[-1]
                yield joined_blocks[1:-1]  # remove -1 and n


def generate_permutations(lst: List[T]) -> Iterable[List[T]]:
    for idx_perm in generate_permutations_idxs(len(lst)):
        yield [lst[i] for i in idx_perm]