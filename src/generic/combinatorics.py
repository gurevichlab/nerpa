from typing import Generator, Iterable, List, Tuple, TypeVar
from itertools import combinations, chain, combinations_with_replacement
from collections import Counter


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
    return (list(xs[i] for i in perm)
            for perm in generate_permutations_idxs(len(xs), max_blocks))


def cyclic_shifts(xs: List[T]) -> Iterable[List[T]]:
    for i in range(len(xs)):
        yield xs[i:] + xs[:i]


def split_into_k_blocks(seq: List[T], k: int) -> Generator[List[List[T]], None, None]:
    """Generate all ways to split a sequence into exactly k non-empty blocks."""
    if k == 1:
        yield [seq]  # Only one way to split the sequence into 1 block (the sequence itself)
        return

    # Iterate through all possible first blocks
    for i in range(1, len(seq) - k + 2):  # Leave room for at least k-1 more blocks
        first_part = seq[:i]
        # Recursively split the remaining sequence into k-1 blocks
        for rest in split_into_k_blocks(seq[i:], k - 1):
            yield [first_part] + rest

def split_sequence_blocks(seq: List[T]) -> Generator[List[List[T]], None, None]:
    """Generate all ways to split a sequence into increasing number of blocks."""
    for k in range(1, len(seq) + 1):  # Iterate through all possible number of blocks
        yield from split_into_k_blocks(seq, k)


def all_subsets(xs: List[T]) -> Generator[List[T], None, None]:
    for k in range(len(xs) + 1):
        yield from combinations(xs, k)


def split_sequence_subseqs_idxs(seq: List[int]) -> Iterable[List[List[int]]]:
    if not seq:
        yield []
        return
    for fst_subset in all_subsets(seq[1:]):  # first subset contains the first element
        fst_subset = [seq[0]] + list(fst_subset)
        rest = [x for x in seq if x not in fst_subset]
        for next_subsets in split_sequence_subseqs_idxs(rest):
            yield [fst_subset] + next_subsets


def split_sequence_subseqs(seq: List[T]) -> Iterable[List[List[T]]]:
    for subseq_idxs_split in split_sequence_subseqs_idxs(list(range(len(seq)))):
        yield [[seq[i] for i in subseq_idxs]
               for subseq_idxs in subseq_idxs_split]


def intersection_with_repeats(xs: Iterable[T], ys: Iterable[T]) -> List[T]:
    '''
    returns a list of elements that are present in both lists (with repeats)
    '''
    return list((Counter(xs) & Counter(ys)).elements())


def is_subsequence(subseq: List[T], seq: List[T]) -> bool:
    '''
    checks if subseq is a subsequence of seq
    '''
    subseq_idx = 0
    seq_idx = 0
    while subseq_idx < len(subseq) and seq_idx < len(seq):
        if subseq[subseq_idx] == seq[seq_idx]:
            subseq_idx += 1
        seq_idx += 1
    return subseq_idx == len(subseq)


def remove_runs_of_equal_elements(xs: Iterable[T]) -> Iterable[T]:
    xs_it = iter(xs)
    try:
        last = next(xs_it)  # Initialize result with the first element of the input list
        yield last

    except StopIteration:
        return []
    while True:
        try:
            while (x := next(xs_it)) == last:
                pass
            yield (last := x)
        except StopIteration:
            break


def filter_unique(xs: Iterable[T]) -> Iterable[T]:
    seen = set()
    xs_it = iter(xs)
    while True:
        try:
            x = next(xs_it)
        except StopIteration:
            break
        if x not in seen:
            seen.add(x)
            yield x