from collections.abc import Hashable
from typing import (
    Callable,
    Generator,
    Iterable,
    List,
    Tuple,
    TypeVar, Optional
)
from itertools import (
    combinations,
    chain,
    combinations_with_replacement,
    groupby, product, pairwise
)
from collections import Counter

from more_itertools import peekable

T = TypeVar('T')
U = TypeVar('U')


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


def all_subsets(xs: List[T],
                subset_len: Optional[int] = None) -> Generator[Tuple[T, ...], None, None]:
    if subset_len is None:
        subset_len_choices = range(len(xs) + 1)
    else:
        subset_len_choices = [subset_len]
    for k in subset_len_choices:
        yield from combinations(xs, k)


def split_sequence_subseqs_idxs(seq: List[int],
                                num_subseqs: Optional[int] = None) -> Iterable[List[List[int]]]:
    """
    Generates all possible ways to partition the sequence into ordered non-empty subsequences.
    Each partition is a list of lists, where each sublist is a subsequence (subset in order) of the original sequence.
    """
    if not seq:
        yield []
        return

    if num_subseqs == 1:
        yield [seq]
        return

    if num_subseqs is None:
        fst_subset_choices = ((seq[0],) + subset
                              for subset in all_subsets(seq[1:]))
        remaining_subseqs = None
    else:
        fst_subset_choices = ((seq[0],) + subset
                              for subset_len in range(len(seq) - num_subseqs + 1)
                              for subset in all_subsets(seq[1:], subset_len))
        remaining_subseqs = num_subseqs - 1

    for fst_subset in fst_subset_choices:  # first subset contains the first element
        rest = [x for x in seq if x not in fst_subset]
        for next_subsets in split_sequence_subseqs_idxs(rest, remaining_subseqs):
            yield [fst_subset] + next_subsets


def split_sequence_subseqs(seq: List[T],
                           num_subseqs: Optional[int] = None) -> Iterable[List[List[T]]]:
    """
    Generates all possible ways to partition the sequence into ordered non-empty subsequences.
    Each partition is a list of lists, where each sublist is a subsequence of the original sequence.
    """
    assert num_subseqs is None or 0 < num_subseqs <= len(seq), \
        "num_subseqs must be None or between 1 and the length of the sequence inclusive."

    for subseq_idxs_split in split_sequence_subseqs_idxs(list(range(len(seq))), num_subseqs):
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


def subseq_occurences(subseq: List[T], seq: List[T]) -> Iterable[List[int]]:
    '''
    returns starting indices of all occurrences of subseq in seq
    '''
    if not is_subsequence(subseq, seq):
        return

    item_occurences = [
        [i for i, x in enumerate(seq) if x == item]
        for item in subseq
    ]

    for occurence_idxs in product(*item_occurences):
        if all(idx1 < idx2
               for idx1, idx2 in pairwise(occurence_idxs)):
            yield list(occurence_idxs)


def remove_runs_of_equal_elements(xs: Iterable[T]) -> Iterable[T]:
    _SENTINEL = object()  # unique sentinel value to distinguish from possible None elements
    xs_it = iter(xs)
    if (last := next(xs_it, _SENTINEL)) != _SENTINEL:
        yield last
    while True:
        while (x := next(xs_it, _SENTINEL)) == last:
            pass
        if x == _SENTINEL:
            break
        yield (last := x)


def filter_unique(xs: Iterable[T],
                  key: Callable[[T], Hashable] = lambda x: x) -> Iterable[T]:
    _SENTINEL = object()
    seen_keys = set()
    xs_it = iter(xs)
    while True:
        if (x := next(xs_it, _SENTINEL)) == _SENTINEL:
            break
        x_key = key(x)
        if x_key not in seen_keys:
            seen_keys.add(x_key)
            yield x


def sort_groupby(items: Iterable[T],
                 key: Callable[[T], U],
                 reverse: bool=False) -> Iterable[Tuple[U, Iterable[T]]]:
    return groupby(sorted(items, key=key, reverse=reverse), key=key)


def split_into_iterations(items: Iterable[T],
                         key: Callable[[T], Hashable] = lambda x: x) \
        -> Generator[List[T], None, None]:
    """
    Extracts "iterations" from an iterable.

    An iteration is a maximal contiguous list, such that
    if one collapses all key repeat runs in it, the remaining keys will be unique.

    Parameters:
    - items: iterable of elements to split (order is preserved).
    - key: function mapping each element to a hashable key (default: identity).

    Yields:
    Each iteration as a List[T]. For an empty input yields nothing.

    Example:
    - For input keys [A, A, B, B, A, C] the result will be [[A, A, B, B], [A, C]]
      because the second occurrence of A starts a new iteration after B.
    """
    _SENTINEL = object()

    def extract_iteration(items_iter: 'peekable[T]') -> Optional[List[T]]:
        """Consume and return a single iteration from the peekable iterator or None."""
        cur = next(items_iter, _SENTINEL)
        if cur == _SENTINEL:
            return None
        cur_key = key(cur)

        seen_keys = {cur_key}
        extracted = [cur]
        while True:
            nxt = items_iter.peek(_SENTINEL)
            if nxt == _SENTINEL:
                break
            nxt_key = key(nxt)
            if nxt_key in seen_keys and nxt_key != cur_key:
                break
            extracted.append(nxt)
            seen_keys.add(nxt_key)
            cur, cur_key = nxt, nxt_key
            next(items_iter)

        return extracted

    items_it = peekable(items)
    while iteration := extract_iteration(items_it):
        yield iteration
