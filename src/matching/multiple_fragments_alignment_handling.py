from typing import (
    Dict,
    Iterable,
    List,
    Optional,
    Tuple
)
from src.data_types import BGC_Fragment, NRP_Monomer, NRP_Fragment, BGC_Module
from src.generic.combinatorics import cyclic_shifts, split_sequence_subseqs
from src.matching.alignment_types import Alignment, combined_alignments_score
from src.matching.scoring_helper import ScoringHelper
from src.matching.dp import LogProb
from src.matching.dp import get_alignment
from src.matching.alignment_types import alignment_score, AlignmentStep, AlignmentStepType
from itertools import chain, combinations, product, permutations


def fragment_monomer_sequences(fragment: NRP_Fragment) -> Iterable[List[NRP_Monomer]]:
    if fragment.is_cyclic:
        return cyclic_shifts(fragment.monomers)
    else:
        return [fragment.monomers]


def fragments_joined_monomer_sequences(fragments: List[NRP_Fragment]) -> Iterable[List[NRP_Monomer]]:
    return (list(chain(*monomer_sequences))
            for monomer_sequences in product(*map(fragment_monomer_sequences, fragments)))


def retrieve_alignments(parent: Dict[Tuple[int, int], Tuple[Tuple[int, int], Optional[Alignment]]],
                        bgc_fragments_len: int,
                        nrp_fragments_len: int) -> List[Alignment]:
    bgc_idx, nrp_idx = bgc_fragments_len, nrp_fragments_len
    alignments: List[Alignment] = []
    while (bgc_idx, nrp_idx) in parent:
        (bgc_prev, nrp_prev), alignment = parent[(bgc_idx, nrp_idx)]
        if alignment is not None:
            alignments.append(alignment)
        bgc_idx, nrp_idx = bgc_prev, nrp_prev
    return alignments[::-1]


def skipping_alignment(modules: List[BGC_Module],
                       dp_helper: ScoringHelper) -> Alignment:
    return ([AlignmentStep(bgc_module=modules[0],
                          nrp_monomer=None,
                          score=dp_helper.scoring_config.bgc_fragment_skip_penalty_at_end,
                          action=AlignmentStepType.BGC_MODULE_SKIP)]
    + [AlignmentStep(bgc_module=module,
                     nrp_monomer=None,
                     score=dp_helper.scoring_config.bgc_fragment_skip_penalty_in_middle,
                     action=AlignmentStepType.BGC_MODULE_SKIP)
            for module in modules[1:]])


def get_alignments_against_monomers(bgc_fragments: List[BGC_Fragment],
                                    nrp_monomers: List[NRP_Monomer],
                                    dp_helper: ScoringHelper) -> List[Alignment]:
    def alignments_for_split(fst_bgc_idx: int, last_bgc_idx: int) -> List[Alignment]:
        aligned_modules = list(chain(*bgc_fragments[fst_bgc_idx: last_bgc_idx]))
        left_skipping_alignments = [skipping_alignment(bgc_fragment, dp_helper)
                                    for bgc_fragment in bgc_fragments[:fst_bgc_idx]]
        right_skipping_alignments = [skipping_alignment(bgc_fragment, dp_helper)
                                        for bgc_fragment in bgc_fragments[last_bgc_idx:]]
        return (left_skipping_alignments
                + [get_alignment(aligned_modules, nrp_monomers, dp_helper)]
                + right_skipping_alignments)

    return max((alignments_for_split(fst_bgc_idx, last_bgc_idx)
                for fst_bgc_idx, last_bgc_idx in combinations(range(len(bgc_fragments) + 1), 2)),
               key=combined_alignments_score)



def get_multiple_fragments_alignment(bgc_fragments: List[BGC_Fragment],
                                     nrp_fragments: List[NRP_Fragment],
                                     dp_helper: ScoringHelper,
                                     permute_fragments: bool = False) -> List[Alignment]:
    if permute_fragments:
        return max((get_multiple_fragments_alignment(bgc_fragments, list(nrp_fragments_permutation), dp_helper, False)
                    for nrp_fragments_permutation in permutations(nrp_fragments)),
                   key=combined_alignments_score)
    # to save time, as these are the most common cases
    if len(nrp_fragments) == 1:
        return max((get_alignments_against_monomers(bgc_fragments, nrp_monomers, dp_helper)
                    for nrp_monomers in fragment_monomer_sequences(nrp_fragments[0])),
                   key=combined_alignments_score)

    if len(bgc_fragments) == 1:
        return [max((get_alignment(bgc_fragments[0], nrp_monomers, dp_helper)
                     for nrp_monomers in fragments_joined_monomer_sequences(nrp_fragments)),
                    key=alignment_score)]

    dp: Dict[Tuple[int, int], LogProb] = {(0, 0): 0.0}  # dp[(i, j)] = score of aligning bgc_fragments[:i] and nrp_fragments[:j]
    parent: Dict[Tuple[int, int], Tuple[Tuple[int, int], Optional[Alignment]]] = {}  # parent[(i, j)] = (i_prev, j_prev, alignment)

    for i, bgc_fragment in enumerate(bgc_fragments, start=1):
        dp[(i, 0)] = dp[(i - 1, 0)] + dp_helper.skip_bgc_fragment_score(bgc_fragment)
        parent[(i, 0)] = ((i - 1, 0), None)

    for j, nrp_fragment in enumerate(nrp_fragments, start=1):
        dp[(0, j)] = dp[(0, j - 1)] + dp_helper.skip_nrp_fragment_score(nrp_fragment)
        parent[(0, j)] = ((0, j - 1), None)

    Transition = Tuple[LogProb, Tuple[int, int], Optional[Alignment]]  # (score, parent_idx, alignment)

    def alignment_transition(bgc_indices: Tuple[int, int], nrp_indices: Tuple[int, int]) -> Transition:
        bgc_start, bgc_end = bgc_indices
        nrp_start, nrp_end = nrp_indices
        bgc_modules = list(chain(*bgc_fragments[bgc_start: bgc_end]))
        nrp_monomers_sequences = fragments_joined_monomer_sequences(nrp_fragments[nrp_start: nrp_end])
        alignment = max((get_alignment(bgc_modules, nrp_monomers, dp_helper)
                         for nrp_monomers in nrp_monomers_sequences),
                        key=alignment_score)

        return (dp[(bgc_start, nrp_start)] + alignment_score(alignment), (bgc_start, nrp_start), alignment)

    for i, bgc_fragment in enumerate(bgc_fragments, start=1):
        for j, nrp_fragment in enumerate(nrp_fragments, start=1):
            transitions: List[Transition] = [
                (dp[(i - 1, j)] + dp_helper.skip_bgc_fragment_score(bgc_fragment), (i - 1, j), None),
                (dp[(i, j - 1)] + dp_helper.skip_nrp_fragment_score(nrp_fragment), (i, j - 1), None)
            ]
            transitions.extend(alignment_transition((i_prev, i), (j_prev, j))
                               for i_prev in range(i)
                               for j_prev in range(j))
            best_transition = max(transitions, key=lambda t: t[0])
            dp[(i, j)], parent[(i, j)] = best_transition[0], best_transition[1:]

    return retrieve_alignments(parent, len(bgc_fragments), len(nrp_fragments))


def get_fragments_alignment(bgc_fragment: BGC_Fragment,
                            nrp_fragment: NRP_Fragment,
                            dp_helper) -> Alignment:
    return max((get_alignment(bgc_fragment, monomers_sequence, dp_helper)
                for monomers_sequence in fragment_monomer_sequences(nrp_fragment)),
               key=alignment_score)


def get_iterative_bgc_alignment(bgc_modules: List[BGC_Module],
                                nrp_fragments: List[NRP_Fragment],
                                dp_helper: ScoringHelper) -> List[Alignment]:
    joined_alignments = []
    nrp_fragments_splits = split_sequence_subseqs(nrp_fragments) if len(nrp_fragments) == 3 \
        else [[[nrp_fragment] for nrp_fragment in nrp_fragments]]  # TODO: do smth more meaningful for more than 3 fragments
    for nrp_fragments_split in nrp_fragments_splits:
        joined_alignment = []
        for nrp_fragments in nrp_fragments_split:
            monomer_sequences = fragments_joined_monomer_sequences(nrp_fragments)
            best_alignment = max((get_alignment(bgc_modules, monomers_sequence, dp_helper)
                                  for monomers_sequence in monomer_sequences),
                                 key=alignment_score)
            joined_alignment.append(best_alignment)
        joined_alignments.append(joined_alignment)

    return max(joined_alignments, key=combined_alignments_score)