from typing import (
    Dict,
    Iterable,
    List,
    Optional,
    Tuple
)
from src.data_types import BGC_Fragment, NRP_Monomer, NRP_Fragment
from src.matching.alignment_types import Alignment, combined_alignments_score
from src.matching.scoring_helper import ScoringHelper
from src.matching.dp import LogProb
from src.matching.dp import get_alignment
from src.matching.alignment_types import alignment_score
from itertools import chain, product


def get_fragments_monomer_sequences(fragments: List[NRP_Fragment]) -> Iterable[List[NRP_Monomer]]:
    fragment_sequences = [[fragment.monomers[i:] + fragment.monomers[:i] for i in range(len(fragment.monomers))]
                          if fragment.is_cyclic else [fragment.monomers]
                          for fragment in fragments]
    return (list(chain(*monomer_sequences))
            for monomer_sequences in product(*fragment_sequences))


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


def get_multiple_fragments_alignment(bgc_fragments: List[BGC_Fragment],
                                     nrp_fragments: List[NRP_Fragment],
                                     dp_helper: ScoringHelper) -> List[Alignment]:
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
        nrp_monomers_sequences = get_fragments_monomer_sequences(nrp_fragments[nrp_start: nrp_end])
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


def get_iterative_nrp_alignment(bgc_fragments: List[BGC_Fragment],
                                nrp_fragments: List[NRP_Fragment],
                                dp_helper: ScoringHelper) -> List[Alignment]:
    def get_fragment_alignment(bgc_fragment: BGC_Fragment,
                               nrp_fragment: NRP_Fragment) -> Alignment:
        return max((get_alignment(bgc_fragment, monomers_sequence, dp_helper)
                    for monomers_sequence in get_fragments_monomer_sequences([nrp_fragment])),
                   key=alignment_score)

    return max(([get_fragment_alignment(bgc_fragment, nrp_fragment) for nrp_fragment in nrp_fragments]
               for bgc_fragment in bgc_fragments),
               key=combined_alignments_score)