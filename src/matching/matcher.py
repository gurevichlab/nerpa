from typing import Iterable, List, Union
from src.data_types import (
    BGC_Module,
    BGC_Variant,
    LogProb,
    NRP_Fragment,
    NRP_Variant
)
from src.matching.scoring_helper import ScoringHelper
from src.matching.alignment_types import Alignment, alignment_score, Match, show_alignment
from src.matching.dp import get_alignment

from itertools import chain, islice, takewhile
from joblib import delayed, Parallel


def get_alignment_fragment(assembly_line: List[BGC_Module],
                           nrp_fragment: NRP_Fragment,
                           dp_helper: ScoringHelper) -> Alignment:

    nrp_cyclic_shifts = [nrp_fragment.monomers] if not nrp_fragment.is_cyclic \
            else [nrp_fragment.monomers[i:] + nrp_fragment.monomers[:i]
                  for i in range(len(nrp_fragment.monomers))]

    alignments = [get_alignment(assembly_line, nrp_cyclic_shift, dp_helper)
                  for nrp_cyclic_shift in nrp_cyclic_shifts]
    with open('alignments.txt', 'a') as f:
        for alignment in alignments:
            f.write(show_alignment(alignment) + '\n\n')

    return max(alignments, key=alignment_score)


def null_hypothesis_score(nrp_fragment: NRP_Fragment,
                          dp_helper: ScoringHelper) -> LogProb:
    return sum(dp_helper.null_hypothesis_score(mon)
               for mon in nrp_fragment.monomers)


def get_match(bgc_variant: BGC_Variant,
              nrp_variant: NRP_Variant,
              dp_helper: ScoringHelper) -> Match:
    fragment_alignments = [get_alignment_fragment(bgc_variant.tentative_assembly_line,
                                                  bgc_fragment,
                                                  dp_helper)
                           for bgc_fragment in nrp_variant.fragments]
    final_score = sum(alignment_score(alignment) - null_hypothesis_score(fragment, dp_helper)
                      for fragment, alignment in zip(nrp_variant.fragments, fragment_alignments))
    return Match(bgc_variant,
                 nrp_variant,
                 fragment_alignments,
                 final_score)


def get_matches_for_bgc_variant(bgc_variant: BGC_Variant,
                                nrp_variants: List[NRP_Variant],
                                scoring_helper: ScoringHelper,
                                min_score: float = 0,
                                max_num_matches: Union[int, None] = None,
                                log=None) -> List[Match]:
    if log is not None:
        log.info(f'Processing BGC variant {bgc_variant.variant_idx}')

    matches = sorted((get_match(bgc_variant, nrp_variant, scoring_helper)
                     for nrp_variant in nrp_variants),
                     key=lambda m: m.normalized_score, reverse=True)

    return list(takewhile(lambda m: m.normalized_score > min_score,
                          islice(matches, max_num_matches)))


def get_matches(bgc_variants: List[BGC_Variant],
                nrp_variants: List[NRP_Variant],
                scoring_helper: ScoringHelper,
                min_score: float = 0,
                max_num_matches: Union[int, None]=None,
                num_threads: int = 1,
                log=None) -> List[Match]:
    if log is not None:
        log.info(f'Matching {len(bgc_variants)} BGC variants against {len(nrp_variants)} NRP variants')
    matches = chain(*Parallel(n_jobs=num_threads)(delayed(get_matches_for_bgc_variant)(bgc_variant, nrp_variants,
                                                                                       scoring_helper,
                                                                                       min_score, max_num_matches,
                                                                                       log)
                                                  for bgc_variant in bgc_variants))
    return sorted(matches, key=lambda m: m.normalized_score, reverse=True)
