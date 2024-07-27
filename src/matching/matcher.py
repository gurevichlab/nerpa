from typing import (
    Iterable,
    List,
    Optional,
    Union
)
from src.data_types import (
    BGC_Module,
    BGC_Variant,
    LogProb,
    NRP_Fragment,
    NRP_Variant
)
from src.matching.scoring_helper import ScoringHelper
from src.matching.alignment_types import (
    Alignment,
    alignment_score,
    combined_alignments_score,
    Match,
    show_alignment
)
from src.matching.dp import get_alignment
from src.matching.multiple_fragments_alignment_handling import get_multiple_fragments_alignment, get_iterative_nrp_alignment

from itertools import chain, islice, takewhile, permutations
from joblib import delayed, Parallel


def get_alignment_fragment(assembly_line: List[BGC_Module],
                           nrp_fragment: NRP_Fragment,
                           dp_helper: ScoringHelper) -> Alignment:

    nrp_cyclic_shifts = [nrp_fragment.monomers] if not nrp_fragment.is_cyclic \
            else [nrp_fragment.monomers[i:] + nrp_fragment.monomers[:i]
                  for i in range(len(nrp_fragment.monomers))]

    alignments = [get_alignment(assembly_line, nrp_cyclic_shift, dp_helper)
                  for nrp_cyclic_shift in nrp_cyclic_shifts]

    return max(alignments, key=alignment_score)


def null_hypothesis_score(nrp_fragment: NRP_Fragment,
                          dp_helper: ScoringHelper) -> LogProb:
    return sum(dp_helper.null_hypothesis_score(mon)
               for mon in nrp_fragment.monomers)


def get_normalized_score(bgc_variant: BGC_Variant,
                         nrp_variant: NRP_Variant,
                         alignments: List[Alignment],
                         dp_helper: ScoringHelper) -> LogProb:
    combined_score = combined_alignments_score(alignments)

    match dp_helper.scoring_config.normalization:
        case 'RANDOM_NULL_MODEL':
            null_score = sum(null_hypothesis_score(fragment, dp_helper) for fragment in nrp_variant.fragments)
            normalized_score = combined_score - null_score
        case 'AVERAGING':
            total_length = sum(len(fragment.monomers) for fragment in nrp_variant.fragments)
            normalized_score = combined_score / total_length
        case 'NONE':
            normalized_score = combined_score
        case _:
            raise ValueError(f'Unknown normalization method: {dp_helper.scoring_config.normalization}')
    return normalized_score


def get_match(bgc_variant: BGC_Variant,
              nrp_variant: NRP_Variant,
              dp_helper: ScoringHelper) -> Match:
    dp_helper.set_pks_domains_in_bgc(bgc_variant.has_pks_domains)  # I don't really like that dp_helper carries some state

    multiple_fragments_alignments = (get_multiple_fragments_alignment(bgc_variant.fragments,
                                                                      list(nrp_fragments),
                                                                      dp_helper)
                                     for nrp_fragments in permutations(nrp_variant.fragments))
    iterative_nrp_alignments = get_iterative_nrp_alignment(bgc_variant.fragments,
                                                           nrp_variant.fragments,
                                                           dp_helper)

    best_fragments_alignment = max(chain(multiple_fragments_alignments, [iterative_nrp_alignments]),
                                   key=combined_alignments_score)
    normalized_score = get_normalized_score(bgc_variant, nrp_variant,
                                            best_fragments_alignment,
                                            dp_helper)
    return Match(bgc_variant,
                 nrp_variant,
                 best_fragments_alignment,
                 normalized_score)


def get_matches_for_bgc_variant(bgc_variant: BGC_Variant,
                                nrp_variants: List[NRP_Variant],
                                scoring_helper: ScoringHelper,
                                min_score: Optional[float] = None,
                                max_num_matches: Optional[int] = None,
                                log=None) -> List[Match]:
    if log is not None:
        log.info(f'Processing BGC variant {bgc_variant.variant_idx}')

    matches = sorted((get_match(bgc_variant, nrp_variant, scoring_helper)
                     for nrp_variant in nrp_variants),
                     key=lambda m: m.normalized_score, reverse=True)

    return list(takewhile(lambda m: (m.normalized_score > min_score) if min_score is not None else True,
                          islice(matches, max_num_matches)))


def get_matches(bgc_variants: List[BGC_Variant],
                nrp_variants: List[NRP_Variant],
                scoring_helper: ScoringHelper,
                min_score: Optional[float] = None,
                max_num_matches: Optional[int] = None,
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
