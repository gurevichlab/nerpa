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
from src.matching.dp import get_alignment
from src.matching.scoring_helper import ScoringHelper
from src.matching.matching_types_alignment import (
    Alignment,
    alignment_score,
    combined_alignments_score,
    show_alignment
)
from src.matching.matching_types_match import Match, Match_BGC_Variant_Info, Match_NRP_Variant_Info
from src.matching.multiple_fragments_alignment_handling import (
    get_multiple_fragments_alignment,
    get_iterative_bgc_alignment,
    get_fragments_alignment,
    fragments_joined_monomer_sequences
)
from src.matching.heuristic_matching import filter_out_bad_nrp_candidates, get_heuristic_discard, DiscardVerdict
from itertools import chain, islice, takewhile, permutations, product
from joblib import delayed, Parallel
from datetime import datetime


def null_hypothesis_score(nrp_fragment: NRP_Fragment,
                          dp_helper: ScoringHelper) -> LogProb:
    return sum(dp_helper.null_hypothesis_score(mon)
               for mon in nrp_fragment.monomers)


def get_normalized_score(bgc_variant: BGC_Variant,
                         nrp_variant: NRP_Variant,
                         alignments: List[Alignment],
                         dp_helper: ScoringHelper) -> LogProb:
    combined_score = combined_alignments_score(alignments)

    '''
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
    '''
    return combined_score


def get_match(bgc_variant: BGC_Variant,
              nrp_variant: NRP_Variant,
              dp_helper: ScoringHelper) -> Optional[Match]:  # match can be discarded based on heuristics
    dp_helper.set_pks_domains_in_bgc(bgc_variant.has_pks_domains())  # I don't really like that dp_helper carries some state

    fragments_alignments: List[List[Alignment]] = []

    heuristic_discard = get_heuristic_discard(bgc_variant, nrp_variant, dp_helper.scoring_config.heuristic_matching_cfg) \
        if dp_helper.heuristic_discard_on else DiscardVerdict(False, False)

    if not heuristic_discard.joined:
        fragments_permutations = permutations(nrp_variant.fragments) \
            if len(nrp_variant.fragments) <= 3 else [nrp_variant.fragments]
        monomer_sequences = chain(*(fragments_joined_monomer_sequences(fragments_permutation)
                                    for fragments_permutation in fragments_permutations))
        fragments_alignments.extend([get_alignment(bgc_variant.modules,
                                                   monomer_sequence,
                                                   dp_helper)]
                                    for monomer_sequence in monomer_sequences)
        # TODO: this is a stub, we should implement a more reasonable way to handle many fragments
        #  instead of just cutting off at some point
        # (e.g. by first finding the best bgc_fragment for each nrp_fragment to retrieve their order)
    if not heuristic_discard.iterative:
        fragments_alignments.append(get_iterative_bgc_alignment(bgc_variant.modules,
                                                                nrp_variant.fragments,
                                                                dp_helper))

    if not fragments_alignments:
        return None
    best_fragments_alignment = max(fragments_alignments, key=combined_alignments_score)
    normalized_score = get_normalized_score(bgc_variant, nrp_variant,
                                            best_fragments_alignment,
                                            dp_helper)
    return Match(bgc_variant_info=Match_BGC_Variant_Info.from_bgc_variant(bgc_variant),
                 nrp_variant_info=Match_NRP_Variant_Info.from_nrp_variant(nrp_variant),
                 alignments=best_fragments_alignment,
                 normalized_score=normalized_score)


def get_matches_for_bgc_variant(bgc_variant: BGC_Variant,
                                nrp_variants: List[NRP_Variant],
                                scoring_helper: ScoringHelper,
                                min_score: Optional[float] = None,
                                max_num_matches: Optional[int] = None,
                                log=None) -> List[Match]:
    if log is not None:
        log.info(f'Processing BGC {bgc_variant.genome_id} variant {bgc_variant.variant_idx}')
        log.info(f'Current time: {datetime.now()}')


    matches = sorted(filter(None, (get_match(bgc_variant, nrp_variant, scoring_helper)
                                   for nrp_variant in nrp_variants)),
                            key=lambda m: m.normalized_score, reverse=True)

    if log is not None:
        log.info(f'Processing complete. '
                 f'Discarded {len(nrp_variants) - len(matches)} out of {len(nrp_variants)} NRP variants '
                 f'based on heuristics')

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
    # q: sort matches by normalized score
    return sorted(matches, key=lambda m: m.normalized_score, reverse=True)

    #return sorted(matches, key=lambda m: m.normalized_score, reverse=True)

