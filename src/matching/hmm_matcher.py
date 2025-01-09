from typing import List, Optional
from src.data_types import BGC_Variant, NRP_Variant, NRP_Fragment
from src.matching.matcher_viterbi_detailed_hmm import DetailedHMM
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.matching.matching_types_match import Match, Match_BGC_Variant_Info, Match_NRP_Variant_Info
from src.matching.matching_types_alignment import Alignment, alignment_score, combined_alignments_score
from src.generic.combinatorics import split_sequence_subseqs
from itertools import chain, permutations, product
from joblib import delayed, Parallel


def fragments_linearizations(fragments: List[NRP_Fragment],
                             max_num_fragments_to_permute: int = 3) -> List[List[rBAN_Monomer]]:
    def fragment_linearizations(fragment: NRP_Fragment) -> List[List[rBAN_Monomer]]:
        if fragment.is_cyclic:
            return [fragment.monomers[i:] + fragment.monomers[:i]
                    for i in range(len(fragment.monomers))]
        else:
            return [fragment.monomers]

    fragments_permutations = permutations(fragments) \
        if len(fragments) <= max_num_fragments_to_permute else [fragments]

    linearizations = []
    for fragments in fragments_permutations:
        for fragments_linearization in product(*(fragment_linearizations(fragment)
                                                 for fragment in fragments)):
            linearizations.append(list(chain(*fragments_linearization)))
    return linearizations


def get_iterative_bgc_alignment(hmm: DetailedHMM,
                                nrp_fragments: List[NRP_Fragment]) -> List[Alignment]:
    joined_alignments = []
    nrp_fragments_splits = split_sequence_subseqs(nrp_fragments) if len(nrp_fragments) == 3 \
        else [[[nrp_fragment] for nrp_fragment in nrp_fragments]]  # TODO: do smth more meaningful for more than 3 fragments
    for nrp_fragments_split in nrp_fragments_splits:
        joined_alignment = []
        for nrp_fragments in nrp_fragments_split:
            linearizations = fragments_linearizations(nrp_fragments)
            best_alignment = max((hmm.get_alignment(linearization) for linearization in linearizations),
                                 key=alignment_score)
            joined_alignment.append(best_alignment)
        joined_alignments.append(joined_alignment)

    return max(joined_alignments, key=combined_alignments_score)


def get_match(hmm: DetailedHMM,
              nrp_variant: NRP_Variant) -> Optional[Match]:
    non_iterative_alignment = max((hmm.get_alignment(nrp_linearization)
                                   for nrp_linearization in fragments_linearizations(nrp_variant.fragments)),
                                  key=alignment_score)
    iterative_alignment = get_iterative_bgc_alignment(hmm, nrp_variant.fragments)
    best_alignments = max([non_iterative_alignment], iterative_alignment, key=combined_alignments_score)
    score = combined_alignments_score(best_alignments)
    return Match(bgc_variant_info=Match_BGC_Variant_Info.from_bgc_variant(hmm.bgc_variant),
                 nrp_variant_info=Match_NRP_Variant_Info.from_nrp_variant(nrp_variant),
                 alignments=best_alignments,
                 normalized_score=score)


def get_matches_for_hmm(hmm: DetailedHMM,
                        nrp_variants: List[NRP_Variant],
                        max_num_matches_per_bgc_variant: Optional[int],
                        log=None) -> List[Match]:
    if log is not None:
        log.info(f'Processing BGC {hmm.bgc_variant.genome_id} variant {hmm.bgc_variant.variant_idx}')
    matches = sorted(filter(None, (get_match(hmm, nrp_variant)
                                   for nrp_variant in nrp_variants)),
                     key=lambda m: m.normalized_score, reverse=True)
    return matches[:max_num_matches_per_bgc_variant]



def get_matches(hmms: List[DetailedHMM],
                nrp_variants: List[NRP_Variant],
                max_num_matches_per_bgc_variant: Optional[int] = None,
                num_threads: int = 1,
                log=None) -> List[Match]:
    if log is not None:
        log.info(f'Matching {len(hmms)} BGC variants against {len(nrp_variants)} NRP variants')
    matches = chain(*Parallel(n_jobs=num_threads)(delayed(get_matches_for_hmm)(hmm, nrp_variants,
                                                                               max_num_matches_per_bgc_variant,
                                                                               log)
                                                  for hmm in hmms))
    # q: sort matches by normalized score
    return sorted(matches, key=lambda m: m.normalized_score, reverse=True)
