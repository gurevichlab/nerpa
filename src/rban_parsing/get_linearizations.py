from typing import (
    Dict,
    List,
    NamedTuple,
)
from src.data_types import NRP_Variant, NRP_Fragment
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.generic.combinatorics import split_sequence_subseqs
from src.matching.match_type import Match_NRP_Variant_Info
from itertools import chain, permutations, product


Linearization = List[rBAN_Monomer]

class NRP_Linearizations(NamedTuple):
    non_iterative: List[Linearization]
    iterative: List[List[List[Linearization]]]
    # fragments are split into groups which are then linearized (and aligned) separately
    # linearizations of one group -> List[Linearization]
    # linearizations of all group members -> List[List[Linearization]]
    # linearizations of all groups -> List[List[List[Linearization]]]

def non_iterative_linearizations(fragments: List[NRP_Fragment],
                                 max_num_fragments_to_permute: int = 3) -> List[Linearization]:
    def fragment_linearizations(fragment: NRP_Fragment) -> List[Linearization]:
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


# fragments are split into groups which are then linearized (and aligned) separately
def iterative_fragments_linearizations(nrp_fragments: List[NRP_Fragment]) -> List[List[List[Linearization]]]:
    splits_into_groups = split_sequence_subseqs(nrp_fragments) if len(nrp_fragments) == 3 \
        else [[[nrp_fragment] for nrp_fragment in nrp_fragments]]  # TODO: do smth more meaningful for more than 3 fragments
    all_splits_linearizations = []
    for nrp_fragments_split in splits_into_groups:
        one_split_linearizations: List[List[Linearization]] = []
        for nrp_fragments_group in nrp_fragments_split:
            group_linearizations = non_iterative_linearizations(nrp_fragments_group)
            one_split_linearizations.append(group_linearizations)
        all_splits_linearizations.append(one_split_linearizations)

    return all_splits_linearizations


def get_nrp_linearizations(nrp_variant: NRP_Variant) -> NRP_Linearizations:
    non_iterative = non_iterative_linearizations(nrp_variant.fragments)
    iterative = iterative_fragments_linearizations(nrp_variant.fragments)
    return NRP_Linearizations(non_iterative, iterative)


def get_all_nrp_linearizations(nrp_variants: List[NRP_Variant]) -> Dict[Match_NRP_Variant_Info, NRP_Linearizations]:  # nrp_id -> List[NRP_Linearizations]
    return {Match_NRP_Variant_Info.from_nrp_variant(nrp_variant): get_nrp_linearizations(nrp_variant)
            for nrp_variant in nrp_variants}
    '''
    nrp_linearizations = {nrp_id: NRP_Linearizations([], [])
                          for nrp_id in set(nrp_variant.nrp_id for nrp_variant in nrp_variants)}
    for nrp_variant in nrp_variants:
        variant_linearizations = get_nrp_linearizations(nrp_variant)
        nrp_linearizations[nrp_variant.nrp_id].non_iterative.extend(variant_linearizations.non_iterative)
        nrp_linearizations[nrp_variant.nrp_id].iterative.extend(variant_linearizations.iterative)
    return nrp_linearizations
    '''