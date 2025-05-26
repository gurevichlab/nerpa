from typing import (
    Dict,
    List,
    NamedTuple,
    Tuple
)
from src.data_types import NRP_Variant, NRP_Fragment
from src.monomer_names_helper import MonomerNamesHelper, MonCode
from src.rban_parsing.rban_monomer import rBAN_Monomer, rBAN_idx
from src.generic.combinatorics import split_sequence_subseqs
from src.matching.match_type import NRP_Variant_ID
from itertools import chain, permutations, product
from io import StringIO


Linearization = List[rBAN_Monomer]
LinearizationLight = Tuple[List[MonCode], List[rBAN_idx]]

def to_mon_codes(linearization: Linearization,
                 monomer_names_helper: MonomerNamesHelper) -> LinearizationLight:
    mon_codes = [monomer_names_helper.mon_to_int[mon.to_base_mon()]
                 for mon in linearization]
    rban_idxs = [mon.rban_idx for mon in linearization]
    return mon_codes, rban_idxs

class NRP_Linearizations(NamedTuple):
    nrp_id: str
    non_iterative: List[Linearization]
    iterative: List[List[List[Linearization]]]
    # fragments are split into groups which are then linearized (and aligned) separately
    # linearizations of one group -> List[Linearization]
    # linearizations of all group members -> List[List[Linearization]]
    # linearizations of all groups -> List[List[List[Linearization]]]

    def to_mon_codes_json(self, monomer_names_helper: MonomerNamesHelper) -> dict:
        return {
            'nrp_id': self.nrp_id,
            'non_iterative': [to_mon_codes(linearization, monomer_names_helper)
                              for linearization in self.non_iterative],
            'iterative': [[[to_mon_codes(linearization, monomer_names_helper)
                           for linearization in group]
                           for group in split]
                          for split in self.iterative]
        }


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
    if len(nrp_fragments) <= 1:
        return []
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
    return NRP_Linearizations(nrp_variant.nrp_variant_id.nrp_id, non_iterative, iterative)


def get_all_nrp_linearizations(nrp_variants: List[NRP_Variant]) -> List[NRP_Linearizations]:
    return list(filter(lambda lis: lis.iterative or lis.non_iterative,
                       (get_nrp_linearizations(nrp_variant) for nrp_variant in nrp_variants)))


def num_linearizations(nrp_linearizations: NRP_Linearizations) -> int:
    return len(nrp_linearizations.non_iterative) + sum(len(group_linearizations)
                                                       for group_linearizations in nrp_linearizations.iterative)