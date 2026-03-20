from typing import (
    List,
    NamedTuple,
    Tuple,
    Optional,
    Iterable,
    TYPE_CHECKING
)

from src.rban_parsing.nrp_variant_types import NRP_Variant, NRP_Fragment
if TYPE_CHECKING:
    from src.hmm.detailed_hmm import DetailedHMM
from src.monomer_names_helper import MonomerNamesHelper, MonCode, NRP_Monomer
from src.rban_parsing.rban_monomer import rBAN_Monomer, rBAN_idx
from src.generic.combinatorics import split_sequence_subseqs, filter_unique
from itertools import chain, permutations, product

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

    def to_mon_codes_json(self,
                          monomer_names_helper: MonomerNamesHelper,
                          any_hmm: Optional['DetailedHMM'] = None) -> dict:  # any_hmm is a stub to compute score_vs_avg_bgc. TODO: refactor
        return {
            'nrp_id': self.nrp_id,
            'non_iterative': [to_mon_codes(linearization, monomer_names_helper)
                              for linearization in self.non_iterative],
            'iterative': [[[to_mon_codes(linearization, monomer_names_helper)
                           for linearization in group]
                           for group in split]
                          for split in self.iterative],
            'score_vs_avg_bgc': any_hmm.score_vs_avg_bgc(self.non_iterative[0])
        }


def fragment_linearizations(fragment: NRP_Fragment) -> List[Linearization]:
    if fragment.is_cyclic:
        return [fragment.monomers[i:] + fragment.monomers[:i]
                for i in range(len(fragment.monomers))]
    else:
        return [fragment.monomers]


def joined_linearizations(fragments: List[NRP_Fragment],
                          max_num_fragments_to_permute: int = 5) -> Iterable[Linearization]:

    fragments_permutations = permutations(fragments) \
        if len(fragments) <= max_num_fragments_to_permute else [fragments]

    for fragments in fragments_permutations:
        for fragments_linearization in product(*(fragment_linearizations(fragment)
                                                 for fragment in fragments)):
            yield list(chain(*fragments_linearization))


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
            group_linearizations = joined_linearizations(nrp_fragments_group)
            one_split_linearizations.append(list(group_linearizations))
        all_splits_linearizations.append(one_split_linearizations)

    return all_splits_linearizations


def get_nrp_linearizations(nrp_variant: NRP_Variant) -> NRP_Linearizations:
    non_iterative = list(joined_linearizations(nrp_variant.fragments))
    iterative = iterative_fragments_linearizations(nrp_variant.fragments)
    return NRP_Linearizations(nrp_variant.nrp_variant_id.nrp_id, non_iterative, iterative)


def get_all_nrp_linearizations(nrp_variants: List[NRP_Variant]) -> List[NRP_Linearizations]:
    return list(filter(lambda lis: lis.iterative or lis.non_iterative,
                       (get_nrp_linearizations(nrp_variant) for nrp_variant in nrp_variants)))


def num_linearizations(nrp_linearizations: NRP_Linearizations) -> int:
    return len(nrp_linearizations.non_iterative) + sum(len(group_linearizations)
                                                       for group_linearizations in nrp_linearizations.iterative)


def _generate_linearizations_new(nrp_variant: NRP_Variant,
                                num_fragment_groups: int) -> Iterable[Tuple[Iterable[Linearization], ...]]:
    '''
    Returns: exactly num_fragment_groups iterators of linearizations
    '''
    assert 1 <= num_fragment_groups <= len(nrp_variant.fragments), \
        f'num_fragments should be in [1, {len(nrp_variant.fragments)}], got {num_fragment_groups}'

    for fragments_grouping in split_sequence_subseqs(nrp_variant.fragments,
                                                     num_subseqs=num_fragment_groups):
         yield tuple(joined_linearizations(group)
                     for group in fragments_grouping)


def generate_linearizations_new(nrp_variant: NRP_Variant,
                                num_fragment_groups: int) -> Iterable[Tuple[Linearization, ...]]:
    '''
    Returns: exactly num_fragment_groups iterators of linearizations
    '''
    assert 1 <= num_fragment_groups <= len(nrp_variant.fragments), \
        f'num_fragments should be in [1, {len(nrp_variant.fragments)}], got {num_fragment_groups}'

    for linearizations_for_grouping in _generate_linearizations_new(nrp_variant,
                                                                    num_fragment_groups=num_fragment_groups):
        for linearizations_tuple in product(*linearizations_for_grouping):
            yield linearizations_tuple

def generate_linearizations_uniq(nrp_variant: NRP_Variant,
                                 num_fragment_groups: int) -> Iterable[Tuple[Linearization, ...]]:

    def linearizations_key(linearizations_tuple: Tuple[Linearization, ...]) -> Tuple[Tuple[NRP_Monomer, ...], ...]:
        return tuple(sorted(tuple(mon.to_base_mon() for mon in linearization)
                            for linearization in linearizations_tuple))

    return filter_unique(generate_linearizations_new(nrp_variant,
                                                      num_fragment_groups=num_fragment_groups),
                         key=linearizations_key)




def generate_linearizations_no_fragmentation(nrp_variant: NRP_Variant) -> Iterable[Linearization]:
    return (linearizations[0]
            for linearizations in generate_linearizations_new(nrp_variant,
                                                              num_fragment_groups=1))
