from typing import List
from src.antismash_parsing.antismash_parser_types import (
    DomainType,
    Gene
)
from itertools import pairwise, dropwhile
from more_itertools import split_at, grouper

def is_a_domain(domain_type: DomainType) -> bool:  # TODO: rename
    return domain_type in (DomainType.A, DomainType.PKS)


def is_iterative_gene(gene: Gene) -> bool:
    joined_domains = [domain_type for module in gene.modules
                      for domain_type in module.domains_sequence]
    last_a_index = next((i for i in reversed(range(len(joined_domains)))
                         if is_a_domain(joined_domains[i])),
                        None)
    if last_a_index is None:
        return False
    return any(domain1_type == DomainType.PCP and domain2_type == DomainType.PCP
               for domain1_type, domain2_type in pairwise(joined_domains[last_a_index:]))


def has_pcp_condensation_pcp_subsequence(interior_domains_types: List[DomainType]) -> bool:
    domains_iter = iter(interior_domains_types)
    try:
        while next(domains_iter) != DomainType.PCP: pass
        while not next(domains_iter).in_c_domain_group(): pass
        while next(domains_iter) != DomainType.PCP: pass
    except StopIteration:
        return False
    return True


# TODO: quite messy with these indexes
def get_iterative_modules_idxs(gene: Gene) -> List[int]:
    joined_domains = [(domain_type, module_idx) for module_idx, module in enumerate(gene.modules)
                      for domain_type in module.domains_sequence]

    modules_idxs = []
    split_domains = list(split_at(joined_domains, lambda p: is_a_domain(p[0]),
                                  keep_separator=True))[1:]  # skip the first (possibly empty) group before the first A domain
    for a_domain_group, interior_domains in grouper(split_domains, 2, fillvalue=[]):
        interior_domains_types = [domain_type for domain_type, module_idx in interior_domains]
        if has_pcp_condensation_pcp_subsequence(interior_domains_types):
            modules_idxs.append(a_domain_group[0][1])
    return modules_idxs
