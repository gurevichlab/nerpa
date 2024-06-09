#!/usr/bin/env python

import src.nerpa_pipeline.handle_helper as handle_helper
from itertools import pairwise, dropwhile
from more_itertools import split_at

def is_a_domain(domain_type: DomainType) -> bool:  # TODO: rename
    return domain_type in (DomainType.A, DomainType.PKS)


def is_iterative_gene(gene: Gene) -> bool:
    joined_domains = [domain_type for module in gene.modules
                      for domain_type in module.domains_sequence]
    last_a_index = next(i for i in reversed(range(len(joined_domains)))
                        if is_a_domain(joined_domains[i]))
    return any(domain1_type == DomainType.PCP and domain2_type == DomainType.PCP
               for domain1_type, domain2_type in pairwise(joined_domains[last_a_index:]))


def has_pcp_condensation_pcp_subsequence(interior_domains_types: List[DomainType]) -> bool:
    domains_iter = iter(interior_domains_types)
    try:
        dropwhile(domains_iter,
                  lambda domain_type: domain_type != DomainType.PCP)
        dropwhile(domains_iter,
                  lambda domain_type: domain_type != DomainType.C)
        dropwhile(domains_iter,
                  lambda domain_type: domain_type != DomainType.PCP)
    except StopIteration:
        return False
    return True


# TODO: quite messy with these indexes
def iterative_modules_idxs(gene: Gene) -> List[int]:
    joined_domains = [(domain_type, module_idx) for module_idx, module in enumerate(gene.modules)
                      for domain_type in module.domains_sequence]
    fst_a_index = next(i for i, (domain_type, module_idx) in enumerate(joined_domains)
                       if is_a_domain(domain_type))

    modules_idxs = []
    for a_domain_group, interior_domains in split_at(joined_domains[fst_a_index:],
                                                     lambda p: is_a_domain(p[0])):
        interior_domains_types = [domain_type for domain_type, module_idx in interior_domains]
        if has_pcp_condensation_pcp_subsequence(interior_domains_types):
            modules_idxs.append(a_domain_group[0][1])
    return modules_idxs
