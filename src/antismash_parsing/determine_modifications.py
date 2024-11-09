from itertools import pairwise
from typing import (
    Dict,
    List,
    Tuple,
    Set
)

from more_itertools import split_at, grouper
from src.antismash_parsing.antismash_parser_types import (
    DomainType,
    Gene,
    GeneId,
    STRAND
)
from src.data_types import BGC_Module_Modification


def is_a_domain(domain_type: DomainType) -> bool:  # TODO: rename
    return domain_type in (DomainType.A, DomainType.PKS)


def ends_with_pcp_pcp(gene: Gene) -> bool:
    if gene.is_iterative:
        return True
    joined_domains = [domain_type for module in gene.modules
                      for domain_type in module.domains_sequence]
    last_a_index = next((i for i in reversed(range(len(joined_domains)))
                         if is_a_domain(joined_domains[i])),
                        None)
    if last_a_index is None:
        return False
    return any(domain1_type == DomainType.PCP and domain2_type == DomainType.PCP
               for domain1_type, domain2_type in pairwise(joined_domains[last_a_index:]))


def get_iterative_genes_orphan_c(genes: List[Gene],  # should be sorted
                                 orphan_c_domains_per_gene: Dict[GeneId, Tuple[bool, bool]]) -> Set[GeneId]:
    iterative_genes_ids = set()
    for gene1, gene2 in pairwise(genes):
        if gene1.coords.strand != gene2.coords.strand:
            continue
        if gene1.coords.strand == gene2.coords.strand == STRAND.REVERSE:
            gene1, gene2 = gene2, gene1
        if orphan_c_domains_per_gene[gene2.gene_id][0]:
            iterative_genes_ids.add(gene1.gene_id)
        if orphan_c_domains_per_gene[gene1.gene_id][1]:
            #and (not gene2.modules or any(domain.in_c_domain_group()
            #                              for domain in gene2.modules[0].domains_sequence)):
            iterative_genes_ids.add(gene1.gene_id)
    if orphan_c_domains_per_gene[genes[-1].gene_id][1]:
        iterative_genes_ids.add(genes[-1].gene_id)
    if genes[0].coords.strand == STRAND.REVERSE and orphan_c_domains_per_gene[genes[0].gene_id][1]:
        iterative_genes_ids.add(genes[0].gene_id)
    return iterative_genes_ids


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


def get_modules_modifications(gene: Gene) -> Dict[int, Tuple[BGC_Module_Modification, ...]]:
    mods_by_module_idx = {}
    for module_idx, module in enumerate(gene.modules):
        next_module_domains = gene.modules[module_idx + 1].domains_sequence \
            if module_idx < len(gene.modules) - 1 else []
        mods = []
        if DomainType.MT in module.domains_sequence:
            mods.append(BGC_Module_Modification.METHYLATION)
        if any([DomainType.E in module.domains_sequence,
                DomainType.C_DCL in next_module_domains,
                DomainType.C_DUAL in next_module_domains]):
            mods.append(BGC_Module_Modification.EPIMERIZATION)

        mods_by_module_idx[module_idx] = tuple(mods)

    return mods_by_module_idx

