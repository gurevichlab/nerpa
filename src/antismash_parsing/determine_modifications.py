from itertools import pairwise, chain
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
    STRAND,
    BGC_Cluster,
    BGC_Module_ID
)
from src.data_types import BGC_Module_Modification


def is_a_domain(domain_type: DomainType) -> bool:  # TODO: rename
    return domain_type in (DomainType.A, DomainType.PKS)


def ends_with_pcp_pcp(gene: Gene) -> bool:
    joined_domains = [domain_type for module in gene.modules
                      for domain_type in module.domains_sequence]
    last_a_index = next((i for i in reversed(range(len(joined_domains)))
                         if is_a_domain(joined_domains[i])),
                        None)
    if last_a_index is None:
        return False
    return any(domain1_type == DomainType.PCP and domain2_type == DomainType.PCP
               for domain1_type, domain2_type in pairwise(joined_domains[last_a_index:]))


def get_iterative_genes_orphan_c(genes: List[Gene]) -> Set[GeneId]:
    """ TODO: remove --- I don't need sorted assumption anymore
    genes_are_sorted = all(gene1.coords.start < gene2.coords.start
                           for gene1, gene2 in pairwise(genes))
    assert genes_are_sorted, "genes should be sorted by coordinates"
    """

    iterative_genes_ids = set()
    for gene1, gene2 in pairwise(genes):
        if gene1.coords.strand != gene2.coords.strand:
            continue
        if gene1.coords.strand == gene2.coords.strand == STRAND.REVERSE:
            gene1, gene2 = gene2, gene1
        if gene2.orphan_c_at_start:
            iterative_genes_ids.add(gene1.gene_id)
        if gene1.orphan_c_at_end:
            #and (not gene2.modules or any(domain.in_c_domain_group()
            #                              for domain in gene2.modules[0].domains_sequence)):
            iterative_genes_ids.add(gene1.gene_id)
    if genes[-1].orphan_c_at_end:
        iterative_genes_ids.add(genes[-1].gene_id)
    if genes[0].coords.strand == STRAND.REVERSE and genes[0].orphan_c_at_end:
        iterative_genes_ids.add(genes[0].gene_id)
    return iterative_genes_ids


def get_iterative_genes(genes: List[Gene]) -> Set[GeneId]:
    #genes = _genes[:]  # make a copy to avoid modifying the original list
    #genes.sort(key=lambda gene: (gene.coords.start, gene.coords.strand))

    iterative_genes_orphan_c = get_iterative_genes_orphan_c(genes)
    iterative_genes_pcp = {gene.gene_id for gene in genes if ends_with_pcp_pcp(gene)}

    return iterative_genes_orphan_c | iterative_genes_pcp



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
def get_iterative_modules_ids_gene(gene: Gene) -> Set[BGC_Module_ID]:
    joined_domains = [(domain_type, module_idx) for module_idx, module in enumerate(gene.modules)
                      for domain_type in module.domains_sequence]

    modules_idxs = set()
    split_domains = list(split_at(joined_domains, lambda p: is_a_domain(p[0]),
                                  keep_separator=True))[1:]  # skip the first (possibly empty) group before the first A domain
    for a_domain_group, interior_domains in grouper(split_domains, 2, fillvalue=[]):
        interior_domains_types = [domain_type for domain_type, module_idx in interior_domains]
        if has_pcp_condensation_pcp_subsequence(interior_domains_types):
            modules_idxs.add(a_domain_group[0][1])
    return {BGC_Module_ID(gene.gene_id, module_idx)
            for module_idx in modules_idxs}


def get_iterative_modules_ids(genes: List[Gene]) -> Set[BGC_Module_ID]:
    return {bgc_module_id
            for gene in genes
            for bgc_module_id in get_iterative_modules_ids_gene(gene)}


def get_modules_modifications(genes: List[Gene]) -> Dict[BGC_Module_ID, Tuple[BGC_Module_Modification, ...]]:
    modules_with_ids = [(BGC_Module_ID(gene.gene_id, module_idx), module)
                        for gene in genes
                        for module_idx, module in enumerate(gene.modules)]

    mods_by_module_id = {}
    for i in range(len(modules_with_ids)):
        module_id, module = modules_with_ids[i]
        module_domains = module.domains_sequence
        next_module_domains = modules_with_ids[i + 1][1].domains_sequence \
            if i < len(modules_with_ids) - 1 else None

        mods = []
        if DomainType.MT in module_domains:
            mods.append(BGC_Module_Modification.METHYLATION)
        if DomainType.E in module_domains or \
            (next_module_domains is not None and any([DomainType.C_DCL in next_module_domains,
                                                      DomainType.C_DUAL in next_module_domains])):
            mods.append(BGC_Module_Modification.EPIMERIZATION)

        mods_by_module_id[module_id] = tuple(mods)

    return mods_by_module_id

