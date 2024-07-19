from typing import List
from src.antismash_parsing.antismash_parser_types import (
    BGC_Cluster,
    Coords,
    Gene,
    Module,
    DomainType,
    STRAND
)
from src.config import antiSMASH_Parsing_Config
from src.generic.generate_permutations import generate_permutations
from src.generic.functional import list_monad_compose
from functools import partial
from itertools import chain, islice, pairwise
from more_itertools import split_before, split_at


def split_by_dist(bgc_cluster: BGC_Cluster,
                  config: antiSMASH_Parsing_Config) -> List[BGC_Cluster]:
    dummy_gene = Gene('', Coords(start=-1, end=-1, strand=STRAND.FORWARD), [])  # dummy gene to start the pairwise
    gene_groups = list(split_before(pairwise([dummy_gene] + bgc_cluster.genes),
                                    lambda p: p[1].coords.start - p[0].coords.end > config.MAX_DISTANCE_BETWEEN_GENES))
    return [BGC_Cluster(genome_id=bgc_cluster.genome_id,
                        contig_id=bgc_cluster.contig_id,
                        bgc_idx=bgc_cluster.bgc_idx,
                        genes=[gene for _, gene in group])
            for group in gene_groups]


def split_by_single_gene_Starter_TE(bgc_cluster: BGC_Cluster,
                                    config: antiSMASH_Parsing_Config) -> List[BGC_Cluster]:  # config is not used but added for consistency
    gene_groups = split_at(bgc_cluster.genes,
                           lambda gene: DomainType.C_STARTER in gene.modules[0].domains_sequence and
                                        DomainType.TE_TD in gene.modules[-1].domains_sequence,
                           keep_separator=True)
    gene_groups = [list(group) for group in gene_groups]
    return [BGC_Cluster(genome_id=bgc_cluster.genome_id,
                        contig_id=bgc_cluster.contig_id,
                        bgc_idx=bgc_cluster.bgc_idx,
                        genes=group)
            for group in gene_groups
            if len(group) > 0]

def a_pcp_module(module: Module) -> bool:
    domains_set = set(module.domains_sequence)
    return domains_set == {DomainType.PKS, DomainType.PCP} or \
        {DomainType.A, DomainType.PCP}.issubset(domains_set) and all(not DomainType.in_c_domain_group(domain)
                                                                     for domain in domains_set)

def genes_sequence_consistent(genes: List[Gene]) -> bool:
    joined_modules = [module for gene in genes for module in gene.modules]

    at_least_one_A_domain = any(DomainType.A in module.domains_sequence for module in joined_modules)
    nc_terms_consistent = DomainType.NTERM not in joined_modules[0].domains_sequence and \
                          DomainType.CTERM not in joined_modules[-1].domains_sequence

    c_starter_consistent = all(DomainType.C_STARTER not in module.domains_sequence  # TODO: refactor: create one function is_starting_module for all cases
                               for module in joined_modules[1:])
    a_pcp_consistent = all(not a_pcp_module(module) for module in joined_modules[1:])
    te_td_consistent = all(DomainType.TE_TD not in module.domains_sequence
                           for module in joined_modules[:-1])

    return all([at_least_one_A_domain,
                a_pcp_consistent,
                c_starter_consistent,
                nc_terms_consistent,
                te_td_consistent])


def reverse_if_all_neg(genes: List[Gene]) -> List[Gene]:
    if all(gene.coords.strand == STRAND.REVERSE for gene in genes):
        return genes[::-1]
    else:
        return genes[:]


def get_genes_rearrangements(_genes: List[Gene], config: antiSMASH_Parsing_Config) -> List[List[Gene]]:
    genes = reverse_if_all_neg(_genes)
    if genes_sequence_consistent(genes):
        return [genes]


    starting_gene = [gene for gene in genes
                     if DomainType.C_STARTER in gene.modules[0].domains_sequence or a_pcp_module(gene.modules[0])]
    terminal_gene = [gene for gene in genes
                     if DomainType.TE_TD in gene.modules[-1].domains_sequence]
    interior_genes = [gene for gene in genes
                      if gene not in starting_gene + terminal_gene]
    return [starting_gene + permuted_interior_genes + terminal_gene
            for permuted_interior_genes in islice(generate_permutations(interior_genes),
                                                  config.MAX_PERMUTATIONS_PER_BGC)
            if genes_sequence_consistent(starting_gene + permuted_interior_genes + terminal_gene)]


def split_and_reorder_inconsistent(bgc: BGC_Cluster, config: antiSMASH_Parsing_Config) -> List[BGC_Cluster]:
    genes = reverse_if_all_neg(bgc.genes)
    if genes_sequence_consistent(genes):
        return [BGC_Cluster(genome_id=bgc.genome_id,
                            contig_id=bgc.contig_id,
                            bgc_idx=bgc.bgc_idx,
                            genes=genes)]

    # remove genes that do not have A domains because they are not relevant for the permutation
    # (the decision to apply all possible rearrangements is already made at this point
    # be careful with this if you ever want to consider other scheme of gene rearrangements
    genes = [gene for gene in genes
             if any(DomainType.A in module.domains_sequence for module in gene.modules)]

    genes_substrings = (genes[i:j]
                        for i in range(len(genes))
                        for j in range(i + 1, len(genes) + 1))
    genes_rearrangements = chain(*map(partial(get_genes_rearrangements, config=config),
                                      genes_substrings))

    return [BGC_Cluster(genome_id=bgc.genome_id,
                        contig_id=bgc.contig_id,
                        bgc_idx=bgc.bgc_idx,
                        genes=rearranged_genes)
            for rearranged_genes in genes_rearrangements]


def split_and_reorder(bgc: BGC_Cluster,
                      config: antiSMASH_Parsing_Config) -> List[BGC_Cluster]:
    return list_monad_compose(*(partial(split_func, config=config)
                                 for split_func in [split_by_dist,
                                                    split_by_single_gene_Starter_TE,
                                                    split_and_reorder_inconsistent]))([bgc])
