from typing import List
from src.antismash_parsing.antismash_parser_types import (
    BGC_Cluster,
    Gene,
    Module,
    DomainType,
    STRAND
)
from config import antiSMASH_Parsing_Config
from src.generic_algorithms import list_monad_compose
from functools import partial
from itertools import pairwise, permutations, chain
from more_itertools import split_after, split_at


def split_by_dist(bgc_cluster: BGC_Cluster,
                  config: antiSMASH_Parsing_Config) -> List[BGC_Cluster]:
    gene_groups = split_after(pairwise(bgc_cluster.genes),
                              lambda p: p[1].coords.start - p[0].coords.end > config.MAX_DISTANCE_BETWEEN_GENES)
    return [BGC_Cluster(genome_id=bgc_cluster.genome_id,
                        contig_id=bgc_cluster.contig_id,
                        bgc_id=f'{bgc_cluster.bgc_id}_{i}',
                        genes=[gene for gene, _ in group])
            for i, group in enumerate(gene_groups)]


def split_by_single_gene_Starter_TE(bgc_cluster: BGC_Cluster) -> List[BGC_Cluster]:
    gene_groups = split_at(bgc_cluster.genes,
                           lambda gene: DomainType.C_STARTER in gene.modules[0].domains_sequence and
                                        DomainType.TE_TD in gene.modules[-1].domains_sequence,
                           keep_separator=True)
    return [BGC_Cluster(genome_id=bgc_cluster.genome_id,
                        contig_id=bgc_cluster.contig_id,
                        bgc_id=f'{bgc_cluster.bgc_id}_{i}',
                        genes=group)
            for i, group in enumerate(gene_groups)]

def a_pcp_module(module: Module) -> bool:
    return set(module.domains_sequence) in ({DomainType.A, DomainType.PCP},
                                            {DomainType.PKS, DomainType.PCP})

def genes_sequence_consistent(genes: List[Gene]) -> bool:
    joined_modules = [module for gene in genes for module in gene.modules]

    at_least_one_A_domain = any(DomainType.A in module.domains_sequence for module in joined_modules)
    nc_terms_consistent = DomainType.NTERM not in joined_modules[0].domains_sequence and \
                          DomainType.CTERM not in joined_modules[-1].domains_sequence

    c_starter_consistent = all(DomainType.C_STARTER not in module.domains_sequence
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


def get_genes_rearrangements(_genes: List[Gene]) -> List[List[Gene]]:
    genes = reverse_if_all_neg(_genes)
    starting_gene = [gene for gene in genes
                     if DomainType.C_STARTER in gene.modules[0].domains_sequence or a_pcp_module(gene.modules[0])]
    terminal_gene = [gene for gene in genes
                     if DomainType.TE_TD in gene.modules[-1].domains_sequence]
    interior_genes = [gene for gene in genes
                      if gene not in starting_gene + terminal_gene]
    genes = starting_gene + interior_genes + terminal_gene
    if genes_sequence_consistent(genes):
        return [genes]
    else:   # TODO: optimize, keep only several top permutations
        return [starting_gene + permuted_interior_genes + terminal_gene
                for permuted_interior_genes in map(list, permutations(interior_genes))
                if genes_sequence_consistent(starting_gene + permuted_interior_genes + terminal_gene)]

def split_and_reorder_inconsistent(bgc: BGC_Cluster) -> List[BGC_Cluster]:
    genes = reverse_if_all_neg(bgc.genes)
    if genes_sequence_consistent(genes):
        return [BGC_Cluster(genome_id=bgc.genome_id,
                            contig_id=bgc.contig_id,
                            bgc_id=bgc.bgc_id,
                            genes=genes)]

    genes_substrings = (genes[i:j]
                        for i in range(len(genes))
                        for j in range(i + 1, len(genes)))
    genes_rearrangements = chain(*map(get_genes_rearrangements, genes_substrings))

    return [BGC_Cluster(genome_id=bgc.genome_id,
                        contig_id=bgc.contig_id,
                        bgc_id=f'{bgc.bgc_id}_{i}',
                        genes=rearranged_genes)
            for i, rearranged_genes in enumerate(genes_rearrangements)]


def split_and_reorder(bgc: BGC_Cluster,
                      config: antiSMASH_Parsing_Config) -> List[BGC_Cluster]:
    return list_monad_compose(partial(split_by_dist, config=config),
                              split_by_single_gene_Starter_TE,
                              split_and_reorder_inconsistent)([bgc])
