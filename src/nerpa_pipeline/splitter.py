#!/usr/bin/env python

import os
import csv

from itertools import pairwise, permutations
from more_itertools import split_after, split_at


def split_by_dist(bgc_cluster: BGC_Cluster) -> List[BGC_Cluster]:
    gene_groups = split_after(pairwise(bgc_cluster.genes),
                              lambda p: p[1].coords.start - p[0].coords.end > config.max_dist_between_genes)
    return [BGC_Cluster(contig_id=bgc_cluster.contig_id,
                        genes=[gene for gene, _ in group])
            for group in gene_groups]


def split_by_single_gene_Starter_TE(bgc_cluster: BGC_Cluster) -> List[BGC_Cluster]:
    gene_groups = split_at(bgc_cluster.genes,
                           lambda gene: gene.modules[0].c_domain == C_Domain.C_STARTER and
                                        gene.modules[-1].terminal_domain,
                           keep_separator=True)
    return [BGC_Cluster(contig_id=bgc_cluster.contig_id,
                        genes=group)
            for group in gene_groups]


def genes_sequence_consistent(genes: List[Gene]) -> bool:
    joined_modules = [module for gene in genes for module in gene.modules]

    at_least_one_A_domain = any(module.a_domain == A_Domain.A for module in joined_modules)
    nc_terms_consistent = joined_modules[0].connecting_domain != ConnectingDomain.NTERM and \
                          joined_modules[-1].connecting_domain != ConnectingDomain.CTERM

    c_starter_consistent = all(module.c_domain != C_Domain.C_STARTER
                               for module in joined_modules[1:])
    te_td_consistent = all(module.terminal_domain == False
                           for module in joined_modules[:-1])

    return all([at_least_one_A_domain,
                nc_terms_consistent,
                c_starter_consistent,
                te_td_consistent])


def reverse_if_all_neg(genes: List[Gene]) -> List[Gene]:
    if all(gene.coords.stand == STRAND.REVERSE for gene in genes):
        return genes[::-1]
    else:
        return genes[:]


def genes_rearrangements(_genes: List[Gene]) -> List[List[Gene]]:
    genes = reverse_if_all_neg(_genes)
    starting_gene = [gene for gene in enumerate(genes)
                     if gene.modules[0].c_domain == C_Domain.C_STARTER]
    terminal_gene = [gene for gene in genes
                     if gene.modules[-1].terminal_domain]
    interior_genes = [gene for gene in genes
                      if gene not in starting_gene + terminal_gene]
    genes = starting_gene + interior_genes + terminal_gene
    if genes_sequence_consistent(genes):
        return [genes]
    else:   # TODO: optimize, keep only several top permutations
        return [starting_gene + permuted_interior_genes + terminal_gene
                for permuted_interior_genes in permutations(interior_genes)
                if genes_sequence_consistent(starting_gene + permuted_interior_genes + terminal_gene)]

def split_inconsistent(bgc: BGC_Cluster) -> List[BGC_Cluster]:
    genes = reverse_if_all_neg(bgc.genes)
    if genes_sequence_consistent(genes):
        return [BGC_Cluster(contig_id=bgc.contig_id, genes=genes)]

    genes_substrings = (genes[i:j]
                    for i in range(len(genes))
                    for j in range(i + 1, len(genes)))

    return [BGC_Cluster(contig_id=bgc.contig_id,
                        genes=rearranged_genes)
            for genes_substring in genes_substrings
            for rearranged_genes in genes_rearrangements(genes_substring)]
