from typing import (
    Iterable,
    List,
    Tuple
)
from src.antismash_parsing.antismash_parser_types import (
    BGC_Cluster,
    Coords,
    Gene,
    GeneId,
    Module,
    DomainType,
    STRAND, Fragmented_BGC_Cluster
)
from src.config import antiSMASH_Processing_Config
from src.generic.combinatorics import generate_permutations, split_sequence_blocks
from functools import partial
from itertools import chain, islice, pairwise, product, groupby
from more_itertools import split_before, split_at
from src.pipeline.logger import NerpaLogger


def split_by_dist(bgc_cluster: BGC_Cluster,
                  config: antiSMASH_Processing_Config) -> List[BGC_Cluster]:
    dummy_gene = Gene(GeneId(''), Coords(start=-1, end=-1, strand=STRAND.FORWARD), [])  # dummy gene to start the pairwise
    gene_groups = list(split_before(pairwise([dummy_gene] + bgc_cluster.genes),
                                    lambda p: p[1].coords.start - p[0].coords.end > config.MAX_DISTANCE_BETWEEN_GENES))
    return [BGC_Cluster(genome_id=bgc_cluster.genome_id,
                        contig_idx=bgc_cluster.contig_idx,
                        bgc_idx=bgc_cluster.bgc_idx,
                        genes=[gene for _, gene in group])
            for group in gene_groups]


def split_by_single_gene_Starter_TE(bgc_cluster: BGC_Cluster,
                                    config: antiSMASH_Processing_Config) -> List[BGC_Cluster]:  # config is not used but added for consistency
    gene_groups = split_at(bgc_cluster.genes,
                           lambda gene: DomainType.C_STARTER in gene.modules[0].domains_sequence and
                                        DomainType.TE_TD in gene.modules[-1].domains_sequence,
                           keep_separator=True)
    gene_groups = [list(group) for group in gene_groups]
    return [BGC_Cluster(genome_id=bgc_cluster.genome_id,
                        contig_idx=bgc_cluster.contig_idx,
                        bgc_idx=bgc_cluster.bgc_idx,
                        genes=group)
            for group in gene_groups
            if len(group) > 0]


def a_pcp_module(module: Module) -> bool:
    domains_set = set(module.domains_sequence)
    return {DomainType.A, DomainType.PCP}.issubset(domains_set) and all(not DomainType.in_c_domain_group(domain)
                                                                        for domain in domains_set)


def genes_contents_consistent(genes: List[Gene]) -> bool:
    if len(genes) == 1:
        return True

    joined_modules = [module for gene in genes for module in gene.modules]

    at_least_one_A_domain = any(DomainType.A in module.domains_sequence for module in joined_modules)
    c_starter_consistent = sum((DomainType.C_STARTER in module.domains_sequence) for module in joined_modules) <= 1
    a_pcp_consistent = sum(a_pcp_module(module) for module in joined_modules) <= 1
    te_td_consistent = sum((DomainType.TE_TD in module.domains_sequence) for module in joined_modules) <= 1

    return all([at_least_one_A_domain,
                a_pcp_consistent,
                c_starter_consistent,
                te_td_consistent])


def genes_sequence_consistent(genes: List[Gene]) -> bool:
    if len(genes) == 1:
        return True

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


# assembly line usually starts with C_STARTER or A-PCP and ends with TE
# however, it can contain C_STARTER/A-PCP/TE (aka border markers) in the middle
# BGC_Fragment is a sequence of genes with the same orientation with no border markers in the middle (but may have them at the ends)
# Thus the actual assembly line is a concatination of several (maybe one) BGC_Fragment
BGC_Fragment = List[Gene]

def get_fragment_rearrangements(genes: BGC_Fragment,
                                config: antiSMASH_Processing_Config) -> Iterable[BGC_Fragment]:
    if genes_sequence_consistent(genes):  # no need to permute
        return [genes]

    # this function is called only after consistency check
    # so there is at most one starting gene and one terminal gene
    starting_gene = next(([gene] for gene in genes
                          if DomainType.C_STARTER in gene.modules[0].domains_sequence or a_pcp_module(gene.modules[0])),
                         [])
    terminal_gene = next(([gene] for gene in genes
                          if DomainType.TE_TD in gene.modules[-1].domains_sequence),
                         [])
    interior_genes = [gene for gene in genes
                      if gene not in starting_gene + terminal_gene]
    result = [starting_gene + permuted_interior_genes + terminal_gene
              for permuted_interior_genes in islice(generate_permutations(interior_genes),
                                                    config.MAX_PERMUTATIONS_PER_BGC)
              if genes_sequence_consistent(starting_gene + permuted_interior_genes + terminal_gene)]
    return result if result else [genes]



def split_genes_into_fragments(genes: List[Gene]) -> List[BGC_Fragment]:
    # it is assumed that all genes are on the same strand

    fragments = []
    fragment: BGC_Fragment = []
    for gene in genes:
        is_start_gene = DomainType.C_STARTER in gene.modules[0].domains_sequence or a_pcp_module(gene.modules[0])
        is_end_gene = DomainType.TE_TD in gene.modules[-1].domains_sequence

        if is_start_gene:
            if fragment:
                fragments.append(fragment)
            fragment = []

        fragment.append(gene)

        if is_end_gene:
            fragments.append(fragment)
            fragment = []

    if fragment:
        fragments.append(fragment)

    return fragments


GeneBlock = List[Gene]  # a continuous sequence of genes from the same strand
def split_gene_block_into_fragments(genes: List[Gene]) -> Iterable[List[BGC_Fragment]]:
    """returns all the ways to split a gene block into fragments"""
    splits = [split_genes_into_fragments(genes)]
    if len(genes) == 2 or genes[0].coords.strand == STRAND.REVERSE:
        splits.append(split_genes_into_fragments(genes[::-1]))
    return splits


def generate_fragmented_bgcs(bgc: BGC_Cluster) -> Iterable[Fragmented_BGC_Cluster]:
    def gene_blocks_to_fragment_sequences(gene_blocks: List[GeneBlock]) -> Iterable[Tuple[List[BGC_Fragment], ...]]:
        return product(*(split_gene_block_into_fragments(gene_block)
                         for gene_block in gene_blocks))

    def build_fragmented_bgc(genes_fragments: Iterable[List[Gene]]) -> Fragmented_BGC_Cluster:
        return Fragmented_BGC_Cluster(bgc_id=bgc.bgc_id,
                                      fragments=[genes_fragment
                                                 for genes_fragment in genes_fragments
                                                 if any(module.a_domain is not None
                                                        for gene in genes_fragment
                                                        for module in gene.modules)])

    # split genes by strand
    gene_blocks: List[GeneBlock] = [list(genes_group)
                                    for _, genes_group in groupby(bgc.genes,
                                                                  key=lambda gene: gene.coords.strand)]
    gene_blocks_rearrangements = [gene_blocks, gene_blocks[::-1]] \
        if len(gene_blocks) == 2 else [gene_blocks]

    for gene_blocks in gene_blocks_rearrangements:
        for fragments_seqs in gene_blocks_to_fragment_sequences(gene_blocks):
            joined_fragments = chain(*fragments_seqs)
            yield build_fragmented_bgc(joined_fragments)


def split_and_reorder(bgc: BGC_Cluster,
                      config: antiSMASH_Processing_Config,
                      log: NerpaLogger) -> List[Fragmented_BGC_Cluster]:
    result = list(islice((fragmented_bgc
                          for fragmented_bgc in generate_fragmented_bgcs(bgc)),
                         config.MAX_VARIANTS_PER_BGC + 1))
    if len(result) > config.MAX_VARIANTS_PER_BGC:
        log.info(f'WARNING: Too many BGC variants. Keeping first {config.MAX_VARIANTS_PER_BGC} of them.')
        # del result[-1]
        result = result[:config.MAX_VARIANTS_PER_BGC]
    return result
