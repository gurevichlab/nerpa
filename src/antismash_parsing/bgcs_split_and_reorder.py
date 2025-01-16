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
    STRAND
)
from src.config import antiSMASH_Parsing_Config
from src.generic.combinatorics import generate_permutations, split_sequence_blocks
from functools import partial
from itertools import chain, islice, pairwise, product, groupby
from more_itertools import split_before, split_at
from src.pipeline.logger import NerpaLogger

Fragmented_BGC_Cluster = List[BGC_Cluster]


def split_by_dist(bgc_cluster: BGC_Cluster,
                  config: antiSMASH_Parsing_Config) -> List[BGC_Cluster]:
    dummy_gene = Gene(GeneId(''), Coords(start=-1, end=-1, strand=STRAND.FORWARD), [])  # dummy gene to start the pairwise
    gene_groups = list(split_before(pairwise([dummy_gene] + bgc_cluster.genes),
                                    lambda p: p[1].coords.start - p[0].coords.end > config.MAX_DISTANCE_BETWEEN_GENES))
    return [BGC_Cluster(genome_id=bgc_cluster.genome_id,
                        contig_idx=bgc_cluster.contig_idx,
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


def reverse_if_all_neg(genes: List[Gene]) -> List[Gene]:
    if all(gene.coords.strand == STRAND.REVERSE for gene in genes):
        reversed_genes = genes[::-1]
        for gene in reversed_genes:
            gene.coords = gene.coords._replace(strand=STRAND.FORWARD)
        return reversed_genes
    else:
        return genes[:]


def get_genes_rearrangements(genes_: Iterable[Gene], config: antiSMASH_Parsing_Config) -> List[List[Gene]]:
    genes = list(genes_)
    if genes_sequence_consistent(genes):
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


BGC_Fragment = List[Gene]  # an assembly line (usually starts with C_STARTER or A-PCP and ends with TE)
Gene_Block = List[BGC_Fragment]  # sequence of BGC fragments with the same orientation

def get_bgc_fragments(block_: Iterable[Gene]) -> List[List[Gene]]:
    # it is assumed that all genes are on the same strand
    block = list(block_)
    if block[0].coords.strand == STRAND.REVERSE:
        block = block[::-1]

    # sometimes genes work in reverse order. This is a stub for the case of just 2 genes
    # sometimes even making fragments inconsistent (see BGC0000437)
    if all([len(block) == 2,
            not genes_sequence_consistent(block),
            genes_sequence_consistent(block[::-1])]):
        return [block[::-1]]

    fragments = []
    fragment = []
    for gene in block:
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


# TODO: handle more complex rearrangements
def get_block_fragments_rearrangements(block: List[BGC_Fragment]) -> Iterable[List[BGC_Fragment]]:
    if len(block) == 1 or block[0][0].coords.strand == STRAND.FORWARD:
        return [block]
    else:
        return [block, block[::-1]]


def generate_fragmented_bgcs(bgc: BGC_Cluster, config: antiSMASH_Parsing_Config) -> Iterable[Fragmented_BGC_Cluster]:
    def get_rearranged_fragments(genes_fragments: List[List[Gene]]) -> Iterable[Tuple[List[Gene]]]:
        return islice(product(*(get_genes_rearrangements(genes_fragment, config=config)
                                for genes_fragment in genes_fragments)),
                      config.MAX_PERMUTATIONS_PER_BGC)

    def build_fragmented_bgc(genes_fragments: Iterable[List[Gene]]) -> Fragmented_BGC_Cluster:
        return [BGC_Cluster(genome_id=bgc.genome_id,
                            contig_idx=bgc.contig_idx,
                            bgc_idx=bgc.bgc_idx,
                            genes=genes_fragment)
                for genes_fragment in genes_fragments
                if any(module.a_domain is not None
                       for gene in genes_fragment
                       for module in gene.modules)]

    gene_blocks: List[Gene_Block] = [get_bgc_fragments(genes_group)
                                     for _, genes_group in groupby(bgc.genes, key=lambda gene: gene.coords.strand)]  # split genes by strand
    block_rearrangements_ = product(*(get_block_fragments_rearrangements(block_fragments)
                                      for block_fragments in gene_blocks))

    if len(gene_blocks) == 2:  # a stub for the cases like BGC0000289 TODO: handle more general cases
        block_rearrangements = chain(block_rearrangements_,
                                     product(*(get_block_fragments_rearrangements(block_fragments)
                                               for block_fragments in gene_blocks[::-1])))
    else:
        block_rearrangements = block_rearrangements_

    fragments_sequences = [[fragment for block in rearranged_blocks for fragment in block]
                            for rearranged_blocks in block_rearrangements]

    return (build_fragmented_bgc(rearranged_fragments)
            for fragments_sequence in fragments_sequences
            for rearranged_fragments in get_rearranged_fragments(fragments_sequence))


def split_and_reorder(bgc_: BGC_Cluster,
                      config: antiSMASH_Parsing_Config,
                      log: NerpaLogger) -> List[Fragmented_BGC_Cluster]:
    result = list(islice((fragmented_bgc
                          for bgc in [bgc_]  # split_by_dist(bgc_, config=config)
                          for fragmented_bgc in generate_fragmented_bgcs(bgc, config=config)),
                         config.MAX_VARIANTS_PER_BGC + 1))
    if len(result) > config.MAX_VARIANTS_PER_BGC:
        log.info(f'WARNING: Too many BGC variants. Keeping first {config.MAX_VARIANTS_PER_BGC} of them.')
        # del result[-1]  # what is the purpose of this line?
        result = result[:config.MAX_VARIANTS_PER_BGC]
    return result
