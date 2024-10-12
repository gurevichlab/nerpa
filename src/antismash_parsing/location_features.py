from typing import (
    List,
    NewType,
    Tuple
)
from enum import Enum, auto
from itertools import (
    chain,
    compress,
    takewhile
)
from src.antismash_parsing.antismash_parser_types import (
    BGC_Cluster,
    DomainType,
    Gene
)

class ModuleLocFeature(Enum):
    START_OF_BGC = auto()
    END_OF_BGC = auto()
    START_OF_FRAGMENT = auto()
    END_OF_FRAGMENT = auto()
    START_OF_GENE = auto()
    END_OF_GENE = auto()
    PKS_UPSTREAM_PREV_GENE = auto()
    PKS_UPSTREAM_SAME_GENE = auto()
    PKS_DOWNSTREAM_NEXT_GENE = auto()
    PKS_DOWNSTREAM_SAME_GENE = auto()


class GeneLocFeature(Enum):
    START_OF_BGC = auto()
    END_OF_BGC = auto()
    START_OF_FRAGMENT = auto()
    END_OF_FRAGMENT = auto()
    PKS_UPSTREAM = auto()
    PKS_DOWNSTREAM = auto()


class BGC_Fragment_Loc_Feature(Enum):
    START_OF_BGC = auto()
    END_OF_BGC = auto()
    PKS_UPSTREAM = auto()
    PKS_DOWNSTREAM = auto()

ModuleLocFeatures = Tuple[ModuleLocFeature, ...]
GeneLocFeatures = Tuple[GeneLocFeature, ...]
BGC_Fragment_Loc_Features = Tuple[BGC_Fragment_Loc_Feature, ...]


def find_pks_upstream(gene: Gene, module_idx: int) -> bool:
    domains_upstream = chain(*(module.domains_sequence[::-1]
                               for module in reversed(gene.modules[:module_idx])))
    return DomainType.PKS in takewhile(lambda domain: domain != domain.A, domains_upstream)


def find_pks_downstream(gene: Gene, module_idx: int) -> bool:
    domains_downstream = chain(*(module.domains_sequence
                                 for module in gene.modules[module_idx + 1:]))
    return DomainType.PKS in takewhile(lambda domain: domain != domain.A, domains_downstream)


def get_bgc_fragment_loc_features(fragment_idx: int,
                                  fragments: List[BGC_Cluster]) -> BGC_Fragment_Loc_Features:
    pairs = [
        (BGC_Fragment_Loc_Feature.START_OF_BGC, fragment_idx == 0),
        (BGC_Fragment_Loc_Feature.END_OF_BGC, fragment_idx == len(fragments) - 1),
        (BGC_Fragment_Loc_Feature.PKS_UPSTREAM,
         fragment_idx > 0 and find_pks_upstream(fragments[fragment_idx - 1].genes[-1], len(fragments[fragment_idx - 1].genes[-1].modules))),
        (BGC_Fragment_Loc_Feature.PKS_DOWNSTREAM,
         fragment_idx < len(fragments) - 1 and find_pks_downstream(fragments[fragment_idx + 1].genes[0], -1))
    ]
    return tuple(compress(*zip(*pairs)))


def get_gene_loc_features(gene_idx: int,
                          fragment: List[Gene],
                          fragment_features: BGC_Fragment_Loc_Features) -> GeneLocFeatures:
    fst_gene_with_a_domains = not any(DomainType.A in module.domains_sequence
                                      for gene in fragment[:gene_idx]
                                      for module in gene.modules)
    last_gene_with_a_domains = not any(DomainType.A in module.domains_sequence
                                       for gene in fragment[gene_idx + 1:]
                                       for module in gene.modules)
    pairs = [
        (GeneLocFeature.START_OF_BGC,
         BGC_Fragment_Loc_Feature.START_OF_BGC in fragment_features and fst_gene_with_a_domains),
        (GeneLocFeature.END_OF_BGC,
         BGC_Fragment_Loc_Feature.END_OF_BGC in fragment_features and last_gene_with_a_domains),
        (GeneLocFeature.START_OF_FRAGMENT, fst_gene_with_a_domains),
        (GeneLocFeature.END_OF_FRAGMENT, last_gene_with_a_domains),
        (GeneLocFeature.PKS_UPSTREAM, any([gene_idx == 0 and BGC_Fragment_Loc_Feature.PKS_UPSTREAM in fragment_features,
                                           gene_idx > 0 and find_pks_upstream(fragment[gene_idx - 1], len(fragment[gene_idx - 1].modules))])),
        (GeneLocFeature.PKS_DOWNSTREAM, any([gene_idx == len(fragment) - 1 and BGC_Fragment_Loc_Feature.PKS_DOWNSTREAM in fragment_features,
                                             gene_idx < len(fragment) - 1 and find_pks_downstream(fragment[gene_idx + 1], -1)]))
    ]
    return tuple(compress(*zip(*pairs)))


def get_module_loc_features(module_idx: int,
                            gene: Gene,
                            gene_features: GeneLocFeatures) -> ModuleLocFeatures:
    fst_module_with_a_domain = not any(DomainType.A in module.domains_sequence
                                       for module in gene.modules[:module_idx])
    last_module_with_a_domain = not any(DomainType.A in module.domains_sequence
                                        for module in gene.modules[module_idx + 1:])
    pairs = [
        (ModuleLocFeature.START_OF_BGC,
         GeneLocFeature.START_OF_BGC in gene_features and fst_module_with_a_domain),
        (ModuleLocFeature.END_OF_BGC,
         GeneLocFeature.END_OF_BGC in gene_features and last_module_with_a_domain),
        (ModuleLocFeature.START_OF_FRAGMENT,
         GeneLocFeature.START_OF_FRAGMENT in gene_features and fst_module_with_a_domain),
        (ModuleLocFeature.END_OF_FRAGMENT,
         GeneLocFeature.END_OF_FRAGMENT in gene_features and last_module_with_a_domain),
        (ModuleLocFeature.START_OF_GENE, fst_module_with_a_domain),
        (ModuleLocFeature.END_OF_GENE, last_module_with_a_domain),
        (ModuleLocFeature.PKS_UPSTREAM_PREV_GENE,
         GeneLocFeature.PKS_UPSTREAM in gene_features and fst_module_with_a_domain),
        (ModuleLocFeature.PKS_UPSTREAM_SAME_GENE,
         find_pks_upstream(gene, module_idx)),
        (ModuleLocFeature.PKS_DOWNSTREAM_NEXT_GENE,
         GeneLocFeature.PKS_DOWNSTREAM in gene_features and last_module_with_a_domain),
        (ModuleLocFeature.PKS_DOWNSTREAM_SAME_GENE,
         find_pks_downstream(gene, module_idx))
    ]
    return tuple(compress(*zip(*pairs)))