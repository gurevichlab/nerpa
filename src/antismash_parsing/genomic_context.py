from __future__ import annotations
from typing import (
    List,
    NewType,
    Optional,
    Tuple,
    Union, Dict
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
    Gene, GeneId, BGC_Module_ID
)
from src.monomer_names_helper import enum_representer
import yaml

class ModuleGenomicContextFeature(Enum):
    START_OF_BGC = auto()
    START_OF_GENE = auto()
    END_OF_GENE = auto()
    ASSEMBLY_LINE_START_MARKER = auto()
    ASSEMBLY_LINE_FINISH_MARKER = auto()
    PKS_UPSTREAM = auto()
    PKS_DOWNSTREAM = auto()
    ONLY_A_DOMAIN = auto()  # the module has only A domain, no other domains

    def __lt__(self, other):
        if not isinstance(other, ModuleGenomicContextFeature):
            return NotImplemented
        return self.name < other.name

yaml.add_representer(ModuleGenomicContextFeature, enum_representer)

ModuleGenomicContext = Tuple[ModuleGenomicContextFeature, ...]


def get_modules_genomic_context(fragments: List[List[Gene]]) -> Dict[BGC_Module_ID, ModuleGenomicContext]:
    MGCF = ModuleGenomicContextFeature  # for brevity

    module_id_to_genomic_context = {}
    for fragment_idx, fragment_genes in enumerate(fragments):
        for gene_idx, gene in enumerate(fragment_genes):
            for module_idx, module in enumerate(gene.modules):
                if module.a_domain is None:
                    continue

                fst_module_with_a_domain_in_gene = not any(DomainType.A in module.domains_sequence
                                                   for module in gene.modules[:module_idx])
                last_module_with_a_domain_in_gene = not any(DomainType.A in module.domains_sequence
                                                    for module in gene.modules[module_idx + 1:])
                fst_module_with_a_domain_in_fragment = (
                    fst_module_with_a_domain_in_gene and
                    not any(
                        DomainType.A in prev_module.domains_sequence
                        for prev_gene in fragment_genes[:gene_idx]
                        for prev_module in prev_gene.modules
                    )
                )
                last_module_with_a_domain_in_fragment = (
                    last_module_with_a_domain_in_gene and
                    not any(
                        DomainType.A in next_module.domains_sequence
                        for next_gene in fragment_genes[gene_idx + 1:]
                        for next_module in next_gene.modules)
                )

                upstream_modules = chain(chain(*(prev_gene.modules for prev_gene in fragment_genes[:gene_idx])),
                                         gene.modules[:module_idx])
                upstream_domains = chain(*(prev_module.domains_sequence for prev_module in upstream_modules))
                upstream_domains_before_A = takewhile(lambda domain: domain != DomainType.A,
                                                      reversed(list(upstream_domains)))

                downstream_modules = chain(gene.modules[module_idx + 1:],
                                           *(gene.modules for gene in fragment_genes[gene_idx + 1:]))
                downstream_domains = chain(*(module.domains_sequence for module in downstream_modules))
                downstream_domains_before_A = takewhile(lambda domain: domain != DomainType.A,
                                                        downstream_domains)

                pks_upstream = DomainType.PKS in upstream_domains_before_A
                pks_downstream = DomainType.PKS in downstream_domains_before_A

                features = (
                    feature
                    for feature, is_present in [
                    (MGCF.START_OF_BGC, fragment_idx == 0 and gene_idx == 0 and fst_module_with_a_domain_in_gene),
                    (MGCF.START_OF_GENE, fst_module_with_a_domain_in_gene),
                    (MGCF.END_OF_GENE, last_module_with_a_domain_in_gene),
                    (MGCF.ASSEMBLY_LINE_START_MARKER, fst_module_with_a_domain_in_fragment),
                    (MGCF.ASSEMBLY_LINE_FINISH_MARKER, last_module_with_a_domain_in_fragment),
                    (MGCF.PKS_UPSTREAM, pks_upstream),
                    (MGCF.PKS_DOWNSTREAM, pks_downstream),
                    (MGCF.ONLY_A_DOMAIN, module.domains_sequence == [DomainType.A])
                ]
                    if is_present
                )
                module_id_to_genomic_context[BGC_Module_ID(gene.gene_id, module_idx)] = tuple(sorted(features, key=lambda x: x.name))

    return module_id_to_genomic_context