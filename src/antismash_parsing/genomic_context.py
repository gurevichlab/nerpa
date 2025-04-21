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
    START_OF_GENE = auto()
    END_OF_GENE = auto()
    PKS_UPSTREAM = auto()
    PKS_DOWNSTREAM = auto()
    ONLY_A_DOMAIN = auto()  # the module has only A domain, no other domains

yaml.add_representer(ModuleGenomicContextFeature, enum_representer)

ModuleGenomicContext = Tuple[ModuleGenomicContextFeature, ...]


def get_modules_genomic_context(genes: List[Gene]) -> Dict[BGC_Module_ID, ModuleGenomicContext]:
    module_id_to_genomic_context = {}
    for gene_idx, gene in enumerate(genes):
        for module_idx, module in enumerate(gene.modules):
            if module.a_domain is None:
                continue

            fst_module_with_a_domain = not any(DomainType.A in module.domains_sequence
                                               for module in gene.modules[:module_idx])
            last_module_with_a_domain = not any(DomainType.A in module.domains_sequence
                                                for module in gene.modules[module_idx + 1:])

            if gene_idx > 0:
                upstream_modules = genes[gene_idx - 1].modules
                pks_upstream = (
                        any(DomainType.PKS in m.domains_sequence for m in upstream_modules)
                        and not any(DomainType.A in m.domains_sequence for m in upstream_modules)
                )
            else:
                pks_upstream = False

            downstream_modules_same_gene = takewhile(lambda module: module.a_domain is not None,
                                                    gene.modules[module_idx + 1:])
            if gene_idx < len(genes) - 1:
                downstream_modules_next_gene = takewhile(lambda module: module.a_domain is not None,
                                                        genes[gene_idx + 1].modules)
            else:
                downstream_modules_next_gene = []
            pks_downstream = any(DomainType.PKS in m.domains_sequence
                                 for m in chain(downstream_modules_same_gene,
                                                downstream_modules_next_gene))

            features = (feature
                        for feature, is_present in [
                            (ModuleGenomicContextFeature.START_OF_GENE, fst_module_with_a_domain),
                            (ModuleGenomicContextFeature.END_OF_GENE, last_module_with_a_domain),
                            (ModuleGenomicContextFeature.PKS_UPSTREAM, pks_upstream),
                            (ModuleGenomicContextFeature.PKS_DOWNSTREAM, pks_downstream),
                            (ModuleGenomicContextFeature.ONLY_A_DOMAIN, module.domains_sequence == [DomainType.A])
                        ] if is_present)
            module_id_to_genomic_context[BGC_Module_ID(gene.gene_id, module_idx)] = tuple(sorted(features, key=lambda x: x.name))

    return module_id_to_genomic_context