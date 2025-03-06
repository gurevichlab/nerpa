from __future__ import annotations
from typing import (
    Dict,
    List,
    Union,
    NamedTuple,
    NewType,
)
from dataclasses import dataclass
from enum import Enum, auto
from src.monomer_names_helper import (
    antiSMASH_MonomerName,
    AA10,
    AA34
)


antiSMASH_record = NewType('antiSMASH_record', dict)
GeneId = NewType('GeneId', str)


class SVM_LEVEL(Enum):
    SINGLE_AMINO = auto()
    SMALL_CLUSTER = auto()
    LARGE_CLUSTER = auto()
    PHYSIOCHEMICAL_CLASS = auto()


class SVM_Prediction(NamedTuple):
    score: float
    substrates: List[antiSMASH_MonomerName]


SVM_Predictions = Dict[SVM_LEVEL, SVM_Prediction]


@dataclass
class A_Domain:
    aa10: AA10
    aa34: AA34
    svm: SVM_Predictions


class DomainType(Enum):
    A = auto()
    PKS = auto()

    PCP = auto()

    C = auto()
    C_STARTER = auto()
    C_LCL = auto()
    C_DCL = auto()
    C_DUAL = auto()


    E = auto()
    MT = auto()

    TE_TD = auto()

    CTERM = auto()
    NTERM = auto()

    def in_c_domain_group(self) -> bool:
        return self in {DomainType.C, DomainType.C_STARTER, DomainType.C_LCL, DomainType.C_DCL, DomainType.C_DUAL}


@dataclass
class Module:
    a_domain: A_Domain = None
    domains_sequence: List[DomainType] = None

    def __post_init__(self):
        if self.domains_sequence is None:
            self.domains_sequence = []

class STRAND(Enum):
    FORWARD = auto()
    REVERSE = auto()
    

class Coords(NamedTuple):
    start: int
    end: int
    strand: STRAND


@dataclass
class Gene:
    gene_id: GeneId
    coords: Coords
    modules: List[Module]  # modules are in the order of appearance in the gene
    is_iterative: bool = False


@dataclass
class BGC_Cluster:
    genome_id: str
    contig_idx: int
    bgc_idx: int
    genes: List[Gene]

    def has_pks_domains(self) -> bool:
        return any(DomainType.PKS in module.domains_sequence
                   for gene in self.genes
                   for module in gene.modules)

    def has_a_domains(self) -> bool:
        return any(module.a_domain is not None
                   for gene in self.genes
                   for module in gene.modules)