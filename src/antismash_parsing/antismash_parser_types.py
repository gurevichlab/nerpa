from __future__ import annotations
from typing import (
    Dict,
    List,
    Union,
    NamedTuple,
    NewType,
    Optional,
)
from dataclasses import dataclass
from enum import Enum, auto
from src.monomer_names_helper import (
    antiSMASH_MonomerName,
    AA10,
    AA34
)
from src.data_types import BGC_ID


antiSMASH_record = NewType('antiSMASH_record', dict)
GeneId = NewType('GeneId', str)

class STRAND(Enum):
    FORWARD = auto()
    REVERSE = auto()


class Coords(NamedTuple):
    start: int
    end: int
    strand: STRAND

    # TODO: I am not sure if this is correct
    @classmethod
    def from_hmm_hit(cls, query_start: int, query_end: int) -> Coords:
        if query_start < query_end:
            return cls(query_start, query_end, STRAND.FORWARD)
        else:
            return cls(query_end, query_start, STRAND.REVERSE)

    def is_in(self, other: Coords) -> bool:
        return self.start >= other.start and self.end <= other.end

    def left_from(self, other: Coords) -> bool:
        return self.end < other.start

    def right_from(self, other: Coords) -> bool:
        return self.start > other.end


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


class BGC_Module_ID(NamedTuple):
    gene_id: GeneId
    module_idx: int


@dataclass
class Module:
    a_domain: Optional[A_Domain]
    domains_sequence: List[DomainType]


@dataclass
class Gene:
    gene_id: GeneId
    coords: Coords
    modules: List[Module]  # modules are in the order of appearance in the gene
    orphan_c_at_start: bool = False
    orphan_c_at_end: bool = False


@dataclass
class BGC_Cluster:
    bgc_id: BGC_ID
    genes: List[Gene]

    def has_pks_domains(self) -> bool:
        return any(DomainType.PKS in module.domains_sequence
                   for gene in self.genes
                   for module in gene.modules)

    def has_a_domains(self) -> bool:
        return any(module.a_domain is not None
                   for gene in self.genes
                   for module in gene.modules)


@dataclass
class Fragmented_BGC_Cluster:
    bgc_id: BGC_ID
    fragments: List[List[Gene]]

