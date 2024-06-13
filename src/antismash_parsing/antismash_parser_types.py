from __future__ import annotations
from dataclasses import dataclass
from enum import Enum, auto
from typing import Dict, List, Union, NamedTuple
from src.data_types import MonomerResidue

class SVM_LEVEL(Enum):
    SINGLE_AMINO = auto()
    SMALL_CLUSTER = auto()
    LARGE_CLUSTER = auto()
    PHYSIOCHEMICAL_CLASS = auto()


class SVM_Prediction(NamedTuple):
    score: float
    monomer_residues: List[MonomerResidue]


SVM_Predictions = Dict[SVM_LEVEL, SVM_Prediction]


@dataclass
class A_Domain:
    aa10: str
    aa34: str
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
    gene_id: str
    coords: Coords
    modules: List[Module]  # modules are in the order of appearance in the gene


@dataclass
class BGC_Cluster:
    genome_id: str
    contig_id: str
    bgc_id: str
    genes: List[Gene]