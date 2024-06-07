from dataclasses import dataclass
from enum import Enum, auto
from typing import Dict, List, Union, NamedTuple



class STRAND(Enum):
    FORWARD = auto()
    REVERSE = auto()


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


class C_Domain(Enum):
    C = auto()
    C_STARTER = auto()
    C_LCL = auto()
    C_DCL = auto()
    C_DUAL = auto()


class ModifyingDomain(Enum):
    MT = auto()
    E = auto()


class ConnectingDomain(Enum):
    CTERM = auto()
    NTERM = auto()

'''
class DomainType(Enum):
    A = auto()
    C = auto()
    Connecting = auto()
    Modifying = auto()
    Terminal = auto()


Domain = Union[A_Domain, C_DOMAIN, ConnectingDomain, ModifyingDomain, TerminalDomain]
Module = Dict[DomainType, Union[Domain, List[Domain]]]
'''

@dataclass
class Module:
    a_domain: A_Domain = None
    c_domain: C_Domain = None
    connecting_domain: ConnectingDomain = None
    modifying_domains: List[ModifyingDomain] = None
    terminal_domain: bool = False

    def __post_init__(self):
        if self.modifying_domains is None:
            self.modifying_domains = []

class Coords(NamedTuple):
    start: int
    end: int
    stand: STRAND


@dataclass
class Gene:
    gene_id: str
    coords: Coords
    modules: List[Module]  # modules are in the order of appearance in the gene

    def modules_ctg_order(self):
        return self.modules[:] if self.coords.stand == STRAND.FORWARD else self.modules[::-1]


@dataclass
class BGC_Cluster:
    contig_id: str
    genes: List[Gene]