from typing import Dict, List, Tuple
from enum import Enum, auto
from functools import total_ordering
from dataclasses import dataclass
from collections import defaultdict

from src.rban_parsing.rban_parser import (
    AtomId,
    AtomInfo,
    AtomicEdge,
    AtomicEdgeInfo
)

@total_ordering
class BindingSiteType(Enum):
    AMINO_C = auto()
    AMINO_N = auto()

    ESTER_C = auto()
    ESTER_O = auto()

    UNKNOWN = auto()

    def __lt__(self, other):
        return self.value < other.value

    def atom(self):
        match self:
            case BindingSiteType.AMINO_C: return 'C'
            case BindingSiteType.AMINO_N: return 'N'
            case BindingSiteType.ESTER_C: return 'C'
            case BindingSiteType.ESTER_O: return 'O'


# Basically Dict[BindingSiteType, int], just hashable
# TODO: dict may be an overkill and inefficient
class BindingSitesFgpt(dict):
    def __init__(self, bs_fgpt: Dict[BindingSiteType, int]):
        super().__init__(bs_fgpt)

    def __key(self):
        return tuple((k, self[k]) for k in sorted(self))

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        return self.__key() == other.__key()


AMINO2_FGPT = BindingSitesFgpt({BindingSiteType.AMINO_C: 1,
                                BindingSiteType.AMINO_N: 1})


@dataclass
class BondType:
    binding_sites: Tuple[BindingSiteType, BindingSiteType]
    name: str


AMINO_BOND = BondType(binding_sites=(BindingSiteType.AMINO_C,
                                     BindingSiteType.AMINO_N),
                      name='AMINO')

ESTER_BOND = BondType(binding_sites=(BindingSiteType.ESTER_C,
                                     BindingSiteType.ESTER_O),
                      name='ESTER')

BONDS = (AMINO_BOND, ESTER_BOND)


def get_binding_sites(monomer_atoms_idxs: List[AtomId],
                      boundary_edges: Dict[AtomicEdge, AtomicEdgeInfo],
                      atoms: Dict[AtomId, AtomInfo]) -> Dict[BindingSiteType, List[AtomId]]:
    binding_sites = defaultdict(list)
    for (u, v), bond_info in boundary_edges.items():
        if u not in monomer_atoms_idxs:
            u, v = v, u  # u is inside monomer, v is outside

        bond_type = next((bond_type for bond_type in BONDS
                          if bond_type.name == bond_info.bondType), None)
        if bond_type is not None:
            site_type = bond_type.binding_sites[0] if atoms[u].name == bond_type.binding_sites[0].atom() else \
                bond_type.binding_sites[1]
        else:
            site_type = BindingSiteType.UNKNOWN

        binding_sites[site_type].append(u)

    return binding_sites
