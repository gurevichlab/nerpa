from typing import (
    Any,
    Dict,
    Iterable,
    NamedTuple,
    NewType,
    Tuple,
    Union
)
from src.generic.functional import compose
from src.nerpa_ms.monomer_graph.monomer import Monomer
from src.nerpa_ms.monomer_graph.monomer_graph import MonomerGraph, MonomerId
from src.rban_parsing.rban_parser import AtomId

from rdkit import Chem
from itertools import chain


AtomExtId = Tuple[MonomerId, AtomId]
AtomicEdgeExt = Tuple[AtomExtId, AtomExtId]

AtomAnyId = Any
AtomicEdge = Any
AtomInfoDict = dict
AtomicEdgeInfoDict = dict
MonomerEdgeInfoDict = dict

SMILES = NewType('SMILES', str)


class AtomicGraphRecord(NamedTuple):
    atoms: Dict[AtomAnyId, AtomInfoDict]  # AtomInfoDict should have 'name' attribute
    bonds: Dict[AtomicEdge, AtomicEdgeInfoDict]  # AtmoicEdgeInfoDict should have 'arity' attribute


class ChemWrap(NamedTuple):
    mol: Chem.rdchem.Mol
    atom_id_to_mol_idx: Dict[AtomExtId, int]


def get_monomer_ext_atoms(G: MonomerGraph,
                          monomer_idx: MonomerId) -> Dict[AtomExtId,
                                                          AtomInfoDict]:
    atomic_graph = G.monomer(monomer_idx).atomic_graph
    return {(monomer_idx, atom_idx): atomic_graph.nodes[atom_idx]
            for atom_idx in atomic_graph}


def get_monomer_ext_edges(G: MonomerGraph,
                          monomer_idx: MonomerId) -> Dict[AtomicEdgeExt,
                                                          AtomicEdgeInfoDict]:
    atomic_graph = G.monomer(monomer_idx).atomic_graph
    return {((monomer_idx, edge_ends[0]),
             (monomer_idx, edge_ends[1])): edge_info
            for edge_ends, edge_info in atomic_graph.edges.items()}


def monomer_edge_to_ext_atomic_edge(edge_info: MonomerEdgeInfoDict) -> Tuple[AtomicEdgeExt,
                                                                             AtomicEdgeInfoDict]:
    atomic_edge_ends_ext = tuple(edge_info['monomer_to_atom'].items())
    atomic_edge_info = {key: value
                        for key, value in edge_info.items()
                        if key != 'monomer_to_atom'}
    return (atomic_edge_ends_ext, atomic_edge_info)


# ext == extended. This can be confused with external, need to change
def to_atomic_graph(G: MonomerGraph) -> AtomicGraphRecord:
    def join_dicts(dicts: Iterable[dict]) -> dict:
        return dict(chain(*(d.items() for d in dicts)))

    atoms = join_dicts(get_monomer_ext_atoms(G, monomer_idx)
                       for monomer_idx in G.nodes)
    interior_edges = join_dicts(get_monomer_ext_edges(G, monomer_idx)
                                for monomer_idx in G.nodes)

    boundary_edges = dict(monomer_edge_to_ext_atomic_edge(edge_info)
                          for (u, v), edge_info in G.edges.items())
    return AtomicGraphRecord(atoms=atoms,
                             bonds=join_dicts([interior_edges, boundary_edges]))


def atomic_graph_to_chem(graph_record: AtomicGraphRecord) -> ChemWrap:
    # create empty editable mol object
    mol = Chem.RWMol()

    # add atoms to mol and keep track of index
    atom_id_to_index = {}
    for atom_id, atom_info in graph_record.atoms.items():
        atom_id_to_index[atom_id] = mol.AddAtom(Chem.Atom(atom_info['name']))

    # add bonds between adjacent atoms
    for (u, v), edge_attr in graph_record.bonds.items():
        # add relevant bond type (there are many more of these)
        if edge_attr['arity'] == 1:
            bond_type = Chem.rdchem.BondType.SINGLE
        elif edge_attr['arity'] == 2:
            bond_type = Chem.rdchem.BondType.DOUBLE
        else:
            raise ValueError('Unknown arity')
        mol.AddBond(atom_id_to_index[u], atom_id_to_index[v], bond_type)

    # Convert RWMol to Mol object
    mol = mol.GetMol()
    return ChemWrap(mol=mol,
                    atom_id_to_mol_idx=atom_id_to_index)


def graph_to_mol(G: MonomerGraph) -> Chem.rdchem.Mol:
    return compose(atomic_graph_to_chem,
                   to_atomic_graph,
                   G).mol


def graph_to_canon_smiles(G: MonomerGraph) -> SMILES:
    return compose(Chem.CanonSmiles,
                   Chem.MolToSmiles,
                   graph_to_mol,
                   G)


def isomorphic(G1: MonomerGraph,
               G2: MonomerGraph) -> bool:
    return graph_to_canon_smiles(G1) == graph_to_canon_smiles(G2)
