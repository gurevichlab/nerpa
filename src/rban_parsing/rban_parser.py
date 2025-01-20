from typing import (
    Dict,
    List,
    NamedTuple,
    Tuple,
    Union
)
from src.rban_parsing import handle_monomers
from src.data_types import (
    Chirality,
    SMILES
)
from src.monomer_names_helper import NorineMonomerName
from networkx import DiGraph, is_isomorphic

from collections import defaultdict
from dataclasses import dataclass

Raw_rBAN_Record = dict

AtomId = int
AtomicEdge = Tuple[AtomId, AtomId]


class AtomInfo(NamedTuple):
    name: str
    hydrogens: int


class AtomicEdgeInfo(NamedTuple):
    arity: int
    bondType: str


MonomerIdx = int
MonomerEdge = Tuple[MonomerIdx, MonomerIdx]


class MonomerInfo(NamedTuple):
    name: NorineMonomerName
    atoms: List[AtomId]
    chirality: Chirality
    is_pks_hybrid: bool = False


class MonomerEdgeInfo(NamedTuple):
    monomer_to_atom: Dict[MonomerIdx, AtomId]
    arity: int  # I heard that there exist fractional arities (e.g. 1.5)
    bondType: str


def parsed_chiralities(chiralities: Dict[MonomerIdx, Union[bool, None]]) -> Dict[MonomerIdx, Chirality]:
    def parsed_chirality(ch: Union[bool, None]):
        if ch is None:
            return Chirality.UNKNOWN
        else:
            return Chirality.D if ch else Chirality.L
    return defaultdict(lambda: Chirality.UNKNOWN,
                       {mon_id: parsed_chirality(ch)
                        for mon_id, ch in chiralities.items()})


@dataclass
class Parsed_rBAN_Record:
    compound_id: str
    monomers: Dict[MonomerIdx, MonomerInfo]
    monomer_bonds: Dict[MonomerEdge, MonomerEdgeInfo]
    atoms: Dict[AtomId, AtomInfo]
    atomic_bonds: Dict[AtomicEdge, AtomicEdgeInfo]

    def __init__(self, rban_record: Raw_rBAN_Record,
                 hybrid_monomers: Dict[MonomerIdx, NorineMonomerName],
                 chiralities: Dict[MonomerIdx, Chirality]):

        def get_monomer_name(idx: MonomerIdx) -> NorineMonomerName:
            if idx in hybrid_monomers:
                name = hybrid_monomers[idx]
            else:
                monomers = rban_record['monomericGraph']['monomericGraph']['monomers']
                name = next(monomer['monomer']['monomer']['monomer']
                            for monomer in monomers
                            if monomer['monomer']['index'] == idx)
            name.replace('C10:0-NH2(2)-Ep(9)-oxo(8)', 'Aeo')  # I have no idea what this is, I just copied
            return NorineMonomerName(name)

        self.compound_id = rban_record['id']
        self.atoms = {atom['cdk_idx']: AtomInfo(name=atom['name'],
                                                hydrogens=atom['hydrogens'])
                      for atom in rban_record['atomicGraph']['atomicGraph']['atoms']}

        self.monomers = {MonomerIdx(idx := monomer['monomer']['index']):
                             MonomerInfo(name=get_monomer_name(idx),
                                         atoms=[AtomId(atom) for atom in monomer['monomer']['atoms']],
                                         chirality=chiralities[idx],
                                         is_pks_hybrid=idx in hybrid_monomers)
                         for monomer in rban_record['monomericGraph']['monomericGraph']['monomers']}

        self.atomic_bonds = {(AtomId(bond['atoms'][0]), AtomId(bond['atoms'][1])):
                                 AtomicEdgeInfo(arity=bond['arity'],
                                                bondType=bond['bondType'])
                             for bond in rban_record['atomicGraph']['atomicGraph']['bonds']}

        self.monomer_bonds = {}
        rban_atomic_bonds = rban_record['atomicGraph']['atomicGraph']['bonds']
        for bond in rban_record['monomericGraph']['monomericGraph']['bonds']:
            bond = bond['bond']
            mon1, mon2 = bond['monomers']
            atomic_bond_idx = bond['atomicIndexes'][0]
            atom1, atom2 = rban_atomic_bonds[atomic_bond_idx]['atoms']
            monomer_to_atom = {mon1: atom1, mon2: atom2} if atom1 in self.monomers[mon1].atoms \
                else {mon1: atom2, mon2: atom1}
            self.monomer_bonds[(mon1, mon2)] = MonomerEdgeInfo(monomer_to_atom=monomer_to_atom,
                                                               arity=self.atomic_bonds[(atom1, atom2)].arity,
                                                               bondType=self.atomic_bonds[(atom1, atom2)].bondType)

    def to_compact_dict(self) -> dict:
        return {'compound_id': self.compound_id,
                'nodes': {mon_idx: mon_info.name
                          for mon_idx, mon_info in self.monomers.items()},
                'edges': [(u, v, edge_info.bondType)
                          for (u, v), edge_info in self.monomer_bonds.items()]}

    def to_nx_monomer_graph(self) -> DiGraph:
        G = DiGraph()
        for mon_idx, mon_info in self.monomers.items():
            G.add_node(mon_idx, name=f'{mon_info.name}_{mon_info.chirality.name}_{mon_info.is_pks_hybrid}')
        for (u, v), edge_info in self.monomer_bonds.items():
            G.add_edge(u, v, bond_type=edge_info.bondType)
        return G

    @classmethod
    def graphs_isomorphic(cls, graph1: DiGraph, graph2: DiGraph) -> bool:
        def node_match(n1, n2):
            name1, name2 = n1.get('name'), n2.get('name')
            return any([name1 == name2,  # same amino acid
                        name1.startswith('X') and name2.startswith('X'),  # both unknown
                        ':' in name1 and ':'])  # both lipid tails

        def edge_match(e1, e2):
            return e1.get('bond_type') == e2.get('bond_type')

        return is_isomorphic(graph1, graph2, node_match=node_match, edge_match=edge_match)


def get_hybrid_monomers_smiles(rban_record: Raw_rBAN_Record) -> List[Tuple[MonomerIdx, SMILES]]:
    hybrid_monomers = []
    for monomer in rban_record["monomericGraph"]["monomericGraph"]['monomers']:
        if monomer['monomer']['monomer']['monomer'].startswith('X'):  # monomer was not recognized
            smi = monomer['monomer']['monomer']['smiles']
            monomer_id = monomer['monomer']['index']
            try:
                aa_smi, pk_smi, _ = handle_monomers.split_aa_pk_hybrid(smi)
                hybrid_monomers.append((monomer_id, aa_smi))
            except handle_monomers.PKError:
                pass  # it's okay

    return hybrid_monomers
