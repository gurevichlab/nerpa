from __future__ import annotations
from typing import (
    Dict,
    List,
    NamedTuple,
    Tuple,
    Union, Optional
)
from src.rban_parsing import handle_monomers
from src.general_type_aliases import (
    SMILES,
)
from src.monomer_names_helper import (
    Chirality,
    NorineMonomerName,
)
from networkx import DiGraph, is_isomorphic

from collections import defaultdict
from dataclasses import dataclass, asdict

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

    @classmethod
    def from_dict(cls, data: dict) -> MonomerInfo:
        return cls(name=NorineMonomerName(data['name']),
                   atoms=[AtomId(atom_id) for atom_id in data['atoms']],
                   chirality=Chirality[data['chirality']],
                   is_pks_hybrid=data.get('is_pks_hybrid', False))

    def to_dict(self) -> dict:
        d = self._asdict()
        d['chirality'] = self.chirality.name
        return d


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
class NRP_metadata:
    name: Optional[str]
    smiles: Optional[SMILES]
    origin: Optional[str]
    inchikey: Optional[str]
    source: Optional[str]

    @classmethod
    def from_dict(cls, data: dict) -> NRP_metadata:
        fields = cls.__dataclass_fields__.keys()
        class_data = {field: None for field in fields}
        try:
            for key, value in data.items():
                if key.lower() in fields:
                    class_data[key.lower()] = data[key]
        except:
            raise ValueError(f"Error parsing metadata:\n {data}") from None
        return cls(**class_data)


@dataclass
class Parsed_rBAN_Record:
    compound_id: str
    monomers: Dict[MonomerIdx, MonomerInfo]
    monomer_bonds: Dict[MonomerEdge, MonomerEdgeInfo]
    atoms: Dict[AtomId, AtomInfo]
    atomic_bonds: Dict[AtomicEdge, AtomicEdgeInfo]
    metadata: NRP_metadata

    def __init__(self, rban_record: Raw_rBAN_Record,
                 hybrid_monomers: Dict[MonomerIdx, NorineMonomerName],
                 chiralities: Dict[MonomerIdx, Chirality],
                 metadata: Dict[str, NRP_metadata]):

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
        self.metadata = metadata[self.compound_id] if self.compound_id in metadata else NRP_metadata(None, None, None, None, None)
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

    def to_dict(self) -> dict:
        return {'compound_id': self.compound_id,
                'monomers': {mon_idx: mon_info.to_dict()
                             for mon_idx, mon_info in self.monomers.items()},
                'monomer_bonds': [[[u, v], edge_info._asdict()]  # saving dict as list of pairs because YAML can't have tuple keys
                                  for (u, v), edge_info in self.monomer_bonds.items()],
                'atoms': {atom_id: atom_info._asdict()
                          for atom_id, atom_info in self.atoms.items()},
                'atomic_bonds': [[[u, v], edge_info._asdict()]  # saving dict as list of pairs because YAML can't have tuple keys
                                 for (u, v), edge_info in self.atomic_bonds.items()],
                'metadata': asdict(self.metadata)}

    @classmethod
    def from_dict(cls, data: dict) -> Parsed_rBAN_Record:
        metadata = NRP_metadata.from_dict(data['metadata'])
        parsed_record = cls.__new__(cls)  # Create an uninitialized instance
        parsed_record.compound_id = data['compound_id']
        parsed_record.monomers = {int(mon_idx): MonomerInfo.from_dict(mon_info)
                                  for mon_idx, mon_info in data['monomers'].items()}
        parsed_record.monomer_bonds = {(int(u), int(v)): MonomerEdgeInfo(**edge_info)  # expecting dict as list of pairs because YAML can't have tuple keys
                                       for (u, v), edge_info in data['monomer_bonds']}
        parsed_record.atoms = {int(atom_id): AtomInfo(**atom_info)
                               for atom_id, atom_info in data['atoms'].items()}
        parsed_record.atomic_bonds = {(int(u), int(v)): AtomicEdgeInfo(**edge_info)
                                      for (u, v), edge_info in data['atomic_bonds']}  # expecting dict as list of pairs because YAML can't have tuple keys
        parsed_record.metadata = metadata
        return parsed_record

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
