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
    NerpaResidue,
    MonCode,
    MonomerNamesHelper
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
    bond_type: str


MonomerIdx = int
MonomerEdge = Tuple[MonomerIdx, MonomerIdx]


class MonomerInfo(NamedTuple):
    name: NorineMonomerName
    atoms: List[AtomId]
    chirality: Chirality
    is_pks_hybrid: bool = False
    nerpa_core: Optional[NerpaResidue] = None
    methylated: Optional[bool] = None
    mon_code: Optional[MonCode] = None


    @classmethod
    def from_dict(cls, data: dict) -> MonomerInfo:
        return cls(name=NorineMonomerName(data['name']),
                   atoms=[AtomId(atom_id) for atom_id in data['atoms']],
                   chirality=Chirality[data['chirality']],
                   is_pks_hybrid=data.get('is_pks_hybrid', False),
                   nerpa_core=NerpaResidue(data['nerpa_core']) if data.get('nerpa_core') else None,
                   methylated=data.get('methylated'),
                   mon_code=MonCode(data['mon_code']) if data.get('mon_code') else None)

    def to_dict(self, monomer_names_helper: Optional[MonomerNamesHelper] = None) -> dict:
        d = self._asdict()
        if monomer_names_helper is not None:
            parsed_name = monomer_names_helper.parsed_name(self.name, name_format='rBAN/Norine')
            d['nerpa_core'] = parsed_name.residue
            d['methylated'] = parsed_name.methylated
            d['mon_code'] = monomer_names_helper.mon_to_int[parsed_name]
        d['chirality'] = self.chirality.name
        return d

# awkward design because I am retroactively adding support for multiple bonds
# between the same monomers, and I don't want to change the existing code too much.
class _MonomerEdgeInfo(NamedTuple):
    monomer_to_atom: Dict[MonomerIdx, AtomId]
    arity: int  # I heard that there exist fractional arities (e.g. 1.5)
    bond_type: str

class MonomerEdgeInfo(NamedTuple):
    monomer_to_atom: Dict[MonomerIdx, AtomId]
    arity: int  # I heard that there exist fractional arities (e.g. 1.5)
    bond_type: str
    all_edges: List[_MonomerEdgeInfo]  # in case of multiple bonds between the same monomers, we need to store all of them


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
                    class_data[key.lower()] = value
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
                                                bond_type=bond['bondType'])
                             for bond in rban_record['atomicGraph']['atomicGraph']['bonds']}

        atom_to_monomer = {
            atom_id: mon_idx
            for mon_idx, mon_info in self.monomers.items()
            for atom_id in mon_info.atoms
        }

        self.monomer_bonds = {}
        for (atom1, atom2), bond_info in self.atomic_bonds.items():
            mon1, mon2 = atom_to_monomer[atom1], atom_to_monomer[atom2]
            if mon1 == mon2:
                continue  # internal bond, not a monomer bond

            _monomer_edge_info = _MonomerEdgeInfo(monomer_to_atom={mon1: atom1, mon2: atom2},
                                                  arity=bond_info.arity,
                                                  bond_type=bond_info.bond_type)
            if (bond_info.bond_type == 'AMINO' and self.atoms[atom1].name == 'N') \
               or (mon2, mon1) in self.monomer_bonds:
                mon1, mon2 = mon2, mon1  # ensure that the amino bond is always directed from C to N

            if (mon1, mon2) not in self.monomer_bonds:
                self.monomer_bonds[(mon1, mon2)] = MonomerEdgeInfo(**{'all_edges': [],
                                                                   **_monomer_edge_info._asdict()})
            self.monomer_bonds[(mon1, mon2)].all_edges.append(_monomer_edge_info)

    def to_compact_dict(self) -> dict:
        return {'compound_id': self.compound_id,
                'nodes': {mon_idx: mon_info.name
                          for mon_idx, mon_info in self.monomers.items()},
                'edges': [(u, v, edge_info.bond_type)
                          for (u, v), edge_info in self.monomer_bonds.items()]}

    def to_dict(self, monomer_names_helper: Optional[MonomerNamesHelper] = None) -> dict:
        return {'compound_id': self.compound_id,
                'monomers': {mon_idx: mon_info.to_dict(monomer_names_helper)
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
            G.add_edge(u, v, bond_type=edge_info.bond_type)
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

