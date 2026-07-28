from __future__ import annotations
from typing import (
    Callable,
    Dict,
    List,
    NamedTuple,
    Tuple,
    Union,
    Optional,
    Hashable,
)
from src.rban_parsing import handle_monomers
from src.general_type_aliases import (
    SMILES,
)
from src.generic.parsing import read_int_pair
from src.monomer_names_helper import (
    Chirality,
    NorineMonomerName,
    NerpaResidue,
    MonCode,
    MonomerNamesHelper
)
import networkx as nx

from collections import defaultdict
from dataclasses import dataclass, asdict

Raw_rBAN_Record = dict

AtomId = int
AtomicEdge = Tuple[AtomId, AtomId]


class AtomInfo(NamedTuple):
    name: str
    hydrogens: int


class AtomicEdgeInfo(NamedTuple):
    arity: str  # I heard that there exist fractional arities (e.g. 1.5), so I use str to be safe
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

class MonomerEdgeInfoSingle(NamedTuple):
    monomer_to_atom: Dict[MonomerIdx, AtomId]  # should have exactly 2 items, but I use dict for convenience of accessing the atoms by monomer id
    atomic_edge: AtomicEdgeInfo

    @classmethod
    def from_dict(cls, data: dict) -> MonomerEdgeInfoSingle:
        try:
            result =  cls(monomer_to_atom={int(mon_id): AtomId(atom_id)
                                           for mon_id, atom_id in data['monomer_to_atom'].items()},
                          atomic_edge=AtomicEdgeInfo(**data['atomic_edge']))
        except Exception as e:
            raise ValueError(f"Error parsing MonomerEdgeInfoSingle from dict:\n {data}") from e
        return result

    def to_dict(self) -> dict:
        return {'monomer_to_atom': self.monomer_to_atom,
                'atomic_edge': self.atomic_edge._asdict()}


MonomerEdgeInfo = List[MonomerEdgeInfoSingle]  # list of dicts because there can be multiple bonds between the same monomers

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
        self.metadata = metadata.get(self.compound_id, NRP_metadata(None, None, None, None, None))
        self.atoms = {atom['cdk_idx']: AtomInfo(name=atom['name'],
                                                hydrogens=atom['hydrogens'])
                      for atom in rban_record['atomicGraph']['atomicGraph']['atoms']}

        self.monomers = {
            MonomerIdx(idx := monomer['monomer']['index']):
            MonomerInfo(name=get_monomer_name(idx),
                        atoms=[AtomId(atom) for atom in monomer['monomer']['atoms']],
                        chirality=chiralities[idx],
                        is_pks_hybrid=idx in hybrid_monomers)
            for monomer in rban_record['monomericGraph']['monomericGraph']['monomers']
        }

        self.atomic_bonds = {(AtomId(bond['atoms'][0]), AtomId(bond['atoms'][1])):
                                 AtomicEdgeInfo(arity=str(bond['arity']),
                                                bond_type=bond['bondType'])
                             for bond in rban_record['atomicGraph']['atomicGraph']['bonds']}

        atom_to_monomer = {
            atom_id: mon_idx
            for mon_idx, mon_info in self.monomers.items()
            for atom_id in mon_info.atoms
        }

        self.monomer_bonds = defaultdict(list)
        for (atom1, atom2), bond_info in self.atomic_bonds.items():
            mon1, mon2 = atom_to_monomer[atom1], atom_to_monomer[atom2]
            if mon1 == mon2:
                continue  # internal bond, not a monomer bond

            if mon1 > mon2:
                mon1, mon2 = mon2, mon1  # canonicalize order of monomers in the bond
                atom1, atom2 = atom2, atom1  # swap atoms accordingly

            self.monomer_bonds[(mon1, mon2)].append(
                MonomerEdgeInfoSingle(monomer_to_atom={mon1: atom1, mon2: atom2},
                                      atomic_edge=bond_info)
            )

    def to_dict(self, monomer_names_helper: Optional[MonomerNamesHelper] = None) -> dict:
        return {'compound_id': self.compound_id,
                'monomers': {mon_idx: mon_info.to_dict(monomer_names_helper)
                             for mon_idx, mon_info in self.monomers.items()},
                'monomer_bonds': [
                    [
                        [u, v],
                        [edge_info.to_dict() for edge_info in edge_info_list] 
                    ]
                    for (u, v), edge_info_list in self.monomer_bonds.items()
                ],
                'atoms': {atom_id: atom_info._asdict()
                          for atom_id, atom_info in self.atoms.items()},
                'atomic_bonds': [[[u, v], edge_info._asdict()]  # saving dict as list of pairs because YAML can't have tuple keys
                                 for (u, v), edge_info in self.atomic_bonds.items()],
                'metadata': asdict(self.metadata)}

    @classmethod
    def from_dict(cls, data: dict) -> Parsed_rBAN_Record:
        metadata = NRP_metadata.from_dict(data['metadata'])
        compound_id = data['compound_id']
        monomers = {int(mon_idx): MonomerInfo.from_dict(mon_info)
                                  for mon_idx, mon_info in data['monomers'].items()}
        atoms = {int(atom_id): AtomInfo(**atom_info)
                               for atom_id, atom_info in data['atoms'].items()}
        monomer_bonds = {
            (int(u), int(v)):
            [
                MonomerEdgeInfoSingle.from_dict(edge_info)
                for edge_info in edge_info_list
            ]
            for (u, v), edge_info_list in data['monomer_bonds']
        }
            
        atomic_bonds = {(int(u), int(v)): AtomicEdgeInfo(**edge_info)
                        for (u, v), edge_info in data['atomic_bonds']}
        parsed_record = Parsed_rBAN_Record.__new__(cls)  # create an uninitialized instance
        parsed_record.compound_id = compound_id
        parsed_record.monomers = monomers
        parsed_record.atoms = atoms
        parsed_record.monomer_bonds = monomer_bonds
        parsed_record.atomic_bonds = atomic_bonds
        parsed_record.metadata = metadata
        return parsed_record

    def to_nx_monomer_graph(self) -> nx.Graph:
        graph = nx.Graph()
        for monomer_idx, monomer_info in self.monomers.items():
            graph.add_node(monomer_idx, data=monomer_info)
        for (mon1_idx, mon2_idx), edges in self.monomer_bonds.items():
            edges_data = [edge_info.atomic_edge for edge_info in edges]
            graph.add_edge(mon1_idx, mon2_idx, data=edges_data)
        return graph

    @classmethod
    def node_key_default(cls, monomer_info: MonomerInfo) -> Hashable:
        old_name = monomer_info.name
        if old_name.startswith('X'):
            new_name = 'X'  # unknown
        elif ':' in old_name:
            new_name = 'LIPID'  # lipid tail
        else:
            new_name = old_name  # regular residue
            
        return (
            new_name,
            monomer_info.chirality,
            monomer_info.is_pks_hybrid
        )

    @classmethod
    def edge_key_default(cls, edges_info: List[AtomicEdgeInfo]) -> Hashable:
        return tuple(sorted(edge_info.bond_type for edge_info in edges_info))

    @classmethod
    def _edge_key_from_mon_edge_info(cls, mon_edge_info: MonomerEdgeInfo) -> Hashable:
        atomic_edges_info = [edge_info.atomic_edge for edge_info in mon_edge_info]
        return cls.edge_key_default(atomic_edges_info)

    def is_isomorphic_to(
            self,
            other: Parsed_rBAN_Record,
            node_key: Optional[Callable[[MonomerInfo], Hashable]] = None,
            edge_key: Optional[Callable[[List[AtomicEdgeInfo]], Hashable]] = None
    ) -> bool:
        if node_key is None:
            node_key = self.node_key_default
        if edge_key is None:
            edge_key = self.edge_key_default

        # compare node and edge keys before committing to full isomorphism check
        self_node_keys = {node_key(mon_info) for mon_info in self.monomers.values()}
        other_node_keys = {node_key(mon_info) for mon_info in other.monomers.values()}
        if self_node_keys != other_node_keys:
            return False

        self_edge_keys = {
            self._edge_key_from_mon_edge_info(edges_info)
            for edges_info in self.monomer_bonds.values()
        }
        other_edge_keys = {
            other._edge_key_from_mon_edge_info(edges_info)
            for edges_info in other.monomer_bonds.values()
        }
        if self_edge_keys != other_edge_keys:
            return False

        # Now perform the full isomorphism check using networkx
        graph1 = self.to_nx_monomer_graph()
        graph2 = other.to_nx_monomer_graph()
        
        return nx.is_isomorphic(
            graph1,
            graph2,
            node_match=lambda n1, n2: node_key(n1['data']) == node_key(n2['data']),
            edge_match=lambda e1, e2: edge_key(e1['data']) == edge_key(e2['data'])
        )


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

