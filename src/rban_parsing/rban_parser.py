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
from src.generic.parsing import read_int_pair
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

MonomerEdgeInfo = List[Dict[MonomerIdx, AtomId]]  # list of dicts because there can be multiple bonds between the same monomers

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
            self.monomer_bonds[(mon1, mon2)].append({mon1: atom1, mon2: atom2})

    def to_dict(self, monomer_names_helper: Optional[MonomerNamesHelper] = None) -> dict:
        return {'compound_id': self.compound_id,
                'monomers': {mon_idx: mon_info.to_dict(monomer_names_helper)
                             for mon_idx, mon_info in self.monomers.items()},
                'monomer_bonds': [[[u, v], edge_info]  # saving dict as list of pairs because YAML can't have tuple keys
                                  for (u, v), edge_info in self.monomer_bonds.items()],
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
        monomer_bonds = {(int(u), int(v)): edge_info
                            for (u, v), edge_info in data['monomer_bonds']}
            
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

