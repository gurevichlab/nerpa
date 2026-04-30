from __future__ import annotations
from typing import (
    Dict,
    NamedTuple,
)
from src.rban_parsing.rban_parser import (
    AtomId,
    Parsed_rBAN_Record,
)
    

from rdkit import Chem


class MolRecord(NamedTuple):
    mol: Chem.rdchem.Mol
    atom_id_to_mol_idx: Dict[AtomId, int]

    @classmethod
    def from_rban_record(cls, record: Parsed_rBAN_Record) -> MolRecord:
        # create empty editable mol object
        mol = Chem.RWMol()

        # add atoms to mol and keep track of index
        atom_id_to_index = {}
        for atom_id, atom_info in record.atoms.items():
            atom_id_to_index[atom_id] = mol.AddAtom(Chem.Atom(atom_info.name))

        # add bonds between adjacent atoms
        for (u, v), edge_info in record.atomic_bonds.items():
            # add relevant bond type (there are many more of these)
            if edge_info.arity == "1":
                bond_type = Chem.rdchem.BondType.SINGLE
            elif edge_info.arity == "2":
                bond_type = Chem.rdchem.BondType.DOUBLE
            else:
                raise ValueError(f'Unknown or not implemented arity: {edge_info.arity}')
            mol.AddBond(atom_id_to_index[u], atom_id_to_index[v], bond_type)

        # Convert RWMol to Mol object
        mol = mol.GetMol()
        return MolRecord(mol=mol,
                         atom_id_to_mol_idx=atom_id_to_index)

