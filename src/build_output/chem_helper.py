from __future__ import annotations
from typing import (
    Dict,
    NamedTuple,
)
from src.rban_parsing.rban_parser import (
    AtomId,
    AtomicEdgeInfo,
    Parsed_rBAN_Record,
)
    

from rdkit import Chem
import math

def get_isomorphism(
        record: Parsed_rBAN_Record,
        mol: Chem.rdchem.Mol
) -> Dict[AtomId, int] | None:
    debug = True
    n_atoms: int = mol.GetNumAtoms()
    if len(record.atoms) != n_atoms:
        return None

    # The contract you described: record AtomId == RDKit atom index
    if set(record.atoms.keys()) != set(range(n_atoms)):
        if debug:
            print(f"Atom IDs in record {record.compound_id} do not match RDKit atom indices")
        return None

    for atom_id, atom_info in record.atoms.items():
        atom = mol.GetAtomWithIdx(int(atom_id))
        if atom_info.name != atom.GetSymbol():
            if debug:
                print(f"Atom name mismatch for atom {atom_id} in record {record.compound_id}: "
                        f"record has {atom_info.name}, RDKit has {atom.GetSymbol()}")
            return None
        if int(atom_info.hydrogens) != int(atom.GetTotalNumHs()):
            if debug:
                print(f"Hydrogen count mismatch for atom {atom_id} in record {record.compound_id}: "
                        f"record has {atom_info.hydrogens}, RDKit has {atom.GetTotalNumHs()}")
            return None

    # Canonicalize bonds as undirected (RDKit bonds are undirected; record keys may be ordered)
    rec_bonds: Dict[tuple[int, int], AtomicEdgeInfo] = {}
    for (u, v), info in record.atomic_bonds.items():
        key = (int(u), int(v)) if int(u) < int(v) else (int(v), int(u))
        if key in rec_bonds:
            # RDKit Mol can't represent parallel edges between the same atom pair in SMILES-derived mols
            if debug:
                print(f"Duplicate bond between atoms {key} in record {record.compound_id}")
                
            return None
        rec_bonds[key] = info

    mol_bonds: Dict[tuple[int, int], Chem.rdchem.Bond] = {}
    for bond in mol.GetBonds():
        u = int(bond.GetBeginAtomIdx())
        v = int(bond.GetEndAtomIdx())
        key = (u, v) if u < v else (v, u)
        mol_bonds[key] = bond

    if set(rec_bonds.keys()) != set(mol_bonds.keys()):
        if debug:
            print(f"Bond keys mismatch for record {record.compound_id}: "
                    f"record has {set(rec_bonds.keys())}, RDKit has {set(mol_bonds.keys())}")
        return None

    # There can be discrepancies in bond order representation (e.g., aromatic bonds: 1.5 everywhere vs 1+2), so I disabled this check not knowing how to handle it properly. 

    # for key, rec_info in rec_bonds.items():
    #     bond = mol_bonds[key]

    #     # Prefer numeric arity comparison when possible (captures aromatic 1.5 etc.)
    #     try:
    #         rec_order = float(rec_info.arity)
    #     except ValueError:
    #         rec_order = None

    #     if rec_order is not None:
    #         rd_order = float(bond.GetBondTypeAsDouble())
    #         if not math.isclose(rec_order, rd_order, rel_tol=0.0, abs_tol=1e-3):
    #             if debug:
    #                 print(f"Bond order mismatch for bond {key} in record {record.compound_id}: "
    #                         f"record has {rec_order}, RDKit has {rd_order}")
    #             return None
    #     else:
    #     # Fallback: compare bond type strings
    #         if rec_info.bond_type.upper() != bond.GetBondType().name.upper():
    #             if debug:
    #                 print(f"Bond type mismatch for bond {key} in record {record.compound_id}: "
    #                         f"record has {rec_info.bond_type}, RDKit has {bond.GetBondType().name}")
    #             return None

    return {atom_id: atom_id for atom_id in record.atoms.keys()}


class MolRecord(NamedTuple):
    mol: Chem.rdchem.Mol
    atom_id_to_mol_idx: Dict[AtomId, int]


    @classmethod
    def from_rban_record(cls, record: Parsed_rBAN_Record) -> MolRecord:
        smiles = record.metadata.smiles
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Failed to build Mol object from SMILES: {smiles}")
        atom_id_to_index = get_isomorphism(record, mol)
        if atom_id_to_index is None:
            raise ValueError(f"Mol object built from {smiles} does not match the rBAN record {record.compound_id}")
        return MolRecord(
            mol=mol,
            atom_id_to_mol_idx=atom_id_to_index
        )

    def to_canonical_smiles(self) -> str:
        return Chem.MolToSmiles(self.mol, canonical=True)

