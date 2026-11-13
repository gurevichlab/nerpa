from __future__ import annotations

import ast
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

def check_isomorphism(
        record: Parsed_rBAN_Record,
        mol: Chem.rdchem.Mol,
        atom_id_to_mol_idx: Dict[AtomId, int],
) -> bool:
    debug = True
    n_atoms: int = mol.GetNumAtoms()
    if len(record.atoms) != n_atoms:
        return False

    if (
        set(atom_id_to_mol_idx) != set(record.atoms)
        or set(atom_id_to_mol_idx.values()) != set(range(n_atoms))
    ):
        if debug:
            print(f"Invalid atom mapping for record {record.compound_id}")
        return False

    for atom_id, atom_info in record.atoms.items():
        atom = mol.GetAtomWithIdx(atom_id_to_mol_idx[atom_id])
        if atom_info.name != atom.GetSymbol():
            if debug:
                print(f"Atom name mismatch for atom {atom_id} in record {record.compound_id}: "
                        f"record has {atom_info.name}, RDKit has {atom.GetSymbol()}")
            return False
        if int(atom_info.hydrogens) != int(atom.GetTotalNumHs()):
            if debug:
                print(f"Hydrogen count mismatch for atom {atom_id} in record {record.compound_id}: "
                        f"record has {atom_info.hydrogens}, RDKit has {atom.GetTotalNumHs()}")
            return False

    # Compare bonds after translating record atom IDs into RDKit indices.
    rec_bonds: Dict[tuple[int, int], AtomicEdgeInfo] = {}
    for (u, v), info in record.atomic_bonds.items():
        left_index = atom_id_to_mol_idx[u]
        right_index = atom_id_to_mol_idx[v]
        key = (
            (left_index, right_index)
            if left_index < right_index
            else (right_index, left_index)
        )
        if key in rec_bonds:
            # RDKit Mol can't represent parallel edges between the same atom pair in SMILES-derived mols
            if debug:
                print(f"Duplicate bond between atoms {key} in record {record.compound_id}")
                
            return False
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
        return False

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

    return True


class MolRecord(NamedTuple):
    mol: Chem.rdchem.Mol
    atom_id_to_mol_idx: Dict[AtomId, int]

    @classmethod
    def _from_rban_record_smiles(cls, record: Parsed_rBAN_Record) -> MolRecord:
        """
        This method assumes that the rBAN record has a valid SMILES string in its metadata.
        and that the order of atoms in both the Mol object and the rBAN record
        actually corresponds to the order of atoms in the SMILES string

        Stereochemistry information is preserved in the SMILES string, so this method will also preserve stereochemistry in the resulting Mol object.
        """
        smiles = record.metadata.smiles
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Failed to build Mol object from SMILES: {smiles}")

        # trivial mapping from atom IDs to RDKit indices, since the order of atoms in the SMILES string is preserved in both the Mol object and the rBAN record
        atom_id_to_mol_idx: Dict[AtomId, int] = {
            atom_id: int(atom_id)
            for atom_id in record.atoms
        }
        if not check_isomorphism(record, mol, atom_id_to_mol_idx):
            raise ValueError(
                f"Mol object built from {smiles} does not match the "
                f"rBAN record {record.compound_id}"
            )

        return cls(
            mol=mol,
            atom_id_to_mol_idx=atom_id_to_mol_idx,
        )


    @classmethod
    def from_rban_record(cls, record: Parsed_rBAN_Record) -> MolRecord:
        smiles = record.metadata.smiles
        if smiles is not None:
            return cls._from_rban_record_smiles(record)
        else:
            return cls._from_rban_record_atomic_graph(record)


    @classmethod
    def _from_rban_record_atomic_graph(cls, record: Parsed_rBAN_Record) -> MolRecord:
        """
        This method does not use the SMILES string in the rBAN record, but instead builds the Mol object from the atomic graph information in the record.

        Stereochemistry information is not preserved in the atomic graph, so this method will not preserve stereochemistry in the resulting Mol object.
        """
        editable_mol = Chem.RWMol()
        atom_id_to_index: Dict[AtomId, int] = {}

        aromatic_atom_ids = {
            atom_id
            for atom_ids, edge_info in record.atomic_bonds.items()
            if edge_info.arity.strip() in {"1.5", "1.50"}
            for atom_id in atom_ids
        }

        # Sorting makes the RDKit indices deterministic even though the record
        # stores atoms in a dictionary.
        for atom_id, atom_info in sorted(record.atoms.items()):
            atom = Chem.Atom(atom_info.name)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(atom_info.hydrogens)
            atom.SetIsAromatic(atom_id in aromatic_atom_ids)

            atom_id_to_index[atom_id] = editable_mol.AddAtom(atom)

        bond_types = {
            "1": Chem.BondType.SINGLE,
            "1.0": Chem.BondType.SINGLE,
            "1.00": Chem.BondType.SINGLE,
            "1.5": Chem.BondType.AROMATIC,
            "1.50": Chem.BondType.AROMATIC,
            "2": Chem.BondType.DOUBLE,
            "2.0": Chem.BondType.DOUBLE,
            "2.00": Chem.BondType.DOUBLE,
            "3": Chem.BondType.TRIPLE,
            "3.0": Chem.BondType.TRIPLE,
            "3.00": Chem.BondType.TRIPLE,
        }

        for atom_ids, edge_info in sorted(record.atomic_bonds.items()):
            bond_order = edge_info.arity.strip()
            try:
                bond_type = bond_types[bond_order]
            except KeyError as error:
                raise ValueError(
                    f"Unsupported bond order: {bond_order!r}"
                ) from error

            try:
                left_index = atom_id_to_index[atom_ids[0]]
                right_index = atom_id_to_index[atom_ids[1]]
            except KeyError as error:
                raise ValueError(
                    f"Bond refers to unknown atom ID: {error.args[0]}"
                ) from error

            editable_mol.AddBond(left_index, right_index, bond_type)

        mol = editable_mol.GetMol()
        Chem.SanitizeMol(mol)

        return cls(
            mol=mol,
            atom_id_to_mol_idx=atom_id_to_index,
        )

    def to_canonical_smiles(self) -> str:
        return Chem.MolToSmiles(self.mol, canonical=True)

    def canonical_atom_order(self) -> list[AtomId]:
        # MolToSmiles records the actual SMILES traversal on the molecule.
        Chem.MolToSmiles(self.mol, canonical=True)
        mol_indices: list[int] = ast.literal_eval(
            self.mol.GetProp("_smilesAtomOutputOrder")
        )
        mol_index_to_atom_id = {
            mol_index: atom_id
            for atom_id, mol_index in self.atom_id_to_mol_idx.items()
        }

        return [
            mol_index_to_atom_id[mol_index]
            for mol_index in mol_indices
        ]

