#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from typing import TypedDict

from rdkit import Chem
from rdkit.Chem.rdchem import BondType


class AtomInput(TypedDict):
    id: int
    element: str
    hydrogens: int


class BondInput(TypedDict):
    left: int
    right: int
    order: str


class GraphInput(TypedDict):
    atoms: list[AtomInput]
    bonds: list[BondInput]


class GraphOutput(TypedDict):
    smiles: str
    atom_order: list[int]


BOND_TYPES: dict[str, BondType] = {
    "1": BondType.SINGLE,
    "1.0": BondType.SINGLE,
    "1.00": BondType.SINGLE,
    "1.5": BondType.AROMATIC,
    "1.50": BondType.AROMATIC,
    "2": BondType.DOUBLE,
    "2.0": BondType.DOUBLE,
    "2.00": BondType.DOUBLE,
    "3": BondType.TRIPLE,
    "3.0": BondType.TRIPLE,
    "3.00": BondType.TRIPLE,
}


def build_molecule(graph: GraphInput) -> tuple[Chem.Mol, list[int]]:
    editable_molecule = Chem.RWMol()
    atom_ids: list[int] = []
    atom_id_to_index: dict[int, int] = {}

    aromatic_atom_ids: set[int] = {
        atom_id
        for bond in graph["bonds"]
        if bond["order"].strip() in {"1.5", "1.50"}
        for atom_id in (bond["left"], bond["right"])
    }

    for atom_data in graph["atoms"]:
        atom_id = atom_data["id"]
        if atom_id in atom_id_to_index:
            raise ValueError(f"Duplicate atom ID: {atom_id}")

        atom = Chem.Atom(atom_data["element"])
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(atom_data["hydrogens"])
        atom.SetIsAromatic(atom_id in aromatic_atom_ids)

        atom_index = editable_molecule.AddAtom(atom)
        atom_ids.append(atom_id)
        atom_id_to_index[atom_id] = atom_index

    for bond_data in graph["bonds"]:
        order = bond_data["order"].strip()

        try:
            bond_type = BOND_TYPES[order]
        except KeyError as error:
            raise ValueError(f"Unsupported bond order: {order!r}") from error

        try:
            left_index = atom_id_to_index[bond_data["left"]]
            right_index = atom_id_to_index[bond_data["right"]]
        except KeyError as error:
            raise ValueError(f"Bond refers to unknown atom ID: {error.args[0]}") from error

        editable_molecule.AddBond(left_index, right_index, bond_type)

        if bond_type == BondType.AROMATIC:
            bond = editable_molecule.GetBondBetweenAtoms(left_index, right_index)
            bond.SetIsAromatic(True)

    molecule = editable_molecule.GetMol()
    Chem.SanitizeMol(molecule)

    return molecule, atom_ids


def analyze_graph(graph: GraphInput) -> GraphOutput:
    molecule, atom_ids = build_molecule(graph)

    smiles = Chem.MolToSmiles(
        molecule,
        canonical=True,
        isomericSmiles=False,
    )
    atom_ranks = list(
        Chem.CanonicalRankAtoms(
            molecule,
            breakTies=True,
            includeChirality=False,
            includeIsotopes=False,
            includeAtomMaps=False,
        )
    )
    atom_order = [
        atom_ids[atom_index]
        for atom_index in sorted(
            range(len(atom_ids)),
            key=atom_ranks.__getitem__,
        )
    ]

    return {
        "smiles": smiles,
        "atom_order": atom_order,
    }


def main() -> None:
    graph: GraphInput = json.load(sys.stdin)
    json.dump(analyze_graph(graph), sys.stdout)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()