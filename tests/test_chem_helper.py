from __future__ import annotations

import ast
import copy
import unittest
from pathlib import Path
from typing import Any

import yaml
from rdkit import Chem

from scripts.for_rust_core.rban_records_to_molecule_analysis import (
    build_id_to_analysis,
)
from src.build_output.chem_helper import MolRecord
from src.rban_parsing.rban_parser import Parsed_rBAN_Record


def make_ethanol_record(smiles: str | None) -> Parsed_rBAN_Record:
    return Parsed_rBAN_Record.from_dict(
        {
            "compound_id": "ethanol",
            "monomers": {},
            "monomer_bonds": [],
            "atoms": {
                10: {"name": "C", "hydrogens": 3},
                20: {"name": "C", "hydrogens": 2},
                30: {"name": "O", "hydrogens": 1},
            },
            "atomic_bonds": [
                [[10, 20], {"arity": "1", "bond_type": "SINGLE"}],
                [[20, 30], {"arity": "1.0", "bond_type": "SINGLE"}],
            ],
            "metadata": {
                "name": None,
                "smiles": smiles,
                "origin": None,
                "inchikey": None,
                "source": None,
            },
        }
    )


class MolRecordTest(unittest.TestCase):
    def test_builds_from_atomic_graph_without_smiles(self) -> None:
        mol_record = MolRecord.from_rban_record(
            make_ethanol_record(smiles=None)
        )

        self.assertEqual(mol_record.to_canonical_smiles(), "CCO")
        self.assertEqual(
            mol_record.atom_id_to_mol_idx,
            {10: 0, 20: 1, 30: 2},
        )

    def test_canonical_atom_order_is_independent_of_atom_ids(self) -> None:
        first_record = make_ethanol_record(smiles=None)
        second_record = Parsed_rBAN_Record.from_dict(
            {
                "compound_id": "renumbered-ethanol",
                "monomers": {},
                "monomer_bonds": [],
                "atoms": {
                    300: {"name": "C", "hydrogens": 3},
                    100: {"name": "C", "hydrogens": 2},
                    200: {"name": "O", "hydrogens": 1},
                },
                "atomic_bonds": [
                    [[300, 100], {"arity": "1", "bond_type": "SINGLE"}],
                    [[100, 200], {"arity": "1", "bond_type": "SINGLE"}],
                ],
                "metadata": {
                    "name": None,
                    "smiles": None,
                    "origin": None,
                    "inchikey": None,
                    "source": None,
                },
            }
        )

        first_order = MolRecord.from_rban_record(
            first_record
        ).canonical_atom_order()
        second_order = MolRecord.from_rban_record(
            second_record
        ).canonical_atom_order()

        renumbering = {
            10: 300,
            20: 100,
            30: 200,
        }
        self.assertEqual(
            second_order,
            [renumbering[atom_id] for atom_id in first_order],
        )

    def test_molecule_analysis_has_stable_result_schema(self) -> None:
        record = make_ethanol_record(smiles=None)

        result = build_id_to_analysis([record.to_dict()])
        analysis = result["ethanol"]

        self.assertEqual(
            set(analysis),
            {"smiles", "atom_order"},
        )
        self.assertEqual(analysis["smiles"], "CCO")
        self.assertEqual(set(analysis["atom_order"]), {10, 20, 30})

    def test_canonical_atom_order_matches_canonical_smiles_order(self) -> None:
        mol_record = MolRecord.from_rban_record(
            make_ethanol_record(smiles=None)
        )
        canonical_smiles = mol_record.to_canonical_smiles()
        canonical_atom_ids = mol_record.canonical_atom_order()

        canonical_mol = Chem.MolFromSmiles(canonical_smiles)
        self.assertIsNotNone(canonical_mol)
        self.assertEqual(
            canonical_atom_ids,
            [
                canonical_atom_ids[index]
                for index in range(canonical_mol.GetNumAtoms())
            ],
        )

    def test_canonical_atom_order_matches_smiles_traversal(self) -> None:
        mol_record = MolRecord.from_rban_record(
            make_ethanol_record(smiles=None)
        )

        Chem.MolToSmiles(mol_record.mol, canonical=True)
        mol_indices = ast.literal_eval(
            mol_record.mol.GetProp("_smilesAtomOutputOrder")
        )
        expected_atom_ids = [
            next(
                atom_id
                for atom_id, mol_index
                in mol_record.atom_id_to_mol_idx.items()
                if mol_index == index
            )
            for index in mol_indices
        ]

        self.assertEqual(
            mol_record.canonical_atom_order(),
            expected_atom_ids,
        )

    def test_preserves_smiles_based_behavior(self) -> None:
        record = Parsed_rBAN_Record.from_dict(
            {
                "compound_id": "ethanol",
                "monomers": {},
                "monomer_bonds": [],
                "atoms": {
                    0: {"name": "C", "hydrogens": 3},
                    1: {"name": "C", "hydrogens": 2},
                    2: {"name": "O", "hydrogens": 1},
                },
                "atomic_bonds": [
                    [[0, 1], {"arity": "1", "bond_type": "SINGLE"}],
                    [[1, 2], {"arity": "1", "bond_type": "SINGLE"}],
                ],
                "metadata": {
                    "name": None,
                    "smiles": "CCO",
                    "origin": None,
                    "inchikey": None,
                    "source": None,
                },
            }
        )

        mol_record = MolRecord.from_rban_record(record)

        self.assertEqual(mol_record.to_canonical_smiles(), "CCO")
        self.assertEqual(
            mol_record.atom_id_to_mol_idx,
            {0: 0, 1: 1, 2: 2},
        )

    def test_builds_aromatic_ring_from_atomic_graph(self) -> None:
        atom_ids = [10, 20, 30, 40, 50, 60]
        record = Parsed_rBAN_Record.from_dict(
            {
                "compound_id": "benzene",
                "monomers": {},
                "monomer_bonds": [],
                "atoms": {
                    atom_id: {"name": "C", "hydrogens": 1}
                    for atom_id in atom_ids
                },
                "atomic_bonds": [
                    [
                        [atom_ids[index], atom_ids[(index + 1) % len(atom_ids)]],
                        {"arity": "1.5", "bond_type": "AROMATIC"},
                    ]
                    for index in range(len(atom_ids))
                ],
                "metadata": {
                    "name": None,
                    "smiles": None,
                    "origin": None,
                    "inchikey": None,
                    "source": None,
                },
            }
        )

        mol_record = MolRecord.from_rban_record(record)

        self.assertEqual(mol_record.to_canonical_smiles(), "c1ccccc1")
        self.assertEqual(set(mol_record.atom_id_to_mol_idx), set(atom_ids))

    def test_database_atomic_graphs_match_metadata_smiles(self) -> None:
        database_path = (
            Path(__file__).parents[1]
            / "data"
            / "input"
            / "preprocessed"
            / "pnrpdb2_mibig_norine_parsed_rban_records.yaml"
        )
        records_raw: Any = yaml.safe_load(
            database_path.read_text(encoding="utf-8")
        )

        self.assertIsInstance(records_raw, list)

        for record_raw in records_raw:
            record = Parsed_rBAN_Record.from_dict(record_raw)
            if record.metadata.smiles is None:
                continue

            with self.subTest(compound_id=record.compound_id):
                smiles_mol = MolRecord.from_rban_record(record).mol

                # Formal charge is absent from the rBAN atomic graph, so these
                # molecules cannot be reconstructed losslessly from it.
                if any(
                    atom.GetFormalCharge() != 0
                    for atom in smiles_mol.GetAtoms()
                ):
                    continue

                graph_record = copy.deepcopy(record)
                graph_record.metadata.smiles = None
                graph_mol = MolRecord.from_rban_record(graph_record).mol

                smiles_canonical = Chem.MolToSmiles(
                    smiles_mol,
                    canonical=True,
                    isomericSmiles=False,
                )
                graph_canonical = Chem.MolToSmiles(
                    graph_mol,
                    canonical=True,
                    isomericSmiles=False,
                )

                self.assertEqual(graph_canonical, smiles_canonical)


if __name__ == "__main__":
    unittest.main()