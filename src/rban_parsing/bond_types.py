from __future__ import annotations
from enum import StrEnum
from typing import List, Optional, Tuple

from src.pipeline.logging.logger import NerpaLogger
from src.monomer_names_helper.monomer_names_helper import (
    MonomerNamesHelper,
    MonomerNamesTableSchema,
    MonomerTag,
    MonomerNameFormat,
)
from src.rban_parsing.rban_parser import (
    MonomerEdgeInfoSingle,
    MonomerIdx,
    Parsed_rBAN_Record,
)


class MonomerBondType(StrEnum):
    """Bond types for chemical bonds."""
    
    AMINO = "AMINO"
    AMINO_TRANS = "AMINO_TRANS"
    OXAZOLE_CYCLE = "OXAZOLE_CYCLE"
    THIAZOLE_CYCLE = "THIAZOLE_CYCLE"
    PYRIMIDINE_CYCLE = "PYRIMIDINE_CYCLE"
    HETEROCYCLE = "HETEROCYCLE"

    IN_CHAIN_ESTER = "IN_CHAIN_ESTER"

    ESTER = "ESTER"

    UNKNOWN = "UNKNOWN"

    def is_backbone_bond(self) -> bool:
        """Check if the bond type is a backbone bond."""
        return self in {
            MonomerBondType.AMINO,
            MonomerBondType.AMINO_TRANS,
            MonomerBondType.OXAZOLE_CYCLE,
            MonomerBondType.THIAZOLE_CYCLE,
            MonomerBondType.PYRIMIDINE_CYCLE,
            MonomerBondType.HETEROCYCLE,
            MonomerBondType.IN_CHAIN_ESTER,
        }

    @classmethod
    def classify_bond(
            cls,
            rban_record: Parsed_rBAN_Record,
            monomer_bond: Tuple[MonomerIdx, MonomerIdx],
            monomer_names_helper: MonomerNamesHelper,
            logger: Optional[NerpaLogger] = None
    ) -> tuple[MonomerBondType, Tuple[MonomerIdx, MonomerIdx] | None]:
        """
        Returns the bond type and its direction (if can be inferred).
        The direction here means the order in the biosynthetic pathway.
        Direction idx1 -> idx2 means that idx1 was incorporated before idx2 in the NRP assembly line.
        The direction is either monomer_bond or reversed monomer_bond, or None if it cannot be inferred.
        """
        
        try:
            edge_info: List[MonomerEdgeInfoSingle] = rban_record.monomer_bonds[monomer_bond]
        except KeyError:
            raise ValueError(f"Bond {monomer_bond} not found in record")

        mon1, mon2 = [
            rban_record.monomers[mon_idx]
            for mon_idx in monomer_bond
        ]
        tags_per_mon = []
        for mon in (mon1, mon2):
            entry = monomer_names_helper.get_table_entry(
                mon.name, name_format=MonomerNameFormat.RBAN_NORINE
            )
            tags_per_mon.append(
                entry[MonomerNamesTableSchema.TAGS]
                if entry is not None
                else []
            )

        rban_bond_types: list[str] = [
            edge.atomic_edge.bond_type
            for edge in edge_info
        ]
        if len(set(rban_bond_types)) > 1 and logger is not None:
            logger.warning(
                f"Ambiguous bond type for bond {monomer_bond} "
                f"in record {rban_record.compound_id}: {rban_bond_types}."
                "Treating as UNKNOWN type and direction."
            )
            return MonomerBondType.UNKNOWN, None

        rban_bond_type = rban_bond_types[0]
        match rban_bond_type:
            case "AMINO":
                if len(edge_info) > 1:
                    # multiple atomic bonds between the same monomer pair
                    # currently unsupported for direction inference, so we return None
                    return MonomerBondType.AMINO, None

                atom_id1, atom_id2 = [
                    edge_info[0].monomer_to_atom[mon_idx]
                    for mon_idx in monomer_bond
                ]
                atom1, atom2 = [
                    rban_record.atoms[atom_id].name
                    for atom_id in (atom_id1, atom_id2)
                ]
                    
                match atom1, atom2:
                    case ('C', 'N'):
                        return MonomerBondType.AMINO, monomer_bond
                    case ('N', 'C'):
                        return MonomerBondType.AMINO, (monomer_bond[1], monomer_bond[0])
                    case _:
                        raise ValueError(
                            "Unexpected atom names for amino bond between atoms "
                            f"{atom_id1} and {atom_id2}: "
                            f"{atom1}, {atom2}"
                        )

            case "ESTER":
                if any(MonomerTag.ALPHA_HYDROXY_AMINO_ACID in tags
                        for tags in tags_per_mon):
                    return MonomerBondType.IN_CHAIN_ESTER, None
                else:
                    return MonomerBondType.ESTER, None
            case _:
                try:
                    bond_type = cls(rban_bond_type)
                except ValueError:
                    return MonomerBondType.UNKNOWN, None
                return bond_type, None

