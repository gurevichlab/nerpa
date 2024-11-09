from typing import Dict
from dataclasses import dataclass
from src.monomer_names_helper import MonomerResidue



@dataclass
class NorineStats:
    total_monomers: int
    methylated: int
    d_chirality: int
    residue_frequencies: Dict[MonomerResidue, float]