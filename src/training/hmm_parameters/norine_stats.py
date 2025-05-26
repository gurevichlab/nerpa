from typing import Dict
from dataclasses import dataclass
from src.monomer_names_helper import MonomerResidue
import yaml
import dacite
from pathlib import Path

@dataclass
class NorineStats:
    total_monomers: int
    methylated: int
    d_chirality: int
    monomer_frequencies: Dict[MonomerResidue, float]


def load_norine_stats(norine_stats_file: Path) -> NorineStats:
    return dacite.from_dict(NorineStats, yaml.safe_load(norine_stats_file.read_text()))