from __future__ import annotations
from typing import NamedTuple
from src.monomer_names_helper import (
    Chirality,
    MonomerResidue,
    NRP_Monomer,
    NorineMonomerName
)
#from monomer_features_types import MonomerFeature, MonomerFeatures
rBAN_idx = int

class rBAN_Monomer(NamedTuple):
    residue: MonomerResidue
    methylated: bool
    chirality: Chirality
    is_pks_hybrid: bool
    rban_name: NorineMonomerName
    rban_idx: rBAN_idx

    @classmethod
    def from_yaml_dict(cls, data: dict) -> rBAN_Monomer:
        if 'methylated' in data:
            methylated = data['methylated']
        else:
            methylated = 'METHYLATION' in data['modifications']
        if 'is_pks_hybrid' in data:
            is_pks_hybrid = data['is_hybrid']
        else:
            is_pks_hybrid = False
        return cls(residue=MonomerResidue(data['residue']),
                   methylated=methylated,
                   chirality=Chirality[data['chirality']],
                   is_pks_hybrid=is_pks_hybrid,
                   rban_name=NorineMonomerName(data['rban_name']),
                   rban_idx=int(data['rban_idx']))

    def to_base_mon(self) -> NRP_Monomer:
        return NRP_Monomer(residue=self.residue,
                           methylated=self.methylated,
                           chirality=self.chirality,
                           is_pks_hybrid=self.is_pks_hybrid)



rBAN_MONOMER_DUMMY = rBAN_Monomer(residue=MonomerResidue(''),
                                  methylated=False,
                                  chirality=Chirality.UNKNOWN,
                                  is_pks_hybrid=False,
                                  rban_name=NorineMonomerName(''),
                                  rban_idx=-1)

