from __future__ import annotations
from typing import NamedTuple
from src.monomer_names_helper import (
    Chirality,
    NerpaResidue,
    NRP_Monomer,
    NorineMonomerName
)
#from monomer_features_types import MonomerFeature, MonomerFeatures
rBAN_idx = int

def rban_name_no_unk_idx(rban_name: NorineMonomerName) -> NorineMonomerName:
    return (
        NorineMonomerName('X') if rban_name.startswith('X')
        else rban_name
    )

class rBAN_Monomer(NamedTuple):
    residue: NerpaResidue
    methylated: bool
    chirality: Chirality
    is_pks_hybrid: bool
    rban_name: NorineMonomerName
    rban_idx: rBAN_idx

    def to_dict(self) -> dict:
        return {'residue': self.residue,
                'methylated': self.methylated,
                'chirality': self.chirality.name,
                'is_pks_hybrid': self.is_pks_hybrid,
                'rban_name': self.rban_name,
                'rban_idx': self.rban_idx,}

    @classmethod
    def from_yaml_dict(cls, data: dict | list) -> rBAN_Monomer:
        if isinstance(data, list):
            return cls.from_list(data)

        if 'methylated' in data:
            methylated = data['methylated']
        else:
            methylated = 'METHYLATION' in data.get('modifications', [])
        if 'is_pks_hybrid' in data:
            is_pks_hybrid = data['is_pks_hybrid']
        else:
            is_pks_hybrid = False
        return cls(residue=NerpaResidue(data['residue']),
                   methylated=methylated,
                   chirality=Chirality[data['chirality']],
                   is_pks_hybrid=is_pks_hybrid,
                   rban_name=NorineMonomerName(data['rban_name']),
                   rban_idx=int(data['rban_idx']))

    @classmethod
    def from_list(cls, data: list) -> rBAN_Monomer:
        return cls(residue=NerpaResidue(data[0]),
                   methylated=data[1],
                   chirality=Chirality[data[2]],
                   is_pks_hybrid=data[3],
                   rban_name=NorineMonomerName(data[4]),
                   rban_idx=int(data[5]))

    def to_base_mon(self) -> NRP_Monomer:
        return NRP_Monomer(residue=self.residue,
                           methylated=self.methylated,
                           chirality=self.chirality,
                           is_pks_hybrid=self.is_pks_hybrid)



rBAN_MONOMER_DUMMY = rBAN_Monomer(residue=NerpaResidue(''),
                                  methylated=False,
                                  chirality=Chirality.UNKNOWN,
                                  is_pks_hybrid=False,
                                  rban_name=NorineMonomerName(''),
                                  rban_idx=-1)
