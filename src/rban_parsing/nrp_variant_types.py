from __future__ import annotations

from itertools import permutations
from typing import List, Optional, NamedTuple, Dict
from src.rban_parsing.rban_monomer import rBAN_Monomer
from dataclasses import dataclass
from src.general_type_aliases import SMILES
from src.rban_parsing.rban_parser import Parsed_rBAN_Record, NRP_metadata


class NRP_Variant_ID(NamedTuple):
    nrp_id: str
    variant_idx: int

    def to_dict(self) -> dict:
        return {'nrp_id': self.nrp_id,
                'variant_idx': self.variant_idx}


@dataclass
class NRP_Fragment:
    monomers: List[rBAN_Monomer]
    is_cyclic: bool

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Fragment:
        return cls(is_cyclic=data['is_cyclic'],
                   monomers=[rBAN_Monomer.from_list(mon_data_lst)
                             for mon_data_lst in data['monomers']])

    def __eq__(self, other):
        if not isinstance(other, NRP_Fragment):
            return NotImplemented
        if len(self.monomers) != len(other.monomers) or self.is_cyclic != other.is_cyclic:
            return False
        mons1 = [mon.to_base_mon() for mon in self.monomers]
        mons2 = [mon.to_base_mon() for mon in other.monomers]
        if self.is_cyclic:
            # Check all cyclic permutations
            return any(mons1 == (mons2[i:] + mons2[:i])
                       for i in range(len(mons2)))
        else:
            return mons1 == mons2




@dataclass
class NRP_Variant:
    nrp_variant_id: NRP_Variant_ID
    fragments: List[NRP_Fragment]
    isolated_unknown_monomers: List[rBAN_Monomer]
    metadata: NRP_metadata

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Variant:
        if 'metadata' in data:
            metadata = NRP_metadata.from_dict(data['metadata'])
        else:
            metadata = NRP_metadata(name=None, smiles=None, origin=None, inchikey=None, source=None)
        return cls(nrp_variant_id=NRP_Variant_ID(*data['nrp_variant_id']),
                   fragments=list(map(NRP_Fragment.from_yaml_dict, data['fragments'])),
                   isolated_unknown_monomers=list(map(rBAN_Monomer.from_yaml_dict, data['isolated_monomers']))
                   if data.get('isolated_monomers') else [],
                   metadata=metadata)

    def __eq__(self, other):
        if not isinstance(other, NRP_Variant):
            return NotImplemented
        if len(self.fragments) != len(other.fragments):
            return False

        return any(self.fragments == list(permuted_other_frags)
                   for permuted_other_frags in permutations(other.fragments))


class NRP_Variants_Info(NamedTuple):
    nrp_variants: List[NRP_Variant]
    rban_records: List[Parsed_rBAN_Record]
    nrp_id_to_repr_id: Dict[NRP_Variant_ID, NRP_Variant_ID]

    def get_representative_nrp_variants(self) -> List[NRP_Variant]:
        repr_ids = set(self.nrp_id_to_repr_id.values())
        return [nrp for nrp in self.nrp_variants if nrp.nrp_variant_id in repr_ids]
