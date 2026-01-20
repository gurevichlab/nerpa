from __future__ import annotations

from itertools import permutations
from typing import List, Optional, NamedTuple, Dict, Callable, Iterable

import networkx as nx

from src.monomer_names_helper import NRP_Monomer
from src.rban_parsing.rban_monomer import rBAN_Monomer
from dataclasses import dataclass, asdict
from src.general_type_aliases import SMILES
from src.rban_parsing.rban_parser import Parsed_rBAN_Record, NRP_metadata


class NRP_Variant_ID(NamedTuple):
    nrp_id: str
    variant_idx: int


@dataclass
class NRP_Fragment:
    monomers: List[rBAN_Monomer]
    is_cyclic: bool

    def to_dict(self) -> Dict:
        return {'monomers': [mon.to_dict() for mon in self.monomers],
                'is_cyclic': self.is_cyclic,}

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Fragment:
        return cls(is_cyclic=data['is_cyclic'],
                   monomers=[rBAN_Monomer.from_yaml_dict(mon_data)
                             for mon_data in data['monomers']])

    def is_isomorphic_to(self,
                         other: NRP_Fragment,
                         monomers_comparator: Callable[[NRP_Monomer, NRP_Monomer], bool] = None) -> bool:
        if not isinstance(other, NRP_Fragment):
            return NotImplemented
        if len(self.monomers) != len(other.monomers) or self.is_cyclic != other.is_cyclic:
            return False

        if monomers_comparator is None:
            monomers_comparator = lambda m1, m2: m1 == m2

        mons1 = [mon.to_base_mon() for mon in self.monomers]
        mons2 = [mon.to_base_mon() for mon in other.monomers]

        def mon_lists_equal(list1: List[NRP_Monomer], list2: List[NRP_Monomer]) -> bool:
            return all(monomers_comparator(mon1, mon2) for mon1, mon2 in zip(list1, list2))

        if self.is_cyclic:
            # Check all cyclic permutations
            return any(mon_lists_equal(mons1, mons2[i:] + mons2[:i])
                       for i in range(len(mons2)))
        else:
            return mon_lists_equal(mons1, mons2)

    def __eq__(self, other):
        return self.is_isomorphic_to(other)


@dataclass
class NRP_Variant:
    nrp_variant_id: NRP_Variant_ID
    fragments: List[NRP_Fragment]
    isolated_unknown_monomers: List[rBAN_Monomer]
    metadata: NRP_metadata

    def to_dict(self) -> Dict:
        return {'nrp_variant_id': self.nrp_variant_id._asdict(),
                'fragments': [frag.to_dict() for frag in self.fragments],
                'isolated_monomers': [mon.to_dict() for mon in self.isolated_unknown_monomers],
                'metadata': asdict(self.metadata),}

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Variant:
        if 'metadata' in data:
            metadata = NRP_metadata.from_dict(data['metadata'])
        else:
            metadata = NRP_metadata(name=None, smiles=None, origin=None, inchikey=None, source=None)

        if isinstance(data['nrp_variant_id'], dict):
            nrp_variant_id = NRP_Variant_ID(**data['nrp_variant_id'])
        elif isinstance(data['nrp_variant_id'], list):
            nrp_variant_id = NRP_Variant_ID(*data['nrp_variant_id'])
        else:
            raise ValueError(f"Unexpected format for nrp_variant_id: {data['nrp_variant_id']}")

        return cls(nrp_variant_id=nrp_variant_id,
                   fragments=list(map(NRP_Fragment.from_yaml_dict, data['fragments'])),
                   isolated_unknown_monomers=[rBAN_Monomer.from_yaml_dict(mon_data)
                                              for mon_data in data.get('isolated_monomers', [])],
                   metadata=metadata)

    def is_isomorphic_to(self, other,
                         monomers_comparator: Callable[[NRP_Monomer, NRP_Monomer], bool] = None) -> bool:
        if not isinstance(other, NRP_Variant):
            return NotImplemented
        if len(self.fragments) != len(other.fragments):
            return False

        def fragment_lists_equal(list1: Iterable[NRP_Fragment],
                                 list2: Iterable[NRP_Fragment]) -> bool:
            return all(frag1.is_isomorphic_to(frag2, monomers_comparator)
                       for frag1, frag2 in zip(list1, list2))

        return any(fragment_lists_equal(self.fragments, permuted_other_frags)
                   for permuted_other_frags in permutations(other.fragments))

    def __eq__(self, other):
        return self.is_isomorphic_to(other)

    def to_str_compact(self) -> str:
        frag_strs = []
        for frag in self.fragments:
            mons_str = ', '.join(mon.rban_name for mon in frag.monomers)
            if frag.is_cyclic:
                mons_str = f'({mons_str})'
            frag_strs.append(mons_str)
        return '\n'.join([self.nrp_variant_id.nrp_id] + frag_strs)

    def to_nx_digraph(self, node_label_key: str = 'monomer') -> nx.DiGraph:
        G = nx.DiGraph()
        for frag_idx, fragment in enumerate(self.fragments):
            for mon_idx, monomer in enumerate(fragment.monomers):
                node_id = (frag_idx, mon_idx)
                G.add_node(node_id, **{node_label_key: monomer})
                if mon_idx > 0:
                    G.add_edge((frag_idx, mon_idx - 1), node_id)
            if fragment.is_cyclic and len(fragment.monomers) > 1:
                G.add_edge((frag_idx, len(fragment.monomers) - 1), (frag_idx, 0))
        return G


class NRP_Variants_Info(NamedTuple):
    nrp_variants: List[NRP_Variant]
    rban_records: List[Parsed_rBAN_Record]
    nrp_id_to_repr_id: Dict[NRP_Variant_ID, NRP_Variant_ID]

    def get_representative_nrp_variants(self) -> List[NRP_Variant]:
        repr_ids = set(self.nrp_id_to_repr_id.values())
        return [nrp for nrp in self.nrp_variants if nrp.nrp_variant_id in repr_ids]
