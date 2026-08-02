from __future__ import annotations

from itertools import permutations
from typing import (
    List,
    Optional,
    NamedTuple,
    Dict,
    Callable,
    Iterable,
    Hashable,
    Tuple,
)

import networkx as nx

from src.monomer_names_helper import NRP_Monomer
from src.rban_parsing.rban_monomer import rBAN_Monomer
from dataclasses import dataclass, asdict
from src.general_type_aliases import SMILES
from src.rban_parsing.rban_parser import Parsed_rBAN_Record, NRP_metadata

NRP_ID = str

def is_mibig_nrp(nrp_id: NRP_ID) -> bool:
    return nrp_id.startswith('BGC')

def is_norine_nrp(nrp_id: NRP_ID) -> bool:
    return nrp_id.startswith('NOR')

class NRP_Variant_ID(NamedTuple):
    nrp_id: NRP_ID
    variant_idx: int = 0


@dataclass
class NRP_Fragment:
    monomers: List[rBAN_Monomer]
    is_cyclic: bool

    def default_sort_key(self) -> Tuple[bool, List[NRP_Monomer]]:
        return (
            self.is_cyclic,
            [mon.to_base_mon() for mon in self.monomers]
        )

    def _canonize(self):
        if self.is_cyclic:
            # For cyclic fragments, choose the lexicographically smallest rotation of the monomer list
            self.monomers = min(
                (
                    self.monomers[i:] + self.monomers[:i]
                    for i in range(len(self.monomers))
                ),
                key=lambda mons: [
                    mon.to_base_mon() for mon in mons
                ]
            )

    def __init__(self, monomers: List[rBAN_Monomer], is_cyclic: bool):
        self.monomers = monomers
        self.is_cyclic = is_cyclic
        self._canonize()

    def to_dict(self) -> Dict:
        return {
            'monomers': [mon.to_dict() for mon in self.monomers],
            'is_cyclic': self.is_cyclic,
        }

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Fragment:
        return cls(
            is_cyclic=data['is_cyclic'],
            monomers=[
                rBAN_Monomer.from_yaml_dict(mon_data)
                for mon_data in data['monomers']
            ])

    def is_isomorphic_to(
            self,
            other: NRP_Fragment,
            node_key: Optional[Callable[[rBAN_Monomer], Hashable]] = None,
    ) -> bool:
        if not isinstance(other, NRP_Fragment):
            raise NotImplemented

        if len(self.monomers) != len(other.monomers) or self.is_cyclic != other.is_cyclic:
            return False

        if node_key is None:
            node_key = lambda mon: mon.to_base_mon()
            default_key = True
        else:
            default_key = False
            
        self_keys = [node_key(mon) for mon in self.monomers]
        other_keys = [node_key(mon) for mon in other.monomers]

        if not self.is_cyclic or default_key:
            # For default key, we can rely on canonization to check equality
            return self_keys == other_keys
        else:
            # For non-default keys, we need to check all rotations for cyclic fragments
            return any(
                self_keys == other_keys[i:] + other_keys[:i]
                for i in range(len(other_keys))
            )



@dataclass
class NRP_Variant:
    nrp_variant_id: NRP_Variant_ID
    fragments: List[NRP_Fragment]
    isolated_unknown_monomers: List[rBAN_Monomer]
    metadata: NRP_metadata

    def _canonize(self):
        self.fragments = sorted(self.fragments, key=NRP_Fragment.default_sort_key)
        self.isolated_unknown_monomers = sorted(self.isolated_unknown_monomers)

    def __init__(
            self,
            nrp_variant_id: NRP_Variant_ID,
            fragments: List[NRP_Fragment],
            isolated_unknown_monomers: List[rBAN_Monomer],
            metadata: NRP_metadata
    ):
        self.nrp_variant_id = nrp_variant_id
        self.fragments = fragments
        self.isolated_unknown_monomers = isolated_unknown_monomers
        self.metadata = metadata
        self._canonize()

    def to_dict(self) -> Dict:
        return {
            'nrp_variant_id': self.nrp_variant_id._asdict(),
            'fragments': [frag.to_dict() for frag in self.fragments],
            'isolated_monomers': [mon.to_dict() for mon in self.isolated_unknown_monomers],
            'metadata': asdict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict) -> NRP_Variant:
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

        return cls(
            nrp_variant_id=nrp_variant_id,
            fragments=list(map(NRP_Fragment.from_yaml_dict, data['fragments'])),
            isolated_unknown_monomers=[
                rBAN_Monomer.from_yaml_dict(mon_data)
                for mon_data in data.get('isolated_monomers', [])
            ],
            metadata=metadata
        )

    def is_isomorphic_to(
            self,
            other,
            node_key: Optional[Callable[[rBAN_Monomer], Hashable]] = None,
    ) -> bool:
        if not isinstance(other, NRP_Variant):
            raise NotImplemented

        if node_key is None:
            node_key = lambda mon: mon.to_base_mon()
            default_key = True
        else:
            default_key = False

        # heuristic: compare sets of keys first to quickly rule out non-isomorphic variants
        self_keys = {node_key(mon) for frag in self.fragments for mon in frag.monomers}
        other_keys = {node_key(mon) for frag in other.fragments for mon in frag.monomers}
        if self_keys != other_keys:
            return False

        def fragment_lists_equal(
                list1: Iterable[NRP_Fragment],
                list2: Iterable[NRP_Fragment],
                node_key: Optional[Callable[[rBAN_Monomer], Hashable]] = None,
        ) -> bool:
            return all(
                frag1.is_isomorphic_to(frag2, node_key=node_key)
                for frag1, frag2 in zip(list1, list2)
            )

        if len(self.fragments) != len(other.fragments):
            return False

        if default_key:
            # It's important to pass node_key=None so that Fragment.is_isomorphic_to can rely on canonization
            return fragment_lists_equal(self.fragments, other.fragments)
        else:
            return any(
                fragment_lists_equal(
                    self.fragments,
                    permuted_other_frags,
                    node_key=node_key
                )
                for permuted_other_frags in permutations(other.fragments)
            )

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
        return [
            nrp
            for nrp in self.nrp_variants
            if nrp.nrp_variant_id in repr_ids
        ]
