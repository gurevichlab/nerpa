from __future__ import annotations
from typing import (
    Dict,
    List,
    NewType
)

from src.generic.graphs import get_boundary_edges, get_interior_edges
from src.nerpa_ms.monomer_graph.bonds import (
    BindingSiteType,
    BindingSitesFgpt,
    get_binding_sites
)
from src.monomer_names_helper import NorineMonomerName
from src.rban_parsing.rban_parser import (
    AtomId,
    Chirality,
    MonomerIdx,
    Parsed_rBAN_Record,
)

import networkx as nx
from dataclasses import dataclass


AtomicGraph = NewType('AtomicGraph', nx.Graph)
Monomer_JSON_Dict = NewType('Monomer_JSON_Dict', dict)  # custom json representation for Monomer


def map_preserves_binding_sites(mon1: Monomer,
                                mon2: Monomer,
                                fst_to_snd_map: Dict[AtomId, AtomId]) -> bool:
    return all(sorted([fst_to_snd_map[u] for u in mon1.binding_sites[binding_site_type]])
               ==
               sorted(mon2.binding_sites[binding_site_type])
               for binding_site_type in mon1.binding_sites.keys())


@dataclass
class Monomer:
    name: NorineMonomerName
    atomic_graph: AtomicGraph
    chirality: Chirality
    binding_sites: Dict[BindingSiteType, List[AtomId]]

    def binding_sites_fgpt(self) -> BindingSitesFgpt:
        return BindingSitesFgpt({binding_site_type: len(self.binding_sites[binding_site_type])
                                 for binding_site_type in self.binding_sites.keys()})

    def __eq__(self, other):
        if any([len(self.atomic_graph.nodes) != len(other.atomic_graph.nodes),
                len(self.atomic_graph.edges) != len(other.atomic_graph.edges),
                self.binding_sites_fgpt() != other.binding_sites_fgpt()]):
            return False

        # check that there exists an isomorphism which maps binding sites to each other
        return any(map_preserves_binding_sites(self, other, mapping)
                   for mapping in nx.vf2pp_all_isomorphisms(self.atomic_graph,
                                                            other.atomic_graph,
                                                            node_label='name'))

    # retrieve monomer with index rban_idx from the rban_record
    @classmethod
    def from_rban(cls,
                  rban_idx: MonomerIdx,
                  rban_record: Parsed_rBAN_Record) -> Monomer:

        monomer_atoms_idxs = rban_record.monomers[rban_idx].atoms
        monomer_edges = get_interior_edges(monomer_atoms_idxs,
                                           rban_record.atomic_bonds.keys())

        atomic_graph = nx.Graph()
        atomic_graph.add_nodes_from([(idx, rban_record.atoms[idx]._asdict())
                                     for idx in monomer_atoms_idxs])
        atomic_graph.add_edges_from((u, v, rban_record.atomic_bonds[(u, v)]._asdict())
                                    for u, v in monomer_edges)

        boundary_edges = {e: rban_record.atomic_bonds[e]
                          for e in get_boundary_edges(monomer_atoms_idxs,
                                                      rban_record.atomic_bonds.keys())}
        binding_sites = get_binding_sites(monomer_atoms_idxs,
                                          boundary_edges,
                                          rban_record.atoms)
        return cls(name=rban_record.monomers[rban_idx].name,
                   atomic_graph=AtomicGraph(atomic_graph),
                   binding_sites=binding_sites,
                   chirality=rban_record.monomers[rban_idx].chirality)

    @classmethod
    def from_json_dict(cls, record: Monomer_JSON_Dict) -> Monomer:
        atomic_graph = nx.Graph()
        atomic_graph.add_nodes_from([(AtomId(int(node_id)), node_attr)
                                     for node_id, node_attr in record['atomic_graph']['nodes'].items()])
        atomic_graph.add_edges_from([(AtomId(edge['endpoints'][0]),
                                      AtomId(edge['endpoints'][1]),
                                      edge['attributes'])
                                     for edge in record['atomic_graph']['edges']])
        return cls(name=record['name'],
                   atomic_graph=AtomicGraph(atomic_graph),
                   binding_sites={BindingSiteType[bs_type_str]: binding_sites
                                  for bs_type_str, binding_sites in record['binding_sites'].items()},
                   chirality=Chirality[record['chirality']])

    def to_json_dict(self) -> Monomer_JSON_Dict:
        # test that from_json(to_json(self)) == self
        nodes = {atom_id: atom_info
                 for atom_id, atom_info in self.atomic_graph.nodes.items()}
        edges = [{'endpoints': edge_endpts,
                  'attributes': edge_attr}
                 for edge_endpts, edge_attr in self.atomic_graph.edges.items()]

        binding_sites_dict = {bs_type.name: atoms
                              for bs_type, atoms in self.binding_sites.items()}
        json_dict = Monomer_JSON_Dict({'name': self.name,
                                       'binding_sites': binding_sites_dict,
                                       'atomic_graph': {'nodes': nodes,
                                                        'edges': edges},
                                       'chirality': self.chirality.name})
        assert self.from_json_dict(json_dict) == self  # for testing
        return json_dict


def get_binding_sites_map(mon1: Monomer,
                          mon2: Monomer) -> Dict[AtomId, AtomId]:  # fst to snd
    assert mon1.binding_sites_fgpt() == mon2.binding_sites_fgpt()
    return {mon1_atom_idx: mon2_atom_idx
            for bs_type in mon1.binding_sites
            for mon1_atom_idx, mon2_atom_idx in
            zip(mon1.binding_sites[bs_type],
                mon2.binding_sites[bs_type])}
