from __future__ import annotations
from typing import NewType, Tuple

from src.nerpa_ms.monomer_graph.bonds import AMINO2_FGPT, BindingSiteType
from src.nerpa_ms.monomer_graph.monomer import Monomer, get_binding_sites_map
from src.rban_parsing.rban_parser import Parsed_rBAN_Record

from networkx import Graph
from copy import deepcopy

MonomerId = int  # NewType('MonomerId', int)
Edge = NewType('Edge', Tuple[MonomerId, MonomerId])


class MonomerGraph(Graph):
    def __init__(self, incoming_graph_data=None):
        super().__init__(incoming_graph_data)

    def monomer(self, mon_id: MonomerId) -> Monomer:
        return self.nodes[mon_id]['monomer']

    @classmethod
    def from_rban_record(cls, rban_record: Parsed_rBAN_Record) -> MonomerGraph:
        monomer_graph = cls()  # empty Graph
        monomer_graph.add_nodes_from([(idx, {'monomer': Monomer.from_rban(idx, rban_record)})
                                      for idx in rban_record.monomers])
        monomer_graph.add_edges_from([(u, v, edge_attr._asdict())
                                      for (u, v), edge_attr in rban_record.monomer_bonds.items()])
        return monomer_graph

    def can_substitute_monomer(self, mon_id: MonomerId,
                               new_monomer: Monomer) -> bool:
        return mon_id in self.nodes and \
            self.monomer(mon_id).binding_sites_fgpt() == new_monomer.binding_sites_fgpt()

    # for now I consider only one variant of binding site mapping
    def substitute_monomer(self,
                           monomer_idx: MonomerId,
                           new_monomer: Monomer):
        binding_site_map = get_binding_sites_map(self.monomer(monomer_idx),
                                                 new_monomer)  # old to new

        self.nodes[monomer_idx]['monomer'] = new_monomer
        for other_monomer_idx, edge_attr in self.adj[monomer_idx].items():
            old_atom = edge_attr['monomer_to_atom'].pop(monomer_idx)
            edge_attr['monomer_to_atom'][monomer_idx] = binding_site_map[old_atom]

    def can_remove_monomer(self, mon_id: MonomerId) -> bool:
        return mon_id in self.nodes and self.monomer(mon_id).binding_sites_fgpt() == AMINO2_FGPT

    def remove_monomer(self,
                       monomer_idx: MonomerId):
        if not self.can_remove_monomer(monomer_idx):
            raise ValueError('Invalid arguments')

        (mon1, atom1), (mon2, atom2) = [(neighb_monomer, edge_attr['monomer_to_atom'][neighb_monomer])
                                        for neighb_monomer, edge_attr in self.adj[monomer_idx].items()]

        self.remove_node(monomer_idx)
        self.add_edge(mon1, mon2,
                      monomer_to_atom={mon1: atom1, mon2: atom2},
                      bondType='AMINO',
                      arity=1)

    def can_insert_monomer(self,
                           e: Edge,
                           new_monomer: Monomer) -> bool:
        return e in self.edges and \
            self.edges[e]['bondType'] == 'AMINO' and \
            new_monomer.binding_sites_fgpt() == AMINO2_FGPT

    def insert_monomer(self,
                       e: Edge,
                       new_monomer: Monomer) -> MonomerId:  # returns the index of inserted monomer
        if not self.can_insert_monomer(e, new_monomer):
            raise ValueError('Invalid arguments')

        new_idx = max(self.nodes) + 1
        self.add_node(new_idx, monomer=new_monomer)

        new_mon_c_site = new_monomer.binding_sites[BindingSiteType.AMINO_C][0]
        new_mon_n_site = new_monomer.binding_sites[BindingSiteType.AMINO_N][0]

        for neighb_mon_idx in e:
            neighb_mon = self.monomer(neighb_mon_idx)
            neighb_mon_atom = self.edges[e]['monomer_to_atom'][neighb_mon_idx]
            new_mon_atom = new_mon_n_site \
                if neighb_mon_atom in neighb_mon.binding_sites[BindingSiteType.AMINO_C] \
                else new_mon_c_site

            self.add_edge(new_idx, neighb_mon_idx,
                          monomer_to_atom={new_idx: new_mon_atom,
                                           neighb_mon_idx: neighb_mon_atom},
                          bondType='AMINO',
                          arity=1)

        self.remove_edge(*e)
        return new_idx

    def to_json_dict(self) -> dict:
        # test that from_json(to_json(self)) == self
        nodes = {monomer_id: self.monomer(monomer_id).to_json_dict()
                 for monomer_id in self}
        edges = [{'endpoints': edge_endpts,
                  'attributes': edge_attr}
                 for edge_endpts, edge_attr in self.edges.items()]
        json_dict = {'nodes': nodes,
                     'edges': edges}
        return json_dict

    @classmethod
    def from_json_dict(cls, record):
        graph = Graph()
        graph.add_nodes_from([(MonomerId(int(node_id)),
                               {'monomer': Monomer.from_json_dict(node_attr)})
                              for node_id, node_attr in record['nodes'].items()])

        def parsed_attributes(edge_attributes):
            new_attributes = deepcopy(edge_attributes)
            new_attributes['monomer_to_atom'] = dict((int(monomer_id), atom_id)
                                                     for monomer_id, atom_id in edge_attributes['monomer_to_atom'].items())
            return new_attributes

        graph.add_edges_from([(edge['endpoints'][0], edge['endpoints'][1], parsed_attributes(edge['attributes']))
                              for edge in record['edges']])
        return cls(graph)
