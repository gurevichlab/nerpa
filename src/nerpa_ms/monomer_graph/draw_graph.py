from typing import Dict, Tuple
from itertools import product, chain
from collections import defaultdict
from pathlib import Path
from rdkit.Chem.Draw import rdMolDraw2D
from src.nerpa_ms.monomer_graph.monomer_graph import MonomerGraph
from src.nerpa_ms.monomer_graph.chem_helper import (
    AtomicGraphRecord,
    atomic_graph_to_chem,
    to_atomic_graph
)
import networkx as nx
import matplotlib.pyplot as plt


def monomer_color_dict() -> Dict[str, Tuple[float, float, float]]:
    standard_aas = ['ala', 'arg', 'asn', 'asp', 'cys',
                    'gln', 'glu', 'gly', 'his', 'ile',
                    'leu', 'lys', 'met', 'phe', 'pro',
                    'ser', 'thr', 'trp', 'tyr', 'val']
    colors_pool = [(0.7 + r/10, 0.7 + g/10, 0.7 + b/10)
                   for r, g, b in product(range(4), repeat=3)
                   if not r == g == b]
    not_supported_color = (0.8, 0.8, 0.8)
    aa_color = defaultdict(lambda: not_supported_color,
                           dict(zip(standard_aas, colors_pool)))
    return aa_color


def draw_molecule(G: MonomerGraph, output_file: Path):
    aa_color = monomer_color_dict()

    atomic_graph_atoms, atomic_graph_edges = to_atomic_graph(G)
    mol, atom_id_to_index = atomic_graph_to_chem(AtomicGraphRecord(atomic_graph_atoms, atomic_graph_edges))

    atom_colors = dict()
    bonds_to_highlight = list()
    atom_labels = dict()
    for monomer_id in G:
        monomer = G.nodes[monomer_id]['monomer']
        monomer_color = aa_color[monomer.name.lower()]
        # monomer_color = aa_color[monomer.name.lower()] if monomer_id not in G.core_nodes\
        #    else core_node_color

        monomer_repr = next((atom_id_to_index[(monomer_id, atom_id)]  # atom where monomers name will be displayed
                             for atom_id in monomer.atomic_graph
                             if monomer.atomic_graph.nodes[atom_id]['name'] == 'C' and \
                             atom_id not in chain(*monomer.binding_sites.values())),
                            None)
    
        if monomer_repr:
            atom_labels[monomer_repr] = f'{monomer.name}_{monomer_id}'

        for atom_id in monomer.atomic_graph:  # set atoms colors
            atom_index = atom_id_to_index[(monomer_id, atom_id)]
            atom_colors[atom_index] = monomer_color

        for atom1_id, atom2_id in monomer.atomic_graph.edges:  # set bonds to highlight (color is deduced automatically)
            index1, index2 = [atom_id_to_index[(monomer_id, atom_id)]
                              for atom_id in (atom1_id, atom2_id)]
            bond_idx = mol.GetBondBetweenAtoms(index1, index2).GetIdx()
            bonds_to_highlight.append(bond_idx)

        drawer = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
        
        opts = drawer.drawOptions()
        for atom_index, label in atom_labels.items():  # write names of monomers
            opts.atomLabels[atom_index] = label
        drawer.DrawMolecule(mol,
                            highlightAtoms=atom_colors.keys(),
                            highlightAtomColors=atom_colors,
                            highlightBonds=bonds_to_highlight)

        drawer.FinishDrawing()
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_bytes(drawer.GetDrawingText())


def make_graph_to_draw(G: MonomerGraph) -> nx.DiGraph:
    H = nx.DiGraph()
    H.add_nodes_from([(u, {'name': f'{G.nodes[u]["monomer"].name}_{u}'}) for u in G.nodes])
    for u, v in G.edges:
        edge_attr = G.edges[u, v]
        if edge_attr['bondType'] == 'AMINO':
            u_atom_idx = edge_attr['monomer_to_atom'][u]
            u_atom = G.nodes[u]['monomer'].atomic_graph.nodes[u_atom_idx]['name']
            if u_atom == 'C':
                H.add_edge(u, v)
            else:
                H.add_edge(v, u)
        else:
            H.add_edge(u, v)
            H.add_edge(v, u)

    return H


def draw_monomer_graph(G: MonomerGraph, output_file: Path,
                       with_rban_indexes: bool=True):
    aa_color = monomer_color_dict()
    H = make_graph_to_draw(G)
    pos = nx.nx_agraph.graphviz_layout(H)
    edge_colors = []
    for u, v in H.edges:
        if (v, u) in H.edges:
            edge_colors.append('red')
        else:
            edge_colors.append('blue')
    def get_residue_name(u):
        full_name = H.nodes[u]['name']
        full_name_wo_index = full_name.split('_')[0]
        return full_name_wo_index.split('-')[-1].lower()

    def get_node_name(u):
        full_name = H.nodes[u]['name']
        if with_rban_indexes:
            return full_name
        else:
            return full_name.split('_')[0]

    nx.draw_networkx(H,
                     pos=pos,
                     node_color=[aa_color[get_residue_name(u)] for u in H.nodes],
                     edge_color=edge_colors,
                     arrowsize=20,
                     labels={u: get_node_name(u) for u in H.nodes},
                     font_size=8,
                     node_size=500,
                     width=3)
    output_file.parent.mkdir(exist_ok=True, parents=True)
    plt.savefig(output_file)
    plt.clf()

