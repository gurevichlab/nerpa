import io
from typing import Dict, Tuple, List, NewType, Optional
from itertools import product, chain
from collections import defaultdict
from pathlib import Path

import graphviz
from rdkit.Chem.Draw import rdMolDraw2D

from src.monomer_names_helper import MonomerNamesHelper
from src.nerpa_ms.monomer_graph.monomer_graph import MonomerGraph
from src.nerpa_ms.monomer_graph.chem_helper import (
    AtomicGraphRecord,
    atomic_graph_to_chem,
    to_atomic_graph
)
import networkx as nx
import matplotlib.pyplot as plt
import colorsys
import cairosvg
from PIL import Image


RGB = NewType('RGB', Tuple[float, float, float])
ColorHex = NewType('ColorHex', str)

def get_n_distinct_rgb_colors(N: int,
                              saturation: float = 0.3,
                              value = 0.8) -> List[RGB]:
    """
    Returns a list of N RGB triplets (floats in [0,1]) that are
    maximally distinct by evenly stepping through the hue wheel.
    """
    colors = []
    for i in range(N):
        hue = i / N
        rgb = RGB(colorsys.hsv_to_rgb(hue, saturation, value))
        colors.append(rgb)
    return colors


def rgb_to_hex(rgb: RGB) -> ColorHex:
    """(r, g, b) with 0≤r,g,b≤1 → “#rrggbb”"""
    return ColorHex("#{:02x}{:02x}{:02x}".format(
        int(rgb[0] * 255),
        int(rgb[1] * 255),
        int(rgb[2] * 255),
    ))


STANDARD_AAS = ['ala', 'arg', 'asn', 'asp', 'cys',
                'gln', 'glu', 'gly', 'his', 'ile',
                'leu', 'lys', 'met', 'phe', 'pro',
                'ser', 'thr', 'trp', 'tyr', 'val']

def monomer_color_dict(monomers: Optional[List[str]] = None,
                       not_supported_color: RGB = (0.8, 0.8, 0.8)) -> Dict[str, RGB]:
    if monomers is None:
        monomers = STANDARD_AAS
    colors = get_n_distinct_rgb_colors(len(monomers))
    return defaultdict(lambda: not_supported_color,
                       dict(zip(monomers, colors)))


def draw_molecule(G: MonomerGraph, output_file: Path):
    aa_color = monomer_color_dict()

    atomic_graph_atoms, atomic_graph_edges = to_atomic_graph(G)
    mol, atom_id_to_index = atomic_graph_to_chem(AtomicGraphRecord(atomic_graph_atoms, atomic_graph_edges))

    atom_colors = dict()
    bonds_to_highlight = list()
    atom_labels = dict()
    for monomer_id in G:
        monomer = G.nodes[monomer_id]['monomer']
        residue = monomer.name.split('-')[-1].lower()  # e.g. 'ala' from 'NMe-Ala'
        monomer_color = aa_color[residue]
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


def draw_monomer_graph(G: MonomerGraph,
                       output_path: Path,
                       with_rban_indexes: bool = True,
                       format: str = 'png',
                       size: Tuple[int, int] = (1000, 1000),
                       dpi: int = 300) -> graphviz.Digraph:
    if format != 'png':
        raise ValueError(f"Unsupported format: {format}. Currently only 'png' is supported.")
    
    aa_color = monomer_color_dict()
    H = make_graph_to_draw(G)

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

    w_inches, h_inches = size[0] / dpi, size[1] / dpi
    fig = graphviz.Digraph(format=format,
                           engine='neato',
                           graph_attr={
                               'splines': 'true',  # tell Graphviz to eliminate overlaps/crossings where possible
                               'overlap': 'false',  # splines='true' gives you smooth curved edges instead of straight lines
                               'size': f'{w_inches},{h_inches}!',
                               'ratio': 'fill',  # fill the whole image
                               'margin': '0',
                               'dpi': str(dpi),
                           })

    # make every edge use a 2‑point pen and 1.5× bigger arrowheads
    fig.attr('edge', penwidth='2', arrowsize='1.5')

    for u in H.nodes:
        color = rgb_to_hex(aa_color[get_residue_name(u)])
        fig.node(str(u),
                 label=get_node_name(u),
                 color=color,
                 style='filled')

    for u, v in H.edges:
        if (v, u) in H.edges:  # undirected edge labels a non-peptide bond
            if u < v: # to avoid duplicates
                fig.edge(str(u), str(v), color='red', dir='none', arrowhead='none')
        else:  # directed edge --> peptide bond
            fig.edge(str(u), str(v), color='blue')

    if output_path:
        output_path.parent.mkdir(exist_ok=True, parents=True)

        if output_path.name.endswith('.png'):
            fname = str(output_path).rsplit('.', 1)[0]
        else:
            fname = str(output_path)

        # some hacks to ensure the image is saved in the correct size
        # 1) render to a temporary PNG bytestring
        png_bytes = fig.pipe(format='png')
        # 2) open the bytestring as an image
        img = Image.open(io.BytesIO(png_bytes))

        # 3) create a blank canvas of exactly the target size
        w_tgt, h_tgt = size
        canvas = Image.new('RGBA', (w_tgt, h_tgt), (255, 255, 255, 255))

        # 4) compute offsets to center
        w_cur, h_cur = img.size
        x_off = max((w_tgt - w_cur) // 2, 0)
        y_off = max((h_tgt - h_cur) // 2, 0)

        # 5) paste the graph onto the canvas
        canvas.paste(img, (x_off, y_off))

        # 6) save out
        output_path.parent.mkdir(exist_ok=True, parents=True)
        canvas.save(str(output_path), format='PNG')

    return fig

