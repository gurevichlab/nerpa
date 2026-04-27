import io
from typing import Dict, Tuple, List, NewType, Optional, Literal
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

import re
import json

def force_svg_pixel_size(svg_bytes: bytes, width_px: int, height_px: int,
                         stretch: bool = False) -> bytes:
    """
    Force SVG to have exact pixel width/height while preserving Graphviz output.
    This edits only the opening <svg ...> tag and avoids XML reserialization.
    """
    svg = svg_bytes.decode("utf-8", errors="strict")

    # Find opening <svg ...> tag
    m = re.search(r"<svg\b[^>]*>", svg, flags=re.IGNORECASE)
    if not m:
        return svg_bytes  # not an SVG? leave it alone

    tag = m.group(0)

    # Replace or insert width/height
    def upsert_attr(tag: str, name: str, value: str) -> str:
        if re.search(rf"\b{name}\s*=", tag):
            return re.sub(rf'\b{name}\s*=\s*"[^"]*"', f'{name}="{value}"', tag)
        else:
            # insert before the closing '>'
            return tag[:-1] + f' {name}="{value}">'

    tag2 = tag
    tag2 = upsert_attr(tag2, "width", f"{width_px}px")
    tag2 = upsert_attr(tag2, "height", f"{height_px}px")

    par = "none" if stretch else "xMidYMid meet"
    tag2 = upsert_attr(tag2, "preserveAspectRatio", par)

    # Replace the tag in the document
    svg2 = svg[:m.start()] + tag2 + svg[m.end():]
    return svg2.encode("utf-8")


def ensure_image_ext(p: Path, fmt: Literal['png', 'svg']) -> Path:
    """Ensure path ends with .svg or .png without stripping mid-name parts like '.1'."""
    s = str(p)
    fmt = fmt.lower()
    if s.lower().endswith('.png') or s.lower().endswith('.svg'):
        s = s[:s.rfind('.')] + f'.{fmt}'
    else:
        s = s + f'.{fmt}'
    return Path(s)

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
                       not_supported_color: RGB = (0.8, 0.8, 0.8),) -> Dict[str, RGB]:
    if monomers is None:
        monomers = STANDARD_AAS
    colors = get_n_distinct_rgb_colors(len(monomers))
    return defaultdict(lambda: not_supported_color,
                       dict(zip(monomers, colors)))


def draw_molecule(G: MonomerGraph,
                  output_file: Path,
                  rban_indexes: bool = True,
                  monomer_labels: bool = True,
                  size: Tuple[int, int] = (1000, 1000),
                  monomer_names_helper: Optional[MonomerNamesHelper] = None) -> None:
    format = output_file.suffix[1:].lower()
    if format not in ('svg', 'png'):
        raise ValueError(f"Unsupported format: {format}. Use 'svg' or 'png'.")
    residues_with_colors = (
        monomer_names_helper.supported_residues
        if monomer_names_helper else
        STANDARD_AAS
    )
    aa_color = monomer_color_dict(monomers=residues_with_colors)

    atomic_graph_atoms, atomic_graph_edges = to_atomic_graph(G)
    mol, atom_id_to_index = atomic_graph_to_chem(AtomicGraphRecord(atomic_graph_atoms, atomic_graph_edges))

    atom_colors = dict()
    bonds_to_highlight = list()
    atom_labels = dict()
    for monomer_id in G:
        monomer = G.nodes[monomer_id]['monomer']
        if monomer_names_helper:
            residue = monomer_names_helper.parsed_name(name=monomer.name,
                                                       name_format='rBAN/Norine').residue
        else:
            residue = monomer.name.split('-')[-1].lower()  # e.g. 'ala' from 'NMe-Ala'
        monomer_color = aa_color[residue]
        # monomer_color = aa_color[monomer.name.lower()] if monomer_id not in G.core_nodes\
        #    else core_node_color

        monomer_repr = next((atom_id_to_index[(monomer_id, atom_id)]  # atom where monomers name will be displayed
                             for atom_id in monomer.atomic_graph
                             if monomer.atomic_graph.nodes[atom_id]['name'] == 'C' and \
                             atom_id not in chain(*monomer.binding_sites.values())),
                            None)

        if monomer_repr and monomer_labels:
            atom_labels[monomer_repr] = (
                f'{monomer.name}_{monomer_id}'
                if rban_indexes else
                monomer.name
            )

        for atom_id in monomer.atomic_graph:  # set atoms colors
            atom_index = atom_id_to_index[(monomer_id, atom_id)]
            atom_colors[atom_index] = monomer_color

        for atom1_id, atom2_id in monomer.atomic_graph.edges:  # set bonds to highlight (color is deduced automatically)
            index1, index2 = [atom_id_to_index[(monomer_id, atom_id)]
                              for atom_id in (atom1_id, atom2_id)]
            # if monomer_id == 1:
            #     print(f'Highlighting bond between atom {atom1_id} and {atom2_id} in monomer {monomer_id} ({monomer.name})')
            bond_idx = mol.GetBondBetweenAtoms(index1, index2).GetIdx()
            bonds_to_highlight.append(bond_idx)

    w, h = size
    if format == "svg":
        drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    elif format == "png":
        drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    else:
        raise ValueError(f"Unsupported format: {format}. Use 'svg' or 'png'.")

    opts = drawer.drawOptions()
    for atom_index, label in atom_labels.items():  # write names of monomers
        opts.atomLabels[atom_index] = label
    drawer.DrawMolecule(mol,
                        highlightAtoms=atom_colors.keys(),
                        highlightAtomColors=atom_colors,
                        highlightBonds=bonds_to_highlight)

    drawer.FinishDrawing()
    data = drawer.GetDrawingText()
    if isinstance(data, str):  # SVG is text
        data = data.encode("utf-8")

    output_file = ensure_image_ext(output_file, format)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_bytes(data)


def make_graph_to_draw(G: MonomerGraph) -> nx.DiGraph:
    H = nx.DiGraph()
    H.add_nodes_from([(u, {'name': f'{G.nodes[u]["monomer"].name}_{u}'}) for u in G.nodes])
    for u, v in G.edges:
        edge_attr = G.edges[u, v]
        if edge_attr['bond_type'] == 'AMINO':
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
                       size: Tuple[int, int] = (1000, 1000),
                       dpi: int = 300,
                       monomer_names_helper: Optional[MonomerNamesHelper] = None) -> graphviz.Digraph:
    format = output_path.suffix[1:].lower()
    # if format not in ('png',):
    #     # SVG
    #     raise ValueError(f"Unsupported format: {format}. "
    #                      f"Currently only '.png' is supported.")

    residues_with_colors = (
        monomer_names_helper.supported_residues
        if monomer_names_helper else
        STANDARD_AAS
    )
    aa_color = monomer_color_dict(monomers=residues_with_colors)
    H = make_graph_to_draw(G)

    def get_residue_name(u):
        full_name = H.nodes[u]['name']
        full_name_wo_index = full_name.split('_')[0]
        if monomer_names_helper:
            parsed = monomer_names_helper.parsed_name(name=full_name_wo_index,
                                                      name_format='rBAN/Norine')
            return parsed.residue
        return full_name_wo_index.split('-')[-1].lower()

    def get_node_name(u):
        full_name = H.nodes[u]['name']
        if with_rban_indexes:
            return full_name
        else:
            return full_name.split('_')[0]

    w_inches, h_inches = size[0] / dpi, size[1] / dpi
    fig = graphviz.Digraph(format=format,
                           #engine='neato',
                           engine='dot',
                           graph_attr={
                               'splines': 'true',  # tell Graphviz to eliminate overlaps/crossings where possible
                               'overlap': 'false',  # splines='true' gives you smooth curved edges instead of straight lines
                               'size': f'{w_inches},{h_inches}',
                               'ratio': 'fill',  # fill the whole image
                               'margin': '0',
                               'dpi': str(dpi),
                               #'K': '4.0',  # increase this to make the graph more spread out
                           })

    # make every edge use a 2‑point pen and 1.5× bigger arrowheads
    fig.attr('edge', penwidth='2', arrowsize='1.5')

    for u in H.nodes:
        color = rgb_to_hex(aa_color[get_residue_name(u)])
        fig.node(str(u),
                 label=get_node_name(u),
                 color=color,
                 style='filled',
                 fontsize='26')

    for u, v in H.edges:
        if (v, u) in H.edges:  # undirected edge labels a non-peptide bond
            if u < v: # to avoid duplicates
                fig.edge(str(u), str(v), color='red', dir='none', arrowhead='none')
        else:  # directed edge --> peptide bond
            fig.edge(str(u), str(v), color='blue')

    output_path = ensure_image_ext(output_path, format)
    output_path.parent.mkdir(exist_ok=True, parents=True)

    if format == 'svg':
        svg_bytes = fig.pipe(format='svg')
        svg_bytes = force_svg_pixel_size(svg_bytes, size[0], size[1], stretch=False)
        output_path.write_bytes(svg_bytes)
    else:  # PNG path stays identical to your original (center onto exact canvas)
        png_bytes = fig.pipe(format='png')
        img = Image.open(io.BytesIO(png_bytes))
        w_tgt, h_tgt = size
        canvas = Image.new('RGBA', (w_tgt, h_tgt), (255, 255, 255, 255))
        w_cur, h_cur = img.size
        x_off = max((w_tgt - w_cur) // 2, 0)
        y_off = max((h_tgt - h_cur) // 2, 0)
        canvas.paste(img, (x_off, y_off))
        canvas.save(output_path, format='PNG')


    json_output = fig.pipe(format='dot_json').decode('utf-8')
    with open( output_path.with_suffix('.json'), 'w') as f: 
       f.write(json_output)

    return fig

'''
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
'''
