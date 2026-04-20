import io
from typing import Dict, Tuple, List, NewType, Optional, Literal
from collections import defaultdict
from pathlib import Path

import graphviz
from rdkit.Chem.Draw import rdMolDraw2D

from src.monomer_names_helper import MonomerNamesHelper
from src.build_output.chem_helper import MolRecord
from src.rban_parsing.rban_parser import (
    AtomId,
    MonomerIdx,
    Parsed_rBAN_Record
)
from src.generic.svg import (
    ensure_image_ext,
    force_svg_pixel_size,
)
import colorsys
from PIL import Image

import re

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

def make_color_dict(labels: List[str],
                    default_color: RGB = RGB((0.8, 0.8, 0.8)),) -> Dict[str, RGB]:
    colors = get_n_distinct_rgb_colors(len(labels))
    return defaultdict(lambda: default_color,
                       dict(zip(labels, colors)))


def get_atom_id_to_display_name(record: Parsed_rBAN_Record,
                                mon_idx: MonomerIdx) -> AtomId:
    atoms_on_border: List[AtomId] = [
        mon_to_atom[mon_idx]
        for _mon_edge, mon_to_atom_lst in record.monomer_bonds.items()
        for mon_to_atom in mon_to_atom_lst
        if mon_idx in mon_to_atom
    ]

    # Choose any carbon atom which is not on the border
    repr_id = next((atom_id
                    for atom_id in record.monomers[mon_idx].atoms
                    if record.atoms[atom_id].name == 'C' and atom_id not in atoms_on_border
                    ),
                   None)

    if repr_id is not None:
        return repr_id
    else:
        # Fallback: just choose any atom
        return record.monomers[mon_idx].atoms[0]

def get_node_colors(record: Parsed_rBAN_Record,
                    monomer_names_helper: Optional[MonomerNamesHelper] = None) -> Dict[MonomerIdx, RGB]:
    residues_with_colors = (
        monomer_names_helper.supported_residues
        if monomer_names_helper else
        STANDARD_AAS
    )
    aa_color = make_color_dict(labels=residues_with_colors)
    monomer_colors = dict()
    for mon_idx, mon_info in record.monomers.items():
        name = mon_info.name
        residue = (
           monomer_names_helper.parsed_name(name=name,
                                            name_format='rBAN/Norine').residue
            if monomer_names_helper else
            name.split('-')[-1].lower()  # e.g. 'ala' from 'NMe-Ala'
        )
 
        monomer_colors[mon_idx] = aa_color[residue]

    return monomer_colors

def get_node_labels(record: Parsed_rBAN_Record,
                    with_rban_indexes: bool = True) -> Dict[MonomerIdx, str]:
    labels = dict()
    for mon_idx, mon_info in record.monomers.items():
        labels[mon_idx] = (
            f'{mon_info.name}_{mon_idx}'
            if with_rban_indexes else
            mon_info.name
        )
    return labels


def draw_molecule(record: Parsed_rBAN_Record,
                  output_file: Path,
                  rban_indexes: bool = True,
                  monomer_labels: bool = True,
                  size: Tuple[int, int] = (1000, 1000),
                  monomer_names_helper: Optional[MonomerNamesHelper] = None) -> None:
    format = output_file.suffix[1:].lower()

    mon_colors = get_node_colors(record, monomer_names_helper)
    mon_labels = get_node_labels(record, with_rban_indexes=rban_indexes)
    mol, atom_id_to_index = MolRecord.from_rban_record(record)

    atom_colors = dict()
    bonds_to_highlight = list()
    atom_labels = dict()
    for mon_idx, mon_info in record.monomers.items():
        if monomer_labels:
            monomer_repr = get_atom_id_to_display_name(record, mon_idx)
            repr_atom_index = atom_id_to_index[monomer_repr]
            atom_labels[repr_atom_index] = mon_labels[mon_idx]

        for atom_id in mon_info.atoms:  # set atoms colors
            atom_index = atom_id_to_index[atom_id]
            atom_colors[atom_index] = mon_colors[mon_idx]

        for atom1_id, atom2_id in record.atomic_bonds:  # set bonds to highlight (color is deduced automatically)
            if atom1_id not in mon_info.atoms or atom2_id not in mon_info.atoms:
                continue  # only highlight bonds internal to the monomer

            index1, index2 = [atom_id_to_index[atom_id]
                              for atom_id in (atom1_id, atom2_id)]

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


def amino_bond_direction(bond: Tuple[MonomerIdx, MonomerIdx],
                         record: Parsed_rBAN_Record) -> Optional[Tuple[MonomerIdx, MonomerIdx]]:
    if bond not in record.monomer_bonds:
        raise ValueError(f'Bond {bond} not found in record')

    mon_to_atom_lst = record.monomer_bonds[bond]
    if len(mon_to_atom_lst) != 1:
        return None  # if multiple atomic bonds between the same monomers, we can't be sure which one is the "main" bond, so we won't classify it as amino

    atom_id1, atom_id2 = mon_to_atom_lst[0].values()
    if (atom_id1, atom_id2) in record.atomic_bonds:
        atomic_bond_info = record.atomic_bonds[(atom_id1, atom_id2)]
    elif (atom_id2, atom_id1) in record.atomic_bonds:
        atomic_bond_info = record.atomic_bonds[(atom_id2, atom_id1)]
    else:
        raise ValueError(f'No atomic bond found between atoms {atom_id1} and {atom_id2} for monomer bond {bond}')

    if atomic_bond_info.bond_type != 'AMINO':
        return None

    if record.atoms[atom_id1].name == 'C' and record.atoms[atom_id2].name == 'N':
        return (bond[0], bond[1])  # direction from C to N
    elif record.atoms[atom_id1].name == 'N' and record.atoms[atom_id2].name == 'C':
        return (bond[1], bond[0])  # direction from C to N
    else:
        raise ValueError(f'Unexpected atom names for amino bond between atoms {atom_id1} and {atom_id2}: {record.atoms[atom_id1].name}, {record.atoms[atom_id2].name}')
    

def draw_monomer_graph(record: Parsed_rBAN_Record,
                       output_path: Path,
                       with_rban_indexes: bool = True,
                       size: Tuple[int, int] = (1000, 1000),
                       dpi: int = 300,
                       monomer_names_helper: Optional[MonomerNamesHelper] = None) -> graphviz.Digraph:
    format = output_path.suffix[1:].lower()
    mon_labels = get_node_labels(record, with_rban_indexes=with_rban_indexes)
    mon_colors = get_node_colors(record, monomer_names_helper)

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

    for u in record.monomers.keys():
        color = rgb_to_hex(mon_colors[u])
        fig.node(str(u),
                 label=mon_labels[u],
                 color=color,
                 style='filled',
                 fontsize='26')

    # print(f'Drawing {len(record.monomer_bonds)} monomer bonds...')
    # print(record.monomer_bonds)

    for u, v in record.monomer_bonds.keys():
        amino_bond_dir = amino_bond_direction((u, v), record)
        if amino_bond_dir is not None:
            if amino_bond_dir == (u, v):
                fig.edge(str(u), str(v), color='blue', dir='forward', arrowhead='normal')
            else:
                fig.edge(str(v), str(u), color='blue', dir='forward', arrowhead='normal')
        else:
            fig.edge(str(u), str(v), color='red', dir='none', arrowhead='none')

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

    return fig
