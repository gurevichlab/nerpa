import io
from typing import (
    Dict,
    Tuple,
    List,
    NewType,
    Optional,
    Literal,
    NamedTuple,
)
from collections import defaultdict
from pathlib import Path

import graphviz
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from rdkit import Chem

from src.monomer_names_helper import MonomerNamesHelper
from src.build_output.chem_helper import MolRecord
from src.rban_parsing.rban_parser import (
    AtomId,
    MonomerIdx,
    Parsed_rBAN_Record
)
from src.rban_parsing.retrieve_nrp_variants import get_bond_direction
from src.generic.svg import (
    ensure_image_ext,
    force_svg_pixel_size,
    join_svgs_in_rectangle,
    svg_with_label,
)
import glasbey
from PIL import Image

import re

RGB = NewType('RGB', Tuple[float, float, float])
ColorHex = NewType('ColorHex', str)

def get_n_distinct_rgb_colors(N: int,
                              saturation: float = 0.3,
                              value = 0.8) -> List[RGB]:
    """
    Returns a list of N RGB triplets (floats in [0,1]) that are
    maximally distinct by using a glasbey palette.
    """
    
    colors = []
    bounds = hsv_to_glasbey_bounds(saturation, value) 
    color_hex = glasbey.create_palette(palette_size=N, lightness_bounds = bounds[0], chroma_bounds=bounds[1])

    for hex in color_hex:
        rgb = hex_to_rgb(ColorHex(hex.replace("#", "")))
        colors.append(rgb)
    return colors

def hsv_to_glasbey_bounds(saturation, value):
    """
    saturation, value in [0,1]
    """

    # perceived lightness
    L = 15 + value * 70

    # perceived chroma
    C = 10 + saturation * value * 80

    return (
        (max(0, L - 10), min(100, L + 10)),
        (max(0, C - 10), C + 10),
    )

def hex_to_rgb(hex: ColorHex) -> RGB:
    return RGB((int(hex[0:2], 16) / 255.0,
               int(hex[2:4], 16) / 255.0,
               int(hex[4:6], 16) / 255.0
    ))

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

def get_atom_id_to_display_name(
        record: Parsed_rBAN_Record,
        mon_idx: MonomerIdx
) -> AtomId:
    atoms_on_border: List[AtomId] = [
        mon_to_atom[mon_idx]
        for _mon_edge, mon_to_atom_lst in record.monomer_bonds.items()
        for mon_to_atom in mon_to_atom_lst
        if mon_idx in mon_to_atom
    ]

    # Choose any carbon atom which is not on the border
    repr_id = next(
        (atom_id
         for atom_id in record.monomers[mon_idx].atoms
         if record.atoms[atom_id].name == 'C' and atom_id not in atoms_on_border
         ),
        None
    )

    if repr_id is not None:
        return repr_id
    else:
        # Fallback: just choose any atom
        return record.monomers[mon_idx].atoms[0]

def get_node_colors(
        record: Parsed_rBAN_Record,
        monomer_names_helper: Optional[MonomerNamesHelper] = None
) -> Dict[MonomerIdx, RGB]:
    residues_with_colors = (
        monomer_names_helper.supported_residues
        if monomer_names_helper else
        STANDARD_AAS
    )
    aa_color = make_color_dict(labels=residues_with_colors)
    monomer_colors = dict()
    for mon_idx, mon_info in record.monomers.items():
        residue = (
            monomer_names_helper.parsed_name(
                mon_info.name,
                name_format='rBAN/Norine'
            ).residue
            if monomer_names_helper
            else mon_info.name.split('-')[-1].lower()
        )
        monomer_colors[mon_idx] = aa_color[residue]

    return monomer_colors


def get_node_labels(
        record: Parsed_rBAN_Record,
        with_rban_indexes: bool = True
) -> Dict[MonomerIdx, str]:
    labels = dict()
    for mon_idx, mon_info in record.monomers.items():
        labels[mon_idx] = (
            f'{mon_info.name}_{mon_idx}'
            if with_rban_indexes else
            mon_info.name
        )
    return labels


class MoleculeDrawingHelper:
    # required
    record: Parsed_rBAN_Record

    # computed
    mol: Chem.rdchem.Mol
    mon_colors: Dict[MonomerIdx, RGB]
    mon_labels: Dict[MonomerIdx, str]
    atom_labels: Dict[AtomId, str]
    atom_colors: Dict[AtomId, RGB]
    bonds_to_highlight: List[int]

    def __init__(
            self,
            record: Parsed_rBAN_Record,
            mon_colors: Optional[Dict[MonomerIdx, RGB]] = None,
            with_rban_indexes: bool = True,
            with_monomer_labels: bool = True,
            monomer_names_helper: Optional[MonomerNamesHelper] = None,
    ) -> None:
        self.record = record
        self.mol: Chem.rdchem.Mol = MolRecord.from_rban_record(record).mol
        self.mon_colors = (
            mon_colors
            if mon_colors is not None
            else get_node_colors(record, monomer_names_helper)
        )
        self.mon_labels: Dict[MonomerIdx, str] = get_node_labels(
            record,
            with_rban_indexes=with_rban_indexes,
        )
        self.atom_labels: Dict[AtomId, str] = {}
        self.atom_colors: Dict[AtomId, RGB] = {}
        self.bonds_to_highlight: List[int] = []

        for mon_idx, mon_info in record.monomers.items():
            if with_monomer_labels:
                repr_atom = get_atom_id_to_display_name(record, mon_idx)
                self.atom_labels[repr_atom] = self.mon_labels[mon_idx]

            for atom_id in mon_info.atoms:
                self.atom_colors[atom_id] = self.mon_colors[mon_idx]

            for atom1_id, atom2_id in record.atomic_bonds:
                if atom1_id not in mon_info.atoms or atom2_id not in mon_info.atoms:
                    continue

                bond = self.mol.GetBondBetweenAtoms(atom1_id, atom2_id)
                self.bonds_to_highlight.append(bond.GetIdx())

    def render(
            self,
            ext: Literal["svg", "png"] | None = "svg",
            output_path: Optional[Path] = None,
            size: Tuple[int, int] = (1000, 1000),
    ) -> str | bytes:
        if ext is None:
            if output_path is None:
                raise ValueError(
                    "Either ext or output_path must be specified."
                )
            ext = output_path.suffix[1:].lower()

        if ext not in ["svg", "png"]:
            raise ValueError(
                f"Unsupported format: {ext}. Use 'svg' or 'png'."
            )
        if (output_path is not None
            and output_path.suffix[1:].lower() != ext):
            raise ValueError(
                f"Output path extension {output_path.suffix} "
                f"does not match the specified format {ext}."
            )

        width, height = size

        if ext == "svg":
            drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        else:
            drawer = rdMolDraw2D.MolDraw2DCairo(width, height)

        options = drawer.drawOptions()
        for atom_id, label in self.atom_labels.items():
            options.atomLabels[atom_id] = label

        drawer.DrawMolecule(
            self.mol,
            highlightAtoms=list(self.atom_colors),
            highlightAtomColors=self.atom_colors,
            highlightBonds=self.bonds_to_highlight,
        )
        drawer.FinishDrawing()

        drawing: str | bytes = drawer.GetDrawingText()

        if output_path is not None:
            # output_path = ensure_image_ext(output_path, ext)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            if isinstance(drawing, str):
                output_path.write_text(drawing, encoding="utf-8")
            else:
                output_path.write_bytes(drawing)

        return drawing

    def get_drawing_data(self) -> dict:
        """Return drawing data with the following structure:

        {
            "a": [
                {
                    "i": str,    # atom index
                    "l": str,    # atom label
                    "x": float,
                    "y": float,
                    "z": float,
                },
                ...
            ],
            "b": [
                {
                    "i": int,    # bond index
                    "b": int,    # beginning atom index
                    "e": int,    # ending atom index
                    "o": float,  # bond order
                },
                ...
            ],
            "monomers": {
                str: list[int],  # monomer label -> atom indices
            },
            "highlightAtomColors": {
                int: tuple[float, float, float],  # RGB
            },
            "highlightBonds": list[int],
        }
        """
        drawing_data = {
            "a": [],
            "b": [],
            "monomers": {},
            "highlightAtomColors": self.atom_colors,
            "highlightBonds": self.bonds_to_highlight,
        }

        rdDepictor.Compute2DCoords(self.mol)
        conformer = self.mol.GetConformer()

        for atom in self.mol.GetAtoms():
            atom_id = atom.GetIdx()
            position = conformer.GetAtomPosition(atom_id)

            drawing_data["a"].append({
                "i": str(atom_id),
                "l": self.atom_labels.get(atom_id, atom.GetSymbol()),
                "x": position.x,
                "y": -position.y,
                "z": position.z,
            })

        for bond in self.mol.GetBonds():
            drawing_data["b"].append({
                "i": bond.GetIdx(),
                "b": bond.GetBeginAtomIdx(),
                "e": bond.GetEndAtomIdx(),
                "o": bond.GetBondTypeAsDouble(),
            })

        for mon_idx, mon_info in self.record.monomers.items():
            monomer_label = self.mon_labels[mon_idx]
            drawing_data["monomers"][monomer_label] = mon_info.atoms

        return drawing_data


class GraphDrawingHelper:
    record: Parsed_rBAN_Record
    mon_colors: Dict[MonomerIdx, RGB]
    mon_labels: Dict[MonomerIdx, str]

    def __init__(
            self,
            record: Parsed_rBAN_Record,
            mon_colors: Optional[Dict[MonomerIdx, RGB]] = None,
            with_rban_indexes: bool = True,
            monomer_names_helper: Optional[MonomerNamesHelper] = None,
    ) -> None:
        self.record = record
        self.mon_colors = (
            mon_colors
            if mon_colors is not None
            else get_node_colors(record, monomer_names_helper)
        )
        self.mon_labels = get_node_labels(
            record,
            with_rban_indexes=with_rban_indexes,
        )

    def render(
            self,
            size: Tuple[int, int] = (1000, 1000),
            dpi: int = 300,
    ) -> graphviz.Digraph:
        w_inches, h_inches = size[0] / dpi, size[1] / dpi
        fig = graphviz.Digraph(
            format='svg',
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
            }
        )

        # make every edge use a 2‑point pen and 1.5× bigger arrowheads
        fig.attr('edge', penwidth='2', arrowsize='1.5')

        for u in sorted(self.record.monomers.keys()):
            color = rgb_to_hex(self.mon_colors[u])
            fig.node(str(u),
                    label=self.mon_labels[u],
                    color=color,
                    style='filled',
                    fontsize='26')

        # print(f'Drawing {len(record.monomer_bonds)} monomer bonds...')
        # print(record.monomer_bonds)

        for u, v in sorted(self.record.monomer_bonds.keys()):
            amino_bond_dir = get_bond_direction((u, v), self.record)
            if amino_bond_dir is not None:
                if amino_bond_dir == (u, v):
                    fig.edge(str(u), str(v), color='blue', dir='forward', arrowhead='normal')
                else:
                    fig.edge(str(v), str(u), color='blue', dir='forward', arrowhead='normal')
            else:
                fig.edge(str(u), str(v), color='red', dir='none', arrowhead='none')

        return fig

    def save_fig(
            self,
            output_path: Path,
            size: Tuple[int, int] = (1000, 1000),
            dpi: int = 300,
    ) -> graphviz.Digraph:
        ext = output_path.suffix[1:].lower()
        output_path.parent.mkdir(exist_ok=True, parents=True)

        fig = self.render(size=size, dpi=dpi)

        match ext:
            case 'svg':
                svg_bytes = fig.pipe(format='svg')
                svg_bytes = force_svg_pixel_size(svg_bytes, size[0], size[1], stretch=False)
                output_path.write_bytes(svg_bytes)
            case 'png':
                png_bytes = fig.pipe(format='png')
                img = Image.open(io.BytesIO(png_bytes))
                w_tgt, h_tgt = size
                canvas = Image.new('RGBA', (w_tgt, h_tgt), (255, 255, 255, 255))
                w_cur, h_cur = img.size
                x_off = max((w_tgt - w_cur) // 2, 0)
                y_off = max((h_tgt - h_cur) // 2, 0)
                canvas.paste(img, (x_off, y_off))
                canvas.save(output_path, format='PNG')
            case 'json':
                output_path.write_text(fig.pipe(format='dot_json').decode('utf-8'),
                                       encoding='utf-8')
            case _:
                raise ValueError(f'Unsupported format: {ext}. Use "svg", "png", or "json".')

        return fig


class GraphDiffColors:
    SAME = RGB((0.8, 0.8, 0.8))  # light gray
    DELETED = RGB((1.0, 0.6, 0.6))  # light red
    INSERTED = RGB((0.6, 1.0, 0.6))  # light green
    SUBSTITUTED = RGB((1.0, 1.0, 0.6))  # light yellow

    original: Dict[MonomerIdx, RGB]
    modified: Dict[MonomerIdx, RGB]

def get_diff_colors(
        original: Parsed_rBAN_Record,
        modified: Parsed_rBAN_Record,
        old_to_new_map: List[Tuple[Optional[MonomerIdx], Optional[MonomerIdx]]]
) -> GraphDiffColors:
    # default: everything same unless proven otherwise
    original_colors: Dict[MonomerIdx, RGB] = {
        mon_idx: GraphDiffColors.SAME for mon_idx in original.monomers.keys()
    }
    modified_colors: Dict[MonomerIdx, RGB] = {
        mon_idx: GraphDiffColors.SAME for mon_idx in modified.monomers.keys()
    }

    for old_idx, new_idx in old_to_new_map:
        if old_idx is None and new_idx is None:
            continue

        if old_idx is None:
            # inserted in modified
            if new_idx in modified_colors:
                modified_colors[new_idx] = GraphDiffColors.INSERTED
            continue

        if new_idx is None:
            # deleted from original
            if old_idx in original_colors:
                original_colors[old_idx] = GraphDiffColors.DELETED
            continue

        # matched: substituted
        original_colors[old_idx] = GraphDiffColors.SUBSTITUTED
        modified_colors[new_idx] = GraphDiffColors.SUBSTITUTED

    out = GraphDiffColors()
    out.original = original_colors
    out.modified = modified_colors
    return out


class GraphDiffOut:
    original_diff_fig: graphviz.Digraph
    modified_diff_fig: graphviz.Digraph


def get_monomer_graph_diff(
        original: Parsed_rBAN_Record,
        modified: Parsed_rBAN_Record,
        old_to_new_map: List[Tuple[Optional[MonomerIdx], Optional[MonomerIdx]]],
        with_rban_indexes: bool = True,
        size: Tuple[int, int] = (1000, 1000),
        dpi: int = 300,
) -> GraphDiffOut:
    colors = get_diff_colors(original, modified, old_to_new_map)

    original_fig = GraphDrawingHelper(
        record=original,
        mon_colors=colors.original,
        with_rban_indexes=with_rban_indexes,
    ).render(
        size=size,
        dpi=dpi,
    )

    modified_fig = GraphDrawingHelper(
        record=modified,
        mon_colors=colors.modified,
        with_rban_indexes=with_rban_indexes,
    ).render(
        size=size,
        dpi=dpi,
    )

    out = GraphDiffOut()
    out.original_diff_fig = original_fig
    out.modified_diff_fig = modified_fig
    return out

class MoleculeDiffOut(NamedTuple):
    ext: Literal['svg', 'png']
    original_diff_data: str | bytes
    modified_diff_data: str | bytes


def get_molecule_diff(
        original: Parsed_rBAN_Record,
        modified: Parsed_rBAN_Record,
        old_to_new_map: List[Tuple[Optional[MonomerIdx], Optional[MonomerIdx]]],
        rban_indexes: bool = True,
        monomer_labels: bool = True,
        size: Tuple[int, int] = (1000, 1000),
        ext: Literal['svg', 'png'] = 'svg',
) -> MoleculeDiffOut:
    colors = get_diff_colors(original, modified, old_to_new_map)

    original_diff_data = MoleculeDrawingHelper(
        record=original,
        mon_colors=colors.original,
        with_rban_indexes=rban_indexes,
        with_monomer_labels=monomer_labels,
    ).render(ext=ext, size=size)

    modified_diff_data = MoleculeDrawingHelper(
        record=modified,
        mon_colors=colors.modified,
        with_rban_indexes=rban_indexes,
        with_monomer_labels=monomer_labels,
    ).render(ext=ext, size=size)

    return MoleculeDiffOut(
        ext=ext,
        original_diff_data=original_diff_data,
        modified_diff_data=modified_diff_data,
    )


def draw_molecule_diff(
        original: Parsed_rBAN_Record,
        modified: Parsed_rBAN_Record,
        original_score: float,
        modified_score: float,
        old_to_new_map: List[Tuple[Optional[MonomerIdx], Optional[MonomerIdx]]],
        output: Path,
        monomer_names_helper: Optional[MonomerNamesHelper] = None,
        rban_indexes: bool = True,
        monomer_labels: bool = True,
        size: Tuple[int, int] = (4000, 4000),
) -> None:
    diff_out = get_molecule_diff(
        original=original,
        modified=modified,
        old_to_new_map=old_to_new_map,
        rban_indexes=rban_indexes,
        monomer_labels=monomer_labels,
        size=(size[0]//2, size[1]//2),  # each molecule gets half of the total size
        ext='svg',
    )
    original_fig = MoleculeDrawingHelper(
        record=original,
        with_rban_indexes=rban_indexes,
        with_monomer_labels=monomer_labels,
        monomer_names_helper=monomer_names_helper,
    ).render(
        ext="svg",
        size=(size[0] // 2, size[1] // 2),
    )

    modified_fig = MoleculeDrawingHelper(
        record=modified,
        with_rban_indexes=rban_indexes,
        with_monomer_labels=monomer_labels,
        monomer_names_helper=monomer_names_helper,
    ).render(
        ext="svg",
        size=(size[0] // 2, size[1] // 2),
    )

    left_side_svg = join_svgs_in_rectangle([[original_fig], [diff_out.original_diff_data]])
    left_side_svg = svg_with_label(svg=left_side_svg,
                                   label=f"Original (score={original_score:.2f})",
                                   position="bottom")

    right_side_svg = join_svgs_in_rectangle([[modified_fig], [diff_out.modified_diff_data]])
    right_side_svg = svg_with_label(svg=right_side_svg,
                                    label=f"Modified (score={modified_score:.2f})",
                                    position="bottom")

    joined_svg = join_svgs_in_rectangle(
        [[left_side_svg, right_side_svg]],
    )
    output.write_bytes(joined_svg.encode('utf-8'))


def draw_monomer_graph_diff(
        original: Parsed_rBAN_Record,
        modified: Parsed_rBAN_Record,
        original_score: float,
        modified_score: float,
        old_to_new_map: List[Tuple[Optional[MonomerIdx], Optional[MonomerIdx]]],
        output: Path,
        monomer_names_helper: Optional[MonomerNamesHelper] = None,
        with_rban_indexes: bool = True,
        size: Tuple[int, int] = (4000, 4000),
        dpi: int = 300,
) -> None:
    cell_size = (size[0] // 2, size[1] // 2)

    diff_out = get_monomer_graph_diff(
        original=original,
        modified=modified,
        old_to_new_map=old_to_new_map,
        with_rban_indexes=with_rban_indexes,
        size=cell_size,
        dpi=dpi,
    )

    original_fig = GraphDrawingHelper(
        record=original,
        with_rban_indexes=with_rban_indexes,
    ).render(
        size=cell_size,
        dpi=dpi,
    )
    original_svg = force_svg_pixel_size(
        original_fig.pipe(format="svg"),
        cell_size[0],
        cell_size[1],
        stretch=False,
    ).decode("utf-8")

    modified_fig = GraphDrawingHelper(
        record=modified,
        with_rban_indexes=with_rban_indexes,
    ).render(
        size=cell_size,
        dpi=dpi,
    )
    modified_svg = force_svg_pixel_size(
        modified_fig.pipe(format='svg'),
        cell_size[0],
        cell_size[1],
        stretch=False
    ).decode('utf-8')

    original_diff_svg = force_svg_pixel_size(
        diff_out.original_diff_fig.pipe(format='svg'),
        cell_size[0],
        cell_size[1],
        stretch=False
    ).decode('utf-8')

    modified_diff_svg = force_svg_pixel_size(
        diff_out.modified_diff_fig.pipe(format='svg'),
        cell_size[0],
        cell_size[1],
        stretch=False
    ).decode('utf-8')

    left_side_svg = join_svgs_in_rectangle([[original_svg], [original_diff_svg]])
    left_side_svg = svg_with_label(svg=left_side_svg,
                                   label=f"Original (score={original_score:.2f})",
                                   position="bottom")

    right_side_svg = join_svgs_in_rectangle([[modified_svg], [modified_diff_svg]])
    right_side_svg = svg_with_label(svg=right_side_svg,
                                    label=f"Modified (score={modified_score:.2f})",
                                    position="bottom")

    joined_svg = join_svgs_in_rectangle(
        [[left_side_svg, right_side_svg]],
    )
    output.write_bytes(joined_svg.encode('utf-8'))
