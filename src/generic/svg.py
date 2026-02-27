from __future__ import annotations

import copy
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable


_SVG_NS = 'http://www.w3.org/2000/svg'
_XLINK_NS = 'http://www.w3.org/1999/xlink'

ET.register_namespace('', _SVG_NS)
ET.register_namespace('xlink', _XLINK_NS)


def _parse_length(value: str) -> float:
    """Parse an SVG length like '123', '123px', '123.4pt' -> 123.4 (units ignored)."""
    m = re.match(r'^\s*([+-]?\d+(?:\.\d+)?)', value)
    if not m:
        raise ValueError(f'Unsupported SVG length: {value!r}')
    return float(m.group(1))


def _svg_dims(root: ET.Element) -> tuple[float, float]:
    """Return (width, height) in user units, using width/height or falling back to viewBox."""
    w_attr = root.get('width')
    h_attr = root.get('height')
    if w_attr and h_attr:
        return (_parse_length(w_attr), _parse_length(h_attr))

    view_box = root.get('viewBox')
    if view_box:
        parts = view_box.replace(',', ' ').split()
        if len(parts) == 4:
            _, _, vw, vh = (float(p) for p in parts)
            return (vw, vh)

    raise ValueError('SVG is missing width/height and has no usable viewBox')


def join_svgs_side_by_side(svg_paths: list[Path], output_path: Path) -> None:
    """
    Join multiple SVGs into a single SVG laid out left-to-right.

    Assumptions:
      - All SVGs have the same height (within a tiny tolerance).
      - Each SVG has width/height, or at least a viewBox.
    """
    if not svg_paths:
        raise ValueError('No SVG paths provided')

    parsed: list[tuple[ET.Element, float, float]] = []
    for p in svg_paths:
        root = ET.parse(p).getroot()
        w, h = _svg_dims(root)
        parsed.append((root, w, h))

    base_h = parsed[0][2]
    if any(abs(h - base_h) > 1e-6 for _, _, h in parsed[1:]):
        raise ValueError('SVG heights do not match')

    total_w = sum(w for _, w, _ in parsed)

    out_root = ET.Element(f'{{{_SVG_NS}}}svg', attrib={
        'width': f'{total_w}',
        'height': f'{base_h}',
        'viewBox': f'0 0 {total_w} {base_h}',
    })

    x = 0.0
    for root, w, h in parsed:
        nested_attrib: dict[str, str] = {
            'x': f'{x}',
            'y': '0',
            'width': f'{w}',
            'height': f'{base_h}',
        }

        # Preserve viewBox + aspect ratio behavior if present.
        vb = root.get('viewBox')
        if vb:
            nested_attrib['viewBox'] = vb
        par = root.get('preserveAspectRatio')
        if par:
            nested_attrib['preserveAspectRatio'] = par

        nested = ET.SubElement(out_root, f'{{{_SVG_NS}}}svg', attrib=nested_attrib)

        # Copy over children (defs, g, paths, etc.) into the nested SVG.
        for child in list(root):
            nested.append(copy.deepcopy(child))

        x += w

    output_path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(out_root).write(output_path, encoding='utf-8', xml_declaration=True)
