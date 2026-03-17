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


def join_svgs_side_by_side(
    svg_paths: list[Path],
    output_path: Path,
    force_same_heights: bool = False,
    padding: float = 0.0,
) -> None:
    """
    Join multiple SVGs into a single SVG laid out left-to-right.

    If force_same_heights=True, each SVG is scaled (by width/height) so its height
    matches the first SVG's height before joining. Widths are adjusted accordingly.

    padding is an optional horizontal gap (in SVG user units) inserted between figures.
    """
    from typing import NamedTuple

    if not svg_paths:
        raise ValueError('No SVG paths provided')

    class ParsedSvg(NamedTuple):
        root: ET.Element
        w: float          # original width
        h: float          # original height
        scaled_w: float   # width after optional scaling to base height

    parsed: list[ParsedSvg] = []
    for p in svg_paths:
        root = ET.parse(p).getroot()
        w, h = _svg_dims(root)
        parsed.append(ParsedSvg(root=root, w=w, h=h, scaled_w=w))

    base_h: float = parsed[0].h

    if force_same_heights:
        # Scale each SVG so its height becomes base_h; adjust width proportionally.
        parsed = [
            ParsedSvg(
                root=ps.root,
                w=ps.w,
                h=ps.h,
                scaled_w=ps.w * (base_h / ps.h),
            )
            for ps in parsed
        ]
    else:
        if any(abs(ps.h - base_h) > 1e-6 for ps in parsed[1:]):
            raise ValueError('SVG heights do not match')

    total_w: float = sum(ps.scaled_w for ps in parsed) + padding * max(0, len(parsed) - 1)

    out_root = ET.Element(f'{{{_SVG_NS}}}svg', attrib={
        'width': f'{total_w}',
        'height': f'{base_h}',
        'viewBox': f'0 0 {total_w} {base_h}',
    })

    x: float = 0.0
    for i, ps in enumerate(parsed):
        nested_attrib: dict[str, str] = {
            'x': f'{x}',
            'y': '0',
            'width': f'{ps.scaled_w}',
            'height': f'{base_h}',
        }

        # Ensure there's a viewBox so changing width/height actually scales the content.
        vb: str | None = ps.root.get('viewBox')
        if vb is None:
            vb = f'0 0 {ps.w} {ps.h}'
        nested_attrib['viewBox'] = vb

        # Preserve viewBox + aspect ratio behavior if present.
        par: str | None = ps.root.get('preserveAspectRatio')
        if par:
            nested_attrib['preserveAspectRatio'] = par

        nested = ET.SubElement(out_root, f'{{{_SVG_NS}}}svg', attrib=nested_attrib)

        # Copy over children (defs, g, paths, etc.) into the nested SVG.
        for child in list(ps.root):
            nested.append(copy.deepcopy(child))

        x += ps.scaled_w
        if i != len(parsed) - 1:
            x += padding

    output_path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(out_root).write(output_path, encoding='utf-8', xml_declaration=True)


# q: make this script runnable from the command line, accepting multiple SVG paths and an output path. Use argparse
def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description='Join multiple SVGs side-by-side into a single SVG.')
    parser.add_argument(
        'svgs',
        nargs='+',
        type=Path,
        help='Input SVG paths (one or more).',
    )
    parser.add_argument(
        '-o',
        '--output',
        required=True,
        type=Path,
        help='Output SVG path.',
    )

    # q: add an optional flag --force-same-heights that forces all SVGs to be scaled to the same height as the first one before joining
    parser.add_argument(
        '--force-same-heights',
        action='store_true',
        help='Scale all SVGs so their height matches the first SVG before joining.',
    )
    parser.add_argument(
        '--padding',
        type=float,
        default=0.0,
        help='Horizontal padding (gap) between SVGs, in user units. Default: 0.',
    )

    args = parser.parse_args(argv)

    join_svgs_side_by_side(args.svgs,
                           args.output,
                           force_same_heights=args.force_same_heights,
                           padding=args.padding)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
