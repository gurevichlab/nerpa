from __future__ import annotations

import copy
import io
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Literal, NamedTuple


_SVG_NS = 'http://www.w3.org/2000/svg'
_XLINK_NS = 'http://www.w3.org/1999/xlink'

ET.register_namespace('', _SVG_NS)
ET.register_namespace('xlink', _XLINK_NS)


class _ParsedSvg(NamedTuple):
    root: ET.Element
    w: float
    h: float
    view_box: str
    preserve_aspect_ratio: str | None


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

    raise ValueError("SVG is missing width/height and has no usable viewBox")


def _parse_svg_root(root: ET.Element) -> _ParsedSvg:
    w, h = _svg_dims(root)

    vb: str | None = root.get('viewBox')
    if vb is None:
        # Without a viewBox, setting width/height on a nested <svg> won't scale its contents.
        vb = f'0 0 {w} {h}'

    par: str | None = root.get('preserveAspectRatio')
    return _ParsedSvg(root=root, w=w, h=h, view_box=vb, preserve_aspect_ratio=par)


def _append_nested_svg(
    parent: ET.Element,
    *,
    parsed: _ParsedSvg,
    x: float,
    y: float,
    width: float,
    height: float,
) -> None:
    # Intent: embed an SVG under another SVG while scaling it using width/height + viewBox.
    nested_attrib: dict[str, str] = {
        'x': f'{x}',
        'y': f'{y}',
        'width': f'{width}',
        'height': f'{height}',
        'viewBox': parsed.view_box,
    }
    if parsed.preserve_aspect_ratio:
        nested_attrib['preserveAspectRatio'] = parsed.preserve_aspect_ratio

    nested = ET.SubElement(parent, f'{{{_SVG_NS}}}svg', attrib=nested_attrib)

    # Copy over children (defs, g, paths, etc.) into the nested SVG.
    for child in list(parsed.root):
        nested.append(copy.deepcopy(child))


def _svg_to_string(root: ET.Element) -> str:
    # ElementTree.tostring() doesn't support xml_declaration=True; write() does.
    buf = io.StringIO()
    ET.ElementTree(root).write(buf, encoding='unicode', xml_declaration=True)
    return buf.getvalue()


def join_svgs_in_rectangle(svgs: list[list[str]]) -> str:
    """
    Given a rectangle of SVG strings [row][col], produce one SVG:
      1) In each row, scale all SVGS to match the first SVG's height (row base height).
      2) Join each row horizontally.
      3) Scale rows so they all have the same total width (max row width).
      4) Join rows vertically.
    """
    if not svgs:
        raise ValueError('No rows provided')
    if any(not row for row in svgs):
        raise ValueError('All rows must be non-empty')

    # Build each row as its own SVG first.
    row_svgs: list[_ParsedSvg] = []
    row_widths: list[float] = []
    row_heights: list[float] = []

    for row in svgs:
        parsed_items: list[_ParsedSvg] = [_parse_svg_root(ET.fromstring(s)) for s in row]
        base_h: float = parsed_items[0].h

        # Intent: normalize all images in the row to the first image's height.
        scaled_ws: list[float] = [p.w * (base_h / p.h) for p in parsed_items]
        row_w: float = sum(scaled_ws)
        row_h: float = base_h

        row_root = ET.Element(f'{{{_SVG_NS}}}svg', attrib={
            'width': f'{row_w}',
            'height': f'{row_h}',
            'viewBox': f'0 0 {row_w} {row_h}',
        })

        # Intent: lay out the normalized row left-to-right.
        x: float = 0.0
        for p, scaled_w in zip(parsed_items, scaled_ws, strict=True):
            _append_nested_svg(
                row_root,
                parsed=p,
                x=x,
                y=0.0,
                width=scaled_w,
                height=base_h,
            )
            x += scaled_w

        row_svgs.append(_parse_svg_root(row_root))
        row_widths.append(row_w)
        row_heights.append(row_h)

    # Intent: make all rows the same width by scaling them up to the widest row.
    target_w: float = max(row_widths)

    scaled_row_heights: list[float] = []
    for row_w, row_h in zip(row_widths, row_heights, strict=True):
        scale: float = (target_w / row_w) if row_w else 1.0
        scaled_row_heights.append(row_h * scale)

    # Intent: vertically stack the scaled rows.
    total_h: float = sum(scaled_row_heights)

    out_root = ET.Element(f'{{{_SVG_NS}}}svg', attrib={
        'width': f'{target_w}',
        'height': f'{total_h}',
        'viewBox': f'0 0 {target_w} {total_h}',
    })

    y: float = 0.0
    for row_svg, scaled_h in zip(row_svgs, scaled_row_heights, strict=True):
        _append_nested_svg(
            out_root,
            parsed=row_svg,
            x=0.0,
            y=y,
            width=target_w,
            height=scaled_h,
        )
        y += scaled_h

    return _svg_to_string(out_root)


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
    if not svg_paths:
        raise ValueError('No SVG paths provided')

    parsed: list[_ParsedSvg] = [_parse_svg_root(ET.parse(p).getroot()) for p in svg_paths]

    base_h: float = parsed[0].h

    if (not force_same_heights) and any(abs(ps.h - base_h) > 1e-6 for ps in parsed[1:]):
        raise ValueError('SVG heights do not match')

    # Intent: optionally normalize heights to base_h so they can be tiled cleanly.
    scaled_ws: list[float] = [
        (ps.w * (base_h / ps.h)) if force_same_heights else ps.w
        for ps in parsed
    ]

    total_w: float = sum(scaled_ws) + padding * max(0, len(parsed) - 1)

    out_root = ET.Element(f'{{{_SVG_NS}}}svg', attrib={
        'width': f'{total_w}',
        'height': f'{base_h}',
        'viewBox': f'0 0 {total_w} {base_h}',
    })

    x: float = 0.0
    for i, (ps, scaled_w) in enumerate(zip(parsed, scaled_ws, strict=True)):
        _append_nested_svg(
            out_root,
            parsed=ps,
            x=x,
            y=0.0,
            width=scaled_w,
            height=base_h,
        )

        x += scaled_w
        if i != len(parsed) - 1:
            x += padding

    output_path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(out_root).write(output_path, encoding='utf-8', xml_declaration=True)


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
    parser.add_argument(
        '--force-same-heights',
        action='store_true',
        help="Scale all SVGs so their height matches the first SVG before joining.",
    )
    parser.add_argument(
        '--padding',
        type=float,
        default=0.0,
        help='Horizontal padding (gap) between SVGs, in user units. Default: 0.',
    )

    args = parser.parse_args(argv)

    join_svgs_side_by_side(
        args.svgs,
        args.output,
        force_same_heights=args.force_same_heights,
        padding=args.padding,
    )
    return 0


def force_svg_pixel_size(
    svg_bytes: bytes,
    width_px: int,
    height_px: int,
    stretch: bool = False,
) -> bytes:
    """
    Force SVG to have exact pixel width/height while preserving Graphviz output.
    This edits only the opening <svg ...> tag and avoids XML reserialization.
    """
    svg = svg_bytes.decode('utf-8', errors='strict')

    # Find opening <svg ...> tag
    m = re.search(r'<svg\b[^>]*>', svg, flags=re.IGNORECASE)
    if not m:
        return svg_bytes  # not an SVG? leave it alone

    tag = m.group(0)

    # Replace or insert width/height
    def upsert_attr(tag: str, name: str, value: str) -> str:
        if re.search(rf'\b{name}\s*=', tag):
            return re.sub(rf'\b{name}\s*=\s*"[^"]*"', f'{name}="{value}"', tag)
        # insert before the closing '>'
        return tag[:-1] + f' {name}="{value}">'

    tag2 = tag
    tag2 = upsert_attr(tag2, 'width', f'{width_px}px')
    tag2 = upsert_attr(tag2, 'height', f'{height_px}px')

    par = 'none' if stretch else 'xMidYMid meet'
    tag2 = upsert_attr(tag2, 'preserveAspectRatio', par)

    # Replace the tag in the document
    svg2 = svg[:m.start()] + tag2 + svg[m.end():]
    return svg2.encode('utf-8')


def ensure_image_ext(p: Path, fmt: Literal['png', 'svg']) -> Path:
    """Ensure path ends with .svg or .png without stripping mid-name parts like '.1'."""
    s = str(p)
    fmt2 = fmt.lower()
    if s.lower().endswith('.png') or s.lower().endswith('.svg'):
        s = s[:s.rfind('.')] + f'.{fmt2}'
    else:
        s = s + f'.{fmt2}'
    return Path(s)


if __name__ == '__main__':
    raise SystemExit(main())
