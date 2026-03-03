from __future__ import annotations

from pathlib import Path
from typing import Any, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes


def my_new_figure(
    base_font_size: int = 18,
    inner_w_px: int = 360,
    inner_h_px: int = 300,
    outer_w_px: int = 450,
    outer_h_px: int = 370,
    margin_x_px: int = 15,
    margin_y_px: int = 10,
    dpi: int = 96,

) -> Tuple[Figure, Axes]:
    plt.rcParams.update({
        'font.size': base_font_size,  # base size
        'axes.titlesize': base_font_size + 2,
        'axes.labelsize': base_font_size - 2,
        'xtick.labelsize': base_font_size - 4,
        'ytick.labelsize': base_font_size - 4,
        'legend.fontsize': base_font_size - 4,
        'figure.titlesize': base_font_size,
    })

    # SVG is vector, but if you want "exact pixels" when viewed in a browser,
    # use dpi=96 so Matplotlib's pt units map to CSS px (1px = 1/96 inch).

    fig = plt.figure(figsize=(outer_w_px / dpi, outer_h_px / dpi), dpi=dpi)


    # Put the axes at the TR corner with a (magin_x, margin_y)
    # margin from the top and right edges of the figure -- so that the labels and title don't get cut off.
    # The axes will be inner_w_px wide and inner_h_px high.
    left_px = outer_w_px - inner_w_px - margin_x_px
    bottom_px = outer_h_px - inner_h_px - margin_y_px
    ax = fig.add_axes([
        left_px / outer_w_px,
        bottom_px / outer_h_px,
        inner_w_px / outer_w_px,
        inner_h_px / outer_h_px,
        ])

    return fig, ax
