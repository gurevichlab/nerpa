from __future__ import annotations
from typing import Optional

from matplotlib.axes import Axes
from src.benchmarking.plots_data_helper import NerpaReport


def highlight_point(
        ax: Axes,
        x: float,
        y: float,
        label: Optional[str] = None,
        color: str = 'red'
) -> None:
    """
    Highlight a point (x, y) on the given Axes object with:
      - a scatter marker
      - dashed lines connecting to x and y axes
      - an optional label with arrow

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The Axes object to draw on.
    x : float
        X-coordinate of the point to highlight.
    y : float
        Y-coordinate of the point to highlight.
    label : Optional[str], default=None
        Text label to annotate the point.
    color : str, default='red'
        Color of the highlighted point and lines.

    Returns
    -------
    None
    """
    # Draw the point
    ax.scatter(x, y, color=color, s=50, zorder=5, marker='o', edgecolors='black', linewidths=1)

    # Draw dashed lines to axes (stopping at the point)
    ax.plot([x, x], [ax.get_ylim()[0], y], color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax.plot([ax.get_xlim()[0], x], [y, y], color='gray', linestyle='--', linewidth=1, alpha=0.7)

    # Add optional annotation
    if label:
        ax.annotate(
            label,
            xy=(x, y),
            xytext=(10, 10),
            textcoords='offset points',
            fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7, edgecolor='black'),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', lw=1)
        )

    # Refresh the plot
    ax.figure.canvas.draw_idle()
