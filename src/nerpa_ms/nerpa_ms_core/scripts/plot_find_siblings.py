# plot_clusters.py
#
# Usage:
#   python3 plot_clusters.py input.json output.svg
#
# input.json format:
# {
#   "title": "...",
#   "bucket_labels": ["0","1","2","3","4+"],
#   "bar_labels": ["total", "..."],
#   "heights": [
#       [10, 8, 7, ...],
#       [ 5, 4, 2, ...]
#   ]
# }

import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def read_payload(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def validate_payload(p: dict) -> None:
    bucket_labels = p["bucket_labels"]
    bar_labels = p["bar_labels"]
    heights = p["heights"]

    if len(bucket_labels) == 0:
        raise ValueError("bucket_labels is empty")
    if len(bar_labels) == 0:
        raise ValueError("bar_labels is empty")
    if len(heights) != len(bucket_labels):
        raise ValueError("heights must have same length as bucket_labels")

    for i, row in enumerate(heights):
        if len(row) != len(bar_labels):
            raise ValueError(
                f"heights[{i}] has {len(row)} values, expected {len(bar_labels)}"
            )


def plot_grouped_bars(payload: dict, out_svg: str) -> None:
    title = payload.get("title", "")
    bucket_labels = payload["bucket_labels"]
    bar_labels = payload["bar_labels"]
    heights = payload["heights"]

    num_buckets = len(bucket_labels)
    num_bars = len(bar_labels)

    # X positions for bucket centers.
    x = list(range(num_buckets))

    # Each bar series gets an offset around the bucket center.
    total_group_width = 0.8
    bar_width = total_group_width / num_bars

    fig_w = max(8.0, 1.3 * num_buckets)
    fig_h = 5.5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for bar_idx, bar_label in enumerate(bar_labels):
        offsets = [xi - total_group_width / 2 + (bar_idx + 0.5) * bar_width for xi in x]
        series = [heights[bucket_idx][bar_idx] for bucket_idx in range(num_buckets)]
        ax.bar(offsets, series, width=bar_width, label=bar_label)

    ax.set_title(title)
    ax.set_ylabel("Count")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xticks(x)
    ax.set_xticklabels(bucket_labels)
    ax.legend(loc="upper right", fontsize="small")
    ax.margins(x=0.02)

    fig.tight_layout()
    fig.savefig(out_svg, format="svg")
    plt.close(fig)


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: python3 plot_clusters.py input.json output.svg", file=sys.stderr)
        return 2

    input_json = sys.argv[1]
    out_svg = sys.argv[2]

    payload = read_payload(input_json)
    validate_payload(payload)
    plot_grouped_bars(payload, out_svg)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
