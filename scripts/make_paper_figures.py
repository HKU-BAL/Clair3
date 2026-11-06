#!/usr/bin/env python3
"""Generate the two figures for the Bioinformatics submission on the Clair3 rewrite.

Figure 1  differential-validation ladder with the defect class detected at each layer
Figure 2  performance attribution for the HG002 Revio chr20 end-to-end benchmark

Every number in this script is taken from the project verification report
(docs/rust_optimization.md) so that the figures can be regenerated from the
archived artifact rather than edited by hand.

Usage:
    python3 scripts/make_paper_figures.py [--output-dir docs/figures]

Requires matplotlib. Output formats: SVG and PDF (vector, for submission) and
PNG at 600 dpi (colour) for inspection and presentation. Figures are rendered
at exactly 178 mm (double column) width; no bounding-box trimming is applied so
that the physical size matches the journal layout.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.patches as patches
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Journal geometry: 178 mm double column.
# ----------------------------------------------------------------------------
MM = 1 / 25.4
DOUBLE_COLUMN_MM = 178.0
DOUBLE_COLUMN_WIDTH = DOUBLE_COLUMN_MM * MM

# Okabe-Ito colour-blind-safe palette.
INK = "#1a1a1a"
MUTED = "#6b6b6b"
PANEL = "#f4f6f8"
CATEGORY_COLOURS = {
    "Argument model": "#0072B2",
    "Launcher model": "#56B4E9",
    "Value model": "#E69F00",
    "Data format": "#009E73",
    "ABI": "#D55E00",
    "Evaluation order": "#CC79A7",
    "Numerical": "#999999",
}

# ----------------------------------------------------------------------------
# Figure 1 data: the five validation layers and the defect classes each one
# first detected. Chip labels are wrapped over two lines so that they stay
# inside their boxes at 178 mm width.
# ----------------------------------------------------------------------------
LAYERS = [
    ("L1", "Component", "tensors, decoded records, selectors,\nLSTM layers, checkpoints"),
    ("L2", "Backend golden", "original-Python outputs regenerated\nper numerical backend"),
    ("L3", "Full workflow", "complete calling workflow\non real GIAB data"),
    ("L4", "Artifact", "training archives and checkpoints read\nby original Python / PyTorch"),
    ("L5", "Consumption", "artifacts used in downstream\ncalling and inference"),
]

DEFECTS = {
    "L1": [],
    "L2": [
        ("Numerical", "ONT LSTM golden mismatch\n(backend difference)"),
        ("Evaluation order", "Quality written into position\n(compiler-dependent order)"),
    ],
    "L3": [
        ("Argument model", "CreateTrainingTensorDirect\nalias absent"),
        ("Argument model", "call_fn / bin_fn share\ninternal ID 'output'"),
        ("Launcher model", "python clair3.py\nunsupported"),
        ("Value model", "--add_indel_length True\nrejected"),
        ("Data format", "gzip detected by\nextension only"),
        ("ABI", "x86_64 vs ARM64\nC++ string ABI"),
    ],
    "L4": [
        ("Data format", "Blosc filter unreadable\nby original h5py"),
    ],
    "L5": [],
}

# ----------------------------------------------------------------------------
# Figure 2 data: runtime attribution (docs/rust_optimization.md, Benchmark table).
# ----------------------------------------------------------------------------
RUNTIMES = [
    ("Original Python / CFFI implementation", 309.10, "#999999", False),
    ("First native Rust port", 301.89, "#56B4E9", False),
    ("Pristine dependencies, stock DP solver", 191.09, "#0072B2", False),
    ("Final configuration: stock DP solver,\n4 pipeline optimisations", 72.01, "#009E73", False),
    ("Excluded: custom DP solver", 42.62, "#D55E00", True),
]

# Line placed directly above the bar, stating the comparison it expresses.
DELTA_NOTES = {
    1: ("−2.3% versus the original (language change only)", "#0072B2"),
    3: ("−62.3% versus the pristine baseline (profiling and pipeline restructuring)", "#009E73"),
}


def _rounded_box(ax, x, y, width, height, face, edge, radius=0.012, linewidth=0.9, zorder=2):
    box = patches.FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=linewidth,
        zorder=zorder,
    )
    ax.add_patch(box)
    return box


def figure_one():
    """Differential-validation ladder and the defects each layer detects."""
    fig, ax = plt.subplots(figsize=(DOUBLE_COLUMN_WIDTH, 4.62))
    fig.subplots_adjust(left=0.010, right=0.990, top=0.990, bottom=0.010)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.0,
        0.990,
        "Figure 1. Differential-validation ladder for a scientific software rewrite",
        fontsize=8.4,
        fontweight="bold",
        color=INK,
        va="top",
    )
    ax.text(
        0.0,
        0.952,
        "Nine compatibility defect classes, seven first detected only at L3 or beyond; none at L1.",
        fontsize=7.0,
        color=MUTED,
        va="top",
    )
    ax.text(
        0.0,
        0.928,
        "Layers increase in scope from unit-level checks to scientific outputs.",
        fontsize=7.0,
        color=MUTED,
        va="top",
    )

    top = 0.868
    bottom = 0.142
    gaps = 0.018
    band = (top - bottom - gaps * (len(LAYERS) - 1)) / len(LAYERS)

    left_x, left_w = 0.005, 0.212
    chip_x = 0.268
    chip_w = 0.727

    for index, (code, name, description) in enumerate(LAYERS):
        band_top = top - index * (band + gaps)
        band_bottom = band_top - band

        _rounded_box(ax, left_x, band_bottom, left_w, band, PANEL, "#c9d2d9", radius=0.010)
        ax.text(
            left_x + 0.011,
            band_top - 0.030,
            code,
            fontsize=8.6,
            fontweight="bold",
            color="#0072B2",
            va="top",
        )
        ax.text(
            left_x + 0.045,
            band_top - 0.030,
            name,
            fontsize=7.5,
            fontweight="bold",
            color=INK,
            va="top",
        )
        ax.text(
            left_x + 0.011,
            band_top - 0.056,
            description,
            fontsize=6.0,
            color=MUTED,
            va="top",
            linespacing=1.25,
        )

        defects = DEFECTS[code]
        if not defects:
            ax.text(
                chip_x,
                band_top - band / 2,
                "no defect class first detected at this layer",
                fontsize=6.3,
                color=MUTED,
                style="italic",
                va="center",
            )
            continue

        columns = 3 if len(defects) > 4 else 2
        rows = (len(defects) + columns - 1) // columns
        chip_h = 0.040
        chip_gap_x = 0.008
        chip_gap_y = 0.006
        cell_w = (chip_w - chip_gap_x * (columns - 1)) / columns
        for position, (category, label) in enumerate(defects):
            row = position // columns
            column = position % columns
            chip_bottom = band_top - (row + 1) * chip_h - row * chip_gap_y
            colour = CATEGORY_COLOURS[category]
            cell_x = chip_x + column * (cell_w + chip_gap_x)
            _rounded_box(
                ax,
                cell_x,
                chip_bottom,
                cell_w,
                chip_h,
                "#ffffff",
                colour,
                radius=0.006,
                linewidth=0.8,
            )
            ax.add_patch(
                patches.Rectangle(
                    (cell_x + 0.004, chip_bottom + 0.005),
                    0.0035,
                    chip_h - 0.010,
                    facecolor=colour,
                    edgecolor="none",
                    zorder=3,
                )
            )
            ax.text(
                cell_x + 0.012,
                chip_bottom + chip_h / 2,
                label,
                fontsize=5.8,
                color=INK,
                va="center",
                ha="left",
                linespacing=1.25,
                zorder=4,
            )

    arrow_x = left_x + left_w + 0.011
    ax.add_patch(
        patches.FancyArrowPatch(
            (arrow_x, top - 0.015),
            (arrow_x, bottom + 0.015),
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=0.9,
            color="#9aa7b0",
            zorder=1,
        )
    )

    legend_y = 0.058
    column_width = 0.240
    for index, (category, colour) in enumerate(CATEGORY_COLOURS.items()):
        column = index % 4
        row = index // 4
        x = 0.005 + column * column_width
        y = legend_y - row * 0.036
        ax.add_patch(
            patches.Rectangle(
                (x, y - 0.013), 0.014, 0.014, facecolor=colour, edgecolor="none"
            )
        )
        ax.text(
            x + 0.020,
            y - 0.006,
            category,
            fontsize=6.4,
            color=INK,
            va="center",
            ha="left",
        )

    return fig


def figure_two():
    """Performance attribution for the end-to-end benchmark."""
    fig, ax = plt.subplots(figsize=(DOUBLE_COLUMN_WIDTH, 4.35))
    fig.subplots_adjust(left=0.012, right=0.988, top=0.855, bottom=0.190)
    ax.set_xlim(-12, 470)
    ax.set_ylim(-1.30, 5.55)
    ax.set_yticks([])

    fig.text(
        0.012,
        0.980,
        "Figure 2. Performance attribution for the HG002 Revio chr20 end-to-end benchmark (8 threads)",
        fontsize=8.4,
        fontweight="bold",
        color=INK,
        va="top",
        ha="left",
    )
    fig.text(
        0.012,
        0.932,
        "Lower is better. The language change alone contributed little; profiling and pipeline restructuring produced the gain.",
        fontsize=7.0,
        color=MUTED,
        va="top",
        ha="left",
    )
    fig.text(
        0.012,
        0.072,
        "One measurement campaign on the platform of Table 2; peak RSS 1.97 GB (original) → 3.93 GB (final).\n"
        "Emulated and CUDA timings are functional measurements, not native benchmarks.",
        fontsize=6.2,
        color=MUTED,
        va="top",
        ha="left",
        linespacing=1.35,
    )

    positions = [4.05, 2.90, 1.75, 0.60, -1.00]

    for index, ((label, value, colour, is_excluded), position) in enumerate(zip(RUNTIMES, positions)):
        bar = ax.barh(
            position,
            value,
            height=0.30,
            color=colour,
            edgecolor="#ffffff",
            linewidth=0.6,
            zorder=3,
        )
        if is_excluded:
            bar[0].set_hatch("///")
            bar[0].set_facecolor("#ffffff")
            bar[0].set_edgecolor(colour)
            bar[0].set_linewidth(1.1)
        ax.text(
            value + 6,
            position,
            f"{value:.2f} s",
            fontsize=7.2,
            va="center",
            ha="left",
            color=INK,
            fontweight="normal" if is_excluded else "bold",
        )
        label_lines = label.count("\n") + 1
        ax.text(
            0,
            position + 0.16,
            label,
            fontsize=6.9,
            va="bottom",
            ha="left",
            color=INK,
            style="italic" if is_excluded else "normal",
            linespacing=1.2,
        )
        if index in DELTA_NOTES:
            note, note_colour = DELTA_NOTES[index]
            ax.text(
                0,
                position + 0.30 + label_lines * 0.105 + 0.14,
                note,
                fontsize=6.0,
                va="bottom",
                ha="left",
                color=note_colour,
            )

    ax.set_xticks([0, 100, 200, 300])
    ax.set_xlabel("wall-clock time (s)", fontsize=7.4, color=INK)
    ax.tick_params(axis="x", labelsize=6.8, colors=MUTED, length=3)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#c9d2d9")
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#e6ebef", linewidth=0.7, zorder=0)

    ax.text(
        RUNTIMES[4][1] + 42,
        positions[4] + 0.02,
        "not shipped: changes upstream solver behaviour",
        fontsize=6.3,
        color="#D55E00",
        va="center",
        ha="left",
    )

    return fig


def save(fig, output_dir: Path, stem: str):
    output_dir.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams["svg.fonttype"] = "path"
    outputs = [
        (f"{stem}.svg", {}),
        (f"{stem}.pdf", {}),
        (f"{stem}.png", {"dpi": 600}),
        # Journal raster for line art is requested at 1,200 dpi; the 600 dpi
        # PNG above is kept for slides and quick inspection.
        (f"{stem}_1200dpi.png", {"dpi": 1200}),
    ]
    for name, kwargs in outputs:
        path = output_dir / name
        fig.savefig(path, facecolor="white", **kwargs)
        if path.suffix == ".svg":
            # Matplotlib emits trailing spaces in SVG path data; strip them so
            # that the committed figure files pass git's whitespace checks.
            lines = path.read_text().splitlines()
            path.write_text("\n".join(line.rstrip() for line in lines) + "\n")
        print(f"wrote {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "docs" / "figures",
    )
    arguments = parser.parse_args()

    save(figure_one(), arguments.output_dir, "figure1_validation_ladder")
    save(figure_two(), arguments.output_dir, "figure2_performance_attribution")


if __name__ == "__main__":
    main()
