"""Redraw manuscript Figure 4 directly from the packaged model-output CSV.

The numerical profiles are treated as immutable inputs. This script does not
rerun or refit the transport model. It recreates the current composite and the
seven separately editable panels in PNG/SVG (optional --pdf).
"""

from pathlib import Path
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from draw_model_schematic import draw_model_schematic


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "Figure4_source_data_all_sensitivities.csv"
OUT = ROOT / "redrawn"
OUT.mkdir(exist_ok=True)

MM = 1 / 25.4
SPLIT = 0.43
LIQUID_VIEW_UM = 150
PRIMARY_DELTA_UM = 10
RATIOS = [0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005]
COLORS = [
    "#A17B00",
    "#626262",
    "#8B4513",
    "#0072B2",
    "#009E73",
    "#C23878",
    "#D55E00",
    "#009BA8",
    "#151515",
    "#7841A8",
]
STYLES = [":", "--", (0, (5, 1)), "-", "-.", (0, (3, 1)), ":", (0, (5, 1, 1, 1)), (0, (6, 2, 1, 2)), "--"]

plt.rcParams.update(
    {
        "font.family": ["Arial", "DejaVu Sans"],
        "font.size": 9.6,
        "axes.labelsize": 9.6,
        "xtick.labelsize": 9.6,
        "ytick.labelsize": 9.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.bbox": None,
    }
)

TABLE = pd.read_csv(DATA_PATH)
DATA = {}
for ratio in RATIOS:
    subset = TABLE[
        np.isclose(TABLE["delta_CO2_um"], PRIMARY_DELTA_UM)
        & np.isclose(TABLE["chi"], ratio)
    ].copy()
    if len(subset) != 144:
        raise RuntimeError(f"Expected 144 rows for chi={ratio}; found {len(subset)}")
    DATA[ratio] = subset

CURVE_CHECKS = []


def save(fig, stem):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
            outside.append(text.get_text())
    if outside:
        raise RuntimeError(f"Text outside canvas for {stem}: {outside}")
    fig.savefig(OUT / f"{stem}.png", dpi=350, facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", facecolor="white")
    if "--pdf" in sys.argv:
        fig.savefig(OUT / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def legend_handles():
    return [
        Line2D([], [], color=color, ls=style, lw=1.6, label=f"{ratio:g}")
        for ratio, color, style in zip(RATIOS, COLORS, STYLES)
    ]


def profile(ax, column, title, letter, ylim=None):
    ax.axvspan(0, SPLIT, color="#F8EDE3", zorder=-5)
    ax.axvspan(SPLIT, 1, color="#EBF5F9", zorder=-5)
    ax.axvline(SPLIT, color="#7C919A", lw=0.8, ls="--", zorder=-3)
    is_isotope = column.startswith("f13")
    factor = 100 if is_isotope else 1

    for ratio, color, style in zip(RATIOS, COLORS, STYLES):
        source = DATA[ratio]
        x = np.where(
            source["region"].eq("catalyst"),
            SPLIT * source["z_nm"] / 200,
            SPLIT + (1 - SPLIT) * source["liquid_y_um"] / LIQUID_VIEW_UM,
        )
        y = source[column].to_numpy() * factor

        gas_boundary_added = column in {
            "CO2_mM",
            "f13_CO2_13gas",
            "f13_CO2_13liquid",
        }
        if gas_boundary_added:
            prescribed = {
                "CO2_mM": 33.0,
                "f13_CO2_13gas": 0.99,
                "f13_CO2_13liquid": 0.0107,
            }[column] * factor
            x = np.r_[0.0, x]
            y = np.r_[prescribed, y]

        line, = ax.plot(x, y, color=color, ls=style, lw=1.4)
        if not np.array_equal(line.get_xdata(), x) or not np.array_equal(line.get_ydata(), y):
            raise RuntimeError(f"Curve mismatch for {column}, chi={ratio}")
        CURVE_CHECKS.append(
            {
                "panel_field": column,
                "chi": ratio,
                "points": len(y),
                "gas_boundary_added": gas_boundary_added,
            }
        )

    maximum = max(DATA[ratio][column].max() * factor for ratio in RATIOS)
    ax.set_ylim(ylim if ylim is not None else (0, 100 if is_isotope else maximum * 1.08))
    ax.set_xlim(0, 1)
    ax.set_xticks(
        [0, SPLIT / 2, SPLIT, SPLIT + (1 - SPLIT) / 2, 1],
        ["0", "100", "200 | 0", "75", "150"],
    )
    ax.yaxis.set_major_locator(MaxNLocator(4))
    low, high = ax.get_ylim()
    ax.set_yticks([value for value in ax.get_yticks() if low <= value <= high])
    ax.tick_params(length=3, pad=3)
    ax.set_ylabel(r"$^{13}$C in CO$_2$ (%)" if is_isotope else "Concentration (mM)", labelpad=3)
    ax.text(SPLIT / 2, -0.21, "Cu depth (nm)", transform=ax.transAxes, ha="center", va="top", fontsize=9.6)
    ax.text((1 + SPLIT) / 2, -0.21, "Liquid distance (µm)", transform=ax.transAxes, ha="center", va="top", fontsize=9.6)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", pad=9)
    ax.text(-0.16, 1.105, letter, transform=ax.transAxes, fontweight="bold", fontsize=10)


PANELS = [
    ("CO2_mM", r"Molecular CO$_2$"),
    ("DIC_mM", "Total dissolved inorganic carbon"),
    ("OH_mM", r"Hydroxide, OH$^-$"),
    ("Cs_mM", r"Caesium, Cs$^+$"),
    ("f13_CO2_13gas", r"$^{13}$CO$_2$ gas + $^{12}$CO$_2$ liquid"),
    ("f13_CO2_13liquid", r"$^{12}$CO$_2$ gas + $^{13}$CO$_2$ liquid"),
]


def main():
    figure = plt.figure(figsize=(183 * MM, 238 * MM))
    schematic = figure.add_axes([0.025, 0.807, 0.95, 0.185])
    draw_model_schematic(schematic)
    for text in schematic.texts:
        text.set_fontsize(max(text.get_fontsize(), 8.7))
        if text.get_text() == "Flow, x":
            text.set_text("Effective flow, x")

    figure.legend(
        handles=legend_handles(),
        title=rf"$D_{{\mathrm{{eff}}}}/D_{{\mathrm{{aq}}}}$ in Cu  |  CO$_2$ transfer length: {PRIMARY_DELTA_UM} µm",
        loc="upper center",
        bbox_to_anchor=(0.53, 0.81),
        ncol=5,
        fontsize=9.6,
        title_fontsize=10,
        handlelength=2.7,
        columnspacing=1.7,
        labelspacing=0.4,
    )

    for index, (column, title) in enumerate(PANELS):
        axis = figure.add_axes(
            [0.11 if index % 2 == 0 else 0.615, 0.555 - (index // 2) * 0.24, 0.35, 0.163]
        )
        profile(axis, column, title, chr(98 + index))
    save(figure, "Figure4_redrawn")

    panel_a = plt.figure(figsize=(168 * MM, 43 * MM))
    schematic = panel_a.add_axes([0, 0, 1, 1])
    draw_model_schematic(schematic)
    for text in schematic.texts:
        text.set_fontsize(max(text.get_fontsize(), 8.7))
        if text.get_text() == "Flow, x":
            text.set_text("Effective flow, x")
    save(panel_a, "Figure4_panel_a_redrawn")

    for index, (column, title) in enumerate(PANELS):
        panel = plt.figure(figsize=(88 * MM, 61 * MM))
        profile(panel.add_axes([0.16, 0.27, 0.81, 0.54]), column, title, chr(98 + index))
        save(panel, f"Figure4_panel_{chr(98 + index)}_redrawn")

    qa = {
        "source_csv": DATA_PATH.name,
        "delta_CO2_um": PRIMARY_DELTA_UM,
        "ratios": RATIOS,
        "curve_checks": CURVE_CHECKS,
        "low_CO2_masked": False,
        "liquid_view_um": LIQUID_VIEW_UM,
    }
    (OUT / "redraw_QA.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    print(f"Redrew Figure 4 and seven panels from {DATA_PATH.name}")


if __name__ == "__main__":
    main()
