"""Reproduce the compact spatial-profile figures from saved screen results.

Usage from the repository root::

    python plot_results.py --results results --out runs/redrawn_profiles

Only the 30 production cases (10 diffusivity ratios × 3 CO2 transfer
thicknesses) are loaded.  Fine-grid and enlarged-domain verification cases are
identified from their ``config.json`` values and are intentionally excluded
from the comparison figures.  The plotting code uses the same colors, line
styles, axes and physical coordinate mapping as the archived compact plotting
style; it does not modify the saved CSV files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


RATIOS = [0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005]
COLORS = ["#A17B00", "#626262", "#8B4513", "#0072B2", "#009E73", "#C23878", "#D55E00", "#009BA8", "#151515", "#7841A8"]
STYLES = [":", "--", (0, (5, 1)), "-", "-.", (0, (3, 1)), ":", (0, (5, 1, 1, 1)), (0, (6, 2, 1, 2)), "--"]
FIELDS = [
    ("CO2_mM", r"CO$_2$"),
    ("OH_mM", r"OH$^-$"),
    ("HCO3_mM", r"HCO$_3^-$"),
    ("CO3_mM", r"CO$_3^{2-}$"),
    ("Cs_mM", r"Cs$^+$"),
    ("SO4_mM", r"SO$_4^{2-}$"),
    ("f13_CO2_13gas", r"$^{13}$CO$_2$ gas + $^{12}$CO$_2$ liquid"),
    ("f13_CO2_13liquid", r"$^{12}$CO$_2$ gas + $^{13}$CO$_2$ liquid"),
]
EXPECTED_COLUMNS = {field for field, _ in FIELDS} | {"z_nm", "liquid_y_um", "region"}
DELTAS = (10, 30, 50)


def _ratio_key(value: float) -> float:
    for ratio in RATIOS:
        if np.isclose(float(value), ratio, rtol=0, atol=1e-12):
            return ratio
    raise ValueError(f"Unexpected diffusivity ratio in production case: {value!r}")


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path.name}")
    return value


def _source_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_production_profiles(results_root: Path) -> tuple[dict[tuple[int, float], pd.DataFrame], dict[str, str]]:
    """Load exactly the production grid from config/result/profile files."""

    cases_root = results_root / "cases"
    if not cases_root.is_dir():
        raise FileNotFoundError(f"Expected {cases_root} containing case directories")
    data: dict[tuple[int, float], pd.DataFrame] = {}
    hashes: dict[str, str] = {}
    for case_dir in sorted(p for p in cases_root.iterdir() if p.is_dir()):
        config_path = case_dir / "config.json"
        if not config_path.is_file():
            continue
        config = _read_json(config_path)
        nc = int(config.get("n_catalyst", -1))
        ny = int(config.get("n_external", -1))
        nx = int(config.get("n_axial", -1))
        if (nc, ny, nx) != (48, 96, 24):
            continue
        delta = int(round(float(config.get("delta_CO2_um", -1))))
        if delta not in DELTAS:
            continue
        chi = _ratio_key(float(config.get("chi", float("nan"))))
        result_path = case_dir / "result.json"
        if not result_path.is_file() or _read_json(result_path).get("status") != "passed":
            raise ValueError(f"Production case is not marked passed: {case_dir.name}")
        profile_path = case_dir / "case_mean_profile.csv"
        if not profile_path.is_file():
            raise FileNotFoundError(f"Missing case_mean_profile.csv in {case_dir.name}")
        key = (delta, chi)
        if key in data:
            raise ValueError(f"Duplicate production configuration for d{delta}, chi={chi:g}")
        frame = pd.read_csv(profile_path)
        missing = sorted(EXPECTED_COLUMNS.difference(frame.columns))
        if missing:
            raise ValueError(f"{case_dir.name}/case_mean_profile.csv is missing columns: {', '.join(missing)}")
        numeric_columns = sorted(EXPECTED_COLUMNS.difference({"region"}))
        if frame.empty or not np.isfinite(frame[numeric_columns].to_numpy(dtype=float, na_value=np.nan)).all():
            raise ValueError(f"Non-finite or empty values in {case_dir.name}/case_mean_profile.csv")
        if set(frame["region"].unique()) != {"catalyst", "liquid"}:
            raise ValueError(f"Unexpected region labels in {case_dir.name}/case_mean_profile.csv")
        data[key] = frame
        hashes[profile_path.relative_to(results_root).as_posix()] = _source_hash(profile_path)
    expected = {(delta, chi) for delta in DELTAS for chi in RATIOS}
    missing = sorted(expected.difference(data), key=lambda item: (item[0], -item[1]))
    if missing:
        formatted = ", ".join(f"d{delta}/chi{chi:g}" for delta, chi in missing)
        raise ValueError(f"Expected all 30 production cases; missing: {formatted}")
    return data, hashes


def build(results_root: Path, out: Path) -> dict[str, Any]:
    """Write the four compact figures and a source-integrity QA record."""

    data, source_hashes = load_production_profiles(results_root)
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 11,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "svg.fonttype": "none",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
        }
    )
    qa: list[dict[str, Any]] = []
    curves: list[dict[str, Any]] = []
    edges: list[dict[str, Any]] = []

    def handles(ratios: list[float]) -> list[Line2D]:
        return [
            Line2D([0], [0], color=COLORS[RATIOS.index(ratio)], ls=STYLES[RATIOS.index(ratio)], lw=2, label=f"{ratio:g}")
            for ratio in ratios
        ]

    def export(fig: Any, stem: str) -> None:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        outside: list[str] = []
        for text in fig.findobj(matplotlib.text.Text):
            if text.get_visible() and text.get_text():
                box = text.get_window_extent(renderer)
                if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
                    outside.append(text.get_text())
        if outside:
            raise ValueError(f"Text outside canvas in {stem}: {outside}")
        fig.savefig(out / f"{stem}.png", dpi=180)
        fig.savefig(out / f"{stem}.svg")
        plt.close(fig)
        qa.append({"figure": stem, "outside_canvas_text": outside})

    def trace(ax: Any, delta: int, chi: float, field: str, cu_only: bool) -> None:
        frame = data[(delta, chi)]
        selected = frame[frame.region == "catalyst"] if cu_only else frame
        x = selected.z_nm.to_numpy() if cu_only else np.where(selected.region == "catalyst", selected.z_nm / 200, 1 + 1.6 * selected.liquid_y_um / 150)
        y = selected[field].to_numpy() * (100 if field.startswith("f13") else 1)
        index = RATIOS.index(chi)
        line, = ax.plot(x, y, color=COLORS[index], ls=STYLES[index], lw=1.75)
        if not np.array_equal(line.get_xdata(), x) or not np.array_equal(line.get_ydata(), y):
            raise AssertionError("Plotted data differ from the source CSV")
        curves.append({"delta_CO2_um": delta, "chi": chi, "field": field, "Cu_only": cu_only, "n_points": len(y), "exact_CSV_values": True})

    # Cu-only isotope comparison, with each transfer scenario on the same axes.
    fig, axes = plt.subplots(3, 2, figsize=(12, 11.5))
    fig.subplots_adjust(left=0.085, right=0.975, bottom=0.075, top=0.815, hspace=0.62, wspace=0.30)
    fig.suptitle("Isotope distribution within the 200 nm Cu layer", x=0.085, ha="left", y=0.987, fontsize=19)
    fig.text(0.085, 0.951, r"All calculated diffusivities  |  $^{12}$C fraction = 100% − $^{13}$C fraction", fontsize=12)
    fig.legend(handles=handles(RATIOS), loc="upper center", bbox_to_anchor=(0.54, 0.930), ncol=5, title=r"$D_{\mathrm{eff}}/D_{\mathrm{aq}}$ in Cu", fontsize=11, title_fontsize=12)
    for row, delta in enumerate(DELTAS):
        for col, (field, title) in enumerate(FIELDS[-2:]):
            ax = axes[row, col]
            for ratio in RATIOS:
                trace(ax, delta, ratio, field, True)
            ax.set_xlim(0, 200)
            ax.set_ylim(0, 100)
            ax.set_xticks([0, 50, 100, 150, 200])
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.set_xlabel("Cu depth from gas contact (nm)")
            ax.set_ylabel(r"$^{13}$C in CO$_2$ (%)")
            ax.set_title(title, loc="left", fontsize=12, weight="bold", pad=22)
            ax.text(0.98, 1.03, rf"CO$_2$ transfer length: {delta} µm", transform=ax.transAxes, ha="right", fontsize=10)
            position = ax.get_position()
            fig.text(position.x0 - 0.045, position.y1 + 0.03, chr(97 + 2 * row + col), weight="bold", fontsize=14)
    export(fig, "Cu_isotope_comparison")

    for delta in DELTAS:
        fig, axes = plt.subplots(4, 2, figsize=(12, 14))
        fig.subplots_adjust(left=0.085, right=0.975, bottom=0.06, top=0.83, hspace=0.80, wspace=0.30)
        fig.suptitle(rf"CO$_2$ equivalent transfer thickness: {delta} µm", x=0.085, ha="left", y=0.986, fontsize=19)
        fig.text(0.085, 0.951, r"300 mA cm$^{-2}$  |  200 nm Cu  |  CO$_2$-dependent CO$_2$RR current", fontsize=12)
        fig.legend(handles=handles(RATIOS), loc="upper center", bbox_to_anchor=(0.54, 0.934), ncol=5, title=r"$D_{\mathrm{eff}}/D_{\mathrm{aq}}$ in Cu", fontsize=11, title_fontsize=12)
        for index, (field, title) in enumerate(FIELDS):
            ax = axes.flat[index]
            isotope = field.startswith("f13")
            ax.axvspan(0, 1, color="#F3E9DC")
            ax.axvspan(1, 2.6, color="#EDF4F8")
            ax.axvline(1, color="#777", lw=0.8)
            for ratio in RATIOS:
                trace(ax, delta, ratio, field, False)
            maximum = max(float(data[(delta, ratio)][field].max()) for ratio in RATIOS)
            ax.set_xlim(0, 2.6)
            ax.set_xticks([0, 0.5, 1, 1.8, 2.6], ["0", "100", "200 | 0", "75", "150"])
            ax.set_ylim(0, 100 if isotope else maximum * 1.08)
            ax.set_ylabel(r"$^{13}$C in CO$_2$ (%)" if isotope else "Concentration (mM)")
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.set_title(title, loc="left", fontsize=12, weight="bold", pad=29)
            ax.text(0.19, 1.04, "Porous Cu", transform=ax.transAxes, ha="center", fontsize=10, color="#704D2C")
            ax.text(0.70, 1.04, "Flowing liquid", transform=ax.transAxes, ha="center", fontsize=10, color="#3B6178")
            ax.text(0.19, -0.23, "Cu depth (nm)", transform=ax.transAxes, ha="center", fontsize=10)
            ax.text(0.70, -0.23, "Distance from Cu (µm)", transform=ax.transAxes, ha="center", fontsize=10)
            box = ax.get_position()
            fig.text(box.x0 - 0.048, box.y1 + 0.029, chr(97 + index), fontsize=14, weight="bold")
        export(fig, f"spatial_profiles_d{delta}")
        for ratio in RATIOS:
            frame = data[(delta, ratio)]
            copper = frame[frame.region == "catalyst"]
            edges.append(
                {
                    "delta_CO2_um": delta,
                    "chi": ratio,
                    "first_Cu_cell_nm": float(copper.z_nm.iloc[0]),
                    "last_Cu_cell_nm": float(copper.z_nm.iloc[-1]),
                    "first_Cu_13C_percent_13gas": float(100 * copper.f13_CO2_13gas.iloc[0]),
                    "last_Cu_13C_percent_13gas": float(100 * copper.f13_CO2_13gas.iloc[-1]),
                }
            )

    after_hashes = {relative: _source_hash(results_root / relative) for relative in source_hashes}
    if after_hashes != source_hashes:
        raise RuntimeError("A source profile changed while plotting")
    qa_record = {
        "schema": "co2rr-spatial-plot-qa-v1",
        "production_case_count": len(data),
        "source_profiles_unchanged": True,
        "source_hashes_sha256": source_hashes,
        "figures": qa,
        "curve_checks": curves,
    }
    (out / "plot_qa.json").write_text(json.dumps(qa_record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    pd.DataFrame(edges).to_csv(out / "Cu_edge_values.csv", index=False)
    return qa_record


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot the 30 production spatial-model cases from saved CSV outputs.")
    parser.add_argument("--results", required=True, help="Results directory containing cases/*/config.json and CSV outputs")
    parser.add_argument("--out", required=True, help="New or empty figure output directory")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    results_root = Path(args.results).resolve()
    out = Path(args.out)
    if not out.is_absolute():
        out = (Path.cwd() / out).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"Output directory is not empty; choose a new directory: {out}")
    out.mkdir(parents=True, exist_ok=True)
    try:
        build(results_root, out)
    except Exception as exc:
        raise SystemExit(f"plot_results.py: {exc}") from exc
