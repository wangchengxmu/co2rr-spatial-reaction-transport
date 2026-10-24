"""Fit isotope-sensitive proton-NMR reporter templates.

Run ``python -m isotope_analysis.nmr_templates --out PATH`` after placing the
portable FIDs and ``acquisition.json`` in ``isotope_analysis/inputs/nmr``.
The fit is the frozen template workflow used for formate, acetate, ethanol
methyl and propanol methyl reporters.  It writes CSV tables and PNG/SVG plots;
no files are created merely by importing this module.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
from scipy import sparse
from scipy.optimize import least_squares, lsq_linear
from scipy.sparse.linalg import spsolve

from .nmr_processing import asls, process_1h


CONDITIONS = ("12gas_12liquid", "13gas_12liquid", "12gas_13liquid")
REPORTER_CENTERS = {
    "formate": (8.3963, 0.025),
    "acetate": (1.8466, 0.022),
    "ethanol_methyl": (1.1127, 0.055),
    "propanol_methyl": (0.8213, 0.055),
}


def ascending(ppm: np.ndarray, signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(ppm)
    return ppm[order], signal[order]


def empirical_template(
    ppm: np.ndarray,
    signal: np.ndarray,
    center: float,
    half_width: float,
) -> tuple[np.ndarray, np.ndarray]:
    x, y = ascending(ppm, signal)
    selected = (x >= center - half_width) & (x <= center + half_width)
    rel = x[selected] - center
    profile = np.clip(y[selected], 0, None)
    edge = max(4, len(profile) // 15)
    profile -= np.linspace(np.median(profile[:edge]), np.median(profile[-edge:]), len(profile))
    profile = np.clip(profile, 0, None)
    area = np.trapezoid(profile, rel)
    return rel, profile / max(area, 1e-15)


def shifted_profile(x: np.ndarray, rel: np.ndarray, profile: np.ndarray, center: float) -> np.ndarray:
    return np.interp(x - center, rel, profile, left=0.0, right=0.0)


def lorentzian(x: np.ndarray, center: float, half_width: float) -> np.ndarray:
    y = half_width**2 / ((x - center) ** 2 + half_width**2)
    return y / np.trapezoid(y, x)


def area_lorentzian(x: np.ndarray, center: float, half_width: float) -> np.ndarray:
    profile = (half_width / np.pi) / ((x - center) ** 2 + half_width**2)
    return profile / np.trapezoid(profile, x)


def bounded_fit(
    y: np.ndarray,
    bases: list[np.ndarray],
    baseline_x: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    matrix = np.column_stack([*bases, np.ones_like(y), baseline_x - baseline_x.mean()])
    lower = np.array([0.0] * len(bases) + [-np.inf, -np.inf])
    upper = np.array([np.inf] * len(matrix.T))
    fit = lsq_linear(matrix, y, bounds=(lower, upper), lsmr_tol="auto", max_iter=500)
    predicted = matrix @ fit.x
    relative = np.linalg.norm(y - predicted) / max(np.linalg.norm(y), 1e-15)
    return fit.x, predicted, float(relative)


def fit_single_group(
    spectra: dict[str, tuple[np.ndarray, np.ndarray, dict]],
    templates: dict[str, tuple[np.ndarray, np.ndarray]],
    condition: str,
    group: str,
    region: tuple[float, float],
    center: float,
    delta: float,
    nuisances: list[float] | None = None,
) -> dict:
    ppm, signal, _ = spectra[condition]
    x, y = ascending(ppm, signal)
    selected = (x >= region[0]) & (x <= region[1])
    x = x[selected]
    y = y[selected]
    rel, profile = templates[group]
    unlabeled = shifted_profile(x, rel, profile, center)
    labeled = 0.5 * (
        shifted_profile(x, rel, profile, center - delta)
        + shifted_profile(x, rel, profile, center + delta)
    )
    nuisance_bases = [lorentzian(x, item, 0.0022) for item in nuisances or []]
    coefficients, predicted, relative = bounded_fit(y, [unlabeled, labeled, *nuisance_bases], x)
    total = coefficients[0] + coefficients[1]
    return {
        "condition": condition,
        "group": group,
        "unlabeled_area": float(coefficients[0]),
        "direct_13C_area": float(coefficients[1]),
        "direct_13C_fraction": float(coefficients[1] / total) if total > 0 else float("nan"),
        "relative_residual": relative,
        "x": x,
        "y": y,
        "predicted": predicted,
        "components": [unlabeled * coefficients[0], labeled * coefficients[1]],
        "nuisance_areas": [float(value) for value in coefficients[2 : 2 + len(nuisance_bases)]],
    }


def fit_alcohol_methyl(
    spectra: dict[str, tuple[np.ndarray, np.ndarray, dict]],
    templates: dict[str, tuple[np.ndarray, np.ndarray]],
    condition: str,
) -> dict:
    ppm, signal, _ = spectra[condition]
    x, y = ascending(ppm, signal)
    selected = (x >= 0.62) & (x <= 1.32)
    x = x[selected]
    y = y[selected]
    e_rel, e_profile = templates["ethanol_methyl"]
    p_rel, p_profile = templates["propanol_methyl"]
    e_center, p_center = 1.1115, 0.8205
    e_delta, p_delta = 0.1255, 0.1300
    e_u = shifted_profile(x, e_rel, e_profile, e_center)
    e_l = 0.5 * (shifted_profile(x, e_rel, e_profile, e_center - e_delta) + shifted_profile(x, e_rel, e_profile, e_center + e_delta))
    p_u = shifted_profile(x, p_rel, p_profile, p_center)
    p_l = 0.5 * (shifted_profile(x, p_rel, p_profile, p_center - p_delta) + shifted_profile(x, p_rel, p_profile, p_center + p_delta))
    nuisance_centers = [1.1753, 1.2623]
    nuisance = [lorentzian(x, center, 0.0022) for center in nuisance_centers]
    coefficients, predicted, relative = bounded_fit(y, [e_u, e_l, p_u, p_l, *nuisance], x)
    design = np.column_stack([e_u, e_l, p_u, p_l, *nuisance, np.ones_like(x), x - x.mean()])
    quiet = ((x >= 0.62) & (x <= 0.67)) | ((x >= 0.72) & (x <= 0.76)) | ((x >= 1.28) & (x <= 1.32))
    residual = y - predicted
    sigma = np.median(np.abs(residual[quiet] - np.median(residual[quiet]))) / 0.6745
    covariance = sigma**2 * np.linalg.pinv(design.T @ design)
    ethanol_upper_area = coefficients[1] + 3 * np.sqrt(max(covariance[1, 1], 0.0))
    propanol_upper_area = coefficients[3] + 3 * np.sqrt(max(covariance[3, 3], 0.0))
    e_total, p_total = coefficients[0] + coefficients[1], coefficients[2] + coefficients[3]
    return {
        "condition": condition,
        "group": "alcohol_methyl_joint",
        "ethanol_unlabeled_area": float(coefficients[0]),
        "ethanol_direct_13C_area": float(coefficients[1]),
        "ethanol_direct_13C_fraction": float(coefficients[1] / e_total) if e_total > 0 else float("nan"),
        "ethanol_direct_13C_fraction_noise_only_3sigma_upper": float(ethanol_upper_area / (coefficients[0] + ethanol_upper_area)),
        "propanol_unlabeled_area": float(coefficients[2]),
        "propanol_direct_13C_area": float(coefficients[3]),
        "propanol_direct_13C_fraction": float(coefficients[3] / p_total) if p_total > 0 else float("nan"),
        "propanol_direct_13C_fraction_noise_only_3sigma_upper": float(propanol_upper_area / (coefficients[2] + propanol_upper_area)),
        "nuisance_1p175_area": float(coefficients[4]),
        "nuisance_1p262_area": float(coefficients[5]),
        "relative_residual": relative,
        "x": x,
        "y": y,
        "predicted": predicted,
        "components": [e_u * coefficients[0], e_l * coefficients[1], p_u * coefficients[2], p_l * coefficients[3], nuisance[0] * coefficients[4], nuisance[1] * coefficients[5]],
    }


def fit_symmetric_satellites(
    spectra: dict[str, tuple[np.ndarray, np.ndarray, dict]],
    condition: str,
    group: str,
    region: tuple[float, float],
    center_guess: float,
    delta_guess: float,
    nuisance_centers: list[float] | None = None,
) -> dict:
    ppm, signal, _ = spectra[condition]
    x, y = ascending(ppm, signal)
    selected = (x >= region[0]) & (x <= region[1])
    x, y = x[selected], y[selected]
    nuisance_centers = nuisance_centers or []

    def model(parameters: np.ndarray) -> np.ndarray:
        area_u, area_l, center, delta, width_u, width_l, baseline, slope = parameters[:8]
        predicted = area_u * area_lorentzian(x, center, width_u)
        predicted += 0.5 * area_l * (area_lorentzian(x, center - delta, width_l) + area_lorentzian(x, center + delta, width_l))
        for amplitude, nuisance_center in zip(parameters[8:], nuisance_centers):
            predicted += amplitude * area_lorentzian(x, nuisance_center, 0.0022)
        return predicted + baseline + slope * (x - x.mean())

    positive_area = np.trapezoid(np.clip(y, 0, None), x)
    initial = np.array([0.9 * positive_area, 0.1 * positive_area, center_guess, delta_guess, 0.003, 0.002, 0.0, 0.0] + [0.005 * positive_area] * len(nuisance_centers))
    lower = np.array([0.0, 0.0, center_guess - 0.008, delta_guess - 0.012, 0.00025, 0.00025, -np.inf, -np.inf] + [0.0] * len(nuisance_centers))
    upper = np.array([np.inf, np.inf, center_guess + 0.008, delta_guess + 0.012, 0.006, 0.004, np.inf, np.inf] + [np.inf] * len(nuisance_centers))
    result = least_squares(lambda values: model(values) - y, initial, bounds=(lower, upper), max_nfev=5000)
    area_u, area_l = result.x[:2]
    return {
        "condition": condition,
        "group": group + "_flexible_symmetric",
        "unlabeled_area": float(area_u),
        "direct_13C_area": float(area_l),
        "direct_13C_fraction": float(area_l / (area_u + area_l)),
        "center_ppm": float(result.x[2]),
        "satellite_half_separation_ppm": float(result.x[3]),
        "unlabeled_half_width_ppm": float(result.x[4]),
        "satellite_half_width_ppm": float(result.x[5]),
        "relative_residual": float(np.linalg.norm(result.fun) / max(np.linalg.norm(y), 1e-15)),
        "x": x,
        "y": y,
        "predicted": model(result.x),
        "components": [],
    }


def _summary_rows(fits: list[dict]) -> list[dict]:
    rows = []
    for fit in fits:
        row = {key: value for key, value in fit.items() if key not in {"x", "y", "predicted", "components", "nuisance_areas"}}
        if "nuisance_areas" in fit:
            row["nuisance_areas"] = json.dumps(fit["nuisance_areas"])
        rows.append(row)
    return rows


def _write_fit_curves(fits: list[dict], path: Path) -> None:
    rows: list[dict[str, float | str]] = []
    for fit in fits:
        components = fit.get("components", [])
        for i, (ppm, observed, predicted) in enumerate(zip(fit["x"], fit["y"], fit["predicted"])):
            row: dict[str, float | str] = {"condition": fit["condition"], "group": fit["group"], "ppm": float(ppm), "observed": float(observed), "predicted": float(predicted)}
            for j, component in enumerate(components):
                row[f"component_{j + 1}"] = float(component[i])
            rows.append(row)
    columns = sorted({key for row in rows for key in row}, key=lambda item: (item not in {"condition", "group", "ppm", "observed", "predicted"}, item))
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _plot_fits(fits: list[dict], out: Path) -> None:
    fig, axes = plt.subplots(4, 4, figsize=(18, 14), constrained_layout=True)
    for ax, fit in zip(axes.flat, fits):
        ax.plot(fit["x"], fit["y"], color="black", lw=0.8, label="processed raw spectrum")
        ax.plot(fit["x"], fit["predicted"], color="#e45756", lw=1.1, label="joint fit")
        for component in fit.get("components", []):
            ax.plot(fit["x"], component, lw=0.65, alpha=0.75)
        ax.set_title(f'{fit["condition"]}: {fit["group"]}')
        ax.invert_xaxis()
        ax.set_yticks([])
    axes[0, 0].legend(frameon=False, fontsize=8)
    for ax in axes.flat[len(fits) :]:
        ax.axis("off")
    for extension in ("png", "svg"):
        fig.savefig(out / f"joint_template_fits.{extension}", dpi=220)
    plt.close(fig)


def analyze(input_dir: Path, out: Path) -> list[dict]:
    """Run all frozen reporter fits and return summary rows."""
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Use an empty output directory: {out}")
    out.mkdir(parents=True, exist_ok=True)
    spectra = {condition: process_1h(condition, input_dir) for condition in CONDITIONS}
    for condition, (ppm, signal, _) in spectra.items():
        with (out / f"{condition}_processed_spectrum.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["ppm", "intensity"])
            writer.writerows(zip(ppm, signal))
    control_ppm, control_signal, _ = spectra["12gas_12liquid"]
    templates = {name: empirical_template(control_ppm, control_signal, center, half_width) for name, (center, half_width) in REPORTER_CENTERS.items()}
    fits: list[dict] = []
    for condition in CONDITIONS:
        fits.append(fit_single_group(spectra, templates, condition, "formate", (8.05, 8.75), 8.3935, 0.1955))
        fits.append(fit_single_group(spectra, templates, condition, "acetate", (1.65, 2.04), 1.8445, 0.1280))
        fits.append(fit_alcohol_methyl(spectra, templates, condition))
    for condition in ("12gas_12liquid", "12gas_13liquid"):
        fits.append(fit_symmetric_satellites(spectra, condition, "formate", (8.08, 8.70), 8.394, 0.1955, nuisance_centers=[8.3523]))
        fits.append(fit_symmetric_satellites(spectra, condition, "acetate", (1.66, 2.03), 1.845, 0.128))
    rows = _summary_rows(fits)
    fields = sorted({key for row in rows for key in row})
    with (out / "template_fit_summary.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    _write_fit_curves(fits, out / "template_fit_curves.csv")
    settings = {
        "conditions": CONDITIONS, "line_broadening_hz": 0.3,
        "zero_fill_size": 32768, "ppm_window": [0.55, 9.1],
        "baseline": {"method": "ASLS", "lambda": 2e7, "p": 0.001, "iterations": 12},
        "curve_components": {
            "formate_and_acetate": {"component_1": "ordinary attached-12C pattern", "component_2": "both attached-13C branches"},
            "alcohol_methyl_joint": {"component_1": "ethanol methyl attached-12C", "component_2": "ethanol methyl attached-13C",
                "component_3": "n-propanol methyl attached-12C", "component_4": "n-propanol methyl attached-13C",
                "component_5": "nuisance at 1.1753 ppm", "component_6": "nuisance at 1.2623 ppm"},
            "flexible_symmetric": "total fitted curve supplied; separate components are not exported"
        },
        "curated_intervals": "Separate input table; not every endpoint is computed by this preferred operator"
    }
    (out / "analysis_settings.json").write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")
    _plot_fits(fits, out)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent / "inputs" / "nmr")
    parser.add_argument("--out", type=Path, default=Path("runs/nmr_templates"))
    args = parser.parse_args()
    rows = analyze(args.input_dir, args.out)
    print(json.dumps(rows, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
