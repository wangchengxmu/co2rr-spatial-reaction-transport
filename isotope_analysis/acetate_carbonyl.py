"""Sensitivity analysis for the acetate carbonyl isotope reporter.

The 144 fits vary window, coupling, line shape and baseline degree.  The
calculation tests the resolved two-bond acetate-methyl/carbonyl reporter only;
it does not add the older whole-molecule fits.

This is the historical, exploratory optimizer. Individual minima and linearized
uncertainties are platform-sensitive; neither a unique optimum nor a calibrated
analytical limit is established. Outputs are not source_fractions inputs.
See docs/ISOTOPE_ANALYSIS.md and docs/VALIDATION.md before interpretation.

Run ``python -m isotope_analysis.acetate_carbonyl --out PATH`` after installing
the pinned dependencies and supplying the portable ``12gas_13liquid`` FID.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import voigt_profile
from scipy.optimize import least_squares

from .nmr_processing import process_1h


SFO_HZ_PER_PPM = 500.17235215
WINDOWS = (0.018, 0.025, 0.035, 0.050)
COUPLINGS_HZ = (5.5, 5.914, 6.8)
LINE_SHAPES = ("gaussian", "lorentzian", "pseudo_voigt", "voigt")
BASELINE_DEGREES = (0, 1, 2)


def normalized_profile(x: np.ndarray, center: float, fwhm: float, eta: float, family: str) -> np.ndarray:
    sigma = fwhm / 2.354820045
    gamma = fwhm / 2.0
    gaussian = np.exp(-0.5 * ((x - center) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))
    lorentzian = gamma / (math.pi * ((x - center) ** 2 + gamma**2))
    if family == "gaussian":
        profile = gaussian
    elif family == "lorentzian":
        profile = lorentzian
    elif family == "voigt":
        profile = voigt_profile(x - center, sigma, gamma)
    else:
        profile = (1 - eta) * gaussian + eta * lorentzian
    return profile / max(np.trapezoid(profile, x), 1e-30)


def fit_profile(
    ppm: np.ndarray,
    signal: np.ndarray,
    center_guess: float,
    window_halfwidth: float,
    j2_hz: float,
    family: str,
    baseline_degree: int,
) -> dict:
    selected = (ppm >= center_guess - window_halfwidth) & (ppm <= center_guess + window_halfwidth)
    x = ppm[selected]
    y = signal[selected]
    xscaled = (x - x.mean()) / window_halfwidth
    d2 = j2_hz / (2 * SFO_HZ_PER_PPM)

    def predict(params: np.ndarray) -> np.ndarray:
        area, q, center, fwhm, eta = params[:5]
        base = params[5:]
        single = normalized_profile(x, center, fwhm, eta, family)
        doublet = 0.5 * (normalized_profile(x, center - d2, fwhm, eta, family) + normalized_profile(x, center + d2, fwhm, eta, family))
        baseline = sum(coefficient * xscaled**power for power, coefficient in enumerate(base))
        return area * ((1 - q) * single + q * doublet) + baseline

    area_guess = float(np.trapezoid(np.clip(y, 0, None), x))
    initial = np.r_[area_guess, 0.005, center_guess, 0.0030, 0.5, np.zeros(baseline_degree + 1)]
    lower = np.r_[0.0, 0.0, center_guess - 0.004, 0.0004, 0.0, np.repeat(-np.inf, baseline_degree + 1)]
    upper = np.r_[np.inf, 0.08, center_guess + 0.004, 0.0200, 1.0, np.repeat(np.inf, baseline_degree + 1)]
    result = least_squares(lambda params: predict(params) - y, initial, bounds=(lower, upper), max_nfev=6000)
    residual = y - predict(result.x)
    dof = max(len(x) - len(result.x), 1)
    sigma2 = float(np.sum(residual**2) / dof)
    covariance = sigma2 * np.linalg.pinv(result.jac.T @ result.jac)
    q_se = math.sqrt(max(float(covariance[1, 1]), 0.0))
    return {
        "window_halfwidth_ppm": window_halfwidth,
        "j2_hz": j2_hz,
        "line_shape": family,
        "baseline_degree": baseline_degree,
        "best_fraction_12CH3_13COO": float(result.x[1]),
        "formal_one_sided_95pct_upper": float(min(0.08, result.x[1] + 1.645 * q_se)),
        "q_standard_error_linearized": q_se,
        "rss": float(np.sum(residual**2)),
        "residual_mad_sigma": float(np.median(np.abs(residual - np.median(residual))) / 0.6745),
        "center_ppm": float(result.x[2]),
        "best_fwhm_ppm": float(result.x[3]),
        "eta": float(result.x[4]),
        "n_points": len(x),
        "n_parameters": len(result.x),
        "relative_residual": float(np.linalg.norm(residual) / max(np.linalg.norm(y), 1e-30)),
    }


def analyze(input_dir: Path, out: Path) -> dict:
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Use an empty output directory: {out}")
    out.mkdir(parents=True, exist_ok=True)
    ppm, signal, _ = process_1h("12gas_13liquid", input_dir, line_broadening_hz=0.25, size=65536, low_ppm=0.35)
    order = np.argsort(ppm)
    ppm, signal = ppm[order], signal[order]
    selected = (ppm >= 1.82) & (ppm <= 1.87)
    if not np.any(selected):
        raise ValueError("Acetate carbonyl search window is absent from the processed spectrum")
    center_guess = float(ppm[selected][np.argmax(signal[selected])])
    all_rows = [fit_profile(ppm, signal, center_guess, window, coupling, family, degree) for window in WINDOWS for coupling in COUPLINGS_HZ for family in LINE_SHAPES for degree in BASELINE_DEGREES]
    import csv

    fields = list(all_rows[0])
    with (out / "acetate_carbonyl_grid.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_rows)
    best_values = np.array([row["best_fraction_12CH3_13COO"] for row in all_rows])
    upper_values = np.array([row["formal_one_sided_95pct_upper"] for row in all_rows])
    residuals = np.array([row["relative_residual"] for row in all_rows])
    accepted = [row for row in all_rows if row["relative_residual"] <= np.quantile(residuals, 0.50)]
    accepted_best = np.array([row["best_fraction_12CH3_13COO"] for row in accepted])
    accepted_upper = np.array([row["formal_one_sided_95pct_upper"] for row in accepted])
    result = {
        "condition": "12gas_13liquid",
        "isotopomer": "12CH3-13COO",
        "spectral_signature": "central acetate methyl doublet from two-bond 2J(H,Ccarbonyl)",
        "j2_hz_measured_from_13gas_12liquid": 5.914,
        "n_nuisance_models": len(all_rows),
        "all_models": {"best_fraction_range": [float(best_values.min()), float(best_values.max())], "formal_95pct_upper_range": [float(upper_values.min()), float(upper_values.max())]},
        "lower_residual_half_of_models": {"n": len(accepted), "median_best_fraction": float(np.median(accepted_best)), "maximum_best_fraction": float(accepted_best.max()), "median_formal_95pct_upper": float(np.median(accepted_upper)), "maximum_formal_95pct_upper": float(accepted_upper.max())},
        "claim_boundary": "model-robustness envelope only; a calibrated 1-13C acetate spike is required for a defensible analytical detection limit",
    }
    (out / "acetate_carbonyl_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    axes[0].hist(100 * best_values, bins=np.arange(-0.025, 8.1, 0.25), color="#4c78a8", alpha=0.85)
    axes[0].set_xlabel("best-fit 12CH3-13COO (%)")
    axes[0].set_ylabel("nuisance models")
    axes[0].set_title("Line-shape sensitivity of the best fit")
    axes[1].hist(100 * upper_values, bins=np.arange(-0.025, 8.1, 0.25), color="#f58518", alpha=0.85)
    axes[1].set_xlabel("formal within-model 95% upper (%)")
    axes[1].set_ylabel("nuisance models")
    axes[1].set_title("Statistical limit varies with nuisance model")
    for axis in axes:
        axis.axvline(100 * np.median(accepted_best if axis is axes[0] else accepted_upper), color="black", ls="--", lw=1.2)
    for extension in ("png", "svg"):
        fig.savefig(out / f"acetate_carbonyl_sensitivity.{extension}", dpi=220)
    plt.close(fig)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent / "inputs" / "nmr")
    parser.add_argument("--out", type=Path, default=Path("runs/acetate_carbonyl"))
    args = parser.parse_args()
    print(json.dumps(analyze(args.input_dir, args.out), indent=2))


if __name__ == "__main__":
    main()
