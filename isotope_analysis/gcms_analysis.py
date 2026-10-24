"""Estimate ethylene isotopologues using a declared fragmentation-response family.

Run: python -m isotope_analysis.gcms_analysis --out runs/gcms
Inputs are numerical mass spectra, not integrated chromatographic peak areas.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import json
from math import comb

import numpy as np
import pandas as pd
from scipy.optimize import nnls

HERE = Path(__file__).resolve().parent
NATURAL_13C = 0.0107


def natural_probabilities(n_carbons: int, p13: float = NATURAL_13C) -> list[float]:
    return [comb(n_carbons, k) * p13**k * (1-p13)**(n_carbons-k)
            for k in range(n_carbons+1)]


def ethylene_response_matrix(kernel: np.ndarray) -> np.ndarray:
    k26, k27, k28 = kernel
    return np.array([[k26, 0., 0.], [k27, k26, 0.], [k28, k27, k26],
                     [0., k28, k27], [0., 0., k28]])


def fit_kernel_from_natural(y: np.ndarray) -> np.ndarray:
    p = np.asarray(natural_probabilities(2))
    design = np.array([[p[0],0.,0.], [p[1],p[0],0.], [p[2],p[1],p[0]],
                       [0.,p[2],p[1]], [0.,0.,p[2]]])
    kernel, _ = nnls(design, y)
    return kernel/kernel[-1]


def fit_ethylene(y: np.ndarray, kernel: np.ndarray) -> dict:
    response = ethylene_response_matrix(kernel)
    amount, residual = nnls(response, y)
    abundance = amount/amount.sum()
    return {"p_12_12":float(abundance[0]), "p_12_13":float(abundance[1]),
            "p_13_13":float(abundance[2]),
            "f13_atom":float((abundance[1]+2*abundance[2])/2),
            "relative_residual":float(residual/np.linalg.norm(y)),
            "predicted":(response@amount).tolist()}


def ethylene_sensitivity(y: np.ndarray, kernels: list[np.ndarray]) -> tuple[dict, list[dict]]:
    """Original 81 by 81 response grid; accept residual <= minimum + 0.02."""
    r26 = [k[0] for k in kernels]
    r27 = [k[1] for k in kernels]
    candidates = []
    for k26 in np.linspace(max(.25,.9*min(r26)),1.1*max(r26),81):
        for k27 in np.linspace(max(.25,.9*min(r27)),1.1*max(r27),81):
            fit = fit_ethylene(y,np.array([k26,k27,1.]))
            candidates.append({**{k:float(fit[k]) for k in
                ["p_12_12","p_12_13","p_13_13","f13_atom","relative_residual"]},
                "kernel_26_over_28":float(k26),"kernel_27_over_28":float(k27)})
    best = min(row["relative_residual"] for row in candidates)
    accepted = [row for row in candidates if row["relative_residual"] <= best+.02]
    interval = {"best_relative_residual":best,"accepted_count":float(len(accepted))}
    for key in ["p_12_12","p_12_13","p_13_13","f13_atom","relative_residual"]:
        values = np.array([row[key] for row in accepted])
        for suffix,value in [("min",values.min()),("max",values.max()),("median",np.median(values))]:
            interval[f"{key}_{suffix}"] = float(value)
    return interval,accepted


def analyze(input_path: Path, out: Path) -> pd.DataFrame:
    data = pd.read_csv(input_path)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Use an empty output directory: {out}")
    out.mkdir(parents=True,exist_ok=True)
    spectra = {}
    for condition,group in data.groupby("condition"):
        group = group.sort_values("m_over_z")
        if group.m_over_z.tolist() != list(range(26,31)):
            raise ValueError(f"Expected masses 26 to 30 for {condition}")
        y = group.intensity.to_numpy(float)
        if not np.isfinite(y).all() or np.any(y<0):
            raise ValueError("Mass intensities must be finite and nonnegative")
        spectra[condition] = y
    kernels = [fit_kernel_from_natural(spectra[c]) for c in
               ["standard","12C_reaction_control"]]
    pd.DataFrame(kernels,columns=["k26_over_28","k27_over_28","k28_over_28"],
                 index=["standard","12C_reaction_control"]).to_csv(out/"reference_responses.csv",index_label="condition")
    envelopes=[];all_accepted=[]
    for condition in ["13gas_12liquid","12gas_13liquid"]:
        interval,accepted = ethylene_sensitivity(spectra[condition],kernels)
        rows = pd.DataFrame(accepted)
        p = rows.f13_atom
        rows["random_mixed"] = 2*p*(1-p)
        rows["mixed_over_random"] = rows.p_12_13/rows.random_mixed
        rows["mixed_depletion"] = rows.random_mixed-rows.p_12_13
        rows.insert(0,"condition",condition)
        assert np.allclose(rows.p_12_12+rows.p_12_13+rows.p_13_13,1,atol=1e-13)
        all_accepted.append(rows)
        envelopes.append({"condition":condition,**interval,
            **{f"{field}_{edge}":float(getattr(rows[field],edge)()) for field in
               ["mixed_over_random","mixed_depletion"] for edge in ["min","max"]}})
    table=pd.DataFrame(envelopes)
    table.to_csv(out/"ethylene_response_envelopes.csv",index=False)
    pd.concat(all_accepted,ignore_index=True).to_csv(out/"accepted_response_fits.csv",index=False)
    (out/"analysis_settings.json").write_text(json.dumps({
        "natural_13C":NATURAL_13C,"grid_points_per_response_ratio":81,
        "acceptance":"relative L2 residual <= best + 0.02",
        "uncertainty":"declared response-model range, not a confidence interval",
        "pair_order":"12C-13C includes both orders; random mixed fraction is 2p(1-p)",
        "scope":"ethylene only; no quantitative CO or methane isotope assignment"
    },indent=2)+"\n",encoding="utf-8")
    return table


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,default=HERE/"inputs/ethylene_mass_spectra.csv")
    parser.add_argument("--out",type=Path,default=Path("runs/gcms"))
    args=parser.parse_args()
    print(analyze(args.input,args.out).to_string(index=False))


if __name__ == "__main__":
    main()
