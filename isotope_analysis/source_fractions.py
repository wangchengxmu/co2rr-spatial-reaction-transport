"""Convert reporter-carbon isotope ranges to gas-/liquid-carbon source ranges.

Run: python -m isotope_analysis.source_fractions --out runs/source_fractions
NMR intervals are curated analytical inputs, not raw data or statistical CIs.
Ethylene intervals are taken from the saved or newly recomputed response family.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent


def gas_fraction(f13,feed: str,natural: float=.0107,enriched: float=.99):
    """Return the *unclipped* two-endmember source fraction."""
    if not 0 <= natural < enriched <= 1:
        raise ValueError("Require 0 <= natural < enriched <= 1")
    f13=np.asarray(f13,dtype=float)
    if not np.isfinite(f13).all() or np.any((f13<0)|(f13>1)):
        raise ValueError("Isotope fractions must lie between zero and one")
    if feed=="13gas_12liquid":return (f13-natural)/(enriched-natural)
    if feed=="12gas_13liquid":return (enriched-f13)/(enriched-natural)
    raise ValueError(f"Unknown isotope-feed condition: {feed}")


def analyze(nmr: Path,gcms: Path,out: Path,natural: float=.0107,enriched: float=.99):
    observations=pd.read_csv(nmr)
    eth=pd.read_csv(gcms)
    rows=observations.to_dict("records")
    for _,row in eth.iterrows():
        rows.append(dict(product="ethylene",reporter="mean of both carbons",feed=row.condition,
            f13_low=row.f13_atom_min,f13_high=row.f13_atom_max,
            evidence="mass-fragment response family",interval_origin="recomputed response envelope"))
    result=[]
    for row in rows:
        if row["f13_low"]>row["f13_high"]:raise ValueError("Reversed isotope interval")
        raw=np.sort(gas_fraction([row["f13_low"],row["f13_high"]],row["feed"],natural,enriched))
        clipped=np.clip(raw,0.,1.)
        result.append({**row,"natural_13C":natural,"enriched_13C":enriched,
            "gas_source_low_raw":raw[0],"gas_source_high_raw":raw[1],
            "gas_source_low":clipped[0],"gas_source_high":clipped[1],
            "liquid_source_low":1-clipped[1],"liquid_source_high":1-clipped[0],
            "physical_boundary_clipped":bool(np.any(raw!=clipped))})
    per_feed=pd.DataFrame(result)
    combined=[]
    for (product,reporter),group in per_feed.groupby(["product","reporter"],sort=False):
        assert len(group)==2
        combined.append(dict(product=product,reporter=reporter,
            gas_source_low=group.gas_source_low.min(),gas_source_high=group.gas_source_high.max(),
            liquid_source_low=group.liquid_source_low.min(),liquid_source_high=group.liquid_source_high.max(),
            feed_midpoint_difference=float(np.ptp((group.gas_source_low+group.gas_source_high)/2)),
            interpretation="span of the two independent feed estimates, not a pooled statistical interval"))
    if out.exists() and any(out.iterdir()):raise FileExistsError(f"Use an empty output directory: {out}")
    out.mkdir(parents=True,exist_ok=True)
    per_feed.to_csv(out/"source_fractions_by_feed.csv",index=False)
    pd.DataFrame(combined).to_csv(out/"product_source_summary.csv",index=False)
    (out/"source_assumptions.json").write_text(json.dumps({"natural_13C":natural,"enriched_13C":enriched,
        "equation":"F13 = g F13_gas + (1-g) F13_liquid",
        "local_endmembers":"feed enrichments used as assumptions; local reacting-reservoir enrichments not measured",
        "out_of_bounds":"raw estimates retained alongside physical [0,1] clipping",
        "NMR_scope":"numeric percentages refer to named reporter carbons; other resolved carbons have an independent assignment ledger",
        "not_a_depth":"source fraction is not a unique reaction depth or a mechanistic assignment"
    },indent=2)+"\n",encoding="utf-8")
    return per_feed,pd.DataFrame(combined)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nmr",type=Path,default=HERE/"inputs/nmr_reporter_intervals.csv")
    parser.add_argument("--gcms",type=Path,default=HERE/"reference_outputs/gcms/ethylene_response_envelopes.csv")
    parser.add_argument("--natural",type=float,default=.0107)
    parser.add_argument("--enriched",type=float,default=.99)
    parser.add_argument("--out",type=Path,default=Path("runs/source_fractions"))
    args=parser.parse_args()
    _,summary=analyze(args.nmr,args.gcms,args.out,args.natural,args.enriched)
    print(summary.to_string(index=False))


if __name__=="__main__":main()
