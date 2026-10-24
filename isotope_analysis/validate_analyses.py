"""Compare freshly generated analytical tables against frozen reference results."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent


def compare(actual: Path,reference: Path,keys: list[str],rtol=2e-5,atol=5e-9):
    left=pd.read_csv(actual).sort_values(keys).reset_index(drop=True)
    right=pd.read_csv(reference).sort_values(keys).reset_index(drop=True)
    assert list(left.columns)==list(right.columns),(actual.name,"columns")
    assert left.shape==right.shape,(actual.name,"shape")
    for column in right:
        if pd.api.types.is_numeric_dtype(right[column]):
            assert np.allclose(left[column],right[column],rtol=rtol,atol=atol,equal_nan=True),(actual.name,column)
        else:
            assert left[column].fillna("").tolist()==right[column].fillna("").tolist(),(actual.name,column)
    return len(left)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acetate",type=Path,required=True,help="Newly computed acetate output directory")
    args=p.parse_args()
    count=compare(args.acetate/"acetate_carbonyl_grid.csv",
        HERE/"reference_outputs/nmr/acetate_carbonyl_grid.csv",
        ["window_halfwidth_ppm","j2_hz","line_shape","baseline_degree"])
    actual=json.loads((args.acetate/"acetate_carbonyl_summary.json").read_text())
    reference=json.loads((HERE/"reference_outputs/nmr/acetate_carbonyl_summary.json").read_text())
    assert actual["n_nuisance_models"]==reference["n_nuisance_models"]==144
    for key in reference["lower_residual_half_of_models"]:
        assert np.isclose(actual["lower_residual_half_of_models"][key],
                          reference["lower_residual_half_of_models"][key],rtol=2e-5,atol=5e-9),key
    print(json.dumps({"acetate_models_compared":count,"regression":"passed",
        "meaning":"numerical reproduction, not an independently calibrated analytical limit"},indent=2))


if __name__=="__main__":main()
