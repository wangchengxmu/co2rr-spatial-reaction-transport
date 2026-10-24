"""Compare freshly generated analytical tables against frozen reference results."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent


def _json_value(value):
    """Convert pandas/numpy scalars to strict-JSON-friendly values."""
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_json_value(item) for item in value]
    missing = pd.isna(value) if value is not None else True
    if isinstance(missing, (bool, np.bool_)) and bool(missing):
        return None
    value = value.item() if isinstance(value, np.generic) else value
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _relative_difference(actual, reference):
    """Return a finite relative difference where possible.

    A zero reference is common for some nuisance terms.  Returning ``None``
    together with the explicit ``relative_difference_infinite`` flag keeps the
    diagnostic valid JSON instead of emitting non-standard ``Infinity``.
    """
    if pd.isna(actual) and pd.isna(reference):
        return 0.0, False
    if pd.isna(actual) or pd.isna(reference):
        return None, True
    absolute = abs(float(actual) - float(reference))
    if float(reference) == 0.0:
        return (0.0, False) if absolute == 0.0 else (None, True)
    return absolute / abs(float(reference)), False


def _difference(scope, key_rows, column, actual, reference, passed):
    absolute = None
    relative = None
    relative_infinite = False
    if isinstance(actual, (int, float, np.number)) and isinstance(reference, (int, float, np.number)):
        if not pd.isna(actual) and not pd.isna(reference):
            absolute = abs(float(actual) - float(reference))
            relative, relative_infinite = _relative_difference(actual, reference)
    return {
        "scope": scope,
        "key_rows": key_rows,
        "column": column,
        "actual": _json_value(actual),
        "reference": _json_value(reference),
        "absolute_difference": absolute,
        "relative_difference": relative,
        "relative_difference_infinite": relative_infinite,
        "pass": bool(passed),
    }


def _compare_table(actual: Path, reference: Path, keys: list[str], rtol=2e-5, atol=5e-9):
    """Compare a CSV and return its row count plus every failed field."""
    left = pd.read_csv(actual)
    right = pd.read_csv(reference)
    differences = []
    if list(left.columns) != list(right.columns):
        differences.append(_difference("csv", None, "columns", list(left.columns), list(right.columns), False))
        return len(left), differences
    if left.shape != right.shape:
        differences.append(_difference("csv", None, "shape", list(left.shape), list(right.shape), False))
    missing_keys = [key for key in keys if key not in left.columns]
    if missing_keys:
        differences.append(_difference("csv", None, "sort_keys", missing_keys, [], False))
        return len(left), differences
    left = left.sort_values(keys).reset_index(drop=True)
    right = right.sort_values(keys).reset_index(drop=True)
    row_count = min(len(left), len(right))
    for row_index in range(row_count):
        key_rows = {key: _json_value(right.iloc[row_index][key]) for key in keys}
        for column in right.columns:
            actual_value = left.iloc[row_index][column]
            reference_value = right.iloc[row_index][column]
            if pd.api.types.is_numeric_dtype(right[column]):
                passed = bool(np.isclose(actual_value, reference_value, rtol=rtol, atol=atol, equal_nan=True))
            else:
                actual_text = "" if pd.isna(actual_value) else str(actual_value)
                reference_text = "" if pd.isna(reference_value) else str(reference_value)
                passed = actual_text == reference_text
            if not passed:
                differences.append(_difference("row", key_rows, column, actual_value, reference_value, False))
    return len(left), differences


def compare(actual: Path, reference: Path, keys: list[str], rtol=2e-5, atol=5e-9):
    """Backward-compatible strict comparison helper.

    The command-line validator uses ``_compare_table`` so it can aggregate all
    mismatches before asserting.  Callers using this helper still receive the
    historical assertion on failure, but with the complete structured report.
    """
    count, differences = _compare_table(actual, reference, keys, rtol=rtol, atol=atol)
    if differences:
        raise AssertionError(json.dumps({"pass": False, "differences": differences}, indent=2))
    return count


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acetate",type=Path,required=True,help="Newly computed acetate output directory")
    args=p.parse_args()
    count, differences = _compare_table(args.acetate/"acetate_carbonyl_grid.csv",
        HERE/"reference_outputs/nmr/acetate_carbonyl_grid.csv",
        ["window_halfwidth_ppm","j2_hz","line_shape","baseline_degree"])
    actual=json.loads((args.acetate/"acetate_carbonyl_summary.json").read_text())
    reference=json.loads((HERE/"reference_outputs/nmr/acetate_carbonyl_summary.json").read_text())
    expected_models = reference["n_nuisance_models"]
    if actual.get("n_nuisance_models") != expected_models:
        differences.append(_difference("summary", {"summary_section": "top_level"},
            "n_nuisance_models", actual.get("n_nuisance_models"), expected_models, False))
    for key in reference["lower_residual_half_of_models"]:
        actual_value = actual.get("lower_residual_half_of_models", {}).get(key)
        reference_value = reference["lower_residual_half_of_models"][key]
        passed = actual_value is not None and bool(np.isclose(actual_value, reference_value,
            rtol=2e-5, atol=5e-9, equal_nan=True))
        if not passed:
            differences.append(_difference("summary", {"summary_section": "lower_residual_half_of_models"},
                key, actual_value, reference_value, False))
    report = {
        "acetate_models_compared": count,
        "regression": "passed" if not differences else "failed",
        "pass": not differences,
        "differences": differences,
        "meaning": "numerical reproduction, not an independently calibrated analytical limit",
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if differences:
        raise AssertionError("Acetate regression failed; complete diagnostics were printed above")


if __name__=="__main__":main()
