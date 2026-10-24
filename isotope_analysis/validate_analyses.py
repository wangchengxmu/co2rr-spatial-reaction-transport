"""Compare freshly generated analytical tables against frozen reference results."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent

GRID_COLUMNS = [
    "window_halfwidth_ppm", "j2_hz", "line_shape", "baseline_degree",
    "best_fraction_12CH3_13COO", "formal_one_sided_95pct_upper",
    "q_standard_error_linearized", "rss", "residual_mad_sigma", "center_ppm",
    "best_fwhm_ppm", "eta", "n_points", "n_parameters", "relative_residual",
]
SETTING_COLUMNS = ["window_halfwidth_ppm", "j2_hz", "line_shape", "baseline_degree"]
EXPECTED_WINDOWS = (0.018, 0.025, 0.035, 0.050)
EXPECTED_COUPLINGS = (5.5, 5.914, 6.8)
EXPECTED_LINE_SHAPES = ("gaussian", "lorentzian", "pseudo_voigt", "voigt")
EXPECTED_BASELINES = (0, 1, 2)
NUMERIC_COLUMNS = [column for column in GRID_COLUMNS if column != "line_shape"]
NONNEGATIVE_COLUMNS = [
    "q_standard_error_linearized", "rss", "residual_mad_sigma", "relative_residual",
]
INTEGER_COLUMNS = ["n_points", "n_parameters"]
CONTRACT_RTOL = 1e-10
CONTRACT_ATOL = 1e-12


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


def _contract_check(name, passed, details=None):
    return {"name": name, "pass": bool(passed), "details": details or {}}


def _violating_rows(frame, mask):
    """Return setting keys for every row violating a contract predicate."""
    if not np.any(mask):
        return []
    return [
        {key: _json_value(row[key]) for key in SETTING_COLUMNS}
        for _, row in frame.loc[mask, SETTING_COLUMNS].iterrows()
    ]


def _contract_report(frame: pd.DataFrame, summary: dict):
    """Check the physical/structural contract independently of fit identity."""
    checks = []
    checks.append(_contract_check(
        "schema",
        list(frame.columns) == GRID_COLUMNS,
        {"actual": list(frame.columns), "expected": GRID_COLUMNS},
    ))
    required = set(GRID_COLUMNS)
    if not required.issubset(frame.columns):
        checks.append(_contract_check(
            "required_columns_present", False,
            {"missing": sorted(required.difference(frame.columns))},
        ))
        return {
            "pass": False,
            "checks": checks,
            "meaning": "structural/physical contract and internal summary arithmetic; not a fit-convergence or global-optimum claim",
        }

    expected_settings = {
        (window, coupling, shape, baseline)
        for window in EXPECTED_WINDOWS
        for coupling in EXPECTED_COUPLINGS
        for shape in EXPECTED_LINE_SHAPES
        for baseline in EXPECTED_BASELINES
    }
    observed_settings = set(frame[SETTING_COLUMNS].itertuples(index=False, name=None))
    missing_settings = sorted(expected_settings.difference(observed_settings), key=str)
    extra_settings = sorted(observed_settings.difference(expected_settings), key=str)
    checks.append(_contract_check(
        "exact_144_unique_settings",
        len(frame) == 144 and len(observed_settings) == 144 and not missing_settings and not extra_settings,
        {
            "rows": int(len(frame)),
            "unique_rows": int(len(observed_settings)),
            "missing": [_json_value(row) for row in missing_settings],
            "extra": [_json_value(row) for row in extra_settings],
        },
    ))

    numeric_ok = True
    nonnumeric = []
    nonfinite = {}
    for column in NUMERIC_COLUMNS:
        if not pd.api.types.is_numeric_dtype(frame[column]):
            numeric_ok = False
            nonnumeric.append(column)
            continue
        values = frame[column].to_numpy(dtype=float)
        finite = np.isfinite(values)
        if not np.all(finite):
            numeric_ok = False
            nonfinite[column] = np.flatnonzero(~finite).tolist()
    checks.append(_contract_check(
        "numeric_values_finite",
        numeric_ok,
        {"nonnumeric_columns": nonnumeric, "nonfinite_rows": nonfinite},
    ))

    def values(column):
        return frame[column].to_numpy(dtype=float)

    q = values("best_fraction_12CH3_13COO")
    upper = values("formal_one_sided_95pct_upper")
    q_bound_mask = (q < -CONTRACT_ATOL) | (q > 0.08 + CONTRACT_ATOL)
    upper_bound_mask = (upper < -CONTRACT_ATOL) | (upper > 0.08 + CONTRACT_ATOL)
    q_upper_mask = q > upper + CONTRACT_ATOL
    checks.append(_contract_check(
        "q_bounds",
        not np.any(q_bound_mask),
        {"range": [float(np.nanmin(q)), float(np.nanmax(q))], "violations": _violating_rows(frame, q_bound_mask)},
    ))
    checks.append(_contract_check(
        "formal_upper_bounds",
        not np.any(upper_bound_mask),
        {"range": [float(np.nanmin(upper)), float(np.nanmax(upper))], "violations": _violating_rows(frame, upper_bound_mask)},
    ))
    checks.append(_contract_check(
        "q_not_above_formal_upper",
        not np.any(q_upper_mask),
        {"violations": _violating_rows(frame, q_upper_mask)},
    ))

    for column in NONNEGATIVE_COLUMNS:
        data = values(column)
        mask = data < 0.0
        checks.append(_contract_check(
            f"nonnegative::{column}",
            not np.any(mask),
            {"minimum": float(np.nanmin(data)), "violations": _violating_rows(frame, mask)},
        ))

    fwhm = values("best_fwhm_ppm")
    fwhm_mask = (fwhm < 0.0004 - CONTRACT_ATOL) | (fwhm > 0.0200 + CONTRACT_ATOL)
    checks.append(_contract_check(
        "best_fwhm_range",
        not np.any(fwhm_mask),
        {"range": [float(np.nanmin(fwhm)), float(np.nanmax(fwhm))], "violations": _violating_rows(frame, fwhm_mask)},
    ))
    eta = values("eta")
    eta_mask = (eta < -CONTRACT_ATOL) | (eta > 1.0 + CONTRACT_ATOL)
    checks.append(_contract_check(
        "eta_range",
        not np.any(eta_mask),
        {"range": [float(np.nanmin(eta)), float(np.nanmax(eta))], "violations": _violating_rows(frame, eta_mask)},
    ))

    for column in INTEGER_COLUMNS:
        data = values(column)
        mask = data != np.floor(data)
        checks.append(_contract_check(
            f"integer::{column}",
            not np.any(mask),
            {"violations": _violating_rows(frame, mask)},
        ))
    npoints = values("n_points")
    nparameters = values("n_parameters")
    order_mask = (nparameters <= 0) | (npoints <= nparameters)
    checks.append(_contract_check(
        "n_points_gt_n_parameters_gt_zero",
        not np.any(order_mask),
        {"violations": _violating_rows(frame, order_mask)},
    ))

    residual_threshold = float(frame["relative_residual"].quantile(0.50))
    lower = frame.loc[frame["relative_residual"] <= residual_threshold]
    computed = {
        "n_nuisance_models": int(len(frame)),
        "all_models": {
            "best_fraction_range": [float(frame["best_fraction_12CH3_13COO"].min()), float(frame["best_fraction_12CH3_13COO"].max())],
            "formal_95pct_upper_range": [float(frame["formal_one_sided_95pct_upper"].min()), float(frame["formal_one_sided_95pct_upper"].max())],
        },
        "lower_residual_half_of_models": {
            "n": int(len(lower)),
            "median_best_fraction": float(lower["best_fraction_12CH3_13COO"].quantile(0.50)),
            "maximum_best_fraction": float(lower["best_fraction_12CH3_13COO"].max()),
            "median_formal_95pct_upper": float(lower["formal_one_sided_95pct_upper"].quantile(0.50)),
            "maximum_formal_95pct_upper": float(lower["formal_one_sided_95pct_upper"].max()),
        },
    }
    summary_fields = [
        ("n_nuisance_models", computed["n_nuisance_models"], summary.get("n_nuisance_models"), True),
        ("all_models.best_fraction_range", computed["all_models"]["best_fraction_range"], summary.get("all_models", {}).get("best_fraction_range"), False),
        ("all_models.formal_95pct_upper_range", computed["all_models"]["formal_95pct_upper_range"], summary.get("all_models", {}).get("formal_95pct_upper_range"), False),
        ("lower_residual_half_of_models.n", computed["lower_residual_half_of_models"]["n"], summary.get("lower_residual_half_of_models", {}).get("n"), True),
        ("lower_residual_half_of_models.median_best_fraction", computed["lower_residual_half_of_models"]["median_best_fraction"], summary.get("lower_residual_half_of_models", {}).get("median_best_fraction"), False),
        ("lower_residual_half_of_models.maximum_best_fraction", computed["lower_residual_half_of_models"]["maximum_best_fraction"], summary.get("lower_residual_half_of_models", {}).get("maximum_best_fraction"), False),
        ("lower_residual_half_of_models.median_formal_95pct_upper", computed["lower_residual_half_of_models"]["median_formal_95pct_upper"], summary.get("lower_residual_half_of_models", {}).get("median_formal_95pct_upper"), False),
        ("lower_residual_half_of_models.maximum_formal_95pct_upper", computed["lower_residual_half_of_models"]["maximum_formal_95pct_upper"], summary.get("lower_residual_half_of_models", {}).get("maximum_formal_95pct_upper"), False),
    ]
    summary_ok = True
    summary_details = {}
    for name, expected, reported, exact in summary_fields:
        if exact:
            passed = reported == expected
        else:
            if reported is None:
                passed = False
            else:
                passed = bool(np.allclose(reported, expected, rtol=CONTRACT_RTOL, atol=CONTRACT_ATOL, equal_nan=True))
        summary_ok = summary_ok and passed
        summary_details[name] = {"computed": expected, "reported": reported, "pass": bool(passed)}
    checks.append(_contract_check("summary_recomputed_from_csv", summary_ok, summary_details))
    return {
        "pass": all(check["pass"] for check in checks),
        "checks": checks,
        "computed_summary": computed,
        "meaning": "structural/physical contract and internal summary arithmetic; not a fit-convergence or global-optimum claim",
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
    p.add_argument("--contract-only",action="store_true",
                   help="Exit based on the physical/structural contract while still reporting strict regression differences")
    args=p.parse_args()
    grid_path = args.acetate/"acetate_carbonyl_grid.csv"
    summary_path = args.acetate/"acetate_carbonyl_summary.json"
    count, differences = _compare_table(grid_path,
        HERE/"reference_outputs/nmr/acetate_carbonyl_grid.csv",
        ["window_halfwidth_ppm","j2_hz","line_shape","baseline_degree"])
    grid = pd.read_csv(grid_path)
    actual=json.loads(summary_path.read_text())
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
    regression_report = {
        "acetate_models_compared": count,
        "regression": "passed" if not differences else "failed",
        "pass": not differences,
        "differences": differences,
        "meaning": "numerical reproduction, not an independently calibrated analytical limit",
    }
    contract = _contract_report(grid, actual)
    report = {
        "mode": "contract-only" if args.contract_only else "strict",
        "acetate_models_compared": count,
        "contract": contract,
        "regression": regression_report,
        "pass": contract["pass"] if args.contract_only else contract["pass"] and regression_report["pass"],
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["pass"]:
        raise AssertionError("Acetate validation failed; complete diagnostics were printed above")


if __name__=="__main__":main()
