"""Focused regression tests for the portable NMR template workflow."""
from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .nmr_templates import analyze


HERE = Path(__file__).resolve().parent
INPUT_DIR = HERE / "inputs" / "nmr"
REFERENCE = HERE / "reference_outputs" / "nmr" / "template_fit_summary.csv"
# Nonlinear minimizers need unit-aware tolerances rather than bit identity.
REGRESSION_RTOL = 2e-6
REGRESSION_ATOL = 2e-10
# 1e-5 ppm is 0.005 Hz at 500 MHz, about 1/50 of the 0.245 Hz FFT spacing.
# Weak-control satellite position is less identifiable than isotope fraction.
PPM_ATOL = 1e-5
class TemplateFitRegression(unittest.TestCase):
    def test_principal_summary_matches_frozen_output(self) -> None:
        with tempfile.TemporaryDirectory(prefix="co2rr_nmr_test_") as directory:
            analyze(INPUT_DIR, Path(directory))
            with (Path(directory) / "template_fit_summary.csv").open(newline="", encoding="utf-8-sig") as handle:
                observed = list(csv.DictReader(handle))
            with REFERENCE.open(newline="", encoding="utf-8-sig") as handle:
                expected = list(csv.DictReader(handle))
        self.assertEqual(len(observed), len(expected))
        self.assertEqual([row.get("condition") for row in observed], [row.get("condition") for row in expected])
        self.assertEqual([row.get("group") for row in observed], [row.get("group") for row in expected])
        self.assertEqual(set(observed[0]), set(expected[0]))
        differences = []
        for actual, reference in zip(observed, expected):
            for key in actual:
                if key == "nuisance_areas":
                    if actual[key] == "" or reference[key] == "":
                        self.assertEqual(actual[key], reference[key], key)
                    else:
                        np.testing.assert_allclose(json.loads(actual[key]), json.loads(reference[key]), rtol=REGRESSION_RTOL, atol=REGRESSION_ATOL)
                elif key in {"condition", "group"}:
                    self.assertEqual(actual[key], reference[key], key)
                else:
                    if actual[key] == "" or reference[key] == "":
                        self.assertEqual(actual[key], reference[key], key)
                        continue
                    atol = PPM_ATOL if key.endswith("_ppm") else REGRESSION_ATOL
                    if not np.isclose(float(actual[key]), float(reference[key]), rtol=REGRESSION_RTOL, atol=atol):
                        differences.append({"condition": actual["condition"], "group": actual["group"],
                            "field": key, "actual": float(actual[key]), "reference": float(reference[key]),
                            "absolute_difference": abs(float(actual[key])-float(reference[key]))})
        self.assertFalse(differences, json.dumps(differences, indent=2))


if __name__ == "__main__":
    unittest.main()
