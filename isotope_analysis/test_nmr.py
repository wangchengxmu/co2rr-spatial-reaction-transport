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
        for actual, reference in zip(observed, expected):
            for key in actual:
                if key == "nuisance_areas":
                    if actual[key] == "" or reference[key] == "":
                        self.assertEqual(actual[key], reference[key], key)
                    else:
                        np.testing.assert_allclose(json.loads(actual[key]), json.loads(reference[key]), rtol=2e-7, atol=2e-10)
                elif key in {"condition", "group"}:
                    self.assertEqual(actual[key], reference[key], key)
                else:
                    if actual[key] == "" or reference[key] == "":
                        self.assertEqual(actual[key], reference[key], key)
                        continue
                    self.assertTrue(np.isclose(float(actual[key]), float(reference[key]), rtol=2e-7, atol=2e-10), key)


if __name__ == "__main__":
    unittest.main()
