"""Contract and strict-regression tests for the acetate validator."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REFERENCE_GRID = HERE / "reference_outputs" / "nmr" / "acetate_carbonyl_grid.csv"
REFERENCE_SUMMARY = HERE / "reference_outputs" / "nmr" / "acetate_carbonyl_summary.json"


class ValidationContractTests(unittest.TestCase):
    def _copy_reference(self, directory: Path) -> None:
        shutil.copy2(REFERENCE_GRID, directory / "acetate_carbonyl_grid.csv")
        shutil.copy2(REFERENCE_SUMMARY, directory / "acetate_carbonyl_summary.json")

    def _run(self, directory: Path, contract_only: bool = True) -> tuple[subprocess.CompletedProcess[str], dict]:
        command = [
            sys.executable, "-m", "isotope_analysis.validate_analyses",
            "--acetate", str(directory),
        ]
        if contract_only:
            command.append("--contract-only")
        result = subprocess.run(
            command,
            cwd=HERE.parent,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertTrue(result.stdout.strip(), result.stderr)
        return result, json.loads(result.stdout)

    def _mutate_grid(self, directory: Path, **updates: object) -> None:
        grid = pd.read_csv(directory / "acetate_carbonyl_grid.csv")
        for column, value in updates.items():
            grid.loc[0, column] = value
        grid.to_csv(directory / "acetate_carbonyl_grid.csv", index=False)

    def test_reference_passes_contract(self) -> None:
        with tempfile.TemporaryDirectory(prefix="co2rr_validator_contract_") as raw:
            directory = Path(raw)
            self._copy_reference(directory)
            result, report = self._run(directory)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(report["contract"]["pass"])
        self.assertTrue(report["pass"])

    def test_negative_q_fails_contract(self) -> None:
        with tempfile.TemporaryDirectory(prefix="co2rr_validator_contract_") as raw:
            directory = Path(raw)
            self._copy_reference(directory)
            self._mutate_grid(directory, best_fraction_12CH3_13COO=-0.001)
            result, report = self._run(directory)
        self.assertNotEqual(result.returncode, 0)
        checks = {check["name"]: check for check in report["contract"]["checks"]}
        self.assertFalse(checks["q_bounds"]["pass"])

    def test_formal_upper_below_q_fails_contract(self) -> None:
        with tempfile.TemporaryDirectory(prefix="co2rr_validator_contract_") as raw:
            directory = Path(raw)
            self._copy_reference(directory)
            self._mutate_grid(
                directory,
                best_fraction_12CH3_13COO=0.05,
                formal_one_sided_95pct_upper=0.01,
            )
            result, report = self._run(directory)
        self.assertNotEqual(result.returncode, 0)
        checks = {check["name"]: check for check in report["contract"]["checks"]}
        self.assertFalse(checks["q_not_above_formal_upper"]["pass"])

    def test_wrong_summary_fails_contract(self) -> None:
        with tempfile.TemporaryDirectory(prefix="co2rr_validator_contract_") as raw:
            directory = Path(raw)
            self._copy_reference(directory)
            summary_path = directory / "acetate_carbonyl_summary.json"
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            summary["lower_residual_half_of_models"]["median_best_fraction"] = 0.123
            summary_path.write_text(json.dumps(summary), encoding="utf-8")
            result, report = self._run(directory)
        self.assertNotEqual(result.returncode, 0)
        checks = {check["name"]: check for check in report["contract"]["checks"]}
        self.assertFalse(checks["summary_recomputed_from_csv"]["pass"])

    def test_unused_eta_change_contract_only_passes_but_strict_fails(self) -> None:
        with tempfile.TemporaryDirectory(prefix="co2rr_validator_contract_") as raw:
            directory = Path(raw)
            self._copy_reference(directory)
            self._mutate_grid(directory, eta=0.75)
            contract_result, contract_report = self._run(directory, contract_only=True)
            strict_result, strict_report = self._run(directory, contract_only=False)
        self.assertEqual(contract_result.returncode, 0, contract_result.stderr)
        self.assertTrue(contract_report["contract"]["pass"])
        self.assertFalse(contract_report["regression"]["pass"])
        self.assertNotEqual(strict_result.returncode, 0)
        self.assertTrue(strict_report["contract"]["pass"])
        self.assertFalse(strict_report["regression"]["pass"])
        self.assertTrue(any(item["column"] == "eta" for item in strict_report["regression"]["differences"]))


if __name__ == "__main__":
    unittest.main()
