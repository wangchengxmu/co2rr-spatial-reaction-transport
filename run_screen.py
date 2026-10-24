"""Run the portable species-transfer sensitivity screen.

The public repository keeps the implementation in ``model/species_transfer_model.py``.
This runner deliberately starts each case in a fresh Python subprocess so that a
failed case cannot contaminate a later case.  It does not change the model or
the archived calculation files.

Example
-------
    python run_screen.py --out runs/new_screen --workers 2

The output directory must be new or empty.  The optional ``--verification``
flag adds the six refinement/domain checks used by the numerical audit.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


DELTAS = (10, 30, 50)
CHIS = (0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005)
TIMEOUT_SECONDS = 1800


def _tag(value: float) -> str:
    """Return the stable case-name representation used by the archived screen."""

    return f"{value:g}".replace(".", "p")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _portable_reason(value: Any) -> str:
    """Keep a short diagnostic without embedding a local absolute path."""

    text = str(value or "model did not report a reason").replace("\r", " ").replace("\n", " ")
    text = re.sub(r"[A-Za-z]:[\\/][^\s,;)]*", "<local-path>", text)
    text = re.sub(r"(?<![A-Za-z0-9])/(?:[^\s,;)]+/?)+", "<local-path>", text)
    return text[:500]


def make_cases(include_verification: bool = False) -> list[dict[str, Any]]:
    """Return the 30 production cases and, optionally, six audit checks."""

    cases = [
        {
            "id": f"d{delta}_chi{_tag(chi)}",
            "delta": delta,
            "chi": chi,
            "nc": 48,
            "ny": 96,
            "nx": 24,
            "outer": 500,
            "kind": "production",
        }
        for delta in DELTAS
        for chi in CHIS
    ]
    if include_verification:
        for delta in DELTAS:
            cases.append(
                {
                    "id": f"d{delta}_chi0p05_fine",
                    "delta": delta,
                    "chi": 0.05,
                    "nc": 96,
                    "ny": 192,
                    "nx": 48,
                    "outer": 500,
                    "kind": "verification_fine",
                }
            )
            cases.append(
                {
                    "id": f"d{delta}_chi0p05_domain",
                    "delta": delta,
                    "chi": 0.05,
                    "nc": 48,
                    "ny": 112,
                    "nx": 24,
                    "outer": 1000,
                    "kind": "verification_domain",
                }
            )
    return cases


def _run_case(case: dict[str, Any], repo_root: Path, out_root: Path) -> dict[str, Any]:
    case_out = out_root / "cases" / str(case["id"])
    case_out.mkdir(parents=True, exist_ok=False)
    model = repo_root / "model" / "species_transfer_model.py"
    command = [
        sys.executable,
        "-X",
        "utf8",
        str(model),
        "--out",
        str(case_out.resolve()),
        "--delta",
        str(case["delta"]),
        "--chi",
        str(case["chi"]),
        "--nc",
        str(case["nc"]),
        "--ny",
        str(case["ny"]),
        "--nx",
        str(case["nx"]),
        "--outer",
        str(case["outer"]),
    ]
    environment = os.environ.copy()
    for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[key] = "1"
    started = time.monotonic()
    return_code: int | None = None
    timed_out = False
    result: dict[str, Any] = {}
    with (case_out / "stdout.log").open("w", encoding="utf-8") as log:
        try:
            process = subprocess.run(
                command,
                cwd=str(repo_root),
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=TIMEOUT_SECONDS,
                check=False,
            )
            return_code = process.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
    result_path = case_out / "result.json"
    if result_path.exists():
        try:
            parsed = json.loads(result_path.read_text(encoding="utf-8"))
            if isinstance(parsed, dict):
                result = parsed
        except (OSError, json.JSONDecodeError):
            result = {}
    accepted = bool(result.get("status") == "passed" and return_code == 0 and not timed_out)
    if timed_out:
        reason = f"worker timeout after {TIMEOUT_SECONDS} s"
    elif return_code != 0:
        reason = f"worker process exited with code {return_code}"
    elif result.get("status") != "passed":
        reason = _portable_reason(result.get("reason"))
    else:
        reason = "accepted"
    status = {
        "id": case["id"],
        "kind": case["kind"],
        "delta_um": case["delta"],
        "chi": case["chi"],
        "grid": {"n_catalyst": case["nc"], "n_external": case["ny"], "n_axial": case["nx"], "outer_um": case["outer"]},
        "status": "passed" if accepted else "not_accepted",
        "model_status": result.get("status", "missing"),
        "return_code": return_code,
        "elapsed_s": round(time.monotonic() - started, 3),
        "reason": reason,
    }
    _write_json(case_out / "portable_status.json", status)
    return status


def _empty_or_missing(path: Path) -> bool:
    return not path.exists() or not any(path.iterdir())


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    model_path = repo_root / "model" / "species_transfer_model.py"
    if not model_path.is_file():
        raise FileNotFoundError("Expected model/species_transfer_model.py in the repository")
    out_root = Path(args.out)
    if not out_root.is_absolute():
        out_root = (Path.cwd() / out_root).resolve()
    if not _empty_or_missing(out_root):
        raise FileExistsError(f"Output directory is not empty; choose a new directory: {out_root}")
    out_root.mkdir(parents=True, exist_ok=True)
    cases = make_cases(args.verification)
    manifest = {
        "schema": "co2rr-spatial-screen-v1",
        "created_utc": _now(),
        "model": "model/species_transfer_model.py",
        "production_case_count": 30,
        "verification_case_count": len(cases) - 30,
        "cases": cases,
        "timeout_s": TIMEOUT_SECONDS,
        "workers": args.workers,
    }
    _write_json(out_root / "screen_manifest.json", manifest)
    pending: list[dict[str, Any]] = cases[: args.limit] if args.limit is not None else list(cases)
    statuses: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        running = {pool.submit(_run_case, case, repo_root, out_root): case for case in pending}
        while running:
            done, _ = wait(running, return_when=FIRST_COMPLETED)
            for future in done:
                case = running.pop(future)
                try:
                    status = future.result()
                except Exception as exc:  # preserve a portable failure record
                    status = {
                        "id": case["id"],
                        "kind": case["kind"],
                        "delta_um": case["delta"],
                        "chi": case["chi"],
                        "status": "not_accepted",
                        "reason": _portable_reason(exc),
                    }
                    _write_json(out_root / "cases" / str(case["id"]) / "portable_status.json", status)
                statuses.append(status)
                print(json.dumps({"id": status["id"], "status": status["status"], "elapsed_s": status.get("elapsed_s")}), flush=True)
    statuses.sort(key=lambda row: str(row["id"]))
    summary = {
        "schema": "co2rr-spatial-screen-summary-v1",
        "created_utc": _now(),
        "requested_cases": len(pending),
        "passed_cases": sum(row.get("status") == "passed" for row in statuses),
        "not_accepted_cases": sum(row.get("status") != "passed" for row in statuses),
        "complete": len(statuses) == len(pending) and all(row.get("status") == "passed" for row in statuses),
        "cases": statuses,
    }
    _write_json(out_root / "screen_summary.json", summary)
    return summary


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run 30 CO2RR species-transfer cases in isolated subprocesses. Output must be new or empty.",
        epilog="The model is expected at model/species_transfer_model.py; use --verification for six refinement/domain checks.",
    )
    parser.add_argument("--out", required=True, help="New or empty output directory (for example runs/new_screen)")
    parser.add_argument("--workers", type=int, default=2, help="Maximum concurrent workers (default: 2)")
    parser.add_argument("--limit", type=int, help="Run only the first N cases, useful for a smoke test")
    parser.add_argument("--verification", action="store_true", help="Append the six fine-grid/domain verification cases")
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be at least 1")
    return args


if __name__ == "__main__":
    try:
        run(parse_args())
    except Exception as exc:
        print(f"run_screen.py: {exc}", file=sys.stderr)
        raise SystemExit(2)
