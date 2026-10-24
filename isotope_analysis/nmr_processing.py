"""Portable one-dimensional proton-NMR preprocessing.

The repository stores a complex FID and the small acquisition dictionary
needed by the Bruker/nmrglue processing steps.  This module deliberately does
not read or write data at import time; callers choose a condition and invoke
``process_1h`` explicitly.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def _nmrglue():
    """Import nmrglue only when processing is requested."""
    try:
        import nmrglue as ng
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            "NMR processing requires nmrglue==0.11; install requirements.txt"
        ) from exc
    return ng


def asls(
    y: np.ndarray,
    lam: float = 2e7,
    p: float = 0.001,
    iterations: int = 12,
) -> np.ndarray:
    """Return an asymmetric least-squares baseline (the original settings)."""
    y = np.asarray(y, dtype=float)
    if y.ndim != 1 or len(y) < 4:
        raise ValueError("ASLS requires a one-dimensional signal with at least 4 points")
    n = len(y)
    diff = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n))
    weights = np.ones(n)
    for _ in range(iterations):
        system = sparse.spdiags(weights, 0, n, n) + lam * (diff.T @ diff)
        baseline = spsolve(system, weights * y)
        weights = p * (y > baseline) + (1 - p) * (y <= baseline)
    return baseline


def _load_acquisition(input_dir: Path) -> dict[str, Any]:
    path = Path(input_dir) / "acquisition.json"
    if not path.is_file():
        raise FileNotFoundError(f"Missing NMR acquisition manifest: {path}")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("acquisition.json must contain a condition mapping")
    return manifest


def load_fid(condition: str, input_dir: Path) -> tuple[dict[str, Any], np.ndarray]:
    """Load one portable condition and its acquisition dictionary.

    ``input_dir/acquisition.json`` maps each condition to a dictionary with
    ``acqus`` and ``procs`` entries.  The corresponding FID is read from
    ``input_dir/{condition}_fid.npz`` under the single array key ``fid``.
    """
    input_dir = Path(input_dir)
    manifest = _load_acquisition(input_dir)
    if condition not in manifest:
        raise KeyError(f"Condition {condition!r} is absent from acquisition.json")
    dic = manifest[condition]
    if not isinstance(dic, dict) or not isinstance(dic.get("acqus"), dict):
        raise ValueError(f"Acquisition entry for {condition!r} must contain an acqus dictionary")
    if not isinstance(dic.get("procs", {}), dict):
        raise ValueError(f"Acquisition entry for {condition!r} must contain a procs dictionary")
    path = input_dir / f"{condition}_fid.npz"
    if not path.is_file():
        raise FileNotFoundError(f"Missing portable NMR FID: {path}")
    with np.load(path, allow_pickle=False) as arrays:
        if set(arrays.files) != {"fid"}:
            raise ValueError(f"{path.name} must contain only the complex array 'fid'")
        fid = np.asarray(arrays["fid"])
    if fid.ndim != 1 or fid.size == 0 or not np.iscomplexobj(fid):
        raise ValueError(f"{path.name}: fid must be a non-empty one-dimensional complex array")
    if not np.isfinite(fid.real).all() or not np.isfinite(fid.imag).all():
        raise ValueError(f"{path.name}: fid contains non-finite values")
    # Make a detached JSON-like dictionary so nmrglue can safely add fields.
    dic = json.loads(json.dumps(dic))
    return dic, fid.astype(np.complex128, copy=False)


def process_1h(
    condition: str,
    input_dir: Path,
    line_broadening_hz: float = 0.3,
    size: int = 32768,
    low_ppm: float = 0.55,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Process a stored FID with the frozen one-dimensional workflow.

    The operation order follows the original analysis: digital-filter removal,
    exponential apodization, zero filling, FFT, ACME autophasing, Bruker
    chemical-shift calibration, and ASLS baseline subtraction.  The returned
    ppm axis keeps nmrglue's original descending ordering.
    """
    if line_broadening_hz < 0:
        raise ValueError("line_broadening_hz must be non-negative")
    if size < 1:
        raise ValueError("size must be positive")
    ng = _nmrglue()
    dic, fid = load_fid(condition, input_dir)
    fid = ng.bruker.remove_digital_filter(dic, fid)
    try:
        spectral_width_hz = float(dic["acqus"]["SW_h"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("acquisition.json must provide numeric acqus.SW_h") from exc
    if spectral_width_hz <= 0:
        raise ValueError("acqus.SW_h must be positive")
    dwell = 1.0 / spectral_width_hz
    fid = fid * np.exp(-np.pi * float(line_broadening_hz) * np.arange(len(fid)) * dwell)
    fid = ng.proc_base.zf_size(fid, int(size))
    spectrum = ng.proc_base.rev(ng.proc_base.fft(fid))
    spectrum = ng.proc_autophase.autops(spectrum, "acme", disp=False)
    real = np.real(spectrum)
    udic = ng.bruker.guess_udic(dic, spectrum)
    udic[0]["size"] = int(size)
    ppm = ng.fileiobase.uc_from_udic(udic, 0).ppm_scale()
    high_ppm = 9.1
    mask = (ppm >= float(low_ppm)) & (ppm <= high_ppm)
    corrected = real[mask] - asls(real[mask])
    return ppm[mask], corrected, dic
