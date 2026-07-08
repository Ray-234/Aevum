"""Equivalence checks for P107 terminal audit outputs.

Performance passes should be able to prove that they did not change generated
terrain.  This module compares the terminal array archive and compact terminal
metrics from two existing P107 runs while ignoring profile-only metadata by
default.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np


P107_EQUIVALENCE_SCHEMA = "aevum.p107_equivalence.v1"
DEFAULT_METRIC_SKIP_KEYS = (
    "terrain_internal_profile",
    "assets",
    "array_archive",
)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def resolve_p107_run_dir(path: str | Path) -> Path:
    """Resolve a P107 root/run/metrics path to a single terminal run directory."""
    source = Path(path)
    if source.is_file():
        if source.name == "p107_terminal_metrics.json":
            return source.parent
        if source.name == "p107_audit_summary.json":
            return resolve_p107_run_dir(source.parent)
        raise ValueError(f"unsupported P107 input file: {source}")
    if not source.exists():
        raise ValueError(f"P107 input does not exist: {source}")
    if (
        (source / "p107_terminal_arrays.npz").is_file()
        and (source / "p107_terminal_metrics.json").is_file()
    ):
        return source
    candidates = sorted(
        candidate.parent
        for candidate in source.glob("*/p107_terminal_arrays.npz")
        if (candidate.parent / "p107_terminal_metrics.json").is_file()
    )
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise ValueError(f"no P107 terminal run found under: {source}")
    raise ValueError(
        f"multiple P107 terminal runs found under {source}; pass one run dir")


def _array_max_abs_diff(a: np.ndarray, b: np.ndarray) -> float | None:
    if a.shape != b.shape:
        return None
    if not (np.issubdtype(a.dtype, np.number) and np.issubdtype(b.dtype, np.number)):
        return None
    try:
        diff = np.abs(a.astype(np.float64) - b.astype(np.float64))
    except (TypeError, ValueError):
        return None
    finite = diff[np.isfinite(diff)]
    if finite.size == 0:
        return 0.0
    return float(np.max(finite))


def compare_p107_arrays(baseline_run_dir: str | Path,
                        candidate_run_dir: str | Path) -> dict[str, Any]:
    """Compare two P107 terminal array archives exactly."""
    baseline_path = Path(baseline_run_dir) / "p107_terminal_arrays.npz"
    candidate_path = Path(candidate_run_dir) / "p107_terminal_arrays.npz"
    if not baseline_path.is_file():
        raise ValueError(f"missing baseline array archive: {baseline_path}")
    if not candidate_path.is_file():
        raise ValueError(f"missing candidate array archive: {candidate_path}")

    changed: list[dict[str, Any]] = []
    dtype_mismatches: list[dict[str, Any]] = []
    with np.load(baseline_path) as baseline, np.load(candidate_path) as candidate:
        baseline_keys = set(baseline.files)
        candidate_keys = set(candidate.files)
        common = sorted(baseline_keys & candidate_keys)
        for key in common:
            left = baseline[key]
            right = candidate[key]
            if left.dtype != right.dtype:
                dtype_mismatches.append({
                    "key": key,
                    "baseline_dtype": str(left.dtype),
                    "candidate_dtype": str(right.dtype),
                })
            if left.shape != right.shape or not np.array_equal(
                left, right, equal_nan=True,
            ):
                changed.append({
                    "key": key,
                    "baseline_shape": list(left.shape),
                    "candidate_shape": list(right.shape),
                    "baseline_dtype": str(left.dtype),
                    "candidate_dtype": str(right.dtype),
                    "max_abs_diff": _array_max_abs_diff(left, right),
                })
    return {
        "common_count": len(common),
        "missing_keys": sorted(baseline_keys - candidate_keys),
        "extra_keys": sorted(candidate_keys - baseline_keys),
        "dtype_mismatches": dtype_mismatches,
        "changed": changed,
        "changed_count": len(changed),
        "equivalent": (
            not changed
            and not dtype_mismatches
            and baseline_keys == candidate_keys
        ),
    }


def _metric_equal(left: Any, right: Any, *, float_atol: float) -> tuple[bool, float]:
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        left_float = float(left)
        right_float = float(right)
        if math.isnan(left_float) and math.isnan(right_float):
            return True, 0.0
        diff = abs(left_float - right_float)
        return diff <= float_atol, diff
    return left == right, 0.0


def _compare_metric_node(left: Any,
                         right: Any,
                         *,
                         path: str,
                         skip_top_keys: set[str],
                         float_atol: float,
                         changed: list[dict[str, Any]],
                         max_float: dict[str, Any]) -> None:
    if path == "" and isinstance(left, dict) and isinstance(right, dict):
        keys = sorted((set(left) | set(right)) - skip_top_keys)
    elif isinstance(left, dict) and isinstance(right, dict):
        keys = sorted(set(left) | set(right))
    else:
        keys = []

    if keys:
        for key in keys:
            child_path = f"{path}.{key}" if path else str(key)
            if key not in left:
                changed.append({"path": child_path, "kind": "missing_baseline"})
                continue
            if key not in right:
                changed.append({"path": child_path, "kind": "missing_candidate"})
                continue
            _compare_metric_node(
                left[key],
                right[key],
                path=child_path,
                skip_top_keys=skip_top_keys,
                float_atol=float_atol,
                changed=changed,
                max_float=max_float,
            )
        return

    if isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            changed.append({
                "path": path,
                "kind": "list_length",
                "baseline": len(left),
                "candidate": len(right),
            })
            return
        for index, (left_item, right_item) in enumerate(zip(left, right)):
            _compare_metric_node(
                left_item,
                right_item,
                path=f"{path}[{index}]",
                skip_top_keys=skip_top_keys,
                float_atol=float_atol,
                changed=changed,
                max_float=max_float,
            )
        return

    equal, diff = _metric_equal(left, right, float_atol=float_atol)
    if diff > float(max_float["diff"]):
        max_float["diff"] = float(diff)
        max_float["path"] = path
    if not equal:
        changed.append({
            "path": path,
            "kind": "value",
            "baseline": left,
            "candidate": right,
            "abs_diff": diff if diff else None,
        })


def compare_p107_metrics(baseline_run_dir: str | Path,
                         candidate_run_dir: str | Path,
                         *,
                         skip_top_keys: tuple[str, ...] = DEFAULT_METRIC_SKIP_KEYS,
                         float_atol: float = 1e-12) -> dict[str, Any]:
    """Compare compact P107 terminal metrics, skipping profile-only top keys."""
    baseline_path = Path(baseline_run_dir) / "p107_terminal_metrics.json"
    candidate_path = Path(candidate_run_dir) / "p107_terminal_metrics.json"
    if not baseline_path.is_file():
        raise ValueError(f"missing baseline metrics: {baseline_path}")
    if not candidate_path.is_file():
        raise ValueError(f"missing candidate metrics: {candidate_path}")
    baseline = json.loads(baseline_path.read_text())
    candidate = json.loads(candidate_path.read_text())
    changed: list[dict[str, Any]] = []
    max_float = {"diff": 0.0, "path": ""}
    _compare_metric_node(
        baseline,
        candidate,
        path="",
        skip_top_keys=set(skip_top_keys),
        float_atol=float(float_atol),
        changed=changed,
        max_float=max_float,
    )
    return {
        "skip_top_keys": list(skip_top_keys),
        "float_atol": float(float_atol),
        "changed": changed,
        "changed_count": len(changed),
        "max_float_diff": max_float,
        "equivalent": not changed,
    }


def compare_p107_outputs(baseline: str | Path,
                         candidate: str | Path,
                         *,
                         skip_top_keys: tuple[str, ...] = DEFAULT_METRIC_SKIP_KEYS,
                         float_atol: float = 1e-12) -> dict[str, Any]:
    """Compare two existing P107 outputs and return a gate report."""
    baseline_run_dir = resolve_p107_run_dir(baseline)
    candidate_run_dir = resolve_p107_run_dir(candidate)
    arrays = compare_p107_arrays(baseline_run_dir, candidate_run_dir)
    metrics = compare_p107_metrics(
        baseline_run_dir,
        candidate_run_dir,
        skip_top_keys=skip_top_keys,
        float_atol=float_atol,
    )
    return {
        "schema": P107_EQUIVALENCE_SCHEMA,
        "baseline": {
            "input": str(baseline),
            "run_dir": str(baseline_run_dir),
        },
        "candidate": {
            "input": str(candidate),
            "run_dir": str(candidate_run_dir),
        },
        "arrays": arrays,
        "metrics": metrics,
        "equivalent": bool(arrays["equivalent"] and metrics["equivalent"]),
    }


def write_p107_equivalence_report(report: dict[str, Any],
                                  out_path: str | Path) -> Path:
    """Write a P107 equivalence report as formatted JSON."""
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, default=_json_default) + "\n")
    return path
