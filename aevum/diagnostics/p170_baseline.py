"""P170 six-world historical geomorphology baseline runner."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
from pathlib import Path
import time
from typing import Any

import numpy as np

from aevum.diagnostics.historical_geomorphology import (
    historical_geomorphology_summary,
    write_historical_geomorphology_audit,
)
from aevum.engine import Engine
from aevum.spec.presets import get_preset


SCHEMA = "aevum.p170_six_world_baseline.v1"

DEFAULT_JOBS: tuple[tuple[str, str, int], ...] = (
    ("waterworld", "waterworld_seed7", 7),
    ("waterworld", "waterworld_seed707", 707),
    ("earthlike", "earthlike_seed42", 42),
    ("earthlike", "earthlike_seed909", 909),
    ("arid", "arid_seed101", 101),
    ("arid", "arid_seed1001", 1001),
)

DEFAULT_ACTIVE_MODULES = frozenset({
    "stellar",
    "interior",
    "impacts",
    "tectonics",
    "terrain",
})

TARGET_TIMES_MYR = (500.0, 1500.0, 2500.0, 3500.0, 4500.0)


def run_p170_six_world_baseline(
    outdir: str | Path,
    *,
    cells: int = 8000,
    frames: int = 90,
    max_workers: int = 5,
    jobs: tuple[tuple[str, str, int], ...] = DEFAULT_JOBS,
    active_modules: frozenset[str] = DEFAULT_ACTIVE_MODULES,
) -> dict[str, Any]:
    """Run the planned P170 six-world baseline and write an aggregate summary."""
    _set_process_thread_env()
    root = Path(outdir)
    root.mkdir(parents=True, exist_ok=True)
    config = {
        "cells": int(cells),
        "frames": int(frames),
        "max_workers": int(max_workers),
        "jobs": [
            {"preset": preset, "label": label, "seed": int(seed)}
            for preset, label, seed in jobs
        ],
        "active_modules": sorted(active_modules),
        "target_times_myr": list(TARGET_TIMES_MYR),
    }
    (root / "p170_baseline_config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    started = time.time()
    rows: list[dict[str, Any]] = []
    worker_count = max(1, min(int(max_workers), len(jobs)))
    with ProcessPoolExecutor(max_workers=worker_count) as pool:
        futures = [
            pool.submit(
                _run_one,
                root,
                preset,
                label,
                int(seed),
                int(cells),
                int(frames),
                tuple(sorted(active_modules)),
            )
            for preset, label, seed in jobs
        ]
        for future in as_completed(futures):
            rows.append(future.result())
            rows.sort(key=lambda row: str(row["label"]))
            (root / "p170_baseline_partial_summary.json").write_text(
                json.dumps(_aggregate_rows(rows, config, started), indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

    summary = _aggregate_rows(rows, config, started)
    (root / "p170_six_world_baseline_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def _run_one(
    root: Path,
    preset_name: str,
    label: str,
    seed: int,
    cells: int,
    frames: int,
    active_modules: tuple[str, ...],
) -> dict[str, Any]:
    _set_process_thread_env()
    started = time.time()
    outdir = root / label
    outdir.mkdir(parents=True, exist_ok=True)

    spec = get_preset(preset_name)
    spec.seed = int(seed)
    spec.grid_cells = int(cells)
    spec.t_end_myr = 4500.0
    eng = Engine.build(spec)
    active = set(active_modules)
    eng.scheduler.modules = [
        sm for sm in eng.scheduler.modules if sm.module.name in active
    ]
    eng.run(n_frames=int(frames), progress=False)

    p170 = write_historical_geomorphology_audit(eng.world, eng.archive, outdir)
    compact = _compact_world_summary(label, preset_name, seed, cells, frames, p170)
    compact["runtime_seconds"] = float(time.time() - started)
    compact["active_modules"] = list(active_modules)
    compact["p170_audit_path"] = str(outdir / "p170_historical_geomorphology_audit.json")
    (outdir / "summary.json").write_text(
        json.dumps(compact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return compact


def _compact_world_summary(
    label: str,
    preset_name: str,
    seed: int,
    cells: int,
    frames: int,
    p170: dict[str, Any],
) -> dict[str, Any]:
    frame_rows = list(p170.get("frame_rows", []) or [])
    flag_summary = dict(p170.get("flag_summary", {}) or {})
    target_rows = [_nearest_frame_row(frame_rows, target) for target in TARGET_TIMES_MYR]
    return {
        "schema": "aevum.p170_world_baseline_row.v1",
        "label": str(label),
        "preset": str(preset_name),
        "seed": int(seed),
        "cells": int(cells),
        "requested_frames": int(frames),
        "frame_count": int(p170.get("frame_count", 0)),
        "usable_frame_count": int(p170.get("usable_frame_count", 0)),
        "required_metric_keys_present": bool(p170.get("required_metric_keys_present", False)),
        "missing_field_counts": dict(p170.get("missing_field_counts", {}) or {}),
        "flag_summary": flag_summary,
        "metric_extremes": p170.get("metric_extremes", {}),
        "target_frame_rows": target_rows,
        "land_baseline_issue": bool(flag_summary.get("ordinary_plateau_frame_count", 0) > 0),
        "ocean_baseline_issue": bool(flag_summary.get("ordinary_deep_ocean_frame_count", 0) > 0),
    }


def _nearest_frame_row(rows: list[dict[str, Any]], target_myr: float) -> dict[str, Any]:
    if not rows:
        return {
            "target_myr": float(target_myr),
            "available": False,
            "time_myr": 0.0,
            "diagnostic_flags": {},
            "land_metrics": {},
            "ocean_metrics": {},
        }
    row = min(rows, key=lambda item: abs(float(item.get("time_myr", 0.0)) - float(target_myr)))
    return {
        "target_myr": float(target_myr),
        "available": True,
        "time_myr": float(row.get("time_myr", 0.0)),
        "diagnostic_flags": dict(row.get("diagnostic_flags", {}) or {}),
        "land_metrics": _pick_metrics(
            row.get("land_metrics", {}),
            (
                "inland_detail_entropy",
                "ordinary_plateau_fraction",
                "broad_flat_inland_component_count",
                "continent_province_count_p50",
                "old_orogen_expression_fraction",
                "rift_basin_expression_fraction",
            ),
        ),
        "ocean_metrics": _pick_metrics(
            row.get("ocean_metrics", {}),
            (
                "ocean_fabric_entropy",
                "ridge_visible_fraction",
                "ridge_age_symmetry_score",
                "fracture_zone_length_fraction",
                "abyssal_plain_fraction",
                "hotspot_track_count",
                "seamount_chain_count",
                "unparented_shoal_fraction",
            ),
        ),
    }


def _pick_metrics(metrics: Any, keys: tuple[str, ...]) -> dict[str, float | int]:
    source = metrics if isinstance(metrics, dict) else {}
    out: dict[str, float | int] = {}
    for key in keys:
        value = source.get(key, 0.0)
        if isinstance(value, (int, np.integer)):
            out[key] = int(value)
        else:
            out[key] = float(value)
    return out


def _aggregate_rows(
    rows: list[dict[str, Any]],
    config: dict[str, Any],
    started: float,
) -> dict[str, Any]:
    rows = sorted(rows, key=lambda row: str(row["label"]))
    land_issue_rows = [row for row in rows if bool(row.get("land_baseline_issue", False))]
    ocean_issue_rows = [row for row in rows if bool(row.get("ocean_baseline_issue", False))]
    terminal_land_jumps = [
        float(row.get("flag_summary", {}).get("land_terminal_quality_jump", 0.0))
        for row in rows
    ]
    terminal_ocean_jumps = [
        float(row.get("flag_summary", {}).get("ocean_terminal_quality_jump", 0.0))
        for row in rows
    ]
    required_ok = all(bool(row.get("required_metric_keys_present", False)) for row in rows)
    return {
        "schema": SCHEMA,
        "config": config,
        "runtime_seconds": float(time.time() - started),
        "world_count": int(len(rows)),
        "completed_labels": [str(row["label"]) for row in rows],
        "worlds": rows,
        "aggregate": {
            "land_issue_world_count": int(len(land_issue_rows)),
            "ocean_issue_world_count": int(len(ocean_issue_rows)),
            "land_issue_labels": [str(row["label"]) for row in land_issue_rows],
            "ocean_issue_labels": [str(row["label"]) for row in ocean_issue_rows],
            "max_land_terminal_quality_jump": (
                float(max(terminal_land_jumps)) if terminal_land_jumps else 0.0
            ),
            "max_ocean_terminal_quality_jump": (
                float(max(terminal_ocean_jumps)) if terminal_ocean_jumps else 0.0
            ),
            "median_land_terminal_quality_jump": (
                float(np.median(np.asarray(terminal_land_jumps, dtype=np.float64)))
                if terminal_land_jumps else 0.0
            ),
            "median_ocean_terminal_quality_jump": (
                float(np.median(np.asarray(terminal_ocean_jumps, dtype=np.float64)))
                if terminal_ocean_jumps else 0.0
            ),
        },
        "acceptance": {
            "baseline_completed": bool(len(rows) == len(config.get("jobs", []))),
            "required_metric_keys_present": bool(required_ok),
            "generation_behavior_changed": False,
            "six_world_8000_baseline_completed": bool(
                len(rows) == 6 and int(config.get("cells", 0)) == 8000
            ),
        },
    }


def _set_process_thread_env() -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="python -m aevum.diagnostics.p170_baseline",
        description="Run the P170 six-world historical geomorphology baseline.",
    )
    parser.add_argument("--out", default="out_p170_six_world_baseline")
    parser.add_argument("--cells", type=int, default=8000)
    parser.add_argument("--frames", type=int, default=90)
    parser.add_argument("--max-workers", type=int, default=5)
    args = parser.parse_args(argv)
    summary = run_p170_six_world_baseline(
        args.out,
        cells=int(args.cells),
        frames=int(args.frames),
        max_workers=int(args.max_workers),
    )
    print("== aevum :: P170 six-world historical geomorphology baseline ==")
    print(f"   worlds: {summary['world_count']}")
    print(f"   land issue worlds: {summary['aggregate']['land_issue_world_count']}")
    print(f"   ocean issue worlds: {summary['aggregate']['ocean_issue_world_count']}")
    print(f"   wrote {Path(args.out) / 'p170_six_world_baseline_summary.json'}")


if __name__ == "__main__":
    main()
