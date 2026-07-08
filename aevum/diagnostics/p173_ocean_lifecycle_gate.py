"""P173 multi-world ocean-floor lifecycle gate."""
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
    write_historical_geomorphology_audit,
)
from aevum.diagnostics.historical_objects import write_historical_object_audit
from aevum.diagnostics.p173_unsupported_shoal_attribution import (
    write_p173_unsupported_shoal_attribution,
)
from aevum.engine import Engine
from aevum.spec.presets import get_preset


SCHEMA = "aevum.p173_ocean_lifecycle_gate.v1"

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

OCEAN_LIFECYCLE_COLLECTIONS = (
    "terrain.ocean_fabric",
    "terrain.margin_landforms",
    "terrain.arc_plume_landforms",
    "terrain.rift_margin_sequences",
)

P173_RESPONSE_KEYS = (
    "terrain.last_p173_age_aware_ocean_response_object_count",
    "terrain.last_p173_age_aware_ocean_response_area_fraction",
    "terrain.last_p173_age_aware_ocean_response_mean_abs_delta_m",
    "terrain.last_p173_age_aware_ocean_response_max_abs_delta_m",
)

P1732_DEPTH_FLOOR_KEYS = (
    "terrain.last_p1732_young_open_ocean_age_depth_used",
    "terrain.last_p1732_young_open_ocean_age_depth_land_mask_preserved",
    "terrain.last_p1732_young_open_ocean_age_depth_candidate_fraction",
    "terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction",
    "terrain.last_p1732_young_open_ocean_age_depth_mean_depth_before_m",
    "terrain.last_p1732_young_open_ocean_age_depth_mean_depth_after_m",
)

P1115_GUARDRAIL_KEYS = (
    "terrain.last_p1115_endpoint_readability_guardrail_mode",
    "terrain.last_p1115_process_relief_enabled",
    "terrain.last_p1115_process_relief_adjusted_area_fraction",
    "terrain.last_p1115_ocean_floor_adjusted_area_fraction",
    "terrain.last_p1115_ocean_floor_land_mask_preserved",
)

OCEAN_EXTREME_KEYS = (
    "ocean_fabric_entropy",
    "ridge_visible_fraction",
    "ridge_age_symmetry_score",
    "fracture_zone_length_fraction",
    "abyssal_plain_fraction",
    "hotspot_track_count",
    "seamount_chain_count",
    "oceanic_plateau_fraction",
    "microcontinent_fraction",
    "unparented_shoal_fraction",
)


def run_p173_ocean_lifecycle_gate(
    outdir: str | Path,
    *,
    cells: int = 8000,
    frames: int = 36,
    max_workers: int = 5,
    jobs: tuple[tuple[str, str, int], ...] = DEFAULT_JOBS,
    active_modules: frozenset[str] = DEFAULT_ACTIVE_MODULES,
) -> dict[str, Any]:
    """Run P173 ocean lifecycle checks and write per-world audits."""
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
        "min_ocean_area_fraction_for_gate": 0.10,
        "max_unparented_shoal_fraction": 0.05,
        "max_p173_response_area_fraction": 0.92,
        "min_p173_response_area_fraction_for_mean_abs_gate": 0.05,
        "max_p173_response_mean_abs_delta_m": 550.0,
    }
    (root / "p173_ocean_lifecycle_gate_config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    rows: list[dict[str, Any]] = []
    started = time.time()
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
                config,
            )
            for preset, label, seed in jobs
        ]
        for future in as_completed(futures):
            rows.append(future.result())
            rows.sort(key=lambda row: str(row["label"]))
            (root / "p173_ocean_lifecycle_partial_summary.json").write_text(
                json.dumps(_aggregate_rows(rows, config, started), indent=2,
                           sort_keys=True) + "\n",
                encoding="utf-8",
            )

    summary = _aggregate_rows(rows, config, started)
    (root / "p173_ocean_lifecycle_gate_summary.json").write_text(
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
    config: dict[str, Any],
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
    p171 = write_historical_object_audit(eng.world, eng.archive, outdir)
    p1731 = write_p173_unsupported_shoal_attribution(
        eng.world, eng.archive, outdir)
    compact = _compact_world_summary(
        label,
        preset_name,
        seed,
        cells,
        frames,
        p170,
        p171,
        p1731,
        world_globals=dict(eng.world.globals),
        config=config,
    )
    compact["runtime_seconds"] = float(time.time() - started)
    compact["active_modules"] = list(active_modules)
    compact["p170_audit_path"] = str(outdir / "p170_historical_geomorphology_audit.json")
    compact["p171_audit_path"] = str(outdir / "p171_historical_object_persistence_audit.json")
    compact["p1731_audit_path"] = str(outdir / "p173_unsupported_shoal_attribution.json")
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
    p171: dict[str, Any],
    p1731: dict[str, Any] | None = None,
    *,
    world_globals: dict[str, Any] | None = None,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    config = config or {}
    world_globals = world_globals or {}
    frame_rows = list(p170.get("frame_rows", []) or [])
    ocean_extremes = dict(
        p170.get("metric_extremes", {}).get("ocean_metrics", {}) or {})
    flag_summary = dict(p170.get("flag_summary", {}) or {})
    response = _pick_metrics(world_globals, P173_RESPONSE_KEYS)
    p1732_depth_floor = _pick_metrics(world_globals, P1732_DEPTH_FLOOR_KEYS)
    p1115_guardrail = _pick_metrics(world_globals, P1115_GUARDRAIL_KEYS)
    response_object_count = float(response[
        "terrain.last_p173_age_aware_ocean_response_object_count"])
    response_area = float(response[
        "terrain.last_p173_age_aware_ocean_response_area_fraction"])
    response_mean_abs = float(response[
        "terrain.last_p173_age_aware_ocean_response_mean_abs_delta_m"])
    ocean_area_max = _frame_metric_max(
        frame_rows, "ocean_area_fraction_of_world", group="ocean_metrics")
    ocean_gate_domain = (
        ocean_area_max >= float(config.get("min_ocean_area_fraction_for_gate", 0.10))
    )
    unparented_max = _extreme_value(
        ocean_extremes, "unparented_shoal_fraction", "max")
    ocean_object_counts = {
        key: int(p171.get("collection_object_counts", {}).get(key, 0))
        for key in OCEAN_LIFECYCLE_COLLECTIONS
    }
    ocean_object_total = int(sum(ocean_object_counts.values()))
    missing_response = bool(ocean_gate_domain and response_object_count <= 0.0)
    missing_objects = bool(ocean_gate_domain and ocean_object_total <= 0)
    overbroad_area = response_area > float(
        config.get("max_p173_response_area_fraction", 0.92))
    broad_high_mean_abs = (
        response_area >= float(
            config.get("min_p173_response_area_fraction_for_mean_abs_gate", 0.05))
        and response_mean_abs > float(
            config.get("max_p173_response_mean_abs_delta_m", 550.0))
    )
    overbroad_response = bool(overbroad_area or broad_high_mean_abs)
    unsupported_shoal = bool(
        unparented_max > float(config.get("max_unparented_shoal_fraction", 0.05))
    )
    p1731_compact = _compact_p1731_summary(p1731 or {})

    return {
        "schema": "aevum.p173_ocean_lifecycle_world_row.v1",
        "label": str(label),
        "preset": str(preset_name),
        "seed": int(seed),
        "cells": int(cells),
        "requested_frames": int(frames),
        "p170": {
            "frame_count": int(p170.get("frame_count", 0)),
            "usable_frame_count": int(p170.get("usable_frame_count", 0)),
            "required_metric_keys_present": bool(
                p170.get("required_metric_keys_present", False)),
            "ordinary_deep_ocean_frame_count": int(
                flag_summary.get("ordinary_deep_ocean_frame_count", 0)),
            "ocean_metric_extremes": _pick_ocean_extremes(ocean_extremes),
            "target_frame_rows": [
                _nearest_frame_row(frame_rows, target)
                for target in TARGET_TIMES_MYR
            ],
        },
        "p171": {
            "frame_count": int(p171.get("frame_count", 0)),
            "total_object_observations": int(
                p171.get("total_object_observations", 0)),
            "unique_object_id_count": int(p171.get("unique_object_id_count", 0)),
            "recurring_object_id_count": int(
                p171.get("recurring_object_id_count", 0)),
            "required_fields_complete": bool(
                p171.get("required_fields_complete", False)),
            "missing_required_field_slot_count": int(
                p171.get("missing_required_field_slot_count", 0)),
            "ocean_lifecycle_collection_counts": ocean_object_counts,
        },
        "p173": {
            "ocean_gate_domain": bool(ocean_gate_domain),
            "ocean_area_fraction_max": float(ocean_area_max),
            "unparented_shoal_fraction_max": float(unparented_max),
            "missing_ocean_lifecycle_objects": bool(missing_objects),
            "missing_ocean_response": bool(missing_response),
            "overbroad_ocean_response": bool(overbroad_response),
            "overbroad_ocean_response_area_gate": bool(overbroad_area),
            "overbroad_ocean_response_mean_abs_gate": bool(broad_high_mean_abs),
            "unsupported_shoal_candidate": bool(unsupported_shoal),
            "age_aware_ocean_response": response,
        },
        "p173_1": p1731_compact,
        "p173_2": {
            "young_open_ocean_depth_floor": p1732_depth_floor,
        },
        "p1115": {
            "terminal_guardrails": p1115_guardrail,
        },
        "acceptance": {
            "audit_completed": bool(frame_rows),
            "required_metric_keys_present": bool(
                p170.get("required_metric_keys_present", False)),
            "object_fields_complete": bool(
                p171.get("required_fields_complete", False)),
            "object_persistence_checked": bool(
                p171.get("acceptance", {}).get("persistence_checked", False)),
            "generation_behavior_changed_by_gate": False,
        },
    }


def _aggregate_rows(
    rows: list[dict[str, Any]],
    config: dict[str, Any],
    started: float,
) -> dict[str, Any]:
    rows = sorted(rows, key=lambda row: str(row["label"]))
    object_incomplete = [
        row for row in rows
        if not bool(row.get("p171", {}).get("required_fields_complete", False))
    ]
    metric_incomplete = [
        row for row in rows
        if not bool(row.get("p170", {}).get("required_metric_keys_present", False))
    ]
    missing_objects = [
        row for row in rows
        if bool(row.get("p173", {}).get("missing_ocean_lifecycle_objects", False))
    ]
    missing_response = [
        row for row in rows
        if bool(row.get("p173", {}).get("missing_ocean_response", False))
    ]
    overbroad_response = [
        row for row in rows
        if bool(row.get("p173", {}).get("overbroad_ocean_response", False))
    ]
    unsupported_shoal = [
        row for row in rows
        if bool(row.get("p173", {}).get("unsupported_shoal_candidate", False))
    ]
    p1731_missing = [
        row for row in rows
        if not bool(row.get("p173_1", {}).get("attribution_completed", False))
    ]
    response_rows = [
        dict(row.get("p173", {}).get("age_aware_ocean_response", {}) or {})
        for row in rows
    ]
    p1115_rows = [
        dict(row.get("p1115", {}).get("terminal_guardrails", {}) or {})
        for row in rows
    ]
    response_object_counts = [
        float(row.get(
            "terrain.last_p173_age_aware_ocean_response_object_count", 0.0))
        for row in response_rows
    ]
    response_areas = [
        float(row.get(
            "terrain.last_p173_age_aware_ocean_response_area_fraction", 0.0))
        for row in response_rows
    ]
    response_mean_abs = [
        float(row.get(
            "terrain.last_p173_age_aware_ocean_response_mean_abs_delta_m", 0.0))
        for row in response_rows
    ]
    p1732_candidate = [
        float(row.get("p173_1", {}).get(
            "max_p1732_young_open_ocean_age_depth_candidate_fraction", 0.0))
        for row in rows
    ]
    p1732_adjusted = [
        float(row.get("p173_1", {}).get(
            "max_p1732_young_open_ocean_age_depth_adjusted_fraction", 0.0))
        for row in rows
    ]
    p1115_endpoint_mode = [
        float(row.get(
            "terrain.last_p1115_endpoint_readability_guardrail_mode", 0.0))
        for row in p1115_rows
    ]
    p1115_process_enabled = [
        float(row.get("terrain.last_p1115_process_relief_enabled", 0.0))
        for row in p1115_rows
    ]
    p1115_process_adjusted = [
        float(row.get(
            "terrain.last_p1115_process_relief_adjusted_area_fraction", 0.0))
        for row in p1115_rows
    ]
    p1115_adjusted = [
        float(row.get("terrain.last_p1115_ocean_floor_adjusted_area_fraction", 0.0))
        for row in p1115_rows
    ]
    unparented = [
        float(row.get("p173", {}).get("unparented_shoal_fraction_max", 0.0))
        for row in rows
    ]
    fracture = [
        _world_ocean_extreme(row, "fracture_zone_length_fraction", "max")
        for row in rows
    ]
    p1731_residual = [
        float(row.get("p173_1", {}).get(
            "max_post_cleanup_residual_fraction_of_ocean", 0.0))
        for row in rows
    ]
    p1731_candidate = [
        float(row.get("p173_1", {}).get(
            "max_cleanup_candidate_fraction_of_ocean", 0.0))
        for row in rows
    ]
    return {
        "schema": SCHEMA,
        "config": config,
        "runtime_seconds": float(time.time() - started),
        "world_count": int(len(rows)),
        "completed_labels": [str(row["label"]) for row in rows],
        "worlds": rows,
        "aggregate": {
            "object_field_incomplete_world_count": int(len(object_incomplete)),
            "metric_incomplete_world_count": int(len(metric_incomplete)),
            "missing_ocean_lifecycle_object_world_count": int(len(missing_objects)),
            "missing_ocean_lifecycle_object_labels": [
                str(row["label"]) for row in missing_objects
            ],
            "missing_ocean_response_world_count": int(len(missing_response)),
            "missing_ocean_response_labels": [
                str(row["label"]) for row in missing_response
            ],
            "overbroad_ocean_response_world_count": int(len(overbroad_response)),
            "overbroad_ocean_response_labels": [
                str(row["label"]) for row in overbroad_response
            ],
            "unsupported_shoal_world_count": int(len(unsupported_shoal)),
            "unsupported_shoal_labels": [
                str(row["label"]) for row in unsupported_shoal
            ],
            "p1731_attribution_missing_world_count": int(len(p1731_missing)),
            "p1731_attribution_missing_labels": [
                str(row["label"]) for row in p1731_missing
            ],
            "ocean_response_world_count": int(sum(
                1 for value in response_object_counts if value > 0.0)),
            "max_ocean_response_area_fraction": (
                float(max(response_areas)) if response_areas else 0.0),
            "max_ocean_response_mean_abs_delta_m": (
                float(max(response_mean_abs)) if response_mean_abs else 0.0),
            "max_unparented_shoal_fraction": (
                float(max(unparented)) if unparented else 0.0),
            "max_fracture_zone_length_fraction": (
                float(max(fracture)) if fracture else 0.0),
            "max_p1731_cleanup_candidate_fraction": (
                float(max(p1731_candidate)) if p1731_candidate else 0.0),
            "max_p1731_post_cleanup_residual_fraction": (
                float(max(p1731_residual)) if p1731_residual else 0.0),
            "max_p1732_young_open_ocean_age_depth_candidate_fraction": (
                float(max(p1732_candidate)) if p1732_candidate else 0.0),
            "max_p1732_young_open_ocean_age_depth_adjusted_fraction": (
                float(max(p1732_adjusted)) if p1732_adjusted else 0.0),
            "p1115_endpoint_guardrail_mode_world_count": int(
                sum(1 for value in p1115_endpoint_mode if value > 0.0)),
            "p1115_process_relief_enabled_world_count": int(
                sum(1 for value in p1115_process_enabled if value > 0.0)),
            "max_p1115_process_relief_adjusted_area_fraction": (
                float(max(p1115_process_adjusted)) if p1115_process_adjusted else 0.0),
            "max_p1115_ocean_floor_adjusted_area_fraction": (
                float(max(p1115_adjusted)) if p1115_adjusted else 0.0),
        },
        "acceptance": {
            "gate_completed": bool(len(rows) == len(config.get("jobs", []))),
            "required_metric_keys_present": bool(not metric_incomplete),
            "object_fields_complete": bool(not object_incomplete),
            "p173_diagnostic_completed": bool(len(rows) > 0),
            "ocean_lifecycle_objects_present": bool(not missing_objects),
            "ocean_response_present": bool(not missing_response),
            "ocean_response_not_overbroad": bool(not overbroad_response),
            "unsupported_shoal_gate_passed": bool(not unsupported_shoal),
            "p1731_attribution_available": bool(not p1731_missing),
            "p1115_terminal_process_relief_contracted": bool(
                not any(value > 0.0 for value in p1115_process_enabled)
            ),
            "generation_behavior_changed_by_gate": False,
            "six_world_8000_gate_completed": bool(
                len(rows) == 6 and int(config.get("cells", 0)) == 8000),
        },
    }


def _nearest_frame_row(rows: list[dict[str, Any]], target_myr: float) -> dict[str, Any]:
    if not rows:
        return {
            "target_myr": float(target_myr),
            "available": False,
            "time_myr": 0.0,
            "ocean_metrics": {},
        }
    row = min(rows, key=lambda item: abs(float(item.get("time_myr", 0.0)) - float(target_myr)))
    return {
        "target_myr": float(target_myr),
        "available": True,
        "time_myr": float(row.get("time_myr", 0.0)),
        "ocean_metrics": _pick_metrics(
            row.get("ocean_metrics", {}),
            OCEAN_EXTREME_KEYS,
        ),
    }


def _compact_p1731_summary(summary: dict[str, Any]) -> dict[str, Any]:
    extremes = dict(summary.get("metric_extremes", {}) or {})
    p1732_extremes = dict(summary.get("p1732_metric_extremes", {}) or {})
    peak = dict(summary.get("peak_residual_frame", {}) or {})
    candidate = dict(extremes.get("cleanup_candidate_fraction_of_ocean", {}) or {})
    residual = dict(extremes.get("post_cleanup_residual_fraction_of_ocean", {}) or {})
    preserve = dict(extremes.get("structural_preserve_fraction_of_candidate", {}) or {})
    object_support = dict(extremes.get("object_support_fraction_of_candidate", {}) or {})
    semantic_support = dict(extremes.get("semantic_support_fraction_of_candidate", {}) or {})
    p1732_candidate = dict(p1732_extremes.get(
        "terrain.last_p1732_young_open_ocean_age_depth_candidate_fraction", {}) or {})
    p1732_adjusted = dict(p1732_extremes.get(
        "terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction", {}) or {})
    return {
        "frame_count": int(summary.get("frame_count", 0)),
        "usable_frame_count": int(summary.get("usable_frame_count", 0)),
        "attribution_completed": bool(
            summary.get("acceptance", {}).get("attribution_completed", False)),
        "max_cleanup_candidate_fraction_of_ocean": float(candidate.get("max", 0.0)),
        "max_post_cleanup_residual_fraction_of_ocean": float(residual.get("max", 0.0)),
        "max_structural_preserve_fraction_of_candidate": float(
            preserve.get("max", 0.0)),
        "max_object_support_fraction_of_candidate": float(
            object_support.get("max", 0.0)),
        "max_semantic_support_fraction_of_candidate": float(
            semantic_support.get("max", 0.0)),
        "max_p1732_young_open_ocean_age_depth_candidate_fraction": float(
            p1732_candidate.get("max", 0.0)),
        "max_p1732_young_open_ocean_age_depth_adjusted_fraction": float(
            p1732_adjusted.get("max", 0.0)),
        "peak_residual_frame": peak,
    }


def _pick_ocean_extremes(ocean_extremes: dict[str, Any]) -> dict[str, Any]:
    return {
        key: dict(ocean_extremes.get(key, {}) or {})
        for key in OCEAN_EXTREME_KEYS
    }


def _world_ocean_extreme(row: dict[str, Any], metric: str, stat: str) -> float:
    extremes = dict(row.get("p170", {}).get("ocean_metric_extremes", {}) or {})
    return _extreme_value(extremes, metric, stat)


def _extreme_value(
    extremes: dict[str, Any],
    metric: str,
    stat: str,
) -> float:
    value = dict(extremes.get(metric, {}) or {}).get(stat, 0.0)
    return float(value)


def _frame_metric_max(
    rows: list[dict[str, Any]],
    metric: str,
    *,
    group: str,
) -> float:
    values = [
        float(row.get(group, {}).get(metric, 0.0))
        for row in rows
    ]
    return float(max(values)) if values else 0.0


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


def _set_process_thread_env() -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="python -m aevum.diagnostics.p173_ocean_lifecycle_gate",
        description="Run the P173 ocean-floor lifecycle diagnostic gate.",
    )
    parser.add_argument("--out", default="out_p173_ocean_lifecycle_gate")
    parser.add_argument("--cells", type=int, default=8000)
    parser.add_argument("--frames", type=int, default=36)
    parser.add_argument("--max-workers", type=int, default=5)
    parser.add_argument(
        "--job",
        action="append",
        default=[],
        help="Optional preset:label:seed job. May be repeated.",
    )
    args = parser.parse_args(argv)
    jobs = _parse_jobs(tuple(args.job)) if args.job else DEFAULT_JOBS
    summary = run_p173_ocean_lifecycle_gate(
        args.out,
        cells=int(args.cells),
        frames=int(args.frames),
        max_workers=int(args.max_workers),
        jobs=jobs,
    )
    print("== aevum :: P173 ocean lifecycle gate ==")
    print(f"   worlds: {summary['world_count']}")
    print(
        "   unsupported shoal worlds: "
        f"{summary['aggregate']['unsupported_shoal_world_count']}"
    )
    print(
        "   overbroad ocean response worlds: "
        f"{summary['aggregate']['overbroad_ocean_response_world_count']}"
    )
    print(
        "   max response area: "
        f"{summary['aggregate']['max_ocean_response_area_fraction']:.3f}"
    )
    print(f"   wrote {Path(args.out) / 'p173_ocean_lifecycle_gate_summary.json'}")


def _parse_jobs(values: tuple[str, ...]) -> tuple[tuple[str, str, int], ...]:
    jobs: list[tuple[str, str, int]] = []
    for value in values:
        parts = str(value).split(":")
        if len(parts) != 3:
            raise SystemExit(
                f"invalid --job {value!r}; expected preset:label:seed")
        preset, label, seed_text = parts
        try:
            seed = int(seed_text)
        except ValueError as exc:
            raise SystemExit(
                f"invalid --job {value!r}; seed must be an integer") from exc
        jobs.append((preset, label, seed))
    return tuple(jobs)


if __name__ == "__main__":
    main()
