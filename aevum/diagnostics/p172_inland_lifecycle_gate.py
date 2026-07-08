"""P172 multi-world inland landform lifecycle gate."""
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
from aevum.diagnostics.historical_objects import write_historical_object_audit
from aevum.engine import Engine
from aevum.spec.presets import get_preset


SCHEMA = "aevum.p172_inland_lifecycle_gate.v1"

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


def run_p172_inland_lifecycle_gate(
    outdir: str | Path,
    *,
    cells: int = 8000,
    frames: int = 36,
    max_workers: int = 5,
    jobs: tuple[tuple[str, str, int], ...] = DEFAULT_JOBS,
    active_modules: frozenset[str] = DEFAULT_ACTIVE_MODULES,
) -> dict[str, Any]:
    """Run P172 inland lifecycle checks and write per-world audits."""
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
        "residual_rift_mature_max_threshold": 0.65,
        "residual_rift_any_max_threshold": 0.75,
        "min_inland_area_fraction_for_rift_gate": 0.02,
    }
    (root / "p172_inland_lifecycle_gate_config.json").write_text(
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
            (root / "p172_inland_lifecycle_partial_summary.json").write_text(
                json.dumps(_aggregate_rows(rows, config, started), indent=2,
                           sort_keys=True) + "\n",
                encoding="utf-8",
            )

    summary = _aggregate_rows(rows, config, started)
    (root / "p172_inland_lifecycle_gate_summary.json").write_text(
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
    compact = _compact_world_summary(
        label,
        preset_name,
        seed,
        cells,
        frames,
        p170,
        p171,
        world_globals=dict(eng.world.globals),
        config=config,
    )
    compact["runtime_seconds"] = float(time.time() - started)
    compact["active_modules"] = list(active_modules)
    compact["p170_audit_path"] = str(outdir / "p170_historical_geomorphology_audit.json")
    compact["p171_audit_path"] = str(outdir / "p171_historical_object_persistence_audit.json")
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
    *,
    world_globals: dict[str, Any] | None = None,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    config = config or {}
    world_globals = world_globals or {}
    frame_rows = list(p170.get("frame_rows", []) or [])
    target_rows = [_nearest_frame_row(frame_rows, target) for target in TARGET_TIMES_MYR]
    land_extremes = dict(p170.get("metric_extremes", {}).get("land_metrics", {}) or {})
    flag_summary = dict(p170.get("flag_summary", {}) or {})
    p174_continuity = dict(p170.get("p174_continuity", {}) or {})
    rift_max = _extreme_value(land_extremes, "rift_basin_expression_fraction", "max")
    rift_median = _extreme_value(land_extremes, "rift_basin_expression_fraction", "median")
    gated_rift_max = _filtered_metric_max(
        frame_rows,
        "rift_basin_expression_fraction",
        min_time_myr=0.0,
        min_inland_area_fraction=float(
            config.get("min_inland_area_fraction_for_rift_gate", 0.02)),
    )
    mature_rift_max = _mature_metric_max(
        frame_rows,
        "rift_basin_expression_fraction",
        min_time_myr=2500.0,
        min_inland_area_fraction=float(
            config.get("min_inland_area_fraction_for_rift_gate", 0.02)),
    )
    residual_rift = bool(
        gated_rift_max >= float(config.get("residual_rift_any_max_threshold", 0.75))
        or mature_rift_max >= float(config.get("residual_rift_mature_max_threshold", 0.65))
    )

    return {
        "schema": "aevum.p172_inland_lifecycle_world_row.v1",
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
            "ordinary_plateau_frame_count": int(
                flag_summary.get("ordinary_plateau_frame_count", 0)),
            "ordinary_deep_ocean_frame_count": int(
                flag_summary.get("ordinary_deep_ocean_frame_count", 0)),
            "lowland_plain_deficient_frame_count": int(
                flag_summary.get("lowland_plain_deficient_frame_count", 0)),
            "p174_continuity": _pick_metrics(
                p174_continuity,
                (
                    "mature_support_frame_count",
                    "mature_lowland_deficient_frame_count",
                    "mature_lowland_deficient_fraction",
                    "mature_lowland_continuity_score",
                    "mature_lowland_residual_area_limited_frame_count",
                    "mature_lowland_residual_parentage_limited_frame_count",
                    "mature_lowland_residual_relief_boundary_frame_count",
                    "mature_lowland_residual_support_to_plain_frame_count",
                    "mature_lowland_residual_upstream_support_frame_count",
                    "mature_lowland_residual_active_exclusion_frame_count",
                    "mature_lowland_residual_component_segmentation_frame_count",
                    "lowland_plain_max_positive_step",
                    "lowland_plain_max_negative_step",
                    "terminal_lowland_plain_jump",
                    "terminal_lowland_pop_candidate",
                    "ocean_fabric_entropy_max_positive_step",
                    "ocean_fabric_entropy_max_negative_step",
                    "terminal_ocean_fabric_entropy_jump",
                ),
            ),
            "land_metric_extremes": _pick_land_extremes(land_extremes),
            "target_frame_rows": target_rows,
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
            "land_lifecycle_collection_counts": {
                key: int(p171.get("collection_object_counts", {}).get(key, 0))
                for key in (
                    "terrain.continental_landforms",
                    "tectonics.continental_provinces",
                    "terrain.mountain_ranges",
                    "terrain.plateau_inventory",
                )
            },
        },
        "p172": {
            "rift_basin_expression_max": float(rift_max),
            "gated_rift_basin_expression_max": float(gated_rift_max),
            "rift_basin_expression_median": float(rift_median),
            "mature_rift_basin_expression_max": float(mature_rift_max),
            "residual_rift_overpaint_candidate": bool(residual_rift),
            "land_lifecycle_objects_present": bool(
                int(p171.get("collection_object_counts", {}).get(
                    "terrain.continental_landforms", 0)) > 0
                or str(preset_name) == "waterworld"
            ),
            "age_aware_inland_response": _pick_metrics(
                world_globals,
                (
                    "terrain.last_p172_age_aware_inland_response_object_count",
                    "terrain.last_p172_age_aware_inland_response_area_fraction",
                    "terrain.last_p172_age_aware_inland_response_mean_abs_delta_m",
                    "terrain.last_p172_age_aware_inland_response_max_abs_delta_m",
                ),
            ),
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
    residual = [
        row for row in rows
        if bool(row.get("p172", {}).get("residual_rift_overpaint_candidate", False))
    ]
    object_incomplete = [
        row for row in rows
        if not bool(row.get("p171", {}).get("required_fields_complete", False))
    ]
    metric_incomplete = [
        row for row in rows
        if not bool(row.get("p170", {}).get("required_metric_keys_present", False))
    ]
    rift_max_values = [
        float(row.get("p172", {}).get("rift_basin_expression_max", 0.0))
        for row in rows
    ]
    mature_rift_max_values = [
        float(row.get("p172", {}).get("mature_rift_basin_expression_max", 0.0))
        for row in rows
    ]
    age_aware_rows = [
        dict(row.get("p172", {}).get("age_aware_inland_response", {}) or {})
        for row in rows
    ]
    age_aware_object_counts = [
        float(row.get(
            "terrain.last_p172_age_aware_inland_response_object_count", 0.0))
        for row in age_aware_rows
    ]
    age_aware_area_fractions = [
        float(row.get(
            "terrain.last_p172_age_aware_inland_response_area_fraction", 0.0))
        for row in age_aware_rows
    ]
    age_aware_mean_abs = [
        float(row.get(
            "terrain.last_p172_age_aware_inland_response_mean_abs_delta_m", 0.0))
        for row in age_aware_rows
    ]
    p174_continuity_rows = [
        dict(row.get("p170", {}).get("p174_continuity", {}) or {})
        for row in rows
    ]
    p174_lowland_risk = [
        row for row in rows
        if bool(row.get("p170", {}).get("p174_continuity", {}).get(
            "terminal_lowland_pop_candidate", False))
        or float(row.get("p170", {}).get("p174_continuity", {}).get(
            "mature_lowland_continuity_score", 1.0)) < 0.70
    ]
    p174_scores = [
        float(row.get("mature_lowland_continuity_score", 1.0))
        for row in p174_continuity_rows
    ]
    p174_deficient_counts = [
        float(row.get("mature_lowland_deficient_frame_count", 0.0))
        for row in p174_continuity_rows
    ]
    return {
        "schema": SCHEMA,
        "config": config,
        "runtime_seconds": float(time.time() - started),
        "world_count": int(len(rows)),
        "completed_labels": [str(row["label"]) for row in rows],
        "worlds": rows,
        "aggregate": {
            "residual_rift_overpaint_world_count": int(len(residual)),
            "residual_rift_overpaint_labels": [str(row["label"]) for row in residual],
            "object_field_incomplete_world_count": int(len(object_incomplete)),
            "metric_incomplete_world_count": int(len(metric_incomplete)),
            "max_rift_basin_expression": (
                float(max(rift_max_values)) if rift_max_values else 0.0),
            "median_rift_basin_expression_max": (
                float(np.median(np.asarray(rift_max_values, dtype=np.float64)))
                if rift_max_values else 0.0),
            "max_mature_rift_basin_expression": (
                float(max(mature_rift_max_values)) if mature_rift_max_values else 0.0),
            "age_aware_inland_response_world_count": int(sum(
                1 for value in age_aware_object_counts if value > 0.0)),
            "max_age_aware_inland_response_area_fraction": (
                float(max(age_aware_area_fractions))
                if age_aware_area_fractions else 0.0),
            "max_age_aware_inland_response_mean_abs_delta_m": (
                float(max(age_aware_mean_abs)) if age_aware_mean_abs else 0.0),
            "p174_lowland_continuity_risk_world_count": int(len(p174_lowland_risk)),
            "p174_lowland_continuity_risk_labels": [
                str(row["label"]) for row in p174_lowland_risk
            ],
            "min_p174_lowland_continuity_score": (
                float(min(p174_scores)) if p174_scores else 1.0),
            "max_p174_mature_lowland_deficient_frame_count": (
                float(max(p174_deficient_counts)) if p174_deficient_counts else 0.0),
        },
        "acceptance": {
            "gate_completed": bool(len(rows) == len(config.get("jobs", []))),
            "required_metric_keys_present": bool(not metric_incomplete),
            "object_fields_complete": bool(not object_incomplete),
            "p172_diagnostic_completed": bool(len(rows) > 0),
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
            "land_metrics": {},
        }
    row = min(rows, key=lambda item: abs(float(item.get("time_myr", 0.0)) - float(target_myr)))
    return {
        "target_myr": float(target_myr),
        "available": True,
        "time_myr": float(row.get("time_myr", 0.0)),
        "land_metrics": _pick_metrics(
            row.get("land_metrics", {}),
            (
                "inland_detail_entropy",
                "ordinary_plateau_fraction",
                "broad_flat_inland_component_count",
                "continent_province_count_p50",
                "old_orogen_expression_fraction",
                "rift_basin_expression_fraction",
                "craton_shield_platform_split_fraction",
                "lowland_plain_fraction",
                "lowland_elevation_parented_fraction",
                "largest_lowland_elevation_parented_component_fraction",
                "lowland_local_relief_blocked_fraction",
                "broad_lowland_plain_component_count",
                "largest_lowland_plain_component_fraction",
                "lowland_plain_parented_fraction",
                "lowland_plain_area_gap_fraction",
                "lowland_plain_component_gap_fraction",
                "lowland_plain_parentage_gap_fraction",
                "lowland_relief_boundary_gap_fraction",
                "lowland_support_to_plain_gap_fraction",
                "lowland_upstream_support_gap_fraction",
                "lowland_active_exclusion_fraction",
                "lowland_near_floor_plain_component_count",
                "lowland_near_floor_parented_component_count",
                "lowland_residual_area_limited",
                "lowland_residual_parentage_limited",
                "lowland_residual_relief_boundary_limited",
                "lowland_residual_support_to_plain_limited",
                "lowland_residual_upstream_support_limited",
                "lowland_residual_active_exclusion_limited",
                "lowland_residual_component_segmentation_limited",
                "lowland_residual_dominant_code",
                "p104f_pre_p174_lowland_prep_area_fraction",
                "p104f_pre_p174_lowland_prep_mean_lowering_m",
                "p104f_pre_p174_lowland_support_area_fraction",
                "p104f_pre_p174_lowland_support_largest_component_fraction",
                "p104f_pre_p174_lowland_memory_seed_area_fraction",
                "p174_lowland_plain_response_area_fraction",
                "p174_lowland_plain_candidate_area_fraction",
                "p174_lowland_plain_parent_area_fraction",
                "p174_lowland_plain_continuity_memory_area_fraction",
                "p174_lowland_plain_continuity_parent_area_fraction",
                "p174_lowland_plain_response_mean_abs_delta_m",
                "p174_lowland_plain_fraction_before",
                "p174_lowland_plain_fraction_after",
                "p174_lowland_plain_largest_component_fraction_after",
                "p174_lowland_plain_parented_fraction_after",
                "p174_lowland_plain_response_stage_code",
                "p174_support_component_response_area_fraction",
                "p174_support_component_response_largest_component_fraction",
            ),
        ),
    }


def _pick_land_extremes(land_extremes: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "inland_detail_entropy",
        "ordinary_plateau_fraction",
        "broad_flat_inland_component_count",
        "continent_province_count_p50",
        "old_orogen_expression_fraction",
        "rift_basin_expression_fraction",
        "craton_shield_platform_split_fraction",
        "lowland_plain_fraction",
        "lowland_elevation_parented_fraction",
        "largest_lowland_elevation_parented_component_fraction",
        "lowland_local_relief_blocked_fraction",
        "broad_lowland_plain_component_count",
        "largest_lowland_plain_component_fraction",
        "lowland_plain_parented_fraction",
        "lowland_plain_area_gap_fraction",
        "lowland_plain_component_gap_fraction",
        "lowland_plain_parentage_gap_fraction",
        "lowland_relief_boundary_gap_fraction",
        "lowland_support_to_plain_gap_fraction",
        "lowland_upstream_support_gap_fraction",
        "lowland_active_exclusion_fraction",
        "lowland_near_floor_plain_component_count",
        "lowland_near_floor_parented_component_count",
        "lowland_residual_area_limited",
        "lowland_residual_parentage_limited",
        "lowland_residual_relief_boundary_limited",
        "lowland_residual_support_to_plain_limited",
        "lowland_residual_upstream_support_limited",
        "lowland_residual_active_exclusion_limited",
        "lowland_residual_component_segmentation_limited",
        "lowland_residual_dominant_code",
        "p104f_pre_p174_lowland_prep_area_fraction",
        "p104f_pre_p174_lowland_prep_mean_lowering_m",
        "p104f_pre_p174_lowland_support_area_fraction",
        "p104f_pre_p174_lowland_support_largest_component_fraction",
        "p104f_pre_p174_lowland_memory_seed_area_fraction",
        "p174_lowland_plain_response_area_fraction",
        "p174_lowland_plain_candidate_area_fraction",
        "p174_lowland_plain_parent_area_fraction",
        "p174_lowland_plain_continuity_memory_area_fraction",
        "p174_lowland_plain_continuity_parent_area_fraction",
        "p174_lowland_plain_response_mean_abs_delta_m",
        "p174_lowland_plain_fraction_before",
        "p174_lowland_plain_fraction_after",
        "p174_lowland_plain_largest_component_fraction_after",
        "p174_lowland_plain_parented_fraction_after",
        "p174_lowland_plain_response_stage_code",
        "p174_support_component_response_area_fraction",
        "p174_support_component_response_largest_component_fraction",
    )
    return {
        key: dict(land_extremes.get(key, {}) or {})
        for key in keys
    }


def _extreme_value(
    extremes: dict[str, Any],
    metric: str,
    stat: str,
) -> float:
    value = dict(extremes.get(metric, {}) or {}).get(stat, 0.0)
    return float(value)


def _mature_metric_max(
    rows: list[dict[str, Any]],
    metric: str,
    *,
    min_time_myr: float,
    min_inland_area_fraction: float = 0.0,
) -> float:
    return _filtered_metric_max(
        rows,
        metric,
        min_time_myr=min_time_myr,
        min_inland_area_fraction=min_inland_area_fraction,
    )


def _filtered_metric_max(
    rows: list[dict[str, Any]],
    metric: str,
    *,
    min_time_myr: float,
    min_inland_area_fraction: float,
) -> float:
    values = [
        float(row.get("land_metrics", {}).get(metric, 0.0))
        for row in rows
        if float(row.get("time_myr", 0.0)) >= float(min_time_myr)
        and float(row.get("land_metrics", {}).get(
            "inland_area_fraction_of_world", 0.0)) >= float(min_inland_area_fraction)
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
        prog="python -m aevum.diagnostics.p172_inland_lifecycle_gate",
        description="Run the P172 inland landform lifecycle diagnostic gate.",
    )
    parser.add_argument("--out", default="out_p172_inland_lifecycle_gate")
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
    summary = run_p172_inland_lifecycle_gate(
        args.out,
        cells=int(args.cells),
        frames=int(args.frames),
        max_workers=int(args.max_workers),
        jobs=jobs,
    )
    print("== aevum :: P172 inland lifecycle gate ==")
    print(f"   worlds: {summary['world_count']}")
    print(
        "   residual rift worlds: "
        f"{summary['aggregate']['residual_rift_overpaint_world_count']}"
    )
    print(
        "   max mature rift expression: "
        f"{summary['aggregate']['max_mature_rift_basin_expression']:.3f}"
    )
    print(f"   wrote {Path(args.out) / 'p172_inland_lifecycle_gate_summary.json'}")


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
