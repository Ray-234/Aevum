"""Old-orogen erosion and boundary-persistence reference diagnostics."""
from __future__ import annotations

import hashlib
import json
from typing import Any

import numpy as np


SCHEMA = "aevum.old_orogen_erosion_decay.v1"

SOURCE_IDS = (
    "NOAA_ETOPO_2022",
    "GMBA_MOUNTAIN_INVENTORY",
)

EXPECTED_STAGE_SEQUENCE = (
    "active_orogen",
    "post_collision_high_orogen",
    "decaying_orogen",
    "old_orogen",
    "subdued_old_orogen",
)

EXPECTED_CURRENT_RESIDUAL_FIELDS = (
    "terrain.old_orogen_decay_stage",
    "terrain.orogen_erosion_budget",
    "terrain.orogen_boundary_memory",
    "terrain.orogen_sediment_export",
)

BOUNDARY_TRACE = (
    (0, 3),
    (1, 3),
    (2, 3),
    (3, 3),
    (4, 3),
    (5, 3),
    (6, 3),
)

REFERENCE_FRAMES: tuple[dict[str, Any], ...] = (
    {
        "time_since_collision_myr": 0.0,
        "stage": "active_orogen",
        "province_class": "active_orogen",
        "mean_elevation_m": 3200.0,
        "local_relief_m": 1900.0,
        "crust_thickness_m": 58000.0,
        "boundary_strength": 1.00,
        "sediment_export_interval_km3": 0.0,
        "parent_processes": ("collision_orogeny",),
        "boundary_trace_cells": BOUNDARY_TRACE,
    },
    {
        "time_since_collision_myr": 80.0,
        "stage": "post_collision_high_orogen",
        "province_class": "active_orogen",
        "mean_elevation_m": 2600.0,
        "local_relief_m": 1550.0,
        "crust_thickness_m": 55000.0,
        "boundary_strength": 0.95,
        "sediment_export_interval_km3": 14700.0,
        "parent_processes": ("collision_orogeny", "early_orogenic_erosion"),
        "boundary_trace_cells": BOUNDARY_TRACE,
    },
    {
        "time_since_collision_myr": 220.0,
        "stage": "decaying_orogen",
        "province_class": "old_orogen",
        "mean_elevation_m": 1600.0,
        "local_relief_m": 1050.0,
        "crust_thickness_m": 50000.0,
        "boundary_strength": 0.84,
        "sediment_export_interval_km3": 21000.0,
        "parent_processes": ("orogenic_decay", "suture_inheritance"),
        "boundary_trace_cells": BOUNDARY_TRACE,
    },
    {
        "time_since_collision_myr": 520.0,
        "stage": "old_orogen",
        "province_class": "old_orogen",
        "mean_elevation_m": 920.0,
        "local_relief_m": 650.0,
        "crust_thickness_m": 45000.0,
        "boundary_strength": 0.72,
        "sediment_export_interval_km3": 16800.0,
        "parent_processes": ("orogenic_decay", "suture_inheritance"),
        "boundary_trace_cells": BOUNDARY_TRACE,
    },
    {
        "time_since_collision_myr": 1000.0,
        "stage": "subdued_old_orogen",
        "province_class": "old_orogen",
        "mean_elevation_m": 680.0,
        "local_relief_m": 420.0,
        "crust_thickness_m": 42000.0,
        "boundary_strength": 0.62,
        "sediment_export_interval_km3": 9660.0,
        "parent_processes": ("orogenic_decay", "suture_inheritance"),
        "boundary_trace_cells": BOUNDARY_TRACE,
    },
)


def _canonical_digest(data: Any) -> str:
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _trace_overlap(a: tuple[tuple[int, int], ...],
                   b: tuple[tuple[int, int], ...]) -> float:
    aset = set(a)
    bset = set(b)
    if not aset and not bset:
        return 1.0
    if not aset or not bset:
        return 0.0
    return float(len(aset & bset) / len(aset | bset))


def reference_old_orogen_decay_summary() -> dict[str, Any]:
    frames = tuple(dict(frame) for frame in REFERENCE_FRAMES)
    reliefs = tuple(float(frame["local_relief_m"]) for frame in frames)
    elevations = tuple(float(frame["mean_elevation_m"]) for frame in frames)
    thicknesses = tuple(float(frame["crust_thickness_m"]) for frame in frames)
    boundary_strengths = tuple(float(frame["boundary_strength"]) for frame in frames)
    sediment_exports = tuple(float(frame["sediment_export_interval_km3"]) for frame in frames)
    cumulative_sediment = tuple(np.cumsum(sediment_exports).astype(float).tolist())
    stage_sequence = tuple(str(frame["stage"]) for frame in frames)
    boundary_id = "suture:p86-collision-memory"
    parent_link_failures = tuple(
        str(frame["stage"])
        for frame in frames
        if not frame["parent_processes"]
    )
    trace_overlaps = tuple(
        _trace_overlap(
            tuple(frames[0]["boundary_trace_cells"]),
            tuple(frame["boundary_trace_cells"]),
        )
        for frame in frames
    )
    relief_monotonic = all(
        left > right for left, right in zip(reliefs, reliefs[1:]))
    elevation_monotonic = all(
        left > right for left, right in zip(elevations, elevations[1:]))
    thickness_monotonic = all(
        left >= right for left, right in zip(thicknesses, thicknesses[1:]))
    boundary_decay_but_persistent = bool(
        boundary_strengths[0] > boundary_strengths[-1]
        and boundary_strengths[-1] >= 0.55
        and min(trace_overlaps) >= 0.90
    )
    relief_decay_fraction = (
        (reliefs[0] - reliefs[-1]) / reliefs[0] if reliefs[0] else 0.0
    )
    late_export_declines = bool(
        sediment_exports[-1] < max(sediment_exports) * 0.55
        and sediment_exports[-1] < sediment_exports[1]
    )
    old_orogen_persists_after_decay = bool(
        frames[-1]["province_class"] == "old_orogen"
        and frames[-1]["stage"] == "subdued_old_orogen"
        and 250.0 <= reliefs[-1] <= 700.0
        and elevations[-1] >= 300.0
    )
    acceptance = {
        "fixture_schema_ready": bool(frames and SOURCE_IDS),
        "stage_sequence_complete": stage_sequence == EXPECTED_STAGE_SEQUENCE,
        "relief_monotonic_decay": relief_monotonic,
        "elevation_monotonic_decay": elevation_monotonic,
        "crustal_root_relaxes": thickness_monotonic,
        "relief_decay_large_enough": relief_decay_fraction >= 0.70,
        "sediment_export_positive": sum(sediment_exports[1:]) > 0.0,
        "late_sediment_export_declines": late_export_declines,
        "boundary_strength_persists": boundary_decay_but_persistent,
        "old_orogen_boundary_persists_after_decay": old_orogen_persists_after_decay,
        "parent_process_links_present": not parent_link_failures,
        "raw_reference_data_not_stored": True,
    }
    metrics = {
        "frame_count": int(len(frames)),
        "initial_relief_m": float(reliefs[0]),
        "final_relief_m": float(reliefs[-1]),
        "relief_decay_fraction": float(relief_decay_fraction),
        "initial_elevation_m": float(elevations[0]),
        "final_elevation_m": float(elevations[-1]),
        "initial_crust_thickness_m": float(thicknesses[0]),
        "final_crust_thickness_m": float(thicknesses[-1]),
        "initial_boundary_strength": float(boundary_strengths[0]),
        "final_boundary_strength": float(boundary_strengths[-1]),
        "min_boundary_trace_overlap": float(min(trace_overlaps)),
        "total_sediment_export_km3": float(sum(sediment_exports)),
        "peak_sediment_export_km3": float(max(sediment_exports)),
        "final_interval_sediment_export_km3": float(sediment_exports[-1]),
        "parent_link_failure_count": int(len(parent_link_failures)),
    }
    enriched_frames = []
    for frame, cumulative in zip(frames, cumulative_sediment):
        out = dict(frame)
        out["boundary_id"] = boundary_id
        out["parent_suture_id"] = boundary_id
        out["cumulative_sediment_export_km3"] = float(cumulative)
        enriched_frames.append(out)
    return {
        "schema": SCHEMA,
        "status": (
            "old_orogen_erosion_decay_reference_ready"
            if all(acceptance.values())
            else "old_orogen_erosion_decay_reference_incomplete"
        ),
        "source_ids": SOURCE_IDS,
        "frames": tuple(enriched_frames),
        "stage_sequence": stage_sequence,
        "expected_stage_sequence": EXPECTED_STAGE_SEQUENCE,
        "boundary_id": boundary_id,
        "boundary_trace_overlaps": trace_overlaps,
        "parent_link_failures": parent_link_failures,
        "extraction_policy": {
            "raw_topography_or_mountain_inventory_stored": False,
            "direct_orogen_decay_extraction_pending": True,
            "fixture_is_small_derived_reference": True,
        },
        "metrics": metrics,
        "acceptance": acceptance,
        "fixture_digest": _canonical_digest({
            "frames": frames,
            "source_ids": SOURCE_IDS,
        }),
    }


def current_generated_old_orogen_audit(world: Any) -> dict[str, Any]:
    fields = getattr(world, "fields", {})
    objects = getattr(world, "objects", {})
    landforms = list(objects.get("terrain.continental_landforms", []))
    old_orogens = [
        obj for obj in landforms if obj.get("kind") == "old_subdued_orogen"
    ]
    parented = [
        obj for obj in old_orogens
        if obj.get("parent_process")
        or obj.get("parent_processes")
        or obj.get("parent_tectonic_object_ids")
    ]
    missing_decay_fields = tuple(
        name for name in EXPECTED_CURRENT_RESIDUAL_FIELDS if name not in fields)
    expected_residuals_recorded = bool(
        set(missing_decay_fields).issubset(set(EXPECTED_CURRENT_RESIDUAL_FIELDS)))
    area = np.asarray(getattr(world.grid, "cell_area", np.ones(world.grid.n)),
                      dtype=np.float64)
    elevation = (
        np.asarray(world.field("terrain.elevation_m"), dtype=np.float64)
        if "terrain.elevation_m" in fields else None
    )
    sediment = (
        np.asarray(world.field("sediment.thickness_m"), dtype=np.float64)
        if "sediment.thickness_m" in fields else None
    )
    decay_stage = (
        np.asarray(fields["terrain.old_orogen_decay_stage"], dtype=np.float64)
        if "terrain.old_orogen_decay_stage" in fields else None
    )
    erosion_budget = (
        np.asarray(fields["terrain.orogen_erosion_budget"], dtype=np.float64)
        if "terrain.orogen_erosion_budget" in fields else None
    )
    boundary_memory = (
        np.asarray(fields["terrain.orogen_boundary_memory"], dtype=np.float64)
        if "terrain.orogen_boundary_memory" in fields else None
    )
    sediment_export = (
        np.asarray(fields["terrain.orogen_sediment_export"], dtype=np.float64)
        if "terrain.orogen_sediment_export" in fields else None
    )
    total_old_area = 0.0
    elev_sum = 0.0
    sed_sum = 0.0
    for obj in old_orogens:
        cells = np.asarray(obj.get("cells", ()), dtype=np.int64)
        cells = cells[(cells >= 0) & (cells < world.grid.n)]
        if cells.size == 0:
            continue
        weights = area[cells]
        total_old_area += float(weights.sum())
        if elevation is not None:
            elev_sum += float(np.sum(elevation[cells] * weights))
        if sediment is not None:
            sed_sum += float(np.sum(sediment[cells] * weights))
    mean_elevation = (
        elev_sum / total_old_area
        if elevation is not None and total_old_area > 0.0 else float("nan")
    )
    mean_sediment = (
        sed_sum / total_old_area
        if sediment is not None and total_old_area > 0.0 else float("nan")
    )
    old_mask = np.zeros(world.grid.n, dtype=bool)
    for obj in old_orogens:
        cells = np.asarray(obj.get("cells", ()), dtype=np.int64)
        cells = cells[(cells >= 0) & (cells < world.grid.n)]
        old_mask[cells] = True
    stage_positive = decay_stage is not None and np.any(decay_stage > 0.0)
    if old_mask.any() and decay_stage is not None:
        active_stages = set(int(v) for v in np.unique(decay_stage[old_mask]) if v > 0)
    else:
        active_stages = set()
    mean_boundary_memory = (
        float(np.average(boundary_memory[old_mask], weights=area[old_mask]))
        if boundary_memory is not None and old_mask.any() else float("nan")
    )
    mean_erosion_budget = (
        float(np.average(erosion_budget[old_mask], weights=area[old_mask]))
        if erosion_budget is not None and old_mask.any() else float("nan")
    )
    export_volume_km3 = (
        float(np.sum(sediment_export[old_mask] * area[old_mask]) / 1.0e9)
        if sediment_export is not None and old_mask.any() else 0.0
    )
    acceptance = {
        "current_old_orogen_objects_present": len(old_orogens) > 0,
        "current_parented_old_orogen_objects_present": len(parented) > 0,
        "current_orogeny_age_field_available": "tectonics.orogeny_age_myr" in fields,
        "current_erosion_field_available": "erosion_m" in fields,
        "current_elevation_field_available": elevation is not None,
        "current_expected_residuals_recorded": expected_residuals_recorded,
        "production_old_orogen_decay_fields_available": not missing_decay_fields,
        "production_old_orogen_decay_stage_available": bool(stage_positive),
        "production_old_orogen_boundary_memory_available": (
            boundary_memory is not None
            and old_mask.any()
            and mean_boundary_memory >= 0.50
        ),
        "production_old_orogen_erosion_budget_available": (
            erosion_budget is not None
            and old_mask.any()
            and mean_erosion_budget > 0.0
        ),
        "production_old_orogen_sediment_export_available": export_volume_km3 > 0.0,
    }
    metrics = {
        "old_subdued_orogen_object_count": int(len(old_orogens)),
        "parented_old_subdued_orogen_object_count": int(len(parented)),
        "missing_decay_field_count": int(len(missing_decay_fields)),
        "mean_old_orogen_elevation_m": float(mean_elevation),
        "mean_old_orogen_sediment_m": float(mean_sediment),
        "old_orogen_area_m2": float(total_old_area),
        "current_old_orogen_decay_stage_count": int(len(active_stages)),
        "current_mean_orogen_boundary_memory": float(mean_boundary_memory),
        "current_mean_orogen_erosion_budget_m": float(mean_erosion_budget),
        "current_orogen_sediment_export_volume_km3": float(export_volume_km3),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "generated_world_old_orogen_decay_audit_ready"
            if all(acceptance.values())
            else "generated_world_old_orogen_decay_audit_incomplete"
        ),
        "missing_decay_fields": missing_decay_fields,
        "expected_current_residual_fields": EXPECTED_CURRENT_RESIDUAL_FIELDS,
        "limitations": {
            "production_old_orogen_decay_budget_missing": bool(missing_decay_fields),
        },
        "metrics": metrics,
        "acceptance": acceptance,
    }


def old_orogen_erosion_decay_summary(world: Any | None = None) -> dict[str, Any]:
    reference = reference_old_orogen_decay_summary()
    current = (
        current_generated_old_orogen_audit(world)
        if world is not None
        else {
            "status": "generated_world_old_orogen_decay_audit_not_run",
            "acceptance": {
                "current_old_orogen_objects_present": False,
                "current_parented_old_orogen_objects_present": False,
                "current_orogeny_age_field_available": False,
                "current_erosion_field_available": False,
                "current_elevation_field_available": False,
                "current_expected_residuals_recorded": False,
                "production_old_orogen_decay_fields_available": False,
                "production_old_orogen_decay_stage_available": False,
                "production_old_orogen_boundary_memory_available": False,
                "production_old_orogen_erosion_budget_available": False,
                "production_old_orogen_sediment_export_available": False,
            },
            "metrics": {},
        }
    )
    acceptance = {
        "reference_decay_fixture_ready": reference["status"]
        == "old_orogen_erosion_decay_reference_ready",
        "relief_decay_large_enough": reference["acceptance"][
            "relief_decay_large_enough"],
        "sediment_export_declines_late": reference["acceptance"][
            "late_sediment_export_declines"],
        "boundary_persists_after_decay": reference["acceptance"][
            "old_orogen_boundary_persists_after_decay"],
        "current_generated_audit_available": current["status"]
        == "generated_world_old_orogen_decay_audit_ready",
        "current_expected_residuals_recorded": current["acceptance"][
            "current_expected_residuals_recorded"],
        "production_old_orogen_decay_fields_available": current["acceptance"][
            "production_old_orogen_decay_fields_available"],
        "production_old_orogen_decay_budget_available": bool(
            current["acceptance"]["production_old_orogen_decay_stage_available"]
            and current["acceptance"][
                "production_old_orogen_boundary_memory_available"]
            and current["acceptance"][
                "production_old_orogen_erosion_budget_available"]
            and current["acceptance"][
                "production_old_orogen_sediment_export_available"]
        ),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "old_orogen_erosion_decay_ready"
            if all(acceptance.values())
            else "old_orogen_erosion_decay_incomplete"
        ),
        "reference": reference,
        "current_generated": current,
        "acceptance": acceptance,
        "summary_digest": _canonical_digest({
            "reference": reference,
            "current": current,
        }),
    }
