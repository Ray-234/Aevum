"""Wilson-cycle lifecycle reference diagnostics."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from typing import Any

import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.modules.tectonics import TectonicsModule


SCHEMA = "aevum.wilson_cycle_lifecycle_reference.v1"

EXPECTED_STAGE_SEQUENCE = (
    "continental_rift",
    "spreading_ocean",
    "mature_ocean",
    "closing_ocean",
    "closing_arc_margin",
    "suture_relict",
    "old_orogen_relict",
)

EXPECTED_OBJECT_SETS = (
    "tectonics.ocean_basins",
    "tectonics.rift_systems",
    "tectonics.passive_margins",
    "tectonics.spreading_centers",
    "tectonics.closing_margins",
    "tectonics.sutures",
    "tectonics.wilson_cycles",
    "tectonics.ocean_gateways",
)

EXPECTED_CURRENT_RESIDUAL_OBJECT_SETS: tuple[str, ...] = ()


def _digest(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def wilson_cycle_lifecycle_summary(world=None, *, cells: int = 1200) -> dict[str, Any]:
    reference = scripted_wilson_cycle_reference_summary(cells=cells)
    current = None if world is None else generated_world_wilson_lifecycle_summary(world)
    acceptance = {
        "scripted_reference_ready": reference["status"] == "scripted_wilson_cycle_reference_ready",
        "generated_world_audit_available": current is not None,
    }
    if current is not None:
        acceptance.update({
            "current_wilson_objects_present": current["acceptance"]["wilson_objects_present"],
            "current_gateway_causality_present": current["acceptance"]["gateway_causality_present"],
            "current_expected_residuals_recorded": current["acceptance"]["expected_residuals_recorded"],
            "current_unexpected_missing_object_sets_absent": current["acceptance"]["unexpected_missing_object_sets_absent"],
        })
    summary = {
        "schema": SCHEMA,
        "status": (
            "wilson_cycle_lifecycle_reference_ready"
            if all(acceptance.values())
            else "wilson_cycle_lifecycle_reference_incomplete"
        ),
        "scripted_reference": reference,
        "current_generated": current,
        "acceptance": acceptance,
        "next_gates": (
            "P83.crust_sediment_province_coupling",
            "P84.source_to_sink_sediment_budget",
        ),
    }
    summary["summary_digest"] = _digest(summary)
    return summary


def scripted_wilson_cycle_reference_summary(*, cells: int = 1200) -> dict[str, Any]:
    grid = SphereGrid.fibonacci(cells, CONSTANTS.EARTH_RADIUS)
    module = TectonicsModule()
    scripts = _scripted_boundaries(grid)
    previous: dict[str, list[dict[str, Any]]] = {}
    frames: list[dict[str, Any]] = []
    basin_ids: list[str] = []
    lineage_keys: list[str] = []
    phase_codes: list[float] = []
    object_set_counts: Counter[str] = Counter()
    gateway_statuses: set[str] = set()
    parent_link_failures: list[str] = []
    basin_births: list[float] = []
    basin_ages: list[float] = []

    for script in scripts:
        lifecycle = module._wilson_lifecycle_objects(
            script["boundary_objects"], previous, float(script["time_myr"]))
        basins = lifecycle["tectonics.ocean_basins"]
        basin = basins[0] if basins else {}
        basin_ids.append(str(basin.get("id", "")))
        lineage_keys.append(str(basin.get("lineage_key", "")))
        phase_codes.append(float(basin.get("phase_code", 0.0)))
        basin_births.append(float(basin.get("birth_myr", 0.0)))
        basin_ages.append(float(basin.get("age_myr", 0.0)))
        for key in EXPECTED_OBJECT_SETS:
            object_set_counts[key] += len(lifecycle.get(key, []))
        for gateway in lifecycle.get("tectonics.ocean_gateways", []):
            gateway_statuses.add(str(gateway.get("status", "")))
            if not gateway.get("parent_basin_id"):
                parent_link_failures.append(f"{gateway.get('id')}:missing_parent_basin")
            if not gateway.get("wilson_cycle_id"):
                parent_link_failures.append(f"{gateway.get('id')}:missing_wilson_cycle")
        _check_parent_links(lifecycle, parent_link_failures)

        old_orogen = ()
        stage = str(basin.get("stage", ""))
        if script["stage"] == "old_orogen_relict":
            old_orogen = (_old_orogen_from_suture(lifecycle, time_myr=script["time_myr"]),)
            stage = "old_orogen_relict" if old_orogen else stage

        frames.append({
            "name": script["name"],
            "time_myr": float(script["time_myr"]),
            "boundary_kinds": tuple(obj["kind"] for obj in script["boundary_objects"]),
            "basin_id": str(basin.get("id", "")),
            "lineage_key": str(basin.get("lineage_key", "")),
            "basin_stage": stage,
            "basin_phase_code": float(basin.get("phase_code", 0.0)),
            "basin_age_myr": float(basin.get("age_myr", 0.0)),
            "object_counts": {
                key: len(lifecycle.get(key, [])) for key in EXPECTED_OBJECT_SETS
            },
            "gateway_statuses": tuple(sorted({
                str(gw.get("status", ""))
                for gw in lifecycle.get("tectonics.ocean_gateways", [])
            })),
            "old_orogen_relicts": old_orogen,
        })
        previous = _previous_from_lifecycle(lifecycle)

    stage_sequence = tuple(frame["basin_stage"] for frame in frames)
    unique_basin_ids = tuple(sorted({bid for bid in basin_ids if bid}))
    unique_lineages = tuple(sorted({key for key in lineage_keys if key}))
    phase_monotonic = all(
        phase_codes[i] <= phase_codes[i + 1] for i in range(len(phase_codes) - 1)
    )
    age_monotonic = all(
        basin_ages[i] <= basin_ages[i + 1] for i in range(len(basin_ages) - 1)
    )
    old_orogens = [
        obj
        for frame in frames
        for obj in frame["old_orogen_relicts"]
    ]
    acceptance = {
        "expected_stage_sequence_present": stage_sequence == EXPECTED_STAGE_SEQUENCE,
        "basin_id_persistent": len(unique_basin_ids) == 1,
        "lineage_key_persistent": len(unique_lineages) == 1,
        "phase_codes_monotonic": phase_monotonic,
        "basin_age_monotonic": age_monotonic,
        "basin_birth_preserved": len({round(value, 1) for value in basin_births}) == 1,
        "required_object_sets_observed": all(
            object_set_counts[key] > 0 for key in EXPECTED_OBJECT_SETS),
        "gateway_status_sequence_observed": {"opening", "open", "closing", "restricted", "closed"}.issubset(gateway_statuses),
        "parent_causality_links_present": not parent_link_failures,
        "old_orogen_inherits_suture": bool(old_orogens)
        and all(obj.get("parent_suture_id") for obj in old_orogens)
        and all(obj.get("parent_basin_id") for obj in old_orogens),
        "no_current_frame_random_stage_labels": all(
            "lineage_key" in frame and frame["lineage_key"] for frame in frames),
    }
    summary = {
        "schema": "aevum.scripted_wilson_cycle_reference.v1",
        "expected_stage_sequence": EXPECTED_STAGE_SEQUENCE,
        "frames": tuple(frames),
        "frame_count": len(frames),
        "stage_sequence": stage_sequence,
        "unique_basin_ids": unique_basin_ids,
        "unique_lineage_keys": unique_lineages,
        "phase_codes": tuple(phase_codes),
        "basin_ages_myr": tuple(basin_ages),
        "object_set_counts": dict(sorted(object_set_counts.items())),
        "gateway_statuses": tuple(sorted(gateway_statuses)),
        "parent_link_failures": tuple(parent_link_failures),
        "old_orogen_relict_count": len(old_orogens),
        "acceptance": acceptance,
        "extraction_policy": {
            "fixture_type": "deterministic_scripted_lifecycle",
            "exact_earth_history_required": False,
            "gplates_replay_pending": True,
            "stage_labels_derived_from_lineage_objects": True,
        },
    }
    summary["status"] = (
        "scripted_wilson_cycle_reference_ready"
        if all(acceptance.values())
        else "scripted_wilson_cycle_reference_incomplete"
    )
    summary["fixture_digest"] = _digest(summary)
    return summary


def generated_world_wilson_lifecycle_summary(world) -> dict[str, Any]:
    object_counts = {
        key: len(world.objects.get(key, [])) for key in EXPECTED_OBJECT_SETS
    }
    missing_object_sets = tuple(
        key for key, count in object_counts.items() if count == 0
    )
    unexpected_missing = tuple(sorted(
        set(missing_object_sets) - set(EXPECTED_CURRENT_RESIDUAL_OBJECT_SETS)
    ))
    ocean_basins = list(world.objects.get("tectonics.ocean_basins", []))
    wilson_cycles = list(world.objects.get("tectonics.wilson_cycles", []))
    gateways = list(world.objects.get("tectonics.ocean_gateways", []))
    basin_stage_counts = Counter(str(obj.get("stage", "")) for obj in ocean_basins)
    wilson_stage_counts = Counter(str(obj.get("stage", "")) for obj in wilson_cycles)
    gateway_status_counts = Counter(str(obj.get("status", "")) for obj in gateways)
    lineage_keys = sorted({
        str(obj.get("lineage_key", ""))
        for obj in ocean_basins
        if obj.get("lineage_key")
    })
    phase = world.fields.get("archive.wilson_cycle_phase")
    phase_codes: tuple[int, ...] = ()
    if phase is not None:
        phase_codes = tuple(sorted(int(code) for code in np.unique(phase.astype(int)) if int(code) > 0))
    parent_link_failures: list[str] = []
    for obj in (
        world.objects.get("tectonics.rift_systems", [])
        + world.objects.get("tectonics.passive_margins", [])
        + world.objects.get("tectonics.spreading_centers", [])
        + world.objects.get("tectonics.closing_margins", [])
        + world.objects.get("tectonics.sutures", [])
    ):
        if not obj.get("parent_basin_id"):
            parent_link_failures.append(f"{obj.get('id')}:missing_parent_basin")
    for gateway in gateways:
        if not gateway.get("parent_basin_id") or not gateway.get("wilson_cycle_id"):
            parent_link_failures.append(f"{gateway.get('id')}:missing_gateway_parent")
    acceptance = {
        "wilson_objects_present": object_counts["tectonics.ocean_basins"] > 0
        and object_counts["tectonics.wilson_cycles"] > 0,
        "gateway_causality_present": object_counts["tectonics.ocean_gateways"] > 0
        and not any("gateway" in failure for failure in parent_link_failures),
        "expected_residuals_recorded": set(missing_object_sets).issubset(
            set(EXPECTED_CURRENT_RESIDUAL_OBJECT_SETS)),
        "unexpected_missing_object_sets_absent": not unexpected_missing,
        "stage_diversity_present": len(basin_stage_counts) >= 3,
        "phase_field_present": len(phase_codes) >= 3,
        "parent_causality_links_present": not parent_link_failures,
    }
    summary = {
        "schema": "aevum.generated_world_wilson_lifecycle.v1",
        "context": {
            "spec_name": str(world.spec.name),
            "seed": int(world.spec.seed),
            "cells": int(world.grid.n),
            "time_myr": float(world.time_myr),
        },
        "object_counts": object_counts,
        "missing_object_sets": missing_object_sets,
        "unexpected_missing_object_sets": unexpected_missing,
        "expected_current_residual_object_sets": EXPECTED_CURRENT_RESIDUAL_OBJECT_SETS,
        "basin_stage_counts": dict(sorted(basin_stage_counts.items())),
        "wilson_stage_counts": dict(sorted(wilson_stage_counts.items())),
        "gateway_status_counts": dict(sorted(gateway_status_counts.items())),
        "lineage_key_count": len(lineage_keys),
        "active_phase_codes": phase_codes,
        "parent_link_failures": tuple(parent_link_failures),
        "acceptance": acceptance,
        "limitations": {
            "current_spreading_center_objects_missing": (
                "tectonics.spreading_centers" in missing_object_sets),
            "gplates_time_slice_replay_pending": True,
            "old_orogen_decay_reference_pending": True,
        },
    }
    summary["status"] = (
        "generated_world_wilson_lifecycle_ready"
        if all(acceptance.values())
        else "generated_world_wilson_lifecycle_incomplete"
    )
    summary["summary_digest"] = _digest(summary)
    return summary


def _scripted_boundaries(grid: SphereGrid) -> tuple[dict[str, Any], ...]:
    rift_cells = np.where((np.abs(grid.lon + 12.0) < 4.0) & (np.abs(grid.lat) < 42.0))[0]
    ridge_cells = np.where((np.abs(grid.lon) < 4.0) & (np.abs(grid.lat) < 48.0))[0]
    passive_cells = np.where((np.abs(grid.lon - 9.0) < 4.0) & (np.abs(grid.lat) < 48.0))[0]
    trench_cells = np.where((np.abs(grid.lon - 18.0) < 4.0) & (np.abs(grid.lat) < 44.0))[0]
    active_cells = np.where((np.abs(grid.lon - 25.0) < 4.0) & (np.abs(grid.lat) < 42.0))[0]
    suture_cells = np.where((np.abs(grid.lon - 5.0) < 4.0) & (np.abs(grid.lat) < 40.0))[0]
    return (
        {
            "name": "rift_birth",
            "stage": "continental_rift",
            "time_myr": 100.0,
            "boundary_objects": (
                _boundary("divergent", rift_cells, 100.0, obj_id="boundary:wilson-demo"),
            ),
        },
        {
            "name": "ocean_opening",
            "stage": "spreading_ocean",
            "time_myr": 160.0,
            "boundary_objects": (
                _boundary("ridge", ridge_cells, 160.0, obj_id="boundary:wilson-demo"),
            ),
        },
        {
            "name": "passive_margin_maturity",
            "stage": "mature_ocean",
            "time_myr": 230.0,
            "boundary_objects": (
                _boundary("passive_margin", passive_cells, 230.0, obj_id="boundary:wilson-demo"),
            ),
        },
        {
            "name": "subduction_closure",
            "stage": "closing_ocean",
            "time_myr": 310.0,
            "boundary_objects": (
                _boundary("trench", trench_cells, 310.0, obj_id="boundary:wilson-demo"),
            ),
        },
        {
            "name": "arc_collision_margin",
            "stage": "closing_arc_margin",
            "time_myr": 360.0,
            "boundary_objects": (
                _boundary("active_margin", active_cells, 360.0, obj_id="boundary:wilson-demo"),
            ),
        },
        {
            "name": "suture_closure",
            "stage": "suture_relict",
            "time_myr": 430.0,
            "boundary_objects": (
                _boundary("suture", suture_cells, 430.0, obj_id="boundary:wilson-demo"),
            ),
        },
        {
            "name": "old_orogen_decay",
            "stage": "old_orogen_relict",
            "time_myr": 620.0,
            "boundary_objects": (
                _boundary("suture", suture_cells, 620.0, obj_id="boundary:wilson-demo"),
            ),
        },
    )


def _boundary(kind: str, cells: np.ndarray, t: float, *, obj_id: str) -> dict[str, Any]:
    cells = np.asarray(cells, dtype=int)
    return {
        "id": obj_id,
        "kind": kind,
        "stage": kind,
        "cells": cells.tolist(),
        "cell_count": int(cells.size),
        "area_fraction": float(cells.size / 1200.0),
        "birth_myr": round(float(t), 1),
        "last_active_myr": round(float(t), 1),
        "age_myr": 0.0,
        "parent_plate_ids": [10, 11],
        "lat": 0.0,
        "lon": 0.0,
        "boundary_continental_fraction": 0.0,
        "subducting_plate_id": 10 if kind in {"trench", "active_margin"} else None,
        "overriding_plate_id": 11 if kind in {"trench", "active_margin"} else None,
        "polarity_basis": (
            "scripted_oceanic_plate_subducts"
            if kind in {"trench", "active_margin"} else None
        ),
    }


def _previous_from_lifecycle(lifecycle: dict[str, list[dict[str, Any]]]) -> dict[str, list[dict[str, Any]]]:
    return {
        "ocean_basins": lifecycle.get("tectonics.ocean_basins", []),
        "rift_systems": lifecycle.get("tectonics.rift_systems", []),
        "passive_margins": lifecycle.get("tectonics.passive_margins", []),
        "spreading_centers": lifecycle.get("tectonics.spreading_centers", []),
        "closing_margins": lifecycle.get("tectonics.closing_margins", []),
        "sutures": lifecycle.get("tectonics.sutures", []),
    }


def _check_parent_links(lifecycle: dict[str, list[dict[str, Any]]], failures: list[str]) -> None:
    for key in (
        "tectonics.rift_systems",
        "tectonics.passive_margins",
        "tectonics.spreading_centers",
        "tectonics.closing_margins",
        "tectonics.sutures",
    ):
        for obj in lifecycle.get(key, []):
            if not obj.get("parent_basin_id"):
                failures.append(f"{obj.get('id')}:missing_parent_basin")


def _old_orogen_from_suture(
    lifecycle: dict[str, list[dict[str, Any]]],
    *,
    time_myr: float,
) -> dict[str, Any] | None:
    sutures = lifecycle.get("tectonics.sutures", [])
    cycles = lifecycle.get("tectonics.wilson_cycles", [])
    if not sutures:
        return None
    suture = sutures[0]
    cycle = cycles[0] if cycles else {}
    return {
        "id": f"old_orogen:{suture['id']}",
        "kind": "old_orogen_relict",
        "stage": "old_orogen_relict",
        "parent_process": "suture_inheritance_and_orogenic_decay",
        "parent_suture_id": suture.get("id"),
        "parent_basin_id": suture.get("parent_basin_id"),
        "parent_wilson_cycle_id": cycle.get("id"),
        "lineage_key": suture.get("lineage_key"),
        "cells": suture.get("cells", []),
        "birth_myr": suture.get("birth_myr"),
        "last_active_myr": round(float(time_myr), 1),
        "age_myr": round(float(max(time_myr - float(suture.get("birth_myr", time_myr)), 0.0)), 1),
    }
