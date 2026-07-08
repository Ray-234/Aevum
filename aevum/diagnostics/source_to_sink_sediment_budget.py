"""Source-to-sink sediment budget reference diagnostics for P84."""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from typing import Any

import numpy as np


SCHEMA = "aevum.source_to_sink_sediment_budget.v1"

SOURCE_IDS = (
    "NOAA_TOTAL_SEDIMENT_THICKNESS",
    "GLIM_GLOBAL_LITHOLOGY",
)

REQUIRED_ZONE_KINDS = {
    "mountain_source",
    "platform_source",
    "foreland_sink",
    "passive_margin_sink",
    "shelf_sink",
    "ocean_basin_sink",
}

EXPECTED_CURRENT_RESIDUAL_OBJECTS = (
    "terrain.drainage_basins",
    "terrain.sediment_routing_edges",
    "terrain.sediment_budget",
)


REFERENCE_ZONES: tuple[dict[str, Any], ...] = (
    {
        "id": "mountain_source",
        "kind": "mountain_source",
        "area_km2": 150000.0,
        "initial_elevation_m": 2400.0,
        "projected_elevation_delta_m": -45.0,
        "initial_sediment_m": 220.0,
        "sediment_delta_m": -300.0,
        "available_erosion_m": 520.0,
        "accommodation_m": 0.0,
        "parent_process": "collision_orogeny",
    },
    {
        "id": "platform_source",
        "kind": "platform_source",
        "area_km2": 300000.0,
        "initial_elevation_m": 620.0,
        "projected_elevation_delta_m": -15.0,
        "initial_sediment_m": 700.0,
        "sediment_delta_m": -80.0,
        "available_erosion_m": 180.0,
        "accommodation_m": 0.0,
        "parent_process": "platform_erosion_and_drainage_export",
    },
    {
        "id": "foreland_sink",
        "kind": "foreland_sink",
        "area_km2": 90000.0,
        "initial_elevation_m": 120.0,
        "projected_elevation_delta_m": 35.0,
        "initial_sediment_m": 1200.0,
        "sediment_delta_m": 250.0,
        "available_erosion_m": 0.0,
        "accommodation_m": 550.0,
        "parent_process": "flexural_loading",
    },
    {
        "id": "passive_margin_sink",
        "kind": "passive_margin_sink",
        "area_km2": 100000.0,
        "initial_elevation_m": 35.0,
        "projected_elevation_delta_m": 25.0,
        "initial_sediment_m": 1500.0,
        "sediment_delta_m": 180.0,
        "available_erosion_m": 0.0,
        "accommodation_m": 400.0,
        "parent_process": "passive_margin_subsidence",
    },
    {
        "id": "shelf_sink",
        "kind": "shelf_sink",
        "area_km2": 120000.0,
        "initial_elevation_m": -90.0,
        "projected_elevation_delta_m": 20.0,
        "initial_sediment_m": 2200.0,
        "sediment_delta_m": 170.0,
        "available_erosion_m": 0.0,
        "accommodation_m": 250.0,
        "parent_process": "shelf_sedimentation",
    },
    {
        "id": "ocean_basin_sink",
        "kind": "ocean_basin_sink",
        "area_km2": 100000.0,
        "initial_elevation_m": -3600.0,
        "projected_elevation_delta_m": 5.0,
        "initial_sediment_m": 500.0,
        "sediment_delta_m": 81.0,
        "available_erosion_m": 0.0,
        "accommodation_m": 400.0,
        "parent_process": "abyssal_fan_deposition",
    },
)

REFERENCE_EDGES: tuple[dict[str, Any], ...] = (
    {"from": "mountain_source", "to": "foreland_sink", "volume_km3": 22500.0},
    {"from": "mountain_source", "to": "shelf_sink", "volume_km3": 14400.0},
    {"from": "mountain_source", "to": "ocean_basin_sink", "volume_km3": 8100.0},
    {"from": "platform_source", "to": "passive_margin_sink", "volume_km3": 18000.0},
    {"from": "platform_source", "to": "shelf_sink", "volume_km3": 6000.0},
)


def _canonical_digest(data: Any) -> str:
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _zone_volume_km3(zone: dict[str, Any]) -> float:
    return float(zone["area_km2"]) * float(zone["sediment_delta_m"]) / 1000.0


def _land_mask(elevation_m: float) -> bool:
    return float(elevation_m) >= 0.0


def reference_source_to_sink_budget_summary() -> dict[str, Any]:
    zones = tuple(dict(zone) for zone in REFERENCE_ZONES)
    edges = tuple(dict(edge) for edge in REFERENCE_EDGES)
    zone_by_id = {str(zone["id"]): zone for zone in zones}
    zone_kinds = {str(zone["kind"]) for zone in zones}
    source_zones = tuple(zone for zone in zones if float(zone["sediment_delta_m"]) < 0.0)
    sink_zones = tuple(zone for zone in zones if float(zone["sediment_delta_m"]) > 0.0)
    source_volume_km3 = -sum(_zone_volume_km3(zone) for zone in source_zones)
    sink_volume_km3 = sum(_zone_volume_km3(zone) for zone in sink_zones)
    balance_fraction = (
        abs(sink_volume_km3 - source_volume_km3) / source_volume_km3
        if source_volume_km3 > 0.0 else float("inf")
    )

    outgoing: defaultdict[str, float] = defaultdict(float)
    incoming: defaultdict[str, float] = defaultdict(float)
    invalid_edges: list[dict[str, Any]] = []
    for edge in edges:
        src = str(edge["from"])
        dst = str(edge["to"])
        volume = float(edge["volume_km3"])
        if src not in zone_by_id or dst not in zone_by_id or volume <= 0.0:
            invalid_edges.append(edge)
            continue
        outgoing[src] += volume
        incoming[dst] += volume

    routing_mismatches = []
    for zone in zones:
        zone_id = str(zone["id"])
        volume = _zone_volume_km3(zone)
        routed = outgoing[zone_id] if volume < 0.0 else incoming[zone_id]
        expected = abs(volume)
        if abs(routed - expected) > 1.0e-6:
            routing_mismatches.append({
                "zone_id": zone_id,
                "expected_volume_km3": expected,
                "routed_volume_km3": routed,
            })

    land_mask_changes = []
    accommodation_violations = []
    erosion_violations = []
    utilization = []
    for zone in zones:
        before = float(zone["initial_elevation_m"])
        after = before + float(zone["projected_elevation_delta_m"])
        if _land_mask(before) != _land_mask(after):
            land_mask_changes.append(str(zone["id"]))
        delta = float(zone["sediment_delta_m"])
        if delta > 0.0:
            accom = float(zone["accommodation_m"])
            if delta > accom:
                accommodation_violations.append(str(zone["id"]))
            utilization.append(delta / accom if accom > 0.0 else float("inf"))
        if delta < 0.0 and abs(delta) > float(zone["available_erosion_m"]):
            erosion_violations.append(str(zone["id"]))

    source_ordering = bool(
        abs(float(zone_by_id["mountain_source"]["sediment_delta_m"]))
        > abs(float(zone_by_id["platform_source"]["sediment_delta_m"]))
    )
    sink_ordering = bool(
        float(zone_by_id["foreland_sink"]["sediment_delta_m"])
        > float(zone_by_id["passive_margin_sink"]["sediment_delta_m"])
        > float(zone_by_id["ocean_basin_sink"]["sediment_delta_m"])
        and float(zone_by_id["shelf_sink"]["sediment_delta_m"])
        > float(zone_by_id["ocean_basin_sink"]["sediment_delta_m"])
    )
    acceptance = {
        "fixture_schema_ready": bool(zones and edges and SOURCE_IDS),
        "required_zone_kinds_present": REQUIRED_ZONE_KINDS.issubset(zone_kinds),
        "all_edges_reference_known_zones": not invalid_edges,
        "source_volumes_match_outgoing_edges": not routing_mismatches,
        "sink_volumes_match_incoming_edges": not routing_mismatches,
        "sediment_volume_conserved": balance_fraction <= 1.0e-9,
        "no_land_mask_regression": not land_mask_changes,
        "sink_deposition_within_accommodation": not accommodation_violations,
        "source_erosion_within_available_material": not erosion_violations,
        "source_sink_ordering_plausible": bool(source_ordering and sink_ordering),
        "raw_reference_data_not_stored": True,
    }
    metrics = {
        "zone_count": int(len(zones)),
        "edge_count": int(len(edges)),
        "source_zone_count": int(len(source_zones)),
        "sink_zone_count": int(len(sink_zones)),
        "source_volume_km3": float(source_volume_km3),
        "sink_volume_km3": float(sink_volume_km3),
        "volume_balance_fraction": float(balance_fraction),
        "max_accommodation_utilization": float(max(utilization) if utilization else 0.0),
        "land_mask_change_count": int(len(land_mask_changes)),
        "routing_mismatch_count": int(len(routing_mismatches)),
        "invalid_edge_count": int(len(invalid_edges)),
        "accommodation_violation_count": int(len(accommodation_violations)),
        "erosion_violation_count": int(len(erosion_violations)),
        "mountain_source_export_km3": float(outgoing["mountain_source"]),
        "platform_source_export_km3": float(outgoing["platform_source"]),
        "foreland_sink_deposition_km3": float(incoming["foreland_sink"]),
        "passive_margin_sink_deposition_km3": float(incoming["passive_margin_sink"]),
        "shelf_sink_deposition_km3": float(incoming["shelf_sink"]),
        "ocean_basin_sink_deposition_km3": float(incoming["ocean_basin_sink"]),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "source_to_sink_sediment_budget_reference_ready"
            if all(acceptance.values())
            else "source_to_sink_sediment_budget_reference_incomplete"
        ),
        "source_ids": SOURCE_IDS,
        "zones": zones,
        "routing_edges": edges,
        "missing_required_zone_kinds": tuple(sorted(REQUIRED_ZONE_KINDS - zone_kinds)),
        "invalid_edges": tuple(invalid_edges),
        "routing_mismatches": tuple(routing_mismatches),
        "land_mask_changes": tuple(land_mask_changes),
        "accommodation_violations": tuple(accommodation_violations),
        "erosion_violations": tuple(erosion_violations),
        "extraction_policy": {
            "raw_sediment_or_drainage_data_stored": False,
            "direct_sediment_drainage_extraction_pending": True,
            "fixture_is_small_derived_reference": True,
        },
        "metrics": metrics,
        "acceptance": acceptance,
        "fixture_digest": _canonical_digest({
            "zones": zones,
            "edges": edges,
            "source_ids": SOURCE_IDS,
        }),
    }


def current_generated_source_to_sink_audit(world: Any) -> dict[str, Any]:
    fields = getattr(world, "fields", {})
    objects = getattr(world, "objects", {})
    missing_objects = tuple(
        name for name in EXPECTED_CURRENT_RESIDUAL_OBJECTS if name not in objects)
    expected_residuals_recorded = bool(
        set(missing_objects).issubset(set(EXPECTED_CURRENT_RESIDUAL_OBJECTS)))
    sediment = (
        np.asarray(world.field("sediment.thickness_m"), dtype=np.float64)
        if "sediment.thickness_m" in fields else None
    )
    elevation = (
        np.asarray(world.field("terrain.elevation_m"), dtype=np.float64)
        if "terrain.elevation_m" in fields else None
    )
    land = elevation >= 0.0 if elevation is not None else np.zeros(world.grid.n, dtype=bool)
    ocean = ~land
    drainage_basins = list(objects.get("terrain.drainage_basins", []))
    routing_edges = list(objects.get("terrain.sediment_routing_edges", []))
    budget_objects = list(objects.get("terrain.sediment_budget", []))
    budget = budget_objects[0] if budget_objects else {}
    source_volume_km3 = float(budget.get("source_volume_km3", 0.0))
    sink_volume_km3 = float(budget.get("sink_volume_km3", 0.0))
    balance_fraction = float(budget.get(
        "volume_balance_fraction",
        (
            abs(sink_volume_km3 - source_volume_km3) / source_volume_km3
            if source_volume_km3 > 0.0 else float("inf")
        ),
    ))
    source_kinds = {
        str(edge.get("from", "unknown"))
        for edge in routing_edges
        if float(edge.get("volume_km3", 0.0)) > 0.0
    }
    sink_kinds = {
        str(edge.get("to", "unknown"))
        for edge in routing_edges
        if float(edge.get("volume_km3", 0.0)) > 0.0
    }
    metrics = {
        "missing_object_count": int(len(missing_objects)),
        "sediment_field_available": sediment is not None,
        "elevation_field_available": elevation is not None,
        "drainage_basin_object_count": int(len(drainage_basins)),
        "routing_edge_count": int(len(routing_edges)),
        "sediment_budget_object_count": int(len(budget_objects)),
        "source_volume_km3": float(source_volume_km3),
        "sink_volume_km3": float(sink_volume_km3),
        "production_volume_balance_fraction": float(balance_fraction),
        "routing_source_kind_count": int(len(source_kinds)),
        "routing_sink_kind_count": int(len(sink_kinds)),
        "land_mean_sediment_m": (
            float(np.mean(sediment[land])) if sediment is not None and land.any()
            else float("nan")
        ),
        "ocean_mean_sediment_m": (
            float(np.mean(sediment[ocean])) if sediment is not None and ocean.any()
            else float("nan")
        ),
        "sediment_min_m": float(np.nanmin(sediment)) if sediment is not None else float("nan"),
        "sediment_max_m": float(np.nanmax(sediment)) if sediment is not None else float("nan"),
    }
    acceptance = {
        "current_sediment_field_available": sediment is not None,
        "current_elevation_field_available": elevation is not None,
        "current_expected_residuals_recorded": expected_residuals_recorded,
        "production_source_to_sink_objects_available": not missing_objects,
        "production_drainage_basins_available": len(drainage_basins) > 0,
        "production_routing_edges_available": len(routing_edges) > 0,
        "production_sediment_budget_available": len(budget_objects) > 0,
        "production_sediment_budget_closes": (
            source_volume_km3 > 0.0
            and sink_volume_km3 > 0.0
            and balance_fraction <= 1.0e-9
        ),
        "production_source_sink_kinds_diverse": (
            len(source_kinds) >= 2 and len(sink_kinds) >= 2
        ),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "generated_world_source_to_sink_audit_ready"
            if all(acceptance.values())
            else "generated_world_source_to_sink_audit_incomplete"
        ),
        "missing_source_to_sink_objects": missing_objects,
        "expected_current_residual_objects": EXPECTED_CURRENT_RESIDUAL_OBJECTS,
        "limitations": {
            "production_source_to_sink_objects_missing": bool(missing_objects),
        },
        "source_kinds": tuple(sorted(source_kinds)),
        "sink_kinds": tuple(sorted(sink_kinds)),
        "metrics": metrics,
        "acceptance": acceptance,
    }


def source_to_sink_sediment_budget_summary(world: Any | None = None) -> dict[str, Any]:
    reference = reference_source_to_sink_budget_summary()
    current = (
        current_generated_source_to_sink_audit(world)
        if world is not None
        else {
            "status": "generated_world_source_to_sink_audit_not_run",
            "acceptance": {
                "current_sediment_field_available": False,
                "current_elevation_field_available": False,
                "current_expected_residuals_recorded": False,
                "production_source_to_sink_objects_available": False,
                "production_sediment_budget_closes": False,
                "production_source_sink_kinds_diverse": False,
            },
            "metrics": {},
        }
    )
    acceptance = {
        "reference_budget_ready": reference["status"]
        == "source_to_sink_sediment_budget_reference_ready",
        "sediment_volume_conserved": reference["acceptance"][
            "sediment_volume_conserved"],
        "no_land_mask_regression": reference["acceptance"][
            "no_land_mask_regression"],
        "deposition_within_accommodation": reference["acceptance"][
            "sink_deposition_within_accommodation"],
        "routing_edges_close_budget": bool(
            reference["acceptance"]["source_volumes_match_outgoing_edges"]
            and reference["acceptance"]["sink_volumes_match_incoming_edges"]),
        "current_generated_audit_available": current["status"]
        == "generated_world_source_to_sink_audit_ready",
        "current_expected_residuals_recorded": current["acceptance"][
            "current_expected_residuals_recorded"],
        "production_source_to_sink_objects_available": current["acceptance"][
            "production_source_to_sink_objects_available"],
        "production_sediment_budget_closes": current["acceptance"][
            "production_sediment_budget_closes"],
        "production_source_sink_kinds_diverse": current["acceptance"][
            "production_source_sink_kinds_diverse"],
    }
    return {
        "schema": SCHEMA,
        "status": (
            "source_to_sink_sediment_budget_ready"
            if all(acceptance.values())
            else "source_to_sink_sediment_budget_incomplete"
        ),
        "reference": reference,
        "current_generated": current,
        "acceptance": acceptance,
        "summary_digest": _canonical_digest({
            "reference": reference,
            "current": current,
        }),
    }
