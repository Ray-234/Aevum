"""GMBA-style mountain inventory expression diagnostics."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from typing import Any

import numpy as np


SCHEMA = "aevum.mountain_inventory_expression.v1"

SOURCE_IDS = (
    "GMBA_MOUNTAIN_INVENTORY",
    "NOAA_ETOPO_2022",
    "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
)

MOUNTAIN_LANDFORM_KINDS = (
    "orogen",
    "old_subdued_orogen",
    "plateau",
    "arc_microcontinent",
)

EXPECTED_CURRENT_RESIDUAL_FIELDS = (
    "terrain.mountain_ranges",
    "terrain.mountain_inventory",
    "terrain.mountain_hierarchy_level",
    "tectonics.mountain_belt_id",
    "tectonics.mountain_parent_process_id",
)

EXPECTED_CURRENT_RESIDUAL_KINDS = (
    "orogen",
    "plateau",
)

REFERENCE_RANGES: tuple[dict[str, Any], ...] = (
    {
        "id": "andes_system",
        "name": "Andes-style continental margin system",
        "mountain_class": "active_margin_orogen",
        "hierarchy_level": 1,
        "parent_range_id": None,
        "parent_process": "ocean_continent_subduction",
        "area_fraction_world": 0.010,
        "elongation_ratio": 9.0,
        "relief_p90_p10_m": 3200.0,
        "mean_elevation_m": 2100.0,
        "object_backed": True,
        "threshold_only": False,
    },
    {
        "id": "central_andes_plateau",
        "name": "Central Andes plateau core",
        "mountain_class": "collision_plateau",
        "hierarchy_level": 2,
        "parent_range_id": "andes_system",
        "parent_process": "crustal_shortening_and_magmatism",
        "area_fraction_world": 0.005,
        "elongation_ratio": 2.6,
        "relief_p90_p10_m": 1800.0,
        "mean_elevation_m": 3700.0,
        "object_backed": True,
        "threshold_only": False,
    },
    {
        "id": "northern_andes_subrange",
        "name": "Northern Andes subrange",
        "mountain_class": "active_margin_orogen",
        "hierarchy_level": 3,
        "parent_range_id": "andes_system",
        "parent_process": "ocean_continent_subduction",
        "area_fraction_world": 0.003,
        "elongation_ratio": 5.4,
        "relief_p90_p10_m": 2600.0,
        "mean_elevation_m": 2300.0,
        "object_backed": True,
        "threshold_only": False,
    },
    {
        "id": "alps_himalaya_system",
        "name": "Alps-Himalaya-style collision system",
        "mountain_class": "active_collision_orogen",
        "hierarchy_level": 1,
        "parent_range_id": None,
        "parent_process": "continent_continent_collision",
        "area_fraction_world": 0.015,
        "elongation_ratio": 8.4,
        "relief_p90_p10_m": 4300.0,
        "mean_elevation_m": 2600.0,
        "object_backed": True,
        "threshold_only": False,
    },
    {
        "id": "tibetan_plateau_core",
        "name": "Tibetan-style collision plateau core",
        "mountain_class": "collision_plateau",
        "hierarchy_level": 2,
        "parent_range_id": "alps_himalaya_system",
        "parent_process": "crustal_thickening",
        "area_fraction_world": 0.008,
        "elongation_ratio": 2.2,
        "relief_p90_p10_m": 1500.0,
        "mean_elevation_m": 4500.0,
        "object_backed": True,
        "threshold_only": False,
    },
    {
        "id": "high_himalaya_subrange",
        "name": "High Himalaya subrange",
        "mountain_class": "active_collision_orogen",
        "hierarchy_level": 3,
        "parent_range_id": "alps_himalaya_system",
        "parent_process": "continent_continent_collision",
        "area_fraction_world": 0.002,
        "elongation_ratio": 6.8,
        "relief_p90_p10_m": 5000.0,
        "mean_elevation_m": 4100.0,
        "object_backed": True,
        "threshold_only": False,
    },
    {
        "id": "appalachian_caledonian_relict",
        "name": "Appalachian-Caledonian-style relict belt",
        "mountain_class": "old_subdued_orogen",
        "hierarchy_level": 1,
        "parent_range_id": None,
        "parent_process": "orogenic_decay",
        "area_fraction_world": 0.006,
        "elongation_ratio": 7.0,
        "relief_p90_p10_m": 900.0,
        "mean_elevation_m": 900.0,
        "object_backed": True,
        "threshold_only": False,
    },
    {
        "id": "ural_old_orogen",
        "name": "Ural-style old suture mountain belt",
        "mountain_class": "old_subdued_orogen",
        "hierarchy_level": 1,
        "parent_range_id": None,
        "parent_process": "suture_inheritance",
        "area_fraction_world": 0.004,
        "elongation_ratio": 6.4,
        "relief_p90_p10_m": 700.0,
        "mean_elevation_m": 850.0,
        "object_backed": True,
        "threshold_only": False,
    },
    {
        "id": "east_african_rift_shoulders",
        "name": "East-African-style rift shoulder ranges",
        "mountain_class": "rift_shoulder_range",
        "hierarchy_level": 1,
        "parent_range_id": None,
        "parent_process": "continental_extension",
        "area_fraction_world": 0.004,
        "elongation_ratio": 5.6,
        "relief_p90_p10_m": 1200.0,
        "mean_elevation_m": 1700.0,
        "object_backed": True,
        "threshold_only": False,
    },
    {
        "id": "west_pacific_arc_chain",
        "name": "West-Pacific-style volcanic arc chain",
        "mountain_class": "volcanic_arc_chain",
        "hierarchy_level": 1,
        "parent_range_id": None,
        "parent_process": "subduction_arc_magmatism",
        "area_fraction_world": 0.004,
        "elongation_ratio": 10.2,
        "relief_p90_p10_m": 2500.0,
        "mean_elevation_m": 1200.0,
        "object_backed": True,
        "threshold_only": False,
    },
    {
        "id": "basin_range_extensional_highlands",
        "name": "Basin-and-Range-style extensional highlands",
        "mountain_class": "extensional_range",
        "hierarchy_level": 1,
        "parent_range_id": None,
        "parent_process": "distributed_extension",
        "area_fraction_world": 0.005,
        "elongation_ratio": 3.8,
        "relief_p90_p10_m": 1100.0,
        "mean_elevation_m": 1500.0,
        "object_backed": True,
        "threshold_only": False,
    },
)


def _canonical_digest(data: Any) -> str:
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _percentile(values: list[float] | tuple[float, ...], q: float) -> float:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size == 0:
        return 0.0
    return float(np.percentile(arr, q))


def _circular_longitude_span_deg(lons: np.ndarray) -> float:
    lons = np.asarray(lons, dtype=np.float64)
    if lons.size <= 1:
        return 0.0
    wrapped = np.sort((lons + 360.0) % 360.0)
    gaps = np.diff(np.r_[wrapped, wrapped[0] + 360.0])
    return float(360.0 - gaps.max(initial=0.0))


def _object_cells(obj: dict[str, Any], n: int) -> np.ndarray:
    cells = np.asarray(obj.get("cells", ()), dtype=np.int64)
    if cells.size == 0:
        return cells
    return cells[(cells >= 0) & (cells < n)]


def _object_elongation(grid: Any, cells: np.ndarray) -> float:
    cells = np.asarray(cells, dtype=np.int64)
    if cells.size <= 1:
        return 1.0
    lon_span = _circular_longitude_span_deg(np.asarray(grid.lon[cells], dtype=np.float64))
    lat_span = float(np.ptp(np.asarray(grid.lat[cells], dtype=np.float64)))
    mean_lat = float(np.mean(np.asarray(grid.lat[cells], dtype=np.float64)))
    zonal = lon_span * max(float(np.cos(np.deg2rad(mean_lat))), 0.25)
    meridional = max(lat_span, 1.0e-6)
    major = max(zonal, meridional, 1.0e-6)
    minor = max(min(zonal, meridional), 1.0e-6)
    return float(major / minor)


def reference_mountain_inventory_summary() -> dict[str, Any]:
    ranges = tuple(dict(item) for item in REFERENCE_RANGES)
    ids = {str(item["id"]) for item in ranges}
    class_counts = Counter(str(item["mountain_class"]) for item in ranges)
    level_counts = Counter(int(item["hierarchy_level"]) for item in ranges)
    process_counts = Counter(str(item["parent_process"]) for item in ranges)
    area_values = tuple(float(item["area_fraction_world"]) for item in ranges)
    elongation_values = tuple(float(item["elongation_ratio"]) for item in ranges)
    relief_values = tuple(float(item["relief_p90_p10_m"]) for item in ranges)
    parent_link_failures = tuple(
        str(item["id"])
        for item in ranges
        if int(item["hierarchy_level"]) > 1
        and str(item.get("parent_range_id")) not in ids
    )
    parent_process_failures = tuple(
        str(item["id"]) for item in ranges if not item.get("parent_process")
    )
    threshold_only_ranges = tuple(
        str(item["id"]) for item in ranges if bool(item.get("threshold_only"))
    )
    object_backing_failures = tuple(
        str(item["id"]) for item in ranges if not bool(item.get("object_backed"))
    )
    acceptance = {
        "fixture_schema_ready": bool(ranges and SOURCE_IDS),
        "range_count_sufficient": len(ranges) >= 10,
        "mountain_class_diversity": len(class_counts) >= 6,
        "hierarchy_levels_present": set(level_counts) >= {1, 2, 3},
        "hierarchy_parent_links_valid": not parent_link_failures,
        "parent_processes_present": not parent_process_failures,
        "object_backed_not_threshold_only": (
            not threshold_only_ranges and not object_backing_failures
        ),
        "area_distribution_finite": (
            0.04 <= sum(area_values) <= 0.10
            and min(area_values, default=0.0) >= 0.001
            and max(area_values, default=1.0) <= 0.020
        ),
        "elongation_distribution_plausible": (
            _percentile(elongation_values, 50.0) >= 3.0
            and max(elongation_values, default=0.0) <= 12.0
            and sum(1 for value in elongation_values if value >= 5.0) >= 6
        ),
        "relief_distribution_plausible": (
            min(relief_values, default=0.0) >= 500.0
            and max(relief_values, default=0.0) >= 4000.0
        ),
        "raw_reference_data_not_stored": True,
    }
    metrics = {
        "range_count": int(len(ranges)),
        "mountain_class_count": int(len(class_counts)),
        "hierarchy_level_count": int(len(level_counts)),
        "parent_process_count": int(len(process_counts)),
        "hierarchy_parent_link_failure_count": int(len(parent_link_failures)),
        "parent_process_failure_count": int(len(parent_process_failures)),
        "threshold_only_range_count": int(len(threshold_only_ranges)),
        "object_backing_failure_count": int(len(object_backing_failures)),
        "total_mountain_area_fraction_world": float(sum(area_values)),
        "max_range_area_fraction_world": float(max(area_values, default=0.0)),
        "min_range_area_fraction_world": float(min(area_values, default=0.0)),
        "median_elongation_ratio": _percentile(elongation_values, 50.0),
        "max_elongation_ratio": float(max(elongation_values, default=0.0)),
        "elongated_range_count": int(sum(
            1 for value in elongation_values if value >= 5.0)),
        "min_relief_p90_p10_m": float(min(relief_values, default=0.0)),
        "max_relief_p90_p10_m": float(max(relief_values, default=0.0)),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "mountain_inventory_reference_ready"
            if all(acceptance.values())
            else "mountain_inventory_reference_incomplete"
        ),
        "source_ids": SOURCE_IDS,
        "ranges": ranges,
        "mountain_class_counts": dict(sorted(class_counts.items())),
        "hierarchy_level_counts": dict(sorted(level_counts.items())),
        "parent_process_counts": dict(sorted(process_counts.items())),
        "parent_link_failures": parent_link_failures,
        "parent_process_failures": parent_process_failures,
        "threshold_only_ranges": threshold_only_ranges,
        "object_backing_failures": object_backing_failures,
        "extraction_policy": {
            "raw_mountain_inventory_stored": False,
            "direct_gmba_extraction_pending": True,
            "fixture_is_small_derived_reference": True,
        },
        "metrics": metrics,
        "acceptance": acceptance,
        "fixture_digest": _canonical_digest({
            "ranges": ranges,
            "source_ids": SOURCE_IDS,
        }),
    }


def current_generated_mountain_inventory_audit(world: Any) -> dict[str, Any]:
    fields = getattr(world, "fields", {})
    objects = getattr(world, "objects", {})
    grid = world.grid
    continental_landforms = list(objects.get("terrain.continental_landforms", []))
    production_ranges = list(objects.get("terrain.mountain_ranges", []))
    landforms = production_ranges if production_ranges else continental_landforms
    mountain_objects = [
        obj for obj in landforms
        if str(obj.get("kind", "")) in set(MOUNTAIN_LANDFORM_KINDS)
    ]
    kind_counts = Counter(str(obj.get("kind", "")) for obj in mountain_objects)
    missing_kinds = tuple(
        kind for kind in EXPECTED_CURRENT_RESIDUAL_KINDS
        if kind_counts.get(kind, 0) == 0
    )
    missing_fields = tuple(
        name for name in EXPECTED_CURRENT_RESIDUAL_FIELDS if name not in fields)
    parented = [
        obj for obj in mountain_objects
        if obj.get("parent_process")
    ]
    parent_linked = [
        obj for obj in mountain_objects
        if obj.get("parent_tectonic_object_ids")
        or obj.get("parent_range_id")
        or obj.get("parent_process")
    ]
    expressed = [
        obj for obj in mountain_objects
        if float(obj.get("mean_elevation_m", -1.0e9)) > 0.0
    ]
    submerged = [
        obj for obj in mountain_objects
        if float(obj.get("mean_elevation_m", 0.0)) < -250.0
    ]
    tiny = [
        obj for obj in mountain_objects
        if int(obj.get("cell_count", 0)) <= 2
    ]
    area_values = tuple(float(obj.get("area_fraction", 0.0)) for obj in mountain_objects)
    expressed_area_values = tuple(
        float(obj.get("area_fraction", 0.0)) for obj in expressed)
    relief_values = tuple(
        float(obj.get("relief_p90_minus_p10_m", 0.0)) for obj in mountain_objects)
    elongations: list[float] = []
    object_summaries: list[dict[str, Any]] = []
    for obj in mountain_objects:
        cells = _object_cells(obj, grid.n)
        elongation = _object_elongation(grid, cells) if cells.size else 0.0
        elongations.append(elongation)
        object_summaries.append({
            "id": str(obj.get("id", "")),
            "kind": str(obj.get("kind", "")),
            "parent_process": str(obj.get("parent_process", "")),
            "parent_tectonic_object_count": int(
                len(obj.get("parent_tectonic_object_ids", ()))),
            "cell_count": int(obj.get("cell_count", 0)),
            "area_fraction": float(obj.get("area_fraction", 0.0)),
            "mean_elevation_m": float(obj.get("mean_elevation_m", 0.0)),
            "relief_p90_minus_p10_m": float(
                obj.get("relief_p90_minus_p10_m", 0.0)),
            "elongation_ratio": float(elongation),
        })
    parent_process_coverage = (
        len(parented) / len(mountain_objects) if mountain_objects else 0.0
    )
    parent_link_coverage = (
        len(parent_linked) / len(mountain_objects) if mountain_objects else 0.0
    )
    expected_residuals_recorded = (
        set(missing_fields).issubset(set(EXPECTED_CURRENT_RESIDUAL_FIELDS))
        and set(missing_kinds).issubset(set(EXPECTED_CURRENT_RESIDUAL_KINDS))
    )
    production_fields_available = not missing_fields
    production_expected_kinds_available = not missing_kinds
    total_area = float(sum(area_values))
    max_area = float(max(area_values, default=0.0))
    acceptance = {
        "current_mountain_objects_present": len(mountain_objects) > 0,
        "current_expressed_mountain_objects_present": len(expressed) > 0,
        "current_parent_process_coverage_present": parent_process_coverage >= 0.95,
        "current_parent_object_context_present": parent_link_coverage >= 0.95,
        "current_mountain_area_finite_and_capped": (
            0.0 < total_area <= 0.20 and max_area <= 0.08
        ),
        "current_shape_metrics_available": bool(elongations),
        "current_expected_residuals_recorded": expected_residuals_recorded,
        "production_mountain_inventory_fields_available": production_fields_available,
        "production_expected_mountain_kinds_available": (
            production_expected_kinds_available),
    }
    metrics = {
        "continental_landform_object_count": int(len(continental_landforms)),
        "production_mountain_range_object_count": int(len(production_ranges)),
        "mountain_candidate_object_count": int(len(mountain_objects)),
        "expressed_mountain_object_count": int(len(expressed)),
        "parented_mountain_object_count": int(len(parented)),
        "parent_linked_mountain_object_count": int(len(parent_linked)),
        "mountain_kind_count": int(len(kind_counts)),
        "missing_expected_kind_count": int(len(missing_kinds)),
        "missing_inventory_field_count": int(len(missing_fields)),
        "submerged_mountain_candidate_count": int(len(submerged)),
        "tiny_mountain_candidate_count": int(len(tiny)),
        "total_mountain_area_fraction_world": total_area,
        "expressed_mountain_area_fraction_world": float(sum(expressed_area_values)),
        "max_mountain_object_area_fraction_world": max_area,
        "median_mountain_object_area_fraction_world": _percentile(area_values, 50.0),
        "parent_process_coverage_fraction": float(parent_process_coverage),
        "parent_object_context_coverage_fraction": float(parent_link_coverage),
        "median_elongation_ratio": _percentile(elongations, 50.0),
        "max_elongation_ratio": float(max(elongations, default=0.0)),
        "elongated_mountain_object_count": int(sum(
            1 for value in elongations if value >= 3.0)),
        "max_relief_p90_minus_p10_m": float(max(relief_values, default=0.0)),
        "median_relief_p90_minus_p10_m": _percentile(relief_values, 50.0),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "generated_world_mountain_inventory_audit_ready"
            if all(acceptance.values())
            else "generated_world_mountain_inventory_audit_incomplete"
        ),
        "mountain_landform_kinds": MOUNTAIN_LANDFORM_KINDS,
        "mountain_kind_counts": dict(sorted(kind_counts.items())),
        "missing_expected_mountain_kinds": missing_kinds,
        "missing_inventory_fields": missing_fields,
        "expected_current_residual_fields": EXPECTED_CURRENT_RESIDUAL_FIELDS,
        "expected_current_residual_kinds": EXPECTED_CURRENT_RESIDUAL_KINDS,
        "mountain_objects": tuple(object_summaries),
        "limitations": {
            "first_class_mountain_inventory_missing": bool(missing_fields),
            "active_orogen_or_plateau_expression_missing": bool(missing_kinds),
            "elongated_range_expression_underdeveloped": (
                (
                    metrics["max_elongation_ratio"] < 1.65
                )
                if production_ranges else (
                    metrics["elongated_mountain_object_count"] < 2
                    or metrics["median_elongation_ratio"] < 3.0
                )
            ),
            "submerged_arc_or_microcontinent_candidates_present": bool(submerged),
            "tiny_mountain_components_present": bool(tiny),
        },
        "metrics": metrics,
        "acceptance": acceptance,
    }


def mountain_inventory_expression_summary(world: Any | None = None) -> dict[str, Any]:
    reference = reference_mountain_inventory_summary()
    current = (
        current_generated_mountain_inventory_audit(world)
        if world is not None
        else {
            "status": "generated_world_mountain_inventory_audit_not_run",
            "acceptance": {
                "current_mountain_objects_present": False,
                "current_expressed_mountain_objects_present": False,
                "current_parent_process_coverage_present": False,
                "current_parent_object_context_present": False,
                "current_mountain_area_finite_and_capped": False,
                "current_shape_metrics_available": False,
                "current_expected_residuals_recorded": False,
                "production_mountain_inventory_fields_available": False,
                "production_expected_mountain_kinds_available": False,
            },
            "metrics": {},
        }
    )
    acceptance = {
        "reference_inventory_ready": (
            reference["status"] == "mountain_inventory_reference_ready"),
        "reference_hierarchy_and_parentage_ready": (
            reference["acceptance"]["hierarchy_parent_links_valid"]
            and reference["acceptance"]["parent_processes_present"]
            and reference["acceptance"]["object_backed_not_threshold_only"]
        ),
        "current_generated_audit_available": (
            current["status"] == "generated_world_mountain_inventory_audit_ready"),
        "current_object_backed_mountains_present": current["acceptance"][
            "current_mountain_objects_present"],
        "current_parentage_context_present": (
            current["acceptance"]["current_parent_process_coverage_present"]
            and current["acceptance"]["current_parent_object_context_present"]
        ),
        "current_area_and_shape_metrics_available": (
            current["acceptance"]["current_mountain_area_finite_and_capped"]
            and current["acceptance"]["current_shape_metrics_available"]
        ),
        "current_expected_residuals_recorded": current["acceptance"][
            "current_expected_residuals_recorded"],
        "production_mountain_inventory_fields_available": current["acceptance"][
            "production_mountain_inventory_fields_available"],
        "production_expected_mountain_kinds_available": current["acceptance"][
            "production_expected_mountain_kinds_available"],
    }
    return {
        "schema": SCHEMA,
        "status": (
            "mountain_inventory_expression_ready"
            if all(acceptance.values())
            else "mountain_inventory_expression_incomplete"
        ),
        "reference": reference,
        "current_generated": current,
        "acceptance": acceptance,
        "summary_digest": _canonical_digest({
            "reference": reference,
            "current": current,
        }),
    }
