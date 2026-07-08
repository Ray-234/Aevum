"""Rift-to-passive-margin escarpment sequence diagnostics."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from typing import Any

import numpy as np

from aevum.modules.terrain import (
    OCEAN_DEPTH_ABYSS,
    OCEAN_DEPTH_RISE,
    OCEAN_DEPTH_SHELF,
    OCEAN_DEPTH_SLOPE,
    OCEAN_MARGIN_PASSIVE,
    RIFT_MARGIN_STAGE_ABYSS,
    RIFT_MARGIN_STAGE_ESCARPMENT,
    RIFT_MARGIN_STAGE_PASSIVE_LOWLAND,
    RIFT_MARGIN_STAGE_RIFT_BASIN,
    RIFT_MARGIN_STAGE_RISE,
    RIFT_MARGIN_STAGE_SHELF,
    RIFT_MARGIN_STAGE_SHOULDER,
    RIFT_MARGIN_STAGE_SLOPE,
)


SCHEMA = "aevum.rift_margin_escarpment_sequence.v1"

SOURCE_IDS = (
    "NOAA_ETOPO_2022",
    "GEBCO_GRIDDED_BATHYMETRY",
    "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
    "NOAA_TOTAL_SEDIMENT_THICKNESS",
)

EXPECTED_STAGE_SEQUENCE = (
    "stable_platform",
    "rift_shoulder",
    "rift_basin",
    "rift_axis",
    "opposite_rift_shoulder",
    "passive_margin_escarpment",
    "passive_margin_lowland",
    "continental_shelf",
    "continental_slope",
    "continental_rise",
    "abyssal_plain",
)

EXPECTED_CURRENT_RESIDUAL_ITEMS = (
    "terrain.rift_shoulders",
    "terrain.escarpments",
    "terrain.rift_margin_sequence_id",
    "terrain.rift_margin_stage",
    "tectonics.rift_margin_lineage_id",
)

REFERENCE_ZONES: tuple[dict[str, Any], ...] = (
    {
        "id": "platform",
        "stage": "stable_platform",
        "zone_class": "platform",
        "parent_process": "stable_continental_interior",
        "elevation_m": 620.0,
        "water_depth_m": 0.0,
        "sediment_m": 520.0,
        "relief_m": 180.0,
        "land": True,
    },
    {
        "id": "rift_shoulder_west",
        "stage": "rift_shoulder",
        "zone_class": "rift_shoulder",
        "parent_process": "rift_shoulder_uplift",
        "elevation_m": 1080.0,
        "water_depth_m": 0.0,
        "sediment_m": 280.0,
        "relief_m": 520.0,
        "land": True,
    },
    {
        "id": "rift_basin",
        "stage": "rift_basin",
        "zone_class": "rift_basin",
        "parent_process": "rift_basin_subsidence",
        "elevation_m": -80.0,
        "water_depth_m": 0.0,
        "sediment_m": 2800.0,
        "relief_m": 260.0,
        "land": True,
    },
    {
        "id": "rift_axis",
        "stage": "rift_axis",
        "zone_class": "rift_axis",
        "parent_process": "continental_extension",
        "elevation_m": 40.0,
        "water_depth_m": 0.0,
        "sediment_m": 2400.0,
        "relief_m": 320.0,
        "land": True,
    },
    {
        "id": "rift_shoulder_east",
        "stage": "opposite_rift_shoulder",
        "zone_class": "rift_shoulder",
        "parent_process": "rift_shoulder_uplift",
        "elevation_m": 980.0,
        "water_depth_m": 0.0,
        "sediment_m": 340.0,
        "relief_m": 480.0,
        "land": True,
    },
    {
        "id": "escarpment",
        "stage": "passive_margin_escarpment",
        "zone_class": "escarpment",
        "parent_process": "passive_margin_uplift_and_backwearing",
        "elevation_m": 760.0,
        "water_depth_m": 0.0,
        "sediment_m": 460.0,
        "relief_m": 780.0,
        "land": True,
    },
    {
        "id": "coastal_lowland",
        "stage": "passive_margin_lowland",
        "zone_class": "passive_margin_lowland",
        "parent_process": "passive_margin_subsidence",
        "elevation_m": 55.0,
        "water_depth_m": 0.0,
        "sediment_m": 1850.0,
        "relief_m": 150.0,
        "land": True,
    },
    {
        "id": "shelf",
        "stage": "continental_shelf",
        "zone_class": "continental_shelf",
        "parent_process": "shelf_sedimentation",
        "elevation_m": -120.0,
        "water_depth_m": 120.0,
        "sediment_m": 3200.0,
        "relief_m": 80.0,
        "land": False,
    },
    {
        "id": "slope",
        "stage": "continental_slope",
        "zone_class": "continental_slope",
        "parent_process": "margin_slope_transition",
        "elevation_m": -1700.0,
        "water_depth_m": 1700.0,
        "sediment_m": 1350.0,
        "relief_m": 900.0,
        "land": False,
    },
    {
        "id": "rise",
        "stage": "continental_rise",
        "zone_class": "continental_rise",
        "parent_process": "slope_apron_sedimentation",
        "elevation_m": -3100.0,
        "water_depth_m": 3100.0,
        "sediment_m": 820.0,
        "relief_m": 360.0,
        "land": False,
    },
    {
        "id": "abyss",
        "stage": "abyssal_plain",
        "zone_class": "abyssal_plain",
        "parent_process": "old_oceanic_crust_sediment_smoothing",
        "elevation_m": -4300.0,
        "water_depth_m": 4300.0,
        "sediment_m": 460.0,
        "relief_m": 220.0,
        "land": False,
    },
)

REFERENCE_EDGES = tuple(
    (REFERENCE_ZONES[idx]["zone_class"], REFERENCE_ZONES[idx + 1]["zone_class"])
    for idx in range(len(REFERENCE_ZONES) - 1)
)


def _canonical_digest(data: Any) -> str:
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _field_or_empty(world: Any, name: str, fill: float = 0.0) -> np.ndarray:
    if name in getattr(world, "fields", {}):
        return np.asarray(world.field(name), dtype=np.float64)
    return np.full(world.grid.n, fill, dtype=np.float64)


def _item_present(world: Any, name: str) -> bool:
    return (
        name in getattr(world, "fields", {})
        or name in getattr(world, "objects", {})
        or name in getattr(world, "networks", {})
        or name in getattr(world, "globals", {})
    )


def _object_mask(world: Any, object_set: str, kind: str | None = None) -> np.ndarray:
    mask = np.zeros(world.grid.n, dtype=bool)
    for obj in getattr(world, "objects", {}).get(object_set, []):
        if kind is not None and str(obj.get("kind", "")) != kind:
            continue
        cells = np.asarray(obj.get("cells", ()), dtype=np.int64)
        cells = cells[(cells >= 0) & (cells < world.grid.n)]
        if cells.size:
            mask[cells] = True
    return mask


def _dilate(grid: Any, mask: np.ndarray, passes: int) -> np.ndarray:
    out = np.asarray(mask, dtype=bool).copy()
    frontier = out.copy()
    for _ in range(int(passes)):
        if not frontier.any():
            break
        nxt = np.zeros(grid.n, dtype=bool)
        for cell in np.where(frontier)[0]:
            nxt[grid.neighbors[int(cell)]] = True
        nxt &= ~out
        out |= nxt
        frontier = nxt
    return out


def _area_fraction(area: np.ndarray, mask: np.ndarray, within: np.ndarray | None = None) -> float:
    mask = np.asarray(mask, dtype=bool)
    if within is not None:
        within = np.asarray(within, dtype=bool)
        denom = float(area[within].sum())
        mask = mask & within
    else:
        denom = float(area.sum())
    if denom <= 0.0:
        return 0.0
    return float(area[mask].sum() / denom)


def _near_fraction(
    grid: Any,
    area: np.ndarray,
    source: np.ndarray,
    target: np.ndarray,
    passes: int,
) -> float:
    source = np.asarray(source, dtype=bool)
    if not source.any():
        return 0.0
    target_near = _dilate(grid, np.asarray(target, dtype=bool), passes=passes)
    return float(area[source & target_near].sum() / max(float(area[source].sum()), 1.0e-12))


def _percentile(values: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    return float(np.percentile(values, q))


def reference_rift_margin_escarpment_sequence_summary() -> dict[str, Any]:
    zones = tuple(dict(zone) for zone in REFERENCE_ZONES)
    stage_sequence = tuple(str(zone["stage"]) for zone in zones)
    classes = tuple(str(zone["zone_class"]) for zone in zones)
    edge_set = set(tuple(edge) for edge in REFERENCE_EDGES)
    required_edges = {
        ("platform", "rift_shoulder"),
        ("rift_shoulder", "rift_basin"),
        ("rift_basin", "rift_axis"),
        ("rift_shoulder", "escarpment"),
        ("escarpment", "passive_margin_lowland"),
        ("passive_margin_lowland", "continental_shelf"),
        ("continental_shelf", "continental_slope"),
        ("continental_slope", "continental_rise"),
        ("continental_rise", "abyssal_plain"),
    }
    by_class = {str(zone["zone_class"]): zone for zone in zones}
    missing_edges = tuple(sorted(required_edges - edge_set))
    parent_process_failures = tuple(
        str(zone["id"]) for zone in zones if not zone.get("parent_process")
    )
    bathymetry_ordered = (
        by_class["continental_shelf"]["water_depth_m"]
        < by_class["continental_slope"]["water_depth_m"]
        < by_class["continental_rise"]["water_depth_m"]
        < by_class["abyssal_plain"]["water_depth_m"]
    )
    lowland_shelf_coupled = (
        by_class["passive_margin_lowland"]["elevation_m"]
        > by_class["continental_shelf"]["elevation_m"]
        and by_class["continental_shelf"]["water_depth_m"] <= 250.0
    )
    rift_relief_ordered = (
        by_class["rift_shoulder"]["elevation_m"]
        > by_class["platform"]["elevation_m"]
        > by_class["rift_basin"]["elevation_m"]
        and by_class["rift_shoulder"]["elevation_m"]
        > by_class["passive_margin_lowland"]["elevation_m"]
    )
    escarpment_expressed = (
        by_class["escarpment"]["relief_m"] >= 600.0
        and by_class["escarpment"]["elevation_m"]
        > by_class["passive_margin_lowland"]["elevation_m"] + 500.0
    )
    sediment_ordered = (
        by_class["continental_shelf"]["sediment_m"]
        > by_class["passive_margin_lowland"]["sediment_m"]
        > by_class["platform"]["sediment_m"]
        and by_class["rift_basin"]["sediment_m"]
        > by_class["rift_shoulder"]["sediment_m"]
    )
    acceptance = {
        "fixture_schema_ready": bool(zones and SOURCE_IDS),
        "stage_sequence_complete": stage_sequence == EXPECTED_STAGE_SEQUENCE,
        "required_classes_present": set(classes) >= {
            "platform", "rift_shoulder", "rift_basin", "rift_axis",
            "escarpment", "passive_margin_lowland", "continental_shelf",
            "continental_slope", "continental_rise", "abyssal_plain",
        },
        "required_adjacency_edges_present": not missing_edges,
        "rift_relief_ordering": rift_relief_ordered,
        "escarpment_expression_present": escarpment_expressed,
        "lowland_shelf_coupled": lowland_shelf_coupled,
        "shelf_slope_rise_abyss_ordered": bathymetry_ordered,
        "sediment_ordering_plausible": sediment_ordered,
        "parent_processes_present": not parent_process_failures,
        "raw_reference_data_not_stored": True,
    }
    metrics = {
        "zone_count": int(len(zones)),
        "class_count": int(len(set(classes))),
        "edge_count": int(len(edge_set)),
        "missing_required_edge_count": int(len(missing_edges)),
        "parent_process_failure_count": int(len(parent_process_failures)),
        "rift_shoulder_elevation_m": float(by_class["rift_shoulder"]["elevation_m"]),
        "rift_basin_elevation_m": float(by_class["rift_basin"]["elevation_m"]),
        "escarpment_relief_m": float(by_class["escarpment"]["relief_m"]),
        "passive_margin_lowland_elevation_m": float(
            by_class["passive_margin_lowland"]["elevation_m"]),
        "shelf_depth_m": float(by_class["continental_shelf"]["water_depth_m"]),
        "slope_depth_m": float(by_class["continental_slope"]["water_depth_m"]),
        "rise_depth_m": float(by_class["continental_rise"]["water_depth_m"]),
        "abyss_depth_m": float(by_class["abyssal_plain"]["water_depth_m"]),
        "shelf_sediment_m": float(by_class["continental_shelf"]["sediment_m"]),
        "passive_margin_lowland_sediment_m": float(
            by_class["passive_margin_lowland"]["sediment_m"]),
        "rift_basin_sediment_m": float(by_class["rift_basin"]["sediment_m"]),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "rift_margin_escarpment_reference_ready"
            if all(acceptance.values())
            else "rift_margin_escarpment_reference_incomplete"
        ),
        "source_ids": SOURCE_IDS,
        "zones": zones,
        "stage_sequence": stage_sequence,
        "expected_stage_sequence": EXPECTED_STAGE_SEQUENCE,
        "adjacency_edges": REFERENCE_EDGES,
        "missing_required_edges": missing_edges,
        "parent_process_failures": parent_process_failures,
        "extraction_policy": {
            "raw_topography_bathymetry_or_geologic_vectors_stored": False,
            "direct_rift_margin_extraction_pending": True,
            "fixture_is_small_derived_reference": True,
        },
        "metrics": metrics,
        "acceptance": acceptance,
        "fixture_digest": _canonical_digest({
            "zones": zones,
            "source_ids": SOURCE_IDS,
        }),
    }


def current_generated_rift_margin_audit(world: Any) -> dict[str, Any]:
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    fields = getattr(world, "fields", {})
    objects = getattr(world, "objects", {})
    elevation = _field_or_empty(world, "terrain.elevation_m")
    sea = float(getattr(world, "sea_level", world.g("ocean.sea_level_m", 0.0)))
    depth = np.maximum(sea - elevation, 0.0)
    ocean = _field_or_empty(world, "ocean.mask").astype(bool)
    depth_province = _field_or_empty(world, "ocean.depth_province").astype(int)
    margin_type = _field_or_empty(world, "ocean.margin_type").astype(int)
    shelf_width = _field_or_empty(world, "ocean.shelf_width")
    passive_lowland = _field_or_empty(world, "terrain.passive_margin_lowland").astype(bool)
    landforms = list(objects.get("terrain.continental_landforms", []))
    margin_landforms = list(objects.get("terrain.margin_landforms", []))
    landform_counts = Counter(str(obj.get("kind", "")) for obj in landforms)
    margin_counts = Counter(str(obj.get("kind", "")) for obj in margin_landforms)
    rift_basin = _object_mask(world, "terrain.continental_landforms", "rift_basin")
    passive_lowland_objects = _object_mask(
        world, "terrain.continental_landforms", "passive_margin_lowland")
    passive_lowland_all = passive_lowland | passive_lowland_objects
    passive_wedge = _object_mask(world, "terrain.margin_landforms", "passive_margin_wedge")
    production_shoulder = _field_or_empty(world, "terrain.rift_shoulders")
    production_escarpment = _field_or_empty(world, "terrain.escarpments")
    production_sequence = _field_or_empty(world, "terrain.rift_margin_sequence_id")
    production_stage = _field_or_empty(world, "terrain.rift_margin_stage").astype(int)
    production_lineage = _field_or_empty(world, "tectonics.rift_margin_lineage_id")
    shelf = ocean & (depth_province == OCEAN_DEPTH_SHELF)
    slope = ocean & (depth_province == OCEAN_DEPTH_SLOPE)
    rise = ocean & (depth_province == OCEAN_DEPTH_RISE)
    abyss = ocean & (depth_province == OCEAN_DEPTH_ABYSS)
    passive_margin_ocean = ocean & (margin_type == OCEAN_MARGIN_PASSIVE)
    missing_items = tuple(
        item for item in EXPECTED_CURRENT_RESIDUAL_ITEMS if not _item_present(world, item)
    )
    expected_residuals_recorded = set(missing_items).issubset(
        set(EXPECTED_CURRENT_RESIDUAL_ITEMS))
    production_fields_available = not missing_items

    shoulder_cells = production_shoulder > 0.0
    escarpment_cells = production_escarpment > 0.0
    sequence_cells = production_sequence > 0.0
    lineage_cells = production_lineage > 0.0
    stage_masks = {
        "rift_shoulder": production_stage == RIFT_MARGIN_STAGE_SHOULDER,
        "rift_basin": production_stage == RIFT_MARGIN_STAGE_RIFT_BASIN,
        "passive_margin_escarpment": (
            production_stage == RIFT_MARGIN_STAGE_ESCARPMENT),
        "passive_margin_lowland": (
            production_stage == RIFT_MARGIN_STAGE_PASSIVE_LOWLAND),
        "continental_shelf": production_stage == RIFT_MARGIN_STAGE_SHELF,
        "continental_slope": production_stage == RIFT_MARGIN_STAGE_SLOPE,
        "continental_rise": production_stage == RIFT_MARGIN_STAGE_RISE,
        "abyssal_plain": production_stage == RIFT_MARGIN_STAGE_ABYSS,
    }
    passive_lowland_all |= stage_masks["passive_margin_lowland"]
    production_stage_count = int(sum(mask.any() for mask in stage_masks.values()))
    production_sequence_ids = {
        int(value) for value in np.unique(production_sequence) if int(value) > 0
    }
    production_lineage_ids = {
        int(value) for value in np.unique(production_lineage) if int(value) > 0
    }
    production_sequence_objects = list(
        objects.get("terrain.rift_margin_sequences", []))
    production_shoulder_objects = list(objects.get("terrain.rift_shoulders", []))
    production_escarpment_objects = list(objects.get("terrain.escarpments", []))

    production_depths = {
        "shelf_depth_p75_m": _percentile(
            depth[stage_masks["continental_shelf"]], 75.0),
        "slope_depth_p50_m": _percentile(
            depth[stage_masks["continental_slope"]], 50.0),
        "rise_depth_p50_m": _percentile(
            depth[stage_masks["continental_rise"]], 50.0),
        "abyss_depth_p50_m": _percentile(
            depth[stage_masks["abyssal_plain"]], 50.0),
    }
    production_bathymetry_ordered = (
        np.isfinite(production_depths["shelf_depth_p75_m"])
        and np.isfinite(production_depths["slope_depth_p50_m"])
        and np.isfinite(production_depths["rise_depth_p50_m"])
        and np.isfinite(production_depths["abyss_depth_p50_m"])
        and production_depths["shelf_depth_p75_m"]
        < production_depths["slope_depth_p50_m"]
        < production_depths["rise_depth_p50_m"]
        < production_depths["abyss_depth_p50_m"]
    )

    lowland_near_shelf_p1 = _near_fraction(
        grid, area, passive_lowland_all, shelf, passes=1)
    lowland_near_shelf_p2 = _near_fraction(
        grid, area, passive_lowland_all, shelf, passes=2)
    lowland_near_wedge_p2 = _near_fraction(
        grid, area, passive_lowland_all, passive_wedge, passes=2)
    rift_near_passive_p3 = _near_fraction(
        grid, area, rift_basin, passive_margin_ocean, passes=3)
    rift_near_passive_p5 = _near_fraction(
        grid, area, rift_basin, passive_margin_ocean, passes=5)
    wedge_near_lowland_p2 = _near_fraction(
        grid, area, passive_wedge, passive_lowland_all, passes=2)

    passive_lowland_object_area = sum(
        float(obj.get("area_fraction", 0.0))
        for obj in landforms
        if str(obj.get("kind", "")) == "passive_margin_lowland"
    )
    rift_object_area = sum(
        float(obj.get("area_fraction", 0.0))
        for obj in landforms
        if str(obj.get("kind", "")) == "rift_basin"
    )
    wedge_object_area = sum(
        float(obj.get("area_fraction", 0.0))
        for obj in margin_landforms
        if str(obj.get("kind", "")) == "passive_margin_wedge"
    )
    lowland_parented = [
        obj for obj in landforms
        if str(obj.get("kind", "")) == "passive_margin_lowland"
        and obj.get("parent_process")
    ]
    rift_parented = [
        obj for obj in landforms
        if str(obj.get("kind", "")) == "rift_basin" and obj.get("parent_process")
    ]
    wedge_parented = [
        obj for obj in margin_landforms
        if str(obj.get("kind", "")) == "passive_margin_wedge"
        and obj.get("parent_process")
    ]
    metrics = {
        "continental_landform_object_count": int(len(landforms)),
        "margin_landform_object_count": int(len(margin_landforms)),
        "rift_basin_object_count": int(landform_counts.get("rift_basin", 0)),
        "passive_margin_lowland_object_count": int(
            landform_counts.get("passive_margin_lowland", 0)),
        "passive_margin_wedge_object_count": int(
            margin_counts.get("passive_margin_wedge", 0)),
        "delta_fan_object_count": int(margin_counts.get("delta_fan", 0)),
        "missing_sequence_item_count": int(len(missing_items)),
        "rift_basin_area_fraction_world": float(rift_object_area),
        "passive_margin_lowland_object_area_fraction_world": float(
            passive_lowland_object_area),
        "passive_margin_lowland_field_area_fraction_world": _area_fraction(
            area, passive_lowland),
        "passive_margin_wedge_area_fraction_world": float(wedge_object_area),
        "shelf_fraction_of_ocean": _area_fraction(area, shelf, ocean),
        "slope_fraction_of_ocean": _area_fraction(area, slope, ocean),
        "rise_fraction_of_ocean": _area_fraction(area, rise, ocean),
        "abyss_fraction_of_ocean": _area_fraction(area, abyss, ocean),
        "passive_margin_ocean_fraction": _area_fraction(
            area, passive_margin_ocean, ocean),
        "lowland_near_shelf_fraction_p1": float(lowland_near_shelf_p1),
        "lowland_near_shelf_fraction_p2": float(lowland_near_shelf_p2),
        "lowland_near_wedge_fraction_p2": float(lowland_near_wedge_p2),
        "rift_near_passive_margin_fraction_p3": float(rift_near_passive_p3),
        "rift_near_passive_margin_fraction_p5": float(rift_near_passive_p5),
        "wedge_near_lowland_fraction_p2": float(wedge_near_lowland_p2),
        "shelf_depth_p75_m": _percentile(depth[shelf], 75.0),
        "slope_depth_p50_m": _percentile(depth[slope], 50.0),
        "rise_depth_p50_m": _percentile(depth[rise], 50.0),
        "abyss_depth_p50_m": _percentile(depth[abyss], 50.0),
        "mean_shelf_width_steps": float(np.mean(shelf_width[shelf])) if shelf.any() else 0.0,
        "parented_rift_basin_object_count": int(len(rift_parented)),
        "parented_passive_margin_lowland_object_count": int(len(lowland_parented)),
        "parented_passive_margin_wedge_object_count": int(len(wedge_parented)),
        "production_sequence_field_cell_count": int(np.count_nonzero(sequence_cells)),
        "production_sequence_id_count": int(len(production_sequence_ids)),
        "production_lineage_id_count": int(len(production_lineage_ids)),
        "production_stage_count": int(production_stage_count),
        "production_rift_shoulder_cell_count": int(np.count_nonzero(shoulder_cells)),
        "production_escarpment_cell_count": int(np.count_nonzero(escarpment_cells)),
        "production_sequence_object_count": int(len(production_sequence_objects)),
        "production_rift_shoulder_object_count": int(len(production_shoulder_objects)),
        "production_escarpment_object_count": int(len(production_escarpment_objects)),
        "production_passive_lowland_area_fraction_world": _area_fraction(
            area, stage_masks["passive_margin_lowland"]),
        "production_shelf_cell_count": int(np.count_nonzero(
            stage_masks["continental_shelf"])),
        "production_slope_cell_count": int(np.count_nonzero(
            stage_masks["continental_slope"])),
        "production_rise_cell_count": int(np.count_nonzero(
            stage_masks["continental_rise"])),
        "production_abyss_cell_count": int(np.count_nonzero(
            stage_masks["abyssal_plain"])),
        "production_shelf_depth_p75_m": production_depths["shelf_depth_p75_m"],
        "production_slope_depth_p50_m": production_depths["slope_depth_p50_m"],
        "production_rise_depth_p50_m": production_depths["rise_depth_p50_m"],
        "production_abyss_depth_p50_m": production_depths["abyss_depth_p50_m"],
    }
    lowland_shelf_coupled = (
        (
            metrics["passive_margin_lowland_object_count"] > 0
            or metrics["production_passive_lowland_area_fraction_world"] > 0.0
        )
        and metrics["passive_margin_wedge_object_count"] > 0
        and metrics["lowland_near_shelf_fraction_p2"] >= 0.75
        and metrics["lowland_near_wedge_fraction_p2"] >= 0.50
    )
    shelf_profile_available = (
        np.isfinite(metrics["shelf_depth_p75_m"])
        and np.isfinite(metrics["abyss_depth_p50_m"])
        and metrics["shelf_depth_p75_m"] < 800.0
        and metrics["abyss_depth_p50_m"] > metrics["shelf_depth_p75_m"] + 1000.0
    )
    rift_passive_context = (
        metrics["rift_basin_object_count"] > 0
        and metrics["rift_near_passive_margin_fraction_p5"] >= 0.50
    )
    parentage_present = (
        metrics["parented_rift_basin_object_count"] > 0
        and (
            metrics["parented_passive_margin_lowland_object_count"] > 0
            or metrics["production_sequence_object_count"] > 0
        )
        and metrics["parented_passive_margin_wedge_object_count"] > 0
    )
    production_sequence_ready = (
        production_fields_available
        and metrics["production_sequence_field_cell_count"] > 0
        and metrics["production_sequence_id_count"] > 0
        and metrics["production_lineage_id_count"] > 0
        and metrics["production_stage_count"] >= 6
        and metrics["production_rift_shoulder_cell_count"] > 0
        and metrics["production_escarpment_cell_count"] > 0
        and metrics["production_sequence_object_count"] > 0
        and metrics["production_rift_shoulder_object_count"] > 0
        and metrics["production_escarpment_object_count"] > 0
        and metrics["production_passive_lowland_area_fraction_world"] > 0.0
        and metrics["production_shelf_cell_count"] > 0
        and metrics["production_slope_cell_count"] > 0
        and metrics["production_rise_cell_count"] > 0
        and metrics["production_abyss_cell_count"] > 0
        and production_bathymetry_ordered
    )
    acceptance = {
        "current_core_fields_available": all(
            item in fields
            for item in (
                "terrain.elevation_m",
                "terrain.passive_margin_lowland",
                "ocean.depth_province",
                "ocean.margin_type",
                "ocean.shelf_width",
            )
        ),
        "current_required_objects_present": (
            metrics["rift_basin_object_count"] > 0
            and (
                metrics["passive_margin_lowland_object_count"] > 0
                or metrics["production_passive_lowland_area_fraction_world"] > 0.0
            )
            and metrics["passive_margin_wedge_object_count"] > 0
        ),
        "current_lowland_shelf_coupled": lowland_shelf_coupled,
        "current_shelf_profile_not_deep_nearshore": shelf_profile_available,
        "current_rift_passive_context_present": rift_passive_context,
        "current_parentage_present": parentage_present,
        "current_expected_residuals_recorded": expected_residuals_recorded,
        "production_rift_margin_sequence_fields_available": production_fields_available,
        "production_rift_margin_sequence_ready": production_sequence_ready,
        "production_shelf_slope_rise_abyss_ordered": production_bathymetry_ordered,
    }
    return {
        "schema": SCHEMA,
        "status": (
            "generated_world_rift_margin_sequence_audit_ready"
            if all(acceptance.values())
            else "generated_world_rift_margin_sequence_audit_incomplete"
        ),
        "landform_kind_counts": dict(sorted(landform_counts.items())),
        "margin_landform_kind_counts": dict(sorted(margin_counts.items())),
        "missing_sequence_items": missing_items,
        "expected_current_residual_items": EXPECTED_CURRENT_RESIDUAL_ITEMS,
        "limitations": {
            "first_class_rift_margin_sequence_missing": bool(missing_items),
            "rift_shoulder_objects_missing": (
                "terrain.rift_shoulders" in missing_items
                or metrics["production_rift_shoulder_object_count"] == 0),
            "escarpment_objects_missing": (
                "terrain.escarpments" in missing_items
                or metrics["production_escarpment_object_count"] == 0),
            "passive_margin_lowland_objects_tiny": (
                metrics["passive_margin_lowland_object_area_fraction_world"] < 0.01
                and metrics["production_passive_lowland_area_fraction_world"] <= 0.0
            ),
            "rift_to_margin_lineage_missing": (
                "tectonics.rift_margin_lineage_id" in missing_items
                or metrics["production_lineage_id_count"] == 0),
        },
        "metrics": metrics,
        "acceptance": acceptance,
    }


def rift_margin_escarpment_sequence_summary(world: Any | None = None) -> dict[str, Any]:
    reference = reference_rift_margin_escarpment_sequence_summary()
    current = (
        current_generated_rift_margin_audit(world)
        if world is not None
        else {
            "status": "generated_world_rift_margin_sequence_audit_not_run",
            "acceptance": {
                "current_core_fields_available": False,
                "current_required_objects_present": False,
                "current_lowland_shelf_coupled": False,
                "current_shelf_profile_not_deep_nearshore": False,
                "current_rift_passive_context_present": False,
                "current_parentage_present": False,
                "current_expected_residuals_recorded": False,
                "production_rift_margin_sequence_fields_available": False,
                "production_rift_margin_sequence_ready": False,
                "production_shelf_slope_rise_abyss_ordered": False,
            },
            "metrics": {},
        }
    )
    acceptance = {
        "reference_sequence_ready": (
            reference["status"] == "rift_margin_escarpment_reference_ready"),
        "reference_ordering_ready": (
            reference["acceptance"]["rift_relief_ordering"]
            and reference["acceptance"]["escarpment_expression_present"]
            and reference["acceptance"]["shelf_slope_rise_abyss_ordered"]
        ),
        "current_generated_audit_available": (
            current["status"] == "generated_world_rift_margin_sequence_audit_ready"),
        "current_lowland_shelf_coupled": current["acceptance"][
            "current_lowland_shelf_coupled"],
        "current_rift_passive_context_present": current["acceptance"][
            "current_rift_passive_context_present"],
        "current_expected_residuals_recorded": current["acceptance"][
            "current_expected_residuals_recorded"],
        "production_rift_margin_sequence_fields_available": current["acceptance"][
            "production_rift_margin_sequence_fields_available"],
        "production_rift_margin_sequence_ready": current["acceptance"][
            "production_rift_margin_sequence_ready"],
        "production_shelf_slope_rise_abyss_ordered": current["acceptance"][
            "production_shelf_slope_rise_abyss_ordered"],
    }
    return {
        "schema": SCHEMA,
        "status": (
            "rift_margin_escarpment_sequence_ready"
            if all(acceptance.values())
            else "rift_margin_escarpment_sequence_incomplete"
        ),
        "reference": reference,
        "current_generated": current,
        "acceptance": acceptance,
        "summary_digest": _canonical_digest({
            "reference": reference,
            "current": current,
        }),
    }
