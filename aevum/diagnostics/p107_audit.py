"""P107 high-resolution plate and geomorphology audit.

P107 is a diagnostics-first phase: before changing plate identity, boundary
province logic, or terrain response, the terminal world must preserve enough raw
state to measure what changed.  This module writes compact arrays plus stable
JSON metrics so 8000- and 24000-cell runs can be compared without reading PNG
colors.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
import time
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aevum.engine import Engine
from aevum.modules.tectonics import (
    DOMAIN_ACCRETED_TERRANE,
    DOMAIN_CONTINENTAL_MARGIN,
    DOMAIN_CRATON,
    DOMAIN_LIP,
    DOMAIN_SUTURE,
    ORIGIN_ARC,
    ORIGIN_PLUME_IMPACT,
)
from aevum.spec.presets import get_preset


SCHEMA = "aevum.p107_terminal_world_audit.v1"
SUMMARY_SCHEMA = "aevum.p107_audit_ladder.v1"
P107_PLATE_TERRAIN_MODULES = (
    "stellar",
    "interior",
    "impacts",
    "tectonics",
    "terrain",
    "climate",
    "biogeochem",
    "biosphere",
)

MAJOR_PLATE_AREA_FRACTION = 0.05
MINOR_PLATE_AREA_FRACTION = 0.008
MICRO_PLATE_AREA_FRACTION = 0.0005

LAND_ELEVATION_BANDS_M = (
    (0.0, 200.0),
    (200.0, 500.0),
    (500.0, 1000.0),
    (1000.0, 2000.0),
    (2000.0, 3000.0),
    (3000.0, 4500.0),
    (4500.0, float("inf")),
)

P109_LAND_HYPSOMETRY_ENVELOPE = {
    "mean_m": (650.0, 1050.0),
    "p50_m": (250.0, 800.0),
    "p90_m": (1700.0, 2700.0),
    "p98_m": (3000.0, 4700.0),
}

P109_LAND_BAND_ENVELOPE = {
    "0_200m": (0.15, 0.35),
    "200_500m": (0.15, 0.38),
    "500_1000m": (0.10, 0.34),
    "1000_2000m": (0.08, 0.30),
    "2000_3000m": (0.01, 0.15),
    "3000_4500m": (0.015, 0.09),
    "gt4500m": (0.003, 0.035),
}

P110A_MAJOR_LAND_COMPONENT_SHARE = 0.03
P110A_MAJOR_OCEAN_BASIN_SHARE = 0.03
P110A_DISCONNECTED_OCEAN_WARNING_SHARE = 0.10
P110A_MODERN_PLANFORM_APPLICABLE_MYR = 3500.0
P110B_HISTORICAL_SUPERCONTINENT_SHARE = 0.70
P110B_HISTORICAL_SUPERCONTINENT_MAX_FRAME_FRACTION = 0.45
P110B_HISTORICAL_SUPERCONTINENT_MAX_DURATION_MYR = 1200.0
P110B_HISTORICAL_SUPERCONTINENT_MAX_WINDOW_MYR = 1800.0
P110B_ARCHETYPE_NAMES = {
    0: "none",
    1: "modern_multipolar",
    2: "afro_eurasia_analog",
    3: "post_supercontinent_breakup",
    4: "archipelago_active_margin",
    5: "supercontinent_final",
}

P107_FIELD_KEYS = (
    "tectonics.plate_id",
    "tectonics.plate_rank",
    "tectonics.protected_plate_id",
    "tectonics.boundary_province_kind",
    "terrain.elevation_m",
    "crust.type",
    "crust.age_myr",
    "ocean.depth_province",
    "terrain.province",
    "ocean.margin_type",
    "ocean.gateway_id",
    "ocean.gateway_system_id",
    "ocean.shelf_width",
    "terrain.continental_detail",
    "terrain.continental_detail_region_code",
    "terrain.internal_geographic_block_region_code",
    "terrain.inland_geomorphology_region_code",
    "terrain.continental_province_code",
    "terrain.mountain_ranges",
    "terrain.mountain_inventory",
    "terrain.mountain_hierarchy_level",
    "terrain.orogenic_parent_hierarchy",
    "terrain.orogenic_hierarchy_spine",
    "terrain.orogenic_shoulder_halo",
    "terrain.orogenic_highland_apron",
    "tectonics.mountain_belt_id",
    "tectonics.mountain_parent_process_id",
    "tectonics.continent_id",
    "tectonics.terrane_id",
    "tectonics.continental_province_code",
    "tectonics.internal_geographic_block_code",
    "tectonics.orogeny_age_myr",
    "tectonics.volcanism_age_myr",
    "ocean.basin_id",
)

P107_BOUNDARY_KEYS = (
    "divergent",
    "ridge",
    "convergent",
    "collision",
    "subduction",
    "trench",
    "suture",
    "active_margin",
    "passive_margin",
    "transform",
)

P107_OBJECT_SETS = (
    "tectonics.boundary_objects",
    "tectonics.boundary_polylines",
    "tectonics.boundary_provinces",
    "tectonics.rift_systems",
    "tectonics.passive_margins",
    "tectonics.spreading_centers",
    "tectonics.closing_margins",
    "tectonics.sutures",
    "tectonics.ocean_basins",
    "tectonics.ocean_gateways",
    "tectonics.deforming_networks",
    "tectonics.cratons",
    "tectonics.shields",
    "tectonics.platforms",
    "tectonics.interior_basins",
    "tectonics.internal_geographic_blocks",
    "tectonics.continents",
    "tectonics.continent_lineage_splits",
    "tectonics.terranes",
    "ocean.basins",
    "ocean.margins",
    "ocean.gateways",
    "ocean.gateway_systems",
    "terrain.ocean_fabric",
    "terrain.p110b_internal_domain_seaways",
    "terrain.p110b_internal_domain_boundaries",
    "terrain.p111619_late_domain_partition_seaways",
    "terrain.margin_landforms",
    "terrain.arc_plume_landforms",
    "terrain.continental_landforms",
    "terrain.mountain_ranges",
    "terrain.plateau_inventory",
    "terrain.orogenic_spine_objects",
    "terrain.orogenic_axis_polylines",
    "terrain.orogenic_spine_repair_candidates",
    "terrain.landform_inventory",
)

P107_OBJECT_MASK_SPECS = (
    ("boundary_ridge", "tectonics.boundary_objects", {"ridge"}),
    ("boundary_transform", "tectonics.boundary_objects", {"transform"}),
    ("boundary_trench", "tectonics.boundary_objects", {"trench"}),
    ("boundary_convergent_parent", "tectonics.boundary_objects", {"convergent_parent"}),
    ("boundary_subduction_parent", "tectonics.boundary_objects", {"subduction_parent"}),
    ("boundary_suture", "tectonics.boundary_objects", {"suture"}),
    ("boundary_active_margin", "tectonics.boundary_objects", {"active_margin"}),
    ("boundary_passive_margin", "tectonics.boundary_objects", {"passive_margin"}),
    ("boundary_polyline_ridge", "tectonics.boundary_polylines", {"ridge_polyline"}),
    ("boundary_polyline_transform", "tectonics.boundary_polylines", {"transform_polyline"}),
    ("boundary_polyline_trench", "tectonics.boundary_polylines", {"trench_polyline"}),
    ("boundary_polyline_suture", "tectonics.boundary_polylines", {"suture_polyline"}),
    ("boundary_polyline_active_margin", "tectonics.boundary_polylines", {"active_margin_polyline"}),
    ("boundary_polyline_passive_margin", "tectonics.boundary_polylines", {"passive_margin_polyline"}),
    ("boundary_polyline_convergent_parent", "tectonics.boundary_polylines", {"convergent_parent_polyline"}),
    ("boundary_polyline_subduction_parent", "tectonics.boundary_polylines", {"subduction_parent_polyline"}),
    ("province_mid_ocean_ridge", "tectonics.boundary_provinces", {"mid_ocean_ridge"}),
    ("province_continental_arc_margin", "tectonics.boundary_provinces", {"continental_arc_margin"}),
    ("province_island_arc_trench", "tectonics.boundary_provinces", {"island_arc_trench"}),
    ("province_collision", "tectonics.boundary_provinces", {"continent_continent_collision"}),
    ("province_passive_margin", "tectonics.boundary_provinces", {"passive_margin"}),
    ("spreading_center", "terrain.ocean_fabric", {"spreading_center"}),
    ("transform_fault", "terrain.ocean_fabric", {"transform_fault"}),
    ("fracture_zone", "terrain.ocean_fabric", {"fracture_zone"}),
    ("abyssal_plain", "terrain.ocean_fabric", {"abyssal_plain"}),
    ("abyssal_hill", "terrain.ocean_fabric", {"abyssal_hill"}),
    ("margin_trench", "terrain.margin_landforms", {"trench"}),
    ("volcanic_arc", "terrain.margin_landforms", {"volcanic_arc"}),
    ("island_arc", "terrain.arc_plume_landforms", {"island_arc"}),
    ("back_arc_basin", "terrain.arc_plume_landforms", {"back_arc_basin"}),
    ("seamount_chain", "terrain.arc_plume_landforms", {"seamount_chain", "hotspot_track"}),
    ("oceanic_plateau", "terrain.arc_plume_landforms", {"oceanic_plateau"}),
    ("microcontinent", "terrain.arc_plume_landforms", {"microcontinent"}),
    ("mountain_orogen", "terrain.mountain_ranges", {"orogen"}),
    ("mountain_old_orogen", "terrain.mountain_ranges", {"old_subdued_orogen"}),
    ("mountain_plateau", "terrain.mountain_ranges", {"plateau"}),
    ("mountain_arc_chain", "terrain.mountain_ranges", {"arc_microcontinent"}),
    ("parent_orogen_crest", "terrain.continental_landforms", {"orogen_crest"}),
    ("parent_orogen_branch", "terrain.continental_landforms", {"orogen_branch_range"}),
    ("parent_orogen_foreland", "terrain.continental_landforms", {"orogen_foreland_slope"}),
    ("parent_orogen_highland_apron", "terrain.continental_landforms", {"orogen_highland_apron"}),
    ("orogenic_main_crest_axis", "terrain.orogenic_axis_polylines", {"orogenic_main_crest_axis"}),
    ("orogenic_branch_axis", "terrain.orogenic_axis_polylines", {"orogenic_branch_axis"}),
)


@dataclass
class P107AuditRun:
    cells: int
    n_plates: int
    seed: int | None = None
    label: str | None = None


@dataclass
class P107AuditConfig:
    preset: str = "earthlike"
    runs: tuple[P107AuditRun, ...] = (
        P107AuditRun(cells=8000, n_plates=36, label="fast_8000_36p"),
        P107AuditRun(cells=24000, n_plates=60, label="main_24000_60p"),
    )
    t_end_myr: float = 4500.0
    frames: int = 5
    render_world_assets: bool = True
    render_contact_sheet: bool = True
    include_earth_reference: bool = True
    enable_ranked_plate_policy: bool = True
    enable_boundary_province_response: bool = True
    enable_p108_boundary_width_guard: bool = True
    enable_p108_high_mountain_coherence: bool = True
    enable_p109_continental_hypsometry_rebalance: bool = True
    enabled_modules: tuple[str, ...] | None = None
    global_overrides: dict[str, float] = field(default_factory=dict)


def plate_rank_metrics(grid, plate: np.ndarray) -> dict[str, Any]:
    """Rank terminal plates by global surface area fraction."""
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    plate = np.asarray(plate, dtype=np.int64)
    valid = plate >= 0
    ids = sorted(int(x) for x in np.unique(plate[valid])) if valid.any() else []

    rows: list[dict[str, Any]] = []
    for pid in ids:
        mask = valid & (plate == pid)
        frac = float(area[mask].sum() / total)
        if frac >= MAJOR_PLATE_AREA_FRACTION:
            rank = "major"
        elif frac >= MINOR_PLATE_AREA_FRACTION:
            rank = "minor"
        elif frac >= MICRO_PLATE_AREA_FRACTION:
            rank = "micro"
        else:
            rank = "tiny"
        rows.append({
            "plate_id": int(pid),
            "area_fraction": frac,
            "cell_count": int(mask.sum()),
            "rank": rank,
        })
    rows.sort(key=lambda row: (-float(row["area_fraction"]), int(row["plate_id"])))
    area_fracs = np.asarray([float(row["area_fraction"]) for row in rows],
                            dtype=np.float64)
    hist_bins = np.asarray(
        [0.0, MICRO_PLATE_AREA_FRACTION, 0.002, MINOR_PLATE_AREA_FRACTION,
         0.02, MAJOR_PLATE_AREA_FRACTION, 0.10, 1.01],
        dtype=np.float64,
    )
    hist_counts = np.histogram(area_fracs, bins=hist_bins)[0].astype(int).tolist()

    def count(rank: str) -> int:
        return int(sum(1 for row in rows if row["rank"] == rank))

    return {
        "plate_rank_thresholds": {
            "major_min_area_fraction": MAJOR_PLATE_AREA_FRACTION,
            "minor_min_area_fraction": MINOR_PLATE_AREA_FRACTION,
            "micro_min_area_fraction": MICRO_PLATE_AREA_FRACTION,
        },
        "terminal_active_plate_count": int(len(rows)),
        "terminal_major_plate_count": count("major"),
        "terminal_minor_plate_count": count("minor"),
        "terminal_microplate_count": count("micro"),
        "terminal_tiny_plate_count": count("tiny"),
        "plate_rank_counts": {
            "major": count("major"),
            "minor": count("minor"),
            "micro": count("micro"),
            "tiny": count("tiny"),
        },
        "plate_area_fraction_histogram": {
            "bins": [float(x) for x in hist_bins.tolist()],
            "counts": hist_counts,
        },
        "plate_area_fraction_top12": [
            float(row["area_fraction"]) for row in rows[:12]
        ],
        "per_plate": rows,
    }


def _p1146_reject_reason_name(code: float) -> str:
    return {
        0: "not_attempted_or_disabled",
        1: "no_spine",
        2: "already_good",
        3: "no_candidate",
        4: "no_multi_spine_support_component",
        5: "no_path_pairs",
        6: "path_search_failed",
        7: "no_bridge_cells",
        8: "branch_area_guard",
        9: "spine_membership_guard",
        10: "linework_score_guard",
        11: "short_component_guard",
        12: "branch_attachment_guard",
        13: "top3_share_guard",
        14: "component_count_guard",
        15: "no_material_improvement",
        16: "accepted",
        17: "terminal_cleanup_proxy_guard",
        18: "class_profile_guard",
        19: "no_class_path_options",
    }.get(int(round(float(code))), "unknown")


def _terrain_internal_profile_metrics(world: Any) -> dict[str, Any]:
    raw = world.objects.get("terrain.internal_profile", {})
    if not isinstance(raw, dict):
        raw = {}
    stage_raw = raw.get("stage_seconds", {})
    if not isinstance(stage_raw, dict):
        stage_raw = {}
    stage_seconds: dict[str, float] = {}
    for key, value in stage_raw.items():
        seconds = float(value)
        if np.isfinite(seconds) and seconds >= 0.0:
            stage_seconds[str(key)] = seconds
    total = float(
        raw.get(
            "total_seconds",
            world.g("terrain.last_internal_profile_total_seconds", 0.0),
        )
    )
    if not np.isfinite(total) or total < 0.0:
        total = 0.0
    profiled_total = max(float(sum(stage_seconds.values())), 1.0e-12)
    top = [
        {
            "stage": key,
            "seconds": float(value),
            "share": float(value / profiled_total),
        }
        for key, value in sorted(
            stage_seconds.items(),
            key=lambda item: item[1],
            reverse=True,
        )[:8]
    ]
    subprofiles: dict[str, Any] = {}
    subprofiles_raw = raw.get("subprofiles", {})
    if isinstance(subprofiles_raw, dict):
        for profile_name, profile_raw in sorted(subprofiles_raw.items()):
            if not isinstance(profile_raw, dict):
                continue
            sub_stage_raw = profile_raw.get("stage_seconds", {})
            if not isinstance(sub_stage_raw, dict):
                continue
            sub_stage_seconds: dict[str, float] = {}
            for key, value in sub_stage_raw.items():
                seconds = float(value)
                if np.isfinite(seconds) and seconds >= 0.0:
                    sub_stage_seconds[str(key)] = seconds
            sub_total = max(float(sum(sub_stage_seconds.values())), 1.0e-12)
            subprofiles[str(profile_name)] = {
                "stage_seconds": sub_stage_seconds,
                "profiled_stage_seconds_total": float(sum(sub_stage_seconds.values())),
                "stage_count": int(len(sub_stage_seconds)),
                "top": [
                    {
                        "stage": key,
                        "seconds": float(value),
                        "share": float(value / sub_total),
                    }
                    for key, value in sorted(
                        sub_stage_seconds.items(),
                        key=lambda item: item[1],
                        reverse=True,
                    )[:8]
                ],
            }
    return {
        "schema": "aevum.terrain_internal_profile.v1",
        "available": bool(stage_seconds),
        "enabled": bool(raw.get("enabled", bool(stage_seconds))),
        "time_myr": float(raw.get("time_myr", world.time_myr)),
        "cell_count": int(raw.get("cell_count", world.grid.n)),
        "total_seconds": float(total),
        "profiled_stage_seconds_total": float(sum(stage_seconds.values())),
        "stage_count": int(len(stage_seconds)),
        "stage_seconds": stage_seconds,
        "top": top,
        "subprofiles": subprofiles,
    }


def terminal_audit_metrics(
    world: Any,
    events: Iterable[Any] = (),
    archive: Any | None = None,
) -> dict[str, Any]:
    """Compute P107.0 terminal-world metrics from a completed world."""
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    sea_level = float(world.sea_level)
    elev = np.asarray(world.get_field("terrain.elevation_m", 0.0), dtype=np.float64)
    rel = elev - sea_level
    land = rel >= 0.0
    ocean = ~land
    plate = np.asarray(world.get_field("tectonics.plate_id", -1.0), dtype=np.int64)
    depth = np.maximum(sea_level - elev, 0.0)
    depth_province = np.asarray(
        world.get_field("ocean.depth_province", -1.0), dtype=np.int64)
    terrain_province = np.asarray(
        world.get_field("terrain.province", -1.0), dtype=np.int64)
    continental_detail = np.asarray(
        world.get_field("terrain.continental_detail", -1.0), dtype=np.int64)
    mountain_inventory = np.asarray(
        world.get_field("terrain.mountain_inventory", -1.0), dtype=np.int64)
    orogenic_parent_hierarchy = np.asarray(
        world.get_field("terrain.orogenic_parent_hierarchy", np.zeros(grid.n)),
        dtype=np.int64,
    )
    if orogenic_parent_hierarchy.shape != (grid.n,):
        orogenic_parent_hierarchy = np.zeros(grid.n, dtype=np.int64)
    orogenic_shoulder_halo = np.asarray(
        world.get_field("terrain.orogenic_shoulder_halo", np.zeros(grid.n)),
        dtype=np.int64,
    )
    if orogenic_shoulder_halo.shape != (grid.n,):
        orogenic_shoulder_halo = np.zeros(grid.n, dtype=np.int64)
    orogenic_highland_apron = np.asarray(
        world.get_field("terrain.orogenic_highland_apron", np.zeros(grid.n)),
        dtype=np.int64,
    )
    if orogenic_highland_apron.shape != (grid.n,):
        orogenic_highland_apron = np.zeros(grid.n, dtype=np.int64)
    orogenic_hierarchy_spine = np.asarray(
        world.get_field("terrain.orogenic_hierarchy_spine", np.zeros(grid.n)),
        dtype=np.int64,
    )
    if orogenic_hierarchy_spine.shape != (grid.n,):
        orogenic_hierarchy_spine = np.zeros(grid.n, dtype=np.int64)
    crust_type = np.asarray(world.get_field("crust.type", -1.0), dtype=np.int64)

    ranks = plate_rank_metrics(grid, plate)
    boundary_counts = _boundary_counts(world)
    object_counts = _object_counts(world)
    object_kind_counts = _object_kind_counts(world)
    protected_plate_ids = _protected_plate_ids(world)
    micro_ids = {
        int(row["plate_id"]) for row in ranks["per_plate"]
        if row["rank"] == "micro"
    }
    merged_microplates = _merged_microplate_count(events)

    ridge_network = _ridge_transform_network_metrics(world)
    boundary_width = _boundary_width_diagnostics(world)
    collision_area = _collision_plateau_area_fraction(world, rel, land)
    high_mountain = _high_mountain_coherence_metrics(world, rel, land)
    relief_planform = _terminal_relief_planform_metrics(
        world,
        rel,
        land,
        orogenic_parent_hierarchy,
        orogenic_hierarchy_spine,
        orogenic_shoulder_halo,
        orogenic_highland_apron,
    )
    orogenic_belt = _orogenic_belt_morphology_metrics(
        world,
        rel,
        land,
        orogenic_parent_hierarchy,
        orogenic_hierarchy_spine,
        orogenic_shoulder_halo,
        orogenic_highland_apron,
    )
    spine_aligned_orogen = _spine_aligned_elevation_morphology_metrics(
        world,
        rel,
        land,
        orogenic_parent_hierarchy,
        orogenic_hierarchy_spine,
        orogenic_shoulder_halo,
        orogenic_highland_apron,
    )
    spine_linework = _orogenic_spine_linework_metrics(
        world,
        land,
        orogenic_parent_hierarchy,
        orogenic_hierarchy_spine,
    )
    polar_edge_orogen = _polar_edge_orogen_overclassification_metrics(
        world,
        land,
        orogenic_parent_hierarchy,
        orogenic_hierarchy_spine,
    )
    arc_trench = _arc_trench_adjacency_metrics(world)
    ocean_feature = _ocean_feature_metrics(world, depth, ocean, depth_province)
    p1116_linework = _p1116_boundary_linework_metrics(world, depth, ocean, depth_province)
    p1116_crust = _p1116_oceanic_crust_structure_metrics(
        world,
        depth,
        ocean,
        depth_province,
    )
    land_stats = _hypsometry_stats(rel, area, land)
    ocean_stats = _hypsometry_stats(-rel, area, ocean)
    land_bands = _hypsometry_band_fractions(rel, area, land)
    earth_hypsometry = _same_grid_earth_hypsometry(grid)
    earth_planform = _same_grid_earth_planform(grid)
    p109_comparison = _p109_hypsometry_comparison(
        land_stats,
        land_bands,
        earth_hypsometry,
    )
    p109_sea_level = _p109_sea_level_consistency_metrics(
        world,
        sea_level,
        float(area[land].sum() / total),
    )
    terrain_internal_profile = _terrain_internal_profile_metrics(world)
    p120_semantic_geometry = _p120_continental_semantic_geometry_metrics(
        world,
        rel,
        land,
        crust_type,
        orogenic_parent_hierarchy,
        orogenic_hierarchy_spine,
        orogenic_shoulder_halo,
        orogenic_highland_apron,
        mountain_inventory,
    )
    p1134_cleanup = {
        "schema": "aevum.p1134_lowres_false_oceanic_island_cleanup.v1",
        "accepted": bool(
            world.g(
                "terrain.last_p1134_lowres_false_oceanic_island_cleanup_accepted",
                0.0,
            )
            >= 1.0
        ),
        "rejected_by_land_floor": bool(
            world.g("terrain.last_p1134_rejected_by_land_floor", 0.0) >= 1.0
        ),
        "candidate_cell_count": float(
            world.g("terrain.last_p1134_candidate_cell_count", 0.0)
        ),
        "candidate_component_count": float(
            world.g("terrain.last_p1134_candidate_component_count", 0.0)
        ),
        "candidate_area_fraction": float(
            world.g("terrain.last_p1134_candidate_area_fraction", 0.0)
        ),
        "adjusted_cell_count": float(
            world.g("terrain.last_p1134_adjusted_cell_count", 0.0)
        ),
        "adjusted_component_count": float(
            world.g("terrain.last_p1134_adjusted_component_count", 0.0)
        ),
        "adjusted_area_fraction": float(
            world.g("terrain.last_p1134_adjusted_area_fraction", 0.0)
        ),
        "component_count_before": float(
            world.g("terrain.last_p1134_component_count_before", 0.0)
        ),
        "component_count_after": float(
            world.g("terrain.last_p1134_component_count_after", 0.0)
        ),
        "land_fraction_before": float(
            world.g("terrain.last_p1134_land_fraction_before", 0.0)
        ),
        "land_fraction_after_cleanup": float(
            world.g("terrain.last_p1134_land_fraction_after", 0.0)
        ),
        "floor_fraction": float(
            world.g("terrain.last_p1134_floor_fraction", 0.0)
        ),
        "terminal_land_fraction": float(area[land].sum() / total),
    }
    p1134b_neck_widening = {
        "schema": "aevum.p1134b_lowres_articulation_neck_widening.v1",
        "accepted": bool(
            world.g(
                "terrain.last_p1134b_lowres_articulation_neck_widening_accepted",
                0.0,
            )
            >= 1.0
        ),
        "rejected_by_land_ceiling": bool(
            world.g("terrain.last_p1134b_rejected_by_land_ceiling", 0.0) >= 1.0
        ),
        "candidate_cell_count": float(
            world.g("terrain.last_p1134b_candidate_cell_count", 0.0)
        ),
        "adjusted_cell_count": float(
            world.g("terrain.last_p1134b_adjusted_cell_count", 0.0)
        ),
        "adjusted_area_fraction": float(
            world.g("terrain.last_p1134b_adjusted_area_fraction", 0.0)
        ),
        "component_count_before": float(
            world.g("terrain.last_p1134b_component_count_before", 0.0)
        ),
        "component_count_after": float(
            world.g("terrain.last_p1134b_component_count_after", 0.0)
        ),
        "land_fraction_before": float(
            world.g("terrain.last_p1134b_land_fraction_before", 0.0)
        ),
        "land_fraction_after": float(
            world.g("terrain.last_p1134b_land_fraction_after", 0.0)
        ),
        "narrow_neck_cells_before": float(
            world.g("terrain.last_p1134b_narrow_neck_cells_before", 0.0)
        ),
        "narrow_neck_cells_after": float(
            world.g("terrain.last_p1134b_narrow_neck_cells_after", 0.0)
        ),
        "narrow_neck_per_1000_before": float(
            world.g("terrain.last_p1134b_narrow_neck_per_1000_before", 0.0)
        ),
        "narrow_neck_per_1000_after": float(
            world.g("terrain.last_p1134b_narrow_neck_per_1000_after", 0.0)
        ),
        "ceiling_fraction": float(
            world.g("terrain.last_p1134b_ceiling_fraction", 0.0)
        ),
        "terminal_land_fraction": float(area[land].sum() / total),
    }
    p1134c_shelf_depth_cap = {
        "schema": "aevum.p1134c_lowres_shelf_depth_cap.v1",
        "adjusted_cell_count": float(
            world.g(
                "terrain.last_p1134c_lowres_shelf_depth_cap_adjusted_cell_count",
                0.0,
            )
        ),
        "adjusted_area_fraction": float(
            world.g(
                "terrain.last_p1134c_lowres_shelf_depth_cap_adjusted_area_fraction",
                0.0,
            )
        ),
        "terminal_land_fraction": float(area[land].sum() / total),
    }
    province_hypsometry = {
        "terrain_province": _category_hypsometry_metrics(
            terrain_province, rel, area, land),
        "continental_detail": _category_hypsometry_metrics(
            continental_detail, rel, area, land),
        "mountain_inventory": _category_hypsometry_metrics(
            mountain_inventory, rel, area, land),
        "crust_type": _category_hypsometry_metrics(crust_type, rel, area, land),
    }

    metrics = {
        "schema": SCHEMA,
        "cells": int(grid.n),
        "time_myr": float(world.time_myr),
        "sea_level_m": sea_level,
        "land_fraction": float(area[land].sum() / total),
        "ocean_fraction": float(area[ocean].sum() / total),
        "land_hypsometry": land_stats,
        "land_elevation_band_fraction": land_bands,
        "ocean_depth_hypsometry": ocean_stats,
        "earth_reference_hypsometry": earth_hypsometry,
        "earth_reference_planform": earth_planform,
        "p109_hypsometry_comparison": p109_comparison,
        "p109_sea_level_consistency": p109_sea_level,
        "terrain_internal_profile": terrain_internal_profile,
        "p120_continental_semantic_geometry": p120_semantic_geometry,
        "province_hypsometry": province_hypsometry,
        "p110a_modern_planform": _p110a_modern_planform_metrics(world, land, ocean),
        "p110b_historical_supercontinent_trajectory": (
            _p110b_historical_supercontinent_trajectory_metrics(archive)
        ),
        **{key: ranks[key] for key in (
            "terminal_active_plate_count",
            "terminal_major_plate_count",
            "terminal_minor_plate_count",
            "terminal_microplate_count",
            "terminal_tiny_plate_count",
            "plate_rank_counts",
            "plate_area_fraction_histogram",
            "plate_area_fraction_top12",
        )},
        "plate_rank_thresholds": ranks["plate_rank_thresholds"],
        "protected_microplate_count": int(len(micro_ids & protected_plate_ids)),
        "protected_plate_count": int(len(protected_plate_ids)),
        "merged_microplate_count": int(merged_microplates),
        "boundary_cell_count_by_kind": boundary_counts,
        "boundary_province_object_count_by_kind": object_kind_counts.get(
            "tectonics.boundary_provinces",
            object_kind_counts.get("tectonics.boundary_objects", {}),
        ),
        "object_count_by_set": object_counts,
        "object_kind_count_by_set": object_kind_counts,
        "ridge_transform_network_continuity": ridge_network,
        "boundary_width_diagnostics": boundary_width,
        "p1116_boundary_linework": p1116_linework,
        "p1116_oceanic_crust_structure": p1116_crust,
        "p11166_trench_arc_refinement": {
            "schema": "aevum.p11166_trench_arc_refinement.v1",
            "accepted": bool(
                world.g("terrain.last_p11166_trench_arc_refinement_accepted", 0.0)
                >= 1.0
            ),
            "components_before": float(
                world.g("terrain.last_p11166_trench_arc_components_before", 0.0)),
            "components_after": float(
                world.g("terrain.last_p11166_trench_arc_components_after", 0.0)),
            "area_fraction_ocean_before": float(
                world.g(
                    "terrain.last_p11166_trench_arc_area_fraction_ocean_before",
                    0.0,
                )
            ),
            "area_fraction_ocean_after": float(
                world.g(
                    "terrain.last_p11166_trench_arc_area_fraction_ocean_after",
                    0.0,
                )
            ),
        },
        "p111632_parent_trench_linework": {
            "schema": "aevum.p111632_parent_trench_linework.v1",
            "accepted": bool(
                world.g(
                    "terrain.last_p111632_parent_trench_refinement_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "components_before": float(
                world.g(
                    "terrain.last_p111632_parent_trench_components_before",
                    0.0,
                )
            ),
            "components_after": float(
                world.g(
                    "terrain.last_p111632_parent_trench_components_after",
                    0.0,
                )
            ),
            "area_fraction_ocean_before": float(
                world.g(
                    "terrain.last_p111632_parent_trench_area_fraction_ocean_before",
                    0.0,
                )
            ),
            "area_fraction_ocean_after": float(
                world.g(
                    "terrain.last_p111632_parent_trench_area_fraction_ocean_after",
                    0.0,
                )
            ),
            "parent_overlap_before": float(
                world.g(
                    "terrain.last_p111632_parent_trench_parent_overlap_before",
                    0.0,
                )
            ),
            "parent_overlap_after": float(
                world.g(
                    "terrain.last_p111632_parent_trench_parent_overlap_after",
                    0.0,
                )
            ),
            "added_area_fraction_ocean": float(
                world.g(
                    "terrain.last_p111632_parent_trench_added_area_fraction_ocean",
                    0.0,
                )
            ),
            "removed_area_fraction_ocean": float(
                world.g(
                    "terrain.last_p111632_parent_trench_removed_area_fraction_ocean",
                    0.0,
                )
            ),
        },
        "p11167_convergent_parent_linework": {
            "schema": "aevum.p11167_convergent_parent_linework.v1",
            "convergent_parent_accepted": bool(
                world.g(
                    "tectonics.last_p11167_convergent_parent_refinement_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "subduction_parent_accepted": bool(
                world.g(
                    "tectonics.last_p11167_subduction_parent_refinement_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "convergent_components_before": float(
                world.g("tectonics.last_p11167_convergent_parent_components_before", 0.0)
            ),
            "convergent_components_after": float(
                world.g("tectonics.last_p11167_convergent_parent_components_after", 0.0)
            ),
            "subduction_components_before": float(
                world.g("tectonics.last_p11167_subduction_parent_components_before", 0.0)
            ),
            "subduction_components_after": float(
                world.g("tectonics.last_p11167_subduction_parent_components_after", 0.0)
            ),
            "convergent_cells_before": float(
                world.g("tectonics.last_p11167_convergent_parent_cells_before", 0.0)
            ),
            "convergent_cells_after": float(
                world.g("tectonics.last_p11167_convergent_parent_cells_after", 0.0)
            ),
            "subduction_cells_before": float(
                world.g("tectonics.last_p11167_subduction_parent_cells_before", 0.0)
            ),
            "subduction_cells_after": float(
                world.g("tectonics.last_p11167_subduction_parent_cells_after", 0.0)
            ),
        },
        "p11168_orogenic_parent_hierarchy": {
            "schema": "aevum.p11168_orogenic_parent_hierarchy.v1",
            "accepted": bool(
                world.g(
                    "terrain.last_p11168_parent_orogen_hierarchy_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "parent_line_cell_count": float(
                world.g("terrain.last_p11168_parent_line_cell_count", 0.0)),
            "parent_line_area_fraction": float(
                world.g("terrain.last_p11168_parent_line_area_fraction", 0.0)),
            "active_orogen_candidate_cell_count": float(
                world.g(
                    "terrain.last_p11168_active_orogen_candidate_cell_count",
                    0.0,
                )
            ),
            "parent_overlap_fraction": float(
                world.g("terrain.last_p11168_parent_overlap_fraction", 0.0)),
            "p146_parent_orogen_corridor_accepted": bool(
                world.g(
                    "terrain.last_p146_parent_orogen_corridor_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p146_parent_orogen_corridor_object_count": float(
                world.g(
                    "terrain.last_p146_parent_orogen_corridor_object_count",
                    0.0,
                )
            ),
            "p146_parent_orogen_corridor_component_count": float(
                world.g(
                    "terrain.last_p146_parent_orogen_corridor_component_count",
                    0.0,
                )
            ),
            "p146_candidate_corridor_component_count": float(
                world.g(
                    "terrain.last_p146_candidate_corridor_component_count",
                    0.0,
                )
            ),
            "p146_candidate_corridor_cell_count": float(
                world.g("terrain.last_p146_candidate_corridor_cell_count", 0.0)
            ),
            "p146_retained_corridor_component_count": float(
                world.g(
                    "terrain.last_p146_retained_corridor_component_count",
                    0.0,
                )
            ),
            "p146_retained_corridor_cell_count": float(
                world.g("terrain.last_p146_retained_corridor_cell_count", 0.0)
            ),
            "p146_pruned_corridor_component_count": float(
                world.g(
                    "terrain.last_p146_pruned_corridor_component_count",
                    0.0,
                )
            ),
            "p146_trunk_parent_overlap_fraction": float(
                world.g(
                    "terrain.last_p146_trunk_parent_overlap_fraction",
                    0.0,
                )
            ),
            "p146_reject_code": float(
                world.g("terrain.last_p146_reject_code", 0.0)
            ),
            "p146_corridor_trunk_cell_count": float(
                world.g("terrain.last_p146_corridor_trunk_cell_count", 0.0)
            ),
            "p146_corridor_branch_cell_count": float(
                world.g("terrain.last_p146_corridor_branch_cell_count", 0.0)
            ),
            "p146_corridor_foreland_cell_count": float(
                world.g("terrain.last_p146_corridor_foreland_cell_count", 0.0)
            ),
            "p146_corridor_bridge_cell_count": float(
                world.g("terrain.last_p146_corridor_bridge_cell_count", 0.0)
            ),
            "p146_corridor_bridge_path_count": float(
                world.g("terrain.last_p146_corridor_bridge_path_count", 0.0)
            ),
            "p147_boundary_group_corridor_accepted": bool(
                world.g(
                    "terrain.last_p147_boundary_group_corridor_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p147_parent_group_count": float(
                world.g("terrain.last_p147_parent_group_count", 0.0)
            ),
            "p147_polyline_group_count": float(
                world.g("terrain.last_p147_polyline_group_count", 0.0)
            ),
            "p147_object_group_count": float(
                world.g("terrain.last_p147_object_group_count", 0.0)
            ),
            "p147_attempted_group_count": float(
                world.g("terrain.last_p147_attempted_group_count", 0.0)
            ),
            "p147_accepted_group_count": float(
                world.g("terrain.last_p147_accepted_group_count", 0.0)
            ),
            "p147_rejected_group_count": float(
                world.g("terrain.last_p147_rejected_group_count", 0.0)
            ),
            "p147_low_coverage_rejected_group_count": float(
                world.g(
                    "terrain.last_p147_low_coverage_rejected_group_count",
                    0.0,
                )
            ),
            "p147_candidate_corridor_component_count": float(
                world.g(
                    "terrain.last_p147_candidate_corridor_component_count",
                    0.0,
                )
            ),
            "p147_candidate_corridor_cell_count": float(
                world.g("terrain.last_p147_candidate_corridor_cell_count", 0.0)
            ),
            "p147_retained_corridor_component_count": float(
                world.g(
                    "terrain.last_p147_retained_corridor_component_count",
                    0.0,
                )
            ),
            "p147_retained_corridor_cell_count": float(
                world.g("terrain.last_p147_retained_corridor_cell_count", 0.0)
            ),
            "p147_trunk_parent_overlap_mean": float(
                world.g("terrain.last_p147_trunk_parent_overlap_mean", 0.0)
            ),
            "p147_trunk_parent_overlap_min": float(
                world.g("terrain.last_p147_trunk_parent_overlap_min", 0.0)
            ),
            "p147_side_to_trunk_ratio": float(
                world.g("terrain.last_p147_side_to_trunk_ratio", 0.0)
            ),
            "p147_aggregate_guard_rejected": bool(
                world.g(
                    "terrain.last_p147_aggregate_guard_rejected",
                    0.0,
                )
                >= 1.0
            ),
            "p147_corridor_object_count": float(
                world.g("terrain.last_p147_corridor_object_count", 0.0)
            ),
            "p147_corridor_cell_count": float(
                world.g("terrain.last_p147_corridor_cell_count", 0.0)
            ),
            "p148_p147_trial_guard_enabled": bool(
                world.g("terrain.last_p148_p147_trial_guard_enabled", 0.0)
                >= 1.0
            ),
            "p148_p147_trial_guard_accepted": bool(
                world.g("terrain.last_p148_p147_trial_guard_accepted", 0.0)
                >= 1.0
            ),
            "p148_p147_trial_guard_rejected": bool(
                world.g("terrain.last_p148_p147_trial_guard_rejected", 0.0)
                >= 1.0
            ),
            "p148_p147_trial_reject_code": float(
                world.g("terrain.last_p148_p147_trial_reject_code", 0.0)
            ),
            "p148_baseline_class_score": float(
                world.g("terrain.last_p148_baseline_class_score", 0.0)
            ),
            "p148_trial_class_score": float(
                world.g("terrain.last_p148_trial_class_score", 0.0)
            ),
            "p148_baseline_class_small_area_fraction": float(
                world.g(
                    "terrain.last_p148_baseline_class_small_area_fraction",
                    0.0,
                )
            ),
            "p148_trial_class_small_area_fraction": float(
                world.g("terrain.last_p148_trial_class_small_area_fraction", 0.0)
            ),
            "p148_baseline_crest_component_count": float(
                world.g("terrain.last_p148_baseline_crest_component_count", 0.0)
            ),
            "p148_trial_crest_component_count": float(
                world.g("terrain.last_p148_trial_crest_component_count", 0.0)
            ),
            "p148_baseline_branch_component_count": float(
                world.g("terrain.last_p148_baseline_branch_component_count", 0.0)
            ),
            "p148_trial_branch_component_count": float(
                world.g("terrain.last_p148_trial_branch_component_count", 0.0)
            ),
            "p148_trial_added_corridor_cell_count": float(
                world.g("terrain.last_p148_trial_added_corridor_cell_count", 0.0)
            ),
            "p149_staged_promotion_enabled": bool(
                world.g("terrain.last_p149_staged_promotion_enabled", 0.0)
                >= 1.0
            ),
            "p149_staged_promotion_accepted": bool(
                world.g("terrain.last_p149_staged_promotion_accepted", 0.0)
                >= 1.0
            ),
            "p149_trunk_trial_accepted": bool(
                world.g("terrain.last_p149_trunk_trial_accepted", 0.0) >= 1.0
            ),
            "p149_trunk_trial_rejected": bool(
                world.g("terrain.last_p149_trunk_trial_rejected", 0.0) >= 1.0
            ),
            "p149_side_trial_accepted": bool(
                world.g("terrain.last_p149_side_trial_accepted", 0.0) >= 1.0
            ),
            "p149_side_trial_rejected": bool(
                world.g("terrain.last_p149_side_trial_rejected", 0.0) >= 1.0
            ),
            "p149_reject_code": float(
                world.g("terrain.last_p149_reject_code", 0.0)
            ),
            "p149_trunk_trial_class_score": float(
                world.g("terrain.last_p149_trunk_trial_class_score", 0.0)
            ),
            "p149_trunk_trial_class_small_area_fraction": float(
                world.g(
                    "terrain.last_p149_trunk_trial_class_small_area_fraction",
                    0.0,
                )
            ),
            "p149_side_trial_class_score": float(
                world.g("terrain.last_p149_side_trial_class_score", 0.0)
            ),
            "p149_side_trial_class_small_area_fraction": float(
                world.g(
                    "terrain.last_p149_side_trial_class_small_area_fraction",
                    0.0,
                )
            ),
            "p149_committed_trunk_cell_count": float(
                world.g("terrain.last_p149_committed_trunk_cell_count", 0.0)
            ),
            "p149_committed_branch_cell_count": float(
                world.g("terrain.last_p149_committed_branch_cell_count", 0.0)
            ),
            "p149_committed_foreland_cell_count": float(
                world.g("terrain.last_p149_committed_foreland_cell_count", 0.0)
            ),
            "p150_group_trunk_repair_enabled": bool(
                world.g("terrain.last_p150_group_trunk_repair_enabled", 0.0)
                >= 1.0
            ),
            "p150_group_trunk_trial_count": float(
                world.g("terrain.last_p150_group_trunk_trial_count", 0.0)
            ),
            "p150_group_trunk_bridge_accepted_count": float(
                world.g(
                    "terrain.last_p150_group_trunk_bridge_accepted_count",
                    0.0,
                )
            ),
            "p150_group_trunk_small_rejected_count": float(
                world.g(
                    "terrain.last_p150_group_trunk_small_rejected_count",
                    0.0,
                )
            ),
            "p150_group_trunk_bridge_cell_count": float(
                world.g("terrain.last_p150_group_trunk_bridge_cell_count", 0.0)
            ),
            "p150_group_trunk_bridge_path_count": float(
                world.g("terrain.last_p150_group_trunk_bridge_path_count", 0.0)
            ),
            "p150_group_trunk_class_score_before_mean": float(
                world.g(
                    "terrain.last_p150_group_trunk_class_score_before_mean",
                    0.0,
                )
            ),
            "p150_group_trunk_class_score_after_mean": float(
                world.g(
                    "terrain.last_p150_group_trunk_class_score_after_mean",
                    0.0,
                )
            ),
            "p150_group_trunk_class_small_before_max": float(
                world.g(
                    "terrain.last_p150_group_trunk_class_small_before_max",
                    0.0,
                )
            ),
            "p150_group_trunk_class_small_after_max": float(
                world.g(
                    "terrain.last_p150_group_trunk_class_small_after_max",
                    0.0,
                )
            ),
            "p150_component_preserving_cap_enabled": bool(
                world.g(
                    "terrain.last_p150_component_preserving_cap_enabled",
                    0.0,
                ) >= 1.0
            ),
            "p150_aggregate_trunk_bridge_accepted": bool(
                world.g(
                    "terrain.last_p150_aggregate_trunk_bridge_accepted",
                    0.0,
                ) >= 1.0
            ),
            "p150_aggregate_trunk_bridge_cell_count": float(
                world.g(
                    "terrain.last_p150_aggregate_trunk_bridge_cell_count",
                    0.0,
                )
            ),
            "p150_aggregate_trunk_bridge_path_count": float(
                world.g(
                    "terrain.last_p150_aggregate_trunk_bridge_path_count",
                    0.0,
                )
            ),
            "p150_aggregate_trunk_class_score_before": float(
                world.g(
                    "terrain.last_p150_aggregate_trunk_class_score_before",
                    0.0,
                )
            ),
            "p150_aggregate_trunk_class_score_after": float(
                world.g(
                    "terrain.last_p150_aggregate_trunk_class_score_after",
                    0.0,
                )
            ),
            "p150_aggregate_trunk_class_small_before": float(
                world.g(
                    "terrain.last_p150_aggregate_trunk_class_small_before",
                    0.0,
                )
            ),
            "p150_aggregate_trunk_class_small_after": float(
                world.g(
                    "terrain.last_p150_aggregate_trunk_class_small_after",
                    0.0,
                )
            ),
            "p151_guarded_object_cap_enabled": bool(
                world.g("terrain.last_p151_guarded_object_cap_enabled", 0.0)
                >= 1.0
            ),
            "p151_trial_component_cap_evaluated": bool(
                world.g("terrain.last_p151_trial_component_cap_evaluated", 0.0)
                >= 1.0
            ),
            "p151_trial_component_cap_accepted": bool(
                world.g("terrain.last_p151_trial_component_cap_accepted", 0.0)
                >= 1.0
            ),
            "p151_trial_component_cap_rejected": bool(
                world.g("terrain.last_p151_trial_component_cap_rejected", 0.0)
                >= 1.0
            ),
            "p151_trial_selected_component_cap": bool(
                world.g("terrain.last_p151_trial_selected_component_cap", 0.0)
                >= 1.0
            ),
            "p151_trial_old_class_score": float(
                world.g("terrain.last_p151_trial_old_class_score", 0.0)
            ),
            "p151_trial_component_class_score": float(
                world.g("terrain.last_p151_trial_component_class_score", 0.0)
            ),
            "p151_trial_old_class_small_area_fraction": float(
                world.g(
                    "terrain.last_p151_trial_old_class_small_area_fraction",
                    0.0,
                )
            ),
            "p151_trial_component_class_small_area_fraction": float(
                world.g(
                    "terrain.last_p151_trial_component_class_small_area_fraction",
                    0.0,
                )
            ),
            "p151_trial_reject_code": float(
                world.g("terrain.last_p151_trial_reject_code", 0.0)
            ),
            "p151_final_component_cap_evaluated": bool(
                world.g("terrain.last_p151_final_component_cap_evaluated", 0.0)
                >= 1.0
            ),
            "p151_final_component_cap_accepted": bool(
                world.g("terrain.last_p151_final_component_cap_accepted", 0.0)
                >= 1.0
            ),
            "p151_final_component_cap_rejected": bool(
                world.g("terrain.last_p151_final_component_cap_rejected", 0.0)
                >= 1.0
            ),
            "p151_final_old_class_score": float(
                world.g("terrain.last_p151_final_old_class_score", 0.0)
            ),
            "p151_final_component_class_score": float(
                world.g("terrain.last_p151_final_component_class_score", 0.0)
            ),
            "p151_final_old_class_small_area_fraction": float(
                world.g(
                    "terrain.last_p151_final_old_class_small_area_fraction",
                    0.0,
                )
            ),
            "p151_final_component_class_small_area_fraction": float(
                world.g(
                    "terrain.last_p151_final_component_class_small_area_fraction",
                    0.0,
                )
            ),
            "p151_final_reject_code": float(
                world.g("terrain.last_p151_final_reject_code", 0.0)
            ),
            "crest_cell_count": float(
                world.g("terrain.last_p11168_crest_cell_count", 0.0)),
            "branch_cell_count": float(
                world.g("terrain.last_p11168_branch_cell_count", 0.0)),
            "foreland_cell_count": float(
                world.g("terrain.last_p11168_foreland_cell_count", 0.0)),
            "crest_area_fraction": float(
                world.g("terrain.last_p11168_crest_area_fraction", 0.0)),
            "branch_area_fraction": float(
                world.g("terrain.last_p11168_branch_area_fraction", 0.0)),
            "foreland_area_fraction": float(
                world.g("terrain.last_p11168_foreland_area_fraction", 0.0)),
            "crest_component_count": float(
                world.g("terrain.last_p11168_crest_component_count", 0.0)),
            "branch_component_count": float(
                world.g("terrain.last_p11168_branch_component_count", 0.0)),
            "foreland_component_count": float(
                world.g("terrain.last_p11168_foreland_component_count", 0.0)),
            "field_crest_area_fraction": float(
                area[orogenic_parent_hierarchy == 3].sum() / total),
            "field_branch_area_fraction": float(
                area[orogenic_parent_hierarchy == 2].sum() / total),
            "field_foreland_area_fraction": float(
                area[orogenic_parent_hierarchy == 1].sum() / total),
            "field_crest_cell_count": int(
                np.count_nonzero(orogenic_parent_hierarchy == 3)
            ),
            "field_branch_cell_count": int(
                np.count_nonzero(orogenic_parent_hierarchy == 2)
            ),
            "field_foreland_cell_count": int(
                np.count_nonzero(orogenic_parent_hierarchy == 1)
            ),
            "field_crest_component_count": int(
                len(_components(grid, orogenic_parent_hierarchy == 3))
            ),
            "field_branch_component_count": int(
                len(_components(grid, orogenic_parent_hierarchy == 2))
            ),
            "field_foreland_component_count": int(
                len(_components(grid, orogenic_parent_hierarchy == 1))
            ),
            "p111610_refinement_accepted": bool(
                world.g(
                    "terrain.last_p111610_hierarchy_refinement_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p111610_crest_components_before_refine": float(
                world.g("terrain.last_p111610_crest_components_before_refine", 0.0)
            ),
            "p111610_crest_components_after_refine": float(
                world.g("terrain.last_p111610_crest_components_after_refine", 0.0)
            ),
            "p111610_branch_components_before_refine": float(
                world.g("terrain.last_p111610_branch_components_before_refine", 0.0)
            ),
            "p111610_branch_components_after_refine": float(
                world.g("terrain.last_p111610_branch_components_after_refine", 0.0)
            ),
            "p111610_bridge_cell_count": float(
                world.g("terrain.last_p111610_bridge_cell_count", 0.0)),
            "p111610_bridge_area_fraction": float(
                world.g("terrain.last_p111610_bridge_area_fraction", 0.0)),
            "p111613_refinement_accepted": bool(
                world.g(
                    "terrain.last_p111613_hierarchy_geometry_refinement_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p111613_branch_components_before_refine": float(
                world.g("terrain.last_p111613_branch_components_before_refine", 0.0)
            ),
            "p111613_branch_components_after_refine": float(
                world.g("terrain.last_p111613_branch_components_after_refine", 0.0)
            ),
            "p111613_foreland_components_before_refine": float(
                world.g("terrain.last_p111613_foreland_components_before_refine", 0.0)
            ),
            "p111613_foreland_components_after_refine": float(
                world.g("terrain.last_p111613_foreland_components_after_refine", 0.0)
            ),
            "p111613_branch_extension_cell_count": float(
                world.g("terrain.last_p111613_branch_extension_cell_count", 0.0)
            ),
            "p111613_branch_extension_area_fraction": float(
                world.g("terrain.last_p111613_branch_extension_area_fraction", 0.0)
            ),
            "p111613_foreland_removed_cell_count": float(
                world.g("terrain.last_p111613_foreland_removed_cell_count", 0.0)
            ),
            "p111613_foreland_removed_area_fraction": float(
                world.g("terrain.last_p111613_foreland_removed_area_fraction", 0.0)
            ),
            "p111613_foreland_removed_component_count": float(
                world.g("terrain.last_p111613_foreland_removed_component_count", 0.0)
            ),
            "p111623_refinement_accepted": bool(
                world.g(
                    "terrain.last_p111623_branch_spine_continuity_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p111623_branch_components_before_refine": float(
                world.g("terrain.last_p111623_branch_components_before_refine", 0.0)
            ),
            "p111623_branch_components_after_refine": float(
                world.g("terrain.last_p111623_branch_components_after_refine", 0.0)
            ),
            "p111623_peak_components_before_refine": float(
                world.g("terrain.last_p111623_peak_components_before_refine", 0.0)
            ),
            "p111623_peak_components_after_refine": float(
                world.g("terrain.last_p111623_peak_components_after_refine", 0.0)
            ),
            "p111623_bridge_cell_count": float(
                world.g("terrain.last_p111623_bridge_cell_count", 0.0)
            ),
            "p111623_bridge_area_fraction": float(
                world.g("terrain.last_p111623_bridge_area_fraction", 0.0)
            ),
            "p111623_bridge_candidate_cell_count": float(
                world.g("terrain.last_p111623_bridge_candidate_cell_count", 0.0)
            ),
            "p111623_bridge_candidate_area_fraction": float(
                world.g("terrain.last_p111623_bridge_candidate_area_fraction", 0.0)
            ),
            "p111623_path_count": float(
                world.g("terrain.last_p111623_path_count", 0.0)
            ),
            "p111623_spine_components_before_refine": float(
                world.g("terrain.last_p111623_spine_components_before_refine", 0.0)
            ),
            "p111623_spine_components_after_refine": float(
                world.g("terrain.last_p111623_spine_components_after_refine", 0.0)
            ),
            "p111623_integrated_spine_cell_count": float(
                world.g("terrain.last_p111623_integrated_spine_cell_count", 0.0)
            ),
            "p111623_integrated_spine_area_fraction": float(
                world.g("terrain.last_p111623_integrated_spine_area_fraction", 0.0)
            ),
            "p111623_integrated_spine_accepted": bool(
                world.g("terrain.last_p111623_integrated_spine_accepted", 0.0) >= 1.0
            ),
            "p111614_spine_refinement_accepted": bool(
                world.g("terrain.last_p111614_spine_refinement_accepted", 0.0) >= 1.0
            ),
            "p111614_crest_spine_cell_count": float(
                world.g("terrain.last_p111614_crest_spine_cell_count", 0.0)
            ),
            "p111614_branch_spine_cell_count": float(
                world.g("terrain.last_p111614_branch_spine_cell_count", 0.0)
            ),
            "p111614_crest_spine_area_fraction": float(
                world.g("terrain.last_p111614_crest_spine_area_fraction", 0.0)
            ),
            "p111614_branch_spine_area_fraction": float(
                world.g("terrain.last_p111614_branch_spine_area_fraction", 0.0)
            ),
            "p111614_crest_spine_component_count": float(
                world.g("terrain.last_p111614_crest_spine_component_count", 0.0)
            ),
            "p111614_branch_spine_component_count": float(
                world.g("terrain.last_p111614_branch_spine_component_count", 0.0)
            ),
            "p111614_crest_width_ratio": float(
                world.g("terrain.last_p111614_crest_width_ratio", 0.0)
            ),
            "p111614_branch_width_ratio": float(
                world.g("terrain.last_p111614_branch_width_ratio", 0.0)
            ),
            "p111614_spine_overlap_valid": bool(
                world.g("terrain.last_p111614_spine_overlap_valid", 0.0) >= 1.0
            ),
            "p111633c_final_spine_refinement_accepted": bool(
                world.g(
                    "terrain.last_p111633c_final_spine_refinement_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p111633c_peak_components_before_refine": float(
                world.g("terrain.last_p111633c_peak_components_before_refine", 0.0)
            ),
            "p111633c_peak_components_after_refine": float(
                world.g("terrain.last_p111633c_peak_components_after_refine", 0.0)
            ),
            "p111633c_spine_less_components_before_refine": float(
                world.g(
                    "terrain.last_p111633c_spine_less_components_before_refine",
                    0.0,
                )
            ),
            "p111633c_spine_less_components_after_refine": float(
                world.g(
                    "terrain.last_p111633c_spine_less_components_after_refine",
                    0.0,
                )
            ),
            "p111633c_orphan_peak_reclassified_cell_count": float(
                world.g(
                    "terrain.last_p111633c_orphan_peak_reclassified_cell_count",
                    0.0,
                )
            ),
            "p111633c_orphan_high_peak_reclassified_cell_count": float(
                world.g(
                    "terrain.last_p111633c_orphan_high_peak_reclassified_cell_count",
                    0.0,
                )
            ),
            "p111633c_orphan_extreme_peak_reclassified_cell_count": float(
                world.g(
                    "terrain.last_p111633c_orphan_extreme_peak_reclassified_cell_count",
                    0.0,
                )
            ),
            "p111633c_final_spine_cell_count": float(
                world.g("terrain.last_p111633c_final_spine_cell_count", 0.0)
            ),
            "p111633c_final_spine_component_count": float(
                world.g("terrain.last_p111633c_final_spine_component_count", 0.0)
            ),
            "p1141_ranked_orogenic_spine_paths_accepted": bool(
                world.g(
                    "terrain.last_p1141_ranked_orogenic_spine_paths_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p1141_ranked_path_used_count": float(
                world.g("terrain.last_p1141_ranked_path_used_count", 0.0)
            ),
            "p1141_ranked_path_changed_count": float(
                world.g("terrain.last_p1141_ranked_path_changed_count", 0.0)
            ),
            "p1141_ranked_path_extra_step_count": float(
                world.g("terrain.last_p1141_ranked_path_extra_step_count", 0.0)
            ),
            "p1142_ranked_path_guard_rejected_count": float(
                world.g("terrain.last_p1142_ranked_path_guard_rejected_count", 0.0)
            ),
            "p1142_ranked_path_directness_rejected_count": float(
                world.g(
                    "terrain.last_p1142_ranked_path_directness_rejected_count",
                    0.0,
                )
            ),
            "p1144_ranked_path_component_guard_rejected_count": float(
                world.g(
                    "terrain.last_p1144_ranked_path_component_guard_rejected_count",
                    0.0,
                )
            ),
            "p1144_ranked_path_rank_rejected_count": float(
                world.g("terrain.last_p1144_ranked_path_rank_rejected_count", 0.0)
            ),
            "p1144_ranked_path_relief_rejected_count": float(
                world.g("terrain.last_p1144_ranked_path_relief_rejected_count", 0.0)
            ),
            "p1144_ranked_path_support_rejected_count": float(
                world.g("terrain.last_p1144_ranked_path_support_rejected_count", 0.0)
            ),
            "p1144_ranked_path_peer_rejected_count": float(
                world.g("terrain.last_p1144_ranked_path_peer_rejected_count", 0.0)
            ),
            "p1144_ranked_path_endpoint_rejected_count": float(
                world.g("terrain.last_p1144_ranked_path_endpoint_rejected_count", 0.0)
            ),
            "p111635_spine_linework_smoothing_accepted": bool(
                world.g(
                    "terrain.last_p111635_spine_linework_smoothing_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p111635_spine_components_before_refine": float(
                world.g("terrain.last_p111635_spine_components_before_refine", 0.0)
            ),
            "p111635_spine_components_after_refine": float(
                world.g("terrain.last_p111635_spine_components_after_refine", 0.0)
            ),
            "p111635_short_spine_components_before_refine": float(
                world.g(
                    "terrain.last_p111635_short_spine_components_before_refine",
                    0.0,
                )
            ),
            "p111635_short_spine_components_after_refine": float(
                world.g(
                    "terrain.last_p111635_short_spine_components_after_refine",
                    0.0,
                )
            ),
            "p111635_branch_attachment_fraction_before": float(
                world.g(
                    "terrain.last_p111635_branch_attachment_fraction_before",
                    0.0,
                )
            ),
            "p111635_branch_attachment_fraction_after": float(
                world.g(
                    "terrain.last_p111635_branch_attachment_fraction_after",
                    0.0,
                )
            ),
            "p111635_bridge_cell_count": float(
                world.g("terrain.last_p111635_bridge_cell_count", 0.0)
            ),
            "p111635_bridge_area_fraction": float(
                world.g("terrain.last_p111635_bridge_area_fraction", 0.0)
            ),
            "p111635_bridge_candidate_cell_count": float(
                world.g("terrain.last_p111635_bridge_candidate_cell_count", 0.0)
            ),
            "p111635_bridge_candidate_area_fraction": float(
                world.g("terrain.last_p111635_bridge_candidate_area_fraction", 0.0)
            ),
            "p111635_path_count": float(
                world.g("terrain.last_p111635_path_count", 0.0)
            ),
            "p124_orogenic_spine_geometry_regularization_accepted": bool(
                world.g(
                    "terrain.last_p124_orogenic_spine_geometry_regularization_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p124_guard_reverted": bool(
                world.g("terrain.last_p124_guard_reverted", 0.0) >= 1.0
            ),
            "p124_candidate_cell_count": float(
                world.g("terrain.last_p124_candidate_cell_count", 0.0)
            ),
            "p124_added_cell_count": float(
                world.g("terrain.last_p124_added_cell_count", 0.0)
            ),
            "p124_component_count_before": float(
                world.g("terrain.last_p124_component_count_before", 0.0)
            ),
            "p124_component_count_after": float(
                world.g("terrain.last_p124_component_count_after", 0.0)
            ),
            "p124_endpoint_count_before": float(
                world.g("terrain.last_p124_endpoint_count_before", 0.0)
            ),
            "p124_endpoint_count_after": float(
                world.g("terrain.last_p124_endpoint_count_after", 0.0)
            ),
            "p124_linework_score_before": float(
                world.g("terrain.last_p124_linework_score_before", 0.0)
            ),
            "p124_linework_score_after": float(
                world.g("terrain.last_p124_linework_score_after", 0.0)
            ),
            "p111636_polar_edge_refinement_accepted": bool(
                world.g(
                    "terrain.last_p111636_polar_edge_refinement_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p111636_polar_peak_area_fraction_before": float(
                world.g(
                    "terrain.last_p111636_polar_peak_area_fraction_before",
                    0.0,
                )
            ),
            "p111636_polar_peak_area_fraction_after": float(
                world.g(
                    "terrain.last_p111636_polar_peak_area_fraction_after",
                    0.0,
                )
            ),
            "p111636_edge_peak_component_count_before": float(
                world.g(
                    "terrain.last_p111636_edge_peak_component_count_before",
                    0.0,
                )
            ),
            "p111636_edge_peak_component_count_after": float(
                world.g(
                    "terrain.last_p111636_edge_peak_component_count_after",
                    0.0,
                )
            ),
            "p111636_candidate_cell_count": float(
                world.g("terrain.last_p111636_candidate_cell_count", 0.0)
            ),
            "p111636_candidate_area_fraction": float(
                world.g("terrain.last_p111636_candidate_area_fraction", 0.0)
            ),
            "p111636_reclassified_cell_count": float(
                world.g("terrain.last_p111636_reclassified_cell_count", 0.0)
            ),
            "p111636_reclassified_area_fraction": float(
                world.g("terrain.last_p111636_reclassified_area_fraction", 0.0)
            ),
            "p111636_polar_reclassified_cell_count": float(
                world.g("terrain.last_p111636_polar_reclassified_cell_count", 0.0)
            ),
            "p111636_edge_reclassified_cell_count": float(
                world.g("terrain.last_p111636_edge_reclassified_cell_count", 0.0)
            ),
            "p111636_extreme_reclassified_cell_count": float(
                world.g("terrain.last_p111636_extreme_reclassified_cell_count", 0.0)
            ),
            "p125_terminal_orogenic_semantic_land_consistency_accepted": bool(
                world.g(
                    "terrain.last_p125_terminal_orogenic_semantic_land_consistency_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p125_guard_reverted": bool(
                world.g("terrain.last_p125_guard_reverted", 0.0) >= 1.0
            ),
            "p125_candidate_cell_count": float(
                world.g("terrain.last_p125_candidate_cell_count", 0.0)
            ),
            "p125_cleared_hierarchy_cell_count": float(
                world.g("terrain.last_p125_cleared_hierarchy_cell_count", 0.0)
            ),
            "p125_cleared_spine_cell_count": float(
                world.g("terrain.last_p125_cleared_spine_cell_count", 0.0)
            ),
            "p126_terminal_orogenic_fringe_support_regularization_accepted": bool(
                world.g(
                    "terrain.last_p126_terminal_orogenic_fringe_support_regularization_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p126_guard_reverted": bool(
                world.g("terrain.last_p126_guard_reverted", 0.0) >= 1.0
            ),
            "p126_candidate_cell_count": float(
                world.g("terrain.last_p126_candidate_cell_count", 0.0)
            ),
            "p126_candidate_area_fraction": float(
                world.g("terrain.last_p126_candidate_area_fraction", 0.0)
            ),
            "p126_cleared_foreland_cell_count": float(
                world.g("terrain.last_p126_cleared_foreland_cell_count", 0.0)
            ),
            "p126_cleared_halo_cell_count": float(
                world.g("terrain.last_p126_cleared_halo_cell_count", 0.0)
            ),
            "p126_cleared_apron_cell_count": float(
                world.g("terrain.last_p126_cleared_apron_cell_count", 0.0)
            ),
            "p126_foreland_component_count_before": float(
                world.g("terrain.last_p126_foreland_component_count_before", 0.0)
            ),
            "p126_foreland_component_count_after": float(
                world.g("terrain.last_p126_foreland_component_count_after", 0.0)
            ),
            "p126_halo_component_count_before": float(
                world.g("terrain.last_p126_halo_component_count_before", 0.0)
            ),
            "p126_halo_component_count_after": float(
                world.g("terrain.last_p126_halo_component_count_after", 0.0)
            ),
            "p126_apron_component_count_before": float(
                world.g("terrain.last_p126_apron_component_count_before", 0.0)
            ),
            "p126_apron_component_count_after": float(
                world.g("terrain.last_p126_apron_component_count_after", 0.0)
            ),
            "p126_tiny_foreland_component_count_before": float(
                world.g(
                    "terrain.last_p126_tiny_foreland_component_count_before",
                    0.0,
                )
            ),
            "p126_tiny_foreland_component_count_after": float(
                world.g(
                    "terrain.last_p126_tiny_foreland_component_count_after",
                    0.0,
                )
            ),
            "p126_tiny_halo_component_count_before": float(
                world.g("terrain.last_p126_tiny_halo_component_count_before", 0.0)
            ),
            "p126_tiny_halo_component_count_after": float(
                world.g("terrain.last_p126_tiny_halo_component_count_after", 0.0)
            ),
            "p126_tiny_apron_component_count_before": float(
                world.g("terrain.last_p126_tiny_apron_component_count_before", 0.0)
            ),
            "p126_tiny_apron_component_count_after": float(
                world.g("terrain.last_p126_tiny_apron_component_count_after", 0.0)
            ),
            "p126_peak_hierarchy_changed_cell_count": float(
                world.g("terrain.last_p126_peak_hierarchy_changed_cell_count", 0.0)
            ),
            "p126_spine_changed_cell_count": float(
                world.g("terrain.last_p126_spine_changed_cell_count", 0.0)
            ),
            "p133_terminal_orogenic_fringe_band_ordering_accepted": bool(
                world.g(
                    "terrain.last_p133_terminal_orogenic_fringe_band_ordering_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p133_guard_reverted": bool(
                world.g("terrain.last_p133_guard_reverted", 0.0) >= 1.0
            ),
            "p133_candidate_cell_count": float(
                world.g("terrain.last_p133_candidate_cell_count", 0.0)
            ),
            "p133_candidate_area_fraction": float(
                world.g("terrain.last_p133_candidate_area_fraction", 0.0)
            ),
            "p133_cleared_foreland_cell_count": float(
                world.g("terrain.last_p133_cleared_foreland_cell_count", 0.0)
            ),
            "p133_cleared_halo_cell_count": float(
                world.g("terrain.last_p133_cleared_halo_cell_count", 0.0)
            ),
            "p133_cleared_apron_cell_count": float(
                world.g("terrain.last_p133_cleared_apron_cell_count", 0.0)
            ),
            "p133_foreland_component_count_before": float(
                world.g("terrain.last_p133_foreland_component_count_before", 0.0)
            ),
            "p133_foreland_component_count_after": float(
                world.g("terrain.last_p133_foreland_component_count_after", 0.0)
            ),
            "p133_halo_component_count_before": float(
                world.g("terrain.last_p133_halo_component_count_before", 0.0)
            ),
            "p133_halo_component_count_after": float(
                world.g("terrain.last_p133_halo_component_count_after", 0.0)
            ),
            "p133_apron_component_count_before": float(
                world.g("terrain.last_p133_apron_component_count_before", 0.0)
            ),
            "p133_apron_component_count_after": float(
                world.g("terrain.last_p133_apron_component_count_after", 0.0)
            ),
            "p133_tiny_foreland_component_count_before": float(
                world.g(
                    "terrain.last_p133_tiny_foreland_component_count_before",
                    0.0,
                )
            ),
            "p133_tiny_foreland_component_count_after": float(
                world.g(
                    "terrain.last_p133_tiny_foreland_component_count_after",
                    0.0,
                )
            ),
            "p133_tiny_halo_component_count_before": float(
                world.g("terrain.last_p133_tiny_halo_component_count_before", 0.0)
            ),
            "p133_tiny_halo_component_count_after": float(
                world.g("terrain.last_p133_tiny_halo_component_count_after", 0.0)
            ),
            "p133_tiny_apron_component_count_before": float(
                world.g("terrain.last_p133_tiny_apron_component_count_before", 0.0)
            ),
            "p133_tiny_apron_component_count_after": float(
                world.g("terrain.last_p133_tiny_apron_component_count_after", 0.0)
            ),
            "p133_band_violation_cell_count": float(
                world.g("terrain.last_p133_band_violation_cell_count", 0.0)
            ),
            "p133_peak_hierarchy_changed_cell_count": float(
                world.g("terrain.last_p133_peak_hierarchy_changed_cell_count", 0.0)
            ),
            "p133_spine_changed_cell_count": float(
                world.g("terrain.last_p133_spine_changed_cell_count", 0.0)
            ),
            "p137_terminal_orogenic_fringe_gap_consolidation_accepted": bool(
                world.g(
                    "terrain.last_p137_terminal_orogenic_fringe_gap_consolidation_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p137_guard_reverted": bool(
                world.g("terrain.last_p137_guard_reverted", 0.0) >= 1.0
            ),
            "p137_candidate_cell_count": float(
                world.g("terrain.last_p137_candidate_cell_count", 0.0)
            ),
            "p137_candidate_area_fraction": float(
                world.g("terrain.last_p137_candidate_area_fraction", 0.0)
            ),
            "p137_added_foreland_cell_count": float(
                world.g("terrain.last_p137_added_foreland_cell_count", 0.0)
            ),
            "p137_added_foreland_area_fraction": float(
                world.g("terrain.last_p137_added_foreland_area_fraction", 0.0)
            ),
            "p137_added_halo_cell_count": float(
                world.g("terrain.last_p137_added_halo_cell_count", 0.0)
            ),
            "p137_added_halo_area_fraction": float(
                world.g("terrain.last_p137_added_halo_area_fraction", 0.0)
            ),
            "p137_added_apron_cell_count": float(
                world.g("terrain.last_p137_added_apron_cell_count", 0.0)
            ),
            "p137_added_apron_area_fraction": float(
                world.g("terrain.last_p137_added_apron_area_fraction", 0.0)
            ),
            "p137_foreland_gap_component_count": float(
                world.g("terrain.last_p137_foreland_gap_component_count", 0.0)
            ),
            "p137_halo_gap_component_count": float(
                world.g("terrain.last_p137_halo_gap_component_count", 0.0)
            ),
            "p137_apron_gap_component_count": float(
                world.g("terrain.last_p137_apron_gap_component_count", 0.0)
            ),
            "p137_foreland_component_count_before": float(
                world.g("terrain.last_p137_foreland_component_count_before", 0.0)
            ),
            "p137_foreland_component_count_after": float(
                world.g("terrain.last_p137_foreland_component_count_after", 0.0)
            ),
            "p137_halo_component_count_before": float(
                world.g("terrain.last_p137_halo_component_count_before", 0.0)
            ),
            "p137_halo_component_count_after": float(
                world.g("terrain.last_p137_halo_component_count_after", 0.0)
            ),
            "p137_apron_component_count_before": float(
                world.g("terrain.last_p137_apron_component_count_before", 0.0)
            ),
            "p137_apron_component_count_after": float(
                world.g("terrain.last_p137_apron_component_count_after", 0.0)
            ),
            "p137_tiny_foreland_component_count_before": float(
                world.g("terrain.last_p137_tiny_foreland_component_count_before", 0.0)
            ),
            "p137_tiny_foreland_component_count_after": float(
                world.g("terrain.last_p137_tiny_foreland_component_count_after", 0.0)
            ),
            "p137_tiny_halo_component_count_before": float(
                world.g("terrain.last_p137_tiny_halo_component_count_before", 0.0)
            ),
            "p137_tiny_halo_component_count_after": float(
                world.g("terrain.last_p137_tiny_halo_component_count_after", 0.0)
            ),
            "p137_tiny_apron_component_count_before": float(
                world.g("terrain.last_p137_tiny_apron_component_count_before", 0.0)
            ),
            "p137_tiny_apron_component_count_after": float(
                world.g("terrain.last_p137_tiny_apron_component_count_after", 0.0)
            ),
            "p137_peak_hierarchy_changed_cell_count": float(
                world.g("terrain.last_p137_peak_hierarchy_changed_cell_count", 0.0)
            ),
            "p137_spine_changed_cell_count": float(
                world.g("terrain.last_p137_spine_changed_cell_count", 0.0)
            ),
            "p138_terminal_orogenic_fringe_component_thickening_accepted": bool(
                world.g(
                    "terrain.last_p138_terminal_orogenic_fringe_component_thickening_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p138_guard_reverted": bool(
                world.g("terrain.last_p138_guard_reverted", 0.0) >= 1.0
            ),
            "p138_candidate_cell_count": float(
                world.g("terrain.last_p138_candidate_cell_count", 0.0)
            ),
            "p138_candidate_area_fraction": float(
                world.g("terrain.last_p138_candidate_area_fraction", 0.0)
            ),
            "p138_added_foreland_cell_count": float(
                world.g("terrain.last_p138_added_foreland_cell_count", 0.0)
            ),
            "p138_added_halo_cell_count": float(
                world.g("terrain.last_p138_added_halo_cell_count", 0.0)
            ),
            "p138_added_apron_cell_count": float(
                world.g("terrain.last_p138_added_apron_cell_count", 0.0)
            ),
            "p138_grown_foreland_component_count": float(
                world.g("terrain.last_p138_grown_foreland_component_count", 0.0)
            ),
            "p138_grown_halo_component_count": float(
                world.g("terrain.last_p138_grown_halo_component_count", 0.0)
            ),
            "p138_grown_apron_component_count": float(
                world.g("terrain.last_p138_grown_apron_component_count", 0.0)
            ),
            "p138_foreland_component_count_before": float(
                world.g("terrain.last_p138_foreland_component_count_before", 0.0)
            ),
            "p138_foreland_component_count_after": float(
                world.g("terrain.last_p138_foreland_component_count_after", 0.0)
            ),
            "p138_halo_component_count_before": float(
                world.g("terrain.last_p138_halo_component_count_before", 0.0)
            ),
            "p138_halo_component_count_after": float(
                world.g("terrain.last_p138_halo_component_count_after", 0.0)
            ),
            "p138_apron_component_count_before": float(
                world.g("terrain.last_p138_apron_component_count_before", 0.0)
            ),
            "p138_apron_component_count_after": float(
                world.g("terrain.last_p138_apron_component_count_after", 0.0)
            ),
            "p138_foreland_small_area_fraction_before": float(
                world.g("terrain.last_p138_foreland_small_area_fraction_before", 0.0)
            ),
            "p138_foreland_small_area_fraction_after": float(
                world.g("terrain.last_p138_foreland_small_area_fraction_after", 0.0)
            ),
            "p138_halo_small_area_fraction_before": float(
                world.g("terrain.last_p138_halo_small_area_fraction_before", 0.0)
            ),
            "p138_halo_small_area_fraction_after": float(
                world.g("terrain.last_p138_halo_small_area_fraction_after", 0.0)
            ),
            "p138_apron_small_area_fraction_before": float(
                world.g("terrain.last_p138_apron_small_area_fraction_before", 0.0)
            ),
            "p138_apron_small_area_fraction_after": float(
                world.g("terrain.last_p138_apron_small_area_fraction_after", 0.0)
            ),
            "p138_peak_hierarchy_changed_cell_count": float(
                world.g("terrain.last_p138_peak_hierarchy_changed_cell_count", 0.0)
            ),
            "p138_spine_changed_cell_count": float(
                world.g("terrain.last_p138_spine_changed_cell_count", 0.0)
            ),
            "p139_terminal_branch_range_component_thickening_accepted": bool(
                world.g(
                    "terrain.last_p139_terminal_branch_range_component_thickening_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p139_guard_reverted": bool(
                world.g("terrain.last_p139_guard_reverted", 0.0) >= 1.0
            ),
            "p139_candidate_cell_count": float(
                world.g("terrain.last_p139_candidate_cell_count", 0.0)
            ),
            "p139_candidate_area_fraction": float(
                world.g("terrain.last_p139_candidate_area_fraction", 0.0)
            ),
            "p139_added_branch_cell_count": float(
                world.g("terrain.last_p139_added_branch_cell_count", 0.0)
            ),
            "p139_cleared_halo_cell_count": float(
                world.g("terrain.last_p139_cleared_halo_cell_count", 0.0)
            ),
            "p139_cleared_apron_cell_count": float(
                world.g("terrain.last_p139_cleared_apron_cell_count", 0.0)
            ),
            "p139_grown_branch_component_count": float(
                world.g("terrain.last_p139_grown_branch_component_count", 0.0)
            ),
            "p139_branch_cell_count_before": float(
                world.g("terrain.last_p139_branch_cell_count_before", 0.0)
            ),
            "p139_branch_cell_count_after": float(
                world.g("terrain.last_p139_branch_cell_count_after", 0.0)
            ),
            "p139_branch_component_count_before": float(
                world.g("terrain.last_p139_branch_component_count_before", 0.0)
            ),
            "p139_branch_component_count_after": float(
                world.g("terrain.last_p139_branch_component_count_after", 0.0)
            ),
            "p139_branch_small_area_fraction_before": float(
                world.g("terrain.last_p139_branch_small_area_fraction_before", 0.0)
            ),
            "p139_branch_small_area_fraction_after": float(
                world.g("terrain.last_p139_branch_small_area_fraction_after", 0.0)
            ),
            "p139_halo_component_count_before": float(
                world.g("terrain.last_p139_halo_component_count_before", 0.0)
            ),
            "p139_halo_component_count_after": float(
                world.g("terrain.last_p139_halo_component_count_after", 0.0)
            ),
            "p139_apron_component_count_before": float(
                world.g("terrain.last_p139_apron_component_count_before", 0.0)
            ),
            "p139_apron_component_count_after": float(
                world.g("terrain.last_p139_apron_component_count_after", 0.0)
            ),
            "p139_crest_changed_cell_count": float(
                world.g("terrain.last_p139_crest_changed_cell_count", 0.0)
            ),
            "p139_spine_changed_cell_count": float(
                world.g("terrain.last_p139_spine_changed_cell_count", 0.0)
            ),
            "p111637a_spine_object_promotion_accepted": bool(
                world.g(
                    "terrain.last_p111637a_spine_object_promotion_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p111637a_spine_components_before": float(
                world.g("terrain.last_p111637a_spine_components_before", 0.0)
            ),
            "p111637a_spine_components_after": float(
                world.g("terrain.last_p111637a_spine_components_after", 0.0)
            ),
            "p111637a_short_spine_components_before": float(
                world.g("terrain.last_p111637a_short_spine_components_before", 0.0)
            ),
            "p111637a_short_spine_components_after": float(
                world.g("terrain.last_p111637a_short_spine_components_after", 0.0)
            ),
            "p111637a_branch_attachment_fraction_before": float(
                world.g(
                    "terrain.last_p111637a_branch_attachment_fraction_before",
                    0.0,
                )
            ),
            "p111637a_branch_attachment_fraction_after": float(
                world.g(
                    "terrain.last_p111637a_branch_attachment_fraction_after",
                    0.0,
                )
            ),
            "p111637a_spine_top3_share_before": float(
                world.g("terrain.last_p111637a_spine_top3_share_before", 0.0)
            ),
            "p111637a_spine_top3_share_after": float(
                world.g("terrain.last_p111637a_spine_top3_share_after", 0.0)
            ),
            "p111637a_linework_score_before": float(
                world.g("terrain.last_p111637a_linework_score_before", 0.0)
            ),
            "p111637a_linework_score_after": float(
                world.g("terrain.last_p111637a_linework_score_after", 0.0)
            ),
            "p111637a_bridge_cell_count": float(
                world.g("terrain.last_p111637a_bridge_cell_count", 0.0)
            ),
            "p111637a_bridge_area_fraction": float(
                world.g("terrain.last_p111637a_bridge_area_fraction", 0.0)
            ),
            "p111637a_candidate_cell_count": float(
                world.g("terrain.last_p111637a_candidate_cell_count", 0.0)
            ),
            "p111637a_candidate_area_fraction": float(
                world.g("terrain.last_p111637a_candidate_area_fraction", 0.0)
            ),
            "p111637a_path_count": float(
                world.g("terrain.last_p111637a_path_count", 0.0)
            ),
            "p132_parent_anchor_spine_promotion_accepted": bool(
                world.g(
                    "terrain.last_p132_parent_anchor_spine_promotion_accepted",
                    0.0,
                ) >= 1.0
            ),
            "p132_candidate_cell_count": float(
                world.g("terrain.last_p132_candidate_cell_count", 0.0)
            ),
            "p132_candidate_area_fraction": float(
                world.g("terrain.last_p132_candidate_area_fraction", 0.0)
            ),
            "p132_promoted_spine_cell_count": float(
                world.g("terrain.last_p132_promoted_spine_cell_count", 0.0)
            ),
            "p132_promoted_spine_area_fraction": float(
                world.g("terrain.last_p132_promoted_spine_area_fraction", 0.0)
            ),
            "p132_parent_aligned_spine_fraction_before": float(
                world.g(
                    "terrain.last_p132_parent_aligned_spine_fraction_before",
                    0.0,
                )
            ),
            "p132_parent_aligned_spine_fraction_after": float(
                world.g(
                    "terrain.last_p132_parent_aligned_spine_fraction_after",
                    0.0,
                )
            ),
            "p132_linework_score_before": float(
                world.g("terrain.last_p132_linework_score_before", 0.0)
            ),
            "p132_linework_score_after": float(
                world.g("terrain.last_p132_linework_score_after", 0.0)
            ),
            "p132_spine_components_before": float(
                world.g("terrain.last_p132_spine_components_before", 0.0)
            ),
            "p132_spine_components_after": float(
                world.g("terrain.last_p132_spine_components_after", 0.0)
            ),
            "p132_short_spine_components_before": float(
                world.g("terrain.last_p132_short_spine_components_before", 0.0)
            ),
            "p132_short_spine_components_after": float(
                world.g("terrain.last_p132_short_spine_components_after", 0.0)
            ),
            "p1145_whole_mask_spine_planner_accepted": bool(
                world.g(
                    "terrain.last_p1145_whole_mask_spine_planner_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p1145_bridge_cell_count": float(
                world.g("terrain.last_p1145_bridge_cell_count", 0.0)
            ),
            "p1145_candidate_cell_count": float(
                world.g("terrain.last_p1145_candidate_cell_count", 0.0)
            ),
            "p1145_path_count": float(
                world.g("terrain.last_p1145_path_count", 0.0)
            ),
            "p1145_linework_score_before": float(
                world.g("terrain.last_p1145_linework_score_before", 0.0)
            ),
            "p1145_linework_score_after": float(
                world.g("terrain.last_p1145_linework_score_after", 0.0)
            ),
            "p1145_spine_components_before": float(
                world.g("terrain.last_p1145_spine_components_before", 0.0)
            ),
            "p1145_spine_components_after": float(
                world.g("terrain.last_p1145_spine_components_after", 0.0)
            ),
            "p1145_short_spine_components_before": float(
                world.g("terrain.last_p1145_short_spine_components_before", 0.0)
            ),
            "p1145_short_spine_components_after": float(
                world.g("terrain.last_p1145_short_spine_components_after", 0.0)
            ),
            "p1145_branch_attachment_fraction_before": float(
                world.g("terrain.last_p1145_branch_attachment_fraction_before", 0.0)
            ),
            "p1145_branch_attachment_fraction_after": float(
                world.g("terrain.last_p1145_branch_attachment_fraction_after", 0.0)
            ),
            "p1145_spine_top3_share_before": float(
                world.g("terrain.last_p1145_spine_top3_share_before", 0.0)
            ),
            "p1145_spine_top3_share_after": float(
                world.g("terrain.last_p1145_spine_top3_share_after", 0.0)
            ),
            "p142_class_aware_combined_score": float(
                world.g("terrain.last_p142_class_aware_combined_score", 0.0)
            ),
            "p142_class_aware_class_score": float(
                world.g("terrain.last_p142_class_aware_class_score", 0.0)
            ),
            "p142_class_aware_crest_small_area_fraction": float(
                world.g(
                    "terrain.last_p142_class_aware_crest_small_area_fraction",
                    0.0,
                )
            ),
            "p142_class_aware_branch_small_area_fraction": float(
                world.g(
                    "terrain.last_p142_class_aware_branch_small_area_fraction",
                    0.0,
                )
            ),
            "p142_class_aware_class_small_area_fraction": float(
                world.g(
                    "terrain.last_p142_class_aware_class_small_area_fraction",
                    0.0,
                )
            ),
            "p142_class_aware_crest_component_count": float(
                world.g(
                    "terrain.last_p142_class_aware_crest_component_count",
                    0.0,
                )
            ),
            "p142_class_aware_branch_component_count": float(
                world.g(
                    "terrain.last_p142_class_aware_branch_component_count",
                    0.0,
                )
            ),
            "p142_class_aware_blind_spot": bool(
                world.g("terrain.last_p142_class_aware_blind_spot", 0.0) >= 1.0
            ),
            "p143_class_repair_needed": bool(
                world.g("terrain.last_p143_class_repair_needed", 0.0) >= 1.0
            ),
            "p143_class_profile_improved": bool(
                world.g("terrain.last_p143_class_profile_improved", 0.0) >= 1.0
            ),
            "p143_class_score_before": float(
                world.g("terrain.last_p143_class_score_before", 0.0)
            ),
            "p143_class_score_after": float(
                world.g("terrain.last_p143_class_score_after", 0.0)
            ),
            "p143_class_small_area_fraction_before": float(
                world.g("terrain.last_p143_class_small_area_fraction_before", 0.0)
            ),
            "p143_class_small_area_fraction_after": float(
                world.g("terrain.last_p143_class_small_area_fraction_after", 0.0)
            ),
            "p143_crest_small_area_fraction_before": float(
                world.g("terrain.last_p143_crest_small_area_fraction_before", 0.0)
            ),
            "p143_crest_small_area_fraction_after": float(
                world.g("terrain.last_p143_crest_small_area_fraction_after", 0.0)
            ),
            "p143_branch_small_area_fraction_before": float(
                world.g("terrain.last_p143_branch_small_area_fraction_before", 0.0)
            ),
            "p143_branch_small_area_fraction_after": float(
                world.g("terrain.last_p143_branch_small_area_fraction_after", 0.0)
            ),
            "p144_class_path_option_count": float(
                world.g("terrain.last_p144_class_path_option_count", 0.0)
            ),
            "p144_class_path_selected_count": float(
                world.g("terrain.last_p144_class_path_selected_count", 0.0)
            ),
            "p144_crest_path_option_count": float(
                world.g("terrain.last_p144_crest_path_option_count", 0.0)
            ),
            "p144_branch_path_option_count": float(
                world.g("terrain.last_p144_branch_path_option_count", 0.0)
            ),
            "p144_class_attempted_path_count": float(
                world.g("terrain.last_p144_class_attempted_path_count", 0.0)
            ),
            "p144_class_found_path_count": float(
                world.g("terrain.last_p144_class_found_path_count", 0.0)
            ),
            "p144_crest_promoted_spine_cell_count": float(
                world.g("terrain.last_p144_crest_promoted_spine_cell_count", 0.0)
            ),
            "p144_branch_promoted_spine_cell_count": float(
                world.g("terrain.last_p144_branch_promoted_spine_cell_count", 0.0)
            ),
            "p144_class_promoted_spine_cell_count": float(
                world.g("terrain.last_p144_class_promoted_spine_cell_count", 0.0)
            ),
            "p144_class_path_profile_rejected_count": float(
                world.g("terrain.last_p144_class_path_profile_rejected_count", 0.0)
            ),
            "p145_class_hierarchy_component_consolidation_accepted": bool(
                world.g(
                    "terrain.last_p145_class_hierarchy_component_consolidation_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p145_removed_crest_component_count": float(
                world.g("terrain.last_p145_removed_crest_component_count", 0.0)
            ),
            "p145_removed_branch_component_count": float(
                world.g("terrain.last_p145_removed_branch_component_count", 0.0)
            ),
            "p145_removed_crest_cell_count": float(
                world.g("terrain.last_p145_removed_crest_cell_count", 0.0)
            ),
            "p145_removed_branch_cell_count": float(
                world.g("terrain.last_p145_removed_branch_cell_count", 0.0)
            ),
            "p145_removed_crest_area_fraction": float(
                world.g("terrain.last_p145_removed_crest_area_fraction", 0.0)
            ),
            "p145_removed_branch_area_fraction": float(
                world.g("terrain.last_p145_removed_branch_area_fraction", 0.0)
            ),
            "p145_added_crest_cell_count": float(
                world.g("terrain.last_p145_added_crest_cell_count", 0.0)
            ),
            "p145_added_branch_cell_count": float(
                world.g("terrain.last_p145_added_branch_cell_count", 0.0)
            ),
            "p145_added_crest_area_fraction": float(
                world.g("terrain.last_p145_added_crest_area_fraction", 0.0)
            ),
            "p145_added_branch_area_fraction": float(
                world.g("terrain.last_p145_added_branch_area_fraction", 0.0)
            ),
            "p145_bridge_path_count": float(
                world.g("terrain.last_p145_bridge_path_count", 0.0)
            ),
            "p145_class_score_before": float(
                world.g("terrain.last_p145_class_score_before", 0.0)
            ),
            "p145_class_score_after": float(
                world.g("terrain.last_p145_class_score_after", 0.0)
            ),
            "p145_class_small_area_fraction_before": float(
                world.g("terrain.last_p145_class_small_area_fraction_before", 0.0)
            ),
            "p145_class_small_area_fraction_after": float(
                world.g("terrain.last_p145_class_small_area_fraction_after", 0.0)
            ),
            "p145_crest_small_area_fraction_before": float(
                world.g("terrain.last_p145_crest_small_area_fraction_before", 0.0)
            ),
            "p145_crest_small_area_fraction_after": float(
                world.g("terrain.last_p145_crest_small_area_fraction_after", 0.0)
            ),
            "p145_branch_small_area_fraction_before": float(
                world.g("terrain.last_p145_branch_small_area_fraction_before", 0.0)
            ),
            "p145_branch_small_area_fraction_after": float(
                world.g("terrain.last_p145_branch_small_area_fraction_after", 0.0)
            ),
            "p1146_reject_code": int(
                round(float(world.g("terrain.last_p1146_reject_code", 0.0)))
            ),
            "p1146_reject_reason": _p1146_reject_reason_name(
                world.g("terrain.last_p1146_reject_code", 0.0)
            ),
            "p1146_support_component_count": float(
                world.g("terrain.last_p1146_support_component_count", 0.0)
            ),
            "p1146_multi_spine_support_component_count": float(
                world.g("terrain.last_p1146_multi_spine_support_component_count", 0.0)
            ),
            "p1146_attempted_path_count": float(
                world.g("terrain.last_p1146_attempted_path_count", 0.0)
            ),
            "p1146_found_path_count": float(
                world.g("terrain.last_p1146_found_path_count", 0.0)
            ),
            "p1146_bridge_path_count": float(
                world.g("terrain.last_p1146_bridge_path_count", 0.0)
            ),
            "p1146_attempted_bridge_cell_count": float(
                world.g("terrain.last_p1146_attempted_bridge_cell_count", 0.0)
            ),
            "p1146_trial_linework_score": float(
                world.g("terrain.last_p1146_trial_linework_score", 0.0)
            ),
            "p1146_trial_spine_components": float(
                world.g("terrain.last_p1146_trial_spine_components", 0.0)
            ),
            "p1146_trial_short_spine_components": float(
                world.g("terrain.last_p1146_trial_short_spine_components", 0.0)
            ),
            "p1146_trial_branch_attachment_fraction": float(
                world.g("terrain.last_p1146_trial_branch_attachment_fraction", 0.0)
            ),
            "p1146_trial_spine_top3_share": float(
                world.g("terrain.last_p1146_trial_spine_top3_share", 0.0)
            ),
            "p1147_pair_option_count": float(
                world.g("terrain.last_p1147_pair_option_count", 0.0)
            ),
            "p1147_pair_selected_count": float(
                world.g("terrain.last_p1147_pair_selected_count", 0.0)
            ),
            "p1147_pair_balance_rejected_count": float(
                world.g("terrain.last_p1147_pair_balance_rejected_count", 0.0)
            ),
            "p1147_pair_profile_rejected_count": float(
                world.g("terrain.last_p1147_pair_profile_rejected_count", 0.0)
            ),
            "p1148_terminal_proxy_enabled": bool(
                world.g("terrain.last_p1148_terminal_proxy_enabled", 0.0) >= 1.0
            ),
            "p1148_terminal_proxy_score_before": float(
                world.g("terrain.last_p1148_terminal_proxy_score_before", 0.0)
            ),
            "p1148_terminal_proxy_score_after": float(
                world.g("terrain.last_p1148_terminal_proxy_score_after", 0.0)
            ),
            "p1148_terminal_proxy_component_count_before": float(
                world.g(
                    "terrain.last_p1148_terminal_proxy_component_count_before",
                    0.0,
                )
            ),
            "p1148_terminal_proxy_component_count_after": float(
                world.g(
                    "terrain.last_p1148_terminal_proxy_component_count_after",
                    0.0,
                )
            ),
            "p1148_terminal_proxy_short_count_before": float(
                world.g("terrain.last_p1148_terminal_proxy_short_count_before", 0.0)
            ),
            "p1148_terminal_proxy_short_count_after": float(
                world.g("terrain.last_p1148_terminal_proxy_short_count_after", 0.0)
            ),
            "p111615_refinement_accepted": bool(
                world.g(
                    "terrain.last_p111615_belt_morphology_refinement_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p111615_shoulder_halo_cell_count": float(
                world.g("terrain.last_p111615_shoulder_halo_cell_count", 0.0)
            ),
            "p111615_shoulder_halo_area_fraction": float(
                world.g("terrain.last_p111615_shoulder_halo_area_fraction", 0.0)
            ),
            "p111615_shoulder_halo_component_count": float(
                world.g("terrain.last_p111615_shoulder_halo_component_count", 0.0)
            ),
            "p111615_field_shoulder_halo_area_fraction": float(
                area[orogenic_shoulder_halo > 0].sum() / total
            ),
            "p111615_crest_pruned_cell_count": float(
                world.g("terrain.last_p111615_crest_pruned_cell_count", 0.0)
            ),
            "p111615_branch_pruned_cell_count": float(
                world.g("terrain.last_p111615_branch_pruned_cell_count", 0.0)
            ),
            "p111615_crest_width_ratio_before": float(
                world.g("terrain.last_p111615_crest_width_ratio_before", 0.0)
            ),
            "p111615_branch_width_ratio_before": float(
                world.g("terrain.last_p111615_branch_width_ratio_before", 0.0)
            ),
            "p111615_peak_hierarchy_cell_count_before": float(
                world.g("terrain.last_p111615_peak_hierarchy_cell_count_before", 0.0)
            ),
            "p111615_peak_hierarchy_cell_count_after": float(
                world.g("terrain.last_p111615_peak_hierarchy_cell_count_after", 0.0)
            ),
            "p111615_high_peak_removed_cell_count": float(
                world.g("terrain.last_p111615_high_peak_removed_cell_count", 0.0)
            ),
            "p111615_halo_hierarchy_overlap_valid": bool(
                world.g("terrain.last_p111615_halo_hierarchy_overlap_valid", 0.0)
                >= 1.0
            ),
            "p111616_refinement_accepted": bool(
                world.g(
                    "terrain.last_p111616_highlat_morphology_refinement_accepted",
                    0.0,
                )
                >= 1.0
            ),
            "p111616_highlat_candidate_cell_count": float(
                world.g("terrain.last_p111616_highlat_candidate_cell_count", 0.0)
            ),
            "p111616_highland_apron_cell_count": float(
                world.g("terrain.last_p111616_highland_apron_cell_count", 0.0)
            ),
            "p111616_highland_apron_area_fraction": float(
                world.g("terrain.last_p111616_highland_apron_area_fraction", 0.0)
            ),
            "p111616_highland_apron_component_count": float(
                world.g("terrain.last_p111616_highland_apron_component_count", 0.0)
            ),
            "p111616_field_highland_apron_area_fraction": float(
                area[orogenic_highland_apron > 0].sum() / total
            ),
            "p111616_highlat_halo_extension_cell_count": float(
                world.g("terrain.last_p111616_highlat_halo_extension_cell_count", 0.0)
            ),
            "p111616_highlat_halo_extension_area_fraction": float(
                world.g("terrain.last_p111616_highlat_halo_extension_area_fraction", 0.0)
            ),
            "p111616_crest_reclassified_cell_count": float(
                world.g("terrain.last_p111616_crest_reclassified_cell_count", 0.0)
            ),
            "p111616_branch_reclassified_cell_count": float(
                world.g("terrain.last_p111616_branch_reclassified_cell_count", 0.0)
            ),
            "p111616_peak_hierarchy_cell_count_before": float(
                world.g("terrain.last_p111616_peak_hierarchy_cell_count_before", 0.0)
            ),
            "p111616_peak_hierarchy_cell_count_after": float(
                world.g("terrain.last_p111616_peak_hierarchy_cell_count_after", 0.0)
            ),
            "p111616_extreme_peak_reclassified_cell_count": float(
                world.g("terrain.last_p111616_extreme_peak_reclassified_cell_count", 0.0)
            ),
            "p111616_highland_apron_overlap_valid": bool(
                world.g("terrain.last_p111616_highland_apron_overlap_valid", 0.0)
                >= 1.0
            ),
        },
        "collision_plateau_area_fraction": float(collision_area),
        "high_mountain_coherence": high_mountain,
        "p166167_terminal_relief_planform": relief_planform,
        "p111633_orogenic_belt_morphology": orogenic_belt,
        "p111634_spine_aligned_elevation_morphology": spine_aligned_orogen,
        "p111635_orogenic_spine_linework": spine_linework,
        "p111636_polar_edge_orogen_overclassification": polar_edge_orogen,
        "p1134_lowres_false_oceanic_island_cleanup": p1134_cleanup,
        "p1134b_lowres_articulation_neck_widening": p1134b_neck_widening,
        "p1134c_lowres_shelf_depth_cap": p1134c_shelf_depth_cap,
        "active_margin_arc_trench_adjacency_fraction": float(
            arc_trench["active_margin_arc_trench_adjacency_fraction"]),
        "active_margin_arc_trench_adjacency": arc_trench,
        "island_arc_chain_count": int(
            object_kind_counts.get("terrain.arc_plume_landforms", {}).get("island_arc", 0)),
        "microcontinent_object_count": int(
            object_kind_counts.get("terrain.arc_plume_landforms", {}).get("microcontinent", 0)),
        "parented_oceanic_island_chain_count": int(
            ocean_feature["parented_oceanic_island_chain_count"]),
        "deep_trench_fraction_below_6000m": float(
            ocean_feature["deep_trench_fraction_below_6000m"]),
        "p1115_unsupported_open_ocean_shoal_fraction": float(
            ocean_feature["p1115_unsupported_open_ocean_shoal_fraction"]),
        "p1115_object_backed_open_ocean_shoal_fraction": float(
            ocean_feature["p1115_object_backed_open_ocean_shoal_fraction"]),
        "p1115_object_backed_ocean_relief_fraction": float(
            ocean_feature["p1115_object_backed_ocean_relief_fraction"]),
        "p1115_emerged_object_backed_island_fraction_world": float(
            ocean_feature["p1115_emerged_object_backed_island_fraction_world"]),
        "ocean_feature_metrics": ocean_feature,
        "terrain_province_area_fraction": _code_area_fraction(terrain_province, area, total),
        "depth_province_area_fraction": _code_area_fraction(depth_province, area, total),
        "per_plate": ranks["per_plate"],
    }
    metrics["acceptance"] = {
        "p107_0_observability_complete": True,
        "raw_arrays_required": True,
        "quality_improvement_not_claimed": True,
        "has_plate_rank_metrics": metrics["terminal_active_plate_count"] >= 0,
        "has_boundary_object_counts": isinstance(
            metrics["boundary_province_object_count_by_kind"], dict),
        "has_ocean_feature_metrics": isinstance(metrics["ocean_feature_metrics"], dict),
        "has_p108_boundary_width_metrics": isinstance(
            metrics["boundary_width_diagnostics"], dict),
        "has_p1116_boundary_linework_metrics": isinstance(
            metrics["p1116_boundary_linework"], dict),
        "has_p1116_oceanic_crust_structure_metrics": isinstance(
            metrics["p1116_oceanic_crust_structure"], dict),
        "has_p108_high_mountain_metrics": isinstance(
            metrics["high_mountain_coherence"], dict),
        "has_p166167_terminal_relief_planform_metrics": isinstance(
            metrics["p166167_terminal_relief_planform"], dict),
        "has_p111633_orogenic_belt_morphology_metrics": isinstance(
            metrics["p111633_orogenic_belt_morphology"], dict),
        "has_p111634_spine_aligned_elevation_morphology_metrics": isinstance(
            metrics["p111634_spine_aligned_elevation_morphology"], dict),
        "has_p111635_orogenic_spine_linework_metrics": isinstance(
            metrics["p111635_orogenic_spine_linework"], dict),
        "has_p111636_polar_edge_orogen_overclassification_metrics": isinstance(
            metrics["p111636_polar_edge_orogen_overclassification"], dict),
        "has_p120_continental_semantic_geometry_metrics": isinstance(
            metrics["p120_continental_semantic_geometry"], dict),
        "has_p109_hypsometry_metrics": isinstance(
            metrics["p109_hypsometry_comparison"], dict),
        "has_p109_sea_level_consistency_metrics": isinstance(
            metrics["p109_sea_level_consistency"], dict),
        "has_p110a_planform_metrics": isinstance(
            metrics["p110a_modern_planform"], dict),
        "has_p110b_historical_supercontinent_trajectory": isinstance(
            metrics["p110b_historical_supercontinent_trajectory"], dict),
        "within_p110a_modern_planform_envelope": bool(
            metrics["p110a_modern_planform"].get(
                "within_p110a_modern_planform_envelope", False)
        ),
        "p1115_unsupported_open_ocean_shoal_fraction": float(
            metrics["p1115_unsupported_open_ocean_shoal_fraction"]),
        "p1115_object_backed_open_ocean_shoal_fraction": float(
            metrics["p1115_object_backed_open_ocean_shoal_fraction"]),
        "p1115_object_backed_ocean_relief_fraction": float(
            metrics["p1115_object_backed_ocean_relief_fraction"]),
        "p1115_emerged_object_backed_island_fraction_world": float(
            metrics["p1115_emerged_object_backed_island_fraction_world"]),
    }
    return metrics


def write_terminal_audit(
    world: Any,
    outdir: Path,
    *,
    events: Iterable[Any] = (),
    archive: Any | None = None,
    render_assets: bool = True,
    include_earth_reference: bool = True,
    contact_sheet: bool = True,
) -> dict[str, Any]:
    """Write P107 terminal audit files for an already completed world."""
    outdir.mkdir(parents=True, exist_ok=True)
    metrics = terminal_audit_metrics(world, events, archive=archive)
    arrays, array_manifest = _audit_arrays(world)
    npz_path = outdir / "p107_terminal_arrays.npz"
    np.savez_compressed(npz_path, **arrays)
    metrics["array_archive"] = {
        "path": str(npz_path),
        "key_count": int(len(arrays)),
        "manifest": array_manifest,
    }
    assets: dict[str, str] = {}
    if render_assets:
        from aevum import render

        render.render_world(world, outdir)
        for path in sorted(outdir.glob("*.png")):
            assets[path.name] = str(path)
    if include_earth_reference:
        earth_path = _render_earth_reference_same_grid(world, outdir)
        if earth_path is not None:
            assets[earth_path.name] = str(earth_path)
            metrics["earth_reference"] = {
                "available": True,
                "path": str(earth_path),
                "source": "data/reference/etopo5/ETOPO5.DAT",
            }
        else:
            metrics["earth_reference"] = {
                "available": False,
                "path": "",
                "source": "data/reference/etopo5/ETOPO5.DAT",
            }
    if contact_sheet:
        sheet = _render_contact_sheet(world, outdir, assets)
        if sheet is not None:
            assets[sheet.name] = str(sheet)
    metrics["assets"] = assets
    metrics_path = outdir / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, default=_json_default) + "\n")
    return metrics


def run_p107_audit(config: P107AuditConfig, outdir: Path | None = None) -> dict[str, Any]:
    """Run the configured P107 audit ladder and write per-run artefacts."""
    if outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
    entries: list[dict[str, Any]] = []
    t0 = time.time()
    for idx, run in enumerate(config.runs):
        stage_seconds: dict[str, float] = {}
        stage_start = time.perf_counter()
        label = run.label or f"{int(run.cells)}cells_{int(run.n_plates)}p"
        spec = get_preset(config.preset)
        spec.grid_cells = int(run.cells)
        spec.n_plates = int(run.n_plates)
        spec.t_end_myr = float(config.t_end_myr)
        if run.seed is not None:
            spec.seed = int(run.seed)
        stage_seconds["configure"] = float(time.perf_counter() - stage_start)

        stage_start = time.perf_counter()
        engine = Engine.build(spec)
        stage_seconds["build"] = float(time.perf_counter() - stage_start)

        stage_start = time.perf_counter()
        scheduler_modules = _select_scheduler_modules(
            engine.scheduler,
            config.enabled_modules,
        )
        stage_seconds["select_scheduler_modules"] = float(
            time.perf_counter() - stage_start
        )

        stage_start = time.perf_counter()
        if config.enable_ranked_plate_policy:
            engine.world.set_g("tectonics.enable_p107_ranked_plate_policy", 1.0)
        if config.enable_boundary_province_response:
            engine.world.set_g("terrain.enable_p107_boundary_province_response", 1.0)
        if config.enable_p108_boundary_width_guard:
            engine.world.set_g("tectonics.enable_p108_boundary_width_guard", 1.0)
            engine.world.set_g("terrain.enable_p108_boundary_width_guard", 1.0)
        if config.enable_p108_high_mountain_coherence:
            engine.world.set_g("terrain.enable_p108_high_mountain_coherence", 1.0)
        if config.enable_p109_continental_hypsometry_rebalance:
            engine.world.set_g(
                "terrain.enable_p109_continental_hypsometry_rebalance", 1.0)
        for key, value in sorted(config.global_overrides.items()):
            engine.world.set_g(str(key), float(value))
        stage_seconds["configure_world"] = float(time.perf_counter() - stage_start)

        stage_start = time.perf_counter()
        engine.run(n_frames=int(config.frames))
        run_seconds = float(time.perf_counter() - stage_start)
        stage_seconds["engine_run"] = run_seconds
        run_dir = None if outdir is None else outdir / f"{idx:02d}_{label}"

        stage_start = time.perf_counter()
        metrics = (
            terminal_audit_metrics(
                engine.world,
                engine.bus.events,
                archive=engine.archive,
            )
            if run_dir is None
            else write_terminal_audit(
                engine.world,
                run_dir,
                events=engine.bus.events,
                archive=engine.archive,
                render_assets=config.render_world_assets,
                include_earth_reference=config.include_earth_reference,
                contact_sheet=config.render_contact_sheet,
            )
        )
        stage_seconds["terminal_audit_write"] = float(
            time.perf_counter() - stage_start
        )
        entry = {
            "label": label,
            "preset": config.preset,
            "cells": int(run.cells),
            "n_plates": int(run.n_plates),
            "seed": int(spec.seed),
            "t_end_myr": float(config.t_end_myr),
            "frames": int(config.frames),
            "run_seconds": float(run_seconds),
            "stage_seconds": stage_seconds,
            "total_seconds": float(sum(stage_seconds.values())),
            "module_seconds": _scheduler_module_seconds(engine.scheduler.history),
            "terrain_internal_profile": _scheduler_terrain_internal_profile(
                engine.scheduler.history),
            "scheduler_modules": scheduler_modules,
            "outdir": "" if run_dir is None else str(run_dir),
            "metrics": _compact_entry_metrics(metrics),
        }
        entries.append(entry)
    summary = {
        "schema": SUMMARY_SCHEMA,
        "config": {
            "preset": config.preset,
            "t_end_myr": float(config.t_end_myr),
            "frames": int(config.frames),
            "runs": [
                {
                    "cells": int(run.cells),
                    "n_plates": int(run.n_plates),
                    "seed": run.seed,
                    "label": run.label,
                }
                for run in config.runs
            ],
            "global_overrides": {
                str(key): float(value)
                for key, value in sorted(config.global_overrides.items())
            },
            "enable_ranked_plate_policy": bool(config.enable_ranked_plate_policy),
            "enable_boundary_province_response": bool(
                config.enable_boundary_province_response),
            "enable_p108_boundary_width_guard": bool(
                config.enable_p108_boundary_width_guard),
            "enable_p108_high_mountain_coherence": bool(
                config.enable_p108_high_mountain_coherence),
            "enable_p109_continental_hypsometry_rebalance": bool(
                config.enable_p109_continental_hypsometry_rebalance),
            "enabled_modules": (
                None
                if config.enabled_modules is None
                else [str(name) for name in config.enabled_modules]
            ),
            "render_world_assets": bool(config.render_world_assets),
            "render_contact_sheet": bool(config.render_contact_sheet),
            "include_earth_reference": bool(config.include_earth_reference),
        },
        "runtime_seconds": float(time.time() - t0),
        "entries": entries,
        "acceptance": {
            "p107_0_observability_complete": all(
                entry["metrics"]["array_archive_present"] for entry in entries
            ) if outdir is not None else True,
            "has_8000_tier": any(entry["cells"] >= 8000 for entry in entries),
            "has_24000_tier": any(entry["cells"] >= 24000 for entry in entries),
            "quality_improvement_not_claimed": True,
        },
    }
    if outdir is not None:
        (outdir / "p107_audit_summary.json").write_text(
            json.dumps(summary, indent=2, default=_json_default) + "\n")
    return summary


def _select_scheduler_modules(
    scheduler: Any,
    enabled_modules: tuple[str, ...] | None,
) -> dict[str, Any]:
    configured = [str(sm.module.name) for sm in scheduler.modules]
    if enabled_modules is None:
        return {
            "mode": "full",
            "configured": configured,
            "requested": None,
            "enabled": configured,
            "disabled": [],
        }

    requested = tuple(dict.fromkeys(str(name) for name in enabled_modules))
    if not requested:
        raise ValueError("enabled_modules cannot be empty")
    configured_set = set(configured)
    unknown = [name for name in requested if name not in configured_set]
    if unknown:
        raise ValueError(
            "enabled_modules contains unknown scheduler module(s): "
            + ", ".join(unknown)
        )
    requested_set = set(requested)
    scheduler.modules[:] = [
        sm for sm in scheduler.modules if str(sm.module.name) in requested_set
    ]
    enabled = [str(sm.module.name) for sm in scheduler.modules]
    disabled = [name for name in configured if name not in requested_set]
    return {
        "mode": "allowlist",
        "configured": configured,
        "requested": list(requested),
        "enabled": enabled,
        "disabled": disabled,
    }


def _scheduler_module_seconds(history: Iterable[Any]) -> dict[str, Any]:
    totals: dict[str, float] = {}
    counts: dict[str, int] = {}
    for rec in history:
        for name, seconds in getattr(rec, "module_seconds", {}).items():
            totals[str(name)] = totals.get(str(name), 0.0) + float(seconds)
            counts[str(name)] = counts.get(str(name), 0) + 1
    top = [
        {"module": name, "seconds": float(seconds), "runs": int(counts.get(name, 0))}
        for name, seconds in sorted(totals.items(), key=lambda item: item[1], reverse=True)
    ]
    return {
        "totals": {name: float(seconds) for name, seconds in sorted(totals.items())},
        "runs": {name: int(count) for name, count in sorted(counts.items())},
        "top": top,
    }


def _scheduler_terrain_internal_profile(history: Iterable[Any]) -> dict[str, Any]:
    stage_totals: dict[str, float] = {}
    subprofile_totals: dict[str, dict[str, float]] = {}
    terrain_module_seconds = 0.0
    profile_total_seconds = 0.0
    profiled_step_count = 0
    terrain_run_count = 0
    latest_time_myr = 0.0
    latest_cell_count = 0
    for rec in history:
        module_seconds = getattr(rec, "module_seconds", {}) or {}
        if "terrain" in module_seconds:
            terrain_run_count += 1
            terrain_module_seconds += float(module_seconds.get("terrain", 0.0))
        diagnostics = getattr(rec, "diagnostics", {}) or {}
        terrain_diag = diagnostics.get("terrain", {}) if isinstance(diagnostics, dict) else {}
        if not isinstance(terrain_diag, dict):
            continue
        profile = terrain_diag.get("terrain_internal_profile", {})
        if not isinstance(profile, dict):
            continue
        stage_seconds = profile.get("stage_seconds", {})
        if not isinstance(stage_seconds, dict) or not stage_seconds:
            continue
        profiled_step_count += 1
        latest_time_myr = float(profile.get("time_myr", latest_time_myr))
        latest_cell_count = int(profile.get("cell_count", latest_cell_count))
        profile_total_seconds += float(profile.get("total_seconds", 0.0))
        for stage, seconds_raw in stage_seconds.items():
            seconds = float(seconds_raw)
            if np.isfinite(seconds) and seconds >= 0.0:
                stage_totals[str(stage)] = stage_totals.get(str(stage), 0.0) + seconds
        subprofiles = profile.get("subprofiles", {})
        if isinstance(subprofiles, dict):
            for profile_name, profile_raw in subprofiles.items():
                if not isinstance(profile_raw, dict):
                    continue
                sub_stage_seconds = profile_raw.get("stage_seconds", {})
                if not isinstance(sub_stage_seconds, dict):
                    continue
                target = subprofile_totals.setdefault(str(profile_name), {})
                for stage, seconds_raw in sub_stage_seconds.items():
                    seconds = float(seconds_raw)
                    if np.isfinite(seconds) and seconds >= 0.0:
                        target[str(stage)] = target.get(str(stage), 0.0) + seconds

    profiled_stage_total = float(sum(stage_totals.values()))
    share_total = max(profiled_stage_total, 1.0e-12)
    top = [
        {
            "stage": stage,
            "seconds": float(seconds),
            "share": float(seconds / share_total),
        }
        for stage, seconds in sorted(
            stage_totals.items(),
            key=lambda item: item[1],
            reverse=True,
        )[:8]
    ]
    coverage = (
        profiled_stage_total / max(terrain_module_seconds, 1.0e-12)
        if terrain_module_seconds > 0.0
        else 0.0
    )
    subprofile_summaries: dict[str, Any] = {}
    for profile_name, stage_map in sorted(subprofile_totals.items()):
        sub_total = max(float(sum(stage_map.values())), 1.0e-12)
        subprofile_summaries[profile_name] = {
            "profiled_stage_seconds_total": float(sum(stage_map.values())),
            "stage_count": int(len(stage_map)),
            "stage_seconds": {
                key: float(value)
                for key, value in sorted(stage_map.items())
            },
            "top": [
                {
                    "stage": stage,
                    "seconds": float(seconds),
                    "share": float(seconds / sub_total),
                }
                for stage, seconds in sorted(
                    stage_map.items(),
                    key=lambda item: item[1],
                    reverse=True,
                )[:8]
            ],
        }
    return {
        "schema": "aevum.terrain_internal_profile.aggregate.v1",
        "available": bool(stage_totals),
        "profiled_step_count": int(profiled_step_count),
        "terrain_run_count": int(terrain_run_count),
        "time_myr": float(latest_time_myr),
        "cell_count": int(latest_cell_count),
        "terrain_module_seconds": float(terrain_module_seconds),
        "profile_total_seconds": float(profile_total_seconds),
        "profiled_stage_seconds_total": float(profiled_stage_total),
        "profile_coverage_fraction_of_terrain_module_seconds": float(coverage),
        "stage_count": int(len(stage_totals)),
        "stage_seconds": {
            key: float(value)
            for key, value in sorted(stage_totals.items())
        },
        "top": top,
        "subprofiles": subprofile_summaries,
    }


def _audit_arrays(world: Any) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    grid = world.grid
    arrays: dict[str, np.ndarray] = {
        "grid_lat": np.asarray(grid.lat, dtype=np.float32),
        "grid_lon": np.asarray(grid.lon, dtype=np.float32),
        "grid_cell_area_m2": np.asarray(grid.cell_area, dtype=np.float64),
    }
    manifest: dict[str, Any] = {
        "fields": {},
        "boundaries": {},
        "object_masks": {},
    }
    for field_name in P107_FIELD_KEYS:
        default = -1.0 if _compact_dtype_name(field_name) is np.int32 else np.nan
        arr = np.asarray(world.get_field(field_name, default), dtype=np.float64)
        key = _array_key("field", field_name)
        arrays[key] = arr.astype(_compact_dtype_name(field_name))
        manifest["fields"][field_name] = key

    boundaries = world.networks.get("tectonics.boundaries", {})
    for kind in P107_BOUNDARY_KEYS:
        mask = np.zeros(grid.n, dtype=np.uint8)
        cells = np.asarray(boundaries.get(kind, []), dtype=np.int64)
        cells = cells[(0 <= cells) & (cells < grid.n)]
        if cells.size:
            mask[np.unique(cells)] = 1
        key = _array_key("boundary", kind)
        arrays[key] = mask
        manifest["boundaries"][kind] = key

    for label, object_set, kinds in P107_OBJECT_MASK_SPECS:
        mask = _object_mask(world, object_set, kinds)
        key = _array_key("object", label)
        arrays[key] = mask.astype(np.uint8)
        manifest["object_masks"][label] = {
            "key": key,
            "object_set": object_set,
            "kinds": sorted(kinds),
        }
    return arrays, manifest


def _compact_dtype_name(field_name: str):
    if field_name.endswith(".plate_id") or field_name in {
        "tectonics.plate_rank",
        "tectonics.boundary_province_kind",
        "crust.type",
        "ocean.depth_province",
        "ocean.basin_id",
        "ocean.margin_type",
        "ocean.gateway_id",
        "ocean.gateway_system_id",
        "terrain.province",
        "terrain.continental_detail",
        "terrain.mountain_ranges",
        "terrain.mountain_inventory",
        "terrain.mountain_hierarchy_level",
        "terrain.orogenic_parent_hierarchy",
        "terrain.orogenic_hierarchy_spine",
        "terrain.orogenic_shoulder_halo",
        "terrain.orogenic_highland_apron",
        "tectonics.mountain_belt_id",
        "tectonics.mountain_parent_process_id",
        "tectonics.continent_id",
        "tectonics.terrane_id",
        "tectonics.continental_province_code",
        "tectonics.internal_geographic_block_code",
    }:
        return np.int32
    return np.float32


def _array_key(prefix: str, name: str) -> str:
    return f"{prefix}__{name.replace('.', '_').replace('-', '_')}"


def _boundary_counts(world: Any) -> dict[str, int]:
    grid = world.grid
    boundaries = world.networks.get("tectonics.boundaries", {})
    out: dict[str, int] = {}
    for kind in P107_BOUNDARY_KEYS:
        cells = np.asarray(boundaries.get(kind, []), dtype=np.int64)
        cells = cells[(0 <= cells) & (cells < grid.n)]
        out[kind] = int(np.unique(cells).size)
    return out


def _object_counts(world: Any) -> dict[str, int]:
    return {
        object_set: int(len(world.objects.get(object_set, [])))
        for object_set in P107_OBJECT_SETS
    }


def _object_kind_counts(world: Any) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = {}
    for object_set in P107_OBJECT_SETS:
        counts: dict[str, int] = {}
        for obj in world.objects.get(object_set, []):
            kind = str(obj.get("kind", obj.get("type", "unknown")))
            counts[kind] = counts.get(kind, 0) + 1
        out[object_set] = counts
    return out


def _object_mask(world: Any, object_set: str, kinds: set[str]) -> np.ndarray:
    grid = world.grid
    mask = np.zeros(grid.n, dtype=bool)
    for obj in world.objects.get(object_set, []):
        if str(obj.get("kind", obj.get("type", ""))) not in kinds:
            continue
        cells = np.asarray(obj.get("cells", []), dtype=np.int64)
        cells = cells[(0 <= cells) & (cells < grid.n)]
        if cells.size:
            mask[cells] = True
        cell = obj.get("cell")
        if cell is not None and 0 <= int(cell) < grid.n:
            mask[int(cell)] = True
    return mask


def _protected_plate_ids(world: Any) -> set[int]:
    out: set[int] = set()
    field = world.fields.get("tectonics.protected_plate_id")
    if field is not None:
        arr = np.asarray(field, dtype=np.int64)
        out.update(int(x) for x in np.unique(arr) if int(x) >= 0)
    for object_set in ("tectonics.microplates", "tectonics.protected_microplates"):
        for obj in world.objects.get(object_set, []):
            for key in ("plate_id", "protected_plate_id"):
                if key in obj and obj[key] is not None:
                    out.add(int(obj[key]))
            for pid in obj.get("plate_ids", []) or []:
                out.add(int(pid))
    return out


def _merged_microplate_count(events: Iterable[Any]) -> int:
    count = 0
    for event in events:
        params = getattr(event, "params", {}) or {}
        if not isinstance(params, dict):
            continue
        for merge in params.get("merged", []) or []:
            frac = float(merge.get("area_fraction", 1.0))
            if frac < MINOR_PLATE_AREA_FRACTION:
                count += 1
    return int(count)


def _ridge_transform_network_metrics(world: Any) -> dict[str, Any]:
    grid = world.grid
    ridge = _boundary_mask(world, "ridge") | _object_mask(
        world, "terrain.ocean_fabric", {"spreading_center"})
    transform = _boundary_mask(world, "transform") | _object_mask(
        world, "terrain.ocean_fabric", {"transform_fault", "fracture_zone"})
    ridge_components = _component_sizes(grid, ridge)
    transform_components = _component_sizes(grid, transform)
    ridge_area = float(grid.cell_area[ridge].sum())
    largest_ridge = max(ridge_components, default=0.0)
    return {
        "ridge_cell_count": int(ridge.sum()),
        "transform_cell_count": int(transform.sum()),
        "ridge_component_count": int(len(ridge_components)),
        "transform_component_count": int(len(transform_components)),
        "ridge_largest_component_share": float(largest_ridge / max(ridge_area, 1.0e-12)),
        "transform_to_ridge_cell_ratio": float(
            transform.sum() / max(int(ridge.sum()), 1)),
        "has_ridge_transform_network": bool(ridge.any() and transform.any()),
    }


def _boundary_width_diagnostics(world: Any) -> dict[str, Any]:
    grid = world.grid
    total = max(float(grid.cell_area.sum()), 1.0e-12)
    masks = {
        "boundary_ridge": _boundary_mask(world, "ridge"),
        "boundary_trench": _boundary_mask(world, "trench"),
        "boundary_transform": _boundary_mask(world, "transform"),
        "province_mid_ocean_ridge": _object_mask(
            world, "tectonics.boundary_provinces", {"mid_ocean_ridge"}),
        "province_trench": (
            _object_mask(world, "tectonics.boundary_provinces", {"ocean_ocean_subduction_trench"})
            | _object_mask(world, "tectonics.boundary_provinces", {"island_arc_trench"})
        ),
        "ocean_fabric_spreading_center": _object_mask(
            world, "terrain.ocean_fabric", {"spreading_center"}),
        "ocean_fabric_transform_fracture": _object_mask(
            world, "terrain.ocean_fabric", {"transform_fault", "fracture_zone"}),
    }
    diagnostics = {
        name: _mask_width_stats(grid, mask, total)
        for name, mask in masks.items()
    }
    ridge_area = diagnostics["boundary_ridge"]["area_fraction_world"]
    trench_area = diagnostics["boundary_trench"]["area_fraction_world"]
    transform_area = diagnostics["boundary_transform"]["area_fraction_world"]
    diagnostics["summary"] = {
        "ridge_trench_transform_area_fraction_world": float(
            ridge_area + trench_area + transform_area),
        "max_p90_width_steps": float(max(
            (diagnostics[name]["p90_width_steps"] for name in diagnostics),
            default=0.0,
        )),
        "max_fraction_width_gt2": float(max(
            (diagnostics[name]["fraction_width_gt2"] for name in diagnostics),
            default=0.0,
        )),
    }
    return diagnostics


def _mask_width_stats(grid: Any, mask: np.ndarray, total_area: float) -> dict[str, Any]:
    mask = np.asarray(mask, dtype=bool)
    area = np.asarray(grid.cell_area, dtype=np.float64)
    if not mask.any():
        return {
            "cell_count": 0,
            "area_fraction_world": 0.0,
            "component_count": 0,
            "median_width_steps": 0.0,
            "p90_width_steps": 0.0,
            "max_width_steps": 0.0,
            "fraction_width_gt1": 0.0,
            "fraction_width_gt2": 0.0,
        }
    width = _mask_interior_width_steps(grid, mask)
    values = width[mask]
    weights = area[mask]
    mask_area = max(float(weights.sum()), 1.0e-12)
    return {
        "cell_count": int(mask.sum()),
        "area_fraction_world": float(mask_area / max(total_area, 1.0e-12)),
        "component_count": int(len(_components(grid, mask))),
        "median_width_steps": float(np.percentile(values, 50)),
        "p90_width_steps": float(np.percentile(values, 90)),
        "max_width_steps": float(np.max(values)),
        "fraction_width_gt1": float(area[mask & (width > 1.0)].sum() / mask_area),
        "fraction_width_gt2": float(area[mask & (width > 2.0)].sum() / mask_area),
    }


def _mask_interior_width_steps(grid: Any, mask: np.ndarray) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    width = np.zeros(grid.n, dtype=np.float64)
    if not mask.any():
        return width
    boundary = np.zeros(grid.n, dtype=bool)
    for c in np.where(mask)[0]:
        if np.any(~mask[grid.neighbors[int(c)]]):
            boundary[int(c)] = True
    if not boundary.any():
        width[mask] = 1.0
        return width
    queue = [int(c) for c in np.where(boundary)[0]]
    seen = np.zeros(grid.n, dtype=bool)
    seen[queue] = True
    width[queue] = 1.0
    head = 0
    while head < len(queue):
        c = queue[head]
        head += 1
        for nb in grid.neighbors[c]:
            nb = int(nb)
            if not mask[nb] or seen[nb]:
                continue
            seen[nb] = True
            width[nb] = width[c] + 1.0
            queue.append(nb)
    return width


def _collision_plateau_area_fraction(world: Any, rel: np.ndarray,
                                     land: np.ndarray) -> float:
    grid = world.grid
    total = max(float(grid.cell_area.sum()), 1.0e-12)
    collision = (
        _boundary_mask(world, "collision")
        | _boundary_mask(world, "suture")
        | _object_mask(world, "tectonics.boundary_objects", {"suture"})
    )
    if collision.any():
        collision = _dilate(grid, collision, passes=4)
    terrain_province = np.asarray(
        world.get_field("terrain.province", -1.0), dtype=np.int64)
    plateau = land & (rel >= 2000.0) & (
        collision | np.isin(terrain_province, [5, 7])
    )
    return float(grid.cell_area[plateau].sum() / total)


def _high_mountain_coherence_metrics(world: Any, rel: np.ndarray,
                                     land: np.ndarray) -> dict[str, Any]:
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    crust_type = np.asarray(world.get_field("crust.type", 0.0), dtype=np.float64)
    continental = crust_type >= 0.5
    high = land & continental & (rel >= 3000.0)
    extreme = land & continental & (rel >= 4500.0)
    parent = _high_mountain_parent_mask(world)
    near_parent = _dilate(grid, parent, passes=3) if parent.any() else parent

    high_components = _component_area_rows(grid, high)
    extreme_components = _component_area_rows(grid, extreme)
    high_area = max(float(area[high].sum()), 1.0e-12)
    extreme_area = max(float(area[extreme].sum()), 1.0e-12)
    largest_high = max((row["area"] for row in high_components), default=0.0)
    mountain_ranges = np.asarray(
        world.get_field("terrain.mountain_ranges", np.zeros(grid.n)),
        dtype=np.float64,
    )
    if mountain_ranges.shape != (grid.n,):
        mountain_ranges = np.zeros(grid.n, dtype=np.float64)
    mountain_mask = mountain_ranges > 0.0
    mountain_near = _dilate(grid, mountain_mask, passes=2) if mountain_mask.any() else mountain_mask
    small_high_component_count = 0
    small_unparented_high_component_count = 0
    small_high_component_area = 0.0
    isolated_extreme_area = 0.0
    for comp in _components(grid, high):
        comp_area = float(area[comp].sum())
        if comp.size <= 8:
            small_high_component_count += 1
            small_high_component_area += comp_area
            support = float(
                area[comp[mountain_near[comp] | near_parent[comp]]].sum()
                / max(comp_area, 1.0e-12))
            if support < 0.25:
                small_unparented_high_component_count += 1
    for comp in _components(grid, extreme):
        comp_area = float(area[comp].sum())
        support_cells = comp[near_parent[comp]]
        if comp.size <= 6 and float(area[support_cells].sum() / max(comp_area, 1.0e-12)) < 0.25:
            isolated_extreme_area += comp_area
    parent_overlap = float(area[high & near_parent].sum() / high_area) if high.any() else 0.0
    top3_high_area = sum(float(row["area"]) for row in high_components[:3])
    top3_high_share = (
        float(np.clip(top3_high_area / high_area, 0.0, 1.0))
        if high.any() else 0.0
    )
    small_high_area_share = (
        float(small_high_component_area / high_area) if high.any() else 0.0
    )
    fragmentation_index = (
        float(np.clip(1.0 - top3_high_share + 0.5 * small_high_area_share, 0.0, 1.0))
        if high.any() else 0.0
    )
    return {
        "high_mountain_threshold_m": 3000.0,
        "extreme_mountain_threshold_m": 4500.0,
        "high_mountain_cell_count": int(high.sum()),
        "extreme_mountain_cell_count": int(extreme.sum()),
        "high_mountain_area_fraction_world": float(area[high].sum() / total),
        "extreme_mountain_area_fraction_world": float(area[extreme].sum() / total),
        "high_mountain_component_count": int(len(high_components)),
        "extreme_mountain_component_count": int(len(extreme_components)),
        "largest_high_mountain_component_share": float(largest_high / high_area),
        "top3_high_mountain_component_share": top3_high_share,
        "small_high_mountain_component_area_fraction_of_high": small_high_area_share,
        "high_mountain_fragmentation_index": fragmentation_index,
        "isolated_extreme_peak_area_fraction_of_extreme": float(
            isolated_extreme_area / extreme_area if extreme.any() else 0.0),
        "high_mountain_parent_overlap_fraction": parent_overlap,
        "high_mountain_mountain_range_overlap_fraction": float(
            area[high & mountain_near].sum() / high_area) if high.any() else 0.0,
        "small_high_mountain_component_count": int(small_high_component_count),
        "small_unparented_high_mountain_component_count": int(
            small_unparented_high_component_count),
        "top_high_mountain_component_area_fractions_world": [
            float(row["area"] / total) for row in high_components[:10]
        ],
        "p1113_bridge_area_fraction": float(
            world.g("terrain.last_p1113_bridge_area_fraction", 0.0)),
        "p1113_high_component_count_before": float(
            world.g("terrain.last_p1113_high_component_count_before", 0.0)),
        "p1113_high_component_count_after": float(
            world.g("terrain.last_p1113_high_component_count_after", 0.0)),
        "p1113_top3_high_share_before": float(
            world.g("terrain.last_p1113_top3_high_share_before", 0.0)),
        "p1113_top3_high_share_after": float(
            world.g("terrain.last_p1113_top3_high_share_after", 0.0)),
        "p1113_fragmentation_index_before": float(
            world.g("terrain.last_p1113_fragmentation_index_before", 0.0)),
        "p1113_fragmentation_index_after": float(
            world.g("terrain.last_p1113_fragmentation_index_after", 0.0)),
        "p1113_range_count_considered": float(
            world.g("terrain.last_p1113_range_count_considered", 0.0)),
        "p1113_land_mask_preserved": bool(
            world.g("terrain.last_p1113_land_mask_preserved", 0.0) >= 1.0),
        "p11169_hierarchy_used": bool(
            world.g("terrain.last_p11169_hierarchy_used", 0.0) >= 1.0),
        "p11169_hierarchy_bridge_area_fraction": float(
            world.g("terrain.last_p11169_hierarchy_bridge_area_fraction", 0.0)),
        "p11169_hierarchy_bridge_cell_count": float(
            world.g("terrain.last_p11169_hierarchy_bridge_cell_count", 0.0)),
        "p11169_legacy_bridge_area_fraction": float(
            world.g("terrain.last_p11169_legacy_bridge_area_fraction", 0.0)),
        "p11169_foreland_peak_cell_count": float(
            world.g("terrain.last_p11169_foreland_peak_cell_count", 0.0)),
        "p11169_fragmentation_delta": float(
            world.g("terrain.last_p11169_fragmentation_delta", 0.0)),
        "p11169_crest_high_overlap_fraction": float(
            world.g("terrain.last_p11169_crest_high_overlap_fraction", 0.0)),
        "p11169_bridge_restricted_to_peak_hierarchy": bool(
            world.g(
                "terrain.last_p11169_bridge_restricted_to_peak_hierarchy",
                0.0,
            ) >= 1.0
        ),
        "p11169_hypsometry_guard_reverted": bool(
            world.g("terrain.last_p11169_hypsometry_guard_reverted", 0.0) >= 1.0
        ),
        "p153_terminal_high_fleck_cleanup_accepted": bool(
            world.g(
                "terrain.last_p153_terminal_high_fleck_cleanup_accepted",
                0.0,
            ) >= 1.0
        ),
        "p153_candidate_cell_count": float(
            world.g("terrain.last_p153_candidate_cell_count", 0.0)),
        "p153_candidate_area_fraction": float(
            world.g("terrain.last_p153_candidate_area_fraction", 0.0)),
        "p153_softened_cell_count": float(
            world.g("terrain.last_p153_softened_cell_count", 0.0)),
        "p153_softened_area_fraction": float(
            world.g("terrain.last_p153_softened_area_fraction", 0.0)),
        "p153_high_component_count_before": float(
            world.g("terrain.last_p153_high_component_count_before", 0.0)),
        "p153_high_component_count_after": float(
            world.g("terrain.last_p153_high_component_count_after", 0.0)),
        "p153_top3_high_share_before": float(
            world.g("terrain.last_p153_top3_high_share_before", 0.0)),
        "p153_top3_high_share_after": float(
            world.g("terrain.last_p153_top3_high_share_after", 0.0)),
        "p153_fragmentation_index_before": float(
            world.g("terrain.last_p153_fragmentation_index_before", 0.0)),
        "p153_fragmentation_index_after": float(
            world.g("terrain.last_p153_fragmentation_index_after", 0.0)),
        "p153_guard_reverted": bool(
            world.g("terrain.last_p153_guard_reverted", 0.0) >= 1.0
        ),
        "p155_terminal_high_relief_consistency_gate_accepted": bool(
            world.g(
                "terrain.last_p155_terminal_high_relief_consistency_gate_accepted",
                0.0,
            ) >= 1.0
        ),
        "p155_guard_reverted": bool(
            world.g("terrain.last_p155_guard_reverted", 0.0) >= 1.0
        ),
        "p155_land_mask_preserved": bool(
            world.g("terrain.last_p155_land_mask_preserved", 0.0) >= 1.0
        ),
        "p155_candidate_cell_count": float(
            world.g("terrain.last_p155_candidate_cell_count", 0.0)),
        "p155_candidate_area_fraction": float(
            world.g("terrain.last_p155_candidate_area_fraction", 0.0)),
        "p155_raised_cell_count": float(
            world.g("terrain.last_p155_raised_cell_count", 0.0)),
        "p155_raised_area_fraction": float(
            world.g("terrain.last_p155_raised_area_fraction", 0.0)),
        "p155_softened_cell_count": float(
            world.g("terrain.last_p155_softened_cell_count", 0.0)),
        "p155_softened_area_fraction": float(
            world.g("terrain.last_p155_softened_area_fraction", 0.0)),
        "p155_high_component_count_before": float(
            world.g("terrain.last_p155_high_component_count_before", 0.0)),
        "p155_high_component_count_after": float(
            world.g("terrain.last_p155_high_component_count_after", 0.0)),
        "p155_spine_3000_component_count_before": float(
            world.g("terrain.last_p155_spine_3000_component_count_before", 0.0)),
        "p155_spine_3000_component_count_after": float(
            world.g("terrain.last_p155_spine_3000_component_count_after", 0.0)),
        "p155_spine_3000_coverage_before": float(
            world.g("terrain.last_p155_spine_3000_coverage_before", 0.0)),
        "p155_spine_3000_coverage_after": float(
            world.g("terrain.last_p155_spine_3000_coverage_after", 0.0)),
        "p155_top3_high_share_before": float(
            world.g("terrain.last_p155_top3_high_share_before", 0.0)),
        "p155_top3_high_share_after": float(
            world.g("terrain.last_p155_top3_high_share_after", 0.0)),
        "p155_fragmentation_index_before": float(
            world.g("terrain.last_p155_fragmentation_index_before", 0.0)),
        "p155_fragmentation_index_after": float(
            world.g("terrain.last_p155_fragmentation_index_after", 0.0)),
        "p155_parent_overlap_before": float(
            world.g("terrain.last_p155_parent_overlap_before", 0.0)),
        "p155_parent_overlap_after": float(
            world.g("terrain.last_p155_parent_overlap_after", 0.0)),
        "p155_high_area_fraction_before": float(
            world.g("terrain.last_p155_high_area_fraction_before", 0.0)),
        "p155_high_area_fraction_after": float(
            world.g("terrain.last_p155_high_area_fraction_after", 0.0)),
        "p155_p90_land_relief_before_m": float(
            world.g("terrain.last_p155_p90_land_relief_before_m", 0.0)),
        "p155_p90_land_relief_after_m": float(
            world.g("terrain.last_p155_p90_land_relief_after_m", 0.0)),
        "p155_p98_land_relief_before_m": float(
            world.g("terrain.last_p155_p98_land_relief_before_m", 0.0)),
        "p155_p98_land_relief_after_m": float(
            world.g("terrain.last_p155_p98_land_relief_after_m", 0.0)),
        "p155_reject_code": float(
            world.g("terrain.last_p155_reject_code", 0.0)),
        "p162_terminal_high_mountain_fragment_cleanup_accepted": bool(
            world.g(
                "terrain.last_p162_terminal_high_mountain_fragment_cleanup_accepted",
                0.0,
            ) >= 1.0
        ),
        "p162_guard_reverted": bool(
            world.g("terrain.last_p162_guard_reverted", 0.0) >= 1.0
        ),
        "p162_land_mask_preserved": bool(
            world.g("terrain.last_p162_land_mask_preserved", 0.0) >= 1.0
        ),
        "p162_candidate_cell_count": float(
            world.g("terrain.last_p162_candidate_cell_count", 0.0)),
        "p162_candidate_area_fraction": float(
            world.g("terrain.last_p162_candidate_area_fraction", 0.0)),
        "p162_softened_cell_count": float(
            world.g("terrain.last_p162_softened_cell_count", 0.0)),
        "p162_softened_area_fraction": float(
            world.g("terrain.last_p162_softened_area_fraction", 0.0)),
        "p162_high_component_count_before": float(
            world.g("terrain.last_p162_high_component_count_before", 0.0)),
        "p162_high_component_count_after": float(
            world.g("terrain.last_p162_high_component_count_after", 0.0)),
        "p162_small_high_component_count_before": float(
            world.g("terrain.last_p162_small_high_component_count_before", 0.0)),
        "p162_small_high_component_count_after": float(
            world.g("terrain.last_p162_small_high_component_count_after", 0.0)),
        "p162_spine_3000_component_count_before": float(
            world.g("terrain.last_p162_spine_3000_component_count_before", 0.0)),
        "p162_spine_3000_component_count_after": float(
            world.g("terrain.last_p162_spine_3000_component_count_after", 0.0)),
        "p162_top3_high_share_before": float(
            world.g("terrain.last_p162_top3_high_share_before", 0.0)),
        "p162_top3_high_share_after": float(
            world.g("terrain.last_p162_top3_high_share_after", 0.0)),
        "p162_fragmentation_index_before": float(
            world.g("terrain.last_p162_fragmentation_index_before", 0.0)),
        "p162_fragmentation_index_after": float(
            world.g("terrain.last_p162_fragmentation_index_after", 0.0)),
        "p162_high_area_fraction_before": float(
            world.g("terrain.last_p162_high_area_fraction_before", 0.0)),
        "p162_high_area_fraction_after": float(
            world.g("terrain.last_p162_high_area_fraction_after", 0.0)),
        "p162_p90_land_relief_before_m": float(
            world.g("terrain.last_p162_p90_land_relief_before_m", 0.0)),
        "p162_p90_land_relief_after_m": float(
            world.g("terrain.last_p162_p90_land_relief_after_m", 0.0)),
        "p162_p98_land_relief_before_m": float(
            world.g("terrain.last_p162_p98_land_relief_before_m", 0.0)),
        "p162_p98_land_relief_after_m": float(
            world.g("terrain.last_p162_p98_land_relief_after_m", 0.0)),
        "p162_mean_lower_m": float(
            world.g("terrain.last_p162_mean_lower_m", 0.0)),
        "p162_max_lower_m": float(
            world.g("terrain.last_p162_max_lower_m", 0.0)),
        "p162_extreme_softened_cell_count": float(
            world.g("terrain.last_p162_extreme_softened_cell_count", 0.0)),
        "p162_reject_code": float(
            world.g("terrain.last_p162_reject_code", 0.0)),
        "p156_path_candidate_count": float(
            world.g("terrain.last_p155_p156_path_candidate_count", 0.0)),
        "p156_path_candidate_cell_count": float(
            world.g("terrain.last_p155_p156_path_candidate_cell_count", 0.0)),
        "p156_selected_path_count": float(
            world.g("terrain.last_p155_p156_selected_path_count", 0.0)),
        "p156_selected_path_cell_count": float(
            world.g("terrain.last_p155_p156_selected_path_cell_count", 0.0)),
        "p157_seed_path_candidate_count": float(
            world.g("terrain.last_p155_p157_seed_path_candidate_count", 0.0)),
        "p157_seed_path_candidate_cell_count": float(
            world.g("terrain.last_p155_p157_seed_path_candidate_cell_count", 0.0)),
        "p157_selected_seed_path_count": float(
            world.g("terrain.last_p155_p157_selected_seed_path_count", 0.0)),
        "p157_selected_seed_path_cell_count": float(
            world.g("terrain.last_p155_p157_selected_seed_path_cell_count", 0.0)),
        "p111611_bridge_candidate_cell_count": float(
            world.g("terrain.last_p111611_bridge_candidate_cell_count", 0.0)),
        "p111611_bridge_candidate_area_fraction": float(
            world.g("terrain.last_p111611_bridge_candidate_area_fraction", 0.0)),
        "p111611_peak_hierarchy_shoulder_cell_count": float(
            world.g("terrain.last_p111611_peak_hierarchy_shoulder_cell_count", 0.0)),
        "p111611_high_pair_count": float(
            world.g("terrain.last_p111611_high_pair_count", 0.0)),
        "p111611_safe_path_count": float(
            world.g("terrain.last_p111611_safe_path_count", 0.0)),
        "p111611_blocked_high_pair_count": float(
            world.g("terrain.last_p111611_blocked_high_pair_count", 0.0)),
        "p111611_diagnostic_system_count": float(
            world.g("terrain.last_p111611_diagnostic_system_count", 0.0)),
        "p111611_bridge_deferred_no_safe_path": bool(
            world.g("terrain.last_p111611_bridge_deferred_no_safe_path", 0.0) >= 1.0
        ),
        "p111622_spine_guided_repair_used": bool(
            world.g("terrain.last_p111622_spine_guided_repair_used", 0.0) >= 1.0
        ),
        "p111622_spine_guided_candidate_cell_count": float(
            world.g("terrain.last_p111622_spine_guided_candidate_cell_count", 0.0)),
        "p111622_spine_guided_candidate_area_fraction": float(
            world.g("terrain.last_p111622_spine_guided_candidate_area_fraction", 0.0)),
        "p111622_spine_guided_bridge_cell_count": float(
            world.g("terrain.last_p111622_spine_guided_bridge_cell_count", 0.0)),
        "p111622_spine_guided_bridge_area_fraction": float(
            world.g("terrain.last_p111622_spine_guided_bridge_area_fraction", 0.0)),
        "p111622_high_component_delta": float(
            world.g("terrain.last_p111622_high_component_delta", 0.0)),
        "p111622_top3_high_share_delta": float(
            world.g("terrain.last_p111622_top3_high_share_delta", 0.0)),
        "p111622_fragmentation_delta": float(
            world.g("terrain.last_p111622_fragmentation_delta", 0.0)),
        "p111624_spine_to_terrain_response_used": bool(
            world.g(
                "terrain.last_p111624_spine_to_terrain_response_used",
                0.0,
            ) >= 1.0
        ),
        "p111624_spine_response_candidate_cell_count": float(
            world.g(
                "terrain.last_p111624_spine_response_candidate_cell_count",
                0.0,
            )
        ),
        "p111624_spine_response_candidate_area_fraction": float(
            world.g(
                "terrain.last_p111624_spine_response_candidate_area_fraction",
                0.0,
            )
        ),
        "p111624_spine_response_bridge_cell_count": float(
            world.g(
                "terrain.last_p111624_spine_response_bridge_cell_count",
                0.0,
            )
        ),
        "p111624_spine_response_bridge_area_fraction": float(
            world.g(
                "terrain.last_p111624_spine_response_bridge_area_fraction",
                0.0,
            )
        ),
        "p111624_spine_saddle_response_cell_count": float(
            world.g(
                "terrain.last_p111624_spine_saddle_response_cell_count",
                0.0,
            )
        ),
        "p111624_spine_saddle_response_area_fraction": float(
            world.g(
                "terrain.last_p111624_spine_saddle_response_area_fraction",
                0.0,
            )
        ),
        "p111624_spine_saddle_guard_reverted": bool(
            world.g(
                "terrain.last_p111624_spine_saddle_guard_reverted",
                0.0,
            ) >= 1.0
        ),
        "p111624_high_component_delta": float(
            world.g("terrain.last_p111624_high_component_delta", 0.0)),
        "p111624_top3_high_share_delta": float(
            world.g("terrain.last_p111624_top3_high_share_delta", 0.0)),
        "p111624_fragmentation_delta": float(
            world.g("terrain.last_p111624_fragmentation_delta", 0.0)),
        "p111624_guard_accepted_small_spine_response": bool(
            world.g(
                "terrain.last_p111624_guard_accepted_small_spine_response",
                0.0,
            ) >= 1.0
        ),
    }


def _orogenic_belt_morphology_metrics(
    world: Any,
    rel: np.ndarray,
    land: np.ndarray,
    hierarchy: np.ndarray,
    spine: np.ndarray,
    shoulder_halo: np.ndarray,
    highland_apron: np.ndarray,
) -> dict[str, Any]:
    """Measure continuous mountain-belt expression below the peak mask.

    High mountain pixels can be naturally discontinuous.  This diagnostic
    separates peak fragmentation from the more important question of whether
    the parent orogenic belt is expressed as a connected ridge/branch/shoulder
    system across lower relief tiers.
    """
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    rel = np.asarray(rel, dtype=np.float64)
    land = np.asarray(land, dtype=bool)
    hierarchy = np.asarray(hierarchy, dtype=np.int64)
    spine = np.asarray(spine, dtype=np.int64)
    shoulder_halo = np.asarray(shoulder_halo, dtype=np.int64)
    highland_apron = np.asarray(highland_apron, dtype=np.int64)
    crust_type = np.asarray(world.get_field("crust.type", 0.0), dtype=np.float64)
    if crust_type.shape != (grid.n,):
        crust_type = np.zeros(grid.n, dtype=np.float64)
    continental_land = land & (crust_type >= 0.5)

    crest = continental_land & (hierarchy >= 3)
    branch = continental_land & (hierarchy == 2)
    foreland = continental_land & (hierarchy == 1)
    peak_hierarchy = crest | branch
    spine_crest = continental_land & (spine >= 3)
    spine_branch = continental_land & (spine == 2)
    spine_peak = spine_crest | spine_branch
    halo = continental_land & (shoulder_halo > 0)
    apron = continental_land & (highland_apron > 0)

    def stats(mask: np.ndarray) -> dict[str, Any]:
        mask = np.asarray(mask, dtype=bool)
        rows = _component_area_rows(grid, mask)
        mask_area = float(area[mask].sum())
        if mask_area <= 0.0:
            return {
                "cell_count": 0,
                "area_fraction_world": 0.0,
                "component_count": 0,
                "largest_component_share": 0.0,
                "top3_component_share": 0.0,
                "small_component_area_fraction": 0.0,
                "fragmentation_index": 0.0,
            }
        largest = max((row["area"] for row in rows), default=0.0)
        top3 = float(np.clip(
            sum(float(row["area"]) for row in rows[:3]) / mask_area,
            0.0,
            1.0,
        ))
        small_area = 0.0
        for comp in _components(grid, mask):
            if comp.size <= 8:
                small_area += float(area[comp].sum())
        small_share = float(small_area / mask_area)
        fragmentation = float(np.clip(1.0 - top3 + 0.5 * small_share, 0.0, 1.0))
        return {
            "cell_count": int(mask.sum()),
            "area_fraction_world": float(mask_area / total),
            "component_count": int(len(rows)),
            "largest_component_share": float(largest / mask_area),
            "top3_component_share": top3,
            "small_component_area_fraction": small_share,
            "fragmentation_index": fragmentation,
        }

    peak_stats = stats(peak_hierarchy)
    spine_stats = stats(spine_peak)
    crest_stats = stats(crest)
    branch_stats = stats(branch)
    halo_stats = stats(halo)
    apron_stats = stats(apron)
    peak_area = max(float(area[peak_hierarchy].sum()), 1.0e-12)
    spine_area = max(float(area[spine_peak].sum()), 1.0e-12)
    crest_area = max(float(area[crest].sum()), 1.0e-12)
    branch_area = max(float(area[branch].sum()), 1.0e-12)
    peak_component_count = max(int(peak_stats["component_count"]), 1)

    def relief_tier(threshold: float) -> dict[str, Any]:
        tier = peak_hierarchy & (rel >= float(threshold))
        tier_stats = stats(tier)
        tier_area = float(area[tier].sum())
        tier_stats.update({
            "threshold_m": float(threshold),
            "coverage_fraction_of_peak_hierarchy": float(tier_area / peak_area),
            "coverage_fraction_of_spine": (
                float(area[tier & spine_peak].sum() / spine_area)
                if spine_peak.any() else 0.0
            ),
            "spine_overlap_fraction_of_tier": (
                float(area[tier & spine_peak].sum() / tier_area)
                if tier_area > 0.0 else 0.0
            ),
            "crest_coverage_fraction": float(area[tier & crest].sum() / crest_area),
            "branch_coverage_fraction": float(area[tier & branch].sum() / branch_area),
            "components_per_peak_hierarchy_component": float(
                int(tier_stats["component_count"]) / peak_component_count
            ),
        })
        continuity = (
            float(tier_stats["top3_component_share"])
            * (1.0 - 0.5 * float(tier_stats["small_component_area_fraction"]))
        )
        tier_stats["continuity_score"] = float(np.clip(continuity, 0.0, 1.0))
        return tier_stats

    relief_1800 = relief_tier(1800.0)
    relief_2400 = relief_tier(2400.0)
    relief_3000 = relief_tier(3000.0)
    relief_4500 = relief_tier(4500.0)
    graded_score = float(np.clip(
        0.40 * float(relief_1800["continuity_score"])
        + 0.35 * float(relief_2400["continuity_score"])
        + 0.25 * float(relief_3000["continuity_score"]),
        0.0,
        1.0,
    ))

    worst_peak_subcomponents = 0
    worst_peak_relief_3000_cells = 0
    for comp in _components(grid, peak_hierarchy):
        comp_mask = np.zeros(grid.n, dtype=bool)
        comp_mask[comp] = True
        high_parts = _components(grid, comp_mask & (rel >= 3000.0))
        if len(high_parts) > worst_peak_subcomponents:
            worst_peak_subcomponents = len(high_parts)
            worst_peak_relief_3000_cells = int(np.count_nonzero(comp_mask & (rel >= 3000.0)))

    near_high_saddles = peak_hierarchy & (rel >= 2400.0) & (rel < 3000.0)
    return {
        "schema": "aevum.p111633_orogenic_belt_morphology.v1",
        "hierarchy": peak_stats,
        "crest": crest_stats,
        "branch": branch_stats,
        "foreland": stats(foreland),
        "spine": spine_stats,
        "spine_crest": stats(spine_crest),
        "spine_branch": stats(spine_branch),
        "shoulder_halo": halo_stats,
        "highland_apron": apron_stats,
        "relief_1800": relief_1800,
        "relief_2400": relief_2400,
        "relief_3000": relief_3000,
        "relief_4500": relief_4500,
        "graded_belt_continuity_score": graded_score,
        "spine_components_per_hierarchy_component": float(
            int(spine_stats["component_count"]) / peak_component_count
        ),
        "high_peak_fragmentation_pressure": float(
            int(relief_3000["component_count"]) / peak_component_count
        ),
        "worst_peak_hierarchy_high_subcomponent_count": int(worst_peak_subcomponents),
        "worst_peak_hierarchy_high_cell_count": int(worst_peak_relief_3000_cells),
        "near_high_saddle_cell_count": int(np.count_nonzero(near_high_saddles)),
        "near_high_saddle_area_fraction_world": float(
            area[near_high_saddles].sum() / total
        ),
        "high_relief_spine_gap_pressure": float(
            max(
                0.0,
                float(relief_2400["components_per_peak_hierarchy_component"])
                - float(relief_3000["components_per_peak_hierarchy_component"]),
            )
        ),
        "p111633_belt_relief_response_used": bool(
            world.g("terrain.last_p111633_belt_relief_response_used", 0.0) >= 1.0
        ),
        "p111633_belt_relief_candidate_cell_count": float(
            world.g("terrain.last_p111633_belt_relief_candidate_cell_count", 0.0)
        ),
        "p111633_belt_relief_candidate_area_fraction": float(
            world.g("terrain.last_p111633_belt_relief_candidate_area_fraction", 0.0)
        ),
        "p111633_belt_relief_bridge_cell_count": float(
            world.g("terrain.last_p111633_belt_relief_bridge_cell_count", 0.0)
        ),
        "p111633_belt_relief_bridge_area_fraction": float(
            world.g("terrain.last_p111633_belt_relief_bridge_area_fraction", 0.0)
        ),
        "p111633_belt_relief_1800_component_count_before": float(
            world.g(
                "terrain.last_p111633_belt_relief_1800_component_count_before",
                0.0,
            )
        ),
        "p111633_belt_relief_1800_component_count_after": float(
            world.g(
                "terrain.last_p111633_belt_relief_1800_component_count_after",
                0.0,
            )
        ),
        "p111633_belt_relief_2400_component_count_before": float(
            world.g(
                "terrain.last_p111633_belt_relief_2400_component_count_before",
                0.0,
            )
        ),
        "p111633_belt_relief_2400_component_count_after": float(
            world.g(
                "terrain.last_p111633_belt_relief_2400_component_count_after",
                0.0,
            )
        ),
        "p111633_belt_relief_guard_reverted": bool(
            world.g("terrain.last_p111633_belt_relief_guard_reverted", 0.0) >= 1.0
        ),
    }


def _spine_aligned_elevation_morphology_metrics(
    world: Any,
    rel: np.ndarray,
    land: np.ndarray,
    hierarchy: np.ndarray,
    spine: np.ndarray,
    shoulder_halo: np.ndarray,
    highland_apron: np.ndarray,
) -> dict[str, Any]:
    """Measure whether high relief is organized by the orogenic spine graph."""
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    rel = np.asarray(rel, dtype=np.float64)
    land = np.asarray(land, dtype=bool)
    hierarchy = np.asarray(hierarchy, dtype=np.int64)
    spine = np.asarray(spine, dtype=np.int64)
    shoulder_halo = np.asarray(shoulder_halo, dtype=np.int64)
    highland_apron = np.asarray(highland_apron, dtype=np.int64)
    crust_type = np.asarray(world.get_field("crust.type", 0.0), dtype=np.float64)
    if crust_type.shape != (grid.n,):
        crust_type = np.zeros(grid.n, dtype=np.float64)
    continental_land = land & (crust_type >= 0.5)

    crest = continental_land & (hierarchy >= 3)
    branch = continental_land & (hierarchy == 2)
    foreland = continental_land & (hierarchy == 1)
    peak_hierarchy = crest | branch
    spine_crest = continental_land & (spine >= 3)
    spine_branch = continental_land & (spine == 2)
    spine_peak = spine_crest | spine_branch
    halo = continental_land & (shoulder_halo > 0)
    apron = continental_land & (highland_apron > 0)
    high = peak_hierarchy & (rel >= 3000.0)
    spine_near_d1 = peak_hierarchy & _dilate(grid, spine_peak, passes=1)
    spine_near_d2 = peak_hierarchy & _dilate(grid, spine_peak, passes=2)
    high_area = float(area[high].sum())
    peak_component_count = len(_components(grid, peak_hierarchy))
    high_component_count = len(_components(grid, high))
    spine_component_count = len(_components(grid, spine_peak))

    def fraction(mask: np.ndarray, denom_area: float) -> float:
        if denom_area <= 0.0:
            return 0.0
        return float(np.clip(float(area[mask].sum()) / denom_area, 0.0, 1.0))

    def top3_share(mask: np.ndarray) -> float:
        mask_area = float(area[mask].sum())
        if mask_area <= 0.0:
            return 0.0
        rows = _component_area_rows(grid, mask)
        return float(np.clip(
            sum(float(row["area"]) for row in rows[:3]) / mask_area,
            0.0,
            1.0,
        ))

    def median_relief(mask: np.ndarray) -> float:
        mask = np.asarray(mask, dtype=bool)
        if not mask.any():
            return 0.0
        return float(np.median(rel[mask]))

    spine_relief_2400 = spine_peak & (rel >= 2400.0)
    spine_relief_3000 = spine_peak & (rel >= 3000.0)
    spine_area = float(area[spine_peak].sum())
    spine_2400_coverage = fraction(spine_relief_2400, spine_area)
    high_near_d1 = fraction(high & spine_near_d1, high_area)
    high_near_d2 = fraction(high & spine_near_d2, high_area)
    high_far = fraction(high & ~spine_near_d2, high_area)
    ridge_axis_relief_continuity_score = float(np.clip(
        spine_2400_coverage
        * top3_share(spine_relief_2400)
        * (0.5 + 0.5 * high_near_d2 if high_area > 0.0 else 1.0),
        0.0,
        1.0,
    ))

    medians = {
        "crest_spine_median_relief_m": median_relief(spine_crest),
        "branch_spine_median_relief_m": median_relief(spine_branch),
        "shoulder_halo_median_relief_m": median_relief(halo),
        "foreland_median_relief_m": median_relief(foreland),
        "highland_apron_median_relief_m": median_relief(apron),
    }
    comparisons: list[bool] = []
    if spine_crest.any() and spine_branch.any():
        comparisons.append(
            medians["crest_spine_median_relief_m"]
            >= medians["branch_spine_median_relief_m"])
    if spine_branch.any() and halo.any():
        comparisons.append(
            medians["branch_spine_median_relief_m"]
            >= medians["shoulder_halo_median_relief_m"])
    if halo.any() and foreland.any():
        comparisons.append(
            medians["shoulder_halo_median_relief_m"]
            >= medians["foreland_median_relief_m"])
    if spine_crest.any() and apron.any():
        comparisons.append(
            medians["crest_spine_median_relief_m"]
            >= medians["highland_apron_median_relief_m"])
    gradient_order_score = (
        float(sum(1 for value in comparisons if value) / len(comparisons))
        if comparisons else 0.0
    )

    return {
        "schema": "aevum.p111634_spine_aligned_elevation_morphology.v1",
        "peak_hierarchy_cell_count": int(np.count_nonzero(peak_hierarchy)),
        "peak_hierarchy_component_count": int(peak_component_count),
        "spine_cell_count": int(np.count_nonzero(spine_peak)),
        "spine_component_count": int(spine_component_count),
        "spine_relief_2400_cell_count": int(np.count_nonzero(spine_relief_2400)),
        "spine_relief_2400_component_count": int(
            len(_components(grid, spine_relief_2400))),
        "spine_relief_3000_cell_count": int(np.count_nonzero(spine_relief_3000)),
        "spine_relief_3000_component_count": int(
            len(_components(grid, spine_relief_3000))),
        "spine_relief_2400_coverage_fraction": float(spine_2400_coverage),
        "high_cell_count": int(np.count_nonzero(high)),
        "high_area_fraction_world": float(high_area / total),
        "high_component_count": int(high_component_count),
        "high_components_per_peak_hierarchy_component": float(
            high_component_count / max(float(peak_component_count), 1.0)
        ),
        "high_near_spine_fraction_d1": float(high_near_d1),
        "high_near_spine_fraction_d2": float(high_near_d2),
        "high_far_from_spine_fraction": float(high_far),
        "ridge_axis_relief_continuity_score": ridge_axis_relief_continuity_score,
        "gradient_order_score": float(np.clip(gradient_order_score, 0.0, 1.0)),
        **medians,
        "p111634_spine_aligned_response_used": bool(
            world.g("terrain.last_p111634_spine_aligned_response_used", 0.0) >= 1.0
        ),
        "p111634_axis_raise_candidate_cell_count": float(
            world.g("terrain.last_p111634_axis_raise_candidate_cell_count", 0.0)
        ),
        "p111634_axis_raise_candidate_area_fraction": float(
            world.g("terrain.last_p111634_axis_raise_candidate_area_fraction", 0.0)
        ),
        "p111634_axis_raised_cell_count": float(
            world.g("terrain.last_p111634_axis_raised_cell_count", 0.0)
        ),
        "p111634_axis_raised_area_fraction": float(
            world.g("terrain.last_p111634_axis_raised_area_fraction", 0.0)
        ),
        "p111634_offaxis_high_candidate_cell_count": float(
            world.g("terrain.last_p111634_offaxis_high_candidate_cell_count", 0.0)
        ),
        "p111634_offaxis_high_candidate_area_fraction": float(
            world.g("terrain.last_p111634_offaxis_high_candidate_area_fraction", 0.0)
        ),
        "p111634_offaxis_high_softened_cell_count": float(
            world.g("terrain.last_p111634_offaxis_high_softened_cell_count", 0.0)
        ),
        "p111634_offaxis_high_softened_area_fraction": float(
            world.g("terrain.last_p111634_offaxis_high_softened_area_fraction", 0.0)
        ),
        "p111634_high_component_count_before": float(
            world.g("terrain.last_p111634_high_component_count_before", 0.0)
        ),
        "p111634_high_component_count_after": float(
            world.g("terrain.last_p111634_high_component_count_after", 0.0)
        ),
        "p111634_high_near_spine_fraction_before": float(
            world.g("terrain.last_p111634_high_near_spine_fraction_before", 0.0)
        ),
        "p111634_high_near_spine_fraction_after": float(
            world.g("terrain.last_p111634_high_near_spine_fraction_after", 0.0)
        ),
        "p111634_guard_reverted": bool(
            world.g("terrain.last_p111634_guard_reverted", 0.0) >= 1.0
        ),
        "p111637c_anti_raster_orogenic_expression_accepted": bool(
            world.g(
                "terrain.last_p111637c_anti_raster_orogenic_expression_accepted",
                0.0,
            ) >= 1.0
        ),
        "p111637c_guard_reverted": bool(
            world.g("terrain.last_p111637c_guard_reverted", 0.0) >= 1.0
        ),
        "p111637c_land_mask_preserved": bool(
            world.g("terrain.last_p111637c_land_mask_preserved", 1.0) >= 1.0
        ),
        "p111637c_candidate_cell_count": float(
            world.g("terrain.last_p111637c_candidate_cell_count", 0.0)
        ),
        "p111637c_candidate_area_fraction": float(
            world.g("terrain.last_p111637c_candidate_area_fraction", 0.0)
        ),
        "p111637c_adjusted_cell_count": float(
            world.g("terrain.last_p111637c_adjusted_cell_count", 0.0)
        ),
        "p111637c_adjusted_area_fraction": float(
            world.g("terrain.last_p111637c_adjusted_area_fraction", 0.0)
        ),
        "p111637c_roughness_before_m": float(
            world.g("terrain.last_p111637c_roughness_before_m", 0.0)
        ),
        "p111637c_roughness_after_m": float(
            world.g("terrain.last_p111637c_roughness_after_m", 0.0)
        ),
        "p111637c_mean_abs_delta_m": float(
            world.g("terrain.last_p111637c_mean_abs_delta_m", 0.0)
        ),
        "p111637c_max_abs_delta_m": float(
            world.g("terrain.last_p111637c_max_abs_delta_m", 0.0)
        ),
        "p111637c_high_component_count_before": float(
            world.g("terrain.last_p111637c_high_component_count_before", 0.0)
        ),
        "p111637c_high_component_count_after": float(
            world.g("terrain.last_p111637c_high_component_count_after", 0.0)
        ),
        "p123_spine_high_relief_continuity_expression_accepted": bool(
            world.g(
                "terrain.last_p123_spine_high_relief_continuity_expression_accepted",
                0.0,
            ) >= 1.0
        ),
        "p123_guard_reverted": bool(
            world.g("terrain.last_p123_guard_reverted", 0.0) >= 1.0
        ),
        "p123_land_mask_preserved": bool(
            world.g("terrain.last_p123_land_mask_preserved", 1.0) >= 1.0
        ),
        "p123_candidate_cell_count": float(
            world.g("terrain.last_p123_candidate_cell_count", 0.0)
        ),
        "p123_candidate_area_fraction": float(
            world.g("terrain.last_p123_candidate_area_fraction", 0.0)
        ),
        "p123_adjusted_cell_count": float(
            world.g("terrain.last_p123_adjusted_cell_count", 0.0)
        ),
        "p123_adjusted_area_fraction": float(
            world.g("terrain.last_p123_adjusted_area_fraction", 0.0)
        ),
        "p123_spine_2400_component_count_before": float(
            world.g("terrain.last_p123_spine_2400_component_count_before", 0.0)
        ),
        "p123_spine_2400_component_count_after": float(
            world.g("terrain.last_p123_spine_2400_component_count_after", 0.0)
        ),
        "p123_spine_3000_component_count_before": float(
            world.g("terrain.last_p123_spine_3000_component_count_before", 0.0)
        ),
        "p123_spine_3000_component_count_after": float(
            world.g("terrain.last_p123_spine_3000_component_count_after", 0.0)
        ),
        "p123_spine_2400_coverage_before": float(
            world.g("terrain.last_p123_spine_2400_coverage_before", 0.0)
        ),
        "p123_spine_2400_coverage_after": float(
            world.g("terrain.last_p123_spine_2400_coverage_after", 0.0)
        ),
        "p123_high_component_count_before": float(
            world.g("terrain.last_p123_high_component_count_before", 0.0)
        ),
        "p123_high_component_count_after": float(
            world.g("terrain.last_p123_high_component_count_after", 0.0)
        ),
        "p123_mean_lift_m": float(
            world.g("terrain.last_p123_mean_lift_m", 0.0)
        ),
        "p123_max_lift_m": float(
            world.g("terrain.last_p123_max_lift_m", 0.0)
        ),
        "p123_p90_land_relief_before_m": float(
            world.g("terrain.last_p123_p90_land_relief_before_m", 0.0)
        ),
        "p123_p90_land_relief_after_m": float(
            world.g("terrain.last_p123_p90_land_relief_after_m", 0.0)
        ),
        "p123_p98_land_relief_before_m": float(
            world.g("terrain.last_p123_p98_land_relief_before_m", 0.0)
        ),
        "p123_p98_land_relief_after_m": float(
            world.g("terrain.last_p123_p98_land_relief_after_m", 0.0)
        ),
        "p134_terminal_spine_relief_gap_bridging_accepted": bool(
            world.g(
                "terrain.last_p134_terminal_spine_relief_gap_bridging_accepted",
                0.0,
            ) >= 1.0
        ),
        "p134_guard_reverted": bool(
            world.g("terrain.last_p134_guard_reverted", 0.0) >= 1.0
        ),
        "p134_land_mask_preserved": bool(
            world.g("terrain.last_p134_land_mask_preserved", 1.0) >= 1.0
        ),
        "p134_candidate_cell_count": float(
            world.g("terrain.last_p134_candidate_cell_count", 0.0)
        ),
        "p134_candidate_area_fraction": float(
            world.g("terrain.last_p134_candidate_area_fraction", 0.0)
        ),
        "p134_adjusted_cell_count": float(
            world.g("terrain.last_p134_adjusted_cell_count", 0.0)
        ),
        "p134_adjusted_area_fraction": float(
            world.g("terrain.last_p134_adjusted_area_fraction", 0.0)
        ),
        "p134_bridge_gap_component_count_2400": float(
            world.g("terrain.last_p134_bridge_gap_component_count_2400", 0.0)
        ),
        "p134_bridge_gap_component_count_3000": float(
            world.g("terrain.last_p134_bridge_gap_component_count_3000", 0.0)
        ),
        "p134_bridge_gap_cell_count_2400": float(
            world.g("terrain.last_p134_bridge_gap_cell_count_2400", 0.0)
        ),
        "p134_bridge_gap_cell_count_3000": float(
            world.g("terrain.last_p134_bridge_gap_cell_count_3000", 0.0)
        ),
        "p134_spine_2400_component_count_before": float(
            world.g("terrain.last_p134_spine_2400_component_count_before", 0.0)
        ),
        "p134_spine_2400_component_count_after": float(
            world.g("terrain.last_p134_spine_2400_component_count_after", 0.0)
        ),
        "p134_spine_3000_component_count_before": float(
            world.g("terrain.last_p134_spine_3000_component_count_before", 0.0)
        ),
        "p134_spine_3000_component_count_after": float(
            world.g("terrain.last_p134_spine_3000_component_count_after", 0.0)
        ),
        "p134_spine_2400_coverage_before": float(
            world.g("terrain.last_p134_spine_2400_coverage_before", 0.0)
        ),
        "p134_spine_2400_coverage_after": float(
            world.g("terrain.last_p134_spine_2400_coverage_after", 0.0)
        ),
        "p134_spine_3000_coverage_before": float(
            world.g("terrain.last_p134_spine_3000_coverage_before", 0.0)
        ),
        "p134_spine_3000_coverage_after": float(
            world.g("terrain.last_p134_spine_3000_coverage_after", 0.0)
        ),
        "p134_high_component_count_before": float(
            world.g("terrain.last_p134_high_component_count_before", 0.0)
        ),
        "p134_high_component_count_after": float(
            world.g("terrain.last_p134_high_component_count_after", 0.0)
        ),
        "p134_mean_lift_m": float(
            world.g("terrain.last_p134_mean_lift_m", 0.0)
        ),
        "p134_max_lift_m": float(
            world.g("terrain.last_p134_max_lift_m", 0.0)
        ),
        "p134_p90_land_relief_before_m": float(
            world.g("terrain.last_p134_p90_land_relief_before_m", 0.0)
        ),
        "p134_p90_land_relief_after_m": float(
            world.g("terrain.last_p134_p90_land_relief_after_m", 0.0)
        ),
        "p134_p98_land_relief_before_m": float(
            world.g("terrain.last_p134_p98_land_relief_before_m", 0.0)
        ),
        "p134_p98_land_relief_after_m": float(
            world.g("terrain.last_p134_p98_land_relief_after_m", 0.0)
        ),
        "p135_terminal_crest_core_relief_gap_bridging_accepted": bool(
            world.g(
                "terrain.last_p135_terminal_crest_core_relief_gap_bridging_accepted",
                0.0,
            ) >= 1.0
        ),
        "p135_guard_reverted": bool(
            world.g("terrain.last_p135_guard_reverted", 0.0) >= 1.0
        ),
        "p135_land_mask_preserved": bool(
            world.g("terrain.last_p135_land_mask_preserved", 1.0) >= 1.0
        ),
        "p135_candidate_cell_count": float(
            world.g("terrain.last_p135_candidate_cell_count", 0.0)
        ),
        "p135_candidate_area_fraction": float(
            world.g("terrain.last_p135_candidate_area_fraction", 0.0)
        ),
        "p135_adjusted_cell_count": float(
            world.g("terrain.last_p135_adjusted_cell_count", 0.0)
        ),
        "p135_adjusted_area_fraction": float(
            world.g("terrain.last_p135_adjusted_area_fraction", 0.0)
        ),
        "p135_gap_component_count": float(
            world.g("terrain.last_p135_gap_component_count", 0.0)
        ),
        "p135_gap_cell_count": float(
            world.g("terrain.last_p135_gap_cell_count", 0.0)
        ),
        "p135_crest_3000_component_count_before": float(
            world.g("terrain.last_p135_crest_3000_component_count_before", 0.0)
        ),
        "p135_crest_3000_component_count_after": float(
            world.g("terrain.last_p135_crest_3000_component_count_after", 0.0)
        ),
        "p135_spine_3000_component_count_before": float(
            world.g("terrain.last_p135_spine_3000_component_count_before", 0.0)
        ),
        "p135_spine_3000_component_count_after": float(
            world.g("terrain.last_p135_spine_3000_component_count_after", 0.0)
        ),
        "p135_high_component_count_before": float(
            world.g("terrain.last_p135_high_component_count_before", 0.0)
        ),
        "p135_high_component_count_after": float(
            world.g("terrain.last_p135_high_component_count_after", 0.0)
        ),
        "p135_crest_3000_coverage_before": float(
            world.g("terrain.last_p135_crest_3000_coverage_before", 0.0)
        ),
        "p135_crest_3000_coverage_after": float(
            world.g("terrain.last_p135_crest_3000_coverage_after", 0.0)
        ),
        "p135_spine_3000_coverage_before": float(
            world.g("terrain.last_p135_spine_3000_coverage_before", 0.0)
        ),
        "p135_spine_3000_coverage_after": float(
            world.g("terrain.last_p135_spine_3000_coverage_after", 0.0)
        ),
        "p135_mean_lift_m": float(
            world.g("terrain.last_p135_mean_lift_m", 0.0)
        ),
        "p135_max_lift_m": float(
            world.g("terrain.last_p135_max_lift_m", 0.0)
        ),
        "p135_p90_land_relief_before_m": float(
            world.g("terrain.last_p135_p90_land_relief_before_m", 0.0)
        ),
        "p135_p90_land_relief_after_m": float(
            world.g("terrain.last_p135_p90_land_relief_after_m", 0.0)
        ),
        "p135_p98_land_relief_before_m": float(
            world.g("terrain.last_p135_p98_land_relief_before_m", 0.0)
        ),
        "p135_p98_land_relief_after_m": float(
            world.g("terrain.last_p135_p98_land_relief_after_m", 0.0)
        ),
        "p136_terminal_high_peak_speckle_rebalance_accepted": bool(
            world.g(
                "terrain.last_p136_terminal_high_peak_speckle_rebalance_accepted",
                0.0,
            ) >= 1.0
        ),
        "p136_guard_reverted": bool(
            world.g("terrain.last_p136_guard_reverted", 0.0) >= 1.0
        ),
        "p136_land_mask_preserved": bool(
            world.g("terrain.last_p136_land_mask_preserved", 1.0) >= 1.0
        ),
        "p136_candidate_cell_count": float(
            world.g("terrain.last_p136_candidate_cell_count", 0.0)
        ),
        "p136_candidate_area_fraction": float(
            world.g("terrain.last_p136_candidate_area_fraction", 0.0)
        ),
        "p136_adjusted_cell_count": float(
            world.g("terrain.last_p136_adjusted_cell_count", 0.0)
        ),
        "p136_adjusted_area_fraction": float(
            world.g("terrain.last_p136_adjusted_area_fraction", 0.0)
        ),
        "p136_candidate_component_count": float(
            world.g("terrain.last_p136_candidate_component_count", 0.0)
        ),
        "p136_lowered_component_count": float(
            world.g("terrain.last_p136_lowered_component_count", 0.0)
        ),
        "p136_high_component_count_before": float(
            world.g("terrain.last_p136_high_component_count_before", 0.0)
        ),
        "p136_high_component_count_after": float(
            world.g("terrain.last_p136_high_component_count_after", 0.0)
        ),
        "p136_spine_3000_component_count_before": float(
            world.g("terrain.last_p136_spine_3000_component_count_before", 0.0)
        ),
        "p136_spine_3000_component_count_after": float(
            world.g("terrain.last_p136_spine_3000_component_count_after", 0.0)
        ),
        "p136_crest_3000_component_count_before": float(
            world.g("terrain.last_p136_crest_3000_component_count_before", 0.0)
        ),
        "p136_crest_3000_component_count_after": float(
            world.g("terrain.last_p136_crest_3000_component_count_after", 0.0)
        ),
        "p136_spine_2400_component_count_before": float(
            world.g("terrain.last_p136_spine_2400_component_count_before", 0.0)
        ),
        "p136_spine_2400_component_count_after": float(
            world.g("terrain.last_p136_spine_2400_component_count_after", 0.0)
        ),
        "p136_high_cell_count_before": float(
            world.g("terrain.last_p136_high_cell_count_before", 0.0)
        ),
        "p136_high_cell_count_after": float(
            world.g("terrain.last_p136_high_cell_count_after", 0.0)
        ),
        "p136_mean_lower_m": float(
            world.g("terrain.last_p136_mean_lower_m", 0.0)
        ),
        "p136_max_lower_m": float(
            world.g("terrain.last_p136_max_lower_m", 0.0)
        ),
        "p136_p90_land_relief_before_m": float(
            world.g("terrain.last_p136_p90_land_relief_before_m", 0.0)
        ),
        "p136_p90_land_relief_after_m": float(
            world.g("terrain.last_p136_p90_land_relief_after_m", 0.0)
        ),
        "p136_p98_land_relief_before_m": float(
            world.g("terrain.last_p136_p98_land_relief_before_m", 0.0)
        ),
        "p136_p98_land_relief_after_m": float(
            world.g("terrain.last_p136_p98_land_relief_after_m", 0.0)
        ),
    }


def _p120_continental_semantic_geometry_metrics(
    world: Any,
    rel: np.ndarray,
    land: np.ndarray,
    crust_type: np.ndarray,
    hierarchy: np.ndarray,
    spine: np.ndarray,
    shoulder_halo: np.ndarray,
    highland_apron: np.ndarray,
    mountain_inventory: np.ndarray,
) -> dict[str, Any]:
    """Measure raw vs regional continental semantic geometry for P120."""
    grid = world.grid
    n = grid.n
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    land = np.asarray(land, dtype=bool)
    crust = np.asarray(crust_type, dtype=np.int64)
    if crust.shape != (n,):
        crust = np.zeros(n, dtype=np.int64)
    cont_land = land & (crust >= 1)
    hierarchy = np.asarray(hierarchy, dtype=np.int64)
    spine = np.asarray(spine, dtype=np.int64)
    shoulder_halo = np.asarray(shoulder_halo, dtype=np.int64)
    highland_apron = np.asarray(highland_apron, dtype=np.int64)
    mountain_inventory = np.asarray(mountain_inventory, dtype=np.int64)
    orogenic_support = (
        (hierarchy >= 1)
        | (spine >= 2)
        | (shoulder_halo > 0)
        | (highland_apron > 0)
        | (mountain_inventory > 0)
    )
    unsupported_interior = cont_land & ~orogenic_support
    rel = np.asarray(rel, dtype=np.float64)

    def code_field(name: str) -> np.ndarray:
        arr = np.asarray(world.get_field(name, np.zeros(n)), dtype=np.float64)
        if arr.shape != (n,):
            return np.zeros(n, dtype=np.int64)
        return arr.astype(np.int64)

    fields = {
        "terrain_province": code_field("terrain.province"),
        "continental_detail_raw": code_field("terrain.continental_detail"),
        "continental_detail_region": code_field(
            "terrain.continental_detail_region_code"),
        "internal_block_raw": code_field(
            "tectonics.internal_geographic_block_code"),
        "internal_block_region": code_field(
            "terrain.internal_geographic_block_region_code"),
        "inland_geomorphology_region": code_field(
            "terrain.inland_geomorphology_region_code"),
        "terrain_continental_province": code_field(
            "terrain.continental_province_code"),
    }

    def component_elongation(comp: np.ndarray) -> float:
        if comp.size < 3:
            return 1.0
        xyz = np.asarray(grid.xyz[comp], dtype=np.float64)
        weights = area[comp]
        center = np.average(xyz, axis=0, weights=weights)
        centered = xyz - center
        cov = (centered * weights[:, None]).T @ centered / max(
            float(weights.sum()), 1.0)
        values = np.sort(np.maximum(np.linalg.eigvalsh(cov), 1.0e-18))[::-1]
        return float(np.sqrt(values[0] / max(values[1], 1.0e-18)))

    def field_metrics(code: np.ndarray, mask: np.ndarray) -> dict[str, Any]:
        code = np.asarray(code, dtype=np.int64)
        mask = np.asarray(mask, dtype=bool) & (code > 0)
        mask_area = max(float(area[mask].sum()), 1.0e-12)
        if not mask.any():
            return {
                "cell_count": 0,
                "area_fraction_world": 0.0,
                "same_neighbor_fraction": 1.0,
                "tiny_component_area_fraction": 0.0,
                "elongated_component_area_fraction": 0.0,
                "component_count": 0,
                "largest_component_area_fraction": 0.0,
                "largest_class_area_fraction": 0.0,
                "class_count": 0,
            }
        edges = np.asarray(grid.edges, dtype=np.int64)
        i = edges[:, 0]
        j = edges[:, 1]
        valid_edges = mask[i] & mask[j]
        same_neighbor = (
            float(np.count_nonzero(code[i[valid_edges]] == code[j[valid_edges]])
                  / max(np.count_nonzero(valid_edges), 1))
            if valid_edges.any() else 1.0
        )
        tiny_area = 0.0
        elongated_area = 0.0
        component_count = 0
        largest_component = 0.0
        class_areas: list[float] = []
        for value in sorted(int(x) for x in np.unique(code[mask]) if int(x) > 0):
            class_mask = mask & (code == value)
            class_area = float(area[class_mask].sum())
            if class_area > 0.0:
                class_areas.append(class_area)
            for comp in _components(grid, class_mask):
                component_count += 1
                comp_area = float(area[comp].sum())
                largest_component = max(largest_component, comp_area / mask_area)
                comp_share = comp_area / mask_area
                if comp.size <= 2 or comp_share <= 0.0025:
                    tiny_area += comp_area
                elongation = component_elongation(comp)
                if comp.size >= 8 and comp_share >= 0.006 and elongation >= 3.4:
                    elongated_area += comp_area
        return {
            "cell_count": int(np.count_nonzero(mask)),
            "area_fraction_world": float(area[mask].sum() / total),
            "same_neighbor_fraction": same_neighbor,
            "tiny_component_area_fraction": float(tiny_area / mask_area),
            "elongated_component_area_fraction": float(elongated_area / mask_area),
            "component_count": int(component_count),
            "largest_component_area_fraction": float(largest_component),
            "largest_class_area_fraction": float(
                max(class_areas) / mask_area if class_areas else 0.0),
            "class_count": int(len(class_areas)),
        }

    metrics: dict[str, Any] = {
        "schema": "aevum.p120_continental_semantic_geometry.v1",
        "continental_land_cell_count": int(np.count_nonzero(cont_land)),
        "continental_land_area_fraction_world": float(area[cont_land].sum() / total),
        "unsupported_interior_cell_count": int(np.count_nonzero(unsupported_interior)),
        "unsupported_interior_area_fraction_world": float(
            area[unsupported_interior].sum() / total),
        "region_fields_available": bool(
            np.any(fields["continental_detail_region"] > 0)
            or np.any(fields["internal_block_region"] > 0)
            or np.any(fields["inland_geomorphology_region"] > 0)
        ),
        "field_metrics": {},
    }
    for name, code in fields.items():
        metrics["field_metrics"][name] = {
            "continental_land": field_metrics(code, cont_land),
            "unsupported_interior": field_metrics(code, unsupported_interior),
        }

    raw_detail = metrics["field_metrics"]["continental_detail_raw"][
        "unsupported_interior"]
    region_detail = metrics["field_metrics"]["continental_detail_region"][
        "unsupported_interior"]
    raw_internal = metrics["field_metrics"]["internal_block_raw"][
        "unsupported_interior"]
    region_internal = metrics["field_metrics"]["internal_block_region"][
        "unsupported_interior"]
    terrain_province = metrics["field_metrics"]["terrain_province"][
        "unsupported_interior"]
    metrics["unsupported_detail_region_elongated_delta"] = float(
        region_detail["elongated_component_area_fraction"]
        - raw_detail["elongated_component_area_fraction"])
    metrics["unsupported_internal_region_elongated_delta"] = float(
        region_internal["elongated_component_area_fraction"]
        - raw_internal["elongated_component_area_fraction"])
    metrics["unsupported_detail_region_tiny_delta"] = float(
        region_detail["tiny_component_area_fraction"]
        - raw_detail["tiny_component_area_fraction"])
    metrics["unsupported_internal_region_tiny_delta"] = float(
        region_internal["tiny_component_area_fraction"]
        - raw_internal["tiny_component_area_fraction"])
    metrics["unsupported_terrain_province_elongated_area_fraction"] = float(
        terrain_province["elongated_component_area_fraction"])
    metrics["unsupported_continental_detail_raw_elongated_area_fraction"] = float(
        raw_detail["elongated_component_area_fraction"])
    metrics["unsupported_continental_detail_region_elongated_area_fraction"] = float(
        region_detail["elongated_component_area_fraction"])
    metrics["unsupported_internal_block_raw_elongated_area_fraction"] = float(
        raw_internal["elongated_component_area_fraction"])
    metrics["unsupported_internal_block_region_elongated_area_fraction"] = float(
        region_internal["elongated_component_area_fraction"])
    return metrics


def _orogenic_spine_linework_metrics(
    world: Any,
    land: np.ndarray,
    hierarchy: np.ndarray,
    spine: np.ndarray,
) -> dict[str, Any]:
    """Measure whether orogenic spines behave like continuous linework."""
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    land = np.asarray(land, dtype=bool)
    hierarchy = np.asarray(hierarchy, dtype=np.int64)
    spine = np.asarray(spine, dtype=np.int64)
    crust_type = np.asarray(world.get_field("crust.type", 0.0), dtype=np.float64)
    if crust_type.shape != (grid.n,):
        crust_type = np.zeros(grid.n, dtype=np.float64)
    continental_land = land & (crust_type >= 0.5)
    peak = continental_land & (hierarchy >= 2)
    crest = continental_land & (hierarchy >= 3)
    branch = continental_land & (hierarchy == 2)
    crest_spine = continental_land & crest & (spine >= 3)
    branch_spine = continental_land & branch & (spine == 2)
    spine_peak = crest_spine | branch_spine
    rel = (
        np.asarray(world.get_field("terrain.elevation_m", 0.0), dtype=np.float64)
        - float(getattr(world, "sea_level", 0.0))
    )
    if rel.shape != (grid.n,):
        rel = np.zeros(grid.n, dtype=np.float64)
    shoulder = np.asarray(
        world.get_field("terrain.orogenic_shoulder_halo", np.zeros(grid.n)),
        dtype=np.float64,
    )
    if shoulder.shape != (grid.n,):
        shoulder = np.zeros(grid.n, dtype=np.float64)
    apron = np.asarray(
        world.get_field("terrain.orogenic_highland_apron", np.zeros(grid.n)),
        dtype=np.float64,
    )
    if apron.shape != (grid.n,):
        apron = np.zeros(grid.n, dtype=np.float64)
    mountain_inventory = np.asarray(
        world.get_field("terrain.mountain_inventory", np.zeros(grid.n)),
        dtype=np.float64,
    )
    if mountain_inventory.shape != (grid.n,):
        mountain_inventory = np.zeros(grid.n, dtype=np.float64)
    spine_area = float(area[spine_peak].sum())
    peak_components = _components(grid, peak)
    spine_components = _components(grid, spine_peak)
    spine_component_count = len(spine_components)
    peak_component_count = len(peak_components)
    short_components = [comp for comp in spine_components if comp.size <= 2]
    spine_neighbor_count = _same_kind_neighbor_count(grid, spine_peak)
    endpoint_mask = spine_peak & (spine_neighbor_count <= 1)
    endpoint_count = int(np.count_nonzero(endpoint_mask))
    junction_count = int(np.count_nonzero(spine_peak & (spine_neighbor_count >= 3)))
    high_degree_count = int(np.count_nonzero(
        spine_peak & (spine_neighbor_count >= 4)))
    junction_pressure = float(
        high_degree_count / max(float(np.count_nonzero(spine_peak)), 1.0))
    component_areas = sorted(
        (float(area[comp].sum()) for comp in spine_components),
        reverse=True,
    )
    top3_share = (
        float(np.clip(sum(component_areas[:3]) / spine_area, 0.0, 1.0))
        if spine_area > 0.0 else 0.0
    )

    branch_components = _components(grid, branch_spine)
    if branch_components:
        crest_near = _dilate(grid, crest_spine, passes=1)
        attached = sum(1 for comp in branch_components if np.any(crest_near[comp]))
        branch_attachment_fraction = float(attached / max(float(len(branch_components)), 1.0))
    else:
        branch_attachment_fraction = 1.0

    short_fraction = float(
        len(short_components) / max(float(spine_component_count), 1.0))
    endpoints_per_component = float(
        endpoint_count / max(float(spine_component_count), 1.0))
    components_per_peak_component = float(
        spine_component_count / max(float(peak_component_count), 1.0))
    polar = np.abs(np.asarray(grid.lat, dtype=np.float64)) >= 60.0
    polar_spine_fraction = (
        float(area[spine_peak & polar].sum() / spine_area)
        if spine_area > 0.0 else 0.0
    )
    continuity_score = float(np.clip(
        top3_share
        * (1.0 - 0.65 * short_fraction)
        * (0.65 + 0.35 * branch_attachment_fraction)
        / max(components_per_peak_component, 1.0),
        0.0,
        1.0,
    ))
    high_500 = continental_land & (rel >= 500.0)
    high_1000 = continental_land & (rel >= 1000.0)
    orogenic_support = (
        (hierarchy >= 1)
        | (spine >= 2)
        | (shoulder > 0.0)
        | (apron > 0.0)
    )
    peak_support = hierarchy >= 2
    mountain_support = mountain_inventory > 0.0

    def area_fraction(mask: np.ndarray, support: np.ndarray) -> float:
        mask_area = float(area[mask].sum())
        if mask_area <= 0.0:
            return 1.0
        return float(np.clip(float(area[mask & support].sum()) / mask_area, 0.0, 1.0))

    non_orogenic_500 = high_500 & ~orogenic_support
    non_orogenic_1000 = high_1000 & ~orogenic_support

    return {
        "schema": "aevum.p111635_orogenic_spine_linework.v1",
        "peak_hierarchy_cell_count": int(np.count_nonzero(peak)),
        "peak_hierarchy_component_count": int(peak_component_count),
        "spine_cell_count": int(np.count_nonzero(spine_peak)),
        "spine_component_count": int(spine_component_count),
        "spine_components_per_peak_hierarchy_component": components_per_peak_component,
        "short_spine_component_count": int(len(short_components)),
        "short_spine_component_fraction": short_fraction,
        "spine_endpoint_count": endpoint_count,
        "spine_endpoints_per_component": endpoints_per_component,
        "spine_junction_count": junction_count,
        "spine_high_degree_count": high_degree_count,
        "spine_high_degree_fraction": junction_pressure,
        "spine_top3_component_share": top3_share,
        "branch_spine_component_count": int(len(branch_components)),
        "branch_attachment_fraction": branch_attachment_fraction,
        "polar_spine_area_fraction": polar_spine_fraction,
        "linework_continuity_score": continuity_score,
        "p117_highland_500_cell_count": int(np.count_nonzero(high_500)),
        "p117_highland_1000_cell_count": int(np.count_nonzero(high_1000)),
        "p117_highland_500_orogenic_support_fraction": area_fraction(
            high_500, orogenic_support),
        "p117_highland_1000_orogenic_support_fraction": area_fraction(
            high_1000, orogenic_support),
        "p117_highland_500_peak_support_fraction": area_fraction(
            high_500, peak_support),
        "p117_highland_1000_peak_support_fraction": area_fraction(
            high_1000, peak_support),
        "p117_highland_500_mountain_inventory_fraction": area_fraction(
            high_500, mountain_support),
        "p117_highland_1000_mountain_inventory_fraction": area_fraction(
            high_1000, mountain_support),
        "p117_non_orogenic_highland_500_cell_count": int(
            np.count_nonzero(non_orogenic_500)),
        "p117_non_orogenic_highland_1000_cell_count": int(
            np.count_nonzero(non_orogenic_1000)),
        "p117_non_orogenic_highland_500_area_fraction": float(
            area[non_orogenic_500].sum() / total),
        "p117_non_orogenic_highland_1000_area_fraction": float(
            area[non_orogenic_1000].sum() / total),
        "p118_non_orogenic_midland_stripe_suppression_accepted": bool(
            world.g(
                "terrain.last_p118_non_orogenic_midland_stripe_suppression_accepted",
                0.0,
            ) >= 1.0
        ),
        "p118_guard_reverted": bool(
            world.g("terrain.last_p118_guard_reverted", 0.0) >= 1.0
        ),
        "p118_land_mask_preserved": bool(
            world.g("terrain.last_p118_land_mask_preserved", 0.0) >= 1.0
        ),
        "p118_candidate_cell_count": float(
            world.g("terrain.last_p118_candidate_cell_count", 0.0)
        ),
        "p118_candidate_area_fraction": float(
            world.g("terrain.last_p118_candidate_area_fraction", 0.0)
        ),
        "p118_adjusted_cell_count": float(
            world.g("terrain.last_p118_adjusted_cell_count", 0.0)
        ),
        "p118_adjusted_area_fraction": float(
            world.g("terrain.last_p118_adjusted_area_fraction", 0.0)
        ),
        "p118_unsupported_midland_area_fraction_before": float(
            world.g(
                "terrain.last_p118_unsupported_midland_area_fraction_before",
                0.0,
            )
        ),
        "p118_unsupported_midland_area_fraction_after": float(
            world.g(
                "terrain.last_p118_unsupported_midland_area_fraction_after",
                0.0,
            )
        ),
        "p118_non_orogenic_highland_500_area_fraction_before": float(
            world.g(
                "terrain.last_p118_non_orogenic_highland_500_area_fraction_before",
                0.0,
            )
        ),
        "p118_non_orogenic_highland_500_area_fraction_after": float(
            world.g(
                "terrain.last_p118_non_orogenic_highland_500_area_fraction_after",
                0.0,
            )
        ),
        "p118_mean_lowering_m": float(
            world.g("terrain.last_p118_mean_lowering_m", 0.0)
        ),
        "p118_max_lowering_m": float(
            world.g("terrain.last_p118_max_lowering_m", 0.0)
        ),
        "p118_mean_land_relief_before_m": float(
            world.g("terrain.last_p118_mean_land_relief_before_m", 0.0)
        ),
        "p118_mean_land_relief_after_m": float(
            world.g("terrain.last_p118_mean_land_relief_after_m", 0.0)
        ),
        "p118_p50_land_relief_before_m": float(
            world.g("terrain.last_p118_p50_land_relief_before_m", 0.0)
        ),
        "p118_p50_land_relief_after_m": float(
            world.g("terrain.last_p118_p50_land_relief_after_m", 0.0)
        ),
        "p118_p90_land_relief_before_m": float(
            world.g("terrain.last_p118_p90_land_relief_before_m", 0.0)
        ),
        "p118_p90_land_relief_after_m": float(
            world.g("terrain.last_p118_p90_land_relief_after_m", 0.0)
        ),
        "p119_non_orogenic_interior_anti_raster_smoothing_accepted": bool(
            world.g(
                "terrain.last_p119_non_orogenic_interior_anti_raster_smoothing_accepted",
                0.0,
            ) >= 1.0
        ),
        "p119_guard_reverted": bool(
            world.g("terrain.last_p119_guard_reverted", 0.0) >= 1.0
        ),
        "p119_land_mask_preserved": bool(
            world.g("terrain.last_p119_land_mask_preserved", 0.0) >= 1.0
        ),
        "p119_candidate_cell_count": float(
            world.g("terrain.last_p119_candidate_cell_count", 0.0)
        ),
        "p119_candidate_area_fraction": float(
            world.g("terrain.last_p119_candidate_area_fraction", 0.0)
        ),
        "p119_selected_cell_count": float(
            world.g("terrain.last_p119_selected_cell_count", 0.0)
        ),
        "p119_selected_area_fraction": float(
            world.g("terrain.last_p119_selected_area_fraction", 0.0)
        ),
        "p119_adjusted_cell_count": float(
            world.g("terrain.last_p119_adjusted_cell_count", 0.0)
        ),
        "p119_adjusted_area_fraction": float(
            world.g("terrain.last_p119_adjusted_area_fraction", 0.0)
        ),
        "p119_compatible_edge_count": float(
            world.g("terrain.last_p119_compatible_edge_count", 0.0)
        ),
        "p119_roughness_before_m": float(
            world.g("terrain.last_p119_roughness_before_m", 0.0)
        ),
        "p119_roughness_after_m": float(
            world.g("terrain.last_p119_roughness_after_m", 0.0)
        ),
        "p119_non_orogenic_highland_500_area_fraction_before": float(
            world.g(
                "terrain.last_p119_non_orogenic_highland_500_area_fraction_before",
                0.0,
            )
        ),
        "p119_non_orogenic_highland_500_area_fraction_after": float(
            world.g(
                "terrain.last_p119_non_orogenic_highland_500_area_fraction_after",
                0.0,
            )
        ),
        "p119_mean_abs_delta_m": float(
            world.g("terrain.last_p119_mean_abs_delta_m", 0.0)
        ),
        "p119_max_abs_delta_m": float(
            world.g("terrain.last_p119_max_abs_delta_m", 0.0)
        ),
        "p119_p50_land_relief_before_m": float(
            world.g("terrain.last_p119_p50_land_relief_before_m", 0.0)
        ),
        "p119_p50_land_relief_after_m": float(
            world.g("terrain.last_p119_p50_land_relief_after_m", 0.0)
        ),
        "p119_p90_land_relief_before_m": float(
            world.g("terrain.last_p119_p90_land_relief_before_m", 0.0)
        ),
        "p119_p90_land_relief_after_m": float(
            world.g("terrain.last_p119_p90_land_relief_after_m", 0.0)
        ),
        "p120_continental_semantic_geometry_repair_accepted": bool(
            world.g(
                "terrain.last_p120_continental_semantic_geometry_repair_accepted",
                0.0,
            ) >= 1.0
        ),
        "p120_guard_reverted": bool(
            world.g("terrain.last_p120_guard_reverted", 0.0) >= 1.0
        ),
        "p120_candidate_cell_count": float(
            world.g("terrain.last_p120_candidate_cell_count", 0.0)
        ),
        "p120_candidate_area_fraction": float(
            world.g("terrain.last_p120_candidate_area_fraction", 0.0)
        ),
        "p120_terrain_changed_area_fraction": float(
            world.g("terrain.last_p120_terrain_changed_area_fraction", 0.0)
        ),
        "p120_internal_changed_area_fraction": float(
            world.g("terrain.last_p120_internal_changed_area_fraction", 0.0)
        ),
        "p120_detail_changed_area_fraction": float(
            world.g("terrain.last_p120_detail_changed_area_fraction", 0.0)
        ),
        "p120_terrain_elongated_before": float(
            world.g("terrain.last_p120_terrain_elongated_before", 0.0)
        ),
        "p120_terrain_elongated_after": float(
            world.g("terrain.last_p120_terrain_elongated_after", 0.0)
        ),
        "p120_internal_elongated_before": float(
            world.g("terrain.last_p120_internal_elongated_before", 0.0)
        ),
        "p120_internal_elongated_after": float(
            world.g("terrain.last_p120_internal_elongated_after", 0.0)
        ),
        "p120_detail_elongated_before": float(
            world.g("terrain.last_p120_detail_elongated_before", 0.0)
        ),
        "p120_detail_elongated_after": float(
            world.g("terrain.last_p120_detail_elongated_after", 0.0)
        ),
        "p120_terrain_tiny_before": float(
            world.g("terrain.last_p120_terrain_tiny_before", 0.0)
        ),
        "p120_terrain_tiny_after": float(
            world.g("terrain.last_p120_terrain_tiny_after", 0.0)
        ),
        "p120_internal_tiny_before": float(
            world.g("terrain.last_p120_internal_tiny_before", 0.0)
        ),
        "p120_internal_tiny_after": float(
            world.g("terrain.last_p120_internal_tiny_after", 0.0)
        ),
        "p120_detail_tiny_before": float(
            world.g("terrain.last_p120_detail_tiny_before", 0.0)
        ),
        "p120_detail_tiny_after": float(
            world.g("terrain.last_p120_detail_tiny_after", 0.0)
        ),
        "p122_orogenic_detail_footprint_compaction_accepted": bool(
            world.g(
                "terrain.last_p122_orogenic_detail_footprint_compaction_accepted",
                0.0,
            ) >= 1.0
        ),
        "p122_guard_reverted": bool(
            world.g("terrain.last_p122_guard_reverted", 0.0) >= 1.0
        ),
        "p122_candidate_cell_count": float(
            world.g("terrain.last_p122_candidate_cell_count", 0.0)
        ),
        "p122_candidate_area_fraction": float(
            world.g("terrain.last_p122_candidate_area_fraction", 0.0)
        ),
        "p122_detail_changed_area_fraction": float(
            world.g("terrain.last_p122_detail_changed_area_fraction", 0.0)
        ),
        "p122_internal_changed_area_fraction": float(
            world.g("terrain.last_p122_internal_changed_area_fraction", 0.0)
        ),
        "p122_terrain_changed_area_fraction": float(
            world.g("terrain.last_p122_terrain_changed_area_fraction", 0.0)
        ),
        "p122_detail_orogen_area_fraction_before": float(
            world.g("terrain.last_p122_detail_orogen_area_fraction_before", 0.0)
        ),
        "p122_detail_orogen_area_fraction_after": float(
            world.g("terrain.last_p122_detail_orogen_area_fraction_after", 0.0)
        ),
        "p122_low_relief_orogen_area_fraction_before": float(
            world.g(
                "terrain.last_p122_low_relief_orogen_area_fraction_before",
                0.0,
            )
        ),
        "p122_low_relief_orogen_area_fraction_after": float(
            world.g(
                "terrain.last_p122_low_relief_orogen_area_fraction_after",
                0.0,
            )
        ),
        "p122_noncore_orogen_area_fraction_before": float(
            world.g("terrain.last_p122_noncore_orogen_area_fraction_before", 0.0)
        ),
        "p122_noncore_orogen_area_fraction_after": float(
            world.g("terrain.last_p122_noncore_orogen_area_fraction_after", 0.0)
        ),
        "p122_preserved_core_orogen_cell_count": float(
            world.g("terrain.last_p122_preserved_core_orogen_cell_count", 0.0)
        ),
        "p121_semantic_region_elevation_response_accepted": bool(
            world.g(
                "terrain.last_p121_semantic_region_elevation_response_accepted",
                0.0,
            ) >= 1.0
        ),
        "p121_guard_reverted": bool(
            world.g("terrain.last_p121_guard_reverted", 0.0) >= 1.0
        ),
        "p121_land_mask_preserved": bool(
            world.g("terrain.last_p121_land_mask_preserved", 0.0) >= 1.0
        ),
        "p121_candidate_cell_count": float(
            world.g("terrain.last_p121_candidate_cell_count", 0.0)
        ),
        "p121_candidate_area_fraction": float(
            world.g("terrain.last_p121_candidate_area_fraction", 0.0)
        ),
        "p121_selected_cell_count": float(
            world.g("terrain.last_p121_selected_cell_count", 0.0)
        ),
        "p121_selected_area_fraction": float(
            world.g("terrain.last_p121_selected_area_fraction", 0.0)
        ),
        "p121_adjusted_cell_count": float(
            world.g("terrain.last_p121_adjusted_cell_count", 0.0)
        ),
        "p121_adjusted_area_fraction": float(
            world.g("terrain.last_p121_adjusted_area_fraction", 0.0)
        ),
        "p121_compatible_edge_count": float(
            world.g("terrain.last_p121_compatible_edge_count", 0.0)
        ),
        "p121_roughness_before_m": float(
            world.g("terrain.last_p121_roughness_before_m", 0.0)
        ),
        "p121_roughness_after_m": float(
            world.g("terrain.last_p121_roughness_after_m", 0.0)
        ),
        "p121_non_orogenic_highland_500_area_fraction_before": float(
            world.g(
                "terrain.last_p121_non_orogenic_highland_500_area_fraction_before",
                0.0,
            )
        ),
        "p121_non_orogenic_highland_500_area_fraction_after": float(
            world.g(
                "terrain.last_p121_non_orogenic_highland_500_area_fraction_after",
                0.0,
            )
        ),
        "p121_non_orogenic_highland_1000_area_fraction_before": float(
            world.g(
                "terrain.last_p121_non_orogenic_highland_1000_area_fraction_before",
                0.0,
            )
        ),
        "p121_non_orogenic_highland_1000_area_fraction_after": float(
            world.g(
                "terrain.last_p121_non_orogenic_highland_1000_area_fraction_after",
                0.0,
            )
        ),
        "p121_mean_abs_delta_m": float(
            world.g("terrain.last_p121_mean_abs_delta_m", 0.0)
        ),
        "p121_max_abs_delta_m": float(
            world.g("terrain.last_p121_max_abs_delta_m", 0.0)
        ),
        "p121_p50_land_relief_before_m": float(
            world.g("terrain.last_p121_p50_land_relief_before_m", 0.0)
        ),
        "p121_p50_land_relief_after_m": float(
            world.g("terrain.last_p121_p50_land_relief_after_m", 0.0)
        ),
        "p121_p90_land_relief_before_m": float(
            world.g("terrain.last_p121_p90_land_relief_before_m", 0.0)
        ),
        "p121_p90_land_relief_after_m": float(
            world.g("terrain.last_p121_p90_land_relief_after_m", 0.0)
        ),
        "p111635_spine_linework_smoothing_accepted": bool(
            world.g("terrain.last_p111635_spine_linework_smoothing_accepted", 0.0)
            >= 1.0
        ),
        "p111635_spine_components_before_refine": float(
            world.g("terrain.last_p111635_spine_components_before_refine", 0.0)
        ),
        "p111635_spine_components_after_refine": float(
            world.g("terrain.last_p111635_spine_components_after_refine", 0.0)
        ),
        "p111635_short_spine_components_before_refine": float(
            world.g("terrain.last_p111635_short_spine_components_before_refine", 0.0)
        ),
        "p111635_short_spine_components_after_refine": float(
            world.g("terrain.last_p111635_short_spine_components_after_refine", 0.0)
        ),
        "p111635_branch_attachment_fraction_before": float(
            world.g("terrain.last_p111635_branch_attachment_fraction_before", 0.0)
        ),
        "p111635_branch_attachment_fraction_after": float(
            world.g("terrain.last_p111635_branch_attachment_fraction_after", 0.0)
        ),
        "p111635_bridge_cell_count": float(
            world.g("terrain.last_p111635_bridge_cell_count", 0.0)
        ),
        "p111635_bridge_area_fraction": float(
            world.g("terrain.last_p111635_bridge_area_fraction", 0.0)
        ),
        "p111635_bridge_candidate_cell_count": float(
            world.g("terrain.last_p111635_bridge_candidate_cell_count", 0.0)
        ),
        "p111635_bridge_candidate_area_fraction": float(
            world.g("terrain.last_p111635_bridge_candidate_area_fraction", 0.0)
        ),
        "p111635_path_count": float(
            world.g("terrain.last_p111635_path_count", 0.0)
        ),
        "p111637a_spine_object_promotion_accepted": bool(
            world.g(
                "terrain.last_p111637a_spine_object_promotion_accepted",
                0.0,
            ) >= 1.0
        ),
        "p111637a_spine_components_before": float(
            world.g("terrain.last_p111637a_spine_components_before", 0.0)
        ),
        "p111637a_spine_components_after": float(
            world.g("terrain.last_p111637a_spine_components_after", 0.0)
        ),
        "p111637a_short_spine_components_before": float(
            world.g("terrain.last_p111637a_short_spine_components_before", 0.0)
        ),
        "p111637a_short_spine_components_after": float(
            world.g("terrain.last_p111637a_short_spine_components_after", 0.0)
        ),
        "p111637a_branch_attachment_fraction_before": float(
            world.g(
                "terrain.last_p111637a_branch_attachment_fraction_before",
                0.0,
            )
        ),
        "p111637a_branch_attachment_fraction_after": float(
            world.g(
                "terrain.last_p111637a_branch_attachment_fraction_after",
                0.0,
            )
        ),
        "p111637a_spine_top3_share_before": float(
            world.g("terrain.last_p111637a_spine_top3_share_before", 0.0)
        ),
        "p111637a_spine_top3_share_after": float(
            world.g("terrain.last_p111637a_spine_top3_share_after", 0.0)
        ),
        "p111637a_linework_score_before": float(
            world.g("terrain.last_p111637a_linework_score_before", 0.0)
        ),
        "p111637a_linework_score_after": float(
            world.g("terrain.last_p111637a_linework_score_after", 0.0)
        ),
        "p111637a_bridge_cell_count": float(
            world.g("terrain.last_p111637a_bridge_cell_count", 0.0)
        ),
        "p111637a_bridge_area_fraction": float(
            world.g("terrain.last_p111637a_bridge_area_fraction", 0.0)
        ),
        "p111637a_candidate_cell_count": float(
            world.g("terrain.last_p111637a_candidate_cell_count", 0.0)
        ),
        "p111637a_candidate_area_fraction": float(
            world.g("terrain.last_p111637a_candidate_area_fraction", 0.0)
        ),
        "p111637a_path_count": float(
            world.g("terrain.last_p111637a_path_count", 0.0)
        ),
        "p132_parent_anchor_spine_promotion_accepted": bool(
            world.g(
                "terrain.last_p132_parent_anchor_spine_promotion_accepted",
                0.0,
            ) >= 1.0
        ),
        "p132_candidate_cell_count": float(
            world.g("terrain.last_p132_candidate_cell_count", 0.0)
        ),
        "p132_candidate_area_fraction": float(
            world.g("terrain.last_p132_candidate_area_fraction", 0.0)
        ),
        "p132_promoted_spine_cell_count": float(
            world.g("terrain.last_p132_promoted_spine_cell_count", 0.0)
        ),
        "p132_promoted_spine_area_fraction": float(
            world.g("terrain.last_p132_promoted_spine_area_fraction", 0.0)
        ),
        "p132_parent_aligned_spine_fraction_before": float(
            world.g(
                "terrain.last_p132_parent_aligned_spine_fraction_before",
                0.0,
            )
        ),
        "p132_parent_aligned_spine_fraction_after": float(
            world.g(
                "terrain.last_p132_parent_aligned_spine_fraction_after",
                0.0,
            )
        ),
        "p132_linework_score_before": float(
            world.g("terrain.last_p132_linework_score_before", 0.0)
        ),
        "p132_linework_score_after": float(
            world.g("terrain.last_p132_linework_score_after", 0.0)
        ),
        "p132_spine_components_before": float(
            world.g("terrain.last_p132_spine_components_before", 0.0)
        ),
        "p132_spine_components_after": float(
            world.g("terrain.last_p132_spine_components_after", 0.0)
        ),
        "p132_short_spine_components_before": float(
            world.g("terrain.last_p132_short_spine_components_before", 0.0)
        ),
        "p132_short_spine_components_after": float(
            world.g("terrain.last_p132_short_spine_components_after", 0.0)
        ),
        "p1145_whole_mask_spine_planner_accepted": bool(
            world.g(
                "terrain.last_p1145_whole_mask_spine_planner_accepted",
                0.0,
            ) >= 1.0
        ),
        "p1145_bridge_cell_count": float(
            world.g("terrain.last_p1145_bridge_cell_count", 0.0)
        ),
        "p1145_bridge_area_fraction": float(
            world.g("terrain.last_p1145_bridge_area_fraction", 0.0)
        ),
        "p1145_candidate_cell_count": float(
            world.g("terrain.last_p1145_candidate_cell_count", 0.0)
        ),
        "p1145_candidate_area_fraction": float(
            world.g("terrain.last_p1145_candidate_area_fraction", 0.0)
        ),
        "p1145_path_count": float(
            world.g("terrain.last_p1145_path_count", 0.0)
        ),
        "p1145_linework_score_before": float(
            world.g("terrain.last_p1145_linework_score_before", 0.0)
        ),
        "p1145_linework_score_after": float(
            world.g("terrain.last_p1145_linework_score_after", 0.0)
        ),
        "p1145_spine_components_before": float(
            world.g("terrain.last_p1145_spine_components_before", 0.0)
        ),
        "p1145_spine_components_after": float(
            world.g("terrain.last_p1145_spine_components_after", 0.0)
        ),
        "p1145_short_spine_components_before": float(
            world.g("terrain.last_p1145_short_spine_components_before", 0.0)
        ),
        "p1145_short_spine_components_after": float(
            world.g("terrain.last_p1145_short_spine_components_after", 0.0)
        ),
        "p1145_branch_attachment_fraction_before": float(
            world.g("terrain.last_p1145_branch_attachment_fraction_before", 0.0)
        ),
        "p1145_branch_attachment_fraction_after": float(
            world.g("terrain.last_p1145_branch_attachment_fraction_after", 0.0)
        ),
        "p1145_spine_top3_share_before": float(
            world.g("terrain.last_p1145_spine_top3_share_before", 0.0)
        ),
        "p1145_spine_top3_share_after": float(
            world.g("terrain.last_p1145_spine_top3_share_after", 0.0)
        ),
        "p142_class_aware_combined_score": float(
            world.g("terrain.last_p142_class_aware_combined_score", 0.0)
        ),
        "p142_class_aware_class_score": float(
            world.g("terrain.last_p142_class_aware_class_score", 0.0)
        ),
        "p142_class_aware_crest_small_area_fraction": float(
            world.g(
                "terrain.last_p142_class_aware_crest_small_area_fraction",
                0.0,
            )
        ),
        "p142_class_aware_branch_small_area_fraction": float(
            world.g(
                "terrain.last_p142_class_aware_branch_small_area_fraction",
                0.0,
            )
        ),
        "p142_class_aware_class_small_area_fraction": float(
            world.g(
                "terrain.last_p142_class_aware_class_small_area_fraction",
                0.0,
            )
        ),
        "p142_class_aware_crest_component_count": float(
            world.g("terrain.last_p142_class_aware_crest_component_count", 0.0)
        ),
        "p142_class_aware_branch_component_count": float(
            world.g("terrain.last_p142_class_aware_branch_component_count", 0.0)
        ),
        "p142_class_aware_blind_spot": bool(
            world.g("terrain.last_p142_class_aware_blind_spot", 0.0) >= 1.0
        ),
        "p143_class_repair_needed": bool(
            world.g("terrain.last_p143_class_repair_needed", 0.0) >= 1.0
        ),
        "p143_class_profile_improved": bool(
            world.g("terrain.last_p143_class_profile_improved", 0.0) >= 1.0
        ),
        "p143_class_score_before": float(
            world.g("terrain.last_p143_class_score_before", 0.0)
        ),
        "p143_class_score_after": float(
            world.g("terrain.last_p143_class_score_after", 0.0)
        ),
        "p143_class_small_area_fraction_before": float(
            world.g("terrain.last_p143_class_small_area_fraction_before", 0.0)
        ),
        "p143_class_small_area_fraction_after": float(
            world.g("terrain.last_p143_class_small_area_fraction_after", 0.0)
        ),
        "p143_crest_small_area_fraction_before": float(
            world.g("terrain.last_p143_crest_small_area_fraction_before", 0.0)
        ),
        "p143_crest_small_area_fraction_after": float(
            world.g("terrain.last_p143_crest_small_area_fraction_after", 0.0)
        ),
        "p143_branch_small_area_fraction_before": float(
            world.g("terrain.last_p143_branch_small_area_fraction_before", 0.0)
        ),
        "p143_branch_small_area_fraction_after": float(
            world.g("terrain.last_p143_branch_small_area_fraction_after", 0.0)
        ),
        "p144_class_path_option_count": float(
            world.g("terrain.last_p144_class_path_option_count", 0.0)
        ),
        "p144_class_path_selected_count": float(
            world.g("terrain.last_p144_class_path_selected_count", 0.0)
        ),
        "p144_crest_path_option_count": float(
            world.g("terrain.last_p144_crest_path_option_count", 0.0)
        ),
        "p144_branch_path_option_count": float(
            world.g("terrain.last_p144_branch_path_option_count", 0.0)
        ),
        "p144_class_attempted_path_count": float(
            world.g("terrain.last_p144_class_attempted_path_count", 0.0)
        ),
        "p144_class_found_path_count": float(
            world.g("terrain.last_p144_class_found_path_count", 0.0)
        ),
        "p144_crest_promoted_spine_cell_count": float(
            world.g("terrain.last_p144_crest_promoted_spine_cell_count", 0.0)
        ),
        "p144_branch_promoted_spine_cell_count": float(
            world.g("terrain.last_p144_branch_promoted_spine_cell_count", 0.0)
        ),
        "p144_class_promoted_spine_cell_count": float(
            world.g("terrain.last_p144_class_promoted_spine_cell_count", 0.0)
        ),
        "p144_class_path_profile_rejected_count": float(
            world.g("terrain.last_p144_class_path_profile_rejected_count", 0.0)
        ),
        "p145_class_hierarchy_component_consolidation_accepted": bool(
            world.g(
                "terrain.last_p145_class_hierarchy_component_consolidation_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p145_removed_crest_component_count": float(
            world.g("terrain.last_p145_removed_crest_component_count", 0.0)
        ),
        "p145_removed_branch_component_count": float(
            world.g("terrain.last_p145_removed_branch_component_count", 0.0)
        ),
        "p145_removed_crest_cell_count": float(
            world.g("terrain.last_p145_removed_crest_cell_count", 0.0)
        ),
        "p145_removed_branch_cell_count": float(
            world.g("terrain.last_p145_removed_branch_cell_count", 0.0)
        ),
        "p145_removed_crest_area_fraction": float(
            world.g("terrain.last_p145_removed_crest_area_fraction", 0.0)
        ),
        "p145_removed_branch_area_fraction": float(
            world.g("terrain.last_p145_removed_branch_area_fraction", 0.0)
        ),
        "p145_added_crest_cell_count": float(
            world.g("terrain.last_p145_added_crest_cell_count", 0.0)
        ),
        "p145_added_branch_cell_count": float(
            world.g("terrain.last_p145_added_branch_cell_count", 0.0)
        ),
        "p145_added_crest_area_fraction": float(
            world.g("terrain.last_p145_added_crest_area_fraction", 0.0)
        ),
        "p145_added_branch_area_fraction": float(
            world.g("terrain.last_p145_added_branch_area_fraction", 0.0)
        ),
        "p145_bridge_path_count": float(
            world.g("terrain.last_p145_bridge_path_count", 0.0)
        ),
        "p145_class_score_before": float(
            world.g("terrain.last_p145_class_score_before", 0.0)
        ),
        "p145_class_score_after": float(
            world.g("terrain.last_p145_class_score_after", 0.0)
        ),
        "p145_class_small_area_fraction_before": float(
            world.g("terrain.last_p145_class_small_area_fraction_before", 0.0)
        ),
        "p145_class_small_area_fraction_after": float(
            world.g("terrain.last_p145_class_small_area_fraction_after", 0.0)
        ),
        "p145_crest_small_area_fraction_before": float(
            world.g("terrain.last_p145_crest_small_area_fraction_before", 0.0)
        ),
        "p145_crest_small_area_fraction_after": float(
            world.g("terrain.last_p145_crest_small_area_fraction_after", 0.0)
        ),
        "p145_branch_small_area_fraction_before": float(
            world.g("terrain.last_p145_branch_small_area_fraction_before", 0.0)
        ),
        "p145_branch_small_area_fraction_after": float(
            world.g("terrain.last_p145_branch_small_area_fraction_after", 0.0)
        ),
        "p1146_reject_code": int(
            round(float(world.g("terrain.last_p1146_reject_code", 0.0)))
        ),
        "p1146_reject_reason": _p1146_reject_reason_name(
            world.g("terrain.last_p1146_reject_code", 0.0)
        ),
        "p1146_support_component_count": float(
            world.g("terrain.last_p1146_support_component_count", 0.0)
        ),
        "p1146_multi_spine_support_component_count": float(
            world.g("terrain.last_p1146_multi_spine_support_component_count", 0.0)
        ),
        "p1146_attempted_path_count": float(
            world.g("terrain.last_p1146_attempted_path_count", 0.0)
        ),
        "p1146_found_path_count": float(
            world.g("terrain.last_p1146_found_path_count", 0.0)
        ),
        "p1146_bridge_path_count": float(
            world.g("terrain.last_p1146_bridge_path_count", 0.0)
        ),
        "p1146_attempted_bridge_cell_count": float(
            world.g("terrain.last_p1146_attempted_bridge_cell_count", 0.0)
        ),
        "p1146_trial_linework_score": float(
            world.g("terrain.last_p1146_trial_linework_score", 0.0)
        ),
        "p1146_trial_spine_components": float(
            world.g("terrain.last_p1146_trial_spine_components", 0.0)
        ),
        "p1146_trial_short_spine_components": float(
            world.g("terrain.last_p1146_trial_short_spine_components", 0.0)
        ),
        "p1146_trial_branch_attachment_fraction": float(
            world.g("terrain.last_p1146_trial_branch_attachment_fraction", 0.0)
        ),
        "p1146_trial_spine_top3_share": float(
            world.g("terrain.last_p1146_trial_spine_top3_share", 0.0)
        ),
        "p1147_pair_option_count": float(
            world.g("terrain.last_p1147_pair_option_count", 0.0)
        ),
        "p1147_pair_selected_count": float(
            world.g("terrain.last_p1147_pair_selected_count", 0.0)
        ),
        "p1147_pair_balance_rejected_count": float(
            world.g("terrain.last_p1147_pair_balance_rejected_count", 0.0)
        ),
        "p1147_pair_profile_rejected_count": float(
            world.g("terrain.last_p1147_pair_profile_rejected_count", 0.0)
        ),
        "p1148_terminal_proxy_enabled": bool(
            world.g("terrain.last_p1148_terminal_proxy_enabled", 0.0) >= 1.0
        ),
        "p1148_terminal_proxy_score_before": float(
            world.g("terrain.last_p1148_terminal_proxy_score_before", 0.0)
        ),
        "p1148_terminal_proxy_score_after": float(
            world.g("terrain.last_p1148_terminal_proxy_score_after", 0.0)
        ),
        "p1148_terminal_proxy_component_count_before": float(
            world.g(
                "terrain.last_p1148_terminal_proxy_component_count_before",
                0.0,
            )
        ),
        "p1148_terminal_proxy_component_count_after": float(
            world.g(
                "terrain.last_p1148_terminal_proxy_component_count_after",
                0.0,
            )
        ),
        "p1148_terminal_proxy_short_count_before": float(
            world.g("terrain.last_p1148_terminal_proxy_short_count_before", 0.0)
        ),
        "p1148_terminal_proxy_short_count_after": float(
            world.g("terrain.last_p1148_terminal_proxy_short_count_after", 0.0)
        ),
        "p1149_orogenic_spine_object_count": float(
            world.g("terrain.last_p1149_object_count", 0.0)
        ),
        "p1149_orogenic_spine_system_count": float(
            world.g("terrain.last_p1149_system_count", 0.0)
        ),
        "p1149_orogenic_spine_trunk_count": float(
            world.g("terrain.last_p1149_trunk_count", 0.0)
        ),
        "p1149_orogenic_spine_branch_count": float(
            world.g("terrain.last_p1149_branch_count", 0.0)
        ),
        "p1149_fallback_trunk_count": float(
            world.g("terrain.last_p1149_fallback_trunk_count", 0.0)
        ),
        "p1149_orphan_branch_count": float(
            world.g("terrain.last_p1149_orphan_branch_count", 0.0)
        ),
        "p1149_attached_branch_fraction": float(
            world.g("terrain.last_p1149_attached_branch_fraction", 0.0)
        ),
        "p1149_mean_branch_count_per_system": float(
            world.g("terrain.last_p1149_mean_branch_count_per_system", 0.0)
        ),
        "p1149_spine_cell_count": float(
            world.g("terrain.last_p1149_spine_cell_count", 0.0)
        ),
        "p1149_trunk_cell_count": float(
            world.g("terrain.last_p1149_trunk_cell_count", 0.0)
        ),
        "p1149_branch_cell_count": float(
            world.g("terrain.last_p1149_branch_cell_count", 0.0)
        ),
        "p1149_endpoint_count": float(
            world.g("terrain.last_p1149_endpoint_count", 0.0)
        ),
        "p128_orogenic_axis_polyline_object_count": float(
            world.g("terrain.last_p128_object_count", 0.0)
        ),
        "p128_main_axis_count": float(
            world.g("terrain.last_p128_main_axis_count", 0.0)
        ),
        "p128_branch_axis_count": float(
            world.g("terrain.last_p128_branch_axis_count", 0.0)
        ),
        "p128_fallback_main_axis_count": float(
            world.g("terrain.last_p128_fallback_main_axis_count", 0.0)
        ),
        "p128_attached_branch_axis_fraction": float(
            world.g("terrain.last_p128_attached_branch_axis_fraction", 0.0)
        ),
        "p128_source_spine_cell_count": float(
            world.g("terrain.last_p128_source_spine_cell_count", 0.0)
        ),
        "p128_axis_cell_count": float(
            world.g("terrain.last_p128_axis_cell_count", 0.0)
        ),
        "p128_axis_source_coverage_fraction": float(
            world.g("terrain.last_p128_axis_source_coverage_fraction", 0.0)
        ),
        "p128_mean_path_coverage_fraction": float(
            world.g("terrain.last_p128_mean_path_coverage_fraction", 0.0)
        ),
        "p128_mean_directness": float(
            world.g("terrain.last_p128_mean_directness", 0.0)
        ),
        "p128_mean_sinuosity": float(
            world.g("terrain.last_p128_mean_sinuosity", 1.0)
        ),
        "p128_max_sinuosity": float(
            world.g("terrain.last_p128_max_sinuosity", 1.0)
        ),
        "p128_source_junction_cell_count": float(
            world.g("terrain.last_p128_source_junction_cell_count", 0.0)
        ),
        "p128_source_high_degree_cell_count": float(
            world.g("terrain.last_p128_source_high_degree_cell_count", 0.0)
        ),
        "p128_polyline_ready": bool(
            world.g("terrain.last_p128_polyline_ready", 0.0) >= 1.0
        ),
        "p115_spine_repair_candidate_count": float(
            world.g("terrain.last_p115_candidate_count", 0.0)
        ),
        "p115_viable_spine_repair_candidate_count": float(
            world.g("terrain.last_p115_viable_candidate_count", 0.0)
        ),
        "p115_trunk_bridge_attempt_count": float(
            world.g("terrain.last_p115_trunk_bridge_attempt_count", 0.0)
        ),
        "p115_branch_attachment_attempt_count": float(
            world.g("terrain.last_p115_branch_attachment_attempt_count", 0.0)
        ),
        "p115_trunk_bridge_candidate_count": float(
            world.g("terrain.last_p115_trunk_bridge_candidate_count", 0.0)
        ),
        "p115_branch_attachment_candidate_count": float(
            world.g("terrain.last_p115_branch_attachment_candidate_count", 0.0)
        ),
        "p115_rejected_proxy_count": float(
            world.g("terrain.last_p115_rejected_proxy_count", 0.0)
        ),
        "p115_multi_trunk_system_count": float(
            world.g("terrain.last_p115_multi_trunk_system_count", 0.0)
        ),
        "p115_detached_branch_component_count": float(
            world.g("terrain.last_p115_detached_branch_component_count", 0.0)
        ),
        "p115_best_proxy_score_delta": float(
            world.g("terrain.last_p115_best_proxy_score_delta", 0.0)
        ),
        "p115_best_component_delta": float(
            world.g("terrain.last_p115_best_component_delta", 0.0)
        ),
        "p115_candidate_cell_count": float(
            world.g("terrain.last_p115_candidate_cell_count", 0.0)
        ),
        "p115_viable_candidate_cell_count": float(
            world.g("terrain.last_p115_viable_candidate_cell_count", 0.0)
        ),
        "p116_candidate_promotion_enabled": bool(
            world.g("terrain.last_p116_enabled", 0.0) >= 1.0
        ),
        "p116_candidate_promotion_accepted": bool(
            world.g("terrain.last_p116_accepted", 0.0) >= 1.0
        ),
        "p116_guard_reverted": bool(
            world.g("terrain.last_p116_guard_reverted", 0.0) >= 1.0
        ),
        "p116_input_candidate_count": float(
            world.g("terrain.last_p116_input_candidate_count", 0.0)
        ),
        "p116_input_viable_candidate_count": float(
            world.g("terrain.last_p116_input_viable_candidate_count", 0.0)
        ),
        "p116_considered_candidate_count": float(
            world.g("terrain.last_p116_considered_candidate_count", 0.0)
        ),
        "p116_selected_candidate_count": float(
            world.g("terrain.last_p116_selected_candidate_count", 0.0)
        ),
        "p116_applied_candidate_count": float(
            world.g("terrain.last_p116_applied_candidate_count", 0.0)
        ),
        "p116_rejected_candidate_count": float(
            world.g("terrain.last_p116_rejected_candidate_count", 0.0)
        ),
        "p116_rejected_overlap_count": float(
            world.g("terrain.last_p116_rejected_overlap_count", 0.0)
        ),
        "p116_rejected_support_count": float(
            world.g("terrain.last_p116_rejected_support_count", 0.0)
        ),
        "p116_rejected_profile_count": float(
            world.g("terrain.last_p116_rejected_profile_count", 0.0)
        ),
        "p116_applied_cell_count": float(
            world.g("terrain.last_p116_applied_cell_count", 0.0)
        ),
        "p116_applied_area_fraction": float(
            world.g("terrain.last_p116_applied_area_fraction", 0.0)
        ),
        "p116_area_budget_fraction": float(
            world.g("terrain.last_p116_area_budget_fraction", 0.0)
        ),
        "p116_linework_score_before": float(
            world.g("terrain.last_p116_linework_score_before", 0.0)
        ),
        "p116_linework_score_after": float(
            world.g("terrain.last_p116_linework_score_after", 0.0)
        ),
        "p116_component_count_before": float(
            world.g("terrain.last_p116_component_count_before", 0.0)
        ),
        "p116_component_count_after": float(
            world.g("terrain.last_p116_component_count_after", 0.0)
        ),
        "p116_short_count_before": float(
            world.g("terrain.last_p116_short_count_before", 0.0)
        ),
        "p116_short_count_after": float(
            world.g("terrain.last_p116_short_count_after", 0.0)
        ),
        "p116_branch_attachment_fraction_before": float(
            world.g("terrain.last_p116_branch_attachment_fraction_before", 0.0)
        ),
        "p116_branch_attachment_fraction_after": float(
            world.g("terrain.last_p116_branch_attachment_fraction_after", 0.0)
        ),
        "p116_spine_top3_share_before": float(
            world.g("terrain.last_p116_spine_top3_share_before", 0.0)
        ),
        "p116_spine_top3_share_after": float(
            world.g("terrain.last_p116_spine_top3_share_after", 0.0)
        ),
        "p124_orogenic_spine_geometry_regularization_accepted": bool(
            world.g(
                "terrain.last_p124_orogenic_spine_geometry_regularization_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p124_guard_reverted": bool(
            world.g("terrain.last_p124_guard_reverted", 0.0) >= 1.0
        ),
        "p124_candidate_cell_count": float(
            world.g("terrain.last_p124_candidate_cell_count", 0.0)
        ),
        "p124_candidate_area_fraction": float(
            world.g("terrain.last_p124_candidate_area_fraction", 0.0)
        ),
        "p124_added_cell_count": float(
            world.g("terrain.last_p124_added_cell_count", 0.0)
        ),
        "p124_added_area_fraction": float(
            world.g("terrain.last_p124_added_area_fraction", 0.0)
        ),
        "p124_component_count_before": float(
            world.g("terrain.last_p124_component_count_before", 0.0)
        ),
        "p124_component_count_after": float(
            world.g("terrain.last_p124_component_count_after", 0.0)
        ),
        "p124_endpoint_count_before": float(
            world.g("terrain.last_p124_endpoint_count_before", 0.0)
        ),
        "p124_endpoint_count_after": float(
            world.g("terrain.last_p124_endpoint_count_after", 0.0)
        ),
        "p124_short_component_count_before": float(
            world.g("terrain.last_p124_short_component_count_before", 0.0)
        ),
        "p124_short_component_count_after": float(
            world.g("terrain.last_p124_short_component_count_after", 0.0)
        ),
        "p124_branch_attachment_fraction_before": float(
            world.g("terrain.last_p124_branch_attachment_fraction_before", 0.0)
        ),
        "p124_branch_attachment_fraction_after": float(
            world.g("terrain.last_p124_branch_attachment_fraction_after", 0.0)
        ),
        "p124_linework_score_before": float(
            world.g("terrain.last_p124_linework_score_before", 0.0)
        ),
        "p124_linework_score_after": float(
            world.g("terrain.last_p124_linework_score_after", 0.0)
        ),
        "p127_terminal_orogenic_spine_node_thinning_accepted": bool(
            world.g(
                "terrain.last_p127_terminal_orogenic_spine_node_thinning_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p127_guard_reverted": bool(
            world.g("terrain.last_p127_guard_reverted", 0.0) >= 1.0
        ),
        "p127_candidate_cell_count": float(
            world.g("terrain.last_p127_candidate_cell_count", 0.0)
        ),
        "p127_candidate_area_fraction": float(
            world.g("terrain.last_p127_candidate_area_fraction", 0.0)
        ),
        "p127_demoted_cell_count": float(
            world.g("terrain.last_p127_demoted_cell_count", 0.0)
        ),
        "p127_demoted_area_fraction": float(
            world.g("terrain.last_p127_demoted_area_fraction", 0.0)
        ),
        "p127_spine_cell_count_before": float(
            world.g("terrain.last_p127_spine_cell_count_before", 0.0)
        ),
        "p127_spine_cell_count_after": float(
            world.g("terrain.last_p127_spine_cell_count_after", 0.0)
        ),
        "p127_component_count_before": float(
            world.g("terrain.last_p127_component_count_before", 0.0)
        ),
        "p127_component_count_after": float(
            world.g("terrain.last_p127_component_count_after", 0.0)
        ),
        "p127_endpoint_count_before": float(
            world.g("terrain.last_p127_endpoint_count_before", 0.0)
        ),
        "p127_endpoint_count_after": float(
            world.g("terrain.last_p127_endpoint_count_after", 0.0)
        ),
        "p127_junction_count_before": float(
            world.g("terrain.last_p127_junction_count_before", 0.0)
        ),
        "p127_junction_count_after": float(
            world.g("terrain.last_p127_junction_count_after", 0.0)
        ),
        "p127_high_degree_count_before": float(
            world.g("terrain.last_p127_high_degree_count_before", 0.0)
        ),
        "p127_high_degree_count_after": float(
            world.g("terrain.last_p127_high_degree_count_after", 0.0)
        ),
        "p127_short_component_count_before": float(
            world.g("terrain.last_p127_short_component_count_before", 0.0)
        ),
        "p127_short_component_count_after": float(
            world.g("terrain.last_p127_short_component_count_after", 0.0)
        ),
        "p127_branch_attachment_fraction_before": float(
            world.g("terrain.last_p127_branch_attachment_fraction_before", 0.0)
        ),
        "p127_branch_attachment_fraction_after": float(
            world.g("terrain.last_p127_branch_attachment_fraction_after", 0.0)
        ),
        "p127_linework_score_before": float(
            world.g("terrain.last_p127_linework_score_before", 0.0)
        ),
        "p127_linework_score_after": float(
            world.g("terrain.last_p127_linework_score_after", 0.0)
        ),
        "p127_junction_pressure_before": float(
            world.g("terrain.last_p127_junction_pressure_before", 0.0)
        ),
        "p127_junction_pressure_after": float(
            world.g("terrain.last_p127_junction_pressure_after", 0.0)
        ),
        "p140_terminal_branch_spine_component_thickening_accepted": bool(
            world.g(
                "terrain.last_p140_terminal_branch_spine_component_thickening_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p140_guard_reverted": bool(
            world.g("terrain.last_p140_guard_reverted", 0.0) >= 1.0
        ),
        "p140_candidate_cell_count": float(
            world.g("terrain.last_p140_candidate_cell_count", 0.0)
        ),
        "p140_candidate_area_fraction": float(
            world.g("terrain.last_p140_candidate_area_fraction", 0.0)
        ),
        "p140_promoted_branch_spine_cell_count": float(
            world.g("terrain.last_p140_promoted_branch_spine_cell_count", 0.0)
        ),
        "p140_promoted_branch_spine_area_fraction": float(
            world.g("terrain.last_p140_promoted_branch_spine_area_fraction", 0.0)
        ),
        "p140_grown_branch_spine_component_count": float(
            world.g("terrain.last_p140_grown_branch_spine_component_count", 0.0)
        ),
        "p140_branch_spine_cell_count_before": float(
            world.g("terrain.last_p140_branch_spine_cell_count_before", 0.0)
        ),
        "p140_branch_spine_cell_count_after": float(
            world.g("terrain.last_p140_branch_spine_cell_count_after", 0.0)
        ),
        "p140_branch_spine_component_count_before": float(
            world.g("terrain.last_p140_branch_spine_component_count_before", 0.0)
        ),
        "p140_branch_spine_component_count_after": float(
            world.g("terrain.last_p140_branch_spine_component_count_after", 0.0)
        ),
        "p140_branch_spine_small_area_fraction_before": float(
            world.g(
                "terrain.last_p140_branch_spine_small_area_fraction_before",
                0.0,
            )
        ),
        "p140_branch_spine_small_area_fraction_after": float(
            world.g(
                "terrain.last_p140_branch_spine_small_area_fraction_after",
                0.0,
            )
        ),
        "p140_all_spine_component_count_before": float(
            world.g("terrain.last_p140_all_spine_component_count_before", 0.0)
        ),
        "p140_all_spine_component_count_after": float(
            world.g("terrain.last_p140_all_spine_component_count_after", 0.0)
        ),
        "p140_junction_count_before": float(
            world.g("terrain.last_p140_junction_count_before", 0.0)
        ),
        "p140_junction_count_after": float(
            world.g("terrain.last_p140_junction_count_after", 0.0)
        ),
        "p140_high_degree_count_before": float(
            world.g("terrain.last_p140_high_degree_count_before", 0.0)
        ),
        "p140_high_degree_count_after": float(
            world.g("terrain.last_p140_high_degree_count_after", 0.0)
        ),
        "p140_junction_fraction_before": float(
            world.g("terrain.last_p140_junction_fraction_before", 0.0)
        ),
        "p140_junction_fraction_after": float(
            world.g("terrain.last_p140_junction_fraction_after", 0.0)
        ),
        "p140_high_degree_fraction_before": float(
            world.g("terrain.last_p140_high_degree_fraction_before", 0.0)
        ),
        "p140_high_degree_fraction_after": float(
            world.g("terrain.last_p140_high_degree_fraction_after", 0.0)
        ),
        "p140_branch_attachment_fraction_before": float(
            world.g("terrain.last_p140_branch_attachment_fraction_before", 0.0)
        ),
        "p140_branch_attachment_fraction_after": float(
            world.g("terrain.last_p140_branch_attachment_fraction_after", 0.0)
        ),
        "p140_hierarchy_changed_cell_count": float(
            world.g("terrain.last_p140_hierarchy_changed_cell_count", 0.0)
        ),
        "p140_crest_spine_changed_cell_count": float(
            world.g("terrain.last_p140_crest_spine_changed_cell_count", 0.0)
        ),
        "p141_terminal_crest_spine_component_thickening_accepted": bool(
            world.g(
                "terrain.last_p141_terminal_crest_spine_component_thickening_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p141_guard_reverted": bool(
            world.g("terrain.last_p141_guard_reverted", 0.0) >= 1.0
        ),
        "p141_candidate_cell_count": float(
            world.g("terrain.last_p141_candidate_cell_count", 0.0)
        ),
        "p141_candidate_area_fraction": float(
            world.g("terrain.last_p141_candidate_area_fraction", 0.0)
        ),
        "p141_promoted_crest_spine_cell_count": float(
            world.g("terrain.last_p141_promoted_crest_spine_cell_count", 0.0)
        ),
        "p141_promoted_crest_spine_area_fraction": float(
            world.g("terrain.last_p141_promoted_crest_spine_area_fraction", 0.0)
        ),
        "p141_connected_crest_spine_component_count": float(
            world.g("terrain.last_p141_connected_crest_spine_component_count", 0.0)
        ),
        "p141_grown_crest_spine_component_count": float(
            world.g("terrain.last_p141_grown_crest_spine_component_count", 0.0)
        ),
        "p141_crest_spine_cell_count_before": float(
            world.g("terrain.last_p141_crest_spine_cell_count_before", 0.0)
        ),
        "p141_crest_spine_cell_count_after": float(
            world.g("terrain.last_p141_crest_spine_cell_count_after", 0.0)
        ),
        "p141_crest_spine_component_count_before": float(
            world.g("terrain.last_p141_crest_spine_component_count_before", 0.0)
        ),
        "p141_crest_spine_component_count_after": float(
            world.g("terrain.last_p141_crest_spine_component_count_after", 0.0)
        ),
        "p141_crest_spine_small_area_fraction_before": float(
            world.g(
                "terrain.last_p141_crest_spine_small_area_fraction_before",
                0.0,
            )
        ),
        "p141_crest_spine_small_area_fraction_after": float(
            world.g(
                "terrain.last_p141_crest_spine_small_area_fraction_after",
                0.0,
            )
        ),
        "p141_all_spine_component_count_before": float(
            world.g("terrain.last_p141_all_spine_component_count_before", 0.0)
        ),
        "p141_all_spine_component_count_after": float(
            world.g("terrain.last_p141_all_spine_component_count_after", 0.0)
        ),
        "p141_junction_count_before": float(
            world.g("terrain.last_p141_junction_count_before", 0.0)
        ),
        "p141_junction_count_after": float(
            world.g("terrain.last_p141_junction_count_after", 0.0)
        ),
        "p141_high_degree_count_before": float(
            world.g("terrain.last_p141_high_degree_count_before", 0.0)
        ),
        "p141_high_degree_count_after": float(
            world.g("terrain.last_p141_high_degree_count_after", 0.0)
        ),
        "p141_junction_fraction_before": float(
            world.g("terrain.last_p141_junction_fraction_before", 0.0)
        ),
        "p141_junction_fraction_after": float(
            world.g("terrain.last_p141_junction_fraction_after", 0.0)
        ),
        "p141_high_degree_fraction_before": float(
            world.g("terrain.last_p141_high_degree_fraction_before", 0.0)
        ),
        "p141_high_degree_fraction_after": float(
            world.g("terrain.last_p141_high_degree_fraction_after", 0.0)
        ),
        "p141_branch_attachment_fraction_before": float(
            world.g("terrain.last_p141_branch_attachment_fraction_before", 0.0)
        ),
        "p141_branch_attachment_fraction_after": float(
            world.g("terrain.last_p141_branch_attachment_fraction_after", 0.0)
        ),
        "p141_hierarchy_changed_cell_count": float(
            world.g("terrain.last_p141_hierarchy_changed_cell_count", 0.0)
        ),
        "p141_branch_spine_changed_cell_count": float(
            world.g("terrain.last_p141_branch_spine_changed_cell_count", 0.0)
        ),
        "p164_terminal_orogenic_axis_skeleton_simplification_accepted": bool(
            world.g(
                "terrain.last_p164_terminal_orogenic_axis_skeleton_simplification_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p164_guard_reverted": bool(
            world.g("terrain.last_p164_guard_reverted", 0.0) >= 1.0
        ),
        "p164_candidate_cell_count": float(
            world.g("terrain.last_p164_candidate_cell_count", 0.0)
        ),
        "p164_candidate_area_fraction": float(
            world.g("terrain.last_p164_candidate_area_fraction", 0.0)
        ),
        "p164_demoted_cell_count": float(
            world.g("terrain.last_p164_demoted_cell_count", 0.0)
        ),
        "p164_demoted_area_fraction": float(
            world.g("terrain.last_p164_demoted_area_fraction", 0.0)
        ),
        "p164_off_axis_candidate_cell_count": float(
            world.g("terrain.last_p164_off_axis_candidate_cell_count", 0.0)
        ),
        "p164_spine_cell_count_before": float(
            world.g("terrain.last_p164_spine_cell_count_before", 0.0)
        ),
        "p164_spine_cell_count_after": float(
            world.g("terrain.last_p164_spine_cell_count_after", 0.0)
        ),
        "p164_component_count_before": float(
            world.g("terrain.last_p164_component_count_before", 0.0)
        ),
        "p164_component_count_after": float(
            world.g("terrain.last_p164_component_count_after", 0.0)
        ),
        "p164_short_component_count_before": float(
            world.g("terrain.last_p164_short_component_count_before", 0.0)
        ),
        "p164_short_component_count_after": float(
            world.g("terrain.last_p164_short_component_count_after", 0.0)
        ),
        "p164_endpoint_count_before": float(
            world.g("terrain.last_p164_endpoint_count_before", 0.0)
        ),
        "p164_endpoint_count_after": float(
            world.g("terrain.last_p164_endpoint_count_after", 0.0)
        ),
        "p164_junction_count_before": float(
            world.g("terrain.last_p164_junction_count_before", 0.0)
        ),
        "p164_junction_count_after": float(
            world.g("terrain.last_p164_junction_count_after", 0.0)
        ),
        "p164_high_degree_count_before": float(
            world.g("terrain.last_p164_high_degree_count_before", 0.0)
        ),
        "p164_high_degree_count_after": float(
            world.g("terrain.last_p164_high_degree_count_after", 0.0)
        ),
        "p164_junction_fraction_before": float(
            world.g("terrain.last_p164_junction_fraction_before", 0.0)
        ),
        "p164_junction_fraction_after": float(
            world.g("terrain.last_p164_junction_fraction_after", 0.0)
        ),
        "p164_high_degree_fraction_before": float(
            world.g("terrain.last_p164_high_degree_fraction_before", 0.0)
        ),
        "p164_high_degree_fraction_after": float(
            world.g("terrain.last_p164_high_degree_fraction_after", 0.0)
        ),
        "p164_branch_attachment_fraction_before": float(
            world.g("terrain.last_p164_branch_attachment_fraction_before", 0.0)
        ),
        "p164_branch_attachment_fraction_after": float(
            world.g("terrain.last_p164_branch_attachment_fraction_after", 0.0)
        ),
        "p164_linework_score_before": float(
            world.g("terrain.last_p164_linework_score_before", 0.0)
        ),
        "p164_linework_score_after": float(
            world.g("terrain.last_p164_linework_score_after", 0.0)
        ),
        "p164_axis_score_before": float(
            world.g("terrain.last_p164_axis_score_before", 0.0)
        ),
        "p164_axis_score_after": float(
            world.g("terrain.last_p164_axis_score_after", 0.0)
        ),
        "p164_protected_extreme_demoted_cell_count": float(
            world.g("terrain.last_p164_protected_extreme_demoted_cell_count", 0.0)
        ),
        "p164_reject_code": float(
            world.g("terrain.last_p164_reject_code", 0.0)
        ),
        "p142_class_aware_combined_score": float(
            world.g("terrain.last_p142_class_aware_combined_score", 0.0)
        ),
        "p142_class_aware_class_score": float(
            world.g("terrain.last_p142_class_aware_class_score", 0.0)
        ),
        "p142_class_aware_crest_small_area_fraction": float(
            world.g(
                "terrain.last_p142_class_aware_crest_small_area_fraction",
                0.0,
            )
        ),
        "p142_class_aware_branch_small_area_fraction": float(
            world.g(
                "terrain.last_p142_class_aware_branch_small_area_fraction",
                0.0,
            )
        ),
        "p142_class_aware_class_small_area_fraction": float(
            world.g(
                "terrain.last_p142_class_aware_class_small_area_fraction",
                0.0,
            )
        ),
        "p142_class_aware_crest_component_count": float(
            world.g("terrain.last_p142_class_aware_crest_component_count", 0.0)
        ),
        "p142_class_aware_branch_component_count": float(
            world.g("terrain.last_p142_class_aware_branch_component_count", 0.0)
        ),
        "p142_class_aware_blind_spot": bool(
            world.g("terrain.last_p142_class_aware_blind_spot", 0.0) >= 1.0
        ),
        "p143_class_repair_needed": bool(
            world.g("terrain.last_p143_class_repair_needed", 0.0) >= 1.0
        ),
        "p143_class_profile_improved": bool(
            world.g("terrain.last_p143_class_profile_improved", 0.0) >= 1.0
        ),
        "p143_class_score_before": float(
            world.g("terrain.last_p143_class_score_before", 0.0)
        ),
        "p143_class_score_after": float(
            world.g("terrain.last_p143_class_score_after", 0.0)
        ),
        "p143_class_small_area_fraction_before": float(
            world.g("terrain.last_p143_class_small_area_fraction_before", 0.0)
        ),
        "p143_class_small_area_fraction_after": float(
            world.g("terrain.last_p143_class_small_area_fraction_after", 0.0)
        ),
        "p143_crest_small_area_fraction_before": float(
            world.g("terrain.last_p143_crest_small_area_fraction_before", 0.0)
        ),
        "p143_crest_small_area_fraction_after": float(
            world.g("terrain.last_p143_crest_small_area_fraction_after", 0.0)
        ),
        "p143_branch_small_area_fraction_before": float(
            world.g("terrain.last_p143_branch_small_area_fraction_before", 0.0)
        ),
        "p143_branch_small_area_fraction_after": float(
            world.g("terrain.last_p143_branch_small_area_fraction_after", 0.0)
        ),
        "p144_class_path_option_count": float(
            world.g("terrain.last_p144_class_path_option_count", 0.0)
        ),
        "p144_class_path_selected_count": float(
            world.g("terrain.last_p144_class_path_selected_count", 0.0)
        ),
        "p144_crest_path_option_count": float(
            world.g("terrain.last_p144_crest_path_option_count", 0.0)
        ),
        "p144_branch_path_option_count": float(
            world.g("terrain.last_p144_branch_path_option_count", 0.0)
        ),
        "p144_class_attempted_path_count": float(
            world.g("terrain.last_p144_class_attempted_path_count", 0.0)
        ),
        "p144_class_found_path_count": float(
            world.g("terrain.last_p144_class_found_path_count", 0.0)
        ),
        "p144_crest_promoted_spine_cell_count": float(
            world.g("terrain.last_p144_crest_promoted_spine_cell_count", 0.0)
        ),
        "p144_branch_promoted_spine_cell_count": float(
            world.g("terrain.last_p144_branch_promoted_spine_cell_count", 0.0)
        ),
        "p144_class_promoted_spine_cell_count": float(
            world.g("terrain.last_p144_class_promoted_spine_cell_count", 0.0)
        ),
        "p144_class_path_profile_rejected_count": float(
            world.g("terrain.last_p144_class_path_profile_rejected_count", 0.0)
        ),
        "p145_class_hierarchy_component_consolidation_accepted": bool(
            world.g(
                "terrain.last_p145_class_hierarchy_component_consolidation_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p145_removed_crest_component_count": float(
            world.g("terrain.last_p145_removed_crest_component_count", 0.0)
        ),
        "p145_removed_branch_component_count": float(
            world.g("terrain.last_p145_removed_branch_component_count", 0.0)
        ),
        "p145_removed_crest_cell_count": float(
            world.g("terrain.last_p145_removed_crest_cell_count", 0.0)
        ),
        "p145_removed_branch_cell_count": float(
            world.g("terrain.last_p145_removed_branch_cell_count", 0.0)
        ),
        "p145_removed_crest_area_fraction": float(
            world.g("terrain.last_p145_removed_crest_area_fraction", 0.0)
        ),
        "p145_removed_branch_area_fraction": float(
            world.g("terrain.last_p145_removed_branch_area_fraction", 0.0)
        ),
        "p145_added_crest_cell_count": float(
            world.g("terrain.last_p145_added_crest_cell_count", 0.0)
        ),
        "p145_added_branch_cell_count": float(
            world.g("terrain.last_p145_added_branch_cell_count", 0.0)
        ),
        "p145_added_crest_area_fraction": float(
            world.g("terrain.last_p145_added_crest_area_fraction", 0.0)
        ),
        "p145_added_branch_area_fraction": float(
            world.g("terrain.last_p145_added_branch_area_fraction", 0.0)
        ),
        "p145_bridge_path_count": float(
            world.g("terrain.last_p145_bridge_path_count", 0.0)
        ),
        "p145_class_score_before": float(
            world.g("terrain.last_p145_class_score_before", 0.0)
        ),
        "p145_class_score_after": float(
            world.g("terrain.last_p145_class_score_after", 0.0)
        ),
        "p145_class_small_area_fraction_before": float(
            world.g("terrain.last_p145_class_small_area_fraction_before", 0.0)
        ),
        "p145_class_small_area_fraction_after": float(
            world.g("terrain.last_p145_class_small_area_fraction_after", 0.0)
        ),
        "p145_crest_small_area_fraction_before": float(
            world.g("terrain.last_p145_crest_small_area_fraction_before", 0.0)
        ),
        "p145_crest_small_area_fraction_after": float(
            world.g("terrain.last_p145_crest_small_area_fraction_after", 0.0)
        ),
        "p145_branch_small_area_fraction_before": float(
            world.g("terrain.last_p145_branch_small_area_fraction_before", 0.0)
        ),
        "p145_branch_small_area_fraction_after": float(
            world.g("terrain.last_p145_branch_small_area_fraction_after", 0.0)
        ),
        "p111637a_terminal_spine_land_consistency_accepted": bool(
            world.g(
                "terrain.last_p111637a_terminal_spine_land_consistency_accepted",
                0.0,
            ) >= 1.0
        ),
        "p111637a_terminal_submerged_spine_cell_count_before": float(
            world.g(
                "terrain.last_p111637a_terminal_submerged_spine_cell_count_before",
                0.0,
            )
        ),
        "p111637a_terminal_submerged_spine_cell_count_after": float(
            world.g(
                "terrain.last_p111637a_terminal_submerged_spine_cell_count_after",
                0.0,
            )
        ),
        "p111637a_terminal_spine_lift_cell_count": float(
            world.g("terrain.last_p111637a_terminal_spine_lift_cell_count", 0.0)
        ),
        "p111637a_terminal_spine_lift_area_fraction": float(
            world.g("terrain.last_p111637a_terminal_spine_lift_area_fraction", 0.0)
        ),
        "p111637b_final_spine_elevation_consistency_accepted": bool(
            world.g(
                "terrain.last_p111637b_final_spine_elevation_consistency_accepted",
                0.0,
            ) >= 1.0
        ),
        "p111637b_guard_reverted": bool(
            world.g("terrain.last_p111637b_guard_reverted", 0.0) >= 1.0
        ),
        "p111637b_candidate_cell_count": float(
            world.g("terrain.last_p111637b_candidate_cell_count", 0.0)
        ),
        "p111637b_candidate_area_fraction": float(
            world.g("terrain.last_p111637b_candidate_area_fraction", 0.0)
        ),
        "p111637b_adjusted_cell_count": float(
            world.g("terrain.last_p111637b_adjusted_cell_count", 0.0)
        ),
        "p111637b_adjusted_area_fraction": float(
            world.g("terrain.last_p111637b_adjusted_area_fraction", 0.0)
        ),
        "p111637b_submerged_spine_cell_count_before": float(
            world.g("terrain.last_p111637b_submerged_spine_cell_count_before", 0.0)
        ),
        "p111637b_submerged_spine_cell_count_after": float(
            world.g("terrain.last_p111637b_submerged_spine_cell_count_after", 0.0)
        ),
        "p111637b_bridge_blocked_submerged_spine_cell_count": float(
            world.g(
                "terrain.last_p111637b_bridge_blocked_submerged_spine_cell_count",
                0.0,
            )
        ),
        "p111637b_bridge_blocked_submerged_spine_area_fraction": float(
            world.g(
                "terrain.last_p111637b_bridge_blocked_submerged_spine_area_fraction",
                0.0,
            )
        ),
        "p111637b_underexpressed_spine_cell_count_before": float(
            world.g(
                "terrain.last_p111637b_underexpressed_spine_cell_count_before",
                0.0,
            )
        ),
        "p111637b_underexpressed_spine_cell_count_after": float(
            world.g(
                "terrain.last_p111637b_underexpressed_spine_cell_count_after",
                0.0,
            )
        ),
        "p111637b_land_delta_area_fraction": float(
            world.g("terrain.last_p111637b_land_delta_area_fraction", 0.0)
        ),
        "p111637b_mean_lift_m": float(
            world.g("terrain.last_p111637b_mean_lift_m", 0.0)
        ),
        "p111637b_max_lift_m": float(
            world.g("terrain.last_p111637b_max_lift_m", 0.0)
        ),
        "p111637b_linework_score_before": float(
            world.g("terrain.last_p111637b_linework_score_before", 0.0)
        ),
        "p111637b_linework_score_after": float(
            world.g("terrain.last_p111637b_linework_score_after", 0.0)
        ),
        "p1138_terminal_spine_stub_cleanup_accepted": bool(
            world.g(
                "terrain.last_p1138_terminal_spine_stub_cleanup_accepted",
                0.0,
            ) >= 1.0
        ),
        "p1138_guard_reverted": bool(
            world.g("terrain.last_p1138_guard_reverted", 0.0) >= 1.0
        ),
        "p1138_candidate_cell_count": float(
            world.g("terrain.last_p1138_candidate_cell_count", 0.0)
        ),
        "p1138_demoted_cell_count": float(
            world.g("terrain.last_p1138_demoted_cell_count", 0.0)
        ),
        "p1138_component_count_before": float(
            world.g("terrain.last_p1138_component_count_before", 0.0)
        ),
        "p1138_component_count_after": float(
            world.g("terrain.last_p1138_component_count_after", 0.0)
        ),
        "p1138_short_component_count_before": float(
            world.g("terrain.last_p1138_short_component_count_before", 0.0)
        ),
        "p1138_short_component_count_after": float(
            world.g("terrain.last_p1138_short_component_count_after", 0.0)
        ),
        "p1138_linework_score_before": float(
            world.g("terrain.last_p1138_linework_score_before", 0.0)
        ),
        "p1138_linework_score_after": float(
            world.g("terrain.last_p1138_linework_score_after", 0.0)
        ),
        "p1139_terminal_spine_leaf_tip_cleanup_accepted": bool(
            world.g(
                "terrain.last_p1139_terminal_spine_leaf_tip_cleanup_accepted",
                0.0,
            ) >= 1.0
        ),
        "p1139_guard_reverted": bool(
            world.g("terrain.last_p1139_guard_reverted", 0.0) >= 1.0
        ),
        "p1139_candidate_cell_count": float(
            world.g("terrain.last_p1139_candidate_cell_count", 0.0)
        ),
        "p1139_demoted_cell_count": float(
            world.g("terrain.last_p1139_demoted_cell_count", 0.0)
        ),
        "p1139_endpoint_count_before": float(
            world.g("terrain.last_p1139_endpoint_count_before", 0.0)
        ),
        "p1139_endpoint_count_after": float(
            world.g("terrain.last_p1139_endpoint_count_after", 0.0)
        ),
        "p1139_component_count_before": float(
            world.g("terrain.last_p1139_component_count_before", 0.0)
        ),
        "p1139_component_count_after": float(
            world.g("terrain.last_p1139_component_count_after", 0.0)
        ),
        "p1139_short_component_count_before": float(
            world.g("terrain.last_p1139_short_component_count_before", 0.0)
        ),
        "p1139_short_component_count_after": float(
            world.g("terrain.last_p1139_short_component_count_after", 0.0)
        ),
        "p1139_linework_score_before": float(
            world.g("terrain.last_p1139_linework_score_before", 0.0)
        ),
        "p1139_linework_score_after": float(
            world.g("terrain.last_p1139_linework_score_after", 0.0)
        ),
        "p1140_terminal_spine_triangular_kink_cleanup_accepted": bool(
            world.g(
                "terrain.last_p1140_terminal_spine_triangular_kink_cleanup_accepted",
                0.0,
            ) >= 1.0
        ),
        "p1140_guard_reverted": bool(
            world.g("terrain.last_p1140_guard_reverted", 0.0) >= 1.0
        ),
        "p1140_candidate_cell_count": float(
            world.g("terrain.last_p1140_candidate_cell_count", 0.0)
        ),
        "p1140_demoted_cell_count": float(
            world.g("terrain.last_p1140_demoted_cell_count", 0.0)
        ),
        "p1140_redundant_kink_count_before": float(
            world.g("terrain.last_p1140_redundant_kink_count_before", 0.0)
        ),
        "p1140_redundant_kink_count_after": float(
            world.g("terrain.last_p1140_redundant_kink_count_after", 0.0)
        ),
        "p1140_endpoint_count_before": float(
            world.g("terrain.last_p1140_endpoint_count_before", 0.0)
        ),
        "p1140_endpoint_count_after": float(
            world.g("terrain.last_p1140_endpoint_count_after", 0.0)
        ),
        "p1140_component_count_before": float(
            world.g("terrain.last_p1140_component_count_before", 0.0)
        ),
        "p1140_component_count_after": float(
            world.g("terrain.last_p1140_component_count_after", 0.0)
        ),
        "p1140_short_component_count_before": float(
            world.g("terrain.last_p1140_short_component_count_before", 0.0)
        ),
        "p1140_short_component_count_after": float(
            world.g("terrain.last_p1140_short_component_count_after", 0.0)
        ),
        "p1140_linework_score_before": float(
            world.g("terrain.last_p1140_linework_score_before", 0.0)
        ),
        "p1140_linework_score_after": float(
            world.g("terrain.last_p1140_linework_score_after", 0.0)
        ),
    }


def _terminal_relief_planform_metrics(
    world: Any,
    rel: np.ndarray,
    land: np.ndarray,
    hierarchy: np.ndarray,
    spine: np.ndarray,
    shoulder: np.ndarray,
    apron: np.ndarray,
) -> dict[str, Any]:
    """Measure final straight-highland and inland-plateau pressure."""
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    rel = np.asarray(rel, dtype=np.float64)
    land = np.asarray(land, dtype=bool)
    hierarchy = np.asarray(hierarchy, dtype=np.float64)
    spine = np.asarray(spine, dtype=np.float64)
    shoulder = np.asarray(shoulder, dtype=np.float64)
    apron = np.asarray(apron, dtype=np.float64)
    crust_type = np.asarray(world.get_field("crust.type", 0.0), dtype=np.float64)
    if crust_type.shape != (grid.n,):
        crust_type = np.zeros(grid.n, dtype=np.float64)
    if hierarchy.shape != (grid.n,):
        hierarchy = np.zeros(grid.n, dtype=np.float64)
    if spine.shape != (grid.n,):
        spine = np.zeros(grid.n, dtype=np.float64)
    if shoulder.shape != (grid.n,):
        shoulder = np.zeros(grid.n, dtype=np.float64)
    if apron.shape != (grid.n,):
        apron = np.zeros(grid.n, dtype=np.float64)

    continental = land & (crust_type >= 0.5)
    spine_core = continental & (spine >= 2.0)
    spine_near_1 = _dilate(grid, spine_core, passes=1) if spine_core.any() else spine_core
    spine_near_2 = _dilate(grid, spine_core, passes=2) if spine_core.any() else spine_core
    semantic_support = (
        (continental & (hierarchy >= 2.0))
        | spine_core
        | (continental & (shoulder > 0.0))
        | (continental & (apron > 0.0))
    )
    protected_core = spine_near_2 | (continental & (hierarchy >= 3.0))

    xyz = np.asarray(getattr(grid, "xyz", np.zeros((grid.n, 3))), dtype=np.float64)
    if xyz.shape != (grid.n, 3):
        xyz = np.zeros((grid.n, 3), dtype=np.float64)

    def component_shape(comp: np.ndarray) -> tuple[float, float]:
        comp = np.asarray(comp, dtype=np.int64)
        if comp.size < 3:
            return 1.0, 0.0
        pts = xyz[comp]
        pts = pts - np.mean(pts, axis=0)
        cov = (pts.T @ pts) / max(float(comp.size), 1.0)
        eig = np.sort(np.linalg.eigvalsh(cov))
        linearity = float(eig[-1] / max(float(eig[-2]), 1.0e-12))
        comp_mask = np.zeros(grid.n, dtype=bool)
        comp_mask[comp] = True

        def farthest(start: int) -> tuple[int, int]:
            seen = np.zeros(grid.n, dtype=bool)
            queue = [int(start)]
            seen[int(start)] = True
            dist = {int(start): 0}
            far = int(start)
            head = 0
            while head < len(queue):
                cell = queue[head]
                head += 1
                far = cell
                for nb in grid.neighbors[cell]:
                    nb = int(nb)
                    if not comp_mask[nb] or seen[nb]:
                        continue
                    seen[nb] = True
                    dist[nb] = dist[cell] + 1
                    queue.append(nb)
            return far, int(dist.get(far, 0))

        endpoint, _ = farthest(int(comp[0]))
        _, diameter = farthest(endpoint)
        elongation = float(diameter / max(np.sqrt(float(comp.size)), 1.0))
        return linearity, elongation

    highland = continental & (rel >= 900.0) & (rel < 4300.0)
    straight_count = 0
    straight_area = 0.0
    straight_score = 0.0
    for comp in _components(grid, highland):
        comp = np.asarray(comp, dtype=np.int64)
        if comp.size < 8:
            continue
        linearity, elongation = component_shape(comp)
        straight_like = (
            (linearity >= 4.2 and elongation >= 1.72)
            or elongation >= 3.80
        )
        if not straight_like:
            continue
        comp_area = float(area[comp].sum())
        support_fraction = float(
            area[comp[semantic_support[comp]]].sum() / max(comp_area, 1.0e-12)
        )
        straight_count += 1
        straight_area += comp_area
        straight_score += (
            comp_area
            / total
            * min(linearity / 8.0, 2.0)
            * (1.0 - 0.55 * min(support_fraction, 1.0))
        )

    def neighbor_count(mask: np.ndarray) -> np.ndarray:
        count = np.zeros(grid.n, dtype=np.int16)
        for cell in np.where(mask)[0]:
            count[np.asarray(grid.neighbors[int(cell)], dtype=np.int64)] += 1
        return count

    ocean = ~land
    coastal_land = continental & (neighbor_count(ocean) > 0)
    near_coast = _dilate(grid, coastal_land, passes=5) if coastal_land.any() else coastal_land
    inland = continental & ~near_coast
    plateau = (
        inland
        & (rel >= 1400.0)
        & (rel < 3800.0)
        & ~protected_core
    )
    broad_min_cells = max(12, int(np.ceil(0.0020 * float(grid.n))))
    plateau_components = _components(grid, plateau)
    broad_components = [
        comp for comp in plateau_components
        if np.asarray(comp, dtype=np.int64).size >= broad_min_cells
    ]
    plateau_area = float(area[plateau].sum())
    broad_areas = sorted(
        (float(area[np.asarray(comp, dtype=np.int64)].sum())
         for comp in broad_components),
        reverse=True,
    )
    largest_share = (
        float(broad_areas[0] / max(plateau_area, 1.0e-12))
        if broad_areas else 0.0
    )
    if plateau.any():
        q25 = _weighted_percentile(rel[plateau], area[plateau], 25.0)
        q75 = _weighted_percentile(rel[plateau], area[plateau], 75.0)
        plateau_iqr = float(q75 - q25)
    else:
        plateau_iqr = 0.0

    return {
        "schema": "aevum.p166167_terminal_relief_planform.v1",
        "straight_highland_component_count": int(straight_count),
        "straight_highland_area_fraction_world": float(straight_area / total),
        "straight_highland_score": float(straight_score),
        "inland_plateau_area_fraction_world": float(plateau_area / total),
        "broad_inland_plateau_component_count": int(len(broad_components)),
        "largest_inland_plateau_component_share": largest_share,
        "inland_plateau_relief_iqr_m": plateau_iqr,
        "p166_terminal_straight_highland_relief_softening_accepted": bool(
            world.g(
                "terrain.last_p166_terminal_straight_highland_relief_softening_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p166_guard_reverted": bool(
            world.g("terrain.last_p166_guard_reverted", 0.0) >= 1.0
        ),
        "p166_land_mask_preserved": bool(
            world.g("terrain.last_p166_land_mask_preserved", 0.0) >= 1.0
        ),
        "p166_candidate_cell_count": float(
            world.g("terrain.last_p166_candidate_cell_count", 0.0)
        ),
        "p166_adjusted_cell_count": float(
            world.g("terrain.last_p166_adjusted_cell_count", 0.0)
        ),
        "p166_straightness_score_before": float(
            world.g("terrain.last_p166_straightness_score_before", 0.0)
        ),
        "p166_straightness_score_after": float(
            world.g("terrain.last_p166_straightness_score_after", 0.0)
        ),
        "p166_straight_highland_area_fraction_before": float(
            world.g(
                "terrain.last_p166_straight_highland_area_fraction_before",
                0.0,
            )
        ),
        "p166_straight_highland_area_fraction_after": float(
            world.g(
                "terrain.last_p166_straight_highland_area_fraction_after",
                0.0,
            )
        ),
        "p166_protected_core_lowered_cell_count": float(
            world.g("terrain.last_p166_protected_core_lowered_cell_count", 0.0)
        ),
        "p166_extreme_softened_cell_count": float(
            world.g("terrain.last_p166_extreme_softened_cell_count", 0.0)
        ),
        "p166_reject_code": float(
            world.g("terrain.last_p166_reject_code", 0.0)
        ),
        "p167_terminal_inland_plateau_diversity_repair_accepted": bool(
            world.g(
                "terrain.last_p167_terminal_inland_plateau_diversity_repair_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p167_guard_reverted": bool(
            world.g("terrain.last_p167_guard_reverted", 0.0) >= 1.0
        ),
        "p167_land_mask_preserved": bool(
            world.g("terrain.last_p167_land_mask_preserved", 0.0) >= 1.0
        ),
        "p167_candidate_cell_count": float(
            world.g("terrain.last_p167_candidate_cell_count", 0.0)
        ),
        "p167_adjusted_cell_count": float(
            world.g("terrain.last_p167_adjusted_cell_count", 0.0)
        ),
        "p167_inland_plateau_area_fraction_before": float(
            world.g("terrain.last_p167_inland_plateau_area_fraction_before", 0.0)
        ),
        "p167_inland_plateau_area_fraction_after": float(
            world.g("terrain.last_p167_inland_plateau_area_fraction_after", 0.0)
        ),
        "p167_largest_plateau_component_share_before": float(
            world.g("terrain.last_p167_largest_plateau_component_share_before", 0.0)
        ),
        "p167_largest_plateau_component_share_after": float(
            world.g("terrain.last_p167_largest_plateau_component_share_after", 0.0)
        ),
        "p167_plateau_relief_iqr_before_m": float(
            world.g("terrain.last_p167_plateau_relief_iqr_before_m", 0.0)
        ),
        "p167_plateau_relief_iqr_after_m": float(
            world.g("terrain.last_p167_plateau_relief_iqr_after_m", 0.0)
        ),
        "p167_protected_core_lowered_cell_count": float(
            world.g("terrain.last_p167_protected_core_lowered_cell_count", 0.0)
        ),
        "p167_extreme_softened_cell_count": float(
            world.g("terrain.last_p167_extreme_softened_cell_count", 0.0)
        ),
        "p167_reject_code": float(
            world.g("terrain.last_p167_reject_code", 0.0)
        ),
        "p168_terminal_surface_derasterization_and_inland_diversity_accepted": bool(
            world.g(
                "terrain.last_p168_terminal_surface_derasterization_and_inland_diversity_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p168_guard_reverted": bool(
            world.g("terrain.last_p168_guard_reverted", 0.0) >= 1.0
        ),
        "p168_land_mask_preserved": bool(
            world.g("terrain.last_p168_land_mask_preserved", 0.0) >= 1.0
        ),
        "p168_candidate_cell_count": float(
            world.g("terrain.last_p168_candidate_cell_count", 0.0)
        ),
        "p168_adjusted_cell_count": float(
            world.g("terrain.last_p168_adjusted_cell_count", 0.0)
        ),
        "p168_land_adjusted_cell_count": float(
            world.g("terrain.last_p168_land_adjusted_cell_count", 0.0)
        ),
        "p168_ocean_adjusted_cell_count": float(
            world.g("terrain.last_p168_ocean_adjusted_cell_count", 0.0)
        ),
        "p168_weak_land_iqr_before_m": float(
            world.g("terrain.last_p168_weak_land_iqr_before_m", 0.0)
        ),
        "p168_weak_land_iqr_after_m": float(
            world.g("terrain.last_p168_weak_land_iqr_after_m", 0.0)
        ),
        "p168_weak_land_edge_delta_p95_before_m": float(
            world.g("terrain.last_p168_weak_land_edge_delta_p95_before_m", 0.0)
        ),
        "p168_weak_land_edge_delta_p95_after_m": float(
            world.g("terrain.last_p168_weak_land_edge_delta_p95_after_m", 0.0)
        ),
        "p168_weak_ocean_edge_delta_p95_before_m": float(
            world.g("terrain.last_p168_weak_ocean_edge_delta_p95_before_m", 0.0)
        ),
        "p168_weak_ocean_edge_delta_p95_after_m": float(
            world.g("terrain.last_p168_weak_ocean_edge_delta_p95_after_m", 0.0)
        ),
        "p168_broad_plateau_area_fraction_before": float(
            world.g("terrain.last_p168_broad_plateau_area_fraction_before", 0.0)
        ),
        "p168_broad_plateau_area_fraction_after": float(
            world.g("terrain.last_p168_broad_plateau_area_fraction_after", 0.0)
        ),
        "p168_protected_core_changed_cell_count": float(
            world.g("terrain.last_p168_protected_core_changed_cell_count", 0.0)
        ),
        "p168_ridge_trench_changed_cell_count": float(
            world.g("terrain.last_p168_ridge_trench_changed_cell_count", 0.0)
        ),
        "p168_extreme_changed_cell_count": float(
            world.g("terrain.last_p168_extreme_changed_cell_count", 0.0)
        ),
        "p168_reject_code": float(
            world.g("terrain.last_p168_reject_code", 0.0)
        ),
    }


def _polar_edge_orogen_overclassification_metrics(
    world: Any,
    land: np.ndarray,
    hierarchy: np.ndarray,
    spine: np.ndarray,
) -> dict[str, Any]:
    """Measure residual polar and map-edge peak-hierarchy overclassification."""
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    land = np.asarray(land, dtype=bool)
    hierarchy = np.asarray(hierarchy, dtype=np.int64)
    spine = np.asarray(spine, dtype=np.int64)
    halo = np.asarray(
        world.get_field("terrain.orogenic_shoulder_halo", np.zeros(grid.n)),
        dtype=np.int64,
    )
    apron = np.asarray(
        world.get_field("terrain.orogenic_highland_apron", np.zeros(grid.n)),
        dtype=np.int64,
    )
    if halo.shape != (grid.n,):
        halo = np.zeros(grid.n, dtype=np.int64)
    if apron.shape != (grid.n,):
        apron = np.zeros(grid.n, dtype=np.int64)
    crust_type = np.asarray(world.get_field("crust.type", 0.0), dtype=np.float64)
    if crust_type.shape != (grid.n,):
        crust_type = np.zeros(grid.n, dtype=np.float64)
    continental_land = land & (crust_type >= 0.5)
    submerged = ~land
    peak = continental_land & (hierarchy >= 2)
    spine_peak = continental_land & (spine >= 2)
    submerged_hierarchy = submerged & (hierarchy > 0)
    submerged_spine = submerged & (spine > 0)
    submerged_halo = submerged & (halo > 0)
    submerged_apron = submerged & (apron > 0)
    submerged_any_semantic = (
        submerged_hierarchy | submerged_spine | submerged_halo | submerged_apron
    )
    polar_band = np.abs(np.asarray(grid.lat, dtype=np.float64)) >= 60.0
    edge_band = np.abs(np.asarray(grid.lon, dtype=np.float64)) >= 172.0
    polar_peak = peak & polar_band
    polar_spine = spine_peak & polar_band
    edge_peak = peak & edge_band
    edge_components = _components(grid, edge_peak)
    polar_peak_area = float(area[polar_peak].sum())
    polar_spine_area = float(area[polar_spine].sum())
    polar_width_ratio = (
        float(polar_peak_area / polar_spine_area)
        if polar_spine_area > 0.0
        else (0.0 if polar_peak_area <= 0.0 else float("inf"))
    )
    if not np.isfinite(polar_width_ratio):
        polar_width_ratio = 1.0e9
    edge_short_components = sum(1 for comp in edge_components if comp.size <= 4)
    before_polar = float(
        world.g("terrain.last_p111636_polar_peak_area_fraction_before", 0.0)
    )
    after_polar = float(
        world.g("terrain.last_p111636_polar_peak_area_fraction_after", 0.0)
    )
    before_edge = float(
        world.g("terrain.last_p111636_edge_peak_component_count_before", 0.0)
    )
    after_edge = float(
        world.g("terrain.last_p111636_edge_peak_component_count_after", 0.0)
    )
    extreme_reclassified = float(
        world.g("terrain.last_p111636_extreme_reclassified_cell_count", 0.0)
    )
    p160_polar_before = float(
        world.g("terrain.last_p160_polar_spine_area_fraction_before", 0.0)
    )
    p160_polar_after = float(
        world.g("terrain.last_p160_polar_spine_area_fraction_after", 0.0)
    )
    p160_edge_before = float(
        world.g("terrain.last_p160_edge_spine_component_count_before", 0.0)
    )
    p160_edge_after = float(
        world.g("terrain.last_p160_edge_spine_component_count_after", 0.0)
    )
    p160_spine_components_before = float(
        world.g("terrain.last_p160_spine_component_count_before", 0.0)
    )
    p160_spine_components_after = float(
        world.g("terrain.last_p160_spine_component_count_after", 0.0)
    )
    p160_short_before = float(
        world.g("terrain.last_p160_short_spine_component_count_before", 0.0)
    )
    p160_short_after = float(
        world.g("terrain.last_p160_short_spine_component_count_after", 0.0)
    )
    p160_extreme_demoted = float(
        world.g("terrain.last_p160_extreme_demoted_cell_count", 0.0)
    )
    p165_polar_semantic_before = float(
        world.g("terrain.last_p165_polar_semantic_area_fraction_before", 0.0)
    )
    p165_polar_semantic_after = float(
        world.g("terrain.last_p165_polar_semantic_area_fraction_after", 0.0)
    )
    p165_edge_semantic_before = float(
        world.g("terrain.last_p165_edge_semantic_area_fraction_before", 0.0)
    )
    p165_edge_semantic_after = float(
        world.g("terrain.last_p165_edge_semantic_area_fraction_after", 0.0)
    )
    p165_polar_peak_before = float(
        world.g("terrain.last_p165_polar_peak_area_fraction_before", 0.0)
    )
    p165_polar_peak_after = float(
        world.g("terrain.last_p165_polar_peak_area_fraction_after", 0.0)
    )
    p165_edge_components_before = float(
        world.g("terrain.last_p165_edge_peak_component_count_before", 0.0)
    )
    p165_edge_components_after = float(
        world.g("terrain.last_p165_edge_peak_component_count_after", 0.0)
    )
    p165_extreme_reclassified = float(
        world.g("terrain.last_p165_protected_extreme_reclassified_cell_count", 0.0)
    )
    p165_spine_changed = float(
        world.g("terrain.last_p165_protected_spine_changed_cell_count", 0.0)
    )
    return {
        "schema": "aevum.p111636_polar_edge_orogen_overclassification.v3",
        "polar_peak_cell_count": int(np.count_nonzero(polar_peak)),
        "polar_peak_area_fraction_world": float(polar_peak_area / total),
        "polar_spine_cell_count": int(np.count_nonzero(polar_spine)),
        "polar_spine_area_fraction_world": float(polar_spine_area / total),
        "polar_peak_to_spine_width_ratio": float(polar_width_ratio),
        "edge_peak_cell_count": int(np.count_nonzero(edge_peak)),
        "edge_peak_component_count": int(len(edge_components)),
        "edge_short_peak_component_count": int(edge_short_components),
        "submerged_orogenic_semantic_cell_count": int(
            np.count_nonzero(submerged_any_semantic)
        ),
        "submerged_hierarchy_cell_count": int(np.count_nonzero(submerged_hierarchy)),
        "submerged_spine_cell_count": int(np.count_nonzero(submerged_spine)),
        "submerged_halo_cell_count": int(np.count_nonzero(submerged_halo)),
        "submerged_apron_cell_count": int(np.count_nonzero(submerged_apron)),
        "p111636_polar_edge_refinement_accepted": bool(
            world.g(
                "terrain.last_p111636_polar_edge_refinement_accepted",
                0.0,
            ) >= 1.0
        ),
        "p111636_polar_peak_area_fraction_before": before_polar,
        "p111636_polar_peak_area_fraction_after": after_polar,
        "p111636_polar_peak_area_fraction_delta": float(before_polar - after_polar),
        "p111636_edge_peak_component_count_before": before_edge,
        "p111636_edge_peak_component_count_after": after_edge,
        "p111636_edge_peak_component_count_delta": float(before_edge - after_edge),
        "p111636_candidate_cell_count": float(
            world.g("terrain.last_p111636_candidate_cell_count", 0.0)
        ),
        "p111636_candidate_area_fraction": float(
            world.g("terrain.last_p111636_candidate_area_fraction", 0.0)
        ),
        "p111636_reclassified_cell_count": float(
            world.g("terrain.last_p111636_reclassified_cell_count", 0.0)
        ),
        "p111636_reclassified_area_fraction": float(
            world.g("terrain.last_p111636_reclassified_area_fraction", 0.0)
        ),
        "p111636_polar_reclassified_cell_count": float(
            world.g("terrain.last_p111636_polar_reclassified_cell_count", 0.0)
        ),
        "p111636_edge_reclassified_cell_count": float(
            world.g("terrain.last_p111636_edge_reclassified_cell_count", 0.0)
        ),
        "p111636_extreme_reclassified_cell_count": extreme_reclassified,
        "p111636_no_extreme_reclassification": bool(extreme_reclassified == 0.0),
        "p111636_nonexpansive_polar_edge": bool(
            after_polar <= before_polar + 1.0e-12
            and after_edge <= before_edge + 1.0e-12
        ),
        "p160_polar_edge_spine_thinning_accepted": bool(
            world.g(
                "terrain.last_p160_polar_edge_spine_thinning_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p160_guard_reverted": bool(
            world.g("terrain.last_p160_guard_reverted", 0.0) >= 1.0
        ),
        "p160_candidate_cell_count": float(
            world.g("terrain.last_p160_candidate_cell_count", 0.0)
        ),
        "p160_candidate_area_fraction": float(
            world.g("terrain.last_p160_candidate_area_fraction", 0.0)
        ),
        "p160_demoted_cell_count": float(
            world.g("terrain.last_p160_demoted_cell_count", 0.0)
        ),
        "p160_demoted_area_fraction": float(
            world.g("terrain.last_p160_demoted_area_fraction", 0.0)
        ),
        "p160_polar_demoted_cell_count": float(
            world.g("terrain.last_p160_polar_demoted_cell_count", 0.0)
        ),
        "p160_edge_demoted_cell_count": float(
            world.g("terrain.last_p160_edge_demoted_cell_count", 0.0)
        ),
        "p160_extreme_demoted_cell_count": p160_extreme_demoted,
        "p160_polar_spine_area_fraction_before": p160_polar_before,
        "p160_polar_spine_area_fraction_after": p160_polar_after,
        "p160_polar_spine_area_fraction_delta": float(
            p160_polar_before - p160_polar_after
        ),
        "p160_edge_spine_component_count_before": p160_edge_before,
        "p160_edge_spine_component_count_after": p160_edge_after,
        "p160_edge_spine_component_count_delta": float(
            p160_edge_before - p160_edge_after
        ),
        "p160_spine_component_count_before": p160_spine_components_before,
        "p160_spine_component_count_after": p160_spine_components_after,
        "p160_short_spine_component_count_before": p160_short_before,
        "p160_short_spine_component_count_after": p160_short_after,
        "p160_branch_attachment_fraction_before": float(
            world.g("terrain.last_p160_branch_attachment_fraction_before", 0.0)
        ),
        "p160_branch_attachment_fraction_after": float(
            world.g("terrain.last_p160_branch_attachment_fraction_after", 0.0)
        ),
        "p160_linework_score_before": float(
            world.g("terrain.last_p160_linework_score_before", 0.0)
        ),
        "p160_linework_score_after": float(
            world.g("terrain.last_p160_linework_score_after", 0.0)
        ),
        "p160_no_extreme_demotion": bool(p160_extreme_demoted == 0.0),
        "p160_nonexpansive_polar_edge_spine": bool(
            p160_polar_after <= p160_polar_before + 1.0e-12
            and p160_edge_after <= p160_edge_before + 1.0e-12
            and p160_spine_components_after <= (
                p160_spine_components_before + 1.0e-12)
            and p160_short_after <= p160_short_before + 1.0e-12
        ),
        "p165_terminal_polar_edge_orogenic_semantic_compaction_accepted": bool(
            world.g(
                "terrain.last_p165_terminal_polar_edge_orogenic_semantic_compaction_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p165_guard_reverted": bool(
            world.g("terrain.last_p165_guard_reverted", 0.0) >= 1.0
        ),
        "p165_candidate_cell_count": float(
            world.g("terrain.last_p165_candidate_cell_count", 0.0)
        ),
        "p165_candidate_area_fraction": float(
            world.g("terrain.last_p165_candidate_area_fraction", 0.0)
        ),
        "p165_peak_reclassified_cell_count": float(
            world.g("terrain.last_p165_peak_reclassified_cell_count", 0.0)
        ),
        "p165_peak_reclassified_area_fraction": float(
            world.g("terrain.last_p165_peak_reclassified_area_fraction", 0.0)
        ),
        "p165_fringe_cleared_cell_count": float(
            world.g("terrain.last_p165_fringe_cleared_cell_count", 0.0)
        ),
        "p165_fringe_cleared_area_fraction": float(
            world.g("terrain.last_p165_fringe_cleared_area_fraction", 0.0)
        ),
        "p165_polar_semantic_area_fraction_before": p165_polar_semantic_before,
        "p165_polar_semantic_area_fraction_after": p165_polar_semantic_after,
        "p165_polar_semantic_area_fraction_delta": float(
            p165_polar_semantic_before - p165_polar_semantic_after
        ),
        "p165_edge_semantic_area_fraction_before": p165_edge_semantic_before,
        "p165_edge_semantic_area_fraction_after": p165_edge_semantic_after,
        "p165_edge_semantic_area_fraction_delta": float(
            p165_edge_semantic_before - p165_edge_semantic_after
        ),
        "p165_polar_peak_area_fraction_before": p165_polar_peak_before,
        "p165_polar_peak_area_fraction_after": p165_polar_peak_after,
        "p165_polar_peak_area_fraction_delta": float(
            p165_polar_peak_before - p165_polar_peak_after
        ),
        "p165_edge_peak_component_count_before": p165_edge_components_before,
        "p165_edge_peak_component_count_after": p165_edge_components_after,
        "p165_edge_peak_component_count_delta": float(
            p165_edge_components_before - p165_edge_components_after
        ),
        "p165_peak_component_count_before": float(
            world.g("terrain.last_p165_peak_component_count_before", 0.0)
        ),
        "p165_peak_component_count_after": float(
            world.g("terrain.last_p165_peak_component_count_after", 0.0)
        ),
        "p165_hierarchy_changed_cell_count": float(
            world.g("terrain.last_p165_hierarchy_changed_cell_count", 0.0)
        ),
        "p165_halo_changed_cell_count": float(
            world.g("terrain.last_p165_halo_changed_cell_count", 0.0)
        ),
        "p165_apron_changed_cell_count": float(
            world.g("terrain.last_p165_apron_changed_cell_count", 0.0)
        ),
        "p165_protected_spine_changed_cell_count": p165_spine_changed,
        "p165_protected_extreme_reclassified_cell_count": p165_extreme_reclassified,
        "p165_reject_code": float(
            world.g("terrain.last_p165_reject_code", 0.0)
        ),
        "p165_no_protected_spine_change": bool(p165_spine_changed == 0.0),
        "p165_no_extreme_reclassification": bool(p165_extreme_reclassified == 0.0),
        "p165_nonexpansive_polar_edge_semantics": bool(
            p165_polar_semantic_after <= p165_polar_semantic_before + 1.0e-12
            and p165_edge_semantic_after <= p165_edge_semantic_before + 1.0e-12
            and p165_polar_peak_after <= p165_polar_peak_before + 1.0e-12
            and p165_edge_components_after <= p165_edge_components_before + 2.0
        ),
        "p125_terminal_orogenic_semantic_land_consistency_accepted": bool(
            world.g(
                "terrain.last_p125_terminal_orogenic_semantic_land_consistency_accepted",
                0.0,
            )
            >= 1.0
        ),
        "p125_guard_reverted": bool(
            world.g("terrain.last_p125_guard_reverted", 0.0) >= 1.0
        ),
        "p125_candidate_cell_count": float(
            world.g("terrain.last_p125_candidate_cell_count", 0.0)
        ),
        "p125_candidate_area_fraction": float(
            world.g("terrain.last_p125_candidate_area_fraction", 0.0)
        ),
        "p125_cleared_hierarchy_cell_count": float(
            world.g("terrain.last_p125_cleared_hierarchy_cell_count", 0.0)
        ),
        "p125_cleared_hierarchy_area_fraction": float(
            world.g("terrain.last_p125_cleared_hierarchy_area_fraction", 0.0)
        ),
        "p125_cleared_spine_cell_count": float(
            world.g("terrain.last_p125_cleared_spine_cell_count", 0.0)
        ),
        "p125_cleared_spine_area_fraction": float(
            world.g("terrain.last_p125_cleared_spine_area_fraction", 0.0)
        ),
        "p125_cleared_halo_cell_count": float(
            world.g("terrain.last_p125_cleared_halo_cell_count", 0.0)
        ),
        "p125_cleared_halo_area_fraction": float(
            world.g("terrain.last_p125_cleared_halo_area_fraction", 0.0)
        ),
        "p125_cleared_apron_cell_count": float(
            world.g("terrain.last_p125_cleared_apron_cell_count", 0.0)
        ),
        "p125_cleared_apron_area_fraction": float(
            world.g("terrain.last_p125_cleared_apron_area_fraction", 0.0)
        ),
        "p125_submerged_hierarchy_cell_count_before": float(
            world.g("terrain.last_p125_submerged_hierarchy_cell_count_before", 0.0)
        ),
        "p125_submerged_hierarchy_cell_count_after": float(
            world.g("terrain.last_p125_submerged_hierarchy_cell_count_after", 0.0)
        ),
        "p125_submerged_spine_cell_count_before": float(
            world.g("terrain.last_p125_submerged_spine_cell_count_before", 0.0)
        ),
        "p125_submerged_spine_cell_count_after": float(
            world.g("terrain.last_p125_submerged_spine_cell_count_after", 0.0)
        ),
        "p125_submerged_halo_cell_count_before": float(
            world.g("terrain.last_p125_submerged_halo_cell_count_before", 0.0)
        ),
        "p125_submerged_halo_cell_count_after": float(
            world.g("terrain.last_p125_submerged_halo_cell_count_after", 0.0)
        ),
        "p125_submerged_apron_cell_count_before": float(
            world.g("terrain.last_p125_submerged_apron_cell_count_before", 0.0)
        ),
        "p125_submerged_apron_cell_count_after": float(
            world.g("terrain.last_p125_submerged_apron_cell_count_after", 0.0)
        ),
        "p125_polar_cleared_cell_count": float(
            world.g("terrain.last_p125_polar_cleared_cell_count", 0.0)
        ),
        "p125_edge_cleared_cell_count": float(
            world.g("terrain.last_p125_edge_cleared_cell_count", 0.0)
        ),
        "p125_land_semantic_cleared_cell_count": float(
            world.g("terrain.last_p125_land_semantic_cleared_cell_count", 0.0)
        ),
        "p125_no_submerged_terminal_orogenic_semantics": bool(
            np.count_nonzero(submerged_any_semantic) == 0
        ),
        "p125_land_semantics_preserved": bool(
            world.g("terrain.last_p125_land_semantic_cleared_cell_count", 0.0)
            == 0.0
        ),
    }


def _high_mountain_parent_mask(world: Any) -> np.ndarray:
    grid = world.grid
    parent = (
        _boundary_mask(world, "collision")
        | _boundary_mask(world, "suture")
        | _boundary_mask(world, "active_margin")
        | _object_mask(world, "tectonics.boundary_provinces", {
            "continent_continent_collision",
            "suture_zone",
            "continental_arc_margin",
            "continental_rift",
        })
        | _object_mask(world, "terrain.mountain_ranges", {
            "orogen",
            "old_subdued_orogen",
            "plateau",
            "arc_microcontinent",
            "active_orogen",
            "continental_arc",
            "old_orogen",
            "rift_shoulder",
            "plateau_margin",
        })
        | _object_mask(world, "terrain.continental_landforms", {
            "active_orogen",
            "old_orogen",
            "volcanic_arc",
            "rift_shoulder",
            "plateau_margin",
        })
    )
    orogenic_load = np.asarray(
        world.get_field("terrain.orogenic_load", np.zeros(grid.n)),
        dtype=np.float64,
    )
    if orogenic_load.shape == (grid.n,):
        parent |= orogenic_load >= 0.25
    orog_age = np.asarray(
        world.get_field("tectonics.orogeny_age_myr", np.full(grid.n, -1.0)),
        dtype=np.float64,
    )
    if orog_age.shape == (grid.n,):
        parent |= orog_age >= max(float(world.time_myr) - 900.0, 0.0)
    return parent


def _component_area_rows(grid: Any, mask: np.ndarray) -> list[dict[str, Any]]:
    rows = [
        {"cell_count": int(comp.size), "area": float(grid.cell_area[comp].sum())}
        for comp in _components(grid, mask)
    ]
    rows.sort(key=lambda row: (-float(row["area"]), -int(row["cell_count"])))
    return rows


def _arc_trench_adjacency_metrics(world: Any) -> dict[str, Any]:
    grid = world.grid
    trench = (
        _boundary_mask(world, "trench")
        | _object_mask(world, "terrain.margin_landforms", {"trench"})
    )
    arcs = (
        _object_mask(world, "terrain.arc_plume_landforms", {"island_arc"})
        | _object_mask(world, "terrain.margin_landforms", {"volcanic_arc"})
    )
    if trench.any():
        near_trench = _dilate(grid, trench, passes=5)
    else:
        near_trench = np.zeros(grid.n, dtype=bool)
    arc_area = float(grid.cell_area[arcs].sum())
    adjacent_area = float(grid.cell_area[arcs & near_trench].sum())
    return {
        "arc_cell_count": int(arcs.sum()),
        "trench_cell_count": int(trench.sum()),
        "arc_area_fraction_world": float(arc_area / max(float(grid.cell_area.sum()), 1.0e-12)),
        "active_margin_arc_trench_adjacency_fraction": float(
            adjacent_area / max(arc_area, 1.0e-12)),
    }


def _p1116_boundary_linework_metrics(
    world: Any,
    depth: np.ndarray,
    ocean: np.ndarray,
    depth_province: np.ndarray,
) -> dict[str, Any]:
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    raw_ridge = _boundary_mask(world, "ridge")
    object_ridge = _object_mask(world, "terrain.ocean_fabric", {"spreading_center"})
    depth_ridge = (np.asarray(depth_province, dtype=np.int64) == 5)
    if object_ridge.any():
        ridge = object_ridge
        ridge_source = "terrain_object"
    elif depth_ridge.any():
        ridge = depth_ridge
        ridge_source = "depth_province"
    else:
        ridge = raw_ridge
        ridge_source = "raw_boundary"

    raw_transform = _boundary_mask(world, "transform")
    object_transform = _object_mask(
        world,
        "terrain.ocean_fabric",
        {"transform_fault", "fracture_zone"},
    )
    if object_transform.any():
        transform = object_transform
        transform_source = "terrain_object"
    else:
        transform = raw_transform
        transform_source = "raw_boundary"

    raw_trench = _boundary_mask(world, "trench")
    object_trench = _object_mask(world, "terrain.margin_landforms", {"trench"})
    depth_trench = (np.asarray(depth_province, dtype=np.int64) == 6)
    if object_trench.any():
        trench = object_trench
        trench_source = "terrain_object"
    elif depth_trench.any():
        trench = depth_trench
        trench_source = "depth_province"
    else:
        trench = raw_trench
        trench_source = "raw_boundary"
    raw_convergent = (
        _boundary_mask(world, "convergent")
        | _boundary_mask(world, "subduction")
        | _boundary_mask(world, "active_margin")
        | _object_mask(world, "tectonics.boundary_provinces", {
            "ocean_ocean_subduction_trench",
            "island_arc_trench",
            "continental_arc_margin",
        })
        | _object_mask(world, "tectonics.boundary_objects", {
            "subduction",
            "active_margin",
        })
    )
    object_convergent_parent = _object_mask(
        world,
        "tectonics.boundary_objects",
        {"convergent_parent", "subduction_parent"},
    )
    if object_convergent_parent.any():
        convergent = object_convergent_parent
        convergent_source = "tectonic_parent_object"
    else:
        convergent = raw_convergent
        convergent_source = "raw_boundary"
    collision = (
        _boundary_mask(world, "collision")
        | _boundary_mask(world, "suture")
        | _object_mask(world, "tectonics.boundary_provinces", {
            "continent_continent_collision",
            "suture_zone",
        })
    )
    boundary_polyline_all = _object_mask(
        world,
        "tectonics.boundary_polylines",
        {
            "ridge_polyline",
            "transform_polyline",
            "trench_polyline",
            "suture_polyline",
            "active_margin_polyline",
            "passive_margin_polyline",
            "convergent_parent_polyline",
            "subduction_parent_polyline",
        },
    )
    boundary_polyline_ridge = _object_mask(
        world, "tectonics.boundary_polylines", {"ridge_polyline"})
    boundary_polyline_transform = _object_mask(
        world, "tectonics.boundary_polylines", {"transform_polyline"})
    boundary_polyline_trench = _object_mask(
        world, "tectonics.boundary_polylines", {"trench_polyline"})
    boundary_polyline_convergent_parent = _object_mask(
        world,
        "tectonics.boundary_polylines",
        {"convergent_parent_polyline", "subduction_parent_polyline"},
    )

    trench_near_subduction = (
        _dilate(grid, convergent | raw_convergent, passes=2)
        if (convergent | raw_convergent).any()
        else np.zeros(grid.n, dtype=bool)
    )
    ridge_near_transform = (
        _dilate(grid, transform, passes=2) if transform.any()
        else np.zeros(grid.n, dtype=bool)
    )
    deep_trench = trench & np.asarray(ocean, dtype=bool) & (np.asarray(depth) >= 5200.0)

    return {
        "schema": "aevum.p1116_boundary_linework.v3",
        "ridge_source": ridge_source,
        "transform_source": transform_source,
        "trench_source": trench_source,
        "ridge": _linework_stats(grid, ridge, total),
        "raw_ridge_boundary": _linework_stats(grid, raw_ridge, total),
        "object_ridge": _linework_stats(grid, object_ridge, total),
        "depth_ridge_province": _linework_stats(grid, depth_ridge, total),
        "transform": _linework_stats(grid, transform, total),
        "raw_transform_boundary": _linework_stats(grid, raw_transform, total),
        "object_transform": _linework_stats(grid, object_transform, total),
        "trench": _linework_stats(grid, trench, total),
        "raw_trench_boundary": _linework_stats(grid, raw_trench, total),
        "object_trench": _linework_stats(grid, object_trench, total),
        "depth_trench_province": _linework_stats(grid, depth_trench, total),
        "convergent_source": convergent_source,
        "convergent": _linework_stats(grid, convergent, total),
        "raw_convergent_boundary": _linework_stats(grid, raw_convergent, total),
        "object_convergent_parent": _linework_stats(
            grid, object_convergent_parent, total),
        "collision_suture": _linework_stats(grid, collision, total),
        "p129_boundary_polyline_all": _linework_stats(
            grid, boundary_polyline_all, total),
        "p129_boundary_polyline_ridge": _linework_stats(
            grid, boundary_polyline_ridge, total),
        "p129_boundary_polyline_transform": _linework_stats(
            grid, boundary_polyline_transform, total),
        "p129_boundary_polyline_trench": _linework_stats(
            grid, boundary_polyline_trench, total),
        "p129_boundary_polyline_convergent_parent": _linework_stats(
            grid, boundary_polyline_convergent_parent, total),
        "p129_boundary_polyline_object_count": float(
            world.g("tectonics.last_p129_object_count", 0.0)
        ),
        "p129_boundary_polyline_source_cell_count": float(
            world.g("tectonics.last_p129_source_cell_count", 0.0)
        ),
        "p129_boundary_polyline_axis_cell_count": float(
            world.g("tectonics.last_p129_axis_cell_count", 0.0)
        ),
        "p129_boundary_polyline_axis_source_coverage_fraction": float(
            world.g("tectonics.last_p129_axis_source_coverage_fraction", 0.0)
        ),
        "p129_boundary_polyline_mean_path_coverage_fraction": float(
            world.g("tectonics.last_p129_mean_path_coverage_fraction", 0.0)
        ),
        "p129_boundary_polyline_mean_directness": float(
            world.g("tectonics.last_p129_mean_directness", 0.0)
        ),
        "p129_boundary_polyline_mean_sinuosity": float(
            world.g("tectonics.last_p129_mean_sinuosity", 1.0)
        ),
        "p129_boundary_polyline_max_sinuosity": float(
            world.g("tectonics.last_p129_max_sinuosity", 1.0)
        ),
        "p129_boundary_polyline_source_junction_cell_count": float(
            world.g("tectonics.last_p129_source_junction_cell_count", 0.0)
        ),
        "p129_boundary_polyline_source_high_degree_cell_count": float(
            world.g("tectonics.last_p129_source_high_degree_cell_count", 0.0)
        ),
        "p129_boundary_polyline_ready": bool(
            world.g("tectonics.last_p129_polyline_ready", 0.0) >= 1.0
        ),
        "p131_gap_bridge_cell_count": float(
            world.g("tectonics.last_p131_gap_bridge_cell_count", 0.0)
        ),
        "p131_gap_bridge_object_count": float(
            world.g("tectonics.last_p131_gap_bridge_object_count", 0.0)
        ),
        "p131_gap_bridge_used": bool(
            world.g("tectonics.last_p131_gap_bridge_used", 0.0) >= 1.0
        ),
        "p129_ridge_polyline_count": float(
            world.g("tectonics.last_p129_ridge_polyline_count", 0.0)
        ),
        "p129_transform_polyline_count": float(
            world.g("tectonics.last_p129_transform_polyline_count", 0.0)
        ),
        "p129_trench_polyline_count": float(
            world.g("tectonics.last_p129_trench_polyline_count", 0.0)
        ),
        "p129_suture_polyline_count": float(
            world.g("tectonics.last_p129_suture_polyline_count", 0.0)
        ),
        "p129_active_margin_polyline_count": float(
            world.g("tectonics.last_p129_active_margin_polyline_count", 0.0)
        ),
        "p129_passive_margin_polyline_count": float(
            world.g("tectonics.last_p129_passive_margin_polyline_count", 0.0)
        ),
        "p129_convergent_parent_polyline_count": float(
            world.g("tectonics.last_p129_convergent_parent_polyline_count", 0.0)
        ),
        "p129_subduction_parent_polyline_count": float(
            world.g("tectonics.last_p129_subduction_parent_polyline_count", 0.0)
        ),
        "p131_gap_bridge_cell_count_by_kind": {
            kind: float(
                world.g(f"tectonics.last_p131_{kind}_gap_bridge_cell_count", 0.0)
            )
            for kind in (
                "ridge",
                "transform",
                "trench",
                "suture",
                "active_margin",
                "passive_margin",
                "convergent_parent",
                "subduction_parent",
            )
        },
        "p130_terrain_polyline_consumption": {
            name: {
                "polyline_cell_count": float(
                    world.g(f"terrain.last_p130_{name}_polyline_cell_count", 0.0)
                ),
                "polyline_added_cell_count": float(
                    world.g(
                        f"terrain.last_p130_{name}_polyline_added_cell_count",
                        0.0,
                    )
                ),
                "process_cell_count": float(
                    world.g(f"terrain.last_p130_{name}_process_cell_count", 0.0)
                ),
                "polyline_used": bool(
                    world.g(f"terrain.last_p130_{name}_polyline_used", 0.0) >= 1.0
                ),
            }
            for name in (
                "ridge",
                "trench",
                "active_margin",
                "passive_margin",
                "suture",
                "transform",
            )
        },
        "ridge_transform_adjacency_fraction": float(
            area[(ridge & ridge_near_transform)].sum() / max(float(area[ridge].sum()), 1.0e-12)
        ) if ridge.any() else 0.0,
        "trench_subduction_attachment_fraction": float(
            area[(trench & trench_near_subduction)].sum() / max(float(area[trench].sum()), 1.0e-12)
        ) if trench.any() else 0.0,
        "deep_trench_subduction_attachment_fraction": float(
            area[(deep_trench & trench_near_subduction)].sum()
            / max(float(area[deep_trench].sum()), 1.0e-12)
        ) if deep_trench.any() else 0.0,
    }


def _linework_stats(grid: Any, mask: np.ndarray, total_area: float) -> dict[str, Any]:
    mask = np.asarray(mask, dtype=bool)
    area = np.asarray(grid.cell_area, dtype=np.float64)
    mask_area = float(area[mask].sum())
    if not mask.any():
        return {
            "cell_count": 0,
            "area_fraction_world": 0.0,
            "component_count": 0,
            "largest_component_share": 0.0,
            "top3_component_share": 0.0,
            "isolated_cell_fraction": 0.0,
            "endpoint_cell_fraction": 0.0,
            "branch_cell_fraction": 0.0,
            "median_same_kind_neighbor_count": 0.0,
            "line_coherence_score": 0.0,
        }
    components = _components(grid, mask)
    comp_areas = sorted((float(area[comp].sum()) for comp in components), reverse=True)
    same = _same_kind_neighbor_count(grid, mask)
    isolated = mask & (same <= 0)
    endpoints = mask & (same == 1)
    branches = mask & (same >= 3)
    largest_share = float(comp_areas[0] / max(mask_area, 1.0e-12))
    top3_share = float(sum(comp_areas[:3]) / max(mask_area, 1.0e-12))
    isolated_fraction = float(area[isolated].sum() / max(mask_area, 1.0e-12))
    branch_fraction = float(area[branches].sum() / max(mask_area, 1.0e-12))
    line_coherence = float(np.clip(
        top3_share * (1.0 - isolated_fraction) * (1.0 - 0.5 * branch_fraction),
        0.0,
        1.0,
    ))
    return {
        "cell_count": int(mask.sum()),
        "area_fraction_world": float(mask_area / max(total_area, 1.0e-12)),
        "component_count": int(len(components)),
        "largest_component_share": largest_share,
        "top3_component_share": top3_share,
        "isolated_cell_fraction": isolated_fraction,
        "endpoint_cell_fraction": float(area[endpoints].sum() / max(mask_area, 1.0e-12)),
        "branch_cell_fraction": branch_fraction,
        "median_same_kind_neighbor_count": float(np.median(same[mask])),
        "line_coherence_score": line_coherence,
    }


def _p1116_oceanic_crust_structure_metrics(
    world: Any,
    depth: np.ndarray,
    ocean: np.ndarray,
    depth_province: np.ndarray,
) -> dict[str, Any]:
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    depth = np.asarray(depth, dtype=np.float64)
    ocean = np.asarray(ocean, dtype=bool)
    age = np.asarray(world.get_field("crust.age_myr", np.full(grid.n, np.nan)),
                     dtype=np.float64)
    crust_type = np.asarray(world.get_field("crust.type", np.zeros(grid.n)),
                            dtype=np.float64)
    if age.shape != (grid.n,):
        age = np.full(grid.n, np.nan, dtype=np.float64)
    if crust_type.shape != (grid.n,):
        crust_type = np.zeros(grid.n, dtype=np.float64)
    ridge = (
        _boundary_mask(world, "ridge")
        | _object_mask(world, "terrain.ocean_fabric", {"spreading_center"})
        | (np.asarray(depth_province, dtype=np.int64) == 5)
    ) & ocean
    trench = (
        _boundary_mask(world, "trench")
        | _object_mask(world, "terrain.margin_landforms", {"trench"})
        | (np.asarray(depth_province, dtype=np.int64) == 6)
    ) & ocean
    basin_ocean = ocean & np.isin(
        np.asarray(depth_province, dtype=np.int64),
        [4, 5],
    )
    open_ocean = basin_ocean & ~trench
    oceanic = open_ocean & (crust_type < 0.5) & np.isfinite(age)
    if not oceanic.any():
        return {
            "schema": "aevum.p1116_oceanic_crust_structure.v1",
            "sample_cell_count": 0,
            "age_depth_correlation": 0.0,
            "age_distance_from_ridge_correlation": 0.0,
            "old_minus_young_median_depth_m": 0.0,
            "far_minus_near_ridge_median_age_myr": 0.0,
            "young_near_ridge_fraction": 0.0,
            "old_far_from_ridge_fraction": 0.0,
            "ridge_median_depth_m": 0.0,
            "nonridge_open_ocean_median_depth_m": 0.0,
            "nonridge_minus_ridge_median_depth_m": 0.0,
            "trench_median_depth_m": 0.0,
            "nontrench_ocean_median_depth_m": 0.0,
            "trench_minus_nontrench_median_depth_m": 0.0,
        }

    distance = _distance_steps_to_mask(grid, ridge, max_steps=16)
    finite_distance = oceanic & np.isfinite(distance)
    values_age = age[oceanic]
    values_depth = depth[oceanic]
    age_depth_corr = _finite_corr(values_age, values_depth)
    age_distance_corr = _finite_corr(age[finite_distance], distance[finite_distance])

    young_cut = float(np.percentile(values_age, 25))
    old_cut = float(np.percentile(values_age, 75))
    young = oceanic & (age <= young_cut)
    old = oceanic & (age >= old_cut)
    near = finite_distance & (distance <= 2.0)
    far = finite_distance & (distance >= 5.0)
    ridge_depth = depth[ridge & ocean]
    nonridge = open_ocean & ~ridge
    trench_depth = depth[trench]
    nontrench_depth = depth[ocean & ~trench]

    return {
        "schema": "aevum.p1116_oceanic_crust_structure.v1",
        "sample_cell_count": int(oceanic.sum()),
        "age_depth_correlation": age_depth_corr,
        "age_distance_from_ridge_correlation": age_distance_corr,
        "old_minus_young_median_depth_m": float(
            np.median(depth[old]) - np.median(depth[young])
        ) if old.any() and young.any() else 0.0,
        "far_minus_near_ridge_median_age_myr": float(
            np.median(age[far]) - np.median(age[near])
        ) if far.any() and near.any() else 0.0,
        "young_near_ridge_fraction": float(
            area[young & near].sum() / max(float(area[young].sum()), 1.0e-12)
        ) if young.any() else 0.0,
        "old_far_from_ridge_fraction": float(
            area[old & far].sum() / max(float(area[old].sum()), 1.0e-12)
        ) if old.any() else 0.0,
        "ridge_median_depth_m": float(np.median(ridge_depth)) if ridge_depth.size else 0.0,
        "nonridge_open_ocean_median_depth_m": (
            float(np.median(depth[nonridge])) if nonridge.any() else 0.0
        ),
        "nonridge_minus_ridge_median_depth_m": (
            float(np.median(depth[nonridge]) - np.median(ridge_depth))
            if nonridge.any() and ridge_depth.size else 0.0
        ),
        "trench_median_depth_m": float(np.median(trench_depth)) if trench_depth.size else 0.0,
        "nontrench_ocean_median_depth_m": (
            float(np.median(nontrench_depth)) if nontrench_depth.size else 0.0
        ),
        "trench_minus_nontrench_median_depth_m": (
            float(np.median(trench_depth) - np.median(nontrench_depth))
            if trench_depth.size and nontrench_depth.size else 0.0
        ),
    }


def _ocean_feature_metrics(world: Any, depth: np.ndarray, ocean: np.ndarray,
                           depth_province: np.ndarray) -> dict[str, Any]:
    grid = world.grid
    ocean_area = max(float(grid.cell_area[ocean].sum()), 1.0e-12)
    total_area = max(float(grid.cell_area.sum()), 1.0e-12)
    shelf_width = np.asarray(
        world.get_field("ocean.shelf_width", np.zeros(grid.n)),
        dtype=np.float64,
    )
    if shelf_width.shape != (grid.n,):
        shelf_width = np.zeros(grid.n, dtype=np.float64)
    trench_context = (
        _boundary_mask(world, "trench")
        | _object_mask(world, "terrain.margin_landforms", {"trench"})
        | (depth_province == 6)
    )
    ridge_province = depth_province == 5
    raw_ridge = _boundary_mask(world, "ridge")
    object_ridge = _object_mask(world, "terrain.ocean_fabric", {"spreading_center"})
    ridge_core = object_ridge if object_ridge.any() else raw_ridge
    transform_fracture = _object_mask(
        world,
        "terrain.ocean_fabric",
        {"transform_fault", "fracture_zone"},
    )
    abyssal_plain = _object_mask(world, "terrain.ocean_fabric", {"abyssal_plain"})
    abyssal_hill = _object_mask(world, "terrain.ocean_fabric", {"abyssal_hill"})
    back_arc_basin = _object_mask(
        world,
        "terrain.arc_plume_landforms",
        {"back_arc_basin"},
    )
    seamount_chain = _object_mask(
        world,
        "terrain.arc_plume_landforms",
        {"seamount_chain", "hotspot_track"},
    )
    oceanic_plateau = _object_mask(
        world,
        "terrain.arc_plume_landforms",
        {"oceanic_plateau", "large_igneous_province"},
    )
    microcontinent = _object_mask(
        world,
        "terrain.arc_plume_landforms",
        {"microcontinent"},
    )
    island_arc = _object_mask(world, "terrain.arc_plume_landforms", {"island_arc"})
    narrow_relief_core = (
        ridge_core
        | transform_fracture
        | trench_context
        | seamount_chain
        | island_arc
    ) & ocean
    narrow_relief_halo = (
        _dilate(grid, narrow_relief_core, passes=1) & ocean
        if narrow_relief_core.any() else narrow_relief_core
    )
    broad_relief_core = (oceanic_plateau | microcontinent) & ocean
    object_backed_ocean_relief = (narrow_relief_core | broad_relief_core) & ocean
    object_backed_ocean_support = (narrow_relief_halo | broad_relief_core) & ocean

    land = ~ocean
    emergent_object_core = land & (
        seamount_chain | oceanic_plateau | microcontinent | island_arc
    )
    domain = np.asarray(world.get_field("crust.domain", np.zeros(grid.n)),
                        dtype=np.float64)
    origin = np.asarray(world.get_field("crust.origin", np.zeros(grid.n)),
                        dtype=np.float64)
    terrane_id = np.asarray(
        world.get_field("tectonics.terrane_id", np.full(grid.n, -1.0)),
        dtype=np.float64,
    )
    if domain.shape == (grid.n,) and origin.shape == (grid.n,) and terrane_id.shape == (grid.n,):
        emergent_object_core &= (
            (terrane_id >= 0.0)
            | (domain == DOMAIN_LIP)
            | (domain == DOMAIN_ACCRETED_TERRANE)
            | (origin == ORIGIN_ARC)
            | (origin == ORIGIN_PLUME_IMPACT)
        )
    emerged_object_backed_islands = np.zeros(grid.n, dtype=bool)
    if emergent_object_core.any():
        for comp in _components(grid, land):
            comp_mask = np.zeros(grid.n, dtype=bool)
            comp_mask[comp] = True
            comp_area = max(float(grid.cell_area[comp].sum()), 1.0e-12)
            comp_world_fraction = float(comp_area / total_area)
            object_share = float(
                grid.cell_area[comp_mask & emergent_object_core].sum() / comp_area
            )
            if comp_world_fraction < 0.03 or (
                comp_world_fraction < 0.06 and object_share >= 0.50
            ):
                emerged_object_backed_islands |= comp_mask & emergent_object_core
    parented_chain_count = 0
    for obj in world.objects.get("terrain.arc_plume_landforms", []):
        if str(obj.get("kind", "")) not in {"seamount_chain", "hotspot_track", "island_arc"}:
            continue
        parented = bool(obj.get("parent_tectonic_object_ids") or obj.get("parent_continent_ids"))
        ocean_fraction = float(obj.get("ocean_fraction", 0.0))
        if parented and ocean_fraction >= 0.5:
            parented_chain_count += 1
    deep = ocean & (depth >= 6000.0)
    deep_trench = deep & (_dilate(grid, trench_context, passes=1) if trench_context.any()
                          else trench_context)
    open_ocean = ocean & ((shelf_width <= 0.0) | (shelf_width >= 5.0))
    shallow_open_ocean = open_ocean & (depth > 0.0) & (depth < 2600.0)
    unsupported_shallow = shallow_open_ocean & ~object_backed_ocean_support
    object_backed_shallow = shallow_open_ocean & object_backed_ocean_support

    def area_fraction(mask: np.ndarray, denom: float = ocean_area) -> float:
        return float(grid.cell_area[mask].sum() / max(denom, 1.0e-12))

    return {
        "parented_oceanic_island_chain_count": int(parented_chain_count),
        "deep_ocean_fraction_below_6000m": float(grid.cell_area[deep].sum() / ocean_area),
        "deep_trench_fraction_below_6000m": float(
            grid.cell_area[deep_trench].sum() / ocean_area),
        "deep_nontrench_fraction_below_6000m": float(
            grid.cell_area[deep & ~deep_trench].sum() / ocean_area),
        "p1115_open_ocean_shoal_fraction": area_fraction(shallow_open_ocean),
        "p1115_object_backed_open_ocean_shoal_fraction": area_fraction(
            object_backed_shallow),
        "p1115_unsupported_open_ocean_shoal_fraction": area_fraction(
            unsupported_shallow),
        "p1115_object_backed_ocean_relief_fraction": area_fraction(
            object_backed_ocean_relief),
        "p1115_object_backed_ocean_support_fraction": area_fraction(
            object_backed_ocean_support),
        "p1115_emerged_object_backed_island_fraction_world": area_fraction(
            emerged_object_backed_islands,
            total_area,
        ),
        "p1115_ridge_area_fraction_ocean": area_fraction(ridge_core & ocean),
        "p1115_ridge_province_area_fraction_ocean": area_fraction(
            ridge_province & ocean),
        "p1115_transform_fracture_area_fraction_ocean": area_fraction(
            transform_fracture & ocean),
        "p1115_trench_area_fraction_ocean": area_fraction(trench_context & ocean),
        "p1115_abyssal_plain_area_fraction_ocean": area_fraction(
            abyssal_plain & ocean),
        "p1115_abyssal_hill_area_fraction_ocean": area_fraction(
            abyssal_hill & ocean),
        "p1115_back_arc_basin_area_fraction_ocean": area_fraction(
            back_arc_basin & ocean),
        "p1115_seamount_chain_area_fraction_ocean": area_fraction(
            seamount_chain & ocean),
        "p1115_oceanic_plateau_area_fraction_ocean": area_fraction(
            oceanic_plateau & ocean),
        "p1115_microcontinent_area_fraction_ocean": area_fraction(
            microcontinent & ocean),
        "p1115_island_arc_area_fraction_ocean": area_fraction(island_arc & ocean),
        "p1115_terrain_unsupported_shoal_deepened_fraction": float(
            world.g(
                "terrain.last_p1115_unsupported_open_ocean_shoal_deepened_fraction",
                0.0,
            )
        ),
        "p1115_terrain_object_backed_shoal_preserved_fraction": float(
            world.g(
                "terrain.last_p1115_object_backed_open_ocean_shoal_area_fraction",
                0.0,
            )
        ),
        "p1115_final_ocean_floor_adjusted_area_fraction_world": float(
            world.g("terrain.last_p1115_ocean_floor_adjusted_area_fraction", 0.0)
        ),
        "p1115_final_ocean_floor_land_mask_preserved": bool(
            world.g("terrain.last_p1115_ocean_floor_land_mask_preserved", 1.0) >= 1.0
        ),
        "p111625_abyssal_depth_calibration_used": bool(
            world.g("terrain.last_p111625_abyssal_depth_calibration_used", 0.0)
            >= 1.0
        ),
        "p111625_abyssal_depth_calibration_candidate_area_fraction": float(
            world.g(
                "terrain.last_p111625_abyssal_depth_calibration_candidate_area_fraction",
                0.0,
            )
        ),
        "p111625_abyssal_depth_calibration_adjusted_area_fraction": float(
            world.g(
                "terrain.last_p111625_abyssal_depth_calibration_adjusted_area_fraction",
                0.0,
            )
        ),
        "p111625_slope_min_depth_adjusted_area_fraction": float(
            world.g(
                "terrain.last_p111625_slope_min_depth_adjusted_area_fraction",
                0.0,
            )
        ),
        "p111625_rise_min_depth_adjusted_area_fraction": float(
            world.g(
                "terrain.last_p111625_rise_min_depth_adjusted_area_fraction",
                0.0,
            )
        ),
        "p111625_far_ocean_age_depth_adjusted_area_fraction": float(
            world.g(
                "terrain.last_p111625_far_ocean_age_depth_adjusted_area_fraction",
                0.0,
            )
        ),
        "p111625_ocean_mean_depth_before_m": float(
            world.g("terrain.last_p111625_ocean_mean_depth_before_m", 0.0)
        ),
        "p111625_ocean_mean_depth_after_m": float(
            world.g("terrain.last_p111625_ocean_mean_depth_after_m", 0.0)
        ),
        "p111625_ocean_p50_depth_before_m": float(
            world.g("terrain.last_p111625_ocean_p50_depth_before_m", 0.0)
        ),
        "p111625_ocean_p50_depth_after_m": float(
            world.g("terrain.last_p111625_ocean_p50_depth_after_m", 0.0)
        ),
        "p111626_ridge_flank_narrowing_used": bool(
            world.g("terrain.last_p111626_ridge_flank_narrowing_used", 0.0)
            >= 1.0
        ),
        "p111626_ridge_axis_protected_area_fraction": float(
            world.g(
                "terrain.last_p111626_ridge_axis_protected_area_fraction",
                0.0,
            )
        ),
        "p111626_ridge_flank_candidate_area_fraction": float(
            world.g(
                "terrain.last_p111626_ridge_flank_candidate_area_fraction",
                0.0,
            )
        ),
        "p111626_ridge_flank_adjusted_area_fraction": float(
            world.g(
                "terrain.last_p111626_ridge_flank_adjusted_area_fraction",
                0.0,
            )
        ),
        "p111626_ridge_open_shoal_before_fraction": float(
            world.g(
                "terrain.last_p111626_ridge_open_shoal_before_fraction",
                0.0,
            )
        ),
        "p111626_ridge_open_shoal_after_fraction": float(
            world.g(
                "terrain.last_p111626_ridge_open_shoal_after_fraction",
                0.0,
            )
        ),
        "p111627_depth_ridge_province_area_fraction_before": float(
            world.g(
                "terrain.last_p111627_depth_ridge_province_area_fraction_before",
                0.0,
            )
        ),
        "p111627_depth_ridge_axis_area_fraction": float(
            world.g(
                "terrain.last_p111627_depth_ridge_axis_area_fraction",
                0.0,
            )
        ),
        "p111627_depth_ridge_province_area_fraction_after": float(
            world.g(
                "terrain.last_p111627_depth_ridge_province_area_fraction_after",
                0.0,
            )
        ),
        "p111627_spreading_center_source_area_fraction": float(
            world.g(
                "terrain.last_p111627_spreading_center_source_area_fraction",
                0.0,
            )
        ),
        "p111627_spreading_center_axis_area_fraction": float(
            world.g(
                "terrain.last_p111627_spreading_center_axis_area_fraction",
                0.0,
            )
        ),
        "p111628_final_shoal_floor_used": bool(
            world.g("terrain.last_p111628_final_shoal_floor_used", 0.0) >= 1.0
        ),
        "p111628_final_shoal_candidate_area_fraction": float(
            world.g(
                "terrain.last_p111628_final_shoal_candidate_area_fraction",
                0.0,
            )
        ),
        "p111628_final_shoal_adjusted_area_fraction": float(
            world.g(
                "terrain.last_p111628_final_shoal_adjusted_area_fraction",
                0.0,
            )
        ),
        "p111628_ridge_shoulder_candidate_area_fraction": float(
            world.g(
                "terrain.last_p111628_ridge_shoulder_candidate_area_fraction",
                0.0,
            )
        ),
        "p111628_ridge_shoulder_adjusted_area_fraction": float(
            world.g(
                "terrain.last_p111628_ridge_shoulder_adjusted_area_fraction",
                0.0,
            )
        ),
        "p111628_final_semantic_cleanup_used": bool(
            world.g("terrain.last_p111628_final_semantic_cleanup_used", 0.0) >= 1.0
        ),
        "p111628_final_semantic_cleanup_candidate_area_fraction": float(
            world.g(
                "terrain.last_p111628_final_semantic_cleanup_candidate_area_fraction",
                0.0,
            )
        ),
        "p111628_final_semantic_cleanup_adjusted_area_fraction": float(
            world.g(
                "terrain.last_p111628_final_semantic_cleanup_adjusted_area_fraction",
                0.0,
            )
        ),
        "p111628_final_semantic_cleanup_ridge_area_fraction": float(
            world.g(
                "terrain.last_p111628_final_semantic_cleanup_ridge_area_fraction",
                0.0,
            )
        ),
        "p111628_final_semantic_cleanup_abyss_area_fraction": float(
            world.g(
                "terrain.last_p111628_final_semantic_cleanup_abyss_area_fraction",
                0.0,
            )
        ),
        "p111630_object_halo_shoal_cleanup_used": bool(
            world.g("terrain.last_p111630_object_halo_shoal_cleanup_used", 0.0)
            >= 1.0
        ),
        "p111630_object_halo_shoal_cleanup_candidate_area_fraction": float(
            world.g(
                "terrain.last_p111630_object_halo_shoal_cleanup_candidate_area_fraction",
                0.0,
            )
        ),
        "p111630_object_halo_shoal_cleanup_adjusted_area_fraction": float(
            world.g(
                "terrain.last_p111630_object_halo_shoal_cleanup_adjusted_area_fraction",
                0.0,
            )
        ),
        "p111630_object_backed_shoal_narrow_core_fraction": area_fraction(
            shallow_open_ocean & narrow_relief_core),
        "p111630_object_backed_shoal_halo_only_fraction": area_fraction(
            shallow_open_ocean
            & narrow_relief_halo
            & ~narrow_relief_core
            & ~broad_relief_core),
        "p111630_object_backed_shoal_broad_core_fraction": area_fraction(
            shallow_open_ocean & broad_relief_core),
        "p111630_shoal_microcontinent_fraction": area_fraction(
            shallow_open_ocean & microcontinent),
        "p111630_shoal_oceanic_plateau_fraction": area_fraction(
            shallow_open_ocean & oceanic_plateau),
        "p111630_shoal_seamount_chain_fraction": area_fraction(
            shallow_open_ocean & seamount_chain),
        "p111630_shoal_island_arc_fraction": area_fraction(
            shallow_open_ocean & island_arc),
        "p111631_microcontinent_raw_area_fraction_ocean": float(
            world.g("terrain.last_p111631_microcontinent_raw_area_fraction_ocean", 0.0)
        ),
        "p111631_microcontinent_area_fraction_ocean": float(
            world.g("terrain.last_p111631_microcontinent_area_fraction_ocean", 0.0)
        ),
        "p111631_microcontinent_removed_area_fraction_ocean": float(
            world.g(
                "terrain.last_p111631_microcontinent_removed_area_fraction_ocean",
                0.0,
            )
        ),
        "p111631_microcontinent_raw_component_count": float(
            world.g("terrain.last_p111631_microcontinent_raw_component_count", 0.0)
        ),
        "p111631_microcontinent_component_count": float(
            world.g("terrain.last_p111631_microcontinent_component_count", 0.0)
        ),
        "p111631_microcontinent_relief_raw_area_fraction_ocean": float(
            world.g(
                "terrain.last_p111631_microcontinent_relief_raw_area_fraction_ocean",
                0.0,
            )
        ),
        "p111631_microcontinent_relief_area_fraction_ocean": float(
            world.g(
                "terrain.last_p111631_microcontinent_relief_area_fraction_ocean",
                0.0,
            )
        ),
        "p111631_microcontinent_relief_removed_area_fraction_ocean": float(
            world.g(
                "terrain.last_p111631_microcontinent_relief_removed_area_fraction_ocean",
                0.0,
            )
        ),
        "p111631_rejected_microcontinent_shoal_deepened_fraction_ocean": float(
            world.g(
                "terrain.last_p111631_rejected_microcontinent_shoal_deepened_fraction_ocean",
                0.0,
            )
        ),
        "p111628_shoal_depth_abyss_fraction": area_fraction(
            shallow_open_ocean & (depth_province == 4)),
        "p111628_unsupported_shoal_depth_abyss_fraction": area_fraction(
            unsupported_shallow & (depth_province == 4)),
        "p111628_shoal_depth_ridge_fraction": area_fraction(
            shallow_open_ocean & ridge_province),
        "p111628_unsupported_shoal_depth_ridge_fraction": area_fraction(
            unsupported_shallow & ridge_province),
        "p111628_shoal_back_arc_basin_fraction": area_fraction(
            shallow_open_ocean & back_arc_basin),
        "p111628_unsupported_shoal_back_arc_basin_fraction": area_fraction(
            unsupported_shallow & back_arc_basin),
        "p111628_shoal_abyssal_plain_fraction": area_fraction(
            shallow_open_ocean & abyssal_plain),
        "p111628_unsupported_shoal_abyssal_plain_fraction": area_fraction(
            unsupported_shallow & abyssal_plain),
        "p111628_shoal_abyssal_hill_fraction": area_fraction(
            shallow_open_ocean & abyssal_hill),
        "p111628_unsupported_shoal_abyssal_hill_fraction": area_fraction(
            unsupported_shallow & abyssal_hill),
    }


def _p110a_modern_planform_metrics(world: Any, land: np.ndarray,
                                   ocean: np.ndarray) -> dict[str, Any]:
    """Measure modern Earth-like land/ocean planform on the sphere."""
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    time_myr = float(getattr(world, "time_myr", 0.0))
    modern_planform_applicable = bool(
        time_myr >= P110A_MODERN_PLANFORM_APPLICABLE_MYR
    )
    land_topology = _mask_component_topology(
        grid,
        land,
        major_share_threshold=P110A_MAJOR_LAND_COMPONENT_SHARE,
    )
    ocean_connectivity = _mask_component_topology(
        grid,
        ocean,
        major_share_threshold=0.01,
    )
    ocean_basins = _ocean_basin_planform_metrics(world, ocean)
    ocean_gateway_topology = _ocean_gateway_topology_metrics(world, ocean)
    archetype = _p110b_final_state_archetype_metrics(world)
    seaway = _seaway_effectiveness_metrics(world)
    payback = _land_payback_bias_metrics(world)
    lineage = _p110b_lineage_survival_metrics(
        world,
        land,
        major_share_threshold=P110A_MAJOR_LAND_COMPONENT_SHARE,
    )
    terminal_supercontinent = _p110b_terminal_supercontinent_diagnostics(
        world,
        land,
        major_share_threshold=P110A_MAJOR_LAND_COMPONENT_SHARE,
    )

    land_fraction = float(area[land].sum() / total)
    largest_land_share = float(land_topology["largest_component_share_of_mask"])
    major_land_count = int(land_topology["major_component_count"])
    top_land = [
        float(x) for x in land_topology["top_component_shares_of_mask"]
    ]
    second_share = top_land[1] if len(top_land) > 1 else 0.0
    third_share = top_land[2] if len(top_land) > 2 else 0.0
    fourth_share = top_land[3] if len(top_land) > 3 else 0.0
    major_basin_count = int(ocean_basins["major_basin_count"])
    disconnected_ocean_share = max(
        0.0,
        1.0 - float(ocean_connectivity["largest_component_share_of_mask"]),
    )
    island_arc_count = sum(
        1 for obj in world.objects.get("terrain.arc_plume_landforms", [])
        if isinstance(obj, dict) and str(obj.get("kind", "")) == "island_arc"
    )
    microcontinent_count = sum(
        1 for obj in world.objects.get("terrain.arc_plume_landforms", [])
        if isinstance(obj, dict) and str(obj.get("kind", "")) == "microcontinent"
    )
    archetype_name = str(archetype.get("name", ""))
    post_breakup_three_major = (
        archetype_name == "post_supercontinent_breakup"
        and major_land_count >= 3
        and second_share >= 0.16
        and third_share >= 0.16
    )
    archipelago_two_major = (
        archetype_name == "archipelago_active_margin"
        and major_land_count >= 2
        and second_share >= 0.20
        and island_arc_count >= 12
        and microcontinent_count >= 6
    )
    modern_balanced_three_major = (
        archetype_name == "modern_multipolar"
        and major_land_count >= 3
        and largest_land_share <= 0.54
        and second_share >= 0.20
        and third_share >= 0.16
        and not terminal_supercontinent["terminal_supercontinent_like"]
        and not terminal_supercontinent["modern_multipolar_overconnected"]
        and (
            int(terminal_supercontinent[
                "largest_land_significant_continent_domain_count"
            ]) >= 2
            or int(lineage["independent_primary_lineage_count"]) >= 2
        )
    )
    split_secondary_pair = (
        major_land_count >= 4
        and second_share >= 0.18
        and third_share >= 0.040
        and fourth_share >= 0.030
        and (third_share + fourth_share) >= 0.080
    )

    out_of_envelope: list[str] = []
    deferred_out_of_envelope: list[str] = []
    warnings: list[str] = []

    def add_modern_planform_failure(flag: str) -> None:
        if modern_planform_applicable:
            out_of_envelope.append(flag)
        else:
            deferred_out_of_envelope.append(flag)
            warnings.append(f"immature_modern_planform_{flag}")

    if land_fraction < 0.22 or land_fraction > 0.32:
        add_modern_planform_failure("land_fraction")
    if major_land_count < 4:
        if post_breakup_three_major:
            warnings.append("major_land_component_count_low_p110b_post_breakup")
        elif archipelago_two_major:
            warnings.append("major_land_component_count_low_p110b_archipelago")
        elif modern_balanced_three_major:
            warnings.append(
                "major_land_component_count_low_p110b_modern_balanced_three_major"
            )
        else:
            add_modern_planform_failure("major_land_component_count_low")
    elif major_land_count > 7:
        warnings.append("major_land_component_count_high")
    if largest_land_share > 0.65:
        add_modern_planform_failure("largest_land_component_share_hard")
    elif largest_land_share > 0.60:
        warnings.append("largest_land_component_share_soft")
    elif largest_land_share > 0.56:
        warnings.append("largest_land_component_share_above_preferred")
    if land_fraction >= 0.20 and sum(1 for share in top_land[1:] if share >= 0.08) < 2:
        if split_secondary_pair:
            warnings.append("nonlargest_major_component_pair_split")
        elif archipelago_two_major:
            warnings.append("nonlargest_major_component_count_low_p110b_archipelago")
        else:
            add_modern_planform_failure("nonlargest_major_component_count_low")
    if major_basin_count < 2:
        add_modern_planform_failure("major_ocean_basin_count_low")
    if disconnected_ocean_share > 0.20:
        add_modern_planform_failure("closed_ocean_ring_score_hard")
    elif disconnected_ocean_share > P110A_DISCONNECTED_OCEAN_WARNING_SHARE:
        warnings.append("closed_ocean_ring_score_soft")
    if deferred_out_of_envelope:
        warnings.insert(0, "modern_planform_gate_deferred_until_3500myr")
    if (
        lineage["continental_major_component_count"] > 0
        and lineage["continental_lineage_supported_major_component_count"]
        < lineage["continental_major_component_count"]
    ):
        warnings.append("p110b_lineage_support_incomplete")
    if (
        lineage["major_component_count"] > 0
        and lineage["provenance_supported_major_component_count"]
        < lineage["major_component_count"]
    ):
        warnings.append("p110b_landform_provenance_incomplete")
    if (
        lineage["continental_major_component_count"] >= 4
        and lineage["independent_primary_lineage_count"] < 3
    ):
        warnings.append("p110b_independent_lineage_count_low")
    if (
        archipelago_two_major
        and lineage["independent_primary_lineage_count"] < 2
    ):
        warnings.append("p110b_archipelago_independent_lineage_count_low")
    if terminal_supercontinent["terminal_supercontinent_like"]:
        warnings.append("p110b_terminal_supercontinent_like")
    if terminal_supercontinent["modern_multipolar_overconnected"]:
        warnings.append("p110b_modern_multipolar_overconnected")

    return {
        "schema": "aevum.p110a_modern_planform.v1",
        "time_myr": time_myr,
        "modern_planform_applicable": modern_planform_applicable,
        "modern_planform_applicable_after_myr":
            P110A_MODERN_PLANFORM_APPLICABLE_MYR,
        "major_land_component_share_threshold": P110A_MAJOR_LAND_COMPONENT_SHARE,
        "major_ocean_basin_share_threshold": P110A_MAJOR_OCEAN_BASIN_SHARE,
        "land": land_topology,
        "ocean_connectivity": {
            **ocean_connectivity,
            "closed_ocean_ring_score": float(disconnected_ocean_share),
            "disconnected_ocean_component_share": float(disconnected_ocean_share),
        },
        "ocean_basins": ocean_basins,
        "ocean_gateway_topology": ocean_gateway_topology,
        "p110b_final_state_archetype": archetype,
        "p110b_lineage_survival": lineage,
        "p110b_terminal_supercontinent_diagnostics": terminal_supercontinent,
        "seaway_cut_effectiveness": seaway,
        "land_payback_largest_component_bias": payback,
        "summary": {
            "land_fraction": land_fraction,
            "time_myr": time_myr,
            "modern_planform_applicable": modern_planform_applicable,
            "major_land_component_count": major_land_count,
            "largest_land_component_share": largest_land_share,
            "second_land_component_share": float(second_share),
            "third_land_component_share": float(third_share),
            "major_ocean_basin_count": major_basin_count,
            "terminal_ocean_gateway_count": int(
                ocean_gateway_topology["terminal_gateway_count"]),
            "terminal_interbasin_ocean_gateway_count": int(
                ocean_gateway_topology["terminal_interbasin_gateway_count"]),
            "terminal_phase_backed_ocean_gateway_count": int(
                ocean_gateway_topology["terminal_phase_backed_gateway_count"]),
            "terminal_ocean_gateway_system_count": int(
                ocean_gateway_topology["terminal_gateway_system_count"]),
            "terminal_interbasin_ocean_gateway_system_count": int(
                ocean_gateway_topology[
                    "terminal_interbasin_gateway_system_count"]),
            "terminal_phase_backed_ocean_gateway_system_count": int(
                ocean_gateway_topology[
                    "terminal_phase_backed_gateway_system_count"]),
            "ocean_gateway_fragment_to_system_ratio": float(
                ocean_gateway_topology["gateway_fragment_to_system_ratio"]),
            "tectonic_ocean_gateway_count": int(
                ocean_gateway_topology["tectonic_gateway_count"]),
            "restricted_ocean_fraction": float(
                ocean_gateway_topology["restricted_ocean_fraction"]),
            "unbacked_major_disconnected_ocean_component_count": int(
                ocean_gateway_topology[
                    "unbacked_major_disconnected_ocean_component_count"]),
            "largest_ocean_basin_share": float(
                ocean_basins["largest_basin_share_of_ocean"]),
            "closed_ocean_ring_score": float(disconnected_ocean_share),
            "lineage_supported_major_component_count": int(
                lineage["lineage_supported_major_component_count"]),
            "continental_major_component_count": int(
                lineage["continental_major_component_count"]),
            "continental_lineage_supported_major_component_count": int(
                lineage["continental_lineage_supported_major_component_count"]),
            "oceanic_landform_major_component_count": int(
                lineage["oceanic_landform_major_component_count"]),
            "provenance_supported_major_component_count": int(
                lineage["provenance_supported_major_component_count"]),
            "independent_primary_lineage_count": int(
                lineage["independent_primary_lineage_count"]),
            "terminal_supercontinent_score": float(
                terminal_supercontinent["terminal_supercontinent_score"]),
            "terminal_supercontinent_like": bool(
                terminal_supercontinent["terminal_supercontinent_like"]),
            "largest_land_significant_continent_domain_count": int(
                terminal_supercontinent[
                    "largest_land_significant_continent_domain_count"]),
            "largest_land_effective_continent_domain_count": float(
                terminal_supercontinent[
                    "largest_land_effective_continent_domain_count"]),
            "largest_land_robust_piece_count_after_neck_removal": int(
                terminal_supercontinent[
                    "largest_land_robust_piece_count_after_neck_removal"]),
            "largest_land_bridge_candidate_fraction": float(
                terminal_supercontinent[
                    "largest_land_bridge_candidate_fraction"]),
            "internal_domain_seaway_opening_count": int(
                terminal_supercontinent["internal_domain_seaway_opening_count"]),
            "internal_domain_seaway_area_fraction_world": float(
                terminal_supercontinent[
                    "internal_domain_seaway_area_fraction_world"]),
            "internal_domain_boundary_count": int(
                terminal_supercontinent["internal_domain_boundary_count"]),
            "internal_domain_boundary_area_fraction_world": float(
                terminal_supercontinent[
                    "internal_domain_boundary_area_fraction_world"]),
        },
        "out_of_envelope": out_of_envelope,
        "deferred_out_of_envelope": deferred_out_of_envelope,
        "warning_flags": warnings,
        "out_of_envelope_count": int(len(out_of_envelope)),
        "deferred_out_of_envelope_count": int(len(deferred_out_of_envelope)),
        "warning_count": int(len(warnings)),
        "within_p110a_modern_planform_envelope": bool(len(out_of_envelope) == 0),
    }


def _p110b_terminal_supercontinent_diagnostics(
    world: Any,
    land: np.ndarray,
    *,
    major_share_threshold: float,
) -> dict[str, Any]:
    """Diagnose visually over-connected terminal landmasses.

    This intentionally separates "large connected landmass" from "bad
    supercontinent".  Modern Earth has a very large Afro-Eurasian connected
    component, so the diagnostic also looks for inherited continent domains and
    whether low, narrow land bridges are carrying most of the connectivity.
    """
    grid = world.grid
    n = int(grid.n)
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    land = np.asarray(land, dtype=bool)
    land_area = max(float(area[land].sum()), 1.0e-12)

    def int_field(name: str, default: float) -> np.ndarray:
        arr = np.asarray(world.get_field(name, default), dtype=np.float64)
        if arr.shape != (n,):
            arr = np.full(n, default, dtype=np.float64)
        return arr.astype(int)

    def float_field(name: str, default: float) -> np.ndarray:
        arr = np.asarray(world.get_field(name, default), dtype=np.float64)
        if arr.shape != (n,):
            arr = np.full(n, default, dtype=np.float64)
        return np.nan_to_num(arr, nan=default, posinf=default, neginf=default)

    components = _components(grid, land)
    components.sort(key=lambda comp: float(area[comp].sum()), reverse=True)
    if not components:
        return {
            "schema": "aevum.p110b_terminal_supercontinent_diagnostics.v1",
            "major_share_threshold": float(major_share_threshold),
            "terminal_supercontinent_score": 0.0,
            "terminal_supercontinent_like": False,
            "modern_multipolar_overconnected": False,
            "land_component_count": 0,
            "major_land_component_count": 0,
            "nonlargest_major_land_component_count": 0,
            "largest_land_share": 0.0,
            "second_land_share": 0.0,
            "third_land_share": 0.0,
            "largest_land_dominant_continent_share": 0.0,
            "largest_land_significant_continent_domain_count": 0,
            "largest_land_effective_continent_domain_count": 0.0,
            "largest_land_bridge_candidate_fraction": 0.0,
            "largest_land_robust_piece_count_after_neck_removal": 0,
            "largest_land_largest_piece_share_after_neck_removal": 0.0,
            "largest_land_neck_split_gain": 0,
            "all_land_robust_major_piece_count_after_neck_removal": 0,
            "internal_domain_seaway_opening_count": 0,
            "internal_domain_seaway_area_fraction_world": 0.0,
            "internal_domain_seaway_object_count": 0,
            "internal_domain_boundary_count": 0,
            "internal_domain_boundary_area_fraction_world": 0.0,
            "internal_domain_boundary_object_count": 0,
            "component_rows_top5": [],
            "continent_domain_rows_top8": [],
        }

    elevation = float_field("terrain.elevation_m", 0.0)
    sea_level = float(world.sea_level)
    rel = elevation - sea_level
    continent_id = int_field("tectonics.continent_id", -1.0)
    terrain_detail = int_field("terrain.continental_detail", -1.0)

    component_rows: list[dict[str, Any]] = []
    shares: list[float] = []
    for idx, comp in enumerate(components):
        comp_area = float(area[comp].sum())
        share = float(comp_area / land_area)
        shares.append(share)
        component_rows.append({
            "component_index": int(idx),
            "cell_count": int(comp.size),
            "area_fraction_world": float(comp_area / total),
            "share_of_land": share,
        })

    largest = components[0]
    largest_mask = np.zeros(n, dtype=bool)
    largest_mask[largest] = True
    largest_area = max(float(area[largest].sum()), 1.0e-12)
    largest_share = float(largest_area / land_area)
    second_share = shares[1] if len(shares) > 1 else 0.0
    third_share = shares[2] if len(shares) > 2 else 0.0
    major_count = int(sum(share >= major_share_threshold for share in shares))
    nonlargest_major_count = int(
        sum(share >= 0.08 for share in shares[1:])
    )

    domain_rows: list[dict[str, Any]] = []
    valid = largest[continent_id[largest] >= 0]
    if valid.size:
        for cid in sorted(int(x) for x in np.unique(continent_id[valid])):
            mask = largest[continent_id[largest] == cid]
            dom_area = float(area[mask].sum())
            domain_rows.append({
                "continent_id": int(cid),
                "cell_count": int(mask.size),
                "share_of_largest_land_component": float(dom_area / largest_area),
            })
    domain_rows.sort(key=lambda row: (
        -float(row["share_of_largest_land_component"]),
        int(row["continent_id"]),
    ))
    domain_shares = np.asarray([
        float(row["share_of_largest_land_component"]) for row in domain_rows
    ], dtype=np.float64)
    dominant_domain_share = float(domain_shares[0]) if domain_shares.size else 0.0
    significant_domain_count = int(np.count_nonzero(domain_shares >= 0.08))
    effective_domain_count = float(
        1.0 / max(float(np.sum(domain_shares ** 2)), 1.0e-12)
    ) if domain_shares.size else 0.0

    land_neighbor_count = np.zeros(n, dtype=np.int16)
    largest_neighbor_count = np.zeros(n, dtype=np.int16)
    for c in np.where(land)[0]:
        land_neighbor_count[int(c)] = int(
            sum(bool(land[int(nb)]) for nb in grid.neighbors[int(c)])
        )
    for c in largest:
        largest_neighbor_count[int(c)] = int(
            sum(bool(largest_mask[int(nb)]) for nb in grid.neighbors[int(c)])
        )

    lowland_or_shelf = rel <= 650.0
    weak_detail = np.isin(terrain_detail, [-1, 0, 1, 2, 3, 4])
    largest_bridge_candidate = (
        largest_mask
        & (largest_neighbor_count <= 3)
        & lowland_or_shelf
        & weak_detail
    )
    candidate_fraction = float(
        area[largest_bridge_candidate].sum() / largest_area
    )
    largest_robust_mask = largest_mask & ~largest_bridge_candidate
    robust_pieces = _components(grid, largest_robust_mask)
    robust_shares = sorted(
        (float(area[comp].sum()) / largest_area for comp in robust_pieces),
        reverse=True,
    )
    robust_piece_count = int(sum(share >= 0.06 for share in robust_shares))
    largest_piece_after = float(robust_shares[0]) if robust_shares else 0.0
    neck_split_gain = max(0, robust_piece_count - 1)

    all_bridge_candidate = (
        land
        & (land_neighbor_count <= 3)
        & lowland_or_shelf
        & weak_detail
    )
    all_robust_pieces = _components(grid, land & ~all_bridge_candidate)
    all_robust_major_count = int(sum(
        float(area[comp].sum()) / land_area >= major_share_threshold
        for comp in all_robust_pieces
    ))
    internal_seaway_openings = int(round(float(world.g(
        "terrain.last_p110b_internal_domain_seaway_openings", 0.0))))
    internal_seaway_area_fraction = float(world.g(
        "terrain.last_p110b_internal_domain_seaway_area_fraction", 0.0))
    internal_seaway_object_count = int(
        len(world.objects.get("terrain.p110b_internal_domain_seaways", []))
    )
    internal_boundary_count = int(round(float(world.g(
        "terrain.last_p110b_internal_domain_boundary_count", 0.0))))
    internal_boundary_area_fraction = float(world.g(
        "terrain.last_p110b_internal_domain_boundary_area_fraction", 0.0))
    internal_boundary_object_count = int(
        len(world.objects.get("terrain.p110b_internal_domain_boundaries", []))
    )

    monolithic_large_component = (
        largest_share >= 0.40
        and significant_domain_count <= 1
        and effective_domain_count < 1.6
    )
    under_separated_large_component = (
        largest_share >= 0.54
        and nonlargest_major_count < 2
    )
    neck_supported_overconnection = (
        largest_share >= 0.48
        and robust_piece_count >= 2
        and candidate_fraction <= 0.12
    )
    no_deep_split_large_component = (
        largest_share >= 0.52
        and robust_piece_count <= 1
        and significant_domain_count <= 2
    )
    terminal_supercontinent_like = bool(
        largest_share >= 0.64
        or monolithic_large_component
        or under_separated_large_component
        or no_deep_split_large_component
    )
    modern_multipolar_overconnected = bool(
        (largest_share >= 0.50 and significant_domain_count <= 2)
        or (largest_share >= 0.42 and significant_domain_count <= 1)
        or (largest_share >= 0.46 and robust_piece_count <= 1
            and nonlargest_major_count < 2)
    )

    score = 0.0
    score += 0.35 * float(np.clip((largest_share - 0.42) / 0.24, 0.0, 1.0))
    score += 0.20 * float(np.clip((dominant_domain_share - 0.58) / 0.34, 0.0, 1.0))
    score += 0.20 * float(np.clip((2 - nonlargest_major_count) / 2, 0.0, 1.0))
    score += 0.15 if robust_piece_count <= 1 and largest_share >= 0.46 else 0.0
    score += 0.10 if significant_domain_count <= 1 and largest_share >= 0.40 else 0.0
    if neck_supported_overconnection:
        score = max(score, 0.45)
    score = float(np.clip(score, 0.0, 1.0))

    return {
        "schema": "aevum.p110b_terminal_supercontinent_diagnostics.v1",
        "major_share_threshold": float(major_share_threshold),
        "terminal_supercontinent_score": score,
        "terminal_supercontinent_like": terminal_supercontinent_like,
        "modern_multipolar_overconnected": modern_multipolar_overconnected,
        "land_component_count": int(len(components)),
        "major_land_component_count": major_count,
        "nonlargest_major_land_component_count": nonlargest_major_count,
        "largest_land_share": largest_share,
        "second_land_share": float(second_share),
        "third_land_share": float(third_share),
        "largest_land_dominant_continent_share": dominant_domain_share,
        "largest_land_significant_continent_domain_count": significant_domain_count,
        "largest_land_effective_continent_domain_count": effective_domain_count,
        "largest_land_bridge_candidate_fraction": candidate_fraction,
        "largest_land_robust_piece_count_after_neck_removal": robust_piece_count,
        "largest_land_largest_piece_share_after_neck_removal": largest_piece_after,
        "largest_land_neck_split_gain": int(neck_split_gain),
        "all_land_robust_major_piece_count_after_neck_removal": (
            all_robust_major_count
        ),
        "internal_domain_seaway_opening_count": internal_seaway_openings,
        "internal_domain_seaway_area_fraction_world": internal_seaway_area_fraction,
        "internal_domain_seaway_object_count": internal_seaway_object_count,
        "internal_domain_boundary_count": internal_boundary_count,
        "internal_domain_boundary_area_fraction_world": internal_boundary_area_fraction,
        "internal_domain_boundary_object_count": internal_boundary_object_count,
        "monolithic_large_component": bool(monolithic_large_component),
        "under_separated_large_component": bool(under_separated_large_component),
        "neck_supported_overconnection": bool(neck_supported_overconnection),
        "no_deep_split_large_component": bool(no_deep_split_large_component),
        "component_rows_top5": component_rows[:5],
        "continent_domain_rows_top8": domain_rows[:8],
    }


def _p110b_lineage_survival_metrics(
    world: Any,
    land: np.ndarray,
    *,
    major_share_threshold: float,
) -> dict[str, Any]:
    """Measure whether major terminal landmasses have inherited lineage support."""
    grid = world.grid
    n = int(grid.n)
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    land = np.asarray(land, dtype=bool)
    land_area = max(float(area[land].sum()), 1.0e-12)

    def int_field(name: str, default: float) -> np.ndarray:
        arr = np.asarray(world.get_field(name, default), dtype=np.float64)
        if arr.shape != (n,):
            arr = np.full(n, default, dtype=np.float64)
        return arr.astype(int)

    def float_field(name: str, default: float) -> np.ndarray:
        arr = np.asarray(world.get_field(name, default), dtype=np.float64)
        if arr.shape != (n,):
            arr = np.full(n, default, dtype=np.float64)
        return np.nan_to_num(arr, nan=default, posinf=default, neginf=default)

    continent_id = int_field("tectonics.continent_id", -1.0)
    terrane_id = int_field("tectonics.terrane_id", -1.0)
    domain = int_field("crust.domain", -1.0)
    crust_type = float_field("crust.type", 0.0)
    stability = float_field("crust.stability", 0.0)

    passive_process = (
        _boundary_mask(world, "passive_margin")
        | _object_mask(world, "tectonics.boundary_objects", {"passive_margin"})
        | _object_mask(world, "tectonics.boundary_provinces", {"passive_margin"})
        | _object_mask(world, "tectonics.passive_margins", {"passive_margin"})
    )
    suture_process = (
        _boundary_mask(world, "suture")
        | _boundary_mask(world, "collision")
        | _object_mask(world, "tectonics.boundary_objects", {"suture", "collision"})
        | _object_mask(world, "tectonics.sutures", {"suture", "collision"})
    )
    craton_process = (
        _object_mask(world, "tectonics.cratons", {"craton"})
        | _object_mask(world, "tectonics.shields", {"shield", "craton"})
        | _object_mask(world, "tectonics.platforms", {"platform"})
    )
    terrane_process = (
        _object_mask(world, "tectonics.terranes", {"terrane", "accreted_terrane"})
        | _object_mask(world, "terrain.arc_plume_landforms", {"microcontinent"})
    )
    oceanic_landform_process = (
        _boundary_mask(world, "ridge")
        | _boundary_mask(world, "trench")
        | _boundary_mask(world, "subduction")
        | _object_mask(world, "tectonics.boundary_objects", {"ridge", "trench", "subduction"})
        | _object_mask(world, "tectonics.boundary_provinces", {
            "mid_ocean_ridge",
            "island_arc_trench",
            "ocean_ocean_subduction_trench",
        })
        | _object_mask(world, "terrain.ocean_fabric", {
            "spreading_center",
            "seamount_chain",
            "hotspot_track",
            "oceanic_plateau",
        })
        | _object_mask(world, "terrain.margin_landforms", {"volcanic_arc", "trench"})
        | _object_mask(world, "terrain.arc_plume_landforms", {
            "island_arc",
            "seamount_chain",
            "hotspot_track",
            "oceanic_plateau",
            "back_arc_basin",
        })
    )
    craton_support = land & (
        (domain == int(DOMAIN_CRATON))
        | (stability >= 0.78)
        | craton_process
    )
    boundary_memory = land & (
        (domain == int(DOMAIN_CONTINENTAL_MARGIN))
        | (domain == int(DOMAIN_SUTURE))
        | (domain == int(DOMAIN_ACCRETED_TERRANE))
        | passive_process
        | suture_process
        | terrane_process
    )
    lineage_mask = land & (continent_id >= 0)
    terrane_mask = land & (terrane_id >= 0)
    geologic_support = craton_support | boundary_memory | terrane_mask

    rows: list[dict[str, Any]] = []
    components = _components(grid, land)
    for comp_index, comp in enumerate(components):
        comp_area = float(area[comp].sum())
        share = comp_area / land_area
        if share < major_share_threshold:
            continue
        valid = comp[continent_id[comp] >= 0]
        dominant_id = -1
        dominant_share = 0.0
        if valid.size:
            ids = sorted(int(x) for x in np.unique(continent_id[valid]))
            best_area = -1.0
            for cid in ids:
                cid_area = float(area[comp[continent_id[comp] == cid]].sum())
                if cid_area > best_area:
                    best_area = cid_area
                    dominant_id = int(cid)
            dominant_share = float(best_area / max(comp_area, 1.0e-12))
        lineage_area_share = float(area[valid].sum() / max(comp_area, 1.0e-12))
        craton_fraction = float(area[comp[craton_support[comp]]].sum()
                                / max(comp_area, 1.0e-12))
        boundary_fraction = float(area[comp[boundary_memory[comp]]].sum()
                                  / max(comp_area, 1.0e-12))
        terrane_fraction = float(area[comp[terrane_mask[comp]]].sum()
                                 / max(comp_area, 1.0e-12))
        support_fraction = float(area[comp[geologic_support[comp]]].sum()
                                 / max(comp_area, 1.0e-12))
        continental_fraction = float(area[comp[crust_type[comp] >= 0.5]].sum()
                                     / max(comp_area, 1.0e-12))
        oceanic_fraction = float(area[comp[crust_type[comp] < 0.5]].sum()
                                 / max(comp_area, 1.0e-12))
        oceanic_process_fraction = float(
            area[comp[oceanic_landform_process[comp]]].sum()
            / max(comp_area, 1.0e-12)
        )
        dominant_supported = bool(dominant_id >= 0 and dominant_share >= 0.35)
        collage_supported = bool(lineage_area_share >= 0.75)
        continental_subject = bool(
            continental_fraction >= 0.25
            or lineage_area_share >= 0.25
            or craton_fraction >= 0.10
            or terrane_fraction >= 0.10
        )
        oceanic_landform_supported = bool(
            not continental_subject
            and oceanic_fraction >= 0.70
            and oceanic_process_fraction >= 0.15
        )
        lineage_supported = bool(
            continental_subject
            and
            (dominant_supported or collage_supported)
            and support_fraction >= 0.03
        )
        provenance_supported = bool(lineage_supported or oceanic_landform_supported)
        support_basis = (
            "dominant_continent"
            if dominant_supported
            else (
                "lineage_collage"
                if collage_supported
                else (
                    "oceanic_landform_provenance"
                    if oceanic_landform_supported
                    else "insufficient_lineage"
                )
            )
        )
        component_class = (
            "continental_lineage_subject"
            if continental_subject
            else (
                "oceanic_exposed_landform"
                if oceanic_landform_supported
                else "unsupported_exposed_landform"
            )
        )
        rows.append({
            "component_index": int(comp_index),
            "cell_count": int(comp.size),
            "area_fraction_world": float(comp_area / total),
            "share_of_land": float(share),
            "component_class": component_class,
            "continental_crust_fraction": float(continental_fraction),
            "oceanic_crust_fraction": float(oceanic_fraction),
            "oceanic_landform_process_fraction": float(oceanic_process_fraction),
            "dominant_continent_id": int(dominant_id),
            "dominant_continent_share": float(dominant_share),
            "continent_lineage_area_share": float(lineage_area_share),
            "craton_core_fraction": float(craton_fraction),
            "boundary_memory_fraction": float(boundary_fraction),
            "terrane_lineage_fraction": float(terrane_fraction),
            "geologic_support_fraction": float(support_fraction),
            "continental_lineage_subject": bool(continental_subject),
            "lineage_supported": lineage_supported,
            "oceanic_landform_supported": bool(oceanic_landform_supported),
            "provenance_supported": provenance_supported,
            "lineage_support_basis": support_basis,
        })

    rows.sort(key=lambda row: (
        -float(row["share_of_land"]),
        int(row["component_index"]),
    ))
    supported_rows = [row for row in rows if bool(row["lineage_supported"])]
    continental_rows = [
        row for row in rows if bool(row["continental_lineage_subject"])
    ]
    continental_supported_rows = [
        row for row in continental_rows if bool(row["lineage_supported"])
    ]
    oceanic_landform_rows = [
        row for row in rows if bool(row["oceanic_landform_supported"])
    ]
    provenance_supported_rows = [
        row for row in rows if bool(row["provenance_supported"])
    ]
    supported_ids = {
        int(row["dominant_continent_id"])
        for row in supported_rows
        if int(row["dominant_continent_id"]) >= 0
    }
    unsupported_rows = [
        int(row["component_index"]) for row in rows
        if not bool(row["provenance_supported"])
    ]
    unsupported_continental_rows = [
        int(row["component_index"]) for row in continental_rows
        if not bool(row["lineage_supported"])
    ]
    min_dominant = min(
        (float(row["dominant_continent_share"]) for row in rows),
        default=0.0,
    )
    min_support = min(
        (float(row["geologic_support_fraction"]) for row in rows),
        default=0.0,
    )
    return {
        "schema": "aevum.p110b_lineage_survival.v2",
        "major_share_threshold": float(major_share_threshold),
        "major_component_count": int(len(rows)),
        "lineage_supported_major_component_count": int(len(supported_rows)),
        "continental_major_component_count": int(len(continental_rows)),
        "continental_lineage_supported_major_component_count": int(
            len(continental_supported_rows)
        ),
        "oceanic_landform_major_component_count": int(len(oceanic_landform_rows)),
        "provenance_supported_major_component_count": int(
            len(provenance_supported_rows)
        ),
        "nonlargest_lineage_supported_component_count": int(
            sum(1 for row in rows[1:] if bool(row["lineage_supported"]))),
        "independent_primary_lineage_count": int(len(supported_ids)),
        "lineage_area_fraction_of_land": float(
            area[lineage_mask].sum() / land_area),
        "terrane_area_fraction_of_land": float(
            area[terrane_mask].sum() / land_area),
        "craton_support_fraction_of_land": float(
            area[craton_support].sum() / land_area),
        "boundary_memory_fraction_of_land": float(
            area[boundary_memory].sum() / land_area),
        "geologic_support_fraction_of_land": float(
            area[geologic_support].sum() / land_area),
        "min_major_dominant_continent_share": float(min_dominant),
        "min_major_continent_lineage_area_share": float(min(
            (float(row["continent_lineage_area_share"]) for row in rows),
            default=0.0,
        )),
        "min_major_geologic_support_fraction": float(min_support),
        "unsupported_major_component_indices": unsupported_rows,
        "unsupported_continental_major_component_indices": unsupported_continental_rows,
        "oceanic_landform_component_indices": [
            int(row["component_index"]) for row in oceanic_landform_rows
        ],
        "component_rows": rows[:10],
    }


def _mask_component_topology(
    grid: Any,
    mask: np.ndarray,
    *,
    major_share_threshold: float,
) -> dict[str, Any]:
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    mask = np.asarray(mask, dtype=bool)
    mask_area = float(area[mask].sum())
    comps = _components(grid, mask)
    rows: list[dict[str, Any]] = []
    for comp in comps:
        comp_area = float(area[comp].sum())
        rows.append({
            "cell_count": int(comp.size),
            "area_fraction_world": float(comp_area / total),
            "share_of_mask": float(comp_area / max(mask_area, 1.0e-12)),
        })
    rows.sort(key=lambda row: (
        -float(row["share_of_mask"]),
        -int(row["cell_count"]),
    ))
    top_shares = [float(row["share_of_mask"]) for row in rows[:10]]
    top_world = [float(row["area_fraction_world"]) for row in rows[:10]]
    major_count = int(
        sum(float(row["share_of_mask"]) >= major_share_threshold for row in rows)
    )
    return {
        "component_count": int(len(rows)),
        "major_component_count": major_count,
        "major_share_threshold": float(major_share_threshold),
        "largest_component_share_of_mask": float(top_shares[0] if top_shares else 0.0),
        "largest_component_area_fraction_world": float(top_world[0] if top_world else 0.0),
        "top_component_shares_of_mask": top_shares,
        "top_component_area_fractions_world": top_world,
        "component_rows_top10": rows[:10],
    }


def _p110b_historical_supercontinent_trajectory_metrics(
    archive: Any | None,
    *,
    share_threshold: float = P110B_HISTORICAL_SUPERCONTINENT_SHARE,
    min_land_fraction: float = 0.10,
) -> dict[str, Any]:
    """Measure whether an archive spends too long in one connected landmass.

    Terminal P110B planform metrics can pass after late coastline or lineage
    edits while the simulated history was still locked into a single continent
    for most of deep time.  This diagnostic intentionally looks only at saved
    archive frames, so its confidence is tied to the requested frame count.
    """
    if archive is None:
        return {
            "schema": "aevum.p110b_historical_supercontinent_trajectory.v1",
            "skipped": True,
            "skip_reason": "archive_not_available",
            "frame_count": 0,
            "usable_frame_count": 0,
            "supercontinent_share_threshold": float(share_threshold),
            "supercontinent_frame_count": 0,
            "supercontinent_frame_fraction": 0.0,
            "max_largest_land_component_share": 0.0,
            "median_largest_land_component_share": 0.0,
            "max_consecutive_supercontinent_duration_myr": 0.0,
            "supercontinent_time_window_myr": 0.0,
            "long_lived_supercontinent_like": False,
            "frame_rows": [],
        }

    world = getattr(archive, "world", None)
    grid = getattr(world, "grid", None)
    frames = list(getattr(archive, "frames", []) or [])
    if grid is None:
        return {
            "schema": "aevum.p110b_historical_supercontinent_trajectory.v1",
            "skipped": True,
            "skip_reason": "archive_grid_not_available",
            "frame_count": int(len(frames)),
            "usable_frame_count": 0,
            "supercontinent_share_threshold": float(share_threshold),
            "supercontinent_frame_count": 0,
            "supercontinent_frame_fraction": 0.0,
            "max_largest_land_component_share": 0.0,
            "median_largest_land_component_share": 0.0,
            "max_consecutive_supercontinent_duration_myr": 0.0,
            "supercontinent_time_window_myr": 0.0,
            "long_lived_supercontinent_like": False,
            "frame_rows": [],
        }

    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    rows: list[dict[str, Any]] = []
    for frame in frames:
        fields = getattr(frame, "fields", {}) or {}
        if "terrain.elevation_m" not in fields:
            continue
        elev = np.asarray(fields["terrain.elevation_m"], dtype=np.float64)
        if elev.shape != (int(grid.n),):
            continue
        globals_ = getattr(frame, "globals", {}) or {}
        sea_level = float(globals_.get("ocean.sea_level_m", 0.0))
        land = elev >= sea_level
        land_area = float(area[land].sum())
        topology = _mask_component_topology(
            grid,
            land,
            major_share_threshold=P110A_MAJOR_LAND_COMPONENT_SHARE,
        )
        shares = list(topology["top_component_shares_of_mask"])
        largest = float(topology["largest_component_share_of_mask"])
        supercontinent_like = bool(
            land_area / total >= min_land_fraction
            and largest >= float(share_threshold)
        )
        rows.append({
            "time_myr": float(getattr(frame, "time_myr", 0.0)),
            "sea_level_m": sea_level,
            "land_fraction": float(land_area / total),
            "land_component_count": int(topology["component_count"]),
            "major_land_component_count": int(topology["major_component_count"]),
            "largest_land_component_share": largest,
            "second_land_component_share": float(shares[1] if len(shares) > 1 else 0.0),
            "third_land_component_share": float(shares[2] if len(shares) > 2 else 0.0),
            "supercontinent_like": supercontinent_like,
        })

    rows.sort(key=lambda row: float(row["time_myr"]))
    if not rows:
        return {
            "schema": "aevum.p110b_historical_supercontinent_trajectory.v1",
            "skipped": True,
            "skip_reason": "no_usable_terrain_frames",
            "frame_count": int(len(frames)),
            "usable_frame_count": 0,
            "supercontinent_share_threshold": float(share_threshold),
            "supercontinent_frame_count": 0,
            "supercontinent_frame_fraction": 0.0,
            "max_largest_land_component_share": 0.0,
            "median_largest_land_component_share": 0.0,
            "max_consecutive_supercontinent_duration_myr": 0.0,
            "supercontinent_time_window_myr": 0.0,
            "long_lived_supercontinent_like": False,
            "frame_rows": [],
        }

    largest_values = np.asarray([
        float(row["largest_land_component_share"]) for row in rows
    ], dtype=np.float64)
    flagged_times = [
        float(row["time_myr"]) for row in rows if bool(row["supercontinent_like"])
    ]
    flagged_count = int(len(flagged_times))
    frame_fraction = float(flagged_count / max(len(rows), 1))
    window_myr = (
        float(max(flagged_times) - min(flagged_times))
        if len(flagged_times) >= 2 else 0.0
    )
    longest_consecutive = 0.0
    current = 0.0
    for prev, cur in zip(rows[:-1], rows[1:]):
        dt = max(float(cur["time_myr"]) - float(prev["time_myr"]), 0.0)
        if bool(prev["supercontinent_like"]) and bool(cur["supercontinent_like"]):
            current += dt
            longest_consecutive = max(longest_consecutive, current)
        else:
            current = 0.0

    wide_recurrence_window = bool(
        window_myr > P110B_HISTORICAL_SUPERCONTINENT_MAX_WINDOW_MYR
        and flagged_count >= 3
    )
    long_lived = bool(
        len(rows) >= 4
        and (
            frame_fraction > P110B_HISTORICAL_SUPERCONTINENT_MAX_FRAME_FRACTION
            or longest_consecutive > P110B_HISTORICAL_SUPERCONTINENT_MAX_DURATION_MYR
        )
    )
    final_row = rows[-1]
    return {
        "schema": "aevum.p110b_historical_supercontinent_trajectory.v1",
        "skipped": False,
        "skip_reason": "",
        "frame_count": int(len(frames)),
        "usable_frame_count": int(len(rows)),
        "supercontinent_share_threshold": float(share_threshold),
        "supercontinent_frame_count": flagged_count,
        "supercontinent_frame_fraction": frame_fraction,
        "max_largest_land_component_share": float(np.max(largest_values)),
        "median_largest_land_component_share": float(np.median(largest_values)),
        "final_largest_land_component_share": float(
            final_row["largest_land_component_share"]),
        "final_supercontinent_like": bool(final_row["supercontinent_like"]),
        "max_consecutive_supercontinent_duration_myr": float(longest_consecutive),
        "supercontinent_time_window_myr": float(window_myr),
        "wide_recurrent_supercontinent_window": wide_recurrence_window,
        "long_lived_supercontinent_like": long_lived,
        "late_recovery_after_long_supercontinent": bool(
            long_lived and not bool(final_row["supercontinent_like"])),
        "thresholds": {
            "max_frame_fraction": P110B_HISTORICAL_SUPERCONTINENT_MAX_FRAME_FRACTION,
            "max_consecutive_duration_myr": (
                P110B_HISTORICAL_SUPERCONTINENT_MAX_DURATION_MYR
            ),
            "max_time_window_myr": P110B_HISTORICAL_SUPERCONTINENT_MAX_WINDOW_MYR,
        },
        "frame_rows": rows,
    }


def _ocean_basin_planform_metrics(world: Any, ocean: np.ndarray) -> dict[str, Any]:
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    ocean = np.asarray(ocean, dtype=bool)
    ocean_area = max(float(area[ocean].sum()), 1.0e-12)
    basin_id_field = world.fields.get("ocean.basin_id")
    rows: list[dict[str, Any]] = []
    source = "ocean_basin_objects"

    object_rows = _ocean_basin_object_rows(world, ocean_area, total)
    if object_rows:
        rows = object_rows
    else:
        source = "connected_ocean_components"

    if not rows and basin_id_field is not None:
        basin_id = np.asarray(basin_id_field, dtype=np.int64)
        valid = ocean & (basin_id >= 0)
        if basin_id.shape == (grid.n,) and valid.any():
            source = "ocean.basin_id"
            for bid in sorted(int(x) for x in np.unique(basin_id[valid])):
                mask = valid & (basin_id == bid)
                barea = float(area[mask].sum())
                rows.append({
                    "basin_id": int(bid),
                    "cell_count": int(mask.sum()),
                    "area_fraction_world": float(barea / total),
                    "share_of_ocean": float(barea / ocean_area),
                    "source": source,
                })

    if not rows:
        for idx, comp in enumerate(_components(grid, ocean)):
            barea = float(area[comp].sum())
            rows.append({
                "basin_id": int(idx),
                "cell_count": int(comp.size),
                "area_fraction_world": float(barea / total),
                "share_of_ocean": float(barea / ocean_area),
                "source": source,
            })

    rows.sort(key=lambda row: (
        -float(row["share_of_ocean"]),
        int(row["basin_id"]),
    ))
    top_shares = [float(row["share_of_ocean"]) for row in rows[:10]]
    major_count = int(
        sum(float(row["share_of_ocean"]) >= P110A_MAJOR_OCEAN_BASIN_SHARE for row in rows)
    )
    return {
        "source": source,
        "basin_count": int(len(rows)),
        "major_basin_count": major_count,
        "largest_basin_share_of_ocean": float(top_shares[0] if top_shares else 0.0),
        "top_basin_shares_of_ocean": top_shares,
        "top_basin_area_fractions_world": [
            float(row["area_fraction_world"]) for row in rows[:10]
        ],
        "basin_rows_top10": rows[:10],
    }


def _ocean_gateway_topology_metrics(world: Any, ocean: np.ndarray) -> dict[str, Any]:
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    ocean = np.asarray(ocean, dtype=bool)
    ocean_area = max(float(area[ocean].sum()), 1.0e-12)

    def int_field(name: str, default: float = -1.0) -> np.ndarray:
        arr = np.asarray(world.get_field(name, default), dtype=np.float64)
        if arr.shape != (grid.n,):
            arr = np.full(grid.n, default, dtype=np.float64)
        return arr.astype(np.int64)

    basin_id = int_field("ocean.basin_id", -1.0)
    gateway_id = int_field("ocean.gateway_id", -1.0)
    gateway_system_id = int_field("ocean.gateway_system_id", -1.0)
    depth_province = int_field("ocean.depth_province", -1.0)
    gateway_cells = ocean & (gateway_id >= 0)
    gateway_system_cells = ocean & (gateway_system_id >= 0)
    restricted = ocean & (depth_province == 7)

    tect_gateway_cells = np.zeros(grid.n, dtype=bool)
    tect_gateway_phase_codes: set[int] = set()
    tect_status_counts: dict[str, int] = {}
    for obj in world.objects.get("tectonics.ocean_gateways", []):
        if str(obj.get("kind", obj.get("type", ""))) not in {"ocean_gateway", "gateway"}:
            continue
        status = str(obj.get("status", "unknown"))
        tect_status_counts[status] = tect_status_counts.get(status, 0) + 1
        phase = obj.get("phase_code")
        if phase is not None:
            tect_gateway_phase_codes.add(int(round(float(phase))))
        cells = np.asarray(obj.get("cells", []), dtype=np.int64)
        cells = cells[(0 <= cells) & (cells < grid.n)]
        if cells.size:
            tect_gateway_cells[cells] = True

    terminal_rows: list[dict[str, Any]] = []
    for obj in world.objects.get("ocean.gateways", []):
        if str(obj.get("kind", obj.get("type", ""))) not in {"ocean_gateway", "gateway"}:
            continue
        basin_ids = sorted(
            int(x) for x in obj.get("basin_ids", []) or []
            if x is not None and int(x) >= 0
        )
        phase_codes = sorted(
            int(round(float(x))) for x in obj.get("wilson_phase_codes", []) or []
            if x is not None
        )
        if obj.get("phase_code") is not None:
            phase_codes = sorted(dict.fromkeys([
                *phase_codes,
                int(round(float(obj.get("phase_code")))),
            ]))
        cells = np.asarray(obj.get("cells", []), dtype=np.int64)
        cells = cells[(0 <= cells) & (cells < grid.n)]
        area_fraction_of_ocean = obj.get("area_fraction_of_ocean")
        if area_fraction_of_ocean is None and cells.size:
            area_fraction_of_ocean = float(area[cells].sum() / ocean_area)
        elif area_fraction_of_ocean is None:
            area_fraction_of_ocean = 0.0
        terminal_rows.append({
            "id": obj.get("id", len(terminal_rows)),
            "cell_count": int(obj.get("cell_count", int(cells.size))),
            "area_fraction_of_ocean": float(area_fraction_of_ocean),
            "basin_ids": basin_ids,
            "basin_id_count": int(len(basin_ids)),
            "interbasin": bool(len(basin_ids) >= 2),
            "wilson_phase_codes": phase_codes,
            "phase_backed": bool(phase_codes),
        })

    terminal_rows.sort(key=lambda row: (
        -float(row["area_fraction_of_ocean"]),
        str(row["id"]),
    ))
    gateway_field_component_count = int(
        len(np.unique(gateway_id[gateway_cells])) if gateway_cells.any() else 0
    )
    terminal_gateway_count = int(len(terminal_rows))
    terminal_interbasin_count = int(
        sum(1 for row in terminal_rows if bool(row["interbasin"]))
    )
    terminal_phase_backed_count = int(
        sum(1 for row in terminal_rows if bool(row["phase_backed"]))
    )
    gateway_system_rows: list[dict[str, Any]] = []
    for obj in world.objects.get("ocean.gateway_systems", []):
        obj_type = str(obj.get("type", ""))
        if obj_type and obj_type != "ocean_gateway_system":
            continue
        basin_ids = sorted(
            int(x) for x in obj.get("basin_ids", []) or []
            if x is not None and int(x) >= 0
        )
        phase_codes = sorted(
            int(round(float(x))) for x in obj.get("wilson_phase_codes", []) or []
            if x is not None
        )
        cells = np.asarray(obj.get("cells", []), dtype=np.int64)
        cells = cells[(0 <= cells) & (cells < grid.n)]
        area_fraction_of_ocean = obj.get("area_fraction_of_ocean")
        if area_fraction_of_ocean is None and cells.size:
            area_fraction_of_ocean = float(area[cells].sum() / ocean_area)
        elif area_fraction_of_ocean is None:
            area_fraction_of_ocean = 0.0
        fragment_count = int(obj.get("fragment_count", 1))
        gateway_system_rows.append({
            "id": obj.get("id", len(gateway_system_rows)),
            "kind": str(obj.get("kind", "unknown")),
            "fragment_count": fragment_count,
            "cell_count": int(obj.get("cell_count", int(cells.size))),
            "area_fraction_of_ocean": float(area_fraction_of_ocean),
            "basin_ids": basin_ids,
            "basin_id_count": int(len(basin_ids)),
            "interbasin": bool(obj.get("interbasin", len(basin_ids) >= 2)),
            "wilson_phase_codes": phase_codes,
            "phase_backed": bool(obj.get("phase_backed", bool(phase_codes))),
            "restricted_fraction": float(obj.get("restricted_fraction", 0.0)),
        })

    if not gateway_system_rows and terminal_rows:
        gateway_system_rows = [
            {
                "id": row["id"],
                "kind": "legacy_terminal_gateway",
                "fragment_count": 1,
                "cell_count": int(row["cell_count"]),
                "area_fraction_of_ocean": float(row["area_fraction_of_ocean"]),
                "basin_ids": list(row["basin_ids"]),
                "basin_id_count": int(row["basin_id_count"]),
                "interbasin": bool(row["interbasin"]),
                "wilson_phase_codes": list(row["wilson_phase_codes"]),
                "phase_backed": bool(row["phase_backed"]),
                "restricted_fraction": 0.0,
            }
            for row in terminal_rows
        ]

    gateway_system_rows.sort(key=lambda row: (
        -float(row["area_fraction_of_ocean"]),
        str(row["id"]),
    ))
    gateway_system_field_component_count = int(
        len(np.unique(gateway_system_id[gateway_system_cells]))
        if gateway_system_cells.any()
        else 0
    )
    terminal_gateway_system_count = int(len(gateway_system_rows))
    terminal_interbasin_gateway_system_count = int(
        sum(1 for row in gateway_system_rows if bool(row["interbasin"]))
    )
    terminal_phase_backed_gateway_system_count = int(
        sum(1 for row in gateway_system_rows if bool(row["phase_backed"]))
    )
    gateway_system_kind_counts: dict[str, int] = {}
    for row in gateway_system_rows:
        kind = str(row["kind"])
        gateway_system_kind_counts[kind] = gateway_system_kind_counts.get(kind, 0) + 1
    gateway_fragment_to_system_ratio = float(
        terminal_gateway_count / max(terminal_gateway_system_count, 1)
    )

    ocean_topology = _mask_component_topology(
        grid,
        ocean,
        major_share_threshold=P110A_MAJOR_OCEAN_BASIN_SHARE,
    )
    major_disconnected = 0
    unbacked_major_disconnected = 0
    for row_index, row in enumerate(ocean_topology["component_rows_top10"]):
        share = float(row["share_of_mask"])
        if row_index == 0 or share < P110A_MAJOR_OCEAN_BASIN_SHARE:
            continue
        major_disconnected += 1

    if major_disconnected:
        comps = _components(grid, ocean)
        comps.sort(key=lambda comp: float(area[comp].sum()), reverse=True)
        for comp in comps[1:]:
            comp_area = float(area[comp].sum())
            if comp_area / ocean_area < P110A_MAJOR_OCEAN_BASIN_SHARE:
                continue
            comp_mask = np.zeros(grid.n, dtype=bool)
            comp_mask[comp] = True
            if not ((gateway_cells | tect_gateway_cells) & comp_mask).any():
                unbacked_major_disconnected += 1

    restricted_components = [
        comp for comp in _components(grid, restricted)
        if float(area[comp].sum()) / ocean_area >= 0.001
    ]
    restricted_basin_count = 0
    for obj in world.objects.get("ocean.basins", []):
        if str(obj.get("kind", obj.get("type", ""))) not in {
            "ocean_basin", "ocean_basin_lifecycle", "basin", ""
        }:
            continue
        if float(obj.get("restricted_fraction", 0.0)) >= 0.15:
            restricted_basin_count += 1

    return {
        "schema": "aevum.p110b_ocean_gateway_topology.v1",
        "terminal_gateway_count": terminal_gateway_count,
        "gateway_field_component_count": gateway_field_component_count,
        "terminal_interbasin_gateway_count": terminal_interbasin_count,
        "terminal_phase_backed_gateway_count": terminal_phase_backed_count,
        "terminal_gateway_system_count": terminal_gateway_system_count,
        "gateway_system_field_component_count": gateway_system_field_component_count,
        "terminal_interbasin_gateway_system_count": (
            terminal_interbasin_gateway_system_count
        ),
        "terminal_phase_backed_gateway_system_count": (
            terminal_phase_backed_gateway_system_count
        ),
        "gateway_fragment_to_system_ratio": gateway_fragment_to_system_ratio,
        "gateway_system_kind_counts": dict(sorted(gateway_system_kind_counts.items())),
        "tectonic_gateway_count": int(sum(tect_status_counts.values())),
        "tectonic_gateway_status_counts": dict(sorted(tect_status_counts.items())),
        "tectonic_gateway_phase_code_count": int(len(tect_gateway_phase_codes)),
        "restricted_ocean_fraction": float(area[restricted].sum() / ocean_area),
        "restricted_component_count": int(len(restricted_components)),
        "restricted_basin_count": int(restricted_basin_count),
        "major_ocean_component_count": int(ocean_topology["major_component_count"]),
        "disconnected_major_ocean_component_count": int(major_disconnected),
        "unbacked_major_disconnected_ocean_component_count": int(
            unbacked_major_disconnected
        ),
        "terminal_gateway_rows_top10": terminal_rows[:10],
        "terminal_gateway_system_rows_top10": gateway_system_rows[:10],
        "area_fraction_world_gateway_cells": float(area[gateway_cells].sum() / total),
        "area_fraction_world_gateway_system_cells": float(
            area[gateway_system_cells].sum() / total
        ),
        "area_fraction_world_tectonic_gateway_cells": float(
            area[tect_gateway_cells & ocean].sum() / total
        ),
    }


def _ocean_basin_object_rows(world: Any, ocean_area: float,
                             total_area: float) -> list[dict[str, Any]]:
    for object_set in ("ocean.basins", "tectonics.ocean_basins"):
        rows: list[dict[str, Any]] = []
        for obj in world.objects.get(object_set, []):
            if not isinstance(obj, dict):
                continue
            kind = str(obj.get("kind", obj.get("type", "")))
            if kind not in {"ocean_basin", "ocean_basin_lifecycle", "basin", ""}:
                continue
            share = obj.get("fraction_of_ocean")
            area_fraction = obj.get("area_fraction")
            if share is None and area_fraction is not None:
                share = float(area_fraction) * total_area / max(ocean_area, 1.0e-12)
            if area_fraction is None and share is not None:
                area_fraction = float(share) * ocean_area / max(total_area, 1.0e-12)
            if share is None:
                cells = np.asarray(obj.get("cells", []), dtype=np.int64)
                cells = cells[(0 <= cells) & (cells < world.grid.n)]
                if cells.size == 0:
                    continue
                cell_area = float(world.grid.cell_area[cells].sum())
                share = cell_area / max(ocean_area, 1.0e-12)
                area_fraction = cell_area / max(total_area, 1.0e-12)
            rows.append({
                "basin_id": _safe_int(obj.get("id", len(rows)), len(rows)),
                "cell_count": int(obj.get("cell_count", 0)),
                "area_fraction_world": float(area_fraction or 0.0),
                "share_of_ocean": float(share or 0.0),
                "source": object_set,
            })
        if rows:
            return rows
    return rows


def _seaway_effectiveness_metrics(world: Any) -> dict[str, Any]:
    attempts = list(world.objects.get("terrain.breakup_seaway_attempts", []))
    modern_endpoint_objects = list(
        world.objects.get("terrain.p1114_modern_endpoint_seaways", [])
    )
    reductions: list[float] = []
    applied_reductions: list[float] = []
    for attempt in attempts:
        before = float(attempt.get("largest_share_before", 0.0))
        after = float(attempt.get("best_new_largest_share", before))
        reduction = max(0.0, before - after)
        reductions.append(reduction)
        if bool(attempt.get("applied", False)):
            applied_reductions.append(reduction)
    return {
        "modern_underpartitioned_reference": bool(
            world.g("terrain.last_modern_seaway_underpartitioned_reference", 0.0) >= 1.0),
        "modern_applied": bool(
            world.g("terrain.last_modern_seaway_applied", 0.0) >= 1.0),
        "modern_best_reduction": float(
            world.g("terrain.last_modern_seaway_best_reduction", 0.0)),
        "modern_best_largest_share": float(
            world.g("terrain.last_modern_seaway_best_largest_share", 0.0)),
        "modern_candidate_count": int(
            world.g("terrain.last_modern_seaway_candidate_count", 0.0)),
        "modern_viable_candidate_count": int(
            world.g("terrain.last_modern_seaway_viable_candidate_count", 0.0)),
        "p110a_final_polish_passes": int(
            world.g("terrain.last_p110a_final_seaway_polish_passes", 0.0)),
        "p110a_final_polish_area_fraction_world": float(
            world.g("terrain.last_p110a_final_seaway_polish_area_fraction", 0.0)),
        "p110a_final_polish_largest_share_before": float(
            world.g("terrain.last_p110a_final_seaway_largest_share_before", 0.0)),
        "p110a_final_polish_largest_share_after": float(
            world.g("terrain.last_p110a_final_seaway_largest_share_after", 0.0)),
        "p154_final_planform_gate_attempted": bool(
            world.g("terrain.last_p154_final_planform_gate_attempted", 0.0) >= 1.0
        ),
        "p154_final_planform_gate_applied": bool(
            world.g("terrain.last_p154_final_planform_gate_applied", 0.0) >= 1.0
        ),
        "p154_final_planform_gate_reverted": bool(
            world.g("terrain.last_p154_final_planform_gate_reverted", 0.0) >= 1.0
        ),
        "p154_largest_share_before": float(
            world.g("terrain.last_p154_largest_share_before", 0.0)
        ),
        "p154_largest_share_after": float(
            world.g("terrain.last_p154_largest_share_after", 0.0)
        ),
        "p154_component_count_before": int(
            world.g("terrain.last_p154_component_count_before", 0.0)
        ),
        "p154_component_count_after": int(
            world.g("terrain.last_p154_component_count_after", 0.0)
        ),
        "p154_major_component_count_before": int(
            world.g("terrain.last_p154_major_component_count_before", 0.0)
        ),
        "p154_major_component_count_after": int(
            world.g("terrain.last_p154_major_component_count_after", 0.0)
        ),
        "p154_land_fraction_before": float(
            world.g("terrain.last_p154_land_fraction_before", 0.0)
        ),
        "p154_land_fraction_after": float(
            world.g("terrain.last_p154_land_fraction_after", 0.0)
        ),
        "p154_opened_area_fraction_world": float(
            world.g("terrain.last_p154_opened_area_fraction", 0.0)
        ),
        "p154_reject_code": int(
            world.g("terrain.last_p154_reject_code", 0.0)
        ),
        "breakup_opening_count": int(
            world.g("terrain.last_breakup_seaway_openings", 0.0)),
        "breakup_area_fraction_world": float(
            world.g("terrain.last_breakup_seaway_area_fraction", 0.0)),
        "breakup_attempt_count": int(len(attempts)),
        "breakup_attempt_max_reduction": float(max(reductions, default=0.0)),
        "breakup_applied_max_reduction": float(max(applied_reductions, default=0.0)),
        "p1114_modern_endpoint_seaway_count": int(
            world.g(
                "terrain.last_p1114_modern_endpoint_seaway_count",
                len(modern_endpoint_objects),
            )
        ),
        "p1114_modern_endpoint_seaway_area_fraction_world": float(
            world.g(
                "terrain.last_p1114_modern_endpoint_seaway_area_fraction",
                sum(
                    float(obj.get("area_fraction", 0.0))
                    for obj in modern_endpoint_objects
                ),
            )
        ),
        "p1114_modern_endpoint_seaway_domain_backed_count": int(
            world.g(
                "terrain.last_p1114_modern_endpoint_seaway_domain_backed_count",
                sum(
                    1 for obj in modern_endpoint_objects
                    if str(obj.get("basis", "")) == "continent_domain_boundary"
                ),
            )
        ),
        "p1114_modern_endpoint_seaway_weak_bridge_count": int(
            world.g(
                "terrain.last_p1114_modern_endpoint_seaway_weak_bridge_count",
                sum(
                    1 for obj in modern_endpoint_objects
                    if str(obj.get("basis", "")) == "weak_lowland_bridge"
                ),
            )
        ),
        "p1114_strict_endpoint_polish_attempted": bool(
            world.g("terrain.last_p1114_strict_endpoint_polish_attempted", 0.0)
            >= 1.0
        ),
        "p1114_strict_endpoint_polish_applied": bool(
            world.g("terrain.last_p1114_strict_endpoint_polish_applied", 0.0)
            >= 1.0
        ),
        "p1114_strict_endpoint_polish_area_fraction_world": float(
            world.g("terrain.last_p1114_strict_endpoint_polish_area_fraction", 0.0)
        ),
        "p1114_strict_endpoint_polish_preferred_ceiling": float(
            world.g("terrain.last_p1114_strict_endpoint_polish_preferred_ceiling", 0.0)
        ),
        "p1114_strict_endpoint_polish_min_largest_after": float(
            world.g("terrain.last_p1114_strict_endpoint_polish_min_largest_after", 0.0)
        ),
        "p1114_strict_endpoint_polish_largest_share_before": float(
            world.g(
                "terrain.last_p1114_strict_endpoint_polish_largest_share_before",
                0.0,
            )
        ),
        "p1114_strict_endpoint_polish_largest_share_after": float(
            world.g(
                "terrain.last_p1114_strict_endpoint_polish_largest_share_after",
                0.0,
            )
        ),
        "p1114_strict_endpoint_no_split_candidate_count": int(
            world.g("terrain.last_p1114_strict_endpoint_no_split_candidate_count", 0.0)
        ),
        "p1114_strict_endpoint_secondary_reject_count": int(
            world.g("terrain.last_p1114_strict_endpoint_secondary_reject_count", 0.0)
        ),
        "p1114_strict_endpoint_reduction_reject_count": int(
            world.g("terrain.last_p1114_strict_endpoint_reduction_reject_count", 0.0)
        ),
        "p1114_strict_endpoint_quality_reject_count": int(
            world.g("terrain.last_p1114_strict_endpoint_quality_reject_count", 0.0)
        ),
        "p111617_late_closed_ocean_attempted": bool(
            world.g("terrain.last_p111617_late_closed_ocean_attempted", 0.0)
            >= 1.0
        ),
        "p111617_late_closed_ocean_applied": bool(
            world.g("terrain.last_p111617_late_closed_ocean_applied", 0.0)
            >= 1.0
        ),
        "p111617_late_closed_ocean_component_count_before": int(
            world.g(
                "terrain.last_p111617_late_closed_ocean_component_count_before",
                0.0,
            )
        ),
        "p111617_late_closed_ocean_component_count_after": int(
            world.g(
                "terrain.last_p111617_late_closed_ocean_component_count_after",
                0.0,
            )
        ),
        "p111617_late_closed_ocean_score_before": float(
            world.g("terrain.last_p111617_late_closed_ocean_score_before", 0.0)
        ),
        "p111617_late_closed_ocean_score_after": float(
            world.g("terrain.last_p111617_late_closed_ocean_score_after", 0.0)
        ),
        "p111617_late_closed_ocean_max_candidate_fraction_before": float(
            world.g(
                "terrain.last_p111617_late_closed_ocean_max_candidate_fraction_before",
                0.0,
            )
        ),
        "p111617_late_closed_ocean_openings": int(
            world.g("terrain.last_p111617_late_closed_ocean_openings", 0.0)
        ),
        "p111617_late_closed_ocean_area_fraction_world": float(
            world.g("terrain.last_p111617_late_closed_ocean_area_fraction", 0.0)
        ),
        "p111619_late_domain_partition_attempted": bool(
            world.g("terrain.last_p111619_late_domain_partition_attempted", 0.0)
            >= 1.0
        ),
        "p111619_late_domain_partition_applied": bool(
            world.g("terrain.last_p111619_late_domain_partition_applied", 0.0)
            >= 1.0
        ),
        "p111619_late_domain_partition_candidate_count": int(
            world.g("terrain.last_p111619_late_domain_partition_candidate_count", 0.0)
        ),
        "p111619_late_domain_partition_area_fraction_world": float(
            world.g("terrain.last_p111619_late_domain_partition_area_fraction", 0.0)
        ),
        "p111619_late_domain_partition_largest_share_before": float(
            world.g("terrain.last_p111619_late_domain_partition_largest_before", 0.0)
        ),
        "p111619_late_domain_partition_largest_share_after": float(
            world.g("terrain.last_p111619_late_domain_partition_largest_after", 0.0)
        ),
        "p111619_late_domain_partition_major_count_before": int(
            world.g("terrain.last_p111619_late_domain_partition_major_count_before", 0.0)
        ),
        "p111619_late_domain_partition_major_count_after": int(
            world.g("terrain.last_p111619_late_domain_partition_major_count_after", 0.0)
        ),
        "p111619_late_domain_partition_domain_count": int(
            world.g("terrain.last_p111619_late_domain_partition_domain_count", 0.0)
        ),
        "p111619_late_domain_partition_domain_floor": float(
            world.g("terrain.last_p111619_late_domain_partition_domain_floor", 0.0)
        ),
        "p111619_late_domain_partition_land_fraction_before": float(
            world.g("terrain.last_p111619_late_domain_partition_land_fraction_before", 0.0)
        ),
        "p111619_late_domain_partition_land_fraction_after": float(
            world.g("terrain.last_p111619_late_domain_partition_land_fraction_after", 0.0)
        ),
        "p1114_modern_endpoint_seaway_rows_top10": [
            {
                "object_id": obj.get("object_id", obj.get("id", "")),
                "basis": obj.get("basis", ""),
                "cell_count": int(obj.get("cell_count", 0)),
                "area_fraction": float(obj.get("area_fraction", 0.0)),
                "largest_share_before": float(
                    obj.get("largest_share_before", 0.0)),
                "best_largest_share_after": float(
                    obj.get("best_largest_share_after", 0.0)),
                "domain_boundary_fraction": float(
                    obj.get("domain_boundary_fraction", 0.0)),
                "margin_suture_fraction": float(
                    obj.get("margin_suture_fraction", 0.0)),
                "weak_lowland_fraction": float(
                    obj.get("weak_lowland_fraction", 0.0)),
            }
            for obj in modern_endpoint_objects[:10]
        ],
    }


def _p110b_final_state_archetype_metrics(world: Any) -> dict[str, Any]:
    code = int(round(float(world.g(
        "terrain.last_p110b_final_state_archetype_code", 0.0))))
    return {
        "code": int(code),
        "name": P110B_ARCHETYPE_NAMES.get(code, "unknown"),
        "largest_share_preferred_ceiling": float(world.g(
            "terrain.last_p110b_largest_share_preferred_ceiling", 0.0)),
        "largest_share_soft_ceiling": float(world.g(
            "terrain.last_p110b_largest_share_soft_ceiling", 0.0)),
        "min_nonlargest_large_components": int(round(float(world.g(
            "terrain.last_p110b_min_nonlargest_large_components", 0.0)))),
    }


def _land_payback_bias_metrics(world: Any) -> dict[str, Any]:
    stages = list(world.objects.get("terrain.land_component_stage_telemetry", []))
    before = _stage_by_name(stages, "after_drop_unsupported")
    after = _stage_by_name(stages, "after_final_drop_payback")
    available = bool(
        before
        and after
        and not bool(before.get("approximate", False))
        and not bool(after.get("approximate", False))
    )
    before_share = float(before.get("largest_land_component_fraction", 0.0)) if before else 0.0
    after_share = float(after.get("largest_land_component_fraction", 0.0)) if after else 0.0
    return {
        "available": available,
        "before_stage": "after_drop_unsupported" if before else "",
        "after_stage": "after_final_drop_payback" if after else "",
        "largest_component_share_before": before_share,
        "largest_component_share_after": after_share,
        "largest_component_share_delta": float(after_share - before_share)
        if available else 0.0,
        "modern_seaway_payback_area_fraction": float(
            world.g("terrain.last_modern_seaway_payback_area_fraction", 0.0)),
        "modern_seaway_payback_candidate_area_fraction": float(
            world.g("terrain.last_modern_seaway_payback_candidate_area_fraction", 0.0)),
        "modern_coastline_payback_candidate_area_fraction": float(
            world.g("terrain.last_modern_coastline_payback_candidate_area_fraction", 0.0)),
        "p110a_component_payback_added_fraction": float(
            world.g("terrain.last_p110a_component_payback_added_fraction", 0.0)),
        "p110a_component_payback_candidate_fraction": float(
            world.g("terrain.last_p110a_component_payback_candidate_fraction", 0.0)),
        "p110a_component_payback_largest_share_before": float(
            world.g("terrain.last_p110a_component_payback_largest_share_before", 0.0)),
        "p110a_component_payback_largest_share_after": float(
            world.g("terrain.last_p110a_component_payback_largest_share_after", 0.0)),
        "p110a_component_payback_major_count_before": int(
            world.g("terrain.last_p110a_component_payback_major_count_before", 0.0)),
        "p110a_component_payback_major_count_after": int(
            world.g("terrain.last_p110a_component_payback_major_count_after", 0.0)),
        "p110a_component_payback_nonlargest_large_before": int(
            world.g("terrain.last_p110a_component_payback_nonlargest_large_before", 0.0)),
        "p110a_component_payback_nonlargest_large_after": int(
            world.g("terrain.last_p110a_component_payback_nonlargest_large_after", 0.0)),
        "p110a_component_payback_reject_code": int(
            world.g("terrain.last_p110a_component_payback_reject_code", 0.0)),
        "p110a_component_payback_trial_added_fraction": float(
            world.g("terrain.last_p110a_component_payback_trial_added_fraction", 0.0)),
        "p110a_component_payback_component_count_before": int(
            world.g("terrain.last_p110a_component_payback_component_count_before", 0.0)),
        "p110a_component_payback_component_count_after": int(
            world.g("terrain.last_p110a_component_payback_component_count_after", 0.0)),
        "p110a_component_payback_ribbon_before": float(
            world.g("terrain.last_p110a_component_payback_ribbon_before", 0.0)),
        "p110a_component_payback_ribbon_after": float(
            world.g("terrain.last_p110a_component_payback_ribbon_after", 0.0)),
        "p110a_component_payback_second_share_after": float(
            world.g("terrain.last_p110a_component_payback_second_share_after", 0.0)),
        "p110a_component_payback_third_share_after": float(
            world.g("terrain.last_p110a_component_payback_third_share_after", 0.0)),
        "p110a_component_payback_fourth_share_after": float(
            world.g("terrain.last_p110a_component_payback_fourth_share_after", 0.0)),
    }


def _stage_by_name(stages: list[Any], name: str) -> dict[str, Any]:
    for stage in stages:
        if isinstance(stage, dict) and str(stage.get("stage", "")) == name:
            return stage
    return {}


def _safe_int(value: Any, default: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(default)


def _boundary_mask(world: Any, kind: str) -> np.ndarray:
    grid = world.grid
    mask = np.zeros(grid.n, dtype=bool)
    cells = np.asarray(
        world.networks.get("tectonics.boundaries", {}).get(kind, []),
        dtype=np.int64,
    )
    cells = cells[(0 <= cells) & (cells < grid.n)]
    if cells.size:
        mask[cells] = True
    return mask


def _component_sizes(grid: Any, mask: np.ndarray) -> list[float]:
    return [float(grid.cell_area[comp].sum()) for comp in _components(grid, mask)]


def _components(grid: Any, mask: np.ndarray) -> list[np.ndarray]:
    mask = np.asarray(mask, dtype=bool)
    nodes = np.where(mask)[0]
    seen = np.zeros(grid.n, dtype=bool)
    out: list[np.ndarray] = []
    for start in nodes:
        if seen[int(start)]:
            continue
        stack = [int(start)]
        seen[int(start)] = True
        comp: list[int] = []
        while stack:
            c = stack.pop()
            comp.append(c)
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if mask[nb] and not seen[nb]:
                    seen[nb] = True
                    stack.append(nb)
        out.append(np.asarray(comp, dtype=np.int64))
    return out


def _dilate(grid: Any, mask: np.ndarray, passes: int = 1) -> np.ndarray:
    out = np.asarray(mask, dtype=bool).copy()
    for _ in range(max(0, int(passes))):
        add = np.zeros(grid.n, dtype=bool)
        for c in np.where(out)[0]:
            add[grid.neighbors[int(c)]] = True
        out |= add
    return out


def _same_kind_neighbor_count(grid: Any, mask: np.ndarray) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    count = np.zeros(grid.n, dtype=np.int16)
    for c in np.where(mask)[0]:
        count[int(c)] = int(np.count_nonzero(mask[grid.neighbors[int(c)]]))
    return count


def _distance_steps_to_mask(grid: Any, source: np.ndarray,
                            max_steps: int = 16) -> np.ndarray:
    source = np.asarray(source, dtype=bool)
    distance = np.full(grid.n, np.inf, dtype=np.float64)
    if not source.any():
        return distance
    queue = [int(c) for c in np.where(source)[0]]
    distance[queue] = 0.0
    head = 0
    while head < len(queue):
        c = queue[head]
        head += 1
        if distance[c] >= float(max_steps):
            continue
        for nb in grid.neighbors[c]:
            nb = int(nb)
            if np.isfinite(distance[nb]):
                continue
            distance[nb] = distance[c] + 1.0
            queue.append(nb)
    return distance


def _finite_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    valid = np.isfinite(x) & np.isfinite(y)
    if int(np.count_nonzero(valid)) < 3:
        return 0.0
    xv = x[valid]
    yv = y[valid]
    if float(np.std(xv)) <= 1.0e-12 or float(np.std(yv)) <= 1.0e-12:
        return 0.0
    return float(np.corrcoef(xv, yv)[0, 1])


def _hypsometry_stats(values: np.ndarray, area: np.ndarray,
                      mask: np.ndarray) -> dict[str, float]:
    if not np.any(mask):
        return {
            "mean_m": 0.0,
            "p02_m": 0.0,
            "p10_m": 0.0,
            "p50_m": 0.0,
            "p90_m": 0.0,
            "p98_m": 0.0,
            "max_m": 0.0,
        }
    vals = np.asarray(values[mask], dtype=np.float64)
    weights = np.asarray(area[mask], dtype=np.float64)
    return {
        "mean_m": float(np.average(vals, weights=weights)),
        "p02_m": float(_weighted_percentile(vals, weights, 2.0)),
        "p10_m": float(_weighted_percentile(vals, weights, 10.0)),
        "p50_m": float(_weighted_percentile(vals, weights, 50.0)),
        "p90_m": float(_weighted_percentile(vals, weights, 90.0)),
        "p98_m": float(_weighted_percentile(vals, weights, 98.0)),
        "max_m": float(np.max(vals)),
    }


def _weighted_percentile(values: np.ndarray, weights: np.ndarray,
                         percentile: float) -> float:
    vals = np.asarray(values, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    finite = np.isfinite(vals) & np.isfinite(w) & (w > 0.0)
    if not np.any(finite):
        return 0.0
    vals = vals[finite]
    w = w[finite]
    order = np.argsort(vals)
    vals = vals[order]
    w = w[order]
    cdf = np.cumsum(w)
    cutoff = float(np.clip(percentile, 0.0, 100.0)) / 100.0 * float(cdf[-1])
    idx = int(np.searchsorted(cdf, cutoff, side="left"))
    return float(vals[np.clip(idx, 0, vals.size - 1)])


def _band_key(lo: float, hi: float) -> str:
    if np.isinf(hi):
        return f"gt{int(lo)}m"
    return f"{int(lo)}_{int(hi)}m"


def _hypsometry_band_fractions(values: np.ndarray, area: np.ndarray,
                               mask: np.ndarray) -> dict[str, float]:
    out: dict[str, float] = {}
    denom = max(float(area[mask].sum()), 1.0e-12)
    vals = np.asarray(values, dtype=np.float64)
    for lo, hi in LAND_ELEVATION_BANDS_M:
        band = mask & (vals >= lo)
        if np.isfinite(hi):
            band &= vals < hi
        out[_band_key(lo, hi)] = float(area[band].sum() / denom)
    return out


def _same_grid_earth_hypsometry(grid: Any) -> dict[str, Any]:
    sampled = _sample_etopo5_to_grid(grid)
    if sampled is None:
        return {
            "available": False,
            "source": "data/reference/etopo5/ETOPO5.DAT",
        }
    area = np.asarray(grid.cell_area, dtype=np.float64)
    elev = np.asarray(sampled, dtype=np.float64)
    land = elev >= 0.0
    return {
        "available": True,
        "source": "data/reference/etopo5/ETOPO5.DAT",
        "land_fraction": float(area[land].sum() / max(float(area.sum()), 1.0e-12)),
        "land_hypsometry": _hypsometry_stats(elev, area, land),
        "land_elevation_band_fraction": _hypsometry_band_fractions(elev, area, land),
    }


def _same_grid_earth_planform(grid: Any) -> dict[str, Any]:
    sampled = _sample_etopo5_to_grid(grid)
    if sampled is None:
        return {
            "available": False,
            "source": "data/reference/etopo5/ETOPO5.DAT",
        }
    elev = np.asarray(sampled, dtype=np.float64)
    land = elev >= 0.0
    ocean = ~land
    land_topology = _mask_component_topology(
        grid,
        land,
        major_share_threshold=P110A_MAJOR_LAND_COMPONENT_SHARE,
    )
    ocean_connectivity = _mask_component_topology(
        grid,
        ocean,
        major_share_threshold=0.01,
    )
    disconnected_ocean_share = max(
        0.0,
        1.0 - float(ocean_connectivity["largest_component_share_of_mask"]),
    )
    return {
        "available": True,
        "source": "data/reference/etopo5/ETOPO5.DAT",
        "land": land_topology,
        "ocean_connectivity": {
            **ocean_connectivity,
            "closed_ocean_ring_score": float(disconnected_ocean_share),
            "disconnected_ocean_component_share": float(disconnected_ocean_share),
        },
        "summary": {
            "major_land_component_count": int(land_topology["major_component_count"]),
            "largest_land_component_share": float(
                land_topology["largest_component_share_of_mask"]),
            "closed_ocean_ring_score": float(disconnected_ocean_share),
        },
    }


def _p109_hypsometry_comparison(
    generated_stats: dict[str, float],
    generated_bands: dict[str, float],
    earth_hypsometry: dict[str, Any],
) -> dict[str, Any]:
    out_of_envelope: list[str] = []
    deltas: dict[str, float] = {}
    earth_stats = earth_hypsometry.get("land_hypsometry", {})
    earth_bands = earth_hypsometry.get("land_elevation_band_fraction", {})
    for key, value in generated_stats.items():
        if isinstance(value, (int, float)):
            deltas[key] = float(value) - float(earth_stats.get(key, 0.0))
    for key, value in generated_bands.items():
        deltas[f"band_{key}"] = float(value) - float(earth_bands.get(key, 0.0))
    for key, (lo, hi) in P109_LAND_HYPSOMETRY_ENVELOPE.items():
        value = float(generated_stats.get(key, 0.0))
        if value < lo or value > hi:
            out_of_envelope.append(key)
    for key, (lo, hi) in P109_LAND_BAND_ENVELOPE.items():
        value = float(generated_bands.get(key, 0.0))
        if value < lo or value > hi:
            out_of_envelope.append(f"band_{key}")
    return {
        "schema": "aevum.p109_hypsometry_comparison.v1",
        "earth_reference_available": bool(earth_hypsometry.get("available", False)),
        "out_of_envelope": out_of_envelope,
        "out_of_envelope_count": int(len(out_of_envelope)),
        "within_p109_envelope": bool(len(out_of_envelope) == 0),
        "generated_minus_earth": deltas,
    }


def _category_hypsometry_metrics(codes: np.ndarray, values: np.ndarray,
                                 area: np.ndarray,
                                 land: np.ndarray) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    total = max(float(area.sum()), 1.0e-12)
    land_area = max(float(area[land].sum()), 1.0e-12)
    codes = np.asarray(codes, dtype=np.int64)
    vals = np.asarray(values, dtype=np.float64)
    for code in sorted(int(x) for x in np.unique(codes[land])):
        mask = land & (codes == code)
        if not mask.any():
            continue
        stats = _hypsometry_stats(vals, area, mask)
        code_area = float(area[mask].sum())
        out[str(code)] = {
            "area_fraction_world": float(code_area / total),
            "area_fraction_land": float(code_area / land_area),
            **stats,
            "fraction_ge_2000m_of_land": float(
                area[mask & (vals >= 2000.0)].sum() / land_area),
            "fraction_ge_3000m_of_land": float(
                area[mask & (vals >= 3000.0)].sum() / land_area),
            "fraction_ge_4500m_of_land": float(
                area[mask & (vals >= 4500.0)].sum() / land_area),
        }
    return out


def _p109_sea_level_consistency_metrics(
    world: Any,
    sea_level_m: float,
    land_fraction: float,
) -> dict[str, Any]:
    target_land = float(getattr(world.spec, "target_land_fraction", 0.0))
    quantile_level = float(world.g("terrain.last_quantile_sea_level_m", sea_level_m))
    water_budget_level = float(
        world.g("terrain.last_water_budget_sea_level_m", sea_level_m))
    water_delta = float(
        world.g(
            "terrain.last_sea_level_water_budget_delta_m",
            quantile_level - water_budget_level,
        )
    )
    p109_mask_preserved = float(
        world.g("terrain.last_p109_land_mask_preserved", 0.0))
    return {
        "schema": "aevum.p109_sea_level_consistency.v1",
        "sea_level_m": float(sea_level_m),
        "target_land_fraction": float(target_land),
        "land_fraction": float(land_fraction),
        "land_fraction_minus_target": float(land_fraction - target_land),
        "quantile_sea_level_m": quantile_level,
        "water_budget_sea_level_m": water_budget_level,
        "quantile_minus_water_budget_m": water_delta,
        "p109_land_mask_preserved": bool(p109_mask_preserved >= 1.0),
        "p109_land_mean_before_m": float(
            world.g("terrain.last_p109_land_mean_before_m", 0.0)),
        "p109_land_mean_after_m": float(
            world.g("terrain.last_p109_land_mean_after_m", 0.0)),
        "p109_land_p50_before_m": float(
            world.g("terrain.last_p109_land_p50_before_m", 0.0)),
        "p109_land_p50_after_m": float(
            world.g("terrain.last_p109_land_p50_after_m", 0.0)),
        "p109_land_p90_before_m": float(
            world.g("terrain.last_p109_land_p90_before_m", 0.0)),
        "p109_land_p90_after_m": float(
            world.g("terrain.last_p109_land_p90_after_m", 0.0)),
        "p109_post_p110a_tail_area_fraction": float(
            world.g("terrain.last_p109_post_p110a_tail_area_fraction", 0.0)),
        "p109_post_p110a_shoulder_area_fraction": float(
            world.g("terrain.last_p109_post_p110a_shoulder_area_fraction", 0.0)),
        "p109_post_p110a_midland_area_fraction": float(
            world.g("terrain.last_p109_post_p110a_midland_area_fraction", 0.0)),
        "p109_post_p110a_lowland_area_fraction": float(
            world.g("terrain.last_p109_post_p110a_lowland_area_fraction", 0.0)),
        "p109_post_p110a_overcap_area_fraction": float(
            world.g("terrain.last_p109_post_p110a_overcap_area_fraction", 0.0)),
        "p109_post_p110a_p90_before_overcap_m": float(
            world.g("terrain.last_p109_post_p110a_p90_before_overcap_m", 0.0)),
        "p109_post_p110a_p90_after_overcap_m": float(
            world.g("terrain.last_p109_post_p110a_p90_after_overcap_m", 0.0)),
        "p109_post_p110a_p90_before_compensation_m": float(
            world.g("terrain.last_p109_post_p110a_p90_before_compensation_m", 0.0)),
        "p109_post_p110a_p90_after_compensation_m": float(
            world.g("terrain.last_p109_post_p110a_p90_after_compensation_m", 0.0)),
        "p109_final_p90_floor_area_fraction": float(
            world.g("terrain.last_p109_final_p90_floor_area_fraction", 0.0)),
        "p109_final_p90_floor_before_m": float(
            world.g("terrain.last_p109_final_p90_floor_before_m", 0.0)),
        "p109_final_p90_floor_after_m": float(
            world.g("terrain.last_p109_final_p90_floor_after_m", 0.0)),
        "p109_final_p90_ceiling_area_fraction": float(
            world.g("terrain.last_p109_final_p90_ceiling_area_fraction", 0.0)),
        "p109_final_p90_ceiling_before_m": float(
            world.g("terrain.last_p109_final_p90_ceiling_before_m", 0.0)),
        "p109_final_p90_ceiling_after_m": float(
            world.g("terrain.last_p109_final_p90_ceiling_after_m", 0.0)),
        "p109_final_lowland_shoulder_area_fraction": float(
            world.g("terrain.last_p109_final_lowland_shoulder_area_fraction", 0.0)),
        "p109_final_lowland_shoulder_before_fraction": float(
            world.g("terrain.last_p109_final_lowland_shoulder_before_fraction", 0.0)),
        "p109_final_lowland_shoulder_after_fraction": float(
            world.g("terrain.last_p109_final_lowland_shoulder_after_fraction", 0.0)),
        "p1112_adjusted_area_fraction": float(
            world.g("terrain.last_p1112_adjusted_area_fraction", 0.0)),
        "p1112_low500_fraction_before": float(
            world.g("terrain.last_p1112_low500_fraction_before", 0.0)),
        "p1112_low500_fraction_after": float(
            world.g("terrain.last_p1112_low500_fraction_after", 0.0)),
        "p1112_mean_before_m": float(
            world.g("terrain.last_p1112_mean_before_m", 0.0)),
        "p1112_mean_after_m": float(
            world.g("terrain.last_p1112_mean_after_m", 0.0)),
        "p1112_p50_before_m": float(
            world.g("terrain.last_p1112_p50_before_m", 0.0)),
        "p1112_p50_after_m": float(
            world.g("terrain.last_p1112_p50_after_m", 0.0)),
        "p1112_p90_before_m": float(
            world.g("terrain.last_p1112_p90_before_m", 0.0)),
        "p1112_p90_after_m": float(
            world.g("terrain.last_p1112_p90_after_m", 0.0)),
        "p1112_highland_preservation_delta_p50_m": float(
            world.g("terrain.last_p1112_highland_preservation_delta_p50_m", 0.0)),
        "p1112_land_mask_preserved": bool(
            world.g("terrain.last_p1112_land_mask_preserved", 0.0) >= 1.0),
        "land_fraction_consistent_with_target": bool(
            abs(float(land_fraction) - target_land) <= 0.05),
    }


def _code_area_fraction(codes: np.ndarray, area: np.ndarray,
                        total: float) -> dict[str, float]:
    out: dict[str, float] = {}
    for code in sorted(int(x) for x in np.unique(codes)):
        out[str(code)] = float(area[codes == code].sum() / max(total, 1.0e-12))
    return out


def _compact_entry_metrics(metrics: dict[str, Any]) -> dict[str, Any]:
    archive = metrics.get("array_archive", {})
    return {
        "cells": int(metrics["cells"]),
        "land_fraction": float(metrics["land_fraction"]),
        "land_hypsometry": metrics.get("land_hypsometry", {}),
        "land_elevation_band_fraction": metrics.get(
            "land_elevation_band_fraction", {}),
        "p109_hypsometry_comparison": {
            key: metrics.get("p109_hypsometry_comparison", {}).get(key)
            for key in (
                "earth_reference_available",
                "out_of_envelope_count",
                "within_p109_envelope",
                "out_of_envelope",
            )
        },
        "p109_sea_level_consistency": {
            key: metrics.get("p109_sea_level_consistency", {}).get(key)
            for key in (
                "land_fraction_minus_target",
                "quantile_minus_water_budget_m",
                "p109_land_mask_preserved",
                "p109_final_p90_ceiling_area_fraction",
                "p109_final_p90_ceiling_before_m",
                "p109_final_p90_ceiling_after_m",
                "p109_final_lowland_shoulder_area_fraction",
                "p109_final_lowland_shoulder_before_fraction",
                "p109_final_lowland_shoulder_after_fraction",
                "land_fraction_consistent_with_target",
            )
        },
        "terrain_internal_profile": {
            key: metrics.get("terrain_internal_profile", {}).get(key)
            for key in (
                "schema",
                "available",
                "enabled",
                "time_myr",
                "cell_count",
                "total_seconds",
                "profiled_stage_seconds_total",
                "stage_count",
                "top",
                "subprofiles",
            )
        },
        "p120_continental_semantic_geometry": {
            key: metrics.get("p120_continental_semantic_geometry", {}).get(key)
            for key in (
                "schema",
                "continental_land_cell_count",
                "continental_land_area_fraction_world",
                "unsupported_interior_cell_count",
                "unsupported_interior_area_fraction_world",
                "region_fields_available",
                "unsupported_detail_region_elongated_delta",
                "unsupported_internal_region_elongated_delta",
                "unsupported_detail_region_tiny_delta",
                "unsupported_internal_region_tiny_delta",
                "unsupported_terrain_province_elongated_area_fraction",
                "unsupported_continental_detail_raw_elongated_area_fraction",
                "unsupported_continental_detail_region_elongated_area_fraction",
                "unsupported_internal_block_raw_elongated_area_fraction",
                "unsupported_internal_block_region_elongated_area_fraction",
            )
        },
        "p110a_modern_planform": {
            key: metrics.get("p110a_modern_planform", {}).get(key)
            for key in (
                "time_myr",
                "modern_planform_applicable",
                "modern_planform_applicable_after_myr",
                "summary",
                "p110b_final_state_archetype",
                "p110b_lineage_survival",
                "p110b_terminal_supercontinent_diagnostics",
                "ocean_gateway_topology",
                "out_of_envelope_count",
                "out_of_envelope",
                "deferred_out_of_envelope_count",
                "deferred_out_of_envelope",
                "warning_flags",
                "within_p110a_modern_planform_envelope",
            )
        },
        "p110b_historical_supercontinent_trajectory": {
            key: metrics.get(
                "p110b_historical_supercontinent_trajectory", {}
            ).get(key)
            for key in (
                "skipped",
                "skip_reason",
                "usable_frame_count",
                "supercontinent_frame_count",
                "supercontinent_frame_fraction",
                "max_largest_land_component_share",
                "median_largest_land_component_share",
                "final_largest_land_component_share",
                "final_supercontinent_like",
                "max_consecutive_supercontinent_duration_myr",
                "supercontinent_time_window_myr",
                "wide_recurrent_supercontinent_window",
                "long_lived_supercontinent_like",
                "late_recovery_after_long_supercontinent",
            )
        },
        "terminal_active_plate_count": int(metrics["terminal_active_plate_count"]),
        "terminal_major_plate_count": int(metrics["terminal_major_plate_count"]),
        "terminal_minor_plate_count": int(metrics["terminal_minor_plate_count"]),
        "terminal_microplate_count": int(metrics["terminal_microplate_count"]),
        "protected_microplate_count": int(metrics["protected_microplate_count"]),
        "merged_microplate_count": int(metrics["merged_microplate_count"]),
        "boundary_province_object_count_by_kind": metrics[
            "boundary_province_object_count_by_kind"],
        "ridge_transform_network_continuity": metrics[
            "ridge_transform_network_continuity"],
        "boundary_width_summary": metrics.get(
            "boundary_width_diagnostics", {}).get("summary", {}),
        "p1116_boundary_linework": {
            key: metrics.get("p1116_boundary_linework", {}).get(key)
            for key in (
                "ridge",
                "ridge_source",
                "raw_ridge_boundary",
                "object_ridge",
                "depth_ridge_province",
                "transform",
                "transform_source",
                "raw_transform_boundary",
                "object_transform",
                "trench",
                "trench_source",
                "raw_trench_boundary",
                "object_trench",
                "depth_trench_province",
                "convergent",
                "convergent_source",
                "raw_convergent_boundary",
                "object_convergent_parent",
                "ridge_transform_adjacency_fraction",
                "trench_subduction_attachment_fraction",
                "deep_trench_subduction_attachment_fraction",
            )
        },
        "p1116_oceanic_crust_structure": {
            key: metrics.get("p1116_oceanic_crust_structure", {}).get(key)
            for key in (
                "age_depth_correlation",
                "age_distance_from_ridge_correlation",
                "old_minus_young_median_depth_m",
                "far_minus_near_ridge_median_age_myr",
                "nonridge_minus_ridge_median_depth_m",
                "trench_minus_nontrench_median_depth_m",
            )
        },
        "p11166_trench_arc_refinement": metrics.get(
            "p11166_trench_arc_refinement", {}),
        "p111632_parent_trench_linework": metrics.get(
            "p111632_parent_trench_linework", {}),
        "p11167_convergent_parent_linework": metrics.get(
            "p11167_convergent_parent_linework", {}),
        "p11168_orogenic_parent_hierarchy": metrics.get(
            "p11168_orogenic_parent_hierarchy", {}),
        "p111633_orogenic_belt_morphology": {
            key: metrics.get("p111633_orogenic_belt_morphology", {}).get(key)
            for key in (
                "schema",
                "graded_belt_continuity_score",
                "spine_components_per_hierarchy_component",
                "high_peak_fragmentation_pressure",
                "worst_peak_hierarchy_high_subcomponent_count",
                "near_high_saddle_cell_count",
                "near_high_saddle_area_fraction_world",
                "p111633_belt_relief_response_used",
                "p111633_belt_relief_candidate_cell_count",
                "p111633_belt_relief_bridge_cell_count",
                "p111633_belt_relief_1800_component_count_before",
                "p111633_belt_relief_1800_component_count_after",
                "p111633_belt_relief_2400_component_count_before",
                "p111633_belt_relief_2400_component_count_after",
                "p111633_belt_relief_guard_reverted",
            )
        },
        "p111634_spine_aligned_elevation_morphology": {
            key: metrics.get(
                "p111634_spine_aligned_elevation_morphology", {}).get(key)
            for key in (
                "schema",
                "ridge_axis_relief_continuity_score",
                "gradient_order_score",
                "high_component_count",
                "high_components_per_peak_hierarchy_component",
                "high_near_spine_fraction_d2",
                "high_far_from_spine_fraction",
                "spine_relief_2400_coverage_fraction",
                "p111634_spine_aligned_response_used",
                "p111634_axis_raise_candidate_cell_count",
                "p111634_axis_raised_cell_count",
                "p111634_offaxis_high_candidate_cell_count",
                "p111634_offaxis_high_softened_cell_count",
                "p111634_high_component_count_before",
                "p111634_high_component_count_after",
                "p111634_high_near_spine_fraction_before",
                "p111634_high_near_spine_fraction_after",
                "p111634_guard_reverted",
                "p111637c_anti_raster_orogenic_expression_accepted",
                "p111637c_guard_reverted",
                "p111637c_land_mask_preserved",
                "p111637c_candidate_cell_count",
                "p111637c_adjusted_cell_count",
                "p111637c_roughness_before_m",
                "p111637c_roughness_after_m",
                "p111637c_mean_abs_delta_m",
                "p111637c_max_abs_delta_m",
                "p111637c_high_component_count_before",
                "p111637c_high_component_count_after",
                "p134_terminal_spine_relief_gap_bridging_accepted",
                "p134_guard_reverted",
                "p134_land_mask_preserved",
                "p134_candidate_cell_count",
                "p134_adjusted_cell_count",
                "p134_spine_2400_component_count_before",
                "p134_spine_2400_component_count_after",
                "p134_spine_3000_component_count_before",
                "p134_spine_3000_component_count_after",
                "p134_spine_2400_coverage_before",
                "p134_spine_2400_coverage_after",
                "p134_spine_3000_coverage_before",
                "p134_spine_3000_coverage_after",
                "p134_high_component_count_before",
                "p134_high_component_count_after",
                "p135_terminal_crest_core_relief_gap_bridging_accepted",
                "p135_guard_reverted",
                "p135_land_mask_preserved",
                "p135_candidate_cell_count",
                "p135_adjusted_cell_count",
                "p135_gap_component_count",
                "p135_gap_cell_count",
                "p135_crest_3000_component_count_before",
                "p135_crest_3000_component_count_after",
                "p135_spine_3000_component_count_before",
                "p135_spine_3000_component_count_after",
                "p135_high_component_count_before",
                "p135_high_component_count_after",
                "p135_crest_3000_coverage_before",
                "p135_crest_3000_coverage_after",
                "p135_spine_3000_coverage_before",
                "p135_spine_3000_coverage_after",
                "p136_terminal_high_peak_speckle_rebalance_accepted",
                "p136_guard_reverted",
                "p136_land_mask_preserved",
                "p136_candidate_cell_count",
                "p136_adjusted_cell_count",
                "p136_candidate_component_count",
                "p136_lowered_component_count",
                "p136_high_component_count_before",
                "p136_high_component_count_after",
                "p136_spine_3000_component_count_before",
                "p136_spine_3000_component_count_after",
                "p136_crest_3000_component_count_before",
                "p136_crest_3000_component_count_after",
                "p136_spine_2400_component_count_before",
                "p136_spine_2400_component_count_after",
                "p136_high_cell_count_before",
                "p136_high_cell_count_after",
            )
        },
        "p111635_orogenic_spine_linework": {
            key: metrics.get("p111635_orogenic_spine_linework", {}).get(key)
            for key in (
                "schema",
                "linework_continuity_score",
                "spine_components_per_peak_hierarchy_component",
                "short_spine_component_fraction",
                "spine_endpoints_per_component",
                "branch_attachment_fraction",
                "polar_spine_area_fraction",
                "p111635_spine_linework_smoothing_accepted",
                "p111635_spine_components_before_refine",
                "p111635_spine_components_after_refine",
                "p111635_short_spine_components_before_refine",
                "p111635_short_spine_components_after_refine",
                "p111635_branch_attachment_fraction_before",
                "p111635_branch_attachment_fraction_after",
                "p111635_bridge_cell_count",
                "p111635_bridge_candidate_cell_count",
                "p111635_path_count",
                "p111637a_spine_object_promotion_accepted",
                "p111637a_spine_components_before",
                "p111637a_spine_components_after",
                "p111637a_short_spine_components_before",
                "p111637a_short_spine_components_after",
                "p111637a_branch_attachment_fraction_before",
                "p111637a_branch_attachment_fraction_after",
                "p111637a_spine_top3_share_before",
                "p111637a_spine_top3_share_after",
                "p111637a_linework_score_before",
                "p111637a_linework_score_after",
                "p111637a_bridge_cell_count",
                "p111637a_candidate_cell_count",
                "p111637a_path_count",
                "p132_parent_anchor_spine_promotion_accepted",
                "p132_candidate_cell_count",
                "p132_promoted_spine_cell_count",
                "p132_parent_aligned_spine_fraction_before",
                "p132_parent_aligned_spine_fraction_after",
                "p132_linework_score_before",
                "p132_linework_score_after",
                "p132_spine_components_before",
                "p132_spine_components_after",
                "p132_short_spine_components_before",
                "p132_short_spine_components_after",
                "p1145_whole_mask_spine_planner_accepted",
                "p1145_bridge_cell_count",
                "p1145_candidate_cell_count",
                "p1145_path_count",
                "p1145_linework_score_before",
                "p1145_linework_score_after",
                "p1145_spine_components_before",
                "p1145_spine_components_after",
                "p1145_short_spine_components_before",
                "p1145_short_spine_components_after",
                "p1145_branch_attachment_fraction_before",
                "p1145_branch_attachment_fraction_after",
                "p1145_spine_top3_share_before",
                "p1145_spine_top3_share_after",
                "p142_class_aware_combined_score",
                "p142_class_aware_class_score",
                "p142_class_aware_crest_small_area_fraction",
                "p142_class_aware_branch_small_area_fraction",
                "p142_class_aware_class_small_area_fraction",
                "p142_class_aware_crest_component_count",
                "p142_class_aware_branch_component_count",
                "p142_class_aware_blind_spot",
                "p143_class_repair_needed",
                "p143_class_profile_improved",
                "p143_class_score_before",
                "p143_class_score_after",
                "p143_class_small_area_fraction_before",
                "p143_class_small_area_fraction_after",
                "p143_crest_small_area_fraction_before",
                "p143_crest_small_area_fraction_after",
                "p143_branch_small_area_fraction_before",
                "p143_branch_small_area_fraction_after",
                "p144_class_path_option_count",
                "p144_class_path_selected_count",
                "p144_crest_path_option_count",
                "p144_branch_path_option_count",
                "p144_class_attempted_path_count",
                "p144_class_found_path_count",
                "p144_crest_promoted_spine_cell_count",
                "p144_branch_promoted_spine_cell_count",
                "p144_class_promoted_spine_cell_count",
                "p144_class_path_profile_rejected_count",
                "p145_class_hierarchy_component_consolidation_accepted",
                "p145_removed_crest_component_count",
                "p145_removed_branch_component_count",
                "p145_removed_crest_cell_count",
                "p145_removed_branch_cell_count",
                "p145_removed_crest_area_fraction",
                "p145_removed_branch_area_fraction",
                "p145_added_crest_cell_count",
                "p145_added_branch_cell_count",
                "p145_added_crest_area_fraction",
                "p145_added_branch_area_fraction",
                "p145_bridge_path_count",
                "p145_class_score_before",
                "p145_class_score_after",
                "p145_class_small_area_fraction_before",
                "p145_class_small_area_fraction_after",
                "p145_crest_small_area_fraction_before",
                "p145_crest_small_area_fraction_after",
                "p145_branch_small_area_fraction_before",
                "p145_branch_small_area_fraction_after",
                "p1146_reject_code",
                "p1146_reject_reason",
                "p1146_support_component_count",
                "p1146_multi_spine_support_component_count",
                "p1146_attempted_path_count",
                "p1146_found_path_count",
                "p1146_bridge_path_count",
                "p1146_attempted_bridge_cell_count",
                "p1146_trial_linework_score",
                "p1146_trial_spine_components",
                "p1146_trial_short_spine_components",
                "p1146_trial_branch_attachment_fraction",
                "p1146_trial_spine_top3_share",
                "p1147_pair_option_count",
                "p1147_pair_selected_count",
                "p1147_pair_balance_rejected_count",
                "p1147_pair_profile_rejected_count",
                "p1148_terminal_proxy_enabled",
                "p1148_terminal_proxy_score_before",
                "p1148_terminal_proxy_score_after",
                "p1148_terminal_proxy_component_count_before",
                "p1148_terminal_proxy_component_count_after",
                "p1148_terminal_proxy_short_count_before",
                "p1148_terminal_proxy_short_count_after",
                "p117_highland_500_cell_count",
                "p117_highland_1000_cell_count",
                "p117_highland_500_orogenic_support_fraction",
                "p117_highland_1000_orogenic_support_fraction",
                "p117_highland_500_peak_support_fraction",
                "p117_highland_1000_peak_support_fraction",
                "p117_highland_500_mountain_inventory_fraction",
                "p117_highland_1000_mountain_inventory_fraction",
                "p117_non_orogenic_highland_500_cell_count",
                "p117_non_orogenic_highland_1000_cell_count",
                "p117_non_orogenic_highland_500_area_fraction",
                "p117_non_orogenic_highland_1000_area_fraction",
                "p118_non_orogenic_midland_stripe_suppression_accepted",
                "p118_guard_reverted",
                "p118_land_mask_preserved",
                "p118_candidate_cell_count",
                "p118_candidate_area_fraction",
                "p118_adjusted_cell_count",
                "p118_adjusted_area_fraction",
                "p118_unsupported_midland_area_fraction_before",
                "p118_unsupported_midland_area_fraction_after",
                "p118_non_orogenic_highland_500_area_fraction_before",
                "p118_non_orogenic_highland_500_area_fraction_after",
                "p118_mean_lowering_m",
                "p118_max_lowering_m",
                "p118_mean_land_relief_before_m",
                "p118_mean_land_relief_after_m",
                "p118_p50_land_relief_before_m",
                "p118_p50_land_relief_after_m",
                "p118_p90_land_relief_before_m",
                "p118_p90_land_relief_after_m",
                "p119_non_orogenic_interior_anti_raster_smoothing_accepted",
                "p119_guard_reverted",
                "p119_land_mask_preserved",
                "p119_candidate_cell_count",
                "p119_candidate_area_fraction",
                "p119_selected_cell_count",
                "p119_selected_area_fraction",
                "p119_adjusted_cell_count",
                "p119_adjusted_area_fraction",
                "p119_compatible_edge_count",
                "p119_roughness_before_m",
                "p119_roughness_after_m",
                "p119_non_orogenic_highland_500_area_fraction_before",
                "p119_non_orogenic_highland_500_area_fraction_after",
                "p119_mean_abs_delta_m",
                "p119_max_abs_delta_m",
                "p119_p50_land_relief_before_m",
                "p119_p50_land_relief_after_m",
                "p119_p90_land_relief_before_m",
                "p119_p90_land_relief_after_m",
                "p120_continental_semantic_geometry_repair_accepted",
                "p120_guard_reverted",
                "p120_candidate_cell_count",
                "p120_candidate_area_fraction",
                "p120_terrain_changed_area_fraction",
                "p120_internal_changed_area_fraction",
                "p120_detail_changed_area_fraction",
                "p120_terrain_elongated_before",
                "p120_terrain_elongated_after",
                "p120_internal_elongated_before",
                "p120_internal_elongated_after",
                "p120_detail_elongated_before",
                "p120_detail_elongated_after",
                "p120_terrain_tiny_before",
                "p120_terrain_tiny_after",
                "p120_internal_tiny_before",
                "p120_internal_tiny_after",
                "p120_detail_tiny_before",
                "p120_detail_tiny_after",
                "p121_semantic_region_elevation_response_accepted",
                "p121_guard_reverted",
                "p121_land_mask_preserved",
                "p121_candidate_cell_count",
                "p121_candidate_area_fraction",
                "p121_selected_cell_count",
                "p121_selected_area_fraction",
                "p121_adjusted_cell_count",
                "p121_adjusted_area_fraction",
                "p121_compatible_edge_count",
                "p121_roughness_before_m",
                "p121_roughness_after_m",
                "p121_non_orogenic_highland_500_area_fraction_before",
                "p121_non_orogenic_highland_500_area_fraction_after",
                "p121_non_orogenic_highland_1000_area_fraction_before",
                "p121_non_orogenic_highland_1000_area_fraction_after",
                "p121_mean_abs_delta_m",
                "p121_max_abs_delta_m",
                "p121_p50_land_relief_before_m",
                "p121_p50_land_relief_after_m",
                "p121_p90_land_relief_before_m",
                "p121_p90_land_relief_after_m",
                "p1149_orogenic_spine_object_count",
                "p1149_orogenic_spine_system_count",
                "p1149_orogenic_spine_trunk_count",
                "p1149_orogenic_spine_branch_count",
                "p1149_fallback_trunk_count",
                "p1149_orphan_branch_count",
                "p1149_attached_branch_fraction",
                "p1149_mean_branch_count_per_system",
                "p1149_spine_cell_count",
                "p1149_trunk_cell_count",
                "p1149_branch_cell_count",
                "p1149_endpoint_count",
                "p115_spine_repair_candidate_count",
                "p115_viable_spine_repair_candidate_count",
                "p115_trunk_bridge_attempt_count",
                "p115_branch_attachment_attempt_count",
                "p115_trunk_bridge_candidate_count",
                "p115_branch_attachment_candidate_count",
                "p115_rejected_proxy_count",
                "p115_multi_trunk_system_count",
                "p115_detached_branch_component_count",
                "p115_best_proxy_score_delta",
                "p115_best_component_delta",
                "p115_candidate_cell_count",
                "p115_viable_candidate_cell_count",
                "p116_candidate_promotion_enabled",
                "p116_candidate_promotion_accepted",
                "p116_guard_reverted",
                "p116_input_candidate_count",
                "p116_input_viable_candidate_count",
                "p116_considered_candidate_count",
                "p116_selected_candidate_count",
                "p116_applied_candidate_count",
                "p116_rejected_candidate_count",
                "p116_rejected_overlap_count",
                "p116_rejected_support_count",
                "p116_rejected_profile_count",
                "p116_applied_cell_count",
                "p116_applied_area_fraction",
                "p116_area_budget_fraction",
                "p116_linework_score_before",
                "p116_linework_score_after",
                "p116_component_count_before",
                "p116_component_count_after",
                "p116_short_count_before",
                "p116_short_count_after",
                "p116_branch_attachment_fraction_before",
                "p116_branch_attachment_fraction_after",
                "p116_spine_top3_share_before",
                "p116_spine_top3_share_after",
                "p124_orogenic_spine_geometry_regularization_accepted",
                "p124_guard_reverted",
                "p124_candidate_cell_count",
                "p124_added_cell_count",
                "p124_component_count_before",
                "p124_component_count_after",
                "p124_endpoint_count_before",
                "p124_endpoint_count_after",
                "p124_short_component_count_before",
                "p124_short_component_count_after",
                "p124_branch_attachment_fraction_before",
                "p124_branch_attachment_fraction_after",
                "p124_linework_score_before",
                "p124_linework_score_after",
                "p140_terminal_branch_spine_component_thickening_accepted",
                "p140_guard_reverted",
                "p140_candidate_cell_count",
                "p140_promoted_branch_spine_cell_count",
                "p140_grown_branch_spine_component_count",
                "p140_branch_spine_cell_count_before",
                "p140_branch_spine_cell_count_after",
                "p140_branch_spine_component_count_before",
                "p140_branch_spine_component_count_after",
                "p140_branch_spine_small_area_fraction_before",
                "p140_branch_spine_small_area_fraction_after",
                "p140_all_spine_component_count_before",
                "p140_all_spine_component_count_after",
                "p140_junction_count_before",
                "p140_junction_count_after",
                "p140_high_degree_count_before",
                "p140_high_degree_count_after",
                "p140_junction_fraction_before",
                "p140_junction_fraction_after",
                "p140_high_degree_fraction_before",
                "p140_high_degree_fraction_after",
                "p140_branch_attachment_fraction_before",
                "p140_branch_attachment_fraction_after",
                "p140_hierarchy_changed_cell_count",
                "p140_crest_spine_changed_cell_count",
                "p141_terminal_crest_spine_component_thickening_accepted",
                "p141_guard_reverted",
                "p141_candidate_cell_count",
                "p141_promoted_crest_spine_cell_count",
                "p141_connected_crest_spine_component_count",
                "p141_grown_crest_spine_component_count",
                "p141_crest_spine_cell_count_before",
                "p141_crest_spine_cell_count_after",
                "p141_crest_spine_component_count_before",
                "p141_crest_spine_component_count_after",
                "p141_crest_spine_small_area_fraction_before",
                "p141_crest_spine_small_area_fraction_after",
                "p141_all_spine_component_count_before",
                "p141_all_spine_component_count_after",
                "p141_junction_count_before",
                "p141_junction_count_after",
                "p141_high_degree_count_before",
                "p141_high_degree_count_after",
                "p141_junction_fraction_before",
                "p141_junction_fraction_after",
                "p141_high_degree_fraction_before",
                "p141_high_degree_fraction_after",
                "p141_branch_attachment_fraction_before",
                "p141_branch_attachment_fraction_after",
                "p141_hierarchy_changed_cell_count",
                "p141_branch_spine_changed_cell_count",
                "p164_terminal_orogenic_axis_skeleton_simplification_accepted",
                "p164_guard_reverted",
                "p164_candidate_cell_count",
                "p164_demoted_cell_count",
                "p164_off_axis_candidate_cell_count",
                "p164_spine_cell_count_before",
                "p164_spine_cell_count_after",
                "p164_component_count_before",
                "p164_component_count_after",
                "p164_short_component_count_before",
                "p164_short_component_count_after",
                "p164_endpoint_count_before",
                "p164_endpoint_count_after",
                "p164_junction_count_before",
                "p164_junction_count_after",
                "p164_high_degree_count_before",
                "p164_high_degree_count_after",
                "p164_junction_fraction_before",
                "p164_junction_fraction_after",
                "p164_high_degree_fraction_before",
                "p164_high_degree_fraction_after",
                "p164_branch_attachment_fraction_before",
                "p164_branch_attachment_fraction_after",
                "p164_linework_score_before",
                "p164_linework_score_after",
                "p164_axis_score_before",
                "p164_axis_score_after",
                "p164_protected_extreme_demoted_cell_count",
                "p164_reject_code",
                "p111637a_terminal_spine_land_consistency_accepted",
                "p111637a_terminal_submerged_spine_cell_count_before",
                "p111637a_terminal_submerged_spine_cell_count_after",
                "p111637a_terminal_spine_lift_cell_count",
                "p111637a_terminal_spine_lift_area_fraction",
                "p111637b_final_spine_elevation_consistency_accepted",
                "p111637b_guard_reverted",
                "p111637b_candidate_cell_count",
                "p111637b_adjusted_cell_count",
                "p111637b_submerged_spine_cell_count_before",
                "p111637b_submerged_spine_cell_count_after",
                "p111637b_bridge_blocked_submerged_spine_cell_count",
                "p111637b_bridge_blocked_submerged_spine_area_fraction",
                "p111637b_underexpressed_spine_cell_count_before",
                "p111637b_underexpressed_spine_cell_count_after",
                "p111637b_land_delta_area_fraction",
                "p111637b_mean_lift_m",
                "p111637b_max_lift_m",
                "p111637b_linework_score_before",
                "p111637b_linework_score_after",
                "p1138_terminal_spine_stub_cleanup_accepted",
                "p1138_guard_reverted",
                "p1138_candidate_cell_count",
                "p1138_demoted_cell_count",
                "p1138_component_count_before",
                "p1138_component_count_after",
                "p1138_short_component_count_before",
                "p1138_short_component_count_after",
                "p1138_linework_score_before",
                "p1138_linework_score_after",
                "p1139_terminal_spine_leaf_tip_cleanup_accepted",
                "p1139_guard_reverted",
                "p1139_candidate_cell_count",
                "p1139_demoted_cell_count",
                "p1139_endpoint_count_before",
                "p1139_endpoint_count_after",
                "p1139_component_count_before",
                "p1139_component_count_after",
                "p1139_short_component_count_before",
                "p1139_short_component_count_after",
                "p1139_linework_score_before",
                "p1139_linework_score_after",
                "p1140_terminal_spine_triangular_kink_cleanup_accepted",
                "p1140_guard_reverted",
                "p1140_candidate_cell_count",
                "p1140_demoted_cell_count",
                "p1140_redundant_kink_count_before",
                "p1140_redundant_kink_count_after",
                "p1140_endpoint_count_before",
                "p1140_endpoint_count_after",
                "p1140_component_count_before",
                "p1140_component_count_after",
                "p1140_short_component_count_before",
                "p1140_short_component_count_after",
                "p1140_linework_score_before",
                "p1140_linework_score_after",
            )
        },
        "p111636_polar_edge_orogen_overclassification": {
            key: metrics.get(
                "p111636_polar_edge_orogen_overclassification", {}).get(key)
            for key in (
                "schema",
                "polar_peak_area_fraction_world",
                "polar_spine_area_fraction_world",
                "polar_peak_to_spine_width_ratio",
                "edge_peak_component_count",
                "edge_short_peak_component_count",
                "p111636_polar_edge_refinement_accepted",
                "p111636_polar_peak_area_fraction_before",
                "p111636_polar_peak_area_fraction_after",
                "p111636_polar_peak_area_fraction_delta",
                "p111636_edge_peak_component_count_before",
                "p111636_edge_peak_component_count_after",
                "p111636_edge_peak_component_count_delta",
                "p111636_candidate_cell_count",
                "p111636_reclassified_cell_count",
                "p111636_polar_reclassified_cell_count",
                "p111636_edge_reclassified_cell_count",
                "p111636_extreme_reclassified_cell_count",
                "p111636_no_extreme_reclassification",
                "p111636_nonexpansive_polar_edge",
                "p160_polar_edge_spine_thinning_accepted",
                "p160_guard_reverted",
                "p160_candidate_cell_count",
                "p160_candidate_area_fraction",
                "p160_demoted_cell_count",
                "p160_demoted_area_fraction",
                "p160_polar_demoted_cell_count",
                "p160_edge_demoted_cell_count",
                "p160_extreme_demoted_cell_count",
                "p160_polar_spine_area_fraction_before",
                "p160_polar_spine_area_fraction_after",
                "p160_polar_spine_area_fraction_delta",
                "p160_edge_spine_component_count_before",
                "p160_edge_spine_component_count_after",
                "p160_edge_spine_component_count_delta",
                "p160_spine_component_count_before",
                "p160_spine_component_count_after",
                "p160_short_spine_component_count_before",
                "p160_short_spine_component_count_after",
                "p160_branch_attachment_fraction_before",
                "p160_branch_attachment_fraction_after",
                "p160_linework_score_before",
                "p160_linework_score_after",
                "p160_no_extreme_demotion",
                "p160_nonexpansive_polar_edge_spine",
                "p165_terminal_polar_edge_orogenic_semantic_compaction_accepted",
                "p165_guard_reverted",
                "p165_candidate_cell_count",
                "p165_candidate_area_fraction",
                "p165_peak_reclassified_cell_count",
                "p165_peak_reclassified_area_fraction",
                "p165_fringe_cleared_cell_count",
                "p165_fringe_cleared_area_fraction",
                "p165_polar_semantic_area_fraction_before",
                "p165_polar_semantic_area_fraction_after",
                "p165_polar_semantic_area_fraction_delta",
                "p165_edge_semantic_area_fraction_before",
                "p165_edge_semantic_area_fraction_after",
                "p165_edge_semantic_area_fraction_delta",
                "p165_polar_peak_area_fraction_before",
                "p165_polar_peak_area_fraction_after",
                "p165_polar_peak_area_fraction_delta",
                "p165_edge_peak_component_count_before",
                "p165_edge_peak_component_count_after",
                "p165_edge_peak_component_count_delta",
                "p165_peak_component_count_before",
                "p165_peak_component_count_after",
                "p165_hierarchy_changed_cell_count",
                "p165_halo_changed_cell_count",
                "p165_apron_changed_cell_count",
                "p165_protected_spine_changed_cell_count",
                "p165_protected_extreme_reclassified_cell_count",
                "p165_reject_code",
                "p165_no_protected_spine_change",
                "p165_no_extreme_reclassification",
                "p165_nonexpansive_polar_edge_semantics",
            )
        },
        "p166167_terminal_relief_planform": {
            key: metrics.get("p166167_terminal_relief_planform", {}).get(key)
            for key in (
                "schema",
                "straight_highland_component_count",
                "straight_highland_area_fraction_world",
                "straight_highland_score",
                "inland_plateau_area_fraction_world",
                "broad_inland_plateau_component_count",
                "largest_inland_plateau_component_share",
                "inland_plateau_relief_iqr_m",
                "p166_terminal_straight_highland_relief_softening_accepted",
                "p166_guard_reverted",
                "p166_land_mask_preserved",
                "p166_candidate_cell_count",
                "p166_adjusted_cell_count",
                "p166_straightness_score_before",
                "p166_straightness_score_after",
                "p166_straight_highland_area_fraction_before",
                "p166_straight_highland_area_fraction_after",
                "p166_protected_core_lowered_cell_count",
                "p166_extreme_softened_cell_count",
                "p166_reject_code",
                "p167_terminal_inland_plateau_diversity_repair_accepted",
                "p167_guard_reverted",
                "p167_land_mask_preserved",
                "p167_candidate_cell_count",
                "p167_adjusted_cell_count",
                "p167_inland_plateau_area_fraction_before",
                "p167_inland_plateau_area_fraction_after",
                "p167_largest_plateau_component_share_before",
                "p167_largest_plateau_component_share_after",
                "p167_plateau_relief_iqr_before_m",
                "p167_plateau_relief_iqr_after_m",
                "p167_protected_core_lowered_cell_count",
                "p167_extreme_softened_cell_count",
                "p167_reject_code",
                "p168_terminal_surface_derasterization_and_inland_diversity_accepted",
                "p168_guard_reverted",
                "p168_land_mask_preserved",
                "p168_candidate_cell_count",
                "p168_adjusted_cell_count",
                "p168_land_adjusted_cell_count",
                "p168_ocean_adjusted_cell_count",
                "p168_weak_land_iqr_before_m",
                "p168_weak_land_iqr_after_m",
                "p168_weak_land_edge_delta_p95_before_m",
                "p168_weak_land_edge_delta_p95_after_m",
                "p168_weak_ocean_edge_delta_p95_before_m",
                "p168_weak_ocean_edge_delta_p95_after_m",
                "p168_broad_plateau_area_fraction_before",
                "p168_broad_plateau_area_fraction_after",
                "p168_protected_core_changed_cell_count",
                "p168_ridge_trench_changed_cell_count",
                "p168_extreme_changed_cell_count",
                "p168_reject_code",
            )
        },
        "collision_plateau_area_fraction": float(
            metrics["collision_plateau_area_fraction"]),
        "high_mountain_coherence": {
            key: metrics.get("high_mountain_coherence", {}).get(key)
            for key in (
                "high_mountain_area_fraction_world",
                "high_mountain_component_count",
                "largest_high_mountain_component_share",
                "top3_high_mountain_component_share",
                "small_high_mountain_component_area_fraction_of_high",
                "high_mountain_fragmentation_index",
                "isolated_extreme_peak_area_fraction_of_extreme",
                "high_mountain_parent_overlap_fraction",
                "high_mountain_mountain_range_overlap_fraction",
                "small_high_mountain_component_count",
                "small_unparented_high_mountain_component_count",
                "p1113_bridge_area_fraction",
                "p1113_high_component_count_before",
                "p1113_high_component_count_after",
                "p1113_top3_high_share_before",
                "p1113_top3_high_share_after",
                "p1113_fragmentation_index_before",
                "p1113_fragmentation_index_after",
                "p1113_range_count_considered",
                "p1113_land_mask_preserved",
                "p11169_hierarchy_used",
                "p11169_hierarchy_bridge_area_fraction",
                "p11169_hierarchy_bridge_cell_count",
                "p11169_legacy_bridge_area_fraction",
                "p11169_foreland_peak_cell_count",
                "p11169_fragmentation_delta",
                "p11169_crest_high_overlap_fraction",
                "p11169_bridge_restricted_to_peak_hierarchy",
                "p11169_hypsometry_guard_reverted",
                "p153_terminal_high_fleck_cleanup_accepted",
                "p153_candidate_cell_count",
                "p153_softened_cell_count",
                "p153_high_component_count_before",
                "p153_high_component_count_after",
                "p153_fragmentation_index_before",
                "p153_fragmentation_index_after",
                "p153_guard_reverted",
                "p155_terminal_high_relief_consistency_gate_accepted",
                "p155_guard_reverted",
                "p155_land_mask_preserved",
                "p155_candidate_cell_count",
                "p155_candidate_area_fraction",
                "p155_raised_cell_count",
                "p155_raised_area_fraction",
                "p155_softened_cell_count",
                "p155_softened_area_fraction",
                "p155_high_component_count_before",
                "p155_high_component_count_after",
                "p155_spine_3000_component_count_before",
                "p155_spine_3000_component_count_after",
                "p155_spine_3000_coverage_before",
                "p155_spine_3000_coverage_after",
                "p155_top3_high_share_before",
                "p155_top3_high_share_after",
                "p155_fragmentation_index_before",
                "p155_fragmentation_index_after",
                "p155_parent_overlap_before",
                "p155_parent_overlap_after",
                "p155_high_area_fraction_before",
                "p155_high_area_fraction_after",
                "p155_p90_land_relief_before_m",
                "p155_p90_land_relief_after_m",
                "p155_p98_land_relief_before_m",
                "p155_p98_land_relief_after_m",
                "p155_reject_code",
                "p162_terminal_high_mountain_fragment_cleanup_accepted",
                "p162_guard_reverted",
                "p162_land_mask_preserved",
                "p162_candidate_cell_count",
                "p162_candidate_area_fraction",
                "p162_softened_cell_count",
                "p162_softened_area_fraction",
                "p162_high_component_count_before",
                "p162_high_component_count_after",
                "p162_small_high_component_count_before",
                "p162_small_high_component_count_after",
                "p162_spine_3000_component_count_before",
                "p162_spine_3000_component_count_after",
                "p162_top3_high_share_before",
                "p162_top3_high_share_after",
                "p162_fragmentation_index_before",
                "p162_fragmentation_index_after",
                "p162_high_area_fraction_before",
                "p162_high_area_fraction_after",
                "p162_p90_land_relief_before_m",
                "p162_p90_land_relief_after_m",
                "p162_p98_land_relief_before_m",
                "p162_p98_land_relief_after_m",
                "p162_mean_lower_m",
                "p162_max_lower_m",
                "p162_extreme_softened_cell_count",
                "p162_reject_code",
                "p156_path_candidate_count",
                "p156_path_candidate_cell_count",
                "p156_selected_path_count",
                "p156_selected_path_cell_count",
                "p157_seed_path_candidate_count",
                "p157_seed_path_candidate_cell_count",
                "p157_selected_seed_path_count",
                "p157_selected_seed_path_cell_count",
                "p111611_bridge_candidate_cell_count",
                "p111611_bridge_candidate_area_fraction",
                "p111611_peak_hierarchy_shoulder_cell_count",
                "p111611_high_pair_count",
                "p111611_safe_path_count",
                "p111611_blocked_high_pair_count",
                "p111611_diagnostic_system_count",
                "p111611_bridge_deferred_no_safe_path",
            )
        },
        "active_margin_arc_trench_adjacency_fraction": float(
            metrics["active_margin_arc_trench_adjacency_fraction"]),
        "island_arc_chain_count": int(metrics["island_arc_chain_count"]),
        "microcontinent_object_count": int(metrics["microcontinent_object_count"]),
        "parented_oceanic_island_chain_count": int(
            metrics["parented_oceanic_island_chain_count"]),
        "deep_trench_fraction_below_6000m": float(
            metrics["deep_trench_fraction_below_6000m"]),
        "p1115_unsupported_open_ocean_shoal_fraction": float(
            metrics["p1115_unsupported_open_ocean_shoal_fraction"]),
        "p1115_object_backed_open_ocean_shoal_fraction": float(
            metrics["p1115_object_backed_open_ocean_shoal_fraction"]),
        "p1115_object_backed_ocean_relief_fraction": float(
            metrics["p1115_object_backed_ocean_relief_fraction"]),
        "p1115_emerged_object_backed_island_fraction_world": float(
            metrics["p1115_emerged_object_backed_island_fraction_world"]),
        "array_archive_present": bool(archive.get("path") and Path(archive["path"]).exists()),
        "asset_count": int(len(metrics.get("assets", {}))),
    }


def _render_contact_sheet(world: Any, outdir: Path,
                          assets: dict[str, str]) -> Path | None:
    from aevum import render

    grid = world.grid
    sea = float(world.sea_level)
    elev = np.asarray(world.get_field("terrain.elevation_m", 0.0), dtype=np.float64)
    rel_r = render.to_raster_continuous(grid, elev - sea, preserve_sign=True)
    plate_r = render.to_raster(grid, world.get_field("tectonics.plate_id", -1.0))
    crust_age_r = render.to_raster(grid, world.get_field("crust.age_myr", 0.0))
    depth = np.maximum(sea - elev, 0.0)
    bathy_r = render.to_raster_continuous(
        grid,
        np.where(elev < sea, depth, np.nan),
    )
    object_mask = np.zeros(grid.n, dtype=float)
    for idx, (_, object_set, kinds) in enumerate(P107_OBJECT_MASK_SPECS, start=1):
        mask = _object_mask(world, object_set, kinds)
        object_mask[mask] = float(idx)
    object_r = render.to_raster(grid, object_mask)
    earth_r = None
    if "earth_reference_same_grid.png" in assets:
        earth = _sample_etopo5_to_grid(grid)
        if earth is not None:
            earth_r = render.to_raster_continuous(
                grid, earth, preserve_sign=True)

    panels = [
        ("Aevum elevation", rel_r, "elevation"),
        ("Aevum plates", plate_r, "tab20"),
        ("Bathymetry depth", bathy_r, "bathy"),
        ("Object masks", object_r, "tab20"),
        ("Crust age", crust_age_r, "viridis"),
        ("Earth ETOPO5 same grid", earth_r, "elevation"),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(12, 9), constrained_layout=True)
    for ax, (title, raster, mode) in zip(axes.ravel(), panels):
        if raster is None:
            ax.axis("off")
            continue
        if mode == "elevation":
            im = render.render_elevation_raster(ax, raster, title=title)
        elif mode == "bathy":
            cmap = plt.get_cmap("cividis_r").copy()
            cmap.set_bad("#efe8d7")
            valid = raster[np.isfinite(raster)]
            vmax = max(float(np.percentile(valid, 98)) if valid.size else 1.0, 1000.0)
            im = ax.imshow(np.ma.masked_invalid(raster), cmap=cmap, vmin=0.0, vmax=vmax,
                           extent=[-180, 180, -90, 90])
            ax.set_title(title)
        else:
            im = ax.imshow(raster, cmap=mode, extent=[-180, 180, -90, 90])
            ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.72)
    path = outdir / "p107_contact_sheet.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_earth_reference_same_grid(world: Any, outdir: Path) -> Path | None:
    from aevum import render

    sampled = _sample_etopo5_to_grid(world.grid)
    if sampled is None:
        return None
    raster = render.to_raster(world.grid, sampled)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = render.render_elevation_raster(
        ax, raster, title=f"Earth ETOPO5 sampled to {world.grid.n} Aevum cells")
    render.add_elevation_colorbar(fig, ax, im)
    path = outdir / "earth_reference_same_grid.png"
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return path


def _sample_etopo5_to_grid(grid: Any) -> np.ndarray | None:
    path = Path(__file__).resolve().parents[2] / "data" / "reference" / "etopo5" / "ETOPO5.DAT"
    if not path.exists():
        return None
    raw = np.fromfile(path, dtype=">i2")
    if raw.size != 2160 * 4320:
        return None
    elev = raw.reshape(2160, 4320).astype(np.float32)
    lat_step = 180.0 / 2160.0
    lon_step = 360.0 / 4320.0
    lat_idx = np.rint((90.0 - lat_step / 2.0 - grid.lat) / lat_step).astype(np.int64)
    lon = (grid.lon + 360.0) % 360.0
    lon_idx = np.rint((lon - lon_step / 2.0) / lon_step).astype(np.int64)
    lat_idx = np.clip(lat_idx, 0, 2159)
    lon_idx %= 4320
    return elev[lat_idx, lon_idx].astype(np.float32)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return str(o)
