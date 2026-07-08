"""P170 historical land and ocean geomorphology diagnostics.

This module is intentionally read-only.  It measures whether archived frames
already contain process-time land and ocean detail, instead of only looking
plausible after terminal terrain and bathymetry polish.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from aevum.modules.tectonics import (
    DOMAIN_ACCRETED_TERRANE,
    DOMAIN_CRATON,
    DOMAIN_LIP,
    ORIGIN_PLUME_IMPACT,
)
from aevum.modules.terrain import (
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATFORM,
    CONT_DETAIL_PLATEAU,
    CONT_DETAIL_RIFT_BASIN,
    CONT_DETAIL_SHIELD,
    CONT_PROVINCE_ACTIVE_OROGEN,
    CONT_PROVINCE_FORELAND_BASIN,
    CONT_PROVINCE_INTRACRATONIC_BASIN,
    CONT_PROVINCE_OLD_OROGEN,
    CONT_PROVINCE_PASSIVE_MARGIN_LOWLAND,
    CONT_PROVINCE_PLATFORM,
    CONT_PROVINCE_VOLCANIC_LIP_PLATEAU,
    INLAND_PROVINCE_OLD_OROGEN,
    INLAND_PROVINCE_PLATFORM,
    INLAND_PROVINCE_PLATFORM_SWELL,
    INLAND_PROVINCE_RIFT,
    INLAND_PROVINCE_SAG_BASIN,
    INLAND_PROVINCE_SHIELD,
    OCEAN_DEPTH_ABYSS,
    OCEAN_DEPTH_RIDGE,
    OCEAN_DEPTH_TRENCH,
    RIFT_MARGIN_STAGE_PASSIVE_LOWLAND,
    RIFT_MARGIN_STAGE_RIFT_BASIN,
)


SCHEMA = "aevum.p170_historical_geomorphology.v1"

LAND_SEMANTIC_METRIC_KEYS = (
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
)

LAND_GLOBAL_TELEMETRY_METRICS = {
    "p104f_pre_p174_lowland_prep_area_fraction": (
        "terrain.last_p104f_pre_p174_lowland_prep_area_fraction"
    ),
    "p104f_pre_p174_lowland_prep_mean_lowering_m": (
        "terrain.last_p104f_pre_p174_lowland_prep_mean_lowering_m"
    ),
    "p104f_pre_p174_lowland_support_area_fraction": (
        "terrain.last_p104f_pre_p174_lowland_support_area_fraction"
    ),
    "p104f_pre_p174_lowland_support_largest_component_fraction": (
        "terrain.last_p104f_pre_p174_lowland_support_largest_component_fraction"
    ),
    "p104f_pre_p174_lowland_memory_seed_area_fraction": (
        "terrain.last_p104f_pre_p174_lowland_memory_seed_area_fraction"
    ),
    "p174_lowland_plain_response_area_fraction": (
        "terrain.last_p174_lowland_plain_response_area_fraction"
    ),
    "p174_lowland_plain_candidate_area_fraction": (
        "terrain.last_p174_lowland_plain_candidate_area_fraction"
    ),
    "p174_lowland_plain_parent_area_fraction": (
        "terrain.last_p174_lowland_plain_parent_area_fraction"
    ),
    "p174_lowland_plain_continuity_memory_area_fraction": (
        "terrain.last_p174_lowland_plain_continuity_memory_area_fraction"
    ),
    "p174_lowland_plain_continuity_parent_area_fraction": (
        "terrain.last_p174_lowland_plain_continuity_parent_area_fraction"
    ),
    "p174_lowland_plain_response_mean_abs_delta_m": (
        "terrain.last_p174_lowland_plain_response_mean_abs_delta_m"
    ),
    "p174_lowland_plain_fraction_before": (
        "terrain.last_p174_lowland_plain_fraction_before"
    ),
    "p174_lowland_plain_fraction_after": (
        "terrain.last_p174_lowland_plain_fraction_after"
    ),
    "p174_lowland_plain_largest_component_fraction_after": (
        "terrain.last_p174_lowland_plain_largest_component_fraction_after"
    ),
    "p174_lowland_plain_parented_fraction_after": (
        "terrain.last_p174_lowland_plain_parented_fraction_after"
    ),
    "p174_lowland_plain_response_stage_code": (
        "terrain.last_p174_lowland_plain_response_stage_code"
    ),
    "p174_support_component_response_area_fraction": (
        "terrain.last_p174_support_component_response_area_fraction"
    ),
    "p174_support_component_response_largest_component_fraction": (
        "terrain.last_p174_support_component_response_largest_component_fraction"
    ),
    "p174_component_stitch_response_area_fraction": (
        "terrain.last_p174_component_stitch_response_area_fraction"
    ),
    "p174_component_stitch_response_largest_component_fraction": (
        "terrain.last_p174_component_stitch_response_largest_component_fraction"
    ),
    "p174_component_stitch_multiring_response_area_fraction": (
        "terrain.last_p174_component_stitch_multiring_response_area_fraction"
    ),
    "p174_component_stitch_multiring_response_largest_component_fraction": (
        "terrain.last_p174_component_stitch_multiring_response_largest_component_fraction"
    ),
    "p174_component_gap_repair_response_area_fraction": (
        "terrain.last_p174_component_gap_repair_response_area_fraction"
    ),
    "p174_component_gap_repair_response_largest_component_fraction": (
        "terrain.last_p174_component_gap_repair_response_largest_component_fraction"
    ),
    "p174_component_gap_repair_gate_active": (
        "terrain.last_p174_component_gap_repair_gate_active"
    ),
    "p174_component_gap_repair_plain_fraction_before": (
        "terrain.last_p174_component_gap_repair_plain_fraction_before"
    ),
    "p174_component_gap_repair_largest_component_fraction_before": (
        "terrain.last_p174_component_gap_repair_largest_component_fraction_before"
    ),
    "p174_component_gap_repair_component_gap_fraction": (
        "terrain.last_p174_component_gap_repair_component_gap_fraction"
    ),
    "p174_component_gap_repair_candidate_area_fraction": (
        "terrain.last_p174_component_gap_repair_candidate_area_fraction"
    ),
    "p174_component_gap_repair_infill_area_fraction": (
        "terrain.last_p174_component_gap_repair_infill_area_fraction"
    ),
    "p174_component_gap_repair_candidate_largest_component_fraction": (
        "terrain.last_p174_component_gap_repair_candidate_largest_component_fraction"
    ),
    "p174_component_gap_repair_diagnostic_connector_domain_area_fraction": (
        "terrain.last_p174_component_gap_repair_diagnostic_connector_domain_area_fraction"
    ),
    "p174_component_gap_repair_diagnostic_connector_response_area_fraction": (
        "terrain.last_p174_component_gap_repair_diagnostic_connector_response_area_fraction"
    ),
    "p174_component_gap_repair_diagnostic_connector_path_count": (
        "terrain.last_p174_component_gap_repair_diagnostic_connector_path_count"
    ),
    "p174_component_gap_repair_diagnostic_connector_target_component_area_fraction": (
        "terrain.last_p174_component_gap_repair_diagnostic_connector_target_component_area_fraction"
    ),
    "p174_component_gap_repair_largest_growth_domain_area_fraction": (
        "terrain.last_p174_component_gap_repair_largest_growth_domain_area_fraction"
    ),
    "p174_component_gap_repair_largest_growth_picked_area_fraction": (
        "terrain.last_p174_component_gap_repair_largest_growth_picked_area_fraction"
    ),
    "p174_component_gap_repair_graph_component_count": (
        "terrain.last_p174_component_gap_repair_graph_component_count"
    ),
    "p174_component_gap_repair_accepted_component_count": (
        "terrain.last_p174_component_gap_repair_accepted_component_count"
    ),
    "p174_component_gap_repair_accepted_candidate_area_fraction": (
        "terrain.last_p174_component_gap_repair_accepted_candidate_area_fraction"
    ),
    "p174_component_gap_repair_budget_area_fraction": (
        "terrain.last_p174_component_gap_repair_budget_area_fraction"
    ),
    "p174_component_gap_repair_reject_empty_count": (
        "terrain.last_p174_component_gap_repair_reject_empty_count"
    ),
    "p174_component_gap_repair_reject_too_small_count": (
        "terrain.last_p174_component_gap_repair_reject_too_small_count"
    ),
    "p174_component_gap_repair_reject_too_large_count": (
        "terrain.last_p174_component_gap_repair_reject_too_large_count"
    ),
    "p174_component_gap_repair_reject_single_label_count": (
        "terrain.last_p174_component_gap_repair_reject_single_label_count"
    ),
    "p174_area_topup_source_area_fraction": (
        "terrain.last_p174_area_topup_source_area_fraction"
    ),
    "p174_area_topup_domain_area_fraction": (
        "terrain.last_p174_area_topup_domain_area_fraction"
    ),
    "p174_area_topup_candidate_area_fraction": (
        "terrain.last_p174_area_topup_candidate_area_fraction"
    ),
    "p174_area_topup_response_area_fraction": (
        "terrain.last_p174_area_topup_response_area_fraction"
    ),
    "p174_area_topup_response_largest_component_fraction": (
        "terrain.last_p174_area_topup_response_largest_component_fraction"
    ),
}

LAND_METRIC_KEYS = (
    LAND_SEMANTIC_METRIC_KEYS + tuple(LAND_GLOBAL_TELEMETRY_METRICS.keys())
)

OCEAN_METRIC_KEYS = (
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


def historical_geomorphology_summary(
    world: Any,
    archive: Any,
    *,
    inland_width_steps: float = 3.0,
    ordinary_plateau_min_m: float = 900.0,
    ordinary_plateau_max_m: float = 2800.0,
) -> dict[str, Any]:
    """Summarize land/ocean geomorphology detail through saved archive frames."""
    grid = getattr(world, "grid", None) or getattr(getattr(archive, "world", None), "grid", None)
    frames = list(getattr(archive, "frames", []) or [])
    if grid is None:
        return _skipped_summary("archive_grid_not_available", len(frames))

    rows: list[dict[str, Any]] = []
    missing_field_counts: dict[str, int] = {}
    for frame in frames:
        row = frame_geomorphology_metrics(
            grid,
            frame,
            default_sea_level=float(getattr(world, "sea_level", 0.0)),
            inland_width_steps=float(inland_width_steps),
            ordinary_plateau_min_m=float(ordinary_plateau_min_m),
            ordinary_plateau_max_m=float(ordinary_plateau_max_m),
        )
        if bool(row.get("usable", False)):
            rows.append(row)
        for field in row.get("missing_fields", ()):
            missing_field_counts[str(field)] = missing_field_counts.get(str(field), 0) + 1

    rows.sort(key=lambda row: float(row["time_myr"]))
    required_present = all(
        all(key in row.get("land_metrics", {}) for key in LAND_METRIC_KEYS)
        and all(key in row.get("ocean_metrics", {}) for key in OCEAN_METRIC_KEYS)
        for row in rows
    )
    land_flags = [row for row in rows if bool(row.get("diagnostic_flags", {}).get("ordinary_plateau_like"))]
    lowland_flags = [row for row in rows if bool(row.get("diagnostic_flags", {}).get("lowland_plain_deficient"))]
    ocean_flags = [row for row in rows if bool(row.get("diagnostic_flags", {}).get("ordinary_deep_ocean_like"))]

    return {
        "schema": SCHEMA,
        "context": {
            "spec_name": str(getattr(getattr(world, "spec", None), "name", "")),
            "seed": int(getattr(getattr(world, "spec", None), "seed", 0)),
            "cells": int(getattr(grid, "n", 0)),
            "final_time_myr": float(getattr(world, "time_myr", 0.0)),
        },
        "parameters": {
            "inland_width_steps": float(inland_width_steps),
            "ordinary_plateau_min_m": float(ordinary_plateau_min_m),
            "ordinary_plateau_max_m": float(ordinary_plateau_max_m),
        },
        "frame_count": int(len(frames)),
        "usable_frame_count": int(len(rows)),
        "required_land_metric_keys": list(LAND_METRIC_KEYS),
        "required_ocean_metric_keys": list(OCEAN_METRIC_KEYS),
        "required_metric_keys_present": bool(required_present),
        "missing_field_counts": dict(sorted(missing_field_counts.items())),
        "metric_extremes": _metric_extremes(rows),
        "p174_continuity": _p174_continuity_summary(rows),
        "flag_summary": {
            "ordinary_plateau_frame_count": int(len(land_flags)),
            "lowland_plain_deficient_frame_count": int(len(lowland_flags)),
            "ordinary_deep_ocean_frame_count": int(len(ocean_flags)),
            "ordinary_plateau_time_windows_myr": _flag_windows(rows, "ordinary_plateau_like"),
            "lowland_plain_deficient_time_windows_myr": _flag_windows(rows, "lowland_plain_deficient"),
            "ordinary_deep_ocean_time_windows_myr": _flag_windows(rows, "ordinary_deep_ocean_like"),
            "land_terminal_quality_jump": _terminal_quality_jump(
                rows, "land_metrics", "inland_detail_entropy"),
            "ocean_terminal_quality_jump": _terminal_quality_jump(
                rows, "ocean_metrics", "ocean_fabric_entropy"),
        },
        "acceptance": {
            "audit_completed": bool(len(rows) > 0),
            "required_metric_keys_present": bool(required_present),
            "generation_behavior_changed": False,
            "identifies_land_time_windows": bool(len(land_flags) > 0),
            "identifies_ocean_time_windows": bool(len(ocean_flags) > 0),
        },
        "frame_rows": rows,
    }


def write_historical_geomorphology_audit(
    world: Any,
    archive: Any,
    outdir: str | Path,
    *,
    filename: str = "p170_historical_geomorphology_audit.json",
) -> dict[str, Any]:
    """Write the P170 summary JSON and return the same dictionary."""
    summary = historical_geomorphology_summary(world, archive)
    path = Path(outdir)
    path.mkdir(parents=True, exist_ok=True)
    (path / filename).write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return summary


def frame_geomorphology_metrics(
    grid: Any,
    frame: Any,
    *,
    default_sea_level: float = 0.0,
    inland_width_steps: float = 3.0,
    ordinary_plateau_min_m: float = 900.0,
    ordinary_plateau_max_m: float = 2800.0,
) -> dict[str, Any]:
    """Return P170 metrics for one archive frame."""
    fields = dict(getattr(frame, "fields", {}) or {})
    globals_ = dict(getattr(frame, "globals", {}) or {})
    n = int(getattr(grid, "n", 0))
    missing: list[str] = []
    elev, elev_ok = _field(fields, "terrain.elevation_m", n, missing)
    if not elev_ok:
        return {
            "time_myr": float(getattr(frame, "time_myr", 0.0)),
            "usable": False,
            "skip_reason": "terrain.elevation_m_missing_or_invalid",
            "missing_fields": missing,
            "land_metrics": _zero_metrics(LAND_METRIC_KEYS),
            "ocean_metrics": _zero_metrics(OCEAN_METRIC_KEYS),
            "diagnostic_flags": {},
        }

    sea_level = float(globals_.get("ocean.sea_level_m", default_sea_level))
    rel = elev - sea_level
    land = rel >= 0.0
    ocean = ~land
    area = np.asarray(getattr(grid, "cell_area", np.ones(n)), dtype=np.float64)
    total_area = max(float(area.sum()), 1.0e-12)
    time_myr = float(getattr(frame, "time_myr", 0.0))

    land_metrics = _land_metrics(
        grid,
        fields,
        globals_,
        missing,
        area,
        rel,
        land,
        time_myr,
        inland_width_steps=float(inland_width_steps),
        ordinary_plateau_min_m=float(ordinary_plateau_min_m),
        ordinary_plateau_max_m=float(ordinary_plateau_max_m),
    )
    ocean_metrics = _ocean_metrics(grid, frame, fields, missing, area, rel, ocean)

    flags = {
        "ordinary_plateau_like": bool(
            land_metrics["inland_area_fraction_of_world"] >= 0.02
            and land_metrics["ordinary_plateau_fraction"] >= 0.30
            and (
                land_metrics["broad_flat_inland_component_count"] > 0
                or land_metrics["inland_detail_entropy"] < 0.35
            )
        ),
        "lowland_plain_deficient": bool(
            land_metrics["continental_land_area_fraction_of_world"] >= 0.12
            and (
                land_metrics["lowland_plain_fraction"] < 0.06
                or land_metrics["broad_lowland_plain_component_count"] <= 0
                or land_metrics["largest_lowland_plain_component_fraction"] < 0.025
                or (
                    land_metrics["lowland_plain_fraction"] >= 0.06
                    and land_metrics["lowland_plain_parented_fraction"] < 0.55
                )
            )
        ),
        "ordinary_deep_ocean_like": bool(
            ocean_metrics["ocean_area_fraction_of_world"] >= 0.25
            and ocean_metrics["ocean_fabric_entropy"] < 0.45
            and ocean_metrics["ridge_visible_fraction"] < 0.02
            and ocean_metrics["fracture_zone_length_fraction"] < 0.015
            and ocean_metrics["hotspot_track_count"] == 0
            and ocean_metrics["seamount_chain_count"] == 0
        ),
    }

    return {
        "time_myr": time_myr,
        "sea_level_m": sea_level,
        "usable": True,
        "missing_fields": sorted(set(missing)),
        "land_fraction": float(area[land].sum() / total_area),
        "ocean_fraction": float(area[ocean].sum() / total_area),
        "land_metrics": land_metrics,
        "ocean_metrics": ocean_metrics,
        "diagnostic_flags": flags,
    }


def _land_metrics(
    grid: Any,
    fields: dict[str, Any],
    globals_: dict[str, Any],
    missing: list[str],
    area: np.ndarray,
    rel: np.ndarray,
    land: np.ndarray,
    time_myr: float,
    *,
    inland_width_steps: float,
    ordinary_plateau_min_m: float,
    ordinary_plateau_max_m: float,
) -> dict[str, Any]:
    n = int(getattr(grid, "n", 0))
    total_area = max(float(area.sum()), 1.0e-12)
    crust_type, _ = _field(fields, "crust.type", n, missing, default=0.0)
    crust_domain, _ = _field(fields, "crust.domain", n, missing, default=0.0)
    detail, has_detail = _field(fields, "terrain.continental_detail", n, missing, default=0.0)
    detail_region, has_detail_region = _field(
        fields, "terrain.continental_detail_region_code", n, missing, default=0.0)
    inland_region, has_inland_region = _field(
        fields, "terrain.inland_geomorphology_region_code", n, missing, default=0.0)
    province_code, has_province_code = _field(
        fields, "terrain.continental_province_code", n, missing, default=0.0)
    orogeny_age, has_orogeny_age = _field(
        fields, "tectonics.orogeny_age_myr", n, missing, default=-1.0)
    rift_stage, has_rift_stage = _field(
        fields, "terrain.rift_margin_stage", n, missing, default=0.0)
    sediment, has_sediment = _field(
        fields, "sediment.thickness_m", n, missing, default=0.0)
    lowland_memory = np.nan_to_num(
        np.asarray(
            fields.get(
                "terrain.p174_lowland_plain_continuity_memory",
                np.zeros(n, dtype=np.float64),
            ),
            dtype=np.float64,
        ),
        nan=0.0,
        posinf=1.0,
        neginf=0.0,
    )
    if lowland_memory.shape != (n,):
        lowland_memory = np.zeros(n, dtype=np.float64)

    continental_land = land & (crust_type >= 0.5)
    width = _mask_width_steps(grid, continental_land)
    inland = continental_land & (width >= float(inland_width_steps))
    if not inland.any():
        inland = continental_land.copy()
    land_area = float(area[land].sum())
    cont_area = float(area[continental_land].sum())
    inland_area = float(area[inland].sum())
    if inland_area <= 0.0:
        out = _zero_metrics(LAND_METRIC_KEYS)
        out.update({
            "land_area_fraction_of_world": float(land_area / total_area),
            "continental_land_area_fraction_of_world": float(cont_area / total_area),
            "inland_area_fraction_of_world": 0.0,
            "inland_area_fraction_of_land": 0.0,
            "inland_cell_count": 0,
        })
        out.update(_land_global_telemetry_metrics(globals_))
        return out

    category = _land_category(
        detail,
        detail_region,
        inland_region,
        province_code,
        crust_domain,
        has_detail=has_detail,
        has_detail_region=has_detail_region,
        has_inland_region=has_inland_region,
        has_province_code=has_province_code,
    )
    local_relief = _neighbor_range(grid, rel)
    continental_local_relief = _neighbor_range_within(grid, rel, continental_land)
    ordinary_semantic = (
        (category <= 0)
        | np.isin(detail.astype(int), [int(CONT_DETAIL_PLATFORM), int(CONT_DETAIL_PLATEAU)])
        | np.isin(
            inland_region.astype(int),
            [
                int(INLAND_PROVINCE_PLATFORM),
                int(INLAND_PROVINCE_PLATFORM_SWELL),
            ],
        )
    )
    ordinary_plateau = (
        inland
        & (rel >= float(ordinary_plateau_min_m))
        & (rel <= float(ordinary_plateau_max_m))
        & (local_relief <= 360.0)
        & ordinary_semantic
    )
    broad_flat_components = 0
    broad_area_floor = max(0.015 * inland_area, 1.0e-12)
    for comp in _component_cell_sets(grid, ordinary_plateau):
        comp_area = float(area[comp].sum())
        if comp_area < broad_area_floor and comp.size < 8:
            continue
        comp_rel = rel[comp]
        if comp_rel.size and float(np.percentile(comp_rel, 90) - np.percentile(comp_rel, 10)) <= 520.0:
            broad_flat_components += 1

    province_counts: list[int] = []
    for comp in _component_cell_sets(grid, continental_land):
        comp_area = float(area[comp].sum())
        if comp_area <= 0.0 or comp_area / max(cont_area, 1.0e-12) < 0.03:
            continue
        comp_categories = category[comp].astype(int)
        count = 0
        for code in np.unique(comp_categories):
            mask = comp_categories == int(code)
            if float(area[comp[mask]].sum()) / comp_area >= 0.02:
                count += 1
        province_counts.append(count)

    old_orogen = inland & (
        (inland_region.astype(int) == int(INLAND_PROVINCE_OLD_OROGEN))
        | (
            (detail.astype(int) == int(CONT_DETAIL_OROGEN))
            & (
                ((time_myr - orogeny_age) >= 250.0)
                if has_orogeny_age else np.ones(n, dtype=bool)
            )
        )
    )
    rift_basin = inland & (
        (detail.astype(int) == int(CONT_DETAIL_RIFT_BASIN))
        | (inland_region.astype(int) == int(INLAND_PROVINCE_RIFT))
        | (
            np.isin(
                rift_stage.astype(int),
                [
                    int(RIFT_MARGIN_STAGE_RIFT_BASIN),
                ],
            )
            if has_rift_stage else np.zeros(n, dtype=bool)
        )
    )
    shield = inland & (
        (detail.astype(int) == int(CONT_DETAIL_SHIELD))
        | (inland_region.astype(int) == int(INLAND_PROVINCE_SHIELD))
        | (crust_domain.astype(int) == int(DOMAIN_CRATON))
    )
    platform = inland & (
        (detail.astype(int) == int(CONT_DETAIL_PLATFORM))
        | (inland_region.astype(int) == int(INLAND_PROVINCE_PLATFORM))
    )
    basin = inland & (
        (detail.astype(int) == int(CONT_DETAIL_BASIN))
        | (inland_region.astype(int) == int(INLAND_PROVINCE_SAG_BASIN))
    )
    shield_share = float(area[shield].sum() / inland_area)
    platform_share = float(area[platform].sum() / inland_area)
    basin_share = float(area[basin].sum() / inland_area)
    split_union = shield | platform | basin
    split_fraction = (
        float(area[split_union].sum() / inland_area)
        if sum(v >= 0.02 for v in (shield_share, platform_share, basin_share)) >= 2
        else 0.0
    )
    lowland_domain = continental_land & (width >= 2.0)
    if not lowland_domain.any():
        lowland_domain = continental_land.copy()
    detail_i = detail.astype(int)
    inland_i = inland_region.astype(int)
    province_i = province_code.astype(int)
    rift_i = rift_stage.astype(int)
    sediment = np.nan_to_num(
        np.asarray(sediment, dtype=np.float64),
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )
    if has_sediment and continental_land.any():
        sed_p35 = float(np.percentile(sediment[continental_land], 35))
        sed_p80 = max(
            sed_p35 + 1.0,
            float(np.percentile(sediment[continental_land], 80)),
        )
        sediment_signal = np.clip(
            (sediment - sed_p35) / max(sed_p80 - sed_p35, 1.0),
            0.0,
            1.0,
        )
        sediment_cut = max(
            80.0,
            float(np.percentile(sediment[continental_land], 65)),
        )
    else:
        sediment_signal = np.zeros(n, dtype=np.float64)
        sediment_cut = np.inf
    passive_lowland = lowland_domain & (
        (province_i == int(CONT_PROVINCE_PASSIVE_MARGIN_LOWLAND))
        | (
            np.isin(rift_i, [int(RIFT_MARGIN_STAGE_PASSIVE_LOWLAND)])
            if has_rift_stage else np.zeros(n, dtype=bool)
        )
    )
    platform_parent = lowland_domain & (
        np.isin(detail_i, [int(CONT_DETAIL_PLATFORM), int(CONT_DETAIL_BASIN)])
        | np.isin(inland_i, [
            int(INLAND_PROVINCE_PLATFORM),
            int(INLAND_PROVINCE_SAG_BASIN),
        ])
        | np.isin(province_i, [
            int(CONT_PROVINCE_PLATFORM),
            int(CONT_PROVINCE_INTRACRATONIC_BASIN),
            int(CONT_PROVINCE_FORELAND_BASIN),
        ])
    )
    old_eroded_margin_parent = lowland_domain & (
        (
            (province_i == int(CONT_PROVINCE_OLD_OROGEN))
            | (inland_i == int(INLAND_PROVINCE_OLD_OROGEN))
        )
        & (continental_local_relief <= 260.0)
        & (rel <= 800.0)
    )
    sediment_parent = lowland_domain & has_sediment & (sediment >= sediment_cut)
    p174_memory_parent = lowland_domain & (lowland_memory >= 0.25)
    sedimentary_lowland_context = (
        lowland_domain
        & (rel >= 0.0)
        & (rel <= 1900.0)
        & (
            platform_parent
            | passive_lowland
            | sediment_parent
            | (sediment_signal >= 0.14)
            | p174_memory_parent
            | np.isin(inland_i, [
                int(INLAND_PROVINCE_SAG_BASIN),
                int(INLAND_PROVINCE_RIFT),
            ])
        )
        & (province_i != int(CONT_PROVINCE_VOLCANIC_LIP_PLATEAU))
        & (crust_domain.astype(int) != int(DOMAIN_LIP))
    )
    accreted_lowland_context = (
        (crust_domain.astype(int) == int(DOMAIN_ACCRETED_TERRANE))
        & sedimentary_lowland_context
    )
    orogenic_lowland_context = (
        np.isin(detail_i, [
            int(CONT_DETAIL_OROGEN),
            int(CONT_DETAIL_PLATEAU),
        ])
        & sedimentary_lowland_context
        & (rel <= 900.0)
    )
    lowland_parent = (
        platform_parent
        | passive_lowland
        | old_eroded_margin_parent
        | sediment_parent
        | p174_memory_parent
        | accreted_lowland_context
        | orogenic_lowland_context
    )
    lowland_area_exclusion_domain = (
        (crust_domain.astype(int) == int(DOMAIN_LIP))
        | (
            (crust_domain.astype(int) == int(DOMAIN_ACCRETED_TERRANE))
            & ~accreted_lowland_context
        )
        | (province_i == int(CONT_PROVINCE_VOLCANIC_LIP_PLATEAU))
    )
    active_orogen_province = province_i == int(CONT_PROVINCE_ACTIVE_OROGEN)
    active_orogen_lowland_context = (
        active_orogen_province
        & (
            ~np.isin(detail_i, [
                int(CONT_DETAIL_OROGEN),
                int(CONT_DETAIL_PLATEAU),
            ])
            | orogenic_lowland_context
        )
        & (
            np.isin(detail_i, [
                int(CONT_DETAIL_BASIN),
                int(CONT_DETAIL_RIFT_BASIN),
            ])
            | np.isin(inland_i, [
                int(INLAND_PROVINCE_SAG_BASIN),
                int(INLAND_PROVINCE_RIFT),
            ])
            | (sediment_signal >= 0.14)
            | (sediment >= sediment_cut)
            | (lowland_memory >= 0.25)
        )
        & (rel <= 1900.0)
        & (crust_domain.astype(int) != int(DOMAIN_LIP))
        & (province_i != int(CONT_PROVINCE_VOLCANIC_LIP_PLATEAU))
    )
    lowland_active_exclusion_domain = (
        (
            np.isin(detail_i, [
                int(CONT_DETAIL_OROGEN),
                int(CONT_DETAIL_PLATEAU),
            ])
            & ~orogenic_lowland_context
        )
        | (active_orogen_province & ~active_orogen_lowland_context)
        | (province_i == int(CONT_PROVINCE_VOLCANIC_LIP_PLATEAU))
        | lowland_area_exclusion_domain
    )
    lowland_elevation_domain = (
        lowland_domain
        & (rel >= 0.0)
        & (rel <= 700.0)
        & ~lowland_area_exclusion_domain
    )
    lowland_relief_buffer_domain = (
        lowland_active_exclusion_domain
        | active_orogen_lowland_context
    )
    lowland_relief_neighbor_domain = (
        lowland_elevation_domain
        & ~lowland_relief_buffer_domain
    )
    lowland_internal_relief = _neighbor_range_within(
        grid,
        rel,
        lowland_relief_neighbor_domain,
    )
    lowland_elevation_parented = (
        lowland_elevation_domain & lowland_parent
    )
    lowland_plain = lowland_elevation_domain & (lowland_internal_relief <= 260.0)
    lowland_local_relief_blocked = (
        lowland_elevation_parented
        & (continental_local_relief > 260.0)
    )
    broad_lowland_components = 0
    largest_lowland_component_area = 0.0
    largest_lowland_elevation_parented_area = 0.0
    broad_lowland_floor = max(0.012 * max(cont_area, 1.0e-12), 1.0e-12)
    near_floor_plain_components = 0
    near_floor_parented_components = 0
    near_floor_area = 0.006 * max(cont_area, 1.0e-12)
    for comp in _component_cell_sets(grid, lowland_elevation_parented):
        comp_area = float(area[comp].sum())
        largest_lowland_elevation_parented_area = max(
            largest_lowland_elevation_parented_area,
            comp_area,
        )
        if comp_area >= near_floor_area or comp.size >= 5:
            near_floor_parented_components += 1
    for comp in _component_cell_sets(grid, lowland_plain):
        comp_area = float(area[comp].sum())
        largest_lowland_component_area = max(largest_lowland_component_area, comp_area)
        if comp_area >= near_floor_area or comp.size >= 5:
            near_floor_plain_components += 1
        if comp_area < broad_lowland_floor and comp.size < 8:
            continue
        comp_measure = comp[~lowland_relief_buffer_domain[comp]]
        if comp_measure.size == 0:
            comp_measure = comp
        comp_rel = rel[comp_measure]
        if comp_rel.size and float(np.percentile(comp_rel, 90) - np.percentile(comp_rel, 10)) <= 420.0:
            broad_lowland_components += 1
    lowland_area = float(area[lowland_plain].sum())
    lowland_parented_area = float(area[lowland_plain & lowland_parent].sum())
    telemetry = _land_global_telemetry_metrics(globals_)
    lowland_plain_fraction = (
        float(lowland_area / cont_area) if cont_area > 0.0 else 0.0
    )
    lowland_elevation_parented_fraction = (
        float(area[lowland_elevation_parented].sum() / cont_area)
        if cont_area > 0.0 else 0.0
    )
    largest_lowland_elevation_parented_fraction = (
        float(largest_lowland_elevation_parented_area / cont_area)
        if cont_area > 0.0 else 0.0
    )
    lowland_local_relief_blocked_fraction = (
        float(area[lowland_local_relief_blocked].sum() / cont_area)
        if cont_area > 0.0 else 0.0
    )
    largest_lowland_plain_fraction = (
        float(largest_lowland_component_area / cont_area) if cont_area > 0.0 else 0.0
    )
    lowland_plain_parented_fraction = (
        float(lowland_parented_area / lowland_area) if lowland_area > 0.0 else 0.0
    )
    lowland_area_floor = 0.06
    lowland_component_floor = 0.025
    lowland_parented_floor = 0.55
    lowland_plain_area_gap = max(lowland_area_floor - lowland_plain_fraction, 0.0)
    lowland_plain_component_gap = max(
        lowland_component_floor - largest_lowland_plain_fraction,
        0.0,
    )
    lowland_plain_parentage_gap = (
        max(lowland_parented_floor - lowland_plain_parented_fraction, 0.0)
        if lowland_plain_fraction >= lowland_area_floor else 0.0
    )
    lowland_relief_boundary_gap = max(
        largest_lowland_elevation_parented_fraction - largest_lowland_plain_fraction,
        0.0,
    )
    support_largest = max(
        float(telemetry.get("p104f_pre_p174_lowland_support_largest_component_fraction", 0.0)),
        float(telemetry.get("p174_support_component_response_largest_component_fraction", 0.0)),
        float(telemetry.get("p174_lowland_plain_largest_component_fraction_after", 0.0)),
    )
    lowland_support_to_plain_gap = max(
        support_largest - largest_lowland_plain_fraction,
        0.0,
    )
    lowland_upstream_support_gap = max(
        lowland_component_floor - max(
            support_largest,
            largest_lowland_elevation_parented_fraction,
        ),
        0.0,
    )
    active_exclusion = (
        lowland_domain
        & (rel >= 0.0)
        & (rel <= 700.0)
        & lowland_active_exclusion_domain
    )
    lowland_active_exclusion_fraction = (
        float(area[active_exclusion].sum() / cont_area) if cont_area > 0.0 else 0.0
    )
    lowland_residual_area_limited = float(
        lowland_plain_fraction < lowland_area_floor
    )
    lowland_residual_parentage_limited = float(
        lowland_plain_fraction >= lowland_area_floor
        and lowland_plain_parented_fraction < lowland_parented_floor
    )
    lowland_residual_relief_boundary_limited = float(
        largest_lowland_plain_fraction < lowland_component_floor
        and largest_lowland_elevation_parented_fraction >= lowland_component_floor
        and lowland_relief_boundary_gap >= 0.003
        and lowland_local_relief_blocked_fraction >= 0.010
    )
    lowland_residual_support_to_plain_limited = float(
        largest_lowland_plain_fraction < lowland_component_floor
        and support_largest >= lowland_component_floor
        and largest_lowland_elevation_parented_fraction < lowland_component_floor
        and lowland_support_to_plain_gap >= 0.003
    )
    lowland_residual_upstream_support_limited = float(
        largest_lowland_plain_fraction < lowland_component_floor
        and lowland_upstream_support_gap > 0.0
        and lowland_plain_fraction >= lowland_area_floor
        and lowland_plain_parented_fraction >= lowland_parented_floor
    )
    lowland_residual_active_exclusion_limited = float(
        largest_lowland_plain_fraction < lowland_component_floor
        and lowland_active_exclusion_fraction >= 0.010
        and largest_lowland_elevation_parented_fraction < lowland_component_floor
    )
    lowland_residual_component_segmentation_limited = float(
        largest_lowland_plain_fraction < lowland_component_floor
        and lowland_plain_fraction >= lowland_area_floor
        and lowland_plain_parented_fraction >= lowland_parented_floor
        and lowland_residual_relief_boundary_limited == 0.0
        and lowland_residual_support_to_plain_limited == 0.0
        and lowland_residual_upstream_support_limited == 0.0
        and lowland_residual_active_exclusion_limited == 0.0
        and (
            near_floor_plain_components >= 2
            or near_floor_parented_components >= 2
        )
    )
    if lowland_residual_area_limited:
        lowland_residual_dominant_code = 1.0
    elif lowland_residual_parentage_limited:
        lowland_residual_dominant_code = 2.0
    elif lowland_residual_active_exclusion_limited:
        lowland_residual_dominant_code = 3.0
    elif lowland_residual_relief_boundary_limited:
        lowland_residual_dominant_code = 4.0
    elif lowland_residual_support_to_plain_limited:
        lowland_residual_dominant_code = 5.0
    elif lowland_residual_upstream_support_limited:
        lowland_residual_dominant_code = 6.0
    elif lowland_residual_component_segmentation_limited:
        lowland_residual_dominant_code = 7.0
    else:
        lowland_residual_dominant_code = 0.0

    out = {
        "land_area_fraction_of_world": float(land_area / total_area),
        "continental_land_area_fraction_of_world": float(cont_area / total_area),
        "inland_area_fraction_of_world": float(inland_area / total_area),
        "inland_area_fraction_of_land": float(inland_area / land_area) if land_area else 0.0,
        "inland_cell_count": int(np.count_nonzero(inland)),
        "inland_detail_entropy": _categorical_entropy(category, area, inland),
        "ordinary_plateau_fraction": float(area[ordinary_plateau].sum() / inland_area),
        "broad_flat_inland_component_count": int(broad_flat_components),
        "continent_province_count_p50": (
            float(np.median(np.asarray(province_counts, dtype=np.float64)))
            if province_counts else 0.0
        ),
        "old_orogen_expression_fraction": float(area[old_orogen].sum() / inland_area),
        "rift_basin_expression_fraction": float(area[rift_basin].sum() / inland_area),
        "craton_shield_platform_split_fraction": float(split_fraction),
        "lowland_plain_fraction": lowland_plain_fraction,
        "lowland_elevation_parented_fraction": lowland_elevation_parented_fraction,
        "largest_lowland_elevation_parented_component_fraction": (
            largest_lowland_elevation_parented_fraction
        ),
        "lowland_local_relief_blocked_fraction": lowland_local_relief_blocked_fraction,
        "broad_lowland_plain_component_count": int(broad_lowland_components),
        "largest_lowland_plain_component_fraction": largest_lowland_plain_fraction,
        "lowland_plain_parented_fraction": lowland_plain_parented_fraction,
        "lowland_plain_area_gap_fraction": float(lowland_plain_area_gap),
        "lowland_plain_component_gap_fraction": float(lowland_plain_component_gap),
        "lowland_plain_parentage_gap_fraction": float(lowland_plain_parentage_gap),
        "lowland_relief_boundary_gap_fraction": float(lowland_relief_boundary_gap),
        "lowland_support_to_plain_gap_fraction": float(lowland_support_to_plain_gap),
        "lowland_upstream_support_gap_fraction": float(lowland_upstream_support_gap),
        "lowland_active_exclusion_fraction": float(lowland_active_exclusion_fraction),
        "lowland_near_floor_plain_component_count": int(near_floor_plain_components),
        "lowland_near_floor_parented_component_count": int(near_floor_parented_components),
        "lowland_residual_area_limited": float(lowland_residual_area_limited),
        "lowland_residual_parentage_limited": float(lowland_residual_parentage_limited),
        "lowland_residual_relief_boundary_limited": float(
            lowland_residual_relief_boundary_limited
        ),
        "lowland_residual_support_to_plain_limited": float(
            lowland_residual_support_to_plain_limited
        ),
        "lowland_residual_upstream_support_limited": float(
            lowland_residual_upstream_support_limited
        ),
        "lowland_residual_active_exclusion_limited": float(
            lowland_residual_active_exclusion_limited
        ),
        "lowland_residual_component_segmentation_limited": float(
            lowland_residual_component_segmentation_limited
        ),
        "lowland_residual_dominant_code": float(lowland_residual_dominant_code),
    }
    out.update(telemetry)
    return out


def _land_global_telemetry_metrics(globals_: dict[str, Any]) -> dict[str, float]:
    """Expose terrain response telemetry as per-frame P170 land metrics."""
    out: dict[str, float] = {}
    for metric_key, global_key in LAND_GLOBAL_TELEMETRY_METRICS.items():
        try:
            out[metric_key] = float(globals_.get(global_key, 0.0))
        except (TypeError, ValueError):
            out[metric_key] = 0.0
    return out


def _ocean_metrics(
    grid: Any,
    frame: Any,
    fields: dict[str, Any],
    missing: list[str],
    area: np.ndarray,
    rel: np.ndarray,
    ocean: np.ndarray,
) -> dict[str, Any]:
    n = int(getattr(grid, "n", 0))
    total_area = max(float(area.sum()), 1.0e-12)
    ocean_area = float(area[ocean].sum())
    depth = np.maximum(-rel, 0.0)
    depth_province, _ = _field(fields, "ocean.depth_province", n, missing, default=0.0)
    crust_age, has_age = _field(fields, "crust.age_myr", n, missing, default=np.nan)
    crust_type, _ = _field(fields, "crust.type", n, missing, default=0.0)
    crust_domain, _ = _field(fields, "crust.domain", n, missing, default=0.0)
    crust_origin, _ = _field(fields, "crust.origin", n, missing, default=0.0)
    if ocean_area <= 0.0:
        out = _zero_metrics(OCEAN_METRIC_KEYS)
        out["ocean_area_fraction_of_world"] = 0.0
        return out

    coast_distance = _distance_from_sources(grid, land_mask=~ocean, domain=ocean, max_passes=8)
    ridge = ocean & (
        (depth_province.astype(int) == int(OCEAN_DEPTH_RIDGE))
        | (
            has_age
            & np.isfinite(crust_age)
            & (crust_age <= 15.0)
            & (depth <= 3300.0)
        )
    )
    trench = ocean & (depth_province.astype(int) == int(OCEAN_DEPTH_TRENCH))
    local_depth_relief = _neighbor_range(grid, depth)
    abyss = ocean & (
        (depth_province.astype(int) == int(OCEAN_DEPTH_ABYSS))
        | ((depth >= 3600.0) & (local_depth_relief <= 360.0) & ~trench)
    )

    age_contrast = _neighbor_range(grid, crust_age) if has_age else np.zeros(n, dtype=np.float64)
    if has_age and np.any(np.isfinite(age_contrast[ocean])):
        valid_contrast = age_contrast[ocean & np.isfinite(age_contrast)]
        contrast_cut = max(25.0, float(np.percentile(valid_contrast, 90))) if valid_contrast.size else 25.0
        fracture_zone = (
            ocean
            & ~ridge
            & ~trench
            & ((coast_distance < 0) | (coast_distance >= 4))
            & np.isfinite(crust_age)
            & (crust_age >= 18.0)
            & (age_contrast >= contrast_cut)
        )
    else:
        fracture_zone = np.zeros(n, dtype=bool)

    plateau = ocean & (
        (crust_domain.astype(int) == int(DOMAIN_LIP))
        | (crust_origin.astype(int) == int(ORIGIN_PLUME_IMPACT))
    )
    microcontinent = (
        ocean
        & (
            (crust_type >= 0.5)
            | (crust_domain.astype(int) == int(DOMAIN_ACCRETED_TERRANE))
        )
        & ((coast_distance < 0) | (coast_distance >= 3))
    )
    terrain_diag = {}
    diagnostics = getattr(frame, "diagnostics", {}) or {}
    if isinstance(diagnostics, dict):
        terrain_diag = diagnostics.get("terrain", {}) or {}
    object_ridge = _frame_object_mask(
        frame,
        "terrain.ocean_fabric",
        {"spreading_center", "ridge_segment"},
        n,
    ) & ocean
    object_fracture = _frame_object_mask(
        frame,
        "terrain.ocean_fabric",
        {"transform_fault", "fracture_zone"},
        n,
    ) & ocean
    object_abyss = _frame_object_mask(
        frame,
        "terrain.ocean_fabric",
        {"abyssal_plain"},
        n,
    ) & ocean
    object_plateau = _frame_object_mask(
        frame,
        "terrain.arc_plume_landforms",
        {"oceanic_plateau", "large_igneous_province"},
        n,
    ) & ocean
    object_microcontinent = _frame_object_mask(
        frame,
        "terrain.arc_plume_landforms",
        {"microcontinent", "accreted_terrane"},
        n,
    ) & ocean
    object_parented_shoal = _frame_object_mask(
        frame,
        "terrain.arc_plume_landforms",
        {"seamount_chain", "hotspot_track", "island_arc"},
        n,
    ) & ocean
    object_margin_parented_shoal = _frame_object_mask(
        frame,
        "terrain.margin_landforms",
        {
            "volcanic_arc",
            "forearc_accretionary_prism",
            "passive_margin_wedge",
            "delta_fan",
        },
        n,
    ) & ocean
    object_backarc_parented_shoal = _frame_object_mask(
        frame,
        "terrain.arc_plume_landforms",
        {"back_arc_basin"},
        n,
    ) & ocean
    object_rift_margin_parented_shoal = _frame_object_mask(
        frame,
        "terrain.rift_margin_sequences",
        {"rift_margin_sequence"},
        n,
    ) & ocean
    object_shoal_support = _dilate_mask(
        grid,
        (
            object_ridge
            | object_fracture
            | object_plateau
            | object_microcontinent
            | object_parented_shoal
            | object_margin_parented_shoal
            | object_backarc_parented_shoal
            | object_rift_margin_parented_shoal
        ),
        passes=1,
    ) & ocean
    ridge |= object_ridge
    fracture_zone |= object_fracture
    abyss |= object_abyss
    plateau |= object_plateau
    microcontinent |= object_microcontinent
    parented_shoal = (
        ridge
        | plateau
        | microcontinent
        | trench
        | fracture_zone
        | object_parented_shoal
        | object_margin_parented_shoal
        | object_backarc_parented_shoal
        | object_rift_margin_parented_shoal
        | object_shoal_support
    )
    far_ocean = ocean & ((coast_distance < 0) | (coast_distance >= 4))
    unparented_shoal = far_ocean & (depth < 1500.0) & ~parented_shoal

    ocean_counts = _frame_object_kind_counts(frame, "terrain.ocean_fabric")
    if not ocean_counts:
        ocean_counts = terrain_diag.get("ocean_fabric_kind_counts", {})
    arc_counts = _frame_object_kind_counts(frame, "terrain.arc_plume_landforms")
    if not arc_counts:
        arc_counts = terrain_diag.get("arc_plume_landform_kind_counts", {})
    ocean_counts = ocean_counts if isinstance(ocean_counts, dict) else {}
    arc_counts = arc_counts if isinstance(arc_counts, dict) else {}

    hotspot_track_count = int(arc_counts.get("hotspot_track", 0))
    seamount_chain_count = int(
        arc_counts.get("seamount_chain", 0) + arc_counts.get("hotspot_track", 0))

    out = {
        "ocean_area_fraction_of_world": float(ocean_area / total_area),
        "ocean_fabric_entropy": _categorical_entropy(depth_province.astype(int), area, ocean),
        "ridge_visible_fraction": float(area[ridge].sum() / ocean_area),
        "ridge_age_symmetry_score": _ridge_age_symmetry_score(grid, ocean, ridge, crust_age),
        "fracture_zone_length_fraction": float(area[fracture_zone].sum() / ocean_area),
        "abyssal_plain_fraction": float(area[abyss].sum() / ocean_area),
        "hotspot_track_count": hotspot_track_count,
        "seamount_chain_count": seamount_chain_count,
        "oceanic_plateau_fraction": float(area[plateau].sum() / ocean_area),
        "microcontinent_fraction": float(area[microcontinent].sum() / ocean_area),
        "unparented_shoal_fraction": float(area[unparented_shoal].sum() / ocean_area),
    }
    return out


def _field(
    fields: dict[str, Any],
    name: str,
    n: int,
    missing: list[str],
    *,
    default: float | int | None = None,
) -> tuple[np.ndarray, bool]:
    if name not in fields:
        missing.append(name)
        fill = 0.0 if default is None else default
        return np.full(n, fill, dtype=np.float64), False
    arr = np.asarray(fields[name], dtype=np.float64)
    if arr.shape != (n,):
        missing.append(name)
        fill = 0.0 if default is None else default
        return np.full(n, fill, dtype=np.float64), False
    return arr, True


def _land_category(
    detail: np.ndarray,
    detail_region: np.ndarray,
    inland_region: np.ndarray,
    province_code: np.ndarray,
    crust_domain: np.ndarray,
    *,
    has_detail: bool,
    has_detail_region: bool,
    has_inland_region: bool,
    has_province_code: bool,
) -> np.ndarray:
    category = np.zeros(detail.shape, dtype=np.int64)
    if has_detail:
        category = detail.astype(int)
    if has_detail_region:
        mask = detail_region > 0
        category[mask] = 100 + detail_region[mask].astype(int)
    if has_inland_region:
        mask = inland_region > 0
        category[mask] = 200 + inland_region[mask].astype(int)
    if has_province_code:
        mask = (category == 0) & (province_code > 0)
        category[mask] = 300 + province_code[mask].astype(int)
    mask = category == 0
    category[mask] = 400 + crust_domain[mask].astype(int)
    return category


def _categorical_entropy(category: np.ndarray, area: np.ndarray, mask: np.ndarray) -> float:
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return 0.0
    total = max(float(area[mask].sum()), 1.0e-12)
    values = np.asarray(category[mask]).astype(int)
    shares = []
    for value in np.unique(values):
        share = float(area[mask & (category.astype(int) == int(value))].sum() / total)
        if share > 0.0:
            shares.append(share)
    if len(shares) <= 1:
        return 0.0
    entropy = -sum(p * float(np.log(p)) for p in shares)
    return float(entropy / max(float(np.log(len(shares))), 1.0e-12))


def _ridge_age_symmetry_score(
    grid: Any,
    ocean: np.ndarray,
    ridge: np.ndarray,
    age: np.ndarray,
) -> float:
    ridge_cells = np.where(ridge & np.isfinite(age))[0]
    if ridge_cells.size == 0:
        return 0.0
    scores = []
    for c in ridge_cells:
        nbs = np.asarray(grid.neighbors[int(c)], dtype=int)
        nbs = nbs[ocean[nbs] & ~ridge[nbs] & np.isfinite(age[nbs])]
        if nbs.size == 0:
            continue
        older = float(np.count_nonzero(age[nbs] >= age[int(c)] + 3.0))
        scores.append(min(older / 2.0, 1.0))
    return float(np.mean(scores)) if scores else 0.0


def _frame_object_kind_counts(frame: Any, key: str) -> dict[str, int]:
    objects_by_key = getattr(frame, "objects", {}) or {}
    objects = objects_by_key.get(key, []) if isinstance(objects_by_key, dict) else []
    if not isinstance(objects, list):
        return {}
    counts: dict[str, int] = {}
    for obj in objects:
        if not isinstance(obj, dict):
            continue
        kind = str(obj.get("kind", ""))
        if not kind:
            continue
        counts[kind] = counts.get(kind, 0) + 1
    return dict(sorted(counts.items()))


def _frame_object_mask(
    frame: Any,
    key: str,
    kinds: set[str],
    n: int,
) -> np.ndarray:
    objects_by_key = getattr(frame, "objects", {}) or {}
    objects = objects_by_key.get(key, []) if isinstance(objects_by_key, dict) else []
    mask = np.zeros(int(n), dtype=bool)
    if not isinstance(objects, list):
        return mask
    for obj in objects:
        if not isinstance(obj, dict):
            continue
        if str(obj.get("kind", "")) not in kinds:
            continue
        cells = obj.get("cells", [])
        try:
            arr = np.asarray(cells, dtype=np.int64)
        except (TypeError, ValueError):
            arr = np.asarray([], dtype=np.int64)
        arr = arr[(arr >= 0) & (arr < int(n))]
        if arr.size:
            mask[arr] = True
            continue
        cell = obj.get("cell")
        if cell is not None:
            try:
                c = int(cell)
            except (TypeError, ValueError):
                continue
            if 0 <= c < int(n):
                mask[c] = True
    return mask


def _dilate_mask(grid: Any, mask: np.ndarray, passes: int = 1) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool).copy()
    if not mask.any() or passes <= 0:
        return mask
    out = mask.copy()
    frontier = mask.copy()
    for _ in range(int(passes)):
        nxt = out.copy()
        for c in np.where(frontier)[0]:
            nbs = np.asarray(grid.neighbors[int(c)], dtype=int)
            nxt[nbs] = True
        frontier = nxt & ~out
        out = nxt
        if not frontier.any():
            break
    return out


def _neighbor_range(grid: Any, values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    out = np.zeros(int(getattr(grid, "n", values.size)), dtype=np.float64)
    for c in range(out.size):
        nbs = np.asarray(grid.neighbors[int(c)], dtype=int)
        samples = values[np.concatenate([np.asarray([c], dtype=int), nbs])]
        samples = samples[np.isfinite(samples)]
        out[c] = float(samples.max() - samples.min()) if samples.size else 0.0
    return out


def _neighbor_range_within(
    grid: Any,
    values: np.ndarray,
    domain: np.ndarray,
) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    domain = np.asarray(domain, dtype=bool)
    out = np.zeros(int(getattr(grid, "n", values.size)), dtype=np.float64)
    for c in range(out.size):
        if not domain[int(c)]:
            continue
        nbs = np.asarray(grid.neighbors[int(c)], dtype=int)
        nbs = nbs[domain[nbs]]
        samples = values[np.concatenate([np.asarray([c], dtype=int), nbs])]
        samples = samples[np.isfinite(samples)]
        out[c] = float(samples.max() - samples.min()) if samples.size else 0.0
    return out


def _mask_width_steps(grid: Any, mask: np.ndarray) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    width = np.zeros(int(getattr(grid, "n", mask.size)), dtype=np.float64)
    if not mask.any():
        return width
    boundary = np.zeros(width.size, dtype=bool)
    for c in np.where(mask)[0]:
        if (~mask[np.asarray(grid.neighbors[int(c)], dtype=int)]).any():
            boundary[int(c)] = True
    if not boundary.any():
        width[mask] = np.sqrt(float(mask.sum()))
        return width
    queue = [int(c) for c in np.where(boundary)[0]]
    seen = np.zeros(width.size, dtype=bool)
    seen[queue] = True
    width[queue] = 1.0
    head = 0
    while head < len(queue):
        c = queue[head]
        head += 1
        for nb in np.asarray(grid.neighbors[c], dtype=int):
            nb = int(nb)
            if not mask[nb] or seen[nb]:
                continue
            width[nb] = width[c] + 1.0
            seen[nb] = True
            queue.append(nb)
    return width


def _distance_from_sources(
    grid: Any,
    *,
    land_mask: np.ndarray,
    domain: np.ndarray,
    max_passes: int,
) -> np.ndarray:
    domain = np.asarray(domain, dtype=bool)
    frontier = np.asarray(land_mask, dtype=bool).copy()
    seen = frontier.copy()
    dist = np.full(int(getattr(grid, "n", domain.size)), -1, dtype=np.int16)
    for p in range(1, int(max_passes) + 1):
        nxt = np.zeros(dist.size, dtype=bool)
        for c in np.where(frontier)[0]:
            nxt[np.asarray(grid.neighbors[int(c)], dtype=int)] = True
        nxt &= domain & ~seen
        if not nxt.any():
            break
        dist[nxt] = int(p)
        seen |= nxt
        frontier = nxt
    return dist


def _component_cell_sets(grid: Any, mask: np.ndarray) -> list[np.ndarray]:
    mask = np.asarray(mask, dtype=bool)
    seen = np.zeros(int(getattr(grid, "n", mask.size)), dtype=bool)
    comps: list[np.ndarray] = []
    for start in np.where(mask)[0]:
        if seen[int(start)]:
            continue
        stack = [int(start)]
        seen[int(start)] = True
        cells: list[int] = []
        while stack:
            c = stack.pop()
            cells.append(c)
            for nb in np.asarray(grid.neighbors[c], dtype=int):
                nb = int(nb)
                if not mask[nb] or seen[nb]:
                    continue
                seen[nb] = True
                stack.append(nb)
        comps.append(np.asarray(cells, dtype=np.int64))
    return comps


def _metric_extremes(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for group_name in ("land_metrics", "ocean_metrics"):
        keys = LAND_METRIC_KEYS if group_name == "land_metrics" else OCEAN_METRIC_KEYS
        out[group_name] = {}
        for key in keys:
            values = np.asarray([
                float(row.get(group_name, {}).get(key, 0.0)) for row in rows
            ], dtype=np.float64)
            out[group_name][key] = {
                "min": float(np.min(values)) if values.size else 0.0,
                "max": float(np.max(values)) if values.size else 0.0,
                "median": float(np.median(values)) if values.size else 0.0,
            }
    return out


def _flag_windows(rows: list[dict[str, Any]], flag: str) -> list[dict[str, float]]:
    windows: list[dict[str, float]] = []
    start: float | None = None
    end: float | None = None
    for row in rows:
        t = float(row["time_myr"])
        active = bool(row.get("diagnostic_flags", {}).get(flag, False))
        if active and start is None:
            start = t
            end = t
        elif active:
            end = t
        elif start is not None and end is not None:
            windows.append({
                "start_myr": float(start),
                "end_myr": float(end),
                "duration_myr": float(max(end - start, 0.0)),
            })
            start = None
            end = None
    if start is not None and end is not None:
        windows.append({
            "start_myr": float(start),
            "end_myr": float(end),
            "duration_myr": float(max(end - start, 0.0)),
        })
    return windows


def _terminal_quality_jump(rows: list[dict[str, Any]], group: str, metric: str) -> float:
    if len(rows) < 3:
        return 0.0
    preterminal = np.asarray([
        float(row.get(group, {}).get(metric, 0.0)) for row in rows[:-1]
    ], dtype=np.float64)
    terminal = float(rows[-1].get(group, {}).get(metric, 0.0))
    return float(terminal - np.median(preterminal)) if preterminal.size else 0.0


def _p174_continuity_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Summarize whether P174 land/ocean detail appears continuously in time."""
    usable = [row for row in rows if bool(row.get("usable", False))]
    mature_rows = []
    for row in usable:
        land = row.get("land_metrics", {})
        cont_fraction = float(land.get("continental_land_area_fraction_of_world", 0.0))
        platform_split = float(land.get("craton_shield_platform_split_fraction", 0.0))
        rift_fraction = float(land.get("rift_basin_expression_fraction", 0.0))
        old_orogen_fraction = float(land.get("old_orogen_expression_fraction", 0.0))
        support_ready = (
            cont_fraction >= 0.12
            and (
                platform_split >= 0.65
                or rift_fraction >= 0.04
                or old_orogen_fraction >= 0.03
            )
        )
        if support_ready:
            mature_rows.append(row)

    lowland_values = np.asarray([
        float(row.get("land_metrics", {}).get("lowland_plain_fraction", 0.0))
        for row in mature_rows
    ], dtype=np.float64)
    ocean_values = np.asarray([
        float(row.get("ocean_metrics", {}).get("ocean_fabric_entropy", 0.0))
        for row in usable
    ], dtype=np.float64)
    mature_deficient = [
        row for row in mature_rows
        if bool(row.get("diagnostic_flags", {}).get("lowland_plain_deficient"))
    ]

    def step_extremes(values: np.ndarray) -> tuple[float, float]:
        if values.size < 2:
            return 0.0, 0.0
        delta = np.diff(values)
        return float(np.max(delta)), float(np.min(delta))

    lowland_max_up, lowland_max_down = step_extremes(lowland_values)
    ocean_max_up, ocean_max_down = step_extremes(ocean_values)
    terminal_lowland_jump = (
        float(lowland_values[-1] - np.median(lowland_values[:-1]))
        if lowland_values.size >= 3 else 0.0
    )
    terminal_ocean_jump = (
        float(ocean_values[-1] - np.median(ocean_values[:-1]))
        if ocean_values.size >= 3 else 0.0
    )
    mature_count = int(len(mature_rows))
    deficient_count = int(len(mature_deficient))
    deficient_fraction = (
        float(deficient_count / mature_count) if mature_count else 0.0
    )
    terminal_lowland_pop = bool(
        lowland_values.size >= 4
        and terminal_lowland_jump > 0.035
        and deficient_fraction > 0.25
        and not bool(mature_rows[-1].get("diagnostic_flags", {}).get(
            "lowland_plain_deficient", False))
    )

    def deficient_metric_count(metric: str) -> int:
        return int(sum(
            float(row.get("land_metrics", {}).get(metric, 0.0)) >= 0.5
            for row in mature_deficient
        ))

    return {
        "mature_support_frame_count": mature_count,
        "mature_lowland_deficient_frame_count": deficient_count,
        "mature_lowland_deficient_fraction": float(deficient_fraction),
        "mature_lowland_continuity_score": float(1.0 - deficient_fraction),
        "mature_lowland_residual_area_limited_frame_count": (
            deficient_metric_count("lowland_residual_area_limited")
        ),
        "mature_lowland_residual_parentage_limited_frame_count": (
            deficient_metric_count("lowland_residual_parentage_limited")
        ),
        "mature_lowland_residual_relief_boundary_frame_count": (
            deficient_metric_count("lowland_residual_relief_boundary_limited")
        ),
        "mature_lowland_residual_support_to_plain_frame_count": (
            deficient_metric_count("lowland_residual_support_to_plain_limited")
        ),
        "mature_lowland_residual_upstream_support_frame_count": (
            deficient_metric_count("lowland_residual_upstream_support_limited")
        ),
        "mature_lowland_residual_active_exclusion_frame_count": (
            deficient_metric_count("lowland_residual_active_exclusion_limited")
        ),
        "mature_lowland_residual_component_segmentation_frame_count": (
            deficient_metric_count("lowland_residual_component_segmentation_limited")
        ),
        "lowland_plain_max_positive_step": float(lowland_max_up),
        "lowland_plain_max_negative_step": float(lowland_max_down),
        "terminal_lowland_plain_jump": float(terminal_lowland_jump),
        "terminal_lowland_pop_candidate": bool(terminal_lowland_pop),
        "ocean_fabric_entropy_max_positive_step": float(ocean_max_up),
        "ocean_fabric_entropy_max_negative_step": float(ocean_max_down),
        "terminal_ocean_fabric_entropy_jump": float(terminal_ocean_jump),
    }


def _zero_metrics(keys: tuple[str, ...]) -> dict[str, float | int]:
    return {key: 0.0 for key in keys}


def _skipped_summary(reason: str, frame_count: int) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "skipped": True,
        "skip_reason": str(reason),
        "frame_count": int(frame_count),
        "usable_frame_count": 0,
        "required_land_metric_keys": list(LAND_METRIC_KEYS),
        "required_ocean_metric_keys": list(OCEAN_METRIC_KEYS),
        "required_metric_keys_present": False,
        "p174_continuity": _p174_continuity_summary([]),
        "acceptance": {
            "audit_completed": False,
            "required_metric_keys_present": False,
            "generation_behavior_changed": False,
            "identifies_land_time_windows": False,
            "identifies_ocean_time_windows": False,
        },
        "frame_rows": [],
    }
