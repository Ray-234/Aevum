"""Current generated-world morphology gap inventory diagnostics."""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
import re
from typing import Any

import numpy as np

from aevum.diagnostics.boundary_process_geometry import (
    boundary_process_geometry_summary,
)
from aevum.diagnostics.crust_sediment_province_coupling import (
    crust_sediment_province_coupling_summary,
)
from aevum.diagnostics.drainage_divide_province_alignment import (
    drainage_divide_province_alignment_summary,
)
from aevum.diagnostics.earth_reference import earth_reference_distribution_metrics
from aevum.diagnostics.generated_province_reference import (
    generated_province_reference_comparison_summary,
)
from aevum.diagnostics.morphology import compute_world_morphology
from aevum.diagnostics.mountain_inventory_expression import (
    mountain_inventory_expression_summary,
)
from aevum.diagnostics.old_orogen_erosion_decay import (
    old_orogen_erosion_decay_summary,
)
from aevum.diagnostics.plateau_area_cap_and_decay import (
    plateau_area_cap_and_decay_summary,
)
from aevum.diagnostics.province_diversity import (
    generated_world_province_diversity_summary,
)
from aevum.diagnostics.province_reference_graph import province_reference_graph_summary
from aevum.diagnostics.real_earth_hypsometry import (
    compare_metrics_to_fixture,
    generated_world_hypsometry_metrics,
    real_earth_hypsometry_fixture_summary,
)
from aevum.diagnostics.reference_source_ledger import source_ledger_summary
from aevum.diagnostics.rift_margin_escarpment_sequence import (
    rift_margin_escarpment_sequence_summary,
)
from aevum.diagnostics.source_to_sink_sediment_budget import (
    source_to_sink_sediment_budget_summary,
)
from aevum.diagnostics.wilson_cycle_lifecycle import (
    wilson_cycle_lifecycle_summary,
)
from aevum.modules.tectonics import CONT, DOMAIN_LIP
from aevum.modules.terrain import (
    CONT_DETAIL_ARC_MICROCONTINENT,
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATEAU,
    CONT_DETAIL_RIFT_BASIN,
)


SCHEMA = "aevum.current_world_morphology_gap_inventory.v1"

OWNER_LAYERS = (
    "planform",
    "province_graph",
    "boundary_lifecycle",
    "crust_sediment",
    "drainage_erosion",
    "landform_expression",
    "bathymetry_margin",
    "compiler_render",
)

GAP_CATEGORIES = (
    "missing_process",
    "missing_field",
    "missing_object",
    "wrong_amplitude",
    "wrong_area_scale",
    "wrong_adjacency",
    "wrong_lifecycle",
    "sediment_crust_coupling",
    "compiler_render_mismatch",
    "asset_review_pending",
)

REQUIRED_REVIEW_ASSETS = (
    "elevation.png",
    "terrain_provinces.png",
    "continental_detail_provinces.png",
    "ocean_depth_provinces.png",
    "crust_age.png",
    "history.png",
    "timeline.png",
    "hexmap.png",
)

FUTURE_STAGE_BY_OWNER = {
    "planform": "R1/P91 planform envelope and high-resolution review",
    "province_graph": "R2/P91 production province graph promotion",
    "boundary_lifecycle": "R3 successors to P81/P82",
    "crust_sediment": "R2/R5 production province-crust-sediment fields",
    "drainage_erosion": "R5 production drainage and source-to-sink routing",
    "landform_expression": "R4 production landform inventories and lifecycle fields",
    "bathymetry_margin": "R1/R4 margin bathymetry profile review",
    "compiler_render": "R7/P91 contact-sheet and compiler consistency audit",
}


def _canonical_digest(data: Any) -> str:
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _slug(text: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_")
    return slug or "unknown"


def _field_or_empty(world: Any, name: str, fill: float = 0.0) -> np.ndarray:
    if name in getattr(world, "fields", {}):
        return np.asarray(world.field(name), dtype=np.float64)
    return np.full(world.grid.n, fill, dtype=np.float64)


def _area_fraction(
    area: np.ndarray,
    mask: np.ndarray,
    within: np.ndarray | None = None,
) -> float:
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


def _object_mask(
    world: Any,
    object_set: str,
    kinds: set[str] | tuple[str, ...] | None = None,
) -> np.ndarray:
    mask = np.zeros(world.grid.n, dtype=bool)
    kind_set = None if kinds is None else {str(kind) for kind in kinds}
    for obj in getattr(world, "objects", {}).get(object_set, []):
        if kind_set is not None and str(obj.get("kind", "")) not in kind_set:
            continue
        cells = np.asarray(obj.get("cells", ()), dtype=np.int64)
        cells = cells[(cells >= 0) & (cells < world.grid.n)]
        if cells.size:
            mask[cells] = True
    return mask


def _one_ring_relief(grid: Any, rel: np.ndarray) -> np.ndarray:
    rel = np.asarray(rel, dtype=np.float64)
    relief = np.zeros(grid.n, dtype=np.float64)
    for cell in range(grid.n):
        cells = np.asarray([cell, *grid.neighbors[cell]], dtype=np.int64)
        vals = rel[cells]
        vals = vals[np.isfinite(vals)]
        if vals.size:
            relief[cell] = float(vals.max() - vals.min())
    return relief


def _add_gap(
    gaps: list[dict[str, Any]],
    *,
    gap_id: str,
    owner_layer: str,
    category: str,
    source_suite: str,
    evidence: str,
    metric_value: Any = None,
    future_stage: str | None = None,
    severity: str = "recorded",
) -> None:
    gaps.append({
        "gap_id": str(gap_id),
        "owner_layer": str(owner_layer),
        "category": str(category),
        "source_suite": str(source_suite),
        "evidence": str(evidence),
        "metric_value": metric_value,
        "future_stage": str(
            future_stage if future_stage is not None
            else FUTURE_STAGE_BY_OWNER.get(owner_layer, "P91 integrated review")
        ),
        "severity": str(severity),
    })


def _direct_surface_metrics(world: Any, morphology: Any) -> dict[str, Any]:
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    rel = _field_or_empty(world, "terrain.elevation_m") - float(world.sea_level)
    detail = _field_or_empty(world, "terrain.continental_detail").astype(int)
    crust = _field_or_empty(world, "crust.type").astype(int)
    domain = _field_or_empty(world, "crust.domain").astype(int)
    land = rel >= 0.0
    continental_land = land & (crust == int(CONT))
    width = np.asarray(
        morphology.fields["tectonics.exposed_continental_width_steps"],
        dtype=np.float64,
    )
    local_relief = _one_ring_relief(grid, rel)
    high_flat = (
        continental_land
        & (width >= 3.0)
        & (rel >= 1000.0)
        & (local_relief <= 450.0)
    )
    highland = continental_land & (rel >= 1600.0)
    parent_mask = (
        np.isin(detail, [
            int(CONT_DETAIL_OROGEN),
            int(CONT_DETAIL_PLATEAU),
            int(CONT_DETAIL_ARC_MICROCONTINENT),
        ])
        | _object_mask(world, "terrain.continental_landforms", {
            "orogen",
            "old_subdued_orogen",
            "plateau",
            "arc_microcontinent",
        })
        | _object_mask(world, "tectonics.lips")
        | (domain == DOMAIN_LIP)
    )
    highland_without_parent = highland & ~parent_mask
    basin_lowland = (
        continental_land
        & (
            (rel < 500.0)
            | np.isin(detail, [int(CONT_DETAIL_BASIN), int(CONT_DETAIL_RIFT_BASIN)])
            | _object_mask(world, "terrain.continental_landforms", {
                "interior_basin",
                "foreland_basin",
                "passive_margin_lowland",
                "rift_basin",
            })
            | (_field_or_empty(world, "terrain.passive_margin_lowland") > 0.0)
        )
    )
    component_labels = np.asarray(
        morphology.fields["tectonics.exposed_continental_component_id"],
        dtype=int,
    )
    major_component_metrics: list[dict[str, Any]] = []
    continental_area = max(float(area[continental_land].sum()), 1.0e-12)
    for component_id in sorted(
        int(value) for value in np.unique(component_labels[continental_land])
        if int(value) >= 0
    ):
        mask = continental_land & (component_labels == component_id)
        comp_area = float(area[mask].sum())
        if comp_area / continental_area < 0.08:
            continue
        major_component_metrics.append({
            "component_id": int(component_id),
            "area_fraction_of_continental_land": float(comp_area / continental_area),
            "lowland_fraction_lt500m": _area_fraction(area, mask & (rel < 500.0), mask),
            "lowland_fraction_lt1000m": _area_fraction(area, mask & (rel < 1000.0), mask),
            "high_flat_fraction": _area_fraction(area, high_flat, mask),
            "basin_lowland_fraction": _area_fraction(area, basin_lowland, mask),
        })

    return {
        "continental_land_fraction_world": _area_fraction(area, continental_land),
        "high_flat_interior_fraction_of_world": _area_fraction(area, high_flat),
        "high_flat_interior_fraction_of_continental_land": _area_fraction(
            area, high_flat, continental_land),
        "highland_fraction_of_continental_land": _area_fraction(
            area, highland, continental_land),
        "highland_without_parent_fraction_of_continental_land": _area_fraction(
            area, highland_without_parent, continental_land),
        "highland_without_parent_fraction_of_highlands": _area_fraction(
            area, highland_without_parent, highland),
        "basin_lowland_fraction_of_continental_land": _area_fraction(
            area, basin_lowland, continental_land),
        "lowland_fraction_lt500m_of_continental_land": _area_fraction(
            area, continental_land & (rel < 500.0), continental_land),
        "lowland_fraction_lt1000m_of_continental_land": _area_fraction(
            area, continental_land & (rel < 1000.0), continental_land),
        "local_relief_p50_m": float(np.percentile(local_relief[continental_land], 50))
        if continental_land.any() else 0.0,
        "local_relief_p90_m": float(np.percentile(local_relief[continental_land], 90))
        if continental_land.any() else 0.0,
        "major_component_metrics": major_component_metrics,
        "major_component_count": int(len(major_component_metrics)),
    }


def current_world_morphology_gap_inventory_summary(
    world: Any,
    *,
    tectonics: dict[str, Any] | None = None,
    compiler_metrics: dict[str, Any] | None = None,
    asset_paths: dict[str, str] | None = None,
) -> dict[str, Any]:
    morphology = compute_world_morphology(world)
    direct = _direct_surface_metrics(world, morphology)
    tectonics = tectonics or {"ocean_geography": {}}
    compiler_metrics = compiler_metrics or {
        "passed_envelope": False,
        "broad_land_to_water_fraction": 0.0,
        "broad_ocean_to_land_fraction": 0.0,
        "shelf_as_deep_ocean_fraction": 0.0,
        "lowland_as_mountain_fraction": 0.0,
        "terrain_elevation_sign_mismatch_fraction": 0.0,
        "invalid_source_array_shape": True,
    }
    asset_paths = asset_paths or {}

    fixture = real_earth_hypsometry_fixture_summary()
    source_ledger = source_ledger_summary()
    generated_hyp = generated_world_hypsometry_metrics(world, tectonics)
    hyp_comparison = compare_metrics_to_fixture(generated_hyp, fixture["metrics"])
    earth_reference = earth_reference_distribution_metrics(world, tectonics=tectonics)
    diversity = generated_world_province_diversity_summary(world)
    reference_graph = province_reference_graph_summary()
    province_comparison = generated_province_reference_comparison_summary(
        world,
        reference_graph=reference_graph,
        diversity=diversity,
    )
    boundary = boundary_process_geometry_summary(world)
    wilson = wilson_cycle_lifecycle_summary(world)
    crust_sediment = crust_sediment_province_coupling_summary(world)
    sediment_budget = source_to_sink_sediment_budget_summary(world)
    drainage = drainage_divide_province_alignment_summary(world)
    old_orogen = old_orogen_erosion_decay_summary(world)
    mountain = mountain_inventory_expression_summary(world)
    rift_margin = rift_margin_escarpment_sequence_summary(world)
    plateau = plateau_area_cap_and_decay_summary(world)

    gaps: list[dict[str, Any]] = []
    for key in earth_reference["out_of_envelope"]:
        _add_gap(
            gaps,
            gap_id=f"planform.{_slug(key)}_out_of_envelope",
            owner_layer="planform",
            category="wrong_area_scale",
            source_suite="P78",
            evidence=f"earth_reference out_of_envelope includes {key}",
            metric_value=earth_reference["metrics"].get(key, {}).get("value"),
        )
    for warning in morphology.metrics["warnings"]:
        _add_gap(
            gaps,
            gap_id=f"planform.{_slug(warning)}",
            owner_layer="planform",
            category="wrong_area_scale",
            source_suite="P78/P90",
            evidence=warning,
        )
    if direct["high_flat_interior_fraction_of_continental_land"] > 0.02:
        _add_gap(
            gaps,
            gap_id="crust_sediment.high_flat_interior_share",
            owner_layer="crust_sediment",
            category="wrong_amplitude",
            source_suite="P90",
            evidence="high interior cells have low local relief without sufficient process partitioning",
            metric_value=direct["high_flat_interior_fraction_of_continental_land"],
        )
    if direct["highland_without_parent_fraction_of_highlands"] > 0.15:
        _add_gap(
            gaps,
            gap_id="landform_expression.unparented_highland_share",
            owner_layer="landform_expression",
            category="missing_process",
            source_suite="P90/P72",
            evidence="highland cells lack orogen/old-orogen/plateau/arc/plume parentage",
            metric_value=direct["highland_without_parent_fraction_of_highlands"],
        )
    if direct["basin_lowland_fraction_of_continental_land"] < 0.20:
        _add_gap(
            gaps,
            gap_id="crust_sediment.basin_lowland_share_low",
            owner_layer="crust_sediment",
            category="sediment_crust_coupling",
            source_suite="P83/P90",
            evidence="continental basin/lowland share is below current calibration floor",
            metric_value=direct["basin_lowland_fraction_of_continental_land"],
        )

    for cls in province_comparison["missing_reference_classes"]:
        _add_gap(
            gaps,
            gap_id=f"province_graph.missing_reference_class.{cls}",
            owner_layer="province_graph",
            category="missing_process",
            source_suite="P80",
            evidence=f"missing reference class {cls}",
            metric_value=cls,
        )
    for edge in province_comparison["missing_required_class_edges"]:
        _add_gap(
            gaps,
            gap_id=f"province_graph.missing_required_edge.{_slug(edge)}",
            owner_layer="province_graph",
            category="wrong_adjacency",
            source_suite="P80",
            evidence=f"missing required class edge {edge}",
            metric_value=edge,
        )
    for process in boundary["current_generated"]["missing_process_types"]:
        _add_gap(
            gaps,
            gap_id=f"boundary_lifecycle.missing_boundary_process.{process}",
            owner_layer="boundary_lifecycle",
            category="missing_process",
            source_suite="P81",
            evidence=f"missing generated boundary process {process}",
            metric_value=process,
        )
    for object_set in wilson["current_generated"]["missing_object_sets"]:
        _add_gap(
            gaps,
            gap_id=f"boundary_lifecycle.missing_object_set.{_slug(object_set)}",
            owner_layer="boundary_lifecycle",
            category="missing_object",
            source_suite="P82",
            evidence=f"missing Wilson-cycle object set {object_set}",
            metric_value=object_set,
        )
    for field in crust_sediment["current_generated"]["missing_production_province_fields"]:
        _add_gap(
            gaps,
            gap_id=f"province_graph.missing_production_field.{_slug(field)}",
            owner_layer="province_graph",
            category="missing_field",
            source_suite="P83",
            evidence=f"missing production province field {field}",
            metric_value=field,
        )
    for obj in sediment_budget["current_generated"]["missing_source_to_sink_objects"]:
        _add_gap(
            gaps,
            gap_id=f"drainage_erosion.missing_source_to_sink.{_slug(obj)}",
            owner_layer="drainage_erosion",
            category="missing_object",
            source_suite="P84",
            evidence=f"missing source-to-sink object/field {obj}",
            metric_value=obj,
        )
    for item in drainage["current_generated"]["missing_drainage_items"]:
        _add_gap(
            gaps,
            gap_id=f"drainage_erosion.missing_drainage_item.{_slug(item)}",
            owner_layer="drainage_erosion",
            category="missing_object",
            source_suite="P85",
            evidence=f"missing drainage item {item}",
            metric_value=item,
        )
    for field in old_orogen["current_generated"]["missing_decay_fields"]:
        _add_gap(
            gaps,
            gap_id=f"drainage_erosion.missing_orogen_decay_field.{_slug(field)}",
            owner_layer="drainage_erosion",
            category="wrong_lifecycle",
            source_suite="P86",
            evidence=f"missing old-orogen decay field {field}",
            metric_value=field,
        )
    for field in mountain["current_generated"]["missing_inventory_fields"]:
        _add_gap(
            gaps,
            gap_id=f"landform_expression.missing_mountain_inventory_field.{_slug(field)}",
            owner_layer="landform_expression",
            category="missing_field",
            source_suite="P87",
            evidence=f"missing mountain inventory field {field}",
            metric_value=field,
        )
    for kind in mountain["current_generated"]["missing_expected_mountain_kinds"]:
        _add_gap(
            gaps,
            gap_id=f"landform_expression.missing_mountain_kind.{kind}",
            owner_layer="landform_expression",
            category="missing_process",
            source_suite="P87",
            evidence=f"missing expected mountain/plateau kind {kind}",
            metric_value=kind,
        )
    if mountain["current_generated"]["limitations"]["elongated_range_expression_underdeveloped"]:
        _add_gap(
            gaps,
            gap_id="landform_expression.elongated_range_expression_underdeveloped",
            owner_layer="landform_expression",
            category="wrong_area_scale",
            source_suite="P87",
            evidence="generated mountain objects lack elongated range expression",
            metric_value=mountain["current_generated"]["metrics"][
                "elongated_mountain_object_count"],
        )
    for item in rift_margin["current_generated"]["missing_sequence_items"]:
        owner = (
            "bathymetry_margin"
            if "rift_margin" in item or "escarpment" in item or "rift_shoulders" in item
            else "landform_expression"
        )
        _add_gap(
            gaps,
            gap_id=f"{owner}.missing_rift_margin_item.{_slug(item)}",
            owner_layer=owner,
            category="wrong_lifecycle",
            source_suite="P88",
            evidence=f"missing rift-margin sequence item {item}",
            metric_value=item,
        )
    if rift_margin["current_generated"]["limitations"]["passive_margin_lowland_objects_tiny"]:
        _add_gap(
            gaps,
            gap_id="bathymetry_margin.passive_margin_lowland_objects_tiny",
            owner_layer="bathymetry_margin",
            category="wrong_area_scale",
            source_suite="P88",
            evidence="passive-margin lowland objects remain too small",
            metric_value=rift_margin["current_generated"]["metrics"][
                "passive_margin_lowland_object_area_fraction_world"],
        )
    for item in plateau["current_generated"]["missing_plateau_items"]:
        _add_gap(
            gaps,
            gap_id=f"landform_expression.missing_plateau_item.{_slug(item)}",
            owner_layer="landform_expression",
            category="missing_field",
            source_suite="P89",
            evidence=f"missing plateau lifecycle/inventory item {item}",
            metric_value=item,
        )
    for kind in plateau["current_generated"]["missing_expected_plateau_kinds"]:
        _add_gap(
            gaps,
            gap_id=f"landform_expression.missing_plateau_kind.{kind}",
            owner_layer="landform_expression",
            category="missing_process",
            source_suite="P89",
            evidence=f"missing expected plateau kind {kind}",
            metric_value=kind,
        )
    if plateau["current_generated"]["limitations"][
        "high_interior_support_needs_p90_gap_audit"
    ]:
        _add_gap(
            gaps,
            gap_id="landform_expression.high_interior_without_plateau_support",
            owner_layer="landform_expression",
            category="wrong_amplitude",
            source_suite="P89/P90",
            evidence="high interior exists without plateau support audit fields",
            metric_value=plateau["current_generated"]["metrics"][
                "high_interior_without_plateau_fraction_of_continental_land"],
        )

    missing_assets = tuple(
        asset for asset in REQUIRED_REVIEW_ASSETS if asset not in asset_paths
    )
    for asset in missing_assets:
        _add_gap(
            gaps,
            gap_id=f"compiler_render.asset_review_pending.{_slug(asset)}",
            owner_layer="compiler_render",
            category="asset_review_pending",
            source_suite="P90/P91",
            evidence=f"required review asset not attached to P90 audit: {asset}",
            metric_value=asset,
        )
    if not bool(compiler_metrics.get("passed_envelope", False)):
        _add_gap(
            gaps,
            gap_id="compiler_render.compiler_consistency_failed",
            owner_layer="compiler_render",
            category="compiler_render_mismatch",
            source_suite="P90",
            evidence="compiled map contradicts source terrain/ocean provinces",
            metric_value={
                key: compiler_metrics.get(key)
                for key in (
                    "broad_land_to_water_fraction",
                    "broad_ocean_to_land_fraction",
                    "shelf_as_deep_ocean_fraction",
                    "lowland_as_mountain_fraction",
                    "terrain_elevation_sign_mismatch_fraction",
                )
            },
            severity="blocking",
        )

    owner_counts = Counter(gap["owner_layer"] for gap in gaps)
    category_counts = Counter(gap["category"] for gap in gaps)
    source_counts = Counter(gap["source_suite"] for gap in gaps)
    layer_summaries = {
        owner: {
            "gap_count": int(owner_counts.get(owner, 0)),
            "future_stage": FUTURE_STAGE_BY_OWNER[owner],
            "status": "gaps_recorded" if owner_counts.get(owner, 0) else "no_gap_recorded",
        }
        for owner in OWNER_LAYERS
    }
    unassigned = tuple(
        gap["gap_id"]
        for gap in gaps
        if gap["owner_layer"] not in OWNER_LAYERS
        or gap["category"] not in GAP_CATEGORIES
        or not gap["future_stage"]
    )
    generic_blockers = tuple(
        gap["gap_id"]
        for gap in gaps
        if "looks wrong" in gap["evidence"].lower()
        or "ugly" in gap["evidence"].lower()
    )
    current_residual_items = tuple(sorted({
        str(gap["metric_value"])
        for gap in gaps
        if isinstance(gap.get("metric_value"), str)
        and (
            "." in str(gap["metric_value"])
            or str(gap["metric_value"]) in {"transform", "plateau", "orogen"}
        )
    }))
    metrics = {
        "cells": int(world.grid.n),
        "world_time_myr": float(world.time_myr),
        "gap_count": int(len(gaps)),
        "owner_layer_count": int(len(owner_counts)),
        "category_count": int(len(category_counts)),
        "source_suite_count": int(len(source_counts)),
        "unassigned_gap_count": int(len(unassigned)),
        "generic_blocker_count": int(len(generic_blockers)),
        "current_residual_item_count": int(len(current_residual_items)),
        "missing_required_asset_count": int(len(missing_assets)),
        "required_review_asset_count": int(len(REQUIRED_REVIEW_ASSETS)),
        "compiler_passed_envelope": bool(compiler_metrics.get("passed_envelope", False)),
        "compiler_broad_land_to_water_fraction": float(
            compiler_metrics.get("broad_land_to_water_fraction", 0.0)),
        "compiler_broad_ocean_to_land_fraction": float(
            compiler_metrics.get("broad_ocean_to_land_fraction", 0.0)),
        "compiler_shelf_as_deep_ocean_fraction": float(
            compiler_metrics.get("shelf_as_deep_ocean_fraction", 0.0)),
        "compiler_lowland_as_mountain_fraction": float(
            compiler_metrics.get("lowland_as_mountain_fraction", 0.0)),
        "compiler_terrain_elevation_sign_mismatch_fraction": float(
            compiler_metrics.get("terrain_elevation_sign_mismatch_fraction", 0.0)),
        "hypsometry_out_of_envelope_count": int(
            len(earth_reference["out_of_envelope"])),
        "hypsometry_core_envelope_pass": bool(
            hyp_comparison["core_hypsometry_envelope_pass"]),
        "province_missing_reference_class_count": int(
            len(province_comparison["missing_reference_classes"])),
        "boundary_missing_process_count": int(
            len(boundary["current_generated"]["missing_process_types"])),
        "wilson_missing_object_set_count": int(
            len(wilson["current_generated"]["missing_object_sets"])),
        "source_to_sink_missing_object_count": int(
            sediment_budget["current_generated"]["metrics"]["missing_object_count"]),
        "drainage_missing_item_count": int(
            drainage["current_generated"]["metrics"]["missing_item_count"]),
        "old_orogen_missing_decay_field_count": int(
            len(old_orogen["current_generated"]["missing_decay_fields"])),
        "mountain_missing_inventory_field_count": int(
            mountain["current_generated"]["metrics"][
                "missing_inventory_field_count"]),
        "plateau_missing_item_count": int(
            plateau["current_generated"]["metrics"]["missing_plateau_item_count"]),
        **direct,
    }
    acceptance = {
        "world_core_fields_available": all(
            name in getattr(world, "fields", {})
            for name in (
                "terrain.elevation_m",
                "terrain.continental_detail",
                "terrain.province",
                "crust.type",
                "crust.domain",
                "ocean.depth_province",
            )
        ),
        "p76_p89_statuses_available": all(
            bool(value) for value in (
                source_ledger["status"],
                fixture["status"],
                earth_reference["status"],
                province_comparison["status"],
                boundary["status"],
                wilson["status"],
                crust_sediment["status"],
                sediment_budget["status"],
                drainage["status"],
                old_orogen["status"],
                mountain["status"],
                rift_margin["status"],
                plateau["status"],
            )
        ),
        "gap_inventory_nonempty": bool(gaps),
        "all_gaps_have_owner_category_and_future_stage": not unassigned,
        "no_generic_visual_blockers": not generic_blockers,
        "required_owner_layers_covered": set(OWNER_LAYERS).issubset(
            set(layer_summaries)),
        "direct_surface_metrics_available": (
            metrics["major_component_count"] >= 1
            and metrics["high_flat_interior_fraction_of_continental_land"] >= 0.0
            and metrics["basin_lowland_fraction_of_continental_land"] >= 0.0
        ),
        "compiler_metrics_available": not bool(
            compiler_metrics.get("invalid_source_array_shape", False)),
        "compiler_consistency_recorded": bool(
            "passed_envelope" in compiler_metrics),
        "asset_review_requirements_defined": (
            metrics["required_review_asset_count"] == len(REQUIRED_REVIEW_ASSETS)
        ),
        "p91_promotion_audit_pending": True,
    }
    return {
        "schema": SCHEMA,
        "status": (
            "current_world_morphology_gap_inventory_ready"
            if all(acceptance.values())
            else "current_world_morphology_gap_inventory_incomplete"
        ),
        "context": {
            "spec_name": str(world.spec.name),
            "seed": int(world.spec.seed),
            "cells": int(world.grid.n),
            "time_myr": float(world.time_myr),
        },
        "owner_layers": OWNER_LAYERS,
        "gap_categories": GAP_CATEGORIES,
        "required_review_assets": REQUIRED_REVIEW_ASSETS,
        "missing_review_assets": missing_assets,
        "current_residual_items": current_residual_items,
        "gaps": tuple(gaps),
        "layer_summaries": layer_summaries,
        "owner_counts": dict(sorted(owner_counts.items())),
        "category_counts": dict(sorted(category_counts.items())),
        "source_suite_counts": dict(sorted(source_counts.items())),
        "unassigned_gaps": unassigned,
        "generic_blockers": generic_blockers,
        "direct_surface_metrics": direct,
        "compiler_metrics": dict(compiler_metrics),
        "p76_p89_statuses": {
            "P76.source_ledger": source_ledger["status"],
            "P77.fixture": fixture["status"],
            "P78.hypsometry": earth_reference["status"],
            "P80.province_graph": province_comparison["status"],
            "P81.boundary_geometry": boundary["status"],
            "P82.wilson_lifecycle": wilson["status"],
            "P83.crust_sediment": crust_sediment["status"],
            "P84.source_to_sink": sediment_budget["status"],
            "P85.drainage": drainage["status"],
            "P86.old_orogen": old_orogen["status"],
            "P87.mountain_inventory": mountain["status"],
            "P88.rift_margin": rift_margin["status"],
            "P89.plateau": plateau["status"],
        },
        "metrics": metrics,
        "acceptance": acceptance,
        "summary_digest": _canonical_digest({
            "gaps": gaps,
            "metrics": metrics,
            "compiler": compiler_metrics,
        }),
    }
