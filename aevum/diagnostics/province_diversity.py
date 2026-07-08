"""Generated-world physiographic province diversity diagnostics."""
from __future__ import annotations

from collections import Counter
from typing import Any

import numpy as np

from aevum.diagnostics.morphology import compute_world_morphology
from aevum.modules.tectonics import (
    CONT,
    DOMAIN_CRATON,
    DOMAIN_LIP,
    DOMAIN_SUTURE,
)
from aevum.modules.terrain import (
    CONT_DETAIL_ARC_MICROCONTINENT,
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATFORM,
    CONT_DETAIL_PLATEAU,
    CONT_DETAIL_RIFT_BASIN,
    CONT_DETAIL_SHIELD,
)


SCHEMA = "aevum.generated_world_province_diversity.v1"

DETAIL_CLASS_NAMES = {
    int(CONT_DETAIL_SHIELD): "shield",
    int(CONT_DETAIL_PLATFORM): "platform",
    int(CONT_DETAIL_BASIN): "basin",
    int(CONT_DETAIL_RIFT_BASIN): "rift_basin",
    int(CONT_DETAIL_OROGEN): "orogen",
    int(CONT_DETAIL_PLATEAU): "plateau",
    int(CONT_DETAIL_ARC_MICROCONTINENT): "arc_microcontinent",
}

PARENTED_HIGHLAND_DETAIL_CLASSES = {
    "shield",
    "rift_basin",
    "orogen",
    "plateau",
    "arc_microcontinent",
}


def generated_world_province_diversity_summary(
    world,
    *,
    major_component_min_fraction: float = 0.08,
    min_detail_share: float = 0.02,
    lowland_cut_m: float = 500.0,
    highland_cut_m: float = 2500.0,
    parented_highland_cut_m: float = 1600.0,
) -> dict[str, Any]:
    morphology = compute_world_morphology(world)
    labels = morphology.fields["tectonics.exposed_continental_component_id"].astype(int)
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total_area = max(float(area.sum()), 1.0e-12)
    elev = np.asarray(world.get_field("terrain.elevation_m", 0.0), dtype=np.float64)
    rel = elev - float(world.sea_level)
    crust_type = np.asarray(world.get_field("crust.type", 0.0), dtype=np.float64).astype(int)
    domain = np.asarray(world.get_field("crust.domain", 0.0), dtype=np.float64).astype(int)
    detail = np.asarray(
        world.get_field("terrain.continental_detail", 0.0),
        dtype=np.float64,
    ).astype(int)
    exposed_continental = (rel >= 0.0) & (crust_type == int(CONT))
    exposed_cont_area = _area(area, exposed_continental)
    detail_global = _detail_fractions(area, detail, exposed_continental, exposed_cont_area)
    landform_objects = list(world.objects.get("terrain.continental_landforms", []))
    landform_kind_counts = Counter(str(obj.get("kind", "unknown")) for obj in landform_objects)
    landform_parent_processes = sorted({
        str(obj.get("parent_process", ""))
        for obj in landform_objects
        if obj.get("parent_process")
    })
    passive_margin_lowland_objects = [
        obj for obj in landform_objects
        if str(obj.get("kind", "")) == "passive_margin_lowland"
    ]
    passive_margin_lowland_field = np.asarray(
        world.get_field("terrain.passive_margin_lowland", np.zeros(grid.n)),
        dtype=np.float64,
    ) > 0.0
    unparented_landform_objects = [
        str(obj.get("id", ""))
        for obj in landform_objects
        if not obj.get("parent_process")
    ]

    components: list[dict[str, Any]] = []
    for cid in sorted(int(x) for x in np.unique(labels[exposed_continental]) if int(x) >= 0):
        mask = exposed_continental & (labels == cid)
        comp_area = _area(area, mask)
        if comp_area <= 0.0:
            continue
        area_fraction = comp_area / exposed_cont_area if exposed_cont_area else 0.0
        detail_fractions = _detail_fractions(area, detail, mask, comp_area)
        component = _component_summary(
            world,
            cid,
            mask,
            area,
            rel,
            domain,
            detail,
            detail_fractions,
            comp_area,
            exposed_cont_area,
            total_area,
            landform_objects,
            area_fraction,
            min_detail_share=min_detail_share,
            lowland_cut_m=lowland_cut_m,
            highland_cut_m=highland_cut_m,
            parented_highland_cut_m=parented_highland_cut_m,
        )
        components.append(component)

    major_components = [
        comp for comp in components
        if comp["area_fraction_of_exposed_continental_land"] >= major_component_min_fraction
    ]
    failing_major_components = [
        comp for comp in major_components
        if not all(comp["acceptance"].values())
    ]
    detail_class_count_global = int(
        sum(1 for value in detail_global.values() if value >= min_detail_share)
    )
    metrics = {
        "cells": int(grid.n),
        "time_myr": float(world.time_myr),
        "exposed_continental_land_fraction_of_world": exposed_cont_area / total_area,
        "component_count": len(components),
        "major_component_count": len(major_components),
        "detail_class_count_global_gt2pct": detail_class_count_global,
        "terrain_landform_object_count": len(landform_objects),
        "terrain_landform_kind_count": len(landform_kind_counts),
        "terrain_landform_parent_process_count": len(landform_parent_processes),
        "passive_margin_lowland_object_count": len(passive_margin_lowland_objects),
        "passive_margin_lowland_object_area_fraction": float(
            sum(float(obj.get("area_fraction", 0.0)) for obj in passive_margin_lowland_objects)),
        "passive_margin_lowland_field_area_fraction": (
            _area(area, passive_margin_lowland_field) / total_area),
        "unparented_landform_object_count": len(unparented_landform_objects),
        "min_province_class_count_per_major": min(
            (comp["province_class_count_gt2pct"] for comp in major_components),
            default=0,
        ),
        "max_largest_internal_province_fraction": max(
            (comp["largest_internal_province_fraction"] for comp in major_components),
            default=0.0,
        ),
        "min_basin_or_lowland_share_per_major": min(
            (comp["basin_or_lowland_share"] for comp in major_components),
            default=0.0,
        ),
        "max_active_highland_or_plateau_fraction": max(
            (comp["active_highland_or_plateau_fraction"] for comp in major_components),
            default=0.0,
        ),
        "max_unparented_highland_fraction": max(
            (comp["unparented_highland_fraction_of_highlands"] for comp in major_components),
            default=0.0,
        ),
        "failing_major_component_count": len(failing_major_components),
    }
    acceptance = {
        "major_continent_count_multiple": metrics["major_component_count"] >= 2,
        "all_major_components_have_three_province_classes": all(
            comp["province_class_count_gt2pct"] >= 3 for comp in major_components),
        "largest_internal_province_fraction_capped": (
            metrics["max_largest_internal_province_fraction"] <= 0.74
        ),
        "basin_or_lowland_present_per_major": all(
            comp["acceptance"]["basin_or_lowland_present"] for comp in major_components),
        "active_highland_area_limited_per_major": all(
            comp["acceptance"]["active_highland_area_limited"] for comp in major_components),
        "highlands_parented_per_major": all(
            comp["acceptance"]["highlands_parented"] for comp in major_components),
        "terrain_landform_parent_processes_present": (
            metrics["unparented_landform_object_count"] == 0
        ),
        "passive_margin_lowland_object_layer_present": (
            metrics["passive_margin_lowland_object_count"] > 0
            and metrics["passive_margin_lowland_field_area_fraction"] > 0.0
        ),
        "generated_world_province_signal_present": (
            metrics["detail_class_count_global_gt2pct"] >= 5
            and metrics["terrain_landform_kind_count"] >= 5
        ),
    }
    return {
        "schema": SCHEMA,
        "status": "generated_world_gate_ready" if all(acceptance.values()) else "generated_world_gate_incomplete",
        "context": {
            "spec_name": str(world.spec.name),
            "seed": int(world.spec.seed),
            "cells": int(grid.n),
            "time_myr": float(world.time_myr),
        },
        "parameters": {
            "major_component_min_fraction": float(major_component_min_fraction),
            "min_detail_share": float(min_detail_share),
            "lowland_cut_m": float(lowland_cut_m),
            "highland_cut_m": float(highland_cut_m),
            "parented_highland_cut_m": float(parented_highland_cut_m),
        },
        "metrics": metrics,
        "acceptance": acceptance,
        "detail_class_fractions_global": detail_global,
        "terrain_landform_kind_counts": dict(sorted(landform_kind_counts.items())),
        "terrain_landform_parent_processes": landform_parent_processes,
        "unparented_landform_objects": unparented_landform_objects,
        "major_components": major_components,
        "failing_major_components": [
            {
                "component_id": comp["component_id"],
                "failed_acceptance": [
                    key for key, value in comp["acceptance"].items() if not value
                ],
            }
            for comp in failing_major_components
        ],
        "limitations": {
            "passive_margin_lowland_detail_class_missing": True,
            "passive_margin_lowland_object_layer_present": (
                metrics["passive_margin_lowland_object_count"] > 0),
            "terrain_rewrite_pending": False,
            "terrain_rewrite_stage": "p74_object_template_active",
        },
        "next_gates": (
            "P73.real_earth_case_study_calibration",
            "P74.terrain_coupling_rewrite",
        ),
    }


def _component_summary(
    world,
    component_id: int,
    mask: np.ndarray,
    area: np.ndarray,
    rel: np.ndarray,
    domain: np.ndarray,
    detail: np.ndarray,
    detail_fractions: dict[str, float],
    comp_area: float,
    exposed_cont_area: float,
    total_area: float,
    landform_objects: list[dict[str, Any]],
    area_fraction: float,
    *,
    min_detail_share: float,
    lowland_cut_m: float,
    highland_cut_m: float,
    parented_highland_cut_m: float,
) -> dict[str, Any]:
    del total_area
    province_class_count = int(
        sum(1 for value in detail_fractions.values() if value >= min_detail_share)
    )
    largest_internal = max(detail_fractions.values(), default=0.0)
    basin_or_lowland_share = (
        detail_fractions.get("basin", 0.0)
        + detail_fractions.get("rift_basin", 0.0)
        + _area(area, mask & (rel < float(lowland_cut_m))) / comp_area
    )
    active_highland_or_plateau = (
        mask
        & (
            np.isin(
                detail,
                [int(CONT_DETAIL_OROGEN), int(CONT_DETAIL_PLATEAU)],
            )
            | (rel > float(highland_cut_m))
        )
    )
    highland = mask & (rel > float(parented_highland_cut_m))
    parented_highland = highland & (
        np.isin(
            detail,
            [
                int(CONT_DETAIL_SHIELD),
                int(CONT_DETAIL_RIFT_BASIN),
                int(CONT_DETAIL_OROGEN),
                int(CONT_DETAIL_PLATEAU),
                int(CONT_DETAIL_ARC_MICROCONTINENT),
            ],
        )
        | np.isin(domain, [int(DOMAIN_CRATON), int(DOMAIN_SUTURE), int(DOMAIN_LIP)])
    )
    highland_area = _area(area, highland)
    unparented_highland_area = max(highland_area - _area(area, parented_highland), 0.0)
    object_kinds, object_parent_processes = _overlapping_landform_context(
        mask, landform_objects)
    lowland_lt500 = _area(area, mask & (rel < 500.0)) / comp_area
    lowland_lt1000 = _area(area, mask & (rel < 1000.0)) / comp_area
    active_highland_fraction = _area(area, active_highland_or_plateau) / comp_area
    unparented_highland_fraction = (
        unparented_highland_area / highland_area if highland_area > 0.0 else 0.0
    )
    acceptance = {
        "minimum_province_class_diversity": province_class_count >= 3,
        "largest_internal_province_fraction_capped": largest_internal <= 0.74,
        "basin_or_lowland_present": (
            basin_or_lowland_share >= 0.18
            or lowland_lt500 >= 0.15
            or lowland_lt1000 >= 0.30
        ),
        "active_highland_area_limited": active_highland_fraction <= 0.45,
        "highlands_parented": unparented_highland_fraction <= 0.15,
        "landform_object_context_present": len(object_kinds) >= 2,
    }
    return {
        "component_id": int(component_id),
        "area_fraction_of_exposed_continental_land": float(area_fraction),
        "area_fraction_of_world": float(comp_area / max(float(area.sum()), 1.0e-12)),
        "cell_count": int(mask.sum()),
        "province_class_count_gt2pct": province_class_count,
        "largest_internal_province_fraction": float(largest_internal),
        "basin_or_lowland_share": float(basin_or_lowland_share),
        "lowland_fraction_lt500m": float(lowland_lt500),
        "lowland_fraction_lt1000m": float(lowland_lt1000),
        "active_highland_or_plateau_fraction": float(active_highland_fraction),
        "highland_fraction_gt2500m": float(
            _area(area, mask & (rel > float(highland_cut_m))) / comp_area),
        "highland_fraction_gt_parentage_cut": float(highland_area / comp_area),
        "unparented_highland_fraction_of_highlands": float(unparented_highland_fraction),
        "relief_p90_p10_m": _percentile(rel[mask], 90.0) - _percentile(rel[mask], 10.0),
        "relief_p95_p05_m": _percentile(rel[mask], 95.0) - _percentile(rel[mask], 5.0),
        "elevation_p50_m": _percentile(rel[mask], 50.0),
        "elevation_p95_m": _percentile(rel[mask], 95.0),
        "detail_class_fractions": detail_fractions,
        "landform_object_kind_count": len(object_kinds),
        "landform_parent_process_count": len(object_parent_processes),
        "landform_object_kinds": object_kinds,
        "landform_parent_processes": object_parent_processes,
        "acceptance": acceptance,
    }


def _detail_fractions(
    area: np.ndarray,
    detail: np.ndarray,
    mask: np.ndarray,
    mask_area: float,
) -> dict[str, float]:
    out = {name: 0.0 for name in DETAIL_CLASS_NAMES.values()}
    if mask_area <= 0.0:
        return out
    for code, name in DETAIL_CLASS_NAMES.items():
        out[name] = _area(area, mask & (detail == int(code))) / mask_area
    return out


def _overlapping_landform_context(
    mask: np.ndarray,
    landform_objects: list[dict[str, Any]],
) -> tuple[list[str], list[str]]:
    kinds: set[str] = set()
    processes: set[str] = set()
    mask_indices = set(int(x) for x in np.where(mask)[0])
    for obj in landform_objects:
        cells = obj.get("cells", [])
        if cells:
            if not mask_indices.intersection(int(cell) for cell in cells):
                continue
        else:
            # Large objects may omit cell lists; keep their global process
            # context visible without assigning them to a specific component.
            continue
        kind = str(obj.get("kind", ""))
        process = str(obj.get("parent_process", ""))
        if kind:
            kinds.add(kind)
        if process:
            processes.add(process)
    return sorted(kinds), sorted(processes)


def _area(area: np.ndarray, mask: np.ndarray) -> float:
    mask = np.asarray(mask, dtype=bool)
    return float(area[mask].sum()) if mask.any() else 0.0


def _percentile(values: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return 0.0
    return float(np.percentile(values, q))
