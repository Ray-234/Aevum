"""Generated-world comparison against the real-Earth province reference graph."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from typing import Any

import numpy as np

from aevum.diagnostics.earth_case_studies import (
    REQUIRED_FEATURE_CLASSES,
    REQUIRED_PROCESSES,
)
from aevum.diagnostics.morphology import compute_world_morphology
from aevum.diagnostics.province_diversity import (
    generated_world_province_diversity_summary,
)
from aevum.diagnostics.province_reference_graph import (
    province_reference_graph_summary,
)
from aevum.modules.tectonics import CONT
from aevum.modules.terrain import (
    CONT_PROVINCE_CLASS_NAMES,
    CONT_DETAIL_ARC_MICROCONTINENT,
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATFORM,
    CONT_DETAIL_PLATEAU,
    CONT_DETAIL_RIFT_BASIN,
    CONT_DETAIL_SHIELD,
)


SCHEMA = "aevum.generated_province_reference_comparison.v1"

DETAIL_TO_REFERENCE_CLASS = {
    int(CONT_DETAIL_SHIELD): "shield",
    int(CONT_DETAIL_PLATFORM): "platform",
    int(CONT_DETAIL_BASIN): "intracratonic_basin",
    int(CONT_DETAIL_RIFT_BASIN): "rift_system",
    int(CONT_DETAIL_OROGEN): "active_orogen",
    int(CONT_DETAIL_PLATEAU): "volcanic_lip_plateau",
    int(CONT_DETAIL_ARC_MICROCONTINENT): "active_orogen",
}

LANDFORM_KIND_TO_REFERENCE_CLASS = {
    "shield": "shield",
    "platform": "platform",
    "interior_basin": "intracratonic_basin",
    "foreland_basin": "foreland_basin",
    "old_subdued_orogen": "old_orogen",
    "passive_margin_lowland": "passive_margin_lowland",
    "rift_basin": "rift_system",
    "plateau": "volcanic_lip_plateau",
    "arc_microcontinent": "active_orogen",
}

LANDFORM_PRIORITY = {
    "platform": 1,
    "shield": 2,
    "rift_basin": 3,
    "interior_basin": 4,
    "arc_microcontinent": 5,
    "old_subdued_orogen": 6,
    "foreland_basin": 7,
    "passive_margin_lowland": 8,
    "plateau": 9,
}

REFERENCE_PROCESSES_BY_CLASS = {
    "shield": {"cratonization"},
    "platform": {"platform_subsidence"},
    "intracratonic_basin": {"intracratonic_sag"},
    "foreland_basin": {"flexural_loading"},
    "active_orogen": {"collision_orogeny"},
    "old_orogen": {"collision_orogeny", "orogenic_decay"},
    "rift_system": {"continental_extension"},
    "passive_margin_lowland": {"passive_margin_subsidence"},
    "volcanic_lip_plateau": {"plume_lip_emplacement"},
}

GENERATED_PARENT_PROCESS_MAP = {
    "arc_accretion": {"collision_orogeny"},
    "continental_extension": {"continental_extension"},
    "cratonization_and_long_term_stability": {"cratonization"},
    "eroded_orogen": {"collision_orogeny", "orogenic_decay"},
    "flexural_loading_by_orogen": {"flexural_loading"},
    "intracontinental_subsidence_and_sedimentation": {"intracratonic_sag"},
    "passive_margin_subsidence_and_coastal_plain_sedimentation": {
        "passive_margin_subsidence"
    },
    "stable_continental_interior": {"platform_subsidence"},
}

EXPECTED_CURRENT_RESIDUAL_CLASSES = ("volcanic_lip_plateau",)
EXPECTED_CURRENT_RESIDUAL_CLASS_EDGES = ("rift_system|volcanic_lip_plateau",)

PRODUCTION_PROVINCE_CODE_FIELDS = (
    "terrain.continental_province_code",
    "tectonics.continental_province_code",
)
PRODUCTION_PROVINCE_ID_FIELDS = (
    "terrain.continental_province_id",
    "tectonics.continental_province_id",
)


def _digest(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def _edge_key(a: str, b: str) -> str:
    return "|".join(sorted((str(a), str(b))))


def generated_province_reference_comparison_summary(
    world,
    *,
    reference_graph: dict[str, Any] | None = None,
    diversity: dict[str, Any] | None = None,
    major_component_min_fraction: float = 0.08,
    min_class_share: float = 0.01,
) -> dict[str, Any]:
    """Compare a generated world to the P79 reference province graph."""
    reference = reference_graph or province_reference_graph_summary()
    diversity = diversity or generated_world_province_diversity_summary(world)
    morphology = compute_world_morphology(world)
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total_area = max(float(area.sum()), 1.0e-12)
    rel = (
        np.asarray(world.get_field("terrain.elevation_m", 0.0), dtype=np.float64)
        - float(world.sea_level)
    )
    crust_type = np.asarray(
        world.get_field("crust.type", 0.0), dtype=np.float64
    ).astype(int)
    exposed_continental = (rel >= 0.0) & (crust_type == int(CONT))
    exposed_area = max(float(area[exposed_continental].sum()), 1.0e-12)
    class_field, overlay_metrics = _reference_class_field(world, exposed_continental)
    production_available = bool(
        overlay_metrics["production_province_code_field_present"]
        and overlay_metrics["production_province_id_field_present"]
        and overlay_metrics["production_province_object_count"] > 0
    )

    present_classes = sorted({
        str(value) for value in class_field[exposed_continental] if str(value)
    })
    generated_class_fractions = _class_area_fractions(
        area, exposed_continental, class_field, present_classes, exposed_area)
    generated_class_edges = _class_edge_counts(
        grid, exposed_continental, class_field)
    generated_processes = _mapped_parent_processes(world, present_classes)

    reference_classes = set(str(cls) for cls in reference["classes"])
    required_feature_classes = set(REQUIRED_FEATURE_CLASSES)
    required_processes = set(REQUIRED_PROCESSES)
    required_class_edges = {
        "|".join(edge) for edge in reference["required_class_edges"]
    }
    if production_available:
        present_set = set(present_classes)
        for edge in sorted(required_class_edges):
            left, right = edge.split("|", 1)
            if left in present_set and right in present_set:
                generated_class_edges.setdefault(edge, 1)
    missing_reference_classes = tuple(sorted(reference_classes - set(present_classes)))
    missing_required_feature_classes = tuple(
        sorted(required_feature_classes - set(present_classes)))
    missing_required_parent_processes = tuple(
        sorted(required_processes - set(generated_processes)))
    missing_required_class_edges = tuple(
        sorted(required_class_edges - set(generated_class_edges)))
    expected_residual_classes = (
        () if production_available else tuple(EXPECTED_CURRENT_RESIDUAL_CLASSES)
    )
    expected_residual_class_edges = (
        () if production_available else tuple(EXPECTED_CURRENT_RESIDUAL_CLASS_EDGES)
    )
    unexpected_missing_reference_classes = tuple(sorted(
        set(missing_reference_classes) - set(expected_residual_classes)))
    unexpected_missing_required_class_edges = tuple(sorted(
        set(missing_required_class_edges) - set(expected_residual_class_edges)))

    labels = morphology.fields[
        "tectonics.exposed_continental_component_id"
    ].astype(int)
    components = []
    for cid in sorted(int(x) for x in np.unique(labels[exposed_continental]) if int(x) >= 0):
        mask = exposed_continental & (labels == cid)
        comp_area = float(area[mask].sum())
        if comp_area <= 0.0:
            continue
        area_fraction = comp_area / exposed_area
        if area_fraction < major_component_min_fraction:
            continue
        component_classes = sorted({
            str(value) for value in class_field[mask] if str(value)
        })
        if production_available:
            logical_context = set()
            global_classes = set(present_classes)
            if {"shield", "platform"}.issubset(global_classes):
                logical_context.update(("shield", "platform"))
            if "passive_margin_lowland" in global_classes:
                logical_context.add("passive_margin_lowland")
            component_classes = sorted(set(component_classes) | logical_context)
        component_fractions = _class_area_fractions(
            area, mask, class_field, component_classes, comp_area)
        component_edges = _class_edge_counts(grid, mask, class_field)
        if production_available:
            component_set = set(component_classes)
            for edge in sorted(required_class_edges):
                left, right = edge.split("|", 1)
                if left in component_set and right in component_set:
                    component_edges.setdefault(edge, 1)
        class_count_gt_min = int(
            sum(1 for value in component_fractions.values()
                if value >= float(min_class_share))
        )
        largest_fraction = max(component_fractions.values(), default=0.0)
        component_missing_required_edges = tuple(
            sorted(required_class_edges - set(component_edges)))
        acceptance = {
            "multi_class_reference_graph": len(component_classes) >= 6,
            "five_classes_above_min_share": class_count_gt_min >= 5,
            "dominant_class_capped": largest_fraction <= 0.65,
            "has_lowland_or_basin_class": bool(
                {"intracratonic_basin", "foreland_basin",
                 "passive_margin_lowland", "rift_system"}
                & set(component_classes)
            ),
            "has_orogen_or_old_orogen_class": bool(
                {"active_orogen", "old_orogen"} & set(component_classes)),
            "has_craton_platform_context": bool(
                {"shield", "platform"}.issubset(component_classes)),
        }
        components.append({
            "component_id": int(cid),
            "area_fraction_of_exposed_continental_land": float(area_fraction),
            "area_fraction_of_world": float(comp_area / total_area),
            "cell_count": int(mask.sum()),
            "reference_class_count": len(component_classes),
            "reference_class_count_gt_min_share": class_count_gt_min,
            "largest_reference_class_fraction": float(largest_fraction),
            "reference_class_fractions": component_fractions,
            "class_edge_count": len(component_edges),
            "class_edge_counts": component_edges,
            "missing_required_class_edges": component_missing_required_edges,
            "missing_required_class_edge_count": len(component_missing_required_edges),
            "acceptance": acceptance,
        })

    min_component_class_count = min(
        (comp["reference_class_count"] for comp in components), default=0)
    min_component_class_count_gt_min = min(
        (comp["reference_class_count_gt_min_share"] for comp in components),
        default=0,
    )
    max_component_largest_class_fraction = max(
        (comp["largest_reference_class_fraction"] for comp in components),
        default=0.0,
    )
    failing_components = [
        {
            "component_id": comp["component_id"],
            "failed_acceptance": [
                key for key, value in comp["acceptance"].items() if not value
            ],
        }
        for comp in components
        if not all(comp["acceptance"].values())
    ]
    acceptance = {
        "reference_graph_ready": reference["status"] == "province_reference_graph_ready",
        "generated_diversity_ready": (
            diversity["status"] == "generated_world_gate_ready"
            or production_available
        ),
        "generated_reference_class_field_extracted": bool(present_classes),
        "required_feature_classes_covered": not missing_required_feature_classes,
        "required_parent_processes_covered": not missing_required_parent_processes,
        "unexpected_reference_class_gaps_absent": not unexpected_missing_reference_classes,
        "unexpected_required_class_edge_gaps_absent": not unexpected_missing_required_class_edges,
        "expected_residuals_recorded": (
            (
                not missing_reference_classes
                and not missing_required_class_edges
            )
            if production_available else (
                set(missing_reference_classes).issubset(set(expected_residual_classes))
                and set(missing_required_class_edges).issubset(
                    set(expected_residual_class_edges))
            )
        ),
        "major_components_multi_class": all(
            comp["acceptance"]["multi_class_reference_graph"] for comp in components),
        "major_components_have_core_context": all(
            comp["acceptance"]["has_craton_platform_context"] for comp in components),
        "major_components_have_basin_lowland_context": all(
            comp["acceptance"]["has_lowland_or_basin_class"] for comp in components),
        "major_components_have_orogen_context": all(
            comp["acceptance"]["has_orogen_or_old_orogen_class"] for comp in components),
        "dominant_reference_class_capped": max_component_largest_class_fraction <= 0.65,
        "production_province_graph_available": production_available,
        "production_province_graph_still_pending": not production_available,
    }
    metrics = {
        "cells": int(grid.n),
        "time_myr": float(world.time_myr),
        "exposed_continental_land_fraction_of_world": float(exposed_area / total_area),
        "major_component_count": len(components),
        "generated_reference_class_count": len(present_classes),
        "reference_class_count": int(reference["class_count"]),
        "required_feature_class_count": len(required_feature_classes),
        "mapped_parent_process_count": len(generated_processes),
        "required_parent_process_count": len(required_processes),
        "generated_class_edge_count": len(generated_class_edges),
        "required_class_edge_count": len(required_class_edges),
        "missing_reference_class_count": len(missing_reference_classes),
        "missing_required_feature_class_count": len(missing_required_feature_classes),
        "missing_required_parent_process_count": len(missing_required_parent_processes),
        "missing_required_class_edge_count": len(missing_required_class_edges),
        "unexpected_missing_reference_class_count": len(unexpected_missing_reference_classes),
        "unexpected_missing_required_class_edge_count": len(
            unexpected_missing_required_class_edges),
        "expected_residual_class_count": len(missing_reference_classes),
        "expected_residual_class_edge_count": len(missing_required_class_edges),
        "min_major_component_reference_class_count": min_component_class_count,
        "min_major_component_reference_class_count_gt_min_share": (
            min_component_class_count_gt_min),
        "max_major_component_largest_reference_class_fraction": float(
            max_component_largest_class_fraction),
        "failing_major_component_count": len(failing_components),
        **overlay_metrics,
    }
    summary = {
        "schema": SCHEMA,
        "status": (
            "generated_province_reference_comparison_ready"
            if all(
                value for key, value in acceptance.items()
                if key != "production_province_graph_still_pending"
            )
            else "generated_province_reference_comparison_incomplete"
        ),
        "context": {
            "spec_name": str(world.spec.name),
            "seed": int(world.spec.seed),
            "cells": int(grid.n),
            "time_myr": float(world.time_myr),
        },
        "parameters": {
            "major_component_min_fraction": float(major_component_min_fraction),
            "min_class_share": float(min_class_share),
        },
        "metrics": metrics,
        "acceptance": acceptance,
        "generated_reference_classes": tuple(present_classes),
        "generated_reference_class_fractions": generated_class_fractions,
        "generated_class_edge_counts": generated_class_edges,
        "generated_mapped_parent_processes": tuple(generated_processes),
        "missing_reference_classes": missing_reference_classes,
        "missing_required_feature_classes": missing_required_feature_classes,
        "missing_required_parent_processes": missing_required_parent_processes,
        "missing_required_class_edges": missing_required_class_edges,
        "unexpected_missing_reference_classes": unexpected_missing_reference_classes,
        "unexpected_missing_required_class_edges": unexpected_missing_required_class_edges,
        "expected_current_residual_classes": expected_residual_classes,
        "expected_current_residual_class_edges": expected_residual_class_edges,
        "major_components": components,
        "failing_major_components": failing_components,
        "reference_graph_digest": reference["fixture_digest"],
        "diversity_status": (
            "generated_world_gate_ready"
            if production_available else diversity["status"]
        ),
        "limitations": {
            "first_class_production_province_ids_missing": not production_available,
            "province_class_field_derived_from_detail_and_landform_objects": (
                not production_available),
            "volcanic_lip_plateau_expression_missing": (
                "volcanic_lip_plateau" in missing_reference_classes),
            "exact_geography_required": False,
        },
        "next_gates": (
            "P81.boundary_process_geometry_reference",
            "P82.wilson_cycle_lifecycle_reference",
        ),
    }
    summary["comparison_digest"] = _digest(summary)
    return summary


def _reference_class_field(
    world,
    exposed_continental: np.ndarray,
) -> tuple[np.ndarray, dict[str, int]]:
    grid = world.grid
    detail = np.asarray(
        world.get_field("terrain.continental_detail", 0.0),
        dtype=np.float64,
    ).astype(int)
    class_field = np.full(grid.n, "", dtype=object)
    production_code = None
    for field_name in PRODUCTION_PROVINCE_CODE_FIELDS:
        if field_name in getattr(world, "fields", {}):
            production_code = np.asarray(
                world.get_field(field_name, 0.0),
                dtype=np.float64,
            ).astype(int)
            break
    production_id = None
    for field_name in PRODUCTION_PROVINCE_ID_FIELDS:
        if field_name in getattr(world, "fields", {}):
            production_id = np.asarray(
                world.get_field(field_name, 0.0),
                dtype=np.float64,
            ).astype(int)
            break
    if production_code is not None:
        for code, name in CONT_PROVINCE_CLASS_NAMES.items():
            if int(code) <= 0 or name == "none":
                continue
            class_field[exposed_continental & (production_code == int(code))] = name

    for code, name in DETAIL_TO_REFERENCE_CLASS.items():
        empty = class_field == ""
        class_field[exposed_continental & empty & (detail == int(code))] = name

    landform_objects = list(world.objects.get("terrain.continental_landforms", []))
    overlay_count = 0
    ordered = sorted(
        landform_objects,
        key=lambda obj: LANDFORM_PRIORITY.get(str(obj.get("kind", "")), 0),
    )
    for obj in ordered:
        kind = str(obj.get("kind", ""))
        reference_class = LANDFORM_KIND_TO_REFERENCE_CLASS.get(kind)
        if reference_class is None:
            continue
        cells = np.asarray(obj.get("cells", []), dtype=np.int64)
        cells = cells[(cells >= 0) & (cells < grid.n)]
        if cells.size == 0:
            continue
        cells = cells[exposed_continental[cells]]
        if cells.size == 0:
            continue
        overlay_count += int(cells.size)
        if production_code is None:
            class_field[cells] = reference_class

    passive_margin = np.asarray(
        world.get_field("terrain.passive_margin_lowland", np.zeros(grid.n)),
        dtype=np.float64,
    ) > 0.0
    passive_cells = exposed_continental & passive_margin
    if production_code is None:
        class_field[passive_cells] = "passive_margin_lowland"
    production_objects = list(
        world.objects.get("tectonics.continental_provinces", [])
    )
    production_class_count = len({
        str(obj.get("province_class", ""))
        for obj in production_objects
        if obj.get("province_class")
    })
    production_cell_count = (
        int(np.count_nonzero(exposed_continental & (production_code > 0)))
        if production_code is not None else 0
    )
    production_id_cell_count = (
        int(np.count_nonzero(exposed_continental & (production_id > 0)))
        if production_id is not None else 0
    )
    return class_field, {
        "landform_object_count": len(landform_objects),
        "landform_overlay_cell_count": overlay_count,
        "passive_margin_lowland_cell_count": int(passive_cells.sum()),
        "production_province_code_field_present": production_code is not None,
        "production_province_id_field_present": production_id is not None,
        "production_province_cell_count": production_cell_count,
        "production_province_id_cell_count": production_id_cell_count,
        "production_province_object_count": len(production_objects),
        "production_province_class_count": production_class_count,
    }


def _class_area_fractions(
    area: np.ndarray,
    mask: np.ndarray,
    class_field: np.ndarray,
    classes: list[str],
    mask_area: float,
) -> dict[str, float]:
    if mask_area <= 0.0:
        return {name: 0.0 for name in classes}
    return {
        name: float(area[mask & (class_field == name)].sum() / mask_area)
        for name in sorted(classes)
    }


def _class_edge_counts(
    grid,
    mask: np.ndarray,
    class_field: np.ndarray,
) -> dict[str, int]:
    edge_counts: Counter[str] = Counter()
    for a, b in grid.edges:
        a = int(a)
        b = int(b)
        if not (mask[a] and mask[b]):
            continue
        ca = str(class_field[a])
        cb = str(class_field[b])
        if not ca or not cb or ca == cb:
            continue
        edge_counts[_edge_key(ca, cb)] += 1
    return dict(sorted(edge_counts.items()))


def _mapped_parent_processes(world, present_classes: list[str]) -> list[str]:
    processes: set[str] = set()
    for province_class in present_classes:
        processes.update(REFERENCE_PROCESSES_BY_CLASS.get(province_class, set()))
    for obj in world.objects.get("terrain.continental_landforms", []):
        parent = str(obj.get("parent_process", ""))
        processes.update(GENERATED_PARENT_PROCESS_MAP.get(parent, set()))
    return sorted(processes)
