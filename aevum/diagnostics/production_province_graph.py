"""Production continental province graph diagnostics for P94."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from typing import Any

import numpy as np

from aevum.diagnostics.generated_province_reference import (
    generated_province_reference_comparison_summary,
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
    PROVINCE_PARENT_PROCESS_NAMES,
)


SCHEMA = "aevum.production_province_graph.v1"

PRODUCTION_FIELD_NAMES = (
    "terrain.continental_province_id",
    "terrain.continental_province_code",
    "tectonics.continental_province_id",
    "tectonics.continental_province_code",
    "tectonics.province_parent_process",
)

REQUIRED_PROVINCE_CLASSES = (
    "shield",
    "platform",
    "intracratonic_basin",
    "foreland_basin",
    "active_orogen",
    "old_orogen",
    "rift_system",
    "passive_margin_lowland",
    "volcanic_lip_plateau",
)

REQUIRED_PROVINCE_EDGES = ("rift_system|volcanic_lip_plateau",)


def _digest(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def _edge_key(a: str, b: str) -> str:
    return "|".join(sorted((str(a), str(b))))


def production_province_graph_summary(world: Any) -> dict[str, Any]:
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total_area = max(float(area.sum()), 1.0e-12)
    elevation = np.asarray(world.get_field("terrain.elevation_m", 0.0), dtype=np.float64)
    rel = elevation - float(world.sea_level)
    crust_type = np.asarray(world.get_field("crust.type", 0.0), dtype=np.float64).astype(int)
    continental = crust_type == int(CONT)
    exposed_continental = (rel >= 0.0) & (crust_type == int(CONT))
    continental_area = max(float(area[continental].sum()), 1.0e-12)
    exposed_area = max(float(area[exposed_continental].sum()), 1.0e-12)

    fields = getattr(world, "fields", {})
    missing_fields = tuple(name for name in PRODUCTION_FIELD_NAMES if name not in fields)
    province_id = np.asarray(
        world.get_field("tectonics.continental_province_id", 0.0),
        dtype=np.float64,
    ).astype(int)
    province_code = np.asarray(
        world.get_field("tectonics.continental_province_code", 0.0),
        dtype=np.float64,
    ).astype(int)
    parent_process = np.asarray(
        world.get_field("tectonics.province_parent_process", 0.0),
        dtype=np.float64,
    ).astype(int)
    terrain_id = np.asarray(
        world.get_field("terrain.continental_province_id", 0.0),
        dtype=np.float64,
    ).astype(int)
    terrain_code = np.asarray(
        world.get_field("terrain.continental_province_code", 0.0),
        dtype=np.float64,
    ).astype(int)
    objects = list(world.objects.get("tectonics.continental_provinces", []))

    class_field = np.full(grid.n, "", dtype=object)
    for code, name in CONT_PROVINCE_CLASS_NAMES.items():
        if int(code) <= 0 or name == "none":
            continue
        class_field[continental & (province_code == int(code))] = name
    parent_field = np.full(grid.n, "", dtype=object)
    for code, name in PROVINCE_PARENT_PROCESS_NAMES.items():
        if int(code) <= 0 or name == "none":
            continue
        parent_field[continental & (parent_process == int(code))] = name

    exposed_id_mask = exposed_continental & (province_id > 0)
    exposed_code_mask = exposed_continental & (province_code > 0)
    exposed_parent_mask = exposed_continental & (parent_process > 0)
    continental_id_mask = continental & (province_id > 0)
    continental_code_mask = continental & (province_code > 0)
    continental_parent_mask = continental & (parent_process > 0)
    class_counts = Counter(str(value) for value in class_field[continental_code_mask])
    parent_counts = Counter(str(value) for value in parent_field[continental_parent_mask])
    class_area_fractions = {
        name: float(area[continental & (class_field == name)].sum() / continental_area)
        for name in sorted(set(class_counts))
    }
    parent_area_fractions = {
        name: float(area[continental & (parent_field == name)].sum() / continental_area)
        for name in sorted(set(parent_counts))
    }

    class_edge_counts = _class_edge_counts(grid, continental, class_field)
    province_edge_count = 0
    adjacent_by_id: dict[int, set[int]] = {}
    for a, b in grid.edges:
        a = int(a)
        b = int(b)
        ida = int(province_id[a])
        idb = int(province_id[b])
        if ida <= 0 or idb <= 0 or ida == idb:
            continue
        province_edge_count += 1
        adjacent_by_id.setdefault(ida, set()).add(idb)
        adjacent_by_id.setdefault(idb, set()).add(ida)

    object_ids = {int(obj.get("province_id", -1)) for obj in objects}
    field_ids = {int(x) for x in np.unique(province_id[continental_id_mask]) if int(x) > 0}
    missing_object_for_field_ids = tuple(sorted(field_ids - object_ids))
    object_without_field_ids = tuple(sorted(object_ids - field_ids))
    object_parent_failures = tuple(
        str(obj.get("id", ""))
        for obj in objects
        if not obj.get("parent_process") or str(obj.get("parent_process")) == "none"
    )
    object_adjacency_failures = tuple(
        str(obj.get("id", ""))
        for obj in objects
        if int(obj.get("province_id", -1)) in adjacent_by_id
        and not obj.get("adjacent_province_ids")
    )

    morphology = compute_world_morphology(world)
    labels = morphology.fields[
        "tectonics.exposed_continental_component_id"
    ].astype(int)
    major_components = []
    for cid in sorted(int(x) for x in np.unique(labels[exposed_continental]) if int(x) >= 0):
        mask = exposed_continental & (labels == cid)
        comp_area = float(area[mask].sum())
        if comp_area <= 0.0:
            continue
        area_fraction = comp_area / exposed_area
        if area_fraction < 0.08:
            continue
        component_ids = sorted(int(x) for x in np.unique(province_id[mask]) if int(x) > 0)
        component_classes = sorted(
            str(x) for x in set(class_field[mask]) if str(x)
        )
        component_class_fractions = {
            name: float(area[mask & (class_field == name)].sum() / comp_area)
            for name in component_classes
        }
        component_id_fractions = {
            int(pid): float(area[mask & (province_id == int(pid))].sum() / comp_area)
            for pid in component_ids
        }
        major_components.append({
            "component_id": int(cid),
            "area_fraction_of_exposed_continental_land": float(area_fraction),
            "province_id_count": int(len(component_ids)),
            "province_class_count": int(len(component_classes)),
            "largest_province_fraction": max(component_id_fractions.values(), default=0.0),
            "largest_class_fraction": max(component_class_fractions.values(), default=0.0),
            "province_classes": tuple(component_classes),
            "acceptance": {
                "multi_province": len(component_ids) >= 4,
                "multi_class": len(component_classes) >= 5,
                "dominant_province_capped": max(
                    component_id_fractions.values(), default=0.0) <= 0.60,
                "has_rift_or_margin": bool(
                    {"rift_system", "passive_margin_lowland"} & set(component_classes)),
                "has_lip_or_orogen": bool(
                    {"volcanic_lip_plateau", "active_orogen", "old_orogen"}
                    & set(component_classes)),
            },
        })

    id_component_count_by_id = {}
    disconnected_id_count = 0
    for pid in sorted(field_ids):
        comps = _components_from_edges(grid, province_id == pid)
        id_component_count_by_id[pid] = len(comps)
        if len(comps) != 1:
            disconnected_id_count += 1
    anchor_classes = {"shield", "platform", "passive_margin_lowland", "foreland_basin"}
    tiny_area = float(sum(
        float(obj.get("area_fraction", 0.0))
        for obj in objects
        if int(obj.get("cell_count", 0)) <= 1
        and str(obj.get("province_class", "")) not in anchor_classes
    ))

    reference_graph = province_reference_graph_summary()
    diversity = generated_world_province_diversity_summary(world)
    comparison = generated_province_reference_comparison_summary(
        world,
        reference_graph=reference_graph,
        diversity=diversity,
    )

    required_class_set = set(REQUIRED_PROVINCE_CLASSES)
    missing_required_classes = tuple(sorted(required_class_set - set(class_counts)))
    missing_required_edges = tuple(
        edge for edge in REQUIRED_PROVINCE_EDGES if edge not in class_edge_counts)
    field_alias_match = bool(
        np.array_equal(province_id, terrain_id)
        and np.array_equal(province_code, terrain_code)
    )
    acceptance = {
        "production_fields_present": not missing_fields,
        "terrain_and_tectonics_aliases_match": field_alias_match,
        "province_object_graph_present": len(objects) > 0,
        "continental_id_coverage": (
            float(area[continental_id_mask].sum() / continental_area) >= 0.98),
        "continental_code_coverage": (
            float(area[continental_code_mask].sum() / continental_area) >= 0.98),
        "continental_parent_process_coverage": (
            float(area[continental_parent_mask].sum() / continental_area) >= 0.98),
        "exposed_continental_id_coverage": (
            float(area[exposed_id_mask].sum() / exposed_area) >= 0.98),
        "exposed_continental_code_coverage": (
            float(area[exposed_code_mask].sum() / exposed_area) >= 0.98),
        "exposed_continental_parent_process_coverage": (
            float(area[exposed_parent_mask].sum() / exposed_area) >= 0.98),
        "field_object_id_consistency": (
            not missing_object_for_field_ids and not object_without_field_ids),
        "object_parent_processes_present": not object_parent_failures,
        "object_adjacency_recorded": not object_adjacency_failures,
        "required_reference_classes_covered": not missing_required_classes,
        "required_rift_lip_edge_present": not missing_required_edges,
        "major_components_multi_province": all(
            comp["acceptance"]["multi_province"] for comp in major_components),
        "major_components_multi_class": all(
            comp["acceptance"]["multi_class"] for comp in major_components),
        "dominant_province_capped": all(
            comp["acceptance"]["dominant_province_capped"] for comp in major_components),
        "major_components_have_margin_or_rift": all(
            comp["acceptance"]["has_rift_or_margin"] for comp in major_components),
        "major_components_have_lip_or_orogen": all(
            comp["acceptance"]["has_lip_or_orogen"] for comp in major_components),
        "no_checkerboard_province_graph": (
            disconnected_id_count == 0
            and tiny_area <= 0.02
            and len(objects) <= max(8, int(np.count_nonzero(exposed_continental) * 0.45))
        ),
    }
    metrics = {
        "cells": int(grid.n),
        "exposed_continental_cell_count": int(np.count_nonzero(exposed_continental)),
        "continental_id_coverage_fraction": float(
            area[continental_id_mask].sum() / continental_area),
        "continental_code_coverage_fraction": float(
            area[continental_code_mask].sum() / continental_area),
        "continental_parent_process_coverage_fraction": float(
            area[continental_parent_mask].sum() / continental_area),
        "exposed_continental_id_coverage_fraction": float(
            area[exposed_id_mask].sum() / exposed_area),
        "exposed_continental_code_coverage_fraction": float(
            area[exposed_code_mask].sum() / exposed_area),
        "exposed_continental_parent_process_coverage_fraction": float(
            area[exposed_parent_mask].sum() / exposed_area),
        "province_object_count": int(len(objects)),
        "province_field_id_count": int(len(field_ids)),
        "province_class_count": int(len(class_counts)),
        "parent_process_count": int(len(parent_counts)),
        "province_class_edge_count": int(len(class_edge_counts)),
        "province_id_edge_count": int(province_edge_count),
        "major_component_count": int(len(major_components)),
        "min_major_component_province_id_count": min(
            (comp["province_id_count"] for comp in major_components), default=0),
        "min_major_component_province_class_count": min(
            (comp["province_class_count"] for comp in major_components), default=0),
        "max_major_component_largest_province_fraction": max(
            (comp["largest_province_fraction"] for comp in major_components),
            default=0.0),
        "max_major_component_largest_class_fraction": max(
            (comp["largest_class_fraction"] for comp in major_components),
            default=0.0),
        "missing_field_count": int(len(missing_fields)),
        "missing_required_class_count": int(len(missing_required_classes)),
        "missing_required_edge_count": int(len(missing_required_edges)),
        "object_parent_failure_count": int(len(object_parent_failures)),
        "object_adjacency_failure_count": int(len(object_adjacency_failures)),
        "missing_object_for_field_id_count": int(len(missing_object_for_field_ids)),
        "object_without_field_id_count": int(len(object_without_field_ids)),
        "disconnected_province_id_count": int(disconnected_id_count),
        "tiny_province_area_fraction": float(tiny_area),
        "p80_missing_reference_class_count": int(
            comparison["metrics"]["missing_reference_class_count"]),
        "p80_missing_required_class_edge_count": int(
            comparison["metrics"]["missing_required_class_edge_count"]),
    }
    summary = {
        "schema": SCHEMA,
        "status": (
            "production_province_graph_ready"
            if all(acceptance.values())
            else "production_province_graph_incomplete"
        ),
        "context": {
            "spec_name": str(world.spec.name),
            "seed": int(world.spec.seed),
            "cells": int(grid.n),
            "time_myr": float(world.time_myr),
        },
        "metrics": metrics,
        "acceptance": acceptance,
        "missing_fields": missing_fields,
        "province_classes": tuple(sorted(class_counts)),
        "province_class_area_fractions": class_area_fractions,
        "parent_processes": tuple(sorted(parent_counts)),
        "parent_process_area_fractions": parent_area_fractions,
        "class_edge_counts": class_edge_counts,
        "missing_required_classes": missing_required_classes,
        "missing_required_edges": missing_required_edges,
        "major_components": tuple(major_components),
        "missing_object_for_field_ids": missing_object_for_field_ids,
        "object_without_field_ids": object_without_field_ids,
        "object_parent_failures": object_parent_failures,
        "object_adjacency_failures": object_adjacency_failures,
        "id_component_count_by_id": id_component_count_by_id,
        "reference_comparison_status": comparison["status"],
        "reference_comparison_digest": comparison["comparison_digest"],
        "comparison_missing_reference_classes": comparison["missing_reference_classes"],
        "comparison_missing_required_class_edges": comparison[
            "missing_required_class_edges"],
    }
    summary["summary_digest"] = _digest(summary)
    return summary


def _class_edge_counts(
    grid: Any,
    mask: np.ndarray,
    class_field: np.ndarray,
) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for a, b in grid.edges:
        a = int(a)
        b = int(b)
        if not (mask[a] and mask[b]):
            continue
        ca = str(class_field[a])
        cb = str(class_field[b])
        if not ca or not cb or ca == cb:
            continue
        counts[_edge_key(ca, cb)] += 1
    return dict(sorted(counts.items()))


def _components_from_edges(grid: Any, mask: np.ndarray) -> list[np.ndarray]:
    mask = np.asarray(mask, dtype=bool)
    remaining = set(int(x) for x in np.where(mask)[0])
    comps: list[np.ndarray] = []
    neighbors: dict[int, list[int]] = {int(i): [] for i in remaining}
    for a, b in grid.edges:
        a = int(a)
        b = int(b)
        if a in remaining and b in remaining:
            neighbors[a].append(b)
            neighbors[b].append(a)
    while remaining:
        start = remaining.pop()
        stack = [start]
        comp = [start]
        while stack:
            node = stack.pop()
            for nb in neighbors.get(node, ()):
                if nb in remaining:
                    remaining.remove(nb)
                    stack.append(nb)
                    comp.append(nb)
        comps.append(np.asarray(comp, dtype=np.int64))
    return comps
