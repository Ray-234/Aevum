"""Continent and landform morphology diagnostics.

These helpers intentionally do not change generation.  They measure whether the
current crust and exposed-land masks contain broad continental interiors or are
dominated by narrow ribbons, necks, and slivers.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
import sys
from typing import Any

import numpy as np


@dataclass
class MaskMorphology:
    labels: np.ndarray
    width_steps: np.ndarray
    ribbon_index: np.ndarray
    complexity_index: np.ndarray
    narrow_neck_index: np.ndarray
    components: list[dict[str, Any]]
    metrics: dict[str, Any]


@dataclass
class WorldMorphology:
    fields: dict[str, np.ndarray]
    objects: dict[str, list[dict[str, Any]]]
    metrics: dict[str, Any]


def analyze_mask_morphology(grid, mask: np.ndarray, *, name: str = "mask") -> MaskMorphology:
    """Return component, width, ribbon, and coastline-complexity diagnostics."""
    mask = np.asarray(mask, dtype=bool)
    labels, cells_by_component = _component_labels(grid, mask)
    width = _width_steps(grid, mask)
    articulation = _articulation_points(grid, mask)
    components, complexity_field, ribbon_field = _component_summaries(
        grid, mask, labels, cells_by_component, width, articulation, name=name)
    narrow_neck = (mask & articulation & (width <= 2.0)).astype(np.float64)
    metrics = _mask_metrics(grid, mask, width, articulation, components, ribbon_field, name)
    return MaskMorphology(
        labels=labels.astype(np.float64),
        width_steps=width,
        ribbon_index=ribbon_field,
        complexity_index=complexity_field,
        narrow_neck_index=narrow_neck,
        components=components,
        metrics=metrics,
    )


def compute_world_morphology(world) -> WorldMorphology:
    """Compute P0 morphology diagnostics for the final world state."""
    grid = world.grid
    sea = world.sea_level
    elev = world.get_field("terrain.elevation_m", 0.0)
    rel = elev - sea
    land = rel >= 0.0
    crust = world.get_field("crust.type", 0.0).astype(int)
    continental = crust == 1
    oceanic = crust == 0
    exposed_continental = land & continental
    exposed_oceanic = land & oceanic

    land_diag = analyze_mask_morphology(grid, land, name="exposed_land")
    cont_diag = analyze_mask_morphology(grid, continental, name="continental_crust")
    exposed_cont_diag = analyze_mask_morphology(
        grid, exposed_continental, name="exposed_continental_land")

    area = grid.cell_area
    land_area = float(area[land].sum())
    cont_area = float(area[continental].sum())
    total_area = float(area.sum())

    high_oceanic_500 = exposed_oceanic & (rel > 500.0)
    high_oceanic_1500 = exposed_oceanic & (rel > 1500.0)
    coupling = {
        "land_area_fraction": land_area / total_area if total_area else 0.0,
        "continental_crust_area_fraction": cont_area / total_area if total_area else 0.0,
        "exposed_continental_land_fraction_of_land": (
            float(area[exposed_continental].sum()) / land_area if land_area else 0.0
        ),
        "exposed_oceanic_land_fraction_of_land": (
            float(area[exposed_oceanic].sum()) / land_area if land_area else 0.0
        ),
        "high_oceanic_land_fraction_gt500m_of_land": (
            float(area[high_oceanic_500].sum()) / land_area if land_area else 0.0
        ),
        "high_oceanic_land_fraction_gt1500m_of_land": (
            float(area[high_oceanic_1500].sum()) / land_area if land_area else 0.0
        ),
        "high_oceanic_land_cells_gt500m": int(high_oceanic_500.sum()),
        "high_oceanic_land_cells_gt1500m": int(high_oceanic_1500.sum()),
    }

    fields = {
        "tectonics.land_component_id": land_diag.labels,
        "tectonics.continental_component_id": cont_diag.labels,
        "tectonics.exposed_continental_component_id": exposed_cont_diag.labels,
        "tectonics.land_width_steps": land_diag.width_steps,
        "tectonics.continent_width_steps": cont_diag.width_steps,
        "tectonics.exposed_continental_width_steps": exposed_cont_diag.width_steps,
        "tectonics.ribbon_index": land_diag.ribbon_index,
        "tectonics.continental_ribbon_index": cont_diag.ribbon_index,
        "tectonics.coastline_complexity": land_diag.complexity_index,
        "tectonics.continental_edge_complexity": cont_diag.complexity_index,
        "tectonics.narrow_neck_index": land_diag.narrow_neck_index,
        "tectonics.continental_narrow_neck_index": cont_diag.narrow_neck_index,
    }
    objects = {
        "tectonics.land_components": land_diag.components,
        "tectonics.continental_components": cont_diag.components,
        "tectonics.exposed_continental_components": exposed_cont_diag.components,
    }
    metrics = {
        "context": {
            "spec_name": world.spec.name,
            "time_myr": float(world.time_myr),
            "n_cells": int(grid.n),
        },
        "exposed_land": land_diag.metrics,
        "continental_crust": cont_diag.metrics,
        "exposed_continental_land": exposed_cont_diag.metrics,
        "crust_land_coupling": coupling,
    }
    metrics["warnings"] = morphology_warnings(metrics)
    return WorldMorphology(fields=fields, objects=objects, metrics=metrics)


def ensure_morphology_fields(world) -> WorldMorphology:
    """Compute diagnostics and attach fields/objects to the world for rendering."""
    result = compute_world_morphology(world)
    for name, values in result.fields.items():
        world.fields[name] = np.asarray(values, dtype=np.float64)
    for name, objects in result.objects.items():
        world.objects[name] = objects
    return result


def morphology_warnings(metrics: dict[str, Any]) -> list[str]:
    warnings: list[str] = []
    context = metrics.get("context", {})
    earthlike = str(context.get("spec_name", "")).startswith("earthlike")
    land = metrics["exposed_land"]
    cont = metrics["continental_crust"]
    coupling = metrics["crust_land_coupling"]

    if earthlike and land["component_count"] > 14:
        warnings.append("Earthlike exposed land is split into many components")
    if earthlike and land["largest_component_area_fraction_of_mask"] < 0.35:
        warnings.append("Earthlike lacks a dominant large landmass")
    if earthlike and land["ribbon_area_fraction_gt_0_5"] > 0.08:
        warnings.append("exposed land contains a large ribbon-like fraction")
    if earthlike and cont["ribbon_area_fraction_gt_0_5"] > 0.08:
        warnings.append("continental crust contains a large ribbon-like fraction")
    if earthlike and land["narrow_neck_cells_per_1000_mask_cells"] > 15.0:
        warnings.append("exposed land has many graph articulation necks")
    if earthlike and cont["narrow_fraction_le2_width"] > 0.42:
        warnings.append("continental crust has too little broad interior area")
    if earthlike and coupling["high_oceanic_land_fraction_gt1500m_of_land"] > 0.005:
        warnings.append("high exposed oceanic crust is common")
    if earthlike and land["coastline_complexity_largest_component"] > 8.0:
        warnings.append("largest landmass coastline complexity is high")
    return warnings


def _component_labels(grid, mask: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    labels = np.full(grid.n, -1, dtype=np.int64)
    components: list[np.ndarray] = []
    for start in np.where(mask)[0]:
        if labels[start] >= 0:
            continue
        cid = len(components)
        stack = [int(start)]
        labels[start] = cid
        cells: list[int] = []
        while stack:
            c = stack.pop()
            cells.append(c)
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if mask[nb] and labels[nb] < 0:
                    labels[nb] = cid
                    stack.append(nb)
        components.append(np.asarray(cells, dtype=np.int64))
    return labels, components


def _width_steps(grid, mask: np.ndarray) -> np.ndarray:
    width = np.zeros(grid.n, dtype=np.float64)
    if not mask.any():
        return width
    outside_degree = np.zeros(grid.n, dtype=np.int64)
    for c in np.where(mask)[0]:
        outside_degree[c] = int((~mask[grid.neighbors[int(c)]]).sum())
    boundary = mask & (outside_degree > 0)
    if not boundary.any():
        width[mask] = math.sqrt(float(mask.sum()))
        return width

    q: deque[int] = deque()
    seen = np.zeros(grid.n, dtype=bool)
    for c in np.where(boundary)[0]:
        c = int(c)
        width[c] = 1.0
        seen[c] = True
        q.append(c)
    while q:
        c = q.popleft()
        for nb in grid.neighbors[c]:
            nb = int(nb)
            if not mask[nb] or seen[nb]:
                continue
            width[nb] = width[c] + 1.0
            seen[nb] = True
            q.append(nb)
    return width


def _articulation_points(grid, mask: np.ndarray) -> np.ndarray:
    art = np.zeros(grid.n, dtype=bool)
    if int(mask.sum()) < 3:
        return art

    old_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old_limit, grid.n + 1000))
    disc = np.full(grid.n, -1, dtype=np.int64)
    low = np.full(grid.n, -1, dtype=np.int64)
    parent = np.full(grid.n, -1, dtype=np.int64)
    time = 0

    def dfs(u: int) -> None:
        nonlocal time
        disc[u] = time
        low[u] = time
        time += 1
        children = 0
        for v in grid.neighbors[u]:
            v = int(v)
            if not mask[v]:
                continue
            if disc[v] < 0:
                parent[v] = u
                children += 1
                dfs(v)
                low[u] = min(low[u], low[v])
                if parent[u] < 0 and children > 1:
                    art[u] = True
                if parent[u] >= 0 and low[v] >= disc[u]:
                    art[u] = True
            elif v != parent[u]:
                low[u] = min(low[u], disc[v])

    try:
        for start in np.where(mask)[0]:
            start = int(start)
            if disc[start] < 0:
                dfs(start)
    finally:
        sys.setrecursionlimit(old_limit)
    return art


def _component_summaries(
    grid,
    mask: np.ndarray,
    labels: np.ndarray,
    cells_by_component: list[np.ndarray],
    width: np.ndarray,
    articulation: np.ndarray,
    *,
    name: str,
) -> tuple[list[dict[str, Any]], np.ndarray, np.ndarray]:
    ncomp = len(cells_by_component)
    area = grid.cell_area
    total_area = float(area.sum())
    mask_area = float(area[mask].sum())
    perim_edges = np.zeros(ncomp, dtype=np.float64)
    perim_length = np.zeros(ncomp, dtype=np.float64)

    if ncomp:
        for edge_idx, (i, j) in enumerate(grid.edges):
            i = int(i)
            j = int(j)
            if mask[i] == mask[j]:
                continue
            edge_len = float(grid.edge_lengths[edge_idx])
            if mask[i]:
                cid = int(labels[i])
                perim_edges[cid] += 1.0
                perim_length[cid] += edge_len
            if mask[j]:
                cid = int(labels[j])
                perim_edges[cid] += 1.0
                perim_length[cid] += edge_len

    components: list[dict[str, Any]] = []
    complexity_field = np.zeros(grid.n, dtype=np.float64)
    ribbon_field = np.zeros(grid.n, dtype=np.float64)

    for cid, cells in enumerate(cells_by_component):
        cells = np.asarray(cells, dtype=np.int64)
        cell_count = int(cells.size)
        comp_area = float(area[cells].sum())
        p_edges = float(perim_edges[cid])
        complexity = p_edges / max(math.sqrt(float(cell_count)), 1.0)
        compactness = (
            min(1.0, 4.0 * math.pi * float(cell_count) / max(p_edges * p_edges, 1.0))
            if cell_count else 0.0
        )
        local_width = width[cells]
        weights = area[cells]
        narrow1 = _weighted_fraction(local_width <= 1.0, weights)
        narrow2 = _weighted_fraction(local_width <= 2.0, weights)
        narrow3 = _weighted_fraction(local_width <= 3.0, weights)
        art = articulation[cells]
        neck_area_fraction = _weighted_fraction(art & (local_width <= 2.0), weights)
        shape_score = float(np.clip((complexity - 4.6) / 4.0, 0.0, 1.0))
        size_score = float(np.clip((cell_count - 3.0) / 12.0, 0.0, 1.0))
        component_ribbon_score = shape_score * size_score * (0.55 * narrow2 + 0.45 * narrow3)
        local_narrow_score = np.clip((3.25 - local_width) / 2.25, 0.0, 1.0)
        ribbon_field[cells] = local_narrow_score * shape_score * size_score
        complexity_field[cells] = np.clip((complexity - 3.0) / 6.0, 0.0, 1.0)
        components.append({
            "id": int(cid),
            "name": name,
            "cell_count": cell_count,
            "area_m2": comp_area,
            "area_fraction_of_total": comp_area / total_area if total_area else 0.0,
            "area_fraction_of_mask": comp_area / mask_area if mask_area else 0.0,
            "perimeter_edges": p_edges,
            "perimeter_length_m": float(perim_length[cid]),
            "coastline_complexity": complexity,
            "compactness_proxy": compactness,
            "width_p10_steps": _percentile(local_width, 10),
            "width_p50_steps": _percentile(local_width, 50),
            "width_p90_steps": _percentile(local_width, 90),
            "width_max_steps": float(local_width.max()) if local_width.size else 0.0,
            "narrow_fraction_le1_width": narrow1,
            "narrow_fraction_le2_width": narrow2,
            "narrow_fraction_le3_width": narrow3,
            "narrow_neck_area_fraction": neck_area_fraction,
            "narrow_neck_cells": int((art & (local_width <= 2.0)).sum()),
            "ribbon_score": component_ribbon_score,
        })

    return components, complexity_field, ribbon_field


def _mask_metrics(
    grid,
    mask: np.ndarray,
    width: np.ndarray,
    articulation: np.ndarray,
    components: list[dict[str, Any]],
    ribbon_field: np.ndarray,
    name: str,
) -> dict[str, Any]:
    area = grid.cell_area
    total_area = float(area.sum())
    mask_area = float(area[mask].sum())
    cells = np.where(mask)[0]
    component_areas = [float(c["area_m2"]) for c in components]
    largest = max(component_areas, default=0.0)
    if cells.size:
        cell_weights = area[cells]
        width_values = width[cells]
        ribbon_values = ribbon_field[cells]
        ribbon_gt = ribbon_values >= 0.5
        narrow_neck = articulation[cells] & (width_values <= 2.0)
    else:
        cell_weights = np.array([], dtype=np.float64)
        width_values = np.array([], dtype=np.float64)
        ribbon_gt = np.array([], dtype=bool)
        narrow_neck = np.array([], dtype=bool)

    small_component_area_fraction = 0.0
    small_component_count = 0
    if mask_area > 0.0:
        for comp in components:
            if comp["area_fraction_of_mask"] <= 0.005:
                small_component_area_fraction += float(comp["area_fraction_of_mask"])
                small_component_count += 1

    largest_component = None
    if components:
        largest_component = max(components, key=lambda c: c["area_m2"])

    return {
        "name": name,
        "component_count": int(len(components)),
        "mask_cell_count": int(mask.sum()),
        "mask_area_fraction_of_total": mask_area / total_area if total_area else 0.0,
        "largest_component_area_fraction_of_mask": largest / mask_area if mask_area else 0.0,
        "small_component_count_le_0_5pct_mask": int(small_component_count),
        "small_component_area_fraction_le_0_5pct_mask": small_component_area_fraction,
        "width_p10_steps": _percentile(width_values, 10),
        "width_p50_steps": _percentile(width_values, 50),
        "width_p90_steps": _percentile(width_values, 90),
        "width_max_steps": float(width_values.max()) if width_values.size else 0.0,
        "narrow_fraction_le1_width": _weighted_fraction(width_values <= 1.0, cell_weights),
        "narrow_fraction_le2_width": _weighted_fraction(width_values <= 2.0, cell_weights),
        "narrow_fraction_le3_width": _weighted_fraction(width_values <= 3.0, cell_weights),
        "ribbon_area_fraction_gt_0_5": _weighted_fraction(ribbon_gt, cell_weights),
        "ribbon_index_area_mean": _weighted_mean(ribbon_values, cell_weights),
        "ribbon_index_p95": _percentile(ribbon_values, 95),
        "narrow_neck_cells": int(narrow_neck.sum()),
        "narrow_neck_cells_per_1000_mask_cells": (
            float(narrow_neck.sum()) * 1000.0 / float(cells.size) if cells.size else 0.0
        ),
        "coastline_complexity_p50_component": _component_percentile(
            components, "coastline_complexity", 50),
        "coastline_complexity_p95_component": _component_percentile(
            components, "coastline_complexity", 95),
        "coastline_complexity_largest_component": (
            float(largest_component["coastline_complexity"]) if largest_component else 0.0
        ),
        "compactness_p50_component": _component_percentile(
            components, "compactness_proxy", 50),
        "compactness_largest_component": (
            float(largest_component["compactness_proxy"]) if largest_component else 0.0
        ),
    }


def _weighted_fraction(mask: np.ndarray, weights: np.ndarray) -> float:
    mask = np.asarray(mask, dtype=bool)
    weights = np.asarray(weights, dtype=np.float64)
    if mask.size == 0 or weights.size == 0:
        return 0.0
    denom = float(weights.sum())
    if denom <= 0.0:
        return 0.0
    return float(weights[mask].sum() / denom)


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if values.size == 0 or weights.size == 0:
        return 0.0
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not valid.any():
        return 0.0
    return float(np.average(values[valid], weights=weights[valid]))


def _percentile(values: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0
    return float(np.percentile(values, q))


def _component_percentile(components: list[dict[str, Any]], key: str, q: float) -> float:
    if not components:
        return 0.0
    values = np.asarray([float(c[key]) for c in components], dtype=np.float64)
    return _percentile(values, q)
