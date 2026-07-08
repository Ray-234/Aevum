"""P29 diagnostics for inland geomorphology complexity.

This read-only summary separates two ideas that are easy to conflate:
large low-elevation continental interiors are Earth-like, but a broad
single-elevation platform with few landform objects is underexpressed.
"""
from __future__ import annotations

from collections import Counter
from typing import Any

import numpy as np

from aevum.modules.tectonics import CONT
from aevum.modules.terrain import (
    CONT_DETAIL_ARC_MICROCONTINENT,
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATFORM,
    CONT_DETAIL_PLATEAU,
    CONT_DETAIL_RIFT_BASIN,
    CONT_DETAIL_SHIELD,
)


SCHEMA = "aevum.p29_inland_geomorphology.v1"

DETAIL_NAMES = {
    int(CONT_DETAIL_SHIELD): "shield",
    int(CONT_DETAIL_PLATFORM): "platform",
    int(CONT_DETAIL_BASIN): "basin",
    int(CONT_DETAIL_RIFT_BASIN): "rift_basin",
    int(CONT_DETAIL_OROGEN): "orogen",
    int(CONT_DETAIL_PLATEAU): "plateau",
    int(CONT_DETAIL_ARC_MICROCONTINENT): "arc_microcontinent",
}

INLAND_OBJECT_KINDS = (
    "shield",
    "platform",
    "interior_basin",
    "old_subdued_orogen",
    "rift_basin",
    "plateau",
    "orogen",
    "foreland_basin",
)


def inland_geomorphology_summary(
    world,
    *,
    inland_width_steps: float = 3.0,
    lowland_min_m: float = 100.0,
    lowland_max_m: float = 1050.0,
    highland_min_m: float = 1600.0,
) -> dict[str, Any]:
    """Summarize broad continental interior relief and landform diversity."""
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total_area = float(area.sum())
    elev = _field(world, "terrain.elevation_m") - float(world.sea_level)
    land = elev >= 0.0
    continental = _field(world, "crust.type").astype(int) == int(CONT)
    continental_land = continental & land
    width = _mask_width_steps(grid, continental_land)
    inland = continental_land & (width >= float(inland_width_steps))
    if not inland.any():
        inland = continental_land.copy()

    land_area = _area(area, land)
    continental_land_area = _area(area, continental_land)
    inland_area = _area(area, inland)
    detail = _field(world, "terrain.continental_detail").astype(int)
    detail_fractions = _detail_fractions(area, detail, inland, inland_area)
    objects = _object_summary(world, area, total_area)
    rel = elev[inland]

    flat_lowland = (
        inland
        & (elev >= float(lowland_min_m))
        & (elev <= float(lowland_max_m))
    )
    highland = inland & (elev >= float(highland_min_m))
    monotone_flat = (
        inland_area > 0.0
        and _percentile(rel, 95.0) - _percentile(rel, 5.0) < 320.0
        and _area(area, flat_lowland) / inland_area > 0.86
        and len([v for v in detail_fractions.values() if v > 0.02]) <= 2
        and objects["inland_object_kind_count"] <= 2
    )
    low_object_diversity = (
        inland_area > 0.0
        and objects["inland_object_kind_count"] < 4
        and len([v for v in detail_fractions.values() if v > 0.02]) < 4
    )

    return {
        "schema": SCHEMA,
        "context": {
            "spec_name": str(world.spec.name),
            "seed": int(world.spec.seed),
            "cells": int(grid.n),
            "time_myr": float(world.time_myr),
        },
        "parameters": {
            "inland_width_steps": float(inland_width_steps),
            "lowland_min_m": float(lowland_min_m),
            "lowland_max_m": float(lowland_max_m),
            "highland_min_m": float(highland_min_m),
        },
        "metrics": {
            "continental_land_area_fraction_of_world": (
                continental_land_area / total_area if total_area else 0.0
            ),
            "inland_area_fraction_of_world": (
                inland_area / total_area if total_area else 0.0
            ),
            "inland_area_fraction_of_land": (
                inland_area / land_area if land_area else 0.0
            ),
            "inland_area_fraction_of_continental_land": (
                inland_area / continental_land_area if continental_land_area else 0.0
            ),
            "inland_cell_count": int(inland.sum()),
            "inland_width_p50_steps": _percentile(width[inland], 50.0),
            "inland_elevation_p05_m": _percentile(rel, 5.0),
            "inland_elevation_p25_m": _percentile(rel, 25.0),
            "inland_elevation_p50_m": _percentile(rel, 50.0),
            "inland_elevation_p75_m": _percentile(rel, 75.0),
            "inland_elevation_p95_m": _percentile(rel, 95.0),
            "inland_relief_p95_p05_m": _percentile(rel, 95.0) - _percentile(rel, 5.0),
            "inland_relief_p90_p10_m": _percentile(rel, 90.0) - _percentile(rel, 10.0),
            "inland_relief_iqr_m": _percentile(rel, 75.0) - _percentile(rel, 25.0),
            "inland_flat_lowland_fraction": (
                _area(area, flat_lowland) / inland_area if inland_area else 0.0
            ),
            "inland_highland_fraction": (
                _area(area, highland) / inland_area if inland_area else 0.0
            ),
            "inland_detail_diversity_gt2pct": int(
                sum(1 for v in detail_fractions.values() if v > 0.02)
            ),
            "inland_object_kind_count": int(objects["inland_object_kind_count"]),
            "inland_landform_object_count": int(objects["inland_landform_object_count"]),
        },
        "detail_fractions": detail_fractions,
        "object_kind_counts": objects["kind_counts"],
        "top_inland_landform_objects": objects["top_objects"],
        "diagnostic_hints": {
            "monotone_flat_inland": bool(monotone_flat),
            "low_inland_object_diversity": bool(low_object_diversity),
            "broad_lowland_mode_present": (
                inland_area > 0.0
                and _area(area, flat_lowland) / inland_area >= 0.35
            ),
            "highland_tail_present": (
                inland_area > 0.0
                and _area(area, highland) / inland_area >= 0.02
            ),
        },
    }


def _field(world, name: str) -> np.ndarray:
    if name in world.fields:
        return np.asarray(world.fields[name], dtype=np.float64)
    return np.zeros(world.grid.n, dtype=np.float64)


def _area(area: np.ndarray, mask: np.ndarray) -> float:
    mask = np.asarray(mask, dtype=bool)
    return float(area[mask].sum()) if mask.any() else 0.0


def _percentile(values: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return 0.0
    return float(np.percentile(values, q))


def _detail_fractions(
    area: np.ndarray, detail: np.ndarray, inland: np.ndarray, inland_area: float
) -> dict[str, float]:
    out = {name: 0.0 for name in DETAIL_NAMES.values()}
    if inland_area <= 0.0:
        return out
    for code, name in DETAIL_NAMES.items():
        out[name] = _area(area, inland & (detail == int(code))) / inland_area
    return out


def _object_summary(world, area: np.ndarray, total_area: float) -> dict[str, Any]:
    objects = list(world.objects.get("terrain.continental_landforms", []))
    kind_counts = Counter(str(obj.get("kind", "unknown")) for obj in objects)
    inland_objects = [
        obj for obj in objects
        if str(obj.get("kind", "")) in set(INLAND_OBJECT_KINDS)
    ]
    top = []
    for obj in sorted(
        inland_objects,
        key=lambda item: float(item.get("area_fraction", 0.0)),
        reverse=True,
    )[:10]:
        top.append({
            "id": str(obj.get("id", "")),
            "kind": str(obj.get("kind", "")),
            "area_fraction": float(obj.get("area_fraction", 0.0)),
            "cell_count": int(obj.get("cell_count", 0)),
            "mean_elevation_m": float(obj.get("mean_elevation_m", 0.0)),
            "mean_sediment_m": float(obj.get("mean_sediment_m", 0.0)),
            "parent_tectonic_object_count": int(
                len(obj.get("parent_tectonic_object_ids", []))
            ),
        })
    del area, total_area
    return {
        "kind_counts": dict(sorted(kind_counts.items())),
        "inland_landform_object_count": int(len(inland_objects)),
        "inland_object_kind_count": int(
            len({str(obj.get("kind", "")) for obj in inland_objects})
        ),
        "top_objects": top,
    }


def _mask_width_steps(grid, mask: np.ndarray) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    width = np.zeros(grid.n, dtype=np.float64)
    if not mask.any():
        return width
    boundary = np.zeros(grid.n, dtype=bool)
    for c in np.where(mask)[0]:
        if (~mask[grid.neighbors[int(c)]]).any():
            boundary[int(c)] = True
    if not boundary.any():
        width[mask] = np.sqrt(float(mask.sum()))
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
            width[nb] = width[c] + 1.0
            seen[nb] = True
            queue.append(nb)
    return width
