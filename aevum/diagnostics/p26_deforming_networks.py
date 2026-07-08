"""P26 diagnostics for explicit deforming-network state.

The P26 production changes should consume active deformation as its own state,
not infer active tectonics from `crust.origin` or `crust.reworked_age_myr`.
This read-only summary records whether the new deformation layer is localized,
which styles dominate it, and how strongly it overlaps ribbon-like land or
continental crust.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from aevum.diagnostics.morphology import compute_world_morphology
from aevum.modules.tectonics import (
    CONT,
    DEFORM_COLLISION_CORE,
    DEFORM_COLLISION_SHOULDER,
    DEFORM_NONE,
    DEFORM_RIFT,
    DEFORM_SUBDUCTION_CORE,
    DEFORM_SUBDUCTION_SHOULDER,
    DEFORM_TRANSFORM,
)


SCHEMA = "aevum.p26_deforming_networks.v1"

STYLE_NAMES = {
    int(DEFORM_NONE): "none",
    int(DEFORM_COLLISION_CORE): "collision_core",
    int(DEFORM_COLLISION_SHOULDER): "collision_shoulder",
    int(DEFORM_SUBDUCTION_CORE): "subduction_core",
    int(DEFORM_SUBDUCTION_SHOULDER): "subduction_shoulder",
    int(DEFORM_RIFT): "rift",
    int(DEFORM_TRANSFORM): "transform",
}

CORE_CODES = {
    int(DEFORM_COLLISION_CORE),
    int(DEFORM_SUBDUCTION_CORE),
}
SHOULDER_CODES = {
    int(DEFORM_COLLISION_SHOULDER),
    int(DEFORM_SUBDUCTION_SHOULDER),
}


def deforming_network_summary(
    world,
    *,
    ribbon_threshold: float = 0.5,
    top_n: int = 8,
) -> dict[str, Any]:
    """Summarize active deformation style, footprint, and ribbon overlap."""
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total_area = float(area.sum())
    intensity = _field(world, "tectonics.deformation_intensity")
    style = _field(world, "tectonics.deformation_style").astype(int)
    active = (intensity > 0.0) | (style != int(DEFORM_NONE))

    crust = _field(world, "crust.type").astype(int)
    continental = crust == int(CONT)
    rel = _field(world, "terrain.elevation_m") - float(world.sea_level)
    land = rel >= 0.0
    morphology = compute_world_morphology(world)
    land_ribbon = morphology.fields["tectonics.ribbon_index"] >= float(ribbon_threshold)
    continental_ribbon = (
        morphology.fields["tectonics.continental_ribbon_index"] >= float(ribbon_threshold)
    )

    active_area = _area(area, active)
    continental_area = _area(area, continental)
    land_area = _area(area, land)
    land_ribbon_area = _area(area, land_ribbon)
    continental_ribbon_area = _area(area, continental_ribbon)
    core = active & np.isin(style, list(CORE_CODES))
    shoulder = active & np.isin(style, list(SHOULDER_CODES))
    core_area = _area(area, core)
    shoulder_area = _area(area, shoulder)

    style_summaries = {
        name: _style_summary(
            area,
            intensity,
            style,
            int(code),
            active_area=active_area,
            total_area=total_area,
            continental=continental,
            land=land,
            land_ribbon=land_ribbon,
            continental_ribbon=continental_ribbon,
        )
        for code, name in STYLE_NAMES.items()
        if int(code) != int(DEFORM_NONE)
    }

    objects = _object_summary(world, top_n=top_n)
    land_ribbon_overlap = active & land_ribbon
    continental_ribbon_overlap = active & continental_ribbon
    return {
        "schema": SCHEMA,
        "context": {
            "spec_name": str(world.spec.name),
            "seed": int(world.spec.seed),
            "cells": int(grid.n),
            "time_myr": float(world.time_myr),
        },
        "parameters": {
            "ribbon_threshold": float(ribbon_threshold),
            "top_n": int(top_n),
        },
        "metrics": {
            "active_deformation_area_fraction_of_world": (
                active_area / total_area if total_area else 0.0
            ),
            "active_deformation_area_fraction_of_continental": (
                _area(area, active & continental) / continental_area
                if continental_area else 0.0
            ),
            "active_deformation_area_fraction_of_land": (
                _area(area, active & land) / land_area if land_area else 0.0
            ),
            "mean_active_intensity": _weighted_mean(area, intensity, active),
            "p95_active_intensity": _percentile(intensity[active], 95.0),
            "core_area_fraction_of_active": (
                core_area / active_area if active_area else 0.0
            ),
            "shoulder_area_fraction_of_active": (
                shoulder_area / active_area if active_area else 0.0
            ),
            "core_to_shoulder_area_ratio": (
                core_area / shoulder_area if shoulder_area else 0.0
            ),
            "land_ribbon_deformation_coverage_fraction": (
                _area(area, land_ribbon_overlap) / land_ribbon_area
                if land_ribbon_area else 0.0
            ),
            "continental_ribbon_deformation_coverage_fraction": (
                _area(area, continental_ribbon_overlap) / continental_ribbon_area
                if continental_ribbon_area else 0.0
            ),
            "active_deformation_overlap_fraction_with_land_ribbon": (
                _area(area, land_ribbon_overlap) / active_area
                if active_area else 0.0
            ),
            "active_deformation_overlap_fraction_with_continental_ribbon": (
                _area(area, continental_ribbon_overlap) / active_area
                if active_area else 0.0
            ),
            "deforming_network_object_count": int(objects["object_count"]),
            "largest_deforming_network_area_fraction": float(
                objects["largest_area_fraction"]),
        },
        "style_summaries": style_summaries,
        "object_kind_counts": objects["kind_counts"],
        "top_deforming_network_objects": objects["top_objects"],
        "diagnostic_hints": _diagnostic_hints(
            active_area=active_area,
            total_area=total_area,
            core_area=core_area,
            shoulder_area=shoulder_area,
            land_ribbon_overlap=_area(area, land_ribbon_overlap),
            land_ribbon_area=land_ribbon_area,
            continental_ribbon_overlap=_area(area, continental_ribbon_overlap),
            continental_ribbon_area=continental_ribbon_area,
        ),
    }


def _field(world, name: str) -> np.ndarray:
    if name in world.fields:
        return np.asarray(world.fields[name], dtype=np.float64)
    return np.zeros(world.grid.n, dtype=np.float64)


def _style_summary(
    area: np.ndarray,
    intensity: np.ndarray,
    style: np.ndarray,
    code: int,
    *,
    active_area: float,
    total_area: float,
    continental: np.ndarray,
    land: np.ndarray,
    land_ribbon: np.ndarray,
    continental_ribbon: np.ndarray,
) -> dict[str, Any]:
    mask = style == int(code)
    style_area = _area(area, mask)
    return {
        "cell_count": int(mask.sum()),
        "area_fraction_of_world": style_area / total_area if total_area else 0.0,
        "area_fraction_of_active": style_area / active_area if active_area else 0.0,
        "mean_intensity": _weighted_mean(area, intensity, mask),
        "continental_fraction": _mask_fraction(mask & continental, mask),
        "land_fraction": _mask_fraction(mask & land, mask),
        "land_ribbon_overlap_fraction": _mask_fraction(mask & land_ribbon, mask),
        "continental_ribbon_overlap_fraction": _mask_fraction(
            mask & continental_ribbon, mask),
    }


def _object_summary(world, *, top_n: int) -> dict[str, Any]:
    objects = list(world.objects.get("tectonics.deforming_networks", []))
    kind_counts: dict[str, int] = {}
    for obj in objects:
        kind = str(obj.get("kind", "unknown"))
        kind_counts[kind] = kind_counts.get(kind, 0) + 1
    top = sorted(
        objects,
        key=lambda obj: (
            -float(obj.get("area_fraction", 0.0)),
            str(obj.get("id", "")),
        ),
    )[:max(0, int(top_n))]
    return {
        "object_count": int(len(objects)),
        "largest_area_fraction": max(
            [float(obj.get("area_fraction", 0.0)) for obj in objects] or [0.0]
        ),
        "kind_counts": kind_counts,
        "top_objects": [
            {
                "id": str(obj.get("id", "")),
                "kind": str(obj.get("kind", "")),
                "style_code": int(obj.get("style_code", 0)),
                "cell_count": int(obj.get("cell_count", 0)),
                "area_fraction": float(obj.get("area_fraction", 0.0)),
                "mean_intensity": float(obj.get("mean_intensity", 0.0)),
                "continental_fraction": float(obj.get("continental_fraction", 0.0)),
                "lat": float(obj.get("lat", 0.0)),
                "lon": float(obj.get("lon", 0.0)),
                "last_active_myr": float(obj.get("last_active_myr", 0.0)),
            }
            for obj in top
        ],
    }


def _diagnostic_hints(
    *,
    active_area: float,
    total_area: float,
    core_area: float,
    shoulder_area: float,
    land_ribbon_overlap: float,
    land_ribbon_area: float,
    continental_ribbon_overlap: float,
    continental_ribbon_area: float,
) -> dict[str, Any]:
    active_fraction = active_area / total_area if total_area else 0.0
    shoulder_share = shoulder_area / active_area if active_area else 0.0
    land_ribbon_coverage = (
        land_ribbon_overlap / land_ribbon_area if land_ribbon_area else 0.0
    )
    continental_ribbon_coverage = (
        continental_ribbon_overlap / continental_ribbon_area
        if continental_ribbon_area else 0.0
    )
    return {
        "active_network_is_world_broad": bool(active_fraction > 0.18),
        "shoulder_dominated_active_network": bool(
            active_area > 0.0 and shoulder_share > 0.65),
        "land_ribbon_is_deformation_coupled": bool(land_ribbon_coverage > 0.35),
        "continental_ribbon_is_deformation_coupled": bool(
            continental_ribbon_coverage > 0.35),
        "core_without_shoulder": bool(core_area > 0.0 and shoulder_area <= 0.0),
    }


def _area(area: np.ndarray, mask: np.ndarray) -> float:
    mask = np.asarray(mask, dtype=bool)
    return float(area[mask].sum()) if mask.any() else 0.0


def _mask_fraction(mask: np.ndarray, denom_mask: np.ndarray) -> float:
    denom = int(np.asarray(denom_mask, dtype=bool).sum())
    if denom == 0:
        return 0.0
    return float(np.asarray(mask, dtype=bool).sum() / denom)


def _weighted_mean(area: np.ndarray, values: np.ndarray, mask: np.ndarray) -> float:
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return 0.0
    return float(np.average(values[mask], weights=area[mask]))


def _percentile(values: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return 0.0
    return float(np.percentile(values, q))
