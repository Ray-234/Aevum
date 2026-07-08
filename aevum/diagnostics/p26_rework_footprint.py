"""P26 diagnostics for recent tectonic rework footprint width.

The current P26 ribbon attribution shows that Earth-like ribbon area is mostly
recently reworked suture/accretionary crust.  This helper measures whether that
recent rework is spatially concentrated near explicit boundary corridors, or
whether it spreads through broad continental interiors.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from aevum.modules.tectonics import (
    CONT,
    DOMAIN_ACCRETED_TERRANE,
    DOMAIN_SUTURE,
    ORIGIN_ARC,
    ORIGIN_SUTURE,
)


SCHEMA = "aevum.p26_rework_footprint.v1"


def recent_rework_footprint_summary(
    world,
    *,
    boundary_mask: np.ndarray | None = None,
    recent_myr: float = 220.0,
    corridor_passes: int = 2,
) -> dict[str, Any]:
    """Measure whether recent arc/suture rework is localized near boundaries."""
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    crust = world.get_field("crust.type", 0.0).astype(int)
    domain = world.get_field("crust.domain", 0.0).astype(int)
    origin = world.get_field("crust.origin", 0.0).astype(int)
    reworked = np.asarray(world.get_field("crust.reworked_age_myr", -1.0), dtype=np.float64)
    orogeny = np.asarray(world.get_field("tectonics.orogeny_age_myr", -1.0), dtype=np.float64)
    volcanism = np.asarray(
        world.get_field("tectonics.volcanism_age_myr", -1.0), dtype=np.float64)

    cont = crust == int(CONT)
    continental_area = float(area[cont].sum())
    if boundary_mask is None:
        boundary = _boundary_mask_from_world(world)
    else:
        boundary = np.asarray(boundary_mask, dtype=bool)
        if boundary.shape != (grid.n,):
            boundary = np.zeros(grid.n, dtype=bool)
    corridor = _dilate_mask(grid, boundary, passes=int(corridor_passes)) & cont

    recent_reworked = _time_since_event(reworked, world.time_myr) < float(recent_myr)
    recent_orogeny = _time_since_event(orogeny, world.time_myr) < float(recent_myr)
    recent_volcanism = _time_since_event(volcanism, world.time_myr) < float(recent_myr)
    active_time = recent_reworked | recent_orogeny | recent_volcanism
    active_kind = (
        np.isin(origin, [int(ORIGIN_ARC), int(ORIGIN_SUTURE)])
        | np.isin(domain, [int(DOMAIN_ACCRETED_TERRANE), int(DOMAIN_SUTURE)])
    )
    active = cont & active_time & active_kind

    active_area = float(area[active].sum())
    corridor_area = float(area[corridor].sum())
    outside = active & ~corridor
    inside = active & corridor
    outside_area = float(area[outside].sum())
    inside_area = float(area[inside].sum())
    boundary_area = float(area[boundary].sum())

    active_fraction_of_cont = (
        active_area / continental_area if continental_area > 0.0 else 0.0
    )
    outside_fraction_of_active = (
        outside_area / active_area if active_area > 0.0 else 0.0
    )
    inside_fraction_of_active = (
        inside_area / active_area if active_area > 0.0 else 0.0
    )
    inside_fraction_of_corridor = (
        inside_area / corridor_area if corridor_area > 0.0 else 0.0
    )
    active_to_corridor_area_ratio = (
        active_area / corridor_area if corridor_area > 0.0 else 0.0
    )
    overbroad = (
        active_fraction_of_cont > 0.04
        and outside_fraction_of_active > 0.35
        and active_to_corridor_area_ratio > 1.45
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
            "recent_myr": float(recent_myr),
            "corridor_passes": int(corridor_passes),
        },
        "metrics": {
            "continental_area_fraction_of_world": (
                continental_area / float(area.sum()) if float(area.sum()) > 0.0 else 0.0
            ),
            "boundary_area_fraction_of_world": (
                boundary_area / float(area.sum()) if float(area.sum()) > 0.0 else 0.0
            ),
            "boundary_corridor_area_fraction_of_continental": (
                corridor_area / continental_area if continental_area > 0.0 else 0.0
            ),
            "active_rework_area_fraction_of_continental": active_fraction_of_cont,
            "active_rework_outside_corridor_fraction_of_active": outside_fraction_of_active,
            "active_rework_inside_corridor_fraction_of_active": inside_fraction_of_active,
            "active_rework_inside_corridor_fraction_of_corridor": inside_fraction_of_corridor,
            "active_to_corridor_area_ratio": active_to_corridor_area_ratio,
            "recent_reworked_fraction_of_active": _fraction(recent_reworked & active, active),
            "recent_orogeny_fraction_of_active": _fraction(recent_orogeny & active, active),
            "recent_volcanism_fraction_of_active": _fraction(recent_volcanism & active, active),
            "active_rework_cell_count": int(active.sum()),
            "boundary_corridor_cell_count": int(corridor.sum()),
            "overbroad_recent_rework": bool(overbroad),
        },
    }


def _boundary_mask_from_world(world) -> np.ndarray:
    grid = world.grid
    out = np.zeros(grid.n, dtype=bool)
    boundaries = world.networks.get("tectonics.boundaries", {})
    if not isinstance(boundaries, dict):
        return out
    for key in ("collision", "suture", "subduction", "trench", "active_margin"):
        cells = boundaries.get(key, [])
        if len(cells) == 0:
            continue
        idx = np.asarray(cells, dtype=int)
        idx = idx[(0 <= idx) & (idx < grid.n)]
        out[idx] = True
    return out


def _time_since_event(event_age_myr: np.ndarray, time_myr: float) -> np.ndarray:
    event_age_myr = np.asarray(event_age_myr, dtype=np.float64)
    out = np.full(event_age_myr.shape, float(time_myr) + 1.0, dtype=np.float64)
    valid = np.isfinite(event_age_myr) & (event_age_myr >= 0.0)
    out[valid] = np.maximum(float(time_myr) - event_age_myr[valid], 0.0)
    return out


def _dilate_mask(grid, mask: np.ndarray, *, passes: int = 1) -> np.ndarray:
    out = np.asarray(mask, dtype=bool).copy()
    for _ in range(max(0, int(passes))):
        add = np.zeros(grid.n, dtype=bool)
        for c in np.where(out)[0]:
            add[grid.neighbors[int(c)]] = True
        out |= add
    return out


def _fraction(mask: np.ndarray, denom_mask: np.ndarray) -> float:
    denom = int(np.asarray(denom_mask, dtype=bool).sum())
    if denom == 0:
        return 0.0
    return float(np.asarray(mask, dtype=bool).sum() / denom)
