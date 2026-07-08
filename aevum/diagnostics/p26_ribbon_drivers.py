"""P26 ribbon-driver attribution diagnostics.

P26 continent/margin geometry work needs more than a world-level fail/pass gate:
when ribbon fraction regresses, the summary should identify which components
and crust histories are responsible.  This module is read-only.  It inspects the
current world morphology and reports the top exposed-land and continental-crust
components contributing to ribbon-like area.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from aevum.diagnostics.morphology import compute_world_morphology
from aevum.modules.tectonics import (
    DOMAIN_ACCRETED_TERRANE,
    DOMAIN_CONTINENTAL_INTERIOR,
    DOMAIN_CONTINENTAL_MARGIN,
    DOMAIN_CRATON,
    DOMAIN_LIP,
    DOMAIN_OCEANIC,
    DOMAIN_SUTURE,
    ORIGIN_ARC,
    ORIGIN_CRATON,
    ORIGIN_PLUME_IMPACT,
    ORIGIN_PRIMORDIAL,
    ORIGIN_RIDGE,
    ORIGIN_SUTURE,
)
from aevum.modules.terrain import (
    CONT_DETAIL_ARC_MICROCONTINENT,
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OCEAN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATFORM,
    CONT_DETAIL_PLATEAU,
    CONT_DETAIL_RIFT_BASIN,
    CONT_DETAIL_SHIELD,
)


SCHEMA = "aevum.p26_ribbon_drivers.v1"

DOMAIN_NAMES = {
    int(DOMAIN_OCEANIC): "oceanic",
    int(DOMAIN_CRATON): "craton",
    int(DOMAIN_CONTINENTAL_INTERIOR): "continental_interior",
    int(DOMAIN_CONTINENTAL_MARGIN): "continental_margin",
    int(DOMAIN_ACCRETED_TERRANE): "accreted_terrane",
    int(DOMAIN_SUTURE): "suture",
    int(DOMAIN_LIP): "lip",
}

ORIGIN_NAMES = {
    int(ORIGIN_RIDGE): "ridge",
    int(ORIGIN_PRIMORDIAL): "primordial",
    int(ORIGIN_ARC): "arc",
    int(ORIGIN_SUTURE): "suture",
    int(ORIGIN_PLUME_IMPACT): "plume_impact",
    int(ORIGIN_CRATON): "craton",
}

DETAIL_NAMES = {
    int(CONT_DETAIL_OCEAN): "ocean",
    int(CONT_DETAIL_SHIELD): "shield",
    int(CONT_DETAIL_PLATFORM): "platform",
    int(CONT_DETAIL_BASIN): "basin",
    int(CONT_DETAIL_RIFT_BASIN): "rift_basin",
    int(CONT_DETAIL_OROGEN): "orogen",
    int(CONT_DETAIL_PLATEAU): "plateau",
    int(CONT_DETAIL_ARC_MICROCONTINENT): "arc_microcontinent",
}


def ribbon_driver_summary(
    world,
    *,
    top_n: int = 6,
    ribbon_threshold: float = 0.5,
) -> dict[str, Any]:
    """Return component-level attribution for exposed-land/continent ribbons."""
    morphology = compute_world_morphology(world)
    area = np.asarray(world.grid.cell_area, dtype=np.float64)
    total_area = float(area.sum())

    land_labels = morphology.fields["tectonics.land_component_id"].astype(int)
    cont_labels = morphology.fields["tectonics.continental_component_id"].astype(int)
    exposed_cont_labels = morphology.fields[
        "tectonics.exposed_continental_component_id"].astype(int)
    land_ribbon = morphology.fields["tectonics.ribbon_index"]
    cont_ribbon = morphology.fields["tectonics.continental_ribbon_index"]

    land_components = _top_components(
        world,
        morphology.objects["tectonics.land_components"],
        land_labels,
        land_ribbon,
        "exposed_land",
        top_n=top_n,
        ribbon_threshold=ribbon_threshold,
    )
    cont_components = _top_components(
        world,
        morphology.objects["tectonics.continental_components"],
        cont_labels,
        cont_ribbon,
        "continental_crust",
        top_n=top_n,
        ribbon_threshold=ribbon_threshold,
    )
    exposed_cont_components = _top_components(
        world,
        morphology.objects["tectonics.exposed_continental_components"],
        exposed_cont_labels,
        land_ribbon,
        "exposed_continental_land",
        top_n=top_n,
        ribbon_threshold=ribbon_threshold,
    )

    land_ribbon_area = _ribbon_area(area, land_labels, land_ribbon, ribbon_threshold)
    cont_ribbon_area = _ribbon_area(area, cont_labels, cont_ribbon, ribbon_threshold)
    land_temporal = _ribbon_temporal_summary(
        world, land_labels, land_ribbon, ribbon_threshold)
    cont_temporal = _ribbon_temporal_summary(
        world, cont_labels, cont_ribbon, ribbon_threshold)
    primary = land_components[0] if land_components else {}
    primary_temporal_hint = _temporal_driver_hint(land_temporal)
    return {
        "schema": SCHEMA,
        "context": {
            "spec_name": str(world.spec.name),
            "seed": int(world.spec.seed),
            "cells": int(world.grid.n),
            "time_myr": float(world.time_myr),
        },
        "ribbon_threshold": float(ribbon_threshold),
        "summary": {
            "exposed_land_component_count": int(
                morphology.metrics["exposed_land"]["component_count"]),
            "continental_component_count": int(
                morphology.metrics["continental_crust"]["component_count"]),
            "exposed_land_ribbon_area_fraction_of_land": float(
                morphology.metrics["exposed_land"]["ribbon_area_fraction_gt_0_5"]),
            "continental_ribbon_area_fraction_of_continental_crust": float(
                morphology.metrics["continental_crust"]["ribbon_area_fraction_gt_0_5"]),
            "exposed_land_ribbon_area_fraction_of_world": (
                land_ribbon_area / total_area if total_area else 0.0
            ),
            "continental_ribbon_area_fraction_of_world": (
                cont_ribbon_area / total_area if total_area else 0.0
            ),
            "top_exposed_land_component_ribbon_share": float(
                primary.get("ribbon_share_of_total_ribbon", 0.0)),
            "primary_driver_hint": str(primary.get("driver_hint", "")),
            "exposed_land_active_rework_ribbon_share_lt220_myr": float(
                land_temporal["active_rework_ribbon_share_lt220_myr"]),
            "exposed_land_quiet_inherited_arc_suture_ribbon_share_gt650_myr": float(
                land_temporal["quiet_inherited_arc_suture_ribbon_share_gt650_myr"]),
            "exposed_land_mean_time_since_rework_myr": float(
                land_temporal["mean_time_since_rework_myr"]),
            "continental_active_rework_ribbon_share_lt220_myr": float(
                cont_temporal["active_rework_ribbon_share_lt220_myr"]),
            "continental_quiet_inherited_arc_suture_ribbon_share_gt650_myr": float(
                cont_temporal["quiet_inherited_arc_suture_ribbon_share_gt650_myr"]),
            "continental_mean_time_since_rework_myr": float(
                cont_temporal["mean_time_since_rework_myr"]),
            "primary_temporal_driver_hint": primary_temporal_hint,
        },
        "top_exposed_land_components": land_components,
        "top_continental_crust_components": cont_components,
        "top_exposed_continental_components": exposed_cont_components,
    }


def _top_components(
    world,
    components: list[dict[str, Any]],
    labels: np.ndarray,
    ribbon: np.ndarray,
    kind: str,
    *,
    top_n: int,
    ribbon_threshold: float,
) -> list[dict[str, Any]]:
    if not components:
        return []
    area = np.asarray(world.grid.cell_area, dtype=np.float64)
    total_area = float(area.sum())
    total_ribbon = _ribbon_area(area, labels, ribbon, ribbon_threshold)
    attributed = [
        _component_attribution(
            world,
            comp,
            labels,
            ribbon,
            kind,
            total_area=total_area,
            total_ribbon=total_ribbon,
            ribbon_threshold=ribbon_threshold,
        )
        for comp in components
    ]
    attributed.sort(
        key=lambda item: (
            -float(item["ribbon_area_fraction_of_world"]),
            -float(item["ribbon_area_fraction_of_component"]),
            -float(item["area_fraction_of_world"]),
            int(item["component_id"]),
        )
    )
    return attributed[:max(0, int(top_n))]


def _component_attribution(
    world,
    component: dict[str, Any],
    labels: np.ndarray,
    ribbon: np.ndarray,
    kind: str,
    *,
    total_area: float,
    total_ribbon: float,
    ribbon_threshold: float,
) -> dict[str, Any]:
    area = np.asarray(world.grid.cell_area, dtype=np.float64)
    cid = int(component["id"])
    cells = np.where(labels == cid)[0]
    if cells.size == 0:
        comp_area = 0.0
        ribbon_cells = cells
        weights = np.array([], dtype=np.float64)
    else:
        weights = area[cells]
        comp_area = float(weights.sum())
        ribbon_cells = cells[np.asarray(ribbon[cells]) >= ribbon_threshold]

    ribbon_area = float(area[ribbon_cells].sum()) if ribbon_cells.size else 0.0
    crust_type = world.get_field("crust.type", 0.0).astype(int)
    domain = world.get_field("crust.domain", DOMAIN_OCEANIC).astype(int)
    origin = world.get_field("crust.origin", ORIGIN_RIDGE).astype(int)
    stability = np.asarray(world.get_field("crust.stability", 0.0), dtype=np.float64)
    age = np.asarray(world.get_field("crust.age_myr", 0.0), dtype=np.float64)
    reworked = np.asarray(
        world.get_field("crust.reworked_age_myr", -1.0), dtype=np.float64)
    orogeny_age = np.asarray(
        world.get_field("tectonics.orogeny_age_myr", -1.0), dtype=np.float64)
    volcanism_age = np.asarray(
        world.get_field("tectonics.volcanism_age_myr", -1.0), dtype=np.float64)
    detail = world.get_field("terrain.continental_detail", CONT_DETAIL_OCEAN).astype(int)
    rel = np.asarray(world.get_field("terrain.elevation_m", 0.0), dtype=np.float64)
    rel = rel - float(world.sea_level)
    time_myr = float(world.time_myr)
    time_since_rework = _time_since_event(reworked[cells], time_myr)
    time_since_orogeny = _time_since_event(orogeny_age[cells], time_myr)
    time_since_volcanism = _time_since_event(volcanism_age[cells], time_myr)
    active_rework = (
        (time_since_rework < 220.0)
        | (time_since_orogeny < 220.0)
        | (time_since_volcanism < 220.0)
    )
    quiet_inherited = (
        np.isin(origin[cells], [int(ORIGIN_ARC), int(ORIGIN_SUTURE)])
        & (time_since_rework > 650.0)
        & (time_since_orogeny > 650.0)
        & (time_since_volcanism > 650.0)
    )

    out = {
        "component_id": cid,
        "component_kind": kind,
        "cell_count": int(component.get("cell_count", cells.size)),
        "area_fraction_of_world": comp_area / total_area if total_area else 0.0,
        "area_fraction_of_parent_mask": float(component.get("area_fraction_of_mask", 0.0)),
        "ribbon_area_fraction_of_world": ribbon_area / total_area if total_area else 0.0,
        "ribbon_area_fraction_of_component": (
            ribbon_area / comp_area if comp_area else 0.0
        ),
        "ribbon_share_of_total_ribbon": (
            ribbon_area / total_ribbon if total_ribbon else 0.0
        ),
        "mean_ribbon_index": _weighted_mean(ribbon[cells], weights),
        "width_p50_steps": float(component.get("width_p50_steps", 0.0)),
        "width_p90_steps": float(component.get("width_p90_steps", 0.0)),
        "width_max_steps": float(component.get("width_max_steps", 0.0)),
        "coastline_complexity": float(component.get("coastline_complexity", 0.0)),
        "compactness_proxy": float(component.get("compactness_proxy", 0.0)),
        "narrow_fraction_le2_width": float(
            component.get("narrow_fraction_le2_width", 0.0)),
        "narrow_neck_area_fraction": float(
            component.get("narrow_neck_area_fraction", 0.0)),
        "continental_crust_share": _weighted_fraction(
            crust_type[cells] == 1, weights),
        "oceanic_crust_share": _weighted_fraction(crust_type[cells] == 0, weights),
        "exposed_land_share": _weighted_fraction(rel[cells] >= 0.0, weights),
        "domain_shares": _code_shares(domain[cells], weights, DOMAIN_NAMES),
        "origin_shares": _code_shares(origin[cells], weights, ORIGIN_NAMES),
        "continental_detail_shares": _code_shares(detail[cells], weights, DETAIL_NAMES),
        "mean_stability": _weighted_mean(stability[cells], weights),
        "low_stability_share_lt030": _weighted_fraction(stability[cells] < 0.30, weights),
        "stable_craton_share_gt075": _weighted_fraction(stability[cells] > 0.75, weights),
        "mean_crust_age_myr": _weighted_mean(age[cells], weights),
        "ancient_share_gt2500_myr": _weighted_fraction(age[cells] > 2500.0, weights),
        "mean_time_since_rework_myr": _finite_weighted_mean(time_since_rework, weights),
        "recent_reworked_share_lt220_myr": _weighted_fraction(
            time_since_rework < 220.0, weights),
        "quiet_reworked_share_gt650_myr": _weighted_fraction(
            time_since_rework > 650.0, weights),
        "recent_orogeny_share_lt220_myr": _weighted_fraction(
            time_since_orogeny < 220.0, weights),
        "recent_volcanism_share_lt220_myr": _weighted_fraction(
            time_since_volcanism < 220.0, weights),
        "active_rework_share_lt220_myr": _weighted_fraction(active_rework, weights),
        "quiet_inherited_arc_suture_share_gt650_myr": _weighted_fraction(
            quiet_inherited, weights),
    }
    out["driver_hint"] = _driver_hint(out)
    return out


def _ribbon_temporal_summary(
    world,
    labels: np.ndarray,
    ribbon: np.ndarray,
    ribbon_threshold: float,
) -> dict[str, float]:
    area = np.asarray(world.grid.cell_area, dtype=np.float64)
    labels = np.asarray(labels, dtype=int)
    ribbon = np.asarray(ribbon, dtype=np.float64)
    cells = np.where((labels >= 0) & (ribbon >= float(ribbon_threshold)))[0]
    if cells.size == 0:
        return {
            "active_rework_ribbon_share_lt220_myr": 0.0,
            "quiet_inherited_arc_suture_ribbon_share_gt650_myr": 0.0,
            "recent_reworked_ribbon_share_lt220_myr": 0.0,
            "recent_orogeny_ribbon_share_lt220_myr": 0.0,
            "recent_volcanism_ribbon_share_lt220_myr": 0.0,
            "mean_time_since_rework_myr": 0.0,
        }
    weights = area[cells]
    time_myr = float(world.time_myr)
    origin = world.get_field("crust.origin", ORIGIN_RIDGE).astype(int)
    reworked = np.asarray(
        world.get_field("crust.reworked_age_myr", -1.0), dtype=np.float64)
    orogeny_age = np.asarray(
        world.get_field("tectonics.orogeny_age_myr", -1.0), dtype=np.float64)
    volcanism_age = np.asarray(
        world.get_field("tectonics.volcanism_age_myr", -1.0), dtype=np.float64)
    time_since_rework = _time_since_event(reworked[cells], time_myr)
    time_since_orogeny = _time_since_event(orogeny_age[cells], time_myr)
    time_since_volcanism = _time_since_event(volcanism_age[cells], time_myr)
    active = (
        (time_since_rework < 220.0)
        | (time_since_orogeny < 220.0)
        | (time_since_volcanism < 220.0)
    )
    quiet_inherited = (
        np.isin(origin[cells], [int(ORIGIN_ARC), int(ORIGIN_SUTURE)])
        & (time_since_rework > 650.0)
        & (time_since_orogeny > 650.0)
        & (time_since_volcanism > 650.0)
    )
    return {
        "active_rework_ribbon_share_lt220_myr": _weighted_fraction(active, weights),
        "quiet_inherited_arc_suture_ribbon_share_gt650_myr": _weighted_fraction(
            quiet_inherited, weights),
        "recent_reworked_ribbon_share_lt220_myr": _weighted_fraction(
            time_since_rework < 220.0, weights),
        "recent_orogeny_ribbon_share_lt220_myr": _weighted_fraction(
            time_since_orogeny < 220.0, weights),
        "recent_volcanism_ribbon_share_lt220_myr": _weighted_fraction(
            time_since_volcanism < 220.0, weights),
        "mean_time_since_rework_myr": _finite_weighted_mean(time_since_rework, weights),
    }


def _temporal_driver_hint(summary: dict[str, float]) -> str:
    active = float(summary.get("active_rework_ribbon_share_lt220_myr", 0.0))
    quiet = float(summary.get("quiet_inherited_arc_suture_ribbon_share_gt650_myr", 0.0))
    recent_orogeny = float(summary.get("recent_orogeny_ribbon_share_lt220_myr", 0.0))
    recent_volcanism = float(summary.get("recent_volcanism_ribbon_share_lt220_myr", 0.0))
    if active > 0.50 and recent_orogeny >= recent_volcanism:
        return "recent collision/suture rework dominates ribbon area"
    if active > 0.50:
        return "recent arc/subduction rework dominates ribbon area"
    if quiet > 0.45:
        return "quiet inherited arc/suture provenance dominates ribbon area"
    if active > quiet:
        return "mixed ribbon, leaning recent active rework"
    if quiet > active:
        return "mixed ribbon, leaning inherited provenance"
    return "weak or mixed temporal ribbon signal"


def _driver_hint(item: dict[str, Any]) -> str:
    domain = item["domain_shares"]
    origin = item["origin_shares"]
    detail = item["continental_detail_shares"]
    accretionary = (
        domain.get("continental_margin", 0.0)
        + domain.get("accreted_terrane", 0.0)
        + domain.get("suture", 0.0)
    )
    stable_core = domain.get("craton", 0.0) + domain.get("continental_interior", 0.0)
    arc_suture_origin = origin.get("arc", 0.0) + origin.get("suture", 0.0)
    active_detail = (
        detail.get("orogen", 0.0)
        + detail.get("plateau", 0.0)
        + detail.get("arc_microcontinent", 0.0)
    )

    if item["oceanic_crust_share"] > 0.30 and item["exposed_land_share"] > 0.50:
        return "exposed oceanic or arc-island chain"
    if accretionary > 0.55 and item["mean_stability"] < 0.45:
        return "young unstable accretionary-margin ribbon"
    if accretionary > 0.55 or arc_suture_origin > 0.55:
        return "narrow accretionary margin or suture collage"
    if stable_core > 0.60 and item["coastline_complexity"] > 8.0:
        return "stable continent with over-complex coastline"
    if item["narrow_neck_area_fraction"] > 0.05:
        return "component split by narrow articulation necks"
    if active_detail > 0.45:
        return "active highland belt dominates component planform"
    if item["ribbon_area_fraction_of_component"] > 0.30:
        return "underwide continent or margin component"
    return "low ribbon contribution"


def _ribbon_area(
    area: np.ndarray,
    labels: np.ndarray,
    ribbon: np.ndarray,
    threshold: float,
) -> float:
    mask = (labels >= 0) & (np.asarray(ribbon, dtype=np.float64) >= float(threshold))
    return float(area[mask].sum()) if mask.any() else 0.0


def _code_shares(
    values: np.ndarray,
    weights: np.ndarray,
    names: dict[int, str],
) -> dict[str, float]:
    values = np.asarray(values, dtype=int)
    weights = np.asarray(weights, dtype=np.float64)
    denom = float(weights.sum()) if weights.size else 0.0
    shares = {name: 0.0 for name in names.values()}
    if denom <= 0.0 or values.size == 0:
        return shares
    for code, name in names.items():
        shares[name] = float(weights[values == int(code)].sum() / denom)
    known = np.isin(values, list(names))
    if not bool(known.all()):
        shares["unknown"] = float(weights[~known].sum() / denom)
    return shares


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
    if not bool(valid.any()):
        return 0.0
    return float(np.average(values[valid], weights=weights[valid]))


def _finite_weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if values.size == 0 or weights.size == 0:
        return 0.0
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not bool(valid.any()):
        return float("inf")
    return float(np.average(values[valid], weights=weights[valid]))


def _time_since_event(event_age_myr: np.ndarray, time_myr: float) -> np.ndarray:
    event_age_myr = np.asarray(event_age_myr, dtype=np.float64)
    out = np.full(event_age_myr.shape, float(time_myr) + 1.0, dtype=np.float64)
    valid = np.isfinite(event_age_myr) & (event_age_myr >= 0.0)
    out[valid] = np.maximum(float(time_myr) - event_age_myr[valid], 0.0)
    return out
