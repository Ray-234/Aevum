"""Validation: the project must not merely "look right".

Five families of automatic checks (project plan, section 10):
  1. conservation  -- mass / energy / water / carbon inventories
  2. topology      -- plate coverage, acyclic rivers, ocean connectivity
  3. causality     -- resources / biomes / landforms have queryable origins
  4. regression    -- a fixed seed reproduces (determinism)
  5. world-combo   -- Earth is only one scenario, not the sole calibration target
"""
from __future__ import annotations

from dataclasses import dataclass, field
import heapq
from typing import Any

import numpy as np

from aevum.core.units import CONSTANTS
from aevum.diagnostics.morphology import compute_world_morphology
from aevum.modules.tectonics import (
    DOMAIN_ACCRETED_TERRANE,
    DOMAIN_LIP,
    ORIGIN_ARC,
    ORIGIN_PLUME_IMPACT,
)


@dataclass
class CheckResult:
    name: str
    passed: bool
    detail: dict[str, Any] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.passed


# ----------------------------------------------------------------------
def _component_areas(grid, mask: np.ndarray) -> list[float]:
    """Connected-component areas for a boolean mask on the sphere graph."""
    nodes = np.where(mask)[0]
    if nodes.size == 0:
        return []
    seen = np.zeros(grid.n, dtype=bool)
    out: list[float] = []
    for start in nodes:
        if seen[start]:
            continue
        stack = [int(start)]
        seen[start] = True
        acc = 0.0
        while stack:
            c = stack.pop()
            acc += float(grid.cell_area[c])
            for nb in grid.neighbors[c]:
                if mask[nb] and not seen[nb]:
                    seen[nb] = True
                    stack.append(int(nb))
        out.append(acc)
    return out


def _component_cell_sets(grid, mask: np.ndarray) -> list[np.ndarray]:
    """Connected-component cell arrays for a boolean mask on the sphere graph."""
    mask = np.asarray(mask, dtype=bool)
    nodes = np.where(mask)[0]
    if nodes.size == 0:
        return []
    seen = np.zeros(grid.n, dtype=bool)
    out: list[np.ndarray] = []
    for start in nodes:
        if seen[start]:
            continue
        stack = [int(start)]
        seen[start] = True
        comp: list[int] = []
        while stack:
            c = stack.pop()
            comp.append(c)
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if mask[nb] and not seen[nb]:
                    seen[nb] = True
                    stack.append(nb)
        out.append(np.asarray(comp, dtype=np.int64))
    return out


def _pct(values: np.ndarray, q: float) -> float:
    if values.size == 0:
        return 0.0
    return float(np.percentile(values, q))


def _area_frac(area: np.ndarray, mask: np.ndarray, denom_mask: np.ndarray | None = None) -> float:
    denom = float(area.sum()) if denom_mask is None else float(area[denom_mask].sum())
    if denom <= 0.0:
        return 0.0
    return float(area[mask].sum() / denom)


def _object_cells_mask(n: int, objects: list[dict[str, Any]],
                       kinds: set[str] | None = None) -> np.ndarray:
    mask = np.zeros(n, dtype=bool)
    for obj in objects:
        if kinds is not None and str(obj.get("kind", "")) not in kinds:
            continue
        cells = np.asarray(obj.get("cells", []), dtype=int)
        cells = cells[(0 <= cells) & (cells < n)]
        if cells.size:
            mask[cells] = True
    return mask


def _plate_fragmentation_from_array(grid, plate: np.ndarray) -> dict[str, Any]:
    area = grid.cell_area
    plate = np.asarray(plate).astype(int)
    valid = plate >= 0
    ids = sorted(int(x) for x in np.unique(plate[valid])) if valid.any() else []

    plate_area_fracs: list[float] = []
    component_counts: list[int] = []
    largest_component_shares: list[float] = []
    for pid in ids:
        mask = valid & (plate == pid)
        p_area = float(area[mask].sum())
        plate_area_fracs.append(p_area / float(area.sum()))
        comps = _component_areas(grid, mask)
        component_counts.append(len(comps))
        largest_component_shares.append(max(comps) / p_area if comps and p_area > 0.0 else 0.0)

    active = len(ids)
    fragmented = sum(1 for n in component_counts if n > 1)
    return {
        "n_active_plates": active,
        "min_plate_area_fraction": min(plate_area_fracs, default=0.0),
        "max_plate_area_fraction": max(plate_area_fracs, default=0.0),
        "mean_plate_area_fraction": float(np.mean(plate_area_fracs)) if plate_area_fracs else 0.0,
        "total_plate_components": int(sum(component_counts)),
        "max_components_per_plate": max(component_counts, default=0),
        "mean_components_per_plate": float(np.mean(component_counts)) if component_counts else 0.0,
        "fragmented_plate_fraction": fragmented / active if active else 0.0,
        "min_largest_component_share": min(largest_component_shares, default=0.0),
        "mean_largest_component_share": float(np.mean(largest_component_shares))
        if largest_component_shares else 0.0,
    }


def _plate_fragmentation_metrics(world) -> dict[str, Any]:
    return _plate_fragmentation_from_array(
        world.grid,
        world.get_field("tectonics.plate_id", -1.0),
    )


def _crust_distribution_metrics(world) -> dict[str, Any]:
    area = world.grid.cell_area
    ctype = world.get_field("crust.type", 0.0).astype(int)
    thick = world.get_field("crust.thickness_m", 0.0)
    age = world.get_field("crust.age_myr", 0.0)
    origin = world.get_field("crust.origin", 0.0).astype(int)
    reworked = world.get_field("crust.reworked_age_myr", -1.0)
    stability = world.get_field("crust.stability", 0.0)
    ocean = ctype == 0
    cont = ctype == 1
    total = float(area.sum())
    origin_counts = {
        f"origin_{int(code)}_fraction": _area_frac(area, origin == code)
        for code in np.unique(origin)
    }

    detail = {
        "oceanic_area_fraction": float(area[ocean].sum() / total) if total else 0.0,
        "continental_area_fraction": float(area[cont].sum() / total) if total else 0.0,
        "oceanic_age_p50_myr": _pct(age[ocean], 50),
        "oceanic_age_p95_myr": _pct(age[ocean], 95),
        "oceanic_age_max_myr": float(age[ocean].max()) if ocean.any() else 0.0,
        "oceanic_young_fraction_lt30_myr": _area_frac(area, ocean & (age < 30.0), ocean),
        "oceanic_old_fraction_gt300_myr": _area_frac(area, ocean & (age > 300.0), ocean),
        "continental_age_p50_myr": _pct(age[cont], 50),
        "continental_age_p95_myr": _pct(age[cont], 95),
        "continental_age_max_myr": float(age[cont].max()) if cont.any() else 0.0,
        "ancient_continental_fraction_gt2500_myr": _area_frac(
            area, cont & (age > 2500.0), cont),
        "continental_age_exceeds_world_cells": int((cont & (age > world.time_myr + 1.0)).sum()),
        "oceanic_thickness_p50_m": _pct(thick[ocean], 50),
        "oceanic_thickness_p90_m": _pct(thick[ocean], 90),
        "continental_thickness_p50_m": _pct(thick[cont], 50),
        "continental_thickness_p90_m": _pct(thick[cont], 90),
        "continental_stability_p50": _pct(stability[cont], 50),
        "continental_stability_p90": _pct(stability[cont], 90),
        "stable_craton_fraction_gt075": _area_frac(area, cont & (stability > 0.75), cont),
        "recently_reworked_fraction": _area_frac(area, reworked >= max(world.time_myr - 300.0, 0.0)),
        "negative_age_cells": int((age < 0.0).sum()),
        "negative_thickness_cells": int((thick < 0.0).sum()),
        "non_binary_crust_cells": int((~np.isin(ctype, [0, 1])).sum()),
        "invalid_origin_cells": int((origin < 0).sum()),
        "stability_out_of_range_cells": int(((stability < 0.0) | (stability > 1.0)).sum()),
        "continental_ridge_origin_cells": int((cont & (origin == 0)).sum()),
    }
    detail.update(origin_counts)
    return detail


def _boundary_metrics(world) -> dict[str, Any]:
    boundaries = world.networks.get("tectonics.boundaries", {})
    counts: dict[str, int] = {}
    all_cells: list[np.ndarray] = []
    for name in ("divergent", "ridge", "convergent", "collision", "subduction",
                 "trench", "suture", "active_margin", "passive_margin", "transform"):
        cells = np.asarray(boundaries.get(name, []), dtype=int)
        counts[f"{name}_cells"] = int(np.unique(cells).size)
        if cells.size:
            all_cells.append(cells)
    unique = np.unique(np.concatenate(all_cells)) if all_cells else np.array([], dtype=int)
    counts["any_boundary_cells"] = int(unique.size)
    counts["boundary_cell_fraction"] = float(unique.size / max(world.grid.n, 1))
    return counts


def _graph_distance_from_sources(grid, sources: np.ndarray, allowed: np.ndarray) -> np.ndarray:
    sources = np.asarray(sources, dtype=bool)
    allowed = np.asarray(allowed, dtype=bool)
    dist = np.full(grid.n, np.inf, dtype=np.float64)
    heap: list[tuple[float, int]] = []
    for c in np.where(sources & allowed)[0]:
        dist[c] = 0.0
        heapq.heappush(heap, (0.0, int(c)))
    while heap:
        d, c = heapq.heappop(heap)
        if d != dist[c]:
            continue
        for nb in grid.neighbors[c]:
            nb = int(nb)
            if not allowed[nb]:
                continue
            nd = d + grid.great_circle_distance(c, nb)
            if nd < dist[nb]:
                dist[nb] = nd
                heapq.heappush(heap, (nd, nb))
    return dist


def _dilate_mask(grid, mask: np.ndarray, allowed: np.ndarray, passes: int = 1) -> np.ndarray:
    out = np.asarray(mask, dtype=bool).copy()
    allowed = np.asarray(allowed, dtype=bool)
    for _ in range(passes):
        add = np.zeros(grid.n, dtype=bool)
        for c in np.where(out)[0]:
            add[grid.neighbors[c]] = True
        out |= add & allowed
    return out


def _weighted_mean(values: np.ndarray, weights: np.ndarray, mask: np.ndarray | None = None,
                   default: float = 0.0) -> float:
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        values = values[mask]
        weights = weights[mask]
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not valid.any():
        return default
    return float(np.average(values[valid], weights=weights[valid]))


def _neighbor_scalar_delta_stats(grid, values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    if grid.edges.size == 0:
        return {
            "neighbor_delta_p95": 0.0,
            "neighbor_delta_p99": 0.0,
            "neighbor_delta_max": 0.0,
        }
    i, j = grid.edges[:, 0], grid.edges[:, 1]
    delta = np.abs(values[i] - values[j])
    delta = delta[np.isfinite(delta)]
    if delta.size == 0:
        return {
            "neighbor_delta_p95": 0.0,
            "neighbor_delta_p99": 0.0,
            "neighbor_delta_max": 0.0,
        }
    return {
        "neighbor_delta_p95": float(np.percentile(delta, 95)),
        "neighbor_delta_p99": float(np.percentile(delta, 99)),
        "neighbor_delta_max": float(delta.max()),
    }


def _lat_band_temperature_metrics(world, temp_c: np.ndarray,
                                  band_width_deg: float = 10.0) -> dict[str, Any]:
    grid = world.grid
    bins = np.arange(-90.0, 90.0 + band_width_deg, band_width_deg)
    means: list[float] = []
    ranges: list[tuple[float, float]] = []
    for idx, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
        if idx == len(bins) - 2:
            mask = (grid.lat >= lo) & (grid.lat <= hi)
        else:
            mask = (grid.lat >= lo) & (grid.lat < hi)
        means.append(_weighted_mean(temp_c, grid.cell_area, mask))
        ranges.append((float(lo), float(hi)))

    deltas = [abs(b - a) for a, b in zip(means[:-1], means[1:])]
    if deltas:
        max_idx = int(np.argmax(deltas))
        max_delta = float(deltas[max_idx])
        pair = [ranges[max_idx], ranges[max_idx + 1]]
    else:
        max_delta = 0.0
        pair = []

    return {
        "lat_band_width_deg": float(band_width_deg),
        "lat_band_means_C": [round(float(x), 3) for x in means],
        "max_adjacent_lat_band_delta_C": max_delta,
        "max_adjacent_lat_band_pair_deg": pair,
    }


def _lat_band_residual_abs_p95(
    grid,
    values: np.ndarray,
    mask: np.ndarray | None = None,
    band_width_deg: float = 10.0,
) -> float:
    values = np.asarray(values, dtype=np.float64)
    active = np.isfinite(values)
    if mask is not None:
        active &= np.asarray(mask, dtype=bool)
    if int(np.count_nonzero(active)) < 8:
        return 0.0
    bins = np.arange(-90.0, 90.0 + band_width_deg, band_width_deg)
    predicted = np.full(values.shape, np.nan, dtype=np.float64)
    for idx, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
        if idx == len(bins) - 2:
            band = active & (grid.lat >= lo) & (grid.lat <= hi)
        else:
            band = active & (grid.lat >= lo) & (grid.lat < hi)
        if band.any():
            predicted[band] = _weighted_mean(values, grid.cell_area, band)
    valid = active & np.isfinite(predicted)
    if int(np.count_nonzero(valid)) < 8:
        return 0.0
    return _pct(np.abs(values[valid] - predicted[valid]), 95)


def _temperature_metrics(world) -> dict[str, Any]:
    grid = world.grid
    area = grid.cell_area
    has_temp = "climate.surface_temperature" in world.fields
    temp_k = world.get_field("climate.surface_temperature", 288.0)
    temp_c = temp_k - 273.15
    land = world.land_mask()
    ocean = ~land
    band = _lat_band_temperature_metrics(world, temp_c)
    neighbor = _neighbor_scalar_delta_stats(grid, temp_c)
    land_mean = _weighted_mean(temp_c, area, land)
    ocean_mean = _weighted_mean(temp_c, area, ocean)
    return {
        "has_surface_temperature": bool(has_temp),
        "nonfinite_temperature_cells": int((~np.isfinite(temp_k)).sum()),
        "mean_temp_C": _weighted_mean(temp_c, area),
        "min_temp_C": float(np.nanmin(temp_c)) if temp_c.size else 0.0,
        "max_temp_C": float(np.nanmax(temp_c)) if temp_c.size else 0.0,
        "land_mean_temp_C": land_mean,
        "ocean_mean_temp_C": ocean_mean,
        "land_ocean_temperature_contrast_C": land_mean - ocean_mean,
        "abs_land_ocean_temperature_contrast_C": abs(land_mean - ocean_mean),
        "max_adjacent_lat_band_delta_C": band["max_adjacent_lat_band_delta_C"],
        "max_adjacent_lat_band_pair_deg": band["max_adjacent_lat_band_pair_deg"],
        "lat_band_means_C": band["lat_band_means_C"],
        "lat_band_residual_abs_p95_C": _lat_band_residual_abs_p95(
            grid, temp_c, np.isfinite(temp_c)),
        "neighbor_temperature_delta_p95_C": neighbor["neighbor_delta_p95"],
        "neighbor_temperature_delta_p99_C": neighbor["neighbor_delta_p99"],
        "neighbor_temperature_delta_max_C": neighbor["neighbor_delta_max"],
    }


def _land_relief_m(grid, elev: np.ndarray) -> np.ndarray:
    relief = np.zeros(grid.n, dtype=np.float64)
    for c, nbs in enumerate(grid.neighbors):
        if nbs.size:
            relief[c] = float(np.max(np.abs(elev[nbs] - elev[c])))
    return relief


def _precipitation_metrics(world) -> dict[str, Any]:
    grid = world.grid
    area = grid.cell_area
    has_precip = "climate.precipitation" in world.fields
    precip = world.get_field("climate.precipitation", 0.0)
    elev = world.get_field("terrain.elevation_m", 0.0) - world.sea_level
    land = world.land_mask()

    detail: dict[str, Any] = {
        "has_precipitation": bool(has_precip),
        "nonfinite_precipitation_cells": int((~np.isfinite(precip)).sum()),
        "negative_precipitation_cells": int((precip < 0.0).sum()),
        "mean_precip_mm_yr": _weighted_mean(precip, area),
        "land_precip_p50_mm_yr": _pct(precip[land], 50),
        "land_precip_p75_mm_yr": _pct(precip[land], 75),
        "land_precip_p95_mm_yr": _pct(precip[land], 95),
        "land_wet_fraction_gt500mm": _area_frac(area, land & (precip > 500.0), land),
        "precip_orographic_concentration": 0.0,
        "orographic_area_fraction": 0.0,
        "orographic_precip_fraction": 0.0,
        "precip_relief_correlation": 0.0,
    }
    if not land.any():
        return detail

    relief = _land_relief_m(grid, elev)
    relief_land = relief[land]
    if relief_land.size == 0:
        return detail

    cutoff = float(np.percentile(relief_land, 85))
    orographic = land & (relief >= cutoff) & (relief > 0.0)
    land_area = float(area[land].sum())
    land_precip_total = float(np.sum(np.maximum(precip[land], 0.0) * area[land]))
    oro_area_fraction = float(area[orographic].sum() / land_area) if land_area > 0.0 else 0.0
    oro_precip_total = float(np.sum(np.maximum(precip[orographic], 0.0) * area[orographic]))
    oro_precip_fraction = (
        oro_precip_total / land_precip_total if land_precip_total > 0.0 else 0.0
    )
    concentration = (
        oro_precip_fraction / max(oro_area_fraction, 1e-12)
        if oro_area_fraction > 0.0 else 0.0
    )
    relief_valid = relief[land]
    precip_valid = precip[land]
    valid = np.isfinite(relief_valid) & np.isfinite(precip_valid)
    corr = 0.0
    if valid.sum() >= 8 and float(np.std(relief_valid[valid])) > 0.0:
        if float(np.std(precip_valid[valid])) > 0.0:
            corr = float(np.corrcoef(relief_valid[valid], precip_valid[valid])[0, 1])

    detail.update({
        "precip_orographic_concentration": concentration,
        "orographic_area_fraction": oro_area_fraction,
        "orographic_precip_fraction": oro_precip_fraction,
        "precip_relief_correlation": corr,
    })
    return detail


def _coastal_temperature_metrics(world) -> dict[str, Any]:
    grid = world.grid
    area = grid.cell_area
    temp_c = world.get_field("climate.surface_temperature", 288.0) - 273.15
    land = world.land_mask()
    ocean = ~land
    east_axis = np.array([0.0, 0.0, 1.0])
    coastal = np.zeros(grid.n, dtype=bool)
    east_component = np.zeros(grid.n, dtype=np.float64)

    for c in np.where(land)[0]:
        nbs = grid.neighbors[int(c)]
        ocean_nbs = nbs[ocean[nbs]]
        if ocean_nbs.size == 0:
            continue
        coastal[c] = True
        ocean_vec = np.mean(grid.xyz[ocean_nbs], axis=0)
        normal = grid.xyz[c]
        tangent = ocean_vec - normal * float(ocean_vec @ normal)
        norm = float(np.linalg.norm(tangent))
        if norm <= 1e-12:
            continue
        tangent /= norm
        east = np.cross(east_axis, normal)
        east_norm = float(np.linalg.norm(east))
        if east_norm <= 1e-12:
            continue
        east /= east_norm
        east_component[c] = float(tangent @ east)

    east_facing = coastal & (east_component > 0.15)
    west_facing = coastal & (east_component < -0.15)
    band_diffs: list[float] = []
    band_width = 15.0
    bins = np.arange(-75.0, 75.0 + band_width, band_width)
    for idx, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
        if idx == len(bins) - 2:
            band = (grid.lat >= lo) & (grid.lat <= hi)
        else:
            band = (grid.lat >= lo) & (grid.lat < hi)
        east_mask = east_facing & band
        west_mask = west_facing & band
        if int(east_mask.sum()) < 3 or int(west_mask.sum()) < 3:
            continue
        east_mean = _weighted_mean(temp_c, area, east_mask)
        west_mean = _weighted_mean(temp_c, area, west_mask)
        band_diffs.append(abs(east_mean - west_mean))

    asymmetry = float(np.mean(band_diffs)) if band_diffs else 0.0
    return {
        "coastal_land_cells": int(coastal.sum()),
        "east_facing_coast_cells": int(east_facing.sum()),
        "west_facing_coast_cells": int(west_facing.sum()),
        "coastal_temperature_asymmetry_index": asymmetry,
        "coastal_temperature_bands_compared": int(len(band_diffs)),
        "coastal_temperature_band_diffs_C": [round(float(x), 3) for x in band_diffs],
    }


def _seasonality_metrics(world) -> dict[str, Any]:
    land = world.land_mask()
    ocean = ~land
    area = world.grid.cell_area
    n = world.grid.n
    seasonal_temp = world.fields.get("climate.seasonal_temperature")
    seasonal_precip = world.fields.get("climate.seasonal_precipitation")
    monsoon_corridor = world.fields.get("climate.monsoon_rainfall_corridor")
    storm_corridor = world.fields.get("climate.storm_track_rainfall_corridor")
    rain_shadow = world.fields.get("climate.rain_shadow_index")
    regional_response = world.fields.get("climate.regional_precipitation_response")
    moisture_flow_source = world.fields.get("atmosphere.moisture_flow_source")
    moisture_flow_pathway = world.fields.get("atmosphere.moisture_flow_pathway")
    moisture_source_basin_id = world.fields.get("atmosphere.moisture_source_basin_id")
    moisture_flow_network_id = world.fields.get("climate.moisture_flow_network_id")
    moisture_flow_precip_response = world.fields.get(
        "climate.moisture_flow_precipitation_response")
    moisture_budget_region_id = world.fields.get("climate.moisture_budget_region_id")
    precipitation_response_region_id = world.fields.get(
        "climate.precipitation_response_region_id")
    receiver_catchment_id = world.fields.get("climate.receiver_catchment_id")
    source_basin_supply_index = world.fields.get("climate.source_basin_supply_index")
    receiver_supply_balance = world.fields.get(
        "climate.receiver_catchment_supply_balance")
    receiver_supply_feedback = world.fields.get(
        "climate.receiver_supply_precipitation_feedback")
    hydro_regions = world.objects.get("climate.hydroclimate_regions", [])
    if not isinstance(hydro_regions, list):
        hydro_regions = []
    moisture_networks = world.objects.get("climate.moisture_flow_networks", [])
    if not isinstance(moisture_networks, list):
        moisture_networks = []
    precipitation_response_regions = world.objects.get(
        "climate.precipitation_response_regions", [])
    if not isinstance(precipitation_response_regions, list):
        precipitation_response_regions = []
    receiver_catchments = world.objects.get("climate.receiver_catchments", [])
    if not isinstance(receiver_catchments, list):
        receiver_catchments = []
    detail: dict[str, Any] = {
        "has_seasonal_temperature": seasonal_temp is not None,
        "has_seasonal_precipitation": seasonal_precip is not None,
        "has_monsoon_rainfall_corridor": monsoon_corridor is not None,
        "has_storm_track_rainfall_corridor": storm_corridor is not None,
        "has_rain_shadow_index": rain_shadow is not None,
        "has_regional_precipitation_response": regional_response is not None,
        "has_hydroclimate_region_objects": bool(hydro_regions),
        "has_moisture_flow_source": moisture_flow_source is not None,
        "has_moisture_flow_pathway": moisture_flow_pathway is not None,
        "has_moisture_source_basin_id": moisture_source_basin_id is not None,
        "has_moisture_flow_network_id": moisture_flow_network_id is not None,
        "has_moisture_flow_precipitation_response": (
            moisture_flow_precip_response is not None),
        "has_moisture_budget_region_id": moisture_budget_region_id is not None,
        "has_moisture_flow_network_objects": bool(moisture_networks),
        "has_precipitation_response_region_id": (
            precipitation_response_region_id is not None),
        "has_precipitation_response_region_objects": bool(
            precipitation_response_regions),
        "has_receiver_catchment_id": receiver_catchment_id is not None,
        "has_source_basin_supply_index": source_basin_supply_index is not None,
        "has_receiver_catchment_supply_balance": receiver_supply_balance is not None,
        "has_receiver_supply_precipitation_feedback": (
            receiver_supply_feedback is not None),
        "has_receiver_catchment_objects": bool(receiver_catchments),
        "invalid_seasonal_temperature_shape": False,
        "invalid_seasonal_precipitation_shape": False,
        "invalid_monsoon_rainfall_corridor_shape": False,
        "invalid_storm_track_rainfall_corridor_shape": False,
        "invalid_rain_shadow_index_shape": False,
        "invalid_regional_precipitation_response_shape": False,
        "invalid_moisture_flow_source_shape": False,
        "invalid_moisture_flow_pathway_shape": False,
        "invalid_moisture_source_basin_id_shape": False,
        "invalid_moisture_flow_network_id_shape": False,
        "invalid_moisture_flow_precipitation_response_shape": False,
        "invalid_moisture_budget_region_id_shape": False,
        "invalid_precipitation_response_region_id_shape": False,
        "invalid_receiver_catchment_id_shape": False,
        "invalid_source_basin_supply_index_shape": False,
        "invalid_receiver_catchment_supply_balance_shape": False,
        "invalid_receiver_supply_precipitation_feedback_shape": False,
        "nonfinite_seasonal_precipitation_cells": 0,
        "negative_seasonal_precipitation_cells": 0,
        "nonfinite_monsoon_rainfall_corridor_cells": 0,
        "nonfinite_storm_track_rainfall_corridor_cells": 0,
        "nonfinite_rain_shadow_index_cells": 0,
        "nonfinite_regional_precipitation_response_cells": 0,
        "nonfinite_moisture_flow_source_cells": 0,
        "nonfinite_moisture_flow_pathway_cells": 0,
        "nonfinite_moisture_source_basin_id_cells": 0,
        "nonfinite_moisture_flow_network_id_cells": 0,
        "nonfinite_moisture_flow_precipitation_response_cells": 0,
        "nonfinite_moisture_budget_region_id_cells": 0,
        "nonfinite_precipitation_response_region_id_cells": 0,
        "nonfinite_receiver_catchment_id_cells": 0,
        "nonfinite_source_basin_supply_index_cells": 0,
        "nonfinite_receiver_catchment_supply_balance_cells": 0,
        "nonfinite_receiver_supply_precipitation_feedback_cells": 0,
        "annual_precip_matches_seasonal_aggregate": False,
        "annual_precip_seasonal_aggregate_max_delta_mm_yr": 0.0,
        "land_seasonal_temp_amplitude_p50_C": 0.0,
        "ocean_seasonal_temp_amplitude_p50_C": 0.0,
        "precip_seasonality_p75": 0.0,
        "monsoon_rainfall_corridor_land_p90": 0.0,
        "storm_track_rainfall_corridor_land_p90": 0.0,
        "rain_shadow_index_land_p90": 0.0,
        "regional_precipitation_response_land_p05": 1.0,
        "regional_precipitation_response_land_p95": 1.0,
        "moisture_flow_pathway_land_p90": 0.0,
        "moisture_flow_source_ocean_p90": 0.0,
        "moisture_source_basin_attributed_land_fraction": 0.0,
        "source_basin_supply_index_land_p50": 0.0,
        "source_basin_supply_attributed_land_fraction": 0.0,
        "receiver_catchment_supply_balance_land_p50": 0.0,
        "receiver_supply_precipitation_feedback_land_p05": 1.0,
        "receiver_supply_precipitation_feedback_land_p95": 1.0,
        "receiver_supply_precipitation_feedback_ocean_abs_p95": 0.0,
        "moisture_flow_precipitation_response_land_p05": 1.0,
        "moisture_flow_precipitation_response_land_p95": 1.0,
        "moisture_budget_region_count_p50": 0.0,
        "hydroclimate_region_object_count": 0,
        "hydroclimate_region_kind_count": 0,
        "hydroclimate_region_season_count": 0,
        "monsoon_region_object_count": 0,
        "storm_track_region_object_count": 0,
        "rain_shadow_region_object_count": 0,
        "wet_response_region_object_count": 0,
        "dry_response_region_object_count": 0,
        "largest_hydroclimate_region_area_fraction": 0.0,
        "hydroclimate_region_mean_intensity_p50": 0.0,
        "moisture_flow_network_object_count": 0,
        "moisture_flow_network_kind_count": 0,
        "moisture_flow_network_season_count": 0,
        "largest_moisture_flow_network_area_fraction": 0.0,
        "moisture_flow_network_mean_pathway_p50": 0.0,
        "precipitation_response_region_count_p50": 0.0,
        "precipitation_response_region_object_count": 0,
        "precipitation_response_region_kind_count": 0,
        "precipitation_response_region_season_count": 0,
        "wet_precipitation_response_region_object_count": 0,
        "dry_precipitation_response_region_object_count": 0,
        "largest_precipitation_response_region_area_fraction": 0.0,
        "precipitation_response_region_mean_abs_anomaly_p50": 0.0,
        "precipitation_response_region_source_basin_attribution_p50": 0.0,
        "precipitation_response_region_budget_attribution_p50": 0.0,
        "precipitation_response_region_flow_network_attribution_p50": 0.0,
        "receiver_catchment_count_p50": 0.0,
        "receiver_catchment_object_count": 0,
        "receiver_catchment_kind_count": 0,
        "receiver_catchment_season_count": 0,
        "source_receiver_catchment_object_count": 0,
        "mixed_receiver_catchment_object_count": 0,
        "largest_receiver_catchment_area_fraction": 0.0,
        "receiver_catchment_source_basin_attribution_p50": 0.0,
        "receiver_catchment_budget_attribution_p50": 0.0,
        "receiver_catchment_precip_response_attribution_p50": 0.0,
        "receiver_catchment_source_supply_p50": 0.0,
        "receiver_catchment_supply_balance_p50": 0.0,
        "receiver_catchment_supported_precip_fraction_p50": 0.0,
    }

    if seasonal_temp is not None:
        arr = np.asarray(seasonal_temp, dtype=np.float64)
        detail["invalid_seasonal_temperature_shape"] = bool(arr.shape != (4, n))
        if arr.shape == (4, n):
            amp = np.nanmax(arr, axis=0) - np.nanmin(arr, axis=0)
            detail["land_seasonal_temp_amplitude_p50_C"] = _pct(amp[land], 50)
            detail["ocean_seasonal_temp_amplitude_p50_C"] = _pct(amp[ocean], 50)
            detail["global_seasonal_temp_amplitude_mean_C"] = _weighted_mean(amp, area)

    if seasonal_precip is not None:
        arr = np.asarray(seasonal_precip, dtype=np.float64)
        detail["invalid_seasonal_precipitation_shape"] = bool(arr.shape != (4, n))
        detail["nonfinite_seasonal_precipitation_cells"] = int((~np.isfinite(arr)).sum())
        detail["negative_seasonal_precipitation_cells"] = int((arr < 0.0).sum())
        if arr.shape == (4, n):
            annual_mean = np.maximum(np.nanmean(arr, axis=0), 1e-9)
            seasonality = np.nanmax(arr, axis=0) / annual_mean
            detail["precip_seasonality_p75"] = _pct(seasonality[land], 75)
            detail["global_precip_seasonality_mean"] = _weighted_mean(seasonality, area)
            annual_precip = world.get_field("climate.precipitation", 0.0)
            delta = np.abs(annual_precip - np.nanmean(arr, axis=0))
            detail["annual_precip_seasonal_aggregate_max_delta_mm_yr"] = _pct(delta, 100)
            detail["annual_precip_matches_seasonal_aggregate"] = bool(
                float(np.nanmax(delta)) < 1e-6)

    for arr_obj, invalid_key, nonfinite_key, p90_key in [
        (monsoon_corridor, "invalid_monsoon_rainfall_corridor_shape",
         "nonfinite_monsoon_rainfall_corridor_cells",
         "monsoon_rainfall_corridor_land_p90"),
        (storm_corridor, "invalid_storm_track_rainfall_corridor_shape",
         "nonfinite_storm_track_rainfall_corridor_cells",
         "storm_track_rainfall_corridor_land_p90"),
        (rain_shadow, "invalid_rain_shadow_index_shape",
         "nonfinite_rain_shadow_index_cells", "rain_shadow_index_land_p90"),
    ]:
        if arr_obj is None:
            continue
        arr = np.asarray(arr_obj, dtype=np.float64)
        detail[invalid_key] = bool(arr.shape != (4, n))
        detail[nonfinite_key] = int((~np.isfinite(arr)).sum())
        if arr.shape == (4, n) and land.any():
            detail[p90_key] = _pct(arr[:, land].ravel(), 90)

    if regional_response is not None:
        arr = np.asarray(regional_response, dtype=np.float64)
        detail["invalid_regional_precipitation_response_shape"] = bool(
            arr.shape != (4, n))
        detail["nonfinite_regional_precipitation_response_cells"] = int(
            (~np.isfinite(arr)).sum())
        if arr.shape == (4, n) and land.any():
            vals = arr[:, land].ravel()
            detail["regional_precipitation_response_land_p05"] = _pct(vals, 5)
            detail["regional_precipitation_response_land_p95"] = _pct(vals, 95)

    for arr_obj, invalid_key, nonfinite_key, p90_key, domain in [
        (moisture_flow_source, "invalid_moisture_flow_source_shape",
         "nonfinite_moisture_flow_source_cells",
         "moisture_flow_source_ocean_p90", ocean),
        (moisture_flow_pathway, "invalid_moisture_flow_pathway_shape",
         "nonfinite_moisture_flow_pathway_cells",
         "moisture_flow_pathway_land_p90", land),
        (moisture_source_basin_id, "invalid_moisture_source_basin_id_shape",
         "nonfinite_moisture_source_basin_id_cells", None, None),
        (moisture_flow_network_id, "invalid_moisture_flow_network_id_shape",
         "nonfinite_moisture_flow_network_id_cells", None, None),
        (moisture_budget_region_id, "invalid_moisture_budget_region_id_shape",
         "nonfinite_moisture_budget_region_id_cells", None, None),
        (precipitation_response_region_id,
         "invalid_precipitation_response_region_id_shape",
         "nonfinite_precipitation_response_region_id_cells", None, None),
        (receiver_catchment_id,
         "invalid_receiver_catchment_id_shape",
         "nonfinite_receiver_catchment_id_cells", None, None),
        (source_basin_supply_index,
         "invalid_source_basin_supply_index_shape",
         "nonfinite_source_basin_supply_index_cells", None, None),
        (receiver_supply_balance,
         "invalid_receiver_catchment_supply_balance_shape",
         "nonfinite_receiver_catchment_supply_balance_cells", None, None),
        (receiver_supply_feedback,
         "invalid_receiver_supply_precipitation_feedback_shape",
         "nonfinite_receiver_supply_precipitation_feedback_cells", None, None),
    ]:
        if arr_obj is None:
            continue
        arr = np.asarray(arr_obj, dtype=np.float64)
        detail[invalid_key] = bool(arr.shape != (4, n))
        detail[nonfinite_key] = int((~np.isfinite(arr)).sum())
        if arr.shape == (4, n) and domain is not None and domain.any():
            detail[p90_key] = _pct(arr[:, domain].ravel(), 90)
        if arr_obj is moisture_budget_region_id and arr.shape == (4, n) and land.any():
            counts = [
                len([
                    x for x in np.unique(arr[season, land])
                    if np.isfinite(x) and int(x) > 0
                ])
                for season in range(4)
            ]
            detail["moisture_budget_region_count_p50"] = (
                float(np.percentile(counts, 50)) if counts else 0.0
            )
        if (arr_obj is precipitation_response_region_id
                and arr.shape == (4, n) and land.any()):
            counts = [
                len([
                    x for x in np.unique(arr[season, land])
                    if np.isfinite(x) and int(x) > 0
                ])
                for season in range(4)
            ]
            detail["precipitation_response_region_count_p50"] = (
                float(np.percentile(counts, 50)) if counts else 0.0
            )
        if arr_obj is receiver_catchment_id and arr.shape == (4, n) and land.any():
            counts = [
                len([
                    x for x in np.unique(arr[season, land])
                    if np.isfinite(x) and int(x) > 0
                ])
                for season in range(4)
            ]
            detail["receiver_catchment_count_p50"] = (
                float(np.percentile(counts, 50)) if counts else 0.0
            )

    if moisture_source_basin_id is not None:
        arr = np.asarray(moisture_source_basin_id, dtype=np.float64)
        if arr.shape == (4, n) and land.any():
            attributed = land[None, :] & np.isfinite(arr) & (arr >= 0.0)
            detail["moisture_source_basin_attributed_land_fraction"] = float(
                attributed.sum() / max(4 * int(land.sum()), 1)
            )

    if source_basin_supply_index is not None:
        arr = np.asarray(source_basin_supply_index, dtype=np.float64)
        if arr.shape == (4, n) and land.any():
            vals = arr[:, land].ravel()
            vals = vals[np.isfinite(vals)]
            detail["source_basin_supply_index_land_p50"] = (
                float(np.percentile(vals, 50)) if vals.size else 0.0
            )
            if moisture_source_basin_id is not None:
                source_arr = np.asarray(moisture_source_basin_id, dtype=np.float64)
                if source_arr.shape == (4, n):
                    attributed = (
                        land[None, :]
                        & np.isfinite(source_arr)
                        & (source_arr >= 0.0)
                        & np.isfinite(arr)
                        & (arr > 0.05)
                    )
                    detail["source_basin_supply_attributed_land_fraction"] = float(
                        np.count_nonzero(attributed)
                        / max(4 * int(np.count_nonzero(land)), 1)
                    )

    if receiver_supply_balance is not None:
        arr = np.asarray(receiver_supply_balance, dtype=np.float64)
        if arr.shape == (4, n) and land.any():
            vals = arr[:, land].ravel()
            vals = vals[np.isfinite(vals)]
            detail["receiver_catchment_supply_balance_land_p50"] = (
                float(np.percentile(vals, 50)) if vals.size else 0.0
            )

    if receiver_supply_feedback is not None:
        arr = np.asarray(receiver_supply_feedback, dtype=np.float64)
        if arr.shape == (4, n):
            if land.any():
                vals = arr[:, land].ravel()
                vals = vals[np.isfinite(vals)]
                if vals.size:
                    detail["receiver_supply_precipitation_feedback_land_p05"] = (
                        float(np.percentile(vals, 5)))
                    detail["receiver_supply_precipitation_feedback_land_p95"] = (
                        float(np.percentile(vals, 95)))
            if ocean.any():
                ocean_vals = np.abs(arr[:, ocean].ravel() - 1.0)
                ocean_vals = ocean_vals[np.isfinite(ocean_vals)]
                detail["receiver_supply_precipitation_feedback_ocean_abs_p95"] = (
                    float(np.percentile(ocean_vals, 95)) if ocean_vals.size else 0.0
                )

    if moisture_flow_precip_response is not None:
        arr = np.asarray(moisture_flow_precip_response, dtype=np.float64)
        detail["invalid_moisture_flow_precipitation_response_shape"] = bool(
            arr.shape != (4, n))
        detail["nonfinite_moisture_flow_precipitation_response_cells"] = int(
            (~np.isfinite(arr)).sum())
        if arr.shape == (4, n) and land.any():
            vals = arr[:, land].ravel()
            detail["moisture_flow_precipitation_response_land_p05"] = _pct(vals, 5)
            detail["moisture_flow_precipitation_response_land_p95"] = _pct(vals, 95)

    if hydro_regions:
        kinds = [str(obj.get("kind", "")) for obj in hydro_regions]
        seasons = [str(obj.get("season", "")) for obj in hydro_regions]
        area_fracs = [
            float(obj.get("area_fraction", 0.0))
            for obj in hydro_regions
            if np.isfinite(float(obj.get("area_fraction", 0.0)))
        ]
        intensities = [
            float(obj.get("mean_intensity", 0.0))
            for obj in hydro_regions
            if np.isfinite(float(obj.get("mean_intensity", 0.0)))
        ]
        detail["hydroclimate_region_object_count"] = int(len(hydro_regions))
        detail["hydroclimate_region_kind_count"] = int(len(set(kinds)))
        detail["hydroclimate_region_season_count"] = int(len(set(seasons)))
        detail["monsoon_region_object_count"] = int(
            sum(1 for kind in kinds if kind == "monsoon_rainfall_corridor"))
        detail["storm_track_region_object_count"] = int(
            sum(1 for kind in kinds if kind == "storm_track_rainfall_corridor"))
        detail["rain_shadow_region_object_count"] = int(
            sum(1 for kind in kinds if kind == "rain_shadow_region"))
        detail["wet_response_region_object_count"] = int(
            sum(1 for kind in kinds if kind == "wet_regional_precipitation_response"))
        detail["dry_response_region_object_count"] = int(
            sum(1 for kind in kinds if kind == "dry_regional_precipitation_response"))
        detail["largest_hydroclimate_region_area_fraction"] = (
            max(area_fracs) if area_fracs else 0.0
        )
        detail["hydroclimate_region_mean_intensity_p50"] = (
            float(np.percentile(intensities, 50)) if intensities else 0.0
        )

    if moisture_networks:
        kinds = [str(obj.get("kind", "")) for obj in moisture_networks]
        seasons = [str(obj.get("season", "")) for obj in moisture_networks]
        area_fracs = [
            float(obj.get("area_fraction", 0.0))
            for obj in moisture_networks
            if np.isfinite(float(obj.get("area_fraction", 0.0)))
        ]
        pathways = [
            float(obj.get("mean_pathway", 0.0))
            for obj in moisture_networks
            if np.isfinite(float(obj.get("mean_pathway", 0.0)))
        ]
        detail["moisture_flow_network_object_count"] = int(len(moisture_networks))
        detail["moisture_flow_network_kind_count"] = int(len(set(kinds)))
        detail["moisture_flow_network_season_count"] = int(len(set(seasons)))
        detail["largest_moisture_flow_network_area_fraction"] = (
            max(area_fracs) if area_fracs else 0.0
        )
        detail["moisture_flow_network_mean_pathway_p50"] = (
            float(np.percentile(pathways, 50)) if pathways else 0.0
        )

    if precipitation_response_regions:
        kinds = [str(obj.get("kind", "")) for obj in precipitation_response_regions]
        seasons = [str(obj.get("season", "")) for obj in precipitation_response_regions]
        area_fracs = [
            float(obj.get("area_fraction", 0.0))
            for obj in precipitation_response_regions
            if np.isfinite(float(obj.get("area_fraction", 0.0)))
        ]
        anomalies = [
            float(obj.get("mean_abs_response_anomaly", 0.0))
            for obj in precipitation_response_regions
            if np.isfinite(float(obj.get("mean_abs_response_anomaly", 0.0)))
        ]
        source_attrs = [
            float(obj.get("source_basin_attributed_fraction", 0.0))
            for obj in precipitation_response_regions
            if np.isfinite(float(obj.get("source_basin_attributed_fraction", 0.0)))
        ]
        budget_attrs = [
            float(obj.get("budget_region_attributed_fraction", 0.0))
            for obj in precipitation_response_regions
            if np.isfinite(float(obj.get("budget_region_attributed_fraction", 0.0)))
        ]
        flow_attrs = [
            float(obj.get("flow_network_attributed_fraction", 0.0))
            for obj in precipitation_response_regions
            if np.isfinite(float(obj.get("flow_network_attributed_fraction", 0.0)))
        ]
        detail["precipitation_response_region_object_count"] = int(
            len(precipitation_response_regions))
        detail["precipitation_response_region_kind_count"] = int(len(set(kinds)))
        detail["precipitation_response_region_season_count"] = int(len(set(seasons)))
        detail["wet_precipitation_response_region_object_count"] = int(
            sum(1 for kind in kinds if kind == "wet_precipitation_response_region"))
        detail["dry_precipitation_response_region_object_count"] = int(
            sum(1 for kind in kinds if kind == "dry_precipitation_response_region"))
        detail["largest_precipitation_response_region_area_fraction"] = (
            max(area_fracs) if area_fracs else 0.0
        )
        detail["precipitation_response_region_mean_abs_anomaly_p50"] = (
            float(np.percentile(anomalies, 50)) if anomalies else 0.0
        )
        detail["precipitation_response_region_source_basin_attribution_p50"] = (
            float(np.percentile(source_attrs, 50)) if source_attrs else 0.0
        )
        detail["precipitation_response_region_budget_attribution_p50"] = (
            float(np.percentile(budget_attrs, 50)) if budget_attrs else 0.0
        )
        detail["precipitation_response_region_flow_network_attribution_p50"] = (
            float(np.percentile(flow_attrs, 50)) if flow_attrs else 0.0
        )

    if receiver_catchments:
        kinds = [str(obj.get("kind", "")) for obj in receiver_catchments]
        seasons = [str(obj.get("season", "")) for obj in receiver_catchments]
        area_fracs = [
            float(obj.get("area_fraction", 0.0))
            for obj in receiver_catchments
            if np.isfinite(float(obj.get("area_fraction", 0.0)))
        ]
        source_attrs = [
            float(obj.get("source_basin_attributed_fraction", 0.0))
            for obj in receiver_catchments
            if np.isfinite(float(obj.get("source_basin_attributed_fraction", 0.0)))
        ]
        budget_attrs = [
            float(obj.get("budget_region_attributed_fraction", 0.0))
            for obj in receiver_catchments
            if np.isfinite(float(obj.get("budget_region_attributed_fraction", 0.0)))
        ]
        precip_attrs = [
            float(obj.get("precipitation_response_attributed_fraction", 0.0))
            for obj in receiver_catchments
            if np.isfinite(float(
                obj.get("precipitation_response_attributed_fraction", 0.0)))
        ]
        supply_values = [
            float(obj.get("mean_source_basin_supply_index", 0.0))
            for obj in receiver_catchments
            if np.isfinite(float(obj.get("mean_source_basin_supply_index", 0.0)))
        ]
        supply_balances = [
            float(obj.get("precipitation_supply_balance", 0.0))
            for obj in receiver_catchments
            if np.isfinite(float(obj.get("precipitation_supply_balance", 0.0)))
        ]
        supported_precip = [
            float(obj.get("supply_supported_precipitation_fraction", 0.0))
            for obj in receiver_catchments
            if np.isfinite(float(
                obj.get("supply_supported_precipitation_fraction", 0.0)))
        ]
        detail["receiver_catchment_object_count"] = int(len(receiver_catchments))
        detail["receiver_catchment_kind_count"] = int(len(set(kinds)))
        detail["receiver_catchment_season_count"] = int(len(set(seasons)))
        detail["source_receiver_catchment_object_count"] = int(
            sum(1 for kind in kinds if kind == "source_receiver_catchment"))
        detail["mixed_receiver_catchment_object_count"] = int(
            sum(1 for kind in kinds if kind == "mixed_receiver_catchment"))
        detail["largest_receiver_catchment_area_fraction"] = (
            max(area_fracs) if area_fracs else 0.0
        )
        detail["receiver_catchment_source_basin_attribution_p50"] = (
            float(np.percentile(source_attrs, 50)) if source_attrs else 0.0
        )
        detail["receiver_catchment_budget_attribution_p50"] = (
            float(np.percentile(budget_attrs, 50)) if budget_attrs else 0.0
        )
        detail["receiver_catchment_precip_response_attribution_p50"] = (
            float(np.percentile(precip_attrs, 50)) if precip_attrs else 0.0
        )
        detail["receiver_catchment_source_supply_p50"] = (
            float(np.percentile(supply_values, 50)) if supply_values else 0.0
        )
        detail["receiver_catchment_supply_balance_p50"] = (
            float(np.percentile(supply_balances, 50)) if supply_balances else 0.0
        )
        detail["receiver_catchment_supported_precip_fraction_p50"] = (
            float(np.percentile(supported_precip, 50)) if supported_precip else 0.0
        )

    return detail


def _cryosphere_feedback_metrics(world) -> dict[str, Any]:
    grid = world.grid
    land = world.land_mask()
    ocean = ~land
    area = grid.cell_area
    n = grid.n
    seasonal_sea_ice = world.fields.get("cryosphere.seasonal_sea_ice")
    sea_ice = world.fields.get("cryosphere.sea_ice")
    seasonal_snow = world.fields.get("cryosphere.seasonal_snow")
    snow_persistence = world.fields.get("cryosphere.snow_persistence")
    seasonal_cloud = world.fields.get("climate.seasonal_cloud_albedo_proxy")
    cloud = world.fields.get("climate.cloud_albedo_proxy")
    vegetation = world.fields.get("biosphere.vegetation_climate_feedback")

    detail: dict[str, Any] = {
        "has_seasonal_sea_ice": seasonal_sea_ice is not None,
        "has_sea_ice": sea_ice is not None,
        "has_seasonal_snow": seasonal_snow is not None,
        "has_snow_persistence": snow_persistence is not None,
        "has_seasonal_cloud_albedo_proxy": seasonal_cloud is not None,
        "has_cloud_albedo_proxy": cloud is not None,
        "has_vegetation_climate_feedback": vegetation is not None,
        "invalid_seasonal_sea_ice_shape": False,
        "invalid_sea_ice_shape": False,
        "invalid_seasonal_snow_shape": False,
        "invalid_snow_persistence_shape": False,
        "invalid_seasonal_cloud_albedo_shape": False,
        "invalid_cloud_albedo_shape": False,
        "invalid_vegetation_feedback_shape": False,
        "nonfinite_seasonal_sea_ice_cells": 0,
        "nonfinite_sea_ice_cells": 0,
        "nonfinite_seasonal_snow_cells": 0,
        "nonfinite_snow_persistence_cells": 0,
        "nonfinite_seasonal_cloud_albedo_cells": 0,
        "nonfinite_cloud_albedo_cells": 0,
        "nonfinite_vegetation_feedback_cells": 0,
        "out_of_bounds_seasonal_sea_ice_cells": 0,
        "out_of_bounds_sea_ice_cells": 0,
        "out_of_bounds_seasonal_snow_cells": 0,
        "out_of_bounds_snow_persistence_cells": 0,
        "out_of_bounds_seasonal_cloud_albedo_cells": 0,
        "out_of_bounds_cloud_albedo_cells": 0,
        "out_of_bounds_vegetation_feedback_cells": 0,
        "seasonal_sea_ice_land_abs_max": 0.0,
        "sea_ice_land_abs_max": 0.0,
        "seasonal_snow_ocean_abs_max": 0.0,
        "snow_persistence_ocean_abs_max": 0.0,
        "vegetation_feedback_ocean_abs_max": 0.0,
        "seasonal_sea_ice_ocean_p95": 0.0,
        "sea_ice_ocean_p95": 0.0,
        "snow_persistence_land_p95": 0.0,
        "seasonal_snow_land_p95": 0.0,
        "cloud_albedo_proxy_p50": 0.0,
        "cloud_albedo_proxy_p95": 0.0,
        "vegetation_feedback_land_p50": 0.0,
        "vegetation_feedback_land_p95": 0.0,
        "sea_ice_adjacent_lat_band_jump_max": 0.0,
    }

    def _bounds_count(arr: np.ndarray) -> int:
        return int(((arr < -1.0e-9) | (arr > 1.0 + 1.0e-9)).sum())

    def _lat_band_jump(values: np.ndarray, mask: np.ndarray) -> float:
        bins = np.linspace(-90.0, 90.0, 19)
        means: list[float] = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            band = mask & (grid.lat >= lo) & (grid.lat < hi)
            if band.any():
                means.append(_weighted_mean(values, area, band))
            else:
                means.append(float("nan"))
        jumps = [
            abs(b - a)
            for a, b in zip(means[:-1], means[1:])
            if np.isfinite(a) and np.isfinite(b)
        ]
        return float(max(jumps, default=0.0))

    if seasonal_sea_ice is not None:
        arr = np.asarray(seasonal_sea_ice, dtype=np.float64)
        detail["invalid_seasonal_sea_ice_shape"] = bool(arr.shape != (4, n))
        detail["nonfinite_seasonal_sea_ice_cells"] = int((~np.isfinite(arr)).sum())
        detail["out_of_bounds_seasonal_sea_ice_cells"] = _bounds_count(arr)
        if arr.shape == (4, n):
            if land.any():
                detail["seasonal_sea_ice_land_abs_max"] = float(
                    np.max(np.abs(arr[:, land]), initial=0.0))
            if ocean.any():
                detail["seasonal_sea_ice_ocean_p95"] = _pct(arr[:, ocean].ravel(), 95)
                detail["sea_ice_adjacent_lat_band_jump_max"] = _lat_band_jump(
                    np.nanmean(arr, axis=0), ocean)

    if sea_ice is not None:
        arr = np.asarray(sea_ice, dtype=np.float64)
        detail["invalid_sea_ice_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_sea_ice_cells"] = int((~np.isfinite(arr)).sum())
        detail["out_of_bounds_sea_ice_cells"] = _bounds_count(arr)
        if arr.shape == (n,):
            detail["sea_ice_land_abs_max"] = float(
                np.max(np.abs(arr[land]), initial=0.0))
            detail["sea_ice_ocean_p95"] = _pct(arr[ocean], 95)

    if seasonal_snow is not None:
        arr = np.asarray(seasonal_snow, dtype=np.float64)
        detail["invalid_seasonal_snow_shape"] = bool(arr.shape != (4, n))
        detail["nonfinite_seasonal_snow_cells"] = int((~np.isfinite(arr)).sum())
        detail["out_of_bounds_seasonal_snow_cells"] = _bounds_count(arr)
        if arr.shape == (4, n):
            if ocean.any():
                detail["seasonal_snow_ocean_abs_max"] = float(
                    np.max(np.abs(arr[:, ocean]), initial=0.0))
            if land.any():
                detail["seasonal_snow_land_p95"] = _pct(arr[:, land].ravel(), 95)

    if snow_persistence is not None:
        arr = np.asarray(snow_persistence, dtype=np.float64)
        detail["invalid_snow_persistence_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_snow_persistence_cells"] = int((~np.isfinite(arr)).sum())
        detail["out_of_bounds_snow_persistence_cells"] = _bounds_count(arr)
        if arr.shape == (n,):
            detail["snow_persistence_ocean_abs_max"] = float(
                np.max(np.abs(arr[ocean]), initial=0.0))
            detail["snow_persistence_land_p95"] = _pct(arr[land], 95)

    if seasonal_cloud is not None:
        arr = np.asarray(seasonal_cloud, dtype=np.float64)
        detail["invalid_seasonal_cloud_albedo_shape"] = bool(arr.shape != (4, n))
        detail["nonfinite_seasonal_cloud_albedo_cells"] = int((~np.isfinite(arr)).sum())
        detail["out_of_bounds_seasonal_cloud_albedo_cells"] = _bounds_count(arr)

    if cloud is not None:
        arr = np.asarray(cloud, dtype=np.float64)
        detail["invalid_cloud_albedo_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_cloud_albedo_cells"] = int((~np.isfinite(arr)).sum())
        detail["out_of_bounds_cloud_albedo_cells"] = _bounds_count(arr)
        if arr.shape == (n,):
            detail["cloud_albedo_proxy_p50"] = _pct(arr, 50)
            detail["cloud_albedo_proxy_p95"] = _pct(arr, 95)

    if vegetation is not None:
        arr = np.asarray(vegetation, dtype=np.float64)
        detail["invalid_vegetation_feedback_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_vegetation_feedback_cells"] = int((~np.isfinite(arr)).sum())
        detail["out_of_bounds_vegetation_feedback_cells"] = _bounds_count(arr)
        if arr.shape == (n,):
            detail["vegetation_feedback_ocean_abs_max"] = float(
                np.max(np.abs(arr[ocean]), initial=0.0))
            detail["vegetation_feedback_land_p50"] = _pct(arr[land], 50)
            detail["vegetation_feedback_land_p95"] = _pct(arr[land], 95)

    return detail


def _circulation_metrics(world) -> dict[str, Any]:
    grid = world.grid
    n = grid.n
    seasonal_wind = world.fields.get("atmosphere.seasonal_wind")
    background_wind = world.fields.get("atmosphere.background_seasonal_wind")
    thermal_wind = world.fields.get("atmosphere.thermal_wind_anomaly")
    orographic_wind = world.fields.get("atmosphere.orographic_wind_anomaly")
    annual_wind = world.fields.get("atmosphere.wind")
    pressure = world.fields.get("atmosphere.land_sea_pressure_proxy")
    seasonal_pressure = world.fields.get("atmosphere.seasonal_pressure_proxy")
    moisture_access = world.fields.get("atmosphere.moisture_access")
    monsoon_potential = world.fields.get("atmosphere.monsoon_potential")
    source_ocean_warmth = world.fields.get("atmosphere.source_ocean_warmth")
    terrain_blocking = world.fields.get("atmosphere.terrain_blocking")
    geo_index = world.fields.get("atmosphere.geographic_circulation_index")
    itcz_lat = world.fields.get("atmosphere.itcz_latitude")
    itcz_intensity = world.fields.get("atmosphere.itcz_intensity")
    storm = world.fields.get("atmosphere.storm_track_intensity")
    detail: dict[str, Any] = {
        "has_seasonal_wind": seasonal_wind is not None,
        "has_background_seasonal_wind": background_wind is not None,
        "has_land_sea_pressure_proxy": pressure is not None,
        "has_seasonal_pressure_proxy": seasonal_pressure is not None,
        "has_moisture_access": moisture_access is not None,
        "has_monsoon_potential": monsoon_potential is not None,
        "has_source_ocean_warmth": source_ocean_warmth is not None,
        "has_terrain_blocking": terrain_blocking is not None,
        "has_thermal_wind_anomaly": thermal_wind is not None,
        "has_orographic_wind_anomaly": orographic_wind is not None,
        "has_geographic_circulation_index": geo_index is not None,
        "has_itcz_latitude": itcz_lat is not None,
        "has_itcz_intensity": itcz_intensity is not None,
        "has_storm_track_intensity": storm is not None,
        "invalid_seasonal_wind_shape": False,
        "invalid_background_seasonal_wind_shape": False,
        "invalid_land_sea_pressure_shape": False,
        "invalid_seasonal_pressure_shape": False,
        "invalid_moisture_access_shape": False,
        "invalid_monsoon_potential_shape": False,
        "invalid_source_ocean_warmth_shape": False,
        "invalid_terrain_blocking_shape": False,
        "invalid_thermal_wind_anomaly_shape": False,
        "invalid_orographic_wind_anomaly_shape": False,
        "invalid_geographic_circulation_index_shape": False,
        "invalid_itcz_latitude_shape": False,
        "invalid_itcz_intensity_shape": False,
        "invalid_storm_track_shape": False,
        "nonfinite_seasonal_wind_cells": 0,
        "nonfinite_moisture_access_cells": 0,
        "nonfinite_monsoon_potential_cells": 0,
        "nonfinite_geographic_circulation_cells": 0,
        "max_wind_normal_component": 0.0,
        "max_background_wind_normal_component": 0.0,
        "max_thermal_wind_normal_component": 0.0,
        "max_orographic_wind_normal_component": 0.0,
        "annual_wind_mean_delta": 0.0,
        "background_final_wind_delta_p95": 0.0,
        "thermal_wind_anomaly_p95_mps": 0.0,
        "orographic_wind_anomaly_p95_mps": 0.0,
        "moisture_access_land_p75": 0.0,
        "moisture_access_ocean_p50": 0.0,
        "monsoon_potential_land_p90": 0.0,
        "monsoon_potential_land_p99": 0.0,
        "monsoon_potential_global_p95": 0.0,
        "source_ocean_warmth_ocean_p75": 0.0,
        "terrain_blocking_land_p75": 0.0,
        "geographic_circulation_index_p50": 0.0,
        "geographic_circulation_index_p90": 0.0,
        "itcz_latitudes_deg": [],
        "itcz_migration_span_deg": 0.0,
        "itcz_DJF_lat_deg": 0.0,
        "itcz_JJA_lat_deg": 0.0,
        "NH_winter_storm_track_ratio": 0.0,
        "SH_winter_storm_track_ratio": 0.0,
    }

    if seasonal_wind is not None:
        arr = np.asarray(seasonal_wind, dtype=np.float64)
        detail["invalid_seasonal_wind_shape"] = bool(arr.shape != (4, n, 3))
        detail["nonfinite_seasonal_wind_cells"] = int((~np.isfinite(arr)).sum())
        if arr.shape == (4, n, 3):
            normal = np.abs(np.einsum("snj,nj->sn", arr, grid.xyz))
            detail["max_wind_normal_component"] = float(np.nanmax(normal))
            if annual_wind is not None:
                annual = np.asarray(annual_wind, dtype=np.float64)
                if annual.shape == (n, 3):
                    detail["annual_wind_mean_delta"] = float(
                        np.nanmax(np.abs(annual - arr.mean(axis=0)))
                    )

    for name, field, shape_key, normal_key, p95_key in [
        ("background", background_wind, "invalid_background_seasonal_wind_shape",
         "max_background_wind_normal_component", None),
        ("thermal", thermal_wind, "invalid_thermal_wind_anomaly_shape",
         "max_thermal_wind_normal_component", "thermal_wind_anomaly_p95_mps"),
        ("orographic", orographic_wind, "invalid_orographic_wind_anomaly_shape",
         "max_orographic_wind_normal_component", "orographic_wind_anomaly_p95_mps"),
    ]:
        if field is None:
            continue
        arr = np.asarray(field, dtype=np.float64)
        detail[shape_key] = bool(arr.shape != (4, n, 3))
        if arr.shape == (4, n, 3):
            normal = np.abs(np.einsum("snj,nj->sn", arr, grid.xyz))
            detail[normal_key] = float(np.nanmax(normal))
            if p95_key is not None:
                detail[p95_key] = float(np.nanpercentile(np.linalg.norm(arr, axis=2), 95))

    if seasonal_wind is not None and background_wind is not None:
        sw = np.asarray(seasonal_wind, dtype=np.float64)
        bg = np.asarray(background_wind, dtype=np.float64)
        if sw.shape == (4, n, 3) and bg.shape == (4, n, 3):
            detail["background_final_wind_delta_p95"] = float(
                np.nanpercentile(np.linalg.norm(sw - bg, axis=2), 95)
            )

    if pressure is not None:
        arr = np.asarray(pressure, dtype=np.float64)
        detail["invalid_land_sea_pressure_shape"] = bool(arr.shape != (4, n))

    for field, shape_key, nonfinite_key in [
        (seasonal_pressure, "invalid_seasonal_pressure_shape", None),
        (moisture_access, "invalid_moisture_access_shape",
         "nonfinite_moisture_access_cells"),
        (monsoon_potential, "invalid_monsoon_potential_shape",
         "nonfinite_monsoon_potential_cells"),
        (source_ocean_warmth, "invalid_source_ocean_warmth_shape", None),
    ]:
        if field is None:
            continue
        arr = np.asarray(field, dtype=np.float64)
        detail[shape_key] = bool(arr.shape != (4, n))
        if nonfinite_key is not None:
            detail[nonfinite_key] = int((~np.isfinite(arr)).sum())

    if terrain_blocking is not None:
        arr = np.asarray(terrain_blocking, dtype=np.float64)
        detail["invalid_terrain_blocking_shape"] = bool(arr.shape != (n,))

    land = world.land_mask()
    ocean = ~land
    if moisture_access is not None:
        arr = np.asarray(moisture_access, dtype=np.float64)
        if arr.shape == (4, n):
            if land.any():
                detail["moisture_access_land_p75"] = _pct(arr[:, land].ravel(), 75)
            if ocean.any():
                detail["moisture_access_ocean_p50"] = _pct(arr[:, ocean].ravel(), 50)
    if monsoon_potential is not None:
        arr = np.asarray(monsoon_potential, dtype=np.float64)
        if arr.shape == (4, n):
            detail["monsoon_potential_global_p95"] = _pct(arr.ravel(), 95)
            if land.any():
                land_vals = arr[:, land].ravel()
                detail["monsoon_potential_land_p90"] = _pct(land_vals, 90)
                detail["monsoon_potential_land_p99"] = _pct(land_vals, 99)
    if source_ocean_warmth is not None:
        arr = np.asarray(source_ocean_warmth, dtype=np.float64)
        if arr.shape == (4, n) and ocean.any():
            detail["source_ocean_warmth_ocean_p75"] = _pct(arr[:, ocean].ravel(), 75)
    if terrain_blocking is not None:
        arr = np.asarray(terrain_blocking, dtype=np.float64)
        if arr.shape == (n,) and land.any():
            detail["terrain_blocking_land_p75"] = _pct(arr[land], 75)

    if geo_index is not None:
        arr = np.asarray(geo_index, dtype=np.float64)
        detail["invalid_geographic_circulation_index_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_geographic_circulation_cells"] = int((~np.isfinite(arr)).sum())
        if arr.shape == (n,):
            detail["geographic_circulation_index_p50"] = _pct(arr[np.isfinite(arr)], 50)
            detail["geographic_circulation_index_p90"] = _pct(arr[np.isfinite(arr)], 90)

    if itcz_lat is not None:
        arr = np.asarray(itcz_lat, dtype=np.float64)
        detail["invalid_itcz_latitude_shape"] = bool(arr.shape != (4,))
        if arr.shape == (4,):
            detail["itcz_latitudes_deg"] = [round(float(x), 3) for x in arr]
            detail["itcz_migration_span_deg"] = float(np.nanmax(arr) - np.nanmin(arr))
            detail["itcz_DJF_lat_deg"] = float(arr[0])
            detail["itcz_JJA_lat_deg"] = float(arr[2])

    if itcz_intensity is not None:
        arr = np.asarray(itcz_intensity, dtype=np.float64)
        detail["invalid_itcz_intensity_shape"] = bool(arr.shape != (4, n))

    if storm is not None:
        arr = np.asarray(storm, dtype=np.float64)
        detail["invalid_storm_track_shape"] = bool(arr.shape != (4, n))
        if arr.shape == (4, n):
            nh = (grid.lat >= 35.0) & (grid.lat <= 60.0)
            sh = (grid.lat <= -35.0) & (grid.lat >= -60.0)
            nh_djf = _weighted_mean(arr[0], grid.cell_area, nh)
            nh_jja = _weighted_mean(arr[2], grid.cell_area, nh)
            sh_djf = _weighted_mean(arr[0], grid.cell_area, sh)
            sh_jja = _weighted_mean(arr[2], grid.cell_area, sh)
            detail["NH_winter_storm_track_ratio"] = nh_djf / max(nh_jja, 1e-9)
            detail["SH_winter_storm_track_ratio"] = sh_jja / max(sh_djf, 1e-9)

    return detail


def _ocean_current_metrics(world) -> dict[str, Any]:
    grid = world.grid
    n = grid.n
    ocean = world.ocean_mask()
    land = ~ocean
    currents = world.fields.get("ocean.currents")
    heat = world.fields.get("ocean.current_heat_transport")
    upwelling = world.fields.get("ocean.upwelling")
    basin_id = world.fields.get("ocean.basin_id")
    seasonal_sst = world.fields.get("climate.seasonal_sst")
    ocean_heat_flux = world.fields.get("climate.ocean_heat_flux")
    coupling_residual = world.fields.get("climate.coupling_residual")
    streamfunction = world.fields.get("ocean.current_streamfunction")
    gyre_id = world.fields.get("ocean.gyre_id")
    boundary_current_type = world.fields.get("ocean.boundary_current_type")
    strait_exchange = world.fields.get("ocean.strait_exchange")
    wind_stress_response = world.fields.get("ocean.wind_stress_current_response")
    sst_anomaly = world.fields.get("ocean.sst_anomaly")
    solved_mask = world.fields.get("ocean.solved_mask")
    solved_ocean = ocean
    if solved_mask is not None:
        arr_mask = np.asarray(solved_mask, dtype=np.float64)
        if arr_mask.shape == (n,):
            solved_ocean = arr_mask > 0.5
    solved_land = ~solved_ocean
    detail: dict[str, Any] = {
        "has_currents": currents is not None,
        "has_current_heat_transport": heat is not None,
        "has_upwelling": upwelling is not None,
        "has_basin_id": basin_id is not None,
        "has_seasonal_sst": seasonal_sst is not None,
        "has_ocean_heat_flux": ocean_heat_flux is not None,
        "has_coupling_residual": coupling_residual is not None,
        "has_current_streamfunction": streamfunction is not None,
        "has_gyre_id": gyre_id is not None,
        "has_boundary_current_type": boundary_current_type is not None,
        "has_strait_exchange": strait_exchange is not None,
        "has_wind_stress_current_response": wind_stress_response is not None,
        "has_sst_anomaly": sst_anomaly is not None,
        "has_solved_mask": solved_mask is not None,
        "invalid_current_shape": False,
        "invalid_current_heat_transport_shape": False,
        "invalid_upwelling_shape": False,
        "invalid_basin_id_shape": False,
        "invalid_seasonal_sst_shape": False,
        "invalid_ocean_heat_flux_shape": False,
        "invalid_coupling_residual_shape": False,
        "invalid_current_streamfunction_shape": False,
        "invalid_gyre_id_shape": False,
        "invalid_boundary_current_type_shape": False,
        "invalid_strait_exchange_shape": False,
        "invalid_wind_stress_current_response_shape": False,
        "invalid_sst_anomaly_shape": False,
        "invalid_solved_mask_shape": False,
        "nonfinite_current_cells": 0,
        "nonfinite_current_heat_transport_cells": 0,
        "nonfinite_upwelling_cells": 0,
        "nonfinite_seasonal_sst_cells": 0,
        "nonfinite_ocean_heat_flux_cells": 0,
        "nonfinite_coupling_residual_cells": 0,
        "nonfinite_current_streamfunction_cells": 0,
        "nonfinite_wind_stress_current_response_cells": 0,
        "nonfinite_sst_anomaly_cells": 0,
        "current_over_solved_land_speed_max_mps": 0.0,
        "current_over_solved_land_cells_gt_1e_8": 0,
        "current_over_final_land_speed_max_mps": 0.0,
        "current_over_final_land_cells_gt_1e_8": 0,
        "solved_final_ocean_mask_mismatch_fraction": 0.0,
        "max_current_normal_component": 0.0,
        "ocean_current_speed_p50_mps": 0.0,
        "ocean_current_speed_p95_mps": 0.0,
        "current_heat_transport_p05_C": 0.0,
        "current_heat_transport_p95_C": 0.0,
        "current_heat_transport_abs_p95_C": 0.0,
        "ocean_heat_transport_abs_p95_C": 0.0,
        "land_heat_transport_abs_p95_C": 0.0,
        "warm_current_ocean_fraction_gt_0_5C": 0.0,
        "cold_current_ocean_fraction_lt_minus_0_5C": 0.0,
        "upwelling_p95": 0.0,
        "upwelling_ocean_fraction_gt_0_25": 0.0,
        "ocean_basin_count": 0,
        "largest_ocean_basin_fraction": 0.0,
        "seasonal_sst_ocean_mean_C": 0.0,
        "seasonal_sst_zonal_residual_abs_p95_C": 0.0,
        "ocean_heat_flux_ocean_mean_C": 0.0,
        "ocean_heat_flux_ocean_abs_p95_C": 0.0,
        "coupling_residual_p95": 0.0,
        "current_streamfunction_land_abs_max": 0.0,
        "current_streamfunction_ocean_abs_p95": 0.0,
        "gyre_count": 0,
        "boundary_current_ocean_fraction": 0.0,
        "strait_exchange_p95": 0.0,
        "wind_stress_response_ocean_p50_mps": 0.0,
        "wind_stress_response_ocean_p95_mps": 0.0,
        "wind_stress_response_land_speed_max_mps": 0.0,
        "wind_stress_response_max_normal_component": 0.0,
        "sst_anomaly_ocean_mean_C": 0.0,
        "sst_anomaly_ocean_abs_p95_C": 0.0,
    }

    if solved_mask is not None:
        arr = np.asarray(solved_mask, dtype=np.float64)
        detail["invalid_solved_mask_shape"] = bool(arr.shape != (n,))
        if arr.shape == (n,):
            detail["solved_final_ocean_mask_mismatch_fraction"] = _area_frac(
                grid.cell_area, solved_ocean != ocean)

    if currents is not None:
        arr = np.asarray(currents, dtype=np.float64)
        detail["invalid_current_shape"] = bool(arr.shape != (n, 3))
        detail["nonfinite_current_cells"] = int((~np.isfinite(arr)).sum())
        if arr.shape == (n, 3):
            speed = np.linalg.norm(arr, axis=1)
            detail["current_over_solved_land_speed_max_mps"] = (
                float(np.nanmax(speed[solved_land])) if solved_land.any() else 0.0
            )
            detail["current_over_solved_land_cells_gt_1e_8"] = int(
                (solved_land & (speed > 1e-8)).sum())
            detail["current_over_final_land_speed_max_mps"] = (
                float(np.nanmax(speed[land])) if land.any() else 0.0
            )
            detail["current_over_final_land_cells_gt_1e_8"] = int(
                (land & (speed > 1e-8)).sum())
            normal = np.abs(np.einsum("nj,nj->n", arr, grid.xyz))
            detail["max_current_normal_component"] = float(np.nanmax(normal))
            if ocean.any():
                detail["ocean_current_speed_p50_mps"] = _pct(speed[ocean], 50)
                detail["ocean_current_speed_p95_mps"] = _pct(speed[ocean], 95)

    if heat is not None:
        arr = np.asarray(heat, dtype=np.float64)
        detail["invalid_current_heat_transport_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_current_heat_transport_cells"] = int((~np.isfinite(arr)).sum())
        if arr.shape == (n,):
            finite = arr[np.isfinite(arr)]
            detail["current_heat_transport_p05_C"] = _pct(finite, 5)
            detail["current_heat_transport_p95_C"] = _pct(finite, 95)
            detail["current_heat_transport_abs_p95_C"] = _pct(np.abs(finite), 95)
            if ocean.any():
                detail["ocean_heat_transport_abs_p95_C"] = _pct(np.abs(arr[ocean]), 95)
                detail["warm_current_ocean_fraction_gt_0_5C"] = _area_frac(
                    grid.cell_area, ocean & (arr > 0.5), ocean)
                detail["cold_current_ocean_fraction_lt_minus_0_5C"] = _area_frac(
                    grid.cell_area, ocean & (arr < -0.5), ocean)
            if land.any():
                detail["land_heat_transport_abs_p95_C"] = _pct(np.abs(arr[land]), 95)

    if upwelling is not None:
        arr = np.asarray(upwelling, dtype=np.float64)
        detail["invalid_upwelling_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_upwelling_cells"] = int((~np.isfinite(arr)).sum())
        if arr.shape == (n,):
            detail["upwelling_p95"] = _pct(arr[ocean], 95) if ocean.any() else 0.0
            detail["upwelling_ocean_fraction_gt_0_25"] = _area_frac(
                grid.cell_area, ocean & (arr > 0.25), ocean)

    if basin_id is not None:
        arr = np.asarray(basin_id)
        detail["invalid_basin_id_shape"] = bool(arr.shape != (n,))
        if arr.shape == (n,) and ocean.any():
            ids = [int(x) for x in np.unique(arr[ocean].astype(int)) if int(x) >= 0]
            detail["ocean_basin_count"] = len(ids)
            if ids:
                areas = [
                    float(grid.cell_area[ocean & (arr.astype(int) == bid)].sum())
                    for bid in ids
                ]
                ocean_area = float(grid.cell_area[ocean].sum())
                detail["largest_ocean_basin_fraction"] = max(areas) / ocean_area if ocean_area else 0.0

    if seasonal_sst is not None:
        arr = np.asarray(seasonal_sst, dtype=np.float64)
        detail["invalid_seasonal_sst_shape"] = bool(arr.shape != (4, n))
        detail["nonfinite_seasonal_sst_cells"] = int((~np.isfinite(arr)).sum())
        if arr.shape == (4, n) and ocean.any():
            sst_c = arr.mean(axis=0) - CONSTANTS.ZERO_C
            detail["seasonal_sst_ocean_mean_C"] = _weighted_mean(
                sst_c, grid.cell_area, ocean)
            detail["seasonal_sst_zonal_residual_abs_p95_C"] = (
                _lat_band_residual_abs_p95(grid, sst_c, ocean)
            )

    if ocean_heat_flux is not None:
        arr = np.asarray(ocean_heat_flux, dtype=np.float64)
        detail["invalid_ocean_heat_flux_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_ocean_heat_flux_cells"] = int((~np.isfinite(arr)).sum())
        if arr.shape == (n,) and ocean.any():
            detail["ocean_heat_flux_ocean_mean_C"] = _weighted_mean(
                arr, grid.cell_area, ocean)
            detail["ocean_heat_flux_ocean_abs_p95_C"] = _pct(np.abs(arr[ocean]), 95)

    if coupling_residual is not None:
        arr = np.asarray(coupling_residual, dtype=np.float64)
        detail["invalid_coupling_residual_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_coupling_residual_cells"] = int((~np.isfinite(arr)).sum())
        if arr.shape == (n,):
            detail["coupling_residual_p95"] = _pct(arr[ocean], 95) if ocean.any() else 0.0

    if streamfunction is not None:
        arr = np.asarray(streamfunction, dtype=np.float64)
        detail["invalid_current_streamfunction_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_current_streamfunction_cells"] = int((~np.isfinite(arr)).sum())
        if arr.shape == (n,):
            detail["current_streamfunction_land_abs_max"] = (
                float(np.nanmax(np.abs(arr[land]))) if land.any() else 0.0
            )
            detail["current_streamfunction_ocean_abs_p95"] = (
                _pct(np.abs(arr[ocean]), 95) if ocean.any() else 0.0
            )

    if gyre_id is not None:
        arr = np.asarray(gyre_id)
        detail["invalid_gyre_id_shape"] = bool(arr.shape != (n,))
        if arr.shape == (n,) and ocean.any():
            ids = [int(x) for x in np.unique(arr[ocean].astype(int)) if int(x) > 0]
            detail["gyre_count"] = len(ids)

    if boundary_current_type is not None:
        arr = np.asarray(boundary_current_type, dtype=np.float64)
        detail["invalid_boundary_current_type_shape"] = bool(arr.shape != (n,))
        if arr.shape == (n,):
            detail["boundary_current_ocean_fraction"] = _area_frac(
                grid.cell_area, ocean & (np.abs(arr) > 0.5), ocean)

    if strait_exchange is not None:
        arr = np.asarray(strait_exchange, dtype=np.float64)
        detail["invalid_strait_exchange_shape"] = bool(arr.shape != (n,))
        if arr.shape == (n,):
            detail["strait_exchange_p95"] = _pct(arr[ocean], 95) if ocean.any() else 0.0

    if wind_stress_response is not None:
        arr = np.asarray(wind_stress_response, dtype=np.float64)
        detail["invalid_wind_stress_current_response_shape"] = bool(
            arr.shape != (n, 3))
        detail["nonfinite_wind_stress_current_response_cells"] = int(
            (~np.isfinite(arr)).sum())
        if arr.shape == (n, 3):
            speed = np.linalg.norm(arr, axis=1)
            detail["wind_stress_response_land_speed_max_mps"] = (
                float(np.nanmax(speed[land])) if land.any() else 0.0
            )
            detail["wind_stress_response_ocean_p50_mps"] = (
                _pct(speed[ocean], 50) if ocean.any() else 0.0
            )
            detail["wind_stress_response_ocean_p95_mps"] = (
                _pct(speed[ocean], 95) if ocean.any() else 0.0
            )
            normal = np.abs(np.einsum("nj,nj->n", arr, grid.xyz))
            detail["wind_stress_response_max_normal_component"] = float(
                np.nanmax(normal))

    if sst_anomaly is not None:
        arr = np.asarray(sst_anomaly, dtype=np.float64)
        detail["invalid_sst_anomaly_shape"] = bool(arr.shape != (n,))
        detail["nonfinite_sst_anomaly_cells"] = int((~np.isfinite(arr)).sum())
        if arr.shape == (n,) and ocean.any():
            detail["sst_anomaly_ocean_mean_C"] = _weighted_mean(
                arr, grid.cell_area, ocean)
            detail["sst_anomaly_ocean_abs_p95_C"] = _pct(np.abs(arr[ocean]), 95)

    return detail


def _geography_primitive_metrics(world) -> dict[str, Any]:
    grid = world.grid
    n = grid.n
    final_ocean = world.ocean_mask()
    solved_mask = world.fields.get("ocean.solved_mask")
    ocean = final_ocean
    if solved_mask is not None:
        solved = np.asarray(solved_mask, dtype=np.float64)
        if solved.shape == (n,):
            ocean = solved > 0.5
    land = ~ocean

    continent_id = world.fields.get("climate.continent_id")
    interiority = world.fields.get("climate.continent_interiority")
    coast_orientation = world.fields.get("climate.coast_orientation")
    coast_distance = world.fields.get("climate.coast_distance")
    coast_strength = world.fields.get("climate.coast_strength")
    coast_facing = world.fields.get("climate.coast_facing_east")
    basin_id = world.fields.get("ocean.basin_id")
    shelf = world.fields.get("ocean.shelf_index")
    strait = world.fields.get("ocean.strait_index")
    barrier = world.fields.get("terrain.barrier_index")
    wind_gap = world.fields.get("terrain.wind_gap_index")

    land_coast = np.zeros(n, dtype=bool)
    ocean_coast = np.zeros(n, dtype=bool)
    for c in range(n):
        nbs = grid.neighbors[c]
        if land[c] and nbs.size and ocean[nbs].any():
            land_coast[c] = True
        if ocean[c] and nbs.size and land[nbs].any():
            ocean_coast[c] = True

    detail: dict[str, Any] = {
        "has_continent_id": continent_id is not None,
        "has_continent_interiority": interiority is not None,
        "has_coast_orientation": coast_orientation is not None,
        "has_coast_distance": coast_distance is not None,
        "has_coast_strength": coast_strength is not None,
        "has_coast_facing_east": coast_facing is not None,
        "has_basin_id": basin_id is not None,
        "has_shelf_index": shelf is not None,
        "has_strait_index": strait is not None,
        "has_barrier_index": barrier is not None,
        "has_wind_gap_index": wind_gap is not None,
        "invalid_continent_id_shape": False,
        "invalid_continent_interiority_shape": False,
        "invalid_coast_orientation_shape": False,
        "invalid_coast_distance_shape": False,
        "invalid_coast_strength_shape": False,
        "invalid_coast_facing_east_shape": False,
        "invalid_basin_id_shape": False,
        "invalid_shelf_index_shape": False,
        "invalid_strait_index_shape": False,
        "invalid_barrier_index_shape": False,
        "invalid_wind_gap_index_shape": False,
        "nonfinite_geography_cells": 0,
        "continent_count": 0,
        "largest_continent_fraction_of_land": 0.0,
        "land_without_continent_id_cells": 0,
        "ocean_with_continent_id_cells": 0,
        "ocean_basin_count": 0,
        "largest_basin_fraction_of_ocean": 0.0,
        "ocean_without_basin_id_cells": 0,
        "land_with_basin_id_cells": 0,
        "coastal_land_cells": int(land_coast.sum()),
        "ocean_coast_cells": int(ocean_coast.sum()),
        "max_coast_orientation_normal_component": 0.0,
        "coast_orientation_coastal_speed_p50": 0.0,
        "coast_orientation_noncoastal_speed_p95": 0.0,
        "coast_strength_land_noncoast_max": 0.0,
        "coast_distance_min": 0.0,
        "coast_distance_max": 0.0,
        "shelf_near_coast_p75": 0.0,
        "shelf_deep_ocean_p75": 0.0,
        "shelf_contrast_near_minus_deep": 0.0,
        "strait_land_max": 0.0,
        "strait_ocean_p95": 0.0,
        "barrier_high_relief_p75": 0.0,
        "barrier_lowland_p75": 0.0,
        "barrier_contrast_high_minus_low": 0.0,
        "wind_gap_p95": 0.0,
        "continent_object_count": len(world.objects.get("climate.continents", [])),
        "ocean_basin_object_count": len(world.objects.get("ocean.basins", [])),
        "coastline_segment_count": len(world.objects.get("climate.coastline_segments", [])),
        "strait_object_count": len(world.objects.get("ocean.straits", [])),
        "barrier_belt_object_count": len(world.objects.get("terrain.barrier_belts", [])),
    }

    scalar_fields = [
        continent_id, interiority, coast_distance, coast_strength, coast_facing,
        basin_id, shelf, strait, barrier, wind_gap,
    ]
    for arr in scalar_fields:
        if arr is not None:
            arr = np.asarray(arr, dtype=np.float64)
            detail["nonfinite_geography_cells"] += int((~np.isfinite(arr)).sum())
    if coast_orientation is not None:
        arr = np.asarray(coast_orientation, dtype=np.float64)
        detail["nonfinite_geography_cells"] += int((~np.isfinite(arr)).sum())

    if continent_id is not None:
        arr = np.asarray(continent_id)
        detail["invalid_continent_id_shape"] = bool(arr.shape != (n,))
        if arr.shape == (n,):
            ids = [int(x) for x in np.unique(arr[land].astype(int)) if int(x) >= 0]
            detail["continent_count"] = len(ids)
            detail["land_without_continent_id_cells"] = int((land & (arr < 0)).sum())
            detail["ocean_with_continent_id_cells"] = int((ocean & (arr >= 0)).sum())
            if ids and land.any():
                areas = [float(grid.cell_area[land & (arr.astype(int) == cid)].sum())
                         for cid in ids]
                land_area = float(grid.cell_area[land].sum())
                detail["largest_continent_fraction_of_land"] = max(areas) / land_area

    if basin_id is not None:
        arr = np.asarray(basin_id)
        detail["invalid_basin_id_shape"] = bool(arr.shape != (n,))
        if arr.shape == (n,):
            ids = [int(x) for x in np.unique(arr[ocean].astype(int)) if int(x) >= 0]
            detail["ocean_basin_count"] = len(ids)
            detail["ocean_without_basin_id_cells"] = int((ocean & (arr < 0)).sum())
            detail["land_with_basin_id_cells"] = int((land & (arr >= 0)).sum())
            if ids and ocean.any():
                areas = [float(grid.cell_area[ocean & (arr.astype(int) == bid)].sum())
                         for bid in ids]
                ocean_area = float(grid.cell_area[ocean].sum())
                detail["largest_basin_fraction_of_ocean"] = max(areas) / ocean_area

    for field, key in [
        (interiority, "invalid_continent_interiority_shape"),
        (coast_distance, "invalid_coast_distance_shape"),
        (coast_strength, "invalid_coast_strength_shape"),
        (coast_facing, "invalid_coast_facing_east_shape"),
        (shelf, "invalid_shelf_index_shape"),
        (strait, "invalid_strait_index_shape"),
        (barrier, "invalid_barrier_index_shape"),
        (wind_gap, "invalid_wind_gap_index_shape"),
    ]:
        if field is not None:
            detail[key] = bool(np.asarray(field).shape != (n,))

    if coast_orientation is not None:
        arr = np.asarray(coast_orientation, dtype=np.float64)
        detail["invalid_coast_orientation_shape"] = bool(arr.shape != (n, 3))
        if arr.shape == (n, 3):
            speed = np.linalg.norm(arr, axis=1)
            normal = np.abs(np.einsum("ij,ij->i", arr, grid.xyz))
            detail["max_coast_orientation_normal_component"] = float(np.nanmax(normal))
            detail["coast_orientation_coastal_speed_p50"] = _pct(speed[land_coast], 50)
            noncoast = ~land_coast
            detail["coast_orientation_noncoastal_speed_p95"] = _pct(speed[noncoast], 95)

    if coast_strength is not None:
        arr = np.asarray(coast_strength, dtype=np.float64)
        if arr.shape == (n,):
            detail["coast_strength_land_noncoast_max"] = (
                float(np.nanmax(arr[land & ~land_coast])) if (land & ~land_coast).any()
                else 0.0
            )

    if coast_distance is not None:
        arr = np.asarray(coast_distance, dtype=np.float64)
        if arr.shape == (n,):
            detail["coast_distance_min"] = float(np.nanmin(arr))
            detail["coast_distance_max"] = float(np.nanmax(arr))

    if shelf is not None:
        arr = np.asarray(shelf, dtype=np.float64)
        if arr.shape == (n,):
            deep = ocean & (coast_distance is not None)
            if coast_distance is not None and np.asarray(coast_distance).shape == (n,):
                cd = np.asarray(coast_distance, dtype=np.float64)
                deep = ocean & (cd > np.percentile(cd[ocean], 65) if ocean.any() else False)
            near = ocean_coast
            detail["shelf_near_coast_p75"] = _pct(arr[near], 75)
            detail["shelf_deep_ocean_p75"] = _pct(arr[deep], 75)
            detail["shelf_contrast_near_minus_deep"] = (
                detail["shelf_near_coast_p75"] - detail["shelf_deep_ocean_p75"]
            )

    if strait is not None:
        arr = np.asarray(strait, dtype=np.float64)
        if arr.shape == (n,):
            detail["strait_land_max"] = float(np.nanmax(arr[land])) if land.any() else 0.0
            detail["strait_ocean_p95"] = _pct(arr[ocean], 95)

    if barrier is not None:
        arr = np.asarray(barrier, dtype=np.float64)
        elev = world.get_field("terrain.elevation_m", 0.0) - world.sea_level
        if arr.shape == (n,) and land.any():
            high_cut = float(np.percentile(elev[land], 80))
            low_cut = float(np.percentile(elev[land], 35))
            high = land & (elev >= high_cut)
            low = land & (elev <= low_cut)
            detail["barrier_high_relief_p75"] = _pct(arr[high], 75)
            detail["barrier_lowland_p75"] = _pct(arr[low], 75)
            detail["barrier_contrast_high_minus_low"] = (
                detail["barrier_high_relief_p75"] - detail["barrier_lowland_p75"]
            )

    if wind_gap is not None:
        arr = np.asarray(wind_gap, dtype=np.float64)
        if arr.shape == (n,):
            detail["wind_gap_p95"] = _pct(arr[land], 95)

    return detail


def _ocean_geography_metrics(world) -> dict[str, Any]:
    grid = world.grid
    n = grid.n
    area = grid.cell_area
    rel = world.get_field("terrain.elevation_m", 0.0) - world.sea_level
    ocean = rel < 0.0
    land = ~ocean
    basin_id = world.fields.get("ocean.basin_id")
    margin_type = world.fields.get("ocean.margin_type")
    depth_province = world.fields.get("ocean.depth_province")
    gateway_id = world.fields.get("ocean.gateway_id")
    shelf_width = world.fields.get("ocean.shelf_width")

    detail: dict[str, Any] = {
        "has_basin_id": basin_id is not None,
        "has_margin_type": margin_type is not None,
        "has_depth_province": depth_province is not None,
        "has_gateway_id": gateway_id is not None,
        "has_shelf_width": shelf_width is not None,
        "invalid_basin_id_shape": False,
        "invalid_margin_type_shape": False,
        "invalid_depth_province_shape": False,
        "invalid_gateway_id_shape": False,
        "invalid_shelf_width_shape": False,
        "nonfinite_ocean_geography_cells": 0,
        "land_with_depth_province_cells": 0,
        "ocean_without_depth_province_cells": 0,
        "land_with_margin_type_cells": 0,
        "ocean_without_margin_type_cells": 0,
        "land_with_gateway_id_cells": 0,
        "land_with_shelf_width_cells": 0,
        "ocean_basin_count": 0,
        "ocean_component_count": 0,
        "largest_ocean_component_fraction_of_ocean": 0.0,
        "max_closed_ocean_candidate_fraction_of_ocean": 0.0,
        "max_closed_ocean_candidate_fraction_world": 0.0,
        "closed_ocean_candidate_fraction_of_ocean": 0.0,
        "closed_ocean_candidate_count_gt2pct_world": 0,
        "closed_ocean_candidate_count_gt5pct_ocean": 0,
        "ocean_basin_object_count": len(world.objects.get("ocean.basins", [])),
        "ocean_margin_object_count": len(world.objects.get("ocean.margins", [])),
        "ocean_gateway_object_count": len(world.objects.get("ocean.gateways", [])),
        "shelf_fraction_of_ocean": 0.0,
        "slope_rise_fraction_of_ocean": 0.0,
        "abyss_fraction_of_ocean": 0.0,
        "ridge_fraction_of_ocean": 0.0,
        "trench_fraction_of_ocean": 0.0,
        "restricted_fraction_of_ocean": 0.0,
        "nearshore_depth_p75_m": 0.0,
        "shelf_depth_p75_m": 0.0,
        "abyss_depth_p50_m": 0.0,
        "trench_depth_p50_m": 0.0,
        "shelf_to_abyss_depth_delta_m": 0.0,
        "nearshore_superdeep_fraction_gt2500m": 0.0,
        "far_ocean_shallow_fraction_lt1500m": 0.0,
        "trench_near_active_margin_fraction": 0.0,
        "ocean_fabric_object_count": 0,
        "ocean_fabric_kind_counts": {},
        "spreading_center_object_count": 0,
        "transform_fault_object_count": 0,
        "fracture_zone_object_count": 0,
        "abyssal_plain_object_count": 0,
        "abyssal_hill_object_count": 0,
        "age_isochron_object_count": 0,
        "ridge_transform_fracture_core_present": False,
        "fracture_zone_combined_span_deg": 0.0,
        "transform_fault_combined_span_deg": 0.0,
        "arc_plume_landform_object_count": 0,
        "arc_plume_landform_kind_counts": {},
        "hotspot_track_object_count": 0,
        "seamount_chain_object_count": 0,
        "oceanic_plateau_object_count": 0,
        "back_arc_basin_object_count": 0,
        "island_arc_object_count": 0,
        "accreted_terrane_object_count": 0,
        "microcontinent_object_count": 0,
        "margin_landform_kind_counts": {},
        "passive_margin_wedge_object_count": 0,
        "delta_fan_object_count": 0,
        "margin_trench_object_count": 0,
        "p106_object_category_count": 0,
        "p106_core_ocean_fabric_complete": False,
        "p106_full_ocean_feature_coverage": False,
        "far_ocean_shallow_component_count_lt1500m": 0,
        "unparented_far_ocean_shallow_component_count_lt1500m": 0,
        "isolated_far_ocean_deep_pit_component_count_gt4100m": 0,
        "parented_archipelago_component_count": 0,
        "parented_archipelago_area_fraction_world": 0.0,
        "unparented_oceanic_island_component_count": 0,
    }

    if ocean.any():
        ocean_area = max(float(area[ocean].sum()), 1.0e-12)
        total_area = max(float(area.sum()), 1.0e-12)
        component_areas = sorted(_component_areas(grid, ocean), reverse=True)
        secondary = component_areas[1:]
        detail["ocean_component_count"] = int(len(component_areas))
        detail["largest_ocean_component_fraction_of_ocean"] = (
            float(component_areas[0] / ocean_area) if component_areas else 0.0
        )
        detail["max_closed_ocean_candidate_fraction_of_ocean"] = (
            float(max(secondary) / ocean_area) if secondary else 0.0
        )
        detail["max_closed_ocean_candidate_fraction_world"] = (
            float(max(secondary) / total_area) if secondary else 0.0
        )
        detail["closed_ocean_candidate_fraction_of_ocean"] = (
            float(sum(secondary) / ocean_area) if secondary else 0.0
        )
        detail["closed_ocean_candidate_count_gt2pct_world"] = int(
            sum(1 for value in secondary if value / total_area >= 0.020)
        )
        detail["closed_ocean_candidate_count_gt5pct_ocean"] = int(
            sum(1 for value in secondary if value / ocean_area >= 0.050)
        )

    for field, key in [
        (basin_id, "invalid_basin_id_shape"),
        (margin_type, "invalid_margin_type_shape"),
        (depth_province, "invalid_depth_province_shape"),
        (gateway_id, "invalid_gateway_id_shape"),
        (shelf_width, "invalid_shelf_width_shape"),
    ]:
        if field is not None:
            arr = np.asarray(field, dtype=np.float64)
            detail[key] = bool(arr.shape != (n,))
            detail["nonfinite_ocean_geography_cells"] += int((~np.isfinite(arr)).sum())

    if basin_id is not None:
        arr = np.asarray(basin_id)
        if arr.shape == (n,) and ocean.any():
            detail["ocean_basin_count"] = len(
                [int(x) for x in np.unique(arr[ocean].astype(int)) if int(x) >= 0])

    if depth_province is not None:
        prov = np.asarray(depth_province, dtype=int)
        if prov.shape == (n,):
            detail["land_with_depth_province_cells"] = int((land & (prov != 0)).sum())
            detail["ocean_without_depth_province_cells"] = int((ocean & (prov == 0)).sum())
            if ocean.any():
                detail["shelf_fraction_of_ocean"] = _area_frac(area, ocean & (prov == 1), ocean)
                detail["slope_rise_fraction_of_ocean"] = _area_frac(
                    area, ocean & np.isin(prov, [2, 3]), ocean)
                detail["abyss_fraction_of_ocean"] = _area_frac(area, ocean & (prov == 4), ocean)
                detail["ridge_fraction_of_ocean"] = _area_frac(area, ocean & (prov == 5), ocean)
                detail["trench_fraction_of_ocean"] = _area_frac(area, ocean & (prov == 6), ocean)
                detail["restricted_fraction_of_ocean"] = _area_frac(
                    area, ocean & (prov == 7), ocean)
            depth = np.maximum(-rel, 0.0)
            shelf = ocean & (prov == 1)
            abyss = ocean & (prov == 4)
            trench = ocean & (prov == 6)
            detail["shelf_depth_p75_m"] = _pct(depth[shelf], 75)
            detail["abyss_depth_p50_m"] = _pct(depth[abyss], 50)
            detail["trench_depth_p50_m"] = _pct(depth[trench], 50)
            detail["shelf_to_abyss_depth_delta_m"] = (
                detail["abyss_depth_p50_m"] - detail["shelf_depth_p75_m"]
            )
            if shelf_width is not None and np.asarray(shelf_width).shape == (n,):
                width = np.asarray(shelf_width, dtype=np.float64)
                nearshore = ocean & (width == 1.0) & (prov != 6)
                far = ocean & ((width <= 0.0) | (width >= 4.0)) & (prov != 6)
                detail["nearshore_depth_p75_m"] = _pct(depth[nearshore], 75)
                detail["nearshore_superdeep_fraction_gt2500m"] = _area_frac(
                    area, nearshore & (depth > 2500.0), nearshore)
                detail["far_ocean_shallow_fraction_lt1500m"] = _area_frac(
                    area, far & (depth < 1500.0), far)

            boundaries = world.networks.get("tectonics.boundaries", {})
            active = np.zeros(n, dtype=bool)
            for name in ("active_margin", "trench", "subduction"):
                cells = np.asarray(boundaries.get(name, []), dtype=int)
                cells = cells[(0 <= cells) & (cells < n)]
                active[cells] = True
            active_zone = _dilate_mask(grid, active, np.ones(n, dtype=bool), passes=2)
            detail["trench_near_active_margin_fraction"] = _area_frac(
                area, trench & active_zone, trench)

    if margin_type is not None:
        arr = np.asarray(margin_type, dtype=int)
        if arr.shape == (n,):
            detail["land_with_margin_type_cells"] = int((land & (arr != 0)).sum())
            detail["ocean_without_margin_type_cells"] = int((ocean & (arr == 0)).sum())

    if gateway_id is not None:
        arr = np.asarray(gateway_id, dtype=np.float64)
        if arr.shape == (n,):
            detail["land_with_gateway_id_cells"] = int((land & (arr >= 0.0)).sum())

    if shelf_width is not None:
        arr = np.asarray(shelf_width, dtype=np.float64)
        if arr.shape == (n,):
            detail["land_with_shelf_width_cells"] = int((land & (arr > 0.0)).sum())

    def kind_counts(objects: list[dict[str, Any]]) -> dict[str, int]:
        counts: dict[str, int] = {}
        for obj in objects:
            kind = str(obj.get("kind", ""))
            if not kind:
                continue
            counts[kind] = counts.get(kind, 0) + 1
        return counts

    def combined_span(objects: list[dict[str, Any]], kinds: set[str]) -> float:
        cells: list[int] = []
        span_from_objects = 0.0
        for obj in objects:
            if str(obj.get("kind", "")) not in kinds:
                continue
            span_from_objects = max(
                span_from_objects,
                float(obj.get("lon_span_deg", 0.0) or 0.0),
                float(obj.get("lat_span_deg", 0.0) or 0.0),
            )
            obj_cells = np.asarray(obj.get("cells", []), dtype=int)
            obj_cells = obj_cells[(0 <= obj_cells) & (obj_cells < n)]
            if obj_cells.size:
                cells.extend(int(c) for c in obj_cells)
        if cells:
            arr = np.asarray(cells, dtype=int)
            return float(max(np.ptp(grid.lon[arr]), np.ptp(grid.lat[arr])))
        return float(span_from_objects)

    ocean_fabric = list(world.objects.get("terrain.ocean_fabric", []))
    ocean_fabric_counts = kind_counts(ocean_fabric)
    detail["ocean_fabric_object_count"] = int(len(ocean_fabric))
    detail["ocean_fabric_kind_counts"] = ocean_fabric_counts
    detail["spreading_center_object_count"] = int(
        ocean_fabric_counts.get("spreading_center", 0))
    detail["transform_fault_object_count"] = int(
        ocean_fabric_counts.get("transform_fault", 0))
    detail["fracture_zone_object_count"] = int(
        ocean_fabric_counts.get("fracture_zone", 0))
    detail["abyssal_plain_object_count"] = int(
        ocean_fabric_counts.get("abyssal_plain", 0))
    detail["abyssal_hill_object_count"] = int(
        ocean_fabric_counts.get("abyssal_hill", 0))
    detail["age_isochron_object_count"] = int(
        ocean_fabric_counts.get("age_isochron", 0))
    detail["ridge_transform_fracture_core_present"] = bool(
        detail["spreading_center_object_count"] > 0
        and detail["transform_fault_object_count"] > 0
        and detail["fracture_zone_object_count"] > 0
    )
    detail["fracture_zone_combined_span_deg"] = combined_span(
        ocean_fabric, {"fracture_zone"})
    detail["transform_fault_combined_span_deg"] = combined_span(
        ocean_fabric, {"transform_fault"})

    arc_plume = list(world.objects.get("terrain.arc_plume_landforms", []))
    arc_counts = kind_counts(arc_plume)
    detail["arc_plume_landform_object_count"] = int(len(arc_plume))
    detail["arc_plume_landform_kind_counts"] = arc_counts
    detail["hotspot_track_object_count"] = int(arc_counts.get("hotspot_track", 0))
    detail["seamount_chain_object_count"] = int(
        arc_counts.get("seamount_chain", 0) + arc_counts.get("hotspot_track", 0))
    detail["oceanic_plateau_object_count"] = int(
        arc_counts.get("oceanic_plateau", 0)
        + sum(
            1 for obj in arc_plume
            if str(obj.get("kind", "")) == "large_igneous_province"
            and float(obj.get("ocean_fraction", 0.0) or 0.0) >= 0.25
        )
    )
    detail["back_arc_basin_object_count"] = int(arc_counts.get("back_arc_basin", 0))
    detail["island_arc_object_count"] = int(arc_counts.get("island_arc", 0))
    detail["accreted_terrane_object_count"] = int(arc_counts.get("accreted_terrane", 0))
    detail["microcontinent_object_count"] = int(arc_counts.get("microcontinent", 0))

    margin_landforms = list(world.objects.get("terrain.margin_landforms", []))
    margin_counts = kind_counts(margin_landforms)
    detail["margin_landform_kind_counts"] = margin_counts
    detail["passive_margin_wedge_object_count"] = int(
        margin_counts.get("passive_margin_wedge", 0))
    detail["delta_fan_object_count"] = int(margin_counts.get("delta_fan", 0))
    detail["margin_trench_object_count"] = int(margin_counts.get("trench", 0))

    p106_categories = {
        "spreading_center": detail["spreading_center_object_count"] > 0,
        "transform_fault": detail["transform_fault_object_count"] > 0,
        "fracture_zone": detail["fracture_zone_object_count"] > 0,
        "abyssal_plain": detail["abyssal_plain_object_count"] > 0,
        "abyssal_hill": detail["abyssal_hill_object_count"] > 0,
        "age_isochron": detail["age_isochron_object_count"] > 0,
        "seamount_chain": detail["seamount_chain_object_count"] > 0,
        "oceanic_plateau": detail["oceanic_plateau_object_count"] > 0,
        "back_arc_basin": detail["back_arc_basin_object_count"] > 0,
        "passive_margin_wedge": detail["passive_margin_wedge_object_count"] > 0,
        "delta_fan": detail["delta_fan_object_count"] > 0,
        "trench_arc": (
            detail["margin_trench_object_count"] > 0
            or detail["island_arc_object_count"] > 0
        ),
        "microcontinent_or_accreted_terrane": (
            detail["accreted_terrane_object_count"] > 0
            or detail["microcontinent_object_count"] > 0
        ),
    }
    detail["p106_object_category_count"] = int(sum(p106_categories.values()))
    detail["p106_core_ocean_fabric_complete"] = bool(
        detail["ridge_transform_fracture_core_present"]
        and detail["abyssal_plain_object_count"] > 0
        and detail["age_isochron_object_count"] > 0
    )

    if (
        ocean.any()
        and depth_province is not None
        and shelf_width is not None
        and np.asarray(depth_province).shape == (n,)
        and np.asarray(shelf_width).shape == (n,)
    ):
        prov = np.asarray(depth_province, dtype=int)
        width = np.asarray(shelf_width, dtype=np.float64)
        depth = np.maximum(-rel, 0.0)
        far = ocean & ((width <= 0.0) | (width >= 4.0)) & (prov != 6)
        shallow = far & (depth < 1500.0)
        detail["far_ocean_shallow_component_count_lt1500m"] = int(
            len(_component_cell_sets(grid, shallow)))

        crust_type = np.asarray(world.get_field("crust.type", 0.0), dtype=np.float64)
        if crust_type.shape != (n,):
            crust_type = np.zeros(n, dtype=np.float64)
        parented_shallow = (
            (prov == 5)
            | (crust_type >= 0.5)
            | _object_cells_mask(
                n, ocean_fabric,
                {"spreading_center", "age_isochron", "fracture_zone", "transform_fault"},
            )
            | _object_cells_mask(
                n, arc_plume,
                {"hotspot_track", "seamount_chain", "large_igneous_province",
                 "oceanic_plateau", "back_arc_basin", "accreted_terrane"},
            )
        )
        detail["unparented_far_ocean_shallow_component_count_lt1500m"] = int(
            len(_component_cell_sets(grid, shallow & ~parented_shallow)))

        ocean_area = max(float(area[ocean].sum()), 1.0e-12)
        deep_nontrench = far & (depth > 4100.0) & (prov != 6)
        parented_deep = (
            _object_cells_mask(
                n,
                ocean_fabric,
                {"fracture_zone", "transform_fault"},
            )
            | _object_cells_mask(n, margin_landforms, {"trench"})
        )
        isolated_deep = 0
        for comp in _component_cell_sets(grid, deep_nontrench):
            comp_area = float(area[comp].sum())
            parent_fraction = float(np.mean(parented_deep[comp])) if comp.size else 0.0
            if parent_fraction >= 0.25:
                continue
            if comp.size <= 10 or comp_area / ocean_area <= 0.0015:
                isolated_deep += 1
        detail["isolated_far_ocean_deep_pit_component_count_gt4100m"] = int(isolated_deep)

    if land.any():
        crust_type = np.asarray(world.get_field("crust.type", 0.0), dtype=np.float64)
        domain = np.asarray(world.get_field("crust.domain", 0.0), dtype=np.float64)
        origin = np.asarray(world.get_field("crust.origin", 0.0), dtype=np.float64)
        terrane_id = np.asarray(world.get_field("tectonics.terrane_id", -1.0),
                                dtype=np.float64)
        volcanism_age = np.asarray(
            world.get_field("tectonics.volcanism_age_myr", -1.0),
            dtype=np.float64,
        )
        if crust_type.shape != (n,):
            crust_type = np.zeros(n, dtype=np.float64)
        if domain.shape != (n,):
            domain = np.zeros(n, dtype=np.float64)
        if origin.shape != (n,):
            origin = np.zeros(n, dtype=np.float64)
        if terrane_id.shape != (n,):
            terrane_id = np.full(n, -1.0, dtype=np.float64)
        if volcanism_age.shape != (n,):
            volcanism_age = np.full(n, -1.0, dtype=np.float64)
        recent_parented_volcanism = (
            volcanism_age >= max(float(world.time_myr) - 250.0, 0.0)
        )
        parent_mask = (
            _object_cells_mask(
                n,
                arc_plume,
                {"hotspot_track", "seamount_chain", "oceanic_plateau",
                 "large_igneous_province", "island_arc", "accreted_terrane",
                 "microcontinent", "back_arc_basin"},
            )
            | (domain == float(DOMAIN_LIP))
            | (domain == float(DOMAIN_ACCRETED_TERRANE))
            | (origin == float(ORIGIN_PLUME_IMPACT))
            | (origin == float(ORIGIN_ARC))
            | (terrane_id >= 0.0)
            | recent_parented_volcanism
        )
        comps = _component_cell_sets(grid, land)
        comps.sort(key=lambda comp: float(area[comp].sum()), reverse=True)
        parented_count = 0
        parented_area = 0.0
        unparented_oceanic_count = 0
        total_area = max(float(area.sum()), 1.0e-12)
        for comp in comps[1:]:
            comp_area = float(area[comp].sum())
            if comp_area / total_area > 0.018:
                continue
            parent_fraction = float(np.mean(parent_mask[comp])) if comp.size else 0.0
            oceanic_fraction = float(np.mean(crust_type[comp] < 0.5)) if comp.size else 0.0
            if parent_fraction >= 0.25:
                parented_count += 1
                parented_area += comp_area
            elif oceanic_fraction >= 0.50:
                unparented_oceanic_count += 1
        detail["parented_archipelago_component_count"] = int(parented_count)
        detail["parented_archipelago_area_fraction_world"] = float(parented_area / total_area)
        detail["unparented_oceanic_island_component_count"] = int(unparented_oceanic_count)

    detail["p106_full_ocean_feature_coverage"] = bool(
        detail["p106_core_ocean_fabric_complete"]
        and detail["seamount_chain_object_count"] > 0
        and detail["oceanic_plateau_object_count"] > 0
        and detail["back_arc_basin_object_count"] > 0
        and detail["passive_margin_wedge_object_count"] > 0
        and detail["parented_archipelago_component_count"] > 0
        and detail["unparented_far_ocean_shallow_component_count_lt1500m"] == 0
        and detail["isolated_far_ocean_deep_pit_component_count_gt4100m"] == 0
        and detail["unparented_oceanic_island_component_count"] == 0
    )

    return detail


def _label_mismatch_fraction_on_edges(labels: np.ndarray, mask: np.ndarray,
                                      edges: np.ndarray) -> float:
    if edges.size == 0:
        return 0.0
    i, j = edges[:, 0], edges[:, 1]
    valid = mask[i] & mask[j]
    if not valid.any():
        return 0.0
    return float((labels[i[valid]].astype(int) != labels[j[valid]].astype(int)).mean())


def _seam_continuity_metrics(world) -> dict[str, Any]:
    """Audit whether spherical objects are broken by the map seam/projection."""
    grid = world.grid
    n = grid.n
    lon = grid.lon
    edges = grid.edges
    seam_edges = edges[np.abs(lon[edges[:, 0]] - lon[edges[:, 1]]) > 180.0]
    rel = world.get_field("terrain.elevation_m", 0.0) - world.sea_level
    land = rel >= 0.0
    ocean = ~land
    plate = world.get_field("tectonics.plate_id", -1.0).astype(int)
    crust = world.get_field("crust.type", 0.0).astype(int)
    continent_id = world.fields.get("climate.continent_id")
    if continent_id is None:
        continent_id = world.fields.get("tectonics.continent_id")
    basin_id = world.fields.get("ocean.basin_id")

    detail: dict[str, Any] = {
        "seam_edge_count": int(seam_edges.shape[0]),
        "seam_land_ocean_mismatch_fraction": 0.0,
        "seam_plate_boundary_fraction": 0.0,
        "seam_exposed_land_component_mismatch_fraction": 0.0,
        "seam_continental_crust_component_mismatch_fraction": 0.0,
        "seam_ocean_basin_mismatch_fraction": 0.0,
        "seam_elevation_jump_p95_m": 0.0,
        "global_elevation_jump_p95_m": 0.0,
        "seam_to_global_elevation_jump_ratio": 0.0,
        "render_duplicate_land_mismatch_fraction": 0.0,
        "render_duplicate_plate_mismatch_fraction": 0.0,
        "render_duplicate_basin_mismatch_fraction": 0.0,
        "render_duplicate_elevation_delta_p95_m": 0.0,
        "edge_band_land_mismatch_fraction": 0.0,
        "edge_band_elevation_delta_p95_m": 0.0,
        "north_polar_land_fraction": 0.0,
        "south_polar_land_fraction": 0.0,
        "nonfinite_seam_values": 0,
    }

    if edges.size:
        i_all, j_all = edges[:, 0], edges[:, 1]
        global_jump = np.abs(rel[i_all] - rel[j_all])
        detail["global_elevation_jump_p95_m"] = _pct(global_jump, 95)
    if seam_edges.size:
        i, j = seam_edges[:, 0], seam_edges[:, 1]
        detail["seam_land_ocean_mismatch_fraction"] = float((land[i] != land[j]).mean())
        detail["seam_plate_boundary_fraction"] = float((plate[i] != plate[j]).mean())
        seam_jump = np.abs(rel[i] - rel[j])
        detail["seam_elevation_jump_p95_m"] = _pct(seam_jump, 95)
        detail["seam_to_global_elevation_jump_ratio"] = (
            detail["seam_elevation_jump_p95_m"]
            / max(detail["global_elevation_jump_p95_m"], 1e-9)
        )
        if continent_id is not None:
            cid = np.asarray(continent_id, dtype=int)
            detail["seam_exposed_land_component_mismatch_fraction"] = (
                _label_mismatch_fraction_on_edges(cid, land, seam_edges))
        if "tectonics.continent_id" in world.fields:
            tect_cid = world.fields["tectonics.continent_id"].astype(int)
            cont_crust = crust == 1
            detail["seam_continental_crust_component_mismatch_fraction"] = (
                _label_mismatch_fraction_on_edges(tect_cid, cont_crust, seam_edges))
        if basin_id is not None:
            bid = np.asarray(basin_id, dtype=int)
            detail["seam_ocean_basin_mismatch_fraction"] = (
                _label_mismatch_fraction_on_edges(bid, ocean, seam_edges))

    lats = np.linspace(-89.5, 89.5, 180)
    left_exact = grid.nearest_latlon(lats, np.full_like(lats, -180.0))
    right_exact = grid.nearest_latlon(lats, np.full_like(lats, 180.0))
    detail["render_duplicate_land_mismatch_fraction"] = float(
        (land[left_exact] != land[right_exact]).mean())
    detail["render_duplicate_plate_mismatch_fraction"] = float(
        (plate[left_exact] != plate[right_exact]).mean())
    if basin_id is not None:
        bid = np.asarray(basin_id, dtype=int)
        exact_ocean = ocean[left_exact] & ocean[right_exact]
        if exact_ocean.any():
            detail["render_duplicate_basin_mismatch_fraction"] = float(
                (bid[left_exact[exact_ocean]] != bid[right_exact[exact_ocean]]).mean())
    detail["render_duplicate_elevation_delta_p95_m"] = _pct(
        np.abs(rel[left_exact] - rel[right_exact]), 95)

    half_step = 180.0 / 360.0
    left_band = grid.nearest_latlon(lats, np.full_like(lats, -180.0 + half_step))
    right_band = grid.nearest_latlon(lats, np.full_like(lats, 180.0 - half_step))
    detail["edge_band_land_mismatch_fraction"] = float(
        (land[left_band] != land[right_band]).mean())
    detail["edge_band_elevation_delta_p95_m"] = _pct(
        np.abs(rel[left_band] - rel[right_band]), 95)

    north = grid.lat >= 80.0
    south = grid.lat <= -80.0
    if north.any():
        detail["north_polar_land_fraction"] = _area_frac(grid.cell_area, land & north, north)
    if south.any():
        detail["south_polar_land_fraction"] = _area_frac(grid.cell_area, land & south, south)
    for arr in (rel, plate.astype(float)):
        detail["nonfinite_seam_values"] += int((~np.isfinite(arr)).sum())
    return detail


def climate_diagnostics(engine) -> dict[str, Any]:
    """Metrics for improving the climate layer without changing it.

    These diagnostics deliberately separate structural hard failures from model
    quality warnings.  C0 is meant to make thermal walls, missing seasons,
    weak ocean heat transport, and orographic precipitation striping measurable
    before the climate model itself is redesigned.
    """
    w = engine.world
    detail = {
        "context": {
            "spec_name": w.spec.name,
            "time_myr": w.time_myr,
            "land_fraction": w.land_fraction(),
            "has_ocean_current_heat_transport": "ocean.current_heat_transport" in w.fields,
        },
        "temperature": _temperature_metrics(w),
        "precipitation": _precipitation_metrics(w),
        "coasts": _coastal_temperature_metrics(w),
        "seasonality": _seasonality_metrics(w),
        "circulation": _circulation_metrics(w),
        "ocean_currents": _ocean_current_metrics(w),
        "cryosphere": _cryosphere_feedback_metrics(w),
        "geography": _geography_primitive_metrics(w),
    }
    detail["warnings"] = _climate_warnings(detail)
    return detail


def _climate_warnings(detail: dict[str, Any]) -> list[str]:
    warnings: list[str] = []
    context = detail.get("context", {})
    temp = detail["temperature"]
    precip = detail["precipitation"]
    coasts = detail["coasts"]
    seasonality = detail["seasonality"]
    circulation = detail["circulation"]
    ocean_currents = detail["ocean_currents"]
    cryosphere = detail["cryosphere"]
    geography = detail["geography"]
    earthlike = str(context.get("spec_name", "")).startswith("earthlike")
    waterworld = str(context.get("spec_name", "")).startswith("waterworld")

    if temp["max_adjacent_lat_band_delta_C"] > (18.0 if earthlike else 24.0):
        warnings.append("annual temperature has a steep adjacent latitude-band jump")
    if earthlike and temp["mean_temp_C"] < 8.0:
        warnings.append("Earthlike annual mean temperature is cold-biased")
    if earthlike and temp["mean_temp_C"] > 30.0:
        warnings.append("Earthlike annual mean temperature is hot-biased")
    if temp["neighbor_temperature_delta_p99_C"] > (20.0 if earthlike else 30.0):
        warnings.append("cell-neighbour temperature jumps suggest local thermal walls")
    if precip["precip_orographic_concentration"] > 2.2:
        warnings.append("land precipitation is concentrated in high-relief cells")
    if precip["precip_relief_correlation"] > 0.55:
        warnings.append("precipitation is strongly tied to relief; orographic striping likely")
    if (earthlike and coasts["coastal_temperature_bands_compared"] >= 2
            and coasts["coastal_temperature_asymmetry_index"] < 0.6):
        warnings.append("same-latitude coastal temperature asymmetry is weak")
    if not seasonality["has_seasonal_temperature"]:
        warnings.append("seasonal temperature field is absent")
    if not seasonality["has_seasonal_precipitation"]:
        warnings.append("seasonal precipitation field is absent")
    if not circulation["has_seasonal_wind"]:
        warnings.append("seasonal wind field is absent")
    if not circulation["has_land_sea_pressure_proxy"]:
        warnings.append("land-sea pressure proxy is absent")
    if not circulation["has_thermal_wind_anomaly"]:
        warnings.append("thermal wind anomaly field is absent")
    if not circulation["has_geographic_circulation_index"]:
        warnings.append("geographic circulation index is absent")
    if not circulation["has_seasonal_pressure_proxy"]:
        warnings.append("seasonal pressure proxy is absent")
    if not circulation["has_moisture_access"]:
        warnings.append("moisture access field is absent")
    if not circulation["has_monsoon_potential"]:
        warnings.append("monsoon potential field is absent")
    if not circulation["has_source_ocean_warmth"]:
        warnings.append("source-ocean warmth field is absent")
    if waterworld and circulation["has_monsoon_potential"] and (
            circulation["monsoon_potential_global_p95"] > 0.18):
        warnings.append("waterworld has strong fake monsoon potential")
    if (earthlike and seasonality["has_seasonal_precipitation"]
            and seasonality["precip_seasonality_p75"] < 1.05):
        warnings.append("seasonal precipitation contrast is weak")
    if (earthlike and circulation["has_geographic_circulation_index"]
            and circulation["geographic_circulation_index_p90"] < 0.10):
        warnings.append("geography-driven circulation anomalies are weak")
    if (earthlike and circulation["has_itcz_latitude"]
            and (circulation["itcz_DJF_lat_deg"] >= 0.0
                 or circulation["itcz_JJA_lat_deg"] <= 0.0)):
        warnings.append("ITCZ does not migrate between solstitial hemispheres")
    if (earthlike and circulation["has_storm_track_intensity"]
            and (circulation["NH_winter_storm_track_ratio"] < 1.05
                 or circulation["SH_winter_storm_track_ratio"] < 1.05)):
        warnings.append("winter storm tracks are not stronger than summer storm tracks")
    if not context["has_ocean_current_heat_transport"]:
        warnings.append("ocean heat transport diagnostic is absent")
    elif earthlike and ocean_currents["ocean_heat_transport_abs_p95_C"] < 0.35:
        warnings.append("ocean heat transport is weak for an Earthlike ocean")
    if not cryosphere["has_seasonal_sea_ice"]:
        warnings.append("seasonal sea-ice field is absent")
    if not cryosphere["has_snow_persistence"]:
        warnings.append("snow persistence field is absent")
    if not cryosphere["has_cloud_albedo_proxy"]:
        warnings.append("cloud albedo proxy is absent")
    if not cryosphere["has_vegetation_climate_feedback"]:
        warnings.append("vegetation climate feedback proxy is absent")
    if (earthlike and cryosphere["has_seasonal_sea_ice"]
            and cryosphere["sea_ice_adjacent_lat_band_jump_max"] > 0.75):
        warnings.append("seasonal sea-ice edge has a hard adjacent latitude-band jump")
    if ocean_currents["solved_final_ocean_mask_mismatch_fraction"] > 0.03:
        warnings.append("climate ocean mask is stale relative to final terrain")
    if not geography["has_continent_id"] or not geography["has_basin_id"]:
        warnings.append("shared geography primitives are absent")
    if (earthlike and geography["has_shelf_index"]
            and geography["shelf_contrast_near_minus_deep"] < 0.20):
        warnings.append("shelf index does not clearly distinguish shelves from deep ocean")
    if (earthlike and geography["has_barrier_index"]
            and geography["barrier_contrast_high_minus_low"] < 0.15):
        warnings.append("terrain barrier index does not clearly distinguish high relief")
    return warnings


def _seafloor_age_geometry_metrics(world) -> dict[str, Any]:
    grid = world.grid
    area = grid.cell_area
    ctype = world.get_field("crust.type", 0.0).astype(int)
    age = world.get_field("crust.age_myr", 0.0)
    ocean = ctype == 0
    boundaries = world.networks.get("tectonics.boundaries", {})
    ridge = np.zeros(grid.n, dtype=bool)
    trench = np.zeros(grid.n, dtype=bool)
    ridge[np.asarray(boundaries.get("ridge", []), dtype=int)] = True
    trench[np.asarray(boundaries.get("trench", []), dtype=int)] = True
    ridge &= ocean
    trench &= ocean | _dilate_mask(grid, ocean, np.ones(grid.n, dtype=bool), passes=1)

    detail: dict[str, Any] = {
        "ridge_ocean_cells": int(ridge.sum()),
        "trench_cells": int(trench.sum()),
        "ridge_ocean_age_p75_myr": _pct(age[ridge], 75),
        "ocean_with_ridge_path_fraction": 0.0,
        "age_distance_correlation": 0.0,
        "old_ocean_near_trench_fraction_gt180_myr": 0.0,
    }
    if not ocean.any() or not ridge.any():
        return detail

    dist = _graph_distance_from_sources(grid, ridge, ocean)
    reached = ocean & np.isfinite(dist)
    detail["ocean_with_ridge_path_fraction"] = _area_frac(area, reached, ocean)
    if int(reached.sum()) >= 8:
        d = dist[reached]
        a = age[reached]
        if float(np.std(d)) > 0.0 and float(np.std(a)) > 0.0:
            detail["age_distance_correlation"] = float(np.corrcoef(d, a)[0, 1])
    trench_zone = _dilate_mask(grid, trench, ocean, passes=3)
    old = ocean & (age > 180.0)
    detail["old_ocean_near_trench_fraction_gt180_myr"] = _area_frac(area, old & trench_zone, old)
    return detail


def _frame_continuity_metrics(engine) -> dict[str, Any]:
    frames = [
        f for f in engine.archive.frames
        if "tectonics.plate_id" in f.fields or "crust.type" in f.fields
    ]
    if len(frames) < 2:
        return {"n_frame_pairs": 0, "skipped": True}

    area = engine.world.grid.cell_area
    grid = engine.world.grid
    total = float(area.sum())
    frame_plate_components: list[int] = []
    frame_max_components_per_plate: list[int] = []
    frame_min_largest_share: list[float] = []
    late_frame_max_components_per_plate: list[int] = []
    late_frame_min_largest_share: list[float] = []
    plate_persistence: list[float] = []
    fixed_grid_plate_persistence: list[float] = []
    active_plate_jaccard: list[float] = []
    plate_area_delta: list[float] = []
    late_active_plate_jaccard: list[float] = []
    late_plate_area_delta: list[float] = []
    crust_change: list[float] = []
    cont_delta: list[float] = []
    late_cont_delta: list[float] = []
    plate_crust_delta: list[float] = []
    late_plate_crust_delta: list[float] = []
    plate_crust_delta_rate: list[float] = []
    late_plate_crust_delta_rate: list[float] = []
    cont_jaccard: list[float] = []
    continent_id_jaccard: list[float] = []
    continent_id_area_delta: list[float] = []
    crust_domain_change: list[float] = []
    late_crust_domain_change: list[float] = []
    terrain_province_change: list[float] = []
    late_terrain_province_change: list[float] = []
    wilson_phase_change: list[float] = []
    late_wilson_phase_change: list[float] = []
    exposed_land_change: list[float] = []
    late_exposed_land_change: list[float] = []
    late_cutoff = max(1000.0, 0.25 * engine.world.spec.t_end_myr)

    def plate_overlap_persistence(pa: np.ndarray, pb: np.ndarray) -> float:
        """Area-weighted plate continuity under dominant overlap remapping.

        Exact fixed-grid label equality is not a good lineage metric for a
        moving-lid planet over ~Gyr archive intervals: a coherent plate can move
        away from its old geographic cells.  This score still uses the same
        fixed grid, but it allows the dominant area overlap between adjacent
        frame labels to identify split/merge lineage continuity.
        """
        pa = np.asarray(pa, dtype=np.int64)
        pb = np.asarray(pb, dtype=np.int64)
        if pa.shape != pb.shape or pa.size == 0:
            return 0.0
        ids_a, inv_a = np.unique(pa, return_inverse=True)
        ids_b, inv_b = np.unique(pb, return_inverse=True)
        if ids_a.size == 0 or ids_b.size == 0:
            return 1.0
        overlap = np.zeros((ids_a.size, ids_b.size), dtype=np.float64)
        np.add.at(overlap, (inv_a, inv_b), area)
        # Forward preserves split children; backward preserves merges.
        forward = float(overlap.max(axis=0).sum() / total)
        backward = float(overlap.max(axis=1).sum() / total)
        return float(np.clip(0.5 * (forward + backward), 0.0, 1.0))

    for fr in frames:
        if "tectonics.plate_id" not in fr.fields:
            continue
        pf = _plate_fragmentation_from_array(grid, fr.fields["tectonics.plate_id"])
        frame_plate_components.append(pf["total_plate_components"])
        frame_max_components_per_plate.append(pf["max_components_per_plate"])
        frame_min_largest_share.append(pf["min_largest_component_share"])
        if fr.time_myr >= late_cutoff:
            late_frame_max_components_per_plate.append(pf["max_components_per_plate"])
            late_frame_min_largest_share.append(pf["min_largest_component_share"])

    for a, b in zip(frames[:-1], frames[1:]):
        is_late_pair = min(a.time_myr, b.time_myr) >= late_cutoff
        dt_myr = max(float(b.time_myr - a.time_myr), 1.0)
        if "tectonics.plate_id" in a.fields and "tectonics.plate_id" in b.fields:
            pa = a.fields["tectonics.plate_id"].astype(int)
            pb = b.fields["tectonics.plate_id"].astype(int)
            exact_persistence = float(area[pa == pb].sum() / total)
            fixed_grid_plate_persistence.append(exact_persistence)
            plate_persistence.append(plate_overlap_persistence(pa, pb))
            ids_a = set(int(x) for x in np.unique(pa))
            ids_b = set(int(x) for x in np.unique(pb))
            union = ids_a | ids_b
            active_plate_jaccard.append(len(ids_a & ids_b) / len(union) if union else 1.0)
            fracs_a = {pid: float(area[pa == pid].sum() / total) for pid in union}
            fracs_b = {pid: float(area[pb == pid].sum() / total) for pid in union}
            delta = 0.5 * sum(abs(fracs_b[pid] - fracs_a[pid]) for pid in union)
            plate_area_delta.append(delta)
            if is_late_pair:
                late_active_plate_jaccard.append(active_plate_jaccard[-1])
                late_plate_area_delta.append(delta)
        if "crust.type" in a.fields and "crust.type" in b.fields:
            ca = a.fields["crust.type"].astype(int)
            cb = b.fields["crust.type"].astype(int)
            crust_change.append(float(area[ca != cb].sum() / total))
            aa = ca == 1
            bb = cb == 1
            fa = float(area[aa].sum() / total)
            fb = float(area[bb].sum() / total)
            cont_delta.append(abs(fb - fa))
            if is_late_pair:
                late_cont_delta.append(cont_delta[-1])
            union = aa | bb
            cont_jaccard.append(float(area[aa & bb].sum() / area[union].sum())
                                if union.any() else 1.0)
            if "tectonics.plate_id" in a.fields and "tectonics.plate_id" in b.fields:
                pa = a.fields["tectonics.plate_id"].astype(int)
                pb = b.fields["tectonics.plate_id"].astype(int)
                shared = sorted(set(int(x) for x in np.unique(pa))
                                & set(int(x) for x in np.unique(pb)))
                weighted = 0.0
                weight_sum = 0.0
                for pid in shared:
                    ma = pa == pid
                    mb = pb == pid
                    wa = float(area[ma].sum())
                    wb = float(area[mb].sum())
                    if wa <= 0.0 or wb <= 0.0:
                        continue
                    fa_plate = float(area[ma & aa].sum() / wa)
                    fb_plate = float(area[mb & bb].sum() / wb)
                    w = 0.5 * (wa + wb)
                    weighted += w * abs(fb_plate - fa_plate)
                    weight_sum += w
                delta = weighted / weight_sum if weight_sum > 0.0 else 0.0
                plate_crust_delta.append(delta)
                rate = delta * 100.0 / dt_myr
                plate_crust_delta_rate.append(rate)
                if is_late_pair:
                    late_plate_crust_delta.append(delta)
                    late_plate_crust_delta_rate.append(rate)
        if "tectonics.continent_id" in a.fields and "tectonics.continent_id" in b.fields:
            ia = a.fields["tectonics.continent_id"].astype(int)
            ib = b.fields["tectonics.continent_id"].astype(int)
            active_a = set(int(x) for x in np.unique(ia) if x >= 0)
            active_b = set(int(x) for x in np.unique(ib) if x >= 0)
            union_ids = active_a | active_b
            continent_id_jaccard.append(
                len(active_a & active_b) / len(union_ids) if union_ids else 1.0
            )
            fracs_a = {cid: float(area[ia == cid].sum() / total) for cid in union_ids}
            fracs_b = {cid: float(area[ib == cid].sum() / total) for cid in union_ids}
            continent_id_area_delta.append(
                0.5 * sum(abs(fracs_b[cid] - fracs_a[cid]) for cid in union_ids)
                if union_ids else 0.0
            )
        if "crust.domain" in a.fields and "crust.domain" in b.fields:
            da = a.fields["crust.domain"].astype(int)
            db = b.fields["crust.domain"].astype(int)
            val = float(area[da != db].sum() / total)
            crust_domain_change.append(val)
            if is_late_pair:
                late_crust_domain_change.append(val)
        if "terrain.province" in a.fields and "terrain.province" in b.fields:
            ta = a.fields["terrain.province"].astype(int)
            tb = b.fields["terrain.province"].astype(int)
            val = float(area[ta != tb].sum() / total)
            terrain_province_change.append(val)
            if is_late_pair:
                late_terrain_province_change.append(val)
        if ("archive.wilson_cycle_phase" in a.fields
                and "archive.wilson_cycle_phase" in b.fields):
            wa = a.fields["archive.wilson_cycle_phase"].astype(int)
            wb = b.fields["archive.wilson_cycle_phase"].astype(int)
            active = (wa > 0) | (wb > 0)
            val = float(area[(wa != wb) & active].sum() / total)
            wilson_phase_change.append(val)
            if is_late_pair:
                late_wilson_phase_change.append(val)
        if "terrain.elevation_m" in a.fields and "terrain.elevation_m" in b.fields:
            sea_a = float(a.globals.get("ocean.sea_level_m", 0.0))
            sea_b = float(b.globals.get("ocean.sea_level_m", 0.0))
            la = a.fields["terrain.elevation_m"] >= sea_a
            lb = b.fields["terrain.elevation_m"] >= sea_b
            val = float(area[la != lb].sum() / total)
            exposed_land_change.append(val)
            if is_late_pair:
                late_exposed_land_change.append(val)

    return {
        "n_frame_pairs": len(frames) - 1,
        "archive_total_plate_components_max": max(frame_plate_components, default=0),
        "archive_max_components_per_plate": max(frame_max_components_per_plate, default=0),
        "archive_min_largest_component_share": min(frame_min_largest_share, default=1.0),
        "late_archive_max_components_per_plate": max(late_frame_max_components_per_plate,
                                                     default=0),
        "late_archive_min_largest_component_share": min(late_frame_min_largest_share,
                                                        default=1.0),
        "plate_label_persistence_min": min(plate_persistence, default=0.0),
        "plate_label_persistence_mean": float(np.mean(plate_persistence))
        if plate_persistence else 0.0,
        "fixed_grid_plate_label_persistence_min": min(
            fixed_grid_plate_persistence, default=0.0),
        "fixed_grid_plate_label_persistence_mean": float(
            np.mean(fixed_grid_plate_persistence))
        if fixed_grid_plate_persistence else 0.0,
        "active_plate_id_jaccard_min": min(active_plate_jaccard, default=1.0),
        "active_plate_id_jaccard_mean": float(np.mean(active_plate_jaccard))
        if active_plate_jaccard else 1.0,
        "max_plate_area_distribution_delta": max(plate_area_delta, default=0.0),
        "mean_plate_area_distribution_delta": float(np.mean(plate_area_delta))
        if plate_area_delta else 0.0,
        "late_active_plate_id_jaccard_min": min(late_active_plate_jaccard, default=1.0),
        "late_max_plate_area_distribution_delta": max(late_plate_area_delta, default=0.0),
        "late_cutoff_myr": late_cutoff,
        "fixed_grid_plate_change_min": min(fixed_grid_plate_persistence, default=0.0),
        "max_crust_type_change_fraction": max(crust_change, default=0.0),
        "mean_crust_type_change_fraction": float(np.mean(crust_change)) if crust_change else 0.0,
        "fixed_grid_crust_type_change_max": max(crust_change, default=0.0),
        "max_continental_fraction_delta": max(cont_delta, default=0.0),
        "late_max_continental_fraction_delta": max(late_cont_delta, default=0.0),
        "max_plate_crust_composition_delta": max(plate_crust_delta, default=0.0),
        "late_max_plate_crust_composition_delta": max(late_plate_crust_delta, default=0.0),
        "max_plate_crust_composition_delta_per_100myr": max(plate_crust_delta_rate,
                                                            default=0.0),
        "late_max_plate_crust_composition_delta_per_100myr": max(
            late_plate_crust_delta_rate, default=0.0),
        "mean_continental_jaccard": float(np.mean(cont_jaccard)) if cont_jaccard else 0.0,
        "min_continental_jaccard": min(cont_jaccard, default=0.0),
        "continent_id_jaccard_min": min(continent_id_jaccard, default=1.0),
        "max_continent_id_area_distribution_delta": max(continent_id_area_delta,
                                                         default=0.0),
        "max_crust_domain_change_fraction": max(crust_domain_change, default=0.0),
        "late_max_crust_domain_change_fraction": max(late_crust_domain_change,
                                                      default=0.0),
        "max_terrain_province_change_fraction": max(terrain_province_change,
                                                    default=0.0),
        "late_max_terrain_province_change_fraction": max(late_terrain_province_change,
                                                         default=0.0),
        "max_wilson_phase_change_fraction": max(wilson_phase_change, default=0.0),
        "late_max_wilson_phase_change_fraction": max(late_wilson_phase_change,
                                                     default=0.0),
        "max_exposed_land_change_fraction": max(exposed_land_change, default=0.0),
        "late_max_exposed_land_change_fraction": max(late_exposed_land_change,
                                                     default=0.0),
    }


def tectonic_diagnostics(engine) -> dict[str, Any]:
    """Metrics for improving the tectonics/crust layer without changing it."""
    morphology = compute_world_morphology(engine.world)
    detail = {
        "context": {
            "spec_name": engine.world.spec.name,
            "regime_code": engine.world.g("tectonics.regime_code", 0.0),
            "time_myr": engine.world.time_myr,
        },
        "plate_fragmentation": _plate_fragmentation_metrics(engine.world),
        "crust_distribution": _crust_distribution_metrics(engine.world),
        "boundaries": _boundary_metrics(engine.world),
        "seafloor_age_geometry": _seafloor_age_geometry_metrics(engine.world),
        "ocean_geography": _ocean_geography_metrics(engine.world),
        "seam_continuity": _seam_continuity_metrics(engine.world),
        "frame_continuity": _frame_continuity_metrics(engine),
        "morphology": morphology.metrics,
    }
    detail["warnings"] = _tectonic_warnings(detail)
    return detail


def _tectonic_warnings(detail: dict[str, Any]) -> list[str]:
    warnings: list[str] = []
    context = detail.get("context", {})
    regime_code = float(context.get("regime_code", 0.0))
    plate = detail["plate_fragmentation"]
    crust = detail["crust_distribution"]
    frames = detail["frame_continuity"]
    boundaries = detail["boundaries"]
    seafloor = detail["seafloor_age_geometry"]
    ocean_geo = detail.get("ocean_geography", {})
    seam = detail.get("seam_continuity", {})
    morphology = detail.get("morphology", {})
    morphology_warnings = morphology.get("warnings", [])

    if plate["max_components_per_plate"] > 8 or plate["min_largest_component_share"] < 0.60:
        warnings.append("plate fragmentation is high; final plates may read as patchwork")
    if (frames.get("late_archive_max_components_per_plate", 0) > 4
            or frames.get("late_archive_min_largest_component_share", 1.0) < 0.70):
        warnings.append("archive frames contain fragmented plate identities")
    if (regime_code >= 2.0
            and (crust["oceanic_age_p95_myr"] > 300.0
                 or crust["oceanic_old_fraction_gt300_myr"] > 0.10)):
        warnings.append("oceanic crust age is too old for a mature mobile-lid Earth analogue")
    if crust["continental_age_exceeds_world_cells"] > 0:
        warnings.append("some continental crust ages exceed planet age; age semantics need splitting")
    if (frames.get("late_active_plate_id_jaccard_min", 1.0) < 0.60
            or frames.get("late_max_plate_area_distribution_delta", 0.0) > 0.55):
        warnings.append("archive frames show abrupt plate label turnover")
    if (frames.get("late_max_continental_fraction_delta", 0.0) > 0.12
            or frames.get("late_max_plate_crust_composition_delta_per_100myr", 0.0) > 0.09
            or frames.get("late_max_crust_domain_change_fraction", 0.0) > 0.55):
        warnings.append("archive frames show large crust-type jumps")
    if boundaries["boundary_cell_fraction"] > 0.30:
        warnings.append("plate boundary coverage is very dense; convergent/divergent masks may be noisy")
    if (regime_code >= 2.0
            and seafloor["ridge_ocean_cells"] > 0
            and seafloor["ridge_ocean_age_p75_myr"] > 35.0):
        warnings.append("ridge-adjacent oceanic crust is not consistently young")
    if (regime_code >= 2.0
            and seafloor["ocean_with_ridge_path_fraction"] > 0.40
            and seafloor["age_distance_correlation"] < 0.45):
        warnings.append("oceanic crust age does not organize around ridge distance")
    if (ocean_geo.get("has_depth_province")
            and ocean_geo.get("shelf_to_abyss_depth_delta_m", 0.0) < 1200.0):
        warnings.append("ocean bathymetry does not separate shelves from abyssal plains")
    if ocean_geo.get("nearshore_superdeep_fraction_gt2500m", 0.0) > 0.06:
        warnings.append("ocean bathymetry still has too much superdeep water next to coasts")
    if ocean_geo.get("far_ocean_shallow_fraction_lt1500m", 0.0) > 0.15:
        warnings.append("far ocean contains too much shallow water")
    if (ocean_geo.get("trench_fraction_of_ocean", 0.0) > 0.0
            and ocean_geo.get("trench_near_active_margin_fraction", 1.0) < 0.75):
        warnings.append("trench provinces are weakly tied to active margins")
    if (
        ocean_geo.get("max_closed_ocean_candidate_fraction_world", 0.0) > 0.030
        or ocean_geo.get("closed_ocean_candidate_count_gt2pct_world", 0) > 0
    ):
        warnings.append("large closed-ocean candidates exist; ocean connectivity needs seaway/gateway repair")
    if seam.get("seam_to_global_elevation_jump_ratio", 0.0) > 2.5:
        warnings.append("antimeridian seam has unusually large elevation jumps")
    if seam.get("edge_band_land_mismatch_fraction", 0.0) > 0.35:
        warnings.append("rendered map edge band has high land/ocean mismatch across seam")
    if seam.get("seam_continental_crust_component_mismatch_fraction", 0.0) > 0.0:
        warnings.append("continental crust ids change across some antimeridian neighbours")
    warnings.extend(f"morphology: {w}" for w in morphology_warnings)
    return warnings


def _nonfinite_numbers(obj: Any, path: str = "") -> list[str]:
    if isinstance(obj, dict):
        out: list[str] = []
        for k, v in obj.items():
            out.extend(_nonfinite_numbers(v, f"{path}.{k}" if path else str(k)))
        return out
    if isinstance(obj, (list, tuple)):
        out = []
        for i, v in enumerate(obj):
            out.extend(_nonfinite_numbers(v, f"{path}[{i}]"))
        return out
    if isinstance(obj, (int, float, np.integer, np.floating)):
        return [] if np.isfinite(float(obj)) else [path]
    return []


# ----------------------------------------------------------------------
def check_conservation(engine) -> CheckResult:
    w = engine.world
    init = w.g("bgc.carbon_total_init", 0.0)
    now = w.g("bgc.carbon_atm_mol", 0.0) + w.g("bgc.carbon_buried_mol", 0.0)
    rel = abs(now - init) / max(init, 1e-30)
    return CheckResult("conservation.carbon", rel < 1e-6,
                       {"init_mol": init, "now_mol": now, "rel_error": rel})


def check_topology(engine) -> CheckResult:
    w = engine.world
    grid = w.grid
    detail: dict[str, Any] = {}
    ok = True

    plate = w.get_field("tectonics.plate_id", -1.0)
    detail["all_cells_have_plate"] = bool((plate >= 0).all())
    detail["n_plates"] = int(np.unique(plate).size)
    ok &= detail["all_cells_have_plate"]

    # rivers acyclic: following receiver must reach a fixed point
    rivers = w.networks.get("hydrology.rivers")
    if rivers is not None:
        recv = rivers["receiver"]
        acyclic = True
        sample = range(0, grid.n, max(1, grid.n // 2000))
        for c in sample:
            seen = set()
            x = int(c)
            for _ in range(grid.n):
                if x in seen:
                    acyclic = False
                    break
                seen.add(x)
                nx = int(recv[x])
                if nx == x:
                    break
                x = nx
            if not acyclic:
                break
        detail["rivers_acyclic"] = acyclic
        ok &= acyclic
    return CheckResult("topology", ok, detail)


def check_causality(engine) -> CheckResult:
    w = engine.world
    detail: dict[str, Any] = {}
    ok = True

    deposits = w.objects.get("resources.deposits", [])
    bad = [d for d in deposits if not d.get("genesis_model")
           or d.get("formation_age_myr") is None]
    detail["deposits_total"] = len(deposits)
    detail["deposits_without_genesis"] = len(bad)
    ok &= (len(bad) == 0)

    # every mountain cell must have crust thickening / orogeny cause
    elev = w.get_field("terrain.elevation_m")
    sea = w.sea_level
    thick = w.get_field("crust.thickness_m")
    land = elev > sea
    if land.any():
        le = elev - sea
        mountains = land & (le > np.percentile(le[land], 95))
        orog = w.get_field("tectonics.orogeny_age_myr", -1.0)
        volc = w.get_field("tectonics.volcanism_age_myr", -1.0)
        # A high cell is "explained" if it traces to collision orogeny, arc
        # volcanism, or simply thick buoyant crust (isostasy) -- valid for both
        # mobile-lid mountains and stagnant-lid volcanic highlands.
        explained = ((orog[mountains] >= 0) | (volc[mountains] >= 0)
                     | (thick[mountains] > 33000))
        detail["mountain_cells"] = int(mountains.sum())
        detail["mountains_with_cause"] = int(explained.sum())
        ok &= bool(explained.all()) if mountains.any() else True
    return CheckResult("causality", ok, detail)


def compiler_consistency_metrics(engine, cm=None) -> dict[str, Any]:
    if cm is None:
        from aevum.compiler.map_compiler import MapCompiler
        cm = MapCompiler(engine.world, engine.archive).compile(width=64, height=32, n_starts=4)
    w = engine.world
    sea = w.sea_level
    terrain = cm.terrain
    source_land = cm.source_land_fraction
    source_shelf = cm.source_shelf_fraction
    source_province = cm.source_terrain_province
    ice_over_land = (terrain == 5) & (
        np.asarray(source_land, dtype=np.float64) >= 0.45
        if source_land is not None else True
    )
    ice_over_ocean = (terrain == 5) & ~ice_over_land
    compiled_land = ((terrain >= 2) & (terrain <= 4)) | ice_over_land
    compiled_water = (terrain <= 1) | ice_over_ocean
    compiled_deep = terrain == 0
    compiled_mountain = terrain == 4
    center_rel = w.get_field("terrain.elevation_m", 0.0)[cm.cell_index] - sea
    center_land = center_rel >= 0.0

    detail: dict[str, Any] = {
        "has_source_land_fraction": source_land is not None,
        "has_source_shelf_fraction": source_shelf is not None,
        "has_source_terrain_province": source_province is not None,
        "compiled_land_fraction": float(compiled_land.mean()),
        "compiled_land_or_coast_fraction": float((terrain >= 1).mean()),
        "compiled_coast_fraction": float((terrain == 1).mean()),
        "center_land_fraction": float(center_land.mean()),
        "source_land_fraction_mean": 0.0,
        "land_fraction_abs_delta_from_source": 0.0,
        "center_land_to_water_fraction": 0.0,
        "center_ocean_to_land_fraction": 0.0,
        "broad_land_to_water_fraction": 0.0,
        "broad_ocean_to_land_fraction": 0.0,
        "shelf_as_deep_ocean_fraction": 0.0,
        "lowland_as_mountain_fraction": 0.0,
        "terrain_elevation_sign_mismatch_fraction": 0.0,
        "invalid_source_array_shape": False,
        "passed_envelope": True,
    }

    expected = terrain.shape
    for arr in (source_land, source_shelf, source_province):
        if arr is None or np.asarray(arr).shape != expected:
            detail["invalid_source_array_shape"] = True
    if detail["invalid_source_array_shape"]:
        detail["passed_envelope"] = False
        return detail

    source_land = np.asarray(source_land, dtype=np.float64)
    source_shelf = np.asarray(source_shelf, dtype=np.float64)
    source_province = np.asarray(source_province, dtype=int)
    detail["source_land_fraction_mean"] = float(source_land.mean())
    detail["land_fraction_abs_delta_from_source"] = abs(
        detail["compiled_land_fraction"] - detail["source_land_fraction_mean"])

    center_land_count = max(int(center_land.sum()), 1)
    center_ocean_count = max(int((~center_land).sum()), 1)
    detail["center_land_to_water_fraction"] = float(
        (center_land & compiled_water).sum() / center_land_count)
    detail["center_ocean_to_land_fraction"] = float(
        ((~center_land) & compiled_land).sum() / center_ocean_count)

    broad_land = source_land >= 0.70
    broad_ocean = source_land <= 0.10
    broad_shelf = (source_shelf >= 0.55) & (source_land < 0.45)
    lowland = broad_land & np.isin(source_province, [1, 2, 3]) & (cm.elevation - sea < 700.0)
    detail["broad_land_to_water_fraction"] = (
        float((broad_land & compiled_water).sum() / max(int(broad_land.sum()), 1)))
    detail["broad_ocean_to_land_fraction"] = (
        float((broad_ocean & compiled_land).sum() / max(int(broad_ocean.sum()), 1)))
    detail["shelf_as_deep_ocean_fraction"] = (
        float((broad_shelf & compiled_deep).sum() / max(int(broad_shelf.sum()), 1)))
    detail["lowland_as_mountain_fraction"] = (
        float((lowland & compiled_mountain).sum() / max(int(lowland.sum()), 1)))
    terrain_land_sign = compiled_land & (cm.elevation < sea)
    terrain_water_sign = compiled_water & (cm.elevation >= sea)
    detail["terrain_elevation_sign_mismatch_fraction"] = float(
        (terrain_land_sign | terrain_water_sign).mean())

    detail["passed_envelope"] = bool(
        detail["broad_land_to_water_fraction"] <= 0.08
        and detail["broad_ocean_to_land_fraction"] <= 0.06
        and detail["shelf_as_deep_ocean_fraction"] <= 0.20
        and detail["lowland_as_mountain_fraction"] <= 0.08
        and detail["terrain_elevation_sign_mismatch_fraction"] <= 0.02
    )
    return detail


def check_compiler_consistency(engine) -> CheckResult:
    detail = compiler_consistency_metrics(engine)
    hard_failures: list[str] = []
    if detail["invalid_source_array_shape"]:
        hard_failures.append("compiled map missing source aggregation arrays")
    if not detail["passed_envelope"]:
        hard_failures.append("compiled map contradicts source terrain/ocean provinces")
    detail["hard_failures"] = hard_failures
    return CheckResult("compiler.consistency", len(hard_failures) == 0, detail)


def check_earthlike_plausibility(engine) -> CheckResult:
    """Envelope checks for the Earth-analogue preset.

    These are not universal planetary rules.  They catch regressions in the
    reference world: thermal walls, all-desert land, high oceanic highlands and
    unusable strategy starts.
    """
    w = engine.world
    if not w.spec.name.startswith("earthlike"):
        return CheckResult("plausibility.earthlike", True, {"skipped": True})

    detail: dict[str, Any] = {}
    ok = True
    area = w.grid.cell_area
    elev = w.get_field("terrain.elevation_m")
    rel = elev - w.sea_level
    land = rel >= 0
    crust = w.get_field("crust.type", 0.0).astype(int)
    temp = w.get_field("climate.surface_temperature", 288.0) - 273.15
    precip = w.get_field("climate.precipitation", 0.0)
    biome = w.get_field("biosphere.biome", 0.0).astype(int)

    land_fraction = w.land_fraction()
    high_oceanic = int(((crust == 0) & land & (rel > 1500)).sum())
    lat_bins = np.arange(-90, 91, 10)
    band_means = []
    for lo, hi in zip(lat_bins[:-1], lat_bins[1:]):
        mask = (w.grid.lat >= lo) & (w.grid.lat < hi)
        band_means.append(float(np.average(temp[mask], weights=area[mask])))
    band_delta = max(abs(a - b) for a, b in zip(band_means[:-1], band_means[1:]))
    wet_land_p75 = float(np.percentile(precip[land], 75)) if land.any() else 0.0
    habitable_biomes = int(np.isin(biome[land], [3, 4, 6]).sum()) if land.any() else 0

    detail.update({
        "land_fraction": land_fraction,
        "high_oceanic_land": high_oceanic,
        "mean_temp_C": float(np.average(temp, weights=area)),
        "max_temp_C": float(temp.max()),
        "lat_band_delta_C": band_delta,
        "land_precip_p75": wet_land_p75,
        "habitable_biome_cells": habitable_biomes,
    })
    ok &= 0.20 <= land_fraction <= 0.38
    ok &= high_oceanic == 0
    ok &= 8.0 <= detail["mean_temp_C"] <= 30.0
    ok &= detail["max_temp_C"] < 60.0
    ok &= band_delta < 22.0
    ok &= wet_land_p75 > 300.0
    ok &= habitable_biomes > 0

    try:
        from aevum.compiler.map_compiler import MapCompiler
        cm = MapCompiler(w, engine.archive).compile(width=64, height=32, n_starts=4)
        ice_fraction = float((cm.terrain == 5).mean())
        n_starts = len(cm.starts)
        max_start_hazard = max((cm.hazard[r, c] for r, c in cm.starts), default=99.0)
        detail.update({
            "map_ice_fraction": ice_fraction,
            "starts": n_starts,
            "max_start_hazard": float(max_start_hazard),
            "local_yield_cv": cm.fairness.get("local_yield_cv"),
        })
        compiler_detail = compiler_consistency_metrics(engine, cm)
        detail["compiler_consistency"] = compiler_detail
        ok &= ice_fraction < 0.15
        ok &= n_starts == 4
        ok &= max_start_hazard <= 1.0
        ok &= cm.fairness.get("local_yield_cv", 1.0) < 0.40
        ok &= compiler_detail["passed_envelope"]
    except Exception as exc:  # pragma: no cover - defensive validation detail
        detail["map_compile_error"] = repr(exc)
        ok = False

    return CheckResult("plausibility.earthlike", bool(ok), detail)


EARTH_GEOMORPHOLOGY_FEATURES: tuple[dict[str, Any], ...] = (
    {"feature": "shield", "suite": "E1", "object_set": "terrain.continental_landforms"},
    {"feature": "platform", "suite": "E1", "object_set": "terrain.continental_landforms"},
    {"feature": "interior_basin", "suite": "E1", "object_set": "terrain.continental_landforms"},
    {"feature": "rift_basin", "suite": "E1", "object_set": "terrain.continental_landforms"},
    {"feature": "orogen", "suite": "E1", "object_set": "terrain.continental_landforms"},
    {"feature": "old_subdued_orogen", "suite": "E1", "object_set": "terrain.continental_landforms"},
    {"feature": "foreland_basin", "suite": "E1", "object_set": "terrain.continental_landforms"},
    {"feature": "plateau", "suite": "E1", "object_set": "terrain.continental_landforms"},
    {"feature": "arc_microcontinent", "suite": "E1", "object_set": "terrain.continental_landforms"},
    {"feature": "passive_margin_wedge", "suite": "E2", "object_set": "terrain.margin_landforms"},
    {"feature": "trench", "suite": "E2", "object_set": "terrain.margin_landforms"},
    {"feature": "forearc_accretionary_prism", "suite": "E2", "object_set": "terrain.margin_landforms"},
    {"feature": "volcanic_arc", "suite": "E2", "object_set": "terrain.margin_landforms"},
    {"feature": "delta_fan", "suite": "E2", "object_set": "terrain.margin_landforms"},
    {"feature": "spreading_center", "suite": "E3", "object_set": "terrain.ocean_fabric"},
    {"feature": "transform_fault", "suite": "E3", "object_set": "terrain.ocean_fabric"},
    {"feature": "fracture_zone", "suite": "E3", "object_set": "terrain.ocean_fabric"},
    {"feature": "abyssal_plain", "suite": "E3", "object_set": "terrain.ocean_fabric"},
    {"feature": "age_isochron", "suite": "E3", "object_set": "terrain.ocean_fabric"},
    {"feature": "island_arc", "suite": "E4", "object_set": "terrain.arc_plume_landforms"},
    {"feature": "back_arc_basin", "suite": "E4", "object_set": "terrain.arc_plume_landforms"},
    {"feature": "accreted_terrane", "suite": "E4", "object_set": "terrain.arc_plume_landforms"},
    {"feature": "hotspot_track", "suite": "E4", "object_set": "terrain.arc_plume_landforms"},
    {"feature": "large_igneous_province", "suite": "E4", "object_set": "terrain.arc_plume_landforms"},
    {"feature": "ice_sheet_loading", "suite": "E5", "object_set": "terrain.cryosphere_landforms"},
    {"feature": "glacial_erosion", "suite": "E5", "object_set": "terrain.cryosphere_landforms"},
    {"feature": "postglacial_rebound", "suite": "E5", "object_set": "terrain.cryosphere_landforms"},
)


def _fixture_suite_passed(fixture_suites: dict[str, Any] | None,
                          suite: str) -> bool:
    if not fixture_suites:
        return False
    summary = fixture_suites.get(suite)
    if not isinstance(summary, dict):
        return False
    acceptance = summary.get("acceptance", {})
    return bool(
        summary.get("status") == "pass"
        and acceptance.get("all_microbenchmarks_pass", summary.get("status") == "pass")
    )


def _object_parented(obj: dict[str, Any]) -> bool:
    for key in (
        "parent_tectonic_object_ids",
        "parent_continent_ids",
        "parent_basin_id",
        "parent_boundary_object_id",
        "parent_margin_ids",
        "parent_rift_system_ids",
        "parent_closing_margin_ids",
        "parent_suture_ids",
    ):
        value = obj.get(key)
        if isinstance(value, (list, tuple, set)) and len(value) > 0:
            return True
        if value not in (None, "", [], {}, ()):
            return True
    return False


def _feature_score(*present: bool) -> float:
    if not present:
        return 0.0
    return float(sum(1 for item in present if item) / len(present))


def earth_geomorphology_coverage_metrics(
    engine,
    cm=None,
    fixture_suites: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Integrated Earth-derived geomorphology coverage.

    Fixture suites prove that a landform class can be generated for the right
    process reason.  The generated-world pass separately records whether that
    class is visible in the current Earth-like reference world.  Ice-related
    features can remain fixture-only while the climate/cryosphere model is
    still intentionally paused.
    """
    del cm
    w = engine.world
    feature_rows: list[dict[str, Any]] = []
    by_feature: dict[str, dict[str, Any]] = {}
    object_counts: dict[str, int] = {}
    object_area: dict[str, float] = {}

    for spec in EARTH_GEOMORPHOLOGY_FEATURES:
        feature = str(spec["feature"])
        suite = str(spec["suite"])
        object_set = str(spec["object_set"])
        objs = [
            obj for obj in w.objects.get(object_set, [])
            if obj.get("kind") == feature
        ]
        count = len(objs)
        area_fraction = float(sum(float(obj.get("area_fraction", 0.0)) for obj in objs))
        fixture_passed = _fixture_suite_passed(fixture_suites, suite)
        if count > 0 and area_fraction >= 5e-4:
            level = "world"
        elif count > 0:
            level = "weak"
        elif fixture_passed:
            level = "partial"
        else:
            level = "none"
        row = {
            "feature": feature,
            "suite": suite,
            "object_set": object_set,
            "level": level,
            "world_object_count": int(count),
            "world_area_fraction": area_fraction,
            "fixture_suite_passed": bool(fixture_passed),
            "ice_related": suite == "E5",
        }
        feature_rows.append(row)
        by_feature[feature] = row
        object_counts[feature] = int(count)
        object_area[feature] = area_fraction

    generated_sets = sorted({str(spec["object_set"]) for spec in EARTH_GEOMORPHOLOGY_FEATURES})
    major_area = 0.0
    parentless_area = 0.0
    parentless_count = 0
    generated_count = 0
    for object_set in generated_sets:
        for obj in w.objects.get(object_set, []):
            generated_count += 1
            area_fraction = float(obj.get("area_fraction", 0.0))
            if area_fraction < 0.002:
                continue
            major_area += area_fraction
            if not _object_parented(obj):
                parentless_count += 1
                parentless_area += area_fraction

    def has(feature: str) -> bool:
        return object_counts.get(feature, 0) > 0

    group_world_counts = {
        "continental_interior": sum(
            1 for row in feature_rows
            if row["suite"] == "E1" and row["world_object_count"] > 0),
        "margin_shelf": sum(
            1 for row in feature_rows
            if row["suite"] == "E2" and row["world_object_count"] > 0),
        "ocean_basin_fabric": sum(
            1 for row in feature_rows
            if row["suite"] == "E3" and row["world_object_count"] > 0),
        "arc_plume": sum(
            1 for row in feature_rows
            if row["suite"] == "E4" and row["world_object_count"] > 0),
        "cryosphere_surface_process": sum(
            1 for row in feature_rows
            if row["suite"] == "E5" and row["world_object_count"] > 0),
    }
    scores = {
        "ridge_transform_continuity_score": _feature_score(
            has("spreading_center"), has("transform_fault"),
            has("fracture_zone"), has("age_isochron")),
        "passive_margin_profile_score": _feature_score(
            has("passive_margin_wedge"), has("delta_fan")),
        "active_margin_profile_score": _feature_score(
            has("trench"), has("forearc_accretionary_prism"), has("volcanic_arc")),
        "ocean_age_isochron_score": _feature_score(
            has("age_isochron"), has("spreading_center")),
        "orogen_lifecycle_score": _feature_score(
            has("orogen"), has("old_subdued_orogen"), has("foreland_basin")),
        "sedimentary_margin_score": _feature_score(
            has("passive_margin_wedge"), has("interior_basin"), has("delta_fan")),
        "hotspot_track_score": _feature_score(
            has("hotspot_track"), has("large_igneous_province")),
        "ice_geomorphology_score": (
            1.0 if group_world_counts["cryosphere_surface_process"] > 0
            else (0.5 if _fixture_suite_passed(fixture_suites, "E5") else 0.0)
        ),
    }

    covered = [row for row in feature_rows if row["level"] == "world"]
    weak = [row for row in feature_rows if row["level"] == "weak"]
    partial = [row for row in feature_rows if row["level"] == "partial"]
    missing = [row for row in feature_rows if row["level"] == "none"]
    non_ice_missing = [row for row in missing if not row["ice_related"]]
    ice_missing = [row for row in missing if row["ice_related"]]
    fixture_only = [row for row in feature_rows if row["level"] == "partial"]
    parentless_major_fraction = (
        parentless_area / major_area if major_area > 0.0 else 0.0
    )
    acceptance = {
        "no_non_ice_feature_none": len(non_ice_missing) == 0,
        "plate_boundary_and_ocean_features_partial": all(
            row["level"] != "none"
            for row in feature_rows
            if row["suite"] in {"E2", "E3", "E4"}),
        "generated_world_has_continental_interior_objects": (
            group_world_counts["continental_interior"] >= 3),
        "generated_world_has_margin_objects": group_world_counts["margin_shelf"] >= 3,
        "generated_world_has_ocean_fabric_objects": (
            group_world_counts["ocean_basin_fabric"] >= 2),
        "generated_world_has_arc_plume_objects": group_world_counts["arc_plume"] >= 3,
        "parentless_major_landform_fraction_ok": parentless_major_fraction <= 0.10,
        "fixture_suites_available": all(
            _fixture_suite_passed(fixture_suites, suite)
            for suite in ("E1", "E2", "E3", "E4", "E5")),
    }

    warnings: list[str] = []
    if parentless_major_fraction > 0.05:
        warnings.append(
            "major generated landforms still need stronger persistent parent-object links")
    if fixture_only:
        names = ", ".join(row["feature"] for row in fixture_only[:10])
        warnings.append(f"fixture-only geomorphology features in this world: {names}")
    if group_world_counts["cryosphere_surface_process"] == 0:
        warnings.append(
            "cryosphere landforms are fixture-backed but absent in the current generated world")
    if scores["ridge_transform_continuity_score"] < 0.50:
        warnings.append("generated ocean fabric lacks transform/fracture-zone expression")
    if scores["orogen_lifecycle_score"] < 0.34:
        warnings.append("generated continent lacks a complete orogen lifecycle expression")

    hard_failures: list[str] = []
    if non_ice_missing:
        hard_failures.append(
            "non-ice geomorphology features have no fixture or world coverage: "
            + ", ".join(row["feature"] for row in non_ice_missing)
        )
    if not acceptance["plate_boundary_and_ocean_features_partial"]:
        hard_failures.append("plate-boundary/ocean-basin features are below partial coverage")
    if not acceptance["parentless_major_landform_fraction_ok"]:
        hard_failures.append("too much major geomorphology area lacks parent objects")
    if not acceptance["generated_world_has_continental_interior_objects"]:
        hard_failures.append("generated world lacks enough continental-interior object expression")
    if not acceptance["generated_world_has_margin_objects"]:
        hard_failures.append("generated world lacks enough margin/shelf object expression")
    if not acceptance["generated_world_has_ocean_fabric_objects"]:
        hard_failures.append("generated world lacks enough ocean-basin fabric object expression")
    if not acceptance["generated_world_has_arc_plume_objects"]:
        hard_failures.append("generated world lacks enough arc/plume object expression")

    return {
        "schema": "aevum.earth_geomorphology_coverage.v1",
        "spec_name": w.spec.name,
        "time_myr": float(w.time_myr),
        "feature_count": len(feature_rows),
        "covered_feature_count": len(covered),
        "partial_feature_count": len(partial),
        "weak_feature_count": len(weak),
        "missing_feature_count": len(missing),
        "non_ice_missing_feature_count": len(non_ice_missing),
        "ice_missing_feature_count": len(ice_missing),
        "generated_landform_object_count": int(generated_count),
        "parentless_major_landform_count": int(parentless_count),
        "parentless_major_landform_fraction": float(parentless_major_fraction),
        "group_world_feature_counts": group_world_counts,
        "scores": scores,
        "acceptance": acceptance,
        "warnings": warnings,
        "hard_failures": hard_failures,
        "features": feature_rows,
    }


def check_earth_geomorphology_coverage(
    engine,
    fixture_suites: dict[str, Any] | None = None,
) -> CheckResult:
    detail = earth_geomorphology_coverage_metrics(
        engine, fixture_suites=fixture_suites)
    return CheckResult(
        "geomorphology.earth_coverage",
        len(detail["hard_failures"]) == 0,
        detail,
    )


def check_tectonic_diagnostics(engine) -> CheckResult:
    """Diagnostic-only quality gate for the tectonics/crust layer.

    The first enhancement step is to make shortcomings measurable before
    changing the generator.  Hard failures here are structural integrity issues;
    model-quality problems are reported as warnings in the detail payload.
    """
    detail = tectonic_diagnostics(engine)
    hard_failures: list[str] = []
    plate = detail["plate_fragmentation"]
    crust = detail["crust_distribution"]
    ocean_geo = detail["ocean_geography"]
    seam = detail["seam_continuity"]

    if plate["n_active_plates"] < 1:
        hard_failures.append("no active plates")
    if crust["negative_age_cells"] > 0:
        hard_failures.append("negative crust ages")
    if crust["negative_thickness_cells"] > 0:
        hard_failures.append("negative crust thickness")
    if crust["non_binary_crust_cells"] > 0:
        hard_failures.append("crust.type contains values outside {0, 1}")
    if crust["invalid_origin_cells"] > 0:
        hard_failures.append("crust.origin contains invalid negative codes")
    if crust["stability_out_of_range_cells"] > 0:
        hard_failures.append("crust.stability outside [0, 1]")
    if crust["continental_ridge_origin_cells"] > 0:
        hard_failures.append("continental crust still marked with ridge-oceanic origin")
    for key, message in [
        ("invalid_basin_id_shape", "ocean.basin_id must have shape (n_cells,)"),
        ("invalid_margin_type_shape", "ocean.margin_type must have shape (n_cells,)"),
        ("invalid_depth_province_shape", "ocean.depth_province must have shape (n_cells,)"),
        ("invalid_gateway_id_shape", "ocean.gateway_id must have shape (n_cells,)"),
        ("invalid_shelf_width_shape", "ocean.shelf_width must have shape (n_cells,)"),
    ]:
        if ocean_geo[key]:
            hard_failures.append(message)
    if ocean_geo["nonfinite_ocean_geography_cells"] > 0:
        hard_failures.append("non-finite ocean geography cells")
    if ocean_geo["land_with_depth_province_cells"] > 0:
        hard_failures.append("land cells labelled with ocean.depth_province")
    if ocean_geo["ocean_without_depth_province_cells"] > 0:
        hard_failures.append("ocean cells missing ocean.depth_province")
    if ocean_geo["land_with_margin_type_cells"] > 0:
        hard_failures.append("land cells labelled with ocean.margin_type")
    if ocean_geo["ocean_without_margin_type_cells"] > 0:
        hard_failures.append("ocean cells missing ocean.margin_type")
    if ocean_geo["land_with_gateway_id_cells"] > 0:
        hard_failures.append("land cells labelled with ocean.gateway_id")
    if ocean_geo["land_with_shelf_width_cells"] > 0:
        hard_failures.append("land cells labelled with ocean.shelf_width")
    if seam["nonfinite_seam_values"] > 0:
        hard_failures.append("non-finite seam continuity values")
    if seam["seam_exposed_land_component_mismatch_fraction"] > 0.0:
        hard_failures.append("exposed land component ids break across antimeridian seam")
    if seam["seam_ocean_basin_mismatch_fraction"] > 0.0:
        hard_failures.append("ocean basin ids break across antimeridian seam")
    nonfinite = _nonfinite_numbers(detail)
    if nonfinite:
        hard_failures.append("non-finite diagnostic values: " + ", ".join(nonfinite[:8]))

    detail["hard_failures"] = hard_failures
    return CheckResult("diagnostics.tectonics", len(hard_failures) == 0, detail)


def check_climate_diagnostics(engine) -> CheckResult:
    """Diagnostic-only quality gate for the climate layer.

    Hard failures are invalid or missing data.  Model-quality issues such as
    thermal walls, absent seasons, weak ocean-current influence, or precipitation
    striping are warnings so C1-C4 can improve them incrementally.
    """
    detail = climate_diagnostics(engine)
    hard_failures: list[str] = []
    temp = detail["temperature"]
    precip = detail["precipitation"]
    seasonality = detail["seasonality"]
    circulation = detail["circulation"]
    ocean_currents = detail["ocean_currents"]
    cryosphere = detail["cryosphere"]
    geography = detail["geography"]

    if not temp["has_surface_temperature"]:
        hard_failures.append("missing climate.surface_temperature")
    if not precip["has_precipitation"]:
        hard_failures.append("missing climate.precipitation")
    if temp["nonfinite_temperature_cells"] > 0:
        hard_failures.append("non-finite climate.surface_temperature cells")
    if precip["nonfinite_precipitation_cells"] > 0:
        hard_failures.append("non-finite climate.precipitation cells")
    if precip["negative_precipitation_cells"] > 0:
        hard_failures.append("negative climate.precipitation cells")
    if seasonality["invalid_seasonal_temperature_shape"]:
        hard_failures.append("climate.seasonal_temperature must have shape (4, n_cells)")
    if seasonality["invalid_seasonal_precipitation_shape"]:
        hard_failures.append("climate.seasonal_precipitation must have shape (4, n_cells)")
    if seasonality["nonfinite_seasonal_precipitation_cells"] > 0:
        hard_failures.append("non-finite climate.seasonal_precipitation cells")
    if seasonality["negative_seasonal_precipitation_cells"] > 0:
        hard_failures.append("negative climate.seasonal_precipitation cells")
    if (seasonality["has_seasonal_precipitation"]
            and not seasonality["annual_precip_matches_seasonal_aggregate"]):
        hard_failures.append("annual precipitation must match seasonal aggregate")
    if not cryosphere["has_seasonal_sea_ice"]:
        hard_failures.append("missing cryosphere.seasonal_sea_ice")
    if not cryosphere["has_seasonal_snow"]:
        hard_failures.append("missing cryosphere.seasonal_snow")
    if not cryosphere["has_snow_persistence"]:
        hard_failures.append("missing cryosphere.snow_persistence")
    if not cryosphere["has_seasonal_cloud_albedo_proxy"]:
        hard_failures.append("missing climate.seasonal_cloud_albedo_proxy")
    if not cryosphere["has_cloud_albedo_proxy"]:
        hard_failures.append("missing climate.cloud_albedo_proxy")
    if not cryosphere["has_vegetation_climate_feedback"]:
        hard_failures.append("missing biosphere.vegetation_climate_feedback")
    for key, message in [
        ("invalid_seasonal_sea_ice_shape",
         "cryosphere.seasonal_sea_ice must have shape (4, n_cells)"),
        ("invalid_sea_ice_shape", "cryosphere.sea_ice must have shape (n_cells,)"),
        ("invalid_seasonal_snow_shape",
         "cryosphere.seasonal_snow must have shape (4, n_cells)"),
        ("invalid_snow_persistence_shape",
         "cryosphere.snow_persistence must have shape (n_cells,)"),
        ("invalid_seasonal_cloud_albedo_shape",
         "climate.seasonal_cloud_albedo_proxy must have shape (4, n_cells)"),
        ("invalid_cloud_albedo_shape",
         "climate.cloud_albedo_proxy must have shape (n_cells,)"),
        ("invalid_vegetation_feedback_shape",
         "biosphere.vegetation_climate_feedback must have shape (n_cells,)"),
    ]:
        if cryosphere[key]:
            hard_failures.append(message)
    for key, message in [
        ("nonfinite_seasonal_sea_ice_cells",
         "non-finite cryosphere.seasonal_sea_ice cells"),
        ("nonfinite_sea_ice_cells", "non-finite cryosphere.sea_ice cells"),
        ("nonfinite_seasonal_snow_cells",
         "non-finite cryosphere.seasonal_snow cells"),
        ("nonfinite_snow_persistence_cells",
         "non-finite cryosphere.snow_persistence cells"),
        ("nonfinite_seasonal_cloud_albedo_cells",
         "non-finite climate.seasonal_cloud_albedo_proxy cells"),
        ("nonfinite_cloud_albedo_cells",
         "non-finite climate.cloud_albedo_proxy cells"),
        ("nonfinite_vegetation_feedback_cells",
         "non-finite biosphere.vegetation_climate_feedback cells"),
        ("out_of_bounds_seasonal_sea_ice_cells",
         "cryosphere.seasonal_sea_ice outside [0, 1]"),
        ("out_of_bounds_sea_ice_cells", "cryosphere.sea_ice outside [0, 1]"),
        ("out_of_bounds_seasonal_snow_cells",
         "cryosphere.seasonal_snow outside [0, 1]"),
        ("out_of_bounds_snow_persistence_cells",
         "cryosphere.snow_persistence outside [0, 1]"),
        ("out_of_bounds_seasonal_cloud_albedo_cells",
         "climate.seasonal_cloud_albedo_proxy outside [0, 1]"),
        ("out_of_bounds_cloud_albedo_cells",
         "climate.cloud_albedo_proxy outside [0, 1]"),
        ("out_of_bounds_vegetation_feedback_cells",
         "biosphere.vegetation_climate_feedback outside [0, 1]"),
    ]:
        if cryosphere[key] > 0:
            hard_failures.append(message)
    if cryosphere["seasonal_sea_ice_land_abs_max"] > 1.0e-8:
        hard_failures.append("cryosphere.seasonal_sea_ice leaks onto land")
    if cryosphere["sea_ice_land_abs_max"] > 1.0e-8:
        hard_failures.append("cryosphere.sea_ice leaks onto land")
    if cryosphere["seasonal_snow_ocean_abs_max"] > 1.0e-8:
        hard_failures.append("cryosphere.seasonal_snow leaks onto ocean")
    if cryosphere["snow_persistence_ocean_abs_max"] > 1.0e-8:
        hard_failures.append("cryosphere.snow_persistence leaks onto ocean")
    if cryosphere["vegetation_feedback_ocean_abs_max"] > 1.0e-8:
        hard_failures.append("biosphere.vegetation_climate_feedback leaks onto ocean")
    if seasonality["invalid_moisture_flow_precipitation_response_shape"]:
        hard_failures.append(
            "climate.moisture_flow_precipitation_response must have shape (4, n_cells)")
    if seasonality["nonfinite_moisture_flow_precipitation_response_cells"] > 0:
        hard_failures.append(
            "non-finite climate.moisture_flow_precipitation_response cells")
    if seasonality["invalid_moisture_source_basin_id_shape"]:
        hard_failures.append(
            "atmosphere.moisture_source_basin_id must have shape (4, n_cells)")
    if seasonality["nonfinite_moisture_source_basin_id_cells"] > 0:
        hard_failures.append("non-finite atmosphere.moisture_source_basin_id cells")
    if seasonality["invalid_moisture_budget_region_id_shape"]:
        hard_failures.append(
            "climate.moisture_budget_region_id must have shape (4, n_cells)")
    if seasonality["nonfinite_moisture_budget_region_id_cells"] > 0:
        hard_failures.append("non-finite climate.moisture_budget_region_id cells")
    if seasonality["invalid_precipitation_response_region_id_shape"]:
        hard_failures.append(
            "climate.precipitation_response_region_id must have shape (4, n_cells)")
    if seasonality["nonfinite_precipitation_response_region_id_cells"] > 0:
        hard_failures.append(
            "non-finite climate.precipitation_response_region_id cells")
    if seasonality["invalid_receiver_catchment_id_shape"]:
        hard_failures.append(
            "climate.receiver_catchment_id must have shape (4, n_cells)")
    if seasonality["nonfinite_receiver_catchment_id_cells"] > 0:
        hard_failures.append("non-finite climate.receiver_catchment_id cells")
    if seasonality["invalid_source_basin_supply_index_shape"]:
        hard_failures.append(
            "climate.source_basin_supply_index must have shape (4, n_cells)")
    if seasonality["nonfinite_source_basin_supply_index_cells"] > 0:
        hard_failures.append("non-finite climate.source_basin_supply_index cells")
    if seasonality["invalid_receiver_catchment_supply_balance_shape"]:
        hard_failures.append(
            "climate.receiver_catchment_supply_balance must have shape (4, n_cells)")
    if seasonality["nonfinite_receiver_catchment_supply_balance_cells"] > 0:
        hard_failures.append(
            "non-finite climate.receiver_catchment_supply_balance cells")
    if seasonality["invalid_receiver_supply_precipitation_feedback_shape"]:
        hard_failures.append(
            "climate.receiver_supply_precipitation_feedback must have shape (4, n_cells)")
    if seasonality["nonfinite_receiver_supply_precipitation_feedback_cells"] > 0:
        hard_failures.append(
            "non-finite climate.receiver_supply_precipitation_feedback cells")
    if circulation["invalid_seasonal_wind_shape"]:
        hard_failures.append("atmosphere.seasonal_wind must have shape (4, n_cells, 3)")
    if circulation["invalid_background_seasonal_wind_shape"]:
        hard_failures.append(
            "atmosphere.background_seasonal_wind must have shape (4, n_cells, 3)")
    if circulation["invalid_land_sea_pressure_shape"]:
        hard_failures.append("atmosphere.land_sea_pressure_proxy must have shape (4, n_cells)")
    if circulation["invalid_seasonal_pressure_shape"]:
        hard_failures.append("atmosphere.seasonal_pressure_proxy must have shape (4, n_cells)")
    if circulation["invalid_moisture_access_shape"]:
        hard_failures.append("atmosphere.moisture_access must have shape (4, n_cells)")
    if circulation["invalid_monsoon_potential_shape"]:
        hard_failures.append("atmosphere.monsoon_potential must have shape (4, n_cells)")
    if circulation["invalid_source_ocean_warmth_shape"]:
        hard_failures.append("atmosphere.source_ocean_warmth must have shape (4, n_cells)")
    if circulation["invalid_terrain_blocking_shape"]:
        hard_failures.append("atmosphere.terrain_blocking must have shape (n_cells,)")
    if circulation["invalid_thermal_wind_anomaly_shape"]:
        hard_failures.append(
            "atmosphere.thermal_wind_anomaly must have shape (4, n_cells, 3)")
    if circulation["invalid_orographic_wind_anomaly_shape"]:
        hard_failures.append(
            "atmosphere.orographic_wind_anomaly must have shape (4, n_cells, 3)")
    if circulation["invalid_geographic_circulation_index_shape"]:
        hard_failures.append(
            "atmosphere.geographic_circulation_index must have shape (n_cells,)")
    if circulation["invalid_itcz_latitude_shape"]:
        hard_failures.append("atmosphere.itcz_latitude must have shape (4,)")
    if circulation["invalid_itcz_intensity_shape"]:
        hard_failures.append("atmosphere.itcz_intensity must have shape (4, n_cells)")
    if circulation["invalid_storm_track_shape"]:
        hard_failures.append("atmosphere.storm_track_intensity must have shape (4, n_cells)")
    if circulation["nonfinite_seasonal_wind_cells"] > 0:
        hard_failures.append("non-finite atmosphere.seasonal_wind cells")
    if circulation["nonfinite_moisture_access_cells"] > 0:
        hard_failures.append("non-finite atmosphere.moisture_access cells")
    if circulation["nonfinite_monsoon_potential_cells"] > 0:
        hard_failures.append("non-finite atmosphere.monsoon_potential cells")
    if circulation["nonfinite_geographic_circulation_cells"] > 0:
        hard_failures.append("non-finite atmosphere.geographic_circulation_index cells")
    if circulation["max_wind_normal_component"] > 1e-6:
        hard_failures.append("atmosphere.seasonal_wind contains non-tangent vectors")
    if circulation["max_background_wind_normal_component"] > 1e-6:
        hard_failures.append("atmosphere.background_seasonal_wind contains non-tangent vectors")
    if circulation["max_thermal_wind_normal_component"] > 1e-6:
        hard_failures.append("atmosphere.thermal_wind_anomaly contains non-tangent vectors")
    if circulation["max_orographic_wind_normal_component"] > 1e-6:
        hard_failures.append("atmosphere.orographic_wind_anomaly contains non-tangent vectors")
    if ocean_currents["invalid_current_shape"]:
        hard_failures.append("ocean.currents must have shape (n_cells, 3)")
    if ocean_currents["invalid_current_heat_transport_shape"]:
        hard_failures.append("ocean.current_heat_transport must have shape (n_cells,)")
    if ocean_currents["invalid_upwelling_shape"]:
        hard_failures.append("ocean.upwelling must have shape (n_cells,)")
    if ocean_currents["invalid_basin_id_shape"]:
        hard_failures.append("ocean.basin_id must have shape (n_cells,)")
    if ocean_currents["invalid_seasonal_sst_shape"]:
        hard_failures.append("climate.seasonal_sst must have shape (4, n_cells)")
    if ocean_currents["invalid_ocean_heat_flux_shape"]:
        hard_failures.append("climate.ocean_heat_flux must have shape (n_cells,)")
    if ocean_currents["invalid_coupling_residual_shape"]:
        hard_failures.append("climate.coupling_residual must have shape (n_cells,)")
    if ocean_currents["invalid_current_streamfunction_shape"]:
        hard_failures.append("ocean.current_streamfunction must have shape (n_cells,)")
    if ocean_currents["invalid_gyre_id_shape"]:
        hard_failures.append("ocean.gyre_id must have shape (n_cells,)")
    if ocean_currents["invalid_boundary_current_type_shape"]:
        hard_failures.append("ocean.boundary_current_type must have shape (n_cells,)")
    if ocean_currents["invalid_strait_exchange_shape"]:
        hard_failures.append("ocean.strait_exchange must have shape (n_cells,)")
    if ocean_currents["invalid_wind_stress_current_response_shape"]:
        hard_failures.append(
            "ocean.wind_stress_current_response must have shape (n_cells, 3)")
    if ocean_currents["invalid_sst_anomaly_shape"]:
        hard_failures.append("ocean.sst_anomaly must have shape (n_cells,)")
    if ocean_currents["invalid_solved_mask_shape"]:
        hard_failures.append("ocean.solved_mask must have shape (n_cells,)")
    if ocean_currents["nonfinite_current_cells"] > 0:
        hard_failures.append("non-finite ocean.currents cells")
    if ocean_currents["nonfinite_current_heat_transport_cells"] > 0:
        hard_failures.append("non-finite ocean.current_heat_transport cells")
    if ocean_currents["nonfinite_upwelling_cells"] > 0:
        hard_failures.append("non-finite ocean.upwelling cells")
    if ocean_currents["nonfinite_seasonal_sst_cells"] > 0:
        hard_failures.append("non-finite climate.seasonal_sst cells")
    if ocean_currents["nonfinite_ocean_heat_flux_cells"] > 0:
        hard_failures.append("non-finite climate.ocean_heat_flux cells")
    if ocean_currents["nonfinite_coupling_residual_cells"] > 0:
        hard_failures.append("non-finite climate.coupling_residual cells")
    if ocean_currents["nonfinite_current_streamfunction_cells"] > 0:
        hard_failures.append("non-finite ocean.current_streamfunction cells")
    if ocean_currents["nonfinite_wind_stress_current_response_cells"] > 0:
        hard_failures.append(
            "non-finite ocean.wind_stress_current_response cells")
    if ocean_currents["nonfinite_sst_anomaly_cells"] > 0:
        hard_failures.append("non-finite ocean.sst_anomaly cells")
    if ocean_currents["current_over_solved_land_speed_max_mps"] > 1e-8:
        hard_failures.append("ocean.currents contains nonzero currents outside solved ocean")
    if ocean_currents["max_current_normal_component"] > 1e-6:
        hard_failures.append("ocean.currents contains non-tangent vectors")
    if ocean_currents["wind_stress_response_land_speed_max_mps"] > 1e-8:
        hard_failures.append(
            "ocean.wind_stress_current_response contains nonzero land response")
    if ocean_currents["wind_stress_response_max_normal_component"] > 1e-6:
        hard_failures.append(
            "ocean.wind_stress_current_response contains non-tangent vectors")
    if ocean_currents["current_streamfunction_land_abs_max"] > 1e-8:
        hard_failures.append("ocean.current_streamfunction must be zero over land")
    for key, message in [
        ("invalid_continent_id_shape", "climate.continent_id must have shape (n_cells,)"),
        ("invalid_continent_interiority_shape",
         "climate.continent_interiority must have shape (n_cells,)"),
        ("invalid_coast_orientation_shape",
         "climate.coast_orientation must have shape (n_cells, 3)"),
        ("invalid_coast_distance_shape", "climate.coast_distance must have shape (n_cells,)"),
        ("invalid_coast_strength_shape", "climate.coast_strength must have shape (n_cells,)"),
        ("invalid_coast_facing_east_shape",
         "climate.coast_facing_east must have shape (n_cells,)"),
        ("invalid_shelf_index_shape", "ocean.shelf_index must have shape (n_cells,)"),
        ("invalid_strait_index_shape", "ocean.strait_index must have shape (n_cells,)"),
        ("invalid_barrier_index_shape", "terrain.barrier_index must have shape (n_cells,)"),
        ("invalid_wind_gap_index_shape", "terrain.wind_gap_index must have shape (n_cells,)"),
    ]:
        if geography[key]:
            hard_failures.append(message)
    if geography["nonfinite_geography_cells"] > 0:
        hard_failures.append("non-finite shared geography primitive cells")
    if geography["land_without_continent_id_cells"] > 0:
        hard_failures.append("land cells missing climate.continent_id")
    if geography["ocean_with_continent_id_cells"] > 0:
        hard_failures.append("ocean cells labelled with climate.continent_id")
    if geography["ocean_without_basin_id_cells"] > 0:
        hard_failures.append("ocean cells missing ocean.basin_id")
    if geography["land_with_basin_id_cells"] > 0:
        hard_failures.append("land cells labelled with ocean.basin_id")
    if geography["max_coast_orientation_normal_component"] > 1e-6:
        hard_failures.append("climate.coast_orientation contains non-tangent vectors")
    if geography["coast_orientation_noncoastal_speed_p95"] > 1e-8:
        hard_failures.append("climate.coast_orientation is nonzero away from coastal land")
    if geography["coast_strength_land_noncoast_max"] > 1e-8:
        hard_failures.append("climate.coast_strength is nonzero away from coastal land")
    if geography["strait_land_max"] > 1e-8:
        hard_failures.append("ocean.strait_index is nonzero on land")

    nonfinite = _nonfinite_numbers(detail)
    if nonfinite:
        hard_failures.append("non-finite diagnostic values: " + ", ".join(nonfinite[:8]))

    detail["hard_failures"] = hard_failures
    return CheckResult("diagnostics.climate", len(hard_failures) == 0, detail)


def check_regression(spec_factory, cells: int = 3000) -> CheckResult:
    """Run a short world twice with the same seed; key fields must match."""
    from aevum.engine import Engine

    def run_once():
        spec = spec_factory()
        spec.grid_cells = cells
        spec.t_end_myr = 1500.0
        eng = Engine.build(spec)
        eng.run(n_frames=4)
        return eng.world.get_field("terrain.elevation_m").copy()

    a = run_once()
    b = run_once()
    same = np.allclose(a, b)
    return CheckResult("regression.determinism", bool(same),
                       {"max_abs_diff": float(np.max(np.abs(a - b)))})


def run_all(engine) -> list[CheckResult]:
    return [check_conservation(engine), check_topology(engine),
            check_causality(engine), check_tectonic_diagnostics(engine),
            check_climate_diagnostics(engine), check_compiler_consistency(engine),
            check_earthlike_plausibility(engine)]
