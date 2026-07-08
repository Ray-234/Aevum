"""Terminal terrain -> climate/biome post-processing.

This runner is for the handoff after plate/terrain generation is accepted.  It
rebuilds a deterministic terminal terrain from preset/seed/cell count using the
existing tectonics stack, then runs climate and a static biome classifier at the
final time.  Tectonics and terrain code are not modified or re-entered after the
terminal state is produced.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import json
import os
import time
from pathlib import Path
from typing import Any

import numpy as np

from aevum import render, validation
from aevum.core.module import StepResult
from aevum.core.rng import RNGKey
from aevum.engine import Engine
from aevum.modules.biosphere import BiosphereModule
from aevum.modules.climate import ClimateModule
from aevum.spec.presets import get_preset


TERRAIN_INPUT_MODULES = {"stellar", "interior", "impacts", "tectonics", "terrain"}


@dataclass(frozen=True)
class TerminalClimateJob:
    preset: str
    label: str
    seed: int


DEFAULT_SIX_WORLD_JOBS = (
    TerminalClimateJob("waterworld", "waterworld_seed7", 7),
    TerminalClimateJob("waterworld", "waterworld_seed707", 707),
    TerminalClimateJob("earthlike", "earthlike_seed42", 42),
    TerminalClimateJob("earthlike", "earthlike_seed909", 909),
    TerminalClimateJob("arid", "arid_seed101", 101),
    TerminalClimateJob("arid", "arid_seed1001", 1001),
)


@dataclass(frozen=True)
class TerminalClimateConfig:
    jobs: tuple[TerminalClimateJob, ...] = DEFAULT_SIX_WORLD_JOBS
    cells: int = 8000
    t_end_myr: float = 4500.0
    frames: int = 4
    max_workers: int = 1
    render_assets: bool = True


def _json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _lat_band_residual_abs_p95(
    lat: np.ndarray,
    area: np.ndarray,
    values: np.ndarray,
    mask: np.ndarray,
    *,
    band_width_deg: float = 10.0,
) -> float | None:
    lat = np.asarray(lat, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values) & np.isfinite(area)
    if int(np.count_nonzero(mask)) < 8:
        return None
    bins = np.arange(-90.0, 90.0 + band_width_deg, band_width_deg)
    predicted = np.full(values.shape, np.nan, dtype=np.float64)
    for idx, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
        if idx == len(bins) - 2:
            band = mask & (lat >= lo) & (lat <= hi)
        else:
            band = mask & (lat >= lo) & (lat < hi)
        if not band.any():
            continue
        weights = area[band]
        predicted[band] = float(np.average(values[band], weights=weights))
    valid = mask & np.isfinite(predicted)
    if int(np.count_nonzero(valid)) < 8:
        return None
    return float(np.percentile(np.abs(values[valid] - predicted[valid]), 95))


def _apply_step_result(world, result: StepResult) -> None:
    delta = result.state_delta
    for key, value in delta.get("fields", {}).items():
        world.fields[key] = np.asarray(value)
    for key, value in delta.get("globals", {}).items():
        world.globals[key] = float(value)
    for key, value in delta.get("networks", {}).items():
        world.networks[key] = value
    for key, value in delta.get("objects", {}).items():
        world.objects[key] = value


def _static_npp(world) -> np.ndarray:
    T = world.get_field("climate.surface_temperature", 288.0)
    precip = world.get_field("climate.precipitation", 500.0)
    ocean = world.ocean_mask()
    tc = T - 273.15
    npp_t = 1.0 / (1.0 + np.exp(1.315 - 0.119 * tc))
    npp_w = 1.0 - np.exp(-0.000664 * np.clip(precip, 0.0, None))
    npp = np.where(ocean, 0.65 * npp_t, np.minimum(npp_t, npp_w))
    return np.clip(1.2 * npp, 0.0, 5.0)


def _run_terminal_climate(world, seed: int) -> dict[str, Any]:
    climate = ClimateModule()
    key = RNGKey(seed, "terminal_climate", world.time_myr, 1)
    climate.init_state(world, key)
    result = climate.step(
        world,
        world.time_myr,
        40.0,
        {"terminal_postprocess": True},
        key,
    )
    _apply_step_result(world, result)

    npp = _static_npp(world)
    biosphere = BiosphereModule()
    biosphere.init_state(world, RNGKey(seed, "terminal_biome", world.time_myr, 2))
    world.fields["biosphere.npp"] = npp
    world.fields["biosphere.biome"] = biosphere._biomes(world, npp)
    return result.diagnostics


def _write_arrays(world, outdir: Path) -> Path:
    path = outdir / "terminal_climate_arrays.npz"
    keys = [
        "terrain.elevation_m",
        "crust.type",
        "tectonics.plate_id",
        "climate.surface_temperature",
        "climate.seasonal_temperature",
        "climate.seasonal_sst",
        "climate.ocean_heat_flux",
        "climate.ocean_evaporation_feedback",
        "climate.coupling_residual",
        "climate.seasonal_insolation_anomaly",
        "climate.surface_heat_capacity_class",
        "climate.land_thermal_anomaly",
        "climate.ocean_mixed_layer_thermal_anomaly",
        "climate.elevation_lapse_cooling",
        "climate.snow_ice_albedo_support",
        "climate.sst_gradient_support",
        "climate.same_latitude_sst_anomaly",
        "climate.land_sea_thermal_contrast",
        "climate.precipitation",
        "climate.evaporation",
        "climate.runoff",
        "climate.seasonal_precipitation",
        "climate.precipitation_seasonality",
        "climate.monsoon_index",
        "climate.dry_season_length",
        "climate.moisture_convergence",
        "climate.orographic_precipitation",
        "climate.monsoon_rainfall_corridor",
        "climate.storm_track_rainfall_corridor",
        "climate.rain_shadow_index",
        "climate.regional_precipitation_response",
        "atmosphere.moisture_flow_source",
        "atmosphere.moisture_flow_pathway",
        "atmosphere.moisture_source_basin_id",
        "climate.moisture_flow_network_id",
        "climate.moisture_flow_precipitation_response",
        "climate.moisture_budget_region_id",
        "climate.precipitation_response_region_id",
        "climate.receiver_catchment_id",
        "climate.source_basin_supply_index",
        "climate.receiver_catchment_supply_balance",
        "climate.receiver_supply_precipitation_feedback",
        "climate.hydro_coupling_residual",
        "climate.hydro_feedback_iteration_delta",
        "climate.wet_season_peak",
        "climate.continent_id",
        "climate.continent_interiority",
        "climate.coast_distance",
        "climate.coast_strength",
        "atmosphere.moisture_access",
        "atmosphere.monsoon_potential",
        "atmosphere.seasonal_wind",
        "atmosphere.background_seasonal_wind",
        "atmosphere.thermal_wind_anomaly",
        "atmosphere.orographic_wind_anomaly",
        "atmosphere.geographic_circulation_index",
        "atmosphere.itcz_intensity",
        "atmosphere.storm_track_intensity",
        "atmosphere.land_sea_pressure_proxy",
        "atmosphere.seasonal_pressure_proxy",
        "atmosphere.pressure_center_support",
        "atmosphere.pressure_center_id",
        "atmosphere.stationary_wave_pressure_support",
        "atmosphere.pressure_genesis_source",
        "atmosphere.pressure_genesis_wave_transfer",
        "atmosphere.ocean_pressure_low_source_support",
        "atmosphere.ocean_pressure_high_source_support",
        "atmosphere.land_pressure_source_support",
        "atmosphere.terrain_pressure_wave_source_support",
        "atmosphere.precipitation_pressure_feedback",
        "atmosphere.hydro_coupled_wind_anomaly",
        "atmosphere.source_ocean_warmth",
        "atmosphere.terrain_blocking",
        "terrain.barrier_index",
        "terrain.wind_gap_index",
        "ocean.basin_id",
        "ocean.shelf_index",
        "ocean.strait_index",
        "ocean.solved_mask",
        "ocean.currents",
        "ocean.current_heat_transport",
        "ocean.upwelling",
        "ocean.gyre_id",
        "ocean.current_streamfunction",
        "ocean.boundary_current_type",
        "ocean.strait_exchange",
        "ocean.wind_stress_current_response",
        "ocean.sst_anomaly",
        "cryosphere.sea_ice",
        "cryosphere.seasonal_sea_ice",
        "cryosphere.seasonal_snow",
        "cryosphere.snow_persistence",
        "cryosphere.ice_sheet",
        "climate.seasonal_cloud_albedo_proxy",
        "climate.cloud_albedo_proxy",
        "biosphere.vegetation_climate_feedback",
        "biosphere.biome",
        "biosphere.npp",
    ]
    arrays: dict[str, np.ndarray] = {
        "lat": world.grid.lat,
        "lon": world.grid.lon,
        "cell_area": world.grid.cell_area,
        "sea_level_m": np.asarray([world.sea_level], dtype=np.float64),
    }
    for key in keys:
        if key in world.fields:
            arrays[key.replace(".", "__")] = np.asarray(world.fields[key])
    np.savez_compressed(path, **arrays)
    return path


def _write_hydroclimate_regions(world, outdir: Path) -> Path | None:
    regions = world.objects.get("climate.hydroclimate_regions", [])
    if not isinstance(regions, list):
        regions = []
    if not regions:
        return None
    kind_counts: dict[str, int] = {}
    season_counts: dict[str, int] = {}
    sanitized: list[dict[str, Any]] = []
    for obj in regions:
        if not isinstance(obj, dict):
            continue
        row = dict(obj)
        kind = str(row.get("kind", "unknown"))
        season = str(row.get("season", "unknown"))
        kind_counts[kind] = kind_counts.get(kind, 0) + 1
        season_counts[season] = season_counts.get(season, 0) + 1
        sanitized.append(row)
    payload = {
        "schema": "aevum.hydroclimate_regions.v1",
        "time_myr": float(world.time_myr),
        "preset": world.spec.name,
        "seed": int(world.spec.seed),
        "region_count": int(len(sanitized)),
        "kind_counts": dict(sorted(kind_counts.items())),
        "season_counts": dict(sorted(season_counts.items())),
        "regions": sanitized,
    }
    path = outdir / "hydroclimate_regions.json"
    path.write_text(json.dumps(payload, indent=2, default=_json_default))
    return path


def _write_moisture_flow_networks(world, outdir: Path) -> Path | None:
    networks = world.objects.get("climate.moisture_flow_networks", [])
    if not isinstance(networks, list):
        networks = []
    if not networks:
        return None
    kind_counts: dict[str, int] = {}
    season_counts: dict[str, int] = {}
    sanitized: list[dict[str, Any]] = []
    for obj in networks:
        if not isinstance(obj, dict):
            continue
        row = dict(obj)
        kind = str(row.get("kind", "unknown"))
        season = str(row.get("season", "unknown"))
        kind_counts[kind] = kind_counts.get(kind, 0) + 1
        season_counts[season] = season_counts.get(season, 0) + 1
        sanitized.append(row)
    payload = {
        "schema": "aevum.moisture_flow_networks.v1",
        "time_myr": float(world.time_myr),
        "preset": world.spec.name,
        "seed": int(world.spec.seed),
        "network_count": int(len(sanitized)),
        "kind_counts": dict(sorted(kind_counts.items())),
        "season_counts": dict(sorted(season_counts.items())),
        "networks": sanitized,
    }
    path = outdir / "moisture_flow_networks.json"
    path.write_text(json.dumps(payload, indent=2, default=_json_default))
    return path


def _write_precipitation_response_regions(world, outdir: Path) -> Path | None:
    regions = world.objects.get("climate.precipitation_response_regions", [])
    if not isinstance(regions, list):
        regions = []
    if not regions:
        return None
    kind_counts: dict[str, int] = {}
    season_counts: dict[str, int] = {}
    sanitized: list[dict[str, Any]] = []
    for obj in regions:
        if not isinstance(obj, dict):
            continue
        row = dict(obj)
        kind = str(row.get("kind", "unknown"))
        season = str(row.get("season", "unknown"))
        kind_counts[kind] = kind_counts.get(kind, 0) + 1
        season_counts[season] = season_counts.get(season, 0) + 1
        sanitized.append(row)
    payload = {
        "schema": "aevum.precipitation_response_regions.v1",
        "time_myr": float(world.time_myr),
        "preset": world.spec.name,
        "seed": int(world.spec.seed),
        "region_count": int(len(sanitized)),
        "kind_counts": dict(sorted(kind_counts.items())),
        "season_counts": dict(sorted(season_counts.items())),
        "regions": sanitized,
    }
    path = outdir / "precipitation_response_regions.json"
    path.write_text(json.dumps(payload, indent=2, default=_json_default))
    return path


def _write_receiver_catchments(world, outdir: Path) -> Path | None:
    catchments = world.objects.get("climate.receiver_catchments", [])
    if not isinstance(catchments, list):
        catchments = []
    if not catchments:
        return None
    kind_counts: dict[str, int] = {}
    season_counts: dict[str, int] = {}
    sanitized: list[dict[str, Any]] = []
    for obj in catchments:
        if not isinstance(obj, dict):
            continue
        row = dict(obj)
        kind = str(row.get("kind", "unknown"))
        season = str(row.get("season", "unknown"))
        kind_counts[kind] = kind_counts.get(kind, 0) + 1
        season_counts[season] = season_counts.get(season, 0) + 1
        sanitized.append(row)
    payload = {
        "schema": "aevum.receiver_catchments.v1",
        "time_myr": float(world.time_myr),
        "preset": world.spec.name,
        "seed": int(world.spec.seed),
        "catchment_count": int(len(sanitized)),
        "kind_counts": dict(sorted(kind_counts.items())),
        "season_counts": dict(sorted(season_counts.items())),
        "catchments": sanitized,
    }
    path = outdir / "receiver_catchments.json"
    path.write_text(json.dumps(payload, indent=2, default=_json_default))
    return path


def _write_pressure_centers(world, outdir: Path) -> Path | None:
    centers = world.objects.get("atmosphere.pressure_centers", [])
    if not isinstance(centers, list):
        centers = []
    if not centers:
        return None
    kind_counts: dict[str, int] = {}
    season_counts: dict[str, int] = {}
    domain_counts: dict[str, int] = {}
    sanitized: list[dict[str, Any]] = []
    for obj in centers:
        if not isinstance(obj, dict):
            continue
        row = dict(obj)
        kind = str(row.get("kind", "unknown"))
        season = str(row.get("season", "unknown"))
        domain = str(row.get("domain", "unknown"))
        kind_counts[kind] = kind_counts.get(kind, 0) + 1
        season_counts[season] = season_counts.get(season, 0) + 1
        domain_counts[domain] = domain_counts.get(domain, 0) + 1
        sanitized.append(row)
    payload = {
        "schema": "aevum.pressure_centers.v1",
        "time_myr": float(world.time_myr),
        "preset": world.spec.name,
        "seed": int(world.spec.seed),
        "center_count": int(len(sanitized)),
        "kind_counts": dict(sorted(kind_counts.items())),
        "season_counts": dict(sorted(season_counts.items())),
        "domain_counts": dict(sorted(domain_counts.items())),
        "centers": sanitized,
    }
    path = outdir / "pressure_centers.json"
    path.write_text(json.dumps(payload, indent=2, default=_json_default))
    return path


def _summarize_world(eng: Engine, climate_diag: dict[str, Any], outdir: Path) -> dict[str, Any]:
    w = eng.world
    land = w.land_mask()
    ocean = ~land
    temp_c = w.get_field("climate.surface_temperature") - 273.15
    precip = w.get_field("climate.precipitation")
    biome = w.get_field("biosphere.biome", 0.0).astype(int)
    seasonal_precip = w.fields.get("climate.seasonal_precipitation")
    seasonal_sst = w.fields.get("climate.seasonal_sst")
    ocean_heat_flux = w.fields.get("climate.ocean_heat_flux")
    current_heat_transport = w.fields.get("ocean.current_heat_transport")
    seasonal_sea_ice = w.fields.get("cryosphere.seasonal_sea_ice")
    snow_persistence = w.fields.get("cryosphere.snow_persistence")
    cloud_albedo = w.fields.get("climate.cloud_albedo_proxy")
    vegetation_feedback = w.fields.get("biosphere.vegetation_climate_feedback")
    flow_response = w.fields.get("climate.moisture_flow_precipitation_response")
    source_basin_id = w.fields.get("atmosphere.moisture_source_basin_id")
    moisture_budget_region = w.fields.get("climate.moisture_budget_region_id")
    precipitation_response_region_id = w.fields.get(
        "climate.precipitation_response_region_id")
    receiver_catchment_id = w.fields.get("climate.receiver_catchment_id")
    source_basin_supply_index = w.fields.get("climate.source_basin_supply_index")
    receiver_supply_balance = w.fields.get(
        "climate.receiver_catchment_supply_balance")
    receiver_feedback = w.fields.get("climate.receiver_supply_precipitation_feedback")
    hydro_regions = w.objects.get("climate.hydroclimate_regions", [])
    if not isinstance(hydro_regions, list):
        hydro_regions = []
    moisture_networks = w.objects.get("climate.moisture_flow_networks", [])
    if not isinstance(moisture_networks, list):
        moisture_networks = []
    precipitation_response_regions = w.objects.get(
        "climate.precipitation_response_regions", [])
    if not isinstance(precipitation_response_regions, list):
        precipitation_response_regions = []
    receiver_catchments = w.objects.get("climate.receiver_catchments", [])
    if not isinstance(receiver_catchments, list):
        receiver_catchments = []
    pressure_centers = w.objects.get("atmosphere.pressure_centers", [])
    if not isinstance(pressure_centers, list):
        pressure_centers = []
    hydro_region_kind_counts: dict[str, int] = {}
    hydro_region_season_counts: dict[str, int] = {}
    for obj in hydro_regions:
        kind = str(obj.get("kind", "unknown"))
        season = str(obj.get("season", "unknown"))
        hydro_region_kind_counts[kind] = hydro_region_kind_counts.get(kind, 0) + 1
        hydro_region_season_counts[season] = hydro_region_season_counts.get(season, 0) + 1
    hydro_regions_path = _write_hydroclimate_regions(w, outdir)
    moisture_networks_path = _write_moisture_flow_networks(w, outdir)
    precipitation_response_regions_path = _write_precipitation_response_regions(
        w, outdir)
    receiver_catchments_path = _write_receiver_catchments(w, outdir)
    pressure_centers_path = _write_pressure_centers(w, outdir)
    moisture_network_kind_counts: dict[str, int] = {}
    moisture_network_season_counts: dict[str, int] = {}
    for obj in moisture_networks:
        kind = str(obj.get("kind", "unknown"))
        season = str(obj.get("season", "unknown"))
        moisture_network_kind_counts[kind] = (
            moisture_network_kind_counts.get(kind, 0) + 1)
        moisture_network_season_counts[season] = (
            moisture_network_season_counts.get(season, 0) + 1)
    precipitation_response_kind_counts: dict[str, int] = {}
    precipitation_response_season_counts: dict[str, int] = {}
    for obj in precipitation_response_regions:
        kind = str(obj.get("kind", "unknown"))
        season = str(obj.get("season", "unknown"))
        precipitation_response_kind_counts[kind] = (
            precipitation_response_kind_counts.get(kind, 0) + 1)
        precipitation_response_season_counts[season] = (
            precipitation_response_season_counts.get(season, 0) + 1)
    receiver_catchment_kind_counts: dict[str, int] = {}
    receiver_catchment_season_counts: dict[str, int] = {}
    for obj in receiver_catchments:
        kind = str(obj.get("kind", "unknown"))
        season = str(obj.get("season", "unknown"))
        receiver_catchment_kind_counts[kind] = (
            receiver_catchment_kind_counts.get(kind, 0) + 1)
        receiver_catchment_season_counts[season] = (
            receiver_catchment_season_counts.get(season, 0) + 1)
    pressure_center_kind_counts: dict[str, int] = {}
    pressure_center_season_counts: dict[str, int] = {}
    pressure_center_domain_counts: dict[str, int] = {}
    for obj in pressure_centers:
        kind = str(obj.get("kind", "unknown"))
        season = str(obj.get("season", "unknown"))
        domain = str(obj.get("domain", "unknown"))
        pressure_center_kind_counts[kind] = pressure_center_kind_counts.get(kind, 0) + 1
        pressure_center_season_counts[season] = (
            pressure_center_season_counts.get(season, 0) + 1)
        pressure_center_domain_counts[domain] = (
            pressure_center_domain_counts.get(domain, 0) + 1)
    response_diag = {}
    if isinstance(climate_diag, dict):
        maybe_response_diag = climate_diag.get("moisture_flow_precipitation_response", {})
        if isinstance(maybe_response_diag, dict):
            response_diag = maybe_response_diag
    summary = {
        "preset": w.spec.name,
        "seed": int(w.spec.seed),
        "cells": int(w.grid.n),
        "time_myr": float(w.time_myr),
        "active_terrain_modules": sorted(TERRAIN_INPUT_MODULES),
        "land_fraction": float(w.land_fraction()),
        "mean_temperature_C": float(np.average(temp_c, weights=w.grid.cell_area)),
        "surface_temperature_zonal_residual_abs_p95_C": (
            _lat_band_residual_abs_p95(
                w.grid.lat, w.grid.cell_area, temp_c,
                np.isfinite(temp_c),
            )
        ),
        "seasonal_sst_zonal_residual_abs_p95_C": (
            _lat_band_residual_abs_p95(
                w.grid.lat,
                w.grid.cell_area,
                np.asarray(seasonal_sst, dtype=np.float64).mean(axis=0) - 273.15,
                ocean,
            )
            if seasonal_sst is not None
            and np.asarray(seasonal_sst).shape == (4, w.grid.n)
            and ocean.any()
            else None
        ),
        "ocean_heat_flux_abs_p95_C": (
            float(np.percentile(
                np.abs(np.asarray(ocean_heat_flux, dtype=np.float64)[ocean]), 95))
            if ocean_heat_flux is not None
            and np.asarray(ocean_heat_flux).shape == (w.grid.n,)
            and ocean.any()
            else None
        ),
        "current_heat_transport_abs_p95_C": (
            float(np.percentile(
                np.abs(np.asarray(current_heat_transport, dtype=np.float64)[ocean]), 95))
            if current_heat_transport is not None
            and np.asarray(current_heat_transport).shape == (w.grid.n,)
            and ocean.any()
            else None
        ),
        "seasonal_sea_ice_ocean_p95": (
            float(np.percentile(
                np.asarray(seasonal_sea_ice, dtype=np.float64)[:, ocean], 95))
            if seasonal_sea_ice is not None
            and np.asarray(seasonal_sea_ice).shape == (4, w.grid.n)
            and ocean.any()
            else None
        ),
        "snow_persistence_land_p95": (
            float(np.percentile(
                np.asarray(snow_persistence, dtype=np.float64)[land], 95))
            if snow_persistence is not None
            and np.asarray(snow_persistence).shape == (w.grid.n,)
            and land.any()
            else None
        ),
        "cloud_albedo_proxy_p50": (
            float(np.percentile(np.asarray(cloud_albedo, dtype=np.float64), 50))
            if cloud_albedo is not None
            and np.asarray(cloud_albedo).shape == (w.grid.n,)
            else None
        ),
        "vegetation_climate_feedback_land_p50": (
            float(np.percentile(
                np.asarray(vegetation_feedback, dtype=np.float64)[land], 50))
            if vegetation_feedback is not None
            and np.asarray(vegetation_feedback).shape == (w.grid.n,)
            and land.any()
            else None
        ),
        "land_mean_temperature_C": (
            float(np.average(temp_c[land], weights=w.grid.cell_area[land]))
            if land.any() else None
        ),
        "ocean_mean_temperature_C": (
            float(np.average(temp_c[ocean], weights=w.grid.cell_area[ocean]))
            if ocean.any() else None
        ),
        "mean_precipitation_mm_yr": float(np.average(precip, weights=w.grid.cell_area)),
        "land_precip_p50_mm_yr": (
            float(np.percentile(precip[land], 50)) if land.any() else None
        ),
        "land_precip_p90_mm_yr": (
            float(np.percentile(precip[land], 90)) if land.any() else None
        ),
        "land_monsoon_index_p90": (
            float(np.percentile(w.fields["climate.monsoon_index"][land], 90))
            if land.any() and "climate.monsoon_index" in w.fields else None
        ),
        "precip_seasonality_land_p75": (
            float(np.percentile(w.fields["climate.precipitation_seasonality"][land], 75))
            if land.any() and "climate.precipitation_seasonality" in w.fields else None
        ),
        "biome_counts": {
            str(code): int((biome == code).sum())
            for code in sorted(int(x) for x in np.unique(biome))
        },
        "climate_step_diagnostics": climate_diag,
        "climate_diagnostics": validation.climate_diagnostics(eng),
        "hydroclimate_region_object_count": int(len(hydro_regions)),
        "hydroclimate_region_kind_counts": dict(sorted(hydro_region_kind_counts.items())),
        "hydroclimate_region_season_counts": dict(sorted(hydro_region_season_counts.items())),
        "hydroclimate_regions_json": str(hydro_regions_path) if hydro_regions_path else None,
        "moisture_flow_network_object_count": int(len(moisture_networks)),
        "moisture_flow_network_kind_counts": dict(
            sorted(moisture_network_kind_counts.items())),
        "moisture_flow_network_season_counts": dict(
            sorted(moisture_network_season_counts.items())),
        "moisture_flow_networks_json": (
            str(moisture_networks_path) if moisture_networks_path else None
        ),
        "precipitation_response_region_object_count": int(
            len(precipitation_response_regions)),
        "precipitation_response_region_kind_counts": dict(
            sorted(precipitation_response_kind_counts.items())),
        "precipitation_response_region_season_counts": dict(
            sorted(precipitation_response_season_counts.items())),
        "precipitation_response_regions_json": (
            str(precipitation_response_regions_path)
            if precipitation_response_regions_path else None
        ),
        "precipitation_response_region_count_p50": (
            float(np.percentile([
                len([
                    x for x in np.unique(
                        np.asarray(precipitation_response_region_id)[season, land])
                    if np.isfinite(x) and int(x) > 0
                ])
                for season in range(4)
            ], 50))
            if precipitation_response_region_id is not None and land.any()
            and np.asarray(precipitation_response_region_id).shape == (4, w.grid.n)
            else None
        ),
        "receiver_catchment_object_count": int(len(receiver_catchments)),
        "receiver_catchment_kind_counts": dict(
            sorted(receiver_catchment_kind_counts.items())),
        "receiver_catchment_season_counts": dict(
            sorted(receiver_catchment_season_counts.items())),
        "receiver_catchments_json": (
            str(receiver_catchments_path) if receiver_catchments_path else None
        ),
        "pressure_center_object_count": int(len(pressure_centers)),
        "pressure_center_kind_counts": dict(
            sorted(pressure_center_kind_counts.items())),
        "pressure_center_season_counts": dict(
            sorted(pressure_center_season_counts.items())),
        "pressure_center_domain_counts": dict(
            sorted(pressure_center_domain_counts.items())),
        "pressure_centers_json": (
            str(pressure_centers_path) if pressure_centers_path else None
        ),
        "pressure_center_count_p50": (
            float(np.percentile([
                len([
                    x for x in np.unique(
                        np.asarray(w.fields["atmosphere.pressure_center_id"])[season])
                    if np.isfinite(x) and int(x) > 0
                ])
                for season in range(4)
            ], 50))
            if "atmosphere.pressure_center_id" in w.fields
            and np.asarray(w.fields["atmosphere.pressure_center_id"]).shape == (4, w.grid.n)
            else None
        ),
        "receiver_catchment_count_p50": (
            float(np.percentile([
                len([
                    x for x in np.unique(
                        np.asarray(receiver_catchment_id)[season, land])
                    if np.isfinite(x) and int(x) > 0
                ])
                for season in range(4)
            ], 50))
            if receiver_catchment_id is not None and land.any()
            and np.asarray(receiver_catchment_id).shape == (4, w.grid.n)
            else None
        ),
        "source_basin_supply_index_land_p50": (
            float(np.percentile(np.asarray(source_basin_supply_index)[:, land], 50))
            if source_basin_supply_index is not None and land.any()
            and np.asarray(source_basin_supply_index).shape == (4, w.grid.n)
            else None
        ),
        "receiver_catchment_supply_balance_land_p50": (
            float(np.percentile(np.asarray(receiver_supply_balance)[:, land], 50))
            if receiver_supply_balance is not None and land.any()
            and np.asarray(receiver_supply_balance).shape == (4, w.grid.n)
            else None
        ),
        "receiver_supply_precipitation_feedback_land_p05": (
            float(np.percentile(np.asarray(receiver_feedback)[:, land].ravel(), 5))
            if receiver_feedback is not None and land.any()
            and np.asarray(receiver_feedback).shape == (4, w.grid.n) else None
        ),
        "receiver_supply_precipitation_feedback_land_p95": (
            float(np.percentile(np.asarray(receiver_feedback)[:, land].ravel(), 95))
            if receiver_feedback is not None and land.any()
            and np.asarray(receiver_feedback).shape == (4, w.grid.n) else None
        ),
        "moisture_flow_precipitation_response_land_p05": (
            float(np.percentile(np.asarray(flow_response)[:, land].ravel(), 5))
            if flow_response is not None and land.any()
            and np.asarray(flow_response).shape == (4, w.grid.n) else None
        ),
        "moisture_flow_precipitation_response_land_p95": (
            float(np.percentile(np.asarray(flow_response)[:, land].ravel(), 95))
            if flow_response is not None and land.any()
            and np.asarray(flow_response).shape == (4, w.grid.n) else None
        ),
        "moisture_source_basin_attributed_land_fraction": (
            float(np.count_nonzero(
                np.asarray(source_basin_id)[:, land] >= 0.0
            ) / max(4 * int(np.count_nonzero(land)), 1))
            if source_basin_id is not None and land.any()
            and np.asarray(source_basin_id).shape == (4, w.grid.n) else None
        ),
        "moisture_budget_region_count_p50": (
            float(np.percentile([
                len([
                    x for x in np.unique(np.asarray(moisture_budget_region)[season, land])
                    if np.isfinite(x) and int(x) > 0
                ])
                for season in range(4)
            ], 50))
            if moisture_budget_region is not None and land.any()
            and np.asarray(moisture_budget_region).shape == (4, w.grid.n) else None
        ),
        "moisture_budget_base_region_count": (
            float(response_diag.get("budget_base_region_count"))
            if "budget_base_region_count" in response_diag else None
        ),
        "moisture_budget_sector_split_count_p50": (
            float(response_diag.get("budget_sector_split_count_p50"))
            if "budget_sector_split_count_p50" in response_diag else None
        ),
        "assets_dir": str(outdir),
    }
    if seasonal_precip is not None:
        summary["seasonal_precip_shape"] = list(np.asarray(seasonal_precip).shape)
    return summary


def run_terminal_climate_job(
    job: TerminalClimateJob,
    outdir: Path,
    *,
    cells: int = 8000,
    t_end_myr: float = 4500.0,
    frames: int = 4,
    render_assets: bool = True,
) -> dict[str, Any]:
    t0 = time.time()
    outdir.mkdir(parents=True, exist_ok=True)

    spec = get_preset(job.preset)
    spec.seed = int(job.seed)
    spec.grid_cells = int(cells)
    spec.t_end_myr = float(t_end_myr)

    eng = Engine.build(spec)
    eng.scheduler.modules = [
        sm for sm in eng.scheduler.modules if sm.module.name in TERRAIN_INPUT_MODULES
    ]
    eng.run(n_frames=frames, progress=False)
    climate_diag = _run_terminal_climate(eng.world, spec.seed)
    if render_assets:
        render.render_world(eng.world, outdir)
    arrays_path = _write_arrays(eng.world, outdir)
    summary = _summarize_world(eng, climate_diag, outdir)
    summary["runtime_seconds"] = float(time.time() - t0)
    summary["arrays"] = str(arrays_path)
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    return summary


def _run_one_from_tuple(args: tuple[TerminalClimateJob, str, int, float, int, bool]):
    job, root, cells, t_end_myr, frames, render_assets = args
    return run_terminal_climate_job(
        job,
        Path(root) / job.label,
        cells=cells,
        t_end_myr=t_end_myr,
        frames=frames,
        render_assets=render_assets,
    )


def run_terminal_climate_biome_batch(
    config: TerminalClimateConfig,
    outdir: Path,
) -> dict[str, Any]:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    outdir.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, Any]] = []
    args = [
        (job, str(outdir), config.cells, config.t_end_myr, config.frames,
         config.render_assets)
        for job in config.jobs
    ]
    if config.max_workers <= 1:
        for item in args:
            summaries.append(_run_one_from_tuple(item))
    else:
        with ProcessPoolExecutor(max_workers=int(config.max_workers)) as pool:
            futures = {pool.submit(_run_one_from_tuple, item): item[0] for item in args}
            for future in as_completed(futures):
                summaries.append(future.result())
    summaries.sort(key=lambda item: str(item["assets_dir"]))
    summary = {
        "schema": "aevum.terminal_climate_biome_batch.v1",
        "cells": int(config.cells),
        "t_end_myr": float(config.t_end_myr),
        "job_count": len(config.jobs),
        "summaries": summaries,
    }
    (outdir / "terminal_climate_biome_summary.json").write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    return summary
