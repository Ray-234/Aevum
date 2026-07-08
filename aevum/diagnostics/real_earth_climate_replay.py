"""Run the current climate/biome stack on frozen real-Earth topography.

This is the Earth-fitting baseline: plate and terrain generation are not
entered.  The real-Earth reference grid supplies the surface, and the current
climate plus static biome layers are replayed on that surface.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from aevum import render, validation
from aevum.core.grid import SphereGrid
from aevum.core.state import WorldState
from aevum.diagnostics.terminal_climate_biome import (
    _summarize_world,
    _write_arrays,
)
from aevum.diagnostics.terminal_climate_replay import (
    _restore_stellar_forcing,
    _run_replayed_climate,
)
from aevum.spec.presets import get_preset


SCHEMA = "aevum.real_earth_climate_replay.v1"


@dataclass(frozen=True)
class RealEarthClimateReplayConfig:
    earth_reference_npz: Path
    outdir: Path
    preset: str = "earthlike"
    seed: int = 20260706
    render_assets: bool = True


def _json_default(value: Any):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _weighted_mean(values: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values) & np.isfinite(weights)
    if not mask.any():
        return float("nan")
    return float(np.average(values[mask], weights=weights[mask]))


def _percentile(values: np.ndarray, q: float, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    if not mask.any():
        return float("nan")
    return float(np.percentile(values[mask], q))


def _rmse(values: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    return float(np.sqrt(_weighted_mean(values * values, weights, mask)))


def _unit_fraction(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if np.nanmax(arr) > 1.5:
        arr = arr / 100.0
    return np.clip(arr, 0.0, 1.0)


def _circular_lon_abs_delta(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.abs((np.asarray(a) - np.asarray(b) + 180.0) % 360.0 - 180.0)


def _world_from_earth_reference(
    reference_npz: Path,
    *,
    preset: str,
    seed: int,
) -> tuple[WorldState, dict[str, Any]]:
    with np.load(reference_npz, allow_pickle=False) as z:
        lat = np.asarray(z["lat"], dtype=np.float64)
        lon = np.asarray(z["lon"], dtype=np.float64)
        elevation = np.asarray(z["earth__elevation_m"], dtype=np.float64)
        land = np.asarray(z["earth__land_mask"], dtype=bool)

    spec = get_preset(preset)
    spec.seed = int(seed)
    spec.grid_cells = int(lat.size)
    spec.t_end_myr = 4500.0
    grid = SphereGrid.fibonacci(int(lat.size), spec.radius_m)
    world = WorldState(grid=grid, spec=spec, time_myr=spec.t_end_myr)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_g("diagnostics.real_earth_replay", 1.0)
    world.fields["terrain.elevation_m"] = elevation
    world.fields["crust.type"] = land.astype(np.float64)
    world.fields["tectonics.plate_id"] = np.zeros(grid.n, dtype=np.float64)

    derived_land = elevation >= 0.0
    area = grid.cell_area
    mismatch = derived_land != land
    diagnostics = {
        "cells": int(grid.n),
        "grid_lat_max_abs_delta": float(np.max(np.abs(grid.lat - lat))),
        "grid_lon_max_abs_delta": float(np.max(_circular_lon_abs_delta(grid.lon, lon))),
        "land_mask_mismatch_area_fraction": (
            float(np.sum(area[mismatch]) / max(float(np.sum(area)), 1.0e-12))
        ),
    }
    return world, diagnostics


def _earth_reference_fields(reference_npz: Path) -> dict[str, np.ndarray]:
    with np.load(reference_npz, allow_pickle=False) as z:
        out = {key: np.asarray(z[key]) for key in z.files}
    return out


def _pattern_corr(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(a) & np.isfinite(b)
    if int(np.count_nonzero(mask)) < 4:
        return float("nan")
    aa = a[mask] - float(np.mean(a[mask]))
    bb = b[mask] - float(np.mean(b[mask]))
    denom = float(np.sqrt(np.sum(aa * aa) * np.sum(bb * bb)))
    if denom <= 1.0e-12:
        return float("nan")
    return float(np.sum(aa * bb) / denom)


def _lat_band_jump(values: np.ndarray, lat: np.ndarray, area: np.ndarray) -> dict[str, Any]:
    values = np.asarray(values, dtype=np.float64)
    lat = np.asarray(lat, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)
    bins = np.linspace(-90.0, 90.0, 19)
    means: list[float] = []
    labels: list[tuple[float, float]] = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (
            (lat >= lo)
            & (lat < hi if hi < 90.0 else lat <= hi)
            & np.isfinite(values)
            & np.isfinite(area)
        )
        labels.append((float(lo), float(hi)))
        means.append(_weighted_mean(values, area, mask) if mask.any() else float("nan"))

    means_arr = np.asarray(means, dtype=np.float64)
    deltas = np.abs(np.diff(means_arr))
    finite = np.isfinite(deltas)
    if not finite.any():
        return {
            "max_adjacent_lat_band_delta_C": float("nan"),
            "max_adjacent_lat_band_pair_deg": [],
            "lat_band_means_C": means,
        }
    finite_indices = np.flatnonzero(finite)
    idx = int(finite_indices[int(np.argmax(deltas[finite]))])
    return {
        "max_adjacent_lat_band_delta_C": float(deltas[idx]),
        "max_adjacent_lat_band_pair_deg": [labels[idx], labels[idx + 1]],
        "lat_band_means_C": means,
    }


def _earth_and_replay_metrics(
    world: WorldState,
    reference: dict[str, np.ndarray],
) -> dict[str, Any]:
    area = np.asarray(world.grid.cell_area, dtype=np.float64)
    lat = np.asarray(world.grid.lat, dtype=np.float64)
    land = np.asarray(reference["earth__land_mask"], dtype=bool)
    ocean = ~land

    earth_land_temp = np.asarray(reference["earth__annual_temperature_C"], dtype=np.float64)
    earth_sst = np.asarray(reference["earth__annual_sst_C"], dtype=np.float64)
    earth_surface = np.where(land, earth_land_temp, earth_sst)
    earth_precip = np.asarray(reference["earth__annual_precip_mm"], dtype=np.float64)
    earth_current_speed = np.asarray(
        reference["earth__annual_surface_current_speed_m_s"], dtype=np.float64)
    earth_sea_ice = _unit_fraction(reference["earth__annual_sea_ice_concentration_pct"])

    replay_surface = (
        np.asarray(world.fields["climate.surface_temperature"], dtype=np.float64)
        - 273.15
    )
    replay_precip = np.asarray(world.fields["climate.precipitation"], dtype=np.float64)
    replay_currents = np.asarray(world.fields["ocean.currents"], dtype=np.float64)
    replay_current_speed = np.linalg.norm(replay_currents, axis=1)
    replay_sea_ice = np.asarray(world.fields["cryosphere.sea_ice"], dtype=np.float64)

    earth = {
        "global_mean_surface_temperature_C": _weighted_mean(
            earth_surface, area, np.isfinite(earth_surface)),
        "land_mean_temperature_C": _weighted_mean(earth_land_temp, area, land),
        "ocean_mean_sst_C": _weighted_mean(earth_sst, area, ocean),
        "land_precip_mean_mm_yr": _weighted_mean(earth_precip, area, land),
        "land_precip_p50_mm_yr": _percentile(earth_precip, 50, land),
        "land_precip_p90_mm_yr": _percentile(earth_precip, 90, land),
        "ocean_current_speed_p90_m_s": _percentile(earth_current_speed, 90, ocean),
        "sea_ice_ocean_p95": _percentile(earth_sea_ice, 95, ocean),
    }
    replay = {
        "global_mean_surface_temperature_C": _weighted_mean(
            replay_surface, area, np.isfinite(replay_surface)),
        "land_mean_temperature_C": _weighted_mean(replay_surface, area, land),
        "ocean_mean_sst_C": _weighted_mean(replay_surface, area, ocean),
        "land_precip_mean_mm_yr": _weighted_mean(replay_precip, area, land),
        "land_precip_p50_mm_yr": _percentile(replay_precip, 50, land),
        "land_precip_p90_mm_yr": _percentile(replay_precip, 90, land),
        "ocean_current_speed_p90_m_s": _percentile(
            replay_current_speed, 90, ocean & np.isfinite(replay_current_speed)),
        "sea_ice_ocean_p95": _percentile(replay_sea_ice, 95, ocean),
    }

    temp_delta = replay_surface - earth_surface
    earth_lat_jump = _lat_band_jump(earth_surface, lat, area)
    replay_lat_jump = _lat_band_jump(replay_surface, lat, area)
    precip_delta = replay_precip - earth_precip
    residuals = {
        "surface_temperature_mae_C": _weighted_mean(np.abs(temp_delta), area, np.isfinite(temp_delta)),
        "surface_temperature_rmse_C": _rmse(temp_delta, area, np.isfinite(temp_delta)),
        "land_temperature_mae_C": _weighted_mean(np.abs(temp_delta), area, land),
        "ocean_sst_mae_C": _weighted_mean(np.abs(temp_delta), area, ocean),
        "earth_surface_temperature_max_adjacent_lat_band_delta_C": (
            earth_lat_jump["max_adjacent_lat_band_delta_C"]
        ),
        "replay_surface_temperature_max_adjacent_lat_band_delta_C": (
            replay_lat_jump["max_adjacent_lat_band_delta_C"]
        ),
        "surface_temperature_lat_band_jump_delta_C": (
            replay_lat_jump["max_adjacent_lat_band_delta_C"]
            - earth_lat_jump["max_adjacent_lat_band_delta_C"]
        ),
        "earth_surface_temperature_max_adjacent_lat_band_pair_deg": (
            earth_lat_jump["max_adjacent_lat_band_pair_deg"]
        ),
        "replay_surface_temperature_max_adjacent_lat_band_pair_deg": (
            replay_lat_jump["max_adjacent_lat_band_pair_deg"]
        ),
        "land_precip_mae_mm_yr": _weighted_mean(np.abs(precip_delta), area, land),
        "land_precip_median_delta_mm_yr": (
            replay["land_precip_p50_mm_yr"] - earth["land_precip_p50_mm_yr"]
        ),
        "land_precip_p90_delta_mm_yr": (
            replay["land_precip_p90_mm_yr"] - earth["land_precip_p90_mm_yr"]
        ),
        "ocean_current_speed_p90_ratio": (
            replay["ocean_current_speed_p90_m_s"]
            / max(earth["ocean_current_speed_p90_m_s"], 1.0e-12)
        ),
        "sea_ice_ocean_p95_delta": (
            replay["sea_ice_ocean_p95"] - earth["sea_ice_ocean_p95"]
        ),
    }

    earth_wind = np.asarray(reference["earth__seasonal_wind_u10_v10"], dtype=np.float64)
    replay_wind = np.asarray(world.fields["atmosphere.seasonal_wind"], dtype=np.float64)
    earth_wind_speed = np.linalg.norm(earth_wind, axis=2)
    replay_wind_speed = np.linalg.norm(replay_wind, axis=2)
    annual_earth_wind_speed = np.mean(earth_wind_speed, axis=0)
    annual_replay_wind_speed = np.mean(replay_wind_speed, axis=0)

    earth_slp = np.asarray(reference["earth__seasonal_slp_anomaly_hPa"], dtype=np.float64)
    replay_pressure = np.asarray(
        world.fields["atmosphere.seasonal_pressure_proxy"], dtype=np.float64)
    seasonal_pressure_corr = [
        _pattern_corr(earth_slp[s], replay_pressure[s], np.ones(world.grid.n, dtype=bool))
        for s in range(4)
    ]
    earth_seasonal_precip = np.asarray(
        reference["earth__seasonal_precip_mm_yr_equiv"], dtype=np.float64)
    replay_seasonal_precip = np.asarray(
        world.fields["climate.seasonal_precipitation"], dtype=np.float64)
    nh_land = land & (lat > 0.0)
    sh_land = land & (lat < 0.0)
    dynamics = {
        "annual_wind_speed_global_mae_m_s": _weighted_mean(
            np.abs(annual_replay_wind_speed - annual_earth_wind_speed),
            area,
            np.isfinite(annual_earth_wind_speed) & np.isfinite(annual_replay_wind_speed),
        ),
        "annual_wind_speed_earth_p90_m_s": _percentile(
            annual_earth_wind_speed, 90, np.isfinite(annual_earth_wind_speed)),
        "annual_wind_speed_replay_p90_m_s": _percentile(
            annual_replay_wind_speed, 90, np.isfinite(annual_replay_wind_speed)),
        "seasonal_pressure_pattern_corr_DJF_MAM_JJA_SON": seasonal_pressure_corr,
        "NH_land_JJA_minus_DJF_precip_earth_mm_yr_equiv": _weighted_mean(
            earth_seasonal_precip[2] - earth_seasonal_precip[0], area, nh_land),
        "NH_land_JJA_minus_DJF_precip_replay_mm_yr_equiv": _weighted_mean(
            replay_seasonal_precip[2] - replay_seasonal_precip[0], area, nh_land),
        "SH_land_DJF_minus_JJA_precip_earth_mm_yr_equiv": _weighted_mean(
            earth_seasonal_precip[0] - earth_seasonal_precip[2], area, sh_land),
        "SH_land_DJF_minus_JJA_precip_replay_mm_yr_equiv": _weighted_mean(
            replay_seasonal_precip[0] - replay_seasonal_precip[2], area, sh_land),
    }
    return {
        "earth": earth,
        "replay": replay,
        "residuals": residuals,
        "dynamics_residuals": dynamics,
    }


def _copy_tile(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    image_path: Path,
    label: str,
    *,
    x: int,
    y: int,
    tile_w: int,
    tile_h: int,
    label_h: int,
    font,
) -> None:
    draw.text((x + 4, y + 4), label, fill=(0, 0, 0), font=font)
    if not image_path.exists():
        draw.rectangle(
            [x, y + label_h, x + tile_w, y + label_h + tile_h],
            outline=(180, 180, 180),
        )
        draw.text((x + 12, y + label_h + 24), "missing", fill=(160, 0, 0), font=font)
        return
    with Image.open(image_path) as raw:
        image = raw.convert("RGB")
        image.thumbnail((tile_w, tile_h), Image.Resampling.LANCZOS)
        bg = Image.new("RGB", (tile_w, tile_h), "white")
        bg.paste(image, ((tile_w - image.width) // 2, (tile_h - image.height) // 2))
        canvas.paste(bg, (x, y + label_h))


def _render_contact_sheet(reference_npz: Path, outdir: Path) -> Path | None:
    earth_dir = reference_npz.parent
    earth_prefix = reference_npz.stem
    columns = [
        ("temperature", f"{earth_prefix}_temperature.png", "temperature.png"),
        ("precip", f"{earth_prefix}_precip.png", "precip.png"),
        ("biome", f"{earth_prefix}_biomes_from_koppen_proxy.png", "biomes.png"),
        ("current", f"{earth_prefix}_current_speed.png", "currents.png"),
    ]
    rows = [
        ("Earth reference", [earth_dir / earth_name for _, earth_name, _ in columns]),
        ("Aevum replay", [outdir / replay_name for _, _, replay_name in columns]),
    ]
    tile_w = 300
    tile_h = 170
    label_h = 22
    left_w = 150
    pad = 10
    font = ImageFont.load_default()
    canvas = Image.new(
        "RGB",
        (
            left_w + len(columns) * (tile_w + pad) + pad,
            len(rows) * (tile_h + label_h + pad) + pad,
        ),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    for row_idx, (row_label, paths) in enumerate(rows):
        y = pad + row_idx * (tile_h + label_h + pad)
        draw.text((pad, y + label_h + tile_h // 2 - 8), row_label, fill=(0, 0, 0), font=font)
        for col_idx, ((col_label, _, _), path) in enumerate(zip(columns, paths)):
            x = left_w + col_idx * (tile_w + pad)
            _copy_tile(
                canvas,
                draw,
                path,
                col_label if row_idx == 0 else "",
                x=x,
                y=y,
                tile_w=tile_w,
                tile_h=tile_h,
                label_h=label_h,
                font=font,
            )
    out_path = outdir / "real_earth_replay_contact_sheet.png"
    canvas.save(out_path)
    return out_path


def run_real_earth_climate_replay(
    config: RealEarthClimateReplayConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    world, grid_diagnostics = _world_from_earth_reference(
        Path(config.earth_reference_npz),
        preset=config.preset,
        seed=int(config.seed),
    )
    reference = _earth_reference_fields(Path(config.earth_reference_npz))
    stellar_diag = _restore_stellar_forcing(world, int(config.seed))
    climate_diag = _run_replayed_climate(world, int(config.seed))
    if config.render_assets:
        render.render_world(world, outdir)

    arrays_path = _write_arrays(world, outdir)
    terminal_summary = _summarize_world(SimpleNamespace(world=world), climate_diag, outdir)
    terminal_summary["arrays"] = str(arrays_path)
    terminal_summary["source_reference_npz"] = str(config.earth_reference_npz)
    terminal_summary["stellar_step_diagnostics"] = stellar_diag

    validation_result = validation.check_climate_diagnostics(
        SimpleNamespace(world=world))
    validation_payload = {
        "name": validation_result.name,
        "passed": bool(validation_result.passed),
        "hard_failures": validation_result.detail.get("hard_failures", []),
        "warnings": validation_result.detail.get("warnings", []),
    }
    metrics = _earth_and_replay_metrics(world, reference)
    contact_sheet = None
    if config.render_assets:
        contact_sheet = _render_contact_sheet(Path(config.earth_reference_npz), outdir)

    summary = {
        "schema": SCHEMA,
        "reference_npz": str(config.earth_reference_npz),
        "preset": config.preset,
        "seed": int(config.seed),
        "arrays": str(arrays_path),
        "assets_dir": str(outdir),
        **grid_diagnostics,
        **metrics,
        "validation": validation_payload,
        "climate_step_diagnostics": climate_diag,
        "stellar_step_diagnostics": stellar_diag,
        "terminal_summary": terminal_summary,
        "contact_sheet": None if contact_sheet is None else str(contact_sheet),
    }
    (outdir / "summary.json").write_text(
        json.dumps(terminal_summary, indent=2, default=_json_default)
    )
    (outdir / "real_earth_climate_replay_summary.json").write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    return summary
