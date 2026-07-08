"""Earth-only R2 wind replay comparison.

This diagnostic compares one real-Earth wind subgraph at a time: NOAA/NCEP
seasonal 10 m wind from the reference package against the current Aevum replay
run on the same real-Earth grid.  It is intentionally not a generated-world
guardrail.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.render import to_raster


SCHEMA = "aevum.real_earth_wind_replay.v1"
SEASON_NAMES = ("DJF", "MAM", "JJA", "SON")


@dataclass(frozen=True)
class RealEarthWindReplayConfig:
    earth_reference_npz: Path
    replay_arrays_npz: Path
    outdir: Path


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


def _xyz_from_lat_lon(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat_r = np.radians(np.asarray(lat, dtype=np.float64))
    lon_r = np.radians(np.asarray(lon, dtype=np.float64))
    return np.column_stack([
        np.cos(lat_r) * np.cos(lon_r),
        np.cos(lat_r) * np.sin(lon_r),
        np.sin(lat_r),
    ])


def _tangent_basis(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xyz = _xyz_from_lat_lon(lat, lon)
    z_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    east = np.cross(z_axis, xyz)
    norm = np.linalg.norm(east, axis=1, keepdims=True)
    east = np.where(norm > 1.0e-9, east / np.maximum(norm, 1.0e-9),
                    np.array([1.0, 0.0, 0.0], dtype=np.float64))
    north = np.cross(xyz, east)
    return east, north


def _uv_to_tangent(lat: np.ndarray, lon: np.ndarray, uv: np.ndarray) -> np.ndarray:
    east, north = _tangent_basis(lat, lon)
    uv = np.asarray(uv, dtype=np.float64)
    return uv[:, :, 0, None] * east[None, :, :] + uv[:, :, 1, None] * north[None, :, :]


def _weighted_mean(values: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values) & np.isfinite(weights)
    if not mask.any():
        return float("nan")
    return float(np.average(values[mask], weights=weights[mask]))


def _percentile(values: np.ndarray, q: float, mask: np.ndarray | None = None) -> float:
    values = np.asarray(values, dtype=np.float64)
    if mask is not None:
        values = values[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    return float(np.percentile(values, q)) if values.size else float("nan")


def _corr(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
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


def _standardize_by_season(field: np.ndarray, area: np.ndarray) -> np.ndarray:
    field = np.asarray(field, dtype=np.float64)
    weights = np.asarray(area, dtype=np.float64)
    out = np.zeros_like(field, dtype=np.float64)
    for season in range(field.shape[0]):
        values = field[season]
        valid = np.isfinite(values) & np.isfinite(weights)
        if not valid.any():
            continue
        mean = float(np.average(values[valid], weights=weights[valid]))
        centered = values - mean
        scale = float(np.percentile(np.abs(centered[valid]), 95))
        if scale <= 1.0e-9:
            scale = 1.0
        out[season] = centered / scale
    return out


def _zonal_anomaly(
    field: np.ndarray,
    lat: np.ndarray,
    area: np.ndarray,
    *,
    bin_width_deg: float = 5.0,
) -> np.ndarray:
    field = np.asarray(field, dtype=np.float64)
    lat = np.asarray(lat, dtype=np.float64)
    weights = np.asarray(area, dtype=np.float64)
    out = np.zeros_like(field, dtype=np.float64)
    band_id = np.floor((lat + 90.0) / bin_width_deg).astype(int)
    for season in range(field.shape[0]):
        values = field[season]
        for band in np.unique(band_id):
            mask = (band_id == int(band)) & np.isfinite(values) & np.isfinite(weights)
            if not mask.any():
                continue
            mean = float(np.average(values[mask], weights=weights[mask]))
            out[season, mask] = values[mask] - mean
    return out


def _wind_metrics(
    *,
    lat: np.ndarray,
    lon: np.ndarray,
    area: np.ndarray,
    land: np.ndarray,
    earth_wind: np.ndarray,
    replay_wind: np.ndarray,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    earth_speed = np.linalg.norm(earth_wind, axis=2)
    replay_speed = np.linalg.norm(replay_wind, axis=2)
    speed_delta = replay_speed - earth_speed
    vector_delta = replay_wind - earth_wind
    vector_error = np.linalg.norm(vector_delta, axis=2)
    dot = np.einsum("snk,snk->sn", earth_wind, replay_wind)
    denom = np.maximum(earth_speed * replay_speed, 1.0e-9)
    direction_cosine = np.clip(dot / denom, -1.0, 1.0)
    valid = np.isfinite(earth_speed) & np.isfinite(replay_speed)
    land4 = np.broadcast_to(land, earth_speed.shape)
    ocean4 = ~land4
    east, _ = _tangent_basis(lat, lon)
    earth_eastward = np.einsum("snk,nk->sn", earth_wind, east)
    replay_eastward = np.einsum("snk,nk->sn", replay_wind, east)
    earth_speed_zonal_anomaly = _zonal_anomaly(earth_speed, lat, area)
    replay_speed_zonal_anomaly = _zonal_anomaly(replay_speed, lat, area)
    earth_eastward_zonal_anomaly = _zonal_anomaly(earth_eastward, lat, area)
    replay_eastward_zonal_anomaly = _zonal_anomaly(replay_eastward, lat, area)
    tropics4 = np.broadcast_to(np.abs(lat) < 23.5, earth_speed.shape)
    midlat4 = np.broadcast_to((np.abs(lat) >= 30.0) & (np.abs(lat) <= 60.0),
                              earth_speed.shape)

    summary = {
        "seasonal_speed_earth_p50_m_s": _percentile(earth_speed, 50, valid),
        "seasonal_speed_earth_p90_m_s": _percentile(earth_speed, 90, valid),
        "seasonal_speed_replay_p50_m_s": _percentile(replay_speed, 50, valid),
        "seasonal_speed_replay_p90_m_s": _percentile(replay_speed, 90, valid),
        "seasonal_speed_mae_m_s": _weighted_mean(
            np.abs(speed_delta), np.broadcast_to(area, earth_speed.shape), valid),
        "seasonal_vector_rmse_m_s": float(np.sqrt(_weighted_mean(
            vector_error * vector_error, np.broadcast_to(area, earth_speed.shape), valid))),
        "direction_cosine_p50": _percentile(direction_cosine, 50, valid),
        "direction_cosine_p10": _percentile(direction_cosine, 10, valid),
        "speed_pattern_corr_all": _corr(earth_speed, replay_speed, valid),
        "speed_pattern_corr_land": _corr(earth_speed, replay_speed, valid & land4),
        "speed_pattern_corr_ocean": _corr(earth_speed, replay_speed, valid & ocean4),
        "speed_zonal_anomaly_corr_all": _corr(
            earth_speed_zonal_anomaly, replay_speed_zonal_anomaly, valid),
        "speed_zonal_anomaly_corr_land": _corr(
            earth_speed_zonal_anomaly, replay_speed_zonal_anomaly, valid & land4),
        "speed_zonal_anomaly_corr_ocean": _corr(
            earth_speed_zonal_anomaly, replay_speed_zonal_anomaly, valid & ocean4),
        "eastward_zonal_anomaly_corr_all": _corr(
            earth_eastward_zonal_anomaly, replay_eastward_zonal_anomaly, valid),
        "eastward_zonal_anomaly_corr_land": _corr(
            earth_eastward_zonal_anomaly, replay_eastward_zonal_anomaly, valid & land4),
        "eastward_zonal_anomaly_corr_ocean": _corr(
            earth_eastward_zonal_anomaly, replay_eastward_zonal_anomaly, valid & ocean4),
        "tropical_speed_mae_m_s": _weighted_mean(
            np.abs(speed_delta), np.broadcast_to(area, earth_speed.shape), valid & tropics4),
        "midlatitude_speed_mae_m_s": _weighted_mean(
            np.abs(speed_delta), np.broadcast_to(area, earth_speed.shape), valid & midlat4),
    }

    rows: list[dict[str, Any]] = []
    weights4 = np.broadcast_to(area, earth_speed.shape)
    for season, name in enumerate(SEASON_NAMES):
        mask = valid[season]
        rows.append({
            "season": name,
            "earth_speed_p50_m_s": _percentile(earth_speed[season], 50, mask),
            "earth_speed_p90_m_s": _percentile(earth_speed[season], 90, mask),
            "replay_speed_p50_m_s": _percentile(replay_speed[season], 50, mask),
            "replay_speed_p90_m_s": _percentile(replay_speed[season], 90, mask),
            "speed_mae_m_s": _weighted_mean(
                np.abs(speed_delta[season]), weights4[season], mask),
            "vector_rmse_m_s": float(np.sqrt(_weighted_mean(
                vector_error[season] * vector_error[season], weights4[season], mask))),
            "direction_cosine_p50": _percentile(direction_cosine[season], 50, mask),
            "speed_pattern_corr": _corr(earth_speed[season], replay_speed[season], mask),
            "speed_zonal_anomaly_corr": _corr(
                earth_speed_zonal_anomaly[season],
                replay_speed_zonal_anomaly[season],
                mask,
            ),
            "eastward_zonal_anomaly_corr": _corr(
                earth_eastward_zonal_anomaly[season],
                replay_eastward_zonal_anomaly[season],
                mask,
            ),
        })
    return summary, rows


def _render_seasonal_panels(
    grid: SphereGrid,
    field: np.ndarray,
    out_path: Path,
    *,
    title: str,
    cmap: str,
    vmin: float,
    vmax: float,
) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(10, 5.6), constrained_layout=True)
    im = None
    for ax, name, values in zip(axes.ravel(), SEASON_NAMES, field):
        raster = to_raster(grid, values)
        im = ax.imshow(raster, cmap=cmap, vmin=vmin, vmax=vmax,
                       extent=[-180, 180, -90, 90])
        ax.set_title(f"{name} {title}")
    if im is not None:
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.72)
    fig.savefig(out_path, dpi=115, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _finite_abs_span(values: np.ndarray, *, q: float = 98.0, minimum: float = 1.0) -> float:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float(minimum)
    return max(float(np.percentile(np.abs(arr), q)), float(minimum))


def _finite_span(
    values: np.ndarray,
    *,
    q_low: float = 2.0,
    q_high: float = 98.0,
    fallback: tuple[float, float] = (0.0, 1.0),
) -> tuple[float, float]:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return fallback
    lo = float(np.percentile(arr, q_low))
    hi = float(np.percentile(arr, q_high))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return fallback
    return lo, hi


def _render_static_panel(
    grid: SphereGrid,
    field: np.ndarray,
    out_path: Path,
    *,
    title: str,
    cmap: str,
    vmin: float,
    vmax: float,
) -> Path:
    fig, ax = plt.subplots(1, 1, figsize=(8.8, 4.6), constrained_layout=True)
    raster = to_raster(grid, field)
    im = ax.imshow(raster, cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=[-180, 180, -90, 90])
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.72)
    fig.savefig(out_path, dpi=115, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _render_contact_sheet(paths: list[tuple[str, Path]], out_path: Path) -> Path:
    tile_w = 420
    tile_h = 240
    label_h = 24
    pad = 12
    cols = 2
    rows = max(1, int(np.ceil(len(paths) / cols)))
    font = ImageFont.load_default()
    canvas = Image.new(
        "RGB",
        (tile_w * cols + pad * (cols + 1),
         (tile_h + label_h + pad) * rows + pad),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    for idx, (label, path) in enumerate(paths):
        row = idx // cols
        col = idx % cols
        x = pad + col * (tile_w + pad)
        y = pad + row * (tile_h + label_h + pad)
        draw.text((x + 4, y + 4), label, fill=(0, 0, 0), font=font)
        with Image.open(path) as raw:
            image = raw.convert("RGB")
            image.thumbnail((tile_w, tile_h), Image.Resampling.LANCZOS)
            bg = Image.new("RGB", (tile_w, tile_h), "white")
            bg.paste(image, ((tile_w - image.width) // 2, (tile_h - image.height) // 2))
            canvas.paste(bg, (x, y + label_h))
    canvas.save(out_path)
    return out_path


def _render_replay_support_assets(
    grid: SphereGrid,
    replay_fields: dict[str, np.ndarray],
    outdir: Path,
) -> dict[str, str]:
    assets: dict[str, str] = {}
    boundary_items: list[tuple[str, Path]] = []
    if "terrain__elevation_m" in replay_fields:
        values = replay_fields["terrain__elevation_m"]
        vmin, vmax = _finite_span(values, fallback=(-6000.0, 4500.0))
        path = _render_static_panel(
            grid, values, outdir / "replay_m0_elevation_m.png",
            title="M0 terrain elevation (m)", cmap="terrain", vmin=vmin, vmax=vmax)
        assets["replay_m0_elevation_m"] = str(path)
        boundary_items.append(("elevation", path))
    for key, name, filename, cmap in [
        ("climate__coast_distance", "coast distance", "replay_m0_coast_distance.png", "viridis"),
        ("terrain__barrier_index", "barrier index", "replay_m0_barrier_index.png", "magma"),
        ("terrain__wind_gap_index", "wind gap index", "replay_m0_wind_gap_index.png", "viridis"),
        ("ocean__shelf_index", "shelf index", "replay_m0_shelf_index.png", "viridis"),
        ("ocean__strait_index", "strait index", "replay_m0_strait_index.png", "magma"),
        ("climate__surface_heat_capacity_class", "heat-capacity class", "replay_m1_heat_capacity_class.png", "viridis"),
        ("climate__elevation_lapse_cooling", "lapse cooling (K)", "replay_m1_lapse_cooling.png", "magma"),
    ]:
        if key not in replay_fields:
            continue
        values = replay_fields[key]
        vmax = 1.0
        if key == "climate__elevation_lapse_cooling":
            vmax = max(float(np.nanpercentile(values, 98)), 1.0)
        path = _render_static_panel(
            grid, values, outdir / filename,
            title=f"M0/M1 {name}",
            cmap=cmap,
            vmin=0.0,
            vmax=vmax,
        )
        assets[filename.removesuffix(".png")] = str(path)
        boundary_items.append((name, path))
    if "ocean__basin_id" in replay_fields:
        basin = replay_fields["ocean__basin_id"]
        vmin, vmax = _finite_span(basin, fallback=(-1.0, 8.0))
        path = _render_static_panel(
            grid, basin, outdir / "replay_m0_ocean_basin_id.png",
            title="M0 ocean basin id", cmap="tab20", vmin=vmin, vmax=vmax)
        assets["replay_m0_ocean_basin_id"] = str(path)
        boundary_items.append(("ocean basin id", path))
    if boundary_items:
        boundary_contact = _render_contact_sheet(
            boundary_items,
            outdir / "replay_m0_boundary_support_contact_sheet.png",
        )
        assets["replay_m0_boundary_support_contact_sheet"] = str(boundary_contact)

    seasonal_specs = [
        (
            "climate__seasonal_insolation_anomaly",
            "M1 seasonal insolation anomaly (W/m2)",
            "replay_m1_seasonal_insolation_anomaly.png",
            "coolwarm",
            "abs",
        ),
        (
            "climate__land_thermal_anomaly",
            "M1 land thermal anomaly (K)",
            "replay_m1_land_thermal_anomaly.png",
            "coolwarm",
            "abs",
        ),
        (
            "climate__ocean_mixed_layer_thermal_anomaly",
            "M1 ocean mixed-layer anomaly (K)",
            "replay_m1_ocean_mixed_layer_thermal_anomaly.png",
            "coolwarm",
            "abs",
        ),
        (
            "climate__land_sea_thermal_contrast",
            "M1 land-sea thermal contrast (K)",
            "replay_m1_land_sea_thermal_contrast.png",
            "coolwarm",
            "abs",
        ),
        (
            "climate__same_latitude_sst_anomaly",
            "M1 same-latitude SST anomaly (K)",
            "replay_m1_same_latitude_sst_anomaly.png",
            "coolwarm",
            "abs",
        ),
        (
            "climate__snow_ice_albedo_support",
            "M1 snow/ice albedo support",
            "replay_m1_snow_ice_albedo_support.png",
            "Blues",
            "unit",
        ),
        (
            "climate__sst_gradient_support",
            "M1 SST-gradient support",
            "replay_m1_sst_gradient_support.png",
            "magma",
            "unit15",
        ),
    ]
    energy_items: list[tuple[str, Path]] = []
    for key, title, filename, cmap, scale_kind in seasonal_specs:
        if key not in replay_fields:
            continue
        field = np.asarray(replay_fields[key], dtype=np.float64)
        if field.shape != (4, grid.n):
            continue
        if scale_kind == "abs":
            span = _finite_abs_span(field)
            vmin, vmax = -span, span
        elif scale_kind == "unit15":
            vmin, vmax = 0.0, 1.5
        else:
            vmin, vmax = 0.0, 1.0
        path = _render_seasonal_panels(
            grid, field, outdir / filename,
            title=title, cmap=cmap, vmin=vmin, vmax=vmax)
        assets[filename.removesuffix(".png")] = str(path)
        energy_items.append((title.replace("M1 ", ""), path))
    if energy_items:
        energy_contact = _render_contact_sheet(
            energy_items,
            outdir / "replay_m1_energy_support_contact_sheet.png",
        )
        assets["replay_m1_energy_support_contact_sheet"] = str(energy_contact)
    return assets


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> Path:
    keys = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def run_real_earth_wind_replay(config: RealEarthWindReplayConfig) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with np.load(config.earth_reference_npz, allow_pickle=False) as ref:
        lat = np.asarray(ref["lat"], dtype=np.float64)
        lon = np.asarray(ref["lon"], dtype=np.float64)
        area = np.asarray(ref["cell_area"], dtype=np.float64)
        land = np.asarray(ref["earth__elevation_m"], dtype=np.float64) >= 0.0
        earth_wind = _uv_to_tangent(
            lat, lon, np.asarray(ref["earth__seasonal_wind_u10_v10"], dtype=np.float64))
        earth_slp_anomaly = (
            np.asarray(ref["earth__seasonal_slp_anomaly_hPa"], dtype=np.float64)
            if "earth__seasonal_slp_anomaly_hPa" in ref.files
            else None
        )
    with np.load(config.replay_arrays_npz, allow_pickle=False) as replay:
        replay_wind = np.asarray(replay["atmosphere__seasonal_wind"], dtype=np.float64)
        replay_pressure = (
            np.asarray(replay["atmosphere__land_sea_pressure_proxy"], dtype=np.float64)
            if "atmosphere__land_sea_pressure_proxy" in replay.files
            else None
        )
        replay_pressure_center_support = (
            np.asarray(replay["atmosphere__pressure_center_support"], dtype=np.float64)
            if "atmosphere__pressure_center_support" in replay.files
            else None
        )
        replay_stationary_wave_support = (
            np.asarray(
                replay["atmosphere__stationary_wave_pressure_support"],
                dtype=np.float64,
            )
            if "atmosphere__stationary_wave_pressure_support" in replay.files
            else None
        )
        replay_pressure_genesis_source = (
            np.asarray(replay["atmosphere__pressure_genesis_source"], dtype=np.float64)
            if "atmosphere__pressure_genesis_source" in replay.files
            else None
        )
        replay_pressure_genesis_wave_transfer = (
            np.asarray(
                replay["atmosphere__pressure_genesis_wave_transfer"],
                dtype=np.float64,
            )
            if "atmosphere__pressure_genesis_wave_transfer" in replay.files
            else None
        )
        replay_ocean_low_source_support = (
            np.asarray(
                replay["atmosphere__ocean_pressure_low_source_support"],
                dtype=np.float64,
            )
            if "atmosphere__ocean_pressure_low_source_support" in replay.files
            else None
        )
        replay_ocean_high_source_support = (
            np.asarray(
                replay["atmosphere__ocean_pressure_high_source_support"],
                dtype=np.float64,
            )
            if "atmosphere__ocean_pressure_high_source_support" in replay.files
            else None
        )
        replay_land_source_support = (
            np.asarray(
                replay["atmosphere__land_pressure_source_support"],
                dtype=np.float64,
            )
            if "atmosphere__land_pressure_source_support" in replay.files
            else None
        )
        replay_terrain_source_support = (
            np.asarray(
                replay["atmosphere__terrain_pressure_wave_source_support"],
                dtype=np.float64,
            )
            if "atmosphere__terrain_pressure_wave_source_support" in replay.files
            else None
        )
        support_field_keys = [
            "terrain__elevation_m",
            "climate__coast_distance",
            "terrain__barrier_index",
            "terrain__wind_gap_index",
            "ocean__basin_id",
            "ocean__shelf_index",
            "ocean__strait_index",
            "climate__seasonal_insolation_anomaly",
            "climate__surface_heat_capacity_class",
            "climate__land_thermal_anomaly",
            "climate__ocean_mixed_layer_thermal_anomaly",
            "climate__elevation_lapse_cooling",
            "climate__snow_ice_albedo_support",
            "climate__sst_gradient_support",
            "climate__same_latitude_sst_anomaly",
            "climate__land_sea_thermal_contrast",
            "atmosphere__pressure_genesis_source",
            "atmosphere__pressure_genesis_wave_transfer",
            "atmosphere__ocean_pressure_low_source_support",
            "atmosphere__ocean_pressure_high_source_support",
            "atmosphere__land_pressure_source_support",
            "atmosphere__terrain_pressure_wave_source_support",
        ]
        replay_support_fields = {
            key: np.asarray(replay[key], dtype=np.float64)
            for key in support_field_keys
            if key in replay.files
        }

    if earth_wind.shape != replay_wind.shape:
        raise ValueError(
            f"wind shape mismatch: earth {earth_wind.shape}, replay {replay_wind.shape}"
        )
    grid = SphereGrid.fibonacci(int(lat.size), CONSTANTS.EARTH_RADIUS)
    if float(np.max(np.abs(grid.lat - lat))) > 1.0e-8:
        raise ValueError("reference grid does not match SphereGrid.fibonacci latitude order")

    earth_speed = np.linalg.norm(earth_wind, axis=2)
    replay_speed = np.linalg.norm(replay_wind, axis=2)
    speed_delta = replay_speed - earth_speed
    vector_error = np.linalg.norm(replay_wind - earth_wind, axis=2)
    east, north = _tangent_basis(lat, np.asarray(lon, dtype=np.float64))
    earth_eastward = np.einsum("snk,nk->sn", earth_wind, east)
    replay_eastward = np.einsum("snk,nk->sn", replay_wind, east)
    earth_northward = np.einsum("snk,nk->sn", earth_wind, north)
    replay_northward = np.einsum("snk,nk->sn", replay_wind, north)
    eastward_delta = replay_eastward - earth_eastward
    northward_delta = replay_northward - earth_northward
    earth_speed_zonal_anomaly = _zonal_anomaly(earth_speed, lat, area)
    replay_speed_zonal_anomaly = _zonal_anomaly(replay_speed, lat, area)
    speed_zonal_delta = replay_speed_zonal_anomaly - earth_speed_zonal_anomaly
    earth_eastward_zonal_anomaly = _zonal_anomaly(earth_eastward, lat, area)
    replay_eastward_zonal_anomaly = _zonal_anomaly(replay_eastward, lat, area)
    eastward_zonal_delta = (
        replay_eastward_zonal_anomaly - earth_eastward_zonal_anomaly)

    vmax = max(
        float(np.nanpercentile(earth_speed, 98)),
        float(np.nanpercentile(replay_speed, 98)),
        1.0,
    )
    delta_span = max(float(np.nanpercentile(np.abs(speed_delta), 98)), 1.0)
    error_vmax = max(float(np.nanpercentile(vector_error, 98)), 1.0)
    earth_east_span = max(float(np.nanpercentile(np.abs(earth_eastward), 98)), 1.0)
    replay_east_span = max(float(np.nanpercentile(np.abs(replay_eastward), 98)), 1.0)
    east_span = max(earth_east_span, replay_east_span)
    east_delta_span = max(float(np.nanpercentile(np.abs(eastward_delta), 98)), 1.0)
    north_delta_span = max(float(np.nanpercentile(np.abs(northward_delta), 98)), 1.0)
    speed_zonal_span = max(
        float(np.nanpercentile(np.abs(earth_speed_zonal_anomaly), 98)),
        float(np.nanpercentile(np.abs(replay_speed_zonal_anomaly), 98)),
        1.0,
    )
    speed_zonal_delta_span = max(
        float(np.nanpercentile(np.abs(speed_zonal_delta), 98)), 1.0)
    east_zonal_span = max(
        float(np.nanpercentile(np.abs(earth_eastward_zonal_anomaly), 98)),
        float(np.nanpercentile(np.abs(replay_eastward_zonal_anomaly), 98)),
        1.0,
    )
    east_zonal_delta_span = max(
        float(np.nanpercentile(np.abs(eastward_zonal_delta), 98)), 1.0)
    earth_png = _render_seasonal_panels(
        grid, earth_speed, outdir / "earth_wind_speed_seasons.png",
        title="Earth 10 m wind speed (m/s)", cmap="viridis", vmin=0.0, vmax=vmax)
    replay_png = _render_seasonal_panels(
        grid, replay_speed, outdir / "replay_wind_speed_seasons.png",
        title="Aevum replay wind speed (m/s)", cmap="viridis", vmin=0.0, vmax=vmax)
    delta_png = _render_seasonal_panels(
        grid, speed_delta, outdir / "wind_speed_delta_seasons.png",
        title="replay minus Earth wind speed (m/s)", cmap="coolwarm",
        vmin=-delta_span, vmax=delta_span)
    error_png = _render_seasonal_panels(
        grid, vector_error, outdir / "wind_vector_error_seasons.png",
        title="wind vector error (m/s)", cmap="magma", vmin=0.0, vmax=error_vmax)
    earth_east_png = _render_seasonal_panels(
        grid, earth_eastward, outdir / "earth_eastward_wind_seasons.png",
        title="Earth eastward wind (m/s)", cmap="coolwarm",
        vmin=-east_span, vmax=east_span)
    replay_east_png = _render_seasonal_panels(
        grid, replay_eastward, outdir / "replay_eastward_wind_seasons.png",
        title="Aevum replay eastward wind (m/s)", cmap="coolwarm",
        vmin=-east_span, vmax=east_span)
    east_delta_png = _render_seasonal_panels(
        grid, eastward_delta, outdir / "eastward_wind_delta_seasons.png",
        title="replay minus Earth eastward wind (m/s)", cmap="coolwarm",
        vmin=-east_delta_span, vmax=east_delta_span)
    north_delta_png = _render_seasonal_panels(
        grid, northward_delta, outdir / "northward_wind_delta_seasons.png",
        title="replay minus Earth northward wind (m/s)", cmap="coolwarm",
        vmin=-north_delta_span, vmax=north_delta_span)
    earth_speed_zonal_png = _render_seasonal_panels(
        grid, earth_speed_zonal_anomaly,
        outdir / "earth_wind_speed_zonal_anomaly_seasons.png",
        title="Earth wind-speed zonal anomaly (m/s)", cmap="coolwarm",
        vmin=-speed_zonal_span, vmax=speed_zonal_span)
    replay_speed_zonal_png = _render_seasonal_panels(
        grid, replay_speed_zonal_anomaly,
        outdir / "replay_wind_speed_zonal_anomaly_seasons.png",
        title="Aevum wind-speed zonal anomaly (m/s)", cmap="coolwarm",
        vmin=-speed_zonal_span, vmax=speed_zonal_span)
    speed_zonal_delta_png = _render_seasonal_panels(
        grid, speed_zonal_delta,
        outdir / "wind_speed_zonal_anomaly_delta_seasons.png",
        title="replay minus Earth wind-speed zonal anomaly (m/s)", cmap="coolwarm",
        vmin=-speed_zonal_delta_span, vmax=speed_zonal_delta_span)
    earth_east_zonal_png = _render_seasonal_panels(
        grid, earth_eastward_zonal_anomaly,
        outdir / "earth_eastward_zonal_anomaly_seasons.png",
        title="Earth eastward-wind zonal anomaly (m/s)", cmap="coolwarm",
        vmin=-east_zonal_span, vmax=east_zonal_span)
    replay_east_zonal_png = _render_seasonal_panels(
        grid, replay_eastward_zonal_anomaly,
        outdir / "replay_eastward_zonal_anomaly_seasons.png",
        title="Aevum eastward-wind zonal anomaly (m/s)", cmap="coolwarm",
        vmin=-east_zonal_span, vmax=east_zonal_span)
    east_zonal_delta_png = _render_seasonal_panels(
        grid, eastward_zonal_delta,
        outdir / "eastward_zonal_anomaly_delta_seasons.png",
        title="replay minus Earth eastward-wind zonal anomaly (m/s)", cmap="coolwarm",
        vmin=-east_zonal_delta_span, vmax=east_zonal_delta_span)
    pressure_assets: dict[str, str] = {}
    pressure_metrics: dict[str, float] = {}
    if (
        earth_slp_anomaly is not None
        and replay_pressure is not None
        and earth_slp_anomaly.shape == replay_pressure.shape == (4, lat.size)
    ):
        earth_pressure_std = _standardize_by_season(earth_slp_anomaly, area)
        replay_pressure_std = _standardize_by_season(replay_pressure, area)
        pressure_delta = replay_pressure_std - earth_pressure_std
        earth_pressure_zonal = _zonal_anomaly(earth_pressure_std, lat, area)
        replay_pressure_zonal = _zonal_anomaly(replay_pressure_std, lat, area)
        pressure_zonal_delta = replay_pressure_zonal - earth_pressure_zonal
        pressure_span = max(
            float(np.nanpercentile(np.abs(earth_pressure_std), 98)),
            float(np.nanpercentile(np.abs(replay_pressure_std), 98)),
            1.0,
        )
        pressure_zonal_span = max(
            float(np.nanpercentile(np.abs(earth_pressure_zonal), 98)),
            float(np.nanpercentile(np.abs(replay_pressure_zonal), 98)),
            1.0,
        )
        pressure_zonal_delta_span = max(
            float(np.nanpercentile(np.abs(pressure_zonal_delta), 98)), 1.0)
        pressure_delta_span = max(
            float(np.nanpercentile(np.abs(pressure_delta), 98)), 1.0)
        earth_pressure_png = _render_seasonal_panels(
            grid, earth_pressure_std,
            outdir / "earth_slp_standardized_anomaly_seasons.png",
            title="Earth standardized SLP anomaly", cmap="coolwarm",
            vmin=-pressure_span, vmax=pressure_span)
        replay_pressure_png = _render_seasonal_panels(
            grid, replay_pressure_std,
            outdir / "replay_pressure_proxy_standardized_anomaly_seasons.png",
            title="Aevum standardized pressure proxy", cmap="coolwarm",
            vmin=-pressure_span, vmax=pressure_span)
        pressure_delta_png = _render_seasonal_panels(
            grid, pressure_delta,
            outdir / "pressure_standardized_delta_seasons.png",
            title="replay minus Earth standardized pressure", cmap="coolwarm",
            vmin=-pressure_delta_span, vmax=pressure_delta_span)
        earth_pressure_zonal_png = _render_seasonal_panels(
            grid, earth_pressure_zonal,
            outdir / "earth_pressure_zonal_anomaly_seasons.png",
            title="Earth pressure zonal anomaly", cmap="coolwarm",
            vmin=-pressure_zonal_span, vmax=pressure_zonal_span)
        replay_pressure_zonal_png = _render_seasonal_panels(
            grid, replay_pressure_zonal,
            outdir / "replay_pressure_zonal_anomaly_seasons.png",
            title="Aevum pressure zonal anomaly", cmap="coolwarm",
            vmin=-pressure_zonal_span, vmax=pressure_zonal_span)
        pressure_zonal_delta_png = _render_seasonal_panels(
            grid, pressure_zonal_delta,
            outdir / "pressure_zonal_anomaly_delta_seasons.png",
            title="replay minus Earth pressure zonal anomaly", cmap="coolwarm",
            vmin=-pressure_zonal_delta_span, vmax=pressure_zonal_delta_span)
        if (
            replay_pressure_center_support is not None
            and replay_pressure_center_support.shape == replay_pressure.shape
        ):
            center_support_png = _render_seasonal_panels(
                grid,
                np.clip(replay_pressure_center_support, 0.0, None),
                outdir / "replay_pressure_center_support_seasons.png",
                title="Aevum pressure-center support",
                cmap="viridis",
                vmin=0.0,
                vmax=max(
                    float(np.nanpercentile(replay_pressure_center_support, 98)),
                    1.0,
                ),
            )
            pressure_assets["replay_pressure_center_support_seasons"] = str(
                center_support_png)
        else:
            center_support_png = None
        if (
            replay_stationary_wave_support is not None
            and replay_stationary_wave_support.shape == replay_pressure.shape
        ):
            stationary_support_png = _render_seasonal_panels(
                grid,
                np.clip(replay_stationary_wave_support, 0.0, None),
                outdir / "replay_stationary_wave_pressure_support_seasons.png",
                title="Aevum stationary-wave pressure support",
                cmap="magma",
                vmin=0.0,
                vmax=max(
                    float(np.nanpercentile(replay_stationary_wave_support, 98)),
                    1.0,
                ),
            )
            pressure_assets["replay_stationary_wave_pressure_support_seasons"] = (
                str(stationary_support_png))
        else:
            stationary_support_png = None
        if (
            replay_pressure_genesis_source is not None
            and replay_pressure_genesis_source.shape == replay_pressure.shape
        ):
            source_span = max(
                float(np.nanpercentile(np.abs(replay_pressure_genesis_source), 98)),
                0.20,
            )
            genesis_source_png = _render_seasonal_panels(
                grid,
                replay_pressure_genesis_source,
                outdir / "replay_pressure_genesis_source_seasons.png",
                title="Aevum M2 pressure-genesis source",
                cmap="coolwarm",
                vmin=-source_span,
                vmax=source_span,
            )
            pressure_assets["replay_pressure_genesis_source_seasons"] = str(
                genesis_source_png)
        else:
            genesis_source_png = None
        if (
            replay_pressure_genesis_wave_transfer is not None
            and replay_pressure_genesis_wave_transfer.shape == replay_pressure.shape
        ):
            transfer_span = max(
                float(np.nanpercentile(
                    np.abs(replay_pressure_genesis_wave_transfer), 98)),
                0.10,
            )
            genesis_transfer_png = _render_seasonal_panels(
                grid,
                replay_pressure_genesis_wave_transfer,
                outdir / "replay_pressure_genesis_wave_transfer_seasons.png",
                title="Aevum M2 pressure-genesis wave transfer",
                cmap="coolwarm",
                vmin=-transfer_span,
                vmax=transfer_span,
            )
            pressure_assets["replay_pressure_genesis_wave_transfer_seasons"] = str(
                genesis_transfer_png)
        else:
            genesis_transfer_png = None
        if (
            replay_ocean_low_source_support is not None
            and replay_ocean_low_source_support.shape == replay_pressure.shape
        ):
            ocean_low_source_png = _render_seasonal_panels(
                grid,
                np.clip(replay_ocean_low_source_support, 0.0, None),
                outdir / "replay_ocean_pressure_low_source_support_seasons.png",
                title="Aevum ocean low source support",
                cmap="viridis",
                vmin=0.0,
                vmax=max(
                    float(np.nanpercentile(replay_ocean_low_source_support, 98)),
                    1.0,
                ),
            )
            pressure_assets["replay_ocean_pressure_low_source_support_seasons"] = (
                str(ocean_low_source_png))
        else:
            ocean_low_source_png = None
        if (
            replay_ocean_high_source_support is not None
            and replay_ocean_high_source_support.shape == replay_pressure.shape
        ):
            ocean_high_source_png = _render_seasonal_panels(
                grid,
                np.clip(replay_ocean_high_source_support, 0.0, None),
                outdir / "replay_ocean_pressure_high_source_support_seasons.png",
                title="Aevum ocean high source support",
                cmap="viridis",
                vmin=0.0,
                vmax=max(
                    float(np.nanpercentile(replay_ocean_high_source_support, 98)),
                    1.0,
                ),
            )
            pressure_assets["replay_ocean_pressure_high_source_support_seasons"] = (
                str(ocean_high_source_png))
        else:
            ocean_high_source_png = None
        if (
            replay_land_source_support is not None
            and replay_land_source_support.shape == replay_pressure.shape
        ):
            land_source_png = _render_seasonal_panels(
                grid,
                np.clip(replay_land_source_support, 0.0, None),
                outdir / "replay_land_pressure_source_support_seasons.png",
                title="Aevum land pressure source support",
                cmap="magma",
                vmin=0.0,
                vmax=max(float(np.nanpercentile(replay_land_source_support, 98)), 1.0),
            )
            pressure_assets["replay_land_pressure_source_support_seasons"] = str(
                land_source_png)
        else:
            land_source_png = None
        if (
            replay_terrain_source_support is not None
            and replay_terrain_source_support.shape == replay_pressure.shape
        ):
            terrain_source_png = _render_seasonal_panels(
                grid,
                np.clip(replay_terrain_source_support, 0.0, None),
                outdir / "replay_terrain_pressure_wave_source_support_seasons.png",
                title="Aevum terrain pressure-wave source support",
                cmap="magma",
                vmin=0.0,
                vmax=max(
                    float(np.nanpercentile(replay_terrain_source_support, 98)),
                    1.0,
                ),
            )
            pressure_assets["replay_terrain_pressure_wave_source_support_seasons"] = (
                str(terrain_source_png))
        else:
            terrain_source_png = None
        contact_items = [
            ("Earth SLP anomaly", earth_pressure_png),
            ("Aevum pressure proxy", replay_pressure_png),
            ("standardized residual", pressure_delta_png),
            ("zonal residual", pressure_zonal_delta_png),
        ]
        if genesis_source_png is not None:
            contact_items.append(("M2 pressure-genesis source", genesis_source_png))
        if genesis_transfer_png is not None:
            contact_items.append(("M2 source-to-pressure transfer", genesis_transfer_png))
        if ocean_low_source_png is not None:
            contact_items.append(("ocean low source support", ocean_low_source_png))
        if center_support_png is not None:
            contact_items.append(("pressure-center support", center_support_png))
        if stationary_support_png is not None:
            contact_items.append(("stationary-wave support", stationary_support_png))
        pressure_contact = _render_contact_sheet(
            contact_items,
            outdir / "real_earth_pressure_replay_contact_sheet.png",
        )
        valid_pressure = np.isfinite(earth_pressure_std) & np.isfinite(replay_pressure_std)
        land4 = np.broadcast_to(land, earth_pressure_std.shape)
        ocean4 = ~land4
        pressure_weights = np.broadcast_to(area, earth_pressure_std.shape)
        pressure_metrics = {
            "pressure_standardized_mae_all": _weighted_mean(
                np.abs(pressure_delta), pressure_weights, valid_pressure),
            "pressure_standardized_mae_land": _weighted_mean(
                np.abs(pressure_delta), pressure_weights, valid_pressure & land4),
            "pressure_standardized_mae_ocean": _weighted_mean(
                np.abs(pressure_delta), pressure_weights, valid_pressure & ocean4),
            "pressure_standardized_corr_all": _corr(
                earth_pressure_std, replay_pressure_std, valid_pressure),
            "pressure_standardized_corr_land": _corr(
                earth_pressure_std, replay_pressure_std, valid_pressure & land4),
            "pressure_standardized_corr_ocean": _corr(
                earth_pressure_std, replay_pressure_std, valid_pressure & ocean4),
            "pressure_zonal_anomaly_corr_all": _corr(
                earth_pressure_zonal, replay_pressure_zonal, valid_pressure),
            "pressure_zonal_anomaly_corr_land": _corr(
                earth_pressure_zonal, replay_pressure_zonal, valid_pressure & land4),
            "pressure_zonal_anomaly_corr_ocean": _corr(
                earth_pressure_zonal, replay_pressure_zonal, valid_pressure & ocean4),
        }
        for season, name in enumerate(SEASON_NAMES):
            valid = valid_pressure[season]
            pressure_metrics[f"pressure_standardized_mae_{name.lower()}"] = (
                _weighted_mean(np.abs(pressure_delta[season]), area, valid)
            )
            pressure_metrics[f"pressure_zonal_anomaly_corr_{name.lower()}"] = _corr(
                earth_pressure_zonal[season],
                replay_pressure_zonal[season],
                valid,
            )
        pressure_assets = {
            "earth_slp_standardized_anomaly_seasons": str(earth_pressure_png),
            "replay_pressure_proxy_standardized_anomaly_seasons": str(replay_pressure_png),
            "pressure_standardized_delta_seasons": str(pressure_delta_png),
            "earth_pressure_zonal_anomaly_seasons": str(earth_pressure_zonal_png),
            "replay_pressure_zonal_anomaly_seasons": str(replay_pressure_zonal_png),
            "pressure_zonal_anomaly_delta_seasons": str(pressure_zonal_delta_png),
            **pressure_assets,
            "pressure_contact_sheet": str(pressure_contact),
        }
    contact = _render_contact_sheet(
        [
            ("Earth reference", earth_png),
            ("Aevum Earth replay", replay_png),
            ("speed residual", delta_png),
            ("vector error", error_png),
        ],
        outdir / "real_earth_wind_replay_contact_sheet.png",
    )
    metrics, seasonal_rows = _wind_metrics(
        lat=lat,
        lon=lon,
        area=area,
        land=land,
        earth_wind=earth_wind,
        replay_wind=replay_wind,
    )
    metrics.update(pressure_metrics)
    support_assets = _render_replay_support_assets(grid, replay_support_fields, outdir)
    _write_csv(outdir / "real_earth_wind_replay_seasonal_metrics.csv", seasonal_rows)
    summary = {
        "schema": SCHEMA,
        "earth_reference_npz": str(config.earth_reference_npz),
        "replay_arrays_npz": str(config.replay_arrays_npz),
        "outdir": str(outdir),
        "metrics": metrics,
        "seasonal_metrics": seasonal_rows,
        "assets": {
            "earth_wind_speed_seasons": str(earth_png),
            "replay_wind_speed_seasons": str(replay_png),
            "wind_speed_delta_seasons": str(delta_png),
            "wind_vector_error_seasons": str(error_png),
            "earth_eastward_wind_seasons": str(earth_east_png),
            "replay_eastward_wind_seasons": str(replay_east_png),
            "eastward_wind_delta_seasons": str(east_delta_png),
            "northward_wind_delta_seasons": str(north_delta_png),
            "earth_wind_speed_zonal_anomaly_seasons": str(earth_speed_zonal_png),
            "replay_wind_speed_zonal_anomaly_seasons": str(replay_speed_zonal_png),
            "wind_speed_zonal_anomaly_delta_seasons": str(speed_zonal_delta_png),
            "earth_eastward_zonal_anomaly_seasons": str(earth_east_zonal_png),
            "replay_eastward_zonal_anomaly_seasons": str(replay_east_zonal_png),
            "eastward_zonal_anomaly_delta_seasons": str(east_zonal_delta_png),
            **pressure_assets,
            **support_assets,
            "contact_sheet": str(contact),
        },
    }
    (outdir / "real_earth_wind_replay_summary.json").write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    return summary
