"""Compare generated terminal climates against same-grid real-Earth references."""
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
from matplotlib.colors import BoundaryNorm, ListedColormap
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import cKDTree

from aevum import render


SCHEMA = "aevum.earth_climate_comparison.v1"
BIOME_LABELS = {
    0: "ocean",
    1: "ice",
    2: "desert",
    3: "grassland",
    4: "forest",
    5: "tundra",
    6: "tropical",
}
LAND_COVER_BROAD_LABELS = {
    1: "water",
    2: "cropland",
    3: "forest",
    4: "grass_shrub",
    5: "wetland",
    6: "urban",
    7: "bare_sparse",
    8: "snow_ice",
}
SCALAR_COMPARISON_SCALES = {
    "land_fraction": 0.12,
    "global_mean_temperature_C": 6.0,
    "land_mean_temperature_C": 6.0,
    "ocean_mean_temperature_C": 5.0,
    "land_precip_mean_mm_yr": 450.0,
    "land_precip_p50_mm_yr": 350.0,
    "land_precip_p90_mm_yr": 900.0,
    "land_seasonal_temp_amp_p50_C": 7.0,
    "ocean_seasonal_temp_amp_p50_C": 3.0,
    "precip_seasonality_land_p75": 0.9,
    "current_speed_p50_m_s": 0.10,
    "current_speed_p90_m_s": 0.18,
    "biome_desert_area_fraction": 0.15,
    "biome_forest_area_fraction": 0.15,
    "biome_tropical_area_fraction": 0.08,
}


@dataclass(frozen=True)
class EarthClimateComparisonConfig:
    earth_reference_npz: Path
    terminal_summary_json: Path
    outdir: Path
    render_contact_sheet: bool = True


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


def _weighted_fraction(weights: np.ndarray, mask: np.ndarray) -> float:
    weights = np.asarray(weights, dtype=np.float64)
    total = float(np.nansum(weights))
    if total <= 0.0:
        return 0.0
    return float(np.nansum(weights[np.asarray(mask, dtype=bool)]) / total)


def _percentile(values: np.ndarray, q: float, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    if not mask.any():
        return float("nan")
    return float(np.nanpercentile(values[mask], q))


def _seasonal_amplitude(seasonal: np.ndarray) -> np.ndarray:
    seasonal = np.asarray(seasonal, dtype=np.float64)
    return _nanmax_axis0(seasonal) - _nanmin_axis0(seasonal)


def _nanmax_axis0(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    all_nan = ~np.isfinite(values).any(axis=0)
    filled = np.where(np.isfinite(values), values, -np.inf)
    out = np.max(filled, axis=0)
    out[all_nan] = np.nan
    return out


def _nanmin_axis0(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    all_nan = ~np.isfinite(values).any(axis=0)
    filled = np.where(np.isfinite(values), values, np.inf)
    out = np.min(filled, axis=0)
    out[all_nan] = np.nan
    return out


def _precip_seasonality_from_annual_equiv(
    seasonal_annual_equiv: np.ndarray,
    annual_precip: np.ndarray,
) -> np.ndarray:
    seasonal_annual_equiv = np.asarray(seasonal_annual_equiv, dtype=np.float64)
    annual_precip = np.asarray(annual_precip, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = _nanmax_axis0(seasonal_annual_equiv) / np.maximum(annual_precip, 1.0e-9)
    out[~np.isfinite(out)] = np.nan
    return out


def _biome_area_fractions(
    biome: np.ndarray,
    area: np.ndarray,
    *,
    prefix: str = "biome",
) -> dict[str, float]:
    biome = np.asarray(biome, dtype=np.int16)
    area = np.asarray(area, dtype=np.float64)
    total = max(float(np.nansum(area)), 1.0e-12)
    out: dict[str, float] = {}
    for code in range(0, 7):
        label = BIOME_LABELS.get(code, str(code))
        out[f"{prefix}_{label}_area_fraction"] = float(
            np.nansum(area[biome == code]) / total
        )
    return out


def _land_cover_area_fractions(
    land_cover: np.ndarray,
    area: np.ndarray,
    *,
    prefix: str = "land_cover",
) -> dict[str, float]:
    land_cover = np.asarray(land_cover, dtype=np.uint8)
    area = np.asarray(area, dtype=np.float64)
    total = max(float(np.nansum(area)), 1.0e-12)
    valid = land_cover > 0
    out: dict[str, float] = {
        f"{prefix}_valid_area_fraction": float(np.nansum(area[valid]) / total),
    }
    for code, label in LAND_COVER_BROAD_LABELS.items():
        out[f"{prefix}_{label}_area_fraction"] = float(
            np.nansum(area[land_cover == code]) / total
        )
    return out


def _lat_band_means(
    lat: np.ndarray,
    values: np.ndarray,
    weights: np.ndarray,
    *,
    width_deg: float = 10.0,
) -> list[float | None]:
    lat = np.asarray(lat, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    means: list[float | None] = []
    for lo in np.arange(-90.0, 90.0, width_deg):
        hi = lo + width_deg
        mask = (lat >= lo) & (lat < hi) & np.isfinite(values)
        if not mask.any():
            means.append(None)
        else:
            means.append(float(np.average(values[mask], weights=weights[mask])))
    return means


def _max_adjacent_band_delta(means: list[float | None]) -> float:
    vals = [m for m in means if m is not None]
    if len(vals) < 2:
        return float("nan")
    return float(max(abs(vals[i + 1] - vals[i]) for i in range(len(vals) - 1)))


def earth_reference_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        area = np.asarray(z["cell_area"], dtype=np.float64)
        lat = np.asarray(z["lat"], dtype=np.float64)
        land = np.asarray(z["earth__land_mask"], dtype=bool)
        ocean = ~land
        land_temp = np.asarray(z["earth__annual_temperature_C"], dtype=np.float64)
        ocean_sst = np.asarray(z["earth__annual_sst_C"], dtype=np.float64)
        surface_temp = np.where(land, land_temp, ocean_sst)
        precip = np.asarray(z["earth__annual_precip_mm"], dtype=np.float64)
        seasonal_land_temp = np.asarray(z["earth__seasonal_temperature_C"], dtype=np.float64)
        seasonal_sst = np.asarray(z["earth__seasonal_sst_C"], dtype=np.float64)
        seasonal_surface_temp = np.where(land[None, :], seasonal_land_temp, seasonal_sst)
        seasonal_precip = np.asarray(
            z["earth__seasonal_precip_mm_yr_equiv"],
            dtype=np.float64,
        )
        biome = np.asarray(z["earth__biome_class_proxy"], dtype=np.int16)
        current_speed = np.asarray(
            z["earth__annual_surface_current_speed_m_s"],
            dtype=np.float64,
        )
        land_cover = (
            np.asarray(z["earth__esa_cci_land_cover_broad_class"], dtype=np.uint8)
            if "earth__esa_cci_land_cover_broad_class" in z.files
            else None
        )

    land_amp = _seasonal_amplitude(seasonal_land_temp)
    ocean_amp = _seasonal_amplitude(seasonal_sst)
    precip_seasonality = _precip_seasonality_from_annual_equiv(seasonal_precip, precip)
    lat_means = _lat_band_means(lat, surface_temp, area)
    metrics: dict[str, Any] = {
        "source_npz": str(path),
        "land_fraction": _weighted_fraction(area, land),
        "global_mean_temperature_C": _weighted_mean(surface_temp, area, np.isfinite(surface_temp)),
        "land_mean_temperature_C": _weighted_mean(land_temp, area, land),
        "ocean_mean_temperature_C": _weighted_mean(ocean_sst, area, ocean),
        "land_precip_mean_mm_yr": _weighted_mean(precip, area, land),
        "land_precip_p50_mm_yr": _percentile(precip, 50, land),
        "land_precip_p90_mm_yr": _percentile(precip, 90, land),
        "land_seasonal_temp_amp_p50_C": _percentile(land_amp, 50, land),
        "ocean_seasonal_temp_amp_p50_C": _percentile(ocean_amp, 50, ocean),
        "precip_seasonality_land_p75": _percentile(precip_seasonality, 75, land),
        "current_speed_p50_m_s": _percentile(current_speed, 50, np.isfinite(current_speed)),
        "current_speed_p90_m_s": _percentile(current_speed, 90, np.isfinite(current_speed)),
        "lat_band_mean_temperature_C": lat_means,
        "max_adjacent_lat_band_delta_C": _max_adjacent_band_delta(lat_means),
    }
    metrics.update(_biome_area_fractions(biome, area))
    metrics["biome_desert_area_fraction"] = metrics["biome_desert_area_fraction"]
    metrics["biome_forest_area_fraction"] = metrics["biome_forest_area_fraction"]
    metrics["biome_tropical_area_fraction"] = metrics["biome_tropical_area_fraction"]
    if land_cover is not None:
        metrics.update(_land_cover_area_fractions(land_cover, area))
    return metrics


def generated_world_metrics(summary: dict[str, Any]) -> dict[str, Any]:
    arrays_path = Path(summary["arrays"])
    with np.load(arrays_path, allow_pickle=False) as z:
        area = np.asarray(z["cell_area"], dtype=np.float64)
        lat = np.asarray(z["lat"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"]).ravel()[0])
        elevation = np.asarray(z["terrain__elevation_m"], dtype=np.float64)
        land = elevation >= sea
        ocean = ~land
        surface_temp = np.asarray(z["climate__surface_temperature"], dtype=np.float64) - 273.15
        seasonal_temp = np.asarray(z["climate__seasonal_temperature"], dtype=np.float64) - 273.15
        precip = np.asarray(z["climate__precipitation"], dtype=np.float64)
        seasonal_precip = np.asarray(z["climate__seasonal_precipitation"], dtype=np.float64)
        precip_seasonality = np.asarray(
            z["climate__precipitation_seasonality"],
            dtype=np.float64,
        )
        biome = np.asarray(z["biosphere__biome"], dtype=np.float64).astype(np.int16)
        current = np.asarray(z["ocean__currents"], dtype=np.float64)

    amp = _seasonal_amplitude(seasonal_temp)
    current_speed = np.linalg.norm(current, axis=1)
    current_valid = ocean & np.isfinite(current_speed)
    lat_means = _lat_band_means(lat, surface_temp, area)
    metrics: dict[str, Any] = {
        "label": Path(summary["assets_dir"]).name,
        "preset": summary.get("preset"),
        "seed": int(summary.get("seed", -1)),
        "arrays": str(arrays_path),
        "assets_dir": str(summary["assets_dir"]),
        "land_fraction": _weighted_fraction(area, land),
        "global_mean_temperature_C": _weighted_mean(surface_temp, area, np.isfinite(surface_temp)),
        "land_mean_temperature_C": _weighted_mean(surface_temp, area, land),
        "ocean_mean_temperature_C": _weighted_mean(surface_temp, area, ocean),
        "land_precip_mean_mm_yr": _weighted_mean(precip, area, land),
        "land_precip_p50_mm_yr": _percentile(precip, 50, land),
        "land_precip_p90_mm_yr": _percentile(precip, 90, land),
        "land_seasonal_temp_amp_p50_C": _percentile(amp, 50, land),
        "ocean_seasonal_temp_amp_p50_C": _percentile(amp, 50, ocean),
        "precip_seasonality_land_p75": _percentile(precip_seasonality, 75, land),
        "current_speed_p50_m_s": _percentile(current_speed, 50, current_valid),
        "current_speed_p90_m_s": _percentile(current_speed, 90, current_valid),
        "lat_band_mean_temperature_C": lat_means,
        "max_adjacent_lat_band_delta_C": _max_adjacent_band_delta(lat_means),
        "seasonal_precip_aggregate_matches": bool(
            np.nanmax(np.abs(np.mean(seasonal_precip, axis=0) - precip)) < 1.0e-6
        ),
    }
    metrics.update(_biome_area_fractions(biome, area))
    metrics["biome_desert_area_fraction"] = metrics["biome_desert_area_fraction"]
    metrics["biome_forest_area_fraction"] = metrics["biome_forest_area_fraction"]
    metrics["biome_tropical_area_fraction"] = metrics["biome_tropical_area_fraction"]
    return metrics


def scalar_comparison(
    generated: dict[str, Any],
    earth: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    out: dict[str, dict[str, float | None]] = {}
    for key, scale in SCALAR_COMPARISON_SCALES.items():
        gv = generated.get(key)
        ev = earth.get(key)
        if gv is None or ev is None or not np.isfinite(gv) or not np.isfinite(ev):
            out[key] = {
                "generated": None if gv is None else float(gv),
                "earth": None if ev is None else float(ev),
                "delta": None,
                "normalized_delta": None,
            }
            continue
        delta = float(gv) - float(ev)
        out[key] = {
            "generated": float(gv),
            "earth": float(ev),
            "delta": delta,
            "normalized_delta": abs(delta) / max(float(scale), 1.0e-12),
        }
    return out


def earthlike_flags(comparison: dict[str, dict[str, float | None]]) -> list[str]:
    flags: list[str] = []

    def delta(key: str) -> float:
        value = comparison.get(key, {}).get("delta")
        return float(value) if value is not None and np.isfinite(value) else 0.0

    def generated(key: str) -> float:
        value = comparison.get(key, {}).get("generated")
        return float(value) if value is not None and np.isfinite(value) else float("nan")

    def earth(key: str) -> float:
        value = comparison.get(key, {}).get("earth")
        return float(value) if value is not None and np.isfinite(value) else float("nan")

    if abs(delta("land_fraction")) > 0.18:
        flags.append("earthlike_land_fraction_far_from_earth")
    if abs(delta("global_mean_temperature_C")) > 6.0:
        flags.append("earthlike_global_temperature_far_from_earth")
    if generated("land_precip_mean_mm_yr") < 0.45 * earth("land_precip_mean_mm_yr"):
        flags.append("earthlike_land_precip_too_dry")
    if generated("land_precip_mean_mm_yr") > 1.9 * earth("land_precip_mean_mm_yr"):
        flags.append("earthlike_land_precip_too_wet")
    if generated("precip_seasonality_land_p75") > earth("precip_seasonality_land_p75") + 1.2:
        flags.append("earthlike_precip_too_seasonally_peaked")
    if generated("land_seasonal_temp_amp_p50_C") > earth("land_seasonal_temp_amp_p50_C") + 8.0:
        flags.append("earthlike_land_temperature_too_seasonal")
    if generated("current_speed_p90_m_s") > 2.2 * earth("current_speed_p90_m_s"):
        flags.append("earthlike_current_speed_too_strong")
    if generated("biome_desert_area_fraction") > earth("biome_desert_area_fraction") + 0.20:
        flags.append("earthlike_desert_fraction_too_high")
    return flags


def _distance_score(comparison: dict[str, dict[str, float | None]]) -> float:
    values = [
        row["normalized_delta"]
        for row in comparison.values()
        if row.get("normalized_delta") is not None
    ]
    if not values:
        return float("nan")
    return float(np.mean(np.clip(np.asarray(values, dtype=np.float64), 0.0, 3.0)))


def _copy_contact_sheet_tile(
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
        draw.rectangle([x, y + label_h, x + tile_w, y + label_h + tile_h], outline=(180, 180, 180))
        draw.text((x + 12, y + label_h + 24), "missing", fill=(160, 0, 0), font=font)
        return
    with Image.open(image_path) as raw:
        image = raw.convert("RGB")
        image.thumbnail((tile_w, tile_h), Image.Resampling.LANCZOS)
        bg = Image.new("RGB", (tile_w, tile_h), "white")
        bg.paste(image, ((tile_w - image.width) // 2, (tile_h - image.height) // 2))
        canvas.paste(bg, (x, y + label_h))


def _latlon_raster_index(
    lat: np.ndarray,
    lon: np.ndarray,
    *,
    width: int = 360,
    height: int = 180,
) -> np.ndarray:
    lat_rad = np.radians(np.asarray(lat, dtype=np.float64))
    lon_rad = np.radians(np.asarray(lon, dtype=np.float64))
    source = np.column_stack([
        np.cos(lat_rad) * np.cos(lon_rad),
        np.cos(lat_rad) * np.sin(lon_rad),
        np.sin(lat_rad),
    ])
    out_lon = np.linspace(-180.0, 180.0, int(width))
    out_lat = np.linspace(90.0, -90.0, int(height))
    lat_grid, lon_grid = np.meshgrid(out_lat, out_lon, indexing="ij")
    q_lat = np.radians(lat_grid.ravel())
    q_lon = np.radians(lon_grid.ravel())
    query = np.column_stack([
        np.cos(q_lat) * np.cos(q_lon),
        np.cos(q_lat) * np.sin(q_lon),
        np.sin(q_lat),
    ])
    _, idx = cKDTree(source).query(query)
    return idx.reshape(int(height), int(width))


def _finite_limits(
    values: np.ndarray,
    *,
    lower_q: float = 2.0,
    upper_q: float = 98.0,
    vmin: float | None = None,
    vmax: float | None = None,
    fallback: tuple[float, float] = (0.0, 1.0),
) -> tuple[float, float]:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    lo = float(vmin) if vmin is not None else float(fallback[0])
    hi = float(vmax) if vmax is not None else float(fallback[1])
    if finite.size:
        if vmin is None:
            lo = float(np.nanpercentile(finite, lower_q))
        if vmax is None:
            hi = float(np.nanpercentile(finite, upper_q))
    if not np.isfinite(lo):
        lo = float(fallback[0])
    if not np.isfinite(hi):
        hi = float(fallback[1])
    if hi <= lo:
        span = max(abs(lo) * 0.05, 1.0e-6)
        lo -= span
        hi += span
    return lo, hi


def _render_generated_scalar_preview(
    raster_index: np.ndarray,
    values: np.ndarray,
    out_path: Path,
    *,
    title: str,
    cmap,
    vmin: float | None = None,
    vmax: float | None = None,
) -> Path:
    values = np.asarray(values, dtype=np.float64)
    raster = values[raster_index]
    lo, hi = _finite_limits(values, vmin=vmin, vmax=vmax)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.imshow(
        np.ma.masked_invalid(raster),
        cmap=cmap,
        vmin=lo,
        vmax=hi,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.7)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _render_generated_biome_preview(
    raster_index: np.ndarray,
    values: np.ndarray,
    out_path: Path,
) -> Path:
    raster = np.asarray(values, dtype=np.int16)[raster_index]
    cmap = ListedColormap(render.BIOME_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, len(render.BIOME_COLORS) + 0.5), cmap.N)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    ax.set_title("Biomes")
    cb = fig.colorbar(im, ax=ax, shrink=0.7, ticks=list(range(len(render.BIOME_COLORS))))
    cb.ax.set_yticklabels(["ocean", "ice", "desert", "grass", "forest", "tundra", "tropical"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _ensure_generated_preview_assets(run_metrics: dict[str, Any]) -> list[Path]:
    assets = Path(str(run_metrics["assets_dir"]))
    paths = [
        assets / "temperature.png",
        assets / "precip.png",
        assets / "biomes.png",
        assets / "currents.png",
    ]
    if all(path.exists() for path in paths):
        return paths

    arrays_path = Path(str(run_metrics.get("arrays", "")))
    if not arrays_path.exists():
        return paths

    with np.load(arrays_path, allow_pickle=False) as z:
        if "lat" not in z.files or "lon" not in z.files:
            return paths
        raster_index = _latlon_raster_index(z["lat"], z["lon"])
        if not paths[0].exists() and "climate__surface_temperature" in z.files:
            _render_generated_scalar_preview(
                raster_index,
                np.asarray(z["climate__surface_temperature"], dtype=np.float64) - 273.15,
                paths[0],
                title="Surface temperature (C)",
                cmap=render.TEMPERATURE_CMAP,
            )
        if not paths[1].exists() and "climate__precipitation" in z.files:
            precip = np.asarray(z["climate__precipitation"], dtype=np.float64)
            _, vmax = _finite_limits(precip, vmin=0.0, fallback=(0.0, 1000.0))
            _render_generated_scalar_preview(
                raster_index,
                precip,
                paths[1],
                title="Precipitation (mm/yr)",
                cmap=render.PRECIP_CMAP,
                vmin=0.0,
                vmax=max(vmax, 100.0),
            )
        if not paths[2].exists() and "biosphere__biome" in z.files:
            _render_generated_biome_preview(
                raster_index,
                np.asarray(z["biosphere__biome"], dtype=np.int16),
                paths[2],
            )
        if not paths[3].exists() and "ocean__currents" in z.files:
            currents = np.asarray(z["ocean__currents"], dtype=np.float64)
            if currents.ndim == 2 and currents.shape[1] == 3:
                speed = np.linalg.norm(currents, axis=1)
                _, vmax = _finite_limits(speed, vmin=0.0, fallback=(0.0, 0.2))
                _render_generated_scalar_preview(
                    raster_index,
                    speed,
                    paths[3],
                    title="Ocean current speed (m/s)",
                    cmap="cividis",
                    vmin=0.0,
                    vmax=max(vmax, 0.1),
                )
    return paths


def _earth_reference_row_label(earth_reference_npz: Path) -> str:
    for token in earth_reference_npz.parent.name.split("_"):
        if len(token) > 1 and token[0].lower() == "r" and token[1:].isdigit():
            return f"Earth {token.upper()}"
    return "Earth reference"


def render_comparison_contact_sheet(
    earth_reference_npz: Path,
    run_metrics: list[dict[str, Any]],
    out_path: Path,
) -> Path:
    earth_dir = earth_reference_npz.parent
    stem = earth_reference_npz.stem.replace(".npz", "")
    if stem.endswith("cells"):
        earth_prefix = stem
    else:
        earth_prefix = earth_reference_npz.stem
    columns = [
        ("temperature", "temperature.png"),
        ("precip", "precip.png"),
        ("biome", "biomes.png"),
        ("current", "currents.png"),
    ]
    rows: list[tuple[str, list[Path]]] = [(
        _earth_reference_row_label(earth_reference_npz),
        [
            earth_dir / f"{earth_prefix}_temperature.png",
            earth_dir / f"{earth_prefix}_precip.png",
            earth_dir / f"{earth_prefix}_biomes_from_koppen_proxy.png",
            earth_dir / f"{earth_prefix}_current_speed.png",
        ],
    )]
    for row in run_metrics:
        paths = _ensure_generated_preview_assets(row)
        rows.append((
            str(row["label"]),
            paths,
        ))

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
    for r, (row_label, paths) in enumerate(rows):
        y = pad + r * (tile_h + label_h + pad)
        draw.text((pad, y + label_h + tile_h // 2 - 8), row_label, fill=(0, 0, 0), font=font)
        for c, ((col_label, _), path) in enumerate(zip(columns, paths)):
            x = left_w + c * (tile_w + pad)
            label = col_label if r == 0 else ""
            _copy_contact_sheet_tile(
                canvas,
                draw,
                path,
                label,
                x=x,
                y=y,
                tile_w=tile_w,
                tile_h=tile_h,
                label_h=label_h,
                font=font,
            )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_path)
    return out_path


def _write_csv(path: Path, rows: list[dict[str, Any]], keys: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in keys})
    return path


def run_earth_climate_comparison(
    config: EarthClimateComparisonConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    earth_metrics = earth_reference_metrics(Path(config.earth_reference_npz))
    terminal = json.loads(Path(config.terminal_summary_json).read_text())
    run_entries: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []
    for summary in terminal["summaries"]:
        metrics = generated_world_metrics(summary)
        comparison = scalar_comparison(metrics, earth_metrics)
        is_earthlike = "earthlike" in str(metrics.get("preset", "")).lower()
        flags = earthlike_flags(comparison) if is_earthlike else []
        score = _distance_score(comparison)
        entry = {
            "label": metrics["label"],
            "preset": metrics["preset"],
            "seed": metrics["seed"],
            "mode": "earthlike_calibration" if is_earthlike else "diagnostic_only",
            "earth_distance_score": score,
            "flags": flags,
            "metrics": metrics,
            "comparison": comparison,
        }
        run_entries.append(entry)
        csv_rows.append({
            "label": metrics["label"],
            "preset": metrics["preset"],
            "seed": metrics["seed"],
            "mode": entry["mode"],
            "earth_distance_score": score,
            "flag_count": len(flags),
            "flags": ";".join(flags),
            "land_fraction": metrics["land_fraction"],
            "global_mean_temperature_C": metrics["global_mean_temperature_C"],
            "land_mean_temperature_C": metrics["land_mean_temperature_C"],
            "land_precip_mean_mm_yr": metrics["land_precip_mean_mm_yr"],
            "land_precip_p50_mm_yr": metrics["land_precip_p50_mm_yr"],
            "land_precip_p90_mm_yr": metrics["land_precip_p90_mm_yr"],
            "land_seasonal_temp_amp_p50_C": metrics["land_seasonal_temp_amp_p50_C"],
            "ocean_seasonal_temp_amp_p50_C": metrics["ocean_seasonal_temp_amp_p50_C"],
            "precip_seasonality_land_p75": metrics["precip_seasonality_land_p75"],
            "current_speed_p90_m_s": metrics["current_speed_p90_m_s"],
            "biome_desert_area_fraction": metrics["biome_desert_area_fraction"],
            "biome_forest_area_fraction": metrics["biome_forest_area_fraction"],
            "biome_tropical_area_fraction": metrics["biome_tropical_area_fraction"],
        })

    run_entries.sort(key=lambda row: (row["mode"] != "earthlike_calibration", row["label"]))
    csv_rows.sort(key=lambda row: (row["mode"] != "earthlike_calibration", row["label"]))
    csv_path = _write_csv(
        outdir / "earth_climate_comparison_metrics.csv",
        csv_rows,
        [
            "label", "preset", "seed", "mode", "earth_distance_score",
            "flag_count", "flags", "land_fraction", "global_mean_temperature_C",
            "land_mean_temperature_C", "land_precip_mean_mm_yr",
            "land_precip_p50_mm_yr", "land_precip_p90_mm_yr",
            "land_seasonal_temp_amp_p50_C", "ocean_seasonal_temp_amp_p50_C",
            "precip_seasonality_land_p75", "current_speed_p90_m_s",
            "biome_desert_area_fraction", "biome_forest_area_fraction",
            "biome_tropical_area_fraction",
        ],
    )
    contact_sheet = None
    if config.render_contact_sheet:
        contact_sheet = render_comparison_contact_sheet(
            Path(config.earth_reference_npz),
            [entry["metrics"] for entry in run_entries],
            outdir / "earth_vs_generated_climate_contact_sheet.png",
        )

    summary = {
        "schema": SCHEMA,
        "earth_reference_npz": str(config.earth_reference_npz),
        "terminal_summary_json": str(config.terminal_summary_json),
        "earth_metrics": earth_metrics,
        "run_count": len(run_entries),
        "earthlike_flagged_count": int(sum(1 for row in run_entries if row["flags"])),
        "entries": run_entries,
        "csv": str(csv_path),
        "contact_sheet": None if contact_sheet is None else str(contact_sheet),
    }
    (outdir / "earth_climate_comparison_summary.json").write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    return summary
