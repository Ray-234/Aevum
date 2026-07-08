"""C4d hydroclimate-region object gate.

This gate checks the seasonal region objects derived from C4d corridor fields.
It intentionally evaluates the object layer, not the precipitation solver:
existing scalar/pattern/biome gates still own broad Earth climate envelopes.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import ConvexHull, QhullError


SCHEMA = "aevum.earth_climate_hydro_region_gate.v1"


@dataclass(frozen=True)
class EarthClimateHydroRegionGateConfig:
    terminal_summary_json: Path
    outdir: Path
    render_contact_sheets: bool = True


def _json_default(value: Any):
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _safe_float(value: Any, default: float = float("nan")) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def _candidate_region_paths(row: dict[str, Any]) -> list[Path]:
    paths: list[Path] = []
    raw = row.get("hydroclimate_regions_json")
    if raw:
        paths.append(Path(str(raw)))
    assets = row.get("assets_dir")
    if assets:
        paths.append(Path(str(assets)) / "hydroclimate_regions.json")
    return paths


def _candidate_array_paths(row: dict[str, Any]) -> list[Path]:
    paths: list[Path] = []
    raw = row.get("arrays")
    if raw:
        paths.append(Path(str(raw)))
    assets = row.get("assets_dir")
    if assets:
        paths.append(Path(str(assets)) / "terminal_climate_arrays.npz")
    return paths


def _load_regions(row: dict[str, Any]) -> tuple[list[dict[str, Any]], Path | None]:
    for path in _candidate_region_paths(row):
        if path.exists():
            payload = json.loads(path.read_text())
            regions = payload.get("regions", [])
            if isinstance(regions, list):
                return [obj for obj in regions if isinstance(obj, dict)], path
            return [], path
    return [], None


def _array_path(row: dict[str, Any]) -> Path | None:
    for path in _candidate_array_paths(row):
        if path.exists():
            return path
    return None


def _values(rows: list[dict[str, Any]], key: str) -> list[float]:
    out: list[float] = []
    for row in rows:
        value = _safe_float(row.get(key, float("nan")))
        if np.isfinite(value):
            out.append(value)
    return out


def _max(rows: list[dict[str, Any]], key: str) -> float:
    vals = _values(rows, key)
    return float(max(vals)) if vals else 0.0


def _sum(rows: list[dict[str, Any]], key: str) -> float:
    vals = _values(rows, key)
    return float(sum(vals)) if vals else 0.0


def _pct(rows: list[dict[str, Any]], key: str, q: float) -> float:
    vals = _values(rows, key)
    return float(np.percentile(vals, q)) if vals else 0.0


def _kind(rows: list[dict[str, Any]], kind: str) -> list[dict[str, Any]]:
    return [row for row in rows if str(row.get("kind", "")) == kind]


def _season_count(rows: list[dict[str, Any]]) -> int:
    return int(len({str(row.get("season", "")) for row in rows if row.get("season")}))


def _weighted_mean(values: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values) & np.isfinite(weights)
    if not mask.any():
        return float("nan")
    denom = float(np.sum(weights[mask]))
    if denom <= 0.0:
        return float("nan")
    return float(np.sum(values[mask] * weights[mask]) / denom)


def _corr(values_a: np.ndarray, values_b: np.ndarray, mask: np.ndarray) -> float:
    a = np.asarray(values_a, dtype=np.float64).ravel()
    b = np.asarray(values_b, dtype=np.float64).ravel()
    m = np.asarray(mask, dtype=bool).ravel()
    m = m & np.isfinite(a) & np.isfinite(b)
    if int(m.sum()) < 3:
        return float("nan")
    aa = a[m] - float(np.mean(a[m]))
    bb = b[m] - float(np.mean(b[m]))
    denom = float(np.sqrt(np.sum(aa * aa) * np.sum(bb * bb)))
    if denom <= 0.0:
        return float("nan")
    return float(np.sum(aa * bb) / denom)


def _fallback_edges(n: int) -> np.ndarray:
    if n <= 1:
        return np.zeros((0, 2), dtype=np.int64)
    edges = [(i, i + 1) for i in range(n - 1)]
    if n > 2:
        edges.append((0, n - 1))
    return np.asarray(edges, dtype=np.int64)


def _edges_from_latlon(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    if lat.ndim != 1 or lon.shape != lat.shape or lat.size < 4:
        return _fallback_edges(int(lat.size))
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    xyz = np.column_stack([
        np.cos(lat_rad) * np.cos(lon_rad),
        np.cos(lat_rad) * np.sin(lon_rad),
        np.sin(lat_rad),
    ])
    try:
        hull = ConvexHull(xyz)
    except QhullError:
        return _fallback_edges(int(lat.size))
    edge_set: set[tuple[int, int]] = set()
    for tri in hull.simplices:
        a, b, c = (int(tri[0]), int(tri[1]), int(tri[2]))
        for u, v in ((a, b), (b, c), (a, c)):
            edge_set.add((u, v) if u < v else (v, u))
    if not edge_set:
        return _fallback_edges(int(lat.size))
    return np.asarray(sorted(edge_set), dtype=np.int64)


def _component_shape_metrics(mask: np.ndarray, area: np.ndarray,
                             edges: np.ndarray) -> dict[str, float]:
    mask = np.asarray(mask, dtype=bool)
    area = np.asarray(area, dtype=np.float64)
    edges = np.asarray(edges, dtype=np.int64)
    active = np.where(mask)[0]
    total_area = max(float(np.sum(area)), 1.0e-12)
    if active.size == 0:
        return {
            "component_count": 0.0,
            "largest_component_share": 0.0,
            "active_world_fraction": 0.0,
            "boundary_per_active_cell": 0.0,
        }
    parent = np.arange(mask.size, dtype=np.int64)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = int(parent[x])
        return x

    def union(a: int, b: int) -> None:
        ra = find(a)
        rb = find(b)
        if ra != rb:
            parent[rb] = ra

    if edges.size:
        in_mask = mask[edges[:, 0]] & mask[edges[:, 1]]
        for i, j in edges[in_mask]:
            union(int(i), int(j))
        boundary_edges = (
            int(np.count_nonzero(mask[edges[:, 0]] ^ mask[edges[:, 1]]))
            if mask.size >= 500 else 0
        )
    else:
        boundary_edges = 0

    component_area: dict[int, float] = {}
    for cell in active:
        root = find(int(cell))
        component_area[root] = component_area.get(root, 0.0) + float(area[cell])
    active_area = float(sum(component_area.values()))
    largest = max(component_area.values(), default=0.0)
    return {
        "component_count": float(len(component_area)),
        "largest_component_share": float(largest / max(active_area, 1.0e-12)),
        "active_world_fraction": float(active_area / total_area),
        "boundary_per_active_cell": float(boundary_edges / max(int(active.size), 1)),
    }


def _threshold_high(values: np.ndarray, mask: np.ndarray,
                    floor: float, percentile: float) -> float:
    vals = np.asarray(values, dtype=np.float64)[np.asarray(mask, dtype=bool)]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float(floor)
    return float(max(floor, np.percentile(vals, percentile) * 0.82))


def _threshold_response_high(values: np.ndarray, mask: np.ndarray) -> float:
    vals = np.asarray(values, dtype=np.float64)[np.asarray(mask, dtype=bool)]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 1.012
    return float(max(1.012, np.percentile(vals, 84)))


def _threshold_response_low(values: np.ndarray, mask: np.ndarray) -> float:
    vals = np.asarray(values, dtype=np.float64)[np.asarray(mask, dtype=bool)]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 0.988
    return float(min(0.988, np.percentile(vals, 16)))


def _seasonal_map_metrics(field: np.ndarray, land: np.ndarray, area: np.ndarray,
                          edges: np.ndarray, *,
                          floor: float | None = None,
                          percentile: float | None = None,
                          response: str | None = None) -> dict[str, float]:
    field = np.asarray(field, dtype=np.float64)
    land = np.asarray(land, dtype=bool)
    if field.ndim != 2 or field.shape[0] != 4 or field.shape[1] != land.size:
        return {
            "component_count_p50": 0.0,
            "largest_component_share_p50": 0.0,
            "active_world_fraction_p50": 0.0,
            "active_land_fraction_p50": 0.0,
            "boundary_per_active_cell_p50": 0.0,
        }
    land_area = max(float(np.sum(np.asarray(area)[land])), 1.0e-12)
    rows: list[dict[str, float]] = []
    for season in range(4):
        values = field[season]
        if response == "wet":
            threshold = _threshold_response_high(values, land)
            active = land & np.isfinite(values) & (values >= threshold)
        elif response == "dry":
            threshold = _threshold_response_low(values, land)
            active = land & np.isfinite(values) & (values <= threshold)
        else:
            threshold = _threshold_high(
                values, land, float(floor or 0.0), float(percentile or 82.0))
            active = land & np.isfinite(values) & (values >= threshold)
        metrics = _component_shape_metrics(active, area, edges)
        metrics["active_land_fraction"] = float(
            np.sum(np.asarray(area)[active]) / land_area)
        rows.append(metrics)

    def median(key: str) -> float:
        values = [row[key] for row in rows if np.isfinite(row[key])]
        return float(np.median(values)) if values else 0.0

    return {
        "component_count_p50": median("component_count"),
        "largest_component_share_p50": median("largest_component_share"),
        "active_world_fraction_p50": median("active_world_fraction"),
        "active_land_fraction_p50": median("active_land_fraction"),
        "boundary_per_active_cell_p50": median("boundary_per_active_cell"),
    }


def _prefixed(prefix: str, metrics: dict[str, float]) -> dict[str, float]:
    return {f"{prefix}_{key}": value for key, value in metrics.items()}


def _contact_sheet_arrays(summary_row: dict[str, Any]) -> dict[str, Any] | None:
    path = _array_path(summary_row)
    if path is None:
        return None
    with np.load(path, allow_pickle=False) as z:
        required = (
            "lat",
            "lon",
            "cell_area",
            "terrain__elevation_m",
            "sea_level_m",
            "climate__seasonal_precipitation",
            "climate__monsoon_rainfall_corridor",
            "climate__storm_track_rainfall_corridor",
            "climate__rain_shadow_index",
            "climate__regional_precipitation_response",
        )
        if any(key not in z for key in required):
            return None
        out = {key: np.asarray(z[key]) for key in required}
    sea = float(np.asarray(out["sea_level_m"], dtype=np.float64).ravel()[0])
    elev = np.asarray(out["terrain__elevation_m"], dtype=np.float64)
    out["land"] = elev >= sea
    out["arrays_path"] = str(path)
    return out


def _active_corridor_mask(
    field: np.ndarray,
    land: np.ndarray,
    *,
    floor: float,
    percentile: float = 82.0,
) -> np.ndarray:
    field = np.asarray(field, dtype=np.float64)
    land = np.asarray(land, dtype=bool)
    out = np.zeros_like(field, dtype=bool)
    if field.ndim != 2 or field.shape[0] != 4:
        return out
    for season in range(4):
        threshold = _threshold_high(field[season], land, floor, percentile)
        out[season] = land & np.isfinite(field[season]) & (field[season] >= threshold)
    return out


def _render_contact_sheet(
    summary_row: dict[str, Any],
    generated_row: dict[str, Any],
    outdir: Path,
) -> Path | None:
    arrays = _contact_sheet_arrays(summary_row)
    if arrays is None:
        return None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from aevum.render import PRECIP_CMAP

    label = str(generated_row.get("label", "world"))
    sheets_dir = outdir / "contact_sheets"
    sheets_dir.mkdir(parents=True, exist_ok=True)
    path = sheets_dir / f"{label}_hydroclimate_region_contact_sheet.png"

    lat = np.asarray(arrays["lat"], dtype=np.float64)
    lon = np.asarray(arrays["lon"], dtype=np.float64)
    land = np.asarray(arrays["land"], dtype=bool)
    seasonal_precip = np.asarray(
        arrays["climate__seasonal_precipitation"], dtype=np.float64)
    monsoon = np.asarray(
        arrays["climate__monsoon_rainfall_corridor"], dtype=np.float64)
    storm = np.asarray(
        arrays["climate__storm_track_rainfall_corridor"], dtype=np.float64)
    shadow = np.asarray(
        arrays["climate__rain_shadow_index"], dtype=np.float64)
    response = np.asarray(
        arrays["climate__regional_precipitation_response"], dtype=np.float64)

    monsoon_mask = _active_corridor_mask(monsoon, land, floor=0.006)
    storm_mask = _active_corridor_mask(storm, land, floor=0.020)
    shadow_mask = _active_corridor_mask(shadow, land, floor=0.055)

    seasons = ("DJF", "MAM", "JJA", "SON")
    n = max(int(lat.size), 1)
    marker_size = float(np.clip(28000.0 / n, 0.45, 10.0))
    land_color = np.where(land, "#dfdac8", "#d9edf4")

    fig, axes = plt.subplots(
        5,
        4,
        figsize=(14.5, 10.5),
        constrained_layout=True,
    )
    fig.suptitle(f"C4d hydroclimate region contact sheet: {label}", fontsize=13)

    def base(ax):
        ax.scatter(lon, lat, c=land_color, s=marker_size, linewidths=0,
                   rasterized=True)
        ax.set_xlim(-180.0, 180.0)
        ax.set_ylim(-90.0, 90.0)
        ax.set_xticks([])
        ax.set_yticks([])

    def panel(ax, values, title, cmap, vmin, vmax, mask=None):
        base(ax)
        vals = np.asarray(values, dtype=np.float64)
        if mask is not None:
            vals = np.where(mask, vals, np.nan)
        finite = np.isfinite(vals)
        if finite.any():
            im = ax.scatter(lon[finite], lat[finite], c=vals[finite],
                            s=marker_size * 1.25, linewidths=0, cmap=cmap,
                            vmin=vmin, vmax=vmax, rasterized=True)
        else:
            im = ax.scatter([], [], c=[], cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=8)
        return im

    precip_land = seasonal_precip[:, land] if land.any() else seasonal_precip
    precip_vmax = max(float(np.nanpercentile(precip_land, 98)), 100.0)
    land4 = np.broadcast_to(land[None, :], monsoon.shape)
    monsoon_vmax = max(float(np.nanpercentile(monsoon[land4], 98))
                       if land.any() else float(np.nanpercentile(monsoon, 98)),
                       0.10)
    storm_vmax = max(float(np.nanpercentile(storm[land4], 98))
                     if land.any() else float(np.nanpercentile(storm, 98)),
                     0.10)
    shadow_vmax = max(float(np.nanpercentile(shadow[land4], 98))
                      if land.any() else float(np.nanpercentile(shadow, 98)),
                      0.10)
    response_anom = response - 1.0
    resp_vmax = max(float(np.nanpercentile(np.abs(response_anom), 98)), 0.08)

    row_maps = [
        (
            "seasonal precipitation",
            seasonal_precip,
            None,
            PRECIP_CMAP,
            0.0,
            precip_vmax,
        ),
        (
            "monsoon corridor mask",
            monsoon,
            monsoon_mask,
            "YlGnBu",
            0.0,
            monsoon_vmax,
        ),
        (
            "storm-track corridor mask",
            storm,
            storm_mask,
            "PuBuGn",
            0.0,
            storm_vmax,
        ),
        (
            "rain-shadow mask",
            shadow,
            shadow_mask,
            "YlOrBr",
            0.0,
            shadow_vmax,
        ),
        (
            "regional response anomaly",
            response_anom,
            None,
            "coolwarm",
            -resp_vmax,
            resp_vmax,
        ),
    ]

    row_images = []
    for row_idx, (row_title, field, mask, cmap, vmin, vmax) in enumerate(row_maps):
        image = None
        for season_idx, season in enumerate(seasons):
            image = panel(
                axes[row_idx, season_idx],
                field[season_idx],
                f"{season} {row_title}",
                cmap,
                vmin,
                vmax,
                None if mask is None else mask[season_idx],
            )
        row_images.append(image)

    for row_idx, image in enumerate(row_images):
        fig.colorbar(image, ax=axes[row_idx, :].ravel().tolist(), shrink=0.72)

    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_contact_sheets(
    summary_rows: list[dict[str, Any]],
    generated_rows: list[dict[str, Any]],
    outdir: Path,
) -> list[dict[str, Any]]:
    sheets: list[dict[str, Any]] = []
    by_assets = {
        str(row.get("assets_dir", "")): row
        for row in generated_rows
    }
    for summary_row in summary_rows:
        assets = str(summary_row.get("assets_dir", ""))
        generated_row = by_assets.get(assets)
        if generated_row is None:
            continue
        path = _render_contact_sheet(summary_row, generated_row, outdir)
        if path is None:
            continue
        sheets.append({
            "label": str(generated_row.get("label", "")),
            "preset": str(generated_row.get("preset", "")),
            "path": str(path),
        })
    if sheets:
        manifest = {
            "schema": "aevum.earth_climate_hydro_region_contact_sheets.v1",
            "sheet_count": int(len(sheets)),
            "sheets": sheets,
        }
        (outdir / "earth_climate_hydro_region_contact_sheets.json").write_text(
            json.dumps(manifest, indent=2, default=_json_default),
        )
    return sheets


def _object_season_weighted_lat(
    rows: list[dict[str, Any]],
    kind: str,
    season: str,
) -> float:
    selected = [
        row for row in rows
        if str(row.get("kind", "")) == kind and str(row.get("season", "")) == season
    ]
    weights = [_safe_float(row.get("area_fraction", 0.0), 0.0) for row in selected]
    lats = [_safe_float(row.get("centroid_lat", float("nan"))) for row in selected]
    pairs = [(lat, weight) for lat, weight in zip(lats, weights)
             if np.isfinite(lat) and np.isfinite(weight) and weight > 0.0]
    denom = float(sum(weight for _, weight in pairs))
    if denom <= 0.0:
        return 0.0
    return float(sum(lat * weight for lat, weight in pairs) / denom)


def _placement_metrics(summary_row: dict[str, Any],
                       regions: list[dict[str, Any]]) -> dict[str, Any]:
    defaults = {
        "arrays_found": 0.0,
        "arrays_path": "",
        "monsoon_field_lat_shift_jja_minus_djf": 0.0,
        "monsoon_object_lat_shift_jja_minus_djf": 0.0,
        "storm_track_abs_lat_weighted_median": 0.0,
        "storm_track_coast_distance_weighted_median": 1.0,
        "rain_shadow_dry_response_corr": 0.0,
        "hydro_map_edges_found": 0.0,
        "monsoon_map_component_count_p50": 0.0,
        "monsoon_map_largest_component_share_p50": 0.0,
        "monsoon_map_active_world_fraction_p50": 0.0,
        "monsoon_map_active_land_fraction_p50": 0.0,
        "monsoon_map_boundary_per_active_cell_p50": 0.0,
        "storm_track_map_component_count_p50": 0.0,
        "storm_track_map_largest_component_share_p50": 0.0,
        "storm_track_map_active_world_fraction_p50": 0.0,
        "storm_track_map_active_land_fraction_p50": 0.0,
        "storm_track_map_boundary_per_active_cell_p50": 0.0,
        "rain_shadow_map_component_count_p50": 0.0,
        "rain_shadow_map_largest_component_share_p50": 0.0,
        "rain_shadow_map_active_world_fraction_p50": 0.0,
        "rain_shadow_map_active_land_fraction_p50": 0.0,
        "rain_shadow_map_boundary_per_active_cell_p50": 0.0,
        "wet_response_map_component_count_p50": 0.0,
        "wet_response_map_largest_component_share_p50": 0.0,
        "wet_response_map_active_world_fraction_p50": 0.0,
        "wet_response_map_active_land_fraction_p50": 0.0,
        "wet_response_map_boundary_per_active_cell_p50": 0.0,
        "dry_response_map_component_count_p50": 0.0,
        "dry_response_map_largest_component_share_p50": 0.0,
        "dry_response_map_active_world_fraction_p50": 0.0,
        "dry_response_map_active_land_fraction_p50": 0.0,
        "dry_response_map_boundary_per_active_cell_p50": 0.0,
    }
    path = _array_path(summary_row)
    if path is None:
        return defaults

    with np.load(path, allow_pickle=False) as z:
        required = (
            "lat",
            "lon",
            "cell_area",
            "terrain__elevation_m",
            "sea_level_m",
            "climate__monsoon_rainfall_corridor",
            "climate__storm_track_rainfall_corridor",
            "climate__rain_shadow_index",
            "climate__regional_precipitation_response",
        )
        if any(key not in z for key in required):
            out = dict(defaults)
            out["arrays_found"] = 1.0
            out["arrays_path"] = str(path)
            return out
        lat = np.asarray(z["lat"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        elev = np.asarray(z["terrain__elevation_m"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        land = elev >= sea
        monsoon = np.asarray(z["climate__monsoon_rainfall_corridor"],
                             dtype=np.float64)
        storm = np.asarray(z["climate__storm_track_rainfall_corridor"],
                           dtype=np.float64)
        shadow = np.asarray(z["climate__rain_shadow_index"], dtype=np.float64)
        response = np.asarray(z["climate__regional_precipitation_response"],
                              dtype=np.float64)
        coast_distance = (
            np.asarray(z["climate__coast_distance"], dtype=np.float64)
            if "climate__coast_distance" in z
            else np.ones_like(lat, dtype=np.float64)
        )
        lon = np.asarray(z["lon"], dtype=np.float64)

    def field_lat(field: np.ndarray, season: int) -> float:
        weights = np.maximum(field[season], 0.0) * area
        return _weighted_mean(lat, weights, land)

    def storm_abs_lat(season: int) -> float:
        weights = np.maximum(storm[season], 0.0) * area
        return _weighted_mean(np.abs(lat), weights, land)

    def storm_coast(season: int) -> float:
        weights = np.maximum(storm[season], 0.0) * area
        return _weighted_mean(coast_distance, weights, land)

    monsoon_lats = [field_lat(monsoon, season) for season in range(4)]
    storm_abs = [storm_abs_lat(season) for season in range(4)]
    storm_coast_values = [storm_coast(season) for season in range(4)]
    dry_response = np.maximum(1.0 - response, 0.0)
    land4 = np.broadcast_to(land[None, :], shadow.shape)
    object_shift = (
        _object_season_weighted_lat(regions, "monsoon_rainfall_corridor", "JJA")
        - _object_season_weighted_lat(regions, "monsoon_rainfall_corridor", "DJF")
    )
    edges = _edges_from_latlon(lat, lon)
    out = {
        "arrays_found": 1.0,
        "arrays_path": str(path),
        "monsoon_field_lat_shift_jja_minus_djf": (
            float(monsoon_lats[2] - monsoon_lats[0])
            if np.isfinite(monsoon_lats[2]) and np.isfinite(monsoon_lats[0])
            else 0.0
        ),
        "monsoon_object_lat_shift_jja_minus_djf": float(object_shift),
        "storm_track_abs_lat_weighted_median": (
            float(np.nanmedian(storm_abs)) if storm_abs else 0.0
        ),
        "storm_track_coast_distance_weighted_median": (
            float(np.nanmedian(storm_coast_values)) if storm_coast_values else 1.0
        ),
        "rain_shadow_dry_response_corr": _corr(shadow, dry_response, land4),
        "hydro_map_edges_found": 1.0 if edges.size else 0.0,
    }
    out.update(_prefixed(
        "monsoon_map",
        _seasonal_map_metrics(monsoon, land, area, edges, floor=0.006,
                              percentile=82.0),
    ))
    out.update(_prefixed(
        "storm_track_map",
        _seasonal_map_metrics(storm, land, area, edges, floor=0.020,
                              percentile=82.0),
    ))
    out.update(_prefixed(
        "rain_shadow_map",
        _seasonal_map_metrics(shadow, land, area, edges, floor=0.055,
                              percentile=82.0),
    ))
    out.update(_prefixed(
        "wet_response_map",
        _seasonal_map_metrics(response, land, area, edges, response="wet"),
    ))
    out.update(_prefixed(
        "dry_response_map",
        _seasonal_map_metrics(response, land, area, edges, response="dry"),
    ))
    return out


def _generated_row(summary_row: dict[str, Any]) -> dict[str, Any]:
    label = Path(str(summary_row.get("assets_dir", summary_row.get("label", "")))).name
    regions, path = _load_regions(summary_row)
    kinds = Counter(str(row.get("kind", "unknown")) for row in regions)
    seasons = Counter(str(row.get("season", "unknown")) for row in regions)
    monsoon = _kind(regions, "monsoon_rainfall_corridor")
    storm = _kind(regions, "storm_track_rainfall_corridor")
    shadow = _kind(regions, "rain_shadow_region")
    wet = _kind(regions, "wet_regional_precipitation_response")
    dry = _kind(regions, "dry_regional_precipitation_response")
    cell_counts = _values(regions, "cell_count")
    mean_intensities = _values(regions, "mean_intensity")
    row = {
        "label": label,
        "preset": str(summary_row.get("preset", "")),
        "seed": int(summary_row.get("seed", -1)),
        "assets_dir": str(summary_row.get("assets_dir", "")),
        "regions_json": str(path) if path else "",
        "object_json_found": 1.0 if path is not None else 0.0,
        "object_count": int(len(regions)),
        "kind_count": int(len(kinds)),
        "season_count": int(len(seasons)),
        "kind_counts": dict(sorted(kinds.items())),
        "season_counts": dict(sorted(seasons.items())),
        "largest_area_fraction": _max(regions, "area_fraction"),
        "top5_area_fraction_sum": float(sum(sorted(
            _values(regions, "area_fraction"), reverse=True)[:5])),
        "cell_count_p50": float(np.percentile(cell_counts, 50)) if cell_counts else 0.0,
        "cell_count_p90": float(np.percentile(cell_counts, 90)) if cell_counts else 0.0,
        "mean_intensity_p50": (
            float(np.percentile(mean_intensities, 50)) if mean_intensities else 0.0
        ),
        "monsoon_region_object_count": int(len(monsoon)),
        "storm_track_region_object_count": int(len(storm)),
        "rain_shadow_region_object_count": int(len(shadow)),
        "wet_response_region_object_count": int(len(wet)),
        "dry_response_region_object_count": int(len(dry)),
        "monsoon_region_season_count": _season_count(monsoon),
        "storm_track_region_season_count": _season_count(storm),
        "rain_shadow_region_season_count": _season_count(shadow),
        "wet_response_region_season_count": _season_count(wet),
        "dry_response_region_season_count": _season_count(dry),
        "largest_monsoon_area_fraction": _max(monsoon, "area_fraction"),
        "largest_storm_track_area_fraction": _max(storm, "area_fraction"),
        "largest_rain_shadow_area_fraction": _max(shadow, "area_fraction"),
        "largest_wet_response_area_fraction": _max(wet, "area_fraction"),
        "largest_dry_response_area_fraction": _max(dry, "area_fraction"),
        "monsoon_area_fraction_sum": _sum(monsoon, "area_fraction"),
        "storm_track_area_fraction_sum": _sum(storm, "area_fraction"),
        "rain_shadow_area_fraction_sum": _sum(shadow, "area_fraction"),
        "wet_response_area_fraction_sum": _sum(wet, "area_fraction"),
        "dry_response_area_fraction_sum": _sum(dry, "area_fraction"),
        "monsoon_mean_intensity_p75": _pct(monsoon, "mean_intensity", 75),
        "storm_track_mean_intensity_p75": _pct(storm, "mean_intensity", 75),
        "rain_shadow_mean_intensity_p75": _pct(shadow, "mean_intensity", 75),
    }
    row.update(_placement_metrics(summary_row, regions))
    return row


def _check(
    *,
    label: str,
    group: str,
    metric: str,
    generated: float,
    operator: str,
    threshold: float,
    severity: str,
    message: str,
) -> dict[str, Any]:
    generated = _safe_float(generated)
    skipped = not np.isfinite(generated)
    if skipped:
        passed = True
    elif operator == ">=":
        passed = generated >= threshold
    elif operator == "<=":
        passed = generated <= threshold
    else:
        raise ValueError(f"unsupported operator {operator!r}")
    return {
        "label": label,
        "group": group,
        "metric": metric,
        "generated": generated,
        "operator": operator,
        "threshold": threshold,
        "severity": severity,
        "passed": bool(passed),
        "skipped": bool(skipped),
        "message": message,
    }


def _earthlike_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(label=label, group="object_archive", metric="object_json_found",
               generated=row["object_json_found"], operator=">=", threshold=1.0,
               severity="fail", message="earthlike runs must archive hydroclimate region objects"),
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="earthlike runs must archive C4d arrays for placement checks"),
        _check(label=label, group="map_archive", metric="hydro_map_edges_found",
               generated=row["hydro_map_edges_found"], operator=">=", threshold=1.0,
               severity="fail", message="earthlike C4d arrays must support map-readability checks"),
        _check(label=label, group="object_coverage", metric="object_count",
               generated=row["object_count"], operator=">=", threshold=80.0,
               severity="fail", message="earthlike hydroclimate should expose enough regional objects"),
        _check(label=label, group="object_coverage", metric="kind_count",
               generated=row["kind_count"], operator=">=", threshold=5.0,
               severity="fail", message="earthlike C4d objects should cover all wet/dry corridor kinds"),
        _check(label=label, group="object_coverage", metric="season_count",
               generated=row["season_count"], operator=">=", threshold=4.0,
               severity="fail", message="earthlike C4d objects should be seasonal, not annual-only"),
        _check(label=label, group="monsoon_regions", metric="monsoon_region_object_count",
               generated=row["monsoon_region_object_count"], operator=">=", threshold=8.0,
               severity="fail", message="earthlike worlds should have monsoon-rainfall corridor objects"),
        _check(label=label, group="storm_track_regions", metric="storm_track_region_object_count",
               generated=row["storm_track_region_object_count"], operator=">=", threshold=12.0,
               severity="fail", message="earthlike worlds should have storm-track wet-region objects"),
        _check(label=label, group="rain_shadow_regions", metric="rain_shadow_region_object_count",
               generated=row["rain_shadow_region_object_count"], operator=">=", threshold=12.0,
               severity="fail", message="earthlike worlds should have leeward rain-shadow objects"),
        _check(label=label, group="wet_dry_response", metric="wet_response_region_object_count",
               generated=row["wet_response_region_object_count"], operator=">=", threshold=8.0,
               severity="fail", message="earthlike worlds should archive wet regional-response objects"),
        _check(label=label, group="wet_dry_response", metric="dry_response_region_object_count",
               generated=row["dry_response_region_object_count"], operator=">=", threshold=8.0,
               severity="fail", message="earthlike worlds should archive dry regional-response objects"),
        _check(label=label, group="seasonal_organization", metric="monsoon_region_season_count",
               generated=row["monsoon_region_season_count"], operator=">=", threshold=2.0,
               severity="fail", message="monsoon objects should not be confined to one season"),
        _check(label=label, group="seasonal_organization", metric="storm_track_region_season_count",
               generated=row["storm_track_region_season_count"], operator=">=", threshold=3.0,
               severity="fail", message="storm-track regions should appear across several seasons"),
        _check(label=label, group="coherence_proxy", metric="largest_area_fraction",
               generated=row["largest_area_fraction"], operator=">=", threshold=0.006,
               severity="fail", message="C4d objects should include coherent regions, not only speckles"),
        _check(label=label, group="coherence_proxy", metric="largest_area_fraction",
               generated=row["largest_area_fraction"], operator="<=", threshold=0.18,
               severity="warn", message="single C4d objects should not swallow continental-scale climate zones"),
        _check(label=label, group="monsoon_seasonal_migration",
               metric="monsoon_field_lat_shift_jja_minus_djf",
               generated=row["monsoon_field_lat_shift_jja_minus_djf"],
               operator=">=", threshold=30.0, severity="fail",
               message="monsoon rainfall corridors should migrate toward the summer hemisphere"),
        _check(label=label, group="monsoon_seasonal_migration",
               metric="monsoon_object_lat_shift_jja_minus_djf",
               generated=row["monsoon_object_lat_shift_jja_minus_djf"],
               operator=">=", threshold=20.0, severity="fail",
               message="monsoon objects should record summer-hemisphere migration"),
        _check(label=label, group="storm_track_placement",
               metric="storm_track_abs_lat_weighted_median",
               generated=row["storm_track_abs_lat_weighted_median"],
               operator=">=", threshold=25.0, severity="fail",
               message="storm-track rainfall corridors should sit outside the deep tropics"),
        _check(label=label, group="storm_track_placement",
               metric="storm_track_abs_lat_weighted_median",
               generated=row["storm_track_abs_lat_weighted_median"],
               operator="<=", threshold=55.0, severity="fail",
               message="storm-track rainfall corridors should not collapse onto polar land"),
        _check(label=label, group="storm_track_placement",
               metric="storm_track_coast_distance_weighted_median",
               generated=row["storm_track_coast_distance_weighted_median"],
               operator="<=", threshold=0.22, severity="fail",
               message="storm-track wet regions should remain tied to coastal/moisture corridors"),
        _check(label=label, group="rain_shadow_placement",
               metric="rain_shadow_dry_response_corr",
               generated=row["rain_shadow_dry_response_corr"],
               operator=">=", threshold=0.70, severity="fail",
               message="rain-shadow objects should align with dry regional precipitation response"),
        _check(label=label, group="monsoon_map_readability",
               metric="monsoon_map_active_land_fraction_p50",
               generated=row["monsoon_map_active_land_fraction_p50"],
               operator=">=", threshold=0.10, severity="fail",
               message="monsoon corridor maps should not be too sparse to read"),
        _check(label=label, group="monsoon_map_readability",
               metric="monsoon_map_active_land_fraction_p50",
               generated=row["monsoon_map_active_land_fraction_p50"],
               operator="<=", threshold=0.25, severity="warn",
               message="monsoon corridor maps should not become broad annual wet zones"),
        _check(label=label, group="monsoon_map_readability",
               metric="monsoon_map_largest_component_share_p50",
               generated=row["monsoon_map_largest_component_share_p50"],
               operator=">=", threshold=0.25, severity="fail",
               message="monsoon corridor maps should include coherent connected belts"),
        _check(label=label, group="monsoon_map_readability",
               metric="monsoon_map_boundary_per_active_cell_p50",
               generated=row["monsoon_map_boundary_per_active_cell_p50"],
               operator="<=", threshold=2.60, severity="fail",
               message="monsoon corridor maps should not read as checkerboard texture"),
        _check(label=label, group="storm_track_map_readability",
               metric="storm_track_map_active_land_fraction_p50",
               generated=row["storm_track_map_active_land_fraction_p50"],
               operator=">=", threshold=0.10, severity="fail",
               message="storm-track maps should expose readable wet-coast belts"),
        _check(label=label, group="storm_track_map_readability",
               metric="storm_track_map_largest_component_share_p50",
               generated=row["storm_track_map_largest_component_share_p50"],
               operator=">=", threshold=0.18, severity="fail",
               message="storm-track maps should not fragment into isolated speckles"),
        _check(label=label, group="storm_track_map_readability",
               metric="storm_track_map_boundary_per_active_cell_p50",
               generated=row["storm_track_map_boundary_per_active_cell_p50"],
               operator="<=", threshold=2.40, severity="fail",
               message="storm-track maps should remain region-like, not pixel-like"),
        _check(label=label, group="response_map_readability",
               metric="wet_response_map_largest_component_share_p50",
               generated=row["wet_response_map_largest_component_share_p50"],
               operator=">=", threshold=0.18, severity="fail",
               message="wet regional-response maps should include coherent areas"),
        _check(label=label, group="response_map_readability",
               metric="dry_response_map_largest_component_share_p50",
               generated=row["dry_response_map_largest_component_share_p50"],
               operator=">=", threshold=0.16, severity="fail",
               message="dry regional-response maps should include coherent areas"),
    ]


def _waterworld_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(label=label, group="object_archive", metric="object_json_found",
               generated=row["object_json_found"], operator=">=", threshold=1.0,
               severity="fail", message="waterworld diagnostics should archive hydroclimate region objects"),
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="waterworld diagnostics should archive C4d arrays for placement checks"),
        _check(label=label, group="map_archive", metric="hydro_map_edges_found",
               generated=row["hydro_map_edges_found"], operator=">=", threshold=1.0,
               severity="fail", message="waterworld C4d arrays must support map-readability checks"),
        _check(label=label, group="waterworld_false_positive", metric="object_count",
               generated=row["object_count"], operator="<=", threshold=90.0,
               severity="fail", message="small waterworld islands should not produce many continental regions"),
        _check(label=label, group="waterworld_false_positive", metric="largest_area_fraction",
               generated=row["largest_area_fraction"], operator="<=", threshold=0.025,
               severity="fail", message="waterworld C4d regions should remain small"),
        _check(label=label, group="waterworld_false_positive", metric="largest_monsoon_area_fraction",
               generated=row["largest_monsoon_area_fraction"], operator="<=", threshold=0.018,
               severity="fail", message="waterworlds should not grow continent-scale monsoon rain regions"),
        _check(label=label, group="waterworld_false_positive",
               metric="monsoon_field_lat_shift_jja_minus_djf",
               generated=row["monsoon_field_lat_shift_jja_minus_djf"],
               operator="<=", threshold=35.0, severity="fail",
               message="waterworld island monsoon corridors should not mimic continental seasonal migration"),
        _check(label=label, group="waterworld_false_positive",
               metric="monsoon_map_active_world_fraction_p50",
               generated=row["monsoon_map_active_world_fraction_p50"],
               operator="<=", threshold=0.012, severity="fail",
               message="waterworld monsoon corridor maps should stay island-scale"),
        _check(label=label, group="waterworld_false_positive",
               metric="storm_track_map_active_world_fraction_p50",
               generated=row["storm_track_map_active_world_fraction_p50"],
               operator="<=", threshold=0.014, severity="fail",
               message="waterworld storm-track wet regions should stay island-scale"),
    ]


def _arid_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(label=label, group="object_archive", metric="object_json_found",
               generated=row["object_json_found"], operator=">=", threshold=1.0,
               severity="fail", message="arid diagnostics should archive hydroclimate region objects"),
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="arid diagnostics should archive C4d arrays for placement checks"),
        _check(label=label, group="map_archive", metric="hydro_map_edges_found",
               generated=row["hydro_map_edges_found"], operator=">=", threshold=1.0,
               severity="fail", message="arid C4d arrays must support map-readability checks"),
        _check(label=label, group="arid_dry_response", metric="dry_response_region_object_count",
               generated=row["dry_response_region_object_count"], operator=">=", threshold=4.0,
               severity="fail", message="arid worlds should expose dry regional-response objects"),
        _check(label=label, group="arid_false_monsoon_guard", metric="largest_monsoon_area_fraction",
               generated=row["largest_monsoon_area_fraction"], operator="<=", threshold=0.12,
               severity="fail", message="arid worlds should not become broadly monsoonal"),
        _check(label=label, group="arid_false_wet_guard", metric="largest_wet_response_area_fraction",
               generated=row["largest_wet_response_area_fraction"], operator="<=", threshold=0.16,
               severity="warn", message="arid wet-response regions should remain bounded"),
        _check(label=label, group="arid_rain_shadow_placement",
               metric="rain_shadow_dry_response_corr",
               generated=row["rain_shadow_dry_response_corr"],
               operator=">=", threshold=0.70, severity="fail",
               message="arid rain-shadow objects should align with dry regional response"),
        _check(label=label, group="arid_false_monsoon_guard",
               metric="monsoon_map_active_world_fraction_p50",
               generated=row["monsoon_map_active_world_fraction_p50"],
               operator="<=", threshold=0.08, severity="fail",
               message="arid monsoon corridor maps should not become broad wet belts"),
        _check(label=label, group="arid_false_wet_guard",
               metric="wet_response_map_active_world_fraction_p50",
               generated=row["wet_response_map_active_world_fraction_p50"],
               operator="<=", threshold=0.05, severity="warn",
               message="arid wet-response maps should remain spatially bounded"),
        _check(label=label, group="arid_dry_response",
               metric="dry_response_map_active_world_fraction_p50",
               generated=row["dry_response_map_active_world_fraction_p50"],
               operator=">=", threshold=0.04, severity="fail",
               message="arid maps should expose spatially meaningful dry-response regions"),
    ]


def _checks_for_row(row: dict[str, Any]) -> list[dict[str, Any]]:
    preset = str(row.get("preset", "")).lower()
    label = str(row.get("label", "")).lower()
    if "earthlike" in preset or "earthlike" in label:
        return _earthlike_checks(row)
    if "waterworld" in preset or "waterworld" in label:
        return _waterworld_checks(row)
    if "arid" in preset or "arid" in label:
        return _arid_checks(row)
    return []


def _write_csv(path: Path, rows: list[dict[str, Any]], keys: list[str]) -> Path:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def _render_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Earth Climate Hydroclimate-Region Gate",
        "",
        f"Schema: `{report['schema']}`",
        f"Verdict: **{report['verdict']}**",
        "",
        (
            f"Failures: `{report['failure_count']}`; warnings: "
            f"`{report['warning_count']}`; skipped: `{report['skipped_count']}`."
        ),
        "",
        "This gate evaluates C4d seasonal hydroclimate region objects.  It checks "
        "that Earthlike runs expose coherent seasonal monsoon, storm-track, "
        "rain-shadow, and wet/dry response regions, while waterworld and arid "
        "presets do not produce broad false-positive monsoon regions.  It also "
        "checks placement proxies from the archived C4d arrays: monsoon seasonal "
        "migration, storm-track latitude/coastal placement, and rain-shadow "
        "alignment with dry response.  The C4d array checks also include a first "
        "map-readability proxy: active-area fraction, largest connected belt "
        "share, and boundary roughness for corridor/response maps.",
        "",
        "## Checks",
        "",
    ]
    if report.get("contact_sheets"):
        lines.extend([
            "## Contact Sheets",
            "",
        ])
        for sheet in report["contact_sheets"]:
            lines.append(
                f"- `{sheet['label']}` `{sheet['preset']}`: `{sheet['path']}`"
            )
        lines.append("")
    for row in report["checks"]:
        status = "pass" if row["passed"] else row["severity"]
        lines.append(
            f"- `{status}` `{row['label']}` `{row['group']}` "
            f"`{row['metric']}` = `{row['generated']:.3f}` "
            f"{row['operator']} `{row['threshold']:.3f}`"
        )
    lines.append("")
    return "\n".join(lines)


def run_earth_climate_hydro_region_gate(
    config: EarthClimateHydroRegionGateConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(Path(config.terminal_summary_json).read_text())
    summary_rows = list(summary.get("summaries", []))
    generated = [_generated_row(row) for row in summary_rows]
    contact_sheets = (
        _render_contact_sheets(summary_rows, generated, outdir)
        if config.render_contact_sheets else []
    )
    checks = [
        check
        for row in generated
        for check in _checks_for_row(row)
    ]
    failures = [
        row for row in checks
        if not row["passed"] and row["severity"] == "fail"
    ]
    warnings = [
        row for row in checks
        if not row["passed"] and row["severity"] == "warn"
    ]
    skipped = [row for row in checks if row.get("skipped")]
    report = {
        "schema": SCHEMA,
        "terminal_summary_json": str(config.terminal_summary_json),
        "generated_metrics": generated,
        "checks": checks,
        "failures": failures,
        "warnings": warnings,
        "failure_count": len(failures),
        "warning_count": len(warnings),
        "skipped_count": len(skipped),
        "contact_sheets": contact_sheets,
        "contact_sheet_count": len(contact_sheets),
        "verdict": "fail" if failures else "pass",
    }
    metric_keys = sorted({
        key for row in generated for key in row.keys()
        if key not in {"kind_counts", "season_counts"}
    })
    check_keys = [
        "label", "group", "metric", "generated", "operator", "threshold",
        "severity", "passed", "skipped", "message",
    ]
    _write_csv(outdir / "earth_climate_hydro_region_metrics.csv",
               generated, metric_keys)
    _write_csv(outdir / "earth_climate_hydro_region_checks.csv",
               checks, check_keys)
    (outdir / "earth_climate_hydro_region_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default),
    )
    (outdir / "earth_climate_hydro_region_gate_report.md").write_text(
        _render_markdown(report),
    )
    return report
