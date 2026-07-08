"""C4f conservative moisture-flow precipitation-response gate.

This gate evaluates the active precipitation-response layer that consumes C4e
moisture-flow networks.  Scalar climate, pattern, and biome gates still own the
Earth-fitting envelope; this gate only checks that the C4f response is archived,
finite, conservative, ocean-preserving, pathway/source-coupled, bounded, and
readable.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np

from aevum.diagnostics.earth_climate_hydro_region_gate import (
    _array_path,
    _check,
    _corr,
    _edges_from_latlon,
    _json_default,
    _prefixed,
    _safe_float,
    _write_csv,
)
from aevum.diagnostics.earth_climate_moisture_flow_gate import (
    _median,
    _seasonal_binary_map_metrics,
    _seasonal_domain_pctl,
    _seasonal_corr_median,
)


SCHEMA = "aevum.earth_climate_moisture_response_gate.v6"
SEASONS = ("DJF", "MAM", "JJA", "SON")


@dataclass(frozen=True)
class EarthClimateMoistureResponseGateConfig:
    terminal_summary_json: Path
    outdir: Path
    render_contact_sheets: bool = True


def _label(summary_row: dict[str, Any]) -> str:
    return Path(str(summary_row.get("assets_dir", summary_row.get("label", "")))).name


def _response_diag(summary_row: dict[str, Any]) -> dict[str, Any]:
    diag = summary_row.get("climate_step_diagnostics", {}) or {}
    if not isinstance(diag, dict):
        return {}
    response = diag.get("moisture_flow_precipitation_response", {}) or {}
    return response if isinstance(response, dict) else {}


def _seasonal_abs_dev_pctl(field: np.ndarray, domain: np.ndarray, q: float) -> float:
    field = np.asarray(field, dtype=np.float64)
    domain = np.asarray(domain, dtype=bool)
    if field.ndim != 2 or field.shape[0] != 4 or field.shape[1] != domain.size:
        return 0.0
    values: list[float] = []
    for season in range(4):
        vals = np.abs(field[season, domain] - 1.0)
        vals = vals[np.isfinite(vals)]
        values.append(float(np.percentile(vals, q)) if vals.size else 0.0)
    return _median(values)


def _seasonal_active_fraction_median(
    active4: np.ndarray,
    domain: np.ndarray,
    area: np.ndarray,
) -> float:
    active4 = np.asarray(active4, dtype=bool)
    domain = np.asarray(domain, dtype=bool)
    area = np.asarray(area, dtype=np.float64)
    if active4.ndim != 2 or active4.shape[0] != 4 or active4.shape[1] != domain.size:
        return 0.0
    denom = max(float(np.sum(area[domain])), 1.0e-12)
    return _median([
        float(np.sum(area[active4[season] & domain]) / denom)
        for season in range(4)
    ])


def _seasonal_within_fraction_median(
    active4: np.ndarray,
    domain4: np.ndarray,
    area: np.ndarray,
) -> float:
    active4 = np.asarray(active4, dtype=bool)
    domain4 = np.asarray(domain4, dtype=bool)
    area = np.asarray(area, dtype=np.float64)
    if active4.shape != domain4.shape or active4.ndim != 2 or active4.shape[0] != 4:
        return 0.0
    values: list[float] = []
    for season in range(4):
        domain = domain4[season]
        denom = float(np.sum(area[domain]))
        if denom <= 1.0e-12:
            values.append(0.0)
        else:
            values.append(float(np.sum(area[active4[season] & domain]) / denom))
    return _median(values)


def _budget_source_purity(
    source_basin: np.ndarray,
    budget_region: np.ndarray | None,
    land: np.ndarray,
    area: np.ndarray,
) -> tuple[float, float]:
    if budget_region is None:
        return 0.0, 0.0
    source_basin = np.asarray(source_basin, dtype=np.float64)
    budget_region = np.asarray(budget_region, dtype=np.float64)
    if (
        source_basin.shape != budget_region.shape
        or source_basin.ndim != 2
        or source_basin.shape[0] != 4
    ):
        return 0.0, 0.0
    purities: list[float] = []
    for season in range(4):
        for rid in [
            int(x) for x in np.unique(budget_region[season, land])
            if np.isfinite(x) and int(x) > 0
        ]:
            region = (
                land
                & (budget_region[season] == float(rid))
                & np.isfinite(source_basin[season])
                & (source_basin[season] >= 0.0)
            )
            total = float(np.sum(area[region]))
            if total <= 1.0e-12:
                continue
            basin_areas = [
                float(np.sum(area[region & (source_basin[season] == float(bid))]))
                for bid in [
                    int(x) for x in np.unique(source_basin[season, region])
                    if np.isfinite(x) and int(x) >= 0
                ]
            ]
            if basin_areas:
                purities.append(max(basin_areas) / max(sum(basin_areas), 1.0e-12))
    if not purities:
        return 0.0, 0.0
    return float(np.percentile(purities, 50)), float(np.percentile(purities, 10))


def _precipitation_response_region_archive_metrics(
    summary_row: dict[str, Any],
) -> dict[str, Any]:
    defaults = {
        "precip_region_archive_found": 0.0,
        "precip_region_object_count": 0.0,
        "precip_region_kind_count": 0.0,
        "precip_region_season_count": 0.0,
        "precip_region_wet_object_count": 0.0,
        "precip_region_dry_object_count": 0.0,
        "precip_region_largest_area_fraction": 0.0,
        "precip_region_mean_abs_anomaly_p50": 0.0,
        "precip_region_source_attribution_p50": 0.0,
        "precip_region_budget_attribution_p50": 0.0,
        "precip_region_wet_flow_attribution_p50": 0.0,
    }
    raw_path = summary_row.get("precipitation_response_regions_json")
    path = Path(str(raw_path)) if raw_path else None
    if path is None or not path.exists():
        assets = summary_row.get("assets_dir")
        candidate = Path(str(assets)) / "precipitation_response_regions.json" if assets else None
        path = candidate if candidate is not None and candidate.exists() else None
    if path is None:
        return defaults
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return defaults
    regions = payload.get("regions", [])
    if not isinstance(regions, list):
        regions = []
    rows = [row for row in regions if isinstance(row, dict)]
    if not rows:
        return defaults
    kinds = [str(row.get("kind", "")) for row in rows]
    seasons = [str(row.get("season", "")) for row in rows]
    area_fracs = [
        _safe_float(row.get("area_fraction", 0.0), 0.0)
        for row in rows
    ]
    anomalies = [
        _safe_float(row.get("mean_abs_response_anomaly", 0.0), 0.0)
        for row in rows
    ]
    source_attrs = [
        _safe_float(row.get("source_basin_attributed_fraction", 0.0), 0.0)
        for row in rows
    ]
    budget_attrs = [
        _safe_float(row.get("budget_region_attributed_fraction", 0.0), 0.0)
        for row in rows
    ]
    wet_flow_attrs = [
        _safe_float(row.get("flow_network_attributed_fraction", 0.0), 0.0)
        for row in rows
        if str(row.get("kind", "")) == "wet_precipitation_response_region"
    ]
    return {
        "precip_region_archive_found": 1.0,
        "precip_region_object_count": float(len(rows)),
        "precip_region_kind_count": float(len(set(kinds))),
        "precip_region_season_count": float(len(set(seasons))),
        "precip_region_wet_object_count": float(sum(
            1 for kind in kinds if kind == "wet_precipitation_response_region")),
        "precip_region_dry_object_count": float(sum(
            1 for kind in kinds if kind == "dry_precipitation_response_region")),
        "precip_region_largest_area_fraction": max(area_fracs) if area_fracs else 0.0,
        "precip_region_mean_abs_anomaly_p50": (
            float(np.percentile(anomalies, 50)) if anomalies else 0.0
        ),
        "precip_region_source_attribution_p50": (
            float(np.percentile(source_attrs, 50)) if source_attrs else 0.0
        ),
        "precip_region_budget_attribution_p50": (
            float(np.percentile(budget_attrs, 50)) if budget_attrs else 0.0
        ),
        "precip_region_wet_flow_attribution_p50": (
            float(np.percentile(wet_flow_attrs, 50)) if wet_flow_attrs else 0.0
        ),
    }


def _array_metrics(summary_row: dict[str, Any]) -> dict[str, Any]:
    defaults = {
        "arrays_found": 0.0,
        "arrays_path": "",
        "response_found": 0.0,
        "response_shape_ok": 0.0,
        "response_finite_fraction": 0.0,
        "map_edges_found": 0.0,
        "response_land_p05": 1.0,
        "response_land_p50": 1.0,
        "response_land_p95": 1.0,
        "response_land_span": 0.0,
        "response_ocean_abs_dev_p99": 0.0,
        "response_pathway_corr_median": 0.0,
        "response_support_corr_median": 0.0,
        "response_shadow_corr_median": 0.0,
        "response_high_land_fraction_p50": 0.0,
        "response_low_land_fraction_p50": 0.0,
        "response_high_map_component_count_p50": 0.0,
        "response_high_map_largest_component_share_p50": 0.0,
        "response_high_map_active_world_fraction_p50": 0.0,
        "response_high_map_active_domain_fraction_p50": 0.0,
        "response_high_map_boundary_per_active_cell_p50": 0.0,
        "response_low_map_component_count_p50": 0.0,
        "response_low_map_largest_component_share_p50": 0.0,
        "response_low_map_active_world_fraction_p50": 0.0,
        "response_low_map_active_domain_fraction_p50": 0.0,
        "response_low_map_boundary_per_active_cell_p50": 0.0,
        "budget_region_found": 0.0,
        "budget_region_shape_ok": 0.0,
        "budget_region_finite_fraction": 0.0,
        "budget_region_count_p50": 0.0,
        "source_basin_found": 0.0,
        "source_basin_shape_ok": 0.0,
        "source_basin_finite_fraction": 0.0,
        "source_basin_attributed_land_fraction_p50": 0.0,
        "source_basin_pathway_attributed_fraction_p50": 0.0,
        "source_basin_wet_response_attributed_fraction_p50": 0.0,
        "source_basin_network_attributed_fraction_p50": 0.0,
        "budget_source_purity_p50": 0.0,
        "budget_source_purity_p10": 0.0,
        "precip_region_id_found": 0.0,
        "precip_region_id_shape_ok": 0.0,
        "precip_region_id_finite_fraction": 0.0,
        "precip_region_id_count_p50": 0.0,
        "diagnostic_enabled": 0.0,
        "diagnostic_max_land_mean_delta_mm_yr": 0.0,
        "diagnostic_budget_base_region_count": 0.0,
        "diagnostic_budget_region_count_p50": 0.0,
        "diagnostic_budget_sector_split_count_p50": 0.0,
        "diagnostic_max_budget_region_mean_delta_mm_yr": 0.0,
        "diagnostic_response_land_p05": 1.0,
        "diagnostic_response_land_p95": 1.0,
    }
    diag = _response_diag(summary_row)
    defaults.update({
        "diagnostic_enabled": 1.0 if bool(diag.get("enabled", False)) else 0.0,
        "diagnostic_max_land_mean_delta_mm_yr": _safe_float(
            diag.get("max_land_mean_delta_mm_yr", 0.0), 0.0),
        "diagnostic_budget_base_region_count": _safe_float(
            diag.get("budget_base_region_count", 0.0), 0.0),
        "diagnostic_budget_region_count_p50": _safe_float(
            diag.get("budget_region_count_p50", 0.0), 0.0),
        "diagnostic_budget_sector_split_count_p50": _safe_float(
            diag.get("budget_sector_split_count_p50", 0.0), 0.0),
        "diagnostic_max_budget_region_mean_delta_mm_yr": _safe_float(
            diag.get("max_budget_region_mean_delta_mm_yr", 0.0), 0.0),
        "diagnostic_response_land_p05": _safe_float(
            diag.get("response_land_p05", 1.0), 1.0),
        "diagnostic_response_land_p95": _safe_float(
            diag.get("response_land_p95", 1.0), 1.0),
    })
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
            "climate__moisture_flow_precipitation_response",
            "atmosphere__moisture_flow_pathway",
            "climate__moisture_flow_network_id",
            "climate__monsoon_rainfall_corridor",
            "climate__storm_track_rainfall_corridor",
            "climate__rain_shadow_index",
            "climate__seasonal_precipitation",
        )
        if any(key not in z for key in required):
            out = dict(defaults)
            out["arrays_found"] = 1.0
            out["arrays_path"] = str(path)
            return out
        lat = np.asarray(z["lat"], dtype=np.float64)
        lon = np.asarray(z["lon"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        elev = np.asarray(z["terrain__elevation_m"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        response = np.asarray(
            z["climate__moisture_flow_precipitation_response"], dtype=np.float64)
        pathway = np.asarray(z["atmosphere__moisture_flow_pathway"], dtype=np.float64)
        network_id = np.asarray(z["climate__moisture_flow_network_id"], dtype=np.float64)
        monsoon = np.asarray(z["climate__monsoon_rainfall_corridor"], dtype=np.float64)
        storm = np.asarray(z["climate__storm_track_rainfall_corridor"], dtype=np.float64)
        shadow = np.asarray(z["climate__rain_shadow_index"], dtype=np.float64)
        source_basin = (
            np.asarray(z["atmosphere__moisture_source_basin_id"], dtype=np.float64)
            if "atmosphere__moisture_source_basin_id" in z else None
        )
        budget_region = (
            np.asarray(z["climate__moisture_budget_region_id"], dtype=np.float64)
            if "climate__moisture_budget_region_id" in z else None
        )
        precip_region_id = (
            np.asarray(z["climate__precipitation_response_region_id"], dtype=np.float64)
            if "climate__precipitation_response_region_id" in z else None
        )

    land = elev >= sea
    ocean = ~land
    expected = (4, lat.size)
    shape_ok = (
        response.shape == expected
        and pathway.shape == expected
        and network_id.shape == expected
        and monsoon.shape == expected
        and storm.shape == expected
        and shadow.shape == expected
    )
    out = dict(defaults)
    out.update({
        "arrays_found": 1.0,
        "arrays_path": str(path),
        "response_found": 1.0,
        "response_shape_ok": 1.0 if shape_ok else 0.0,
        "budget_region_found": 1.0 if budget_region is not None else 0.0,
        "source_basin_found": 1.0 if source_basin is not None else 0.0,
        "precip_region_id_found": 1.0 if precip_region_id is not None else 0.0,
    })
    if not shape_ok:
        return out

    finite = np.isfinite(response)
    out["response_finite_fraction"] = float(np.count_nonzero(finite) / response.size)
    if budget_region is not None:
        budget_shape_ok = budget_region.shape == expected
        out["budget_region_shape_ok"] = 1.0 if budget_shape_ok else 0.0
        if budget_shape_ok:
            budget_finite = np.isfinite(budget_region)
            out["budget_region_finite_fraction"] = float(
                np.count_nonzero(budget_finite) / budget_region.size)
            if land.any():
                counts = []
                for season in range(4):
                    ids = [
                        x for x in np.unique(budget_region[season, land])
                        if np.isfinite(x) and int(x) > 0
                    ]
                    counts.append(len(ids))
                out["budget_region_count_p50"] = _median(counts)
    source_shape_ok = source_basin is not None and source_basin.shape == expected
    out["source_basin_shape_ok"] = 1.0 if source_shape_ok else 0.0
    if source_shape_ok:
        finite_source = np.isfinite(source_basin)
        out["source_basin_finite_fraction"] = float(
            np.count_nonzero(finite_source) / source_basin.size)
    edges = _edges_from_latlon(lat, lon)
    out["map_edges_found"] = 1.0 if edges.size else 0.0
    support = np.maximum(monsoon, storm)
    high = response > 1.025
    low = response < 0.975
    if source_shape_ok:
        attributed = np.isfinite(source_basin) & (source_basin >= 0.0)
        pathway_active = np.zeros_like(attributed, dtype=bool)
        for season in range(4):
            if land.any():
                vals = pathway[season, land]
                threshold = max(0.10, float(np.percentile(vals, 75)))
                pathway_active[season] = land & (pathway[season] >= threshold)
        out["source_basin_attributed_land_fraction_p50"] = (
            _seasonal_active_fraction_median(attributed, land, area)
        )
        out["source_basin_pathway_attributed_fraction_p50"] = (
            _seasonal_within_fraction_median(attributed, pathway_active, area)
        )
        out["source_basin_wet_response_attributed_fraction_p50"] = (
            _seasonal_within_fraction_median(attributed, high & land[None, :], area)
        )
        out["source_basin_network_attributed_fraction_p50"] = (
            _seasonal_within_fraction_median(
                attributed, (network_id > 0.0) & land[None, :], area)
        )
        purity_p50, purity_p10 = _budget_source_purity(
            source_basin, budget_region, land, area)
        out["budget_source_purity_p50"] = purity_p50
        out["budget_source_purity_p10"] = purity_p10
    precip_region_shape_ok = (
        precip_region_id is not None and precip_region_id.shape == expected
    )
    out["precip_region_id_shape_ok"] = 1.0 if precip_region_shape_ok else 0.0
    if precip_region_shape_ok:
        finite_region = np.isfinite(precip_region_id)
        out["precip_region_id_finite_fraction"] = float(
            np.count_nonzero(finite_region) / precip_region_id.size)
        if land.any():
            counts = []
            for season in range(4):
                ids = [
                    x for x in np.unique(precip_region_id[season, land])
                    if np.isfinite(x) and int(x) > 0
                ]
                counts.append(len(ids))
            out["precip_region_id_count_p50"] = _median(counts)
    out.update({
        "response_land_p05": _seasonal_domain_pctl(response, land, 5.0),
        "response_land_p50": _seasonal_domain_pctl(response, land, 50.0),
        "response_land_p95": _seasonal_domain_pctl(response, land, 95.0),
        "response_ocean_abs_dev_p99": _seasonal_abs_dev_pctl(response, ocean, 99.0),
        "response_pathway_corr_median": _seasonal_corr_median(response, pathway, land),
        "response_support_corr_median": _seasonal_corr_median(response, support, land),
        "response_shadow_corr_median": _seasonal_corr_median(response, shadow, land),
        "response_high_land_fraction_p50": _seasonal_active_fraction_median(
            high, land, area),
        "response_low_land_fraction_p50": _seasonal_active_fraction_median(
            low, land, area),
    })
    out["response_land_span"] = float(out["response_land_p95"] - out["response_land_p05"])
    out.update(_prefixed(
        "response_high_map",
        _seasonal_binary_map_metrics(high, land, area, edges),
    ))
    out.update(_prefixed(
        "response_low_map",
        _seasonal_binary_map_metrics(low, land, area, edges),
    ))
    return out


def _generated_row(summary_row: dict[str, Any]) -> dict[str, Any]:
    row = {
        "label": _label(summary_row),
        "preset": str(summary_row.get("preset", "")),
        "seed": int(summary_row.get("seed", -1)),
        "assets_dir": str(summary_row.get("assets_dir", "")),
    }
    row.update(_array_metrics(summary_row))
    row.update(_precipitation_response_region_archive_metrics(summary_row))
    return row


def _earthlike_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="earthlike runs must archive C4f arrays"),
        _check(label=label, group="array_archive", metric="response_found",
               generated=row["response_found"], operator=">=", threshold=1.0,
               severity="fail", message="C4f response field must be archived"),
        _check(label=label, group="array_archive", metric="response_shape_ok",
               generated=row["response_shape_ok"], operator=">=", threshold=1.0,
               severity="fail", message="C4f response must have four seasonal maps"),
        _check(label=label, group="array_archive", metric="response_finite_fraction",
               generated=row["response_finite_fraction"], operator=">=", threshold=1.0,
               severity="fail", message="C4f response must be finite"),
        _check(label=label, group="conservation", metric="diagnostic_enabled",
               generated=row["diagnostic_enabled"], operator=">=", threshold=1.0,
               severity="fail", message="C4f diagnostic should report active response"),
        _check(label=label, group="conservation",
               metric="diagnostic_max_land_mean_delta_mm_yr",
               generated=row["diagnostic_max_land_mean_delta_mm_yr"],
               operator="<=", threshold=1.0e-5, severity="fail",
               message="C4f must preserve seasonal land precipitation mean"),
        _check(label=label, group="local_budget", metric="budget_region_found",
               generated=row["budget_region_found"], operator=">=", threshold=1.0,
               severity="fail", message="C4g/C4h moisture budget region id must be archived"),
        _check(label=label, group="local_budget", metric="budget_region_shape_ok",
               generated=row["budget_region_shape_ok"], operator=">=", threshold=1.0,
               severity="fail", message="C4g/C4h budget regions must have four seasonal maps"),
        _check(label=label, group="local_budget", metric="budget_region_finite_fraction",
               generated=row["budget_region_finite_fraction"], operator=">=", threshold=1.0,
               severity="fail", message="C4g/C4h budget regions must be finite"),
        _check(label=label, group="local_budget", metric="budget_region_count_p50",
               generated=row["budget_region_count_p50"], operator=">=", threshold=1.0,
               severity="fail", message="C4g/C4h should expose at least one land budget region"),
        _check(label=label, group="local_budget",
               metric="diagnostic_budget_sector_split_count_p50",
               generated=row["diagnostic_budget_sector_split_count_p50"],
               operator=">=", threshold=1.0, severity="fail",
               message="C4h earthlike runs should split at least one large moisture-network sector"),
        _check(label=label, group="local_budget",
               metric="diagnostic_max_budget_region_mean_delta_mm_yr",
               generated=row["diagnostic_max_budget_region_mean_delta_mm_yr"],
               operator="<=", threshold=1.0e-5, severity="fail",
               message="C4g/C4h must preserve seasonal precipitation mean inside each budget region"),
        _check(label=label, group="source_basin", metric="source_basin_found",
               generated=row["source_basin_found"], operator=">=", threshold=1.0,
               severity="fail", message="C4i source-basin field must be archived"),
        _check(label=label, group="source_basin", metric="source_basin_shape_ok",
               generated=row["source_basin_shape_ok"], operator=">=", threshold=1.0,
               severity="fail", message="C4i source-basin field must have four seasonal maps"),
        _check(label=label, group="source_basin", metric="source_basin_finite_fraction",
               generated=row["source_basin_finite_fraction"], operator=">=", threshold=1.0,
               severity="fail", message="C4i source-basin field must be finite"),
        _check(label=label, group="source_basin",
               metric="source_basin_attributed_land_fraction_p50",
               generated=row["source_basin_attributed_land_fraction_p50"],
               operator=">=", threshold=0.30, severity="fail",
               message="C4i should attribute a meaningful share of land moisture pathways to source basins"),
        _check(label=label, group="source_basin",
               metric="source_basin_pathway_attributed_fraction_p50",
               generated=row["source_basin_pathway_attributed_fraction_p50"],
               operator=">=", threshold=0.90, severity="fail",
               message="active moisture pathways should carry source-basin labels"),
        _check(label=label, group="source_basin",
               metric="source_basin_wet_response_attributed_fraction_p50",
               generated=row["source_basin_wet_response_attributed_fraction_p50"],
               operator=">=", threshold=0.90, severity="fail",
               message="wet response corridors should carry source-basin labels"),
        _check(label=label, group="source_basin",
               metric="source_basin_network_attributed_fraction_p50",
               generated=row["source_basin_network_attributed_fraction_p50"],
               operator=">=", threshold=0.90, severity="fail",
               message="moisture-flow networks should carry source-basin labels"),
        _check(label=label, group="source_basin",
               metric="budget_source_purity_p50",
               generated=row["budget_source_purity_p50"],
               operator=">=", threshold=0.80, severity="fail",
               message="C4h/C4i budget sectors should be dominated by one source basin"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_id_found",
               generated=row["precip_region_id_found"], operator=">=", threshold=1.0,
               severity="fail", message="C4j precipitation-response region id must be archived"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_id_shape_ok",
               generated=row["precip_region_id_shape_ok"], operator=">=", threshold=1.0,
               severity="fail", message="C4j region ids must have four seasonal maps"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_id_finite_fraction",
               generated=row["precip_region_id_finite_fraction"],
               operator=">=", threshold=1.0, severity="fail",
               message="C4j region ids must be finite"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_id_count_p50",
               generated=row["precip_region_id_count_p50"], operator=">=", threshold=1.0,
               severity="fail", message="C4j should mark at least one seasonal response region"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_archive_found",
               generated=row["precip_region_archive_found"],
               operator=">=", threshold=1.0, severity="fail",
               message="C4j precipitation-response region objects must be archived"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_kind_count",
               generated=row["precip_region_kind_count"], operator=">=", threshold=2.0,
               severity="fail", message="C4j should expose both wet and dry response object kinds"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_wet_object_count",
               generated=row["precip_region_wet_object_count"], operator=">=", threshold=1.0,
               severity="fail", message="C4j should expose wet response objects"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_dry_object_count",
               generated=row["precip_region_dry_object_count"], operator=">=", threshold=1.0,
               severity="fail", message="C4j should expose dry response objects"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_source_attribution_p50",
               generated=row["precip_region_source_attribution_p50"],
               operator=">=", threshold=0.70, severity="fail",
               message="C4j response objects should carry source-basin attribution"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_budget_attribution_p50",
               generated=row["precip_region_budget_attribution_p50"],
               operator=">=", threshold=0.85, severity="fail",
               message="C4j response objects should carry local budget-region attribution"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_wet_flow_attribution_p50",
               generated=row["precip_region_wet_flow_attribution_p50"],
               operator=">=", threshold=0.50, severity="fail",
               message="C4j wet response objects should stay tied to moisture-flow networks"),
        _check(label=label, group="ocean_preservation",
               metric="response_ocean_abs_dev_p99",
               generated=row["response_ocean_abs_dev_p99"], operator="<=", threshold=1.0e-9,
               severity="fail", message="C4f must not alter ocean precipitation"),
        _check(label=label, group="response_strength", metric="response_land_p95",
               generated=row["response_land_p95"], operator=">=", threshold=1.035,
               severity="fail", message="earthlike response should visibly wet routed corridors"),
        _check(label=label, group="response_strength", metric="response_land_p05",
               generated=row["response_land_p05"], operator="<=", threshold=0.920,
               severity="fail", message="earthlike response should visibly dry donor regions"),
        _check(label=label, group="response_bounds", metric="response_land_p95",
               generated=row["response_land_p95"], operator="<=", threshold=1.18,
               severity="fail", message="C4f wetting should remain bounded"),
        _check(label=label, group="response_bounds", metric="response_land_p05",
               generated=row["response_land_p05"], operator=">=", threshold=0.70,
               severity="fail", message="C4f drying should remain bounded"),
        _check(label=label, group="pathway_coupling",
               metric="response_pathway_corr_median",
               generated=row["response_pathway_corr_median"], operator=">=", threshold=0.35,
               severity="fail", message="C4f wetting should align with moisture pathways"),
        _check(label=label, group="support_coupling",
               metric="response_support_corr_median",
               generated=row["response_support_corr_median"], operator=">=", threshold=0.10,
               severity="fail", message="C4f wetting should align with monsoon/storm support"),
        _check(label=label, group="rain_shadow_coupling",
               metric="response_shadow_corr_median",
               generated=row["response_shadow_corr_median"], operator="<=", threshold=0.15,
               severity="fail", message="C4f response should not wet rain-shadow regions"),
        _check(label=label, group="map_readability",
               metric="response_high_land_fraction_p50",
               generated=row["response_high_land_fraction_p50"], operator=">=", threshold=0.04,
               severity="fail", message="wet response should not be invisible"),
        _check(label=label, group="map_readability",
               metric="response_high_map_largest_component_share_p50",
               generated=row["response_high_map_largest_component_share_p50"],
               operator=">=", threshold=0.10, severity="fail",
               message="wet response should form readable corridors"),
        _check(label=label, group="map_readability",
               metric="response_high_map_boundary_per_active_cell_p50",
               generated=row["response_high_map_boundary_per_active_cell_p50"],
               operator="<=", threshold=2.90, severity="fail",
               message="wet response should not be checkerboard texture"),
    ]


def _arid_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    checks = [
        check for check in _earthlike_checks(row)
        if not (
            check["group"] == "source_basin"
            and check["metric"] == "budget_source_purity_p50"
        )
    ]
    label = str(row["label"])
    checks.extend([
        _check(label=label, group="source_basin",
               metric="budget_source_purity_p50",
               generated=row["budget_source_purity_p50"],
               operator=">=", threshold=0.55, severity="fail",
               message="arid active budget sectors should retain coarse source-basin coherence"),
        _check(label=label, group="arid_response_bounds", metric="response_land_p95",
               generated=row["response_land_p95"], operator="<=", threshold=1.14,
               severity="fail", message="arid response should not create extreme wet cores"),
        _check(label=label, group="arid_response_bounds",
               metric="response_high_land_fraction_p50",
               generated=row["response_high_land_fraction_p50"],
               operator="<=", threshold=0.32, severity="fail",
               message="arid wet response should remain bounded"),
        _check(label=label, group="arid_response_bounds",
               metric="response_low_land_fraction_p50",
               generated=row["response_low_land_fraction_p50"],
               operator="<=", threshold=0.95, severity="fail",
               message="arid drying should not cover almost all land"),
    ])
    return checks


def _waterworld_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="waterworld diagnostics must archive C4f arrays"),
        _check(label=label, group="array_archive", metric="response_found",
               generated=row["response_found"], operator=">=", threshold=1.0,
               severity="fail", message="C4f response field must be archived"),
        _check(label=label, group="array_archive", metric="response_finite_fraction",
               generated=row["response_finite_fraction"], operator=">=", threshold=1.0,
               severity="fail", message="C4f response must be finite"),
        _check(label=label, group="conservation",
               metric="diagnostic_max_land_mean_delta_mm_yr",
               generated=row["diagnostic_max_land_mean_delta_mm_yr"],
               operator="<=", threshold=1.0e-5, severity="fail",
               message="C4f must preserve seasonal land precipitation mean"),
        _check(label=label, group="local_budget", metric="budget_region_found",
               generated=row["budget_region_found"], operator=">=", threshold=1.0,
               severity="fail", message="C4g/C4h moisture budget region id must be archived"),
        _check(label=label, group="local_budget", metric="budget_region_shape_ok",
               generated=row["budget_region_shape_ok"], operator=">=", threshold=1.0,
               severity="fail", message="C4g/C4h budget regions must have four seasonal maps"),
        _check(label=label, group="local_budget", metric="budget_region_finite_fraction",
               generated=row["budget_region_finite_fraction"], operator=">=", threshold=1.0,
               severity="fail", message="C4g/C4h budget regions must be finite"),
        _check(label=label, group="local_budget",
               metric="diagnostic_max_budget_region_mean_delta_mm_yr",
               generated=row["diagnostic_max_budget_region_mean_delta_mm_yr"],
               operator="<=", threshold=1.0e-5, severity="fail",
               message="C4g/C4h must preserve seasonal precipitation mean inside each budget region"),
        _check(label=label, group="source_basin", metric="source_basin_found",
               generated=row["source_basin_found"], operator=">=", threshold=1.0,
               severity="fail", message="C4i source-basin field must be archived"),
        _check(label=label, group="source_basin", metric="source_basin_shape_ok",
               generated=row["source_basin_shape_ok"], operator=">=", threshold=1.0,
               severity="fail", message="C4i source-basin field must have four seasonal maps"),
        _check(label=label, group="source_basin", metric="source_basin_finite_fraction",
               generated=row["source_basin_finite_fraction"], operator=">=", threshold=1.0,
               severity="fail", message="C4i source-basin field must be finite"),
        _check(label=label, group="source_basin",
               metric="source_basin_pathway_attributed_fraction_p50",
               generated=row["source_basin_pathway_attributed_fraction_p50"],
               operator=">=", threshold=0.90, severity="fail",
               message="active island moisture pathways should carry source-basin labels"),
        _check(label=label, group="source_basin",
               metric="budget_source_purity_p50",
               generated=row["budget_source_purity_p50"],
               operator=">=", threshold=0.80, severity="fail",
               message="waterworld island budget sectors should be source-basin coherent"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_id_found",
               generated=row["precip_region_id_found"], operator=">=", threshold=1.0,
               severity="fail", message="C4j precipitation-response region id must be archived"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_id_shape_ok",
               generated=row["precip_region_id_shape_ok"], operator=">=", threshold=1.0,
               severity="fail", message="C4j region ids must have four seasonal maps"),
        _check(label=label, group="precipitation_object_continuity",
               metric="precip_region_id_finite_fraction",
               generated=row["precip_region_id_finite_fraction"],
               operator=">=", threshold=1.0, severity="fail",
               message="C4j region ids must be finite"),
        _check(label=label, group="ocean_preservation",
               metric="response_ocean_abs_dev_p99",
               generated=row["response_ocean_abs_dev_p99"], operator="<=", threshold=1.0e-9,
               severity="fail", message="C4f must not alter ocean precipitation"),
        _check(label=label, group="waterworld_false_positive",
               metric="response_land_p95",
               generated=row["response_land_p95"], operator="<=", threshold=1.06,
               severity="fail", message="waterworlds should not grow strong continental wet response"),
        _check(label=label, group="waterworld_false_positive",
               metric="response_land_p05",
               generated=row["response_land_p05"], operator=">=", threshold=0.82,
               severity="fail", message="waterworlds should not grow strong continental dry response"),
        _check(label=label, group="waterworld_false_positive",
               metric="response_high_map_active_world_fraction_p50",
               generated=row["response_high_map_active_world_fraction_p50"],
               operator="<=", threshold=0.020, severity="fail",
               message="waterworld wet response should remain island-scale"),
        _check(label=label, group="waterworld_false_positive",
               metric="response_low_map_active_world_fraction_p50",
               generated=row["response_low_map_active_world_fraction_p50"],
               operator="<=", threshold=0.030, severity="fail",
               message="waterworld dry response should remain island-scale"),
        _check(label=label, group="pathway_coupling",
               metric="response_pathway_corr_median",
               generated=row["response_pathway_corr_median"], operator=">=", threshold=0.20,
               severity="fail", message="waterworld island response should still follow pathways"),
    ]


def _checks_for_row(row: dict[str, Any]) -> list[dict[str, Any]]:
    preset = str(row.get("preset", "")).lower()
    label = str(row.get("label", "")).lower()
    if "waterworld" in preset or "waterworld" in label:
        return _waterworld_checks(row)
    if "arid" in preset or "arid" in label:
        return _arid_checks(row)
    if "earthlike" in preset or "earthlike" in label:
        return _earthlike_checks(row)
    return []


def _contact_sheet_arrays(summary_row: dict[str, Any]) -> dict[str, np.ndarray] | None:
    path = _array_path(summary_row)
    if path is None:
        return None
    with np.load(path, allow_pickle=False) as z:
        required = (
            "lat",
            "lon",
            "terrain__elevation_m",
            "sea_level_m",
            "climate__seasonal_precipitation",
            "atmosphere__moisture_flow_pathway",
            "atmosphere__moisture_source_basin_id",
            "climate__moisture_flow_precipitation_response",
            "climate__moisture_budget_region_id",
            "climate__precipitation_response_region_id",
            "climate__moisture_flow_network_id",
            "climate__monsoon_rainfall_corridor",
            "climate__storm_track_rainfall_corridor",
            "climate__rain_shadow_index",
        )
        if any(key not in z for key in required):
            return None
        arrays = {key: np.asarray(z[key]) for key in required}
    sea = float(np.asarray(arrays["sea_level_m"], dtype=np.float64).ravel()[0])
    arrays["land"] = np.asarray(arrays["terrain__elevation_m"], dtype=np.float64) >= sea
    return arrays


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
    from matplotlib.colors import TwoSlopeNorm

    from aevum.render import PRECIP_CMAP

    label = str(generated_row.get("label", "world"))
    sheets_dir = outdir / "contact_sheets"
    sheets_dir.mkdir(parents=True, exist_ok=True)
    path = sheets_dir / f"{label}_moisture_response_contact_sheet.png"

    lat = np.asarray(arrays["lat"], dtype=np.float64)
    lon = np.asarray(arrays["lon"], dtype=np.float64)
    land = np.asarray(arrays["land"], dtype=bool)
    precip = np.asarray(arrays["climate__seasonal_precipitation"], dtype=np.float64)
    pathway = np.asarray(arrays["atmosphere__moisture_flow_pathway"], dtype=np.float64)
    source_basin_id = np.asarray(
        arrays["atmosphere__moisture_source_basin_id"], dtype=np.float64)
    response = np.asarray(
        arrays["climate__moisture_flow_precipitation_response"], dtype=np.float64)
    budget_region = np.asarray(
        arrays["climate__moisture_budget_region_id"], dtype=np.float64)
    precip_region_id = np.asarray(
        arrays["climate__precipitation_response_region_id"], dtype=np.float64)
    network_id = np.asarray(arrays["climate__moisture_flow_network_id"], dtype=np.float64)
    support = np.maximum(
        np.asarray(arrays["climate__monsoon_rainfall_corridor"], dtype=np.float64),
        np.asarray(arrays["climate__storm_track_rainfall_corridor"], dtype=np.float64),
    )
    shadow = np.asarray(arrays["climate__rain_shadow_index"], dtype=np.float64)

    n = max(int(lat.size), 1)
    marker_size = float(np.clip(28000.0 / n, 0.45, 10.0))
    base_color = np.where(land, "#dfdac8", "#d9edf4")
    fig, axes = plt.subplots(9, 4, figsize=(14.5, 18.0), constrained_layout=True)
    fig.suptitle(f"C4f/C4j moisture-flow precipitation response: {label}", fontsize=13)

    def base(ax):
        ax.scatter(lon, lat, c=base_color, s=marker_size, linewidths=0,
                   rasterized=True)
        ax.set_xlim(-180.0, 180.0)
        ax.set_ylim(-90.0, 90.0)
        ax.set_xticks([])
        ax.set_yticks([])

    def panel(ax, values, title, cmap, vmin=None, vmax=None, mask=None, norm=None):
        base(ax)
        vals = np.asarray(values, dtype=np.float64)
        if mask is not None:
            vals = np.where(mask, vals, np.nan)
        finite = np.isfinite(vals)
        if finite.any():
            image = ax.scatter(lon[finite], lat[finite], c=vals[finite],
                               s=marker_size * 1.25, linewidths=0, cmap=cmap,
                               vmin=None if norm is not None else vmin,
                               vmax=None if norm is not None else vmax,
                               norm=norm, rasterized=True)
        else:
            image = ax.scatter([], [], c=[], cmap=cmap, vmin=vmin, vmax=vmax,
                               norm=norm)
        ax.set_title(title, fontsize=8)
        return image

    precip_land = precip[:, land] if land.any() else precip
    precip_vmax = max(float(np.nanpercentile(precip_land, 98)), 100.0)
    pathway_vmax = max(float(np.nanpercentile(pathway, 98)), 0.10)
    support_vmax = max(float(np.nanpercentile(support, 98)), 0.10)
    shadow_vmax = max(float(np.nanpercentile(shadow, 98)), 0.10)
    network_positive = network_id[network_id > 0]
    network_vmax = max(float(np.nanmax(network_positive)) if network_positive.size else 1.0, 1.0)
    budget_positive = budget_region[budget_region > 0]
    budget_vmax = max(float(np.nanmax(budget_positive)) if budget_positive.size else 1.0, 1.0)
    precip_region_positive = precip_region_id[precip_region_id > 0]
    precip_region_vmax = max(
        float(np.nanmax(precip_region_positive)) if precip_region_positive.size else 1.0,
        1.0,
    )
    source_basin_positive = source_basin_id[source_basin_id >= 0.0]
    source_basin_vmax = max(
        float(np.nanmax(source_basin_positive)) if source_basin_positive.size else 1.0,
        1.0,
    )
    span = max(float(np.nanpercentile(np.abs(response - 1.0), 98)), 0.05)
    response_norm = TwoSlopeNorm(vmin=max(0.0, 1.0 - span), vcenter=1.0,
                                 vmax=1.0 + span)

    row_maps = [
        ("seasonal precipitation", precip, None, PRECIP_CMAP, 0.0, precip_vmax, None),
        ("land moisture pathway", pathway, land, "BuPu", 0.0, pathway_vmax, None),
        ("source ocean basin id", source_basin_id, source_basin_id >= 0.0,
         "tab20", 0.0, source_basin_vmax, None),
        ("C4f precip response", response, land, "coolwarm", None, None, response_norm),
        ("C4h budget region id", budget_region, budget_region > 0.0,
         "tab20", 0.0, budget_vmax, None),
        ("C4j response region id", precip_region_id, precip_region_id > 0.0,
         "tab20", 0.0, precip_region_vmax, None),
        ("flow network id", network_id, network_id > 0.0, "tab20", 0.0, network_vmax, None),
        ("monsoon/storm support", support, land, "PuBuGn", 0.0, support_vmax, None),
        ("rain-shadow index", shadow, land, "YlOrBr", 0.0, shadow_vmax, None),
    ]
    images = []
    for row_idx, (row_title, field, mask, cmap, vmin, vmax, norm) in enumerate(row_maps):
        image = None
        for season_idx, season in enumerate(SEASONS):
            local_mask = (
                None if mask is None
                else mask[season_idx] if getattr(mask, "ndim", 1) == 2
                else mask
            )
            image = panel(
                axes[row_idx, season_idx],
                field[season_idx],
                f"{season} {row_title}",
                cmap,
                vmin,
                vmax,
                local_mask,
                norm,
            )
        images.append(image)
    for row_idx, image in enumerate(images):
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
    by_assets = {str(row.get("assets_dir", "")): row for row in generated_rows}
    for summary_row in summary_rows:
        generated_row = by_assets.get(str(summary_row.get("assets_dir", "")))
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
        (outdir / "earth_climate_moisture_response_contact_sheets.json").write_text(
            json.dumps({
                "schema": "aevum.earth_climate_moisture_response_contact_sheets.v1",
                "sheet_count": int(len(sheets)),
                "sheets": sheets,
            }, indent=2, default=_json_default),
        )
    return sheets


def _render_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Earth Climate Moisture-Flow Precipitation-Response Gate",
        "",
        f"Schema: `{report['schema']}`",
        f"Verdict: **{report['verdict']}**",
        "",
        (
            f"Failures: `{report['failure_count']}`; warnings: "
            f"`{report['warning_count']}`; skipped: `{report['skipped_count']}`."
        ),
        "",
        "This gate evaluates the C4f active precipitation response and the C4g/C4h/C4i "
        "local moisture-budget region layer.  It checks that the response is "
        "archived, conservative globally and locally, ocean-preserving, bounded, "
        "coupled to moisture pathways/source basins, and visibly organized without letting "
        "waterworlds develop false continent-scale monsoons.",
        "",
    ]
    if report.get("contact_sheets"):
        lines.extend(["## Contact Sheets", ""])
        for sheet in report["contact_sheets"]:
            lines.append(
                f"- `{sheet['label']}` `{sheet['preset']}`: `{sheet['path']}`")
        lines.append("")
    lines.extend(["## Checks", ""])
    for row in report["checks"]:
        status = "pass" if row["passed"] else row["severity"]
        lines.append(
            f"- `{status}` `{row['label']}` `{row['group']}` "
            f"`{row['metric']}` = `{row['generated']:.6f}` "
            f"{row['operator']} `{row['threshold']:.6f}`"
        )
    lines.append("")
    return "\n".join(lines)


def run_earth_climate_moisture_response_gate(
    config: EarthClimateMoistureResponseGateConfig,
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
    checks = [check for row in generated for check in _checks_for_row(row)]
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
    metric_keys = sorted({key for row in generated for key in row.keys()})
    check_keys = [
        "label", "group", "metric", "generated", "operator", "threshold",
        "severity", "passed", "skipped", "message",
    ]
    _write_csv(outdir / "earth_climate_moisture_response_metrics.csv",
               generated, metric_keys)
    _write_csv(outdir / "earth_climate_moisture_response_checks.csv",
               checks, check_keys)
    (outdir / "earth_climate_moisture_response_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default),
    )
    (outdir / "earth_climate_moisture_response_gate_report.md").write_text(
        _render_markdown(report),
    )
    return report
