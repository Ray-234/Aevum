"""F2 ocean-current and SST spatial-structure gate.

This gate follows the Earth-fitting order after the circulation-layout gate:
surface-current speed may be correctly scaled while the map is still too zonal
or insufficiently tied to coastlines, basins, and cold/warm-current SST
structure.  The checks here compare generated maps to broad OSCAR/OISST Earth
envelopes without requiring generated worlds to copy Earth's geography.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from aevum.diagnostics.earth_climate_circulation_layout_gate import _uv_to_tangent
from aevum.diagnostics.earth_climate_hydro_region_gate import (
    _array_path,
    _check,
    _json_default,
    _safe_float,
)


SCHEMA = "aevum.earth_climate_ocean_spatial_gate.v1"


@dataclass(frozen=True)
class EarthClimateOceanSpatialGateConfig:
    earth_reference_npz: Path
    terminal_summary_json: Path
    outdir: Path


def _xyz_from_lat_lon(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat_r = np.radians(np.asarray(lat, dtype=np.float64))
    lon_r = np.radians(np.asarray(lon, dtype=np.float64))
    return np.column_stack([
        np.cos(lat_r) * np.cos(lon_r),
        np.cos(lat_r) * np.sin(lon_r),
        np.sin(lat_r),
    ])


def _knn_edges(lat: np.ndarray, lon: np.ndarray, k: int = 6) -> np.ndarray:
    xyz = _xyz_from_lat_lon(lat, lon)
    k = max(1, min(int(k), xyz.shape[0] - 1))
    _, idx = cKDTree(xyz).query(xyz, k=k + 1)
    edges: set[tuple[int, int]] = set()
    for i, row in enumerate(np.asarray(idx, dtype=np.int64)):
        for j in row[1:]:
            a, b = sorted((int(i), int(j)))
            edges.add((a, b))
    return np.asarray(sorted(edges), dtype=np.int64)


def _ocean_distance_to_coast(
    lat: np.ndarray,
    lon: np.ndarray,
    land: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    land = np.asarray(land, dtype=bool)
    ocean = ~land
    n = land.size
    edges = _knn_edges(lat, lon, k=6)
    adjacent: list[list[int]] = [[] for _ in range(n)]
    coast = np.zeros(n, dtype=bool)
    for i, j in edges:
        adjacent[int(i)].append(int(j))
        adjacent[int(j)].append(int(i))
        if land[i] != land[j]:
            if ocean[i]:
                coast[i] = True
            if ocean[j]:
                coast[j] = True

    dist = np.full(n, -1, dtype=np.int32)
    frontier = [int(x) for x in np.where(coast)[0]]
    for cell in frontier:
        dist[cell] = 0
    head = 0
    while head < len(frontier):
        cell = frontier[head]
        head += 1
        for nb in adjacent[cell]:
            if ocean[nb] and dist[nb] < 0:
                dist[nb] = dist[cell] + 1
                frontier.append(nb)
    return dist, coast


def _percentile(values: np.ndarray, q: float, mask: np.ndarray | None = None) -> float:
    values = np.asarray(values, dtype=np.float64)
    if mask is not None:
        values = values[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    return float(np.percentile(values, q)) if values.size else float("nan")


def _area_fraction(area: np.ndarray, mask: np.ndarray, domain: np.ndarray) -> float:
    area = np.asarray(area, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    domain = np.asarray(domain, dtype=bool)
    denom = float(np.sum(area[domain & np.isfinite(area)]))
    if denom <= 1.0e-12:
        return float("nan")
    return float(np.sum(area[mask & domain & np.isfinite(area)]) / denom)


def _zonal_r2(lat: np.ndarray, values: np.ndarray, mask: np.ndarray) -> float:
    lat = np.asarray(lat, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    if int(np.count_nonzero(mask)) < 16:
        return float("nan")
    bins = np.arange(-90.0, 91.0, 10.0)
    predicted = np.full(values.shape, np.nan, dtype=np.float64)
    for lower, upper in zip(bins[:-1], bins[1:]):
        band = mask & (lat >= lower) & (lat < upper)
        if band.any():
            predicted[band] = float(np.mean(values[band]))
    valid = mask & np.isfinite(predicted)
    if int(np.count_nonzero(valid)) < 16:
        return float("nan")
    total = float(np.sum((values[valid] - float(np.mean(values[valid]))) ** 2))
    residual = float(np.sum((values[valid] - predicted[valid]) ** 2))
    return float(1.0 - residual / max(total, 1.0e-12))


def _zonal_residual_abs_p95(
    lat: np.ndarray,
    values: np.ndarray,
    mask: np.ndarray,
) -> float:
    lat = np.asarray(lat, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    if int(np.count_nonzero(mask)) < 16:
        return float("nan")
    bins = np.arange(-90.0, 91.0, 10.0)
    predicted = np.full(values.shape, np.nan, dtype=np.float64)
    for lower, upper in zip(bins[:-1], bins[1:]):
        band = mask & (lat >= lower) & (lat < upper)
        if band.any():
            predicted[band] = float(np.mean(values[band]))
    valid = mask & np.isfinite(predicted)
    if int(np.count_nonzero(valid)) < 16:
        return float("nan")
    return _percentile(np.abs(values - predicted), 95, valid)


def _coastal_current_sst_spread(
    lat: np.ndarray,
    speed: np.ndarray,
    sst_c: np.ndarray,
    ocean: np.ndarray,
    coastal_ocean: np.ndarray,
) -> tuple[float, float, float]:
    speed = np.asarray(speed, dtype=np.float64)
    sst_c = np.asarray(sst_c, dtype=np.float64)
    ocean = np.asarray(ocean, dtype=bool)
    coastal_ocean = np.asarray(coastal_ocean, dtype=bool)
    if int(np.count_nonzero(coastal_ocean)) < 8:
        return float("nan"), float("nan"), float("nan")
    anomaly = np.full(sst_c.shape, np.nan, dtype=np.float64)
    bins = np.arange(-90.0, 91.0, 10.0)
    for lower, upper in zip(bins[:-1], bins[1:]):
        band = ocean & (lat >= lower) & (lat < upper) & np.isfinite(sst_c)
        if band.any():
            anomaly[band] = sst_c[band] - float(np.median(sst_c[band]))
    threshold = max(
        _percentile(speed, 75, coastal_ocean),
        0.66 * _percentile(speed, 90, coastal_ocean),
    )
    swift_coast = coastal_ocean & np.isfinite(anomaly) & (speed >= threshold)
    if int(np.count_nonzero(swift_coast)) < 8:
        return float("nan"), float("nan"), float("nan")
    warm = _percentile(anomaly, 75, swift_coast)
    cold = _percentile(anomaly, 25, swift_coast)
    return warm, cold, float(warm - cold)


def _spatial_metrics(
    *,
    label: str,
    preset: str,
    seed: int,
    lat: np.ndarray,
    lon: np.ndarray,
    area: np.ndarray,
    land: np.ndarray,
    current_speed: np.ndarray,
    sst_c: np.ndarray,
    current_land_leak: float,
    heat_transport: np.ndarray | None = None,
) -> dict[str, Any]:
    land = np.asarray(land, dtype=bool)
    ocean = ~land
    current_speed = np.asarray(current_speed, dtype=np.float64)
    sst_c = np.asarray(sst_c, dtype=np.float64)
    dist, _ = _ocean_distance_to_coast(lat, lon, land)
    coastal_ocean = ocean & (dist >= 0) & (dist <= 3)
    far_ocean = ocean & (dist >= 8)
    valid_ocean = ocean & np.isfinite(current_speed)
    swift = valid_ocean & (
        current_speed >= _percentile(current_speed, 90, valid_ocean)
    )
    warm, cold, spread = _coastal_current_sst_spread(
        lat, current_speed, sst_c, ocean, coastal_ocean)
    heat_mean = 0.0
    heat_p95 = 0.0
    if heat_transport is not None:
        heat = np.asarray(heat_transport, dtype=np.float64)
        if heat.shape == area.shape and valid_ocean.any():
            heat_mean = float(np.average(heat[valid_ocean], weights=area[valid_ocean]))
            heat_p95 = _percentile(np.abs(heat), 95, valid_ocean)
    return {
        "label": label,
        "preset": preset,
        "seed": seed,
        "land_fraction": _area_fraction(area, land, np.ones_like(land, dtype=bool)),
        "ocean_fraction": _area_fraction(area, ocean, np.ones_like(land, dtype=bool)),
        "current_speed_p50_m_s": _percentile(current_speed, 50, valid_ocean),
        "current_speed_p90_m_s": _percentile(current_speed, 90, valid_ocean),
        "current_speed_zonal_r2": _zonal_r2(lat, current_speed, valid_ocean),
        "current_speed_zonal_residual_abs_p95_m_s": _zonal_residual_abs_p95(
            lat, current_speed, valid_ocean),
        "sst_zonal_r2": _zonal_r2(lat, sst_c, ocean & np.isfinite(sst_c)),
        "sst_zonal_residual_abs_p95_C": _zonal_residual_abs_p95(
            lat, sst_c, ocean & np.isfinite(sst_c)),
        "swift_current_near_coast_share": _area_fraction(area, swift & coastal_ocean, swift),
        "swift_current_far_ocean_share": _area_fraction(area, swift & far_ocean, swift),
        "coastal_current_speed_p90_ratio": (
            _percentile(current_speed, 90, coastal_ocean)
            / max(_percentile(current_speed, 90, valid_ocean), 1.0e-12)
            if coastal_ocean.any() else float("nan")
        ),
        "coastal_swift_sst_warm_anomaly_p75_C": warm,
        "coastal_swift_sst_cold_anomaly_p25_C": cold,
        "coastal_swift_sst_anomaly_spread_C": spread,
        "current_land_leak_fraction": current_land_leak,
        "current_heat_transport_ocean_mean_C": heat_mean,
        "current_heat_transport_abs_p95_C": heat_p95,
    }


def _earth_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        lat = np.asarray(z["lat"], dtype=np.float64)
        lon = np.asarray(z["lon"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        land = np.asarray(z["earth__elevation_m"], dtype=np.float64) >= 0.0
        speed = np.asarray(z["earth__annual_surface_current_speed_m_s"],
                           dtype=np.float64)
        if "earth__annual_sst_C" in z.files:
            sst = np.asarray(z["earth__annual_sst_C"], dtype=np.float64)
        else:
            sst = np.asarray(z["earth__annual_temperature_C"], dtype=np.float64)
    return _spatial_metrics(
        label="earth_reference",
        preset="earth_reference",
        seed=0,
        lat=lat,
        lon=lon,
        area=area,
        land=land,
        current_speed=speed,
        sst_c=sst,
        current_land_leak=0.0,
    )


def _generated_metrics(summary_row: dict[str, Any], earth: dict[str, Any]) -> dict[str, Any]:
    label = Path(str(summary_row.get("assets_dir", summary_row.get("label", "")))).name
    path = _array_path(summary_row)
    defaults: dict[str, Any] = {
        "label": label,
        "preset": str(summary_row.get("preset", "")),
        "seed": int(summary_row.get("seed", -1)),
        "arrays_found": 0.0,
        "arrays_path": "",
        "land_fraction": 0.0,
        "ocean_fraction": 0.0,
        "current_speed_p50_m_s": float("nan"),
        "current_speed_p90_m_s": float("nan"),
        "current_speed_p90_ratio_to_earth": float("nan"),
        "current_speed_zonal_r2": float("nan"),
        "current_speed_zonal_r2_delta_to_earth": float("nan"),
        "current_speed_zonal_residual_abs_p95_m_s": float("nan"),
        "current_speed_zonal_residual_abs_p95_ratio_to_earth": float("nan"),
        "sst_zonal_r2": float("nan"),
        "sst_zonal_r2_delta_to_earth": float("nan"),
        "sst_zonal_residual_abs_p95_C": float("nan"),
        "sst_zonal_residual_abs_p95_ratio_to_earth": float("nan"),
        "swift_current_near_coast_share": float("nan"),
        "swift_current_far_ocean_share": float("nan"),
        "coastal_current_speed_p90_ratio": float("nan"),
        "coastal_swift_sst_warm_anomaly_p75_C": float("nan"),
        "coastal_swift_sst_cold_anomaly_p25_C": float("nan"),
        "coastal_swift_sst_anomaly_spread_C": float("nan"),
        "current_land_leak_fraction": float("nan"),
        "current_heat_transport_ocean_mean_C": float("nan"),
        "current_heat_transport_abs_p95_C": float("nan"),
    }
    if path is None:
        return defaults
    defaults["arrays_found"] = 1.0
    defaults["arrays_path"] = str(path)
    with np.load(path, allow_pickle=False) as z:
        required = (
            "lat",
            "lon",
            "cell_area",
            "terrain__elevation_m",
            "sea_level_m",
            "ocean__currents",
            "climate__seasonal_sst",
        )
        if any(key not in z.files for key in required):
            return defaults
        lat = np.asarray(z["lat"], dtype=np.float64)
        lon = np.asarray(z["lon"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        land = np.asarray(z["terrain__elevation_m"], dtype=np.float64) >= sea
        currents = np.asarray(z["ocean__currents"], dtype=np.float64)
        speed = np.linalg.norm(currents, axis=1)
        sst = np.asarray(z["climate__seasonal_sst"], dtype=np.float64).mean(axis=0)
        sst_c = sst - 273.15
        heat_transport = (
            np.asarray(z["ocean__current_heat_transport"], dtype=np.float64)
            if "ocean__current_heat_transport" in z.files else None
        )
    land_speed = speed[land & np.isfinite(speed)]
    current_land_leak = float(np.mean(land_speed > 1.0e-6)) if land_speed.size else 0.0
    row = _spatial_metrics(
        label=label,
        preset=str(summary_row.get("preset", "")),
        seed=int(summary_row.get("seed", -1)),
        lat=lat,
        lon=lon,
        area=area,
        land=land,
        current_speed=speed,
        sst_c=sst_c,
        current_land_leak=current_land_leak,
        heat_transport=heat_transport,
    )
    earth_speed = _safe_float(earth.get("current_speed_p90_m_s"))
    row["arrays_found"] = 1.0
    row["arrays_path"] = str(path)
    row["current_speed_p90_ratio_to_earth"] = (
        row["current_speed_p90_m_s"] / earth_speed
        if earth_speed > 0.0 else float("nan")
    )
    row["current_speed_zonal_r2_delta_to_earth"] = (
        row["current_speed_zonal_r2"] - _safe_float(earth.get("current_speed_zonal_r2"))
    )
    earth_current_resid = _safe_float(
        earth.get("current_speed_zonal_residual_abs_p95_m_s"))
    row["current_speed_zonal_residual_abs_p95_ratio_to_earth"] = (
        row["current_speed_zonal_residual_abs_p95_m_s"] / earth_current_resid
        if earth_current_resid > 0.0 else float("nan")
    )
    row["sst_zonal_r2_delta_to_earth"] = (
        row["sst_zonal_r2"] - _safe_float(earth.get("sst_zonal_r2"))
    )
    earth_sst_resid = _safe_float(earth.get("sst_zonal_residual_abs_p95_C"))
    row["sst_zonal_residual_abs_p95_ratio_to_earth"] = (
        row["sst_zonal_residual_abs_p95_C"] / earth_sst_resid
        if earth_sst_resid > 0.0 else float("nan")
    )
    return row


def _checks_for_row(row: dict[str, Any], earth: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    preset = str(row.get("preset", "")).lower()
    checks = [
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="F2 ocean spatial gate requires archived arrays"),
        _check(label=label, group="current_speed",
               metric="current_speed_p90_ratio_to_earth",
               generated=row["current_speed_p90_ratio_to_earth"], operator=">=",
               threshold=0.45, severity="fail",
               message="surface-current p90 should remain in OSCAR order"),
        _check(label=label, group="current_speed",
               metric="current_speed_p90_ratio_to_earth",
               generated=row["current_speed_p90_ratio_to_earth"], operator="<=",
               threshold=1.80, severity="fail",
               message="surface-current p90 should not exceed OSCAR envelope"),
        _check(label=label, group="land_mask", metric="current_land_leak_fraction",
               generated=row["current_land_leak_fraction"], operator="<=",
               threshold=0.001, severity="fail",
               message="surface currents must remain ocean-confined"),
        _check(label=label, group="current_pattern",
               metric="current_speed_zonal_r2_delta_to_earth",
               generated=row["current_speed_zonal_r2_delta_to_earth"],
               operator="<=", threshold=0.14, severity="fail",
               message="current-speed map should not become mostly latitude bands"),
        _check(label=label, group="sst_pattern",
               metric="sst_zonal_r2_delta_to_earth",
               generated=row["sst_zonal_r2_delta_to_earth"],
               operator="<=", threshold=0.045, severity="fail",
               message="SST should retain non-zonal current/coast structure"),
        _check(label=label, group="sst_current_coupling",
               metric="coastal_swift_sst_anomaly_spread_C",
               generated=row["coastal_swift_sst_anomaly_spread_C"],
               operator=">=", threshold=1.0, severity="fail",
               message="swift coastal currents should expose warm/cold SST contrast"),
        _check(label=label, group="energy_conservation",
               metric="current_heat_transport_ocean_mean_C",
               generated=abs(_safe_float(row["current_heat_transport_ocean_mean_C"], 0.0)),
               operator="<=", threshold=0.20, severity="fail",
               message="current heat transport should not create net ocean energy"),
    ]
    if "earthlike" in preset or "earthlike" in label.lower():
        checks.extend([
            _check(label=label, group="boundary_current_structure",
                   metric="swift_current_far_ocean_share",
                   generated=row["swift_current_far_ocean_share"],
                   operator="<=", threshold=0.43, severity="fail",
                   message="earthlike strongest currents should not be dominated by remote open-ocean bands"),
            _check(label=label, group="boundary_current_structure",
                   metric="swift_current_near_coast_share",
                   generated=row["swift_current_near_coast_share"],
                   operator=">=", threshold=0.34, severity="fail",
                   message="earthlike strongest currents should include boundary-current/coastal structure"),
            _check(label=label, group="sst_pattern",
                   metric="sst_zonal_residual_abs_p95_ratio_to_earth",
                   generated=row["sst_zonal_residual_abs_p95_ratio_to_earth"],
                   operator=">=", threshold=0.58, severity="fail",
                   message="earthlike SST should have enough same-latitude structure, not only latitude bands"),
            _check(label=label, group="heat_transport_structure",
                   metric="current_heat_transport_abs_p95_C",
                   generated=row["current_heat_transport_abs_p95_C"],
                   operator=">=", threshold=0.95, severity="fail",
                   message="earthlike ocean heat transport should be strong enough to visibly break zonal SST bands"),
        ])
    if "waterworld" in preset or "waterworld" in label.lower():
        checks.append(
            _check(label=label, group="waterworld_current_structure",
                   metric="current_speed_zonal_r2_delta_to_earth",
                   generated=row["current_speed_zonal_r2_delta_to_earth"],
                   operator="<=", threshold=0.20, severity="fail",
                   message="waterworld currents may be open-ocean dominated but should not be purely zonal")
        )
    return checks


def _write_csv(path: Path, rows: list[dict[str, Any]], keys: list[str]) -> Path:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def _render_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Earth Climate Ocean Spatial Gate",
        "",
        f"Schema: `{report['schema']}`",
        f"Verdict: **{report['verdict']}**",
        "",
        (
            f"Failures: `{report['failure_count']}`; warnings: "
            f"`{report['warning_count']}`; skipped: `{report['skipped_count']}`."
        ),
        "",
        "This F2 gate evaluates spatial ocean-current and SST structure after "
        "current-speed magnitude has been calibrated against OSCAR.",
        "",
        "## Checks",
        "",
    ]
    for check in report["checks"]:
        mark = "PASS" if check["passed"] else "FAIL"
        if check["skipped"]:
            mark = "SKIP"
        lines.append(
            f"- `{mark}` `{check['label']}` `{check['group']}.{check['metric']}`: "
            f"{check['generated']:.4g} {check['operator']} {check['threshold']:.4g} "
            f"({check['severity']})"
        )
    lines.extend(["", "## Outputs", ""])
    lines.append("- `earth_climate_ocean_spatial_metrics.csv`")
    lines.append("- `earth_climate_ocean_spatial_checks.csv`")
    return "\n".join(lines) + "\n"


def run_earth_climate_ocean_spatial_gate(
    config: EarthClimateOceanSpatialGateConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    earth = _earth_metrics(Path(config.earth_reference_npz))
    summary = json.loads(Path(config.terminal_summary_json).read_text())
    rows = summary.get("summaries", summary if isinstance(summary, list) else [])
    generated = [
        _generated_metrics(row, earth)
        for row in rows
        if isinstance(row, dict)
    ]
    checks = [check for row in generated for check in _checks_for_row(row, earth)]
    failures = [
        check for check in checks
        if not check["passed"] and not check["skipped"] and check["severity"] == "fail"
    ]
    warnings = [
        check for check in checks
        if not check["passed"] and not check["skipped"] and check["severity"] == "warn"
    ]
    skipped = [check for check in checks if check["skipped"]]
    report = {
        "schema": SCHEMA,
        "earth_reference_npz": str(config.earth_reference_npz),
        "terminal_summary_json": str(config.terminal_summary_json),
        "earth_metrics": earth,
        "generated": generated,
        "checks": checks,
        "failures": failures,
        "warnings": warnings,
        "skipped": skipped,
        "failure_count": len(failures),
        "warning_count": len(warnings),
        "skipped_count": len(skipped),
        "verdict": "fail" if failures else "pass",
    }
    metric_keys = [
        "label", "preset", "seed", "arrays_found", "arrays_path",
        "land_fraction", "ocean_fraction",
        "current_speed_p50_m_s", "current_speed_p90_m_s",
        "current_speed_p90_ratio_to_earth",
        "current_speed_zonal_r2", "current_speed_zonal_r2_delta_to_earth",
        "current_speed_zonal_residual_abs_p95_m_s",
        "current_speed_zonal_residual_abs_p95_ratio_to_earth",
        "sst_zonal_r2", "sst_zonal_r2_delta_to_earth",
        "sst_zonal_residual_abs_p95_C",
        "sst_zonal_residual_abs_p95_ratio_to_earth",
        "swift_current_near_coast_share", "swift_current_far_ocean_share",
        "coastal_current_speed_p90_ratio",
        "coastal_swift_sst_warm_anomaly_p75_C",
        "coastal_swift_sst_cold_anomaly_p25_C",
        "coastal_swift_sst_anomaly_spread_C",
        "current_land_leak_fraction",
        "current_heat_transport_ocean_mean_C",
        "current_heat_transport_abs_p95_C",
    ]
    check_keys = [
        "label", "group", "metric", "generated", "operator", "threshold",
        "severity", "passed", "skipped", "message",
    ]
    _write_csv(outdir / "earth_climate_ocean_spatial_metrics.csv",
               generated, metric_keys)
    _write_csv(outdir / "earth_climate_ocean_spatial_checks.csv",
               checks, check_keys)
    (outdir / "earth_climate_ocean_spatial_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default))
    (outdir / "earth_climate_ocean_spatial_gate_report.md").write_text(
        _render_markdown(report))
    return report
