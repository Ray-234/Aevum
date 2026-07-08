"""F2/F3 circulation layout gate against broad Earth envelopes.

This gate evaluates the atmospheric-wind and surface-current layers before
later moisture, precipitation, Koppen, and biome tuning.  It intentionally
checks only coarse Earth-calibrated envelopes and geography-coupling proxies:
generated worlds should not copy Earth's geography, but their winds and
currents should be in Earthlike magnitudes and respond measurably to coastlines,
continents, barriers, and ocean basins.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from aevum.diagnostics.earth_climate_hydro_region_gate import (
    _array_path,
    _check,
    _json_default,
    _safe_float,
)


SCHEMA = "aevum.earth_climate_circulation_layout_gate.v1"


@dataclass(frozen=True)
class EarthClimateCirculationLayoutGateConfig:
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


def _tangent_basis(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xyz = _xyz_from_lat_lon(lat, lon)
    z_axis = np.array([0.0, 0.0, 1.0])
    east = np.cross(z_axis, xyz)
    norm = np.linalg.norm(east, axis=1, keepdims=True)
    east = np.where(norm > 1.0e-9, east / np.maximum(norm, 1.0e-9),
                    np.array([1.0, 0.0, 0.0]))
    north = np.cross(xyz, east)
    return xyz, east, north


def _uv_to_tangent(lat: np.ndarray, lon: np.ndarray, uv: np.ndarray) -> np.ndarray:
    _, east, north = _tangent_basis(lat, lon)
    uv = np.asarray(uv, dtype=np.float64)
    if uv.ndim == 2:
        return uv[:, 0, None] * east + uv[:, 1, None] * north
    return uv[:, :, 0, None] * east[None, :, :] + uv[:, :, 1, None] * north[None, :, :]


def _coast_orientation(lat: np.ndarray, lon: np.ndarray, land: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xyz = _xyz_from_lat_lon(lat, lon)
    land = np.asarray(land, dtype=bool)
    k = max(1, min(6, xyz.shape[0] - 1))
    _, idx = cKDTree(xyz).query(xyz, k=k + 1)
    to_ocean = np.zeros((xyz.shape[0], 3), dtype=np.float64)
    strength = np.zeros(xyz.shape[0], dtype=np.float64)
    for cell, row in enumerate(np.asarray(idx, dtype=np.int64)):
        if not land[cell]:
            continue
        targets = [int(nb) for nb in row[1:] if not land[int(nb)]]
        if not targets:
            continue
        vec = np.mean(xyz[targets], axis=0)
        tangent = vec - float(vec @ xyz[cell]) * xyz[cell]
        norm = float(np.linalg.norm(tangent))
        if norm <= 1.0e-12:
            continue
        to_ocean[cell] = tangent / norm
        strength[cell] = min(1.0, len(targets) / max(k, 1))
    return to_ocean, strength


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


def _coastal_wind_metrics(
    lat: np.ndarray,
    lon: np.ndarray,
    land: np.ndarray,
    seasonal_wind: np.ndarray,
) -> dict[str, float]:
    speed = np.linalg.norm(seasonal_wind, axis=2)
    wind_unit = np.where(
        speed[:, :, None] > 1.0e-9,
        seasonal_wind / np.maximum(speed[:, :, None], 1.0e-9),
        0.0,
    )
    to_ocean, coast_strength = _coast_orientation(lat, lon, land)
    coastal_land = np.asarray(land, dtype=bool) & (coast_strength > 0.0)
    if not coastal_land.any():
        return {
            "coastal_land_cell_count": 0.0,
            "coastal_onshore_abs_p75": float("nan"),
            "coastal_onshore_seasonal_amplitude_p75": float("nan"),
        }
    onshore = (
        np.einsum("snk,nk->sn", wind_unit, -to_ocean)
        * coast_strength[None, :]
    )
    abs_onshore = np.max(np.abs(onshore), axis=0)
    seasonal_amp = np.max(onshore, axis=0) - np.min(onshore, axis=0)
    return {
        "coastal_land_cell_count": float(np.count_nonzero(coastal_land)),
        "coastal_onshore_abs_p75": _percentile(abs_onshore, 75, coastal_land),
        "coastal_onshore_seasonal_amplitude_p75": _percentile(
            seasonal_amp, 75, coastal_land),
    }


def _earth_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        lat = np.asarray(z["lat"], dtype=np.float64)
        lon = np.asarray(z["lon"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        elev = np.asarray(z["earth__elevation_m"], dtype=np.float64)
        land = elev >= 0.0
        seasonal_wind = _uv_to_tangent(
            lat, lon, np.asarray(z["earth__seasonal_wind_u10_v10"], dtype=np.float64))
        current_uv = np.asarray(z["earth__annual_surface_current_u_v"], dtype=np.float64)
        current = _uv_to_tangent(lat, lon, current_uv)
    wind_speed = np.linalg.norm(seasonal_wind, axis=2)
    current_speed = np.linalg.norm(current, axis=1)
    ocean_valid = (~land) & np.isfinite(current_speed)
    all_cells = np.ones_like(land, dtype=bool)
    return {
        "label": "earth_reference",
        "preset": "earth_reference",
        "seed": 0,
        "land_fraction": _area_fraction(area, land, all_cells),
        "wind_speed_p50_m_s": _percentile(wind_speed, 50),
        "wind_speed_p90_m_s": _percentile(wind_speed, 90),
        "current_speed_p50_m_s": _percentile(current_speed, 50, ocean_valid),
        "current_speed_p90_m_s": _percentile(current_speed, 90, ocean_valid),
        **_coastal_wind_metrics(lat, lon, land, seasonal_wind),
    }


def _generated_metrics(summary_row: dict[str, Any], earth: dict[str, Any]) -> dict[str, Any]:
    label = Path(str(summary_row.get("assets_dir", summary_row.get("label", "")))).name
    path = _array_path(summary_row)
    defaults = {
        "label": label,
        "preset": str(summary_row.get("preset", "")),
        "seed": int(summary_row.get("seed", -1)),
        "arrays_found": 0.0,
        "arrays_path": "",
        "land_fraction": 0.0,
        "wind_speed_p50_m_s": float("nan"),
        "wind_speed_p90_m_s": float("nan"),
        "wind_speed_p90_ratio_to_earth": float("nan"),
        "thermal_to_background_wind_p90_ratio": float("nan"),
        "orographic_to_background_wind_p90_ratio": float("nan"),
        "coastal_land_cell_count": 0.0,
        "coastal_onshore_abs_p75": float("nan"),
        "coastal_onshore_seasonal_amplitude_p75": float("nan"),
        "current_speed_p50_m_s": float("nan"),
        "current_speed_p90_m_s": float("nan"),
        "current_speed_p90_ratio_to_earth": float("nan"),
        "current_land_leak_fraction": float("nan"),
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
            "atmosphere__seasonal_wind",
            "atmosphere__background_seasonal_wind",
            "atmosphere__thermal_wind_anomaly",
            "atmosphere__orographic_wind_anomaly",
            "ocean__currents",
        )
        if any(key not in z.files for key in required):
            return defaults
        lat = np.asarray(z["lat"], dtype=np.float64)
        lon = np.asarray(z["lon"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        elev = np.asarray(z["terrain__elevation_m"], dtype=np.float64)
        seasonal_wind = np.asarray(z["atmosphere__seasonal_wind"], dtype=np.float64)
        background = np.asarray(z["atmosphere__background_seasonal_wind"], dtype=np.float64)
        thermal = np.asarray(z["atmosphere__thermal_wind_anomaly"], dtype=np.float64)
        orographic = np.asarray(z["atmosphere__orographic_wind_anomaly"], dtype=np.float64)
        currents = np.asarray(z["ocean__currents"], dtype=np.float64)
    land = elev >= sea
    ocean = ~land
    all_cells = np.ones_like(land, dtype=bool)
    wind_speed = np.linalg.norm(seasonal_wind, axis=2)
    background_speed = np.linalg.norm(background, axis=2)
    thermal_speed = np.linalg.norm(thermal, axis=2)
    orographic_speed = np.linalg.norm(orographic, axis=2)
    current_speed = np.linalg.norm(currents, axis=1)
    bg_p90 = max(_percentile(background_speed, 90), 1.0e-9)
    earth_wind_p90 = _safe_float(earth.get("wind_speed_p90_m_s"))
    earth_current_p90 = _safe_float(earth.get("current_speed_p90_m_s"))
    defaults.update({
        "land_fraction": _area_fraction(area, land, all_cells),
        "wind_speed_p50_m_s": _percentile(wind_speed, 50),
        "wind_speed_p90_m_s": _percentile(wind_speed, 90),
        "wind_speed_p90_ratio_to_earth": (
            _percentile(wind_speed, 90) / earth_wind_p90
            if earth_wind_p90 > 0.0 else float("nan")
        ),
        "thermal_to_background_wind_p90_ratio": _percentile(
            thermal_speed, 90) / bg_p90,
        "orographic_to_background_wind_p90_ratio": _percentile(
            orographic_speed, 95) / bg_p90,
        "current_speed_p50_m_s": _percentile(current_speed, 50, ocean),
        "current_speed_p90_m_s": _percentile(current_speed, 90, ocean),
        "current_speed_p90_ratio_to_earth": (
            _percentile(current_speed, 90, ocean) / earth_current_p90
            if earth_current_p90 > 0.0 else float("nan")
        ),
        "current_land_leak_fraction": float(np.mean(current_speed[land] > 1.0e-6))
        if land.any() else 0.0,
    })
    defaults.update(_coastal_wind_metrics(lat, lon, land, seasonal_wind))
    return defaults


def _checks_for_row(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    preset = str(row.get("preset", "")).lower()
    checks = [
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="circulation gate requires archived arrays"),
        _check(label=label, group="wind_speed", metric="wind_speed_p90_ratio_to_earth",
               generated=row["wind_speed_p90_ratio_to_earth"], operator=">=", threshold=0.45,
               severity="fail", message="wind p90 should remain in Earthlike order of magnitude"),
        _check(label=label, group="wind_speed", metric="wind_speed_p90_ratio_to_earth",
               generated=row["wind_speed_p90_ratio_to_earth"], operator="<=", threshold=1.75,
               severity="fail", message="wind p90 should not remain an over-strong zonal template"),
        _check(label=label, group="coastal_response", metric="coastal_onshore_abs_p75",
               generated=row["coastal_onshore_abs_p75"], operator=">=", threshold=0.22,
               severity="fail", message="coastal winds should retain onshore/offshore structure"),
        _check(label=label, group="coastal_response",
               metric="coastal_onshore_seasonal_amplitude_p75",
               generated=row["coastal_onshore_seasonal_amplitude_p75"],
               operator=">=", threshold=0.32, severity="fail",
               message="coastal winds should have seasonal land-sea response"),
        _check(label=label, group="current_speed", metric="current_speed_p90_ratio_to_earth",
               generated=row["current_speed_p90_ratio_to_earth"], operator=">=", threshold=0.45,
               severity="fail", message="surface currents should not collapse to near-zero"),
        _check(label=label, group="current_speed", metric="current_speed_p90_ratio_to_earth",
               generated=row["current_speed_p90_ratio_to_earth"], operator="<=", threshold=1.80,
               severity="fail", message="surface currents should stay near OSCAR magnitude"),
        _check(label=label, group="land_mask", metric="current_land_leak_fraction",
               generated=row["current_land_leak_fraction"], operator="<=", threshold=0.001,
               severity="fail", message="surface currents must not leak onto land"),
    ]
    if "earthlike" in preset or "earthlike" in label.lower():
        checks.append(
            _check(label=label, group="geographic_wind_anomaly",
                   metric="thermal_to_background_wind_p90_ratio",
                   generated=row["thermal_to_background_wind_p90_ratio"],
                   operator=">=", threshold=0.06, severity="fail",
                   message="earthlike winds should contain measurable land-sea thermal anomalies")
        )
    if "waterworld" in preset or "waterworld" in label.lower():
        checks.append(
            _check(label=label, group="waterworld_guard",
                   metric="thermal_to_background_wind_p90_ratio",
                   generated=row["thermal_to_background_wind_p90_ratio"],
                   operator="<=", threshold=0.08, severity="fail",
                   message="island worlds should not grow continent-scale wind anomalies")
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
        "# Earth Climate Circulation Layout Gate",
        "",
        f"Schema: `{report['schema']}`",
        f"Verdict: **{report['verdict']}**",
        "",
        (
            f"Failures: `{report['failure_count']}`; warnings: "
            f"`{report['warning_count']}`; skipped: `{report['skipped_count']}`."
        ),
        "",
        "This gate evaluates the W/O layers in the documented climate pipeline: "
        "seasonal wind and basin-confined surface currents.  It should be run "
        "before moisture, precipitation, Koppen, or biome tuning.",
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
    lines.append("- `earth_climate_circulation_layout_metrics.csv`")
    lines.append("- `earth_climate_circulation_layout_checks.csv`")
    return "\n".join(lines) + "\n"


def run_earth_climate_circulation_layout_gate(
    config: EarthClimateCirculationLayoutGateConfig,
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
    checks = [check for row in generated for check in _checks_for_row(row)]
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
        "label", "preset", "seed", "arrays_found", "arrays_path", "land_fraction",
        "wind_speed_p50_m_s", "wind_speed_p90_m_s",
        "wind_speed_p90_ratio_to_earth",
        "thermal_to_background_wind_p90_ratio",
        "orographic_to_background_wind_p90_ratio",
        "coastal_land_cell_count", "coastal_onshore_abs_p75",
        "coastal_onshore_seasonal_amplitude_p75",
        "current_speed_p50_m_s", "current_speed_p90_m_s",
        "current_speed_p90_ratio_to_earth", "current_land_leak_fraction",
    ]
    check_keys = [
        "label", "group", "metric", "generated", "operator", "threshold",
        "severity", "passed", "skipped", "message",
    ]
    _write_csv(outdir / "earth_climate_circulation_layout_metrics.csv",
               generated, metric_keys)
    _write_csv(outdir / "earth_climate_circulation_layout_checks.csv",
               checks, check_keys)
    (outdir / "earth_climate_circulation_layout_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default))
    (outdir / "earth_climate_circulation_layout_gate_report.md").write_text(
        _render_markdown(report))
    return report
