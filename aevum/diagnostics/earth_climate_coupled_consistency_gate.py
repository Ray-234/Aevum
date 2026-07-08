"""F3 coupled pressure/wind/moisture consistency gate.

This gate sits after the C5b1 ocean-spatial pass.  It does not try to make
generated worlds copy Earth's geography.  Instead it checks that the reduced
climate fields exchange the right first-order information: seasonal warmth
drives pressure, pressure steers wind, SST/current state drives evaporation and
source humidity, moisture/support fields explain precipitation, and high
monsoon potential is not created in dry interiors.
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


SCHEMA = "aevum.earth_climate_coupled_consistency_gate.v1"


@dataclass(frozen=True)
class EarthClimateCoupledConsistencyGateConfig:
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


def _knn(lat: np.ndarray, lon: np.ndarray, k: int = 6) -> tuple[np.ndarray, np.ndarray]:
    xyz = _xyz_from_lat_lon(lat, lon)
    k = max(1, min(int(k), xyz.shape[0] - 1))
    _, idx = cKDTree(xyz).query(xyz, k=k + 1)
    return np.asarray(idx[:, 1:], dtype=np.int64), xyz


def _graph_gradient(lat: np.ndarray, lon: np.ndarray, values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    idx, xyz = _knn(lat, lon, k=6)
    grad = np.zeros((values.size, 3), dtype=np.float64)
    weight = np.zeros(values.size, dtype=np.float64)
    for i, row in enumerate(idx):
        xi = xyz[i]
        for j in row:
            direction = xyz[int(j)] - float(xyz[int(j)] @ xi) * xi
            norm = float(np.linalg.norm(direction))
            if norm <= 1.0e-12:
                continue
            grad[i] += (values[int(j)] - values[i]) * direction / norm
            weight[i] += 1.0
    return grad / np.maximum(weight[:, None], 1.0)


def _percentile(values: np.ndarray, q: float,
                mask: np.ndarray | None = None) -> float:
    values = np.asarray(values, dtype=np.float64)
    if mask is not None:
        values = values[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    return float(np.percentile(values, q)) if values.size else float("nan")


def _weighted_mean(values: np.ndarray, area: np.ndarray, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values) & np.isfinite(area)
    if not mask.any():
        return float("nan")
    return float(np.average(values[mask], weights=area[mask]))


def _corr(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    mask = np.asarray(mask, dtype=bool).ravel()
    mask &= np.isfinite(a) & np.isfinite(b)
    if int(np.count_nonzero(mask)) < 8:
        return float("nan")
    aa = a[mask] - float(np.mean(a[mask]))
    bb = b[mask] - float(np.mean(b[mask]))
    denom = float(np.sqrt(np.sum(aa * aa) * np.sum(bb * bb)))
    if denom <= 1.0e-12:
        return float("nan")
    return float(np.sum(aa * bb) / denom)


def _area_fraction(area: np.ndarray, mask: np.ndarray, domain: np.ndarray) -> float:
    area = np.asarray(area, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    domain = np.asarray(domain, dtype=bool)
    denom = float(np.sum(area[domain & np.isfinite(area)]))
    if denom <= 1.0e-12:
        return float("nan")
    return float(np.sum(area[mask & domain & np.isfinite(area)]) / denom)


def _wind_pressure_alignment(
    lat: np.ndarray,
    lon: np.ndarray,
    pressure: np.ndarray,
    wind: np.ndarray,
) -> float:
    rows: list[float] = []
    for season in range(4):
        high_to_low = -_graph_gradient(lat, lon, pressure[season])
        grad_speed = np.linalg.norm(high_to_low, axis=1)
        wind_speed = np.linalg.norm(wind[season], axis=1)
        alignment = np.einsum("ij,ij->i", high_to_low, wind[season]) / np.maximum(
            grad_speed * wind_speed, 1.0e-12)
        active = (
            (grad_speed >= _percentile(grad_speed, 60))
            & (wind_speed >= _percentile(wind_speed, 40))
        )
        rows.append(_percentile(alignment, 60, active))
    rows = [x for x in rows if np.isfinite(x)]
    return float(np.median(rows)) if rows else float("nan")


def _earth_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        lat = np.asarray(z["lat"], dtype=np.float64)
        lon = np.asarray(z["lon"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        land = np.asarray(z["earth__land_mask"], dtype=bool)
        pressure = np.asarray(z["earth__seasonal_slp_anomaly_hPa"],
                              dtype=np.float64)
        wind = _uv_to_tangent(
            lat, lon, np.asarray(z["earth__seasonal_wind_u10_v10"],
                                 dtype=np.float64))
        temp = np.asarray(z["earth__seasonal_temperature_C"],
                          dtype=np.float64) + 273.15
    land4 = np.broadcast_to(land, pressure.shape)
    temp_anom = temp - np.mean(temp, axis=0, keepdims=True)
    all_cells = np.ones_like(land, dtype=bool)
    return {
        "label": "earth_reference",
        "preset": "earth_reference",
        "seed": 0,
        "land_fraction": _area_fraction(area, land, all_cells),
        "pressure_temperature_land_corr": _corr(pressure, temp_anom, land4),
        "wind_pressure_alignment_p60": _wind_pressure_alignment(
            lat, lon, pressure, wind),
    }


def _spatial_metrics(
    *,
    label: str,
    preset: str,
    seed: int,
    arrays_path: str,
    lat: np.ndarray,
    lon: np.ndarray,
    area: np.ndarray,
    land: np.ndarray,
    seasonal_temp: np.ndarray,
    seasonal_sst: np.ndarray,
    pressure: np.ndarray,
    wind: np.ndarray,
    source_warmth: np.ndarray,
    moisture: np.ndarray,
    monsoon: np.ndarray,
    seasonal_precip: np.ndarray,
    annual_precip: np.ndarray,
    evaporation: np.ndarray,
    heat_transport: np.ndarray,
    upwelling: np.ndarray,
    ocean_heat_flux: np.ndarray,
    coupling_residual: np.ndarray,
    storm_track: np.ndarray,
    itcz: np.ndarray,
) -> dict[str, Any]:
    land = np.asarray(land, dtype=bool)
    ocean = ~land
    land4 = np.broadcast_to(land, pressure.shape)
    ocean4 = np.broadcast_to(ocean, pressure.shape)
    temp_anom = seasonal_temp - np.mean(seasonal_temp, axis=0, keepdims=True)
    support = (
        np.clip(moisture, 0.0, 1.0)
        + 0.55 * np.clip(monsoon, 0.0, 1.2)
        + 0.35 * np.clip(storm_track, 0.0, 1.2)
        + 0.25 * np.clip(itcz, 0.0, 1.2)
    )

    monsoon_positive = np.clip(monsoon, 0.0, None)
    monsoon_threshold = _percentile(monsoon_positive, 85, land4)
    monsoon_top = land4 & np.isfinite(monsoon_positive) & (
        monsoon_positive >= monsoon_threshold)
    land_moisture_median = _percentile(moisture, 50, land4)
    monsoon_top_moisture_p25 = _percentile(moisture, 25, monsoon_top)
    land_precip_median = _percentile(seasonal_precip, 50, land4)
    monsoon_top_precip_p50 = _percentile(seasonal_precip, 50, monsoon_top)

    cold = ocean & (heat_transport < -0.35) & (upwelling > 0.02)
    warm = ocean & (heat_transport > 0.35)
    cold_evap = _weighted_mean(evaporation, area, cold)
    warm_evap = _weighted_mean(evaporation, area, warm)
    heat_mean_abs = abs(_weighted_mean(ocean_heat_flux, area, ocean))

    all_cells = np.ones_like(land, dtype=bool)
    return {
        "label": label,
        "preset": preset,
        "seed": seed,
        "arrays_found": 1.0,
        "arrays_path": arrays_path,
        "required_fields_found": 1.0,
        "land_fraction": _area_fraction(area, land, all_cells),
        "pressure_temperature_land_corr": _corr(pressure, temp_anom, land4),
        "wind_pressure_alignment_p60": _wind_pressure_alignment(
            lat, lon, pressure, wind),
        "source_sst_ocean_corr": _corr(source_warmth, seasonal_sst, ocean4),
        "evap_source_ocean_corr": _corr(
            np.broadcast_to(evaporation, source_warmth.shape), source_warmth, ocean4),
        "moisture_precip_land_corr": _corr(moisture, seasonal_precip, land4),
        "support_precip_land_corr": _corr(support, seasonal_precip, land4),
        "monsoon_pressure_corr_land": _corr(
            monsoon_positive, -pressure, land4),
        "monsoon_moisture_corr_land": _corr(
            monsoon_positive, moisture, land4),
        "monsoon_top_moisture_p25": monsoon_top_moisture_p25,
        "land_moisture_median": land_moisture_median,
        "monsoon_top_moisture_ratio": (
            monsoon_top_moisture_p25 / max(land_moisture_median, 1.0e-9)
        ),
        "monsoon_top_precip_p50": monsoon_top_precip_p50,
        "land_precip_median": land_precip_median,
        "monsoon_top_precip_ratio": (
            monsoon_top_precip_p50 / max(land_precip_median, 1.0e-9)
        ),
        "cold_warm_evap_ratio": (
            cold_evap / max(warm_evap, 1.0e-9)
            if np.isfinite(cold_evap) and np.isfinite(warm_evap) else float("nan")
        ),
        "cold_current_cell_count": float(np.count_nonzero(cold)),
        "warm_current_cell_count": float(np.count_nonzero(warm)),
        "seasonal_precip_aggregate_max_delta": float(
            np.nanmax(np.abs(np.mean(seasonal_precip, axis=0) - annual_precip))),
        "pressure_abs_p95": _percentile(np.abs(pressure), 95),
        "ocean_heat_flux_ocean_mean_abs_C": heat_mean_abs,
        "coupling_residual_ocean_p95_C": _percentile(coupling_residual, 95, ocean),
        "coupled_support_land_p90": _percentile(support, 90, land4),
        "moisture_access_land_p75": _percentile(moisture, 75, land4),
    }


def _generated_metrics(summary_row: dict[str, Any]) -> dict[str, Any]:
    label = Path(str(summary_row.get("assets_dir", summary_row.get("label", "")))).name
    path = _array_path(summary_row)
    defaults: dict[str, Any] = {
        "label": label,
        "preset": str(summary_row.get("preset", "")),
        "seed": int(summary_row.get("seed", -1)),
        "arrays_found": 0.0,
        "arrays_path": "",
        "required_fields_found": 0.0,
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
            "climate__seasonal_temperature",
            "climate__seasonal_sst",
            "climate__seasonal_precipitation",
            "climate__precipitation",
            "climate__evaporation",
            "climate__ocean_heat_flux",
            "climate__coupling_residual",
            "atmosphere__seasonal_pressure_proxy",
            "atmosphere__seasonal_wind",
            "atmosphere__source_ocean_warmth",
            "atmosphere__moisture_access",
            "atmosphere__monsoon_potential",
            "atmosphere__storm_track_intensity",
            "atmosphere__itcz_intensity",
            "ocean__current_heat_transport",
            "ocean__upwelling",
        )
        if any(key not in z.files for key in required):
            return defaults
        lat = np.asarray(z["lat"], dtype=np.float64)
        lon = np.asarray(z["lon"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        land = np.asarray(z["terrain__elevation_m"], dtype=np.float64) >= sea
        row = _spatial_metrics(
            label=label,
            preset=str(summary_row.get("preset", "")),
            seed=int(summary_row.get("seed", -1)),
            arrays_path=str(path),
            lat=lat,
            lon=lon,
            area=area,
            land=land,
            seasonal_temp=np.asarray(
                z["climate__seasonal_temperature"], dtype=np.float64),
            seasonal_sst=np.asarray(z["climate__seasonal_sst"], dtype=np.float64),
            pressure=np.asarray(
                z["atmosphere__seasonal_pressure_proxy"], dtype=np.float64),
            wind=np.asarray(z["atmosphere__seasonal_wind"], dtype=np.float64),
            source_warmth=np.asarray(
                z["atmosphere__source_ocean_warmth"], dtype=np.float64),
            moisture=np.asarray(z["atmosphere__moisture_access"], dtype=np.float64),
            monsoon=np.asarray(
                z["atmosphere__monsoon_potential"], dtype=np.float64),
            seasonal_precip=np.asarray(
                z["climate__seasonal_precipitation"], dtype=np.float64),
            annual_precip=np.asarray(
                z["climate__precipitation"], dtype=np.float64),
            evaporation=np.asarray(z["climate__evaporation"], dtype=np.float64),
            heat_transport=np.asarray(
                z["ocean__current_heat_transport"], dtype=np.float64),
            upwelling=np.asarray(z["ocean__upwelling"], dtype=np.float64),
            ocean_heat_flux=np.asarray(
                z["climate__ocean_heat_flux"], dtype=np.float64),
            coupling_residual=np.asarray(
                z["climate__coupling_residual"], dtype=np.float64),
            storm_track=np.asarray(
                z["atmosphere__storm_track_intensity"], dtype=np.float64),
            itcz=np.asarray(z["atmosphere__itcz_intensity"], dtype=np.float64),
        )
    return row


def _checks_for_row(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    preset = str(row.get("preset", "")).lower()
    checks = [
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row.get("arrays_found", 0.0), operator=">=",
               threshold=1.0, severity="fail",
               message="F3 coupled consistency gate requires archived arrays"),
        _check(label=label, group="array_archive", metric="required_fields_found",
               generated=row.get("required_fields_found", 0.0), operator=">=",
               threshold=1.0, severity="fail",
               message="F3 fields must include pressure, wind, SST, moisture, evaporation, and precipitation"),
        _check(label=label, group="pressure_temperature",
               metric="pressure_temperature_land_corr",
               generated=row.get("pressure_temperature_land_corr", float("nan")),
               operator="<=", threshold=-0.35, severity="fail",
               message="seasonal warm land should correspond to lower pressure proxy"),
        _check(label=label, group="pressure_wind",
               metric="wind_pressure_alignment_p60",
               generated=row.get("wind_pressure_alignment_p60", float("nan")),
               operator=">=", threshold=0.18, severity="fail",
               message="seasonal winds should broadly follow high-to-low pressure gradients"),
        _check(label=label, group="sst_source_evaporation",
               metric="source_sst_ocean_corr",
               generated=row.get("source_sst_ocean_corr", float("nan")),
               operator=">=", threshold=0.70, severity="fail",
               message="source-ocean warmth should be derived from SST structure"),
        _check(label=label, group="sst_source_evaporation",
               metric="evap_source_ocean_corr",
               generated=row.get("evap_source_ocean_corr", float("nan")),
               operator=">=", threshold=0.45, severity="fail",
               message="ocean evaporation should track diagnosed source warmth"),
        _check(label=label, group="moisture_precipitation",
               metric="moisture_precip_land_corr",
               generated=row.get("moisture_precip_land_corr", float("nan")),
               operator=">=", threshold=0.25, severity="fail",
               message="seasonal land precipitation should increase with moisture access"),
        _check(label=label, group="moisture_precipitation",
               metric="support_precip_land_corr",
               generated=row.get("support_precip_land_corr", float("nan")),
               operator=">=", threshold=0.35, severity="fail",
               message="combined moisture/monsoon/storm support should explain seasonal rainfall"),
        _check(label=label, group="monsoon_support",
               metric="monsoon_pressure_corr_land",
               generated=row.get("monsoon_pressure_corr_land", float("nan")),
               operator=">=", threshold=0.18, severity="fail",
               message="positive monsoon potential should align with thermal lows"),
        _check(label=label, group="monsoon_support",
               metric="monsoon_top_moisture_ratio",
               generated=row.get("monsoon_top_moisture_ratio", float("nan")),
               operator=">=", threshold=0.75, severity="fail",
               message="top monsoon-potential land should not be drier than the land median"),
        _check(label=label, group="monsoon_support",
               metric="monsoon_top_precip_ratio",
               generated=row.get("monsoon_top_precip_ratio", float("nan")),
               operator=">=", threshold=1.20, severity="fail",
               message="top monsoon-potential land should receive enhanced seasonal precipitation"),
        _check(label=label, group="cold_current_evaporation",
               metric="cold_warm_evap_ratio",
               generated=row.get("cold_warm_evap_ratio", float("nan")),
               operator="<=", threshold=0.90, severity="fail",
               message="cold-current/upwelling source regions should evaporate less than warm-current regions"),
        _check(label=label, group="budget_closure",
               metric="seasonal_precip_aggregate_max_delta",
               generated=row.get("seasonal_precip_aggregate_max_delta", float("nan")),
               operator="<=", threshold=1.0e-6, severity="fail",
               message="seasonal precipitation must aggregate exactly to annual precipitation"),
        _check(label=label, group="budget_closure",
               metric="ocean_heat_flux_ocean_mean_abs_C",
               generated=row.get("ocean_heat_flux_ocean_mean_abs_C", float("nan")),
               operator="<=", threshold=0.20, severity="fail",
               message="ocean heat flux should redistribute rather than add net energy"),
        _check(label=label, group="coupling_stability",
               metric="coupling_residual_ocean_p95_C",
               generated=row.get("coupling_residual_ocean_p95_C", float("nan")),
               operator="<=", threshold=1.50, severity="fail",
               message="weak ocean-atmosphere coupling should remain bounded"),
    ]
    if "waterworld" in preset or "waterworld" in label.lower():
        checks.append(
            _check(label=label, group="waterworld_pressure",
                   metric="pressure_abs_p95",
                   generated=row.get("pressure_abs_p95", float("nan")),
                   operator="<=", threshold=0.55, severity="fail",
                   message="waterworlds should not invent continent-scale pressure contrast")
        )
    return checks


FIELDNAMES = [
    "label",
    "preset",
    "seed",
    "arrays_found",
    "arrays_path",
    "required_fields_found",
    "land_fraction",
    "pressure_temperature_land_corr",
    "wind_pressure_alignment_p60",
    "source_sst_ocean_corr",
    "evap_source_ocean_corr",
    "moisture_precip_land_corr",
    "support_precip_land_corr",
    "monsoon_pressure_corr_land",
    "monsoon_moisture_corr_land",
    "monsoon_top_moisture_p25",
    "land_moisture_median",
    "monsoon_top_moisture_ratio",
    "monsoon_top_precip_p50",
    "land_precip_median",
    "monsoon_top_precip_ratio",
    "cold_warm_evap_ratio",
    "cold_current_cell_count",
    "warm_current_cell_count",
    "seasonal_precip_aggregate_max_delta",
    "pressure_abs_p95",
    "ocean_heat_flux_ocean_mean_abs_C",
    "coupling_residual_ocean_p95_C",
    "coupled_support_land_p90",
    "moisture_access_land_p75",
]


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _render_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Earth Climate Coupled-Consistency Gate",
        "",
        f"Schema: `{report['schema']}`",
        f"Verdict: **{report['verdict']}**",
        f"Failures: {report['failure_count']}",
        f"Warnings: {report['warning_count']}",
        f"Skipped checks: {report['skipped_count']}",
        "",
        "## Earth Reference",
        "",
    ]
    earth = report["earth_metrics"]
    for key in (
        "pressure_temperature_land_corr",
        "wind_pressure_alignment_p60",
    ):
        lines.append(f"- `{key}`: `{earth.get(key)}`")
    lines.extend(["", "## Generated Runs", ""])
    for row in report["generated"]:
        lines.append(
            f"- `{row['label']}`: pressure/temp corr "
            f"`{row.get('pressure_temperature_land_corr')}`, "
            f"wind/pressure align `{row.get('wind_pressure_alignment_p60')}`, "
            f"monsoon-top moisture ratio "
            f"`{row.get('monsoon_top_moisture_ratio')}`, "
            f"support/precip corr `{row.get('support_precip_land_corr')}`"
        )
    if report["failures"]:
        lines.extend(["", "## Failures", ""])
        for row in report["failures"]:
            lines.append(
                f"- `{row['label']}` `{row['group']}.{row['metric']}` "
                f"{row['operator']} `{row['threshold']}` got `{row['generated']}`: "
                f"{row['message']}"
            )
    return "\n".join(lines) + "\n"


def run_earth_climate_coupled_consistency_gate(
    config: EarthClimateCoupledConsistencyGateConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    earth = _earth_metrics(Path(config.earth_reference_npz))
    with Path(config.terminal_summary_json).open() as f:
        summary = json.load(f)
    generated = [_generated_metrics(row) for row in summary.get("summaries", [])]
    checks = [check for row in generated for check in _checks_for_row(row)]
    failures = [
        row for row in checks
        if not row["passed"] and not row["skipped"] and row["severity"] == "fail"
    ]
    warnings = [
        row for row in checks
        if not row["passed"] and not row["skipped"] and row["severity"] == "warn"
    ]
    skipped = [row for row in checks if row["skipped"]]
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
    _write_csv(outdir / "earth_climate_coupled_consistency_metrics.csv",
               generated, FIELDNAMES)
    _write_csv(
        outdir / "earth_climate_coupled_consistency_checks.csv",
        checks,
        [
            "label", "group", "metric", "generated", "operator", "threshold",
            "severity", "passed", "skipped", "message",
        ],
    )
    (outdir / "earth_climate_coupled_consistency_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default))
    (outdir / "earth_climate_coupled_consistency_gate_report.md").write_text(
        _render_markdown(report))
    return report
