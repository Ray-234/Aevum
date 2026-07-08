"""Windward/leeward precipitation asymmetry gate against Earth envelopes."""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree


SCHEMA = "aevum.earth_climate_windward_leeward_gate.v1"


@dataclass(frozen=True)
class EarthClimateWindwardLeewardGateConfig:
    earth_reference_npz: Path
    terminal_summary_json: Path
    outdir: Path


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
    east_norm = np.linalg.norm(east, axis=1, keepdims=True)
    east = np.where(east_norm > 1.0e-9, east / np.maximum(east_norm, 1.0e-9),
                    np.array([1.0, 0.0, 0.0]))
    north = np.cross(xyz, east)
    return xyz, east, north


def _knn_edges(lat: np.ndarray, lon: np.ndarray, k: int = 6) -> tuple[np.ndarray, np.ndarray]:
    xyz = _xyz_from_lat_lon(lat, lon)
    k = max(1, min(int(k), xyz.shape[0] - 1))
    _, idx = cKDTree(xyz).query(xyz, k=k + 1)
    edges: set[tuple[int, int]] = set()
    for i, row in enumerate(np.asarray(idx, dtype=np.int64)):
        for j in row[1:]:
            a, b = sorted((int(i), int(j)))
            edges.add((a, b))
    return xyz, np.asarray(sorted(edges), dtype=np.int64)


def _graph_gradient_vectors(
    lat: np.ndarray,
    lon: np.ndarray,
    values: np.ndarray,
    edges: np.ndarray,
) -> np.ndarray:
    xyz = _xyz_from_lat_lon(lat, lon)
    values = np.asarray(values, dtype=np.float64)
    i = edges[:, 0]
    j = edges[:, 1]
    delta_pos = xyz[j] - xyz[i]
    dir_i = delta_pos - np.sum(delta_pos * xyz[i], axis=1, keepdims=True) * xyz[i]
    dir_i_norm = np.linalg.norm(dir_i, axis=1, keepdims=True)
    dir_i = np.where(dir_i_norm > 1.0e-12,
                     dir_i / np.maximum(dir_i_norm, 1.0e-12), 0.0)

    delta_neg = -delta_pos
    dir_j = delta_neg - np.sum(delta_neg * xyz[j], axis=1, keepdims=True) * xyz[j]
    dir_j_norm = np.linalg.norm(dir_j, axis=1, keepdims=True)
    dir_j = np.where(dir_j_norm > 1.0e-12,
                     dir_j / np.maximum(dir_j_norm, 1.0e-12), 0.0)

    grad = np.zeros((values.size, 3), dtype=np.float64)
    degree = np.zeros(values.size, dtype=np.float64)
    delta = values[j] - values[i]
    np.add.at(grad, i, delta[:, None] * dir_i)
    np.add.at(grad, j, -delta[:, None] * dir_j)
    np.add.at(degree, i, 1.0)
    np.add.at(degree, j, 1.0)
    return grad / np.maximum(degree[:, None], 1.0)


def _safe_weighted_mean(values: np.ndarray, area: np.ndarray, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values) & np.isfinite(area)
    if not mask.any():
        return float("nan")
    return float(np.average(values[mask], weights=area[mask]))


def _area_fraction(area: np.ndarray, mask: np.ndarray, denom_mask: np.ndarray) -> float:
    area = np.asarray(area, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    denom_mask = np.asarray(denom_mask, dtype=bool)
    denom = float(np.nansum(area[denom_mask]))
    if denom <= 0.0:
        return float("nan")
    return float(np.nansum(area[mask & denom_mask]) / denom)


def _earth_wind_3d(lat: np.ndarray, lon: np.ndarray, seasonal_uv: np.ndarray) -> np.ndarray:
    _, east, north = _tangent_basis(lat, lon)
    seasonal_uv = np.asarray(seasonal_uv, dtype=np.float64)
    return (
        seasonal_uv[:, :, 0, None] * east[None, :, :]
        + seasonal_uv[:, :, 1, None] * north[None, :, :]
    )


def _windward_metrics(
    label: str,
    preset: str,
    seed: int,
    lat: np.ndarray,
    lon: np.ndarray,
    area: np.ndarray,
    land: np.ndarray,
    relative_elev: np.ndarray,
    precip: np.ndarray,
    seasonal_precip: np.ndarray,
    seasonal_wind: np.ndarray,
    *,
    barrier: np.ndarray | None = None,
    orographic_precip: np.ndarray | None = None,
) -> dict[str, Any]:
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)
    land = np.asarray(land, dtype=bool)
    relative_elev = np.asarray(relative_elev, dtype=np.float64)
    precip = np.asarray(precip, dtype=np.float64)
    seasonal_precip = np.asarray(seasonal_precip, dtype=np.float64)
    seasonal_wind = np.asarray(seasonal_wind, dtype=np.float64)

    _, edges = _knn_edges(lat, lon, k=6)
    topo = np.maximum(relative_elev, 0.0)
    topo_scale = max(float(np.nanpercentile(topo[land], 95)) if land.any() else 1.0, 1.0)
    topo_norm = np.clip(topo / topo_scale, 0.0, 2.0)
    topo_grad = _graph_gradient_vectors(lat, lon, topo_norm, edges)
    topo_mag = np.linalg.norm(topo_grad, axis=1)
    if land.any():
        g70 = float(np.nanpercentile(topo_mag[land], 70))
        g92 = max(float(np.nanpercentile(topo_mag[land], 92)), g70 + 1.0e-9)
    else:
        g70, g92 = 0.0, 1.0
    mountain_slope = land & (relative_elev >= 600.0) & (topo_mag >= g70)
    strong_slope = mountain_slope & (topo_mag >= g92)

    wind_speed = np.linalg.norm(seasonal_wind, axis=2)
    wind_unit = np.where(
        wind_speed[:, :, None] > 1.0e-9,
        seasonal_wind / np.maximum(wind_speed[:, :, None], 1.0e-9),
        0.0,
    )
    grad_unit = np.where(
        topo_mag[:, None] > 1.0e-9,
        topo_grad / np.maximum(topo_mag[:, None], 1.0e-9),
        0.0,
    )
    wind_uphill = np.einsum("snk,nk->sn", wind_unit, grad_unit)

    seasonal_ratios: list[float] = []
    seasonal_windward_means: list[float] = []
    seasonal_leeward_means: list[float] = []
    min_side_cells = max(2, int(0.015 * max(int(mountain_slope.sum()), 1)))
    for s in range(min(4, seasonal_wind.shape[0])):
        windward = mountain_slope & (wind_uphill[s] > 0.25)
        leeward = mountain_slope & (wind_uphill[s] < -0.25)
        if int(windward.sum()) < min_side_cells or int(leeward.sum()) < min_side_cells:
            continue
        windward_mean = _safe_weighted_mean(seasonal_precip[s], area, windward)
        leeward_mean = _safe_weighted_mean(seasonal_precip[s], area, leeward)
        if np.isfinite(windward_mean) and np.isfinite(leeward_mean) and leeward_mean > 0.0:
            seasonal_windward_means.append(windward_mean)
            seasonal_leeward_means.append(leeward_mean)
            seasonal_ratios.append(float(windward_mean / leeward_mean))

    annual_weight = np.maximum(wind_speed, 1.0e-6)
    annual_dot = np.average(wind_uphill, axis=0, weights=annual_weight)
    annual_windward = mountain_slope & (annual_dot > 0.20)
    annual_leeward = mountain_slope & (annual_dot < -0.20)
    annual_windward_mean = _safe_weighted_mean(precip, area, annual_windward)
    annual_leeward_mean = _safe_weighted_mean(precip, area, annual_leeward)
    annual_ratio = (
        float(annual_windward_mean / annual_leeward_mean)
        if np.isfinite(annual_windward_mean)
        and np.isfinite(annual_leeward_mean)
        and annual_leeward_mean > 0.0
        else float("nan")
    )

    row = {
        "label": label,
        "preset": preset,
        "seed": seed,
        "land_fraction": _area_fraction(area, land, np.ones_like(land, dtype=bool)),
        "mountain_slope_land_fraction": _area_fraction(area, mountain_slope, land),
        "strong_slope_land_fraction": _area_fraction(area, strong_slope, land),
        "annual_windward_area_fraction_of_land": _area_fraction(area, annual_windward, land),
        "annual_leeward_area_fraction_of_land": _area_fraction(area, annual_leeward, land),
        "annual_windward_precip_mean_mm_yr": annual_windward_mean,
        "annual_leeward_precip_mean_mm_yr": annual_leeward_mean,
        "annual_windward_leeward_precip_ratio": annual_ratio,
        "seasonal_valid_count": len(seasonal_ratios),
        "seasonal_windward_leeward_precip_ratio_min": (
            float(np.min(seasonal_ratios)) if seasonal_ratios else float("nan")
        ),
        "seasonal_windward_leeward_precip_ratio_median": (
            float(np.median(seasonal_ratios)) if seasonal_ratios else float("nan")
        ),
        "seasonal_windward_leeward_precip_ratio_max": (
            float(np.max(seasonal_ratios)) if seasonal_ratios else float("nan")
        ),
        "seasonal_windward_precip_mean_median_mm_yr": (
            float(np.median(seasonal_windward_means))
            if seasonal_windward_means else float("nan")
        ),
        "seasonal_leeward_precip_mean_median_mm_yr": (
            float(np.median(seasonal_leeward_means))
            if seasonal_leeward_means else float("nan")
        ),
    }
    if barrier is not None:
        row["barrier_p75_mountain_slope"] = float(
            np.nanpercentile(np.asarray(barrier, dtype=np.float64)[mountain_slope], 75)
        ) if mountain_slope.any() else float("nan")
    else:
        row["barrier_p75_mountain_slope"] = float("nan")
    if orographic_precip is not None:
        row["orographic_precip_p75_mountain_slope_mm_yr"] = float(
            np.nanpercentile(
                np.asarray(orographic_precip, dtype=np.float64)[mountain_slope], 75)
        ) if mountain_slope.any() else float("nan")
    else:
        row["orographic_precip_p75_mountain_slope_mm_yr"] = float("nan")
    return row


def _earth_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        lat = np.asarray(z["lat"], dtype=np.float64)
        lon = np.asarray(z["lon"], dtype=np.float64)
        return _windward_metrics(
            "Earth windward/leeward reference",
            "earth_reference",
            -1,
            lat,
            lon,
            np.asarray(z["cell_area"], dtype=np.float64),
            np.asarray(z["earth__land_mask"], dtype=bool),
            np.asarray(z["earth__elevation_m"], dtype=np.float64),
            np.asarray(z["earth__annual_precip_mm"], dtype=np.float64),
            np.asarray(z["earth__seasonal_precip_mm_yr_equiv"], dtype=np.float64),
            _earth_wind_3d(lat, lon, np.asarray(z["earth__seasonal_wind_u10_v10"])),
        )


def _generated_metrics(summary: dict[str, Any]) -> dict[str, Any]:
    arrays = Path(summary["arrays"])
    with np.load(arrays, allow_pickle=False) as z:
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        elev = np.asarray(z["terrain__elevation_m"], dtype=np.float64)
        land = elev >= sea
        seasonal_wind = np.asarray(z["atmosphere__seasonal_wind"], dtype=np.float64)
        row = _windward_metrics(
            Path(str(summary.get("assets_dir", ""))).name,
            str(summary.get("preset", "")),
            int(summary.get("seed", -1)),
            np.asarray(z["lat"], dtype=np.float64),
            np.asarray(z["lon"], dtype=np.float64),
            np.asarray(z["cell_area"], dtype=np.float64),
            land,
            elev - sea,
            np.asarray(z["climate__precipitation"], dtype=np.float64),
            np.asarray(z["climate__seasonal_precipitation"], dtype=np.float64),
            seasonal_wind,
            barrier=(
                np.asarray(z["terrain__barrier_index"], dtype=np.float64)
                if "terrain__barrier_index" in z else None
            ),
            orographic_precip=(
                np.asarray(z["climate__orographic_precipitation"], dtype=np.float64)
                if "climate__orographic_precipitation" in z else None
            ),
        )
    row["arrays"] = str(arrays)
    row["mode"] = (
        "earthlike_calibration"
        if "earthlike" in row["preset"].lower()
        else "diagnostic_only"
    )
    return row


def _ratio(value: float, ref: float) -> float:
    if not np.isfinite(value) or not np.isfinite(ref) or abs(ref) <= 1.0e-12:
        return float("nan")
    return float(value / ref)


def _check(
    label: str,
    metric: str,
    generated: float,
    reference: float,
    operator: str,
    threshold: float,
    severity: str,
    message: str,
) -> dict[str, Any]:
    ratio = _ratio(generated, reference)
    skipped = not np.isfinite(float(generated))
    if skipped:
        passed = True
    elif operator == "ratio>=":
        passed = np.isfinite(ratio) and ratio >= threshold
    elif operator == ">=":
        passed = float(generated) >= threshold
    elif operator == "<=":
        passed = float(generated) <= threshold
    else:
        raise ValueError(f"unsupported operator {operator!r}")
    return {
        "label": label,
        "metric": metric,
        "generated": float(generated),
        "reference": float(reference),
        "ratio_to_reference": ratio,
        "operator": operator,
        "threshold": float(threshold),
        "severity": severity,
        "passed": bool(passed),
        "skipped": bool(skipped),
        "message": message,
    }


def _earthlike_checks(row: dict[str, Any], earth: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(
            label,
            "mountain_slope_land_fraction",
            row["mountain_slope_land_fraction"],
            earth["mountain_slope_land_fraction"],
            "ratio>=",
            0.45,
            "warn",
            "diagnostic mountain-slope sample should not collapse relative to Earth",
        ),
        _check(
            label,
            "annual_windward_leeward_precip_ratio",
            row["annual_windward_leeward_precip_ratio"],
            earth["annual_windward_leeward_precip_ratio"],
            "ratio>=",
            0.62,
            "fail",
            "annual mountain precipitation should be wetter on windward slopes",
        ),
        _check(
            label,
            "seasonal_windward_leeward_precip_ratio_median",
            row["seasonal_windward_leeward_precip_ratio_median"],
            earth["seasonal_windward_leeward_precip_ratio_median"],
            "ratio>=",
            0.62,
            "fail",
            "seasonal windward/leeward precipitation contrast is too weak",
        ),
        _check(
            label,
            "seasonal_windward_leeward_precip_ratio_max",
            row["seasonal_windward_leeward_precip_ratio_max"],
            earth["seasonal_windward_leeward_precip_ratio_max"],
            "ratio>=",
            0.62,
            "warn",
            "no season develops a strong orographic wet/dry side contrast",
        ),
        _check(
            label,
            "seasonal_valid_count",
            row["seasonal_valid_count"],
            earth["seasonal_valid_count"],
            ">=",
            3.0,
            "warn",
            "too few seasons have enough windward and leeward mountain-slope area",
        ),
    ]


def _write_csv(path: Path, rows: list[dict[str, Any]], keys: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in keys})
    return path


def _render_markdown(report: dict[str, Any]) -> str:
    earth = report["earth_reference_metrics"]
    lines = [
        "# Earth Climate Windward/Leeward Gate",
        "",
        f"Verdict: `{report['verdict']}`",
        f"Failures: `{report['failure_count']}`",
        f"Warnings: `{report['warning_count']}`",
        "",
        "## Earth Envelope",
        "",
        f"- Mountain-slope land fraction: `{earth['mountain_slope_land_fraction']:.3f}`",
        f"- Annual windward/leeward precip ratio: `{earth['annual_windward_leeward_precip_ratio']:.3f}`",
        f"- Seasonal median windward/leeward precip ratio: `{earth['seasonal_windward_leeward_precip_ratio_median']:.3f}`",
        f"- Seasonal max windward/leeward precip ratio: `{earth['seasonal_windward_leeward_precip_ratio_max']:.3f}`",
        "",
        "## Failed / Warning Checks",
        "",
    ]
    for row in report["checks"]:
        if row["passed"]:
            continue
        status = "FAIL" if row["severity"] == "fail" else "WARN"
        lines.append(
            f"- {status} `{row['label']}` `{row['metric']}` generated "
            f"`{row['generated']:.3f}`, ref `{row['reference']:.3f}`, "
            f"ratio `{row['ratio_to_reference']:.3f}`: {row['message']}"
        )
    lines.extend(["", "## Generated Earthlike Metrics", ""])
    for row in report["generated_metrics"]:
        if row.get("mode") != "earthlike_calibration":
            continue
        lines.extend([
            f"### {row['label']}",
            "",
            f"- Annual windward/leeward precip ratio: `{row['annual_windward_leeward_precip_ratio']:.3f}`",
            f"- Seasonal median/max ratio: `{row['seasonal_windward_leeward_precip_ratio_median']:.3f}` / `{row['seasonal_windward_leeward_precip_ratio_max']:.3f}`",
            f"- Mountain-slope land fraction: `{row['mountain_slope_land_fraction']:.3f}`",
            f"- Orographic precip p75 on mountain slopes: `{row['orographic_precip_p75_mountain_slope_mm_yr']:.3f}`",
            "",
        ])
    return "\n".join(lines)


def run_earth_climate_windward_leeward_gate(
    config: EarthClimateWindwardLeewardGateConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    terminal = json.loads(Path(config.terminal_summary_json).read_text())
    earth = _earth_metrics(Path(config.earth_reference_npz))
    generated = [_generated_metrics(row) for row in terminal.get("summaries", [])]
    generated.sort(key=lambda row: (
        row.get("mode") != "earthlike_calibration",
        str(row.get("label")),
    ))
    checks: list[dict[str, Any]] = []
    for row in generated:
        if row["mode"] == "earthlike_calibration":
            checks.extend(_earthlike_checks(row, earth))

    failures = [row for row in checks if not row["passed"] and row["severity"] == "fail"]
    warnings = [row for row in checks if not row["passed"] and row["severity"] == "warn"]
    skipped = [row for row in checks if row.get("skipped")]
    if failures:
        verdict = "fail"
    elif warnings:
        verdict = "pass_with_warnings"
    else:
        verdict = "pass"

    metric_keys = [
        "label", "preset", "seed", "mode", "land_fraction",
        "mountain_slope_land_fraction", "strong_slope_land_fraction",
        "annual_windward_area_fraction_of_land",
        "annual_leeward_area_fraction_of_land",
        "annual_windward_precip_mean_mm_yr",
        "annual_leeward_precip_mean_mm_yr",
        "annual_windward_leeward_precip_ratio",
        "seasonal_valid_count",
        "seasonal_windward_leeward_precip_ratio_min",
        "seasonal_windward_leeward_precip_ratio_median",
        "seasonal_windward_leeward_precip_ratio_max",
        "seasonal_windward_precip_mean_median_mm_yr",
        "seasonal_leeward_precip_mean_median_mm_yr",
        "barrier_p75_mountain_slope",
        "orographic_precip_p75_mountain_slope_mm_yr",
    ]
    metrics_csv = _write_csv(
        outdir / "earth_climate_windward_leeward_metrics.csv",
        [earth | {"mode": "earth_reference"}] + generated,
        metric_keys,
    )
    checks_csv = _write_csv(
        outdir / "earth_climate_windward_leeward_checks.csv",
        checks,
        [
            "label", "metric", "generated", "reference", "ratio_to_reference",
            "operator", "threshold", "severity", "passed", "skipped", "message",
        ],
    )
    report = {
        "schema": SCHEMA,
        "earth_reference_npz": str(config.earth_reference_npz),
        "terminal_summary_json": str(config.terminal_summary_json),
        "verdict": verdict,
        "failure_count": len(failures),
        "warning_count": len(warnings),
        "skipped_count": len(skipped),
        "earth_reference_metrics": earth,
        "generated_metrics": generated,
        "checks": checks,
        "failures": failures,
        "warnings": warnings,
        "skipped": skipped,
        "metrics_csv": str(metrics_csv),
        "checks_csv": str(checks_csv),
    }
    md_path = outdir / "earth_climate_windward_leeward_gate_report.md"
    report["markdown"] = str(md_path)
    (outdir / "earth_climate_windward_leeward_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default)
    )
    md_path.write_text(_render_markdown(report))
    return report
