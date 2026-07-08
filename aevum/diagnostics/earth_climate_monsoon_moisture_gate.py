"""C4a monsoon/moisture diagnostic gate against broad Earth envelopes."""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "aevum.earth_climate_monsoon_moisture_gate.v1"


@dataclass(frozen=True)
class EarthClimateMonsoonMoistureGateConfig:
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


def _safe_float(value: Any, default: float = float("nan")) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def _percentile(values: np.ndarray, q: float, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    if not mask.any():
        return float("nan")
    return float(np.nanpercentile(values[mask], q))


def _weighted_mean(values: np.ndarray, area: np.ndarray, mask: np.ndarray) -> float:
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


def _summer(lat: np.ndarray, seasonal: np.ndarray) -> np.ndarray:
    return np.where(np.asarray(lat) >= 0.0, seasonal[2], seasonal[0])


def _winter(lat: np.ndarray, seasonal: np.ndarray) -> np.ndarray:
    return np.where(np.asarray(lat) >= 0.0, seasonal[0], seasonal[2])


def _monsoon_like_mask(
    lat: np.ndarray,
    land: np.ndarray,
    annual_precip: np.ndarray,
    seasonal_precip: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    low_mid = (
        np.asarray(land, dtype=bool)
        & (np.abs(np.asarray(lat, dtype=np.float64)) >= 5.0)
        & (np.abs(np.asarray(lat, dtype=np.float64)) <= 35.0)
    )
    summer = _summer(lat, seasonal_precip)
    winter = _winter(lat, seasonal_precip)
    diff = summer - winter
    ratio = (summer + 25.0) / np.maximum(winter + 25.0, 1.0e-9)
    monsoon_like = (
        low_mid
        & (np.asarray(annual_precip, dtype=np.float64) >= 400.0)
        & (diff >= 250.0)
        & (ratio >= 1.35)
    )
    return low_mid, monsoon_like, diff


def _earth_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        lat = np.asarray(z["lat"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        land = np.asarray(z["earth__land_mask"], dtype=bool)
        annual_precip = np.asarray(z["earth__annual_precip_mm"], dtype=np.float64)
        seasonal_precip = np.asarray(
            z["earth__seasonal_precip_mm_yr_equiv"],
            dtype=np.float64,
        )
        slp = np.asarray(z["earth__seasonal_slp_anomaly_hPa"], dtype=np.float64)

    low_mid, monsoon_like, diff = _monsoon_like_mask(
        lat, land, annual_precip, seasonal_precip)
    all_cells = np.ones_like(land, dtype=bool)
    return {
        "label": "earth_reference",
        "preset": "earth_reference",
        "seed": 0,
        "land_fraction": _area_fraction(area, land, all_cells),
        "low_mid_land_share": _area_fraction(area, low_mid, land),
        "monsoon_like_fraction_low_mid": _area_fraction(area, monsoon_like, low_mid),
        "precip_summer_minus_winter_p75_mm_yr": _percentile(diff, 75, low_mid),
        "precip_summer_minus_winter_p90_mm_yr": _percentile(diff, 90, low_mid),
        "monsoon_like_slp_summer_mean_hPa": _weighted_mean(
            _summer(lat, slp), area, monsoon_like),
        "monsoon_like_slp_winter_mean_hPa": _weighted_mean(
            _winter(lat, slp), area, monsoon_like),
        "dry_low_mid_fraction": _area_fraction(
            area, low_mid & (annual_precip < 400.0), low_mid),
    }


def _generated_row(summary_row: dict[str, Any]) -> dict[str, Any]:
    arrays_path = Path(summary_row["assets_dir"]) / "terminal_climate_arrays.npz"
    label = Path(summary_row["assets_dir"]).name
    with np.load(arrays_path, allow_pickle=False) as z:
        lat = np.asarray(z["lat"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        elev = np.asarray(z["terrain__elevation_m"], dtype=np.float64)
        land = elev >= sea
        annual_precip = np.asarray(z["climate__precipitation"], dtype=np.float64)
        seasonal_precip = np.asarray(
            z["climate__seasonal_precipitation"],
            dtype=np.float64,
        )
        moisture = np.asarray(z["atmosphere__moisture_access"], dtype=np.float64)
        monsoon = np.asarray(z["atmosphere__monsoon_potential"], dtype=np.float64)
        pressure = np.asarray(
            z["atmosphere__seasonal_pressure_proxy"],
            dtype=np.float64,
        )

    low_mid, monsoon_like, diff = _monsoon_like_mask(
        lat, land, annual_precip, seasonal_precip)
    summer_potential = _summer(lat, monsoon)
    winter_potential = _winter(lat, monsoon)
    summer_access = _summer(lat, moisture)
    summer_pressure = _summer(lat, pressure)
    winter_pressure = _winter(lat, pressure)
    all_cells = np.ones_like(land, dtype=bool)
    land4 = np.broadcast_to(land, moisture.shape)
    return {
        "label": label,
        "preset": str(summary_row.get("preset", "")),
        "seed": int(summary_row.get("seed", -1)),
        "arrays": str(arrays_path),
        "land_fraction": _area_fraction(area, land, all_cells),
        "low_mid_land_share": _area_fraction(area, low_mid, land),
        "monsoon_like_fraction_low_mid": _area_fraction(area, monsoon_like, low_mid),
        "precip_summer_minus_winter_p75_mm_yr": _percentile(diff, 75, low_mid),
        "precip_summer_minus_winter_p90_mm_yr": _percentile(diff, 90, low_mid),
        "monsoon_potential_summer_p75_low_mid": _percentile(
            summer_potential, 75, low_mid),
        "monsoon_potential_summer_p90_low_mid": _percentile(
            summer_potential, 90, low_mid),
        "monsoon_potential_winter_p75_low_mid": _percentile(
            winter_potential, 75, low_mid),
        "monsoon_potential_summer_minus_winter_p75_low_mid": _percentile(
            summer_potential - winter_potential, 75, low_mid),
        "monsoon_potential_monsoon_like_p75": _percentile(
            summer_potential, 75, monsoon_like),
        "monsoon_potential_land_p95": _percentile(np.max(monsoon, axis=0), 95, land),
        "moisture_access_summer_p75_low_mid": _percentile(
            summer_access, 75, low_mid),
        "moisture_access_summer_p90_low_mid": _percentile(
            summer_access, 90, low_mid),
        "moisture_access_monsoon_like_p75": _percentile(
            summer_access, 75, monsoon_like),
        "moisture_access_land_allseason_p75": _percentile(moisture, 75, land4),
        "pressure_summer_mean_low_mid": _weighted_mean(
            summer_pressure, area, low_mid),
        "pressure_winter_mean_low_mid": _weighted_mean(
            winter_pressure, area, low_mid),
        "pressure_summer_minus_winter_mean_low_mid": _weighted_mean(
            summer_pressure - winter_pressure, area, low_mid),
    }


def _ratio(value: float, ref: float) -> float:
    if not np.isfinite(value) or not np.isfinite(ref) or abs(ref) <= 1.0e-12:
        return float("nan")
    return float(value / ref)


def _check(
    *,
    label: str,
    group: str,
    metric: str,
    generated: float,
    earth: float = float("nan"),
    operator: str,
    threshold: float,
    severity: str,
    message: str,
) -> dict[str, Any]:
    generated = _safe_float(generated)
    earth = _safe_float(earth)
    ratio = _ratio(generated, earth)
    skipped = not np.isfinite(generated)
    if skipped:
        passed = True
    elif operator == ">=":
        passed = generated >= threshold
    elif operator == "<=":
        passed = generated <= threshold
    elif operator == "ratio>=":
        passed = np.isfinite(ratio) and ratio >= threshold
    elif operator == "ratio<=":
        passed = np.isfinite(ratio) and ratio <= threshold
    else:
        raise ValueError(f"unsupported operator {operator!r}")
    return {
        "label": label,
        "group": group,
        "metric": metric,
        "generated": generated,
        "earth": earth,
        "ratio_to_earth": ratio,
        "operator": operator,
        "threshold": threshold,
        "severity": severity,
        "passed": bool(passed),
        "skipped": bool(skipped),
        "message": message,
    }


def _earthlike_checks(row: dict[str, Any], earth: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(
            label=label,
            group="earthlike_monsoon_envelope",
            metric="monsoon_like_fraction_low_mid",
            generated=row["monsoon_like_fraction_low_mid"],
            earth=earth["monsoon_like_fraction_low_mid"],
            operator="ratio>=",
            threshold=0.45,
            severity="fail",
            message="earthlike low/mid-latitude land needs a broad summer-wet envelope",
        ),
        _check(
            label=label,
            group="earthlike_monsoon_potential",
            metric="monsoon_potential_summer_p90_low_mid",
            generated=row["monsoon_potential_summer_p90_low_mid"],
            operator=">=",
            threshold=0.20,
            severity="fail",
            message="summer heated continents should diagnose strong monsoon potential",
        ),
        _check(
            label=label,
            group="earthlike_monsoon_potential",
            metric="monsoon_potential_summer_minus_winter_p75_low_mid",
            generated=row["monsoon_potential_summer_minus_winter_p75_low_mid"],
            operator=">=",
            threshold=0.18,
            severity="fail",
            message="monsoon potential should be summer-enhanced, not seasonless",
        ),
        _check(
            label=label,
            group="earthlike_moisture_corridors",
            metric="moisture_access_monsoon_like_p75",
            generated=row["moisture_access_monsoon_like_p75"],
            operator=">=",
            threshold=0.55,
            severity="fail",
            message="summer-wet monsoon corridors should have high moisture access",
        ),
        _check(
            label=label,
            group="earthlike_pressure",
            metric="pressure_summer_minus_winter_mean_low_mid",
            generated=row["pressure_summer_minus_winter_mean_low_mid"],
            operator="<=",
            threshold=-0.45,
            severity="fail",
            message="low/mid-latitude continents should have summer heat-low pressure",
        ),
    ]


def _waterworld_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(
            label=label,
            group="waterworld_fake_monsoon_guard",
            metric="monsoon_potential_land_p95",
            generated=row["monsoon_potential_land_p95"],
            operator="<=",
            threshold=0.24,
            severity="fail",
            message="small islands on waterworlds should not become continent-scale monsoons",
        ),
        _check(
            label=label,
            group="waterworld_fake_monsoon_guard",
            metric="monsoon_potential_summer_p90_low_mid",
            generated=row["monsoon_potential_summer_p90_low_mid"],
            operator="<=",
            threshold=0.22,
            severity="fail",
            message="waterworld low/mid-latitude land should not diagnose broad monsoon potential",
        ),
        _check(
            label=label,
            group="waterworld_pressure_guard",
            metric="pressure_summer_minus_winter_mean_low_mid",
            generated=row["pressure_summer_minus_winter_mean_low_mid"],
            operator=">=",
            threshold=-0.55,
            severity="fail",
            message="waterworld islands should have weak continent-driven pressure reversal",
        ),
    ]


def _arid_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(
            label=label,
            group="arid_moisture_guard",
            metric="moisture_access_land_allseason_p75",
            generated=row["moisture_access_land_allseason_p75"],
            operator="<=",
            threshold=0.30,
            severity="fail",
            message="arid large-continent worlds should retain moisture-limited interiors",
        ),
        _check(
            label=label,
            group="arid_thermal_pressure",
            metric="pressure_summer_minus_winter_mean_low_mid",
            generated=row["pressure_summer_minus_winter_mean_low_mid"],
            operator="<=",
            threshold=-0.25,
            severity="fail",
            message="arid large continents should still diagnose seasonal thermal pressure",
        ),
        _check(
            label=label,
            group="arid_fake_monsoon_guard",
            metric="monsoon_potential_land_p95",
            generated=row["monsoon_potential_land_p95"],
            operator="<=",
            threshold=0.32,
            severity="fail",
            message="arid worlds should not become broadly monsoonal despite thermal contrast",
        ),
    ]


def _checks_for_row(row: dict[str, Any], earth: dict[str, Any]) -> list[dict[str, Any]]:
    preset = str(row.get("preset", "")).lower()
    label = str(row.get("label", "")).lower()
    if "earthlike" in preset or "earthlike" in label:
        return _earthlike_checks(row, earth)
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
        "# Earth Climate Monsoon/Moisture Gate",
        "",
        f"Schema: `{report['schema']}`",
        f"Verdict: **{report['verdict']}**",
        "",
        (
            f"Failures: `{report['failure_count']}`; warnings: "
            f"`{report['warning_count']}`; skipped: `{report['skipped_count']}`."
        ),
        "",
        "This gate targets C4a geography-derived seasonal pressure, moisture "
        "access, and monsoon potential.  It uses real Earth seasonal "
        "precipitation and pressure as broad envelopes, then applies preset "
        "guards so waterworld islands do not become artificial continental "
        "monsoons and arid worlds remain moisture-limited.",
        "",
        "## Earth Reference",
        "",
    ]
    for key, value in report["earth"].items():
        if key in {"label", "preset", "seed"}:
            continue
        lines.append(f"- `{key}`: `{_safe_float(value):.3f}`")
    lines.extend(["", "## Checks", ""])
    for row in report["checks"]:
        status = "pass" if row["passed"] else row["severity"]
        lines.append(
            f"- `{status}` `{row['label']}` `{row['group']}` "
            f"`{row['metric']}` = `{row['generated']:.3f}` "
            f"{row['operator']} `{row['threshold']:.3f}`"
        )
    lines.append("")
    return "\n".join(lines)


def run_earth_climate_monsoon_moisture_gate(
    config: EarthClimateMonsoonMoistureGateConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    earth = _earth_metrics(Path(config.earth_reference_npz))
    summary = json.loads(Path(config.terminal_summary_json).read_text())
    generated = [_generated_row(row) for row in summary.get("summaries", [])]
    checks = [
        check
        for row in generated
        for check in _checks_for_row(row, earth)
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
        "earth_reference_npz": str(config.earth_reference_npz),
        "terminal_summary_json": str(config.terminal_summary_json),
        "earth": earth,
        "metrics": generated,
        "checks": checks,
        "failures": failures,
        "warnings": warnings,
        "failure_count": len(failures),
        "warning_count": len(warnings),
        "skipped_count": len(skipped),
        "verdict": "fail" if failures else "pass",
    }
    metric_keys = sorted({key for row in [earth, *generated] for key in row.keys()})
    check_keys = [
        "label", "group", "metric", "generated", "earth", "ratio_to_earth",
        "operator", "threshold", "severity", "passed", "skipped", "message",
    ]
    _write_csv(outdir / "earth_climate_monsoon_moisture_metrics.csv",
               [earth, *generated], metric_keys)
    _write_csv(outdir / "earth_climate_monsoon_moisture_checks.csv",
               checks, check_keys)
    (outdir / "earth_climate_monsoon_moisture_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default),
    )
    (outdir / "earth_climate_monsoon_moisture_gate_report.md").write_text(
        _render_markdown(report),
    )
    return report

