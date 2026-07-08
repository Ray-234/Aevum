"""Koppen-like seasonal subtype gate for generated terminal worlds."""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "aevum.earth_climate_seasonal_subtype_gate.v1"


@dataclass(frozen=True)
class EarthClimateSeasonalSubtypeGateConfig:
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


def _area_fraction(area: np.ndarray, mask: np.ndarray, denom_mask: np.ndarray) -> float:
    area = np.asarray(area, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    denom_mask = np.asarray(denom_mask, dtype=bool)
    denom = float(np.nansum(area[denom_mask]))
    if denom <= 0.0:
        return float("nan")
    return float(np.nansum(area[mask & denom_mask]) / denom)


def _percentile(values: np.ndarray, q: float, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    if not mask.any():
        return float("nan")
    return float(np.nanpercentile(values[mask], q))


def _band(lat: np.ndarray, lo: float, hi: float) -> np.ndarray:
    alat = np.abs(np.asarray(lat, dtype=np.float64))
    return (alat >= lo) & (alat < hi)


def _dry_quarter_count(
    seasonal_precip: np.ndarray,
    annual_precip: np.ndarray,
) -> np.ndarray:
    seasonal_precip = np.asarray(seasonal_precip, dtype=np.float64)
    annual_precip = np.asarray(annual_precip, dtype=np.float64)
    threshold = np.minimum(350.0, 0.45 * np.maximum(annual_precip, 1.0e-9))
    return np.sum(seasonal_precip < threshold[None, :], axis=0)


def _peak_share(area: np.ndarray, mask: np.ndarray, peak: np.ndarray, season: int) -> float:
    return _area_fraction(area, mask & (peak == season), mask)


def _seasonal_metrics(
    label: str,
    preset: str,
    seed: int,
    lat: np.ndarray,
    area: np.ndarray,
    land: np.ndarray,
    annual_precip: np.ndarray,
    seasonal_precip: np.ndarray,
) -> dict[str, Any]:
    lat = np.asarray(lat, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)
    land = np.asarray(land, dtype=bool)
    annual_precip = np.asarray(annual_precip, dtype=np.float64)
    seasonal_precip = np.asarray(seasonal_precip, dtype=np.float64)

    low_tropics = land & _band(lat, 0.0, 15.0)
    subtropics = land & _band(lat, 15.0, 30.0)
    low_mid = land & _band(lat, 0.0, 35.0)
    warm_midlat = land & _band(lat, 30.0, 45.0)
    cool_midlat = land & _band(lat, 45.0, 60.0)
    high_lat = land & _band(lat, 60.0, 90.0)

    dry_quarters = _dry_quarter_count(seasonal_precip, annual_precip)
    valid_seasonal = np.isfinite(seasonal_precip).any(axis=0)
    seasonal_max = np.full(annual_precip.shape, np.nan, dtype=np.float64)
    seasonal_min = np.full(annual_precip.shape, np.nan, dtype=np.float64)
    peak = np.zeros(annual_precip.shape, dtype=np.int16)
    if valid_seasonal.any():
        valid = seasonal_precip[:, valid_seasonal]
        seasonal_max[valid_seasonal] = np.max(
            np.where(np.isfinite(valid), valid, -np.inf), axis=0)
        seasonal_min[valid_seasonal] = np.min(
            np.where(np.isfinite(valid), valid, np.inf), axis=0)
        peak[valid_seasonal] = np.argmax(
            np.where(np.isfinite(valid), valid, -np.inf), axis=0).astype(np.int16)
    amp = (seasonal_max - seasonal_min) / np.maximum(annual_precip, 1.0e-9)

    row: dict[str, Any] = {
        "label": label,
        "preset": preset,
        "seed": seed,
        "land_fraction": _area_fraction(area, land, np.ones_like(land, dtype=bool)),
        "low_tropics_dry_quarter_ge1_fraction": _area_fraction(
            area, low_tropics & (dry_quarters >= 1), low_tropics),
        "low_tropics_dry_quarter_ge2_fraction": _area_fraction(
            area, low_tropics & (dry_quarters >= 2), low_tropics),
        "low_tropics_dry_quarter_ge3_fraction": _area_fraction(
            area, low_tropics & (dry_quarters >= 3), low_tropics),
        "low_tropics_amp_p50": _percentile(amp, 50, low_tropics),
        "low_tropics_peak_dominance": max(
            _peak_share(area, low_tropics, peak, season) for season in range(4)
        ),
        "subtropical_dry_quarter_ge1_fraction": _area_fraction(
            area, subtropics & (dry_quarters >= 1), subtropics),
        "subtropical_dry_quarter_ge2_fraction": _area_fraction(
            area, subtropics & (dry_quarters >= 2), subtropics),
        "low_mid_dry_quarter_ge1_fraction": _area_fraction(
            area, low_mid & (dry_quarters >= 1), low_mid),
        "low_mid_dry_quarter_ge2_fraction": _area_fraction(
            area, low_mid & (dry_quarters >= 2), low_mid),
        "low_mid_amp_p75": _percentile(amp, 75, low_mid),
        "warm_midlat_dry_quarter_ge2_fraction": _area_fraction(
            area, warm_midlat & (dry_quarters >= 2), warm_midlat),
        "cool_midlat_amp_p50": _percentile(amp, 50, cool_midlat),
        "high_lat_dry_quarter_ge3_fraction": _area_fraction(
            area, high_lat & (dry_quarters >= 3), high_lat),
        "high_lat_amp_p50": _percentile(amp, 50, high_lat),
    }
    for group_name, group_mask in [
        ("low_tropics", low_tropics),
        ("subtropics", subtropics),
        ("cool_midlat", cool_midlat),
    ]:
        for season in range(4):
            row[f"{group_name}_peak_s{season}_fraction"] = _peak_share(
                area, group_mask, peak, season)
    return row


def _earth_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        return _seasonal_metrics(
            "Earth seasonal reference",
            "earth_reference",
            -1,
            np.asarray(z["lat"], dtype=np.float64),
            np.asarray(z["cell_area"], dtype=np.float64),
            np.asarray(z["earth__land_mask"], dtype=bool),
            np.asarray(z["earth__annual_precip_mm"], dtype=np.float64),
            np.asarray(z["earth__seasonal_precip_mm_yr_equiv"], dtype=np.float64),
        )


def _generated_metrics(summary: dict[str, Any]) -> dict[str, Any]:
    arrays = Path(summary["arrays"])
    with np.load(arrays, allow_pickle=False) as z:
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        land = np.asarray(z["terrain__elevation_m"], dtype=np.float64) >= sea
        row = _seasonal_metrics(
            Path(str(summary.get("assets_dir", ""))).name,
            str(summary.get("preset", "")),
            int(summary.get("seed", -1)),
            np.asarray(z["lat"], dtype=np.float64),
            np.asarray(z["cell_area"], dtype=np.float64),
            land,
            np.asarray(z["climate__precipitation"], dtype=np.float64),
            np.asarray(z["climate__seasonal_precipitation"], dtype=np.float64),
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
    elif operator == "ratio<=":
        passed = np.isfinite(ratio) and ratio <= threshold
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
            "low_tropics_dry_quarter_ge2_fraction",
            row["low_tropics_dry_quarter_ge2_fraction"],
            earth["low_tropics_dry_quarter_ge2_fraction"],
            "ratio>=",
            0.45,
            "fail",
            "low-tropical land lacks enough seasonal-dry/monsoon subtype area",
        ),
        _check(
            label,
            "low_tropics_dry_quarter_ge3_fraction",
            row["low_tropics_dry_quarter_ge3_fraction"],
            earth["low_tropics_dry_quarter_ge3_fraction"],
            "<=",
            0.14,
            "warn",
            "low-tropical land should not become mostly long-dry-season climate",
        ),
        _check(
            label,
            "low_mid_dry_quarter_ge1_fraction",
            row["low_mid_dry_quarter_ge1_fraction"],
            earth["low_mid_dry_quarter_ge1_fraction"],
            "ratio>=",
            0.75,
            "fail",
            "low/mid-latitude land lacks broad dry-season subtype expression",
        ),
        _check(
            label,
            "low_mid_dry_quarter_ge1_fraction",
            row["low_mid_dry_quarter_ge1_fraction"],
            earth["low_mid_dry_quarter_ge1_fraction"],
            "ratio<=",
            1.55,
            "warn",
            "low/mid-latitude dry-season subtype area is too dominant",
        ),
        _check(
            label,
            "subtropical_dry_quarter_ge1_fraction",
            row["subtropical_dry_quarter_ge1_fraction"],
            earth["subtropical_dry_quarter_ge1_fraction"],
            "ratio>=",
            0.70,
            "fail",
            "subtropical seasonal dryness should be visible",
        ),
        _check(
            label,
            "subtropical_dry_quarter_ge1_fraction",
            row["subtropical_dry_quarter_ge1_fraction"],
            earth["subtropical_dry_quarter_ge1_fraction"],
            "ratio<=",
            1.30,
            "warn",
            "subtropical seasonal dryness is too widespread",
        ),
        _check(
            label,
            "low_mid_amp_p75",
            row["low_mid_amp_p75"],
            earth["low_mid_amp_p75"],
            "ratio>=",
            0.65,
            "fail",
            "low/mid-latitude precipitation seasonality amplitude is too weak",
        ),
        _check(
            label,
            "low_mid_amp_p75",
            row["low_mid_amp_p75"],
            earth["low_mid_amp_p75"],
            "ratio<=",
            1.45,
            "warn",
            "low/mid-latitude precipitation seasonality amplitude is too strong",
        ),
        _check(
            label,
            "cool_midlat_amp_p50",
            row["cool_midlat_amp_p50"],
            earth["cool_midlat_amp_p50"],
            "ratio>=",
            0.55,
            "fail",
            "cool-midlatitude seasonality is too flat",
        ),
        _check(
            label,
            "cool_midlat_amp_p50",
            row["cool_midlat_amp_p50"],
            earth["cool_midlat_amp_p50"],
            "ratio<=",
            1.80,
            "warn",
            "cool-midlatitude seasonality is too peaked",
        ),
        _check(
            label,
            "high_lat_dry_quarter_ge3_fraction",
            row["high_lat_dry_quarter_ge3_fraction"],
            earth["high_lat_dry_quarter_ge3_fraction"],
            "<=",
            0.60,
            "warn",
            "high-latitude precipitation seasonality should not become uniformly long-dry-season",
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
        "# Earth Climate Seasonal Subtype Gate",
        "",
        f"Verdict: `{report['verdict']}`",
        f"Failures: `{report['failure_count']}`",
        f"Warnings: `{report['warning_count']}`",
        "",
        "## Earth Seasonal Envelope",
        "",
        f"- Low-tropics dry-quarter >=2: `{earth['low_tropics_dry_quarter_ge2_fraction']:.3f}`",
        f"- Low/mid dry-quarter >=1: `{earth['low_mid_dry_quarter_ge1_fraction']:.3f}`",
        f"- Subtropical dry-quarter >=1: `{earth['subtropical_dry_quarter_ge1_fraction']:.3f}`",
        f"- Low/mid amplitude p75: `{earth['low_mid_amp_p75']:.3f}`",
        f"- Cool-midlat amplitude p50: `{earth['cool_midlat_amp_p50']:.3f}`",
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
    lines.extend(["", "## Generated Earthlike Seasonal Metrics", ""])
    for row in report["generated_metrics"]:
        if row.get("mode") != "earthlike_calibration":
            continue
        lines.extend([
            f"### {row['label']}",
            "",
            f"- Low-tropics dry-quarter >=2 / >=3: `{row['low_tropics_dry_quarter_ge2_fraction']:.3f}` / `{row['low_tropics_dry_quarter_ge3_fraction']:.3f}`",
            f"- Low/mid dry-quarter >=1 / >=2: `{row['low_mid_dry_quarter_ge1_fraction']:.3f}` / `{row['low_mid_dry_quarter_ge2_fraction']:.3f}`",
            f"- Subtropical dry-quarter >=1 / >=2: `{row['subtropical_dry_quarter_ge1_fraction']:.3f}` / `{row['subtropical_dry_quarter_ge2_fraction']:.3f}`",
            f"- Low/mid amp p75: `{row['low_mid_amp_p75']:.3f}`",
            f"- Cool-midlat amp p50: `{row['cool_midlat_amp_p50']:.3f}`",
            "",
        ])
    return "\n".join(lines)


def run_earth_climate_seasonal_subtype_gate(
    config: EarthClimateSeasonalSubtypeGateConfig,
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
        "low_tropics_dry_quarter_ge1_fraction",
        "low_tropics_dry_quarter_ge2_fraction",
        "low_tropics_dry_quarter_ge3_fraction",
        "low_tropics_amp_p50", "low_tropics_peak_dominance",
        "subtropical_dry_quarter_ge1_fraction",
        "subtropical_dry_quarter_ge2_fraction",
        "low_mid_dry_quarter_ge1_fraction",
        "low_mid_dry_quarter_ge2_fraction",
        "low_mid_amp_p75",
        "warm_midlat_dry_quarter_ge2_fraction",
        "cool_midlat_amp_p50",
        "high_lat_dry_quarter_ge3_fraction",
        "high_lat_amp_p50",
    ]
    metrics_csv = _write_csv(
        outdir / "earth_climate_seasonal_subtype_metrics.csv",
        [earth | {"mode": "earth_reference"}] + generated,
        metric_keys,
    )
    checks_csv = _write_csv(
        outdir / "earth_climate_seasonal_subtype_checks.csv",
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
    md_path = outdir / "earth_climate_seasonal_subtype_gate_report.md"
    report["markdown"] = str(md_path)
    (outdir / "earth_climate_seasonal_subtype_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default)
    )
    md_path.write_text(_render_markdown(report))
    return report
