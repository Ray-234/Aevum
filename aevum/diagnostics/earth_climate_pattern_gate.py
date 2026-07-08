"""Earth-pattern climate gate for generated terminal worlds.

This diagnostic complements scalar Earth fitting.  It compares generated
Earthlike worlds against broad Earth envelopes by process/geographic class
rather than by exact longitude/latitude match.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "aevum.earth_climate_pattern_gate.v1"


@dataclass(frozen=True)
class EarthClimatePatternGateConfig:
    earth_reference_npz: Path
    terminal_summary_json: Path
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


def _safe_float(value: Any, default: float = float("nan")) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def _weighted_mean(values: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values) & np.isfinite(weights)
    if not mask.any():
        return float("nan")
    return float(np.average(values[mask], weights=weights[mask]))


def _percentile(values: np.ndarray, q: float, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    if not mask.any():
        return float("nan")
    return float(np.nanpercentile(values[mask], q))


def _area_fraction(area: np.ndarray, mask: np.ndarray, denom_mask: np.ndarray) -> float:
    area = np.asarray(area, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    denom_mask = np.asarray(denom_mask, dtype=bool)
    denom = float(np.nansum(area[denom_mask]))
    if denom <= 0.0:
        return float("nan")
    return float(np.nansum(area[mask & denom_mask]) / denom)


def _seasonal_monsoon_index(
    lat: np.ndarray,
    seasonal_precip: np.ndarray,
    annual_precip: np.ndarray,
) -> np.ndarray:
    lat = np.asarray(lat, dtype=np.float64)
    seasonal_precip = np.asarray(seasonal_precip, dtype=np.float64)
    annual_precip = np.asarray(annual_precip, dtype=np.float64)
    nh = lat >= 0.0
    summer = np.where(nh, seasonal_precip[2], seasonal_precip[0])
    winter = np.where(nh, seasonal_precip[0], seasonal_precip[2])
    with np.errstate(divide="ignore", invalid="ignore"):
        out = (summer - winter) / np.maximum(annual_precip, 1.0e-9)
    out[~np.isfinite(out)] = np.nan
    return out


def _pattern_metrics(
    *,
    label: str,
    preset: str,
    seed: int,
    arrays_path: Path,
    earth: bool,
) -> dict[str, Any]:
    with np.load(arrays_path, allow_pickle=False) as z:
        lat = np.asarray(z["lat"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        if earth:
            land = np.asarray(z["earth__land_mask"], dtype=bool)
            relative_elev = np.asarray(z["earth__elevation_m"], dtype=np.float64)
            temp_c = np.asarray(z["earth__annual_temperature_C"], dtype=np.float64)
            precip = np.asarray(z["earth__annual_precip_mm"], dtype=np.float64)
            seasonal_precip = np.asarray(
                z["earth__seasonal_precip_mm_yr_equiv"],
                dtype=np.float64,
            )
            biome = np.asarray(z["earth__biome_class_proxy"], dtype=np.int16)
        else:
            sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
            relative_elev = (
                np.asarray(z["terrain__elevation_m"], dtype=np.float64) - sea
            )
            land = relative_elev >= 0.0
            temp_c = (
                np.asarray(z["climate__surface_temperature"], dtype=np.float64)
                - 273.15
            )
            precip = np.asarray(z["climate__precipitation"], dtype=np.float64)
            seasonal_precip = np.asarray(
                z["climate__seasonal_precipitation"],
                dtype=np.float64,
            )
            biome = np.asarray(z["biosphere__biome"], dtype=np.float64).astype(np.int16)

    all_cells = np.ones_like(land, dtype=bool)
    tropical = land & (np.abs(lat) <= 23.5)
    low_mid = land & (np.abs(lat) <= 35.0)
    subtropics = land & (np.abs(lat) >= 18.0) & (np.abs(lat) <= 35.0)
    high_lat = land & (np.abs(lat) >= 55.0)
    mountain = land & (relative_elev >= 1000.0)
    warm_mountain = mountain & (np.abs(lat) <= 35.0)
    monsoon = _seasonal_monsoon_index(lat, seasonal_precip, precip)
    forest_or_tropical = land & ((biome == 4) | (biome == 6))
    ice_or_tundra = land & ((biome == 1) | (biome == 5))

    return {
        "label": label,
        "preset": preset,
        "seed": seed,
        "arrays": str(arrays_path),
        "land_fraction": _area_fraction(area, land, all_cells),
        "tropical_land_share": _area_fraction(area, tropical, land),
        "wet_tropics_precip_mean_mm_yr": _weighted_mean(precip, area, tropical),
        "wet_tropics_precip_p75_mm_yr": _percentile(precip, 75, tropical),
        "wet_tropics_precip_p90_mm_yr": _percentile(precip, 90, tropical),
        "wet_tropics_fraction_gt1000mm": _area_fraction(
            area, tropical & (precip >= 1000.0), tropical),
        "wet_tropics_fraction_gt1500mm": _area_fraction(
            area, tropical & (precip >= 1500.0), tropical),
        "low_mid_monsoon_p75": _percentile(monsoon, 75, low_mid),
        "low_mid_monsoon_p90": _percentile(monsoon, 90, low_mid),
        "low_mid_monsoon_wet_fraction_ge1": _area_fraction(
            area, low_mid & (monsoon >= 1.0) & (precip >= 500.0), low_mid),
        "dry_subtropics_fraction_lt250mm": _area_fraction(
            area, subtropics & (precip < 250.0), subtropics),
        "dry_subtropics_fraction_lt500mm": _area_fraction(
            area, subtropics & (precip < 500.0), subtropics),
        "dry_subtropics_precip_p25_mm_yr": _percentile(precip, 25, subtropics),
        "mountain_land_share": _area_fraction(area, mountain, land),
        "mountain_precip_p75_mm_yr": _percentile(precip, 75, mountain),
        "mountain_wet_fraction_gt1000mm": _area_fraction(
            area, mountain & (precip >= 1000.0), mountain),
        "warm_mountain_precip_p75_mm_yr": _percentile(precip, 75, warm_mountain),
        "high_lat_land_share": _area_fraction(area, high_lat, land),
        "high_lat_temperature_mean_C": _weighted_mean(temp_c, area, high_lat),
        "high_lat_temperature_p75_C": _percentile(temp_c, 75, high_lat),
        "high_lat_cold_fraction_lt0C": _area_fraction(
            area, high_lat & (temp_c < 0.0), high_lat),
        "high_lat_ice_tundra_fraction": _area_fraction(
            area, high_lat & ice_or_tundra, high_lat),
        "forest_tropical_land_fraction": _area_fraction(area, forest_or_tropical, land),
        "tropical_biome_land_fraction": _area_fraction(area, land & (biome == 6), land),
    }


def _ratio(generated: float, earth: float) -> float:
    if not np.isfinite(generated) or not np.isfinite(earth) or abs(earth) <= 1.0e-12:
        return float("nan")
    return float(generated / earth)


def _check(
    *,
    label: str,
    group: str,
    metric: str,
    generated: float,
    earth: float,
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


def _earthlike_pattern_checks(
    row: dict[str, Any],
    earth: dict[str, Any],
) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(
            label=label,
            group="wet_tropics",
            metric="wet_tropics_precip_p90_mm_yr",
            generated=row["wet_tropics_precip_p90_mm_yr"],
            earth=earth["wet_tropics_precip_p90_mm_yr"],
            operator="ratio>=",
            threshold=0.50,
            severity="fail",
            message="earthlike wet tropics need a credible high-rainfall tail",
        ),
        _check(
            label=label,
            group="wet_tropics",
            metric="wet_tropics_fraction_gt1000mm",
            generated=row["wet_tropics_fraction_gt1000mm"],
            earth=earth["wet_tropics_fraction_gt1000mm"],
            operator="ratio>=",
            threshold=0.35,
            severity="fail",
            message="earthlike wet tropics should contain broad >1000 mm/yr areas",
        ),
        _check(
            label=label,
            group="dry_subtropics",
            metric="dry_subtropics_fraction_lt250mm",
            generated=row["dry_subtropics_fraction_lt250mm"],
            earth=earth["dry_subtropics_fraction_lt250mm"],
            operator=">=",
            threshold=0.18,
            severity="fail",
            message="earthlike subtropical dry belts are too weak or absent",
        ),
        _check(
            label=label,
            group="monsoon_margins",
            metric="low_mid_monsoon_wet_fraction_ge1",
            generated=row["low_mid_monsoon_wet_fraction_ge1"],
            earth=earth["low_mid_monsoon_wet_fraction_ge1"],
            operator="<=",
            threshold=0.75,
            severity="warn",
            message="monsoon-season wet area should not dominate all low-mid land",
        ),
        _check(
            label=label,
            group="windward_mountains",
            metric="mountain_wet_fraction_gt1000mm",
            generated=row["mountain_wet_fraction_gt1000mm"],
            earth=earth["mountain_wet_fraction_gt1000mm"],
            operator="ratio>=",
            threshold=0.25,
            severity="warn",
            message="mountain/windward wet tail is weak relative to Earth",
        ),
        _check(
            label=label,
            group="high_latitudes",
            metric="high_lat_cold_fraction_lt0C",
            generated=row["high_lat_cold_fraction_lt0C"],
            earth=earth["high_lat_cold_fraction_lt0C"],
            operator=">=",
            threshold=0.65,
            severity="fail",
            message="high-latitude earthlike land should have a broad cold envelope",
        ),
        _check(
            label=label,
            group="high_latitudes",
            metric="high_lat_ice_tundra_fraction",
            generated=row["high_lat_ice_tundra_fraction"],
            earth=earth["high_lat_ice_tundra_fraction"],
            operator=">=",
            threshold=0.20,
            severity="fail",
            message="cold high-latitude land is not being expressed as tundra/ice biome",
        ),
        _check(
            label=label,
            group="biome_envelope",
            metric="forest_tropical_land_fraction",
            generated=row["forest_tropical_land_fraction"],
            earth=earth["forest_tropical_land_fraction"],
            operator="ratio>=",
            threshold=0.40,
            severity="warn",
            message="forest+tropical biome area remains low versus Earth envelope",
        ),
        _check(
            label=label,
            group="biome_envelope",
            metric="tropical_biome_land_fraction",
            generated=row["tropical_biome_land_fraction"],
            earth=earth["tropical_biome_land_fraction"],
            operator="ratio>=",
            threshold=0.35,
            severity="warn",
            message="tropical biome area remains low versus Earth envelope",
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
    lines = [
        "# Earth Climate Pattern Gate",
        "",
        f"Earth reference: `{report['earth_reference_npz']}`",
        f"Terminal summary: `{report['terminal_summary_json']}`",
        "",
        f"Verdict: `{report['verdict']}`",
        f"Failures: `{report['failure_count']}`",
        f"Warnings: `{report['warning_count']}`",
        "",
        "## Earth Baseline",
        "",
    ]
    earth = report["earth_metrics"]
    for key in [
        "wet_tropics_precip_p90_mm_yr",
        "wet_tropics_fraction_gt1000mm",
        "dry_subtropics_fraction_lt250mm",
        "low_mid_monsoon_wet_fraction_ge1",
        "mountain_wet_fraction_gt1000mm",
        "high_lat_cold_fraction_lt0C",
        "high_lat_ice_tundra_fraction",
        "forest_tropical_land_fraction",
    ]:
        value = earth.get(key, float("nan"))
        lines.append(f"- `{key}`: `{value:.3f}`")
    lines.extend(["", "## Checks", ""])
    for row in report["checks"]:
        if row["passed"]:
            continue
        status = "FAIL" if row["severity"] == "fail" else "WARN"
        lines.append(
            f"- {status} `{row['label']}` `{row['metric']}` "
            f"generated `{row['generated']:.3f}`, Earth `{row['earth']:.3f}`, "
            f"ratio `{row['ratio_to_earth']:.3f}`: {row['message']}"
        )
    lines.extend(["", "## Earthlike Metrics", ""])
    for row in report["generated_metrics"]:
        if row.get("mode") != "earthlike_calibration":
            continue
        lines.extend([
            f"### {row['label']}",
            "",
            f"- Wet tropics p90 / Earth: `{_ratio(row['wet_tropics_precip_p90_mm_yr'], earth['wet_tropics_precip_p90_mm_yr']):.3f}`",
            f"- Wet tropics >1000 mm fraction / Earth: `{_ratio(row['wet_tropics_fraction_gt1000mm'], earth['wet_tropics_fraction_gt1000mm']):.3f}`",
            f"- Dry subtropics <250 mm fraction: `{row['dry_subtropics_fraction_lt250mm']:.3f}`",
            f"- Mountain >1000 mm fraction / Earth: `{_ratio(row['mountain_wet_fraction_gt1000mm'], earth['mountain_wet_fraction_gt1000mm']):.3f}`",
            f"- High-lat cold / tundra-ice fractions: `{row['high_lat_cold_fraction_lt0C']:.3f}` / `{row['high_lat_ice_tundra_fraction']:.3f}`",
            "",
        ])
    return "\n".join(lines)


def run_earth_climate_pattern_gate(
    config: EarthClimatePatternGateConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    earth_path = Path(config.earth_reference_npz)
    terminal_path = Path(config.terminal_summary_json)
    terminal = json.loads(terminal_path.read_text())

    earth_metrics = _pattern_metrics(
        label="Earth",
        preset="earth_reference",
        seed=-1,
        arrays_path=earth_path,
        earth=True,
    )

    generated: list[dict[str, Any]] = []
    checks: list[dict[str, Any]] = []
    for summary in terminal.get("summaries", []):
        label = Path(str(summary.get("assets_dir", ""))).name
        preset = str(summary.get("preset", ""))
        row = _pattern_metrics(
            label=label,
            preset=preset,
            seed=int(summary.get("seed", -1)),
            arrays_path=Path(summary["arrays"]),
            earth=False,
        )
        row["mode"] = (
            "earthlike_calibration"
            if "earthlike" in preset.lower()
            else "diagnostic_only"
        )
        generated.append(row)
        if row["mode"] == "earthlike_calibration":
            checks.extend(_earthlike_pattern_checks(row, earth_metrics))

    generated.sort(key=lambda row: (
        row.get("mode") != "earthlike_calibration",
        str(row.get("label")),
    ))
    checks.sort(key=lambda row: (
        row["passed"],
        row["severity"] != "fail",
        row["label"],
        row["group"],
        row["metric"],
    ))
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
        "tropical_land_share", "wet_tropics_precip_mean_mm_yr",
        "wet_tropics_precip_p75_mm_yr", "wet_tropics_precip_p90_mm_yr",
        "wet_tropics_fraction_gt1000mm", "wet_tropics_fraction_gt1500mm",
        "low_mid_monsoon_p75", "low_mid_monsoon_p90",
        "low_mid_monsoon_wet_fraction_ge1",
        "dry_subtropics_fraction_lt250mm",
        "dry_subtropics_fraction_lt500mm",
        "dry_subtropics_precip_p25_mm_yr", "mountain_land_share",
        "mountain_precip_p75_mm_yr", "mountain_wet_fraction_gt1000mm",
        "warm_mountain_precip_p75_mm_yr", "high_lat_land_share",
        "high_lat_temperature_mean_C", "high_lat_temperature_p75_C",
        "high_lat_cold_fraction_lt0C", "high_lat_ice_tundra_fraction",
        "forest_tropical_land_fraction", "tropical_biome_land_fraction",
    ]
    metrics_csv = _write_csv(
        outdir / "earth_climate_pattern_metrics.csv",
        [earth_metrics | {"mode": "earth_reference"}] + generated,
        metric_keys,
    )
    checks_csv = _write_csv(
        outdir / "earth_climate_pattern_checks.csv",
        checks,
        [
            "label", "group", "metric", "generated", "earth", "ratio_to_earth",
            "operator", "threshold", "severity", "passed", "skipped", "message",
        ],
    )

    report = {
        "schema": SCHEMA,
        "earth_reference_npz": str(earth_path),
        "terminal_summary_json": str(terminal_path),
        "verdict": verdict,
        "failure_count": len(failures),
        "warning_count": len(warnings),
        "skipped_count": len(skipped),
        "earth_metrics": earth_metrics,
        "generated_metrics": generated,
        "checks": checks,
        "failures": failures,
        "warnings": warnings,
        "skipped": skipped,
        "metrics_csv": str(metrics_csv),
        "checks_csv": str(checks_csv),
    }
    md_path = outdir / "earth_climate_pattern_gate_report.md"
    report["markdown"] = str(md_path)
    (outdir / "earth_climate_pattern_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default)
    )
    md_path.write_text(_render_markdown(report))
    return report
