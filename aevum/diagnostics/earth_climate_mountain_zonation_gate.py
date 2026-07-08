"""Mountain ecological zonation gate against broad Earth envelopes."""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "aevum.earth_climate_mountain_zonation_gate.v1"


@dataclass(frozen=True)
class EarthClimateMountainZonationGateConfig:
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


def _weighted_mean(values: np.ndarray, area: np.ndarray, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values) & np.isfinite(area)
    if not mask.any():
        return float("nan")
    return float(np.average(values[mask], weights=area[mask]))


def _percentile(values: np.ndarray, q: float, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    if not mask.any():
        return float("nan")
    return float(np.nanpercentile(values[mask], q))


def _mountain_metrics(
    label: str,
    preset: str,
    seed: int,
    lat: np.ndarray,
    area: np.ndarray,
    land: np.ndarray,
    relative_elev: np.ndarray,
    temp_c: np.ndarray,
    precip: np.ndarray,
    biome: np.ndarray,
) -> dict[str, Any]:
    lat = np.asarray(lat, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)
    land = np.asarray(land, dtype=bool)
    relative_elev = np.asarray(relative_elev, dtype=np.float64)
    temp_c = np.asarray(temp_c, dtype=np.float64)
    precip = np.asarray(precip, dtype=np.float64)
    biome = np.asarray(biome, dtype=np.int16)

    lowland = land & (relative_elev >= 0.0) & (relative_elev < 500.0)
    mountain = land & (relative_elev >= 1000.0)
    high_mountain = land & (relative_elev >= 2000.0)
    midlat = np.abs(lat) < 60.0
    low_midlat = lowland & midlat
    mountain_midlat = mountain & midlat
    high_mountain_midlat = high_mountain & midlat
    alpine_ecology = (biome == 1) | (biome == 3) | (biome == 5)
    forest_tropical = (biome == 4) | (biome == 6)
    desert = biome == 2

    return {
        "label": label,
        "preset": preset,
        "seed": seed,
        "land_fraction": _area_fraction(area, land, np.ones_like(land, dtype=bool)),
        "mountain_land_fraction": _area_fraction(area, mountain, land),
        "high_mountain_land_fraction": _area_fraction(area, high_mountain, land),
        "lowland_land_fraction": _area_fraction(area, lowland, land),
        "mountain_temperature_mean_C": _weighted_mean(temp_c, area, mountain),
        "high_mountain_temperature_mean_C": _weighted_mean(temp_c, area, high_mountain),
        "lowland_temperature_mean_C": _weighted_mean(temp_c, area, lowland),
        "midlat_low_minus_mountain_temperature_C": (
            _weighted_mean(temp_c, area, low_midlat)
            - _weighted_mean(temp_c, area, mountain_midlat)
        ),
        "midlat_low_minus_high_mountain_temperature_C": (
            _weighted_mean(temp_c, area, low_midlat)
            - _weighted_mean(temp_c, area, high_mountain_midlat)
        ),
        "mountain_precip_p50_mm_yr": _percentile(precip, 50, mountain),
        "mountain_precip_p75_mm_yr": _percentile(precip, 75, mountain),
        "mountain_wet_fraction_gt1000mm": _area_fraction(
            area, mountain & (precip >= 1000.0), mountain),
        "lowland_precip_p50_mm_yr": _percentile(precip, 50, lowland),
        "mountain_to_lowland_precip_p50_ratio": (
            _percentile(precip, 50, mountain)
            / max(_percentile(precip, 50, lowland), 1.0e-9)
        ),
        "mountain_alpine_ecology_fraction": _area_fraction(
            area, mountain & alpine_ecology, mountain),
        "mountain_forest_tropical_fraction": _area_fraction(
            area, mountain & forest_tropical, mountain),
        "mountain_desert_fraction": _area_fraction(area, mountain & desert, mountain),
        "high_mountain_alpine_ecology_fraction": _area_fraction(
            area, high_mountain & alpine_ecology, high_mountain),
        "high_mountain_desert_fraction": _area_fraction(
            area, high_mountain & desert, high_mountain),
        "midlat_high_mountain_alpine_ecology_fraction": _area_fraction(
            area, high_mountain_midlat & alpine_ecology, high_mountain_midlat),
        "midlat_high_mountain_desert_fraction": _area_fraction(
            area, high_mountain_midlat & desert, high_mountain_midlat),
    }


def _earth_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        return _mountain_metrics(
            "Earth mountain reference",
            "earth_reference",
            -1,
            np.asarray(z["lat"], dtype=np.float64),
            np.asarray(z["cell_area"], dtype=np.float64),
            np.asarray(z["earth__land_mask"], dtype=bool),
            np.asarray(z["earth__elevation_m"], dtype=np.float64),
            np.asarray(z["earth__annual_temperature_C"], dtype=np.float64),
            np.asarray(z["earth__annual_precip_mm"], dtype=np.float64),
            np.asarray(z["earth__biome_class_proxy"], dtype=np.int16),
        )


def _generated_metrics(summary: dict[str, Any]) -> dict[str, Any]:
    arrays = Path(summary["arrays"])
    with np.load(arrays, allow_pickle=False) as z:
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        elev = np.asarray(z["terrain__elevation_m"], dtype=np.float64)
        land = elev >= sea
        row = _mountain_metrics(
            Path(str(summary.get("assets_dir", ""))).name,
            str(summary.get("preset", "")),
            int(summary.get("seed", -1)),
            np.asarray(z["lat"], dtype=np.float64),
            np.asarray(z["cell_area"], dtype=np.float64),
            land,
            elev - sea,
            np.asarray(z["climate__surface_temperature"], dtype=np.float64) - 273.15,
            np.asarray(z["climate__precipitation"], dtype=np.float64),
            np.asarray(z["biosphere__biome"], dtype=np.float64).astype(np.int16),
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
            "midlat_low_minus_mountain_temperature_C",
            row["midlat_low_minus_mountain_temperature_C"],
            earth["midlat_low_minus_mountain_temperature_C"],
            "ratio>=",
            0.45,
            "fail",
            "mountains should be meaningfully cooler than comparable lowlands",
        ),
        _check(
            label,
            "mountain_wet_fraction_gt1000mm",
            row["mountain_wet_fraction_gt1000mm"],
            earth["mountain_wet_fraction_gt1000mm"],
            "ratio>=",
            0.25,
            "warn",
            "mountain wet-tail is weak relative to Earth",
        ),
        _check(
            label,
            "high_mountain_alpine_ecology_fraction",
            row["high_mountain_alpine_ecology_fraction"],
            earth["high_mountain_alpine_ecology_fraction"],
            "ratio>=",
            0.55,
            "fail",
            "high mountains should express alpine grass/tundra/ice ecology",
        ),
        _check(
            label,
            "high_mountain_desert_fraction",
            row["high_mountain_desert_fraction"],
            earth["high_mountain_desert_fraction"],
            "<=",
            0.25,
            "fail",
            "high mountains are overclassified as desert",
        ),
        _check(
            label,
            "midlat_high_mountain_desert_fraction",
            row["midlat_high_mountain_desert_fraction"],
            earth["midlat_high_mountain_desert_fraction"],
            "<=",
            0.30,
            "warn",
            "mid-latitude high mountains should not mostly be desert",
        ),
        _check(
            label,
            "mountain_desert_fraction",
            row["mountain_desert_fraction"],
            earth["mountain_desert_fraction"],
            "<=",
            0.40,
            "warn",
            "mountain desert fraction is high",
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
        "# Earth Climate Mountain Zonation Gate",
        "",
        f"Verdict: `{report['verdict']}`",
        f"Failures: `{report['failure_count']}`",
        f"Warnings: `{report['warning_count']}`",
        "",
        "## Earth Mountain Envelope",
        "",
        f"- Mountain land fraction: `{earth['mountain_land_fraction']:.3f}`",
        f"- High-mountain alpine ecology fraction: `{earth['high_mountain_alpine_ecology_fraction']:.3f}`",
        f"- High-mountain desert fraction: `{earth['high_mountain_desert_fraction']:.3f}`",
        f"- Midlat low-minus-mountain temperature: `{earth['midlat_low_minus_mountain_temperature_C']:.3f} C`",
        f"- Mountain wet fraction >1000 mm/yr: `{earth['mountain_wet_fraction_gt1000mm']:.3f}`",
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
    lines.extend(["", "## Generated Earthlike Mountain Metrics", ""])
    for row in report["generated_metrics"]:
        if row.get("mode") != "earthlike_calibration":
            continue
        lines.extend([
            f"### {row['label']}",
            "",
            f"- Mountain/high-mountain land: `{row['mountain_land_fraction']:.3f}` / `{row['high_mountain_land_fraction']:.3f}`",
            f"- Low-minus-mountain temp: `{row['midlat_low_minus_mountain_temperature_C']:.3f} C`",
            f"- High-mountain alpine/desert: `{row['high_mountain_alpine_ecology_fraction']:.3f}` / `{row['high_mountain_desert_fraction']:.3f}`",
            f"- Mountain wet fraction >1000 mm/yr: `{row['mountain_wet_fraction_gt1000mm']:.3f}`",
            "",
        ])
    return "\n".join(lines)


def run_earth_climate_mountain_zonation_gate(
    config: EarthClimateMountainZonationGateConfig,
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
        "mountain_land_fraction", "high_mountain_land_fraction",
        "lowland_land_fraction", "mountain_temperature_mean_C",
        "high_mountain_temperature_mean_C", "lowland_temperature_mean_C",
        "midlat_low_minus_mountain_temperature_C",
        "midlat_low_minus_high_mountain_temperature_C",
        "mountain_precip_p50_mm_yr", "mountain_precip_p75_mm_yr",
        "mountain_wet_fraction_gt1000mm", "lowland_precip_p50_mm_yr",
        "mountain_to_lowland_precip_p50_ratio",
        "mountain_alpine_ecology_fraction",
        "mountain_forest_tropical_fraction", "mountain_desert_fraction",
        "high_mountain_alpine_ecology_fraction",
        "high_mountain_desert_fraction",
        "midlat_high_mountain_alpine_ecology_fraction",
        "midlat_high_mountain_desert_fraction",
    ]
    metrics_csv = _write_csv(
        outdir / "earth_climate_mountain_zonation_metrics.csv",
        [earth | {"mode": "earth_reference"}] + generated,
        metric_keys,
    )
    checks_csv = _write_csv(
        outdir / "earth_climate_mountain_zonation_checks.csv",
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
    md_path = outdir / "earth_climate_mountain_zonation_gate_report.md"
    report["markdown"] = str(md_path)
    (outdir / "earth_climate_mountain_zonation_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default)
    )
    md_path.write_text(_render_markdown(report))
    return report
