"""Spatial biome-organization gate against broad Earth latitude envelopes."""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "aevum.earth_climate_spatial_biome_gate.v1"


@dataclass(frozen=True)
class EarthClimateSpatialBiomeGateConfig:
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


def _band(lat: np.ndarray, lo: float, hi: float) -> np.ndarray:
    alat = np.abs(np.asarray(lat, dtype=np.float64))
    return (alat >= lo) & (alat < hi)


def _biome_spatial_metrics(
    label: str,
    preset: str,
    seed: int,
    lat: np.ndarray,
    area: np.ndarray,
    land: np.ndarray,
    biome: np.ndarray,
) -> dict[str, Any]:
    lat = np.asarray(lat, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)
    land = np.asarray(land, dtype=bool)
    biome = np.asarray(biome, dtype=np.int16)

    low_tropics = land & _band(lat, 0.0, 15.0)
    subtropics = land & _band(lat, 15.0, 30.0)
    warm_midlat = land & _band(lat, 30.0, 45.0)
    cool_midlat = land & _band(lat, 45.0, 60.0)
    high_lat = land & _band(lat, 60.0, 90.0)
    all_land = land

    forest = biome == 4
    tropical = biome == 6
    desert = biome == 2
    grassland = biome == 3
    tundra_ice = (biome == 1) | (biome == 5)

    tropical_land = land & tropical
    forest_land = land & forest
    desert_land = land & desert

    return {
        "label": label,
        "preset": preset,
        "seed": seed,
        "land_fraction": _area_fraction(area, land, np.ones_like(land, dtype=bool)),
        "low_tropics_land_fraction": _area_fraction(area, low_tropics, all_land),
        "low_tropics_forest_tropical_fraction": _area_fraction(
            area, low_tropics & (forest | tropical), low_tropics),
        "low_tropics_tropical_fraction": _area_fraction(
            area, low_tropics & tropical, low_tropics),
        "subtropical_desert_fraction": _area_fraction(
            area, subtropics & desert, subtropics),
        "warm_midlat_desert_fraction": _area_fraction(
            area, warm_midlat & desert, warm_midlat),
        "warm_midlat_forest_fraction": _area_fraction(
            area, warm_midlat & forest, warm_midlat),
        "cool_midlat_desert_fraction": _area_fraction(
            area, cool_midlat & desert, cool_midlat),
        "cool_midlat_forest_fraction": _area_fraction(
            area, cool_midlat & forest, cool_midlat),
        "cool_midlat_grassland_fraction": _area_fraction(
            area, cool_midlat & grassland, cool_midlat),
        "high_lat_desert_fraction": _area_fraction(area, high_lat & desert, high_lat),
        "high_lat_tundra_ice_fraction": _area_fraction(
            area, high_lat & tundra_ice, high_lat),
        "tropical_low_lat_share": _area_fraction(
            area, tropical_land & _band(lat, 0.0, 30.0), tropical_land),
        "desert_subtropical_midlat_share": _area_fraction(
            area, desert_land & _band(lat, 15.0, 45.0), desert_land),
        "forest_mid_high_share": _area_fraction(
            area, forest_land & _band(lat, 30.0, 75.0), forest_land),
    }


def _earth_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        return _biome_spatial_metrics(
            "Earth Koppen proxy",
            "earth_reference",
            -1,
            np.asarray(z["lat"], dtype=np.float64),
            np.asarray(z["cell_area"], dtype=np.float64),
            np.asarray(z["earth__land_mask"], dtype=bool),
            np.asarray(z["earth__biome_class_proxy"], dtype=np.int16),
        )


def _generated_metrics(summary: dict[str, Any]) -> dict[str, Any]:
    arrays = Path(summary["arrays"])
    with np.load(arrays, allow_pickle=False) as z:
        lat = np.asarray(z["lat"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        land = np.asarray(z["terrain__elevation_m"], dtype=np.float64) >= sea
        biome = np.asarray(z["biosphere__biome"], dtype=np.float64).astype(np.int16)
    row = _biome_spatial_metrics(
        Path(str(summary.get("assets_dir", ""))).name,
        str(summary.get("preset", "")),
        int(summary.get("seed", -1)),
        lat,
        area,
        land,
        biome,
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
            "low_tropics_forest_tropical_fraction",
            row["low_tropics_forest_tropical_fraction"],
            earth["low_tropics_forest_tropical_fraction"],
            "ratio>=",
            0.65,
            "fail",
            "low-latitude land should retain a broad forest/tropical envelope",
        ),
        _check(
            label,
            "tropical_low_lat_share",
            row["tropical_low_lat_share"],
            earth["tropical_low_lat_share"],
            ">=",
            0.85,
            "fail",
            "tropical biome should be concentrated in low latitudes",
        ),
        _check(
            label,
            "subtropical_desert_fraction",
            row["subtropical_desert_fraction"],
            earth["subtropical_desert_fraction"],
            "ratio>=",
            0.40,
            "fail",
            "earthlike worlds should express a subtropical dry-belt envelope",
        ),
        _check(
            label,
            "cool_midlat_forest_fraction",
            row["cool_midlat_forest_fraction"],
            earth["cool_midlat_forest_fraction"],
            "ratio>=",
            0.20,
            "fail",
            "cool mid-latitudes should not lose the temperate/boreal forest belt",
        ),
        _check(
            label,
            "cool_midlat_desert_fraction",
            row["cool_midlat_desert_fraction"],
            earth["cool_midlat_desert_fraction"],
            "<=",
            0.32,
            "fail",
            "cool mid-latitude land is too often classified as desert",
        ),
        _check(
            label,
            "high_lat_tundra_ice_fraction",
            row["high_lat_tundra_ice_fraction"],
            earth["high_lat_tundra_ice_fraction"],
            ">=",
            0.45,
            "fail",
            "high-latitude land should mostly become tundra/ice rather than dry desert",
        ),
        _check(
            label,
            "high_lat_desert_fraction",
            row["high_lat_desert_fraction"],
            earth["high_lat_desert_fraction"],
            "<=",
            0.20,
            "fail",
            "high-latitude desert extent should stay limited in earthlike calibration",
        ),
        _check(
            label,
            "desert_subtropical_midlat_share",
            row["desert_subtropical_midlat_share"],
            earth["desert_subtropical_midlat_share"],
            ">=",
            0.45,
            "warn",
            "most desert area should be tied to subtropical or warm mid-latitude belts",
        ),
        _check(
            label,
            "forest_mid_high_share",
            row["forest_mid_high_share"],
            earth["forest_mid_high_share"],
            "ratio>=",
            0.45,
            "warn",
            "forest area is too concentrated away from temperate/boreal belts",
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
        "# Earth Climate Spatial Biome Gate",
        "",
        f"Verdict: `{report['verdict']}`",
        f"Failures: `{report['failure_count']}`",
        f"Warnings: `{report['warning_count']}`",
        "",
        "## Earth Reference Envelope",
        "",
        f"- Low-tropics forest+tropical fraction: `{earth['low_tropics_forest_tropical_fraction']:.3f}`",
        f"- Subtropical desert fraction: `{earth['subtropical_desert_fraction']:.3f}`",
        f"- Cool-midlat forest fraction: `{earth['cool_midlat_forest_fraction']:.3f}`",
        f"- Cool-midlat desert fraction: `{earth['cool_midlat_desert_fraction']:.3f}`",
        f"- High-lat tundra/ice fraction: `{earth['high_lat_tundra_ice_fraction']:.3f}`",
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
    lines.extend(["", "## Generated Earthlike Spatial Metrics", ""])
    for row in report["generated_metrics"]:
        if row.get("mode") != "earthlike_calibration":
            continue
        lines.extend([
            f"### {row['label']}",
            "",
            f"- Low-tropics forest+tropical: `{row['low_tropics_forest_tropical_fraction']:.3f}`",
            f"- Subtropical desert: `{row['subtropical_desert_fraction']:.3f}`",
            f"- Cool-midlat forest / desert: `{row['cool_midlat_forest_fraction']:.3f}` / `{row['cool_midlat_desert_fraction']:.3f}`",
            f"- High-lat tundra+ice / desert: `{row['high_lat_tundra_ice_fraction']:.3f}` / `{row['high_lat_desert_fraction']:.3f}`",
            f"- Tropical low-lat share: `{row['tropical_low_lat_share']:.3f}`",
            "",
        ])
    return "\n".join(lines)


def run_earth_climate_spatial_biome_gate(
    config: EarthClimateSpatialBiomeGateConfig,
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
        "low_tropics_land_fraction", "low_tropics_forest_tropical_fraction",
        "low_tropics_tropical_fraction", "subtropical_desert_fraction",
        "warm_midlat_desert_fraction", "warm_midlat_forest_fraction",
        "cool_midlat_desert_fraction", "cool_midlat_forest_fraction",
        "cool_midlat_grassland_fraction", "high_lat_desert_fraction",
        "high_lat_tundra_ice_fraction", "tropical_low_lat_share",
        "desert_subtropical_midlat_share", "forest_mid_high_share",
    ]
    metrics_csv = _write_csv(
        outdir / "earth_climate_spatial_biome_metrics.csv",
        [earth | {"mode": "earth_reference_koppen_proxy"}] + generated,
        metric_keys,
    )
    checks_csv = _write_csv(
        outdir / "earth_climate_spatial_biome_checks.csv",
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
    md_path = outdir / "earth_climate_spatial_biome_gate_report.md"
    report["markdown"] = str(md_path)
    (outdir / "earth_climate_spatial_biome_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default)
    )
    md_path.write_text(_render_markdown(report))
    return report
