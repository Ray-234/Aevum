"""F4 seasonal hydroclimate placement gate.

The scalar and subtype gates check broad Earth envelopes.  This gate checks the
process placement of seasonal precipitation in generated worlds: wet cells
should be supported by moisture, ITCZ, monsoon, storm-track, and wet-response
fields; dry cells should be explainable by low moisture or rain-shadow/dry
response; and the wet season should follow the season of maximum support.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from aevum.diagnostics.earth_climate_hydro_region_gate import (
    _array_path,
    _check,
    _json_default,
    _safe_float,
)


SCHEMA = "aevum.earth_climate_seasonal_hydro_placement_gate.v1"


@dataclass(frozen=True)
class EarthClimateSeasonalHydroPlacementGateConfig:
    earth_reference_npz: Path
    terminal_summary_json: Path
    outdir: Path


def _percentile(values: np.ndarray, q: float,
                mask: np.ndarray | None = None) -> float:
    values = np.asarray(values, dtype=np.float64)
    if mask is not None:
        values = values[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    return float(np.percentile(values, q)) if values.size else float("nan")


def _mean(values: np.ndarray, mask: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    return float(np.mean(values[mask])) if mask.any() else float("nan")


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


def _ratio(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or abs(den) <= 1.0e-12:
        return float("nan")
    return float(num / den)


def _earth_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        lat = np.asarray(z["lat"], dtype=np.float64)
        land = np.asarray(z["earth__land_mask"], dtype=bool)
        precip = np.asarray(z["earth__seasonal_precip_mm_yr_equiv"],
                            dtype=np.float64)
        annual = np.asarray(z["earth__annual_precip_mm"], dtype=np.float64)
    land4 = np.broadcast_to(land, precip.shape)
    wet = land4 & (precip >= _percentile(precip, 85, land4))
    dry = land4 & (precip <= _percentile(precip, 25, land4))
    return {
        "label": "earth_reference",
        "preset": "earth_reference",
        "seed": 0,
        "seasonal_precip_aggregate_max_delta": float(
            np.nanmax(np.abs(np.mean(precip, axis=0) - annual))),
        "wet_land_threshold_mm_yr": _percentile(precip, 85, land4),
        "dry_land_threshold_mm_yr": _percentile(precip, 25, land4),
        "tropical_wet_share": float(np.mean(
            wet[np.broadcast_to(np.abs(lat) <= 25.0, precip.shape) & land4]))
        if np.any(np.broadcast_to(np.abs(lat) <= 25.0, precip.shape) & land4) else float("nan"),
        "dry_cell_fraction": float(np.mean(dry[land4])) if land4.any() else float("nan"),
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
            "terrain__elevation_m",
            "sea_level_m",
            "climate__seasonal_precipitation",
            "climate__precipitation",
            "climate__monsoon_rainfall_corridor",
            "climate__storm_track_rainfall_corridor",
            "climate__rain_shadow_index",
            "climate__regional_precipitation_response",
            "climate__moisture_flow_precipitation_response",
            "atmosphere__moisture_access",
            "atmosphere__itcz_intensity",
        )
        if any(key not in z.files for key in required):
            return defaults
        lat = np.asarray(z["lat"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        land = np.asarray(z["terrain__elevation_m"], dtype=np.float64) >= sea
        precip = np.asarray(z["climate__seasonal_precipitation"], dtype=np.float64)
        annual = np.asarray(z["climate__precipitation"], dtype=np.float64)
        moisture = np.asarray(z["atmosphere__moisture_access"], dtype=np.float64)
        monsoon = np.asarray(z["climate__monsoon_rainfall_corridor"],
                             dtype=np.float64)
        storm = np.asarray(z["climate__storm_track_rainfall_corridor"],
                           dtype=np.float64)
        rain_shadow = np.asarray(z["climate__rain_shadow_index"],
                                 dtype=np.float64)
        response = np.asarray(z["climate__regional_precipitation_response"],
                              dtype=np.float64)
        flow_response = np.asarray(
            z["climate__moisture_flow_precipitation_response"], dtype=np.float64)
        itcz = np.asarray(z["atmosphere__itcz_intensity"], dtype=np.float64)

    land4 = np.broadcast_to(land, precip.shape)
    tropical = land4 & np.broadcast_to(np.abs(lat) <= 25.0, precip.shape)
    low_mid = land4 & np.broadcast_to(np.abs(lat) <= 35.0, precip.shape)
    midlat = land4 & np.broadcast_to(
        (np.abs(lat) >= 30.0) & (np.abs(lat) <= 60.0), precip.shape)
    wet_support = (
        np.clip(moisture, 0.0, 1.0)
        + 0.75 * np.clip(monsoon, 0.0, 1.2)
        + 0.55 * np.clip(storm, 0.0, 1.2)
        + 0.35 * np.clip(itcz, 0.0, 1.2)
        + 0.55 * np.clip(response - 1.0, 0.0, 1.0)
        + 0.45 * np.clip(flow_response - 1.0, 0.0, 1.0)
    )
    dry_support = (
        1.0 - np.clip(moisture, 0.0, 1.0)
        + 0.75 * np.clip(rain_shadow, 0.0, 1.2)
        + 0.55 * np.clip(1.0 - response, 0.0, 1.0)
        + 0.45 * np.clip(1.0 - flow_response, 0.0, 1.0)
    )
    wet = land4 & (precip >= _percentile(precip, 85, land4))
    dry = land4 & (precip <= _percentile(precip, 25, land4))
    wet_support_median = _percentile(wet_support, 50, land4)
    dry_support_threshold = _percentile(dry_support, 55, land4)
    low_moisture_threshold = _percentile(moisture, 35, land4)
    dry_explained = dry & (
        (dry_support >= dry_support_threshold)
        | (moisture <= low_moisture_threshold)
    )
    peak_precip = np.argmax(precip, axis=0)
    peak_support = np.argmax(wet_support, axis=0)
    active_land = land & (annual >= 250.0)
    return {
        **defaults,
        "required_fields_found": 1.0,
        "wet_support_p25_ratio": _ratio(
            _percentile(wet_support, 25, wet), wet_support_median),
        "unsupported_wet_fraction": float(np.mean(
            wet_support[wet] < _percentile(wet_support, 45, land4)))
        if wet.any() else float("nan"),
        "dry_explained_fraction": float(np.mean(dry_explained[dry]))
        if dry.any() else float("nan"),
        "support_precip_land_corr": _corr(wet_support, precip, land4),
        "peak_support_match_fraction": float(np.mean(
            peak_precip[active_land] == peak_support[active_land]))
        if active_land.any() else float("nan"),
        "monsoon_precip_ratio": _ratio(
            _mean(precip, low_mid & (monsoon >= _percentile(monsoon, 80, low_mid))),
            _mean(precip, low_mid & (monsoon <= _percentile(monsoon, 45, low_mid))),
        ),
        "storm_track_precip_ratio": _ratio(
            _mean(precip, midlat & (storm >= _percentile(storm, 80, midlat))),
            _mean(precip, midlat & (storm <= _percentile(storm, 45, midlat))),
        ),
        "itcz_precip_ratio": _ratio(
            _mean(precip, tropical & (itcz >= _percentile(itcz, 80, tropical))),
            _mean(precip, tropical & (itcz <= _percentile(itcz, 45, tropical))),
        ),
        "rain_shadow_precip_ratio": _ratio(
            _mean(precip, land4 & (rain_shadow >= _percentile(rain_shadow, 80, land4))),
            _mean(precip, land4 & (rain_shadow <= _percentile(rain_shadow, 45, land4))),
        ),
        "seasonal_precip_aggregate_max_delta": float(
            np.nanmax(np.abs(np.mean(precip, axis=0) - annual))),
    }


def _checks_for_row(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    preset = str(row.get("preset", "")).lower()
    checks = [
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row.get("arrays_found", 0.0), operator=">=",
               threshold=1.0, severity="fail",
               message="F4 placement gate requires archived arrays"),
        _check(label=label, group="array_archive", metric="required_fields_found",
               generated=row.get("required_fields_found", 0.0), operator=">=",
               threshold=1.0, severity="fail",
               message="F4 placement fields must include precipitation and support diagnostics"),
        _check(label=label, group="wet_support", metric="wet_support_p25_ratio",
               generated=row.get("wet_support_p25_ratio", float("nan")),
               operator=">=", threshold=1.10, severity="fail",
               message="top seasonal wet cells should have above-median process support"),
        _check(label=label, group="wet_support", metric="unsupported_wet_fraction",
               generated=row.get("unsupported_wet_fraction", float("nan")),
               operator="<=", threshold=0.15, severity="fail",
               message="top seasonal wet cells should not be unsupported patches"),
        _check(label=label, group="wet_support", metric="support_precip_land_corr",
               generated=row.get("support_precip_land_corr", float("nan")),
               operator=">=", threshold=0.45, severity="fail",
               message="seasonal precipitation should follow wet-process support"),
        _check(label=label, group="seasonal_phase", metric="peak_support_match_fraction",
               generated=row.get("peak_support_match_fraction", float("nan")),
               operator=">=", threshold=0.55, severity="fail",
               message="wet season should usually align with max support season"),
        _check(label=label, group="budget_closure",
               metric="seasonal_precip_aggregate_max_delta",
               generated=row.get("seasonal_precip_aggregate_max_delta", float("nan")),
               operator="<=", threshold=1.0e-6, severity="fail",
               message="seasonal precipitation must aggregate exactly to annual precipitation"),
    ]
    if "waterworld" in preset or "waterworld" in label.lower():
        checks.append(
            _check(label=label, group="dry_explanation",
                   metric="dry_explained_fraction",
                   generated=row.get("dry_explained_fraction", float("nan")),
                   operator=">=", threshold=0.40, severity="fail",
                   message="small-island dry cells should still be mostly explainable")
        )
    else:
        checks.append(
            _check(label=label, group="dry_explanation",
                   metric="dry_explained_fraction",
                   generated=row.get("dry_explained_fraction", float("nan")),
                   operator=">=", threshold=0.70, severity="fail",
                   message="dry seasonal land should be explained by low moisture or rain shadow")
        )
    if "earthlike" in preset or "earthlike" in label.lower():
        checks.extend([
            _check(label=label, group="process_ratios",
                   metric="monsoon_precip_ratio",
                   generated=row.get("monsoon_precip_ratio", float("nan")),
                   operator=">=", threshold=1.40, severity="fail",
                   message="monsoon corridors should be wetter than low-monsoon low/mid-lat land"),
            _check(label=label, group="process_ratios",
                   metric="storm_track_precip_ratio",
                   generated=row.get("storm_track_precip_ratio", float("nan")),
                   operator=">=", threshold=1.25, severity="fail",
                   message="storm-track corridors should be wetter than weak storm-track midlat land"),
            _check(label=label, group="process_ratios",
                   metric="itcz_precip_ratio",
                   generated=row.get("itcz_precip_ratio", float("nan")),
                   operator=">=", threshold=1.40, severity="fail",
                   message="tropical ITCZ-supported cells should be wetter than weak-ITCZ tropical cells"),
            _check(label=label, group="process_ratios",
                   metric="rain_shadow_precip_ratio",
                   generated=row.get("rain_shadow_precip_ratio", float("nan")),
                   operator="<=", threshold=0.85, severity="fail",
                   message="rain-shadow land should be drier than low-rain-shadow land"),
        ])
    return checks


FIELDNAMES = [
    "label",
    "preset",
    "seed",
    "arrays_found",
    "arrays_path",
    "required_fields_found",
    "wet_support_p25_ratio",
    "unsupported_wet_fraction",
    "dry_explained_fraction",
    "support_precip_land_corr",
    "peak_support_match_fraction",
    "monsoon_precip_ratio",
    "storm_track_precip_ratio",
    "itcz_precip_ratio",
    "rain_shadow_precip_ratio",
    "seasonal_precip_aggregate_max_delta",
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
        "# Earth Climate Seasonal Hydro Placement Gate",
        "",
        f"Schema: `{report['schema']}`",
        f"Verdict: **{report['verdict']}**",
        f"Failures: {report['failure_count']}",
        f"Warnings: {report['warning_count']}",
        f"Skipped checks: {report['skipped_count']}",
        "",
        "## Generated Runs",
        "",
    ]
    for row in report["generated"]:
        lines.append(
            f"- `{row['label']}`: wet support ratio "
            f"`{row.get('wet_support_p25_ratio')}`, dry explained "
            f"`{row.get('dry_explained_fraction')}`, support/precip corr "
            f"`{row.get('support_precip_land_corr')}`, peak match "
            f"`{row.get('peak_support_match_fraction')}`"
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


def run_earth_climate_seasonal_hydro_placement_gate(
    config: EarthClimateSeasonalHydroPlacementGateConfig,
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
    _write_csv(outdir / "earth_climate_seasonal_hydro_placement_metrics.csv",
               generated, FIELDNAMES)
    _write_csv(
        outdir / "earth_climate_seasonal_hydro_placement_checks.csv",
        checks,
        [
            "label", "group", "metric", "generated", "operator", "threshold",
            "severity", "passed", "skipped", "message",
        ],
    )
    (outdir / "earth_climate_seasonal_hydro_placement_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default))
    (outdir / "earth_climate_seasonal_hydro_placement_gate_report.md").write_text(
        _render_markdown(report))
    return report
