"""C5e coupling-convergence gate.

This gate checks that reduced climate feedbacks remain bounded after
hydroclimate starts feeding back into pressure and wind.  It is intentionally
not an Earth map-overlap test; it checks internal convergence and regression
guardrails for the current reduced coupling loop.
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
)


SCHEMA = "aevum.earth_climate_coupling_convergence_gate.v1"


@dataclass(frozen=True)
class EarthClimateCouplingConvergenceGateConfig:
    terminal_summary_json: Path
    outdir: Path


def _percentile(values: np.ndarray, q: float, mask: np.ndarray | None = None) -> float:
    values = np.asarray(values, dtype=np.float64)
    if mask is not None:
        values = values[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    return float(np.percentile(values, q)) if values.size else float("nan")


def _ratio(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or abs(den) <= 1.0e-12:
        return float("nan")
    return float(num / den)


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
            "cell_area",
            "terrain__elevation_m",
            "sea_level_m",
            "climate__precipitation",
            "climate__seasonal_precipitation",
            "climate__coupling_residual",
            "climate__ocean_evaporation_feedback",
            "climate__hydro_coupling_residual",
            "climate__hydro_feedback_iteration_delta",
            "atmosphere__land_sea_pressure_proxy",
            "atmosphere__seasonal_pressure_proxy",
            "atmosphere__precipitation_pressure_feedback",
            "atmosphere__seasonal_wind",
            "atmosphere__hydro_coupled_wind_anomaly",
            "ocean__wind_stress_current_response",
        )
        if any(key not in z.files for key in required):
            return defaults
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        area = np.asarray(z["cell_area"], dtype=np.float64)
        land = np.asarray(z["terrain__elevation_m"], dtype=np.float64) >= sea
        precip = np.asarray(z["climate__precipitation"], dtype=np.float64)
        seasonal_precip = np.asarray(
            z["climate__seasonal_precipitation"], dtype=np.float64)
        ocean_residual = np.asarray(z["climate__coupling_residual"], dtype=np.float64)
        evaporation_feedback = np.asarray(
            z["climate__ocean_evaporation_feedback"], dtype=np.float64)
        hydro_residual = np.asarray(
            z["climate__hydro_coupling_residual"], dtype=np.float64)
        iteration_delta = np.asarray(
            z["climate__hydro_feedback_iteration_delta"], dtype=np.float64)
        base_pressure = np.asarray(
            z["atmosphere__land_sea_pressure_proxy"], dtype=np.float64)
        final_pressure = np.asarray(
            z["atmosphere__seasonal_pressure_proxy"], dtype=np.float64)
        feedback = np.asarray(
            z["atmosphere__precipitation_pressure_feedback"], dtype=np.float64)
        wind = np.asarray(z["atmosphere__seasonal_wind"], dtype=np.float64)
        wind_anomaly = np.asarray(
            z["atmosphere__hydro_coupled_wind_anomaly"], dtype=np.float64)
        wind_stress_response = np.asarray(
            z["ocean__wind_stress_current_response"], dtype=np.float64)

    land4 = np.broadcast_to(land, feedback.shape)
    ocean = ~land
    ocean_evap_mean = (
        float(np.average(evaporation_feedback[ocean], weights=area[ocean]))
        if ocean.any() else 0.0
    )
    wind_speed = np.linalg.norm(wind, axis=2)
    wind_anomaly_speed = np.linalg.norm(wind_anomaly, axis=2)
    annual_wind = wind.mean(axis=0)
    annual_wind_speed = np.linalg.norm(annual_wind, axis=1)
    wind_response_speed = np.linalg.norm(wind_stress_response, axis=1)
    wind_response_active = ocean & (wind_response_speed > 0.01)
    wind_alignment = np.full_like(wind_response_speed, np.nan, dtype=np.float64)
    align_denom = np.maximum(wind_response_speed * annual_wind_speed, 1.0e-12)
    wind_alignment[wind_response_active] = (
        np.sum(wind_stress_response * annual_wind, axis=1)[wind_response_active]
        / align_denom[wind_response_active]
    )
    feedback_abs_p95 = _percentile(np.abs(feedback), 95, land4)
    pressure_abs_p95 = _percentile(np.abs(base_pressure), 95, land4)
    wind_anom_p95 = _percentile(wind_anomaly_speed, 95, land4)
    wind_p95 = _percentile(wind_speed, 95, land4)
    return {
        **defaults,
        "required_fields_found": 1.0,
        "land_fraction": float(np.mean(land)),
        "pressure_feedback_abs_p95": feedback_abs_p95,
        "pressure_feedback_abs_max": float(np.nanmax(np.abs(feedback))),
        "pressure_feedback_to_base_pressure_p95_ratio": _ratio(
            feedback_abs_p95, pressure_abs_p95),
        "final_pressure_abs_p99": _percentile(np.abs(final_pressure), 99, land4),
        "wind_anomaly_p95_m_s": wind_anom_p95,
        "wind_anomaly_p99_m_s": _percentile(wind_anomaly_speed, 99, land4),
        "wind_anomaly_to_wind_p95_ratio": _ratio(wind_anom_p95, wind_p95),
        "hydro_coupling_residual_p95": _percentile(hydro_residual, 95, land),
        "hydro_coupling_residual_max": float(np.nanmax(hydro_residual)),
        "hydro_feedback_iteration_delta_p95": _percentile(
            iteration_delta, 95, land),
        "hydro_feedback_iteration_delta_max": float(np.nanmax(iteration_delta)),
        "ocean_coupling_residual_p95": _percentile(ocean_residual, 95),
        "ocean_evaporation_feedback_abs_p95_C": _percentile(
            np.abs(evaporation_feedback), 95, ocean),
        "ocean_evaporation_feedback_abs_max_C": float(
            np.nanmax(np.abs(evaporation_feedback))),
        "ocean_evaporation_feedback_weighted_mean_abs_C": abs(ocean_evap_mean),
        "wind_stress_response_ocean_p50_m_s": _percentile(
            wind_response_speed, 50, ocean),
        "wind_stress_response_ocean_p95_m_s": _percentile(
            wind_response_speed, 95, ocean),
        "wind_stress_response_land_max_m_s": (
            float(np.nanmax(wind_response_speed[land])) if land.any() else 0.0),
        "wind_stress_response_alignment_p50": _percentile(
            wind_alignment, 50, wind_response_active),
        "wind_stress_response_to_wind_p95_ratio": _ratio(
            _percentile(wind_response_speed, 95, ocean),
            _percentile(annual_wind_speed, 95, ocean),
        ),
        "seasonal_precip_aggregate_max_delta": float(
            np.nanmax(np.abs(np.mean(seasonal_precip, axis=0) - precip))),
    }


def _checks_for_row(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    preset = str(row.get("preset", "")).lower()
    label_l = label.lower()
    checks = [
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row.get("arrays_found", 0.0), operator=">=",
               threshold=1.0, severity="fail",
               message="C5e coupling gate requires archived arrays"),
        _check(label=label, group="array_archive", metric="required_fields_found",
               generated=row.get("required_fields_found", 0.0), operator=">=",
               threshold=1.0, severity="fail",
               message="C5e arrays must include coupling feedback fields"),
        _check(label=label, group="pressure_feedback",
               metric="pressure_feedback_abs_p95",
               generated=row.get("pressure_feedback_abs_p95", float("nan")),
               operator="<=", threshold=0.060, severity="fail",
               message="precipitation-pressure feedback must remain bounded"),
        _check(label=label, group="pressure_feedback",
               metric="pressure_feedback_abs_max",
               generated=row.get("pressure_feedback_abs_max", float("nan")),
               operator="<=", threshold=0.076, severity="fail",
               message="precipitation-pressure feedback hard cap must hold"),
        _check(label=label, group="pressure_feedback",
               metric="pressure_feedback_to_base_pressure_p95_ratio",
               generated=row.get(
                   "pressure_feedback_to_base_pressure_p95_ratio", float("nan")),
               operator="<=", threshold=0.35, severity="fail",
               message="hydro feedback should not dominate base seasonal pressure"),
        _check(label=label, group="wind_feedback",
               metric="wind_anomaly_p95_m_s",
               generated=row.get("wind_anomaly_p95_m_s", float("nan")),
               operator="<=", threshold=0.18, severity="fail",
               message="hydro wind anomaly should stay a small near-surface correction"),
        _check(label=label, group="wind_feedback",
               metric="wind_anomaly_p99_m_s",
               generated=row.get("wind_anomaly_p99_m_s", float("nan")),
               operator="<=", threshold=0.42, severity="fail",
               message="hydro wind anomaly tail should remain capped"),
        _check(label=label, group="wind_feedback",
               metric="wind_anomaly_to_wind_p95_ratio",
               generated=row.get("wind_anomaly_to_wind_p95_ratio", float("nan")),
               operator="<=", threshold=0.08, severity="fail",
               message="hydro wind anomaly should not dominate seasonal winds"),
        _check(label=label, group="convergence",
               metric="hydro_coupling_residual_p95",
               generated=row.get("hydro_coupling_residual_p95", float("nan")),
               operator="<=", threshold=0.055, severity="fail",
               message="hydro coupling residual should remain small"),
        _check(label=label, group="convergence",
               metric="hydro_coupling_residual_max",
               generated=row.get("hydro_coupling_residual_max", float("nan")),
               operator="<=", threshold=0.20, severity="fail",
               message="hydro coupling residual hard cap must hold"),
        _check(label=label, group="convergence",
               metric="hydro_feedback_iteration_delta_p95",
               generated=row.get(
                   "hydro_feedback_iteration_delta_p95", float("nan")),
               operator="<=", threshold=0.030, severity="fail",
               message="hydro feedback should converge across bounded iterations"),
        _check(label=label, group="convergence",
               metric="hydro_feedback_iteration_delta_max",
               generated=row.get(
                   "hydro_feedback_iteration_delta_max", float("nan")),
               operator="<=", threshold=0.12, severity="fail",
               message="hydro feedback iteration delta hard cap must hold"),
        _check(label=label, group="convergence",
               metric="ocean_coupling_residual_p95",
               generated=row.get("ocean_coupling_residual_p95", float("nan")),
               operator="<=", threshold=0.08, severity="fail",
               message="ocean/SST/wind fixed-point residual should remain bounded"),
        _check(label=label, group="evaporation_sst_feedback",
               metric="ocean_evaporation_feedback_abs_p95_C",
               generated=row.get(
                   "ocean_evaporation_feedback_abs_p95_C", float("nan")),
               operator="<=", threshold=0.95, severity="fail",
               message="evaporation-SST feedback should stay a small heat-flux correction"),
        _check(label=label, group="evaporation_sst_feedback",
               metric="ocean_evaporation_feedback_abs_max_C",
               generated=row.get(
                   "ocean_evaporation_feedback_abs_max_C", float("nan")),
               operator="<=", threshold=1.05, severity="fail",
               message="evaporation-SST feedback hard cap must hold"),
        _check(label=label, group="evaporation_sst_feedback",
               metric="ocean_evaporation_feedback_weighted_mean_abs_C",
               generated=row.get(
                   "ocean_evaporation_feedback_weighted_mean_abs_C", float("nan")),
               operator="<=", threshold=1.0e-6, severity="fail",
               message="evaporation-SST feedback should redistribute heat without net forcing"),
        _check(label=label, group="wind_stress_current_response",
               metric="wind_stress_response_ocean_p50_m_s",
               generated=row.get("wind_stress_response_ocean_p50_m_s", float("nan")),
               operator=">=", threshold=0.030, severity="fail",
               message="ocean currents should include a visible wind-stress response"),
        _check(label=label, group="wind_stress_current_response",
               metric="wind_stress_response_ocean_p95_m_s",
               generated=row.get("wind_stress_response_ocean_p95_m_s", float("nan")),
               operator="<=", threshold=0.28, severity="fail",
               message="wind-stress current response should remain bounded"),
        _check(label=label, group="wind_stress_current_response",
               metric="wind_stress_response_land_max_m_s",
               generated=row.get("wind_stress_response_land_max_m_s", float("nan")),
               operator="<=", threshold=1.0e-8, severity="fail",
               message="wind-stress current response must stay off land"),
        _check(label=label, group="wind_stress_current_response",
               metric="wind_stress_response_alignment_p50",
               generated=row.get("wind_stress_response_alignment_p50", float("nan")),
               operator=">=", threshold=0.55, severity="fail",
               message="wind-stress current response should align with seasonal wind stress"),
        _check(label=label, group="wind_stress_current_response",
               metric="wind_stress_response_to_wind_p95_ratio",
               generated=row.get(
                   "wind_stress_response_to_wind_p95_ratio", float("nan")),
               operator="<=", threshold=0.06, severity="fail",
               message="wind-stress response should stay a small fraction of wind speed"),
        _check(label=label, group="budget_closure",
               metric="seasonal_precip_aggregate_max_delta",
               generated=row.get("seasonal_precip_aggregate_max_delta", float("nan")),
               operator="<=", threshold=1.0e-6, severity="fail",
               message="seasonal precipitation must aggregate exactly to annual"),
    ]
    if "waterworld" in preset or "waterworld" in label_l:
        checks.extend([
            _check(label=label, group="waterworld_false_positive",
                   metric="pressure_feedback_abs_p95",
                   generated=row.get("pressure_feedback_abs_p95", float("nan")),
                   operator="<=", threshold=0.006, severity="fail",
                   message="waterworld islands should not create broad hydro pressure feedback"),
            _check(label=label, group="waterworld_false_positive",
                   metric="wind_anomaly_p95_m_s",
                   generated=row.get("wind_anomaly_p95_m_s", float("nan")),
                   operator="<=", threshold=0.015, severity="fail",
                   message="waterworld islands should not create broad hydro wind anomalies"),
            _check(label=label, group="waterworld_false_positive",
                   metric="hydro_feedback_iteration_delta_p95",
                   generated=row.get(
                       "hydro_feedback_iteration_delta_p95", float("nan")),
                   operator="<=", threshold=0.003, severity="fail",
                   message="waterworld hydro feedback should not oscillate"),
        ])
    return checks


FIELDNAMES = [
    "label",
    "preset",
    "seed",
    "arrays_found",
    "arrays_path",
    "required_fields_found",
    "land_fraction",
    "pressure_feedback_abs_p95",
    "pressure_feedback_abs_max",
    "pressure_feedback_to_base_pressure_p95_ratio",
    "final_pressure_abs_p99",
    "wind_anomaly_p95_m_s",
    "wind_anomaly_p99_m_s",
    "wind_anomaly_to_wind_p95_ratio",
    "hydro_coupling_residual_p95",
    "hydro_coupling_residual_max",
    "hydro_feedback_iteration_delta_p95",
    "hydro_feedback_iteration_delta_max",
    "ocean_coupling_residual_p95",
    "ocean_evaporation_feedback_abs_p95_C",
    "ocean_evaporation_feedback_abs_max_C",
    "ocean_evaporation_feedback_weighted_mean_abs_C",
    "wind_stress_response_ocean_p50_m_s",
    "wind_stress_response_ocean_p95_m_s",
    "wind_stress_response_land_max_m_s",
    "wind_stress_response_alignment_p50",
    "wind_stress_response_to_wind_p95_ratio",
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
        "# Earth Climate Coupling Convergence Gate",
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
            f"- `{row['label']}`: pressure feedback p95 "
            f"`{row.get('pressure_feedback_abs_p95')}`, wind anomaly p95 "
            f"`{row.get('wind_anomaly_p95_m_s')}`, hydro residual p95 "
            f"`{row.get('hydro_coupling_residual_p95')}`, iteration delta p95 "
            f"`{row.get('hydro_feedback_iteration_delta_p95')}`, ocean evap "
            f"feedback p95 `{row.get('ocean_evaporation_feedback_abs_p95_C')}`, "
            f"wind-stress response p95 "
            f"`{row.get('wind_stress_response_ocean_p95_m_s')}`"
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


def run_earth_climate_coupling_convergence_gate(
    config: EarthClimateCouplingConvergenceGateConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
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
        "terminal_summary_json": str(config.terminal_summary_json),
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
    _write_csv(outdir / "earth_climate_coupling_convergence_metrics.csv",
               generated, FIELDNAMES)
    _write_csv(
        outdir / "earth_climate_coupling_convergence_checks.csv",
        checks,
        [
            "label", "group", "metric", "generated", "operator", "threshold",
            "severity", "passed", "skipped", "message",
        ],
    )
    (outdir / "earth_climate_coupling_convergence_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default))
    (outdir / "earth_climate_coupling_convergence_gate_report.md").write_text(
        _render_markdown(report))
    return report
