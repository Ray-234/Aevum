"""Earth-based fitting report for the generated climate pipeline.

The report is intentionally diagnostic.  It reads an existing
``earth-climate-compare`` summary and turns the Earth-vs-generated deltas into
phase-level priorities for the climate fitting plan.  It does not rerun
tectonics, terrain, climate, or biomes.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "aevum.earth_climate_fitting_report.v2"

PHASES = (
    "F1_temperature_energy",
    "F2_sst_ocean_currents",
    "F3_circulation_moisture_access",
    "F4_seasonal_hydroclimate",
    "F5_koppen_biomes",
)


@dataclass(frozen=True)
class EarthClimateFittingConfig:
    comparison_summary_json: Path
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


def _ratio(generated: Any, earth: Any) -> float:
    g = _safe_float(generated)
    e = _safe_float(earth)
    if not np.isfinite(g) or not np.isfinite(e) or abs(e) <= 1.0e-12:
        return float("nan")
    return float(g / e)


def _mean(values: list[float]) -> float:
    finite = [v for v in values if np.isfinite(v)]
    if not finite:
        return float("nan")
    return float(np.mean(finite))


def _comparison_value(entry: dict[str, Any], key: str, field: str = "generated") -> float:
    row = entry.get("comparison", {}).get(key, {})
    return _safe_float(row.get(field))


def _metrics_value(entry: dict[str, Any], key: str) -> float:
    return _safe_float(entry.get("metrics", {}).get(key))


def _load_world_summary(entry: dict[str, Any]) -> dict[str, Any]:
    assets_dir = entry.get("metrics", {}).get("assets_dir")
    if not assets_dir:
        return {}
    path = Path(str(assets_dir)) / "summary.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def _nested_float(root: dict[str, Any], path: tuple[str, ...]) -> float:
    cur: Any = root
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return float("nan")
        cur = cur[key]
    return _safe_float(cur)


def _summer(lat: np.ndarray, seasonal: np.ndarray) -> np.ndarray:
    return np.where(np.asarray(lat, dtype=np.float64) >= 0.0, seasonal[2], seasonal[0])


def _winter(lat: np.ndarray, seasonal: np.ndarray) -> np.ndarray:
    return np.where(np.asarray(lat, dtype=np.float64) >= 0.0, seasonal[0], seasonal[2])


def _f3_seasonal_array_metrics(metrics: dict[str, Any]) -> dict[str, float]:
    assets_dir = metrics.get("assets_dir")
    arrays = metrics.get("arrays")
    path = Path(str(arrays)) if arrays else (
        Path(str(assets_dir)) / "terminal_climate_arrays.npz" if assets_dir else None
    )
    if path is None or not path.exists():
        return {}
    try:
        with np.load(path, allow_pickle=False) as z:
            required = (
                "lat",
                "terrain__elevation_m",
                "sea_level_m",
                "atmosphere__moisture_access",
                "atmosphere__monsoon_potential",
            )
            if any(key not in z.files for key in required):
                return {}
            lat = np.asarray(z["lat"], dtype=np.float64)
            sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
            land = np.asarray(z["terrain__elevation_m"], dtype=np.float64) >= sea
            moisture = np.asarray(z["atmosphere__moisture_access"], dtype=np.float64)
            monsoon = np.asarray(z["atmosphere__monsoon_potential"], dtype=np.float64)
    except (OSError, ValueError):
        return {}
    if moisture.shape != monsoon.shape or moisture.shape != (4, lat.size):
        return {}
    low_mid = land & (np.abs(lat) >= 5.0) & (np.abs(lat) <= 35.0)
    if not low_mid.any():
        return {}

    def pctl(values: np.ndarray, q: float) -> float:
        vals = np.asarray(values, dtype=np.float64)[low_mid]
        vals = vals[np.isfinite(vals)]
        return float(np.percentile(vals, q)) if vals.size else float("nan")

    summer_moisture = _summer(lat, moisture)
    summer_monsoon = _summer(lat, monsoon)
    winter_monsoon = _winter(lat, monsoon)
    return {
        "moisture_access_summer_p75_low_mid": pctl(summer_moisture, 75),
        "monsoon_potential_summer_p90_low_mid": pctl(summer_monsoon, 90),
        "monsoon_potential_summer_minus_winter_p75_low_mid": pctl(
            summer_monsoon - winter_monsoon, 75),
    }


def _row_for_entry(entry: dict[str, Any], earth: dict[str, Any]) -> dict[str, Any]:
    metrics = entry.get("metrics", {})
    world_summary = _load_world_summary(entry)
    climate_diag = world_summary.get("climate_diagnostics", {})
    circulation = climate_diag.get("circulation", {})
    precipitation = climate_diag.get("precipitation", {})
    step_diag = world_summary.get("climate_step_diagnostics", {})
    f3_arrays = _f3_seasonal_array_metrics(metrics)

    row = {
        "label": entry.get("label"),
        "preset": entry.get("preset"),
        "mode": entry.get("mode"),
        "seed": metrics.get("seed"),
        "earth_distance_score": entry.get("earth_distance_score"),
        "flags": ";".join(entry.get("flags", [])),
        "land_fraction": metrics.get("land_fraction"),
        "global_mean_temperature_C": metrics.get("global_mean_temperature_C"),
        "global_temp_delta_C": _comparison_value(
            entry, "global_mean_temperature_C", "delta"),
        "land_mean_temperature_C": metrics.get("land_mean_temperature_C"),
        "land_temp_delta_C": _comparison_value(
            entry, "land_mean_temperature_C", "delta"),
        "ocean_mean_temperature_C": metrics.get("ocean_mean_temperature_C"),
        "ocean_temp_delta_C": _comparison_value(
            entry, "ocean_mean_temperature_C", "delta"),
        "land_precip_mean_mm_yr": metrics.get("land_precip_mean_mm_yr"),
        "land_precip_mean_ratio_to_earth": _ratio(
            metrics.get("land_precip_mean_mm_yr"),
            earth.get("land_precip_mean_mm_yr"),
        ),
        "land_precip_p50_mm_yr": metrics.get("land_precip_p50_mm_yr"),
        "land_precip_p50_ratio_to_earth": _ratio(
            metrics.get("land_precip_p50_mm_yr"),
            earth.get("land_precip_p50_mm_yr"),
        ),
        "land_precip_p90_mm_yr": metrics.get("land_precip_p90_mm_yr"),
        "land_precip_p90_ratio_to_earth": _ratio(
            metrics.get("land_precip_p90_mm_yr"),
            earth.get("land_precip_p90_mm_yr"),
        ),
        "precip_seasonality_land_p75": metrics.get("precip_seasonality_land_p75"),
        "current_speed_p90_m_s": metrics.get("current_speed_p90_m_s"),
        "current_speed_p90_ratio_to_earth": _ratio(
            metrics.get("current_speed_p90_m_s"),
            earth.get("current_speed_p90_m_s"),
        ),
        "biome_desert_area_fraction": metrics.get("biome_desert_area_fraction"),
        "biome_forest_area_fraction": metrics.get("biome_forest_area_fraction"),
        "biome_tropical_area_fraction": metrics.get("biome_tropical_area_fraction"),
        "moisture_access_land_p75": _safe_float(
            circulation.get("moisture_access_land_p75"),
            _safe_float(step_diag.get("moisture_access_land_p75")),
        ),
        "monsoon_potential_land_p90": _safe_float(
            circulation.get("monsoon_potential_land_p90")),
        "monsoon_potential_land_p99": _safe_float(
            circulation.get("monsoon_potential_land_p99")),
        "moisture_access_summer_p75_low_mid": _safe_float(
            f3_arrays.get("moisture_access_summer_p75_low_mid")),
        "monsoon_potential_summer_p90_low_mid": _safe_float(
            f3_arrays.get("monsoon_potential_summer_p90_low_mid")),
        "monsoon_potential_summer_minus_winter_p75_low_mid": _safe_float(
            f3_arrays.get("monsoon_potential_summer_minus_winter_p75_low_mid")),
        "source_ocean_warmth_ocean_p75": _safe_float(
            circulation.get("source_ocean_warmth_ocean_p75")),
        "terrain_blocking_land_p75": _safe_float(
            circulation.get("terrain_blocking_land_p75")),
        "land_wet_fraction_gt500mm": _safe_float(
            precipitation.get("land_wet_fraction_gt500mm")),
        "precip_orographic_concentration": _safe_float(
            precipitation.get("precip_orographic_concentration")),
        "land_monsoon_index_p90": _safe_float(
            step_diag.get("land_monsoon_index_p90")),
    }
    return row


def _phase_scores(earthlike_rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    if not earthlike_rows:
        return {
            phase: {
                "score": float("nan"),
                "priority": "none",
                "status": "no_earthlike_runs",
                "evidence": [],
            }
            for phase in PHASES
        }

    global_temp_abs = [
        abs(_safe_float(row.get("global_temp_delta_C"))) / 6.0
        for row in earthlike_rows
    ]
    land_temp_abs = [
        abs(_safe_float(row.get("land_temp_delta_C"))) / 6.0
        for row in earthlike_rows
    ]
    ocean_temp_abs = [
        abs(_safe_float(row.get("ocean_temp_delta_C"))) / 5.0
        for row in earthlike_rows
    ]
    f1_score = _mean(global_temp_abs + land_temp_abs + ocean_temp_abs)

    current_ratio = [
        abs(_safe_float(row.get("current_speed_p90_ratio_to_earth")) - 1.0)
        for row in earthlike_rows
    ]
    f2_score = _mean(current_ratio)

    moisture_p75 = _mean([
        _safe_float(row.get("moisture_access_land_p75"))
        for row in earthlike_rows
    ])
    monsoon_p90 = _mean([
        _safe_float(row.get("monsoon_potential_land_p90"))
        for row in earthlike_rows
    ])
    summer_moisture_p75 = _mean([
        _safe_float(row.get("moisture_access_summer_p75_low_mid"))
        for row in earthlike_rows
    ])
    summer_monsoon_p90 = _mean([
        _safe_float(row.get("monsoon_potential_summer_p90_low_mid"))
        for row in earthlike_rows
    ])
    summer_winter_monsoon_p75 = _mean([
        _safe_float(row.get("monsoon_potential_summer_minus_winter_p75_low_mid"))
        for row in earthlike_rows
    ])
    # F3 is seasonal by construction.  Prefer the same low/mid-latitude summer
    # diagnostics used by the monsoon/moisture gate; fall back to the older
    # all-season land summary only for legacy comparison reports without arrays.
    if (
        np.isfinite(summer_moisture_p75)
        and np.isfinite(summer_monsoon_p90)
        and np.isfinite(summer_winter_monsoon_p75)
    ):
        moisture_gap = max(0.0, 0.65 - summer_moisture_p75) / 0.65
        monsoon_gap = max(0.0, 0.20 - summer_monsoon_p90) / 0.20
        seasonality_gap = max(0.0, 0.18 - summer_winter_monsoon_p75) / 0.18
        f3_score = float(np.clip(
            0.35 * moisture_gap + 0.45 * monsoon_gap + 0.20 * seasonality_gap,
            0.0,
            3.0,
        ))
        f3_evidence = [
            f"mean low/mid-lat summer moisture-access p75 {summer_moisture_p75:.3f}",
            f"mean low/mid-lat summer monsoon-potential p90 {summer_monsoon_p90:.3f}",
            (
                "mean low/mid-lat summer-minus-winter monsoon-potential p75 "
                f"{summer_winter_monsoon_p75:.3f}"
            ),
        ]
    else:
        moisture_gap = max(0.0, 0.62 - moisture_p75) / 0.62 if np.isfinite(moisture_p75) else 0.5
        monsoon_gap = max(0.0, 0.14 - monsoon_p90) / 0.14 if np.isfinite(monsoon_p90) else 0.5
        f3_score = float(np.clip(0.55 * moisture_gap + 0.45 * monsoon_gap, 0.0, 3.0))
        f3_evidence = [
            f"mean land moisture-access p75 {moisture_p75:.3f}",
            f"mean monsoon-potential land p90 {monsoon_p90:.3f}",
        ]

    precip_ratios = [
        _safe_float(row.get("land_precip_mean_ratio_to_earth"))
        for row in earthlike_rows
    ]
    p50_ratios = [
        _safe_float(row.get("land_precip_p50_ratio_to_earth"))
        for row in earthlike_rows
    ]
    p90_ratios = [
        _safe_float(row.get("land_precip_p90_ratio_to_earth"))
        for row in earthlike_rows
    ]
    precip_gap = [
        max(0.0, 0.45 - r) / 0.45 for r in precip_ratios + p50_ratios + p90_ratios
        if np.isfinite(r)
    ]
    f4_score = _mean(precip_gap)

    forest_fraction = _mean([
        _safe_float(row.get("biome_forest_area_fraction"))
        for row in earthlike_rows
    ])
    tropical_fraction = _mean([
        _safe_float(row.get("biome_tropical_area_fraction"))
        for row in earthlike_rows
    ])
    desert_excess = [
        max(0.0, _safe_float(row.get("biome_desert_area_fraction")) - 0.12)
        for row in earthlike_rows
    ]
    f5_score = _mean(desert_excess + [
        max(0.0, 0.02 - forest_fraction) / 0.02,
        max(0.0, 0.01 - tropical_fraction) / 0.01,
    ])

    hydro_severe = bool(_mean(precip_ratios) < 0.45)
    biome_blocked = bool(hydro_severe and (forest_fraction < 0.005 or tropical_fraction < 0.005))

    phases = {
        "F1_temperature_energy": {
            "score": f1_score,
            "priority": _priority(f1_score, low=0.45, high=1.0),
            "status": "watch" if f1_score < 1.0 else "needs_tuning",
            "evidence": [
                f"mean absolute global/land/ocean temperature normalized score {f1_score:.2f}",
            ],
        },
        "F2_sst_ocean_currents": {
            "score": f2_score,
            "priority": _priority(f2_score, low=0.35, high=1.0),
            "status": "needs_gate_before_hydro_tuning" if f2_score >= 0.35 else "watch",
            "evidence": [
                f"mean current p90 ratio deviation from Earth {f2_score:.2f}",
            ],
        },
        "F3_circulation_moisture_access": {
            "score": f3_score,
            "priority": _priority(f3_score, low=0.25, high=0.60),
            "status": "needs_tuning" if f3_score >= 0.25 else "watch",
            "evidence": f3_evidence,
        },
        "F4_seasonal_hydroclimate": {
            "score": f4_score,
            "priority": _priority(f4_score, low=0.20, high=0.45),
            "status": "dominant_blocker" if hydro_severe else "watch",
            "evidence": [
                f"mean land precipitation ratio to Earth {_mean(precip_ratios):.3f}",
                f"mean p50 precipitation ratio to Earth {_mean(p50_ratios):.3f}",
                f"mean p90 precipitation ratio to Earth {_mean(p90_ratios):.3f}",
            ],
        },
        "F5_koppen_biomes": {
            "score": f5_score,
            "priority": "blocked" if biome_blocked else _priority(f5_score, 0.25, 0.65),
            "status": (
                "blocked_by_hydroclimate"
                if biome_blocked
                else "needs_tuning" if f5_score >= 0.25
                else "watch"
            ),
            "evidence": [
                f"mean forest fraction {forest_fraction:.3f}",
                f"mean tropical fraction {tropical_fraction:.3f}",
                "biome thresholds should not be tuned while F4 is severe"
                if biome_blocked
                else (
                    "biome envelopes are within the current guardrail"
                    if f5_score < 0.25
                    else "biome tuning can follow climate fitting"
                ),
            ],
        },
    }
    return phases


def _priority(score: float, low: float, high: float) -> str:
    if not np.isfinite(score):
        return "unknown"
    if score >= high:
        return "high"
    if score >= low:
        return "medium"
    return "low"


def _candidate_levers(phases: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    f5_status = phases.get("F5_koppen_biomes", {}).get("status")
    f3_status = phases.get("F3_circulation_moisture_access", {}).get("status")
    return [
        {
            "phase": "F1_temperature_energy",
            "priority": phases["F1_temperature_energy"]["priority"],
            "code_area": "aevum.modules.climate.ClimateModule.step / _seasonal_temperature",
            "suggested_lever": (
                "gate annual temperature and land/ocean bias before changing "
                "precipitation coefficients"
            ),
            "guardrail": "global mean within about 6 C; do not warm oceans only to hide dry bias",
        },
        {
            "phase": "F2_sst_ocean_currents",
            "priority": phases["F2_sst_ocean_currents"]["priority"],
            "code_area": "ClimateModule._ocean_currents / _apply_ocean_current_hydro_adjustment",
            "suggested_lever": (
                "separate vector current speed caps from heat-transport and "
                "cold-current drying scalings"
            ),
            "guardrail": "currents remain basin-confined; heat anomaly global mean stays near zero",
        },
        {
            "phase": "F3_circulation_moisture_access",
            "priority": phases["F3_circulation_moisture_access"]["priority"],
            "code_area": "ClimateModule._seasonal_pressure_moisture / _advective_moisture_access",
            "suggested_lever": (
                "keep low/mid-latitude summer moisture and monsoon diagnostics as "
                "guardrails; avoid circulation retuning unless downstream gates regress"
                if f3_status == "watch"
                else (
                    "increase physically routed moisture access and monsoon potential "
                    "from warm source oceans, onshore flow, and passable terrain corridors"
                )
            ),
            "guardrail": "waterworld monsoon stays near zero; arid interiors remain dry",
        },
        {
            "phase": "F4_seasonal_hydroclimate",
            "priority": phases["F4_seasonal_hydroclimate"]["priority"],
            "code_area": "ClimateModule._seasonal_hydroclimate",
            "suggested_lever": (
                "retune land precipitation rainout, ITCZ/storm-track contribution, "
                "and monsoon/convergence coefficients after F3 diagnostics improve"
            ),
            "guardrail": "p50 and p90 land precipitation rise together; no narrow ridge-rain artifacts",
        },
        {
            "phase": "F5_koppen_biomes",
            "priority": phases["F5_koppen_biomes"]["priority"],
            "code_area": "aevum.modules.biosphere.BiosphereModule",
            "suggested_lever": (
                "defer biome thresholds until F4 precipitation reaches plausible ranges"
                if f5_status == "blocked_by_hydroclimate"
                else (
                    "keep current dry/cold stress thresholds as guardrails; do not "
                    "retune biomes while coarse, spatial, seasonal, mountain, and "
                    "windward/leeward gates remain green"
                    if f5_status == "watch"
                    else "retune dry/cold stress thresholds against Koppen and RESOLVE proxies"
                )
            ),
            "guardrail": "do not create forest/tropical classes by ignoring climate stress",
        },
    ]


def _guardrail_check(
    group: str,
    label: str,
    metric: str,
    value: float,
    op: str,
    threshold: float,
    severity: str,
    message: str,
) -> dict[str, Any]:
    if not np.isfinite(value):
        passed = True
        skipped = True
    elif op == ">=":
        passed = value >= threshold
        skipped = False
    elif op == "<=":
        passed = value <= threshold
        skipped = False
    elif op == ">":
        passed = value > threshold
        skipped = False
    elif op == "<":
        passed = value < threshold
        skipped = False
    else:
        raise ValueError(f"unsupported guardrail operator {op!r}")
    return {
        "group": group,
        "label": label,
        "metric": metric,
        "value": value,
        "operator": op,
        "threshold": threshold,
        "severity": severity,
        "passed": bool(passed),
        "skipped": bool(skipped),
        "message": message,
    }


def _guardrail_assessment(rows: list[dict[str, Any]]) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    for row in rows:
        label = str(row.get("label", "unknown"))
        preset = str(row.get("preset", "")).lower()
        mode = str(row.get("mode", ""))
        if mode == "earthlike_calibration":
            checks.extend([
                _guardrail_check(
                    "earthlike",
                    label,
                    "land_precip_mean_ratio_to_earth",
                    _safe_float(row.get("land_precip_mean_ratio_to_earth")),
                    ">=",
                    0.45,
                    "fail",
                    "earthlike land precipitation must clear the dry-bias floor",
                ),
                _guardrail_check(
                    "earthlike",
                    label,
                    "land_precip_p90_ratio_to_earth",
                    _safe_float(row.get("land_precip_p90_ratio_to_earth")),
                    ">=",
                    0.50,
                    "warn",
                    "earthlike high-rainfall tail remains thin",
                ),
                _guardrail_check(
                    "earthlike",
                    label,
                    "global_temp_delta_abs_C",
                    abs(_safe_float(row.get("global_temp_delta_C"))),
                    "<=",
                    6.0,
                    "fail",
                    "earthlike global temperature must remain near Earth envelope",
                ),
                _guardrail_check(
                    "earthlike",
                    label,
                    "current_speed_p90_ratio_to_earth",
                    _safe_float(row.get("current_speed_p90_ratio_to_earth")),
                    "<=",
                    1.45,
                    "warn",
                    "surface current p90 should stay close to drifter envelope",
                ),
            ])
        elif "arid" in preset:
            checks.extend([
                _guardrail_check(
                    "arid",
                    label,
                    "biome_desert_area_fraction",
                    _safe_float(row.get("biome_desert_area_fraction")),
                    ">=",
                    0.45,
                    "fail",
                    "arid worlds must remain desert-dominated",
                ),
                _guardrail_check(
                    "arid",
                    label,
                    "land_precip_p50_mm_yr",
                    _safe_float(row.get("land_precip_p50_mm_yr")),
                    "<=",
                    180.0,
                    "fail",
                    "arid median land precipitation must remain low",
                ),
                _guardrail_check(
                    "arid",
                    label,
                    "biome_tropical_area_fraction",
                    _safe_float(row.get("biome_tropical_area_fraction")),
                    "<=",
                    0.02,
                    "warn",
                    "arid worlds should not grow broad tropical biome patches",
                ),
                _guardrail_check(
                    "arid",
                    label,
                    "land_wet_fraction_gt500mm",
                    _safe_float(row.get("land_wet_fraction_gt500mm")),
                    "<=",
                    0.15,
                    "warn",
                    "arid worlds should retain mostly dry land area",
                ),
            ])
        elif "waterworld" in preset:
            checks.extend([
                _guardrail_check(
                    "waterworld",
                    label,
                    "land_fraction",
                    _safe_float(row.get("land_fraction")),
                    "<=",
                    0.12,
                    "fail",
                    "waterworld exposed land should remain small",
                ),
                _guardrail_check(
                    "waterworld",
                    label,
                    "biome_desert_area_fraction",
                    _safe_float(row.get("biome_desert_area_fraction")),
                    "<=",
                    0.05,
                    "fail",
                    "waterworld islands should not become broad deserts",
                ),
                _guardrail_check(
                    "waterworld",
                    label,
                    "monsoon_potential_land_p90",
                    _safe_float(row.get("monsoon_potential_land_p90")),
                    "<=",
                    0.25,
                    "warn",
                    "waterworld should not develop strong continent-scale monsoon potential",
                ),
                _guardrail_check(
                    "waterworld",
                    label,
                    "land_precip_p50_mm_yr",
                    _safe_float(row.get("land_precip_p50_mm_yr")),
                    ">=",
                    250.0,
                    "warn",
                    "waterworld exposed land should stay maritime rather than arid",
                ),
            ])

    failures = [row for row in checks if not row["passed"] and row["severity"] == "fail"]
    warnings = [row for row in checks if not row["passed"] and row["severity"] == "warn"]
    skipped = [row for row in checks if row.get("skipped")]
    if failures:
        verdict = "fail"
    elif warnings:
        verdict = "pass_with_warnings"
    else:
        verdict = "pass"
    return {
        "verdict": verdict,
        "failure_count": len(failures),
        "warning_count": len(warnings),
        "skipped_count": len(skipped),
        "checks": checks,
        "failures": failures,
        "warnings": warnings,
        "skipped": skipped,
    }


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
        "# Earth-Based Climate Fitting Report",
        "",
        f"Comparison source: `{report['comparison_summary_json']}`",
        "",
        "## Dominant Conclusion",
        "",
        report["dominant_conclusion"],
        "",
        "## Phase Priorities",
        "",
    ]
    for phase in PHASES:
        row = report["phase_assessment"][phase]
        lines.extend([
            f"### {phase}",
            "",
            f"- Priority: `{row['priority']}`",
            f"- Status: `{row['status']}`",
            f"- Score: `{row['score']:.3f}`" if np.isfinite(row["score"]) else "- Score: `nan`",
        ])
        for item in row.get("evidence", []):
            lines.append(f"- Evidence: {item}")
        lines.append("")

    guard = report.get("guardrail_assessment", {})
    lines.extend([
        "## Cross-Preset Guardrails",
        "",
        f"- Verdict: `{guard.get('verdict', 'unknown')}`",
        f"- Failures: `{guard.get('failure_count', 0)}`",
        f"- Warnings: `{guard.get('warning_count', 0)}`",
        "",
    ])
    for item in guard.get("failures", []):
        lines.append(
            f"- FAIL `{item['label']}` {item['metric']}={item['value']:.3f}: "
            f"{item['message']}"
        )
    for item in guard.get("warnings", []):
        lines.append(
            f"- WARN `{item['label']}` {item['metric']}={item['value']:.3f}: "
            f"{item['message']}"
        )
    lines.append("")

    lines.extend(["## Earthlike Runs", ""])
    for row in report["generated_rows"]:
        if row["mode"] != "earthlike_calibration":
            continue
        lines.extend([
            f"### {row['label']}",
            "",
            f"- Land precipitation mean ratio to Earth: `{row['land_precip_mean_ratio_to_earth']:.3f}`",
            f"- Land precipitation p50/p90 ratios: "
            f"`{row['land_precip_p50_ratio_to_earth']:.3f}` / "
            f"`{row['land_precip_p90_ratio_to_earth']:.3f}`",
            f"- Global temperature delta: `{row['global_temp_delta_C']:.2f} C`",
            f"- Current p90 ratio to Earth: `{row['current_speed_p90_ratio_to_earth']:.2f}`",
            f"- Desert/forest/tropical fractions: "
            f"`{row['biome_desert_area_fraction']:.3f}` / "
            f"`{row['biome_forest_area_fraction']:.3f}` / "
            f"`{row['biome_tropical_area_fraction']:.3f}`",
            f"- Moisture access p75 / monsoon potential p90: "
            f"`{row['moisture_access_land_p75']:.3f}` / "
            f"`{row['monsoon_potential_land_p90']:.3f}`",
            f"- Low/mid-lat summer moisture p75 / monsoon p90 / monsoon seasonality p75: "
            f"`{row['moisture_access_summer_p75_low_mid']:.3f}` / "
            f"`{row['monsoon_potential_summer_p90_low_mid']:.3f}` / "
            f"`{row['monsoon_potential_summer_minus_winter_p75_low_mid']:.3f}`",
            "",
        ])

    lines.extend(["## Candidate Levers", ""])
    for lever in report["candidate_levers"]:
        lines.extend([
            f"- `{lever['phase']}` ({lever['priority']}): {lever['suggested_lever']}",
            f"  Code: `{lever['code_area']}`",
            f"  Guardrail: {lever['guardrail']}",
        ])
    lines.append("")
    return "\n".join(lines)


def run_earth_climate_fitting_report(
    config: EarthClimateFittingConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    comparison_path = Path(config.comparison_summary_json)
    comparison = json.loads(comparison_path.read_text())
    earth = comparison.get("earth_metrics", {})
    rows = [
        _row_for_entry(entry, earth)
        for entry in comparison.get("entries", [])
    ]
    rows.sort(key=lambda row: (
        row.get("mode") != "earthlike_calibration",
        str(row.get("label")),
    ))
    earthlike_rows = [
        row for row in rows if row.get("mode") == "earthlike_calibration"
    ]
    phases = _phase_scores(earthlike_rows)
    guardrails = _guardrail_assessment(rows)
    levers = _candidate_levers(phases)
    hydro_status = phases["F4_seasonal_hydroclimate"]["status"]
    f5_status = phases["F5_koppen_biomes"]["status"]
    if guardrails["verdict"] == "fail":
        conclusion = (
            "Cross-preset guardrails failed.  Do not tune further until the "
            "listed arid/waterworld/earthlike regressions are fixed."
        )
    elif hydro_status == "dominant_blocker":
        conclusion = (
            "The first climate-side fix should target circulation, moisture "
            "access, and seasonal hydroclimate.  Biome thresholds are blocked "
            "because earthlike land precipitation is far below the Earth envelope."
        )
    elif f5_status == "blocked_by_hydroclimate":
        conclusion = (
            "Hydroclimate remains too weak for biome calibration; finish F3/F4 "
            "before retuning Koppen or biome thresholds."
        )
    else:
        conclusion = (
            "No single severe blocker dominates; proceed in phase order and keep "
            "comparison metrics under versioned reports."
        )

    row_keys = [
        "label", "preset", "mode", "seed", "earth_distance_score", "flags",
        "land_fraction", "global_mean_temperature_C", "global_temp_delta_C",
        "land_mean_temperature_C", "land_temp_delta_C",
        "ocean_mean_temperature_C", "ocean_temp_delta_C",
        "land_precip_mean_mm_yr", "land_precip_mean_ratio_to_earth",
        "land_precip_p50_mm_yr", "land_precip_p50_ratio_to_earth",
        "land_precip_p90_mm_yr", "land_precip_p90_ratio_to_earth",
        "precip_seasonality_land_p75", "current_speed_p90_m_s",
        "current_speed_p90_ratio_to_earth", "biome_desert_area_fraction",
        "biome_forest_area_fraction", "biome_tropical_area_fraction",
        "moisture_access_land_p75", "monsoon_potential_land_p90",
        "monsoon_potential_land_p99",
        "moisture_access_summer_p75_low_mid",
        "monsoon_potential_summer_p90_low_mid",
        "monsoon_potential_summer_minus_winter_p75_low_mid",
        "source_ocean_warmth_ocean_p75",
        "terrain_blocking_land_p75", "land_wet_fraction_gt500mm",
        "precip_orographic_concentration", "land_monsoon_index_p90",
    ]
    csv_path = _write_csv(outdir / "earth_climate_fitting_runs.csv", rows, row_keys)
    levers_path = _write_csv(
        outdir / "earth_climate_fitting_levers.csv",
        levers,
        ["phase", "priority", "code_area", "suggested_lever", "guardrail"],
    )
    guardrails_path = _write_csv(
        outdir / "earth_climate_guardrails.csv",
        guardrails["checks"],
        [
            "group", "label", "metric", "value", "operator", "threshold",
            "severity", "passed", "skipped", "message",
        ],
    )
    report = {
        "schema": SCHEMA,
        "comparison_summary_json": str(comparison_path),
        "earth_reference_npz": comparison.get("earth_reference_npz"),
        "earthlike_run_count": len(earthlike_rows),
        "dominant_conclusion": conclusion,
        "phase_assessment": phases,
        "guardrail_assessment": guardrails,
        "overall_verdict": guardrails["verdict"],
        "candidate_levers": levers,
        "generated_rows": rows,
        "runs_csv": str(csv_path),
        "levers_csv": str(levers_path),
        "guardrails_csv": str(guardrails_path),
    }
    md_path = outdir / "earth_climate_fitting_report.md"
    report["markdown"] = str(md_path)
    (outdir / "earth_climate_fitting_report.json").write_text(
        json.dumps(report, indent=2, default=_json_default)
    )
    md_path.write_text(_render_markdown(report))
    return report
