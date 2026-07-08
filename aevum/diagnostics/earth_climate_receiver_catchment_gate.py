"""C4k receiver-catchment object gate.

This gate checks the diagnostic receiver-side catchment layer that follows C4j
precipitation-response objects.  It does not evaluate the precipitation solver
itself; it verifies that source-basin and local-budget semantics are now stable
enough for a later source-basin -> receiver-catchment water-budget closure.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np

from aevum.diagnostics.earth_climate_hydro_region_gate import (
    _array_path,
    _check,
    _json_default,
    _safe_float,
    _write_csv,
)


SCHEMA = "aevum.earth_climate_receiver_catchment_gate.v3"


@dataclass(frozen=True)
class EarthClimateReceiverCatchmentGateConfig:
    terminal_summary_json: Path
    outdir: Path


def _label(summary_row: dict[str, Any]) -> str:
    return Path(str(summary_row.get("assets_dir", summary_row.get("label", "")))).name


def _receiver_catchment_paths(row: dict[str, Any]) -> list[Path]:
    paths: list[Path] = []
    raw = row.get("receiver_catchments_json")
    if raw:
        paths.append(Path(str(raw)))
    assets = row.get("assets_dir")
    if assets:
        paths.append(Path(str(assets)) / "receiver_catchments.json")
    return paths


def _load_receiver_catchments(row: dict[str, Any]) -> tuple[list[dict[str, Any]], Path | None]:
    for path in _receiver_catchment_paths(row):
        if path.exists():
            payload = json.loads(path.read_text())
            catchments = payload.get("catchments", [])
            if isinstance(catchments, list):
                return [obj for obj in catchments if isinstance(obj, dict)], path
            return [], path
    return [], None


def _area_fraction(
    area: np.ndarray,
    mask: np.ndarray,
    domain: np.ndarray,
) -> float:
    denom = float(np.sum(area[domain]))
    if denom <= 1.0e-12:
        return 0.0
    return float(np.sum(area[mask & domain]) / denom)


def _object_pctl(rows: list[dict[str, Any]], key: str, q: float) -> float:
    vals = [
        _safe_float(row.get(key, float("nan")))
        for row in rows
    ]
    vals = [v for v in vals if np.isfinite(v)]
    return float(np.percentile(vals, q)) if vals else 0.0


def _array_metrics(summary_row: dict[str, Any]) -> dict[str, Any]:
    defaults = {
        "arrays_found": 0.0,
        "arrays_path": "",
        "receiver_id_found": 0.0,
        "receiver_id_shape_ok": 0.0,
        "receiver_id_finite_fraction": 0.0,
        "receiver_id_count_p50": 0.0,
        "receiver_land_coverage_p50": 0.0,
        "receiver_source_attribution_p50": 0.0,
        "receiver_budget_attribution_p50": 0.0,
        "receiver_precip_response_attribution_p50": 0.0,
        "source_supply_found": 0.0,
        "source_supply_shape_ok": 0.0,
        "source_supply_finite_fraction": 0.0,
        "source_supply_land_p50": 0.0,
        "source_supply_attributed_land_p50": 0.0,
        "wet_response_source_supply_p50": 0.0,
        "receiver_supply_balance_found": 0.0,
        "receiver_supply_balance_shape_ok": 0.0,
        "receiver_supply_balance_finite_fraction": 0.0,
        "receiver_supply_balance_land_p10": 0.0,
        "receiver_supply_balance_land_p50": 0.0,
        "receiver_feedback_found": 0.0,
        "receiver_feedback_shape_ok": 0.0,
        "receiver_feedback_finite_fraction": 0.0,
        "receiver_feedback_land_p05": 1.0,
        "receiver_feedback_land_p95": 1.0,
        "receiver_feedback_ocean_abs_p95": 0.0,
    }
    path = _array_path(summary_row)
    if path is None:
        return defaults
    defaults["arrays_found"] = 1.0
    defaults["arrays_path"] = str(path)
    with np.load(path, allow_pickle=False) as z:
        required = ("cell_area", "terrain__elevation_m", "sea_level_m")
        if not all(key in z.files for key in required):
            return defaults
        area = np.asarray(z["cell_area"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"]).ravel()[0])
        land = np.asarray(z["terrain__elevation_m"], dtype=np.float64) >= sea
        if "climate__receiver_catchment_id" not in z.files:
            return defaults
        receiver = np.asarray(z["climate__receiver_catchment_id"], dtype=np.float64)
        defaults["receiver_id_found"] = 1.0
        defaults["receiver_id_shape_ok"] = 1.0 if receiver.shape == (4, area.size) else 0.0
        defaults["receiver_id_finite_fraction"] = (
            float(np.isfinite(receiver).mean()) if receiver.size else 0.0
        )
        if receiver.shape != (4, area.size):
            return defaults
        source = np.asarray(
            z["atmosphere__moisture_source_basin_id"], dtype=np.float64
        ) if "atmosphere__moisture_source_basin_id" in z.files else np.full_like(
            receiver, -1.0)
        budget = np.asarray(
            z["climate__moisture_budget_region_id"], dtype=np.float64
        ) if "climate__moisture_budget_region_id" in z.files else np.full_like(
            receiver, -1.0)
        precip_region = np.asarray(
            z["climate__precipitation_response_region_id"], dtype=np.float64
        ) if "climate__precipitation_response_region_id" in z.files else np.full_like(
            receiver, -1.0)
        response = np.asarray(
            z["climate__moisture_flow_precipitation_response"], dtype=np.float64
        ) if "climate__moisture_flow_precipitation_response" in z.files else None
        source_supply = None
        if "climate__source_basin_supply_index" in z.files:
            source_supply = np.asarray(
                z["climate__source_basin_supply_index"], dtype=np.float64)
            defaults["source_supply_found"] = 1.0
            defaults["source_supply_shape_ok"] = (
                1.0 if source_supply.shape == (4, area.size) else 0.0
            )
            defaults["source_supply_finite_fraction"] = (
                float(np.isfinite(source_supply).mean())
                if source_supply.size else 0.0
            )
        receiver_balance = None
        if "climate__receiver_catchment_supply_balance" in z.files:
            receiver_balance = np.asarray(
                z["climate__receiver_catchment_supply_balance"], dtype=np.float64)
            defaults["receiver_supply_balance_found"] = 1.0
            defaults["receiver_supply_balance_shape_ok"] = (
                1.0 if receiver_balance.shape == (4, area.size) else 0.0
            )
            defaults["receiver_supply_balance_finite_fraction"] = (
                float(np.isfinite(receiver_balance).mean())
                if receiver_balance.size else 0.0
            )
        receiver_feedback = None
        if "climate__receiver_supply_precipitation_feedback" in z.files:
            receiver_feedback = np.asarray(
                z["climate__receiver_supply_precipitation_feedback"],
                dtype=np.float64,
            )
            defaults["receiver_feedback_found"] = 1.0
            defaults["receiver_feedback_shape_ok"] = (
                1.0 if receiver_feedback.shape == (4, area.size) else 0.0
            )
            defaults["receiver_feedback_finite_fraction"] = (
                float(np.isfinite(receiver_feedback).mean())
                if receiver_feedback.size else 0.0
            )
    if not land.any():
        return defaults
    counts: list[int] = []
    coverage: list[float] = []
    source_attrs: list[float] = []
    budget_attrs: list[float] = []
    response_attrs: list[float] = []
    source_supply_values: list[float] = []
    attributed_supply_values: list[float] = []
    wet_supply_values: list[float] = []
    receiver_balance_values: list[float] = []
    receiver_feedback_values: list[float] = []
    receiver_feedback_ocean_abs: list[float] = []
    for season in range(4):
        active = land & np.isfinite(receiver[season]) & (receiver[season] > 0.0)
        counts.append(len({
            int(x) for x in receiver[season, active]
            if np.isfinite(x) and int(x) > 0
        }))
        coverage.append(_area_fraction(area, active, land))
        for rid in [
            int(x) for x in np.unique(receiver[season, active])
            if np.isfinite(x) and int(x) > 0
        ]:
            cells = land & (receiver[season] == float(rid))
            if int(np.count_nonzero(cells)) < 2:
                continue
            source_attrs.append(_area_fraction(
                area,
                cells & np.isfinite(source[season]) & (source[season] >= 0.0),
                cells,
            ))
            budget_attrs.append(_area_fraction(
                area,
                cells & np.isfinite(budget[season]) & (budget[season] > 0.0),
                cells,
            ))
            response_attrs.append(_area_fraction(
                area,
                cells & np.isfinite(precip_region[season])
                & (precip_region[season] > 0.0),
                cells,
            ))
        if source_supply is not None and source_supply.shape == (4, area.size):
            season_supply = source_supply[season]
            finite_land = land & np.isfinite(season_supply)
            source_supply_values.extend([
                float(x) for x in season_supply[finite_land]
            ])
            attributed = (
                finite_land
                & np.isfinite(source[season])
                & (source[season] >= 0.0)
            )
            attributed_supply_values.extend([
                float(x) for x in season_supply[attributed]
            ])
            if response is not None and response.shape == (4, area.size):
                wet = (
                    finite_land
                    & np.isfinite(response[season])
                    & (response[season] > 1.02)
                )
                wet_supply_values.extend([
                    float(x) for x in season_supply[wet]
                ])
        if receiver_balance is not None and receiver_balance.shape == (4, area.size):
            finite_balance = land & np.isfinite(receiver_balance[season])
            receiver_balance_values.extend([
                float(x) for x in receiver_balance[season, finite_balance]
            ])
        if receiver_feedback is not None and receiver_feedback.shape == (4, area.size):
            finite_feedback = land & np.isfinite(receiver_feedback[season])
            receiver_feedback_values.extend([
                float(x) for x in receiver_feedback[season, finite_feedback]
            ])
            ocean = ~land
            finite_ocean_feedback = ocean & np.isfinite(receiver_feedback[season])
            receiver_feedback_ocean_abs.extend([
                float(abs(x - 1.0))
                for x in receiver_feedback[season, finite_ocean_feedback]
            ])
    defaults.update({
        "receiver_id_count_p50": float(np.percentile(counts, 50)) if counts else 0.0,
        "receiver_land_coverage_p50": (
            float(np.percentile(coverage, 50)) if coverage else 0.0
        ),
        "receiver_source_attribution_p50": (
            float(np.percentile(source_attrs, 50)) if source_attrs else 0.0
        ),
        "receiver_budget_attribution_p50": (
            float(np.percentile(budget_attrs, 50)) if budget_attrs else 0.0
        ),
        "receiver_precip_response_attribution_p50": (
            float(np.percentile(response_attrs, 50)) if response_attrs else 0.0
        ),
        "source_supply_land_p50": (
            float(np.percentile(source_supply_values, 50))
            if source_supply_values else 0.0
        ),
        "source_supply_attributed_land_p50": (
            float(np.percentile(attributed_supply_values, 50))
            if attributed_supply_values else 0.0
        ),
        "wet_response_source_supply_p50": (
            float(np.percentile(wet_supply_values, 50))
            if wet_supply_values else 0.0
        ),
        "receiver_supply_balance_land_p10": (
            float(np.percentile(receiver_balance_values, 10))
            if receiver_balance_values else 0.0
        ),
        "receiver_supply_balance_land_p50": (
            float(np.percentile(receiver_balance_values, 50))
            if receiver_balance_values else 0.0
        ),
        "receiver_feedback_land_p05": (
            float(np.percentile(receiver_feedback_values, 5))
            if receiver_feedback_values else 1.0
        ),
        "receiver_feedback_land_p95": (
            float(np.percentile(receiver_feedback_values, 95))
            if receiver_feedback_values else 1.0
        ),
        "receiver_feedback_ocean_abs_p95": (
            float(np.percentile(receiver_feedback_ocean_abs, 95))
            if receiver_feedback_ocean_abs else 0.0
        ),
    })
    return defaults


def _archive_metrics(summary_row: dict[str, Any]) -> dict[str, Any]:
    catchments, path = _load_receiver_catchments(summary_row)
    kinds = [str(row.get("kind", "")) for row in catchments]
    seasons = [str(row.get("season", "")) for row in catchments]
    return {
        "receiver_archive_found": 1.0 if path is not None else 0.0,
        "receiver_archive_path": "" if path is None else str(path),
        "receiver_object_count": float(len(catchments)),
        "receiver_kind_count": float(len(set(kinds))),
        "receiver_season_count": float(len(set(seasons))),
        "source_receiver_object_count": float(sum(
            1 for kind in kinds if kind == "source_receiver_catchment")),
        "mixed_receiver_object_count": float(sum(
            1 for kind in kinds if kind == "mixed_receiver_catchment")),
        "receiver_largest_area_fraction": max(
            [_safe_float(row.get("area_fraction", 0.0), 0.0) for row in catchments],
            default=0.0,
        ),
        "receiver_object_source_attribution_p50": _object_pctl(
            catchments, "source_basin_attributed_fraction", 50),
        "receiver_object_source_purity_p50": _object_pctl(
            catchments, "source_basin_purity", 50),
        "receiver_object_budget_attribution_p50": _object_pctl(
            catchments, "budget_region_attributed_fraction", 50),
        "receiver_object_precip_response_attribution_p50": _object_pctl(
            catchments, "precipitation_response_attributed_fraction", 50),
        "receiver_object_precip_response_attribution_p90": _object_pctl(
            catchments, "precipitation_response_attributed_fraction", 90),
        "receiver_object_source_supply_p50": _object_pctl(
            catchments, "mean_source_basin_supply_index", 50),
        "receiver_object_supply_attribution_p50": _object_pctl(
            catchments, "source_basin_supply_attributed_fraction", 50),
        "receiver_object_supply_balance_p50": _object_pctl(
            catchments, "precipitation_supply_balance", 50),
        "receiver_object_supported_precip_fraction_p50": _object_pctl(
            catchments, "supply_supported_precipitation_fraction", 50),
    }


def _metrics_for_row(summary_row: dict[str, Any]) -> dict[str, Any]:
    row = {
        "label": _label(summary_row),
        "preset": str(summary_row.get("preset", "")),
        "seed": int(summary_row.get("seed", -1)),
    }
    row.update(_array_metrics(summary_row))
    row.update(_archive_metrics(summary_row))
    return row


def _checks_for_row(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    checks = [
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="C4k gate requires archived climate arrays"),
        _check(label=label, group="array_archive", metric="receiver_id_found",
               generated=row["receiver_id_found"], operator=">=", threshold=1.0,
               severity="fail", message="climate.receiver_catchment_id must be archived"),
        _check(label=label, group="array_archive", metric="receiver_id_shape_ok",
               generated=row["receiver_id_shape_ok"], operator=">=", threshold=1.0,
               severity="fail", message="receiver catchment id must have four seasonal maps"),
        _check(label=label, group="array_archive", metric="receiver_id_finite_fraction",
               generated=row["receiver_id_finite_fraction"], operator=">=", threshold=1.0,
               severity="fail", message="receiver catchment id must be finite"),
        _check(label=label, group="source_supply", metric="source_supply_found",
               generated=row["source_supply_found"], operator=">=", threshold=1.0,
               severity="fail", message="C5e7 source-basin supply field must be archived"),
        _check(label=label, group="source_supply", metric="source_supply_shape_ok",
               generated=row["source_supply_shape_ok"], operator=">=", threshold=1.0,
               severity="fail", message="source-basin supply field must have four seasonal maps"),
        _check(label=label, group="source_supply", metric="source_supply_finite_fraction",
               generated=row["source_supply_finite_fraction"], operator=">=", threshold=1.0,
               severity="fail", message="source-basin supply field must be finite"),
        _check(label=label, group="receiver_supply",
               metric="receiver_supply_balance_found",
               generated=row["receiver_supply_balance_found"], operator=">=", threshold=1.0,
               severity="fail", message="C5e7 receiver supply-balance field must be archived"),
        _check(label=label, group="receiver_supply",
               metric="receiver_supply_balance_shape_ok",
               generated=row["receiver_supply_balance_shape_ok"],
               operator=">=", threshold=1.0, severity="fail",
               message="receiver supply-balance field must have four seasonal maps"),
        _check(label=label, group="receiver_supply",
               metric="receiver_supply_balance_finite_fraction",
               generated=row["receiver_supply_balance_finite_fraction"],
               operator=">=", threshold=1.0, severity="fail",
               message="receiver supply-balance field must be finite"),
        _check(label=label, group="receiver_feedback",
               metric="receiver_feedback_found",
               generated=row["receiver_feedback_found"], operator=">=", threshold=1.0,
               severity="fail", message="C5e8 receiver-supply feedback field must be archived"),
        _check(label=label, group="receiver_feedback",
               metric="receiver_feedback_shape_ok",
               generated=row["receiver_feedback_shape_ok"], operator=">=", threshold=1.0,
               severity="fail", message="receiver-supply feedback field must have four seasonal maps"),
        _check(label=label, group="receiver_feedback",
               metric="receiver_feedback_finite_fraction",
               generated=row["receiver_feedback_finite_fraction"], operator=">=", threshold=1.0,
               severity="fail", message="receiver-supply feedback field must be finite"),
        _check(label=label, group="receiver_feedback",
               metric="receiver_feedback_land_p05",
               generated=row["receiver_feedback_land_p05"], operator=">=", threshold=0.88,
               severity="fail", message="receiver-supply feedback must stay conservatively bounded"),
        _check(label=label, group="receiver_feedback",
               metric="receiver_feedback_land_p95",
               generated=row["receiver_feedback_land_p95"], operator="<=", threshold=1.14,
               severity="fail", message="receiver-supply feedback must stay conservatively bounded"),
        _check(label=label, group="receiver_feedback",
               metric="receiver_feedback_ocean_abs_p95",
               generated=row["receiver_feedback_ocean_abs_p95"], operator="<=", threshold=1.0e-9,
               severity="fail", message="receiver-supply feedback must leave ocean precipitation untouched"),
        _check(label=label, group="coverage", metric="receiver_land_coverage_p50",
               generated=row["receiver_land_coverage_p50"], operator=">=", threshold=0.80,
               severity="fail", message="receiver catchments should cover most seasonal land"),
        _check(label=label, group="source_supply",
               metric="source_supply_attributed_land_p50",
               generated=row["source_supply_attributed_land_p50"],
               operator=">=", threshold=0.05, severity="fail",
               message="source-attributed receiving land must have diagnosed basin supply"),
        _check(label=label, group="receiver_supply",
               metric="receiver_supply_balance_land_p50",
               generated=row["receiver_supply_balance_land_p50"],
               operator=">=", threshold=0.25, severity="fail",
               message="receiver catchments need nonzero precipitation-supply consistency"),
        _check(label=label, group="object_archive", metric="receiver_archive_found",
               generated=row["receiver_archive_found"], operator=">=", threshold=1.0,
               severity="fail", message="receiver catchment objects must be archived"),
        _check(label=label, group="object_archive", metric="receiver_object_count",
               generated=row["receiver_object_count"], operator=">=", threshold=1.0,
               severity="fail", message="receiver catchment archive should not be empty"),
        _check(label=label, group="object_archive", metric="receiver_season_count",
               generated=row["receiver_season_count"], operator=">=", threshold=4.0,
               severity="fail", message="receiver catchments should cover all seasons"),
        _check(label=label, group="budget_semantics",
               metric="receiver_object_budget_attribution_p50",
               generated=row["receiver_object_budget_attribution_p50"],
               operator=">=", threshold=0.85, severity="fail",
               message="receiver catchments must retain local budget-region attribution"),
        _check(label=label, group="source_semantics",
               metric="receiver_object_source_attribution_p50",
               generated=row["receiver_object_source_attribution_p50"],
               operator=">=", threshold=0.20, severity="warn",
               message="receiver catchments should increasingly carry source-basin attribution"),
        _check(label=label, group="source_supply",
               metric="receiver_object_source_supply_p50",
               generated=row["receiver_object_source_supply_p50"],
               operator=">=", threshold=0.05, severity="fail",
               message="receiver catchment objects should archive source-basin supply support"),
        _check(label=label, group="receiver_supply",
               metric="receiver_object_supply_balance_p50",
               generated=row["receiver_object_supply_balance_p50"],
               operator=">=", threshold=0.20, severity="fail",
               message="receiver catchment objects should archive precipitation-supply balance"),
    ]
    preset = str(row.get("preset", "")).lower()
    if "earthlike" in preset or "earthlike" in label.lower():
        checks.extend([
            _check(label=label, group="response_semantics",
                   metric="receiver_object_precip_response_attribution_p90",
                   generated=row["receiver_object_precip_response_attribution_p90"],
                   operator=">=", threshold=0.05, severity="fail",
                   message="earthlike receiver catchments should bind back to C4j response regions"),
            _check(label=label, group="earthlike_structure",
                   metric="receiver_id_count_p50",
                   generated=row["receiver_id_count_p50"], operator=">=", threshold=2.0,
                   severity="fail",
                   message="earthlike worlds should have multiple receiver catchments"),
            _check(label=label, group="earthlike_structure",
                   metric="source_receiver_object_count",
                   generated=row["source_receiver_object_count"],
                   operator=">=", threshold=1.0, severity="fail",
                   message="earthlike receiver catchments should include source-dominated objects"),
            _check(label=label, group="earthlike_source_supply",
                   metric="wet_response_source_supply_p50",
                   generated=row["wet_response_source_supply_p50"],
                   operator=">=", threshold=0.05, severity="fail",
                   message="earthlike wet response cells should be backed by diagnosed source-basin supply"),
        ])
    if "arid" in preset or "arid" in label.lower():
        checks.append(
            _check(label=label, group="response_semantics",
                   metric="receiver_object_precip_response_attribution_p90",
                   generated=row["receiver_object_precip_response_attribution_p90"],
                   operator=">=", threshold=0.05, severity="fail",
                   message="arid receiver catchments should bind back to C4j response regions")
        )
    if "waterworld" in preset or "waterworld" in label.lower():
        checks.append(
            _check(label=label, group="waterworld_guard",
                   metric="receiver_largest_area_fraction",
                   generated=row["receiver_largest_area_fraction"],
                   operator="<=", threshold=0.12, severity="fail",
                   message="waterworld receiver catchments should remain island scale")
        )
    return checks


def _render_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Earth Climate Receiver-Catchment Gate",
        "",
        f"Schema: `{report['schema']}`",
        f"Verdict: **{report['verdict']}**",
        "",
        (
            f"Failures: `{report['failure_count']}`; warnings: "
            f"`{report['warning_count']}`; skipped: `{report['skipped_count']}`."
        ),
        "",
        "This gate evaluates C4k/C5e7 receiver-side catchment objects.  It checks "
        "that C4j precipitation-response patches have been merged into stable "
        "seasonal receiving domains with local budget-region and source-basin "
        "attribution, verifies that C5e7 source-basin supply and receiver "
        "precipitation-supply balance diagnostics are archived and nonzero, "
        "and checks that the C5e8 receiver-supply precipitation feedback is "
        "archived, bounded, and ocean-neutral.",
        "",
        "## Checks",
        "",
    ]
    for check in report["checks"]:
        mark = "PASS" if check["passed"] else "FAIL"
        lines.append(
            f"- {mark} `{check['label']}` `{check['group']}` "
            f"`{check['metric']}` = `{check['generated']:.3g}` "
            f"{check['operator']} `{check['threshold']}`"
        )
    return "\n".join(lines) + "\n"


def run_earth_climate_receiver_catchment_gate(
    config: EarthClimateReceiverCatchmentGateConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    terminal = json.loads(Path(config.terminal_summary_json).read_text())
    rows = [_metrics_for_row(row) for row in terminal.get("summaries", [])]
    rows.sort(key=lambda row: str(row["label"]))
    checks: list[dict[str, Any]] = []
    for row in rows:
        checks.extend(_checks_for_row(row))
    failures = [
        check for check in checks
        if not check["passed"] and check["severity"] == "fail"
    ]
    warnings = [
        check for check in checks
        if not check["passed"] and check["severity"] == "warn"
    ]
    skipped = [check for check in checks if check.get("skipped")]
    metrics_csv = _write_csv(
        outdir / "earth_climate_receiver_catchment_metrics.csv",
        rows,
        [
            "label", "preset", "seed", "arrays_found", "receiver_id_found",
            "receiver_id_shape_ok", "receiver_id_finite_fraction",
            "receiver_id_count_p50", "receiver_land_coverage_p50",
            "source_supply_found", "source_supply_shape_ok",
            "source_supply_finite_fraction", "source_supply_land_p50",
            "source_supply_attributed_land_p50",
            "wet_response_source_supply_p50",
            "receiver_supply_balance_found",
            "receiver_supply_balance_shape_ok",
            "receiver_supply_balance_finite_fraction",
            "receiver_supply_balance_land_p10",
            "receiver_supply_balance_land_p50",
            "receiver_feedback_found", "receiver_feedback_shape_ok",
            "receiver_feedback_finite_fraction", "receiver_feedback_land_p05",
            "receiver_feedback_land_p95", "receiver_feedback_ocean_abs_p95",
            "receiver_archive_found", "receiver_object_count",
            "receiver_kind_count", "receiver_season_count",
            "source_receiver_object_count", "mixed_receiver_object_count",
            "receiver_largest_area_fraction",
            "receiver_object_source_attribution_p50",
            "receiver_object_source_purity_p50",
            "receiver_object_budget_attribution_p50",
            "receiver_object_precip_response_attribution_p50",
            "receiver_object_precip_response_attribution_p90",
            "receiver_object_source_supply_p50",
            "receiver_object_supply_attribution_p50",
            "receiver_object_supply_balance_p50",
            "receiver_object_supported_precip_fraction_p50",
            "arrays_path", "receiver_archive_path",
        ],
    )
    checks_csv = _write_csv(
        outdir / "earth_climate_receiver_catchment_checks.csv",
        checks,
        [
            "label", "group", "metric", "generated", "operator", "threshold",
            "severity", "passed", "skipped", "message",
        ],
    )
    report = {
        "schema": SCHEMA,
        "terminal_summary_json": str(config.terminal_summary_json),
        "verdict": "fail" if failures else "pass",
        "failure_count": int(len(failures)),
        "warning_count": int(len(warnings)),
        "skipped_count": int(len(skipped)),
        "metrics": rows,
        "checks": checks,
        "failures": failures,
        "warnings": warnings,
        "metrics_csv": str(metrics_csv),
        "checks_csv": str(checks_csv),
    }
    (outdir / "earth_climate_receiver_catchment_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default)
    )
    (outdir / "earth_climate_receiver_catchment_gate_report.md").write_text(
        _render_markdown(report)
    )
    return report
