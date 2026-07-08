"""Small real-Earth hypsometry fixture validation for P77.

This module validates a derived metric fixture.  It does not download or bundle
raw ETOPO/GEBCO rasters; direct raster extraction is a later Stage A task that
must regenerate the same small schema.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from aevum.diagnostics.earth_reference import REFERENCE_ENVELOPES


SCHEMA = "aevum.real_earth_hypsometry_fixture.v1"
SUMMARY_SCHEMA = "aevum.real_earth_hypsometry_extraction.v1"
FIXTURE_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "reference"
    / "earth_hypsometry_fixture_20260627.json"
)

REQUIRED_METRICS = (
    "land_fraction",
    "ocean_fraction",
    "land_elevation_mean_m",
    "land_elevation_p50_m",
    "land_elevation_p95_m",
    "high_land_fraction_gt2500m",
    "lowland_fraction_lt500m",
    "lowland_fraction_lt1000m",
    "ocean_depth_mean_m",
    "shelf_fraction_of_ocean",
    "slope_rise_fraction_of_ocean",
    "abyss_fraction_of_ocean",
    "trench_and_hadal_fraction_of_ocean",
    "shelf_to_abyss_depth_delta_m",
)


def _load_fixture(path: Path = FIXTURE_PATH) -> dict[str, Any]:
    return json.loads(path.read_text())


def _fixture_digest(fixture: dict[str, Any]) -> str:
    payload = json.dumps(fixture, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _sum_bins(fixture: dict[str, Any], *, land: bool | None = None) -> float:
    total = 0.0
    for entry in fixture["bins"]:
        lo = float(entry["min_m"])
        hi = float(entry["max_m"])
        if land is True and hi <= 0.0:
            continue
        if land is False and lo >= 0.0:
            continue
        total += float(entry["area_fraction"])
    return total


def _envelope_checks(metrics: dict[str, float]) -> dict[str, dict[str, Any]]:
    checked_keys = (
        "land_fraction",
        "land_elevation_mean_m",
        "land_elevation_p95_m",
        "high_land_fraction_gt2500m",
        "shelf_fraction_of_ocean",
        "slope_rise_fraction_of_ocean",
        "abyss_fraction_of_ocean",
        "shelf_to_abyss_depth_delta_m",
    )
    checks: dict[str, dict[str, Any]] = {}
    for key in checked_keys:
        spec = REFERENCE_ENVELOPES[key]
        lo, hi = spec["range"]
        value = float(metrics[key])
        checks[key] = {
            "value": value,
            "reference_min": float(lo),
            "reference_max": float(hi),
            "in_envelope": bool(float(lo) <= value <= float(hi)),
            "basis": spec["basis"],
        }
    return checks


def real_earth_hypsometry_fixture_summary(
    path: Path = FIXTURE_PATH,
) -> dict[str, Any]:
    fixture = _load_fixture(path)
    metrics = {key: float(value) for key, value in fixture["derived_metrics"].items()}
    bin_total = _sum_bins(fixture)
    land_bin_total = _sum_bins(fixture, land=True)
    ocean_bin_total = _sum_bins(fixture, land=False)
    missing_metrics = tuple(key for key in REQUIRED_METRICS if key not in metrics)
    source_ids = tuple(str(source_id) for source_id in fixture.get("source_ids", ()))
    envelope_checks = _envelope_checks(metrics) if not missing_metrics else {}
    lowland_ordering_ok = bool(
        0.0 < metrics.get("lowland_fraction_lt500m", -1.0)
        < metrics.get("lowland_fraction_lt1000m", -1.0)
        < 1.0
    )
    ocean_partition_sum = (
        metrics.get("shelf_fraction_of_ocean", 0.0)
        + metrics.get("slope_rise_fraction_of_ocean", 0.0)
        + metrics.get("abyss_fraction_of_ocean", 0.0)
        + metrics.get("trench_and_hadal_fraction_of_ocean", 0.0)
    )
    acceptance = {
        "schema_valid": fixture.get("schema") == SCHEMA,
        "source_ids_present": {"NOAA_ETOPO_2022", "GEBCO_GRIDDED_BATHYMETRY"}.issubset(source_ids),
        "raw_raster_not_stored": fixture.get("raw_data_stored_in_repo") is False,
        "small_derived_fixture": bool(
            fixture.get("quality_flags", {}).get("small_derived_fixture", False)),
        "direct_raster_extraction_marked_pending": not bool(
            fixture.get("quality_flags", {}).get("direct_raster_extraction", True)),
        "required_metrics_present": not missing_metrics,
        "bin_area_fractions_sum_to_one": abs(bin_total - 1.0) <= 1.0e-9,
        "land_ocean_bins_match_metrics": (
            abs(land_bin_total - metrics.get("land_fraction", -1.0)) <= 1.0e-9
            and abs(ocean_bin_total - metrics.get("ocean_fraction", -1.0)) <= 1.0e-9
        ),
        "lowland_ordering_plausible": lowland_ordering_ok,
        "ocean_partition_plausible": 0.80 <= ocean_partition_sum <= 0.86,
        "envelope_checks_pass": all(
            check["in_envelope"] for check in envelope_checks.values()),
    }
    ready = all(acceptance.values())
    return {
        "schema": SUMMARY_SCHEMA,
        "status": "hypsometry_fixture_ready" if ready else "hypsometry_fixture_incomplete",
        "fixture_path": str(path),
        "fixture_digest": _fixture_digest(fixture),
        "source_ids": source_ids,
        "bin_count": len(fixture["bins"]),
        "bin_area_total": bin_total,
        "land_bin_area_total": land_bin_total,
        "ocean_bin_area_total": ocean_bin_total,
        "ocean_partition_sum": ocean_partition_sum,
        "missing_required_metrics": missing_metrics,
        "metrics": metrics,
        "envelope_checks": envelope_checks,
        "quality_flags": fixture.get("quality_flags", {}),
        "acceptance": acceptance,
        "notes": (
            "This is a small derived metric fixture, not a bundled ETOPO/GEBCO raster.",
            "Future direct raster extraction should regenerate this schema and tighten the envelopes.",
        ),
    }


def generated_world_hypsometry_metrics(
    world: Any,
    tectonics: dict[str, Any],
) -> dict[str, float]:
    """Extract P77-comparable hypsometry metrics from a generated world."""
    import numpy as np

    grid = world.grid
    area = grid.cell_area
    total_area = max(float(area.sum()), 1.0e-12)
    rel = np.asarray(world.get_field("terrain.elevation_m", 0.0), dtype=np.float64) - float(world.sea_level)
    land = rel >= 0.0
    ocean = ~land
    land_area = max(float(area[land].sum()), 1.0e-12)
    ocean_area = max(float(area[ocean].sum()), 1.0e-12)
    ocean_geo = tectonics.get("ocean_geography", {})
    metrics = {
        "land_fraction": float(area[land].sum() / total_area),
        "ocean_fraction": float(area[ocean].sum() / total_area),
        "land_elevation_mean_m": (
            float(np.average(rel[land], weights=area[land])) if land.any() else 0.0
        ),
        "land_elevation_p50_m": float(np.percentile(rel[land], 50)) if land.any() else 0.0,
        "land_elevation_p95_m": float(np.percentile(rel[land], 95)) if land.any() else 0.0,
        "high_land_fraction_gt2500m": float(area[land & (rel > 2500.0)].sum() / land_area),
        "lowland_fraction_lt500m": float(area[land & (rel < 500.0)].sum() / land_area),
        "lowland_fraction_lt1000m": float(area[land & (rel < 1000.0)].sum() / land_area),
        "ocean_depth_mean_m": (
            float(np.average(-rel[ocean], weights=area[ocean])) if ocean.any() else 0.0
        ),
        "shelf_fraction_of_ocean": float(ocean_geo.get("shelf_fraction_of_ocean", 0.0)),
        "slope_rise_fraction_of_ocean": float(ocean_geo.get("slope_rise_fraction_of_ocean", 0.0)),
        "abyss_fraction_of_ocean": float(ocean_geo.get("abyss_fraction_of_ocean", 0.0)),
        "trench_and_hadal_fraction_of_ocean": float(ocean_geo.get("trench_fraction_of_ocean", 0.0)),
        "shelf_to_abyss_depth_delta_m": float(
            ocean_geo.get("shelf_to_abyss_depth_delta_m", 0.0)),
    }
    return metrics


def compare_metrics_to_fixture(
    generated_metrics: dict[str, float],
    fixture_metrics: dict[str, float],
) -> dict[str, Any]:
    comparable = tuple(key for key in REQUIRED_METRICS if key in generated_metrics and key in fixture_metrics)
    deltas: dict[str, dict[str, float]] = {}
    for key in comparable:
        generated = float(generated_metrics[key])
        reference = float(fixture_metrics[key])
        deltas[key] = {
            "generated": generated,
            "reference": reference,
            "absolute_delta": float(generated - reference),
            "absolute_delta_abs": abs(float(generated - reference)),
        }
    envelope_checks = _envelope_checks(generated_metrics)
    core_keys = (
        "land_elevation_mean_m",
        "land_elevation_p95_m",
        "high_land_fraction_gt2500m",
        "shelf_fraction_of_ocean",
        "slope_rise_fraction_of_ocean",
        "abyss_fraction_of_ocean",
        "shelf_to_abyss_depth_delta_m",
    )
    return {
        "comparable_metric_count": len(comparable),
        "comparable_metrics": comparable,
        "deltas": deltas,
        "envelope_checks": envelope_checks,
        "out_of_envelope": tuple(
            key for key, check in envelope_checks.items() if not check["in_envelope"]),
        "core_hypsometry_envelope_pass": all(
            envelope_checks[key]["in_envelope"] for key in core_keys),
        "land_fraction_close_to_fixture": (
            abs(float(generated_metrics["land_fraction"]) - float(fixture_metrics["land_fraction"])) <= 0.08
        ),
        "lowland_ordering_plausible": (
            0.0 < float(generated_metrics["lowland_fraction_lt500m"])
            < float(generated_metrics["lowland_fraction_lt1000m"])
            < 1.0
        ),
    }
