import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_moisture_response_gate
from aevum.diagnostics.earth_climate_moisture_response_gate import (
    SCHEMA,
    EarthClimateMoistureResponseGateConfig,
    run_earth_climate_moisture_response_gate,
)


def _write_summary(path, row):
    path.write_text(json.dumps({
        "schema": "aevum.terminal_climate_replay.v1",
        "summaries": [row],
    }))


def _summary_row(assets, preset="earthlike_mobile_lid"):
    return {
        "preset": preset,
        "seed": 1,
        "assets_dir": str(assets),
        "precipitation_response_regions_json": str(
            assets / "precipitation_response_regions.json"),
        "climate_step_diagnostics": {
            "moisture_flow_precipitation_response": {
                "enabled": True,
                "response_land_p05": 0.84,
                "response_land_p95": 1.09,
                "max_land_mean_delta_mm_yr": 1.0e-12,
                "budget_base_region_count": 2.0,
                "budget_region_count_p50": 2.0,
                "budget_sector_split_count_p50": 1.0,
                "max_budget_region_mean_delta_mm_yr": 1.0e-12,
            }
        },
        "arrays": str(assets / "terminal_climate_arrays.npz"),
    }


def _grid(n=160):
    idx = np.arange(n, dtype=float) + 0.5
    z = 1.0 - 2.0 * idx / n
    lat = np.degrees(np.arcsin(np.clip(z, -1.0, 1.0)))
    lon = ((np.degrees(np.pi * (1.0 + 5.0 ** 0.5) * idx) + 180.0) % 360.0) - 180.0
    return lat, lon


def _write_arrays(
    path,
    *,
    missing_response=False,
    ocean_changed=False,
    waterworld=False,
    strong_waterworld=False,
):
    lat, lon = _grid()
    n = lat.size
    area = np.ones(n, dtype=float)
    if waterworld:
        land = (np.abs(lat) < 28.0) & (lon > -28.0) & (lon < 28.0)
    else:
        land = (np.abs(lat) < 62.0) & (lon > -150.0)
    elevation = np.where(land, 120.0, -160.0)
    shape = (4, n)
    pathway = np.zeros(shape, dtype=float)
    source_basin_id = np.full(shape, -1.0, dtype=float)
    network_id = np.full(shape, -1.0, dtype=float)
    budget_region_id = np.full(shape, -1.0, dtype=float)
    monsoon = np.zeros(shape, dtype=float)
    storm = np.zeros(shape, dtype=float)
    shadow = np.zeros(shape, dtype=float)
    precip = np.full(shape, 500.0, dtype=float)
    response = np.ones(shape, dtype=float)
    precip_region_id = np.full(shape, -1.0, dtype=float)
    west_land = land & (lon < 0.0)
    east_land = land & (lon >= 0.0)
    budget_region_id[:, west_land] = 1.0
    budget_region_id[:, east_land] = 2.0
    centers = (-24.0, -4.0, 28.0, 6.0)
    for season, center in enumerate(centers):
        monsoon_band = np.exp(-((lat - center) / 20.0) ** 2)
        storm_band = np.exp(-((np.abs(lat) - 42.0) / 11.0) ** 2)
        support = np.maximum(0.95 * monsoon_band, 0.70 * storm_band)
        pathway[season, land] = 0.12 + 0.90 * support[land]
        source_basin_id[season, land & (lon < 0.0)] = 1.0
        source_basin_id[season, land & (lon >= 0.0)] = 2.0
        network_id[season, land & (pathway[season] > 0.42)] = 1.0 + season
        monsoon[season] = monsoon_band
        storm[season] = storm_band
        shadow[season] = np.clip((np.abs(lat) - 12.0) / 70.0, 0.0, 0.8)
        precip[season, land] *= 1.0 + 0.40 * support[land]
        if land.any():
            signal = pathway[season].copy()
            signal -= float(np.mean(signal[land]))
            scale = max(float(np.percentile(np.abs(signal[land]), 85)), 1.0e-9)
            amp = 0.20 if strong_waterworld else (0.055 if waterworld else 0.13)
            response[season, land] = np.clip(1.0 + amp * signal[land] / scale, 0.76, 1.16)
            precip_region_id[season, land & (response[season] > 1.025)] = float(
                10 + season)
            precip_region_id[season, land & (response[season] < 0.975)] = float(
                20 + season)
    if ocean_changed:
        response[:, ~land] = 0.94

    payload = {
        "lat": lat,
        "lon": lon,
        "cell_area": area,
        "sea_level_m": np.array([0.0], dtype=float),
        "terrain__elevation_m": elevation,
        "atmosphere__moisture_flow_pathway": pathway,
        "atmosphere__moisture_source_basin_id": source_basin_id,
        "climate__moisture_flow_network_id": network_id,
        "climate__moisture_budget_region_id": budget_region_id,
        "climate__monsoon_rainfall_corridor": monsoon,
        "climate__storm_track_rainfall_corridor": storm,
        "climate__rain_shadow_index": shadow,
        "climate__seasonal_precipitation": precip,
    }
    if not missing_response:
        payload["climate__moisture_flow_precipitation_response"] = response
        payload["climate__precipitation_response_region_id"] = precip_region_id
    np.savez_compressed(path, **payload)


def _write_precipitation_response_regions(path):
    regions = []
    for season, name in enumerate(("DJF", "MAM", "JJA", "SON")):
        regions.extend([
            {
                "id": f"precipitation_response:{name.lower()}:{10 + season}",
                "type": "precipitation_response_region",
                "kind": "wet_precipitation_response_region",
                "season": name,
                "season_index": season,
                "region_index": 10 + season,
                "cell_count": 12,
                "area_fraction": 0.025,
                "mean_response": 1.08,
                "mean_abs_response_anomaly": 0.08,
                "source_basin_ids": [1],
                "dominant_source_basin_id": 1,
                "source_basin_attributed_fraction": 1.0,
                "source_basin_purity": 1.0,
                "budget_region_ids": [1],
                "dominant_budget_region_id": 1,
                "budget_region_attributed_fraction": 1.0,
                "budget_region_purity": 1.0,
                "flow_network_ids": [1 + season],
                "dominant_flow_network_id": 1 + season,
                "flow_network_attributed_fraction": 1.0,
                "flow_network_purity": 1.0,
            },
            {
                "id": f"precipitation_response:{name.lower()}:{20 + season}",
                "type": "precipitation_response_region",
                "kind": "dry_precipitation_response_region",
                "season": name,
                "season_index": season,
                "region_index": 20 + season,
                "cell_count": 10,
                "area_fraction": 0.020,
                "mean_response": 0.88,
                "mean_abs_response_anomaly": 0.12,
                "source_basin_ids": [2],
                "dominant_source_basin_id": 2,
                "source_basin_attributed_fraction": 1.0,
                "source_basin_purity": 1.0,
                "budget_region_ids": [2],
                "dominant_budget_region_id": 2,
                "budget_region_attributed_fraction": 1.0,
                "budget_region_purity": 1.0,
                "flow_network_ids": [],
                "dominant_flow_network_id": -1,
                "flow_network_attributed_fraction": 0.0,
                "flow_network_purity": 0.0,
            },
        ])
    path.write_text(json.dumps({
        "schema": "aevum.precipitation_response_regions.v1",
        "region_count": len(regions),
        "regions": regions,
    }))


def test_moisture_response_gate_passes_complete_earthlike_response(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    _write_arrays(assets / "terminal_climate_arrays.npz")
    _write_precipitation_response_regions(
        assets / "precipitation_response_regions.json")
    summary = tmp_path / "summary.json"
    _write_summary(summary, _summary_row(assets))

    report = run_earth_climate_moisture_response_gate(
        EarthClimateMoistureResponseGateConfig(summary, tmp_path / "out")
    )

    assert report["schema"] == SCHEMA
    assert report["verdict"] == "pass"
    assert report["failure_count"] == 0
    assert report["contact_sheet_count"] == 1
    row = report["generated_metrics"][0]
    assert row["budget_region_found"] == 1.0
    assert row["budget_region_count_p50"] >= 1.0
    assert row["precip_region_archive_found"] == 1.0
    assert row["precip_region_kind_count"] >= 2.0
    assert (tmp_path / "out" / "earth_climate_moisture_response_metrics.csv").exists()
    assert (tmp_path / "out" / "earth_climate_moisture_response_checks.csv").exists()
    assert (
        tmp_path / "out" / "contact_sheets"
        / "earthlike_seed1_moisture_response_contact_sheet.png"
    ).exists()


def test_moisture_response_gate_flags_missing_response_array(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    _write_arrays(assets / "terminal_climate_arrays.npz", missing_response=True)
    summary = tmp_path / "summary.json"
    _write_summary(summary, _summary_row(assets))

    report = run_earth_climate_moisture_response_gate(
        EarthClimateMoistureResponseGateConfig(
            summary, tmp_path / "out", render_contact_sheets=False)
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "array_archive" and row["metric"] == "response_found"
        for row in report["failures"]
    )


def test_moisture_response_gate_flags_ocean_precip_response_change(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    _write_arrays(assets / "terminal_climate_arrays.npz", ocean_changed=True)
    summary = tmp_path / "summary.json"
    _write_summary(summary, _summary_row(assets))

    report = run_earth_climate_moisture_response_gate(
        EarthClimateMoistureResponseGateConfig(
            summary, tmp_path / "out", render_contact_sheets=False)
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "ocean_preservation"
        and row["metric"] == "response_ocean_abs_dev_p99"
        for row in report["failures"]
    )


def test_moisture_response_gate_flags_waterworld_strong_response(tmp_path):
    assets = tmp_path / "waterworld_seed7"
    assets.mkdir()
    _write_arrays(
        assets / "terminal_climate_arrays.npz",
        waterworld=True,
        strong_waterworld=True,
    )
    summary = tmp_path / "summary.json"
    _write_summary(summary, _summary_row(assets, preset="waterworld"))

    report = run_earth_climate_moisture_response_gate(
        EarthClimateMoistureResponseGateConfig(
            summary, tmp_path / "out", render_contact_sheets=False)
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "waterworld_false_positive"
        and row["metric"] == "response_land_p95"
        for row in report["failures"]
    )


def test_moisture_response_gate_cli_can_fail(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    _write_arrays(assets / "terminal_climate_arrays.npz", missing_response=True)
    summary = tmp_path / "summary.json"
    _write_summary(summary, _summary_row(assets))

    args = SimpleNamespace(
        terminal_summary=summary,
        out=tmp_path / "out",
        fail_on_moisture_response=True,
        no_contact_sheet=True,
    )
    with pytest.raises(SystemExit) as exc_info:
        cmd_earth_climate_moisture_response_gate(args)

    assert exc_info.value.code == 2
    assert (tmp_path / "out" / "earth_climate_moisture_response_checks.csv").exists()
