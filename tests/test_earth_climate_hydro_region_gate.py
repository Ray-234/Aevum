import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_hydro_region_gate
from aevum.diagnostics.earth_climate_hydro_region_gate import (
    SCHEMA,
    EarthClimateHydroRegionGateConfig,
    run_earth_climate_hydro_region_gate,
)


KINDS = (
    "monsoon_rainfall_corridor",
    "storm_track_rainfall_corridor",
    "rain_shadow_region",
    "wet_regional_precipitation_response",
    "dry_regional_precipitation_response",
)
SEASONS = ("DJF", "MAM", "JJA", "SON")


def _region(kind, season, idx, area=0.01):
    season_lat = {"DJF": -28.0, "MAM": -2.0, "JJA": 32.0, "SON": 8.0}
    return {
        "id": f"{kind}:{season.lower()}:{idx}",
        "type": "hydroclimate_region",
        "kind": kind,
        "season": season,
        "season_index": SEASONS.index(season),
        "cell_count": 4 + idx % 5,
        "area_fraction": area,
        "centroid_lat": season_lat[season],
        "centroid_lon": float(idx),
        "mean_intensity": 0.08,
        "max_intensity": 0.2,
    }


def _write_regions(path, regions):
    path.write_text(json.dumps({
        "schema": "aevum.hydroclimate_regions.v1",
        "region_count": len(regions),
        "regions": regions,
    }))


def _write_summary(path, row):
    path.write_text(json.dumps({
        "schema": "aevum.terminal_climate_replay.v1",
        "summaries": [row],
    }))


def _write_arrays(path, *, waterworld=False, sparse_monsoon=False):
    n = 80
    idx = np.arange(n, dtype=float) + 0.5
    z = 1.0 - 2.0 * idx / n
    lat = np.degrees(np.arcsin(np.clip(z, -1.0, 1.0)))
    lon = ((np.degrees(np.pi * (1.0 + 5.0 ** 0.5) * idx) + 180.0) % 360.0) - 180.0
    n = lat.size
    area = np.ones(n, dtype=float)
    land = np.ones(n, dtype=bool)
    elevation = np.ones(n, dtype=float)
    monsoon = np.zeros((4, n), dtype=float)
    storm = np.zeros((4, n), dtype=float)
    shadow = np.zeros((4, n), dtype=float)
    response = np.ones((4, n), dtype=float)
    seasonal_precip = np.full((4, n), 450.0, dtype=float)
    if sparse_monsoon:
        monsoon[0, 1] = 1.0
        monsoon[1, 5] = 0.7
        monsoon[2, 10] = 1.0
        monsoon[3, 7] = 0.7
    elif waterworld:
        monsoon[:, :] = 0.02
    else:
        monsoon[0, lat < -10] = 1.0
        monsoon[1, np.abs(lat) < 24] = 0.7
        monsoon[2, lat > 10] = 1.0
        monsoon[3, (lat > -5) & (lat < 25)] = 0.7
    storm_band = np.exp(-((np.abs(lat) - 40.0) / 9.0) ** 2)
    storm[:, :] = storm_band[None, :]
    base_shadow = np.clip((lat + 75.0) / 150.0, 0.05, 1.0)
    shadow[:, :] = base_shadow[None, :]
    wet_bump = np.exp(-((np.abs(lat) - 28.0) / 18.0) ** 2)
    response[:, :] = 1.0 + 0.22 * wet_bump[None, :] - 0.32 * shadow
    seasonal_precip *= np.clip(response, 0.35, 1.8)
    np.savez_compressed(
        path,
        lat=lat,
        lon=lon,
        cell_area=area,
        sea_level_m=np.array([0.0], dtype=float),
        terrain__elevation_m=elevation,
        climate__seasonal_precipitation=seasonal_precip,
        climate__monsoon_rainfall_corridor=monsoon,
        climate__storm_track_rainfall_corridor=storm,
        climate__rain_shadow_index=shadow,
        climate__regional_precipitation_response=response,
        climate__coast_distance=np.full(n, 0.05, dtype=float),
    )


def _valid_earthlike_regions():
    regions = []
    idx = 0
    for kind in KINDS:
        for season in SEASONS:
            for _ in range(5):
                regions.append(_region(kind, season, idx, area=0.009 + 0.0001 * idx))
                idx += 1
    return regions


def test_hydro_region_gate_passes_complete_earthlike_objects(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    _write_regions(assets / "hydroclimate_regions.json", _valid_earthlike_regions())
    _write_arrays(assets / "terminal_climate_arrays.npz")
    summary = tmp_path / "summary.json"
    _write_summary(summary, {
        "preset": "earthlike_mobile_lid",
        "seed": 1,
        "assets_dir": str(assets),
    })

    report = run_earth_climate_hydro_region_gate(
        EarthClimateHydroRegionGateConfig(summary, tmp_path / "out")
    )

    assert report["schema"] == SCHEMA
    assert report["verdict"] == "pass"
    assert report["failure_count"] == 0
    assert report["contact_sheet_count"] == 1
    assert (tmp_path / "out" / "earth_climate_hydro_region_metrics.csv").exists()
    assert (tmp_path / "out" / "earth_climate_hydro_region_checks.csv").exists()
    assert (
        tmp_path / "out" / "contact_sheets"
        / "earthlike_seed1_hydroclimate_region_contact_sheet.png"
    ).exists()
    assert (
        tmp_path / "out" / "earth_climate_hydro_region_contact_sheets.json"
    ).exists()


def test_hydro_region_gate_flags_missing_earthlike_archive(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    summary = tmp_path / "summary.json"
    _write_summary(summary, {
        "preset": "earthlike_mobile_lid",
        "seed": 1,
        "assets_dir": str(assets),
    })

    report = run_earth_climate_hydro_region_gate(
        EarthClimateHydroRegionGateConfig(summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "object_archive"
        and row["metric"] == "object_json_found"
        for row in report["failures"]
    )


def test_hydro_region_gate_flags_waterworld_broad_monsoon_object(tmp_path):
    assets = tmp_path / "waterworld_seed7"
    assets.mkdir()
    regions = [_region("monsoon_rainfall_corridor", "JJA", 0, area=0.08)]
    regions.extend(_region("storm_track_rainfall_corridor", "DJF", i, area=0.004)
                   for i in range(6))
    _write_regions(assets / "hydroclimate_regions.json", regions)
    _write_arrays(assets / "terminal_climate_arrays.npz", waterworld=True)
    summary = tmp_path / "summary.json"
    _write_summary(summary, {
        "preset": "waterworld",
        "seed": 7,
        "assets_dir": str(assets),
    })

    report = run_earth_climate_hydro_region_gate(
        EarthClimateHydroRegionGateConfig(summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "waterworld_false_positive"
        and row["metric"] == "largest_monsoon_area_fraction"
        for row in report["failures"]
    )


def test_hydro_region_gate_flags_sparse_earthlike_monsoon_map(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    _write_regions(assets / "hydroclimate_regions.json", _valid_earthlike_regions())
    _write_arrays(assets / "terminal_climate_arrays.npz", sparse_monsoon=True)
    summary = tmp_path / "summary.json"
    _write_summary(summary, {
        "preset": "earthlike_mobile_lid",
        "seed": 1,
        "assets_dir": str(assets),
    })

    report = run_earth_climate_hydro_region_gate(
        EarthClimateHydroRegionGateConfig(summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "monsoon_map_readability"
        and row["metric"] == "monsoon_map_active_land_fraction_p50"
        for row in report["failures"]
    )


def test_hydro_region_gate_cli_can_fail(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    summary = tmp_path / "summary.json"
    _write_summary(summary, {
        "preset": "earthlike_mobile_lid",
        "seed": 1,
        "assets_dir": str(assets),
    })

    args = SimpleNamespace(
        terminal_summary=summary,
        out=tmp_path / "out",
        fail_on_hydro_regions=True,
    )
    with pytest.raises(SystemExit) as exc_info:
        cmd_earth_climate_hydro_region_gate(args)

    assert exc_info.value.code == 2
    assert (tmp_path / "out" / "earth_climate_hydro_region_checks.csv").exists()
