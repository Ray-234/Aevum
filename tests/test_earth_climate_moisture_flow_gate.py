import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_moisture_flow_gate
from aevum.diagnostics.earth_climate_moisture_flow_gate import (
    SCHEMA,
    EarthClimateMoistureFlowGateConfig,
    run_earth_climate_moisture_flow_gate,
)


KINDS = (
    "monsoon_moisture_flow_network",
    "storm_track_moisture_flow_network",
    "mixed_moisture_flow_network",
)
SEASONS = ("DJF", "MAM", "JJA", "SON")


def _network(kind, season, idx, area=0.012):
    season_lat = {"DJF": -24.0, "MAM": -4.0, "JJA": 28.0, "SON": 6.0}
    return {
        "id": f"{kind}:{season.lower()}:{idx}",
        "type": "moisture_flow_network",
        "kind": kind,
        "season": season,
        "season_index": SEASONS.index(season),
        "cell_count": 8 + idx % 7,
        "area_fraction": area,
        "centroid_lat": season_lat[season],
        "centroid_lon": float(idx),
        "mean_pathway": 0.72,
        "mean_moisture_access": 0.68,
        "mean_monsoon_corridor": 0.11 if "monsoon" in kind else 0.02,
        "mean_storm_track_corridor": 0.16 if "storm" in kind else 0.04,
        "mean_rain_shadow": 0.05,
        "mean_precipitation": 720.0,
    }


def _write_networks(path, networks):
    path.write_text(json.dumps({
        "schema": "aevum.moisture_flow_networks.v1",
        "network_count": len(networks),
        "networks": networks,
    }))


def _write_summary(path, row):
    path.write_text(json.dumps({
        "schema": "aevum.terminal_climate_replay.v1",
        "summaries": [row],
    }))


def _write_arrays(path, *, sparse_pathway=False):
    n = 120
    idx = np.arange(n, dtype=float) + 0.5
    z = 1.0 - 2.0 * idx / n
    lat = np.degrees(np.arcsin(np.clip(z, -1.0, 1.0)))
    lon = ((np.degrees(np.pi * (1.0 + 5.0 ** 0.5) * idx) + 180.0) % 360.0) - 180.0
    area = np.ones(n, dtype=float)
    land = (np.abs(lat) < 60.0) & (lon > -155.0)
    elevation = np.where(land, 120.0, -150.0)

    source = np.zeros((4, n), dtype=float)
    pathway = np.zeros((4, n), dtype=float)
    network_id = np.full((4, n), -1.0, dtype=float)
    monsoon = np.zeros((4, n), dtype=float)
    storm = np.zeros((4, n), dtype=float)
    shadow = np.zeros((4, n), dtype=float)
    precip = np.full((4, n), 500.0, dtype=float)
    centers = (-24.0, -4.0, 28.0, 6.0)
    for season, center in enumerate(centers):
        monsoon_band = np.exp(-((lat - center) / 20.0) ** 2)
        storm_band = np.exp(-((np.abs(lat) - 42.0) / 11.0) ** 2)
        support = np.maximum(0.95 * monsoon_band, 0.70 * storm_band)
        source[season, ~land] = 0.72 + 0.16 * np.cos(np.radians(lat[~land])) ** 2
        source[season, land] = 0.08 + 0.12 * support[land]
        if sparse_pathway:
            active = np.zeros(n, dtype=bool)
            active[season * 3 + 4] = True
            pathway[season, active] = 0.95
        else:
            pathway[season, land] = 0.10 + 0.92 * support[land]
        network_id[season, land & (pathway[season] > 0.42)] = 1.0 + season
        monsoon[season] = monsoon_band
        storm[season] = storm_band
        shadow[season] = np.clip((np.abs(lat) - 15.0) / 80.0, 0.0, 0.6)
        precip[season, land] *= 1.0 + 0.55 * support[land]

    np.savez_compressed(
        path,
        lat=lat,
        lon=lon,
        cell_area=area,
        sea_level_m=np.array([0.0], dtype=float),
        terrain__elevation_m=elevation,
        atmosphere__moisture_flow_source=source,
        atmosphere__moisture_flow_pathway=pathway,
        climate__moisture_flow_network_id=network_id,
        climate__seasonal_precipitation=precip,
        climate__monsoon_rainfall_corridor=monsoon,
        climate__storm_track_rainfall_corridor=storm,
        climate__rain_shadow_index=shadow,
    )


def _valid_earthlike_networks():
    networks = []
    idx = 0
    for season in SEASONS:
        for _ in range(2):
            networks.append(_network("monsoon_moisture_flow_network", season, idx))
            idx += 1
        for _ in range(6):
            networks.append(_network("storm_track_moisture_flow_network", season, idx))
            idx += 1
        networks.append(_network("mixed_moisture_flow_network", season, idx, area=0.010))
        idx += 1
    return networks


def test_moisture_flow_gate_passes_complete_earthlike_objects(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    _write_networks(assets / "moisture_flow_networks.json", _valid_earthlike_networks())
    _write_arrays(assets / "terminal_climate_arrays.npz")
    summary = tmp_path / "summary.json"
    _write_summary(summary, {
        "preset": "earthlike_mobile_lid",
        "seed": 1,
        "assets_dir": str(assets),
    })

    report = run_earth_climate_moisture_flow_gate(
        EarthClimateMoistureFlowGateConfig(summary, tmp_path / "out")
    )

    assert report["schema"] == SCHEMA
    assert report["verdict"] == "pass"
    assert report["failure_count"] == 0
    assert report["contact_sheet_count"] == 1
    assert (tmp_path / "out" / "earth_climate_moisture_flow_metrics.csv").exists()
    assert (tmp_path / "out" / "earth_climate_moisture_flow_checks.csv").exists()
    assert (
        tmp_path / "out" / "contact_sheets"
        / "earthlike_seed1_moisture_flow_network_contact_sheet.png"
    ).exists()
    assert (
        tmp_path / "out" / "earth_climate_moisture_flow_contact_sheets.json"
    ).exists()


def test_moisture_flow_gate_flags_missing_earthlike_archive(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    summary = tmp_path / "summary.json"
    _write_summary(summary, {
        "preset": "earthlike_mobile_lid",
        "seed": 1,
        "assets_dir": str(assets),
    })

    report = run_earth_climate_moisture_flow_gate(
        EarthClimateMoistureFlowGateConfig(summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "object_archive"
        and row["metric"] == "object_json_found"
        for row in report["failures"]
    )


def test_moisture_flow_gate_flags_waterworld_broad_monsoon_network(tmp_path):
    assets = tmp_path / "waterworld_seed7"
    assets.mkdir()
    networks = [_network("monsoon_moisture_flow_network", "JJA", 0, area=0.08)]
    networks.extend(_network("storm_track_moisture_flow_network", "DJF", i, area=0.004)
                    for i in range(6))
    _write_networks(assets / "moisture_flow_networks.json", networks)
    _write_arrays(assets / "terminal_climate_arrays.npz")
    summary = tmp_path / "summary.json"
    _write_summary(summary, {
        "preset": "waterworld",
        "seed": 7,
        "assets_dir": str(assets),
    })

    report = run_earth_climate_moisture_flow_gate(
        EarthClimateMoistureFlowGateConfig(summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "waterworld_false_positive"
        and row["metric"] == "largest_monsoon_flow_area_fraction"
        for row in report["failures"]
    )


def test_moisture_flow_gate_flags_sparse_earthlike_pathway_map(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    _write_networks(assets / "moisture_flow_networks.json", _valid_earthlike_networks())
    _write_arrays(assets / "terminal_climate_arrays.npz", sparse_pathway=True)
    summary = tmp_path / "summary.json"
    _write_summary(summary, {
        "preset": "earthlike_mobile_lid",
        "seed": 1,
        "assets_dir": str(assets),
    })

    report = run_earth_climate_moisture_flow_gate(
        EarthClimateMoistureFlowGateConfig(summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "pathway_map_readability"
        and row["metric"] == "pathway_map_active_land_fraction_p50"
        for row in report["failures"]
    )


def test_moisture_flow_gate_cli_can_fail(tmp_path):
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
        fail_on_moisture_flow=True,
        no_contact_sheet=True,
    )
    with pytest.raises(SystemExit) as exc_info:
        cmd_earth_climate_moisture_flow_gate(args)

    assert exc_info.value.code == 2
    assert (tmp_path / "out" / "earth_climate_moisture_flow_checks.csv").exists()
