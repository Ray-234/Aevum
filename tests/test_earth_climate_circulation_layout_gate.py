import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_circulation_layout_gate
from aevum.diagnostics.earth_climate_circulation_layout_gate import (
    EarthClimateCirculationLayoutGateConfig,
    _coast_orientation,
    run_earth_climate_circulation_layout_gate,
)


def _basis(lat, lon):
    lat_r = np.radians(np.asarray(lat, dtype=float))
    lon_r = np.radians(np.asarray(lon, dtype=float))
    xyz = np.column_stack([
        np.cos(lat_r) * np.cos(lon_r),
        np.cos(lat_r) * np.sin(lon_r),
        np.sin(lat_r),
    ])
    z_axis = np.array([0.0, 0.0, 1.0])
    east = np.cross(z_axis, xyz)
    east /= np.maximum(np.linalg.norm(east, axis=1, keepdims=True), 1.0e-9)
    north = np.cross(xyz, east)
    return east, north


def _write_case(tmp_path, *, wind_scale=8.0):
    n = 36
    lat = np.linspace(-50.0, 50.0, n)
    lon = np.linspace(-170.0, 170.0, n)
    east, north = _basis(lat, lon)
    area = np.ones(n, dtype=float)
    land = lon < 0.0

    earth_uv = np.zeros((4, n, 2), dtype=float)
    earth_uv[:, :, 0] = 8.0
    earth_current = np.zeros((n, 2), dtype=float)
    earth_current[~land, 0] = 0.20
    earth = tmp_path / "earth_reference.npz"
    np.savez(
        earth,
        lat=lat,
        lon=lon,
        cell_area=area,
        earth__elevation_m=np.where(land, 100.0, -1000.0),
        earth__seasonal_wind_u10_v10=earth_uv,
        earth__annual_surface_current_u_v=earth_current,
    )

    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    to_ocean, coast_strength = _coast_orientation(lat, lon, land)
    seasonal_wind = np.repeat((wind_scale * east)[None, :, :], 4, axis=0)
    coast = coast_strength > 0.0
    seasonal_wind[0, coast] = wind_scale * to_ocean[coast]
    seasonal_wind[2, coast] = -wind_scale * to_ocean[coast]
    background = np.repeat((8.0 * east)[None, :, :], 4, axis=0)
    thermal = np.zeros_like(seasonal_wind)
    thermal[:, land] = 0.7 * north[None, land, :]
    orographic = np.zeros_like(seasonal_wind)
    currents = np.zeros((n, 3), dtype=float)
    currents[~land] = 0.20 * east[~land]
    np.savez(
        assets / "terminal_climate_arrays.npz",
        lat=lat,
        lon=lon,
        cell_area=area,
        terrain__elevation_m=np.where(land, 100.0, -1000.0),
        sea_level_m=np.asarray([0.0]),
        atmosphere__seasonal_wind=seasonal_wind,
        atmosphere__background_seasonal_wind=background,
        atmosphere__thermal_wind_anomaly=thermal,
        atmosphere__orographic_wind_anomaly=orographic,
        ocean__currents=currents,
    )
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "summaries": [{
            "label": "earthlike_seed1",
            "preset": "earthlike_mobile_lid",
            "seed": 1,
            "assets_dir": str(assets),
            "arrays": str(assets / "terminal_climate_arrays.npz"),
        }],
    }))
    return earth, summary


def test_circulation_layout_gate_passes_reasonable_wind(tmp_path):
    earth, summary = _write_case(tmp_path, wind_scale=8.0)

    report = run_earth_climate_circulation_layout_gate(
        EarthClimateCirculationLayoutGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["failure_count"] == 0
    assert (tmp_path / "out" / "earth_climate_circulation_layout_metrics.csv").exists()


def test_circulation_layout_gate_flags_overstrong_template_wind(tmp_path):
    earth, summary = _write_case(tmp_path, wind_scale=20.0)

    report = run_earth_climate_circulation_layout_gate(
        EarthClimateCirculationLayoutGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "wind_speed"
        and row["metric"] == "wind_speed_p90_ratio_to_earth"
        for row in report["failures"]
    )


def test_circulation_layout_gate_cli_can_fail(tmp_path):
    earth, summary = _write_case(tmp_path, wind_scale=20.0)
    args = SimpleNamespace(
        earth_reference=str(earth),
        terminal_summary=str(summary),
        out=str(tmp_path / "out"),
        fail_on_circulation_layout=True,
    )

    with pytest.raises(SystemExit):
        cmd_earth_climate_circulation_layout_gate(args)
