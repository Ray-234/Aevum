import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_ocean_spatial_gate
from aevum.diagnostics.earth_climate_ocean_spatial_gate import (
    EarthClimateOceanSpatialGateConfig,
    run_earth_climate_ocean_spatial_gate,
)


def _write_ocean_case(tmp_path, *, far_ocean_swift: bool = False):
    lat_vals = np.linspace(-60.0, 60.0, 13)
    lon_vals = np.linspace(-177.5, 177.5, 72)
    lon_grid, lat_grid = np.meshgrid(lon_vals, lat_vals)
    lat = lat_grid.ravel()
    lon = lon_grid.ravel()
    n = lat.size
    area = np.ones(n, dtype=float)
    land = np.abs(lon) < 30.0
    ocean = ~land
    coastal_ocean = ocean & (np.abs(lon) >= 30.0) & (np.abs(lon) < 65.0)
    far_ocean = ocean & (np.abs(lon) > 125.0)

    earth_speed = np.where(ocean, 0.08, 0.0)
    earth_speed[coastal_ocean] = 0.24
    earth_uv = np.zeros((n, 2), dtype=float)
    earth_uv[:, 0] = earth_speed
    sst = 24.0 - 0.35 * np.abs(lat)
    sst[coastal_ocean & (lat > 0.0)] += 1.4
    sst[coastal_ocean & (lat < 0.0)] -= 1.4
    earth = tmp_path / "earth_reference.npz"
    np.savez(
        earth,
        lat=lat,
        lon=lon,
        cell_area=area,
        earth__elevation_m=np.where(land, 120.0, -1200.0),
        earth__annual_surface_current_u_v=earth_uv,
        earth__annual_surface_current_speed_m_s=earth_speed,
        earth__annual_sst_C=sst,
    )

    generated_speed = np.where(ocean, 0.08, 0.0)
    if far_ocean_swift:
        generated_speed[far_ocean] = 0.24
    else:
        generated_speed[coastal_ocean] = 0.24
    currents = np.zeros((n, 3), dtype=float)
    currents[:, 0] = generated_speed
    heat_transport = np.where(
        ocean,
        generated_speed - generated_speed[ocean].mean(),
        0.0,
    )
    heat_scale = max(float(np.percentile(np.abs(heat_transport[ocean]), 95)), 1.0e-9)
    heat_transport = np.where(ocean, 1.25 * heat_transport / heat_scale, 0.0)
    seasonal_sst = np.repeat((sst + 273.15)[None, :], 4, axis=0)
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    np.savez(
        assets / "terminal_climate_arrays.npz",
        lat=lat,
        lon=lon,
        cell_area=area,
        terrain__elevation_m=np.where(land, 120.0, -1200.0),
        sea_level_m=np.asarray([0.0]),
        ocean__currents=currents,
        ocean__current_heat_transport=heat_transport,
        climate__seasonal_sst=seasonal_sst,
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


def test_ocean_spatial_gate_passes_boundary_current_structure(tmp_path):
    earth, summary = _write_ocean_case(tmp_path)

    report = run_earth_climate_ocean_spatial_gate(
        EarthClimateOceanSpatialGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["failure_count"] == 0
    row = report["generated"][0]
    assert row["sst_zonal_residual_abs_p95_ratio_to_earth"] >= 0.58
    assert row["current_heat_transport_abs_p95_C"] >= 0.95
    assert (tmp_path / "out" / "earth_climate_ocean_spatial_metrics.csv").exists()


def test_ocean_spatial_gate_flags_far_ocean_swift_band(tmp_path):
    earth, summary = _write_ocean_case(tmp_path, far_ocean_swift=True)

    report = run_earth_climate_ocean_spatial_gate(
        EarthClimateOceanSpatialGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "boundary_current_structure"
        and row["metric"] == "swift_current_far_ocean_share"
        for row in report["failures"]
    )


def test_ocean_spatial_gate_cli_can_fail(tmp_path):
    earth, summary = _write_ocean_case(tmp_path, far_ocean_swift=True)
    args = SimpleNamespace(
        earth_reference=str(earth),
        terminal_summary=str(summary),
        out=str(tmp_path / "out"),
        fail_on_ocean_spatial=True,
    )

    with pytest.raises(SystemExit):
        cmd_earth_climate_ocean_spatial_gate(args)
