import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_coupled_consistency_gate
from aevum.diagnostics.earth_climate_coupled_consistency_gate import (
    EarthClimateCoupledConsistencyGateConfig,
    run_earth_climate_coupled_consistency_gate,
)


def _basis(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lat_r = np.radians(lat)
    lon_r = np.radians(lon)
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


def _write_coupled_case(tmp_path, *, dry_monsoon: bool = False):
    lat_vals = np.linspace(-55.0, 55.0, 10)
    lon_vals = np.linspace(-165.0, 165.0, 18)
    lon_grid, lat_grid = np.meshgrid(lon_vals, lat_vals)
    lat = lat_grid.ravel()
    lon = lon_grid.ravel()
    n = lat.size
    area = np.ones(n, dtype=float)
    land = np.abs(lon) < 70.0
    ocean = ~land
    east, _ = _basis(lat, lon)

    lon_norm = np.clip(lon / 165.0, -1.0, 1.0)
    land_moisture_pattern = np.clip((lon_norm + 1.0) / 2.0, 0.0, 1.0)
    source_pattern = np.where(ocean, np.clip((lon_norm + 1.0) / 2.0, 0.0, 1.0), 0.0)
    phases = np.array([-1.0, -0.4, 1.0, 0.4])

    earth_pressure = -phases[:, None] * lon_norm[None, :]
    earth_uv = np.zeros((4, n, 2), dtype=float)
    earth_uv[:, :, 0] = np.sign(phases)[:, None] * 7.0
    earth_temp = 16.0 + phases[:, None] * lon_norm[None, :] * 5.0
    earth = tmp_path / "earth_reference.npz"
    np.savez(
        earth,
        lat=lat,
        lon=lon,
        cell_area=area,
        earth__land_mask=land,
        earth__elevation_m=np.where(land, 100.0, -1500.0),
        earth__seasonal_slp_anomaly_hPa=earth_pressure,
        earth__seasonal_wind_u10_v10=earth_uv,
        earth__seasonal_temperature_C=earth_temp,
    )

    pressure = earth_pressure.copy()
    wind = np.sign(phases)[:, None, None] * east[None, :, :] * 6.0
    seasonal_temp = earth_temp + 273.15
    seasonal_sst = np.repeat((276.15 + 20.0 * source_pattern)[None, :], 4, axis=0)
    source = np.repeat(source_pattern[None, :], 4, axis=0)
    moisture = np.repeat(
        np.where(land, 0.18 + 0.75 * land_moisture_pattern, source_pattern)[None, :],
        4,
        axis=0,
    )
    wet_monsoon = np.clip((-pressure) * moisture, 0.0, 1.2)
    dry_core = np.where(land, 1.0 - land_moisture_pattern, 0.0)
    monsoon = (
        np.clip((-pressure) * dry_core[None, :], 0.0, 1.2)
        if dry_monsoon else wet_monsoon
    )
    storm = np.repeat((0.35 * np.exp(-((np.abs(lat) - 40.0) / 12.0) ** 2))[None, :],
                      4, axis=0)
    itcz = np.repeat((0.25 * np.exp(-(lat / 18.0) ** 2))[None, :], 4, axis=0)
    seasonal_precip = np.where(
        land[None, :],
        120.0 + 650.0 * moisture + 260.0 * np.clip(monsoon, 0.0, 1.0),
        900.0 * source + 100.0,
    )
    annual_precip = seasonal_precip.mean(axis=0)
    evaporation = np.where(ocean, 180.0 + 900.0 * source_pattern, 0.0)
    heat = np.where(ocean, 2.0 * (source_pattern - 0.5), 0.0)
    upwelling = np.where(ocean & (source_pattern < 0.35), 0.08, 0.0)
    ocean_heat_flux = np.where(ocean, heat - np.mean(heat[ocean]), 0.0)
    coupling_residual = np.where(ocean, 0.05, 0.0)

    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    np.savez(
        assets / "terminal_climate_arrays.npz",
        lat=lat,
        lon=lon,
        cell_area=area,
        terrain__elevation_m=np.where(land, 100.0, -1500.0),
        sea_level_m=np.asarray([0.0]),
        climate__seasonal_temperature=seasonal_temp,
        climate__seasonal_sst=seasonal_sst,
        climate__seasonal_precipitation=seasonal_precip,
        climate__precipitation=annual_precip,
        climate__evaporation=evaporation,
        climate__ocean_heat_flux=ocean_heat_flux,
        climate__coupling_residual=coupling_residual,
        atmosphere__seasonal_pressure_proxy=pressure,
        atmosphere__seasonal_wind=wind,
        atmosphere__source_ocean_warmth=source,
        atmosphere__moisture_access=moisture,
        atmosphere__monsoon_potential=monsoon,
        atmosphere__storm_track_intensity=storm,
        atmosphere__itcz_intensity=itcz,
        ocean__current_heat_transport=heat,
        ocean__upwelling=upwelling,
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


def test_coupled_consistency_gate_passes_coherent_fields(tmp_path):
    earth, summary = _write_coupled_case(tmp_path)

    report = run_earth_climate_coupled_consistency_gate(
        EarthClimateCoupledConsistencyGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["failure_count"] == 0
    assert (
        tmp_path / "out" / "earth_climate_coupled_consistency_metrics.csv"
    ).exists()


def test_coupled_consistency_gate_flags_dry_monsoon_potential(tmp_path):
    earth, summary = _write_coupled_case(tmp_path, dry_monsoon=True)

    report = run_earth_climate_coupled_consistency_gate(
        EarthClimateCoupledConsistencyGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "monsoon_support"
        and row["metric"] == "monsoon_top_moisture_ratio"
        for row in report["failures"]
    )


def test_coupled_consistency_gate_cli_can_fail(tmp_path):
    earth, summary = _write_coupled_case(tmp_path, dry_monsoon=True)
    args = SimpleNamespace(
        earth_reference=str(earth),
        terminal_summary=str(summary),
        out=str(tmp_path / "out"),
        fail_on_coupled_consistency=True,
    )

    with pytest.raises(SystemExit):
        cmd_earth_climate_coupled_consistency_gate(args)
