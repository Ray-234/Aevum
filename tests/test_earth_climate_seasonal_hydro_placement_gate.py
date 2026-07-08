import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_seasonal_hydro_placement_gate
from aevum.diagnostics.earth_climate_seasonal_hydro_placement_gate import (
    EarthClimateSeasonalHydroPlacementGateConfig,
    run_earth_climate_seasonal_hydro_placement_gate,
)


def _write_hydro_placement_case(tmp_path, *, unsupported_wet: bool = False):
    lat_vals = np.linspace(-55.0, 55.0, 12)
    lon_vals = np.linspace(-170.0, 170.0, 20)
    lon_grid, lat_grid = np.meshgrid(lon_vals, lat_vals)
    lat = lat_grid.ravel()
    lon = lon_grid.ravel()
    n = lat.size
    land = np.abs(lon) < 95.0
    elev = np.where(land, 180.0, -1800.0)

    tropical = np.exp(-(lat / 16.0) ** 2)
    midlat = np.exp(-((np.abs(lat) - 42.0) / 10.0) ** 2)
    east_wet = np.clip((lon + 95.0) / 190.0, 0.0, 1.0)
    west_shadow = land & (lon < -35.0)
    phases = np.array([-1.0, -0.25, 1.0, 0.25])

    moisture = np.repeat((0.22 + 0.66 * east_wet)[None, :], 4, axis=0)
    moisture[:, ~land] = 0.90
    itcz = np.stack([
        0.95 * np.exp(-((lat - 11.0 * phase) / 17.0) ** 2)
        for phase in phases
    ])
    monsoon = np.stack([
        np.where(land, np.clip(1.05 * east_wet * (0.55 + 0.45 * phase), 0.0, 1.1), 0.0)
        for phase in phases
    ])
    storm = np.stack([
        0.92 * midlat * (1.0 + 0.25 * np.where(lat >= 0.0, -phase, phase))
        for phase in phases
    ])
    rain_shadow = np.repeat(np.where(west_shadow, 0.95, 0.05)[None, :], 4, axis=0)
    response = 0.78 + 0.42 * east_wet[None, :] + 0.10 * tropical[None, :]
    response = np.repeat(response, 4, axis=0)
    flow_response = 0.82 + 0.34 * east_wet[None, :] + 0.08 * midlat[None, :]
    flow_response = np.repeat(flow_response, 4, axis=0)

    support = (
        np.clip(moisture, 0.0, 1.0)
        + 0.75 * np.clip(monsoon, 0.0, 1.2)
        + 0.55 * np.clip(storm, 0.0, 1.2)
        + 0.35 * np.clip(itcz, 0.0, 1.2)
        + 0.55 * np.clip(response - 1.0, 0.0, 1.0)
        + 0.45 * np.clip(flow_response - 1.0, 0.0, 1.0)
    )
    seasonal_precip = np.where(
        land[None, :],
        90.0 + 430.0 * support + 240.0 * itcz - 260.0 * rain_shadow,
        850.0 + 240.0 * itcz,
    )
    seasonal_precip = np.maximum(seasonal_precip, 8.0)
    if unsupported_wet:
        unsupported = land & (lon < -45.0) & (np.abs(lat) < 18.0)
        seasonal_precip[:, unsupported] = 1800.0
    annual_precip = seasonal_precip.mean(axis=0)

    earth = tmp_path / "earth_reference.npz"
    np.savez(
        earth,
        lat=lat,
        lon=lon,
        earth__land_mask=land,
        earth__seasonal_precip_mm_yr_equiv=seasonal_precip,
        earth__annual_precip_mm=annual_precip,
    )

    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    np.savez(
        assets / "terminal_climate_arrays.npz",
        lat=lat,
        lon=lon,
        terrain__elevation_m=elev,
        sea_level_m=np.asarray([0.0]),
        climate__seasonal_precipitation=seasonal_precip,
        climate__precipitation=annual_precip,
        climate__monsoon_rainfall_corridor=monsoon,
        climate__storm_track_rainfall_corridor=storm,
        climate__rain_shadow_index=rain_shadow,
        climate__regional_precipitation_response=response,
        climate__moisture_flow_precipitation_response=flow_response,
        atmosphere__moisture_access=moisture,
        atmosphere__itcz_intensity=itcz,
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


def test_seasonal_hydro_placement_gate_passes_process_supported_rain(tmp_path):
    earth, summary = _write_hydro_placement_case(tmp_path)

    report = run_earth_climate_seasonal_hydro_placement_gate(
        EarthClimateSeasonalHydroPlacementGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["failure_count"] == 0
    assert (
        tmp_path / "out" / "earth_climate_seasonal_hydro_placement_metrics.csv"
    ).exists()


def test_seasonal_hydro_placement_gate_flags_unsupported_wet_patch(tmp_path):
    earth, summary = _write_hydro_placement_case(tmp_path, unsupported_wet=True)

    report = run_earth_climate_seasonal_hydro_placement_gate(
        EarthClimateSeasonalHydroPlacementGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(row["group"] == "wet_support" for row in report["failures"])


def test_seasonal_hydro_placement_gate_cli_can_fail(tmp_path):
    earth, summary = _write_hydro_placement_case(tmp_path, unsupported_wet=True)
    args = SimpleNamespace(
        earth_reference=str(earth),
        terminal_summary=str(summary),
        out=str(tmp_path / "out"),
        fail_on_seasonal_hydro_placement=True,
    )

    with pytest.raises(SystemExit):
        cmd_earth_climate_seasonal_hydro_placement_gate(args)
