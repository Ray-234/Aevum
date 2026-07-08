import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_coupling_convergence_gate
from aevum.diagnostics.earth_climate_coupling_convergence_gate import (
    EarthClimateCouplingConvergenceGateConfig,
    run_earth_climate_coupling_convergence_gate,
)


def _write_coupling_case(tmp_path, *, runaway_feedback: bool = False):
    lat_vals = np.linspace(-55.0, 55.0, 10)
    lon_vals = np.linspace(-170.0, 170.0, 18)
    lon_grid, lat_grid = np.meshgrid(lon_vals, lat_vals)
    lat = lat_grid.ravel()
    lon = lon_grid.ravel()
    n = lat.size
    land = np.abs(lon) < 80.0
    elev = np.where(land, 220.0, -1800.0)

    phases = np.array([-1.0, -0.25, 1.0, 0.25])
    base_pattern = 0.55 * np.sin(np.radians(lon)) + 0.25 * np.sin(np.radians(lat))
    base_pressure = phases[:, None] * base_pattern[None, :]
    feedback = np.repeat((0.018 * np.cos(np.radians(lon)))[None, :], 4, axis=0)
    feedback[:, ~land] = 0.0
    if runaway_feedback:
        feedback[:, land & (lon > 0.0)] = 0.12
    final_pressure = np.clip(base_pressure + feedback, -1.8, 1.8)

    wind = np.zeros((4, n, 3), dtype=float)
    wind[:, :, 0] = 6.0
    wind[:, :, 1] = 1.0 * phases[:, None]
    wind_anom = np.zeros_like(wind)
    wind_anom[:, :, 1] = np.where(land, 0.045, 0.0)
    if runaway_feedback:
        wind_anom[:, land, 1] = 0.55

    seasonal_precip = np.repeat((350.0 + 380.0 * land)[None, :], 4, axis=0)
    seasonal_precip += 60.0 * phases[:, None] * np.cos(np.radians(lat))[None, :]
    seasonal_precip = np.maximum(seasonal_precip, 0.0)
    annual_precip = seasonal_precip.mean(axis=0)
    hydro_residual = np.where(land, 0.016, 0.0)
    iteration_delta = np.where(land, 0.004, 0.0)
    if runaway_feedback:
        iteration_delta[land & (lon > 0.0)] = 0.08
    ocean_residual = np.where(~land, 0.006, 0.0)
    ocean_evap_feedback = np.where(~land, 0.22 * np.sin(np.radians(lon)), 0.0)
    ocean_evap_feedback[~land] -= float(np.mean(ocean_evap_feedback[~land]))
    wind_stress_response = np.zeros((n, 3), dtype=float)
    wind_stress_response[~land, 0] = 0.08

    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    np.savez(
        assets / "terminal_climate_arrays.npz",
        lat=lat,
        lon=lon,
        cell_area=np.ones(n, dtype=float),
        terrain__elevation_m=elev,
        sea_level_m=np.asarray([0.0]),
        climate__precipitation=annual_precip,
        climate__seasonal_precipitation=seasonal_precip,
        climate__coupling_residual=ocean_residual,
        climate__ocean_evaporation_feedback=ocean_evap_feedback,
        climate__hydro_coupling_residual=hydro_residual,
        climate__hydro_feedback_iteration_delta=iteration_delta,
        atmosphere__land_sea_pressure_proxy=base_pressure,
        atmosphere__seasonal_pressure_proxy=final_pressure,
        atmosphere__precipitation_pressure_feedback=feedback,
        atmosphere__seasonal_wind=wind,
        atmosphere__hydro_coupled_wind_anomaly=wind_anom,
        ocean__wind_stress_current_response=wind_stress_response,
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
    return summary


def test_coupling_convergence_gate_passes_bounded_feedback(tmp_path):
    summary = _write_coupling_case(tmp_path)

    report = run_earth_climate_coupling_convergence_gate(
        EarthClimateCouplingConvergenceGateConfig(summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["failure_count"] == 0
    assert (
        tmp_path / "out" / "earth_climate_coupling_convergence_metrics.csv"
    ).exists()


def test_coupling_convergence_gate_flags_runaway_feedback(tmp_path):
    summary = _write_coupling_case(tmp_path, runaway_feedback=True)

    report = run_earth_climate_coupling_convergence_gate(
        EarthClimateCouplingConvergenceGateConfig(summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] in {"pressure_feedback", "wind_feedback"}
        for row in report["failures"]
    )


def test_coupling_convergence_gate_cli_can_fail(tmp_path):
    summary = _write_coupling_case(tmp_path, runaway_feedback=True)
    args = SimpleNamespace(
        terminal_summary=str(summary),
        out=str(tmp_path / "out"),
        fail_on_coupling_convergence=True,
    )

    with pytest.raises(SystemExit):
        cmd_earth_climate_coupling_convergence_gate(args)
