import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_windward_leeward_gate
from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.diagnostics.earth_climate_windward_leeward_gate import (
    SCHEMA,
    EarthClimateWindwardLeewardGateConfig,
    run_earth_climate_windward_leeward_gate,
)


def _tangent_wind_toward_x(grid):
    x_axis = np.array([1.0, 0.0, 0.0])
    wind = x_axis[None, :] - (grid.xyz @ x_axis)[:, None] * grid.xyz
    norm = np.linalg.norm(wind, axis=1, keepdims=True)
    return np.where(norm > 1.0e-9, wind / np.maximum(norm, 1.0e-9), 0.0)


def _east_north_basis(grid):
    z_axis = np.array([0.0, 0.0, 1.0])
    east = np.cross(z_axis, grid.xyz)
    norm = np.linalg.norm(east, axis=1, keepdims=True)
    east = np.where(norm > 1.0e-9, east / np.maximum(norm, 1.0e-9),
                    np.array([1.0, 0.0, 0.0]))
    north = np.cross(grid.xyz, east)
    return east, north


def _synthetic_fields():
    grid = SphereGrid.fibonacci(240, CONSTANTS.EARTH_RADIUS)
    wind = _tangent_wind_toward_x(grid)
    ridge = 1.0 - np.abs(grid.xyz[:, 0])
    elevation = 450.0 + 2400.0 * np.clip(ridge, 0.0, 1.0)
    windward = grid.xyz[:, 0] < -0.18
    leeward = grid.xyz[:, 0] > 0.18
    earth_precip = np.where(windward, 1050.0, np.where(leeward, 320.0, 620.0))
    generated_precip = np.full(grid.n, 560.0, dtype=float)
    seasonal_wind = np.repeat(wind[None, :, :], 4, axis=0)
    return grid, elevation, earth_precip, generated_precip, seasonal_wind


def _write_earth(path):
    grid, elevation, earth_precip, _, seasonal_wind = _synthetic_fields()
    east, north = _east_north_basis(grid)
    seasonal_uv = np.empty((4, grid.n, 2), dtype=float)
    seasonal_uv[:, :, 0] = np.einsum("snk,nk->sn", seasonal_wind, east)
    seasonal_uv[:, :, 1] = np.einsum("snk,nk->sn", seasonal_wind, north)
    np.savez(
        path,
        lat=grid.lat,
        lon=grid.lon,
        cell_area=grid.cell_area,
        earth__land_mask=np.ones(grid.n, dtype=bool),
        earth__elevation_m=elevation,
        earth__annual_precip_mm=earth_precip,
        earth__seasonal_precip_mm_yr_equiv=np.repeat(earth_precip[None, :], 4, axis=0),
        earth__seasonal_wind_u10_v10=seasonal_uv,
    )


def _write_generated(path, precip):
    grid, elevation, _, _, seasonal_wind = _synthetic_fields()
    precip = np.asarray(precip, dtype=float)
    np.savez(
        path,
        lat=grid.lat,
        lon=grid.lon,
        cell_area=grid.cell_area,
        sea_level_m=np.array([0.0], dtype=float),
        terrain__elevation_m=elevation,
        climate__precipitation=precip,
        climate__seasonal_precipitation=np.repeat(precip[None, :], 4, axis=0),
        atmosphere__seasonal_wind=seasonal_wind,
        terrain__barrier_index=np.ones(grid.n, dtype=float),
        climate__orographic_precipitation=np.zeros(grid.n, dtype=float),
    )


def _summary(path, preset="earthlike_mobile_lid"):
    return {
        "summaries": [{
            "preset": preset,
            "seed": 1,
            "arrays": str(path),
            "assets_dir": str(path.parent / "earthlike_seed1"),
        }],
    }


def test_windward_leeward_gate_flags_flat_orographic_contrast(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _, _, _, generated_precip, _ = _synthetic_fields()
    _write_generated(generated, generated_precip)
    summary.write_text(json.dumps(_summary(generated)))

    report = run_earth_climate_windward_leeward_gate(
        EarthClimateWindwardLeewardGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["schema"] == SCHEMA
    assert report["verdict"] == "fail"
    failed = {row["metric"] for row in report["failures"]}
    assert "annual_windward_leeward_precip_ratio" in failed
    assert "seasonal_windward_leeward_precip_ratio_median" in failed
    assert (tmp_path / "out" / "earth_climate_windward_leeward_metrics.csv").exists()
    assert (tmp_path / "out" / "earth_climate_windward_leeward_checks.csv").exists()


def test_windward_leeward_gate_skips_diagnostic_only_world(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _, _, _, generated_precip, _ = _synthetic_fields()
    _write_generated(generated, generated_precip)
    summary.write_text(json.dumps(_summary(generated, preset="waterworld")))

    report = run_earth_climate_windward_leeward_gate(
        EarthClimateWindwardLeewardGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["checks"] == []
    assert report["generated_metrics"][0]["mode"] == "diagnostic_only"


def test_windward_leeward_gate_cli_can_fail(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _, _, _, generated_precip, _ = _synthetic_fields()
    _write_generated(generated, generated_precip)
    summary.write_text(json.dumps(_summary(generated)))

    args = SimpleNamespace(
        earth_reference=earth,
        terminal_summary=summary,
        out=tmp_path / "out",
        fail_on_windward_leeward=True,
    )
    with pytest.raises(SystemExit) as exc_info:
        cmd_earth_climate_windward_leeward_gate(args)

    assert exc_info.value.code == 2
    assert (
        tmp_path / "out" / "earth_climate_windward_leeward_checks.csv"
    ).exists()
