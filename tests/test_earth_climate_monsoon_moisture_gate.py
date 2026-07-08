from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from aevum.cli import cmd_earth_climate_monsoon_moisture_gate
from aevum.diagnostics.earth_climate_monsoon_moisture_gate import (
    EarthClimateMonsoonMoistureGateConfig,
    run_earth_climate_monsoon_moisture_gate,
)


def _write_earth(path: Path) -> None:
    n = 8
    lat = np.array([-30, -20, -10, 10, 20, 30, 45, -45], dtype=float)
    seasonal_precip = np.full((4, n), 120.0)
    seasonal_precip[0, :3] = 1500.0
    seasonal_precip[2, 3:6] = 1500.0
    annual_precip = seasonal_precip.mean(axis=0)
    slp = np.zeros((4, n), dtype=float)
    slp[0, :3] = -3.0
    slp[2, :3] = 3.0
    slp[2, 3:6] = -3.0
    slp[0, 3:6] = 3.0
    np.savez_compressed(
        path,
        lat=lat,
        lon=np.linspace(-160.0, 160.0, n),
        cell_area=np.ones(n, dtype=float),
        earth__land_mask=np.ones(n, dtype=bool),
        earth__annual_precip_mm=annual_precip,
        earth__seasonal_precip_mm_yr_equiv=seasonal_precip,
        earth__seasonal_slp_anomaly_hPa=slp,
    )


def _write_generated(path: Path, *, preset: str, monsoon_value: float) -> dict:
    n = 8
    lat = np.array([-30, -20, -10, 10, 20, 30, 45, -45], dtype=float)
    seasonal_precip = np.full((4, n), 120.0)
    seasonal_precip[0, :3] = 1300.0
    seasonal_precip[2, 3:6] = 1300.0
    monsoon = np.zeros((4, n), dtype=float)
    monsoon[0, :3] = monsoon_value
    monsoon[2, 3:6] = monsoon_value
    moisture = np.full((4, n), 0.82, dtype=float)
    pressure = np.zeros((4, n), dtype=float)
    pressure[0, :3] = -1.1
    pressure[2, :3] = 0.3
    pressure[2, 3:6] = -1.1
    pressure[0, 3:6] = 0.3
    np.savez_compressed(
        path / "terminal_climate_arrays.npz",
        lat=lat,
        lon=np.linspace(-160.0, 160.0, n),
        cell_area=np.ones(n, dtype=float),
        sea_level_m=np.array([0.0]),
        terrain__elevation_m=np.ones(n, dtype=float),
        climate__precipitation=seasonal_precip.mean(axis=0),
        climate__seasonal_precipitation=seasonal_precip,
        atmosphere__moisture_access=moisture,
        atmosphere__monsoon_potential=monsoon,
        atmosphere__seasonal_pressure_proxy=pressure,
    )
    return {"preset": preset, "seed": 7, "assets_dir": str(path)}


def _write_summary(path: Path, rows: list[dict]) -> Path:
    import json

    path.write_text(json.dumps({"summaries": rows}))
    return path


def test_monsoon_moisture_gate_flags_waterworld_fake_monsoon(tmp_path):
    earth = tmp_path / "earth.npz"
    _write_earth(earth)
    world_dir = tmp_path / "waterworld_seed7"
    world_dir.mkdir()
    summary = _write_summary(
        tmp_path / "summary.json",
        [_write_generated(world_dir, preset="waterworld", monsoon_value=0.75)],
    )

    report = run_earth_climate_monsoon_moisture_gate(
        EarthClimateMonsoonMoistureGateConfig(
            earth_reference_npz=earth,
            terminal_summary_json=summary,
            outdir=tmp_path / "out",
        )
    )

    assert report["verdict"] == "fail"
    assert any(
        row["group"] == "waterworld_fake_monsoon_guard"
        and row["metric"] == "monsoon_potential_land_p95"
        for row in report["failures"]
    )
    assert (tmp_path / "out" / "earth_climate_monsoon_moisture_checks.csv").exists()


def test_monsoon_moisture_gate_cli_can_fail(tmp_path):
    earth = tmp_path / "earth.npz"
    _write_earth(earth)
    world_dir = tmp_path / "waterworld_seed7"
    world_dir.mkdir()
    summary = _write_summary(
        tmp_path / "summary.json",
        [_write_generated(world_dir, preset="waterworld", monsoon_value=0.75)],
    )

    args = SimpleNamespace(
        earth_reference=str(earth),
        terminal_summary=str(summary),
        out=str(tmp_path / "out"),
        fail_on_monsoon_moisture=True,
    )
    try:
        cmd_earth_climate_monsoon_moisture_gate(args)
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("expected CLI failure for fake waterworld monsoon")

