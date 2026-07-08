import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_pattern_gate
from aevum.diagnostics.earth_climate_pattern_gate import (
    SCHEMA,
    EarthClimatePatternGateConfig,
    run_earth_climate_pattern_gate,
)


def _write_earth_reference(path):
    lat = np.array([-70, -60, -30, -25, -20, -10, 0, 10, 20, 25, 60, 70], dtype=float)
    n = lat.size
    land = np.ones(n, dtype=bool)
    precip = np.array([300, 400, 80, 120, 180, 1600, 2400, 2600, 1700, 120, 350, 320], dtype=float)
    seasonal = np.vstack([
        0.40 * precip,
        0.80 * precip,
        1.60 * precip,
        1.20 * precip,
    ])
    temp = np.array([-15, -12, 20, 21, 24, 26, 27, 27, 25, 23, -8, -18], dtype=float)
    biome = np.array([5, 5, 2, 2, 2, 6, 6, 6, 4, 2, 5, 1], dtype=np.int16)
    np.savez(
        path,
        lat=lat,
        cell_area=np.ones(n, dtype=float),
        earth__land_mask=land,
        earth__elevation_m=np.array([200, 300, 800, 1200, 1300, 500, 200, 1500, 1100, 800, 400, 300], dtype=float),
        earth__annual_temperature_C=temp,
        earth__annual_precip_mm=precip,
        earth__seasonal_precip_mm_yr_equiv=seasonal,
        earth__biome_class_proxy=biome,
    )


def _write_generated(path, *, wet=True, cold_biome=True):
    lat = np.array([-70, -60, -30, -25, -20, -10, 0, 10, 20, 25, 60, 70], dtype=float)
    n = lat.size
    if wet:
        precip = np.array([260, 280, 90, 140, 180, 1500, 2200, 2300, 1500, 160, 280, 260], dtype=float)
        temp_c = np.array([-12, -10, 20, 21, 24, 26, 27, 27, 25, 23, -7, -15], dtype=float)
        biome = np.array([5, 5, 2, 2, 2, 6, 6, 6, 4, 2, 5, 1], dtype=np.int16)
    else:
        precip = np.array([700, 650, 520, 560, 620, 650, 760, 820, 720, 580, 650, 620], dtype=float)
        temp_c = np.array([4, 3, 20, 21, 24, 26, 27, 27, 25, 23, 4, 2], dtype=float)
        biome = np.array([3, 3, 3, 3, 3, 4, 4, 6, 4, 3, 3, 3], dtype=np.int16)
    if not cold_biome:
        biome[[0, 1, 10, 11]] = 3
    seasonal = np.vstack([
        0.40 * precip,
        0.80 * precip,
        1.60 * precip,
        1.20 * precip,
    ])
    np.savez(
        path,
        lat=lat,
        cell_area=np.ones(n, dtype=float),
        sea_level_m=np.array([0.0], dtype=float),
        terrain__elevation_m=np.array([200, 300, 800, 1200, 1300, 500, 200, 1500, 1100, 800, 400, 300], dtype=float),
        climate__surface_temperature=temp_c + 273.15,
        climate__precipitation=precip,
        climate__seasonal_precipitation=seasonal,
        biosphere__biome=biome,
    )


def _summary(path, *, preset="earthlike_mobile_lid"):
    return {
        "schema": "aevum.terminal_climate_replay.v1",
        "summaries": [{
            "preset": preset,
            "seed": 1,
            "arrays": str(path),
            "assets_dir": str(path.parent / "earthlike_seed1"),
        }],
    }


def test_earth_climate_pattern_gate_flags_spatial_failures(tmp_path):
    earth_path = tmp_path / "earth.npz"
    generated_path = tmp_path / "generated.npz"
    _write_earth_reference(earth_path)
    _write_generated(generated_path, wet=False, cold_biome=False)
    summary_path = tmp_path / "terminal.json"
    summary_path.write_text(json.dumps(_summary(generated_path)))

    report = run_earth_climate_pattern_gate(
        EarthClimatePatternGateConfig(earth_path, summary_path, tmp_path / "out")
    )

    assert report["schema"] == SCHEMA
    assert report["verdict"] == "fail"
    failed_metrics = {row["metric"] for row in report["failures"]}
    assert "wet_tropics_precip_p90_mm_yr" in failed_metrics
    assert "dry_subtropics_fraction_lt250mm" in failed_metrics
    assert "high_lat_cold_fraction_lt0C" in failed_metrics
    assert (tmp_path / "out" / "earth_climate_pattern_metrics.csv").exists()
    assert (tmp_path / "out" / "earth_climate_pattern_checks.csv").exists()
    assert (tmp_path / "out" / "earth_climate_pattern_gate_report.md").exists()


def test_earth_climate_pattern_gate_handles_diagnostic_only_world(tmp_path):
    earth_path = tmp_path / "earth.npz"
    generated_path = tmp_path / "generated.npz"
    _write_earth_reference(earth_path)
    _write_generated(generated_path, wet=False, cold_biome=False)
    summary_path = tmp_path / "terminal.json"
    summary_path.write_text(json.dumps(_summary(generated_path, preset="waterworld")))

    report = run_earth_climate_pattern_gate(
        EarthClimatePatternGateConfig(earth_path, summary_path, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["failure_count"] == 0
    assert report["checks"] == []
    assert report["generated_metrics"][0]["mode"] == "diagnostic_only"


def test_earth_climate_pattern_gate_cli_can_fail_on_pattern(tmp_path):
    earth_path = tmp_path / "earth.npz"
    generated_path = tmp_path / "generated.npz"
    _write_earth_reference(earth_path)
    _write_generated(generated_path, wet=False, cold_biome=False)
    summary_path = tmp_path / "terminal.json"
    summary_path.write_text(json.dumps(_summary(generated_path)))

    args = SimpleNamespace(
        earth_reference=earth_path,
        terminal_summary=summary_path,
        out=tmp_path / "out",
        fail_on_pattern=True,
    )
    with pytest.raises(SystemExit) as exc_info:
        cmd_earth_climate_pattern_gate(args)

    assert exc_info.value.code == 2
    assert (tmp_path / "out" / "earth_climate_pattern_checks.csv").exists()
