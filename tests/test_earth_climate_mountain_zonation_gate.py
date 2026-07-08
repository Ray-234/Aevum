import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_mountain_zonation_gate
from aevum.diagnostics.earth_climate_mountain_zonation_gate import (
    SCHEMA,
    EarthClimateMountainZonationGateConfig,
    run_earth_climate_mountain_zonation_gate,
)


def _write_earth(path):
    lat = np.array([-45, -30, -20, -10, 5, 20, 35, 50, 65, 70], dtype=float)
    elev = np.array([2600, 2300, 1800, 1200, 300, 200, 100, 400, 2100, 2600], dtype=float)
    temp = np.array([-7, -4, 2, 7, 22, 18, 14, 8, -8, -12], dtype=float)
    precip = np.array([550, 650, 900, 1200, 1600, 900, 700, 600, 420, 360], dtype=float)
    biome = np.array([5, 5, 3, 4, 6, 4, 4, 4, 5, 1], dtype=np.int16)
    np.savez(
        path,
        lat=lat,
        cell_area=np.ones(lat.size, dtype=float),
        earth__land_mask=np.ones(lat.size, dtype=bool),
        earth__elevation_m=elev,
        earth__annual_temperature_C=temp,
        earth__annual_precip_mm=precip,
        earth__biome_class_proxy=biome,
    )


def _write_generated(path, biome):
    lat = np.array([-45, -30, -20, -10, 5, 20, 35, 50, 65, 70], dtype=float)
    elev = np.array([2600, 2300, 1800, 1200, 300, 200, 100, 400, 2100, 2600], dtype=float)
    temp = np.array([-5, -2, 4, 8, 23, 19, 15, 9, -5, -8], dtype=float)
    precip = np.array([260, 280, 420, 700, 1200, 750, 600, 520, 230, 220], dtype=float)
    np.savez(
        path,
        lat=lat,
        cell_area=np.ones(lat.size, dtype=float),
        sea_level_m=np.array([0.0], dtype=float),
        terrain__elevation_m=elev,
        climate__surface_temperature=temp + 273.15,
        climate__precipitation=precip,
        biosphere__biome=np.asarray(biome, dtype=np.int16),
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


def test_mountain_zonation_gate_flags_high_mountain_desert(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(generated, [2, 2, 3, 4, 6, 4, 4, 4, 2, 2])
    summary.write_text(json.dumps(_summary(generated)))

    report = run_earth_climate_mountain_zonation_gate(
        EarthClimateMountainZonationGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["schema"] == SCHEMA
    assert report["verdict"] == "fail"
    failed = {row["metric"] for row in report["failures"]}
    assert "high_mountain_alpine_ecology_fraction" in failed
    assert "high_mountain_desert_fraction" in failed
    assert (tmp_path / "out" / "earth_climate_mountain_zonation_metrics.csv").exists()
    assert (tmp_path / "out" / "earth_climate_mountain_zonation_checks.csv").exists()
    assert (
        tmp_path / "out" / "earth_climate_mountain_zonation_gate_report.md"
    ).exists()


def test_mountain_zonation_gate_skips_diagnostic_only_world(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(generated, [2, 2, 3, 4, 6, 4, 4, 4, 2, 2])
    summary.write_text(json.dumps(_summary(generated, preset="waterworld")))

    report = run_earth_climate_mountain_zonation_gate(
        EarthClimateMountainZonationGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["checks"] == []
    assert report["generated_metrics"][0]["mode"] == "diagnostic_only"


def test_mountain_zonation_gate_cli_can_fail(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(generated, [2, 2, 3, 4, 6, 4, 4, 4, 2, 2])
    summary.write_text(json.dumps(_summary(generated)))

    args = SimpleNamespace(
        earth_reference=earth,
        terminal_summary=summary,
        out=tmp_path / "out",
        fail_on_mountain_zonation=True,
    )
    with pytest.raises(SystemExit) as exc_info:
        cmd_earth_climate_mountain_zonation_gate(args)

    assert exc_info.value.code == 2
    assert (
        tmp_path / "out" / "earth_climate_mountain_zonation_checks.csv"
    ).exists()
