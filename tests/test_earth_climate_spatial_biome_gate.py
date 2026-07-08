import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_spatial_biome_gate
from aevum.diagnostics.earth_climate_spatial_biome_gate import (
    SCHEMA,
    EarthClimateSpatialBiomeGateConfig,
    run_earth_climate_spatial_biome_gate,
)


def _write_earth(path):
    lat = np.array([-70, -50, -35, -25, -5, 5, 25, 35, 50, 70], dtype=float)
    biome = np.array([5, 4, 4, 2, 6, 6, 2, 2, 4, 1], dtype=np.int16)
    np.savez(
        path,
        lat=lat,
        cell_area=np.ones(lat.size, dtype=float),
        earth__land_mask=np.ones(lat.size, dtype=bool),
        earth__biome_class_proxy=biome,
    )


def _write_generated(path, biome):
    biome = np.asarray(biome, dtype=np.int16)
    np.savez(
        path,
        lat=np.array([-70, -50, -35, -25, -5, 5, 25, 35, 50, 70], dtype=float),
        cell_area=np.ones(biome.size, dtype=float),
        sea_level_m=np.array([0.0], dtype=float),
        terrain__elevation_m=np.ones(biome.size, dtype=float),
        biosphere__biome=biome,
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


def test_spatial_biome_gate_flags_missing_cool_midlat_forest(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(generated, [5, 2, 2, 2, 6, 6, 2, 2, 2, 1])
    summary.write_text(json.dumps(_summary(generated)))

    report = run_earth_climate_spatial_biome_gate(
        EarthClimateSpatialBiomeGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["schema"] == SCHEMA
    assert report["verdict"] == "fail"
    failed = {row["metric"] for row in report["failures"]}
    assert "cool_midlat_forest_fraction" in failed
    assert "cool_midlat_desert_fraction" in failed
    assert (tmp_path / "out" / "earth_climate_spatial_biome_metrics.csv").exists()
    assert (tmp_path / "out" / "earth_climate_spatial_biome_checks.csv").exists()


def test_spatial_biome_gate_skips_diagnostic_only_world(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(generated, [5, 2, 2, 2, 6, 6, 2, 2, 2, 1])
    summary.write_text(json.dumps(_summary(generated, preset="waterworld")))

    report = run_earth_climate_spatial_biome_gate(
        EarthClimateSpatialBiomeGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["checks"] == []
    assert report["generated_metrics"][0]["mode"] == "diagnostic_only"


def test_spatial_biome_gate_cli_can_fail(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(generated, [5, 2, 2, 2, 6, 6, 2, 2, 2, 1])
    summary.write_text(json.dumps(_summary(generated)))

    args = SimpleNamespace(
        earth_reference=earth,
        terminal_summary=summary,
        out=tmp_path / "out",
        fail_on_spatial_biome=True,
    )
    with pytest.raises(SystemExit) as exc_info:
        cmd_earth_climate_spatial_biome_gate(args)

    assert exc_info.value.code == 2
    assert (
        tmp_path / "out" / "earth_climate_spatial_biome_checks.csv"
    ).exists()
