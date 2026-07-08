import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_biome_gate
from aevum.diagnostics.earth_climate_biome_gate import (
    SCHEMA,
    EarthClimateBiomeGateConfig,
    run_earth_climate_biome_gate,
)


def _write_earth(path):
    biome_proxy = np.array([0, 6, 6, 4, 4, 4, 2, 3, 5, 1], dtype=np.int16)
    resolve = np.array([0, 1, 2, 4, 5, 6, 13, 8, 11, 11], dtype=np.int16)
    land = biome_proxy > 0
    np.savez(
        path,
        cell_area=np.ones(biome_proxy.size, dtype=float),
        earth__land_mask=land,
        earth__biome_class_proxy=biome_proxy,
        earth__resolve_biome_class=resolve,
    )


def _write_generated(path, biome):
    biome = np.asarray(biome, dtype=np.int16)
    np.savez(
        path,
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


def test_earth_climate_biome_gate_flags_low_forest_tropical(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(generated, [3, 3, 3, 3, 2, 3, 2, 3, 5, 1])
    summary.write_text(json.dumps(_summary(generated)))

    report = run_earth_climate_biome_gate(
        EarthClimateBiomeGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["schema"] == SCHEMA
    assert report["verdict"] == "fail"
    assert any(
        row["metric"] == "forest_tropical_land_fraction"
        for row in report["failures"]
    )
    assert (tmp_path / "out" / "earth_climate_biome_metrics.csv").exists()
    assert (tmp_path / "out" / "earth_climate_biome_checks.csv").exists()
    assert (tmp_path / "out" / "earth_climate_biome_gate_report.md").exists()


def test_earth_climate_biome_gate_skips_diagnostic_only_world(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(generated, [3, 3, 3, 3, 2, 3, 2, 3, 5, 1])
    summary.write_text(json.dumps(_summary(generated, preset="waterworld")))

    report = run_earth_climate_biome_gate(
        EarthClimateBiomeGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["checks"] == []
    assert report["generated_metrics"][0]["mode"] == "diagnostic_only"


def test_earth_climate_biome_gate_cli_can_fail_on_biome(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(generated, [3, 3, 3, 3, 2, 3, 2, 3, 5, 1])
    summary.write_text(json.dumps(_summary(generated)))

    args = SimpleNamespace(
        earth_reference=earth,
        terminal_summary=summary,
        out=tmp_path / "out",
        fail_on_biome=True,
    )
    with pytest.raises(SystemExit) as exc_info:
        cmd_earth_climate_biome_gate(args)

    assert exc_info.value.code == 2
    assert (tmp_path / "out" / "earth_climate_biome_checks.csv").exists()
