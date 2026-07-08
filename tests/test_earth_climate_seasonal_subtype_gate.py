import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_seasonal_subtype_gate
from aevum.diagnostics.earth_climate_seasonal_subtype_gate import (
    SCHEMA,
    EarthClimateSeasonalSubtypeGateConfig,
    run_earth_climate_seasonal_subtype_gate,
)


def _write_earth(path):
    lat = np.array([-20, -10, -5, 5, 10, 20, 35, 50], dtype=float)
    annual = np.array([700, 1200, 1200, 1200, 1200, 700, 500, 500], dtype=float)
    seasonal = np.array([
        [1200, 1800, 300, 300, 1700, 1200, 700, 450],
        [300, 300, 1700, 1700, 400, 300, 350, 500],
        [1200, 1800, 300, 300, 1700, 1200, 700, 650],
        [100, 900, 2500, 2500, 1000, 100, 250, 400],
    ], dtype=float)
    np.savez(
        path,
        lat=lat,
        cell_area=np.ones(lat.size, dtype=float),
        earth__land_mask=np.ones(lat.size, dtype=bool),
        earth__annual_precip_mm=annual,
        earth__seasonal_precip_mm_yr_equiv=seasonal,
    )


def _write_generated(path, seasonal):
    seasonal = np.asarray(seasonal, dtype=float)
    np.savez(
        path,
        lat=np.array([-20, -10, -5, 5, 10, 20, 35, 50], dtype=float),
        cell_area=np.ones(seasonal.shape[1], dtype=float),
        sea_level_m=np.array([0.0], dtype=float),
        terrain__elevation_m=np.ones(seasonal.shape[1], dtype=float),
        climate__precipitation=seasonal.mean(axis=0),
        climate__seasonal_precipitation=seasonal,
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


def test_seasonal_subtype_gate_flags_missing_low_tropical_dry_season(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(
        generated,
        np.repeat(np.array([[1200, 1200, 1200, 1200, 1200, 700, 500, 500]], dtype=float), 4, axis=0),
    )
    summary.write_text(json.dumps(_summary(generated)))

    report = run_earth_climate_seasonal_subtype_gate(
        EarthClimateSeasonalSubtypeGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["schema"] == SCHEMA
    assert report["verdict"] == "fail"
    failed = {row["metric"] for row in report["failures"]}
    assert "low_tropics_dry_quarter_ge2_fraction" in failed
    assert (tmp_path / "out" / "earth_climate_seasonal_subtype_metrics.csv").exists()
    assert (tmp_path / "out" / "earth_climate_seasonal_subtype_checks.csv").exists()


def test_seasonal_subtype_gate_skips_diagnostic_only_world(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(
        generated,
        np.repeat(np.array([[1200, 1200, 1200, 1200, 1200, 700, 500, 500]], dtype=float), 4, axis=0),
    )
    summary.write_text(json.dumps(_summary(generated, preset="waterworld")))

    report = run_earth_climate_seasonal_subtype_gate(
        EarthClimateSeasonalSubtypeGateConfig(earth, summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["checks"] == []
    assert report["generated_metrics"][0]["mode"] == "diagnostic_only"


def test_seasonal_subtype_gate_cli_can_fail(tmp_path):
    earth = tmp_path / "earth.npz"
    generated = tmp_path / "generated.npz"
    summary = tmp_path / "summary.json"
    _write_earth(earth)
    _write_generated(
        generated,
        np.repeat(np.array([[1200, 1200, 1200, 1200, 1200, 700, 500, 500]], dtype=float), 4, axis=0),
    )
    summary.write_text(json.dumps(_summary(generated)))

    args = SimpleNamespace(
        earth_reference=earth,
        terminal_summary=summary,
        out=tmp_path / "out",
        fail_on_seasonal_subtype=True,
    )
    with pytest.raises(SystemExit) as exc_info:
        cmd_earth_climate_seasonal_subtype_gate(args)

    assert exc_info.value.code == 2
    assert (
        tmp_path / "out" / "earth_climate_seasonal_subtype_checks.csv"
    ).exists()
