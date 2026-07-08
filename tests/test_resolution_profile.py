import json

from aevum.diagnostics.resolution_profile import (
    ResolutionProfileConfig,
    run_resolution_profile,
)


def test_resolution_profile_writes_timed_summary(tmp_path):
    config = ResolutionProfileConfig(
        cells=(240,),
        t_end_myr=120.0,
        frames=1,
        hex_width=16,
        hex_height=8,
        starts=1,
    )

    summary = run_resolution_profile(config, tmp_path)

    path = tmp_path / "resolution_profile_summary.json"
    assert path.exists()
    persisted = json.loads(path.read_text())
    assert persisted["schema"] == "aevum.resolution_profile.v1"

    entry = summary["entries"][0]
    assert entry["cells"] == 240
    assert entry["stage_seconds"]["build"] >= 0.0
    assert entry["stage_seconds"]["run"] >= 0.0
    assert entry["stage_seconds"]["compile"] >= 0.0
    assert entry["stage_seconds"]["tectonic_diagnostics"] >= 0.0
    assert entry["scheduler"]["macro_steps"] > 0
    assert entry["scheduler"]["module_run_counts"]["tectonics"] > 0
    assert 0.0 <= entry["world"]["land_fraction"] <= 1.0
    assert entry["compiler"]["compiled"]
    assert entry["compiler"]["width"] == 16
    assert entry["compiler"]["height"] == 8
    assert entry["morphology"]["land_component_count"] >= 0
    assert not entry["geomorphology_coverage"]["computed"]
    assert not summary["scaling"]["computed"]
    assert summary["config"]["projection_cells"] == [8000, 24000, 72000]
    assert summary["environment"]["active_array_backend"] == "numpy"
    assert "numba" in summary["environment"]["optional_acceleration"]
    assert summary["high_resolution"]["computed"]
    assert summary["high_resolution"]["readiness"] == "preflight_only"
    assert summary["high_resolution"]["resolution_tiers"][-1]["cells"] == "72000"
