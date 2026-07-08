import json

import numpy as np

from aevum.diagnostics.p107_equivalence import (
    compare_p107_outputs,
    resolve_p107_run_dir,
)


def _write_p107_run(root, label, arrays, metrics):
    run_dir = root / f"00_{label}"
    run_dir.mkdir(parents=True)
    np.savez(run_dir / "p107_terminal_arrays.npz", **arrays)
    (run_dir / "p107_terminal_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n")
    (root / "p107_audit_summary.json").write_text(
        json.dumps({
            "schema": "test",
            "entries": [{"outdir": str(run_dir)}],
        }) + "\n")
    return run_dir


def test_p107_equivalence_passes_when_only_profile_metadata_changes(tmp_path):
    arrays = {
        "elevation": np.asarray([0.0, 1.5, np.nan], dtype=np.float64),
        "plate": np.asarray([1, 2, 2], dtype=np.int32),
    }
    baseline_metrics = {
        "land_fraction": 0.29,
        "plate_counts": {"major": 4, "minor": 12},
        "terrain_internal_profile": {"semantic_object_build": 12.0},
        "assets": {"elevation_png": "baseline/path.png"},
        "array_archive": "baseline/archive.npz",
    }
    candidate_metrics = {
        "land_fraction": 0.29,
        "plate_counts": {"major": 4, "minor": 12},
        "terrain_internal_profile": {"semantic_object_build": 8.0},
        "assets": {"elevation_png": "candidate/path.png"},
        "array_archive": "candidate/archive.npz",
    }
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    baseline_run = _write_p107_run(
        baseline, "run", arrays, baseline_metrics)
    _write_p107_run(candidate, "run", arrays, candidate_metrics)

    assert resolve_p107_run_dir(baseline) == baseline_run
    assert resolve_p107_run_dir(baseline_run / "p107_terminal_metrics.json") == baseline_run

    report = compare_p107_outputs(baseline, candidate)

    assert report["equivalent"]
    assert report["arrays"]["common_count"] == 2
    assert report["arrays"]["changed_count"] == 0
    assert report["metrics"]["changed_count"] == 0


def test_p107_equivalence_fails_on_array_or_core_metric_change(tmp_path):
    baseline_arrays = {
        "elevation": np.asarray([0.0, 1.5, np.nan], dtype=np.float64),
        "plate": np.asarray([1, 2, 2], dtype=np.int32),
    }
    candidate_arrays = {
        "elevation": np.asarray([0.0, 2.5, np.nan], dtype=np.float64),
        "plate": np.asarray([1, 2, 2], dtype=np.int32),
    }
    baseline_metrics = {
        "land_fraction": 0.29,
        "plate_counts": {"major": 4, "minor": 12},
        "terrain_internal_profile": {"semantic_object_build": 12.0},
    }
    candidate_metrics = {
        "land_fraction": 0.31,
        "plate_counts": {"major": 4, "minor": 12},
        "terrain_internal_profile": {"semantic_object_build": 8.0},
    }
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    _write_p107_run(baseline, "run", baseline_arrays, baseline_metrics)
    _write_p107_run(candidate, "run", candidate_arrays, candidate_metrics)

    report = compare_p107_outputs(baseline, candidate)

    assert not report["equivalent"]
    assert [item["key"] for item in report["arrays"]["changed"]] == ["elevation"]
    assert [item["path"] for item in report["metrics"]["changed"]] == [
        "land_fraction"
    ]
