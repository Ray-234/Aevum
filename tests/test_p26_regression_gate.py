from aevum.diagnostics.p26_regression_gate import compare_p12_summaries


def _summary(
    *,
    status="warn",
    failures=None,
    validation_passed=True,
    land=0.253,
    components=4,
    ribbon=0.410,
    continental_ribbon=0.397,
    coastline=21.46,
    basins=14,
    seam_basin=0.0,
    seam_render=0.0,
):
    hard = [] if validation_passed else ["ocean basin ids break across antimeridian seam"]
    failures = [] if failures is None else list(failures)
    return {
        "schema": "aevum.p12_tectonics_release_summary.v1",
        "release_decision": {"status": status, "passed": status != "fail"},
        "entries": [
            {
                "preset": "earthlike",
                "land_fraction": land,
                "validation": {
                    "diagnostics.tectonics": {
                        "passed": validation_passed,
                        "hard_failures": hard,
                    }
                },
                "release_gate": {
                    "passed": not failures,
                    "failures": failures,
                    "warnings": [],
                },
                "morphology": {
                    "land_component_count": components,
                    "land_ribbon_fraction_gt_0_5": ribbon,
                    "continental_ribbon_fraction_gt_0_5": continental_ribbon,
                    "land_coastline_complexity_largest": coastline,
                },
                "ocean_geography": {"basin_count": basins},
                "seam_continuity": {
                    "seam_ocean_basin_mismatch_fraction": seam_basin,
                    "render_duplicate_basin_mismatch_fraction": seam_render,
                },
            }
        ],
    }


def test_p26_regression_gate_accepts_non_regressed_candidate():
    baseline = _summary()
    candidate = _summary(
        land=0.255,
        components=4,
        ribbon=0.405,
        continental_ribbon=0.390,
        coastline=20.5,
        basins=13,
    )

    report = compare_p12_summaries(baseline, candidate)

    assert report["schema"] == "aevum.p26_regression_gate.v1"
    assert report["status"] == "pass"
    assert report["acceptance"]["no_world_level_regressions"]
    assert report["failed_checks"] == []


def test_p26_regression_gate_accepts_candidate_moving_from_underpartitioned_basins():
    baseline = _summary(basins=2)
    candidate = _summary(basins=4)

    report = compare_p12_summaries(baseline, candidate)

    assert report["status"] == "pass"
    assert report["acceptance"]["ocean_partition_not_regressed"]
    assert report["failed_checks"] == []


def test_p26_regression_gate_rejects_failed_margin_spike_shape():
    baseline = _summary()
    candidate = _summary(
        status="fail",
        failures=["Earthlike land is dominated by ribbon-like landforms"],
        land=0.216,
        components=16,
        ribbon=0.565,
        continental_ribbon=0.433,
        coastline=34.53,
        basins=40,
    )

    report = compare_p12_summaries(baseline, candidate)

    assert report["status"] == "fail"
    failed = set(report["failed_checks"])
    assert "release_and_validation_pass" in failed
    assert "land_fraction_no_regression" in failed
    assert "land_component_count_no_regression" in failed
    assert "land_ribbon_no_regression" in failed
    assert "largest_coastline_complexity_no_regression" in failed
    assert "ocean_basin_count_no_regression" in failed


def test_p26_regression_gate_rejects_basin_seam_regression():
    baseline = _summary()
    candidate = _summary(
        status="fail",
        validation_passed=False,
        seam_basin=0.05,
        seam_render=0.02,
    )

    report = compare_p12_summaries(baseline, candidate)

    failed = set(report["failed_checks"])
    assert "release_and_validation_pass" in failed
    assert "ocean_basin_seam_continuity" in failed
    assert not report["acceptance"]["ocean_partition_not_regressed"]
