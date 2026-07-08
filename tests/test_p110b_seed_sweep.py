import json

from aevum.diagnostics.p110b_seed_sweep import (
    P110BSeedSweepThresholds,
    summarize_p110b_seed_sweep,
    write_p110b_seed_sweep_summary,
)


def _terminal_metrics(
    *,
    cells=8000,
    land_fraction=0.28,
    largest=0.52,
    second=0.22,
    third=0.12,
    fourth=0.04,
    major_land=5,
    ocean_basins=3,
    ocean_largest=0.45,
    closed=0.005,
    p109=True,
    p110a=True,
    p110a_out=None,
    warnings=(),
    ocean_source="ocean_basin_objects",
    archetype="modern_multipolar",
    robust_piece_count=3,
    internal_seaway_count=0,
    p1114_modern_endpoint_seaway_count=0,
):
    p110a_out_of_envelope = (
        []
        if p110a
        else list(p110a_out or ["largest_land_component_share_hard"])
    )
    return {
        "schema": "aevum.p107_terminal_world_audit.v1",
        "cells": cells,
        "land_fraction": land_fraction,
        "p109_hypsometry_comparison": {
            "within_p109_envelope": p109,
            "out_of_envelope": [] if p109 else ["p90_m"],
        },
        "p110a_modern_planform": {
            "summary": {
                "land_fraction": land_fraction,
                "major_land_component_count": major_land,
                "largest_land_component_share": largest,
                "second_land_component_share": second,
                "third_land_component_share": third,
                "major_ocean_basin_count": ocean_basins,
                "terminal_ocean_gateway_count": 3,
                "terminal_interbasin_ocean_gateway_count": 2,
                "terminal_phase_backed_ocean_gateway_count": 2,
                "terminal_ocean_gateway_system_count": 2,
                "terminal_interbasin_ocean_gateway_system_count": 1,
                "terminal_phase_backed_ocean_gateway_system_count": 2,
                "ocean_gateway_fragment_to_system_ratio": 1.5,
                "tectonic_ocean_gateway_count": 4,
                "restricted_ocean_fraction": 0.015,
                "unbacked_major_disconnected_ocean_component_count": 0,
                "largest_ocean_basin_share": ocean_largest,
                "closed_ocean_ring_score": closed,
            },
            "warning_flags": list(warnings),
            "out_of_envelope": p110a_out_of_envelope,
            "within_p110a_modern_planform_envelope": p110a,
            "p110b_final_state_archetype": {
                "code": 1 if archetype == "modern_multipolar" else 0,
                "name": archetype,
                "largest_share_preferred_ceiling": 0.54,
                "largest_share_soft_ceiling": 0.59,
                "min_nonlargest_large_components": 2,
            },
            "ocean_basins": {
                "source": ocean_source,
                "major_basin_count": ocean_basins,
                "largest_basin_share_of_ocean": ocean_largest,
            },
            "land": {
                "component_count": major_land,
                "major_component_count": major_land,
                "major_share_threshold": 0.03,
                "largest_component_share_of_mask": largest,
                "top_component_shares_of_mask": [
                    largest,
                    second,
                    third,
                    fourth,
                ],
            },
            "p110b_lineage_survival": {
                "lineage_supported_major_component_count": max(
                    0, min(major_land, 4)),
                "continental_major_component_count": max(0, min(major_land, 4)),
                "continental_lineage_supported_major_component_count": max(
                    0, min(major_land, 4)),
                "oceanic_landform_major_component_count": max(
                    0, major_land - min(major_land, 4)),
                "provenance_supported_major_component_count": major_land,
                "independent_primary_lineage_count": max(0, min(major_land, 4)),
                "lineage_area_fraction_of_land": 0.92,
                "geologic_support_fraction_of_land": 0.34,
                "min_major_dominant_continent_share": 0.72,
                "min_major_continent_lineage_area_share": 0.88,
                "min_major_geologic_support_fraction": 0.08,
            },
            "p110b_terminal_supercontinent_diagnostics": {
                "schema": "aevum.p110b_terminal_supercontinent_diagnostics.v1",
                "terminal_supercontinent_score": 0.12,
                "terminal_supercontinent_like": False,
                "modern_multipolar_overconnected": False,
                "largest_land_significant_continent_domain_count": 3,
                "largest_land_effective_continent_domain_count": 2.6,
                "largest_land_robust_piece_count_after_neck_removal": robust_piece_count,
                "largest_land_bridge_candidate_fraction": 0.03,
                "internal_domain_seaway_opening_count": internal_seaway_count,
                "internal_domain_seaway_area_fraction_world": (
                    0.004 if internal_seaway_count else 0.0
                ),
            },
            "ocean_gateway_topology": {
                "schema": "aevum.p110b_ocean_gateway_topology.v1",
                "terminal_gateway_count": 3,
                "terminal_interbasin_gateway_count": 2,
                "terminal_phase_backed_gateway_count": 2,
                "terminal_gateway_system_count": 2,
                "terminal_interbasin_gateway_system_count": 1,
                "terminal_phase_backed_gateway_system_count": 2,
                "gateway_fragment_to_system_ratio": 1.5,
                "tectonic_gateway_count": 4,
                "restricted_ocean_fraction": 0.015,
                "unbacked_major_disconnected_ocean_component_count": 0,
            },
            "seaway_cut_effectiveness": {
                "p1114_modern_endpoint_seaway_count": p1114_modern_endpoint_seaway_count,
                "p1114_modern_endpoint_seaway_area_fraction_world": (
                    0.004 if p1114_modern_endpoint_seaway_count else 0.0
                ),
                "p1114_modern_endpoint_seaway_domain_backed_count": (
                    p1114_modern_endpoint_seaway_count
                ),
            },
        },
        "island_arc_chain_count": 4,
        "microcontinent_object_count": 2,
        "parented_oceanic_island_chain_count": 3,
        "deep_trench_fraction_below_6000m": 0.012,
        "p110b_historical_supercontinent_trajectory": {
            "schema": "aevum.p110b_historical_supercontinent_trajectory.v1",
            "skipped": False,
            "usable_frame_count": 6,
            "supercontinent_frame_count": 1,
            "supercontinent_frame_fraction": 1.0 / 6.0,
            "max_largest_land_component_share": 0.66,
            "median_largest_land_component_share": 0.48,
            "max_consecutive_supercontinent_duration_myr": 0.0,
            "supercontinent_time_window_myr": 0.0,
            "long_lived_supercontinent_like": False,
        },
    }


def test_p110b_seed_sweep_loads_terminal_metrics_from_ladder_outdirs(tmp_path):
    root = tmp_path / "audit"
    run_a = root / "00_8000cells_36p"
    run_b = root / "01_8000cells_36p"
    run_a.mkdir(parents=True)
    run_b.mkdir(parents=True)
    (run_a / "p107_terminal_metrics.json").write_text(
        json.dumps(_terminal_metrics(largest=0.52)) + "\n")
    (run_b / "p107_terminal_metrics.json").write_text(
        json.dumps(_terminal_metrics(
            largest=0.63,
            second=0.14,
            third=0.09,
            warnings=("largest_land_component_share_soft",),
        )) + "\n")
    summary = {
        "schema": "aevum.p107_audit_ladder.v1",
        "entries": [
            {
                "label": "seed101",
                "preset": "earthlike",
                "cells": 8000,
                "n_plates": 36,
                "seed": 101,
                "outdir": str(run_a),
                "run_seconds": 12.0,
                "metrics": {"p110a_modern_planform": {"summary": {}}},
            },
            {
                "label": "seed202",
                "preset": "earthlike",
                "cells": 8000,
                "n_plates": 36,
                "seed": 202,
                "outdir": str(run_b),
                "run_seconds": 13.0,
                "metrics": {"p110a_modern_planform": {"summary": {}}},
            },
        ],
    }
    (root / "p107_audit_summary.json").write_text(json.dumps(summary) + "\n")

    report = summarize_p110b_seed_sweep(
        [root],
        thresholds=P110BSeedSweepThresholds(min_sample_size=2),
    )

    assert report["run_count"] == 2
    assert report["aggregate"]["threshold_pass_count"] == 2
    assert report["aggregate"]["p110b_visual_candidate_count"] == 1
    assert report["aggregate"]["largest_land_component_soft_warning_count"] == 1
    assert report["aggregate"]["object_backed_ocean_basin_count"] == 2
    assert (
        report["aggregate"]["p110b_independent_primary_lineage_count"]["median"]
        == 4.0
    )
    assert (
        report["aggregate"]["p110b_provenance_supported_major_component_count"]["median"]
        == 5.0
    )
    assert report["aggregate"]["p110b_terminal_ocean_gateway_count"]["median"] == 3.0
    assert (
        report["aggregate"]["p110b_terminal_phase_backed_ocean_gateway_count"]["median"]
        == 2.0
    )
    assert report["aggregate"]["p110b_terminal_ocean_gateway_system_count"]["median"] == 2.0
    assert (
        report["aggregate"]["p110b_terminal_interbasin_ocean_gateway_system_count"][
            "median"
        ]
        == 1.0
    )
    assert (
        report["aggregate"]["p110b_ocean_gateway_fragment_to_system_ratio"]["median"]
        == 1.5
    )
    assert report["aggregate"]["p110b_terminal_supercontinent_like_count"] == 0
    assert (
        report["aggregate"][
            "p110b_largest_land_significant_continent_domain_count"
        ]["median"]
        == 3.0
    )
    assert report["runs"][1]["p110b_status"] == "threshold_pass_needs_p110b_review"
    assert "largest_land_component_share_soft" in report["runs"][1]["p110b_review_flags"]


def test_p110b_seed_sweep_flags_terminal_supercontinent_like_world(tmp_path):
    metrics = _terminal_metrics(largest=0.55, second=0.18, third=0.08)
    diag = metrics["p110a_modern_planform"][
        "p110b_terminal_supercontinent_diagnostics"
    ]
    diag.update({
        "terminal_supercontinent_score": 0.72,
        "terminal_supercontinent_like": True,
        "modern_multipolar_overconnected": True,
        "largest_land_significant_continent_domain_count": 1,
        "largest_land_effective_continent_domain_count": 1.0,
        "largest_land_robust_piece_count_after_neck_removal": 1,
        "largest_land_bridge_candidate_fraction": 0.02,
    })
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(metrics) + "\n")

    report = summarize_p110b_seed_sweep(
        [metrics_path],
        thresholds=P110BSeedSweepThresholds(min_sample_size=1),
    )

    assert report["aggregate"]["p110b_terminal_supercontinent_like_count"] == 1
    assert "terminal_supercontinent_like_present" in (
        report["aggregate"]["distribution_flags"]
    )
    assert report["runs"][0]["p110b_status"] == "threshold_pass_needs_p110b_review"
    assert "terminal_supercontinent_like" in report["runs"][0]["p110b_review_flags"]


def test_p110b_seed_sweep_flags_long_lived_historical_supercontinent(tmp_path):
    metrics = _terminal_metrics(largest=0.50, second=0.20, third=0.12)
    metrics["p110b_historical_supercontinent_trajectory"].update({
        "usable_frame_count": 8,
        "supercontinent_frame_count": 6,
        "supercontinent_frame_fraction": 0.75,
        "max_largest_land_component_share": 1.0,
        "median_largest_land_component_share": 0.86,
        "max_consecutive_supercontinent_duration_myr": 2400.0,
        "supercontinent_time_window_myr": 2800.0,
        "long_lived_supercontinent_like": True,
    })
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(metrics) + "\n")

    report = summarize_p110b_seed_sweep(
        [metrics_path],
        thresholds=P110BSeedSweepThresholds(min_sample_size=1),
    )

    assert (
        report["aggregate"][
            "p110b_historical_long_lived_supercontinent_like_count"
        ]
        == 1
    )
    assert "historical_long_lived_supercontinent_like_present" in (
        report["aggregate"]["distribution_flags"]
    )
    assert report["runs"][0]["p110b_status"] == "threshold_pass_needs_p110b_review"
    assert "historical_long_lived_supercontinent_like" in (
        report["runs"][0]["p110b_review_flags"]
    )


def test_p110b_seed_sweep_allows_legacy_modern_balanced_three_major(tmp_path):
    metrics = _terminal_metrics(
        largest=0.533,
        second=0.224,
        third=0.220,
        major_land=3,
        p110a=False,
        p110a_out=("major_land_component_count_low",),
        archetype="modern_multipolar",
    )
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(metrics) + "\n")

    report = summarize_p110b_seed_sweep(
        [metrics_path],
        thresholds=P110BSeedSweepThresholds(min_sample_size=1),
    )

    run = report["runs"][0]
    assert run["p110a_pass"]
    assert run["p110a_reclassified_modern_balanced_three_major"]
    assert run["p110a_out_of_envelope"] == []
    assert (
        "major_land_component_count_low_p110b_modern_balanced_three_major"
        in run["p110a_warning_flags"]
    )
    assert "p110a_planform_failure" not in run["p110b_review_flags"]
    assert "major_land_component_count_low" not in run["p110b_review_flags"]
    assert report["aggregate"]["threshold_pass_count"] == 1
    assert run["p110b_status"] == "p110b_visual_candidate"


def test_p110b_seed_sweep_reads_historical_diagnostic_from_compact_ladder(tmp_path):
    metrics = _terminal_metrics(largest=0.50, second=0.20, third=0.12)
    metrics["p110b_historical_supercontinent_trajectory"].update({
        "usable_frame_count": 7,
        "supercontinent_frame_fraction": 0.72,
        "max_consecutive_supercontinent_duration_myr": 2100.0,
        "supercontinent_time_window_myr": 2400.0,
        "long_lived_supercontinent_like": True,
    })
    summary = {
        "schema": "aevum.p107_audit_ladder.v1",
        "entries": [{
            "label": "compact",
            "preset": "earthlike",
            "cells": 900,
            "n_plates": 36,
            "seed": 505,
            "outdir": "",
            "run_seconds": 1.0,
            "metrics": metrics,
        }],
    }
    (tmp_path / "p107_audit_summary.json").write_text(json.dumps(summary) + "\n")

    report = summarize_p110b_seed_sweep(
        [tmp_path],
        thresholds=P110BSeedSweepThresholds(min_sample_size=1),
    )

    assert report["runs"][0]["p110b_historical_long_lived_supercontinent_like"]
    assert "historical_long_lived_supercontinent_like_present" in (
        report["aggregate"]["distribution_flags"]
    )


def test_p110b_seed_sweep_flags_missing_object_backed_ocean_provenance(tmp_path):
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(_terminal_metrics(
        ocean_source="ocean.basin_id",
    )) + "\n")

    report = summarize_p110b_seed_sweep(
        [metrics_path],
        thresholds=P110BSeedSweepThresholds(min_sample_size=1),
    )

    assert report["run_count"] == 1
    assert report["aggregate"]["threshold_pass_count"] == 1
    assert report["aggregate"]["p110b_visual_candidate_count"] == 0
    assert report["aggregate"]["object_backed_ocean_basin_count"] == 0
    assert report["runs"][0]["p110b_review_flags"] == [
        "ocean_basin_source_not_object_backed"
    ]
    assert "object_backed_ocean_basin_source_incomplete" in (
        report["aggregate"]["distribution_flags"]
    )


def test_p111_seed_sweep_flags_two_body_endpoint_separately(tmp_path):
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(_terminal_metrics(
        largest=0.51,
        second=0.47,
        third=0.01,
        fourth=0.0,
        major_land=2,
    )) + "\n")

    report = summarize_p110b_seed_sweep(
        [metrics_path],
        thresholds=P110BSeedSweepThresholds(min_sample_size=1),
    )

    run = report["runs"][0]
    assert "p111_two_body_endpoint" in run["p111_modern_planform_review_flags"]
    assert (
        "p111_two_body_endpoint_present"
        in report["aggregate"]["p111_distribution_flags"]
    )
    assert report["aggregate"]["p111_two_body_endpoint_count"] == 1
    assert not report["acceptance"]["p111_modern_planform_distribution_ready"]


def test_p111_seed_sweep_flags_overlarge_weak_tertiary_candidate(tmp_path):
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(_terminal_metrics(
        largest=0.52,
        second=0.32,
        third=0.09,
        fourth=0.02,
        major_land=4,
    )) + "\n")

    report = summarize_p110b_seed_sweep(
        [metrics_path],
        thresholds=P110BSeedSweepThresholds(min_sample_size=1),
    )

    run = report["runs"][0]
    assert run["p110b_status"] == "p110b_visual_candidate"
    assert run["p111_modern_planform_status"] == "p111_modern_planform_review"
    assert "p111_largest_land_component_above_modern_preferred" in (
        run["p111_modern_planform_review_flags"]
    )
    assert "p111_weak_third_fourth_land_component" in (
        run["p111_modern_planform_review_flags"]
    )
    assert report["aggregate"]["p111_modern_planform_candidate_count"] == 0
    assert "p111_median_largest_landmass_above_modern_preferred" in (
        report["aggregate"]["p111_distribution_flags"]
    )


def test_p111_seed_sweep_accepts_p1114_endpoint_seaway_provenance(tmp_path):
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(_terminal_metrics(
        largest=0.57,
        second=0.22,
        third=0.12,
        fourth=0.04,
        major_land=5,
        robust_piece_count=1,
        internal_seaway_count=0,
        p1114_modern_endpoint_seaway_count=1,
    )) + "\n")

    report = summarize_p110b_seed_sweep(
        [metrics_path],
        thresholds=P110BSeedSweepThresholds(min_sample_size=1),
    )

    run = report["runs"][0]
    assert run["p1114_modern_endpoint_seaway_count"] == 1
    assert "p111_overlarge_landmass_lacks_internal_seaway_provenance" not in (
        run["p111_modern_planform_review_flags"]
    )
    assert report["aggregate"]["p1114_modern_endpoint_seaway_count"]["max"] == 1.0


def test_p110b_seed_sweep_writes_summary_json(tmp_path):
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(_terminal_metrics()) + "\n")
    report = summarize_p110b_seed_sweep(
        [metrics_path],
        thresholds=P110BSeedSweepThresholds(min_sample_size=1),
    )

    out_path = write_p110b_seed_sweep_summary(report, tmp_path / "summary_out")

    assert out_path.name == "p110b_seed_sweep_summary.json"
    saved = json.loads(out_path.read_text())
    assert saved["schema"] == "aevum.p110b_seed_sweep_summary.v1"
    assert saved["acceptance"]["p110b_seed_sweep_diagnostics_complete"]
