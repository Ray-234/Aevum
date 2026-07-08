from aevum.diagnostics.p173_ocean_lifecycle_gate import (
    _aggregate_rows,
    _compact_world_summary,
    _parse_jobs,
)


def _p170_summary(ocean_values):
    rows = []
    for idx, (time_myr, unparented, fracture) in enumerate(ocean_values):
        rows.append({
            "time_myr": float(time_myr),
            "ocean_metrics": {
                "ocean_area_fraction_of_world": 0.6,
                "ocean_fabric_entropy": 0.5 + 0.02 * idx,
                "ridge_visible_fraction": 0.1,
                "ridge_age_symmetry_score": 0.7,
                "fracture_zone_length_fraction": float(fracture),
                "abyssal_plain_fraction": 0.35,
                "hotspot_track_count": 1.0,
                "seamount_chain_count": 1.0,
                "oceanic_plateau_fraction": 0.03,
                "microcontinent_fraction": 0.02,
                "unparented_shoal_fraction": float(unparented),
            },
        })
    return {
        "frame_count": len(rows),
        "usable_frame_count": len(rows),
        "required_metric_keys_present": True,
        "flag_summary": {
            "ordinary_plateau_frame_count": 0,
            "ordinary_deep_ocean_frame_count": 0,
        },
        "metric_extremes": {
            "ocean_metrics": {
                "unparented_shoal_fraction": {
                    "max": max(v for _, v, _ in ocean_values),
                    "median": 0.01,
                    "min": min(v for _, v, _ in ocean_values),
                },
                "fracture_zone_length_fraction": {
                    "max": max(v for _, _, v in ocean_values),
                    "median": 0.04,
                    "min": min(v for _, _, v in ocean_values),
                },
                "ocean_fabric_entropy": {
                    "max": 0.7,
                    "median": 0.6,
                    "min": 0.5,
                },
            }
        },
        "frame_rows": rows,
    }


def _p171_summary(collection_counts=None):
    counts = {
        "terrain.ocean_fabric": 4,
        "terrain.margin_landforms": 2,
        "terrain.arc_plume_landforms": 3,
        "terrain.rift_margin_sequences": 1,
    }
    if collection_counts is not None:
        counts.update(collection_counts)
    return {
        "frame_count": 4,
        "total_object_observations": 30,
        "unique_object_id_count": 20,
        "recurring_object_id_count": 6,
        "required_fields_complete": True,
        "missing_required_field_slot_count": 0,
        "collection_object_counts": counts,
        "acceptance": {"persistence_checked": True},
    }


def _response_globals(
    object_count=5,
    area=0.42,
    mean_abs=180.0,
    max_abs=650.0,
    p1732_candidate=0.04,
    p1732_adjusted=0.03,
):
    return {
        "terrain.last_p173_age_aware_ocean_response_object_count": object_count,
        "terrain.last_p173_age_aware_ocean_response_area_fraction": area,
        "terrain.last_p173_age_aware_ocean_response_mean_abs_delta_m": mean_abs,
        "terrain.last_p173_age_aware_ocean_response_max_abs_delta_m": max_abs,
        "terrain.last_p1732_young_open_ocean_age_depth_used": 1.0,
        "terrain.last_p1732_young_open_ocean_age_depth_land_mask_preserved": 1.0,
        "terrain.last_p1732_young_open_ocean_age_depth_candidate_fraction": (
            p1732_candidate
        ),
        "terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction": (
            p1732_adjusted
        ),
        "terrain.last_p1732_young_open_ocean_age_depth_mean_depth_before_m": 220.0,
        "terrain.last_p1732_young_open_ocean_age_depth_mean_depth_after_m": 3100.0,
        "terrain.last_p1115_endpoint_readability_guardrail_mode": 1.0,
        "terrain.last_p1115_process_relief_enabled": 0.0,
        "terrain.last_p1115_process_relief_adjusted_area_fraction": 0.0,
        "terrain.last_p1115_ocean_floor_adjusted_area_fraction": 0.02,
        "terrain.last_p1115_ocean_floor_land_mask_preserved": 1.0,
    }


def _p1731_summary(
    candidate=0.2,
    residual=0.01,
    p1732_candidate=0.0,
    p1732_adjusted=0.0,
):
    return {
        "frame_count": 3,
        "usable_frame_count": 3,
        "metric_extremes": {
            "cleanup_candidate_fraction_of_ocean": {
                "max": candidate,
                "median": candidate * 0.5,
                "min": 0.0,
            },
            "post_cleanup_residual_fraction_of_ocean": {
                "max": residual,
                "median": residual * 0.5,
                "min": 0.0,
            },
            "structural_preserve_fraction_of_candidate": {
                "max": 0.4,
                "median": 0.2,
                "min": 0.0,
            },
            "object_support_fraction_of_candidate": {
                "max": 0.3,
                "median": 0.2,
                "min": 0.0,
            },
            "semantic_support_fraction_of_candidate": {
                "max": 0.1,
                "median": 0.05,
                "min": 0.0,
            },
        },
        "p1732_metric_extremes": {
            "terrain.last_p1732_young_open_ocean_age_depth_candidate_fraction": {
                "max": p1732_candidate,
                "median": p1732_candidate * 0.5,
                "min": 0.0,
            },
            "terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction": {
                "max": p1732_adjusted,
                "median": p1732_adjusted * 0.5,
                "min": 0.0,
            },
        },
        "peak_residual_frame": {
            "time_myr": 1900.0,
            "owner_hint": "open_ocean_young_shallow_abyss_rise",
            "post_cleanup_residual_fraction_of_ocean": residual,
            "p1732_young_open_ocean_depth_floor": {
                "terrain.last_p1732_young_open_ocean_age_depth_candidate_fraction": (
                    p1732_candidate
                ),
                "terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction": (
                    p1732_adjusted
                ),
            },
        },
        "acceptance": {"attribution_completed": True},
    }


def test_p173_compact_world_summary_flags_unsupported_shoal_and_overbroad_response():
    row = _compact_world_summary(
        "earthlike_seed1",
        "earthlike",
        1,
        8000,
        12,
        _p170_summary([(500.0, 0.01, 0.02), (2600.0, 0.08, 0.11)]),
        _p171_summary(),
        _p1731_summary(
            candidate=0.33,
            residual=0.08,
            p1732_candidate=0.12,
            p1732_adjusted=0.10,
        ),
        world_globals=_response_globals(area=0.95, mean_abs=620.0),
        config={
            "min_ocean_area_fraction_for_gate": 0.10,
            "max_unparented_shoal_fraction": 0.05,
            "max_p173_response_area_fraction": 0.92,
            "max_p173_response_mean_abs_delta_m": 550.0,
        },
    )

    assert row["p173"]["ocean_gate_domain"]
    assert row["p173"]["unsupported_shoal_candidate"]
    assert row["p173"]["overbroad_ocean_response"]
    assert row["p173"]["overbroad_ocean_response_area_gate"]
    assert row["p173"]["overbroad_ocean_response_mean_abs_gate"]
    assert not row["p173"]["missing_ocean_lifecycle_objects"]
    assert not row["p173"]["missing_ocean_response"]
    assert row["p173"]["age_aware_ocean_response"][
        "terrain.last_p173_age_aware_ocean_response_object_count"] == 5
    assert row["p173_2"]["young_open_ocean_depth_floor"][
        "terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction"
    ] == 0.03
    assert row["p1115"]["terminal_guardrails"][
        "terrain.last_p1115_endpoint_readability_guardrail_mode"
    ] == 1.0
    assert row["p1115"]["terminal_guardrails"][
        "terrain.last_p1115_process_relief_enabled"
    ] == 0.0
    assert row["p173_1"]["attribution_completed"]
    assert row["p173_1"]["max_post_cleanup_residual_fraction_of_ocean"] == 0.08
    assert (
        row["p173_1"][
            "max_p1732_young_open_ocean_age_depth_adjusted_fraction"
        ]
        == 0.10
    )
    assert row["p173_1"]["peak_residual_frame"]["owner_hint"] == (
        "open_ocean_young_shallow_abyss_rise")


def test_p173_compact_world_summary_flags_missing_objects_and_response():
    row = _compact_world_summary(
        "waterworld_seed7",
        "waterworld",
        7,
        8000,
        12,
        _p170_summary([(500.0, 0.0, 0.0), (4500.0, 0.0, 0.0)]),
        _p171_summary({
            "terrain.ocean_fabric": 0,
            "terrain.margin_landforms": 0,
            "terrain.arc_plume_landforms": 0,
            "terrain.rift_margin_sequences": 0,
        }),
        world_globals=_response_globals(object_count=0, area=0.0, mean_abs=0.0),
        config={"min_ocean_area_fraction_for_gate": 0.10},
    )

    assert row["p173"]["ocean_gate_domain"]
    assert row["p173"]["missing_ocean_lifecycle_objects"]
    assert row["p173"]["missing_ocean_response"]
    assert row["p171"]["ocean_lifecycle_collection_counts"][
        "terrain.ocean_fabric"] == 0


def test_p173_compact_world_summary_allows_localized_high_amplitude_response():
    row = _compact_world_summary(
        "arid_seed101",
        "arid",
        101,
        8000,
        12,
        _p170_summary([(500.0, 0.0, 0.02), (4500.0, 0.0, 0.04)]),
        _p171_summary(),
        world_globals=_response_globals(area=0.01, mean_abs=615.0),
        config={
            "min_p173_response_area_fraction_for_mean_abs_gate": 0.05,
            "max_p173_response_mean_abs_delta_m": 550.0,
        },
    )

    assert not row["p173"]["overbroad_ocean_response"]
    assert not row["p173"]["overbroad_ocean_response_area_gate"]
    assert not row["p173"]["overbroad_ocean_response_mean_abs_gate"]


def test_p173_aggregate_reports_labels_and_acceptance_counts():
    clean = _compact_world_summary(
        "earthlike_seed2",
        "earthlike",
        2,
        8000,
        12,
        _p170_summary([(500.0, 0.01, 0.02), (4500.0, 0.02, 0.04)]),
        _p171_summary(),
        _p1731_summary(
            candidate=0.1,
            residual=0.01,
            p1732_candidate=0.10,
            p1732_adjusted=0.08,
        ),
        world_globals=_response_globals(p1732_candidate=0.10, p1732_adjusted=0.08),
        config={},
    )
    flagged = _compact_world_summary(
        "earthlike_seed3",
        "earthlike",
        3,
        8000,
        12,
        _p170_summary([(500.0, 0.01, 0.02), (4500.0, 0.07, 0.12)]),
        _p171_summary(),
        _p1731_summary(
            candidate=0.4,
            residual=0.07,
            p1732_candidate=0.20,
            p1732_adjusted=0.18,
        ),
        world_globals=_response_globals(
            area=0.97,
            mean_abs=620.0,
            p1732_candidate=0.20,
            p1732_adjusted=0.18,
        ),
        config={},
    )

    summary = _aggregate_rows(
        [clean, flagged],
        {
            "jobs": [
                {"preset": "earthlike", "label": "earthlike_seed2", "seed": 2},
                {"preset": "earthlike", "label": "earthlike_seed3", "seed": 3},
            ],
            "cells": 8000,
        },
        0.0,
    )

    assert summary["world_count"] == 2
    assert summary["aggregate"]["unsupported_shoal_world_count"] == 1
    assert summary["aggregate"]["unsupported_shoal_labels"] == ["earthlike_seed3"]
    assert summary["aggregate"]["overbroad_ocean_response_world_count"] == 1
    assert summary["aggregate"]["overbroad_ocean_response_labels"] == [
        "earthlike_seed3"]
    assert summary["aggregate"]["ocean_response_world_count"] == 2
    assert summary["aggregate"]["p1731_attribution_missing_world_count"] == 0
    assert summary["aggregate"]["max_p1731_cleanup_candidate_fraction"] == 0.4
    assert summary["aggregate"]["max_p1731_post_cleanup_residual_fraction"] == 0.07
    assert (
        summary["aggregate"][
            "max_p1732_young_open_ocean_age_depth_candidate_fraction"
        ]
        == 0.2
    )
    assert (
        summary["aggregate"][
            "max_p1732_young_open_ocean_age_depth_adjusted_fraction"
        ]
        == 0.18
    )
    assert summary["aggregate"]["p1115_endpoint_guardrail_mode_world_count"] == 2
    assert summary["aggregate"]["p1115_process_relief_enabled_world_count"] == 0
    assert summary["aggregate"]["max_p1115_process_relief_adjusted_area_fraction"] == 0.0
    assert summary["acceptance"]["p1115_terminal_process_relief_contracted"]
    assert summary["acceptance"]["gate_completed"]
    assert summary["acceptance"]["required_metric_keys_present"]
    assert summary["acceptance"]["object_fields_complete"]
    assert not summary["acceptance"]["unsupported_shoal_gate_passed"]
    assert not summary["acceptance"]["ocean_response_not_overbroad"]
    assert summary["acceptance"]["p1731_attribution_available"]


def test_p173_parse_cli_jobs():
    assert _parse_jobs(("earthlike:earthlike_seed42:42",)) == (
        ("earthlike", "earthlike_seed42", 42),
    )
