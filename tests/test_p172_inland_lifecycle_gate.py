from aevum.diagnostics.p172_inland_lifecycle_gate import (
    _aggregate_rows,
    _compact_world_summary,
    _parse_jobs,
)


def _p170_summary(rift_values):
    rows = []
    for idx, (time_myr, value) in enumerate(rift_values):
        rows.append({
            "time_myr": float(time_myr),
            "land_metrics": {
                "inland_detail_entropy": 0.6 + 0.01 * idx,
                "inland_area_fraction_of_world": 0.12,
                "ordinary_plateau_fraction": 0.0,
                "broad_flat_inland_component_count": 0,
                "continent_province_count_p50": 4.0,
                "old_orogen_expression_fraction": 0.1,
                "rift_basin_expression_fraction": float(value),
                "craton_shield_platform_split_fraction": 0.8,
                "lowland_plain_fraction": 0.02 + 0.01 * idx,
                "broad_lowland_plain_component_count": int(idx > 1),
                "largest_lowland_plain_component_fraction": 0.008 + 0.004 * idx,
                "lowland_plain_parented_fraction": 0.9,
            },
        })
    mature_deficient = sum(
        1 for row in rows
        if row["land_metrics"]["lowland_plain_fraction"] < 0.06
    )
    return {
        "frame_count": len(rows),
        "usable_frame_count": len(rows),
        "required_metric_keys_present": True,
        "flag_summary": {
            "ordinary_plateau_frame_count": 0,
            "ordinary_deep_ocean_frame_count": 0,
            "lowland_plain_deficient_frame_count": mature_deficient,
        },
        "p174_continuity": {
            "mature_support_frame_count": len(rows),
            "mature_lowland_deficient_frame_count": mature_deficient,
            "mature_lowland_deficient_fraction": mature_deficient / max(len(rows), 1),
            "mature_lowland_continuity_score": 1.0 - mature_deficient / max(len(rows), 1),
            "lowland_plain_max_positive_step": 0.01,
            "lowland_plain_max_negative_step": 0.0,
            "terminal_lowland_plain_jump": 0.02,
            "terminal_lowland_pop_candidate": mature_deficient > 1,
            "ocean_fabric_entropy_max_positive_step": 0.0,
            "ocean_fabric_entropy_max_negative_step": 0.0,
            "terminal_ocean_fabric_entropy_jump": 0.0,
        },
        "metric_extremes": {
            "land_metrics": {
                "rift_basin_expression_fraction": {
                    "max": max(v for _, v in rift_values),
                    "median": 0.5,
                    "min": min(v for _, v in rift_values),
                },
                "inland_detail_entropy": {"max": 0.8, "median": 0.7, "min": 0.6},
                "lowland_plain_fraction": {"max": 0.08, "median": 0.05, "min": 0.02},
            }
        },
        "frame_rows": rows,
    }


def _p171_summary():
    return {
        "frame_count": 4,
        "total_object_observations": 30,
        "unique_object_id_count": 20,
        "recurring_object_id_count": 5,
        "required_fields_complete": True,
        "missing_required_field_slot_count": 0,
        "collection_object_counts": {
            "terrain.continental_landforms": 8,
            "tectonics.continental_provinces": 12,
            "terrain.mountain_ranges": 4,
            "terrain.plateau_inventory": 1,
        },
        "acceptance": {"persistence_checked": True},
    }


def test_p172_compact_world_summary_flags_mature_rift_residual():
    row = _compact_world_summary(
        "earthlike_seed1",
        "earthlike",
        1,
        8000,
        12,
        _p170_summary([(500.0, 0.2), (2600.0, 0.72), (4500.0, 0.4)]),
        _p171_summary(),
        world_globals={
            "terrain.last_p172_age_aware_inland_response_object_count": 4.0,
            "terrain.last_p172_age_aware_inland_response_area_fraction": 0.31,
            "terrain.last_p172_age_aware_inland_response_mean_abs_delta_m": 88.0,
            "terrain.last_p172_age_aware_inland_response_max_abs_delta_m": 220.0,
        },
        config={
            "residual_rift_mature_max_threshold": 0.65,
            "residual_rift_any_max_threshold": 0.75,
        },
    )

    assert row["p172"]["residual_rift_overpaint_candidate"]
    assert row["p172"]["mature_rift_basin_expression_max"] == 0.72
    assert row["p172"]["age_aware_inland_response"][
        "terrain.last_p172_age_aware_inland_response_object_count"] == 4.0
    assert row["p172"]["age_aware_inland_response"][
        "terrain.last_p172_age_aware_inland_response_mean_abs_delta_m"] == 88.0
    assert row["p170"]["p174_continuity"]["mature_support_frame_count"] == 3
    assert row["p170"]["p174_continuity"]["terminal_lowland_pop_candidate"] == 1
    assert row["p171"]["required_fields_complete"]
    assert row["acceptance"]["object_persistence_checked"]


def test_p172_aggregate_reports_residual_and_gate_acceptance_counts():
    clean = _compact_world_summary(
        "earthlike_seed2",
        "earthlike",
        2,
        8000,
        12,
        _p170_summary([(500.0, 0.1), (2600.0, 0.2), (4500.0, 0.3)]),
        _p171_summary(),
        config={},
    )
    flagged = _compact_world_summary(
        "earthlike_seed3",
        "earthlike",
        3,
        8000,
        12,
        _p170_summary([(500.0, 0.1), (2600.0, 0.8), (4500.0, 0.2)]),
        _p171_summary(),
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
    assert summary["aggregate"]["residual_rift_overpaint_world_count"] == 1
    assert summary["aggregate"]["residual_rift_overpaint_labels"] == [
        "earthlike_seed3"]
    assert summary["aggregate"]["age_aware_inland_response_world_count"] == 0
    assert summary["aggregate"]["p174_lowland_continuity_risk_world_count"] == 2
    assert summary["aggregate"]["min_p174_lowland_continuity_score"] < 1.0
    assert summary["acceptance"]["gate_completed"]
    assert summary["acceptance"]["required_metric_keys_present"]
    assert summary["acceptance"]["object_fields_complete"]


def test_p172_parse_cli_jobs():
    assert _parse_jobs(("earthlike:earthlike_seed42:42",)) == (
        ("earthlike", "earthlike_seed42", 42),
    )
