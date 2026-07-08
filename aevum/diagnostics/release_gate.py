"""P12 tectonics release-gate diagnostics.

This module is deliberately downstream of generation.  It runs worlds, collects
the existing P0-P11 diagnostics into a stable JSON summary, and renders a compact
contact sheet for manual review.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aevum.compiler.map_compiler import MapCompiler
from aevum.engine import Engine
from aevum.modules.tectonics import P39_OLD_PLATFORM_FUNNEL_METRICS
from aevum.spec.presets import get_preset


P12_CANONICAL_PRESETS = (
    "earthlike",
    "waterworld",
    "arid",
    "stagnant_lid",
    "tidally_locked",
    "frozen",
)


@dataclass
class P12RunConfig:
    presets: tuple[str, ...] = P12_CANONICAL_PRESETS
    cells: int = 3000
    t_end_myr: float | None = None
    frames: int = 4
    hex_width: int = 64
    hex_height: int = 32
    starts: int = 4
    render_world_assets: bool = False
    global_overrides: dict[str, float] = field(default_factory=dict)


def run_p12_release_gate(config: P12RunConfig, outdir: Path | None = None) -> dict[str, Any]:
    """Run the P12 preset matrix and optionally write summary/assets.

    The return value is JSON-serializable and intentionally compact: detailed
    diagnostics remain available through the normal single-world render outputs.
    """
    from aevum import render, validation

    if outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)

    earth_geomorphology_benchmarks = (
        _earth_geomorphology_fixture_summary()
        if any(str(preset).startswith("earthlike") for preset in config.presets)
        else {}
    )

    entries: list[dict[str, Any]] = []
    contact_panels: list[dict[str, Any]] = []
    for preset in config.presets:
        spec = get_preset(preset)
        spec.grid_cells = int(config.cells)
        if config.t_end_myr is not None:
            spec.t_end_myr = float(config.t_end_myr)

        eng = Engine.build(spec)
        for key, value in sorted(config.global_overrides.items()):
            eng.world.set_g(str(key), float(value))
        eng.run(n_frames=config.frames)
        cm = MapCompiler(eng.world, eng.archive).compile(
            width=config.hex_width,
            height=config.hex_height,
            n_starts=config.starts,
        )

        checks = validation.run_all(eng)
        tectonics = validation.tectonic_diagnostics(eng)
        climate = validation.climate_diagnostics(eng)
        compiler = validation.compiler_consistency_metrics(eng, cm)
        geomorphology = validation.earth_geomorphology_coverage_metrics(
            eng, cm, fixture_suites=earth_geomorphology_benchmarks)
        from aevum.diagnostics.earth_reference import earth_reference_distribution_metrics

        earth_reference = earth_reference_distribution_metrics(
            eng.world, tectonics=tectonics, geomorphology=geomorphology)
        entry = _summarize_entry(
            preset, eng, cm, checks, tectonics, climate, compiler,
            geomorphology, earth_reference)
        entry["release_gate"] = _assess_entry(entry)
        entries.append(entry)
        contact_panels.append({"preset": preset, "world": eng.world, "compiled": cm})

        if outdir is not None and config.render_world_assets:
            preset_dir = outdir / f"{preset}_seed{spec.seed}"
            preset_dir.mkdir(parents=True, exist_ok=True)
            render.render_world(eng.world, preset_dir)
            render.render_hexmap(cm, preset_dir)
            render.render_compiler_consistency(cm, preset_dir)
            render.render_timeline(eng.archive.timeline(), preset_dir)
            render.render_history(eng.archive, preset_dir)
            render.render_archive_continuity(eng, preset_dir)

    summary = {
        "schema": "aevum.p12_tectonics_release_summary.v1",
        "config": {
            "presets": list(config.presets),
            "cells": int(config.cells),
            "t_end_myr": None if config.t_end_myr is None else float(config.t_end_myr),
            "frames": int(config.frames),
            "hex_width": int(config.hex_width),
            "hex_height": int(config.hex_height),
            "starts": int(config.starts),
            "global_overrides": {
                str(key): float(value)
                for key, value in sorted(config.global_overrides.items())
            },
        },
        "earth_geomorphology_benchmarks": _compact_geomorphology_benchmarks(
            earth_geomorphology_benchmarks),
        "entries": entries,
    }
    summary["release_decision"] = _release_decision(entries)

    if outdir is not None:
        (outdir / "p12_tectonics_release_summary.json").write_text(
            json.dumps(summary, indent=2, default=_json_default) + "\n")
        _render_contact_sheet(contact_panels, outdir / "p12_preset_matrix_contact_sheet.png")
    return summary


def _earth_geomorphology_fixture_summary() -> dict[str, Any]:
    from aevum.diagnostics.tectonics_bench import run_suite

    return {suite: run_suite(suite) for suite in ("E1", "E2", "E3", "E4", "E5")}


def _compact_geomorphology_benchmarks(summaries: dict[str, Any]) -> dict[str, Any]:
    suites: dict[str, Any] = {}
    for suite, summary in summaries.items():
        benches = summary.get("benchmarks", [])
        suites[suite] = {
            "status": summary.get("status"),
            "benchmark_count": len(benches),
            "passed_count": sum(1 for bench in benches if bench.get("passed")),
            "schema": summary.get("schema"),
        }
    return {
        "schema": "aevum.earth_geomorphology_benchmarks.v1",
        "suite_count": len(suites),
        "all_passed": bool(suites) and all(
            item["status"] == "pass" for item in suites.values()),
        "suites": suites,
    }


def _summarize_entry(preset: str, eng, cm, checks, tectonics: dict[str, Any],
                     climate: dict[str, Any], compiler: dict[str, Any],
                     geomorphology: dict[str, Any],
                     earth_reference: dict[str, Any]) -> dict[str, Any]:
    w = eng.world
    morphology = tectonics["morphology"]
    ocean = tectonics["ocean_geography"]
    crust = tectonics["crust_distribution"]
    frames = tectonics["frame_continuity"]
    seam = tectonics["seam_continuity"]
    climate_ocean = climate["ocean_currents"]
    climate_geo = climate["geography"]

    validation_summary = {
        c.name: {
            "passed": bool(c.passed),
            "hard_failures": list(c.detail.get("hard_failures", [])),
        }
        for c in checks
    }

    land = morphology["exposed_land"]
    cont = morphology["continental_crust"]
    coupling = morphology["crust_land_coupling"]
    from aevum.diagnostics.p26_deforming_networks import deforming_network_summary
    from aevum.diagnostics.p26_rework_footprint import recent_rework_footprint_summary
    from aevum.diagnostics.p26_ribbon_drivers import ribbon_driver_summary
    from aevum.diagnostics.p29_inland_geomorphology import inland_geomorphology_summary

    return {
        "preset": preset,
        "spec_name": w.spec.name,
        "seed": int(w.spec.seed),
        "cells": int(w.grid.n),
        "time_myr": float(w.time_myr),
        "regime": w.spec.initial_tectonic_regime.value,
        "target_land_fraction": float(w.spec.target_land_fraction),
        "land_fraction": float(w.land_fraction()),
        "validation": validation_summary,
        "all_validation_checks_passed": bool(all(c.passed for c in checks)),
        "tectonic_warning_count": len(tectonics.get("warnings", [])),
        "tectonic_warnings": list(tectonics.get("warnings", [])),
        "climate_warning_count": len(climate.get("warnings", [])),
        "climate_warnings": list(climate.get("warnings", [])),
        "morphology": {
            "land_component_count": int(land["component_count"]),
            "largest_land_component_fraction": float(
                land["largest_component_area_fraction_of_mask"]),
            "land_width_p50_steps": float(land["width_p50_steps"]),
            "land_width_p90_steps": float(land["width_p90_steps"]),
            "land_ribbon_fraction_gt_0_5": float(land["ribbon_area_fraction_gt_0_5"]),
            "land_narrow_necks_per_1000": float(
                land["narrow_neck_cells_per_1000_mask_cells"]),
            "land_coastline_complexity_largest": float(
                land["coastline_complexity_largest_component"]),
            "continental_component_count": int(cont["component_count"]),
            "continental_width_p50_steps": float(cont["width_p50_steps"]),
            "continental_width_p90_steps": float(cont["width_p90_steps"]),
            "continental_ribbon_fraction_gt_0_5": float(
                cont["ribbon_area_fraction_gt_0_5"]),
            "continental_narrow_fraction_le2_width": float(cont["narrow_fraction_le2_width"]),
            "exposed_continental_land_fraction": float(
                coupling["exposed_continental_land_fraction_of_land"]),
            "exposed_oceanic_land_fraction": float(
                coupling["exposed_oceanic_land_fraction_of_land"]),
            "high_oceanic_land_fraction_gt1500m": float(
                coupling["high_oceanic_land_fraction_gt1500m_of_land"]),
        },
        "crust": {
            "continental_area_fraction": float(crust["continental_area_fraction"]),
            "stable_craton_fraction_gt075": float(crust["stable_craton_fraction_gt075"]),
            "ancient_continental_fraction_gt2500_myr": float(
                crust["ancient_continental_fraction_gt2500_myr"]),
            "oceanic_age_p95_myr": float(crust["oceanic_age_p95_myr"]),
            "oceanic_old_fraction_gt300_myr": float(crust["oceanic_old_fraction_gt300_myr"]),
        },
        "terrain_detail": _terrain_detail_summary(w, cm),
        "tectonic_object_telemetry": _tectonic_object_telemetry(w),
        "p26_ribbon_drivers": ribbon_driver_summary(w),
        "p26_rework_footprint": recent_rework_footprint_summary(w),
        "p26_deforming_networks": deforming_network_summary(w),
        "p29_inland_geomorphology": inland_geomorphology_summary(w),
        "earth_geomorphology_coverage": geomorphology,
        "earth_reference_distribution": earth_reference,
        "ocean_geography": {
            "basin_count": int(ocean["ocean_basin_count"]),
            "basin_object_count": int(ocean["ocean_basin_object_count"]),
            "margin_object_count": int(ocean["ocean_margin_object_count"]),
            "gateway_object_count": int(ocean["ocean_gateway_object_count"]),
            "shelf_fraction": float(ocean["shelf_fraction_of_ocean"]),
            "slope_rise_fraction": float(ocean["slope_rise_fraction_of_ocean"]),
            "abyss_fraction": float(ocean["abyss_fraction_of_ocean"]),
            "ridge_fraction": float(ocean["ridge_fraction_of_ocean"]),
            "trench_fraction": float(ocean["trench_fraction_of_ocean"]),
            "restricted_fraction": float(ocean["restricted_fraction_of_ocean"]),
            "nearshore_depth_p75_m": float(ocean["nearshore_depth_p75_m"]),
            "shelf_depth_p75_m": float(ocean["shelf_depth_p75_m"]),
            "abyss_depth_p50_m": float(ocean["abyss_depth_p50_m"]),
            "shelf_to_abyss_depth_delta_m": float(
                ocean["shelf_to_abyss_depth_delta_m"]),
            "nearshore_superdeep_fraction_gt2500m": float(
                ocean["nearshore_superdeep_fraction_gt2500m"]),
            "far_ocean_shallow_fraction_lt1500m": float(
                ocean["far_ocean_shallow_fraction_lt1500m"]),
            "trench_near_active_margin_fraction": float(
                ocean["trench_near_active_margin_fraction"]),
        },
        "archive_continuity": {
            "late_plate_crust_delta_per_100myr": float(
                frames["late_max_plate_crust_composition_delta_per_100myr"]),
            "late_exposed_land_change_fraction": float(
                frames["late_max_exposed_land_change_fraction"]),
            "late_crust_domain_change_fraction": float(
                frames["late_max_crust_domain_change_fraction"]),
            "late_terrain_province_change_fraction": float(
                frames["late_max_terrain_province_change_fraction"]),
            "late_wilson_phase_change_fraction": float(
                frames["late_max_wilson_phase_change_fraction"]),
            "plate_label_persistence_mean": float(frames["plate_label_persistence_mean"]),
        },
        "seam_continuity": {
            "render_duplicate_land_mismatch_fraction": float(
                seam["render_duplicate_land_mismatch_fraction"]),
            "render_duplicate_basin_mismatch_fraction": float(
                seam["render_duplicate_basin_mismatch_fraction"]),
            "render_duplicate_elevation_delta_p95_m": float(
                seam["render_duplicate_elevation_delta_p95_m"]),
            "seam_exposed_land_component_mismatch_fraction": float(
                seam["seam_exposed_land_component_mismatch_fraction"]),
            "seam_ocean_basin_mismatch_fraction": float(
                seam["seam_ocean_basin_mismatch_fraction"]),
            "edge_band_land_mismatch_fraction": float(seam["edge_band_land_mismatch_fraction"]),
        },
        "compiler": {
            "passed_envelope": bool(compiler["passed_envelope"]),
            "compiled_land_fraction": float(compiler["compiled_land_fraction"]),
            "compiled_land_or_coast_fraction": float(
                compiler["compiled_land_or_coast_fraction"]),
            "source_land_fraction_mean": float(compiler["source_land_fraction_mean"]),
            "land_fraction_abs_delta_from_source": float(
                compiler["land_fraction_abs_delta_from_source"]),
            "broad_land_to_water_fraction": float(compiler["broad_land_to_water_fraction"]),
            "broad_ocean_to_land_fraction": float(compiler["broad_ocean_to_land_fraction"]),
            "shelf_as_deep_ocean_fraction": float(compiler["shelf_as_deep_ocean_fraction"]),
            "terrain_elevation_sign_mismatch_fraction": float(
                compiler["terrain_elevation_sign_mismatch_fraction"]),
            "start_count": int(len(cm.starts)),
            "local_yield_cv": float(cm.fairness.get("local_yield_cv", 0.0)),
        },
        "climate_facing_prerequisites": {
            "has_basin_id": bool(climate_geo["has_basin_id"]),
            "has_continent_id": bool(climate_geo["has_continent_id"]),
            "has_shelf_index": bool(climate_geo["has_shelf_index"]),
            "has_barrier_index": bool(climate_geo["has_barrier_index"]),
            "ocean_heat_transport_abs_p95_C": float(
                climate_ocean["ocean_heat_transport_abs_p95_C"]),
            "solved_final_ocean_mask_mismatch_fraction": float(
                climate_ocean["solved_final_ocean_mask_mismatch_fraction"]),
        },
    }


def _tectonic_object_telemetry(world) -> dict[str, Any]:
    """Summarize process-object activity needed to interpret morphology gates."""
    breakup = list(world.objects.get("tectonics.breakup_seaways", []))
    rifts = list(world.objects.get("tectonics.rift_systems", []))
    components = list(world.objects.get("tectonics.breakup_component_telemetry", []))
    attempts = list(world.objects.get("terrain.breakup_seaway_attempts", []))
    stages = list(world.objects.get("terrain.land_component_stage_telemetry", []))
    coastline_stages = list(
        world.objects.get("terrain.modern_coastline_smoothing_stages", []))
    payback_filter_stages = list(
        world.objects.get("terrain.modern_coastline_payback_filter_stages", []))
    opened_corridors = list(
        world.objects.get("terrain.breakup_seaway_opened_corridors", []))

    def mean_metric(key: str) -> float:
        vals = [float(obj.get(key, 0.0)) for obj in breakup if key in obj]
        return float(np.mean(vals)) if vals else 0.0

    def max_metric(key: str) -> float:
        vals = [float(obj.get(key, 0.0)) for obj in breakup if key in obj]
        return float(np.max(vals)) if vals else 0.0

    skip_reasons: dict[str, int] = {}
    for comp in components:
        reason = str(comp.get("skip_reason", ""))
        if reason:
            skip_reasons[reason] = skip_reasons.get(reason, 0) + 1
    attempt_reject_reasons: dict[str, int] = {}
    for attempt in attempts:
        reason = str(attempt.get("reject_reason", ""))
        if reason:
            attempt_reject_reasons[reason] = attempt_reject_reasons.get(reason, 0) + 1

    top_components = sorted(
        components,
        key=lambda item: (
            -float(item.get("component_area_fraction", 0.0)),
            int(item.get("component_index", 0)),
        ),
    )[:8]

    return {
        "rift_system_object_count": int(len(rifts)),
        "rift_system_total_area_fraction": float(
            sum(float(obj.get("area_fraction", 0.0)) for obj in rifts)
        ),
        "rift_system_max_cell_count": int(max(
            [len(obj.get("cells", [])) for obj in rifts] or [0])),
        "breakup_seaway_object_count": int(len(breakup)),
        "breakup_seaway_total_area_fraction": float(
            sum(float(obj.get("area_fraction", 0.0)) for obj in breakup)
        ),
        "breakup_seaway_mean_topology_score": mean_metric("topology_score"),
        "breakup_seaway_max_topology_score": max_metric("topology_score"),
        "breakup_seaway_mean_split_balance": mean_metric("split_balance"),
        "breakup_seaway_mean_core_fraction": mean_metric("core_fraction"),
        "tectonics_continent_shape_pressure": float(
            world.g("tectonics.last_continent_shape_pressure", 0.0)),
        "tectonics_background_continent_shape_maintenance": float(
            world.g("tectonics.last_background_continent_shape_maintenance", 0.0)),
        "tectonics_preconserve_craton_promoted_cells": float(
            world.g("tectonics.last_preconserve_craton_promoted_cells", 0.0)),
        "tectonics_preconserve_craton_anchor_cells": float(
            world.g("tectonics.last_preconserve_craton_anchor_cells", 0.0)),
        "tectonics_passive_margin_progradation_cells": float(
            world.g("tectonics.last_passive_margin_progradation_cells", 0.0)),
        "tectonics_p33_parented_sliver_candidate_cells": float(
            world.g("tectonics.last_p33_parented_sliver_candidate_cells", 0.0)),
        "tectonics_p33_parented_sliver_consolidated_cells": float(
            world.g("tectonics.last_p33_parented_sliver_consolidated_cells", 0.0)),
        "tectonics_p33_parented_sliver_area_fraction": float(
            world.g("tectonics.last_p33_parented_sliver_area_fraction", 0.0)),
        "tectonics_p33_parented_sliver_active_preserved_cells": float(
            world.g("tectonics.last_p33_parented_sliver_active_preserved_cells", 0.0)),
        "tectonics_p33_parented_sliver_detached_preserved_cells": float(
            world.g("tectonics.last_p33_parented_sliver_detached_preserved_cells", 0.0)),
        "tectonics_p34_parented_platform_candidate_cells": float(
            world.g("tectonics.last_p34_parented_platform_candidate_cells", 0.0)),
        "tectonics_p34_parented_platform_restored_cells": float(
            world.g("tectonics.last_p34_parented_platform_restored_cells", 0.0)),
        "tectonics_p34_parented_platform_area_fraction": float(
            world.g("tectonics.last_p34_parented_platform_area_fraction", 0.0)),
        "tectonics_p34_parented_platform_protected_cells": float(
            world.g("tectonics.last_p34_parented_platform_protected_cells", 0.0)),
        "tectonics_p35_continental_debt_candidate_cells": float(
            world.g("tectonics.last_p35_continental_debt_candidate_cells", 0.0)),
        "tectonics_p35_continental_debt_restored_cells": float(
            world.g("tectonics.last_p35_continental_debt_restored_cells", 0.0)),
        "tectonics_p35_continental_debt_area_fraction": float(
            world.g("tectonics.last_p35_continental_debt_area_fraction", 0.0)),
        "tectonics_cumulative_p35_continental_debt_candidate_cells": float(
            world.g("tectonics.cumulative_p35_continental_debt_candidate_cells", 0.0)),
        "tectonics_cumulative_p35_continental_debt_restored_cells": float(
            world.g("tectonics.cumulative_p35_continental_debt_restored_cells", 0.0)),
        "tectonics_cumulative_p35_continental_debt_area_fraction": float(
            world.g("tectonics.cumulative_p35_continental_debt_area_fraction", 0.0)),
        "tectonics_p35_stable_parent_cells": float(
            world.g("tectonics.last_p35_stable_parent_cells", 0.0)),
        "tectonics_p35_old_platform_anchor_cells": float(
            world.g("tectonics.last_p35_old_platform_anchor_cells", 0.0)),
        "tectonics_p35_continental_gap_fraction": float(
            world.g("tectonics.last_p35_continental_gap_fraction", 0.0)),
        "tectonics_p37_shape_guard_input_cells": float(
            world.g("tectonics.last_p37_shape_guard_input_cells", 0.0)),
        "tectonics_p37_shape_guard_accepted_cells": float(
            world.g("tectonics.last_p37_shape_guard_accepted_cells", 0.0)),
        "tectonics_p37_shape_guard_rejected_cells": float(
            world.g("tectonics.last_p37_shape_guard_rejected_cells", 0.0)),
        "tectonics_cumulative_p37_shape_guard_input_cells": float(
            world.g("tectonics.cumulative_p37_shape_guard_input_cells", 0.0)),
        "tectonics_cumulative_p37_shape_guard_accepted_cells": float(
            world.g("tectonics.cumulative_p37_shape_guard_accepted_cells", 0.0)),
        "tectonics_cumulative_p37_shape_guard_rejected_cells": float(
            world.g("tectonics.cumulative_p37_shape_guard_rejected_cells", 0.0)),
        "tectonics_p38_lineage_guard_input_cells": float(
            world.g("tectonics.last_p38_lineage_guard_input_cells", 0.0)),
        "tectonics_p38_lineage_guard_accepted_cells": float(
            world.g("tectonics.last_p38_lineage_guard_accepted_cells", 0.0)),
        "tectonics_p38_lineage_guard_rejected_cells": float(
            world.g("tectonics.last_p38_lineage_guard_rejected_cells", 0.0)),
        "tectonics_cumulative_p38_lineage_guard_input_cells": float(
            world.g("tectonics.cumulative_p38_lineage_guard_input_cells", 0.0)),
        "tectonics_cumulative_p38_lineage_guard_accepted_cells": float(
            world.g("tectonics.cumulative_p38_lineage_guard_accepted_cells", 0.0)),
        "tectonics_cumulative_p38_lineage_guard_rejected_cells": float(
            world.g("tectonics.cumulative_p38_lineage_guard_rejected_cells", 0.0)),
        **{
            f"tectonics_p39_old_platform_funnel_{key}": float(
                world.g(f"tectonics.last_p39_old_platform_funnel_{key}", 0.0)
            )
            for key in P39_OLD_PLATFORM_FUNNEL_METRICS
        },
        "tectonics_p40_parented_restoration_inherited_rework_cells": float(
            world.g("tectonics.last_p40_parented_restoration_inherited_rework_cells", 0.0)),
        "tectonics_p40_parented_restoration_reactivated_cells": float(
            world.g("tectonics.last_p40_parented_restoration_reactivated_cells", 0.0)),
        "tectonics_p40_shape_maintenance_inherited_rework_cells": float(
            world.g("tectonics.last_p40_shape_maintenance_inherited_rework_cells", 0.0)),
        "tectonics_p40_shape_maintenance_reactivated_cells": float(
            world.g("tectonics.last_p40_shape_maintenance_reactivated_cells", 0.0)),
        "tectonics_p41_parented_restoration_priority_cells": float(
            world.g("tectonics.last_p41_parented_restoration_priority_cells", 0.0)),
        "tectonics_p41_active_only_accretion_deferred_cells": float(
            world.g("tectonics.last_p41_active_only_accretion_deferred_cells", 0.0)),
        "tectonics_p43_mature_active_only_candidate_cells": float(
            world.g("tectonics.last_p43_mature_active_only_candidate_cells", 0.0)),
        "tectonics_p43_mature_active_only_accepted_cells": float(
            world.g("tectonics.last_p43_mature_active_only_accepted_cells", 0.0)),
        "tectonics_p43_mature_active_only_rejected_cells": float(
            world.g("tectonics.last_p43_mature_active_only_rejected_cells", 0.0)),
        "tectonics_p43_mature_active_only_accepted_area_fraction": float(
            world.g("tectonics.last_p43_mature_active_only_accepted_area_fraction", 0.0)),
        "tectonics_p43_mature_active_only_cap_area_fraction": float(
            world.g("tectonics.last_p43_mature_active_only_cap_area_fraction", 0.0)),
        "tectonics_p44_parented_recovery_candidate_cells": float(
            world.g("tectonics.last_p44_parented_recovery_candidate_cells", 0.0)),
        "tectonics_p44_parented_recovery_restored_cells": float(
            world.g("tectonics.last_p44_parented_recovery_restored_cells", 0.0)),
        "tectonics_p44_parented_recovery_area_fraction": float(
            world.g("tectonics.last_p44_parented_recovery_area_fraction", 0.0)),
        "tectonics_p44_parented_recovery_protected_cells": float(
            world.g("tectonics.last_p44_parented_recovery_protected_cells", 0.0)),
        "tectonics_p45_collision_process_cells": float(
            world.g("tectonics.last_p45_collision_process_cells", 0.0)),
        "tectonics_p45_collision_provenance_cells": float(
            world.g("tectonics.last_p45_collision_provenance_cells", 0.0)),
        "tectonics_p45_arc_process_cells": float(
            world.g("tectonics.last_p45_arc_process_cells", 0.0)),
        "tectonics_p45_arc_provenance_cells": float(
            world.g("tectonics.last_p45_arc_provenance_cells", 0.0)),
        "tectonics_p46_mature_collage_candidate_cells": float(
            world.g("tectonics.last_p46_mature_collage_candidate_cells", 0.0)),
        "tectonics_p46_mature_collage_consolidated_cells": float(
            world.g("tectonics.last_p46_mature_collage_consolidated_cells", 0.0)),
        "tectonics_p46_mature_collage_area_fraction": float(
            world.g("tectonics.last_p46_mature_collage_area_fraction", 0.0)),
        "tectonics_p46_mature_collage_active_preserved_cells": float(
            world.g("tectonics.last_p46_mature_collage_active_preserved_cells", 0.0)),
        "tectonics_p46_mature_collage_detached_preserved_cells": float(
            world.g("tectonics.last_p46_mature_collage_detached_preserved_cells", 0.0)),
        "tectonics_p47_craton_target_fraction": float(
            world.g("tectonics.last_p47_craton_target_fraction", 0.0)),
        "tectonics_p47_pre_stable_craton_fraction": float(
            world.g("tectonics.last_p47_pre_stable_craton_fraction", 0.0)),
        "tectonics_p47_promoted_craton_cells": float(
            world.g("tectonics.last_p47_promoted_craton_cells", 0.0)),
        "tectonics_p47_promoted_craton_area_fraction": float(
            world.g("tectonics.last_p47_promoted_craton_area_fraction", 0.0)),
        "tectonics_p48_shield_object_count": float(
            world.g("tectonics.last_p48_shield_object_count", 0.0)),
        "tectonics_p48_shield_share_of_continental_crust": float(
            world.g("tectonics.last_p48_shield_share_of_continental_crust", 0.0)),
        "tectonics_p48_platform_object_count": float(
            world.g("tectonics.last_p48_platform_object_count", 0.0)),
        "tectonics_p48_platform_share_of_continental_crust": float(
            world.g("tectonics.last_p48_platform_share_of_continental_crust", 0.0)),
        "tectonics_p48_interior_basin_object_count": float(
            world.g("tectonics.last_p48_interior_basin_object_count", 0.0)),
        "tectonics_p48_interior_basin_share_of_continental_crust": float(
            world.g("tectonics.last_p48_interior_basin_share_of_continental_crust", 0.0)),
        "tectonics_p49_platform_subsidence_area_fraction": float(
            world.g("tectonics.last_p49_platform_subsidence_area_fraction", 0.0)),
        "tectonics_p49_platform_subsidence_mean": float(
            world.g("tectonics.last_p49_platform_subsidence_mean", 0.0)),
        "tectonics_p49_platform_subsidence_max": float(
            world.g("tectonics.last_p49_platform_subsidence_max", 0.0)),
        "tectonics_p50_planform_recycle_candidate_cells": float(
            world.g("tectonics.last_p50_planform_recycle_candidate_cells", 0.0)),
        "tectonics_p50_planform_fill_candidate_cells": float(
            world.g("tectonics.last_p50_planform_fill_candidate_cells", 0.0)),
        "tectonics_p50_planform_recycled_cells": float(
            world.g("tectonics.last_p50_planform_recycled_cells", 0.0)),
        "tectonics_p50_planform_filled_cells": float(
            world.g("tectonics.last_p50_planform_filled_cells", 0.0)),
        "tectonics_p50_planform_area_fraction": float(
            world.g("tectonics.last_p50_planform_area_fraction", 0.0)),
        "tectonics_p50_planform_before_narrow_fraction": float(
            world.g("tectonics.last_p50_planform_before_narrow_fraction", 0.0)),
        "tectonics_p50_planform_after_narrow_fraction": float(
            world.g("tectonics.last_p50_planform_after_narrow_fraction", 0.0)),
        "tectonics_p50_planform_narrow_fraction_delta": float(
            world.g("tectonics.last_p50_planform_narrow_fraction_delta", 0.0)),
        "tectonics_cumulative_p50_planform_recycled_cells": float(
            world.g("tectonics.cumulative_p50_planform_recycled_cells", 0.0)),
        "tectonics_cumulative_p50_planform_filled_cells": float(
            world.g("tectonics.cumulative_p50_planform_filled_cells", 0.0)),
        "tectonics_cumulative_p50_planform_area_fraction": float(
            world.g("tectonics.cumulative_p50_planform_area_fraction", 0.0)),
        "tectonics_continent_gain_cells": float(
            world.g("tectonics.last_continent_gain_cells", 0.0)),
        "tectonics_continent_loss_cells": float(
            world.g("tectonics.last_continent_loss_cells", 0.0)),
        "tectonics_unforced_continent_gain_blocked": float(
            world.g("tectonics.last_unforced_continent_gain_blocked", 0.0)),
        "tectonics_unforced_continent_loss_blocked": float(
            world.g("tectonics.last_unforced_continent_loss_blocked", 0.0)),
        "terrain_breakup_seaway_openings": float(
            world.g("terrain.last_breakup_seaway_openings", 0.0)),
        "terrain_breakup_seaway_area_fraction": float(
            world.g("terrain.last_breakup_seaway_area_fraction", 0.0)),
        "terrain_breakup_seaway_source_reuse": float(
            world.g("terrain.last_breakup_seaway_source_reuse", 0.0)),
        "terrain_breakup_seaway_opened_corridor_count": int(len(opened_corridors)),
        "terrain_breakup_seaway_opened_corridor_area_fraction": float(sum(
            float(corridor.get("area_fraction", 0.0))
            for corridor in opened_corridors
        )),
        "terrain_breakup_seaway_attempt_count": int(len(attempts)),
        "terrain_breakup_seaway_applied_attempt_count": int(sum(
            1 for attempt in attempts if bool(attempt.get("applied", False)))),
        "terrain_breakup_seaway_reject_reasons": attempt_reject_reasons,
        "terrain_breakup_seaway_top_attempts": [
            {
                "object_id": str(attempt.get("object_id", "")),
                "largest_share_before": float(attempt.get("largest_share_before", 0.0)),
                "source_land_fraction": float(attempt.get("source_land_fraction", 0.0)),
                "reused_process_source": bool(
                    attempt.get("reused_process_source", False)),
                "candidate_count": int(attempt.get("candidate_count", 0)),
                "viable_candidate_count": int(
                    attempt.get("viable_candidate_count", 0)),
                "best_reduction": float(attempt.get("best_reduction", 0.0)),
                "best_new_largest_share": float(
                    attempt.get("best_new_largest_share", 0.0)),
                "best_cut_area_fraction": float(
                    attempt.get("best_cut_area_fraction", 0.0)),
                "applied": bool(attempt.get("applied", False)),
                "reject_reason": str(attempt.get("reject_reason", "")),
            }
            for attempt in sorted(
                attempts,
                key=lambda item: (
                    -float(item.get("best_reduction", 0.0)),
                    not bool(item.get("applied", False)),
                    str(item.get("object_id", "")),
                ),
            )[:8]
        ],
        "terrain_land_component_stages": [
            {
                "stage": str(stage.get("stage", "")),
                "land_fraction": float(stage.get("land_fraction", 0.0)),
                "land_component_count": int(stage.get("land_component_count", 0)),
                "largest_land_component_fraction": float(
                    stage.get("largest_land_component_fraction", 0.0)),
            }
            for stage in stages[:16]
        ],
        "terrain_modern_coastline_smoothing_stages": [
            {
                "stage": str(stage.get("stage", "")),
                "land_fraction": float(stage.get("land_fraction", 0.0)),
                "land_component_count": int(stage.get("land_component_count", 0)),
                "largest_land_component_fraction": float(
                    stage.get("largest_land_component_fraction", 0.0)),
            }
            for stage in coastline_stages[:12]
        ],
        "terrain_modern_coastline_payback_filter_stages": [
            {
                "stage": str(stage.get("stage", "")),
                "recovery_need_fraction": float(
                    stage.get("recovery_need_fraction", 0.0)),
                "base_ocean_cont_area_fraction": float(
                    stage.get("base_ocean_cont_area_fraction", 0.0)),
                "same_component_area_fraction": float(
                    stage.get("same_component_area_fraction", 0.0)),
                "broad_supported_area_fraction": float(
                    stage.get("broad_supported_area_fraction", 0.0)),
                "support_ok_area_fraction": float(
                    stage.get("support_ok_area_fraction", 0.0)),
                "domain_ok_area_fraction": float(
                    stage.get("domain_ok_area_fraction", 0.0)),
                "geometry_ok_area_fraction": float(
                    stage.get("geometry_ok_area_fraction", 0.0)),
                "candidate_area_fraction": float(
                    stage.get("candidate_area_fraction", 0.0)),
                "base_ocean_cont_cells": int(stage.get("base_ocean_cont_cells", 0)),
                "same_component_cells": int(stage.get("same_component_cells", 0)),
                "broad_supported_cells": int(stage.get("broad_supported_cells", 0)),
                "support_ok_cells": int(stage.get("support_ok_cells", 0)),
                "domain_ok_cells": int(stage.get("domain_ok_cells", 0)),
                "geometry_ok_cells": int(stage.get("geometry_ok_cells", 0)),
                "candidate_cells": int(stage.get("candidate_cells", 0)),
            }
            for stage in payback_filter_stages[:4]
        ],
        "terrain_largest_landmass_shave_cells": float(
            world.g("terrain.last_largest_landmass_shave_cells", 0.0)),
        "terrain_largest_landmass_shave_area_fraction": float(
            world.g("terrain.last_largest_landmass_shave_area_fraction", 0.0)),
        "terrain_modern_coastline_drowned_area_fraction": float(
            world.g("terrain.last_modern_coastline_drowned_area_fraction", 0.0)),
        "terrain_modern_coastline_fill_area_fraction": float(
            world.g("terrain.last_modern_coastline_fill_area_fraction", 0.0)),
        "terrain_modern_coastline_payback_candidate_area_fraction": float(
            world.g(
                "terrain.last_modern_coastline_payback_candidate_area_fraction",
                0.0,
            )),
        "terrain_modern_coastline_payback_recovery_need_fraction": float(
            world.g(
                "terrain.last_modern_coastline_payback_recovery_need_fraction",
                0.0,
            )),
        "terrain_modern_coastline_payback_area_fraction": float(
            world.g("terrain.last_modern_coastline_payback_area_fraction", 0.0)),
        "terrain_modern_coastline_payback_cells": float(
            world.g("terrain.last_modern_coastline_payback_cells", 0.0)),
        "terrain_p43_payback_shape_rejected_cells": float(
            world.g("terrain.last_p43_payback_shape_rejected_cells", 0.0)),
        "terrain_p32_coastline_complexity_before": float(
            world.g("terrain.last_p32_coastline_complexity_before", 0.0)),
        "terrain_p32_coastline_complexity_after": float(
            world.g("terrain.last_p32_coastline_complexity_after", 0.0)),
        "terrain_p32_coastline_swap_area_fraction": float(
            world.g("terrain.last_p32_coastline_swap_area_fraction", 0.0)),
        "terrain_p32_coastline_drown_candidate_area_fraction": float(
            world.g("terrain.last_p32_coastline_drown_candidate_area_fraction", 0.0)),
        "terrain_p32_coastline_fill_candidate_area_fraction": float(
            world.g("terrain.last_p32_coastline_fill_candidate_area_fraction", 0.0)),
        "terrain_p32_coastline_drowned_cells": float(
            world.g("terrain.last_p32_coastline_drowned_cells", 0.0)),
        "terrain_p32_coastline_filled_cells": float(
            world.g("terrain.last_p32_coastline_filled_cells", 0.0)),
        "terrain_deformation_relief_area_fraction": float(
            world.g("terrain.last_deformation_relief_area_fraction", 0.0)),
        "terrain_deformation_relief_mean_m": float(
            world.g("terrain.last_deformation_relief_mean_m", 0.0)),
        "terrain_inland_geomorphology_relief_area_fraction": float(
            world.g("terrain.last_inland_geomorphology_relief_area_fraction", 0.0)),
        "terrain_inland_geomorphology_relief_uplift_area_fraction": float(
            world.g("terrain.last_inland_geomorphology_relief_uplift_area_fraction", 0.0)),
        "terrain_inland_geomorphology_relief_depression_area_fraction": float(
            world.g("terrain.last_inland_geomorphology_relief_depression_area_fraction", 0.0)),
        "terrain_inland_geomorphology_relief_mean_abs_m": float(
            world.g("terrain.last_inland_geomorphology_relief_mean_abs_m", 0.0)),
        "terrain_inland_state_area_fraction": float(
            world.g("terrain.last_inland_state_area_fraction", 0.0)),
        "terrain_inland_state_relief_mean_abs_m": float(
            world.g("terrain.last_inland_state_relief_mean_abs_m", 0.0)),
        "terrain_p49_subsidence_basin_sediment_mean_m": float(
            world.g("terrain.last_p49_subsidence_basin_sediment_mean_m", 0.0)),
        "terrain_p49_subsidence_platform_sediment_mean_m": float(
            world.g("terrain.last_p49_subsidence_platform_sediment_mean_m", 0.0)),
        "terrain_p49_subsidence_sediment_contrast_m": float(
            world.g("terrain.last_p49_subsidence_sediment_contrast_m", 0.0)),
        "breakup_component_evaluated_count": int(len(components)),
        "breakup_component_eligible_count": int(sum(
            1 for comp in components if bool(comp.get("eligible", False)))),
        "breakup_component_candidate_total": int(sum(
            int(comp.get("candidate_count", 0)) for comp in components)),
        "breakup_component_accepted_total": int(sum(
            int(comp.get("accepted_object_count", 0)) for comp in components)),
        "breakup_component_skip_reasons": skip_reasons,
        "breakup_component_top": [
            {
                "component_index": int(comp.get("component_index", 0)),
                "component_area_fraction": float(comp.get("component_area_fraction", 0.0)),
                "component_share_of_continental_crust": float(
                    comp.get("component_share_of_continental_crust", 0.0)),
                "boundary_hit_count": int(comp.get("boundary_hit_count", 0)),
                "boundary_hit_area_fraction": float(
                    comp.get("boundary_hit_area_fraction", 0.0)),
                "weak_area_fraction": float(comp.get("weak_area_fraction", 0.0)),
                "candidate_count": int(comp.get("candidate_count", 0)),
                "accepted_object_count": int(comp.get("accepted_object_count", 0)),
                "best_topology_score": float(comp.get("best_topology_score", 0.0)),
                "eligible": bool(comp.get("eligible", False)),
                "skip_reason": str(comp.get("skip_reason", "")),
            }
            for comp in top_components
        ],
    }


def _terrain_detail_summary(world, cm) -> dict[str, Any]:
    elev = np.asarray(world.get_field("terrain.elevation_m"), dtype=np.float64)
    rel = elev - float(world.sea_level)
    land = rel >= 0.0
    detail = np.asarray(
        world.get_field("terrain.continental_detail", np.zeros(world.grid.n)),
        dtype=int,
    )
    terrain = np.asarray(
        world.get_field("terrain.province", np.zeros(world.grid.n)),
        dtype=int,
    )

    if land.any():
        land_rel = rel[land]
        detail_land = detail[land]
        terrain_land = terrain[land]
        elevation = {
            "land_elevation_mean_m": float(np.mean(land_rel)),
            "land_elevation_p75_m": float(np.percentile(land_rel, 75)),
            "land_elevation_p90_m": float(np.percentile(land_rel, 90)),
            "land_elevation_p95_m": float(np.percentile(land_rel, 95)),
            "land_elevation_p99_m": float(np.percentile(land_rel, 99)),
            "land_elevation_max_m": float(np.max(land_rel)),
        }
        detail_frac = {
            "shield_fraction_of_land": float(np.mean(detail_land == 1)),
            "platform_fraction_of_land": float(np.mean(detail_land == 2)),
            "basin_or_rift_fraction_of_land": float(np.mean(np.isin(detail_land, [3, 4]))),
            "orogen_or_plateau_fraction_of_land": float(np.mean(np.isin(detail_land, [5, 6]))),
            "arc_microcontinent_fraction_of_land": float(np.mean(detail_land == 7)),
        }
        terrain_frac = {
            "suture_lip_highland_fraction_of_land": float(
                np.mean(np.isin(terrain_land, [5, 6, 7]))
            ),
        }
    else:
        elevation = {
            "land_elevation_mean_m": 0.0,
            "land_elevation_p75_m": 0.0,
            "land_elevation_p90_m": 0.0,
            "land_elevation_p95_m": 0.0,
            "land_elevation_p99_m": 0.0,
            "land_elevation_max_m": 0.0,
        }
        detail_frac = {
            "shield_fraction_of_land": 0.0,
            "platform_fraction_of_land": 0.0,
            "basin_or_rift_fraction_of_land": 0.0,
            "orogen_or_plateau_fraction_of_land": 0.0,
            "arc_microcontinent_fraction_of_land": 0.0,
        }
        terrain_frac = {"suture_lip_highland_fraction_of_land": 0.0}

    compiled_land = (cm.terrain >= 2) & (cm.terrain <= 4)
    compiled = {
        "compiled_hills_fraction": float(np.mean(cm.terrain == 3)),
        "compiled_mountain_fraction": float(np.mean(cm.terrain == 4)),
        "compiled_hills_or_mountains_fraction_of_land": 0.0,
    }
    if compiled_land.any():
        compiled["compiled_hills_or_mountains_fraction_of_land"] = float(
            np.mean(np.isin(cm.terrain[compiled_land], [3, 4]))
        )
    out = {}
    out.update(elevation)
    out.update(detail_frac)
    out.update(terrain_frac)
    out.update(compiled)
    return out


def _assess_entry(entry: dict[str, Any]) -> dict[str, Any]:
    failures: list[str] = []
    warnings: list[str] = []
    preset = entry["preset"]
    low_resolution = int(entry.get("cells", 0)) < 4000
    m = entry["morphology"]
    ocean = entry["ocean_geography"]
    crust = entry["crust"]
    terrain = entry.get("terrain_detail", {})
    seam = entry["seam_continuity"]
    compiler = entry["compiler"]
    archive = entry["archive_continuity"]
    climate = entry["climate_facing_prerequisites"]
    geomorphology = entry.get("earth_geomorphology_coverage", {})
    earth_reference = entry.get("earth_reference_distribution", {})

    if not entry["all_validation_checks_passed"]:
        failed = [name for name, result in entry["validation"].items() if not result["passed"]]
        failures.append("validation checks failed: " + ", ".join(failed))
    if not compiler["passed_envelope"]:
        failures.append("compiled map contradicts source terrain or ocean provinces")
    if seam["render_duplicate_land_mismatch_fraction"] > 0.0:
        failures.append("rendered antimeridian duplicate columns disagree on land/sea")
    if seam["render_duplicate_basin_mismatch_fraction"] > 0.0:
        failures.append("rendered antimeridian duplicate columns disagree on ocean basin id")
    if seam["render_duplicate_elevation_delta_p95_m"] > 0.0:
        failures.append("rendered antimeridian duplicate columns disagree on elevation")
    if ocean["nearshore_superdeep_fraction_gt2500m"] > 0.08:
        failures.append("too much superdeep water immediately offshore")
    if ocean["far_ocean_shallow_fraction_lt1500m"] > 0.20:
        if low_resolution:
            warnings.append("low-resolution ocean smoke: far ocean is too shallow")
        else:
            failures.append("far ocean is too shallow")
    if ocean["shelf_to_abyss_depth_delta_m"] < 1000.0 and entry["land_fraction"] < 0.95:
        failures.append("shelves and abyssal plains are not depth-separated")
    if archive["late_plate_crust_delta_per_100myr"] > 0.12:
        failures.append("late archive continuity has large plate/crust jumps")
    if not all(climate[k] for k in ("has_basin_id", "has_continent_id",
                                    "has_shelf_index", "has_barrier_index")):
        failures.append("climate-facing geography primitives are incomplete")

    if preset == "earthlike":
        def morphology_failure(message: str) -> None:
            if low_resolution:
                warnings.append("low-resolution morphology smoke: " + message)
            else:
                failures.append(message)

        land = entry["land_fraction"]
        if not (0.20 <= land <= 0.38):
            failures.append("Earthlike land fraction outside broad Earth analogue envelope")
        if m["largest_land_component_fraction"] < 0.20:
            morphology_failure("Earthlike land is too fragmented to support broad continents")
        if m["largest_land_component_fraction"] > 0.82:
            morphology_failure("Earthlike land is effectively a single supercontinent")
        elif m["largest_land_component_fraction"] > 0.65:
            warnings.append(
                "Earthlike land is supercontinent-like rather than modern-fragmented")
        if m["land_ribbon_fraction_gt_0_5"] > 0.55:
            morphology_failure("Earthlike land is dominated by ribbon-like landforms")
        if m["continental_narrow_fraction_le2_width"] > 0.78:
            morphology_failure("Earthlike continental crust has too little broad interior")
        if m["high_oceanic_land_fraction_gt1500m"] > 0.005:
            failures.append("high exposed oceanic crust is too common")
        if ocean["basin_count"] < 2:
            failures.append("Earthlike ocean lacks multiple basins")
        if crust["stable_craton_fraction_gt075"] < 0.01:
            failures.append("Earthlike world lacks preserved stable cratonic crust")
        if (entry.get("time_myr", 0.0) >= 3000.0
                and crust["ancient_continental_fraction_gt2500_myr"] < 0.05):
            warnings.append("ancient continental crust fraction is low for an Earth analogue")
        if m["land_component_count"] > 14:
            warnings.append("Earthlike exposed land has more island/land components than modern Earth")
        if m["land_ribbon_fraction_gt_0_5"] > 0.08:
            warnings.append("Earthlike land still contains more narrow/ribbon area than desired")
        if m["land_narrow_necks_per_1000"] > 15.0:
            warnings.append("Earthlike land still has many graph articulation necks")
        if m["land_coastline_complexity_largest"] > 8.0:
            warnings.append("largest Earthlike landmass coastline is over-complex")
        if terrain.get("land_elevation_p95_m", 0.0) < 1500.0:
            warnings.append("Earthlike continental hypsometry still has a weak high tail")
        if terrain.get("orogen_or_plateau_fraction_of_land", 0.0) < 0.02:
            warnings.append("Earthlike continental detail lacks orogen/plateau provinces")
        if terrain.get("orogen_or_plateau_fraction_of_land", 0.0) > 0.45:
            warnings.append("Earthlike continental detail overpaints orogen/plateau provinces")
        if climate["ocean_heat_transport_abs_p95_C"] < 0.35:
            warnings.append("ocean heat-transport prerequisite is weak for later climate C4")
        if geomorphology:
            hard_failures = geomorphology.get("hard_failures", [])
            if hard_failures:
                message = (
                    "Earth geomorphology coverage failures: "
                    + "; ".join(str(item) for item in hard_failures[:3])
                )
                if low_resolution:
                    warnings.append("low-resolution geomorphology smoke: " + message)
                else:
                    failures.append(message)
            warnings.extend(
                f"geomorphology: {warning}"
                for warning in geomorphology.get("warnings", [])
            )
        if earth_reference and earth_reference.get("out_of_envelope"):
            warnings.append(
                "earth-reference distribution needs calibration: "
                + ", ".join(str(x) for x in earth_reference["out_of_envelope"][:5])
            )
    elif preset == "waterworld":
        if not (0.005 <= entry["land_fraction"] <= 0.10):
            failures.append("waterworld should retain only sparse emergent land")
    elif preset == "arid":
        if entry["land_fraction"] < 0.55:
            failures.append("arid preset should remain land-dominated")
    elif preset == "stagnant_lid":
        if entry["regime"] != "stagnant_lid":
            failures.append("stagnant-lid preset reported the wrong tectonic regime")
    elif preset == "frozen":
        if entry["land_fraction"] < 0.15:
            failures.append("frozen preset lost most emergent land")

    status = "pass" if not failures else "fail"
    if status == "pass" and warnings:
        status = "warn"
    return {
        "status": status,
        "passed": not failures,
        "failures": failures,
        "warnings": warnings,
    }


def _release_decision(entries: Iterable[dict[str, Any]]) -> dict[str, Any]:
    entries = list(entries)
    failed = [
        {
            "preset": e["preset"],
            "failures": e["release_gate"]["failures"],
        }
        for e in entries
        if not e["release_gate"]["passed"]
    ]
    warned = [
        {
            "preset": e["preset"],
            "warnings": e["release_gate"]["warnings"],
        }
        for e in entries
        if e["release_gate"]["passed"] and e["release_gate"]["warnings"]
    ]
    return {
        "passed": len(failed) == 0,
        "status": "pass" if not failed and not warned else ("warn" if not failed else "fail"),
        "entry_count": len(entries),
        "failed_entries": failed,
        "warned_entries": warned,
    }


def _render_contact_sheet(panels: list[dict[str, Any]], path: Path) -> Path:
    from aevum import render
    from aevum.diagnostics.morphology import ensure_morphology_fields

    if not panels:
        return path
    n = len(panels)
    fig, axes = plt.subplots(n, 4, figsize=(13, max(2.2, 2.0 * n)),
                             squeeze=False, constrained_layout=True)
    for row, item in enumerate(panels):
        w = item["world"]
        cm = item["compiled"]
        ensure_morphology_fields(w)
        sea = w.sea_level
        elev = render.to_raster_continuous(
            w.grid,
            w.get_field("terrain.elevation_m") - sea,
            width=180,
            height=90,
            preserve_sign=True,
        )
        depth = render.to_raster(w.grid, w.get_field("ocean.depth_province", 0.0),
                                 width=180, height=90)
        ribbon = render.to_raster(w.grid, w.get_field("tectonics.ribbon_index", 0.0),
                                  width=180, height=90)
        terrain = cm.terrain
        axes[row, 0].imshow(elev, cmap=render.ELEVATION_CMAP, norm=render.ELEVATION_NORM,
                            extent=[-180, 180, -90, 90])
        axes[row, 1].imshow(depth, cmap="viridis", vmin=0, vmax=7,
                            extent=[-180, 180, -90, 90])
        axes[row, 2].imshow(ribbon, cmap="magma", vmin=0, vmax=1,
                            extent=[-180, 180, -90, 90])
        axes[row, 3].imshow(terrain, cmap=render.ListedColormap(render.TERRAIN_COLORS),
                            vmin=0, vmax=5, interpolation="nearest")
        axes[row, 0].set_ylabel(item["preset"])
        for col, title in enumerate(("elevation", "ocean depth", "ribbon index", "hex terrain")):
            if row == 0:
                axes[row, col].set_title(title)
        for ax in axes[row]:
            ax.set_xticks([])
            ax.set_yticks([])
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return path


def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)
