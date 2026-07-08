import csv
import json

from aevum.diagnostics.planform_reference import planform_reference_summary
from aevum.diagnostics.tectonics_bench import run_suite


def test_r0_randomness_inventory_is_classified_and_deterministic(tmp_path):
    summary = run_suite("R0", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.r0.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]

    inventory = summary["randomness_inventory"]
    assert inventory["unclassified_rng_calls"] == 0
    assert inventory["forbidden_rng_calls"] > 0
    assert inventory["temporary_rng_calls"] > 0
    assert inventory["allowed_rng_calls"] > 0
    assert len(inventory["uses"]) == sum(inventory["counts"].values())

    forbidden_purposes = {
        use["purpose"]
        for use in inventory["uses"]
        if use["category"] == "forbidden"
    }
    assert any("nuclei" in purpose for purpose in forbidden_purposes)
    assert any("Euler poles" in purpose for purpose in forbidden_purposes)
    assert not any("plate reorganization" in purpose for purpose in forbidden_purposes)
    assert not any("plume locations" in purpose for purpose in forbidden_purposes)

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "randomness_inventory.csv"
    assert summary_path.exists()
    assert csv_path.exists()

    written = json.loads(summary_path.read_text())
    assert written["determinism"]["first_digest"] == summary["determinism"]["first_digest"]
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(inventory["uses"])


def test_r1_deep_interior_potential_microbenchmarks_pass(tmp_path):
    summary = run_suite("R1", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.r1.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "R1.heat_diffusion_stability",
        "R1.continental_insulation",
        "R1.slab_downwelling",
        "R1.plume_trigger",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    assert (
        benchmarks["R1.heat_diffusion_stability"]["metrics"]["roughness_final"]
        < benchmarks["R1.heat_diffusion_stability"]["metrics"]["roughness_initial"]
    )
    assert (
        benchmarks["R1.continental_insulation"]["metrics"]["interior_minus_ocean_heat"]
        > 0.025
    )
    assert (
        benchmarks["R1.slab_downwelling"]["metrics"]["trench_to_far_ratio"]
        > 2.0
    )
    assert (
        benchmarks["R1.plume_trigger"]["metrics"]["plume_max_distance_deg"]
        <= 10.0
    )

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "r1_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()

    written = json.loads(summary_path.read_text())
    assert written["determinism"]["summary_digest"] == summary["determinism"]["summary_digest"]
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_r2_torque_proxy_microbenchmarks_pass(tmp_path):
    summary = run_suite("R2", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.r2.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["random_reorg_jitter_removed"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "R2.single_slab_pull",
        "R2.ridge_push_symmetry",
        "R2.collision_locking",
        "R2.no_random_reorg_jitter",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    slab = benchmarks["R2.single_slab_pull"]["metrics"]
    assert slab["velocity_dot_toward_trench"] > 0.65
    assert slab["old_rate"] > slab["young_rate"]
    assert slab["old_rate"] > slab["short_trench_rate"]

    ridge = benchmarks["R2.ridge_push_symmetry"]["metrics"]
    assert ridge["west_velocity_dot_away_from_ridge"] > 0.55
    assert ridge["east_velocity_dot_away_from_ridge"] > 0.55
    assert ridge["net_torque_ratio"] < 0.22

    collision = benchmarks["R2.collision_locking"]["metrics"]
    assert collision["thick_collision_rate"] < collision["subduction_rate"]
    assert collision["thick_collision_resistance"] > collision["thin_collision_resistance"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "r2_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_r3_persistent_boundary_microbenchmarks_pass(tmp_path):
    summary = run_suite("R3", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.r3.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["persistent_boundary_objects_available"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "R3.subduction_polarity_old_ocean",
        "R3.ocean_continent_margin",
        "R3.transform_from_offset_ridge",
        "R3.boundary_persistence",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    polarity = benchmarks["R3.subduction_polarity_old_ocean"]["metrics"]
    assert polarity["west_old_subducting_plate"] == 0
    assert polarity["east_old_subducting_plate"] == 1

    margin = benchmarks["R3.ocean_continent_margin"]["metrics"]
    assert margin["trench_subducting_plate"] == 0
    assert margin["trench_overriding_plate"] == 1
    assert margin["trench_continental_fraction"] < 0.20
    assert margin["active_margin_continental_fraction"] > 0.80

    persistence = benchmarks["R3.boundary_persistence"]["metrics"]
    assert persistence["first_id"] == persistence["second_id"]
    assert persistence["age_increment_myr"] == 25.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "r3_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_r4_ocean_basin_lifecycle_microbenchmarks_pass(tmp_path):
    summary = run_suite("R4", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.r4.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["persistent_ocean_basin_lifecycle_available"]
    assert summary["acceptance"]["gateway_causality_available"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "R4.rift_to_ocean",
        "R4.basin_maturation",
        "R4.closure_to_suture",
        "R4.gateway_causality",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    rift = benchmarks["R4.rift_to_ocean"]["metrics"]
    assert rift["initial_basin_id"] == rift["later_basin_id"]
    assert rift["basin_parent_margin_count"] >= 1

    closure = benchmarks["R4.closure_to_suture"]["metrics"]
    assert closure["opening_basin_id"] == closure["closing_basin_id"]
    assert closure["closing_basin_id"] == closure["suture_basin_id"]
    assert closure["suture_stage"] == "suture_relict"

    gateway = benchmarks["R4.gateway_causality"]["metrics"]
    assert gateway["gateway_count"] >= 1
    assert gateway["parent_basin_id"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "r4_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_r5_continental_margin_evolution_microbenchmarks_pass(tmp_path):
    summary = run_suite("R5", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.r5.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["no_arbitrary_area_flip"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "R5.craton_survival",
        "R5.arc_accretion",
        "R5.passive_margin_progradation",
        "R5.rifted_margin_split",
        "R5.no_arbitrary_area_flip",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    craton = benchmarks["R5.craton_survival"]["metrics"]
    assert craton["lost_margin_cells"] > 0
    assert craton["lost_craton_cells"] == 0

    arc = benchmarks["R5.arc_accretion"]["metrics"]
    assert arc["gained_cells"] > 0
    assert arc["gained_outside_accretion_zone"] == 0

    progradation = benchmarks["R5.passive_margin_progradation"]["metrics"]
    assert progradation["shape_pressure"] > 0.35
    assert progradation["gained_cells"] > 0
    assert progradation["gained_outside_progradation_zone"] == 0
    assert progradation["gained_arc_origin_cells"] == 0
    assert progradation["recorded_progradation_cells"] > 0

    rift = benchmarks["R5.rifted_margin_split"]["metrics"]
    assert rift["candidate_cells"] >= 8
    assert rift["candidate_cells_in_weak_belt"] > rift["candidate_cells_outside_weak_belt"]

    no_flip = benchmarks["R5.no_arbitrary_area_flip"]["metrics"]
    assert no_flip["unforced_gained_cells"] == 0
    assert no_flip["unforced_lost_cells"] == 0
    assert no_flip["blocked_gain_flag"] == 1.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "r5_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_r6_object_derived_terrain_microbenchmarks_pass(tmp_path):
    summary = run_suite("R6", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.r6.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["object_only_terrain_sources"]
    assert summary["acceptance"]["compiler_consistency"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "R6.passive_margin_profile",
        "R6.active_margin_profile",
        "R6.ridge_transform_bathymetry",
        "R6.orogen_width_and_decay",
        "R6.compiler_consistency",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    passive = benchmarks["R6.passive_margin_profile"]["metrics"]
    assert passive["shelf_depth_p75_m"] < 800.0
    assert passive["far_ocean_depth_p50_m"] > passive["rise_depth_p50_m"]
    assert passive["passive_object_classified"]

    active = benchmarks["R6.active_margin_profile"]["metrics"]
    assert active["trench_depth_median_m"] > 3300.0
    assert active["landward_arc_relief_p75_m"] > active["far_continent_relief_median_m"]

    ridge = benchmarks["R6.ridge_transform_bathymetry"]["metrics"]
    assert ridge["ridge_depth_median_m"] < ridge["far_abyss_depth_median_m"]
    assert ridge["transform_depth_median_m"] > ridge["ridge_depth_median_m"]

    orogen = benchmarks["R6.orogen_width_and_decay"]["metrics"]
    assert orogen["recent_orogen_uplift_median_m"] > orogen["old_orogen_uplift_median_m"]

    compiler = benchmarks["R6.compiler_consistency"]["metrics"]
    assert compiler["passed_envelope"]
    assert compiler["terrain_elevation_sign_mismatch_fraction"] < 0.02

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "r6_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p20_resolved_plate_topology_microbenchmarks_pass(tmp_path):
    summary = run_suite("P20", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p20.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["resolved_plate_topology_objects_available"]
    assert summary["acceptance"]["topology_guided_plate_split_available"]
    assert summary["acceptance"]["topology_guided_plate_capture_available"]
    assert summary["acceptance"]["continent_terrane_lifecycle_available"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P20.resolved_plate_topology_objects",
        "P20.topology_guided_plate_split",
        "P20.topology_guided_plate_capture",
        "P20.continent_terrane_lifecycle_events",
    }
    topology = benchmarks["P20.resolved_plate_topology_objects"]["metrics"]
    assert topology["object_count"] == 3
    assert topology["middle_neighbours"] == [0, 2]
    assert 0 in topology["middle_continent_ids"]
    assert 0 in topology["middle_terrane_ids"]
    assert topology["middle_trench_fraction"] > 0.0
    assert topology["middle_ridge_fraction"] > 0.0
    assert topology["middle_transform_fraction"] > 0.0
    split = benchmarks["P20.topology_guided_plate_split"]["metrics"]
    assert split["split_count"] == 1
    assert split["seed_in_weak_belt"]
    assert split["new_plate_weak_belt_fraction"] > 0.05
    assert split["motion_source"] == "p20_topology_split_inherited"
    assert split["deterministic_plate_equal"]
    capture = benchmarks["P20.topology_guided_plate_capture"]["metrics"]
    assert capture["oceanic_capture_count"] == 1
    assert capture["captured_from"] == 1
    assert capture["captured_to"] == 0
    assert capture["capture_basis"] == "p20_topology_shared_edge_capture"
    assert capture["shared_edge_fraction"] >= 0.95
    assert capture["continental_cargo_capture_count"] == 1
    assert capture["continental_cargo_policy"] == "capture_plate_label_preserve_crust_cargo"
    assert capture["continental_source_fraction"] >= 0.95
    assert capture["continental_cargo_crust_preserved"]
    lineage = benchmarks["P20.continent_terrane_lifecycle_events"]["metrics"]
    assert set(lineage["object_kinds"]) >= {
        "continent_split",
        "continent_merge",
        "terrane_capture",
        "microcontinent_plate_capture",
    }
    assert lineage["split_parent_ids"] == [0]
    assert set(lineage["split_child_ids"]) == {0, 2}
    assert set(lineage["merge_parent_ids"]) == {5, 6}
    assert lineage["merge_child_ids"] == [5]
    assert lineage["terrane_capture_parent"] == 5
    assert lineage["microcontinent_capture_policy"] == "capture_plate_label_preserve_crust_cargo"

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p20_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p21_supercontinent_breakup_microbenchmarks_pass(tmp_path):
    summary = run_suite("P21", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p21.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["breakup_seaway_objects_available"]
    assert summary["acceptance"]["object_driven_terrain_seaway_available"]
    assert summary["acceptance"]["continental_sediment_depocenter_contrast"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P21.breakup_seaway_objects",
        "P21.object_driven_terrain_seaway",
        "P21.continental_sediment_depocenters",
    }
    objects = benchmarks["P21.breakup_seaway_objects"]["metrics"]
    assert objects["object_count"] >= 1
    assert objects["mean_rift_potential"] > 0.45
    assert objects["weak_fraction"] > 0.35
    assert objects["core_fraction"] < 0.25
    assert objects["path_hits_weak_belt"] > 0.45

    terrain = benchmarks["P21.object_driven_terrain_seaway"]["metrics"]
    assert terrain["before_largest_land_fraction"] > 0.95
    assert terrain["after_largest_land_fraction"] < 0.72
    assert terrain["after_component_count"] >= 2
    assert terrain["core_preserved"]
    assert terrain["recorded_openings"] >= 1.0

    sediment = benchmarks["P21.continental_sediment_depocenters"]["metrics"]
    assert sediment["basin_mean_sediment_m"] > sediment["platform_mean_sediment_m"] + 900.0
    assert sediment["basin_mean_sediment_m"] > sediment["highland_mean_sediment_m"] + 1100.0
    assert sediment["craton_mean_sediment_m"] < 1300.0
    assert sediment["recorded_contrast_m"] > 650.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p21_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p22_medium_supercontinent_breakup_microbenchmark_pass(tmp_path):
    summary = run_suite("P22", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p22.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["medium_supercontinent_partition_available"]
    assert not summary["acceptance"]["medium_gap_reproduced"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P22.medium_supercontinent_breakup_partition"}
    metrics = benchmarks["P22.medium_supercontinent_breakup_partition"]["metrics"]
    assert metrics["breakup_object_count"] >= 1
    assert metrics["first_object_topology_score"] > 0.55
    assert metrics["first_object_split_balance"] > 0.45
    assert metrics["recorded_openings"] >= 1.0
    assert metrics["after_largest_land_fraction"] < 0.68
    assert metrics["after_largest_land_fraction"] < metrics["before_largest_land_fraction"]
    assert metrics["split_balance"] > 0.22
    assert metrics["interior_rift_ocean_fraction_after"] > 0.18

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p22_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p23_medium_rifted_component_breakup_microbenchmark_pass(tmp_path):
    summary = run_suite("P23", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p23.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["medium_rifted_component_breakup_objects"]
    assert summary["acceptance"]["divergent_boundary_rift_lifecycle"]
    assert summary["acceptance"]["medium_divergent_component_breakup_objects"]
    assert summary["acceptance"]["multi_corridor_component_breakup_objects"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P23.medium_rifted_component_breakup_objects",
        "P23.divergent_boundary_rift_lifecycle",
        "P23.medium_divergent_component_breakup_objects",
        "P23.multi_corridor_component_breakup_objects",
    }
    metrics = benchmarks["P23.medium_rifted_component_breakup_objects"]["metrics"]
    assert metrics["medium_component_count"] == 3
    assert min(metrics["component_area_fractions"]) > 0.080
    assert max(metrics["component_area_fractions"]) < 0.115
    assert metrics["breakup_object_count"] >= 3
    assert metrics["strong_breakup_object_count"] >= 3
    assert metrics["covered_parent_continent_ids"] == [0, 1, 2]
    assert metrics["mean_strong_topology_score"] > 0.75
    assert metrics["mean_strong_split_balance"] > 0.65

    lifecycle = benchmarks["P23.divergent_boundary_rift_lifecycle"]["metrics"]
    assert lifecycle["rift_system_count"] == 1
    assert lifecycle["spreading_center_count"] == 0
    assert lifecycle["rift_stage"] == "continental_rift"
    assert lifecycle["rift_phase_code"] == 1.0
    assert lifecycle["basin_parent_rift_count"] == 1

    divergent = benchmarks["P23.medium_divergent_component_breakup_objects"]["metrics"]
    assert divergent["breakup_object_count"] >= 3
    assert divergent["strong_breakup_object_count"] >= 3
    assert divergent["covered_parent_continent_ids"] == [0, 1, 2]
    assert divergent["mean_strong_boundary_fraction"] > 0.65
    assert divergent["mean_strong_topology_score"] > 0.75

    multi = benchmarks["P23.multi_corridor_component_breakup_objects"]["metrics"]
    assert multi["candidate_count"] >= 6
    assert multi["accepted_object_count"] >= 3
    assert multi["breakup_object_count"] >= 3
    assert sum(1 for hits in multi["corridor_hits"] if hits > 0) >= 3
    assert multi["best_topology_score"] > 0.75

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p23_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p24_object_to_terrain_seaway_microbenchmark_pass(tmp_path):
    summary = run_suite("P24", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p24.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["submerged_breakup_axis_still_propagates"]
    assert summary["acceptance"]["opened_breakup_seaway_survives_ocean_floor_regionalization"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P24.submerged_breakup_axis_still_propagates",
        "P24.opened_breakup_seaway_survives_ocean_floor_regionalization",
    }
    metrics = benchmarks["P24.submerged_breakup_axis_still_propagates"]["metrics"]
    assert metrics["object_land_fraction_before"] == 0.0
    assert metrics["source_reuse_count"] >= 1.0
    assert metrics["recorded_openings"] >= 1.0
    assert metrics["attempt_count"] >= 1
    assert metrics["applied_attempt_count"] >= 1
    assert metrics["best_attempt_reduction"] > 0.25
    assert metrics["best_attempt_candidate_count"] >= 1
    assert metrics["best_attempt_viable_candidate_count"] >= 1
    assert metrics["best_attempt_reused_process_source"]
    assert metrics["before_largest_land_fraction"] > 0.95
    assert metrics["after_largest_land_fraction"] < 0.72
    assert metrics["split_balance"] > 0.35
    assert metrics["propagated_axis_ocean_fraction"] > 0.82

    regionalized = benchmarks[
        "P24.opened_breakup_seaway_survives_ocean_floor_regionalization"
    ]["metrics"]
    assert regionalized["before_largest_land_fraction"] < 0.58
    assert regionalized["after_largest_land_fraction"] < 0.58
    assert regionalized["protected_ocean_fraction"] == 1.0
    assert regionalized["protected_depth_min_m"] >= 250.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p24_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p25_earth_reference_highland_calibration_microbenchmarks_pass(tmp_path):
    summary = run_suite("P25", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p25.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["broad_mature_highland_relaxation"]
    assert summary["acceptance"]["highland_province_requires_active_or_thick_core"]
    assert summary["acceptance"]["modern_coastline_payback_preserves_breakup_seaways"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P25.broad_mature_highland_relaxation",
        "P25.highland_province_requires_active_or_thick_core",
        "P25.modern_coastline_payback_preserves_breakup_seaways",
    }
    relax = benchmarks["P25.broad_mature_highland_relaxation"]["metrics"]
    assert relax["after_broad_mature_mean_m"] < relax["before_broad_mature_mean_m"] - 1000.0
    assert relax["after_broad_mature_high_fraction_gt2500"] < 0.08
    assert relax["after_active_orogen_p50_m"] > 3600.0
    assert relax["recorded_broad_interior_relaxed_fraction"] > 0.03

    province = benchmarks[
        "P25.highland_province_requires_active_or_thick_core"
    ]["metrics"]
    assert province["mature_broad_highland_province_share"] < 0.20
    assert province["active_orogen_highland_province_share"] > 0.60
    assert province["craton_shield_province_share"] > 0.70
    assert province["global_highland_province_share"] < 0.35

    coastline = benchmarks[
        "P25.modern_coastline_payback_preserves_breakup_seaways"
    ]["metrics"]
    assert coastline["drowned_area_fraction"] > 0.004
    assert coastline["payback_candidate_area_fraction"] > 0.006
    assert coastline["payback_area_fraction"] > 0.004
    assert coastline["after_land_fraction"] >= coastline["before_land_fraction"] - 0.009
    assert coastline["protected_seaway_land_fraction"] == 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p25_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()


def test_p26_recent_rework_footprint_microbenchmarks_pass(tmp_path):
    summary = run_suite("P26", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p26.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["diagnostic_only"]
    assert summary["acceptance"]["recent_rework_footprint_flags_overbroad_swath"]
    assert summary["acceptance"]["recent_rework_footprint_accepts_localized_belt"]
    assert summary["acceptance"]["deforming_network_state_separates_core_and_shoulder"]
    assert summary["acceptance"]["deforming_network_state_localizes_rift_axes"]
    assert summary["acceptance"][
        "terrain_response_uses_deforming_network_not_broad_recent_swath"
    ]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P26.recent_rework_footprint_flags_overbroad_swath",
        "P26.recent_rework_footprint_accepts_localized_belt",
        "P26.deforming_network_state_separates_core_and_shoulder",
        "P26.deforming_network_state_localizes_rift_axes",
        "P26.terrain_response_uses_deforming_network_not_broad_recent_swath",
    }
    broad = benchmarks["P26.recent_rework_footprint_flags_overbroad_swath"]["metrics"]
    assert broad["overbroad_recent_rework"]
    assert broad["active_rework_outside_corridor_fraction_of_active"] > 0.45
    assert broad["active_to_corridor_area_ratio"] > 1.45

    localized = benchmarks["P26.recent_rework_footprint_accepts_localized_belt"]["metrics"]
    assert not localized["overbroad_recent_rework"]
    assert localized["active_rework_outside_corridor_fraction_of_active"] < 0.20
    assert localized["active_rework_inside_corridor_fraction_of_active"] > 0.90
    assert localized["active_to_corridor_area_ratio"] < 1.15
    deforming = benchmarks[
        "P26.deforming_network_state_separates_core_and_shoulder"
    ]["metrics"]
    assert deforming["core_cells"] > 0
    assert deforming["shoulder_cells"] > 0
    assert deforming["origin_unchanged"]
    assert {"collision_core", "collision_shoulder"} <= set(deforming["object_kinds"])
    rift = benchmarks["P26.deforming_network_state_localizes_rift_axes"]["metrics"]
    assert rift["rift_cells"] > 0
    assert rift["noncontact_rift_cells"] == 0
    terrain = benchmarks[
        "P26.terrain_response_uses_deforming_network_not_broad_recent_swath"
    ]["metrics"]
    assert terrain["core_active_highland_detail_fraction"] > 0.50
    assert terrain["nondeforming_recent_active_highland_detail_fraction"] < 0.35
    assert terrain["core_active_highland_province_fraction"] > 0.50
    assert terrain["nondeforming_recent_active_highland_province_fraction"] < 0.35
    assert terrain["core_relief_m"] > terrain["shoulder_relief_m"] + 45.0
    assert terrain["shoulder_relief_m"] > 20.0
    assert terrain["nondeforming_recent_relief_m"] < 1.0
    assert terrain["land_mask_unchanged"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p26_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p27_orogenic_foreland_response_microbenchmark_pass(tmp_path):
    summary = run_suite("P27", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p27.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["orogenic_foreland_response_state_available"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P27.orogenic_foreland_response_state"}
    terrain = benchmarks["P27.orogenic_foreland_response_state"]["metrics"]
    assert terrain["orogen_load_mean"] > terrain["foreland_load_mean"]
    assert terrain["foreland_load_mean"] > terrain["distal_load_mean"]
    assert terrain["foreland_accommodation_mean"] > terrain["orogen_accommodation_mean"]
    assert terrain["foreland_accommodation_mean"] > terrain["distal_accommodation_mean"]
    assert terrain["foreland_accommodation_mean"] > terrain["craton_accommodation_mean"]
    assert terrain["foreland_basin_detail_fraction"] >= 0.35
    assert terrain["foreland_basin_objects"] >= 1
    assert terrain["foreland_mean_sediment_m"] > terrain["distal_platform_mean_sediment_m"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p27_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p28_sediment_coupling_budget_gates_pass(tmp_path):
    summary = run_suite("P28", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p28.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["sediment_coupling_budget_gate_available"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P28.foreland_sediment_budget_gate",
        "P28.passive_margin_sediment_budget_gate",
    }
    foreland = benchmarks["P28.foreland_sediment_budget_gate"]["metrics"]
    assert foreland["foreland_delta_mean_m"] > 250.0
    assert foreland["orogen_delta_mean_m"] < 0.0
    assert foreland["volume_change_fraction"] < 0.002
    assert foreland["projected_land_change_fraction"] < 0.01

    passive = benchmarks["P28.passive_margin_sediment_budget_gate"]["metrics"]
    assert passive["margin_delta_mean_m"] > 200.0
    assert passive["shelf_delta_mean_m"] > 150.0
    assert passive["interior_delta_mean_m"] < 0.0
    assert passive["volume_change_fraction"] < 0.002
    assert passive["projected_land_change_fraction"] < 0.012

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p28_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p29_inland_geomorphology_complexity_gates_pass(tmp_path):
    summary = run_suite("P29", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p29.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["inland_mechanism_diversity_fixture"]
    assert summary["acceptance"]["inland_relief_response_preserves_land_mask"]
    assert summary["acceptance"]["flat_interior_detector_available"]
    assert summary["acceptance"]["unparented_speckle_detector_available"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P29.inland_mechanism_diversity_fixture",
        "P29.inland_relief_response_preserves_land_mask",
        "P29.flat_interior_detector",
        "P29.unparented_speckle_detector",
    }
    diversity = benchmarks["P29.inland_mechanism_diversity_fixture"]["metrics"]
    assert diversity["inland_relief_p95_p05_m"] >= 2200.0
    assert diversity["inland_relief_p90_p10_m"] >= 750.0
    assert diversity["inland_detail_diversity"] >= 5
    assert diversity["object_kind_count"] >= 6
    assert diversity["shield_objects"] >= 1
    assert diversity["interior_basin_objects"] >= 1
    assert diversity["old_subdued_orogen_objects"] >= 1
    assert diversity["rift_basin_objects"] >= 1
    assert diversity["plateau_objects"] >= 1
    assert diversity["parented_shield_objects"] >= 1
    assert diversity["parented_old_subdued_orogen_objects"] >= 1
    assert diversity["parented_rift_basin_objects"] >= 1
    assert diversity["parented_plateau_objects"] >= 1
    assert 0.45 <= diversity["inland_flat_lowland_fraction"] <= 0.90

    relief = benchmarks["P29.inland_relief_response_preserves_land_mask"]["metrics"]
    assert relief["land_mask_unchanged"]
    assert relief["old_orogen_delta_mean_m"] > 65.0
    assert relief["plateau_delta_mean_m"] > 75.0
    assert relief["interior_basin_delta_mean_m"] < -45.0
    assert relief["rift_basin_delta_mean_m"] < -75.0
    assert relief["shield_abs_delta_mean_m"] < 8.0
    assert relief["platform_abs_delta_mean_m"] < 8.0

    flat = benchmarks["P29.flat_interior_detector"]["metrics"]
    assert flat["flatness_detected"]
    assert flat["inland_relief_p90_p10_m"] < 220.0
    assert flat["object_kind_count"] <= 2

    speckle = benchmarks["P29.unparented_speckle_detector"]["metrics"]
    assert speckle["speckle_detected"]
    assert speckle["parented_object_kind_count"] == 0
    assert speckle["speckle_high_component_count"] >= 80

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p29_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p31_planform_reference_microbenchmarks_pass(tmp_path):
    summary = run_suite("P31", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p31.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["planform_reference_gate_available"]
    assert summary["acceptance"]["ribbon_land_detector_available"]
    assert summary["acceptance"]["overcomplex_coastline_detector_available"]
    assert summary["acceptance"]["diagnostic_only"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P31.broad_multicontinent_planform_reference",
        "P31.ribbon_land_planform_detector",
        "P31.overcomplex_coastline_planform_detector",
    }
    broad = benchmarks["P31.broad_multicontinent_planform_reference"]["metrics"]
    assert broad["within_initial_planform_envelope"]
    assert broad["land_fraction_in_envelope"]
    assert broad["land_component_count_in_envelope"]
    assert broad["land_ribbon_fraction_gt_0_5"] < 0.35
    assert broad["land_coastline_complexity_largest"] < 8.0
    assert broad["resolution_ladder_declared"]

    ribbon = benchmarks["P31.ribbon_land_planform_detector"]["metrics"]
    assert ribbon["ribbon_outside_initial_reference"]
    assert ribbon["land_ribbon_fraction_gt_0_5"] > 0.85
    assert ribbon["land_width_p50_steps"] <= 1.0
    assert "land_ribbon_fraction_gt_0_5" in ribbon["out_of_envelope"]

    coastline = benchmarks["P31.overcomplex_coastline_planform_detector"]["metrics"]
    assert coastline["coastline_outside_initial_reference"]
    assert not coastline["ribbon_outside_initial_reference"]
    assert coastline["land_coastline_complexity_largest"] > 11.0
    assert coastline["land_ribbon_fraction_gt_0_5"] < 0.35

    release_entry = {
        "land_fraction": 0.252,
        "morphology": {
            "largest_land_component_fraction": 0.48,
            "land_component_count": 4,
            "land_ribbon_fraction_gt_0_5": 0.424,
            "continental_ribbon_fraction_gt_0_5": 0.397,
            "land_coastline_complexity_largest": 23.27,
        },
    }
    release_summary = planform_reference_summary(morphology_metrics=release_entry)
    assert release_summary["status"] == "needs_planform_calibration"
    assert release_summary["diagnostic_flags"]["ribbon_outside_initial_reference"]
    assert release_summary["diagnostic_flags"]["coastline_outside_initial_reference"]
    assert release_summary["metrics"]["land_fraction"]["in_envelope"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p31_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p32_coastline_simplification_microbenchmarks_pass(tmp_path):
    summary = run_suite("P32", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p32.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["coastline_smoothing_reduces_complexity"]
    assert summary["acceptance"]["object_constrained_coastline_simplification"]
    assert summary["acceptance"]["protected_breakup_seaway_preserved"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P32.coastline_smoothing_reduces_complexity",
        "P32.object_constrained_coastline_simplification",
        "P32.protected_breakup_seaway_not_refilled",
    }
    production = benchmarks["P32.coastline_smoothing_reduces_complexity"]["metrics"]
    assert production["before_coastline_complexity"] > 12.0
    assert production["coastline_complexity_drop"] >= 0.75
    assert (
        production["after_coastline_complexity"]
        <= production["before_coastline_complexity"] * 0.94
    )
    assert production["ribbon_fraction_drop"] >= 0.018
    assert production["after_ribbon_fraction"] <= 0.35
    assert 0.25 <= production["after_land_fraction"] <= 0.31
    assert production["after_component_count"] >= 2
    assert production["protected_seaway_land_fraction"] == 0.0
    assert production["fill_area_fraction"] > 0.008
    assert (
        production["payback_area_fraction"] > 0.004
        or production["after_land_fraction"] >= production["before_land_fraction"] - 0.003
    )

    simplified = benchmarks["P32.object_constrained_coastline_simplification"]["metrics"]
    assert simplified["coastline_complexity_delta"] >= 0.35
    assert simplified["land_fraction_delta"] <= 0.003
    assert simplified["p32_swap_area_fraction"] > 0.00025

    protected = benchmarks["P32.protected_breakup_seaway_not_refilled"]["metrics"]
    assert protected["protected_still_ocean"]
    assert protected["coastline_complexity_delta"] >= 0.35
    assert protected["p32_swap_area_fraction"] > 0.00025

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p32_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p33_parented_margin_sliver_microbenchmarks_pass(tmp_path):
    summary = run_suite("P33", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p33.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["parented_margin_sliver_consolidation"]
    assert summary["acceptance"]["state_only_no_land_mask_mutation"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P33.parented_margin_sliver_consolidation"}
    metrics = benchmarks["P33.parented_margin_sliver_consolidation"]["metrics"]
    assert metrics["consolidated_cells"] > 0
    assert metrics["after_mature_accretionary_cells"] < metrics["before_mature_accretionary_cells"]
    assert metrics["active_changed_cells"] == 0
    assert metrics["detached_changed_cells"] == 0
    assert metrics["land_mask_unchanged"]
    assert metrics["after_shape_pressure"] < metrics["before_shape_pressure"]
    assert metrics["recorded_consolidated_cells"] > 0.0
    assert metrics["recorded_active_preserved_cells"] > 0.0
    assert metrics["recorded_detached_preserved_cells"] > 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p33_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p34_parented_platform_progradation_microbenchmarks_pass(tmp_path):
    summary = run_suite("P34", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p34.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["parented_platform_progradation"]
    assert summary["acceptance"]["no_active_or_detached_or_protected_growth"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P34.parented_platform_progradation"}
    metrics = benchmarks["P34.parented_platform_progradation"]["metrics"]
    assert metrics["after_continental_fraction"] > metrics["before_continental_fraction"]
    assert metrics["gained_cells"] > 0
    assert metrics["parented_gain_cells"] == metrics["gained_cells"]
    assert metrics["arc_gain_cells"] == 0
    assert metrics["old_arc_restored_cells"] > 0
    assert metrics["old_shell_restored_cells"] > 0
    assert metrics["protected_gain_cells"] == 0
    assert metrics["detached_gain_cells"] == 0
    assert metrics["active_gain_cells"] == 0
    assert metrics["recorded_candidate_cells"] > 0.0
    assert metrics["recorded_restored_cells"] > 0.0
    assert metrics["recorded_total_progradation_cells"] > 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p34_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p35_continental_debt_recovery_microbenchmarks_pass(tmp_path):
    summary = run_suite("P35", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p35.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["continental_debt_recovery"]
    assert summary["acceptance"]["no_active_or_detached_or_protected_growth"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P35.continental_debt_recovery"}
    metrics = benchmarks["P35.continental_debt_recovery"]["metrics"]
    assert metrics["after_continental_fraction"] > metrics["before_continental_fraction"]
    assert metrics["gained_cells"] > 0
    assert metrics["pocket_gain_cells"] > 0
    assert metrics["parented_gain_cells"] == metrics["gained_cells"]
    assert metrics["arc_gain_cells"] == 0
    assert metrics["protected_gain_cells"] == 0
    assert metrics["active_gain_cells"] == 0
    assert metrics["detached_gain_cells"] == 0
    assert metrics["recorded_candidate_cells"] > 0.0
    assert metrics["recorded_restored_cells"] > 0.0
    assert metrics["recorded_stable_parent_cells"] > 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p35_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p36_mature_platform_anchor_recovery_microbenchmarks_pass(tmp_path):
    summary = run_suite("P36", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p36.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["mature_platform_anchor_recovery"]
    assert summary["acceptance"]["preconserve_craton_readiness"]
    assert summary["acceptance"]["legacy_anchor_absent_but_process_anchor_present"]
    assert summary["acceptance"]["no_active_or_protected_growth"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P36.mature_platform_anchor_recovery",
        "P36.preconserve_craton_readiness",
    }
    metrics = benchmarks["P36.mature_platform_anchor_recovery"]["metrics"]
    assert metrics["legacy_anchor_cells"] == 0
    assert metrics["process_anchor_cells"] > 0
    assert metrics["after_continental_fraction"] > metrics["before_continental_fraction"]
    assert metrics["gained_cells"] > 0
    assert metrics["pocket_gain_cells"] > 0
    assert metrics["parented_gain_cells"] == metrics["gained_cells"]
    assert metrics["arc_gain_cells"] == 0
    assert metrics["protected_gain_cells"] == 0
    assert metrics["active_gain_cells"] == 0
    assert metrics["recorded_candidate_cells"] > 0.0
    assert metrics["recorded_restored_cells"] > 0.0
    assert metrics["recorded_old_platform_anchor_cells"] > 0.0
    readiness = benchmarks["P36.preconserve_craton_readiness"]["metrics"]
    assert readiness["before_craton_cells"] == 0
    assert readiness["before_stable_cells"] == 0
    assert readiness["promoted_craton_cells"] > 0
    assert readiness["after_stable_cells"] > readiness["before_stable_cells"]
    assert readiness["anchor_after_cells"] >= readiness["promoted_craton_cells"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p36_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p37_shape_guard_restoration_microbenchmarks_pass(tmp_path):
    summary = run_suite("P37", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p37.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["shape_guard_restoration"]
    assert summary["acceptance"]["rejects_ribbon_fragment_restoration"]
    assert summary["acceptance"]["accepts_broad_parent_restoration"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P37.shape_guard_restoration"}
    metrics = benchmarks["P37.shape_guard_restoration"]["metrics"]
    assert metrics["broad_candidate_cells"] > 0
    assert metrics["broad_accepted_cells"] > 0
    assert metrics["ribbon_candidate_cells"] > 0
    assert metrics["ribbon_accepted_cells"] == 0
    assert metrics["accepted_cells"] < metrics["candidate_cells"]
    assert metrics["broad_parent_width_p50"] >= 2.75
    assert metrics["ribbon_parent_width_p50"] < 2.75

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p37_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p38_lineage_guard_restoration_microbenchmarks_pass(tmp_path):
    summary = run_suite("P38", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p38.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["lineage_guard_restoration"]
    assert summary["acceptance"]["rejects_orphan_parent_restoration"]
    assert summary["acceptance"]["accepts_valid_broad_lineage_restoration"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P38.lineage_guard_restoration"}
    metrics = benchmarks["P38.lineage_guard_restoration"]["metrics"]
    assert metrics["shape_broad_accepted_cells"] > 0
    assert metrics["shape_orphan_accepted_cells"] > 0
    assert metrics["lineage_broad_accepted_cells"] > 0
    assert metrics["lineage_orphan_accepted_cells"] == 0
    assert metrics["lineage_ribbon_accepted_cells"] == 0
    assert metrics["lineage_guarded_accepted_cells"] < metrics["shape_only_accepted_cells"]
    assert metrics["recorded_lineage_input_cells"] > 0.0
    assert metrics["recorded_lineage_accepted_cells"] > 0.0
    assert metrics["recorded_lineage_rejected_cells"] > 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p38_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p39_old_platform_anchor_funnel_microbenchmarks_pass(tmp_path):
    summary = run_suite("P39", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p39.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["old_platform_anchor_funnel"]
    assert summary["acceptance"]["records_all_blocker_classes"]
    assert summary["acceptance"]["distinguishes_lineaged_from_orphan_anchor"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P39.old_platform_anchor_funnel"}
    metrics = benchmarks["P39.old_platform_anchor_funnel"]["metrics"]
    assert metrics["legacy_anchor_cells"] == 0
    assert metrics["process_anchor_cells"] > 0
    assert metrics["valid_lineaged_anchor_cells"] > 0
    assert metrics["orphan_process_anchor_cells"] > 0
    assert metrics["lineaged_anchor_cells"] == metrics["valid_lineaged_anchor_cells"]
    assert metrics["process_anchor_lineaged_cells"] > 0
    assert metrics["process_anchor_unlineaged_cells"] > 0
    assert metrics["source_blocked_cells"] > 0
    assert metrics["thick_blocked_cells"] > 0
    assert metrics["recently_reworked_cells"] > 0
    assert metrics["age_old_blocked_cells"] > 0
    assert metrics["stability_ge030_blocked_cells"] > 0
    assert metrics["width_ge3_blocked_cells"] > 0
    assert metrics["recorded_process_anchor_cells"] == metrics["process_anchor_cells"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p39_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p40_rework_semantics_microbenchmarks_pass(tmp_path):
    summary = run_suite("P40", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p40.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["parented_restoration_inherits_rework"]
    assert summary["acceptance"]["shape_maintenance_rework_semantics"]
    assert summary["acceptance"]["active_accretion_still_reworks"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P40.parented_restoration_inherits_rework",
        "P40.shape_maintenance_rework_semantics",
    }
    restoration = benchmarks["P40.parented_restoration_inherits_rework"]["metrics"]
    assert restoration["gained_cells"] > 0
    assert restoration["quiet_inherited_gained_cells"] == restoration["gained_cells"]
    assert restoration["current_reworked_gained_cells"] == 0
    assert (
        restoration["recorded_parented_inherited_rework_cells"]
        == restoration["gained_cells"]
    )
    assert restoration["recorded_parented_reactivated_cells"] == 0

    shape = benchmarks["P40.shape_maintenance_rework_semantics"]["metrics"]
    assert shape["quiet_hole_filled"]
    assert shape["quiet_hole_reworked_age_myr"] < 0.0
    assert shape["quiet_hole_origin"] == 1.0
    assert shape["quiet_recorded_inherited_cells"] == 1.0
    assert shape["quiet_recorded_reactivated_cells"] == 0.0
    assert shape["active_hole_filled"]
    assert shape["active_hole_reworked_age_myr"] == 4500.0
    assert shape["active_hole_origin"] == 2.0
    assert shape["active_recorded_inherited_cells"] == 0.0
    assert shape["active_recorded_reactivated_cells"] == 1.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p40_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p41_mature_debt_priority_microbenchmarks_pass(tmp_path):
    summary = run_suite("P41", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p41.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["mature_debt_prefers_parented_restoration"]
    assert summary["acceptance"]["defers_active_only_accretion"]
    assert summary["acceptance"]["preserves_quiet_parented_rework"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P41.mature_debt_prefers_parented_restoration"}
    metrics = benchmarks["P41.mature_debt_prefers_parented_restoration"]["metrics"]
    assert metrics["gained_cells"] > 0
    assert metrics["pocket_gain_cells"] > 0
    assert metrics["active_gain_cells"] == 0
    assert metrics["arc_gain_cells"] == 0
    assert metrics["quiet_inherited_gain_cells"] == metrics["gained_cells"]
    assert metrics["recorded_parented_priority_cells"] == metrics["gained_cells"]
    assert metrics["recorded_active_only_deferred_cells"] > 0
    assert metrics["recorded_p40_reactivated_cells"] == 0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p41_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p42_shape_supported_craton_microbenchmarks_pass(tmp_path):
    summary = run_suite("P42", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p42.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["shape_supported_craton_domain"]
    assert summary["acceptance"]["kernel_maintenance_rejects_all_ribbon_continent"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P42.shape_supported_craton_domain",
        "P42.kernel_maintenance_rejects_all_ribbon_continent",
    }
    domain = benchmarks["P42.shape_supported_craton_domain"]["metrics"]
    assert domain["broad_craton_fraction"] > 0.60
    assert domain["ribbon_craton_fraction"] == 0.0
    assert domain["raw_ribbon_stability_max"] > 0.78
    assert domain["guarded_ribbon_stability_max"] <= 0.74

    kernel = benchmarks[
        "P42.kernel_maintenance_rejects_all_ribbon_continent"]["metrics"]
    assert kernel["ribbon_cells"] > 0
    assert kernel["promoted_craton_cells"] == 0
    assert kernel["after_ribbon_stability_max"] <= 0.75

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p42_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p43_mature_active_accretion_shape_guard_microbenchmarks_pass(tmp_path):
    summary = run_suite("P43", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p43.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["mature_active_accretion_shape_guard"]
    assert summary["acceptance"]["rejects_unsupported_active_only_cells"]
    assert summary["acceptance"]["caps_mature_active_only_area"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P43.mature_active_accretion_shape_guard",
        "P43.payback_fills_broad_embayment_not_terrane_ribbon",
        "P43.payback_recovers_supported_accreted_collage",
    }
    metrics = benchmarks["P43.mature_active_accretion_shape_guard"]["metrics"]
    assert metrics["gained_cells"] > 0
    assert metrics["unsupported_gain_cells"] == 0
    assert metrics["recorded_candidate_cells"] > 0
    assert metrics["recorded_accepted_cells"] == metrics["gained_cells"]
    assert metrics["recorded_rejected_cells"] > 0
    assert metrics["accepted_area_fraction"] <= metrics["cap_area_fraction"] + (1.0 / 3600.0)

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p43_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p43_coastline_payback_microbenchmarks_pass(tmp_path):
    summary = run_suite("P43", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p43.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["payback_fills_broad_embayment_not_terrane_ribbon"]
    assert summary["acceptance"]["does_not_fill_unsupported_terrane_or_suture"]
    assert summary["acceptance"]["payback_recovers_supported_accreted_collage"]
    assert summary["acceptance"]["collage_recovery_has_support"]
    assert summary["acceptance"]["does_not_fragment_or_worsen_ribbon"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P43.mature_active_accretion_shape_guard",
        "P43.payback_fills_broad_embayment_not_terrane_ribbon",
        "P43.payback_recovers_supported_accreted_collage",
    }
    metrics = benchmarks[
        "P43.payback_fills_broad_embayment_not_terrane_ribbon"]["metrics"]
    assert metrics["filled_cells"] > 0
    assert metrics["filled_supported_cells"] == metrics["filled_cells"]
    assert metrics["filled_unsupported_terrane_or_suture_cells"] == 0
    assert metrics["after_component_count"] <= metrics["before_component_count"]
    assert metrics["after_ribbon_fraction"] <= metrics["before_ribbon_fraction"] + 0.012
    assert metrics["protected_seaway_land_fraction"] == 0.0
    assert metrics["payback_area_fraction"] > 0.0
    collage_metrics = benchmarks[
        "P43.payback_recovers_supported_accreted_collage"]["metrics"]
    assert collage_metrics["filled_collage_cells"] > 0
    assert collage_metrics["filled_target_collage_cells"] > 0
    assert collage_metrics["unsupported_collage_filled_cells"] == 0
    assert collage_metrics["protected_filled_cells"] == 0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p43_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p44_parented_recovery_precedes_active_only_microbenchmarks_pass(tmp_path):
    summary = run_suite("P44", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p44.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["parented_recovery_precedes_active_only"]
    assert summary["acceptance"]["parented_recovery_dominates_active_only"]
    assert summary["acceptance"]["p44_changes_budget_order"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P44.parented_recovery_precedes_active_only"}
    metrics = benchmarks["P44.parented_recovery_precedes_active_only"]["metrics"]
    assert metrics["gained_cells"] > 0
    assert metrics["parented_gain_cells"] > metrics["active_gain_cells"]
    assert metrics["pocket_gain_cells"] > 0
    assert metrics["quiet_inherited_gain_cells"] == metrics["parented_gain_cells"]
    assert metrics["protected_parented_gain_cells"] == 0
    assert metrics["p44_disabled_active_gain_share"] > metrics["active_gain_share"]
    assert metrics["p44_disabled_parented_gain_cells"] < metrics["parented_gain_cells"]
    assert metrics["recorded_candidate_cells"] > 0
    assert metrics["recorded_restored_cells"] == metrics["parented_gain_cells"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p44_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p45_process_provenance_localization_microbenchmarks_pass(tmp_path):
    summary = run_suite("P45", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p45.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["process_provenance_localized"]
    assert summary["acceptance"]["provenance_is_subset_of_process"]
    assert summary["acceptance"]["axis_preserved"]
    assert summary["acceptance"]["broad_swath_not_fully_rewritten"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P45.process_provenance_localizes_broad_deformation"}
    metrics = benchmarks[
        "P45.process_provenance_localizes_broad_deformation"]["metrics"]
    assert metrics["process_cells"] > metrics["provenance_cells"] > metrics["axis_cells"]
    assert metrics["provenance_share_of_process"] <= 0.45
    assert metrics["outside_process_cells"] == 0
    assert metrics["missing_axis_cells"] == 0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p45_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p46_mature_parent_supported_collage_microbenchmarks_pass(tmp_path):
    summary = run_suite("P46", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p46.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["mature_parent_supported_collage_belts"]
    assert summary["acceptance"]["state_only_no_land_mask_mutation"]
    assert summary["acceptance"]["active_and_detached_preserved"]
    assert summary["acceptance"]["accretionary_domain_reduced"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P46.mature_parent_supported_collage_belts"}
    metrics = benchmarks["P46.mature_parent_supported_collage_belts"]["metrics"]
    assert metrics["matured_cells"] > 0
    assert metrics["after_mature_accretionary_cells"] < metrics["before_mature_accretionary_cells"]
    assert metrics["after_mature_terrane_mask_cells"] < metrics["before_mature_terrane_mask_cells"]
    assert metrics["active_changed_cells"] == 0
    assert metrics["detached_changed_cells"] == 0
    assert metrics["land_mask_unchanged"]
    assert metrics["recorded_candidate_cells"] > 0.0
    assert metrics["recorded_consolidated_cells"] > 0.0
    assert metrics["recorded_active_preserved_cells"] > 0.0
    assert metrics["recorded_detached_preserved_cells"] > 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p46_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p47_mature_cratonic_kernel_target_microbenchmarks_pass(tmp_path):
    summary = run_suite("P47", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p47.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["mature_cratonic_kernel_target"]
    assert summary["acceptance"]["target_is_time_dependent"]
    assert summary["acceptance"]["broad_interiors_promoted"]
    assert summary["acceptance"]["ribbon_not_promoted"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P47.mature_cratonic_kernel_target"}
    metrics = benchmarks["P47.mature_cratonic_kernel_target"]["metrics"]
    assert metrics["target_fraction"] >= 0.075
    assert metrics["stable_share_of_continent"] >= 0.06
    assert metrics["ancient_share_of_continent"] >= 0.06
    assert metrics["broad_promoted_cells"] > 0
    assert metrics["ribbon_promoted_cells"] == 0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p47_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p48_cratonic_province_object_layer_microbenchmarks_pass(tmp_path):
    summary = run_suite("P48", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p48.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["cratonic_province_object_layer"]
    assert summary["acceptance"]["object_generation_does_not_mutate_land_mask"]
    assert summary["acceptance"]["active_and_accreted_regions_excluded"]
    assert summary["acceptance"]["terrain_landforms_inherit_province_parents"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P48.cratonic_province_object_layer"}
    metrics = benchmarks["P48.cratonic_province_object_layer"]["metrics"]
    assert metrics["shield_objects"] >= 1
    assert metrics["platform_objects"] >= 1
    assert metrics["interior_basin_objects"] >= 1
    assert metrics["shield_coverage"] >= 0.70
    assert metrics["platform_coverage"] >= 0.62
    assert metrics["interior_basin_coverage"] >= 0.60
    assert metrics["basin_landform_has_basin_parent"]
    assert metrics["basin_landform_has_platform_parent"]
    assert metrics["platform_landform_has_platform_parent"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p48_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p49_mature_platform_subsidence_coupling_microbenchmarks_pass(tmp_path):
    summary = run_suite("P49", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p49.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["mature_platform_subsidence_coupling"]
    assert summary["acceptance"]["subsidence_is_tectonically_localized"]
    assert summary["acceptance"]["terrain_relief_and_sediment_respond"]
    assert summary["acceptance"]["land_mask_preserved"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P49.mature_platform_subsidence_coupling"}
    metrics = benchmarks["P49.mature_platform_subsidence_coupling"]["metrics"]
    assert metrics["subsiding_basin_subsidence_mean"] > 0.46
    assert metrics["shield_subsidence_mean"] < 0.05
    assert metrics["active_suture_subsidence_mean"] < 0.08
    assert metrics["tectonic_basin_object_count"] >= 1
    assert metrics["tectonic_basin_object_coverage"] >= 0.58
    assert metrics["detail_after_basin_fraction"] >= 0.70
    assert metrics["basin_relief_delta_mean_m"] < -85.0
    assert metrics["basin_sediment_mean_m"] > metrics["quiet_platform_sediment_mean_m"] + 550.0
    assert metrics["land_mask_unchanged"]
    assert metrics["basin_landform_has_tectonic_basin_parent"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p49_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p50_planform_rebalance_microbenchmarks_pass(tmp_path):
    summary = run_suite("P50", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p50.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["parent_supported_planform_rebalance"]
    assert summary["acceptance"]["reduces_unstable_ribbon_without_area_drift"]
    assert summary["acceptance"]["inherits_from_broad_parent_not_arc"]
    assert summary["acceptance"]["does_not_bridge_gateway"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P50.parent_supported_planform_rebalance",
        "P50.no_cross_component_gateway_fill",
    }
    rebalance = benchmarks["P50.parent_supported_planform_rebalance"]["metrics"]
    assert rebalance["recycled_cells"] >= 3
    assert rebalance["filled_cells"] == rebalance["recycled_cells"]
    assert rebalance["after_ribbon_fraction"] <= rebalance["before_ribbon_fraction"] - 0.004
    assert rebalance["added_arc_or_suture_cells"] == 0
    assert rebalance["area_relative_error"] < 0.002

    gateway = benchmarks["P50.no_cross_component_gateway_fill"]["metrics"]
    assert gateway["gateway_core_new_land_cells"] == 0
    assert gateway["gateway_core_land_cells"] == 0
    assert gateway["after_left_right_component_count"] >= 2

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p50_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p51_generated_world_trajectory_guard_pass(tmp_path):
    summary = run_suite("P51", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p51.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["generated_world_trajectory_guard"]
    assert summary["acceptance"]["land_topology_within_guard_envelope"]
    assert summary["acceptance"]["mature_province_coverage_preserved"]
    assert summary["acceptance"]["p50_remains_conservative"]
    assert summary["acceptance"]["residual_ribbon_gap_visible"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P51.generated_world_trajectory_guard"}
    metrics = benchmarks["P51.generated_world_trajectory_guard"]["metrics"]
    assert 0.20 <= metrics["land_fraction"] <= 0.34
    assert metrics["land_component_count"] <= 8
    assert metrics["land_ribbon_fraction"] <= 0.72
    assert metrics["stable_craton_continental_fraction"] >= 0.12
    assert metrics["ancient_continental_fraction"] >= 0.20
    assert metrics["p62_supported_cratonic_lithosphere_fraction"] >= 0.16
    assert (
        metrics["p62_supported_cratonic_lithosphere_fraction"]
        >= metrics["stable_craton_continental_fraction"]
    )
    assert metrics["p48_mature_cratonic_share"] >= 0.16
    assert metrics["p48_shield_share"] >= 0.07
    assert metrics["p48_platform_share"] >= 0.05
    assert metrics["p50_planform_area_fraction"] <= 0.0025
    assert metrics["cumulative_p50_planform_area_fraction"] <= 0.0075
    assert metrics["residual_ribbon_gap_visible"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p51_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p52_deterministic_broad_block_initial_cargo_pass(tmp_path):
    summary = run_suite("P52", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p52.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["deterministic_broad_block_initial_cargo"]
    assert summary["acceptance"]["initial_cargo_rng_independent"]
    assert summary["acceptance"]["initial_broad_blocks"]
    assert summary["acceptance"]["initial_cratonic_kernels"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P52.deterministic_broad_block_initial_cargo"}
    metrics = benchmarks["P52.deterministic_broad_block_initial_cargo"]["metrics"]
    assert metrics["rng_independent_cargo"]
    assert metrics["nuclei_count"] == 4.0
    assert metrics["component_count"] >= 4.0
    assert 0.275 <= metrics["area_fraction"] <= 0.305
    assert 0.20 <= metrics["largest_component_fraction"] <= 0.38
    assert metrics["width_p50_steps"] >= 3.0
    assert metrics["width_p90_steps"] >= 5.0
    assert metrics["ribbon_fraction"] <= 0.32
    assert 0.52 <= metrics["craton_fraction"] <= 0.68

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p52_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p53_deterministic_initial_plate_motions_pass(tmp_path):
    summary = run_suite("P53", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p53.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["initial_plate_motion_rng_independent"]
    assert summary["acceptance"]["initial_plate_motion_sources_process_based"]
    assert summary["acceptance"]["initial_plate_motion_rates_physical"]
    assert summary["acceptance"]["initial_boundary_proxy_complete"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P53.deterministic_initial_plate_motions"}
    metrics = benchmarks["P53.deterministic_initial_plate_motions"]["metrics"]
    assert metrics["rng_independent_initial_motion"]
    assert metrics["active_plate_count"] >= 8.0
    assert metrics["random_motion_source_count"] == 0.0
    assert metrics["torque_driven_fraction"] >= 0.85
    assert metrics["fallback_fraction"] <= 0.15
    assert 0.0030 <= metrics["min_rate_rad_per_myr"]
    assert metrics["max_rate_rad_per_myr"] <= 0.0140
    assert metrics["ridge_proxy_cells"] > 0.0
    assert metrics["trench_proxy_cells"] > 0.0
    assert metrics["collision_proxy_cells"] > 0.0
    assert metrics["transform_proxy_cells"] > 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p53_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p54_cratonic_province_cargo_inheritance_pass(tmp_path):
    summary = run_suite("P54", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p54.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["cratonic_province_cargo_inheritance"]
    assert summary["acceptance"]["inherited_province_objects_restored"]
    assert summary["acceptance"]["province_memory_roundtrip"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P54.cratonic_province_cargo_inheritance"}
    metrics = benchmarks["P54.cratonic_province_cargo_inheritance"]["metrics"]
    assert metrics["shield_with_memory_share"] > metrics["shield_without_memory_share"]
    assert metrics["platform_with_memory_share"] > metrics["platform_without_memory_share"]
    assert metrics["basin_with_memory_share"] >= 0.02
    assert metrics["updated_shield_memory_share"] >= metrics["shield_with_memory_share"] * 0.95
    assert metrics["telemetry_total_memory_share"] >= 0.25

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p54_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p55_inherited_province_physical_cargo_pass(tmp_path):
    summary = run_suite("P55", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p55.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["inherited_province_physical_cargo"]
    assert summary["acceptance"]["physical_cargo_restored"]
    assert summary["acceptance"]["active_rework_respected"]
    assert summary["acceptance"]["province_objects_restored_after_physical_cargo"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P55.inherited_province_physical_cargo"}
    metrics = benchmarks["P55.inherited_province_physical_cargo"]["metrics"]
    assert metrics["shield_preserved_share"] >= 0.12
    assert metrics["platform_preserved_share"] >= 0.12
    assert metrics["stable_ancient_share"] >= 0.24
    assert metrics["blocked_expected_cells"] > 0.0
    assert metrics["blocked_young_fraction"] >= 0.90
    assert metrics["shield_object_share_after_physical_cargo"] >= 0.12
    assert metrics["platform_object_share_after_physical_cargo"] >= 0.18

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p55_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p56_persistent_basement_cargo_archive_pass(tmp_path):
    summary = run_suite("P56", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p56.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["persistent_basement_cargo_archive"]
    assert summary["acceptance"]["archive_restored"]
    assert summary["acceptance"]["active_respected"]
    assert summary["acceptance"]["true_platform_objects"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P56.persistent_basement_cargo_archive"}
    metrics = benchmarks["P56.persistent_basement_cargo_archive"]["metrics"]
    assert metrics["archive_shield_share"] >= 0.12
    assert metrics["archive_platform_share"] >= 0.12
    assert metrics["restored_shield_share"] >= 0.12
    assert metrics["restored_platform_share"] >= 0.12
    assert metrics["stable_ancient_share"] >= 0.24
    assert metrics["active_unrestored_fraction"] >= 0.90
    assert metrics["platform_object_share"] >= 0.12

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p56_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p57_conservation_basement_cargo_continuity_pass(tmp_path):
    summary = run_suite("P57", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p57.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["conservation_basement_cargo_continuity"]
    assert summary["acceptance"]["parented_cargo_preserved"]
    assert summary["acceptance"]["maintenance_gap_repaired"]
    assert summary["acceptance"]["true_ocean_recycled"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P57.conservation_basement_cargo_continuity"}
    metrics = benchmarks["P57.conservation_basement_cargo_continuity"]["metrics"]
    assert metrics["platform_gain_cells"] > 0
    assert metrics["parented_old_basement_fraction"] >= 0.95
    assert metrics["helper_repaired_old_fraction"] >= 0.95
    assert metrics["stale_ocean_cleared_fraction"] >= 0.95

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p57_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p58_basement_survival_ledger_pass(tmp_path):
    summary = run_suite("P58", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p58.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["basement_survival_ledger"]
    assert summary["acceptance"]["raster_attribution"]
    assert summary["acceptance"]["reorg_attribution"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P58.basement_survival_ledger"}
    metrics = benchmarks["P58.basement_survival_ledger"]["metrics"]
    assert metrics["raster_unrepresented_old_area_fraction"] > 0.0
    assert metrics["raster_duplicate_old_area_fraction"] > 0.0
    assert metrics["ridge_cleared_old_area_fraction"] > 0.0
    assert metrics["arc_cleared_old_area_fraction"] > 0.0
    assert metrics["repair_gain_old_area_fraction"] > 0.0
    assert metrics["reorg_old_plate_changed_area_fraction"] > 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p58_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p59_old_basement_raster_source_coverage_pass(tmp_path):
    summary = run_suite("P59", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p59.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["old_basement_raster_source_coverage"]
    assert summary["acceptance"]["coverage_restored"]
    assert summary["acceptance"]["ledger_improved"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P59.old_basement_raster_source_coverage"}
    metrics = benchmarks["P59.old_basement_raster_source_coverage"]["metrics"]
    assert metrics["before_missing_old_source_cells"] >= 1
    assert metrics["after_missing_old_source_cells"] == 0
    assert metrics["crust_primary_unchanged"]
    assert metrics["target_reassigned_to_missing_source"]
    assert metrics["old_a_still_represented"]
    assert metrics["p58_unrepresented_after_p59"] == 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p59_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p60_basement_to_province_expression_pass(tmp_path):
    summary = run_suite("P60", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p60.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["basement_to_province_expression"]
    assert summary["acceptance"]["objects_expressed"]
    assert summary["acceptance"]["terrain_expressed"]
    assert summary["acceptance"]["archive_only"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P60.basement_to_province_expression"}
    metrics = benchmarks["P60.basement_to_province_expression"]["metrics"]
    assert metrics["shield_object_coverage"] >= 0.58
    assert metrics["platform_object_coverage"] >= 0.48
    assert metrics["basin_object_coverage"] >= 0.48
    assert metrics["active_old_object_fraction"] <= 0.05
    assert metrics["shield_detail_fraction"] >= 0.70
    assert metrics["platform_detail_fraction"] >= 0.58
    assert metrics["basin_detail_fraction"] >= 0.55
    assert metrics["active_old_detail_shield_fraction"] <= 0.05
    assert metrics["current_crust_arrays_unchanged"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p60_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p61_basement_current_state_continuity_pass(tmp_path):
    summary = run_suite("P61", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p61.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["basement_current_state_continuity"]
    assert summary["acceptance"]["shield_recovered"]
    assert summary["acceptance"]["platform_recovered_not_shield"]
    assert summary["acceptance"]["active_mobile_belt_blocked"]
    assert summary["acceptance"]["localized_and_archive_preserved"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P61.basement_current_state_continuity"}
    metrics = benchmarks["P61.basement_current_state_continuity"]["metrics"]
    assert metrics["shield_origin_craton_fraction"] >= 0.80
    assert metrics["shield_current_stable_fraction"] >= 0.80
    assert metrics["platform_current_stable_fraction"] >= 0.70
    assert metrics["passive_margin_current_stable_fraction"] >= 0.70
    assert metrics["platform_origin_craton_fraction"] <= 0.05
    assert metrics["passive_margin_origin_craton_fraction"] <= 0.05
    assert metrics["active_suture_changed_fraction"] <= 0.05
    assert metrics["background_changed_fraction"] <= 0.05
    assert metrics["stable_craton_fraction_after"] > metrics["stable_craton_fraction_before"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p61_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p62_nonmutating_cratonic_lithosphere_support_pass(tmp_path):
    summary = run_suite("P62", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p62.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["nonmutating_cratonic_lithosphere_support"]
    assert summary["acceptance"]["recognizes_current_and_archive"]
    assert summary["acceptance"]["active_mobile_belt_blocked"]
    assert summary["acceptance"]["crust_arrays_unchanged"]
    assert summary["acceptance"]["telemetry_available"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P62.nonmutating_cratonic_lithosphere_support"}
    metrics = benchmarks["P62.nonmutating_cratonic_lithosphere_support"]["metrics"]
    assert metrics["current_core_support_fraction"] >= 0.95
    assert metrics["archive_shield_support_fraction"] >= 0.95
    assert metrics["covered_platform_support_fraction"] >= 0.85
    assert metrics["sag_basin_support_fraction"] >= 0.80
    assert metrics["passive_margin_support_fraction"] >= 0.80
    assert metrics["active_suture_support_fraction"] <= 0.05
    assert metrics["background_support_fraction"] <= 0.05
    assert metrics["telemetry_supported_fraction"] > metrics["telemetry_current_stable_fraction"]
    assert metrics["telemetry_proxy_only_fraction"] > 0.12
    assert metrics["crust_arrays_unchanged"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p62_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p63_flagged_full_history_support_guard_pass(tmp_path):
    summary = run_suite("P63", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p63.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["flagged_full_history_support_guard"]
    assert summary["acceptance"]["topology_preserved"]
    assert summary["acceptance"]["archive_preserved"]
    assert summary["acceptance"]["province_preserved"]
    assert summary["acceptance"]["p62_support_covers_missing_signal"]
    assert summary["acceptance"]["production_default_still_blocked"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P63.flagged_full_history_support_guard"}
    metrics = benchmarks["P63.flagged_full_history_support_guard"]["metrics"]
    assert metrics["old_basement_archive_fraction"] >= 0.214
    assert metrics["p48_platform_share"] >= 0.095
    assert metrics["p48_mature_cratonic_share"] >= 0.153
    assert metrics["p62_supported_cratonic_lithosphere_fraction"] >= metrics["p60_old_candidate_share"]
    assert (
        metrics["p62_supported_cratonic_lithosphere_fraction"]
        >= metrics["stable_craton_continental_fraction"] + 0.04
    )
    assert metrics["p62_proxy_only_fraction"] >= 0.04
    assert metrics["p59_missing_old_source_cells_after"] == 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p63_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p64_support_guided_current_projection_pass(tmp_path):
    summary = run_suite("P64", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p64.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["support_guided_current_projection"]
    assert summary["acceptance"]["topology_preserved"]
    assert summary["acceptance"]["archive_preserved"]
    assert summary["acceptance"]["platform_preserved"]
    assert summary["acceptance"]["current_projection_improved"]
    assert summary["acceptance"]["support_bounded"]
    assert summary["acceptance"]["production_default_still_blocked"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P64.support_guided_current_projection"}
    metrics = benchmarks["P64.support_guided_current_projection"]["metrics"]
    assert metrics["stable_craton_continental_fraction"] >= 0.120
    assert metrics["stable_craton_continental_fraction"] >= (
        metrics["p63_baseline_stable_craton_continental_fraction"] + 0.035
    )
    assert metrics["stable_craton_continental_fraction"] <= (
        metrics["p62_supported_cratonic_lithosphere_fraction"] + 0.006
    )
    assert metrics["old_basement_archive_fraction"] >= 0.214
    assert metrics["p48_platform_share"] >= 0.095
    assert metrics["p64_platform_projected_share"] <= 0.010
    assert metrics["p64_projected_share"] > 0.030
    assert metrics["p59_missing_old_source_cells_after"] == 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p64_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p65_post_projection_cratonic_refresh_pass(tmp_path):
    summary = run_suite("P65", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p65.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["post_projection_cratonic_refresh"]
    assert summary["acceptance"]["topology_preserved"]
    assert summary["acceptance"]["archive_preserved"]
    assert summary["acceptance"]["ancient_gap_closed"]
    assert summary["acceptance"]["province_refresh_closed"]
    assert summary["acceptance"]["support_bounded"]
    assert summary["acceptance"]["refresh_bounded"]
    assert summary["acceptance"]["flagged_path_ready"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P65.post_projection_cratonic_refresh"}
    metrics = benchmarks["P65.post_projection_cratonic_refresh"]["metrics"]
    assert metrics["stable_craton_continental_fraction"] >= 0.120
    assert metrics["ancient_continental_fraction"] >= 0.200
    assert metrics["p48_mature_cratonic_share"] >= 0.160
    assert metrics["old_basement_archive_fraction"] >= 0.213
    assert metrics["p65_ancient_age_refreshed_share"] <= 0.012
    assert metrics["p65_refreshed_p48_mature_share"] >= 0.160
    assert metrics["p65_added_platform_object_share"] <= 0.012
    assert metrics["p65_added_platform_object_count"] >= 0.0
    assert metrics["p65_province_code_changed_share"] <= 0.030
    assert metrics["p59_missing_old_source_cells_after"] == 0.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p65_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p66_multiseed_flagged_cratonic_audit_pass(tmp_path):
    summary = run_suite("P66", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p66.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["multiseed_flagged_cratonic_audit"]
    assert summary["acceptance"]["audit_completed"]
    assert (
        summary["acceptance"]["promotion_ready"]
        or summary["acceptance"]["keeps_default_off_until_multiseed_ready"]
    )

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P66.multiseed_flagged_cratonic_audit"}
    metrics = benchmarks["P66.multiseed_flagged_cratonic_audit"]["metrics"]
    assert metrics["seed_count"] == 3
    assert metrics["p59_missing_old_source_cells_after_max"] == 0.0
    assert metrics["unique_trajectory_fingerprint_count"] >= 1
    assert metrics["archive_min"] >= 0.18
    assert metrics["land_component_max"] <= 12
    assert metrics["land_ribbon_max"] <= 0.80
    assert metrics["continental_ribbon_max"] <= 0.80
    assert metrics["promotion_ready"] == (
        metrics["continuity_ready_all_seeds"]
        and metrics["terrain_response_ready_all_seeds"]
        and not metrics["seed_invariant_flagged_path"]
    )
    assert metrics["promotion_ready"] != metrics["keeps_default_off_until_multiseed_ready"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p66_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p67_deterministic_physical_ensemble_audit_pass(tmp_path):
    summary = run_suite("P67", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p67.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["deterministic_physical_ensemble_audit"]
    assert summary["acceptance"]["audit_completed"]
    assert summary["acceptance"]["physical_ensemble_differentiated"]
    assert summary["acceptance"]["catastrophic_topology_free_all_members"]
    assert (
        summary["acceptance"]["promotion_ready"]
        or summary["acceptance"]["keeps_default_off_until_physical_ensemble_ready"]
    )

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P67.deterministic_physical_ensemble_audit"}
    metrics = benchmarks["P67.deterministic_physical_ensemble_audit"]["metrics"]
    assert metrics["member_count"] == 3
    assert metrics["unique_trajectory_fingerprint_count"] >= 2
    assert metrics["p59_missing_old_source_cells_after_max"] == 0.0
    assert metrics["archive_min"] >= 0.140
    assert metrics["land_fraction_min"] >= 0.10
    assert metrics["land_fraction_max"] <= 0.45
    assert metrics["land_component_max"] <= 32
    assert metrics["land_ribbon_max"] <= 0.82
    assert metrics["continental_ribbon_max"] <= 0.82
    assert metrics["catastrophic_topology_free_count"] == metrics["member_count"]
    assert metrics["cratonic_balance_ok_count"] < metrics["member_count"]
    assert metrics["failed_member_count"] >= 1
    assert metrics["promotion_ready"] == (
        metrics["audit_completed"]
        and metrics["continuity_ready_all_members"]
        and metrics["terrain_response_ready_all_members"]
        and metrics["catastrophic_topology_free_all_members"]
    )
    assert (
        metrics["promotion_ready"]
        != metrics["keeps_default_off_until_physical_ensemble_ready"]
    )

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p67_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p68_cratonic_balance_physical_ensemble_pass(tmp_path):
    summary = run_suite("P68", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p68.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["cratonic_balance_physical_ensemble"]
    assert summary["acceptance"]["audit_completed"]
    assert summary["acceptance"]["physical_ensemble_differentiated"]
    assert summary["acceptance"]["continuity_ready_all_members"]
    assert summary["acceptance"]["cratonic_balance_ready_all_members"]
    assert summary["acceptance"]["terrain_response_ready_all_members"]
    assert summary["acceptance"]["promotion_ready"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P68.cratonic_balance_physical_ensemble"}
    metrics = benchmarks["P68.cratonic_balance_physical_ensemble"]["metrics"]
    assert metrics["member_count"] == 3
    assert metrics["unique_trajectory_fingerprint_count"] >= 2
    assert metrics["stable_min"] >= 0.120
    assert metrics["ancient_min"] >= 0.200
    assert metrics["mature_province_min"] >= 0.160
    assert metrics["archive_min"] >= 0.195
    assert metrics["ancient_max"] <= 0.650
    assert metrics["archive_max"] <= 0.700
    assert metrics["mature_province_max"] <= 0.600
    assert metrics["platform_min"] >= 0.030
    assert metrics["land_component_max"] <= 28
    assert metrics["land_ribbon_max"] <= 0.76
    assert metrics["continental_ribbon_max"] <= 0.76
    assert metrics["p59_missing_old_source_cells_after_max"] == 0.0
    assert metrics["failed_member_count"] == 0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p68_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_p70_reference_corpus_metric_scaffold_pass(tmp_path):
    summary = run_suite("P70", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p70.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["reference_corpus_metric_scaffold"]
    assert summary["acceptance"]["source_inventory_ready"]
    assert summary["acceptance"]["metric_schema_ready"]
    assert summary["acceptance"]["p69_baseline_mapped"]
    assert summary["acceptance"]["reference_data_download_required"]
    assert summary["acceptance"]["production_province_graph_missing"]

    scaffold = summary["reference_scaffold"]
    assert scaffold["schema"] == "aevum.physiographic_reference.v1"
    assert scaffold["status"] == "scaffold_ready"
    assert len(scaffold["phases"]) == 8
    assert scaffold["source_inventory"]["source_count"] >= 14
    assert scaffold["metric_schema"]["metric_count"] >= 30
    assert scaffold["source_inventory"]["acceptance"]["covers_all_declared_phases"]
    assert scaffold["metric_schema"]["acceptance"]["has_required_metric_groups"]
    assert scaffold["p69_baseline_mapping"]["acceptance"]["explicitly_records_missing_province_graph"]
    assert not scaffold["p69_baseline_mapping"]["unmapped_out_of_envelope_metrics"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P70.reference_corpus_metric_scaffold"}
    metrics = benchmarks["P70.reference_corpus_metric_scaffold"]["metrics"]
    assert metrics["phase_count"] == 8
    assert metrics["source_category_count"] >= 4
    assert metrics["metric_group_count"] >= 6
    assert metrics["p69_out_of_envelope_count"] == 6
    assert metrics["future_architecture_metric_count"] >= 5

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p70_microbenchmarks.csv").exists()


def test_p71_province_graph_fixture_suite_pass(tmp_path):
    summary = run_suite("P71", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p71.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["province_graph_fixture_suite"]
    assert summary["acceptance"]["expected_fixture_suite_complete"]
    assert summary["acceptance"]["all_fixtures_pass"]
    assert summary["acceptance"]["required_province_classes_covered"]
    assert summary["acceptance"]["required_parent_processes_covered"]
    assert summary["acceptance"]["parent_process_links_present"]
    assert summary["acceptance"]["province_ids_contiguous_where_expected"]
    assert summary["acceptance"]["no_random_texture_only_provinces"]
    assert summary["acceptance"]["expected_elevation_ordering"]
    assert summary["acceptance"]["expected_sediment_ordering"]
    assert summary["acceptance"]["production_world_integration_pending"]

    fixture_suite = summary["province_graph_fixtures"]
    assert fixture_suite["schema"] == "aevum.province_graph_fixtures.v1"
    assert fixture_suite["status"] == "fixture_suite_ready"
    assert fixture_suite["fixture_count"] == 7
    assert fixture_suite["province_class_count"] >= 13
    assert fixture_suite["parent_process_count"] >= 13
    assert len(fixture_suite["fixtures"]) == 7
    assert {
        fixture["name"] for fixture in fixture_suite["fixtures"]
    } == {
        "craton_platform_basin_fixture",
        "active_orogen_foreland_fixture",
        "old_suture_orogen_fixture",
        "rift_axis_shoulder_fixture",
        "passive_margin_lowland_fixture",
        "volcanic_lip_plateau_fixture",
        "multi_province_continent_fixture",
    }
    for fixture in fixture_suite["fixtures"]:
        assert fixture["passed"]
        assert fixture["acceptance"]["expected_classes_present"]
        assert fixture["acceptance"]["expected_edges_present"]
        assert fixture["acceptance"]["parent_process_links_present"]
        assert fixture["acceptance"]["province_ids_contiguous_where_expected"]
        assert fixture["acceptance"]["no_random_texture_only_provinces"]
        assert fixture["acceptance"]["expected_elevation_ordering"]
        assert fixture["acceptance"]["expected_sediment_ordering"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P71.province_graph_fixture_suite"}
    metrics = benchmarks["P71.province_graph_fixture_suite"]["metrics"]
    assert metrics["fixture_count"] == 7
    assert metrics["passed_fixture_count"] == 7
    assert metrics["province_class_count"] >= 13
    assert metrics["parent_process_count"] >= 13
    assert metrics["failed_fixture_count"] == 0
    assert metrics["random_texture_only_province_count"] == 0
    assert metrics["unparented_province_count"] == 0
    assert metrics["non_contiguous_province_count"] == 0
    assert metrics["missing_expected_class_count"] == 0
    assert metrics["missing_expected_edge_count"] == 0

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p71_microbenchmarks.csv").exists()


def test_p72_generated_world_province_diversity_gate_pass(tmp_path):
    summary = run_suite("P72", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p72.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["generated_world_province_diversity_gate"]
    assert summary["acceptance"]["major_continent_count_multiple"]
    assert summary["acceptance"]["all_major_components_have_three_province_classes"]
    assert summary["acceptance"]["largest_internal_province_fraction_capped"]
    assert summary["acceptance"]["basin_or_lowland_present_per_major"]
    assert summary["acceptance"]["active_highland_area_limited_per_major"]
    assert summary["acceptance"]["highlands_parented_per_major"]
    assert summary["acceptance"]["terrain_landform_parent_processes_present"]
    assert summary["acceptance"]["passive_margin_lowland_object_layer_present"]
    assert summary["acceptance"]["generated_world_province_signal_present"]
    assert summary["acceptance"]["passive_margin_lowland_detail_class_pending"]
    assert summary["acceptance"]["terrain_rewrite_active"]
    assert not summary["acceptance"]["terrain_rewrite_pending"]

    diversity = summary["province_diversity"]
    assert diversity["schema"] == "aevum.generated_world_province_diversity.v1"
    assert diversity["status"] == "generated_world_gate_ready"
    assert diversity["context"]["cells"] == 900
    assert diversity["metrics"]["major_component_count"] >= 2
    assert diversity["metrics"]["min_province_class_count_per_major"] >= 3
    assert diversity["metrics"]["max_largest_internal_province_fraction"] <= 0.74
    assert diversity["metrics"]["min_basin_or_lowland_share_per_major"] >= 0.18
    assert diversity["metrics"]["max_active_highland_or_plateau_fraction"] <= 0.45
    assert diversity["metrics"]["max_unparented_highland_fraction"] <= 0.15
    assert diversity["metrics"]["terrain_landform_kind_count"] >= 5
    assert diversity["metrics"]["terrain_landform_parent_process_count"] >= 5
    assert diversity["metrics"]["passive_margin_lowland_object_count"] > 0
    assert diversity["metrics"]["passive_margin_lowland_field_area_fraction"] > 0.0
    assert diversity["metrics"]["unparented_landform_object_count"] == 0
    assert len(diversity["major_components"]) == diversity["metrics"]["major_component_count"]
    assert not diversity["failing_major_components"]
    for component in diversity["major_components"]:
        assert component["acceptance"]["minimum_province_class_diversity"]
        assert component["acceptance"]["largest_internal_province_fraction_capped"]
        assert component["acceptance"]["basin_or_lowland_present"]
        assert component["acceptance"]["active_highland_area_limited"]
        assert component["acceptance"]["highlands_parented"]
        assert component["acceptance"]["landform_object_context_present"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P72.generated_world_province_diversity_gate"}
    metrics = benchmarks["P72.generated_world_province_diversity_gate"]["metrics"]
    assert metrics["cells"] == 900
    assert metrics["major_component_count"] >= 2
    assert metrics["failing_major_component_count"] == 0

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p72_microbenchmarks.csv").exists()


def test_p73_real_earth_case_study_calibration_pass(tmp_path):
    summary = run_suite("P73", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p73.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["real_earth_case_study_calibration"]
    assert summary["acceptance"]["required_case_suite_complete"]
    assert summary["acceptance"]["all_cases_have_reference_sketch"]
    assert summary["acceptance"]["all_cases_have_metric_envelope"]
    assert summary["acceptance"]["all_case_sources_known"]
    assert summary["acceptance"]["all_adjacencies_valid"]
    assert summary["acceptance"]["feature_class_not_exact_geography_policy"]
    assert summary["acceptance"]["failure_categories_complete"]
    assert summary["acceptance"]["required_feature_classes_covered"]
    assert summary["acceptance"]["required_parent_processes_covered"]
    assert summary["acceptance"]["production_terrain_rewrite_pending"]
    assert summary["acceptance"]["reference_raster_extraction_pending"]

    cases = summary["case_studies"]
    assert cases["schema"] == "aevum.earth_case_studies.v1"
    assert cases["status"] == "case_study_calibration_ready"
    assert cases["case_count"] == 5
    assert cases["province_count"] >= 25
    assert cases["province_class_count"] >= 8
    assert cases["parent_process_count"] >= 8
    assert cases["adjacency_edge_count"] >= 20
    assert set(cases["failure_categories"]) == {
        "missing_process",
        "wrong_amplitude",
        "wrong_scale",
        "wrong_adjacency",
        "compiler_rendering_mismatch",
    }
    assert {case["id"] for case in cases["cases"]} == {
        "north_america",
        "south_america",
        "africa",
        "eurasia",
        "australia",
    }
    for case in cases["cases"]:
        assert case["province_count"] >= 5
        assert case["province_class_count"] >= 4
        assert case["adjacency_edge_count"] >= 4
        assert not case["missing_source_ids"]
        assert not case["invalid_adjacency_edges"]
        assert not case["missing_envelope_keys"]
        assert not case["invalid_envelope_keys"]
        assert case["reproduction_policy"]["feature_class_required"]
        assert not case["reproduction_policy"]["exact_geography_required"]
        assert set(case["failure_categories"]) == {
            "missing_process",
            "wrong_amplitude",
            "wrong_scale",
            "wrong_adjacency",
            "compiler_rendering_mismatch",
        }
        assert all(case["acceptance"].values())

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P73.real_earth_case_study_calibration"}
    metrics = benchmarks["P73.real_earth_case_study_calibration"]["metrics"]
    assert metrics["case_count"] == 5
    assert metrics["province_count"] >= 25
    assert metrics["failure_category_count"] == 5
    assert metrics["failed_case_count"] == 0
    assert metrics["cases_with_unknown_sources"] == 0
    assert metrics["cases_with_invalid_adjacencies"] == 0
    assert metrics["cases_with_missing_envelopes"] == 0

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p73_microbenchmarks.csv").exists()


def test_p74_terrain_coupling_rewrite_pass(tmp_path):
    summary = run_suite("P74", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p74.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["terrain_coupling_rewrite"]
    assert summary["acceptance"]["p74_template_response_active"]
    assert summary["acceptance"]["passive_margin_lowland_object_layer_active"]
    assert summary["acceptance"]["province_diversity_gate_ready"]
    assert summary["acceptance"]["lowland_and_basin_modes_present"]
    assert summary["acceptance"]["rift_and_foreland_lows_present"]
    assert summary["acceptance"]["old_orogen_expression_present"]
    assert summary["acceptance"]["highlands_parented_or_limited"]
    assert summary["acceptance"]["bathymetry_regression_free"]
    assert summary["acceptance"]["compiler_regression_free"]
    assert summary["acceptance"]["release_promotion_audit_pending"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P74.terrain_coupling_rewrite"}
    metrics = benchmarks["P74.terrain_coupling_rewrite"]["metrics"]
    assert metrics["cells"] == 900
    assert metrics["p74_template_area_fraction"] > 0.0
    assert metrics["p74_template_mean_abs_m"] > 25.0
    assert metrics["p74_land_preserved"]
    assert metrics["p74_passive_margin_lowland_area_fraction"] > 0.0
    assert metrics["passive_margin_lowland_field_area_fraction"] > 0.0
    assert metrics["passive_margin_lowland_object_count"] > 0
    assert metrics["passive_margin_lowland_mean_elevation_m"] < 500.0
    assert metrics["platform_object_count"] > 0
    assert metrics["interior_basin_object_count"] > 0
    assert metrics["rift_basin_object_count"] > 0
    assert metrics["foreland_basin_object_count"] > 0
    assert metrics["old_subdued_orogen_object_count"] > 0
    assert metrics["province_diversity_status"] == "generated_world_gate_ready"
    assert metrics["min_province_class_count_per_major"] >= 3
    assert metrics["max_largest_internal_province_fraction"] <= 0.74
    assert metrics["max_unparented_highland_fraction"] <= 0.15
    assert metrics["bathymetry_regression_free"]
    assert metrics["compiler_regression_free"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p74_microbenchmarks.csv").exists()


def test_p75_release_and_promotion_audit_pass(tmp_path):
    summary = run_suite("P75", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p75.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["release_and_promotion_audit"]
    assert summary["acceptance"]["p69_highres_assets_available"]
    assert summary["acceptance"]["p69_highres_ready_all_members"]
    assert summary["acceptance"]["p70_p74_gates_pass"]
    assert summary["acceptance"]["legacy_or_replacement_gates_pass"]
    assert summary["acceptance"]["promotion_decision_recorded"]
    assert not summary["acceptance"]["promotion_ready"]
    assert summary["acceptance"]["keeps_default_off_until_stage_a_reference_metrics"]
    assert summary["acceptance"]["stage_a_next_step_ready"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P75.release_and_promotion_audit"}
    metrics = benchmarks["P75.release_and_promotion_audit"]["metrics"]
    assert metrics["p69_archived_summary_count"] >= 1
    assert metrics["p69_selected_summary_status"] == "pass"
    assert metrics["p69_cells"] == 8000
    assert metrics["p69_member_count"] == 3
    assert metrics["p69_asset_set_complete_count"] == 3
    assert metrics["p69_highres_member_ready_count"] == 3
    assert metrics["actual_key_asset_files"] == metrics["expected_key_asset_files"]
    assert metrics["contact_sheet_present_count"] >= 1
    assert metrics["p70_p74_gate_pass_count"] == 5
    assert set(metrics["p70_p74_gate_statuses"]) == {"P70", "P71", "P72", "P73", "P74"}
    assert all(status == "pass" for status in metrics["p70_p74_gate_statuses"].values())
    assert metrics["legacy_gate_pass_count"] == 4
    assert set(metrics["legacy_gate_statuses"]) == {"P29", "P48", "P49", "P68"}
    assert all(status == "pass" for status in metrics["legacy_gate_statuses"].values())
    assert metrics["audit_completed"]
    assert metrics["promotion_decision_recorded"]
    assert not metrics["promotion_ready"]
    assert metrics["keeps_default_off_until_stage_a_reference_metrics"]
    assert metrics["reference_data_pending"]
    assert metrics["production_province_graph_pending"]
    assert "reference_data_download_or_raster_extraction_pending" in metrics["promotion_blockers"]
    assert "first_class_production_province_graph_pending" in metrics["promotion_blockers"]
    assert metrics["next_recommended_suites"] == [
        "P76.reference_source_ledger_schema",
        "P77.real_earth_hypsometry_extraction",
    ]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p75_microbenchmarks.csv").exists()


def test_p76_reference_source_ledger_schema_pass(tmp_path):
    summary = run_suite("P76", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p76.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["reference_source_ledger_schema"]
    assert summary["acceptance"]["source_ids_match_p70_inventory"]
    assert summary["acceptance"]["all_required_fields_present"]
    assert summary["acceptance"]["license_status_explicit"]
    assert summary["acceptance"]["source_version_note_explicit"]
    assert summary["acceptance"]["projection_resolution_note_explicit"]
    assert summary["acceptance"]["raw_data_policy_explicit"]
    assert summary["acceptance"]["acquisition_separated_from_benchmarks"]
    assert summary["acceptance"]["derived_metric_targets_present"]
    assert summary["acceptance"]["no_external_raw_data_required"]
    assert summary["acceptance"]["p77_hypsometry_extraction_pending"]

    ledger = summary["source_ledger"]
    assert ledger["schema"] == "aevum.reference_source_ledger.v1"
    assert ledger["status"] == "source_ledger_schema_ready"
    assert ledger["source_count"] >= 14
    assert ledger["external_source_count"] >= 10
    assert ledger["internal_source_count"] >= 1
    assert ledger["category_count"] >= 4
    assert len(ledger["required_fields"]) >= 12
    assert not ledger["missing_required_fields"]
    assert not ledger["invalid_phase_refs"]
    assert not ledger["duplicate_source_ids"]
    assert all(ledger["acceptance"].values())
    for entry in ledger["entries"]:
        assert set(ledger["required_fields"]).issubset(entry)
        assert entry["source_id"]
        assert entry["source_version_note"]
        assert entry["license_status"]
        assert entry["projection_resolution_note"]
        assert entry["raw_data_policy"]
        assert entry["derived_metric_targets"]
        if entry["source_kind"] == "external_reference":
            assert entry["acquisition_status"] == "not_downloaded_by_default"
            assert "outside git" in entry["local_storage_policy"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P76.reference_source_ledger_schema"}
    metrics = benchmarks["P76.reference_source_ledger_schema"]["metrics"]
    assert metrics["source_count"] == ledger["source_count"]
    assert metrics["missing_required_field_source_count"] == 0
    assert metrics["invalid_phase_ref_source_count"] == 0
    assert metrics["duplicate_source_id_count"] == 0

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p76_microbenchmarks.csv").exists()


def test_p77_real_earth_hypsometry_extraction_pass(tmp_path):
    summary = run_suite("P77", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p77.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["real_earth_hypsometry_extraction"]
    assert summary["acceptance"]["fixture_schema_ready"]
    assert summary["acceptance"]["source_ids_present"]
    assert summary["acceptance"]["raw_raster_not_stored"]
    assert summary["acceptance"]["small_derived_fixture"]
    assert summary["acceptance"]["required_metrics_present"]
    assert summary["acceptance"]["bin_area_fractions_sum_to_one"]
    assert summary["acceptance"]["land_ocean_bins_match_metrics"]
    assert summary["acceptance"]["envelope_checks_pass"]
    assert summary["acceptance"]["direct_raster_extraction_marked_pending"]
    assert summary["acceptance"]["p78_generated_hypsometry_envelope_pending"]

    fixture = summary["hypsometry_fixture"]
    assert fixture["schema"] == "aevum.real_earth_hypsometry_extraction.v1"
    assert fixture["status"] == "hypsometry_fixture_ready"
    assert "NOAA_ETOPO_2022" in fixture["source_ids"]
    assert "GEBCO_GRIDDED_BATHYMETRY" in fixture["source_ids"]
    assert fixture["bin_count"] >= 10
    assert abs(fixture["bin_area_total"] - 1.0) <= 1.0e-9
    assert abs(
        fixture["land_bin_area_total"] - fixture["metrics"]["land_fraction"]
    ) <= 1.0e-9
    assert abs(
        fixture["ocean_bin_area_total"] - fixture["metrics"]["ocean_fraction"]
    ) <= 1.0e-9
    assert not fixture["missing_required_metrics"]
    assert all(fixture["acceptance"].values())
    assert not fixture["quality_flags"]["direct_raster_extraction"]
    assert fixture["quality_flags"]["small_derived_fixture"]
    assert fixture["quality_flags"]["requires_future_regeneration_from_raw_raster"]
    assert 0.25 <= fixture["metrics"]["land_fraction"] <= 0.33
    assert 400.0 <= fixture["metrics"]["land_elevation_mean_m"] <= 1200.0
    assert 1500.0 <= fixture["metrics"]["land_elevation_p95_m"] <= 5200.0
    assert 0.04 <= fixture["metrics"]["shelf_fraction_of_ocean"] <= 0.14
    assert 0.30 <= fixture["metrics"]["abyss_fraction_of_ocean"] <= 0.65
    assert all(check["in_envelope"] for check in fixture["envelope_checks"].values())

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P77.real_earth_hypsometry_extraction"}
    metrics = benchmarks["P77.real_earth_hypsometry_extraction"]["metrics"]
    assert metrics["bin_count"] == fixture["bin_count"]
    assert metrics["missing_required_metric_count"] == 0
    assert metrics["envelope_check_pass_count"] == metrics["envelope_check_count"]
    assert metrics["raw_raster_not_stored"]
    assert metrics["direct_raster_extraction_marked_pending"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p77_microbenchmarks.csv").exists()


def test_p78_generated_hypsometry_envelope_pass(tmp_path):
    summary = run_suite("P78", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p78.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["generated_hypsometry_envelope"]
    assert summary["acceptance"]["p77_fixture_ready"]
    assert summary["acceptance"]["core_hypsometry_envelope_pass"]
    assert summary["acceptance"]["land_fraction_close_to_fixture"]
    assert summary["acceptance"]["lowland_ordering_plausible"]
    assert summary["acceptance"]["current_expected_residuals_recorded"]
    assert summary["acceptance"]["archived_highres_evidence_available"]
    assert summary["acceptance"]["promotion_calibration_still_required"]
    assert summary["acceptance"]["p79_province_reference_graph_extraction_pending"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P78.generated_hypsometry_envelope"}
    metrics = benchmarks["P78.generated_hypsometry_envelope"]["metrics"]
    assert metrics["cells"] == 900
    assert metrics["fixture_status"] == "hypsometry_fixture_ready"
    assert metrics["fixture_comparable_metric_count"] >= 12
    assert 0.15 <= metrics["generated_land_fraction"] <= 0.35
    assert 400.0 <= metrics["generated_land_elevation_mean_m"] <= 1200.0
    assert 1500.0 <= metrics["generated_land_elevation_p95_m"] <= 5200.0
    assert 0.0 <= metrics["generated_high_land_fraction_gt2500m"] <= 0.18
    assert 0.04 <= metrics["generated_shelf_fraction_of_ocean"] <= 0.14
    assert 0.10 <= metrics["generated_slope_rise_fraction_of_ocean"] <= 0.35
    assert 0.30 <= metrics["generated_abyss_fraction_of_ocean"] <= 0.65
    assert 1200.0 <= metrics["generated_shelf_to_abyss_depth_delta_m"] <= 5200.0
    assert metrics["core_hypsometry_envelope_pass"]
    assert metrics["land_fraction_close_to_fixture"]
    assert metrics["lowland_ordering_plausible"]
    assert metrics["current_out_of_envelope_count"] > 0
    assert metrics["current_expected_residuals_recorded"]
    assert "land_fraction" in metrics["current_out_of_envelope"]
    assert metrics["archived_highres_available"]
    assert metrics["archived_highres_status"] == "pass"
    assert metrics["archived_highres_cells"] == 8000
    assert metrics["archived_highres_member_ready_count"] == metrics["archived_highres_member_count"] == 3
    assert metrics["archived_highres_out_of_envelope_max"] <= 6
    assert metrics["promotion_calibration_still_required"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p78_microbenchmarks.csv").exists()


def test_p79_province_reference_graph_extraction_pass(tmp_path):
    summary = run_suite("P79", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p79.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["province_reference_graph_extraction"]
    assert summary["acceptance"]["case_study_calibration_ready"]
    assert summary["acceptance"]["required_feature_classes_covered"]
    assert summary["acceptance"]["required_parent_processes_covered"]
    assert summary["acceptance"]["required_class_edges_covered"]
    assert summary["acceptance"]["all_case_graphs_connected"]
    assert summary["acceptance"]["all_case_edges_valid"]
    assert summary["acceptance"]["no_duplicate_case_edges"]
    assert summary["acceptance"]["no_isolated_provinces"]
    assert summary["acceptance"]["small_derived_fixture"]
    assert summary["acceptance"]["raw_vectors_not_stored"]
    assert summary["acceptance"]["direct_vector_extraction_marked_pending"]
    assert summary["acceptance"]["p80_generated_province_graph_reference_comparison_pending"]

    graph = summary["province_reference_graph"]
    assert graph["schema"] == "aevum.province_reference_graph.v1"
    assert graph["status"] == "province_reference_graph_ready"
    assert graph["source_case_study_status"] == "case_study_calibration_ready"
    assert graph["case_count"] == 5
    assert graph["node_count"] >= 25
    assert graph["edge_count"] >= 20
    assert graph["class_count"] >= 8
    assert graph["parent_process_count"] >= 8
    assert graph["class_edge_count"] >= 9
    assert not graph["missing_required_feature_classes"]
    assert not graph["missing_required_parent_processes"]
    assert not graph["missing_required_class_edges"]
    assert all(graph["acceptance"].values())
    assert graph["extraction_policy"]["raw_vector_extraction_pending"]
    assert not graph["extraction_policy"]["exact_geography_required"]
    assert graph["extraction_policy"]["feature_class_required"]
    assert "active_orogen|foreland_basin" in graph["class_edge_counts"]
    assert "platform|shield" in graph["class_edge_counts"]
    assert "rift_system|volcanic_lip_plateau" in graph["class_edge_counts"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P79.province_reference_graph_extraction"}
    metrics = benchmarks["P79.province_reference_graph_extraction"]["metrics"]
    assert metrics["case_count"] == graph["case_count"]
    assert metrics["node_count"] == graph["node_count"]
    assert metrics["edge_count"] == graph["edge_count"]
    assert metrics["missing_required_feature_class_count"] == 0
    assert metrics["missing_required_parent_process_count"] == 0
    assert metrics["missing_required_class_edge_count"] == 0
    assert metrics["case_node_count_min"] >= 5
    assert metrics["case_edge_count_min"] >= 4
    assert metrics["case_largest_class_fraction_max"] <= 0.5
    assert metrics["raw_vectors_not_stored"]
    assert metrics["direct_vector_extraction_marked_pending"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p79_microbenchmarks.csv").exists()


def test_p80_generated_province_graph_reference_comparison_pass(tmp_path):
    summary = run_suite("P80", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p80.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["generated_province_graph_reference_comparison"]
    assert summary["acceptance"]["reference_graph_ready"]
    assert summary["acceptance"]["generated_diversity_ready"]
    assert summary["acceptance"]["required_feature_classes_covered"]
    assert summary["acceptance"]["required_parent_processes_covered"]
    assert summary["acceptance"]["unexpected_reference_class_gaps_absent"]
    assert summary["acceptance"]["unexpected_required_class_edge_gaps_absent"]
    assert summary["acceptance"]["expected_residuals_recorded"]
    assert summary["acceptance"]["major_components_multi_class"]
    assert summary["acceptance"]["major_components_have_core_context"]
    assert summary["acceptance"]["major_components_have_basin_lowland_context"]
    assert summary["acceptance"]["major_components_have_orogen_context"]
    assert summary["acceptance"]["dominant_reference_class_capped"]
    assert summary["acceptance"]["production_province_graph_available"]
    assert not summary["acceptance"]["production_province_graph_still_pending"]
    assert summary["acceptance"]["volcanic_lip_plateau_covered"]
    assert summary["acceptance"]["p81_boundary_process_geometry_reference_pending"]

    comparison = summary["generated_province_reference_comparison"]
    assert comparison["schema"] == "aevum.generated_province_reference_comparison.v1"
    assert comparison["status"] == "generated_province_reference_comparison_ready"
    assert comparison["diversity_status"] == "generated_world_gate_ready"
    assert comparison["metrics"]["cells"] == 900
    assert comparison["metrics"]["major_component_count"] >= 2
    assert comparison["metrics"]["generated_reference_class_count"] >= 9
    assert comparison["metrics"]["reference_class_count"] >= 9
    assert comparison["metrics"]["missing_required_feature_class_count"] == 0
    assert comparison["metrics"]["missing_required_parent_process_count"] == 0
    assert comparison["metrics"]["unexpected_missing_reference_class_count"] == 0
    assert comparison["metrics"]["unexpected_missing_required_class_edge_count"] == 0
    assert comparison["metrics"]["missing_reference_class_count"] == 0
    assert comparison["metrics"]["missing_required_class_edge_count"] == 0
    assert comparison["metrics"]["production_province_class_count"] >= 9
    assert comparison["metrics"]["production_province_object_count"] > 0
    assert comparison["metrics"]["min_major_component_reference_class_count"] >= 6
    assert comparison["metrics"]["min_major_component_reference_class_count_gt_min_share"] >= 5
    assert comparison["metrics"]["max_major_component_largest_reference_class_fraction"] <= 0.65
    assert comparison["metrics"]["failing_major_component_count"] == 0
    assert comparison["metrics"]["landform_object_count"] > 0
    assert comparison["metrics"]["landform_overlay_cell_count"] > 0
    assert comparison["missing_required_feature_classes"] == ()
    assert comparison["missing_required_parent_processes"] == ()
    assert comparison["missing_reference_classes"] == ()
    assert comparison["missing_required_class_edges"] == ()
    assert "active_orogen|foreland_basin" in comparison["generated_class_edge_counts"]
    assert "platform|shield" in comparison["generated_class_edge_counts"]
    assert all(
        value for key, value in comparison["acceptance"].items()
        if key != "production_province_graph_still_pending")
    assert comparison["acceptance"]["production_province_graph_available"]
    assert not comparison["acceptance"]["production_province_graph_still_pending"]
    assert not comparison["limitations"]["first_class_production_province_ids_missing"]
    assert not comparison["limitations"]["province_class_field_derived_from_detail_and_landform_objects"]
    assert not comparison["limitations"]["volcanic_lip_plateau_expression_missing"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P80.generated_province_graph_reference_comparison"}
    metrics = benchmarks["P80.generated_province_graph_reference_comparison"]["metrics"]
    assert metrics["comparison_status"] == "generated_province_reference_comparison_ready"
    assert metrics["missing_reference_classes"] == ()
    assert metrics["missing_required_class_edges"] == ()
    assert metrics["production_province_graph_available"]
    assert not metrics["production_province_graph_still_pending"]
    assert metrics["expected_residuals_recorded"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p80_microbenchmarks.csv").exists()


def test_p81_boundary_process_geometry_reference_pass(tmp_path):
    summary = run_suite("P81", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p81.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["boundary_process_geometry_reference"]
    assert summary["acceptance"]["reference_geometry_ready"]
    assert summary["acceptance"]["generated_world_comparison_available"]
    assert summary["acceptance"]["current_boundary_network_present"]
    assert summary["acceptance"]["current_key_types_present"]
    assert summary["acceptance"]["current_expected_residuals_recorded"]
    assert summary["acceptance"]["current_unexpected_missing_types_absent"]
    assert summary["acceptance"]["transform_residual_recorded"]
    assert summary["acceptance"]["raw_vectors_not_stored"]
    assert summary["acceptance"]["direct_vector_extraction_marked_pending"]
    assert summary["acceptance"]["p82_wilson_cycle_lifecycle_reference_pending"]

    geometry = summary["boundary_process_geometry"]
    assert geometry["schema"] == "aevum.boundary_process_geometry.v1"
    assert geometry["status"] == "boundary_process_geometry_ready"
    assert geometry["reference"]["status"] == "boundary_process_geometry_reference_ready"
    assert geometry["current_generated"]["status"] == "generated_boundary_process_geometry_ready"

    reference = geometry["reference"]
    synthetic = reference["synthetic_network"]
    assert reference["missing_source_ids"] == ()
    assert reference["acceptance"]["source_ids_known"]
    assert reference["acceptance"]["raw_vectors_not_stored"]
    assert reference["acceptance"]["direct_vector_extraction_marked_pending"]
    assert synthetic["process_type_count"] == 7
    assert set(synthetic["present_process_types"]) == set(reference["process_types"])
    assert all(check["in_envelope"] for check in synthetic["length_fraction_checks"].values())
    assert synthetic["transform_offset_count"] >= 2
    assert synthetic["antimeridian_ridge_component_count"] == 1
    assert synthetic["adjacency"]["transform_near_ridge_fraction"] >= 0.65
    assert synthetic["adjacency"]["trench_near_active_margin_fraction"] >= 0.75
    assert synthetic["adjacency"]["collision_near_diffuse_fraction"] >= 0.60
    assert synthetic["mean_component_sinuosity"] >= 1.05
    assert synthetic["trench_longitude_std_deg"] >= 2.0

    current = geometry["current_generated"]
    assert current["missing_process_types"] == ()
    assert current["unexpected_missing_process_types"] == ()
    assert not current["limitations"]["current_transform_boundary_missing"]
    assert current["acceptance"]["expected_residuals_recorded"]
    assert current["metrics"]["process_type_count"] >= 6
    assert current["metrics"]["any_boundary_cell_count"] > 0
    assert current["metrics"]["type_cell_counts"]["ridge"] > 0
    assert current["metrics"]["type_cell_counts"]["transform"] > 0
    assert current["metrics"]["type_cell_counts"]["subduction_trench"] > 0
    assert current["metrics"]["type_cell_counts"]["collision_suture"] > 0
    assert current["metrics"]["type_cell_counts"]["diffuse_deformation"] > 0
    assert current["metrics"]["type_cell_counts"]["passive_margin"] > 0
    assert current["metrics"]["type_cell_counts"]["continental_rift"] > 0
    assert current["metrics"]["adjacency"]["transform_near_ridge_fraction"] >= 0.50

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P81.boundary_process_geometry_reference"}
    metrics = benchmarks["P81.boundary_process_geometry_reference"]["metrics"]
    assert metrics["reference_status"] == "boundary_process_geometry_reference_ready"
    assert metrics["current_status"] == "generated_boundary_process_geometry_ready"
    assert metrics["missing_source_id_count"] == 0
    assert metrics["synthetic_process_type_count"] == 7
    assert metrics["synthetic_length_check_pass_count"] == metrics["synthetic_length_check_count"] == 7
    assert metrics["current_missing_process_types"] == ()
    assert not metrics["current_unexpected_missing_process_types"]
    assert not metrics["current_transform_boundary_missing"]
    assert metrics["current_transform_boundary_available"]
    assert metrics["current_transform_cells"] > 0
    assert metrics["current_transform_near_ridge_fraction"] >= 0.50
    assert metrics["expected_transform_residual_recorded"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p81_microbenchmarks.csv").exists()


def test_p82_wilson_cycle_lifecycle_reference_pass(tmp_path):
    summary = run_suite("P82", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p82.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["wilson_cycle_lifecycle_reference"]
    assert summary["acceptance"]["scripted_reference_ready"]
    assert summary["acceptance"]["generated_world_audit_available"]
    assert summary["acceptance"]["current_wilson_objects_present"]
    assert summary["acceptance"]["current_gateway_causality_present"]
    assert summary["acceptance"]["current_expected_residuals_recorded"]
    assert summary["acceptance"]["current_unexpected_missing_object_sets_absent"]
    assert summary["acceptance"]["spreading_center_residual_recorded"]
    assert summary["acceptance"]["gplates_replay_marked_pending"]
    assert summary["acceptance"]["p83_crust_sediment_province_coupling_pending"]

    lifecycle = summary["wilson_cycle_lifecycle"]
    assert lifecycle["schema"] == "aevum.wilson_cycle_lifecycle_reference.v1"
    assert lifecycle["status"] == "wilson_cycle_lifecycle_reference_ready"
    assert lifecycle["scripted_reference"]["status"] == "scripted_wilson_cycle_reference_ready"
    assert lifecycle["current_generated"]["status"] == "generated_world_wilson_lifecycle_ready"

    scripted = lifecycle["scripted_reference"]
    assert scripted["frame_count"] == 7
    assert scripted["stage_sequence"] == scripted["expected_stage_sequence"]
    assert scripted["stage_sequence"] == (
        "continental_rift",
        "spreading_ocean",
        "mature_ocean",
        "closing_ocean",
        "closing_arc_margin",
        "suture_relict",
        "old_orogen_relict",
    )
    assert len(scripted["unique_basin_ids"]) == 1
    assert len(scripted["unique_lineage_keys"]) == 1
    assert scripted["phase_codes"] == (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0)
    assert scripted["basin_ages_myr"][-1] >= 500.0
    assert scripted["old_orogen_relict_count"] >= 1
    assert not scripted["parent_link_failures"]
    assert all(scripted["acceptance"].values())
    for key, count in scripted["object_set_counts"].items():
        assert count > 0, key
    assert {"opening", "open", "closing", "restricted", "closed"}.issubset(
        set(scripted["gateway_statuses"]))

    old_orogen_frame = scripted["frames"][-1]
    assert old_orogen_frame["basin_stage"] == "old_orogen_relict"
    assert old_orogen_frame["old_orogen_relicts"]
    old_orogen = old_orogen_frame["old_orogen_relicts"][0]
    assert old_orogen["parent_suture_id"]
    assert old_orogen["parent_basin_id"]
    assert old_orogen["parent_wilson_cycle_id"]

    current = lifecycle["current_generated"]
    assert current["missing_object_sets"] == ()
    assert current["unexpected_missing_object_sets"] == ()
    assert not current["limitations"]["current_spreading_center_objects_missing"]
    assert current["object_counts"]["tectonics.ocean_basins"] > 0
    assert current["object_counts"]["tectonics.wilson_cycles"] > 0
    assert current["object_counts"]["tectonics.ocean_gateways"] > 0
    assert current["object_counts"]["tectonics.rift_systems"] > 0
    assert current["object_counts"]["tectonics.passive_margins"] > 0
    assert current["object_counts"]["tectonics.spreading_centers"] > 0
    assert current["object_counts"]["tectonics.closing_margins"] > 0
    assert current["object_counts"]["tectonics.sutures"] > 0
    assert len(current["basin_stage_counts"]) >= 3
    assert len(current["active_phase_codes"]) >= 3
    assert not current["parent_link_failures"]
    assert all(current["acceptance"].values())

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P82.wilson_cycle_lifecycle_reference"}
    metrics = benchmarks["P82.wilson_cycle_lifecycle_reference"]["metrics"]
    assert metrics["reference_status"] == "scripted_wilson_cycle_reference_ready"
    assert metrics["current_status"] == "generated_world_wilson_lifecycle_ready"
    assert metrics["scripted_frame_count"] == 7
    assert metrics["scripted_unique_basin_id_count"] == 1
    assert metrics["scripted_unique_lineage_key_count"] == 1
    assert metrics["scripted_required_object_sets_observed_count"] == metrics["scripted_required_object_set_count"]
    assert metrics["scripted_parent_link_failure_count"] == 0
    assert metrics["scripted_old_orogen_relict_count"] >= 1
    assert metrics["current_missing_object_sets"] == ()
    assert not metrics["current_unexpected_missing_object_sets"]
    assert metrics["current_parent_link_failure_count"] == 0
    assert not metrics["current_spreading_center_residual_recorded"]
    assert metrics["current_spreading_center_objects_available"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p82_microbenchmarks.csv").exists()


def test_p83_crust_sediment_province_coupling_pass(tmp_path):
    summary = run_suite("P83", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p83.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["crust_sediment_province_coupling"]
    assert summary["acceptance"]["reference_fixture_ready"]
    assert summary["acceptance"]["current_generated_audit_available"]
    assert summary["acceptance"]["required_classes_present"]
    assert summary["acceptance"]["required_parent_processes_present"]
    assert summary["acceptance"]["crust_thickness_ordering"]
    assert summary["acceptance"]["sediment_accommodation_ordering"]
    assert summary["acceptance"]["shield_old_stable_not_high_flat"]
    assert summary["acceptance"]["basins_and_passive_margins_low_not_erased"]
    assert summary["acceptance"]["current_core_fields_available"]
    assert summary["acceptance"]["current_object_signal_present"]
    assert summary["acceptance"]["current_expected_residuals_recorded"]
    assert summary["acceptance"]["production_province_graph_available"]
    assert summary["acceptance"]["p84_source_to_sink_sediment_budget_pending"]

    coupling = summary["crust_sediment_province_coupling"]
    assert coupling["schema"] == "aevum.crust_sediment_province_coupling.v1"
    assert coupling["status"] == "crust_sediment_province_coupling_ready"

    reference = coupling["reference"]
    assert reference["status"] == "crust_sediment_province_reference_ready"
    assert reference["source_ids"] == (
        "CRUST1_0",
        "GLIM_GLOBAL_LITHOLOGY",
        "NOAA_TOTAL_SEDIMENT_THICKNESS",
    )
    assert not reference["missing_required_classes"]
    assert not reference["missing_required_parent_processes"]
    assert not reference["unparented_province_ids"]
    assert not reference["random_sourced_province_ids"]
    assert all(reference["acceptance"].values())
    assert set(reference["classes_present"]).issuperset({
        "shield",
        "platform",
        "intracratonic_basin",
        "foreland_basin",
        "active_orogen",
        "old_orogen",
        "rift_basin",
        "passive_margin_lowland",
        "continental_shelf",
        "volcanic_lip_plateau",
    })
    class_means = reference["class_means"]
    assert class_means["elevation_m"]["active_orogen"] > class_means["elevation_m"]["platform"]
    assert class_means["elevation_m"]["intracratonic_basin"] < class_means["elevation_m"]["platform"]
    assert class_means["sediment_m"]["foreland_basin"] > class_means["sediment_m"]["platform"]
    assert class_means["sediment_m"]["continental_shelf"] > class_means["sediment_m"]["passive_margin_lowland"]
    assert class_means["crust_thickness_m"]["active_orogen"] > class_means["crust_thickness_m"]["platform"]
    assert class_means["relief_m"]["shield"] > class_means["relief_m"]["platform"]
    assert reference["extraction_policy"]["raw_crust_lithology_sediment_data_stored"] is False
    assert reference["extraction_policy"]["direct_crust_lithology_sediment_extraction_pending"]

    current = coupling["current_generated"]
    assert current["status"] == "generated_world_crust_sediment_province_audit_ready"
    assert current["missing_fields"] == ()
    assert current["missing_production_province_fields"] == ()
    assert not current["limitations"]["production_continental_province_ids_missing"]
    assert current["kind_groups_present"]["core"]
    assert current["kind_groups_present"]["platform"]
    assert current["kind_groups_present"]["basin_lowland"]
    assert current["kind_groups_present"]["orogen"]
    assert current["metrics"]["object_count"] > 0
    assert current["metrics"]["kind_count"] >= 6
    assert current["metrics"]["basin_lowland_mean_elevation_m"] < current["metrics"]["platform_mean_elevation_m"]
    assert current["metrics"]["basin_lowland_mean_sediment_m"] > current["metrics"]["platform_mean_sediment_m"]
    assert current["metrics"]["orogen_mean_elevation_m"] > current["metrics"]["basin_lowland_mean_elevation_m"]
    assert all(current["acceptance"].values())

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P83.crust_sediment_province_coupling"}
    metrics = benchmarks["P83.crust_sediment_province_coupling"]["metrics"]
    assert metrics["reference_status"] == "crust_sediment_province_reference_ready"
    assert metrics["current_status"] == "generated_world_crust_sediment_province_audit_ready"
    assert metrics["province_class_count"] >= 12
    assert metrics["parent_process_count"] >= 12
    assert metrics["missing_required_class_count"] == 0
    assert metrics["missing_required_parent_process_count"] == 0
    assert metrics["unparented_province_id_count"] == 0
    assert metrics["random_sourced_province_id_count"] == 0
    assert metrics["current_missing_field_count"] == 0
    assert metrics["current_missing_production_province_field_count"] == 0
    assert metrics["current_missing_production_province_fields"] == ()
    assert metrics["p84_source_to_sink_sediment_budget_pending"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p83_microbenchmarks.csv").exists()


def test_p84_source_to_sink_sediment_budget_pass(tmp_path):
    summary = run_suite("P84", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p84.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["source_to_sink_sediment_budget"]
    assert summary["acceptance"]["reference_budget_ready"]
    assert summary["acceptance"]["sediment_volume_conserved"]
    assert summary["acceptance"]["no_land_mask_regression"]
    assert summary["acceptance"]["deposition_within_accommodation"]
    assert summary["acceptance"]["routing_edges_close_budget"]
    assert summary["acceptance"]["current_generated_audit_available"]
    assert summary["acceptance"]["current_expected_residuals_recorded"]
    assert summary["acceptance"]["production_source_to_sink_objects_available"]
    assert summary["acceptance"]["production_sediment_budget_closes"]
    assert summary["acceptance"]["production_source_sink_kinds_diverse"]
    assert not summary["acceptance"]["production_source_to_sink_objects_still_pending"]
    assert summary["acceptance"]["p85_drainage_divide_province_alignment_pending"]

    budget = summary["source_to_sink_sediment_budget"]
    assert budget["schema"] == "aevum.source_to_sink_sediment_budget.v1"
    assert budget["status"] == "source_to_sink_sediment_budget_ready"

    reference = budget["reference"]
    assert reference["status"] == "source_to_sink_sediment_budget_reference_ready"
    assert reference["source_ids"] == (
        "NOAA_TOTAL_SEDIMENT_THICKNESS",
        "GLIM_GLOBAL_LITHOLOGY",
    )
    assert not reference["missing_required_zone_kinds"]
    assert not reference["invalid_edges"]
    assert not reference["routing_mismatches"]
    assert not reference["land_mask_changes"]
    assert not reference["accommodation_violations"]
    assert not reference["erosion_violations"]
    assert all(reference["acceptance"].values())
    assert reference["extraction_policy"]["raw_sediment_or_drainage_data_stored"] is False
    assert reference["extraction_policy"]["direct_sediment_drainage_extraction_pending"]

    ref_metrics = reference["metrics"]
    assert ref_metrics["zone_count"] == 6
    assert ref_metrics["edge_count"] == 5
    assert ref_metrics["source_zone_count"] == 2
    assert ref_metrics["sink_zone_count"] == 4
    assert ref_metrics["source_volume_km3"] == ref_metrics["sink_volume_km3"] == 69000.0
    assert ref_metrics["volume_balance_fraction"] <= 1.0e-9
    assert ref_metrics["max_accommodation_utilization"] < 1.0
    assert ref_metrics["land_mask_change_count"] == 0
    assert ref_metrics["routing_mismatch_count"] == 0
    assert ref_metrics["mountain_source_export_km3"] > ref_metrics["platform_source_export_km3"]
    assert ref_metrics["foreland_sink_deposition_km3"] > ref_metrics["passive_margin_sink_deposition_km3"]
    assert ref_metrics["shelf_sink_deposition_km3"] > ref_metrics["ocean_basin_sink_deposition_km3"]

    current = budget["current_generated"]
    assert current["status"] == "generated_world_source_to_sink_audit_ready"
    assert current["missing_source_to_sink_objects"] == ()
    assert not current["limitations"]["production_source_to_sink_objects_missing"]
    assert current["metrics"]["drainage_basin_object_count"] >= 8
    assert current["metrics"]["routing_edge_count"] >= 12
    assert current["metrics"]["sediment_budget_object_count"] == 1
    assert current["metrics"]["source_volume_km3"] > 0.0
    assert current["metrics"]["sink_volume_km3"] > 0.0
    assert current["metrics"]["production_volume_balance_fraction"] <= 1.0e-9
    assert current["metrics"]["routing_source_kind_count"] >= 2
    assert current["metrics"]["routing_sink_kind_count"] >= 2
    assert current["metrics"]["sediment_field_available"]
    assert current["metrics"]["elevation_field_available"]
    assert current["metrics"]["sediment_max_m"] > current["metrics"]["sediment_min_m"]
    assert all(current["acceptance"].values())

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P84.source_to_sink_sediment_budget"}
    metrics = benchmarks["P84.source_to_sink_sediment_budget"]["metrics"]
    assert metrics["reference_status"] == "source_to_sink_sediment_budget_reference_ready"
    assert metrics["current_status"] == "generated_world_source_to_sink_audit_ready"
    assert metrics["zone_count"] == 6
    assert metrics["edge_count"] == 5
    assert metrics["volume_balance_fraction"] <= 1.0e-9
    assert metrics["land_mask_change_count"] == 0
    assert metrics["routing_mismatch_count"] == 0
    assert metrics["current_missing_object_count"] == 0
    assert metrics["current_drainage_basin_object_count"] >= 8
    assert metrics["current_routing_edge_count"] >= 12
    assert metrics["current_sediment_budget_object_count"] == 1
    assert metrics["current_production_volume_balance_fraction"] <= 1.0e-9
    assert metrics["current_sediment_field_available"]
    assert metrics["current_elevation_field_available"]
    assert metrics["p85_drainage_divide_province_alignment_pending"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p84_microbenchmarks.csv").exists()


def test_p85_drainage_divide_province_alignment_pass(tmp_path):
    summary = run_suite("P85", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p85.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["drainage_divide_province_alignment"]
    assert summary["acceptance"]["reference_divide_fixture_ready"]
    assert summary["acceptance"]["divide_boundary_alignment"]
    assert summary["acceptance"]["flow_to_expected_sinks"]
    assert summary["acceptance"]["basins_not_checkerboarded"]
    assert summary["acceptance"]["current_generated_audit_available"]
    assert summary["acceptance"]["current_expected_residuals_recorded"]
    assert summary["acceptance"]["production_drainage_fields_available"]
    assert summary["acceptance"]["production_drainage_objects_available"]
    assert summary["acceptance"]["production_basins_contiguous"]
    assert summary["acceptance"]["production_divides_aligned"]
    assert summary["acceptance"]["production_flow_paths_stay_in_basins"]
    assert not summary["acceptance"]["production_drainage_objects_still_pending"]
    assert summary["acceptance"]["p86_old_orogen_erosion_decay_pending"]

    alignment = summary["drainage_divide_province_alignment"]
    assert alignment["schema"] == "aevum.drainage_divide_province_alignment.v1"
    assert alignment["status"] == "drainage_divide_province_alignment_ready"

    reference = alignment["reference"]
    assert reference["status"] == "drainage_divide_reference_ready"
    assert reference["source_ids"] == ("HYDROSHEDS_HYDROBASINS_HYDRORIVERS",)
    assert not reference["missing_required_province_classes"]
    assert all(reference["acceptance"].values())
    assert reference["extraction_policy"]["raw_hydrology_data_stored"] is False
    assert reference["extraction_policy"]["direct_hydrology_extraction_pending"]

    ref_metrics = reference["metrics"]
    assert ref_metrics["province_class_count"] >= 10
    assert ref_metrics["basin_count"] == 3
    assert ref_metrics["divide_cell_count"] == 5
    assert 0.05 <= ref_metrics["divide_fraction"] <= 0.20
    assert ref_metrics["divide_alignment_fraction"] >= 0.90
    assert ref_metrics["highland_alignment_fraction"] >= 0.75
    assert ref_metrics["flow_path_count"] == 6
    assert ref_metrics["flow_to_sink_consistency_fraction"] >= 0.95
    assert ref_metrics["downhill_step_fraction"] == 1.0
    assert ref_metrics["uphill_step_count"] == 0
    assert ref_metrics["divide_crossing_count"] == 0
    assert ref_metrics["basin_crossing_count"] == 0
    assert ref_metrics["sink_failure_count"] == 0
    assert ref_metrics["max_basin_component_count"] == 1
    assert set(reference["basin_component_counts"].values()) == {1}

    current = alignment["current_generated"]
    assert current["status"] == "generated_world_drainage_context_audit_ready"
    assert current["missing_drainage_items"] == ()
    assert not current["limitations"]["production_drainage_objects_missing"]
    assert current["metrics"]["drainage_basin_field_count"] >= 8
    assert current["metrics"]["drainage_basin_object_count"] >= 8
    assert current["metrics"]["drainage_divide_object_count"] > 0
    assert current["metrics"]["divide_cell_count"] > 0
    assert current["metrics"]["divide_fraction_of_land"] <= 0.55
    assert current["metrics"]["current_divide_alignment_fraction"] >= 0.70
    assert current["metrics"]["major_basin_component_failure_count"] == 0
    assert current["metrics"]["flow_path_crossing_count"] == 0
    assert current["metrics"]["flow_path_downhill_violation_count"] == 0
    assert current["metrics"]["landform_object_count"] > 0
    assert (
        current["metrics"]["divide_parent_context_count"] > 0
        or current["metrics"]["production_divide_parent_context_count"] > 0
    )
    assert (
        current["metrics"]["basin_sink_context_count"] > 0
        or current["metrics"]["production_basin_sink_context_count"] > 0
    )
    assert (
        current["metrics"]["margin_sink_context_count"] > 0
        or current["metrics"]["production_margin_sink_context_count"] > 0
    )
    assert all(current["acceptance"].values())

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P85.drainage_divide_province_alignment"}
    metrics = benchmarks["P85.drainage_divide_province_alignment"]["metrics"]
    assert metrics["reference_status"] == "drainage_divide_reference_ready"
    assert metrics["current_status"] == "generated_world_drainage_context_audit_ready"
    assert metrics["divide_alignment_fraction"] >= 0.90
    assert metrics["flow_to_sink_consistency_fraction"] >= 0.95
    assert metrics["max_basin_component_count"] == 1
    assert metrics["current_missing_item_count"] == 0
    assert metrics["current_drainage_basin_field_count"] >= 8
    assert metrics["current_drainage_basin_object_count"] >= 8
    assert metrics["current_divide_fraction_of_land"] <= 0.55
    assert metrics["current_major_basin_component_failure_count"] == 0
    assert metrics["current_flow_to_sink_consistency_fraction"] >= 0.98
    assert metrics["p86_old_orogen_erosion_decay_pending"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p85_microbenchmarks.csv").exists()


def test_p86_old_orogen_erosion_decay_pass(tmp_path):
    summary = run_suite("P86", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p86.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["old_orogen_erosion_decay"]
    assert summary["acceptance"]["reference_decay_fixture_ready"]
    assert summary["acceptance"]["relief_decay_large_enough"]
    assert summary["acceptance"]["sediment_export_declines_late"]
    assert summary["acceptance"]["boundary_persists_after_decay"]
    assert summary["acceptance"]["current_generated_audit_available"]
    assert summary["acceptance"]["current_expected_residuals_recorded"]
    assert summary["acceptance"]["production_old_orogen_decay_fields_available"]
    assert summary["acceptance"]["production_old_orogen_decay_budget_available"]
    assert not summary["acceptance"]["production_old_orogen_decay_objects_still_pending"]
    assert summary["acceptance"]["p87_mountain_inventory_expression_pending"]

    decay = summary["old_orogen_erosion_decay"]
    assert decay["schema"] == "aevum.old_orogen_erosion_decay.v1"
    assert decay["status"] == "old_orogen_erosion_decay_ready"

    reference = decay["reference"]
    assert reference["status"] == "old_orogen_erosion_decay_reference_ready"
    assert reference["source_ids"] == ("NOAA_ETOPO_2022", "GMBA_MOUNTAIN_INVENTORY")
    assert reference["stage_sequence"] == reference["expected_stage_sequence"]
    assert all(reference["acceptance"].values())
    assert reference["extraction_policy"]["raw_topography_or_mountain_inventory_stored"] is False
    assert reference["extraction_policy"]["direct_orogen_decay_extraction_pending"]

    ref_metrics = reference["metrics"]
    assert ref_metrics["frame_count"] == 5
    assert ref_metrics["initial_relief_m"] >= 1800.0
    assert 250.0 <= ref_metrics["final_relief_m"] <= 700.0
    assert ref_metrics["relief_decay_fraction"] >= 0.70
    assert ref_metrics["final_boundary_strength"] >= 0.55
    assert ref_metrics["min_boundary_trace_overlap"] >= 0.90
    assert ref_metrics["total_sediment_export_km3"] > 50000.0
    assert (
        ref_metrics["final_interval_sediment_export_km3"]
        < ref_metrics["peak_sediment_export_km3"] * 0.55
    )
    assert ref_metrics["parent_link_failure_count"] == 0

    current = decay["current_generated"]
    assert current["status"] == "generated_world_old_orogen_decay_audit_ready"
    assert current["missing_decay_fields"] == ()
    assert not current["limitations"]["production_old_orogen_decay_budget_missing"]
    assert current["metrics"]["old_subdued_orogen_object_count"] > 0
    assert current["metrics"]["parented_old_subdued_orogen_object_count"] > 0
    assert current["metrics"]["current_old_orogen_decay_stage_count"] > 0
    assert current["metrics"]["current_mean_orogen_boundary_memory"] >= 0.50
    assert current["metrics"]["current_mean_orogen_erosion_budget_m"] > 0.0
    assert current["metrics"]["current_orogen_sediment_export_volume_km3"] > 0.0
    assert current["metrics"]["mean_old_orogen_elevation_m"] > 0.0
    assert current["metrics"]["old_orogen_area_m2"] > 0.0
    assert all(current["acceptance"].values())

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P86.old_orogen_erosion_decay"}
    metrics = benchmarks["P86.old_orogen_erosion_decay"]["metrics"]
    assert metrics["reference_status"] == "old_orogen_erosion_decay_reference_ready"
    assert metrics["current_status"] == "generated_world_old_orogen_decay_audit_ready"
    assert metrics["stage_sequence"] == reference["expected_stage_sequence"]
    assert metrics["relief_decay_fraction"] >= 0.70
    assert metrics["final_boundary_strength"] >= 0.55
    assert metrics["min_boundary_trace_overlap"] >= 0.90
    assert metrics["current_old_subdued_orogen_object_count"] > 0
    assert metrics["current_parented_old_subdued_orogen_object_count"] > 0
    assert metrics["current_missing_decay_field_count"] == 0
    assert metrics["current_old_orogen_decay_stage_count"] > 0
    assert metrics["current_mean_orogen_boundary_memory"] >= 0.50
    assert metrics["current_orogen_sediment_export_volume_km3"] > 0.0
    assert metrics["p87_mountain_inventory_expression_pending"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p86_microbenchmarks.csv").exists()


def test_p87_mountain_inventory_expression_pass(tmp_path):
    summary = run_suite("P87", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p87.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["mountain_inventory_expression"]
    assert summary["acceptance"]["reference_inventory_ready"]
    assert summary["acceptance"]["reference_hierarchy_and_parentage_ready"]
    assert summary["acceptance"]["current_generated_audit_available"]
    assert summary["acceptance"]["current_object_backed_mountains_present"]
    assert summary["acceptance"]["current_parentage_context_present"]
    assert summary["acceptance"]["current_area_and_shape_metrics_available"]
    assert summary["acceptance"]["current_expected_residuals_recorded"]
    assert summary["acceptance"]["production_mountain_inventory_fields_available"]
    assert summary["acceptance"]["production_expected_mountain_kinds_available"]
    assert not summary["acceptance"]["production_mountain_inventory_still_pending"]
    assert not summary["acceptance"]["production_active_orogen_plateau_expression_still_pending"]
    assert not summary["acceptance"]["production_range_elongation_still_underdeveloped"]
    assert summary["acceptance"]["p88_rift_margin_escarpment_sequence_pending"]

    inventory = summary["mountain_inventory_expression"]
    assert inventory["schema"] == "aevum.mountain_inventory_expression.v1"
    assert inventory["status"] == "mountain_inventory_expression_ready"

    reference = inventory["reference"]
    assert reference["status"] == "mountain_inventory_reference_ready"
    assert reference["source_ids"] == (
        "GMBA_MOUNTAIN_INVENTORY",
        "NOAA_ETOPO_2022",
        "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
    )
    assert all(reference["acceptance"].values())
    assert reference["extraction_policy"]["raw_mountain_inventory_stored"] is False
    assert reference["extraction_policy"]["direct_gmba_extraction_pending"]

    ref_metrics = reference["metrics"]
    assert ref_metrics["range_count"] >= 10
    assert ref_metrics["mountain_class_count"] >= 6
    assert ref_metrics["hierarchy_level_count"] >= 3
    assert ref_metrics["hierarchy_parent_link_failure_count"] == 0
    assert ref_metrics["parent_process_failure_count"] == 0
    assert ref_metrics["threshold_only_range_count"] == 0
    assert ref_metrics["object_backing_failure_count"] == 0
    assert 0.04 <= ref_metrics["total_mountain_area_fraction_world"] <= 0.10
    assert ref_metrics["max_range_area_fraction_world"] <= 0.020
    assert ref_metrics["median_elongation_ratio"] >= 3.0
    assert ref_metrics["elongated_range_count"] >= 6
    assert ref_metrics["max_relief_p90_p10_m"] >= 4000.0

    current = inventory["current_generated"]
    assert current["status"] == "generated_world_mountain_inventory_audit_ready"
    assert current["missing_expected_mountain_kinds"] == ()
    assert current["missing_inventory_fields"] == ()
    assert not current["limitations"]["first_class_mountain_inventory_missing"]
    assert not current["limitations"]["active_orogen_or_plateau_expression_missing"]
    assert not current["limitations"]["elongated_range_expression_underdeveloped"]
    assert current["metrics"]["production_mountain_range_object_count"] > 0
    assert current["metrics"]["mountain_candidate_object_count"] > 0
    assert current["metrics"]["expressed_mountain_object_count"] > 0
    assert current["metrics"]["parent_process_coverage_fraction"] >= 0.95
    assert current["metrics"]["parent_object_context_coverage_fraction"] >= 0.95
    assert current["metrics"]["total_mountain_area_fraction_world"] <= 0.20
    assert all(current["acceptance"].values())

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P87.mountain_inventory_expression"}
    metrics = benchmarks["P87.mountain_inventory_expression"]["metrics"]
    assert metrics["reference_status"] == "mountain_inventory_reference_ready"
    assert metrics["current_status"] == "generated_world_mountain_inventory_audit_ready"
    assert metrics["reference_range_count"] >= 10
    assert metrics["current_mountain_candidate_object_count"] > 0
    assert metrics["current_missing_expected_kind_count"] == 0
    assert metrics["current_missing_inventory_field_count"] == 0
    assert metrics["production_mountain_inventory_fields_available"]
    assert metrics["production_expected_mountain_kinds_available"]
    assert not metrics["current_active_orogen_or_plateau_residual_recorded"]
    assert not metrics["current_elongation_residual_recorded"]
    assert metrics["p88_rift_margin_escarpment_sequence_pending"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p87_microbenchmarks.csv").exists()


def test_p88_rift_margin_escarpment_sequence_pass(tmp_path):
    summary = run_suite("P88", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p88.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["rift_margin_escarpment_sequence"]
    assert summary["acceptance"]["reference_sequence_ready"]
    assert summary["acceptance"]["reference_ordering_ready"]
    assert summary["acceptance"]["current_generated_audit_available"]
    assert summary["acceptance"]["current_lowland_shelf_coupled"]
    assert summary["acceptance"]["current_rift_passive_context_present"]
    assert summary["acceptance"]["current_expected_residuals_recorded"]
    assert summary["acceptance"]["production_rift_margin_sequence_fields_available"]
    assert summary["acceptance"]["production_rift_margin_sequence_ready"]
    assert summary["acceptance"]["production_shelf_slope_rise_abyss_ordered"]
    assert not summary["acceptance"]["production_rift_margin_sequence_still_pending"]
    assert not summary["acceptance"]["production_rift_shoulder_escarpment_still_pending"]
    assert summary["acceptance"]["p89_plateau_area_cap_and_decay_pending"]

    sequence = summary["rift_margin_escarpment_sequence"]
    assert sequence["schema"] == "aevum.rift_margin_escarpment_sequence.v1"
    assert sequence["status"] == "rift_margin_escarpment_sequence_ready"

    reference = sequence["reference"]
    assert reference["status"] == "rift_margin_escarpment_reference_ready"
    assert reference["source_ids"] == (
        "NOAA_ETOPO_2022",
        "GEBCO_GRIDDED_BATHYMETRY",
        "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
        "NOAA_TOTAL_SEDIMENT_THICKNESS",
    )
    assert reference["stage_sequence"] == reference["expected_stage_sequence"]
    assert all(reference["acceptance"].values())
    assert reference["extraction_policy"]["raw_topography_bathymetry_or_geologic_vectors_stored"] is False
    assert reference["extraction_policy"]["direct_rift_margin_extraction_pending"]

    ref_metrics = reference["metrics"]
    assert ref_metrics["zone_count"] == 11
    assert ref_metrics["class_count"] >= 10
    assert ref_metrics["edge_count"] >= 10
    assert ref_metrics["missing_required_edge_count"] == 0
    assert ref_metrics["parent_process_failure_count"] == 0
    assert ref_metrics["rift_shoulder_elevation_m"] > ref_metrics["rift_basin_elevation_m"] + 900.0
    assert ref_metrics["escarpment_relief_m"] >= 600.0
    assert ref_metrics["shelf_depth_m"] < ref_metrics["slope_depth_m"] < ref_metrics["rise_depth_m"] < ref_metrics["abyss_depth_m"]
    assert ref_metrics["shelf_sediment_m"] > ref_metrics["passive_margin_lowland_sediment_m"]
    assert ref_metrics["rift_basin_sediment_m"] > 2000.0

    current = sequence["current_generated"]
    assert current["status"] == "generated_world_rift_margin_sequence_audit_ready"
    assert current["missing_sequence_items"] == ()
    assert not current["limitations"]["first_class_rift_margin_sequence_missing"]
    assert not current["limitations"]["rift_shoulder_objects_missing"]
    assert not current["limitations"]["escarpment_objects_missing"]
    assert not current["limitations"]["passive_margin_lowland_objects_tiny"]
    assert not current["limitations"]["rift_to_margin_lineage_missing"]
    assert current["metrics"]["rift_basin_object_count"] > 0
    assert (
        current["metrics"]["passive_margin_lowland_object_count"] > 0
        or current["metrics"]["production_passive_lowland_area_fraction_world"] > 0.0
    )
    assert current["metrics"]["passive_margin_wedge_object_count"] > 0
    assert current["metrics"]["lowland_near_shelf_fraction_p2"] >= 0.75
    assert current["metrics"]["lowland_near_wedge_fraction_p2"] >= 0.50
    assert current["metrics"]["rift_near_passive_margin_fraction_p5"] >= 0.50
    assert current["metrics"]["shelf_depth_p75_m"] < 800.0
    assert current["metrics"]["abyss_depth_p50_m"] > current["metrics"]["shelf_depth_p75_m"] + 1000.0
    assert current["metrics"]["parented_rift_basin_object_count"] > 0
    assert (
        current["metrics"]["parented_passive_margin_lowland_object_count"] > 0
        or current["metrics"]["production_sequence_object_count"] > 0
    )
    assert current["metrics"]["parented_passive_margin_wedge_object_count"] > 0
    assert current["metrics"]["production_sequence_field_cell_count"] > 0
    assert current["metrics"]["production_sequence_id_count"] > 0
    assert current["metrics"]["production_lineage_id_count"] > 0
    assert current["metrics"]["production_stage_count"] >= 6
    assert current["metrics"]["production_rift_shoulder_cell_count"] > 0
    assert current["metrics"]["production_escarpment_cell_count"] > 0
    assert current["metrics"]["production_sequence_object_count"] > 0
    assert current["metrics"]["production_rift_shoulder_object_count"] > 0
    assert current["metrics"]["production_escarpment_object_count"] > 0
    assert current["metrics"]["production_shelf_cell_count"] > 0
    assert current["metrics"]["production_slope_cell_count"] > 0
    assert current["metrics"]["production_rise_cell_count"] > 0
    assert current["metrics"]["production_abyss_cell_count"] > 0
    assert (
        current["metrics"]["production_shelf_depth_p75_m"]
        < current["metrics"]["production_slope_depth_p50_m"]
        < current["metrics"]["production_rise_depth_p50_m"]
        < current["metrics"]["production_abyss_depth_p50_m"]
    )
    assert all(current["acceptance"].values())

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P88.rift_margin_escarpment_sequence"}
    metrics = benchmarks["P88.rift_margin_escarpment_sequence"]["metrics"]
    assert metrics["reference_status"] == "rift_margin_escarpment_reference_ready"
    assert metrics["current_status"] == "generated_world_rift_margin_sequence_audit_ready"
    assert metrics["current_rift_basin_object_count"] > 0
    assert (
        metrics["current_passive_margin_lowland_object_count"] > 0
        or metrics["current_production_passive_lowland_area_fraction_world"] > 0.0
    )
    assert metrics["current_passive_margin_wedge_object_count"] > 0
    assert metrics["current_missing_sequence_item_count"] == 0
    assert metrics["current_production_sequence_field_cell_count"] > 0
    assert metrics["current_production_sequence_id_count"] > 0
    assert metrics["current_production_lineage_id_count"] > 0
    assert metrics["current_production_stage_count"] >= 6
    assert metrics["current_production_rift_shoulder_cell_count"] > 0
    assert metrics["current_production_escarpment_cell_count"] > 0
    assert metrics["current_production_sequence_object_count"] > 0
    assert metrics["current_production_rift_shoulder_object_count"] > 0
    assert metrics["current_production_escarpment_object_count"] > 0
    assert metrics["current_production_shelf_depth_p75_m"] < metrics["current_production_slope_depth_p50_m"]
    assert metrics["current_production_slope_depth_p50_m"] < metrics["current_production_rise_depth_p50_m"]
    assert metrics["current_production_rise_depth_p50_m"] < metrics["current_production_abyss_depth_p50_m"]
    assert not metrics["current_rift_shoulder_residual_recorded"]
    assert not metrics["current_escarpment_residual_recorded"]
    assert not metrics["current_sequence_lineage_residual_recorded"]
    assert not metrics["current_lowland_tiny_residual_recorded"]
    assert metrics["p89_plateau_area_cap_and_decay_pending"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p88_microbenchmarks.csv").exists()


def test_p89_plateau_area_cap_and_decay_pass(tmp_path):
    summary = run_suite("P89", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p89.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["plateau_area_cap_and_decay"]
    assert summary["acceptance"]["reference_plateau_fixture_ready"]
    assert summary["acceptance"]["reference_area_caps_and_decay_ready"]
    assert summary["acceptance"]["current_generated_audit_available"]
    assert summary["acceptance"]["current_plateau_area_cap_audited"]
    assert summary["acceptance"]["current_expected_residuals_recorded"]
    assert summary["acceptance"]["production_plateau_inventory_fields_available"]
    assert summary["acceptance"]["production_plateau_expression_available"]
    assert not summary["acceptance"]["production_plateau_inventory_still_pending"]
    assert not summary["acceptance"]["production_plateau_expression_still_pending"]
    assert not summary["acceptance"]["production_volcanic_lip_plateau_expression_still_pending"]
    assert summary["acceptance"]["p90_current_world_morphology_gap_inventory_pending"]

    plateau = summary["plateau_area_cap_and_decay"]
    assert plateau["schema"] == "aevum.plateau_area_cap_and_decay.v1"
    assert plateau["status"] == "plateau_area_cap_and_decay_ready"

    reference = plateau["reference"]
    assert reference["status"] == "plateau_area_cap_decay_reference_ready"
    assert reference["source_ids"] == (
        "NOAA_ETOPO_2022",
        "GMBA_MOUNTAIN_INVENTORY",
        "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
        "CRUST1_0",
    )
    assert reference["variants"] == reference["expected_variants"]
    assert all(reference["acceptance"].values())
    assert reference["extraction_policy"]["raw_topography_mountain_inventory_or_crustal_grids_stored"] is False
    assert reference["extraction_policy"]["direct_plateau_extraction_pending"]

    ref_metrics = reference["metrics"]
    assert ref_metrics["variant_count"] == 2
    assert ref_metrics["frame_count"] >= 8
    assert ref_metrics["stage_count"] >= 8
    assert ref_metrics["collision_frame_count"] >= 4
    assert ref_metrics["volcanic_frame_count"] >= 4
    assert ref_metrics["parent_process_failure_count"] == 0
    assert ref_metrics["parent_object_failure_count"] == 0
    assert 0.015 <= ref_metrics["max_collision_plateau_area_fraction_world"] <= 0.035
    assert 0.010 <= ref_metrics["max_volcanic_plateau_area_fraction_world"] <= 0.025
    assert ref_metrics["combined_peak_plateau_area_fraction_world"] <= 0.060
    assert ref_metrics["min_background_to_plateau_area_ratio"] >= 6.0
    assert ref_metrics["collision_peak_elevation_m"] > ref_metrics["collision_final_elevation_m"]
    assert ref_metrics["collision_elevation_decay_m"] >= 1500.0
    assert ref_metrics["collision_area_decay_fraction_world"] >= 0.008
    assert ref_metrics["volcanic_peak_elevation_m"] > ref_metrics["volcanic_final_elevation_m"]
    assert ref_metrics["volcanic_elevation_decay_m"] >= 800.0
    assert ref_metrics["volcanic_area_decay_fraction_world"] >= 0.006
    assert ref_metrics["max_plateau_delta_above_platform_m"] >= 1200.0

    current = plateau["current_generated"]
    assert current["status"] == "generated_world_plateau_area_cap_audit_ready"
    assert current["missing_plateau_items"] == ()
    assert current["missing_expected_plateau_kinds"] == ()
    assert not current["limitations"]["first_class_plateau_inventory_missing"]
    assert not current["limitations"]["plateau_landform_expression_missing"]
    assert not current["limitations"]["volcanic_lip_plateau_expression_missing"]
    assert not current["limitations"]["plateau_decay_fields_missing"]
    assert not current["limitations"]["plateau_parent_lineage_missing"]
    assert current["metrics"]["tectonics_lip_object_count"] > 0
    assert current["metrics"]["plateau_object_count"] > 0
    assert current["metrics"]["production_plateau_inventory_cell_count"] > 0
    assert current["metrics"]["production_volcanic_plateau_cell_count"] > 0
    assert current["metrics"]["plateau_detail_area_fraction_world"] <= 0.060
    assert current["metrics"]["plateau_object_area_fraction_world"] <= 0.060
    assert current["metrics"]["max_plateau_object_area_fraction_world"] <= 0.035
    assert all(current["acceptance"].values())

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P89.plateau_area_cap_and_decay"}
    metrics = benchmarks["P89.plateau_area_cap_and_decay"]["metrics"]
    assert metrics["reference_status"] == "plateau_area_cap_decay_reference_ready"
    assert metrics["current_status"] == "generated_world_plateau_area_cap_audit_ready"
    assert metrics["current_lip_object_count"] > 0
    assert metrics["current_plateau_object_count"] > 0
    assert metrics["current_production_plateau_inventory_cell_count"] > 0
    assert metrics["current_production_volcanic_plateau_cell_count"] > 0
    assert metrics["current_missing_plateau_item_count"] == 0
    assert metrics["current_missing_plateau_kind_count"] == 0
    assert metrics["production_plateau_inventory_fields_available"]
    assert metrics["production_plateau_expression_available"]
    assert not metrics["current_first_class_plateau_inventory_residual_recorded"]
    assert not metrics["current_plateau_expression_residual_recorded"]
    assert not metrics["current_volcanic_lip_plateau_residual_recorded"]
    assert not metrics["current_plateau_decay_residual_recorded"]
    assert not metrics["current_plateau_lineage_residual_recorded"]
    assert metrics["p90_current_world_morphology_gap_inventory_pending"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p89_microbenchmarks.csv").exists()


def test_p90_current_world_morphology_gap_inventory_pass(tmp_path):
    summary = run_suite("P90", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p90.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["current_world_morphology_gap_inventory"]
    assert summary["acceptance"]["p76_p89_statuses_available"]
    assert summary["acceptance"]["all_gaps_have_owner_category_and_future_stage"]
    assert summary["acceptance"]["no_generic_visual_blockers"]
    assert summary["acceptance"]["direct_surface_metrics_available"]
    assert summary["acceptance"]["compiler_metrics_available"]
    assert summary["acceptance"]["compiler_consistency_recorded"]
    assert summary["acceptance"]["asset_review_requirements_defined"]
    assert summary["acceptance"]["p91_integrated_real_earth_morphology_promotion_audit_pending"]

    inventory = summary["current_world_morphology_gap_inventory"]
    assert inventory["schema"] == "aevum.current_world_morphology_gap_inventory.v1"
    assert inventory["status"] == "current_world_morphology_gap_inventory_ready"
    assert set(inventory["owner_layers"]) == {
        "planform",
        "province_graph",
        "boundary_lifecycle",
        "crust_sediment",
        "drainage_erosion",
        "landform_expression",
        "bathymetry_margin",
        "compiler_render",
    }
    assert set(inventory["gap_categories"]).issuperset({
        "missing_process",
        "missing_field",
        "wrong_amplitude",
        "wrong_area_scale",
        "wrong_lifecycle",
        "asset_review_pending",
    })
    assert inventory["missing_review_assets"] == (
        "elevation.png",
        "terrain_provinces.png",
        "continental_detail_provinces.png",
        "ocean_depth_provinces.png",
        "crust_age.png",
        "history.png",
        "timeline.png",
        "hexmap.png",
    )
    assert inventory["unassigned_gaps"] == ()
    assert inventory["generic_blockers"] == ()
    assert inventory["p76_p89_statuses"]["P76.source_ledger"] == "source_ledger_schema_ready"
    assert inventory["p76_p89_statuses"]["P89.plateau"] == "plateau_area_cap_and_decay_ready"

    metrics = inventory["metrics"]
    assert metrics["cells"] == 900
    assert metrics["gap_count"] >= 19
    assert metrics["owner_layer_count"] == 3
    assert metrics["category_count"] >= 3
    assert metrics["source_suite_count"] >= 4
    assert metrics["unassigned_gap_count"] == 0
    assert metrics["generic_blocker_count"] == 0
    assert metrics["current_residual_item_count"] >= 8
    assert metrics["missing_required_asset_count"] == metrics["required_review_asset_count"] == 8
    assert metrics["compiler_passed_envelope"]
    assert metrics["compiler_broad_land_to_water_fraction"] <= 0.08
    assert metrics["compiler_broad_ocean_to_land_fraction"] <= 0.06
    assert metrics["compiler_shelf_as_deep_ocean_fraction"] <= 0.20
    assert metrics["compiler_lowland_as_mountain_fraction"] <= 0.08
    assert metrics["compiler_terrain_elevation_sign_mismatch_fraction"] <= 0.02
    assert metrics["hypsometry_out_of_envelope_count"] >= 1
    assert metrics["province_missing_reference_class_count"] == 0
    assert metrics["boundary_missing_process_count"] == 0
    assert metrics["wilson_missing_object_set_count"] == 0
    assert metrics["source_to_sink_missing_object_count"] == 0
    assert metrics["drainage_missing_item_count"] == 0
    assert metrics["old_orogen_missing_decay_field_count"] == 0
    assert metrics["mountain_missing_inventory_field_count"] == 0
    assert metrics["plateau_missing_item_count"] == 0
    assert metrics["high_flat_interior_fraction_of_continental_land"] > 0.05
    assert 0.0 <= metrics["highland_without_parent_fraction_of_highlands"] <= 0.20
    assert metrics["basin_lowland_fraction_of_continental_land"] > 0.20
    assert metrics["lowland_fraction_lt500m_of_continental_land"] > 0.20
    assert metrics["lowland_fraction_lt1000m_of_continental_land"] > metrics["lowland_fraction_lt500m_of_continental_land"]
    assert metrics["major_component_count"] >= 2

    owner_counts = inventory["owner_counts"]
    assert owner_counts["planform"] >= 1
    assert owner_counts["crust_sediment"] >= 1
    assert "landform_expression" not in owner_counts
    assert "bathymetry_margin" not in owner_counts
    assert owner_counts["compiler_render"] >= 1
    assert inventory["category_counts"]["asset_review_pending"] == 8
    assert "terrain.plateau_inventory" not in inventory["current_residual_items"]
    assert "terrain.mountain_ranges" not in inventory["current_residual_items"]
    assert "terrain.drainage_basins" not in inventory["current_residual_items"]
    assert "terrain.rift_shoulders" not in inventory["current_residual_items"]
    assert "terrain.escarpments" not in inventory["current_residual_items"]
    assert "terrain.rift_margin_sequence_id" not in inventory["current_residual_items"]
    assert "terrain.rift_margin_stage" not in inventory["current_residual_items"]
    assert "tectonics.rift_margin_lineage_id" not in inventory["current_residual_items"]
    assert "tectonics.spreading_centers" not in inventory["current_residual_items"]
    assert "transform" not in inventory["current_residual_items"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P90.current_world_morphology_gap_inventory"}
    bench_metrics = benchmarks["P90.current_world_morphology_gap_inventory"]["metrics"]
    assert bench_metrics["inventory_status"] == "current_world_morphology_gap_inventory_ready"
    assert bench_metrics["gap_count"] == metrics["gap_count"]
    assert bench_metrics["owner_layer_count"] == 3
    assert bench_metrics["source_to_sink_missing_object_count"] == 0
    assert bench_metrics["drainage_missing_item_count"] == 0
    assert bench_metrics["old_orogen_missing_decay_field_count"] == 0
    assert bench_metrics["compiler_passed_envelope"]
    assert bench_metrics["p91_integrated_real_earth_morphology_promotion_audit_pending"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p90_microbenchmarks.csv").exists()


def test_p91_integrated_real_earth_morphology_promotion_audit_pass(tmp_path):
    summary = run_suite("P91", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p91.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["integrated_real_earth_morphology_promotion_audit"]
    assert summary["acceptance"]["p76_p90_stage_summaries_available"]
    assert summary["acceptance"]["p76_p90_stage_gates_pass"]
    assert summary["acceptance"]["p69_highres_review_assets_available"]
    assert summary["acceptance"]["p69_highres_contact_sheet_available"]
    assert summary["acceptance"]["ci_900_and_2500_worlds_audited"]
    assert summary["acceptance"]["ci_required_assets_available"]
    assert summary["acceptance"]["ci_contact_sheet_available"]
    assert summary["acceptance"]["ci_compiler_consistency_passed"]
    assert summary["acceptance"]["ci_gap_inventories_ready"]
    assert summary["acceptance"]["p90_residuals_mapped"]
    assert summary["acceptance"]["promotion_decision_recorded"]
    assert not summary["acceptance"]["promotion_ready"]
    assert summary["acceptance"]["keeps_default_off_until_named_residuals_resolved"]
    assert summary["acceptance"]["next_recommended_entry_defined"]

    audit = summary["integrated_real_earth_morphology_promotion_audit"]
    assert audit["schema"] == "aevum.integrated_real_earth_morphology_promotion_audit.v1"
    assert audit["status"] == "integrated_real_earth_morphology_promotion_audit_ready"

    metrics = audit["metrics"]
    assert metrics["stage_suite_count"] == 15
    assert metrics["stage_suite_pass_count"] == 15
    assert metrics["stage_missing_count"] == 0
    assert metrics["stage_failing_count"] == 0
    assert metrics["highres_member_count"] >= 3
    assert metrics["highres_required_asset_files"] == metrics["highres_actual_required_asset_files"]
    assert metrics["highres_contact_sheet_present_count"] >= 1
    assert metrics["ci_world_count"] == 2
    assert metrics["ci_required_asset_files"] == metrics["ci_actual_required_asset_files"] == 16
    assert metrics["ci_asset_set_complete_count"] == 2
    assert metrics["ci_compiler_passed_count"] == 2
    assert metrics["ci_inventory_ready_count"] == 2
    assert metrics["ci_contact_sheet_present"]
    assert metrics["p90_non_asset_gap_count"] >= 10
    assert metrics["p90_owner_layer_count"] >= 2
    assert metrics["p90_unassigned_gap_count"] == 0
    assert metrics["p90_generic_blocker_count"] == 0
    assert metrics["p90_compiler_passed_envelope"]
    assert metrics["promotion_blocker_count"] == 4
    assert metrics["audit_completed"]
    assert metrics["promotion_decision_recorded"]
    assert not metrics["promotion_ready"]

    decision = audit["promotion_decision"]
    assert decision["keeps_default_off_until_named_residuals_resolved"]
    assert decision["next_recommended_entry"] == "P92.production_residual_owner_repair_plan"
    assert "p90_current_world_residuals_unresolved" in decision["promotion_blockers"]
    assert "p69_earthlike_reference_needs_calibration" in decision["promotion_blockers"]
    assert "bathymetry_margin_residuals_unresolved" not in decision["promotion_blockers"]
    assert "landform_expression_residuals_unresolved" not in decision["promotion_blockers"]
    assert "drainage_erosion_residuals_unresolved" not in decision["promotion_blockers"]

    stage_rows = audit["stage_matrix"]["rows"]
    assert set(stage_rows) == {f"P{idx}" for idx in range(76, 91)}
    assert all(row["status"] == "pass" for row in stage_rows.values())
    assert audit["ci_review"]["cells"] == (900, 2500)
    assert audit["ci_review"]["contact_sheet_present"]
    assert len(audit["ci_worlds"]) == 2
    assert all(not row["missing_review_assets"] for row in audit["ci_worlds"])
    assert all(row["asset_set_complete"] for row in audit["ci_worlds"])

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P91.integrated_real_earth_morphology_promotion_audit"}
    bench_metrics = benchmarks["P91.integrated_real_earth_morphology_promotion_audit"]["metrics"]
    assert bench_metrics["audit_status"] == "integrated_real_earth_morphology_promotion_audit_ready"
    assert bench_metrics["stage_suite_pass_count"] == 15
    assert bench_metrics["ci_contact_sheet_present"]
    assert not bench_metrics["promotion_ready"]

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p91_microbenchmarks.csv").exists()
    assert (tmp_path / "p91_ci_world_contact_sheet.png").exists()
    assert (tmp_path / "ci_world_assets" / "earthlike_900cells" / "elevation.png").exists()
    assert (tmp_path / "ci_world_assets" / "earthlike_2500cells" / "hexmap.png").exists()


def test_p92_production_residual_owner_repair_plan_pass(tmp_path):
    summary = run_suite("P92", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p92.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["production_residual_owner_repair_plan"]
    assert summary["acceptance"]["p91_audit_available"]
    assert summary["acceptance"]["p91_promotion_blocked_not_ignored"]
    assert summary["acceptance"]["all_p91_blockers_assigned"]
    assert summary["acceptance"]["all_p90_owner_layers_assigned"]
    assert summary["acceptance"]["all_residual_items_assigned"]
    assert summary["acceptance"]["repair_packets_have_microbenchmarks"]
    assert summary["acceptance"]["repair_packets_have_acceptance_targets"]
    assert summary["acceptance"]["repair_packets_have_validation_suites"]
    assert summary["acceptance"]["repair_packets_have_implementation_targets"]
    assert summary["acceptance"]["dependencies_are_ordered"]
    assert summary["acceptance"]["climate_ocean_monsoon_work_excluded"]
    assert summary["acceptance"]["final_p91_reaudit_defined"]

    plan = summary["production_residual_owner_repair_plan"]
    assert plan["schema"] == "aevum.production_residual_owner_repair_plan.v1"
    assert plan["status"] == "production_residual_owner_repair_plan_ready"
    assert plan["unassigned_blockers"] == ()
    assert plan["unassigned_owner_layers"] == ()
    assert plan["unassigned_residual_items"] == ()
    assert plan["climate_targets"] == ()

    metrics = plan["metrics"]
    assert metrics["p91_summary_available"]
    assert metrics["p91_status"] == "pass"
    assert metrics["p91_audit_completed"]
    assert not metrics["p91_promotion_ready"]
    assert metrics["promotion_blocker_count"] == 9
    assert metrics["assigned_blocker_count"] == 9
    assert metrics["owner_layer_count"] == 7
    assert metrics["assigned_owner_layer_count"] == 7
    assert metrics["residual_item_count"] == 32
    assert metrics["assigned_residual_item_count"] == 32
    assert metrics["repair_packet_count"] == 8
    assert metrics["dependency_order_valid"]
    assert metrics["climate_target_count"] == 0
    assert metrics["packets_with_microbenchmarks"] == 8
    assert metrics["packets_with_acceptance_targets"] == 8
    assert metrics["packets_with_validation_suites"] == 8
    assert metrics["packets_with_implementation_targets"] == 8
    assert metrics["final_reaudit_packet_defined"]
    assert metrics["next_implementation_packet"] == "P92.1_planform_and_reference_calibration"
    assert metrics["final_validation_suite"] == "P91"

    packet_ids = tuple(packet["packet_id"] for packet in plan["repair_packets"])
    assert packet_ids == (
        "P92.1_planform_and_reference_calibration",
        "P92.2_production_province_graph_fields",
        "P92.3_boundary_lifecycle_objects",
        "P92.4_crust_sediment_interior_relief_coupling",
        "P92.5_drainage_source_to_sink_fields",
        "P92.6_landform_inventory_lifecycle",
        "P92.7_bathymetry_margin_sequence",
        "P92.8_integrated_reaudit_and_promotion_gate",
    )
    packets = {packet["packet_id"]: packet for packet in plan["repair_packets"]}
    assert packets["P92.1_planform_and_reference_calibration"]["depends_on"] == ()
    assert packets["P92.8_integrated_reaudit_and_promotion_gate"]["depends_on"] == packet_ids[:7]
    assert "P93.planform_reference_calibration" in packets[
        "P92.1_planform_and_reference_calibration"]["microbenchmarks_to_add"]
    assert "P97.production_drainage_source_to_sink_fields" in packets[
        "P92.5_drainage_source_to_sink_fields"]["microbenchmarks_to_add"]
    assert "P100.default_promotion_decision_gate" in packets[
        "P92.8_integrated_reaudit_and_promotion_gate"]["microbenchmarks_to_add"]

    blockers = set(plan["promotion_blockers"])
    assert blockers == {
        "p69_earthlike_reference_needs_calibration",
        "p90_current_world_residuals_unresolved",
        "bathymetry_margin_residuals_unresolved",
        "boundary_lifecycle_residuals_unresolved",
        "crust_sediment_residuals_unresolved",
        "drainage_erosion_residuals_unresolved",
        "landform_expression_residuals_unresolved",
        "planform_residuals_unresolved",
        "province_graph_residuals_unresolved",
    }
    assert plan["blocker_assignments"]["planform_residuals_unresolved"] == (
        "P92.1_planform_and_reference_calibration",)
    assert plan["blocker_assignments"]["p90_current_world_residuals_unresolved"] == (
        "P92.8_integrated_reaudit_and_promotion_gate",)
    assert plan["residual_item_assignments"]["terrain.drainage_basins"] == (
        "P92.5_drainage_source_to_sink_fields")
    assert plan["residual_item_assignments"]["terrain.plateau_inventory"] == (
        "P92.6_landform_inventory_lifecycle")
    assert plan["residual_item_assignments"]["terrain.rift_margin_sequence_id"] == (
        "P92.7_bathymetry_margin_sequence")

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P92.production_residual_owner_repair_plan"}
    bench_metrics = benchmarks["P92.production_residual_owner_repair_plan"]["metrics"]
    assert bench_metrics["plan_status"] == "production_residual_owner_repair_plan_ready"
    assert bench_metrics["repair_packet_count"] == 8
    assert bench_metrics["next_implementation_packet"] == (
        "P92.1_planform_and_reference_calibration")

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p92_microbenchmarks.csv").exists()


def test_p93_planform_reference_calibration_pass(tmp_path):
    summary = run_suite("P93", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p93.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["planform_reference_calibration"]
    assert summary["acceptance"]["generated_component_ribbon_envelope"]
    assert summary["acceptance"]["p92_planform_packet_available"]
    assert summary["acceptance"]["p69_p78_p90_evidence_available"]
    assert summary["acceptance"]["p91_blockers_preserved"]
    assert summary["acceptance"]["all_planform_reference_ranges_available"]
    assert summary["acceptance"]["p90_planform_gaps_covered"]
    assert summary["acceptance"]["calibration_directions_defined"]
    assert summary["acceptance"]["generated_900_and_8000_envelope_defined"]
    assert summary["acceptance"]["candidate_microbenchmarks_declared"]
    assert summary["acceptance"]["trench_gap_deferred_to_bathymetry_packet"]

    calibration = summary["planform_reference_calibration"]
    assert calibration["schema"] == "aevum.planform_reference_calibration.v1"
    assert calibration["status"] == "planform_reference_calibration_ready"
    assert calibration["p92_planform_packet"]["packet_id"] == (
        "P92.1_planform_and_reference_calibration")
    assert calibration["microbenchmarks_declared"] == (
        "P93.planform_reference_calibration",
        "P93.generated_component_ribbon_envelope",
    )

    metrics = calibration["metrics"]
    assert metrics["p92_packet_id"] == "P92.1_planform_and_reference_calibration"
    assert not metrics["p91_promotion_ready"]
    assert metrics["p91_blocker_count"] == 9
    assert metrics["p90_planform_gap_count"] == 9
    assert metrics["calibration_target_count"] == 5
    assert metrics["covered_planform_metric_count"] >= 4
    assert metrics["cross_owner_target_count"] == 1
    assert metrics["p69_member_count"] == 3
    assert metrics["p69_earthlike_out_of_envelope_count"] == 6
    assert metrics["p78_current_out_of_envelope_count"] == 5
    assert metrics["unresolved_primary_planform_metric_count"] >= 4
    assert metrics["generated_envelope_metric_count"] == 5
    assert metrics["next_packet_after_p93"] == "P92.2_production_province_graph_fields"

    targets = {target["metric"]: target for target in calibration["calibration_targets"]}
    assert set(targets) == {
        "land_fraction",
        "largest_land_component_fraction",
        "land_component_count",
        "land_ribbon_fraction_gt_0_5",
        "land_coastline_complexity_largest",
    }
    assert targets["land_fraction"]["current_900_direction"] == "increase"
    assert targets["land_component_count"]["current_900_direction"] == "increase"
    assert targets["land_ribbon_fraction_gt_0_5"]["current_900_direction"] == "decrease"
    assert targets["land_coastline_complexity_largest"]["current_900_direction"] == "decrease"
    assert targets["largest_land_component_fraction"]["current_900_direction"] == "unknown"
    assert targets["largest_land_component_fraction"]["p69_earthlike_direction"] == "decrease"

    cross_owner = {
        target["metric"]: target for target in calibration["cross_owner_targets"]
    }
    assert cross_owner["trench_fraction_of_ocean"]["repair_packet"] == (
        "P92.7_bathymetry_margin_sequence")

    envelope = calibration["generated_envelope"]
    assert envelope["land_component_count"]["current_900_value"] == 2.0
    assert envelope["land_ribbon_fraction_gt_0_5"]["current_900_value"] > 0.35
    assert envelope["land_coastline_complexity_largest"]["current_900_value"] > 8.0

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P93.planform_reference_calibration",
        "P93.generated_component_ribbon_envelope",
    }
    assert benchmarks["P93.planform_reference_calibration"]["metrics"][
        "calibration_status"] == "planform_reference_calibration_ready"
    assert benchmarks["P93.generated_component_ribbon_envelope"]["metrics"][
        "current_900_component_direction"] == "increase"

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p93_microbenchmarks.csv").exists()


def test_p94_production_province_graph_fields_pass(tmp_path):
    summary = run_suite("P94", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p94.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["production_province_graph_fields"]
    assert summary["acceptance"]["volcanic_lip_plateau_edge_coverage"]
    assert summary["acceptance"]["production_fields_present"]
    assert summary["acceptance"]["field_object_consistency"]
    assert summary["acceptance"]["parent_process_coverage"]
    assert summary["acceptance"]["major_continents_multi_province_900_2500"]
    assert summary["acceptance"]["required_reference_classes_covered"]
    assert summary["acceptance"]["volcanic_lip_plateau_and_rift_edge_covered"]
    assert summary["acceptance"]["p80_residuals_cleared"]
    assert summary["acceptance"]["no_checkerboard_province_graph"]
    assert summary["acceptance"]["p95_boundary_lifecycle_objects_pending"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P94.production_province_graph_fields",
        "P94.volcanic_lip_plateau_edge_coverage",
    }
    graph_metrics = benchmarks["P94.production_province_graph_fields"]["metrics"]
    assert graph_metrics["world_count"] == 2
    assert graph_metrics["ready_world_count"] == 2
    assert graph_metrics["min_object_count"] >= 8
    assert graph_metrics["min_field_id_count"] >= 8
    assert graph_metrics["min_class_count"] >= 9
    assert graph_metrics["min_parent_process_count"] >= 8
    assert graph_metrics["min_id_coverage"] >= 0.98
    assert graph_metrics["min_code_coverage"] >= 0.98
    assert graph_metrics["min_parent_process_coverage"] >= 0.98
    assert graph_metrics["max_missing_field_count"] == 0
    assert graph_metrics["max_missing_required_class_count"] == 0
    assert graph_metrics["max_missing_required_edge_count"] == 0
    assert graph_metrics["max_field_object_mismatch_count"] == 0
    assert graph_metrics["max_object_parent_failure_count"] == 0
    assert graph_metrics["max_disconnected_province_id_count"] == 0
    assert graph_metrics["max_tiny_province_area_fraction"] <= 0.02
    assert graph_metrics["min_major_component_count"] >= 2
    assert graph_metrics["min_major_component_province_id_count"] >= 4
    assert graph_metrics["min_major_component_province_class_count"] >= 5
    assert graph_metrics["max_major_component_largest_province_fraction"] <= 0.60
    assert graph_metrics["max_p80_missing_reference_class_count"] == 0
    assert graph_metrics["max_p80_missing_required_class_edge_count"] == 0

    lip_metrics = benchmarks["P94.volcanic_lip_plateau_edge_coverage"]["metrics"]
    assert lip_metrics["worlds_with_volcanic_lip_plateau"] == 2
    assert lip_metrics["worlds_with_required_rift_lip_edge"] == 2
    assert lip_metrics["min_required_edge_count"] > 0
    assert lip_metrics["min_volcanic_lip_area_fraction"] > 0.0

    graphs = summary["production_province_graphs"]
    assert set(graphs) == {"900", "2500"}
    for cells, graph in graphs.items():
        assert graph["status"] == "production_province_graph_ready", cells
        assert graph["missing_fields"] == ()
        assert graph["missing_required_classes"] == ()
        assert graph["missing_required_edges"] == ()
        assert "passive_margin_lowland" in graph["province_classes"]
        assert "volcanic_lip_plateau" in graph["province_classes"]
        assert int(graph["class_edge_counts"]["rift_system|volcanic_lip_plateau"]) > 0
        metrics = graph["metrics"]
        assert metrics["province_object_count"] >= 8
        assert metrics["province_class_count"] >= 9
        assert metrics["continental_id_coverage_fraction"] >= 0.98
        assert metrics["continental_code_coverage_fraction"] >= 0.98
        assert metrics["continental_parent_process_coverage_fraction"] >= 0.98
        assert metrics["missing_field_count"] == 0
        assert metrics["missing_required_class_count"] == 0
        assert metrics["missing_required_edge_count"] == 0
        assert metrics["missing_object_for_field_id_count"] == 0
        assert metrics["object_without_field_id_count"] == 0
        assert metrics["object_parent_failure_count"] == 0
        assert metrics["disconnected_province_id_count"] == 0
        assert metrics["min_major_component_province_id_count"] >= 4
        assert metrics["min_major_component_province_class_count"] >= 5
        assert metrics["tiny_province_area_fraction"] <= 0.02

    assert graphs["900"]["comparison_missing_reference_classes"] == ()
    assert graphs["900"]["comparison_missing_required_class_edges"] == ()
    assert graphs["900"]["metrics"]["p80_missing_reference_class_count"] == 0
    assert graphs["900"]["metrics"]["p80_missing_required_class_edge_count"] == 0

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p94_microbenchmarks.csv").exists()


def test_p95_boundary_lifecycle_objects_pass(tmp_path):
    summary = run_suite("P95", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p95.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["transform_and_spreading_center_objects"]
    assert summary["acceptance"]["boundary_lifecycle_current_world_audit"]
    assert summary["acceptance"]["transform_process_geometry_available"]
    assert summary["acceptance"]["ridge_and_transform_boundary_objects_available"]
    assert summary["acceptance"]["spreading_center_lifecycle_objects_available"]
    assert summary["acceptance"]["boundary_residuals_cleared"]
    assert summary["acceptance"]["wilson_lifecycle_residuals_cleared"]
    assert summary["acceptance"]["parent_lifecycle_links_preserved"]
    assert summary["acceptance"]["p96_crust_sediment_interior_relief_coupling_pending"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P95.transform_and_spreading_center_objects",
        "P95.boundary_lifecycle_current_world_audit",
    }
    transform = benchmarks["P95.transform_and_spreading_center_objects"]["metrics"]
    assert transform["world_count"] == 2
    assert transform["worlds_with_transform_process"] == 2
    assert transform["worlds_with_transform_boundary_objects"] == 2
    assert transform["worlds_with_ridge_boundary_objects"] == 2
    assert transform["worlds_with_spreading_centers"] == 2
    assert transform["min_transform_cell_count"] >= 4
    assert transform["min_ridge_cell_count"] > 0
    assert transform["min_transform_boundary_object_count"] > 0
    assert transform["min_ridge_boundary_object_count"] > 0
    assert transform["min_spreading_center_count"] > 0
    assert transform["min_transform_near_ridge_fraction"] >= 0.50
    assert transform["max_transform_length_fraction"] <= 0.25
    assert transform["max_missing_process_type_count"] == 0
    assert transform["max_missing_wilson_object_set_count"] == 0
    assert transform["max_parent_link_failure_count"] == 0

    audit = benchmarks["P95.boundary_lifecycle_current_world_audit"]["metrics"]
    assert audit["world_count"] == 2
    assert audit["boundary_ready_world_count"] == 2
    assert audit["wilson_ready_world_count"] == 2
    assert audit["min_boundary_process_type_count"] >= 7
    assert audit["min_active_phase_code_count"] >= 3
    assert audit["min_basin_stage_count"] >= 3
    assert audit["min_gateway_count"] > 0
    assert audit["min_wilson_cycle_count"] > 0
    assert audit["max_missing_process_type_count"] == 0
    assert audit["max_missing_wilson_object_set_count"] == 0
    assert audit["max_parent_link_failure_count"] == 0

    worlds = summary["boundary_lifecycle_worlds"]
    assert set(worlds) == {"900", "2500"}
    for cells, world in worlds.items():
        boundary = world["boundary"]["current_generated"]
        wilson = world["wilson"]["current_generated"]
        assert boundary["status"] == "generated_boundary_process_geometry_ready", cells
        assert wilson["status"] == "generated_world_wilson_lifecycle_ready", cells
        assert boundary["missing_process_types"] == ()
        assert wilson["missing_object_sets"] == ()
        assert int(world["boundary_object_kind_counts"]["ridge"]) > 0
        assert int(world["boundary_object_kind_counts"]["transform"]) > 0
        assert boundary["metrics"]["type_cell_counts"]["transform"] >= 4
        assert boundary["metrics"]["adjacency"]["transform_near_ridge_fraction"] >= 0.50
        assert wilson["object_counts"]["tectonics.spreading_centers"] > 0
        assert wilson["parent_link_failures"] == ()

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p95_microbenchmarks.csv").exists()


def test_p96_crust_sediment_surface_ordering_pass(tmp_path):
    summary = run_suite("P96", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p96.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["high_flat_interior_owner_reduction"]
    assert summary["acceptance"]["province_crust_sediment_surface_ordering"]
    assert summary["acceptance"]["high_flat_interiors_below_current_floor"]
    assert summary["acceptance"]["basin_lowland_area_preserved"]
    assert summary["acceptance"]["production_province_graphs_ready"]
    assert summary["acceptance"]["production_coupling_audits_ready"]
    assert summary["acceptance"]["platform_basin_orogen_surface_ordering"]
    assert summary["acceptance"]["sediment_accommodation_signal_preserved"]
    assert summary["acceptance"]["p97_drainage_source_to_sink_fields_pending"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P96.high_flat_interior_owner_reduction",
        "P96.province_crust_sediment_surface_ordering",
    }
    flat = benchmarks["P96.high_flat_interior_owner_reduction"]["metrics"]
    assert flat["world_count"] == 2
    assert flat["max_high_flat_interior_fraction_of_continental_land"] <= 0.080
    assert flat["max_p96_high_flat_fraction_after"] <= 0.020
    assert flat["min_basin_lowland_fraction_of_continental_land"] >= 0.20
    assert flat["min_lowland_fraction_lt500m_of_continental_land"] >= 0.20
    assert flat["min_lowland_fraction_lt1000m_of_continental_land"] >= 0.45
    assert flat["max_highland_without_parent_fraction_of_highlands"] <= 0.35
    assert flat["worlds_with_surface_ordering_telemetry"] == 2

    ordering = benchmarks["P96.province_crust_sediment_surface_ordering"]["metrics"]
    assert ordering["world_count"] == 2
    assert ordering["production_graph_ready_world_count"] == 2
    assert ordering["coupling_ready_world_count"] == 2
    assert ordering["min_province_class_count"] >= 9
    assert ordering["min_parent_process_count"] >= 9
    assert ordering["max_missing_required_class_count"] == 0
    assert ordering["max_missing_required_edge_count"] == 0
    assert ordering["max_p80_missing_reference_class_count"] == 0
    assert ordering["max_p80_missing_required_class_edge_count"] == 0
    assert ordering["min_platform_area_fraction"] >= 0.003
    assert ordering["min_intracratonic_basin_area_fraction"] > 0.0
    assert ordering["min_foreland_basin_area_fraction"] > 0.0
    assert ordering["max_basin_minus_platform_elevation_m"] < 0.0
    assert ordering["min_basin_minus_platform_sediment_m"] >= 500.0
    assert ordering["min_orogen_minus_basin_elevation_m"] >= 400.0
    assert ordering["production_aggregation_world_count"] == 2
    assert ordering["max_missing_production_field_count"] == 0

    worlds = summary["crust_sediment_surface_worlds"]
    assert set(worlds) == {"900", "2500"}
    for cells, world in worlds.items():
        graph = world["production_graph"]
        current = world["crust_sediment_coupling"]["current_generated"]
        direct = world["direct_surface_metrics"]
        telemetry = world["p96_terrain_telemetry"]

        assert graph["status"] == "production_province_graph_ready", cells
        assert graph["missing_required_classes"] == ()
        assert graph["missing_required_edges"] == ()
        assert graph["metrics"]["p80_missing_reference_class_count"] == 0
        assert graph["metrics"]["p80_missing_required_class_edge_count"] == 0
        assert graph["metrics"]["tiny_province_area_fraction"] <= 0.02
        assert current["status"] == "generated_world_crust_sediment_province_audit_ready"
        assert current["aggregation_source"] == "tectonics.continental_provinces"
        assert current["missing_production_province_fields"] == ()
        assert direct["high_flat_interior_fraction_of_continental_land"] <= 0.080
        assert direct["basin_lowland_fraction_of_continental_land"] >= 0.20
        assert "terrain.last_p96_high_flat_fraction_after" in telemetry

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p96_microbenchmarks.csv").exists()


def test_p97_drainage_source_to_sink_fields_pass(tmp_path):
    summary = run_suite("P97", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p97.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["production_drainage_source_to_sink_fields"]
    assert summary["acceptance"]["old_orogen_decay_budget_current_world"]
    assert summary["acceptance"]["source_to_sink_worlds_ready"]
    assert summary["acceptance"]["drainage_worlds_ready"]
    assert summary["acceptance"]["old_orogen_worlds_ready"]
    assert summary["acceptance"]["source_to_sink_objects_present"]
    assert summary["acceptance"]["drainage_fields_present"]
    assert summary["acceptance"]["sediment_budget_closed"]
    assert summary["acceptance"]["regional_drainage_basins_not_fragmented"]
    assert summary["acceptance"]["old_orogen_decay_fields_present"]
    assert summary["acceptance"]["p98_landform_inventory_lifecycle_pending"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P97.production_drainage_source_to_sink_fields",
        "P97.old_orogen_decay_budget_current_world",
    }
    drainage = benchmarks["P97.production_drainage_source_to_sink_fields"]["metrics"]
    assert drainage["world_count"] == 2
    assert drainage["source_to_sink_ready_world_count"] == 2
    assert drainage["drainage_ready_world_count"] == 2
    assert drainage["max_missing_source_to_sink_object_count"] == 0
    assert drainage["max_missing_drainage_item_count"] == 0
    assert drainage["min_drainage_basin_object_count"] >= 8
    assert drainage["min_routing_edge_count"] >= 12
    assert drainage["min_routing_source_kind_count"] >= 2
    assert drainage["min_routing_sink_kind_count"] >= 2
    assert drainage["max_source_sink_balance_fraction"] <= 1.0e-9
    assert drainage["min_drainage_basin_field_count"] >= 8
    assert drainage["max_divide_fraction_of_land"] <= 0.55
    assert drainage["min_divide_alignment_fraction"] >= 0.70
    assert drainage["max_major_basin_component_failure_count"] == 0
    assert drainage["min_flow_to_sink_consistency_fraction"] >= 0.98
    assert drainage["min_downhill_path_fraction"] >= 0.98

    old_orogen = benchmarks["P97.old_orogen_decay_budget_current_world"]["metrics"]
    assert old_orogen["world_count"] == 2
    assert old_orogen["old_orogen_ready_world_count"] == 2
    assert old_orogen["max_missing_decay_field_count"] == 0
    assert old_orogen["min_old_subdued_orogen_object_count"] > 0
    assert old_orogen["min_parented_old_subdued_orogen_object_count"] > 0
    assert old_orogen["min_old_orogen_decay_stage_count"] > 0
    assert old_orogen["min_mean_orogen_boundary_memory"] >= 0.50
    assert old_orogen["min_mean_orogen_erosion_budget_m"] > 0.0
    assert old_orogen["min_orogen_sediment_export_volume_km3"] > 0.0

    worlds = summary["drainage_source_to_sink_worlds"]
    assert set(worlds) == {"900", "2500"}
    for cells, world in worlds.items():
        source = world["source_to_sink"]["current_generated"]
        drainage_current = world["drainage"]["current_generated"]
        old_current = world["old_orogen"]["current_generated"]
        assert world["source_to_sink"]["status"] == "source_to_sink_sediment_budget_ready", cells
        assert world["drainage"]["status"] == "drainage_divide_province_alignment_ready", cells
        assert world["old_orogen"]["status"] == "old_orogen_erosion_decay_ready", cells
        assert source["missing_source_to_sink_objects"] == ()
        assert drainage_current["missing_drainage_items"] == ()
        assert old_current["missing_decay_fields"] == ()
        assert source["metrics"]["production_volume_balance_fraction"] <= 1.0e-9
        assert drainage_current["metrics"]["major_basin_component_failure_count"] == 0
        assert drainage_current["metrics"]["current_flow_to_sink_consistency_fraction"] >= 0.98
        assert old_current["metrics"]["current_orogen_sediment_export_volume_km3"] > 0.0
        assert world["terrain_globals"]["terrain.last_p97_drainage_basin_count"] >= 8.0

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p97_microbenchmarks.csv").exists()


def test_p98_landform_inventory_lifecycle_pass(tmp_path):
    summary = run_suite("P98", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p98.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["production_mountain_inventory_fields"]
    assert summary["acceptance"]["production_plateau_lifecycle_fields"]
    assert summary["acceptance"]["mountain_worlds_ready"]
    assert summary["acceptance"]["plateau_worlds_ready"]
    assert summary["acceptance"]["mountain_inventory_fields_present"]
    assert summary["acceptance"]["mountain_expected_kinds_present"]
    assert summary["acceptance"]["mountain_ranges_finite_and_shaped"]
    assert summary["acceptance"]["plateau_lifecycle_fields_present"]
    assert summary["acceptance"]["plateau_expected_kinds_present"]
    assert summary["acceptance"]["volcanic_lip_plateau_present"]
    assert summary["acceptance"]["plateaus_area_capped"]
    assert summary["acceptance"]["p99_bathymetry_margin_sequence_pending"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P98.production_mountain_inventory_fields",
        "P98.production_plateau_lifecycle_fields",
    }
    mountain = benchmarks["P98.production_mountain_inventory_fields"]["metrics"]
    assert mountain["world_count"] == 2
    assert mountain["mountain_ready_world_count"] == 2
    assert mountain["max_missing_inventory_field_count"] == 0
    assert mountain["max_missing_expected_kind_count"] == 0
    assert mountain["min_production_mountain_range_object_count"] >= 4
    assert mountain["min_mountain_field_id_count"] >= 4
    assert mountain["min_mountain_inventory_class_count"] >= 3
    assert mountain["min_max_mountain_elongation_ratio"] >= 1.65
    assert mountain["max_total_mountain_area_fraction_world"] <= 0.20
    assert mountain["max_mountain_object_area_fraction_world"] <= 0.08

    plateau = benchmarks["P98.production_plateau_lifecycle_fields"]["metrics"]
    assert plateau["world_count"] == 2
    assert plateau["plateau_ready_world_count"] == 2
    assert plateau["max_missing_plateau_item_count"] == 0
    assert plateau["max_missing_plateau_kind_count"] == 0
    assert plateau["min_plateau_inventory_cell_count"] > 0
    assert plateau["min_volcanic_plateau_cell_count"] > 0
    assert plateau["min_plateau_decay_stage_count"] > 0
    assert plateau["max_plateau_area_fraction_world"] <= 0.060
    assert plateau["max_plateau_object_area_fraction_world"] <= 0.060
    assert plateau["max_single_plateau_object_area_fraction_world"] <= 0.035
    assert plateau["max_high_interior_without_plateau_fraction"] <= 0.10

    worlds = summary["landform_inventory_lifecycle_worlds"]
    assert set(worlds) == {"900", "2500"}
    for cells, world in worlds.items():
        assert world["mountain"]["status"] == "mountain_inventory_expression_ready", cells
        assert world["plateau"]["status"] == "plateau_area_cap_and_decay_ready", cells
        assert world["mountain"]["current_generated"]["missing_inventory_fields"] == ()
        assert world["mountain"]["current_generated"]["missing_expected_mountain_kinds"] == ()
        assert world["plateau"]["current_generated"]["missing_plateau_items"] == ()
        assert world["plateau"]["current_generated"]["missing_expected_plateau_kinds"] == ()
        assert world["field_metrics"]["mountain_field_id_count"] >= 4
        assert world["field_metrics"]["plateau_inventory_cell_count"] > 0

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p98_microbenchmarks.csv").exists()


def test_p99_bathymetry_margin_sequence_pass(tmp_path):
    summary = run_suite("P99", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p99.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["production_rift_margin_sequence_fields"]
    assert summary["acceptance"]["production_bathymetry_margin_ordering"]
    assert summary["acceptance"]["sequence_worlds_ready"]
    assert summary["acceptance"]["sequence_fields_present"]
    assert summary["acceptance"]["rift_shoulders_and_escarpments_present"]
    assert summary["acceptance"]["margin_sequence_objects_present"]
    assert summary["acceptance"]["shelf_slope_rise_abyss_ordered"]
    assert summary["acceptance"]["p100_integrated_reaudit_and_promotion_gate_pending"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P99.production_rift_margin_sequence_fields",
        "P99.production_bathymetry_margin_ordering",
    }

    sequence = benchmarks["P99.production_rift_margin_sequence_fields"]["metrics"]
    assert sequence["world_count"] == 2
    assert sequence["sequence_ready_world_count"] == 2
    assert sequence["max_missing_sequence_item_count"] == 0
    assert sequence["min_sequence_id_count"] > 0
    assert sequence["min_lineage_id_count"] > 0
    assert sequence["min_stage_count"] >= 6
    assert sequence["min_rift_shoulder_cell_count"] > 0
    assert sequence["min_escarpment_cell_count"] > 0
    assert sequence["min_sequence_object_count"] > 0
    assert sequence["min_rift_shoulder_object_count"] > 0
    assert sequence["min_escarpment_object_count"] > 0
    assert sequence["min_passive_lowland_stage_cell_count"] > 0

    bathymetry = benchmarks["P99.production_bathymetry_margin_ordering"]["metrics"]
    assert bathymetry["world_count"] == 2
    assert bathymetry["ordered_world_count"] == 2
    assert bathymetry["min_shelf_stage_cell_count"] > 0
    assert bathymetry["min_slope_stage_cell_count"] > 0
    assert bathymetry["min_rise_stage_cell_count"] > 0
    assert bathymetry["min_abyss_stage_cell_count"] > 0
    assert bathymetry["max_shelf_depth_p75_m"] < 800.0
    assert bathymetry["min_shelf_to_abyss_depth_delta_m"] > 1000.0

    worlds = summary["bathymetry_margin_sequence_worlds"]
    assert set(worlds) == {"900", "2500"}
    for cells, world in worlds.items():
        sequence_world = world["rift_margin_escarpment_sequence"]
        current = sequence_world["current_generated"]
        field = world["field_metrics"]
        assert sequence_world["status"] == "rift_margin_escarpment_sequence_ready", cells
        assert current["missing_sequence_items"] == (), cells
        assert current["acceptance"]["production_rift_margin_sequence_fields_available"], cells
        assert current["acceptance"]["production_rift_margin_sequence_ready"], cells
        assert current["acceptance"]["production_shelf_slope_rise_abyss_ordered"], cells
        assert field["sequence_id_count"] > 0
        assert field["lineage_id_count"] > 0
        assert field["stage_count"] >= 6
        assert field["rift_shoulder_cell_count"] > 0
        assert field["escarpment_cell_count"] > 0
        assert field["sequence_object_count"] > 0
        assert field["rift_shoulder_object_count"] > 0
        assert field["escarpment_object_count"] > 0
        assert field["shelf_stage_cell_count"] > 0
        assert field["slope_stage_cell_count"] > 0
        assert field["rise_stage_cell_count"] > 0
        assert field["abyss_stage_cell_count"] > 0
        assert (
            field["shelf_depth_p75_m"]
            < field["slope_depth_p50_m"]
            < field["rise_depth_p50_m"]
            < field["abyss_depth_p50_m"]
        )

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p99_microbenchmarks.csv").exists()


def test_p100_integrated_reaudit_and_promotion_gate_pass(tmp_path):
    summary = run_suite("P100", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p100.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["integrated_owner_repair_reaudit"]
    assert summary["acceptance"]["default_promotion_decision_gate"]
    assert summary["acceptance"]["p93_p99_repair_suites_available"]
    assert summary["acceptance"]["p93_p99_repair_suites_pass"]
    assert summary["acceptance"]["p91_after_p99_reaudit_pass"]
    assert summary["acceptance"]["completed_owner_blockers_cleared_from_root"]
    assert summary["acceptance"]["remaining_blockers_are_named"]
    assert summary["acceptance"]["default_promotion_blocked"]
    assert summary["acceptance"]["next_action_targets_planform_crust_sediment"]

    gate = summary["integrated_reaudit_and_promotion_gate"]
    assert gate["schema"] == "aevum.p100_integrated_reaudit_and_promotion_gate.v1"
    assert gate["status"] == "integrated_reaudit_and_promotion_gate_ready"
    assert all(gate["acceptance"].values())

    matrix = gate["repair_suite_matrix"]
    assert matrix["schema"] == "aevum.p100_repair_suite_matrix.v1"
    rows = {row["suite"]: row for row in matrix["rows"]}
    assert set(rows) == {"P93", "P94", "P95", "P96", "P97", "P98", "P99"}
    assert matrix["missing_suites"] == ()
    assert matrix["failing_suites"] == ()
    assert all(row["status"] == "pass" for row in rows.values())
    assert all(row["all_microbenchmarks_pass"] for row in rows.values())
    assert rows["P93"]["packet_id"] == "P92.1_planform_and_reference_calibration"
    assert rows["P99"]["packet_id"] == "P92.7_bathymetry_margin_sequence"

    reaudit = gate["after_p99_p91_reaudit"]
    assert reaudit["status"] == "pass"
    assert reaudit["promotion_blockers"] == (
        "p69_earthlike_reference_needs_calibration",
        "p90_current_world_residuals_unresolved",
        "crust_sediment_residuals_unresolved",
        "planform_residuals_unresolved",
    )
    assert reaudit["remaining_owner_blockers"] == (
        "crust_sediment_residuals_unresolved",
        "planform_residuals_unresolved",
    )
    assert reaudit["cleared_root_owner_blockers"] == (
        "bathymetry_margin_residuals_unresolved",
        "boundary_lifecycle_residuals_unresolved",
        "drainage_erosion_residuals_unresolved",
        "landform_expression_residuals_unresolved",
        "province_graph_residuals_unresolved",
    )
    assert set(reaudit["root_p90_owner_counts"]) == {"crust_sediment", "planform"}
    assert set(reaudit["ci_owner_counts"][900]) == {"crust_sediment", "planform"}
    assert set(reaudit["ci_owner_counts"][2500]).issuperset({
        "crust_sediment",
        "planform",
    })

    decision = gate["promotion_decision"]
    assert not decision["release_gate_allowed"]
    assert not decision["default_promotion_ready"]
    assert decision["default_promotion_blocked"]
    assert decision["next_recommended_action"] == (
        "P101.planform_crust_sediment_residual_repair")

    metrics = gate["metrics"]
    assert metrics["repair_suite_count"] == 7
    assert metrics["repair_suite_pass_count"] == 7
    assert metrics["missing_repair_suite_count"] == 0
    assert metrics["failing_repair_suite_count"] == 0
    assert metrics["p91_after_p99_status"] == "pass"
    assert metrics["p91_after_p99_deterministic"]
    assert metrics["p91_stage_suite_count"] == 15
    assert metrics["p91_stage_suite_pass_count"] == 15
    assert metrics["p91_ci_world_count"] == 2
    assert metrics["p91_ci_asset_set_complete_count"] == 2
    assert metrics["p91_ci_compiler_passed_count"] == 2
    assert metrics["promotion_blocker_count"] == 4
    assert metrics["expected_blocker_set_matched"]
    assert metrics["remaining_owner_blocker_count"] == 2
    assert metrics["cleared_root_owner_blocker_count"] == 5
    assert metrics["root_p90_non_asset_gap_count"] == 11
    assert metrics["root_p90_owner_layer_count"] == 2
    assert metrics["root_p90_residual_item_count"] == 0
    assert not metrics["release_gate_allowed"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P100.integrated_owner_repair_reaudit",
        "P100.default_promotion_decision_gate",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p100_microbenchmarks.csv").exists()
    assert (tmp_path / "p91_reaudit" / "tectonics_bench_summary.json").exists()


def test_p101_current_residual_attribution_phase0_pass(tmp_path):
    summary = run_suite("P101", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p101.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["planform_residual_baseline"]
    assert summary["acceptance"]["crust_sediment_high_flat_repair"]
    assert summary["acceptance"]["phase0_baseline_reproduced"]
    assert summary["acceptance"]["remaining_owner_layers_match_p100"]
    assert summary["acceptance"]["all_expected_non_asset_gaps_present"]
    assert summary["acceptance"]["all_non_asset_gaps_exactly_attributed"]
    assert summary["acceptance"]["multi_resolution_audit_recorded"]
    assert summary["acceptance"]["repair_targets_recorded_without_threshold_relaxation"]
    assert summary["acceptance"]["phase1_reference_evidence_expansion_pending"]
    assert summary["acceptance"]["phase2_planform_mechanism_repair_pending"]
    assert summary["acceptance"]["phase3_crust_sediment_interior_elevation_repair_pending"]

    attribution = summary["current_residual_attribution"]
    assert attribution["schema"] == "aevum.p101_current_residual_attribution.v1"
    assert attribution["status"] == "current_residual_attribution_ready"
    assert attribution["phase"] == "Phase 0. Baseline Reproduction and Failure Attribution"
    assert attribution["evidence_packet_template"]["packet_id"] == (
        "P101_phase0_current_residual_attribution")
    assert set(attribution["per_resolution"]) == {900, 2500}

    metrics = attribution["metrics"]
    assert metrics["root_gap_count"] >= 19
    assert metrics["root_non_asset_gap_count"] == 11
    assert metrics["root_asset_review_gap_count"] == 8
    assert metrics["root_planform_gap_count"] == 10
    assert metrics["root_crust_sediment_gap_count"] == 1
    assert metrics["expected_gap_id_match_count"] == metrics["expected_non_asset_gap_count"] == 11
    assert metrics["missing_expected_gap_count"] == 0
    assert metrics["unexpected_non_asset_gap_count"] == 0
    assert metrics["attributed_non_asset_gap_count"] == 11
    assert metrics["exact_attribution_count"] == 11
    assert metrics["fallback_attribution_count"] == 0
    assert metrics["root_continental_land_fraction_world"] < metrics["land_fraction_target_min"]
    assert metrics["root_major_component_count"] < metrics["major_component_target_min"]
    assert (
        metrics["root_high_flat_interior_fraction_of_continental_land"]
        > metrics["high_flat_repair_threshold"]
    )
    assert metrics["high_flat_repair_threshold"] == 0.02
    assert metrics["land_fraction_target_min"] == 0.25
    assert metrics["land_fraction_target_max"] == 0.33
    assert metrics["major_component_target_min"] == 4
    assert metrics["largest_component_target_max"] == 0.60
    assert metrics["land_ribbon_target_max"] == 0.35

    rows = attribution["attribution_rows"]
    assert len(rows) == 11
    assert all(row["exact_attribution"] for row in rows)
    assert {row["owner_layer"] for row in rows} == {"planform", "crust_sediment"}
    assert sum(1 for row in rows if row["owner_layer"] == "planform") == 10
    assert sum(1 for row in rows if row["owner_layer"] == "crust_sediment") == 1
    assert all(row["code_targets"] for row in rows)
    assert any(
        "aevum/modules/tectonics.py::_balanced_planform_reshape" in row["code_targets"]
        for row in rows
        if row["owner_layer"] == "planform"
    )
    assert any(
        "aevum/modules/terrain.py::_province_crust_sediment_surface_ordering"
        in row["code_targets"]
        for row in rows
        if row["owner_layer"] == "crust_sediment"
    )

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P101.planform_residual_baseline",
        "P101.crust_sediment_high_flat_repair",
    }
    assert all(bench["passed"] for bench in benchmarks.values())
    assert benchmarks["P101.planform_residual_baseline"]["metrics"][
        "root_planform_gap_count"] == 10
    assert benchmarks["P101.crust_sediment_high_flat_repair"]["metrics"][
        "root_crust_sediment_gap_count"] == 1

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p101_microbenchmarks.csv").exists()


def test_p102_reference_evidence_packet_matrix_pass(tmp_path):
    summary = run_suite("P102", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p102.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["reference_evidence_packet_matrix"]
    assert summary["acceptance"]["packet_count_sufficient"]
    assert summary["acceptance"]["all_source_ids_valid"]
    assert summary["acceptance"]["all_packets_complete"]
    assert summary["acceptance"]["reference_and_generated_tracks_present"]
    assert summary["acceptance"]["raw_data_policy_explicit"]
    assert summary["acceptance"]["core_metric_groups_covered"]
    assert summary["acceptance"]["p101_residual_owners_covered"]
    assert summary["acceptance"]["phase1_source_coverage_expanded"]
    assert summary["acceptance"]["phase2_planform_mechanism_repair_ready_to_plan"]

    packets = summary["reference_evidence_packets"]
    assert packets["schema"] == "aevum.reference_evidence_packets.v1"
    assert packets["status"] == "reference_evidence_packets_ready"
    assert packets["source_ledger_status"] == "source_ledger_schema_ready"
    assert packets["hypsometry_fixture_status"] == "hypsometry_fixture_ready"
    assert packets["province_reference_graph_status"] == "province_reference_graph_ready"
    assert packets["packet_count"] == 6
    assert packets["source_id_count"] >= 12
    assert packets["metric_key_count"] >= 24
    assert set(packets["covered_metric_groups"]).issuperset({
        "planform",
        "hypsometry",
        "province_architecture",
        "process_parentage",
        "drainage",
        "ocean_bathymetry",
    })
    assert {"planform", "crust_sediment"}.issubset(
        set(packets["residual_owner_layers"]))
    assert packets["duplicate_packet_ids"] == ()
    assert packets["invalid_source_refs"] == {}
    assert packets["missing_required_fields"] == {}
    assert packets["incomplete_packets"] == ()
    assert all(packets["acceptance"].values())

    packet_rows = {packet["packet_id"]: packet for packet in packets["packets"]}
    assert set(packet_rows) == {
        "R1_global_hypsometry_planform",
        "R2_province_crust_sediment_basement",
        "R3_boundary_wilson_deeptime",
        "R4_drainage_erosion_source_to_sink",
        "R5_landform_margins_mountains_plateaus",
        "R6_case_study_feature_catalog",
    }
    for packet in packet_rows.values():
        assert packet["schema"] == "aevum.reference_evidence_packet.v1"
        assert packet["source_ids"]
        assert packet["theory_claims"]
        assert packet["derived_metrics"]
        assert packet["reference_fixture"]
        assert packet["generated_world_audit"]
        assert packet["optimization_targets"]
        assert packet["residual_policy"]
        assert packet["residual_owner_layers"]
        assert packet["missing_required_fields"] == ()
        assert packet["invalid_source_ids"] == ()
        assert all(packet["acceptance"].values())

    assert "P101.planform_residual_baseline" in packet_rows[
        "R1_global_hypsometry_planform"]["generated_world_audit"]
    assert "P101.crust_sediment_high_flat_repair" in packet_rows[
        "R2_province_crust_sediment_basement"]["generated_world_audit"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P102.reference_evidence_packet_matrix"}
    assert benchmarks["P102.reference_evidence_packet_matrix"]["passed"]
    metrics = benchmarks["P102.reference_evidence_packet_matrix"]["metrics"]
    assert metrics["packet_status"] == "reference_evidence_packets_ready"
    assert metrics["packet_count"] == 6
    assert metrics["invalid_source_ref_packet_count"] == 0
    assert metrics["missing_required_field_packet_count"] == 0
    assert metrics["incomplete_packet_count"] == 0

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p102_microbenchmarks.csv").exists()


def test_p103_planform_mechanism_repair_pass(tmp_path):
    summary = run_suite("P103", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p103.v1"
    assert summary["status"] == "pass"
    assert summary["phase"] == "Phase 2. Planform Mechanism Repair"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["reference_land_floor_guard"]
    assert summary["acceptance"]["planform_repair_progress"]
    assert summary["acceptance"]["p102_reference_evidence_consumed"]
    assert summary["acceptance"]["multi_resolution_audit_recorded"]
    assert summary["acceptance"]["land_fraction_residual_cleared"]
    assert summary["acceptance"]["no_new_root_unexpected_gaps"]
    assert summary["acceptance"]["compiler_consistency_preserved"]
    assert summary["acceptance"]["phase2_remaining_planform_residuals_recorded"]
    assert summary["acceptance"]["phase3_crust_sediment_interior_elevation_repair_pending"]

    rows = {int(row["cells"]): row for row in summary["planform_repair_rows"]}
    assert set(rows) == {900, 2500}
    root = rows[900]
    assert root["continental_land_fraction_world"] >= 0.25
    assert "planform.land_fraction_out_of_envelope" not in root["non_asset_gap_ids"]
    assert root["planform_gap_count"] < 10
    assert root["new_unexpected_gap_ids"] == ()
    assert root["compiler_passed_envelope"]
    assert root["largest_shave_floor_fraction"] >= 0.25
    assert (
        root["largest_shave_land_fraction_after"] >= 0.25
        or root["largest_shave_rejected_by_land_floor"]
        or root["largest_shave_area_fraction"] <= 0.0
    )
    assert root["post_planform_solve_component_count"] >= 4
    assert root["post_coast_smoothing_land_fraction"] >= 0.25

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "P103.reference_land_floor_guard",
        "P103.planform_repair_progress",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    assert (tmp_path / "tectonics_bench_summary.json").exists()
    assert (tmp_path / "p103_microbenchmarks.csv").exists()


def test_p104a_continental_mosaic_expression_fixture():
    import numpy as np

    from aevum.core.grid import SphereGrid
    from aevum.core.state import WorldState
    from aevum.modules.terrain import (
        CONT_DETAIL_BASIN,
        CONT_DETAIL_OROGEN,
        CONT_DETAIL_PLATFORM,
        CONT_DETAIL_PLATEAU,
        CONT_DETAIL_RIFT_BASIN,
        CONT_DETAIL_SHIELD,
        CONT_PROVINCE_ACTIVE_OROGEN,
        CONT_PROVINCE_INTRACRATONIC_BASIN,
        CONT_PROVINCE_OLD_OROGEN,
        CONT_PROVINCE_PLATFORM,
        CONT_PROVINCE_RIFT_SYSTEM,
        CONT_PROVINCE_SHIELD,
        CONT_PROVINCE_VOLCANIC_LIP_PLATEAU,
        TerrainModule,
    )
    from aevum.spec.presets import get_preset

    grid = SphereGrid.fibonacci(900, 6.371e6)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    continent = (np.abs(grid.lat) < 48.0) & (grid.lon > -155.0) & (grid.lon < 155.0)

    def cap(lat, lon, radius):
        target = np.array([
            np.cos(np.radians(lat)) * np.cos(np.radians(lon)),
            np.cos(np.radians(lat)) * np.sin(np.radians(lon)),
            np.sin(np.radians(lat)),
        ])
        dots = np.clip(grid.xyz @ target, -1.0, 1.0)
        return np.degrees(np.arccos(dots)) <= radius

    shield = cap(18.0, -112.0, 18.0) & continent
    basin = cap(-18.0, -66.0, 18.0) & continent & ~shield
    rift = (np.abs(grid.lon + 18.0) < 9.0) & continent & ~(shield | basin)
    old_orogen = (np.abs(grid.lon - 28.0) < 10.0) & continent & ~(shield | basin | rift)
    active_orogen = (np.abs(grid.lon - 78.0) < 10.0) & continent & ~(shield | basin | rift | old_orogen)
    plateau = cap(20.0, 118.0, 18.0) & continent & ~(shield | basin | rift | old_orogen | active_orogen)
    platform = continent & ~(shield | basin | rift | old_orogen | active_orogen | plateau)

    surface = np.where(continent, 780.0, -4200.0).astype(float)
    surface[basin] = 520.0
    surface[rift] = 620.0
    surface[old_orogen] = 1320.0
    surface[active_orogen] = 1640.0
    surface[plateau] = 1780.0
    crust_type = np.where(continent, 1.0, 0.0)
    terrain_province = np.where(continent, 2.0, 0.0)
    detail = np.where(continent, float(CONT_DETAIL_PLATFORM), 0.0)
    province_code = np.zeros(grid.n, dtype=float)
    province_code[shield] = CONT_PROVINCE_SHIELD
    province_code[basin] = CONT_PROVINCE_INTRACRATONIC_BASIN
    province_code[rift] = CONT_PROVINCE_RIFT_SYSTEM
    province_code[old_orogen] = CONT_PROVINCE_OLD_OROGEN
    province_code[active_orogen] = CONT_PROVINCE_ACTIVE_OROGEN
    province_code[plateau] = CONT_PROVINCE_VOLCANIC_LIP_PLATEAU
    province_code[platform] = CONT_PROVINCE_PLATFORM

    world.set_field("crust.thickness_m", np.where(continent, 70000.0, 7000.0))
    thick = world.field("crust.thickness_m")
    thick[basin] = 61000.0
    thick[rift] = 64000.0
    thick[old_orogen] = 77500.0
    thick[active_orogen] = 83500.0
    thick[plateau] = 87500.0
    world.set_field("crust.stability", np.where(continent, 0.62, 0.0))
    stability = world.field("crust.stability")
    stability[shield] = 0.92
    stability[basin] = 0.48
    stability[rift] = 0.36
    stability[active_orogen] = 0.42
    world.set_field("sediment.thickness_m", np.where(continent, 900.0, 200.0))
    sediment = world.field("sediment.thickness_m")
    sediment[basin] = 3400.0
    sediment[rift] = 2200.0
    world.set_field("tectonics.rift_potential", np.where(rift, 0.92, 0.0))

    out = terrain._apply_continental_mosaic_expression(
        world,
        surface,
        0.0,
        crust_type,
        terrain_province,
        detail,
        {"province_code": province_code},
    )
    detail_out = out["continental_detail"].astype(int)
    surface_out = out["surface"]

    assert np.mean(detail_out[shield] == CONT_DETAIL_SHIELD) > 0.85
    assert np.mean(detail_out[basin] == CONT_DETAIL_BASIN) > 0.85
    assert np.mean(detail_out[rift] == CONT_DETAIL_RIFT_BASIN) > 0.85
    assert np.mean(detail_out[old_orogen] == CONT_DETAIL_OROGEN) > 0.85
    assert np.mean(detail_out[active_orogen] == CONT_DETAIL_OROGEN) > 0.85
    assert np.mean(detail_out[plateau] == CONT_DETAIL_PLATEAU) > 0.85
    platform_detail_fraction = float(np.mean(detail_out[platform] == CONT_DETAIL_PLATFORM))
    platform_split_fraction = float(np.mean(
        np.isin(detail_out[platform], [CONT_DETAIL_SHIELD, CONT_DETAIL_BASIN])
    ))
    assert platform_detail_fraction > 0.65
    assert platform_split_fraction > 0.05
    assert platform_detail_fraction + platform_split_fraction > 0.95

    platform_median = float(np.median(surface_out[platform]))
    assert float(np.median(surface_out[basin])) < platform_median - 90.0
    assert float(np.median(surface_out[rift])) < platform_median - 40.0
    assert float(np.median(surface_out[old_orogen])) > platform_median + 120.0
    assert float(np.median(surface_out[active_orogen])) > platform_median + 220.0
    assert float(np.median(surface_out[plateau])) > platform_median + 250.0
    assert np.all(surface_out[continent] >= 30.0)
    assert world.g("terrain.last_p104a_detail_reclassified_area_fraction") > 0.10
    assert world.g("terrain.last_p104a_mean_abs_surface_delta_m") > 50.0


def test_p104a_generated_world_continental_mosaic_gate_pass(tmp_path):
    summary = run_suite("P104A", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p104a.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["generated_world_continental_mosaic_gate"]
    assert summary["acceptance"]["legacy_gate_separation"]
    assert summary["acceptance"]["per_major_component_internal_mosaic_pass"]

    rows = summary["generated_world_continental_mosaic_rows"]
    assert {row["seed"] for row in rows} == {42, 31415, 16180}
    assert all(row["cells"] == 900 for row in rows)
    assert all(row["p104a_interior_mosaic_pass"] for row in rows)
    assert min(row["min_province_class_count_per_major"] for row in rows) >= 3
    assert max(row["max_largest_internal_province_fraction"] for row in rows) <= 0.74
    assert min(row["p104a_detail_class_count_exposed_continental"] for row in rows) >= 5
    assert min(row["p104a_province_class_count_exposed_continental"] for row in rows) >= 6

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p104a_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    csv_rows = list(csv.DictReader(csv_path.open()))
    assert len(csv_rows) == len(summary["benchmarks"])


def test_p104b_generated_world_interior_elevation_gate_pass(tmp_path):
    summary = run_suite("P104B", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p104b.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["high_flat_land_floor_gate"]
    assert summary["acceptance"]["mosaic_preservation_gate"]
    assert summary["acceptance"]["high_flat_interiors_below_direct_threshold"]
    assert summary["acceptance"]["continental_land_floor_preserved"]
    assert summary["acceptance"]["p104a_mosaic_preserved"]

    rows = summary["generated_world_interior_elevation_rows"]
    assert {row["seed"] for row in rows} == {42, 31415, 16180}
    assert all(row["cells"] == 900 for row in rows)
    assert all(row["p104b_interior_elevation_pass"] for row in rows)
    assert min(row["continental_land_fraction_world"] for row in rows) >= 0.25
    assert max(
        row["high_flat_interior_fraction_of_continental_land"] for row in rows
    ) <= 0.020
    assert min(row["p104a_min_province_class_count_per_major"] for row in rows) >= 3
    assert max(row["p104a_failing_major_component_count"] for row in rows) == 0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p104b_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    csv_rows = list(csv.DictReader(csv_path.open()))
    assert len(csv_rows) == len(summary["benchmarks"])


def test_p104c_internal_geographic_block_initialization_gate_pass(tmp_path):
    summary = run_suite("P104C", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p104c.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["internal_geographic_block_gate"]
    assert summary["acceptance"]["mosaic_and_highflat_preservation_gate"]
    assert summary["acceptance"]["late_platform_fallback_dependency_reduced"]
    assert summary["acceptance"]["p104a_mosaic_preserved"]
    assert summary["acceptance"]["p104b_highflat_preserved"]

    rows = summary["generated_world_internal_block_rows"]
    assert {row["seed"] for row in rows} == {42, 31415, 16180}
    assert all(row["cells"] == 900 for row in rows)
    assert all(row["p104c_internal_block_pass"] for row in rows)
    assert min(
        row["min_internal_block_class_count_per_major_continent"]
        for row in rows
    ) >= 4
    assert max(
        row["max_largest_internal_block_class_fraction"] for row in rows
    ) <= 0.50
    assert max(row["p104a_platform_split_area_fraction"] for row in rows) <= 0.010
    assert max(row["p104a_failing_major_component_count"] for row in rows) == 0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p104c_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    csv_rows = list(csv.DictReader(csv_path.open()))
    assert len(csv_rows) == len(summary["benchmarks"])


def test_p104d_internal_block_region_expression_gate_pass(tmp_path):
    summary = run_suite("P104D", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p104d.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["internal_block_region_expression_gate"]
    assert summary["acceptance"]["p104b_p104c_preservation_gate"]
    assert summary["acceptance"]["region_expression_reduces_local_striping"]
    assert summary["acceptance"]["region_expression_preserves_block_classes"]
    assert summary["acceptance"]["p104b_highflat_preserved"]

    rows = summary["generated_world_internal_block_region_rows"]
    assert {row["seed"] for row in rows} == {42, 31415, 16180}
    assert all(row["cells"] == 900 for row in rows)
    assert all(row["p104d_region_expression_pass"] for row in rows)
    assert all(row["p104d_preservation_pass"] for row in rows)
    assert min(row["p104d_region_class_count"] for row in rows) >= 4
    assert max(
        row["p104d_region_max_largest_block_class_fraction"] for row in rows
    ) <= 0.62
    assert min(row["p104d_same_neighbor_improvement"] for row in rows) >= 0.015
    assert max(row["p104a_failing_major_component_count"] for row in rows) == 0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p104d_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    csv_rows = list(csv.DictReader(csv_path.open()))
    assert len(csv_rows) == len(summary["benchmarks"])


def test_p104e_continental_detail_region_expression_gate_pass(tmp_path):
    summary = run_suite("P104E", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p104e.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["detail_region_expression_gate"]
    assert summary["acceptance"]["compiler_and_preservation_gate"]
    assert summary["acceptance"]["compiler_consumes_detail_region"]
    assert summary["acceptance"]["detail_region_reduces_local_striping"]
    assert summary["acceptance"]["detail_region_preserves_classes"]

    rows = summary["generated_world_detail_region_rows"]
    assert {row["seed"] for row in rows} == {42, 31415, 16180}
    assert all(row["cells"] == 900 for row in rows)
    assert all(row["p104e_detail_region_expression_pass"] for row in rows)
    assert all(row["p104e_compiler_expression_pass"] for row in rows)
    assert all(row["p104e_preservation_pass"] for row in rows)
    assert {
        row["compiler_detail_field"] for row in rows
    } == {"terrain.continental_detail_region_code"}
    assert min(row["p104e_detail_region_class_count"] for row in rows) >= 5
    assert max(
        row["p104e_detail_region_max_largest_class_fraction"] for row in rows
    ) <= 0.66
    assert min(row["p104e_same_neighbor_improvement"] for row in rows) >= 0.060

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p104e_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    csv_rows = list(csv.DictReader(csv_path.open()))
    assert len(csv_rows) == len(summary["benchmarks"])


def test_p104f_inland_elevation_region_response_gate_pass(tmp_path):
    summary = run_suite("P104F", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p104f.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["inland_elevation_region_gate"]
    assert summary["acceptance"]["p104b_p104e_preservation_gate"]
    assert summary["acceptance"]["inland_region_field_available"]
    assert summary["acceptance"]["land_mask_preserved"]
    assert summary["acceptance"]["inland_elevation_bands_preserved_or_improved"]
    assert summary["acceptance"]["p104e_compiler_detail_preserved"]

    rows = summary["generated_world_inland_elevation_region_rows"]
    assert {row["seed"] for row in rows} == {42, 31415, 16180}
    assert all(row["cells"] == 900 for row in rows)
    assert all(row["p104f_inland_elevation_region_pass"] for row in rows)
    assert all(row["p104f_preservation_pass"] for row in rows)
    assert min(row["p104f_inland_region_class_count"] for row in rows) >= 3
    assert min(row["p104f_applied_continental_inland_fraction"] for row in rows) >= 0.20
    assert max(
        row["p104f_max_major_largest_elevation_band_fraction_after"]
        for row in rows
    ) <= 0.72
    assert max(
        row["high_flat_interior_fraction_of_continental_land"] for row in rows
    ) <= 0.020

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p104f_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    csv_rows = list(csv.DictReader(csv_path.open()))
    assert len(csv_rows) == len(summary["benchmarks"])


def test_p30_inland_state_regional_surface_microbenchmark_pass(tmp_path):
    summary = run_suite("P30", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.p30.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["inland_state_candidate_regional_relief"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {"P30.inland_state_candidate_regional_relief"}
    metrics = benchmarks["P30.inland_state_candidate_regional_relief"]["metrics"]
    assert metrics["shield_state_fraction"] > 0.80
    assert metrics["sag_state_fraction"] > 0.55
    assert metrics["old_orogen_state_fraction"] > 0.55
    assert metrics["rift_state_fraction"] > 0.55
    assert metrics["swell_state_fraction"] > 0.55
    assert metrics["plateau_state_fraction"] > 0.55
    assert metrics["land_mask_unchanged"]
    assert metrics["after_inland_p90_p10_m"] > metrics["before_inland_p90_p10_m"] + 90.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "p30_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_e1_real_earth_continental_landform_microbenchmarks_pass(tmp_path):
    summary = run_suite("E1", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.e1.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["continental_landform_objects"]
    assert summary["acceptance"]["real_earth_interior_landform_coverage"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "E1.craton_platform_basin",
        "E1.old_orogen_decay",
        "E1.collision_plateau",
        "E1.foreland_basin",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    craton = benchmarks["E1.craton_platform_basin"]["metrics"]
    assert craton["shield_objects"] >= 1
    assert craton["platform_objects"] >= 1
    assert craton["interior_basin_objects"] >= 1
    assert craton["basin_mean_elevation_m"] < craton["platform_mean_elevation_m"]
    assert craton["basin_mean_sediment_m"] > craton["platform_mean_sediment_m"]

    old = benchmarks["E1.old_orogen_decay"]["metrics"]
    assert old["recent_orogen_objects"] >= 1
    assert old["old_subdued_orogen_objects"] >= 1
    assert old["recent_orogen_uplift_median_m"] > old["old_orogen_uplift_median_m"]

    plateau = benchmarks["E1.collision_plateau"]["metrics"]
    assert plateau["plateau_objects"] >= 1
    assert plateau["plateau_core_detail_fraction"] >= 0.55
    assert plateau["plateau_mean_elevation_m"] > plateau["platform_mean_elevation_m"]
    assert plateau["plateau_parent_object_count"] >= 1

    foreland = benchmarks["E1.foreland_basin"]["metrics"]
    assert foreland["foreland_basin_objects"] >= 1
    assert foreland["orogen_objects"] >= 1
    assert foreland["foreland_mean_elevation_m"] < foreland["orogen_mean_elevation_m"]
    assert foreland["foreland_mean_sediment_m"] > foreland["platform_mean_sediment_m"]

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "e1_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_e2_real_earth_margin_landform_microbenchmarks_pass(tmp_path):
    summary = run_suite("E2", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.e2.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["margin_landform_objects"]
    assert summary["acceptance"]["real_earth_margin_landform_coverage"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "E2.passive_margin_shelf_wedge",
        "E2.active_margin_trench_arc",
        "E2.delta_fan",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    passive = benchmarks["E2.passive_margin_shelf_wedge"]["metrics"]
    assert passive["passive_margin_wedge_objects"] >= 1
    assert passive["shelf_depth_p75_m"] < passive["slope_depth_p50_m"]
    assert passive["slope_depth_p50_m"] < passive["rise_depth_p50_m"]
    assert passive["wedge_shelf_mean_sediment_m"] > passive["wedge_rise_mean_sediment_m"]

    active = benchmarks["E2.active_margin_trench_arc"]["metrics"]
    assert active["trench_objects"] >= 1
    assert active["forearc_accretionary_prism_objects"] >= 1
    assert active["volcanic_arc_objects"] >= 1
    assert active["trench_mean_depth_m"] > 4200.0
    assert active["arc_mean_elevation_m"] > active["forearc_mean_elevation_m"]
    assert active["trench_centroid_lon"] < active["arc_centroid_lon"]

    delta = benchmarks["E2.delta_fan"]["metrics"]
    assert delta["delta_fan_objects"] >= 1
    assert delta["delta_mean_sediment_m"] > 1800.0
    assert delta["delta_parent_object_count"] > 0
    assert delta["delta_centroid_distance_from_mouth_deg"] < 10.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "e2_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_e3_real_earth_ocean_basin_fabric_microbenchmarks_pass(tmp_path):
    summary = run_suite("E3", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.e3.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["ocean_fabric_objects"]
    assert summary["acceptance"]["real_earth_ocean_basin_fabric_coverage"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "E3.ridge_transform_fracture_zone",
        "E3.abyssal_plain_sedimentation",
        "E3.ocean_age_isochrons",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    fabric = benchmarks["E3.ridge_transform_fracture_zone"]["metrics"]
    assert fabric["spreading_center_objects"] >= 1
    assert fabric["transform_fault_objects"] >= 1
    assert fabric["fracture_zone_objects"] >= 1
    assert fabric["ridge_mean_age_myr"] < fabric["fracture_mean_age_myr"]
    assert fabric["fracture_combined_lon_span_deg"] > 60.0

    abyss = benchmarks["E3.abyssal_plain_sedimentation"]["metrics"]
    assert abyss["abyssal_plain_objects"] >= 1
    assert abyss["abyssal_mean_age_myr"] > 70.0
    assert abyss["abyssal_mean_sediment_m"] > 1200.0
    assert abyss["abyssal_depth_relief_p90_minus_p10_m"] < 350.0

    isochron = benchmarks["E3.ocean_age_isochrons"]["metrics"]
    assert isochron["age_isochron_band_count"] >= 5
    assert isochron["monotonic_age"]
    assert isochron["monotonic_ridge_distance"]
    assert isochron["west_east_median_age_delta_myr"] < 8.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "e3_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_e4_real_earth_arc_plume_microbenchmarks_pass(tmp_path):
    summary = run_suite("E4", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.e4.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["arc_plume_landform_objects"]
    assert summary["acceptance"]["real_earth_arc_plume_coverage"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "E4.island_arc_accretion",
        "E4.back_arc_basin",
        "E4.hotspot_track",
        "E4.large_igneous_province",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    arc = benchmarks["E4.island_arc_accretion"]["metrics"]
    assert arc["island_arc_objects"] >= 1
    assert arc["accreted_terrane_objects"] >= 1
    assert arc["island_arc_mean_elevation_m"] > 450.0
    assert 0 in arc["terrane_parent_continent_ids"]

    backarc = benchmarks["E4.back_arc_basin"]["metrics"]
    assert backarc["back_arc_basin_objects"] >= 1
    assert backarc["back_arc_mean_age_myr"] < 35.0
    assert backarc["back_arc_centroid_lon"] > backarc["arc_centroid_lon"]

    hotspot = benchmarks["E4.hotspot_track"]["metrics"]
    assert hotspot["hotspot_track_objects"] >= 1
    assert hotspot["distance_age_correlation"] > 0.78
    assert hotspot["oldest_distance_deg"] > hotspot["youngest_distance_deg"]

    lip = benchmarks["E4.large_igneous_province"]["metrics"]
    assert lip["large_igneous_province_objects"] >= 1
    assert lip["lip_mean_thickness_m"] > 42000.0
    assert lip["lip_parent_object_count"] >= 2

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "e4_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)


def test_e5_real_earth_cryosphere_surface_process_microbenchmarks_pass(tmp_path):
    summary = run_suite("E5", tmp_path)

    assert summary["schema"] == "aevum.tectonics_bench.e5.v1"
    assert summary["status"] == "pass"
    assert summary["determinism"]["passed"]
    assert summary["acceptance"]["all_microbenchmarks_pass"]
    assert summary["acceptance"]["cryosphere_landform_objects"]
    assert summary["acceptance"]["real_earth_cryosphere_landform_coverage"]

    benchmarks = {bench["name"]: bench for bench in summary["benchmarks"]}
    assert set(benchmarks) == {
        "E5.ice_sheet_loading",
        "E5.glacial_erosion",
        "E5.postglacial_rebound",
    }
    assert all(bench["passed"] for bench in benchmarks.values())

    loading = benchmarks["E5.ice_sheet_loading"]["metrics"]
    assert loading["ice_sheet_loading_objects"] >= 1
    assert loading["loading_mean_ice_thickness_m"] > 1200.0
    assert loading["estimated_bed_depression_m"] > 330.0

    glacial = benchmarks["E5.glacial_erosion"]["metrics"]
    assert glacial["glacial_erosion_objects"] >= 1
    assert glacial["glacial_mean_erosion_m"] > 800.0
    assert glacial["glacial_relief_p90_minus_p10_m"] > 350.0

    rebound = benchmarks["E5.postglacial_rebound"]["metrics"]
    assert rebound["postglacial_rebound_objects"] >= 1
    assert rebound["rebound_mean_unloaded_ice_m"] > 1100.0
    assert rebound["estimated_rebound_potential_m"] > 330.0

    summary_path = tmp_path / "tectonics_bench_summary.json"
    csv_path = tmp_path / "e5_microbenchmarks.csv"
    assert summary_path.exists()
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == len(benchmarks)
