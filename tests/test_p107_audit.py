import json
import math
from types import SimpleNamespace

import numpy as np

from aevum.core.events import Event
from aevum.core.grid import SphereGrid
from aevum.core.state import WorldState
from aevum.core.units import CONSTANTS
from aevum.diagnostics.p107_audit import (
    _ocean_feature_metrics,
    _p110b_historical_supercontinent_trajectory_metrics,
    P107_PLATE_TERRAIN_MODULES,
    P107AuditConfig,
    P107AuditRun,
    plate_rank_metrics,
    run_p107_audit,
    terminal_audit_metrics,
    write_terminal_audit,
)
from aevum.modules.tectonics import (
    CONT,
    CONT_THICK,
    DOMAIN_ACCRETED_TERRANE,
    INTERNAL_BLOCK_NONE,
    INTERNAL_BLOCK_INTRACRATONIC_BASIN,
    INTERNAL_BLOCK_RIFTED_MARGIN,
    INTERNAL_BLOCK_STABLE_PLATFORM,
    OCEAN,
    OCEAN_THICK,
    ORIGIN_ARC,
    ORIGIN_PRIMORDIAL,
    ORIGIN_RIDGE,
    PROVINCE_PLATFORM,
    TectonicsModule,
)
from aevum.modules.tectonics import (
    DOMAIN_CONTINENTAL_MARGIN,
    DOMAIN_CONTINENTAL_INTERIOR,
    DOMAIN_CRATON,
    DOMAIN_LIP,
    DOMAIN_OCEANIC,
    DOMAIN_SUTURE,
    INTERNAL_BLOCK_MOBILE_BELT,
    ORIGIN_PLUME_IMPACT,
)
from aevum.modules.terrain import (
    CONT_PROVINCE_ACTIVE_OROGEN,
    CONT_PROVINCE_INTRACRATONIC_BASIN,
    CONT_PROVINCE_RIFT_SYSTEM,
    CONT_PROVINCE_PLATFORM,
    CONT_DETAIL_BASIN,
    CONT_DETAIL_PLATFORM,
    CONT_DETAIL_PLATEAU,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_RIFT_BASIN,
    CONT_DETAIL_SHIELD,
    INLAND_PROVINCE_PLATFORM,
    INLAND_PROVINCE_RIFT,
    INLAND_PROVINCE_SAG_BASIN,
    OCEAN_DEPTH_ABYSS,
    OCEAN_DEPTH_RIDGE,
    OCEAN_DEPTH_TRENCH,
    TerrainModule,
)
from aevum.spec.presets import get_preset


class _DummyGrid:
    def __init__(self, n):
        self.n = int(n)
        self.cell_area = np.ones(self.n, dtype=np.float64)


def test_p107_plate_rank_metrics_separates_major_minor_micro_and_tiny():
    plate = np.full(10000, 0, dtype=np.int64)
    plate[:600] = 1
    plate[600:700] = 2
    plate[700:710] = 3
    plate[710:712] = 4

    metrics = plate_rank_metrics(_DummyGrid(10000), plate)

    assert metrics["terminal_active_plate_count"] == 5
    assert metrics["plate_rank_counts"]["major"] == 2
    assert metrics["plate_rank_counts"]["minor"] == 1
    assert metrics["plate_rank_counts"]["micro"] == 1
    assert metrics["plate_rank_counts"]["tiny"] == 1
    assert metrics["per_plate"][0]["rank"] == "major"


def test_p107_terminal_audit_writes_raw_arrays_and_required_metrics(tmp_path):
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    plate = np.arange(grid.n, dtype=np.int64) % 6
    elev = np.full(grid.n, -3800.0, dtype=np.float64)
    land = grid.lat > 18.0
    elev[land] = 450.0 + 12.0 * grid.lat[land]
    ridge = np.where((np.abs(grid.lat) < 5.0) & (np.abs(grid.lon) < 70.0))[0]
    transform = np.where((np.abs(grid.lat - 10.0) < 5.0) & (np.abs(grid.lon) < 60.0))[0]
    trench = np.where((np.abs(grid.lat + 25.0) < 5.0) & (np.abs(grid.lon - 80.0) < 70.0))[0]
    elev[trench] = -6500.0

    world.fields.update({
        "tectonics.plate_id": plate.astype(np.float64),
        "terrain.elevation_m": elev,
        "crust.type": land.astype(np.float64),
        "crust.age_myr": np.where(land, 1500.0, 80.0),
        "ocean.depth_province": np.where(elev < 0.0, 4.0, 0.0),
        "ocean.basin_id": np.where(elev < 0.0, (grid.lon > 0.0).astype(float), -1.0),
        "terrain.province": np.where(land, 2.0, 0.0),
        "terrain.continental_detail": np.where(land, 1.0, 0.0),
        "terrain.mountain_ranges": np.where(land, 1.0, 0.0),
        "terrain.mountain_inventory": np.where(land, 1.0, 0.0),
        "terrain.mountain_hierarchy_level": np.where(land, 1.0, 0.0),
        "tectonics.mountain_belt_id": np.where(land, 1.0, 0.0),
        "tectonics.mountain_parent_process_id": np.where(land, 5.0, 0.0),
    })
    world.networks["tectonics.boundaries"] = {
        "ridge": ridge.astype(np.int64),
        "transform": transform.astype(np.int64),
        "trench": trench.astype(np.int64),
        "suture": np.where((grid.lat > 35.0) & (np.abs(grid.lon) < 80.0))[0].astype(np.int64),
    }
    world.objects["tectonics.boundary_objects"] = [
        {"kind": "ridge", "cells": ridge.astype(int).tolist()},
        {"kind": "transform", "cells": transform.astype(int).tolist()},
        {"kind": "trench", "cells": trench.astype(int).tolist()},
    ]
    world.objects["terrain.ocean_fabric"] = [
        {"kind": "spreading_center", "cells": ridge.astype(int).tolist()},
        {"kind": "transform_fault", "cells": transform.astype(int).tolist()},
    ]
    world.objects["terrain.arc_plume_landforms"] = [
        {
            "kind": "seamount_chain",
            "cells": transform.astype(int).tolist(),
            "ocean_fraction": 1.0,
            "parent_tectonic_object_ids": ["plume:test"],
        }
    ]
    world.objects["terrain.mountain_ranges"] = [
        {"kind": "orogen", "cells": np.where(land)[0].astype(int).tolist()}
    ]
    world.objects["terrain.internal_profile"] = {
        "schema": "aevum.terrain_internal_profile.v1",
        "enabled": True,
        "time_myr": 4500.0,
        "cell_count": int(grid.n),
        "stage_seconds": {
            "initial_isostasy_hydrology": 0.12,
            "semantic_object_build": 0.34,
            "final_margin_drainage": 0.21,
        },
        "subprofiles": {
            "semantic_object_build": {
                "stage_seconds": {
                    "continental_classes_landforms": 0.08,
                    "initial_ocean_geography": 0.17,
                    "landform_inventory_cleanup": 0.09,
                },
            },
        },
        "total_seconds": 0.75,
    }
    world.set_g("terrain.last_quantile_sea_level_m", 0.0)
    world.set_g("terrain.last_water_budget_sea_level_m", -120.0)
    world.set_g("terrain.last_sea_level_water_budget_delta_m", 120.0)
    world.set_g("terrain.last_p109_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p109_land_mean_before_m", 2200.0)
    world.set_g("terrain.last_p109_land_mean_after_m", 850.0)
    world.set_g("terrain.last_p110b_final_state_archetype_code", 1.0)
    world.set_g("terrain.last_p110b_largest_share_preferred_ceiling", 0.54)
    world.set_g("terrain.last_p110b_largest_share_soft_ceiling", 0.59)
    world.set_g("terrain.last_p110b_min_nonlargest_large_components", 2.0)
    world.set_g("terrain.last_p111619_late_domain_partition_attempted", 1.0)
    world.set_g("terrain.last_p111619_late_domain_partition_applied", 1.0)
    world.set_g("terrain.last_p111619_late_domain_partition_candidate_count", 2.0)
    world.set_g("terrain.last_p111619_late_domain_partition_area_fraction", 0.012)
    world.set_g("terrain.last_p111619_late_domain_partition_largest_before", 0.48)
    world.set_g("terrain.last_p111619_late_domain_partition_largest_after", 0.38)
    world.set_g("terrain.last_p111619_late_domain_partition_major_count_before", 4.0)
    world.set_g("terrain.last_p111619_late_domain_partition_major_count_after", 5.0)
    world.set_g("terrain.last_p111619_late_domain_partition_domain_count", 10.0)
    world.set_g("terrain.last_p111619_late_domain_partition_domain_floor", 0.06)
    world.set_g("terrain.last_p111619_late_domain_partition_land_fraction_before", 0.278)
    world.set_g("terrain.last_p111619_late_domain_partition_land_fraction_after", 0.266)
    world.set_g("terrain.last_p154_final_planform_gate_attempted", 1.0)
    world.set_g("terrain.last_p154_final_planform_gate_applied", 1.0)
    world.set_g("terrain.last_p154_final_planform_gate_reverted", 0.0)
    world.set_g("terrain.last_p154_largest_share_before", 0.974)
    world.set_g("terrain.last_p154_largest_share_after", 0.514)
    world.set_g("terrain.last_p154_component_count_before", 11.0)
    world.set_g("terrain.last_p154_component_count_after", 5.0)
    world.set_g("terrain.last_p154_major_component_count_before", 1.0)
    world.set_g("terrain.last_p154_major_component_count_after", 4.0)
    world.set_g("terrain.last_p154_land_fraction_before", 0.295)
    world.set_g("terrain.last_p154_land_fraction_after", 0.267)
    world.set_g("terrain.last_p154_opened_area_fraction", 0.028)
    world.set_g("terrain.last_p154_reject_code", 0.0)
    world.set_g(
        "terrain.last_p155_terminal_high_relief_consistency_gate_accepted",
        1.0,
    )
    world.set_g("terrain.last_p155_guard_reverted", 0.0)
    world.set_g("terrain.last_p155_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p155_candidate_cell_count", 9.0)
    world.set_g("terrain.last_p155_candidate_area_fraction", 0.014)
    world.set_g("terrain.last_p155_raised_cell_count", 5.0)
    world.set_g("terrain.last_p155_raised_area_fraction", 0.008)
    world.set_g("terrain.last_p155_softened_cell_count", 4.0)
    world.set_g("terrain.last_p155_softened_area_fraction", 0.006)
    world.set_g("terrain.last_p155_high_component_count_before", 12.0)
    world.set_g("terrain.last_p155_high_component_count_after", 7.0)
    world.set_g("terrain.last_p155_spine_3000_component_count_before", 5.0)
    world.set_g("terrain.last_p155_spine_3000_component_count_after", 3.0)
    world.set_g("terrain.last_p155_spine_3000_coverage_before", 0.36)
    world.set_g("terrain.last_p155_spine_3000_coverage_after", 0.49)
    world.set_g("terrain.last_p155_top3_high_share_before", 0.58)
    world.set_g("terrain.last_p155_top3_high_share_after", 0.76)
    world.set_g("terrain.last_p155_fragmentation_index_before", 0.68)
    world.set_g("terrain.last_p155_fragmentation_index_after", 0.41)
    world.set_g("terrain.last_p155_parent_overlap_before", 0.52)
    world.set_g("terrain.last_p155_parent_overlap_after", 0.71)
    world.set_g("terrain.last_p155_high_area_fraction_before", 0.021)
    world.set_g("terrain.last_p155_high_area_fraction_after", 0.024)
    world.set_g("terrain.last_p155_p90_land_relief_before_m", 2520.0)
    world.set_g("terrain.last_p155_p90_land_relief_after_m", 2540.0)
    world.set_g("terrain.last_p155_p98_land_relief_before_m", 3940.0)
    world.set_g("terrain.last_p155_p98_land_relief_after_m", 4020.0)
    world.set_g("terrain.last_p155_reject_code", 0.0)
    world.set_g(
        "terrain.last_p162_terminal_high_mountain_fragment_cleanup_accepted",
        1.0,
    )
    world.set_g("terrain.last_p162_guard_reverted", 0.0)
    world.set_g("terrain.last_p162_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p162_candidate_cell_count", 12.0)
    world.set_g("terrain.last_p162_candidate_area_fraction", 0.018)
    world.set_g("terrain.last_p162_softened_cell_count", 9.0)
    world.set_g("terrain.last_p162_softened_area_fraction", 0.014)
    world.set_g("terrain.last_p162_high_component_count_before", 28.0)
    world.set_g("terrain.last_p162_high_component_count_after", 12.0)
    world.set_g("terrain.last_p162_small_high_component_count_before", 22.0)
    world.set_g("terrain.last_p162_small_high_component_count_after", 8.0)
    world.set_g("terrain.last_p162_spine_3000_component_count_before", 15.0)
    world.set_g("terrain.last_p162_spine_3000_component_count_after", 12.0)
    world.set_g("terrain.last_p162_top3_high_share_before", 0.29)
    world.set_g("terrain.last_p162_top3_high_share_after", 0.52)
    world.set_g("terrain.last_p162_fragmentation_index_before", 0.95)
    world.set_g("terrain.last_p162_fragmentation_index_after", 0.55)
    world.set_g("terrain.last_p162_high_area_fraction_before", 0.019)
    world.set_g("terrain.last_p162_high_area_fraction_after", 0.014)
    world.set_g("terrain.last_p162_p90_land_relief_before_m", 2480.0)
    world.set_g("terrain.last_p162_p90_land_relief_after_m", 2440.0)
    world.set_g("terrain.last_p162_p98_land_relief_before_m", 3900.0)
    world.set_g("terrain.last_p162_p98_land_relief_after_m", 3760.0)
    world.set_g("terrain.last_p162_mean_lower_m", 320.0)
    world.set_g("terrain.last_p162_max_lower_m", 610.0)
    world.set_g("terrain.last_p162_extreme_softened_cell_count", 0.0)
    world.set_g("terrain.last_p162_reject_code", 0.0)
    world.set_g("terrain.last_p155_p156_path_candidate_count", 3.0)
    world.set_g("terrain.last_p155_p156_path_candidate_cell_count", 11.0)
    world.set_g("terrain.last_p155_p156_selected_path_count", 2.0)
    world.set_g("terrain.last_p155_p156_selected_path_cell_count", 7.0)
    world.set_g("terrain.last_p155_p157_seed_path_candidate_count", 4.0)
    world.set_g("terrain.last_p155_p157_seed_path_candidate_cell_count", 13.0)
    world.set_g("terrain.last_p155_p157_selected_seed_path_count", 2.0)
    world.set_g("terrain.last_p155_p157_selected_seed_path_cell_count", 8.0)
    event = Event(
        "plate_reorganization",
        4500.0,
        "tectonics",
        params={"merged": [{"from": 7, "to": 1, "area_fraction": 0.001}]},
    )

    metrics = terminal_audit_metrics(world, [event])
    assert metrics["schema"] == "aevum.p107_terminal_world_audit.v1"
    assert metrics["merged_microplate_count"] == 1
    assert metrics["ridge_transform_network_continuity"]["has_ridge_transform_network"]
    assert metrics["parented_oceanic_island_chain_count"] == 1
    assert "p1115_unsupported_open_ocean_shoal_fraction" in metrics
    assert "p1115_object_backed_ocean_relief_fraction" in metrics
    assert "p1115_ridge_area_fraction_ocean" in metrics["ocean_feature_metrics"]
    assert "land_elevation_band_fraction" in metrics
    assert "p109_hypsometry_comparison" in metrics
    assert "p109_sea_level_consistency" in metrics
    assert "p110a_modern_planform" in metrics
    assert "p1134_lowres_false_oceanic_island_cleanup" in metrics
    assert (
        metrics["p1134_lowres_false_oceanic_island_cleanup"]["terminal_land_fraction"]
        == metrics["land_fraction"]
    )
    assert "p1134b_lowres_articulation_neck_widening" in metrics
    assert (
        metrics["p1134b_lowres_articulation_neck_widening"]["terminal_land_fraction"]
        == metrics["land_fraction"]
    )
    assert "p1134c_lowres_shelf_depth_cap" in metrics
    assert (
        metrics["p1134c_lowres_shelf_depth_cap"]["terminal_land_fraction"]
        == metrics["land_fraction"]
    )
    assert (
        metrics["p110a_modern_planform"]["p110b_final_state_archetype"]["name"]
        == "modern_multipolar"
    )
    assert "p110b_lineage_survival" in metrics["p110a_modern_planform"]
    assert "province_hypsometry" in metrics
    assert metrics["terrain_internal_profile"]["available"]
    assert metrics["terrain_internal_profile"]["top"][0]["stage"] == "semantic_object_build"
    semantic_terminal = metrics["terrain_internal_profile"]["subprofiles"][
        "semantic_object_build"
    ]
    assert semantic_terminal["top"][0]["stage"] == "initial_ocean_geography"
    assert metrics["acceptance"]["has_p109_hypsometry_metrics"]
    assert "p1115_unsupported_open_ocean_shoal_fraction" in metrics["acceptance"]
    assert metrics["acceptance"]["has_p109_sea_level_consistency_metrics"]
    assert metrics["acceptance"]["has_p110a_planform_metrics"]
    assert metrics["p109_sea_level_consistency"]["p109_land_mask_preserved"]
    assert metrics["p110a_modern_planform"]["land"]["component_count"] >= 1
    seaway = metrics["p110a_modern_planform"]["seaway_cut_effectiveness"]
    assert seaway["p111619_late_domain_partition_attempted"]
    assert seaway["p111619_late_domain_partition_applied"]
    assert seaway["p111619_late_domain_partition_candidate_count"] == 2
    assert seaway["p111619_late_domain_partition_largest_share_before"] == 0.48
    assert seaway["p111619_late_domain_partition_largest_share_after"] == 0.38
    assert seaway["p111619_late_domain_partition_major_count_before"] == 4
    assert seaway["p111619_late_domain_partition_major_count_after"] == 5
    assert seaway["p111619_late_domain_partition_domain_count"] == 10
    assert seaway["p111619_late_domain_partition_domain_floor"] == 0.06
    assert seaway["p154_final_planform_gate_attempted"]
    assert seaway["p154_final_planform_gate_applied"]
    assert not seaway["p154_final_planform_gate_reverted"]
    assert seaway["p154_largest_share_before"] == 0.974
    assert seaway["p154_largest_share_after"] == 0.514
    assert seaway["p154_component_count_before"] == 11
    assert seaway["p154_component_count_after"] == 5
    assert seaway["p154_major_component_count_before"] == 1
    assert seaway["p154_major_component_count_after"] == 4
    assert seaway["p154_land_fraction_before"] == 0.295
    assert seaway["p154_land_fraction_after"] == 0.267
    assert seaway["p154_opened_area_fraction_world"] == 0.028
    assert seaway["p154_reject_code"] == 0
    mountains = metrics["high_mountain_coherence"]
    assert mountains["p155_terminal_high_relief_consistency_gate_accepted"]
    assert not mountains["p155_guard_reverted"]
    assert mountains["p155_land_mask_preserved"]
    assert mountains["p155_candidate_cell_count"] == 9.0
    assert mountains["p155_candidate_area_fraction"] == 0.014
    assert mountains["p155_raised_cell_count"] == 5.0
    assert mountains["p155_raised_area_fraction"] == 0.008
    assert mountains["p155_softened_cell_count"] == 4.0
    assert mountains["p155_softened_area_fraction"] == 0.006
    assert mountains["p155_high_component_count_before"] == 12.0
    assert mountains["p155_high_component_count_after"] == 7.0
    assert mountains["p155_spine_3000_component_count_before"] == 5.0
    assert mountains["p155_spine_3000_component_count_after"] == 3.0
    assert mountains["p155_spine_3000_coverage_before"] == 0.36
    assert mountains["p155_spine_3000_coverage_after"] == 0.49
    assert mountains["p155_top3_high_share_before"] == 0.58
    assert mountains["p155_top3_high_share_after"] == 0.76
    assert mountains["p155_fragmentation_index_before"] == 0.68
    assert mountains["p155_fragmentation_index_after"] == 0.41
    assert mountains["p155_parent_overlap_before"] == 0.52
    assert mountains["p155_parent_overlap_after"] == 0.71
    assert mountains["p155_high_area_fraction_before"] == 0.021
    assert mountains["p155_high_area_fraction_after"] == 0.024
    assert mountains["p155_p90_land_relief_before_m"] == 2520.0
    assert mountains["p155_p90_land_relief_after_m"] == 2540.0
    assert mountains["p155_p98_land_relief_before_m"] == 3940.0
    assert mountains["p155_p98_land_relief_after_m"] == 4020.0
    assert mountains["p155_reject_code"] == 0.0
    assert mountains["p162_terminal_high_mountain_fragment_cleanup_accepted"]
    assert not mountains["p162_guard_reverted"]
    assert mountains["p162_land_mask_preserved"]
    assert mountains["p162_candidate_cell_count"] == 12.0
    assert mountains["p162_candidate_area_fraction"] == 0.018
    assert mountains["p162_softened_cell_count"] == 9.0
    assert mountains["p162_softened_area_fraction"] == 0.014
    assert mountains["p162_high_component_count_before"] == 28.0
    assert mountains["p162_high_component_count_after"] == 12.0
    assert mountains["p162_small_high_component_count_before"] == 22.0
    assert mountains["p162_small_high_component_count_after"] == 8.0
    assert mountains["p162_spine_3000_component_count_before"] == 15.0
    assert mountains["p162_spine_3000_component_count_after"] == 12.0
    assert mountains["p162_top3_high_share_before"] == 0.29
    assert mountains["p162_top3_high_share_after"] == 0.52
    assert mountains["p162_fragmentation_index_before"] == 0.95
    assert mountains["p162_fragmentation_index_after"] == 0.55
    assert mountains["p162_high_area_fraction_before"] == 0.019
    assert mountains["p162_high_area_fraction_after"] == 0.014
    assert mountains["p162_p90_land_relief_before_m"] == 2480.0
    assert mountains["p162_p90_land_relief_after_m"] == 2440.0
    assert mountains["p162_p98_land_relief_before_m"] == 3900.0
    assert mountains["p162_p98_land_relief_after_m"] == 3760.0
    assert mountains["p162_mean_lower_m"] == 320.0
    assert mountains["p162_max_lower_m"] == 610.0
    assert mountains["p162_extreme_softened_cell_count"] == 0.0
    assert mountains["p162_reject_code"] == 0.0
    assert mountains["p156_path_candidate_count"] == 3.0
    assert mountains["p156_path_candidate_cell_count"] == 11.0
    assert mountains["p156_selected_path_count"] == 2.0
    assert mountains["p156_selected_path_cell_count"] == 7.0
    assert mountains["p157_seed_path_candidate_count"] == 4.0
    assert mountains["p157_seed_path_candidate_cell_count"] == 13.0
    assert mountains["p157_selected_seed_path_count"] == 2.0
    assert mountains["p157_selected_seed_path_cell_count"] == 8.0
    assert "terrain_province" in metrics["province_hypsometry"]

    written = write_terminal_audit(
        world,
        tmp_path,
        events=[event],
        render_assets=False,
        include_earth_reference=False,
        contact_sheet=False,
    )
    assert written["array_archive"]["key_count"] >= 10
    assert (tmp_path / "p107_terminal_arrays.npz").exists()
    assert (tmp_path / "p107_terminal_metrics.json").exists()
    saved = json.loads((tmp_path / "p107_terminal_metrics.json").read_text())
    assert saved["acceptance"]["p107_0_observability_complete"]
    with np.load(tmp_path / "p107_terminal_arrays.npz") as data:
        assert "field__tectonics_plate_id" in data.files
        assert "field__terrain_elevation_m" in data.files
        assert "field__ocean_basin_id" in data.files
        assert "field__ocean_gateway_id" in data.files
        assert "field__ocean_gateway_system_id" in data.files
        assert "field__ocean_margin_type" in data.files
        assert "field__ocean_shelf_width" in data.files
        assert "field__terrain_mountain_ranges" in data.files
        assert "field__terrain_orogenic_shoulder_halo" in data.files
        assert "field__terrain_orogenic_highland_apron" in data.files
        assert "field__tectonics_mountain_belt_id" in data.files
        assert "boundary__ridge" in data.files
        assert "object__spreading_center" in data.files
        assert "object__mountain_orogen" in data.files


def test_p107_audit_summary_records_stage_timings(tmp_path):
    config = P107AuditConfig(
        runs=(P107AuditRun(cells=90, n_plates=8, seed=1707, label="timing_smoke"),),
        t_end_myr=40.0,
        frames=1,
        render_world_assets=False,
        render_contact_sheet=False,
        include_earth_reference=False,
    )

    summary = run_p107_audit(config, tmp_path)
    entry = summary["entries"][0]
    stages = entry["stage_seconds"]

    assert entry["run_seconds"] == stages["engine_run"]
    assert stages["build"] >= 0.0
    assert stages["configure_world"] >= 0.0
    assert stages["terminal_audit_write"] >= 0.0
    assert entry["total_seconds"] >= entry["run_seconds"]
    assert entry["module_seconds"]["totals"]["terrain"] >= 0.0
    assert entry["module_seconds"]["top"][0]["seconds"] >= 0.0
    aggregate_profile = entry["terrain_internal_profile"]
    assert aggregate_profile["available"]
    assert aggregate_profile["profiled_step_count"] >= 1
    assert aggregate_profile["terrain_run_count"] >= aggregate_profile["profiled_step_count"]
    assert aggregate_profile["profiled_stage_seconds_total"] >= 0.0
    assert aggregate_profile["profile_coverage_fraction_of_terrain_module_seconds"] > 0.0
    assert aggregate_profile["top"]
    semantic_aggregate = aggregate_profile["subprofiles"]["semantic_object_build"]
    assert semantic_aggregate["stage_count"] >= 1
    assert semantic_aggregate["top"]
    terminal_bathymetry_aggregate = aggregate_profile["subprofiles"][
        "terminal_bathymetry_polish"
    ]
    assert terminal_bathymetry_aggregate["stage_count"] >= 1
    assert terminal_bathymetry_aggregate["top"]
    final_margin_aggregate = aggregate_profile["subprofiles"]["final_margin_drainage"]
    assert final_margin_aggregate["stage_count"] >= 1
    assert final_margin_aggregate["top"]
    source_to_sink_aggregate = aggregate_profile["subprofiles"]["source_to_sink"]
    assert source_to_sink_aggregate["stage_count"] >= 1
    assert source_to_sink_aggregate["top"]
    profile = entry["metrics"]["terrain_internal_profile"]
    assert profile["available"]
    assert profile["stage_count"] >= 1
    assert profile["top"]
    assert "semantic_object_build" in profile["subprofiles"]
    assert "terminal_bathymetry_polish" in profile["subprofiles"]
    assert profile["top"][0]["seconds"] >= 0.0
    saved = json.loads((tmp_path / "p107_audit_summary.json").read_text())
    assert "stage_seconds" in saved["entries"][0]
    assert "module_seconds" in saved["entries"][0]
    assert saved["entries"][0]["terrain_internal_profile"]["available"]
    assert (
        saved["entries"][0]["terrain_internal_profile"]["subprofiles"][
            "semantic_object_build"
        ]["top"]
    )
    assert saved["entries"][0]["metrics"]["terrain_internal_profile"]["available"]
    compact_mountains = saved["entries"][0]["metrics"]["high_mountain_coherence"]
    terminal_metrics = json.loads(
        (tmp_path / saved["entries"][0]["outdir"] / "p107_terminal_metrics.json")
        .read_text()
    )
    terminal_mountains = terminal_metrics["high_mountain_coherence"]
    for key in (
        "p155_land_mask_preserved",
        "p155_candidate_area_fraction",
        "p155_raised_area_fraction",
        "p155_softened_area_fraction",
        "p155_spine_3000_component_count_before",
        "p155_spine_3000_component_count_after",
        "p155_spine_3000_coverage_before",
        "p155_spine_3000_coverage_after",
        "p155_top3_high_share_before",
        "p155_top3_high_share_after",
        "p155_high_area_fraction_before",
        "p155_high_area_fraction_after",
        "p155_p90_land_relief_before_m",
        "p155_p90_land_relief_after_m",
        "p155_p98_land_relief_before_m",
        "p155_p98_land_relief_after_m",
        "p155_reject_code",
    ):
        assert key in compact_mountains
        assert compact_mountains[key] == terminal_mountains[key]


def test_p107_plate_terrain_only_preview_keeps_climate_feedbacks(tmp_path):
    config = P107AuditConfig(
        runs=(P107AuditRun(cells=90, n_plates=8, seed=1708, label="preview"),),
        t_end_myr=40.0,
        frames=1,
        render_world_assets=False,
        render_contact_sheet=False,
        include_earth_reference=False,
        enabled_modules=P107_PLATE_TERRAIN_MODULES,
    )

    summary = run_p107_audit(config, tmp_path)
    entry = summary["entries"][0]
    selection = entry["scheduler_modules"]

    assert summary["config"]["enabled_modules"] == list(P107_PLATE_TERRAIN_MODULES)
    assert selection["mode"] == "allowlist"
    assert selection["enabled"] == list(P107_PLATE_TERRAIN_MODULES)
    assert "resources" in selection["disabled"]
    assert "climate" in selection["enabled"]
    assert "biogeochem" in selection["enabled"]
    assert "biosphere" in selection["enabled"]
    assert set(entry["module_seconds"]["totals"]) == set(P107_PLATE_TERRAIN_MODULES)
    assert entry["terrain_internal_profile"]["available"]


def test_terrain_smooth_field_graph_cache_matches_uncached_formula():
    grid = SphereGrid.fibonacci(96, CONSTANTS.EARTH_RADIUS)
    terrain = TerrainModule()
    values = np.sin(np.radians(grid.lat)) + 0.25 * np.cos(np.radians(grid.lon))

    def reference_smooth(input_values, *, passes: int, alpha: float):
        out = input_values.astype(np.float64).copy()
        edges = grid.edges
        i, j = edges[:, 0], edges[:, 1]
        deg = np.zeros(grid.n, dtype=np.float64)
        np.add.at(deg, i, 1.0)
        np.add.at(deg, j, 1.0)
        deg = np.maximum(deg, 1.0)
        for _ in range(passes):
            acc = np.zeros_like(out)
            np.add.at(acc, i, out[j])
            np.add.at(acc, j, out[i])
            out = (1.0 - alpha) * out + alpha * (acc / deg)
        return out

    smoothed = terrain._smooth_field(grid, values, passes=3, alpha=0.21)
    cached_smoothed = terrain._smooth_field(grid, values + 0.1, passes=2, alpha=0.17)

    np.testing.assert_allclose(
        smoothed,
        reference_smooth(values, passes=3, alpha=0.21),
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        cached_smoothed,
        reference_smooth(values + 0.1, passes=2, alpha=0.17),
        rtol=0.0,
        atol=0.0,
    )
    assert getattr(terrain, "_smooth_field_graph_cache")[0][1] == grid.n


def test_terrain_frontier_expansion_helpers_match_reference_loops():
    grid = SphereGrid.fibonacci(180, CONSTANTS.EARTH_RADIUS)
    terrain = TerrainModule()
    seed_cells = np.asarray([3, 17, 48, 93, 144], dtype=np.int64)
    source = np.zeros(grid.n, dtype=bool)
    source[seed_cells] = True

    def reference_dilated(cells, *, passes: int):
        out = np.zeros(grid.n, dtype=bool)
        out[cells] = True
        frontier = out.copy()
        for _ in range(passes):
            nxt = np.zeros(grid.n, dtype=bool)
            for c in np.where(frontier)[0]:
                nxt[grid.neighbors[c]] = True
            nxt &= ~out
            if not nxt.any():
                break
            out |= nxt
            frontier = nxt
        return out

    def reference_band(mask, *, passes: int):
        influence = np.zeros(grid.n, dtype=np.float64)
        frontier = np.asarray(mask, dtype=bool).copy()
        if not frontier.any():
            return influence
        seen = frontier.copy()
        influence[frontier] = 1.0
        weight = 0.68
        for p in range(1, passes + 1):
            nxt = np.zeros(grid.n, dtype=bool)
            for c in np.where(frontier)[0]:
                nxt[grid.neighbors[c]] = True
            nxt &= ~seen
            if not nxt.any():
                break
            influence[nxt] = weight ** p
            seen |= nxt
            frontier = nxt
        return terrain._smooth_field(grid, influence, passes=1, alpha=0.25)

    assert np.array_equal(
        terrain._dilated_bool(grid, seed_cells, passes=4),
        reference_dilated(seed_cells, passes=4),
    )
    np.testing.assert_allclose(
        terrain._band_influence(grid, source, passes=4),
        reference_band(source, passes=4),
        rtol=0.0,
        atol=0.0,
    )


def test_terrain_ocean_age_contrast_matches_reference_neighbor_loop():
    grid = SphereGrid.fibonacci(192, CONSTANTS.EARTH_RADIUS)
    terrain = TerrainModule()
    age = (
        80.0
        + 35.0 * np.sin(np.radians(grid.lon))
        + 12.0 * np.cos(np.radians(2.0 * grid.lat))
    )
    age[::29] = np.nan
    mask = (grid.lat < 62.0) & (np.abs(grid.lon) > 18.0)

    reference = np.zeros(grid.n, dtype=np.float64)
    for c in np.where(mask)[0]:
        c = int(c)
        nbs = grid.neighbors[c]
        valid = nbs[mask[nbs]]
        if valid.size:
            reference[c] = float(np.max(np.abs(age[c] - age[valid])))

    optimized = terrain._ocean_age_contrast(grid, age, mask)
    np.testing.assert_allclose(optimized, reference, equal_nan=True)


def test_terrain_dominant_basin_parent_matches_reference_loop():
    grid = SphereGrid.fibonacci(192, CONSTANTS.EARTH_RADIUS)
    terrain = TerrainModule()
    basin_id = np.full(grid.n, -1, dtype=int)
    basin_id[::2] = 4
    basin_id[1::5] = 9
    basin_id[7::11] = 2
    basin_id[13::17] = 31
    mask = (
        (grid.lat > -52.0)
        & (grid.lat < 66.0)
        & ((grid.lon < -25.0) | (grid.lon > 18.0))
    )

    def reference_dominant_parent():
        valid = mask & (basin_id >= 0)
        if not valid.any():
            return None
        best_id = None
        best_area = -1.0
        for bid in np.unique(basin_id[valid]):
            basin_area = float(grid.cell_area[valid & (basin_id == bid)].sum())
            if basin_area > best_area:
                best_area = basin_area
                best_id = int(bid)
        return best_id

    assert (
        terrain._dominant_basin_parent(basin_id, mask, grid.cell_area)
        == reference_dominant_parent()
    )
    assert terrain._dominant_basin_parent(basin_id, basin_id < 0, grid.cell_area) is None


def test_terrain_ocean_semantic_context_helper_matches_manual_chain():
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 48.0
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[land] = 650.0
    ridge = (~land) & (np.abs(grid.lon) < 8.0)
    surface[ridge] = -2500.0
    crust_type = land.astype(np.float64)
    crust_domain = np.where(land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC)
    crust_origin = np.zeros(grid.n, dtype=np.float64)
    crust_thick = np.where(land, CONT_THICK, 7000.0).astype(np.float64)
    sediment = np.full(grid.n, 180.0, dtype=np.float64)

    world.fields["crust.type"] = crust_type
    world.fields["crust.domain"] = crust_domain.astype(np.float64)
    world.fields["crust.origin"] = crust_origin
    world.fields["crust.thickness_m"] = crust_thick
    world.fields["crust.age_myr"] = np.where(
        ridge, 4.0, np.where(land, 1800.0, 90.0)).astype(np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.continent_id"] = np.where(land, 0.0, -1.0)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(ridge)[0].astype(int).tolist(),
    }

    fields, objects, diag = terrain._ocean_geography(world, surface, sea_level)
    fabric = terrain._ocean_fabric_objects(
        world, surface, sea_level, crust_type, sediment, fields)
    margins = terrain._margin_landform_objects(
        world, surface, sea_level, crust_type, sediment, fields)
    arcs = terrain._arc_plume_landform_objects(
        world, surface, sea_level, crust_type, crust_domain, crust_origin,
        crust_thick, sediment, fields)

    context_profile = {}
    child_profile = {}
    parent_lookup_profile = {}
    object_summary_profile = {}
    context_cache = {}
    helper = terrain._rebuild_ocean_semantic_context(
        world, surface, sea_level, crust_type, crust_domain, crust_origin,
        crust_thick, sediment, profile_stage_seconds=context_profile,
        child_profile_stage_seconds=child_profile,
        parent_lookup_profile_stage_seconds=parent_lookup_profile,
        object_summary_profile_stage_seconds=object_summary_profile,
        profile_prefix="test", cache=context_cache)

    def normalized(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, dict):
            return {str(k): normalized(v) for k, v in sorted(value.items())}
        if isinstance(value, list):
            return [normalized(v) for v in value]
        if isinstance(value, tuple):
            return tuple(normalized(v) for v in value)
        return value

    helper_fields, helper_objects, helper_diag, helper_fabric, helper_margins, helper_arcs = helper
    for key, value in fields.items():
        np.testing.assert_array_equal(helper_fields[key], value)
    assert normalized(helper_objects) == normalized(objects)
    assert normalized(helper_diag) == normalized(diag)
    assert normalized(helper_fabric) == normalized(fabric)
    assert normalized(helper_margins) == normalized(margins)
    assert normalized(helper_arcs) == normalized(arcs)
    assert set(context_profile) == {
        "test.ocean_geography",
        "test.ocean_fabric",
        "test.margin_landforms",
        "test.arc_plume_landforms",
    }
    assert all(value >= 0.0 for value in context_profile.values())
    assert {
        "test.ocean_fabric.prepare_sources",
        "test.ocean_fabric.component_object_summaries",
        "test.arc_plume.source_influence",
        "test.arc_plume.component_object_summaries",
    }.issubset(set(child_profile))
    assert all(value >= 0.0 for value in child_profile.values())
    assert any(
        key.startswith("test.ocean_fabric.parent_lookup.")
        or key.startswith("test.arc_plume.parent_lookup.")
        for key in parent_lookup_profile
    )
    assert all(value >= 0.0 for value in parent_lookup_profile.values())
    assert {
        "test.ocean_fabric.object_summary.component_extract_sort",
        "test.ocean_fabric.object_summary.object_fields",
    }.issubset(set(object_summary_profile))
    assert all(value >= 0.0 for value in object_summary_profile.values())

    cache_profile = {}
    cache_child_profile = {}
    cache_parent_lookup_profile = {}
    cache_object_summary_profile = {}
    cached_helper = terrain._rebuild_ocean_semantic_context(
        world, surface, sea_level, crust_type, crust_domain, crust_origin,
        crust_thick, sediment, profile_stage_seconds=cache_profile,
        child_profile_stage_seconds=cache_child_profile,
        parent_lookup_profile_stage_seconds=cache_parent_lookup_profile,
        object_summary_profile_stage_seconds=cache_object_summary_profile,
        profile_prefix="cached", cache=context_cache)
    assert normalized(cached_helper) == normalized(helper)
    assert set(cache_profile) == {"cached.context_cache_hit"}
    assert cache_child_profile == {}
    assert cache_parent_lookup_profile == {}
    assert cache_object_summary_profile == {}
    assert all(value >= 0.0 for value in cache_profile.values())


def test_terrain_component_summary_direct_masks_match_full_mask_reference():
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    terrain = TerrainModule()

    comp = np.where((np.abs(grid.lat) < 18.0) & (np.abs(grid.lon) < 42.0))[0]
    assert comp.size >= 2
    comp = comp.astype(np.int64)
    full_mask = np.zeros(grid.n, dtype=bool)
    full_mask[comp] = True

    for passes in (1, 2, 3):
        direct = terrain._dilated_bool(grid, comp, passes=passes)
        reference = terrain._dilate_mask(grid, full_mask, passes=passes)
        np.testing.assert_array_equal(direct, reference)

    flag = (grid.lat > 0.0) | (np.abs(grid.lon) < 12.0)
    comp_area = max(float(grid.cell_area[comp].sum()), 1e-12)
    reference_fraction = float(grid.cell_area[full_mask & flag].sum() / comp_area)
    direct_fraction = float(grid.cell_area[comp[flag[comp]]].sum() / comp_area)
    assert direct_fraction == reference_fraction


def test_terrain_nearest_seed_distance_matches_reference_loop():
    spec = get_preset("earthlike")
    spec.grid_cells = 180
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    terrain = TerrainModule()

    cells = np.asarray([3, 17, 24, 51, 88, 123, 149], dtype=np.int64)
    seed_cells = np.asarray([5, 42, 77, 121], dtype=np.int64)

    reference_distance = np.full(cells.size, np.nan, dtype=np.float64)
    reference_seed = None
    reference_seed_distance = float("inf")
    for idx, c in enumerate(cells):
        dots = np.clip(grid.xyz[seed_cells] @ grid.xyz[int(c)], -1.0, 1.0)
        best = int(np.argmax(dots))
        dist = float(np.degrees(np.arccos(dots[best])))
        reference_distance[idx] = dist
        if dist < reference_seed_distance:
            reference_seed_distance = dist
            reference_seed = int(seed_cells[best])

    distance, nearest_seed = terrain._nearest_seed_distance_deg(
        grid,
        cells,
        seed_cells,
    )
    np.testing.assert_allclose(distance, reference_distance, rtol=0.0, atol=1e-12)
    assert nearest_seed == reference_seed


def test_p11222_arc_plume_vectorized_plume_distance_preserves_objects(monkeypatch):
    spec = get_preset("earthlike")
    spec.grid_cells = 420
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 76.0
    ocean = ~land
    surface = np.full(grid.n, -2400.0, dtype=np.float64)
    surface[land] = 450.0
    chain = ocean & (np.abs(grid.lat) < 12.0) & (np.abs(grid.lon) < 36.0)
    assert chain.sum() >= 3

    crust_type = np.zeros(grid.n, dtype=np.float64)
    crust_domain = np.full(grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    crust_origin = np.zeros(grid.n, dtype=np.float64)
    crust_thick = np.full(grid.n, 7000.0, dtype=np.float64)
    sediment = np.full(grid.n, 420.0, dtype=np.float64)
    world.fields["crust.age_myr"] = np.where(ocean, 80.0, 1600.0)
    volcanism_age = np.full(grid.n, -1.0, dtype=np.float64)
    chain_cells = np.where(chain)[0].astype(np.int64)
    volcanism_age[chain_cells] = np.linspace(4380.0, 4490.0, chain_cells.size)
    world.fields["tectonics.volcanism_age_myr"] = volcanism_age
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.continent_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.rift_potential"] = np.zeros(grid.n, dtype=np.float64)

    seed_a = int(np.argmin((grid.lat - 0.0) ** 2 + (grid.lon - 0.0) ** 2))
    seed_b = int(np.argmin((grid.lat - 4.0) ** 2 + (grid.lon - 70.0) ** 2))
    world.objects["tectonics.plumes"] = [
        {"id": "plume:a", "cell": seed_a, "cells": [seed_a]},
        {"id": "plume:b", "cell": seed_b, "cells": [seed_b]},
    ]
    world.objects["tectonics.volcanoes"] = [{
        "id": "volcano:chain",
        "cells": chain_cells.astype(int).tolist(),
    }]
    ocean_geo_fields = {
        "ocean.depth_province": np.where(
            ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.float64),
    }

    vector_objects = terrain._arc_plume_landform_objects(
        world,
        surface,
        sea_level,
        crust_type,
        crust_domain,
        crust_origin,
        crust_thick,
        sediment,
        ocean_geo_fields,
    )
    assert any(
        obj.get("kind") in {"seamount_chain", "hotspot_track"}
        for obj in vector_objects
    )

    def reference_nearest_seed_distance(self, grid_arg, cells_arg, seed_cells_arg):
        cells_arg = np.asarray(cells_arg, dtype=np.int64)
        seed_cells_arg = np.asarray(seed_cells_arg, dtype=np.int64)
        distance = np.full(cells_arg.size, np.nan, dtype=np.float64)
        if cells_arg.size == 0 or seed_cells_arg.size == 0:
            return distance, None
        nearest_seed = None
        nearest_seed_distance = float("inf")
        for idx, c in enumerate(cells_arg):
            dots = np.clip(
                grid_arg.xyz[seed_cells_arg] @ grid_arg.xyz[int(c)],
                -1.0,
                1.0,
            )
            best = int(np.argmax(dots))
            dist = float(np.degrees(np.arccos(dots[best])))
            distance[idx] = dist
            if dist < nearest_seed_distance:
                nearest_seed_distance = dist
                nearest_seed = int(seed_cells_arg[best])
        return distance, nearest_seed

    monkeypatch.setattr(
        TerrainModule,
        "_nearest_seed_distance_deg",
        reference_nearest_seed_distance,
    )
    reference_objects = terrain._arc_plume_landform_objects(
        world,
        surface,
        sea_level,
        crust_type,
        crust_domain,
        crust_origin,
        crust_thick,
        sediment,
        ocean_geo_fields,
    )

    def assert_same(left, right, path="root"):
        if isinstance(left, dict):
            assert set(left) == set(right), path
            for key in sorted(left):
                assert_same(left[key], right[key], f"{path}.{key}")
            return
        if isinstance(left, list):
            assert len(left) == len(right), path
            for index, (l_item, r_item) in enumerate(zip(left, right)):
                assert_same(l_item, r_item, f"{path}[{index}]")
            return
        if isinstance(left, float) or isinstance(right, float):
            l_val = float(left)
            r_val = float(right)
            if math.isnan(l_val) and math.isnan(r_val):
                return
            np.testing.assert_allclose(l_val, r_val, rtol=0.0, atol=1e-12)
            return
        assert left == right, path

    assert_same(vector_objects, reference_objects)


def test_terrain_related_object_ids_profile_matches_unprofiled_lookup():
    spec = get_preset("earthlike")
    spec.grid_cells = 160
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    parent_cell = 17
    nearby_cell = int(grid.neighbors[parent_cell][0])
    world.objects["tectonics.boundary_objects"] = [{
        "id": "boundary:test",
        "cells": [parent_cell],
    }]

    direct_mask = np.zeros(grid.n, dtype=bool)
    direct_mask[parent_cell] = True
    direct_profile = {}
    assert terrain._related_object_ids(
        world,
        direct_mask,
        ("tectonics.boundary_objects",),
    ) == terrain._related_object_ids(
        world,
        direct_mask,
        ("tectonics.boundary_objects",),
        profile_stage_seconds=direct_profile,
        profile_prefix="direct",
    )
    assert {
        "direct.overlap_index",
        "direct.direct_overlap",
    }.issubset(set(direct_profile))

    expanded_mask = np.zeros(grid.n, dtype=bool)
    expanded_mask[nearby_cell] = True
    expanded_profile = {}
    assert terrain._related_object_ids(
        world,
        expanded_mask,
        ("tectonics.boundary_objects",),
        expand_passes=1,
    ) == terrain._related_object_ids(
        world,
        expanded_mask,
        ("tectonics.boundary_objects",),
        expand_passes=1,
        profile_stage_seconds=expanded_profile,
        profile_prefix="expanded",
    )
    assert {
        "expanded.overlap_index",
        "expanded.direct_overlap",
        "expanded.expanded_mask",
        "expanded.expanded_overlap",
    }.issubset(set(expanded_profile))
    assert all(value >= 0.0 for value in direct_profile.values())
    assert all(value >= 0.0 for value in expanded_profile.values())


def test_terrain_tectonic_process_mask_cache_is_copy_safe_and_invalidates():
    spec = get_preset("earthlike")
    spec.grid_cells = 96
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    world.networks["tectonics.boundaries"] = {
        "ridge": np.asarray([0, 1], dtype=np.int64),
    }
    world.objects["tectonics.boundary_objects"] = [{
        "kind": "ridge",
        "cells": [2],
    }]
    world.objects["tectonics.spreading_centers"] = [{
        "id": "ridge:existing",
        "cells": [3],
    }]
    world.objects["tectonics.boundary_polylines"] = [{
        "id": "ridge:polyline",
        "kind": "ridge_polyline",
        "boundary_kind": "ridge",
        "cells": [5],
    }]

    first = terrain._tectonic_process_mask(world, "ridge")
    assert first[[0, 1, 2, 3, 5]].all()
    assert world.g("terrain.last_p130_ridge_polyline_used") == 1.0
    assert world.g("terrain.last_p130_ridge_polyline_cell_count") == 1.0

    first[:] = False
    second = terrain._tectonic_process_mask(world, "ridge")
    assert second[[0, 1, 2, 3, 5]].all()
    assert len(getattr(terrain, "_tectonic_process_mask_cache")) == 1

    world.objects["tectonics.spreading_centers"].append({
        "id": "ridge:new",
        "cells": [4],
    })
    third = terrain._tectonic_process_mask(world, "ridge")
    assert third[[0, 1, 2, 3, 4, 5]].all()
    assert len(getattr(terrain, "_tectonic_process_mask_cache")) == 2
    world.objects["tectonics.boundary_polylines"].append({
        "id": "ridge:polyline:new",
        "kind": "ridge_polyline",
        "boundary_kind": "ridge",
        "cells": [6],
    })
    fourth = terrain._tectonic_process_mask(world, "ridge")
    assert fourth[[0, 1, 2, 3, 4, 5, 6]].all()
    assert len(getattr(terrain, "_tectonic_process_mask_cache")) == 3

    rows = terrain._related_object_context_rows(
        world, ("tectonics.spreading_centers",)
    )
    assert terrain._related_object_context_rows(
        world, ("tectonics.spreading_centers",)
    ) is rows

    world.objects["tectonics.spreading_centers"].append({
        "id": "ridge:newer",
        "cells": [5],
    })
    invalidated_rows = terrain._related_object_context_rows(
        world, ("tectonics.spreading_centers",)
    )
    assert invalidated_rows is not rows
    assert any(row["id"] == "ridge:newer" for row in invalidated_rows)

    overlap_mask = np.zeros(grid.n, dtype=bool)
    overlap_mask[[4, 5]] = True
    assert terrain._overlapping_object_ids(
        world,
        overlap_mask,
        ("tectonics.spreading_centers",),
    ) == ["ridge:new", "ridge:newer"]


def test_terrain_components_large_mask_matches_reference_dfs():
    grid = SphereGrid.fibonacci(900, CONSTANTS.EARTH_RADIUS)
    terrain = TerrainModule()

    def reference_components(mask):
        mask = np.asarray(mask, dtype=bool)
        nodes = np.flatnonzero(mask)
        seen = np.zeros(grid.n, dtype=bool)
        comps = []
        for start in nodes:
            if seen[start]:
                continue
            stack = [int(start)]
            seen[start] = True
            comp = []
            while stack:
                c = stack.pop()
                comp.append(c)
                for nb in grid.neighbors[c]:
                    nb = int(nb)
                    if mask[nb] and not seen[nb]:
                        seen[nb] = True
                        stack.append(nb)
            comps.append(np.asarray(comp, dtype=np.int64))
        return sorted([tuple(sorted(comp.tolist())) for comp in comps])

    rng = np.random.default_rng(123)
    masks = [
        (grid.lat > -12.0) & (grid.lat < 35.0),
        (np.abs(grid.lon) < 55.0) | (grid.lat > 52.0),
    ]
    random_mask = np.zeros(grid.n, dtype=bool)
    random_mask[rng.choice(grid.n, size=260, replace=False)] = True
    masks.append(random_mask)
    threshold_mask = np.zeros(grid.n, dtype=bool)
    threshold_mask[rng.choice(grid.n, size=150, replace=False)] = True
    masks.append(threshold_mask)

    for mask in masks:
        actual = sorted(
            tuple(sorted(comp.tolist()))
            for comp in terrain._components(grid, mask)
        )
        assert actual == reference_components(mask)
    assert getattr(terrain, "_component_adjacency_graph_cache")[0][1] == grid.n


def test_terrain_merge_small_drainage_basins_matches_reference_formula():
    grid = SphereGrid.fibonacci(240, CONSTANTS.EARTH_RADIUS)
    terrain = TerrainModule()
    land = grid.lat > -58.0
    basin_id = np.full(grid.n, -1, dtype=np.int32)
    land_idx = np.where(land)[0]
    raw = (
        ((grid.lon[land] + 180.0) // 24.0).astype(int)
        + 17 * ((grid.lat[land] + 90.0) // 18.0).astype(int)
    )
    basin_id[land_idx] = raw.astype(np.int32)
    area = np.asarray(grid.cell_area, dtype=np.float64)

    def reference_merge():
        labels = basin_id.copy()
        land_area = max(float(area[land].sum()), 1.0)
        min_area = 0.025 * land_area
        for _ in range(256):
            ids = np.array(sorted(int(v) for v in np.unique(labels[land])
                                  if int(v) >= 0), dtype=np.int32)
            if ids.size <= 1:
                break
            basin_area = {
                int(bid): float(area[land & (labels == int(bid))].sum())
                for bid in ids
            }
            candidate_ids = [
                int(bid) for bid in sorted(ids, key=lambda b: basin_area[int(b)])
                if basin_area[int(bid)] < min_area or ids.size > 12
            ]
            if not candidate_ids:
                break
            merged = False
            for bid in candidate_ids:
                neighbor_counts = {}
                cells = np.where(land & (labels == bid))[0]
                for c in cells:
                    for nb in grid.neighbors[int(c)]:
                        nb = int(nb)
                        other = int(labels[nb])
                        if other >= 0 and other != bid:
                            neighbor_counts[other] = neighbor_counts.get(other, 0) + 1
                if not neighbor_counts:
                    continue
                target = max(
                    neighbor_counts,
                    key=lambda other: (
                        neighbor_counts[other],
                        basin_area.get(other, 0.0),
                        -other,
                    ),
                )
                labels[labels == bid] = int(target)
                merged = True
                break
            if not merged:
                break

        compact = np.full(grid.n, -1, dtype=np.int32)
        final_ids = [
            int(v) for v in np.unique(labels[land])
            if int(v) >= 0
        ]
        final_ids.sort(key=lambda bid: -float(area[land & (labels == bid)].sum()))
        for new_id, old_id in enumerate(final_ids):
            compact[land & (labels == old_id)] = int(new_id)
        return compact

    optimized = terrain._merge_small_drainage_basins(
        grid, basin_id, land, area, target_count=12, min_area_fraction=0.025
    )
    np.testing.assert_array_equal(optimized, reference_merge())


def test_p110b_lineage_survival_metrics_detect_supported_major_continents():
    spec = get_preset("earthlike")
    spec.grid_cells = 4000
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    land_a = (
        (grid.lat > 20.0)
        & (grid.lat < 66.0)
        & (grid.lon > -165.0)
        & (grid.lon < -85.0)
    )
    land_b = (
        (grid.lat > 12.0)
        & (grid.lat < 60.0)
        & (grid.lon > -20.0)
        & (grid.lon < 65.0)
    )
    land_c = (
        (grid.lat > -66.0)
        & (grid.lat < -20.0)
        & (grid.lon > 80.0)
        & (grid.lon < 165.0)
    )
    land = land_a | land_b | land_c
    assert int(land_a.sum()) > 20
    assert int(land_b.sum()) > 20
    assert int(land_c.sum()) > 20

    elevation = np.full(grid.n, -3800.0, dtype=np.float64)
    elevation[land] = 520.0
    domain = np.full(grid.n, 0.0, dtype=np.float64)
    domain[land] = DOMAIN_CONTINENTAL_INTERIOR
    stability = np.full(grid.n, 0.15, dtype=np.float64)
    stability[land] = 0.55
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[land_a] = 11.0
    continent_id[land_b] = 22.0
    continent_id[land_c] = 33.0
    for mask in (land_a, land_b, land_c):
        core = np.where(mask)[0][: max(3, int(mask.sum() * 0.20))]
        domain[core] = DOMAIN_CRATON
        stability[core] = 0.86

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 8).astype(np.float64),
        "terrain.elevation_m": elevation,
        "crust.type": land.astype(np.float64),
        "crust.domain": domain,
        "crust.stability": stability,
        "crust.age_myr": np.where(land, 1500.0, 80.0),
        "ocean.depth_province": np.where(land, 0.0, 4.0),
        "terrain.province": np.where(land, 2.0, 0.0),
        "terrain.continental_detail": np.where(land, 2.0, 0.0),
        "terrain.mountain_inventory": np.zeros(grid.n, dtype=np.float64),
        "tectonics.continent_id": continent_id,
        "tectonics.terrane_id": np.full(grid.n, -1.0, dtype=np.float64),
    })
    world.set_g("terrain.last_p110b_final_state_archetype_code", 1.0)
    world.set_g("terrain.last_p110b_largest_share_preferred_ceiling", 0.54)
    world.set_g("terrain.last_p110b_largest_share_soft_ceiling", 0.59)
    world.set_g("terrain.last_p110b_min_nonlargest_large_components", 2.0)

    metrics = terminal_audit_metrics(world)
    lineage = metrics["p110a_modern_planform"]["p110b_lineage_survival"]

    assert lineage["major_component_count"] == 3
    assert lineage["lineage_supported_major_component_count"] == 3
    assert lineage["independent_primary_lineage_count"] == 3
    assert lineage["min_major_dominant_continent_share"] >= 0.95
    assert lineage["min_major_geologic_support_fraction"] > 0.05
    assert lineage["unsupported_major_component_indices"] == []


def test_p110a_modern_planform_allows_balanced_three_major_continents():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    def cap(lat_deg: float, lon_deg: float, radius_deg: float) -> np.ndarray:
        lat = np.deg2rad(lat_deg)
        lon = np.deg2rad(lon_deg)
        center = np.array([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ])
        return (grid.xyz @ center) >= np.cos(np.deg2rad(radius_deg))

    land_a = cap(5.0, -90.0, 40.0)
    land_b = cap(20.0, 45.0, 29.0)
    land_c = cap(-35.0, 80.0, 29.0)
    land = land_a | land_b | land_c
    ocean = ~land

    elevation = np.where(land, 500.0, -4200.0).astype(np.float64)
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[land_a & (grid.lon < -90.0)] = 10.0
    continent_id[land_a & (grid.lon >= -90.0)] = 11.0
    continent_id[land_b] = 20.0
    continent_id[land_c] = 30.0
    domain = np.where(
        land,
        DOMAIN_CONTINENTAL_INTERIOR,
        DOMAIN_OCEANIC,
    ).astype(np.float64)
    stability = np.where(land, 0.58, 0.12).astype(np.float64)
    for mask in (land_a, land_b, land_c):
        core = np.where(mask)[0][: max(3, int(mask.sum() * 0.10))]
        domain[core] = DOMAIN_CRATON
        stability[core] = 0.86

    basin_id = np.full(grid.n, -1.0, dtype=np.float64)
    basin_id[ocean & (grid.lon < -60.0)] = 0.0
    basin_id[ocean & (grid.lon >= -60.0) & (grid.lon < 60.0)] = 1.0
    basin_id[ocean & (grid.lon >= 60.0)] = 2.0

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 12).astype(np.float64),
        "terrain.elevation_m": elevation,
        "crust.type": land.astype(np.float64),
        "crust.domain": domain,
        "crust.stability": stability,
        "crust.age_myr": np.where(land, 1600.0, 90.0),
        "ocean.depth_province": np.where(ocean, 4.0, 0.0),
        "ocean.basin_id": basin_id,
        "terrain.province": np.where(land, 2.0, 0.0),
        "terrain.continental_detail": np.where(land, 2.0, 0.0),
        "terrain.mountain_inventory": np.zeros(grid.n, dtype=np.float64),
        "tectonics.continent_id": continent_id,
        "tectonics.terrane_id": np.full(grid.n, -1.0, dtype=np.float64),
    })
    world.set_g("terrain.last_p110b_final_state_archetype_code", 1.0)
    world.set_g("terrain.last_p110b_largest_share_preferred_ceiling", 0.54)
    world.set_g("terrain.last_p110b_largest_share_soft_ceiling", 0.59)
    world.set_g("terrain.last_p110b_min_nonlargest_large_components", 2.0)

    metrics = terminal_audit_metrics(world, [])
    plan = metrics["p110a_modern_planform"]
    diag = plan["p110b_terminal_supercontinent_diagnostics"]

    assert plan["summary"]["major_land_component_count"] == 3
    assert plan["summary"]["largest_land_component_share"] <= 0.54
    assert plan["summary"]["second_land_component_share"] >= 0.20
    assert plan["summary"]["third_land_component_share"] >= 0.16
    assert not diag["terminal_supercontinent_like"]
    assert not diag["modern_multipolar_overconnected"]
    assert diag["largest_land_significant_continent_domain_count"] >= 2
    assert plan["out_of_envelope"] == []
    assert (
        "major_land_component_count_low_p110b_modern_balanced_three_major"
        in plan["warning_flags"]
    )
    assert plan["within_p110a_modern_planform_envelope"]


def test_p110b_terminal_supercontinent_diagnostics_flags_monolithic_landmass():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    land = (
        (np.abs(grid.lat) < 58.0)
        & (grid.lon > -125.0)
        & (grid.lon < 125.0)
    )
    assert int(land.sum()) > 200
    elevation = np.full(grid.n, -3800.0, dtype=np.float64)
    elevation[land] = 420.0
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[land] = 7.0

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 8).astype(np.float64),
        "terrain.elevation_m": elevation,
        "crust.type": land.astype(np.float64),
        "crust.domain": np.where(land, DOMAIN_CONTINENTAL_INTERIOR, 0.0),
        "crust.stability": np.where(land, 0.52, 0.1),
        "crust.age_myr": np.where(land, 1600.0, 80.0),
        "ocean.depth_province": np.where(land, 0.0, 4.0),
        "terrain.province": np.where(land, 2.0, 0.0),
        "terrain.continental_detail": np.where(land, 2.0, 0.0),
        "terrain.mountain_inventory": np.zeros(grid.n, dtype=np.float64),
        "tectonics.continent_id": continent_id,
        "tectonics.terrane_id": np.full(grid.n, -1.0, dtype=np.float64),
    })
    world.set_g("terrain.last_p110b_final_state_archetype_code", 1.0)
    world.set_g("terrain.last_p110b_largest_share_preferred_ceiling", 0.54)
    world.set_g("terrain.last_p110b_largest_share_soft_ceiling", 0.59)
    world.set_g("terrain.last_p110b_min_nonlargest_large_components", 2.0)

    metrics = terminal_audit_metrics(world)
    plan = metrics["p110a_modern_planform"]
    diag = plan["p110b_terminal_supercontinent_diagnostics"]

    assert diag["terminal_supercontinent_like"]
    assert diag["monolithic_large_component"]
    assert diag["largest_land_significant_continent_domain_count"] == 1
    assert diag["largest_land_effective_continent_domain_count"] == 1.0
    assert diag["terminal_supercontinent_score"] >= 0.5
    assert plan["modern_planform_applicable"]
    assert "largest_land_component_share_hard" in plan["out_of_envelope"]
    assert "p110b_terminal_supercontinent_like" in plan["warning_flags"]


def test_p110a_modern_planform_defers_hard_gate_before_modern_window():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=1200.0)

    land = (
        (np.abs(grid.lat) < 58.0)
        & (grid.lon > -125.0)
        & (grid.lon < 125.0)
    )
    elevation = np.full(grid.n, -3800.0, dtype=np.float64)
    elevation[land] = 420.0
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[land] = 7.0

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 8).astype(np.float64),
        "terrain.elevation_m": elevation,
        "crust.type": land.astype(np.float64),
        "crust.domain": np.where(land, DOMAIN_CONTINENTAL_INTERIOR, 0.0),
        "crust.stability": np.where(land, 0.52, 0.1),
        "crust.age_myr": np.where(land, 1600.0, 80.0),
        "ocean.depth_province": np.where(land, 0.0, 4.0),
        "terrain.province": np.where(land, 2.0, 0.0),
        "terrain.continental_detail": np.where(land, 2.0, 0.0),
        "terrain.mountain_inventory": np.zeros(grid.n, dtype=np.float64),
        "tectonics.continent_id": continent_id,
        "tectonics.terrane_id": np.full(grid.n, -1.0, dtype=np.float64),
    })
    world.set_g("terrain.last_p110b_final_state_archetype_code", 1.0)
    world.set_g("terrain.last_p110b_largest_share_preferred_ceiling", 0.54)
    world.set_g("terrain.last_p110b_largest_share_soft_ceiling", 0.59)
    world.set_g("terrain.last_p110b_min_nonlargest_large_components", 2.0)

    metrics = terminal_audit_metrics(world, [])
    plan = metrics["p110a_modern_planform"]

    assert not plan["modern_planform_applicable"]
    assert plan["modern_planform_applicable_after_myr"] == 3500.0
    assert plan["out_of_envelope"] == []
    assert "land_fraction" in plan["deferred_out_of_envelope"]
    assert "largest_land_component_share_hard" in plan["deferred_out_of_envelope"]
    assert "modern_planform_gate_deferred_until_3500myr" in plan["warning_flags"]
    assert (
        "immature_modern_planform_largest_land_component_share_hard"
        in plan["warning_flags"]
    )
    assert plan["within_p110a_modern_planform_envelope"]


def test_p110b_historical_supercontinent_trajectory_flags_long_lived_lock():
    spec = get_preset("earthlike")
    spec.grid_cells = 360
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    frames = []
    for t in (0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0):
        elev = np.full(grid.n, 250.0, dtype=np.float64)
        frames.append(SimpleNamespace(
            time_myr=t,
            globals={"ocean.sea_level_m": 0.0},
            fields={"terrain.elevation_m": elev},
        ))
    archive = SimpleNamespace(world=world, frames=frames)

    metrics = _p110b_historical_supercontinent_trajectory_metrics(archive)

    assert not metrics["skipped"]
    assert metrics["supercontinent_frame_count"] == 6
    assert metrics["supercontinent_frame_fraction"] == 1.0
    assert metrics["max_consecutive_supercontinent_duration_myr"] == 2500.0
    assert metrics["long_lived_supercontinent_like"]


def test_p110b_historical_supercontinent_trajectory_flags_wide_recurrence_window():
    spec = get_preset("earthlike")
    spec.grid_cells = 360
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    super_land = (
        (np.abs(grid.lat) < 54.0)
        & (grid.lon > -145.0)
        & (grid.lon < 145.0)
    )
    split_land = (
        ((grid.lat > 18.0) & (grid.lat < 68.0) & (grid.lon > -155.0) & (grid.lon < -60.0))
        | ((grid.lat > -68.0) & (grid.lat < -18.0) & (grid.lon > 60.0) & (grid.lon < 155.0))
        | ((np.abs(grid.lat) < 20.0) & (grid.lon > -25.0) & (grid.lon < 25.0))
    )
    super_elev = np.where(super_land, 250.0, -2400.0).astype(np.float64)
    split_elev = np.where(split_land, 250.0, -2400.0).astype(np.float64)

    frames = []
    for t in (0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4500.0):
        elev = super_elev if t in {500.0, 1000.0, 3000.0} else split_elev
        frames.append(SimpleNamespace(
            time_myr=t,
            globals={"ocean.sea_level_m": 0.0},
            fields={"terrain.elevation_m": elev},
        ))
    archive = SimpleNamespace(world=world, frames=frames)

    metrics = _p110b_historical_supercontinent_trajectory_metrics(archive)

    assert metrics["supercontinent_frame_count"] == 3
    assert metrics["supercontinent_frame_fraction"] < 0.45
    assert metrics["max_consecutive_supercontinent_duration_myr"] < 1200.0
    assert metrics["supercontinent_time_window_myr"] > 1800.0
    assert metrics["wide_recurrent_supercontinent_window"]
    assert not metrics["long_lived_supercontinent_like"]


def test_p110b_historical_breakup_pressure_opens_mid_history_seaway():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=1000.0)
    tectonics = TectonicsModule()
    terrain = TerrainModule()

    continent = (
        (np.abs(grid.lat) < 55.0)
        & (grid.lon > -145.0)
        & (grid.lon < 145.0)
    )
    crust = np.where(continent, CONT, OCEAN).astype(np.float64)
    domain = np.where(
        continent, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    stability = np.where(continent, 0.42, 0.10).astype(np.float64)
    core = continent & (
        ((grid.lon < -80.0) & (grid.lat > 15.0))
        | ((grid.lon > 80.0) & (grid.lat < -15.0))
    )
    domain[core] = DOMAIN_CRATON
    stability[core] = 0.86
    rift = np.where(
        continent,
        0.18 + 0.20 * np.exp(-(grid.lat / 18.0) ** 2),
        0.0,
    )
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[continent] = 0.0
    plate = (np.arange(grid.n) % 12).astype(np.int64)
    world.set_g(
        "tectonics.last_p110b_historical_supercontinent_pressure_time_myr",
        700.0,
    )
    world.set_g(
        "tectonics.last_p110b_historical_supercontinent_consecutive_myr",
        600.0,
    )
    world.set_g(
        "tectonics.last_p110b_historical_supercontinent_cumulative_myr",
        900.0,
    )

    objects = tectonics._breakup_seaway_objects(
        world,
        grid,
        crust,
        domain,
        stability,
        rift,
        continent_id,
        plate,
        {},
        [],
        [],
        1000.0,
    )
    assert world.g("tectonics.last_p110b_historical_breakup_active") == 1.0
    assert len(objects) >= 1
    assert objects[0]["basis"] == "p110b_historical_supercontinent_breakup_pressure"

    world.objects["tectonics.breakup_seaways"] = objects
    surface = np.where(continent, 520.0, -3600.0).astype(np.float64)
    shaped = terrain._open_breakup_seaway_objects(
        world,
        surface,
        0.0,
        crust,
        domain,
        stability,
    )
    area = grid.cell_area
    before_land = surface >= 0.0
    after_land = shaped >= 0.0
    before_components = terrain._components(grid, before_land)
    after_components = terrain._components(grid, after_land)
    before_largest = max(float(area[c].sum()) for c in before_components) / float(
        area[before_land].sum())
    after_largest = max(float(area[c].sum()) for c in after_components) / float(
        area[after_land].sum())

    assert world.g("terrain.last_breakup_seaway_openings") >= 1.0
    assert len(after_components) > len(before_components)
    assert before_largest > 0.95
    assert after_largest < 0.70


def test_p110b_historical_breakup_pressure_builds_late_low_rift_fallback():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=2800.0)
    tectonics = TectonicsModule()

    continent = (
        (np.abs(grid.lat) < 58.0)
        & (grid.lon > -150.0)
        & (grid.lon < 150.0)
    )
    crust = np.where(continent, CONT, OCEAN).astype(np.float64)
    domain = np.where(
        continent, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    stability = np.where(continent, 0.60, 0.10).astype(np.float64)
    core = continent & (
        ((grid.lon < -72.0) & (grid.lat > 8.0))
        | ((grid.lon > 72.0) & (grid.lat < -8.0))
    )
    domain[core] = DOMAIN_CRATON
    stability[core] = 0.88
    rift = np.where(continent, 0.08, 0.0).astype(np.float64)
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[continent] = 0.0
    plate = (np.arange(grid.n) % 14).astype(np.int64)
    world.set_g(
        "tectonics.last_p110b_historical_supercontinent_pressure_time_myr",
        2400.0,
    )
    world.set_g(
        "tectonics.last_p110b_historical_supercontinent_consecutive_myr",
        900.0,
    )
    world.set_g(
        "tectonics.last_p110b_historical_supercontinent_cumulative_myr",
        1600.0,
    )

    objects = tectonics._breakup_seaway_objects(
        world,
        grid,
        crust,
        domain,
        stability,
        rift,
        continent_id,
        plate,
        {},
        [],
        [],
        2800.0,
    )

    assert world.g("tectonics.last_p110b_historical_breakup_active") == 1.0
    assert (
        world.g("tectonics.last_p110b_historical_supercontinent_residence_debt_myr")
        >= 399.0
    )
    assert world.g("tectonics.last_p110b_historical_residence_controller_active") == 1.0
    assert len(objects) >= 1
    assert objects[0]["basis"] == "p110b_historical_supercontinent_breakup_pressure"
    assert objects[0]["historical_residence_debt_myr"] >= 399.0
    assert objects[0]["historical_residence_controller_active"]
    assert objects[0]["topology_score"] >= 0.18
    assert objects[0]["weak_fraction"] >= 0.12


def test_p110b_residence_controller_opens_stronger_terrain_breakup():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=3000.0)
    terrain = TerrainModule()

    continent = (
        (np.abs(grid.lat) < 58.0)
        & (grid.lon > -150.0)
        & (grid.lon < 150.0)
    )
    surface = np.where(continent, 620.0, -3600.0).astype(np.float64)
    crust = np.where(continent, CONT, OCEAN).astype(np.float64)
    domain = np.where(
        continent, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    stability = np.where(continent, 0.68, 0.10).astype(np.float64)
    source = continent & (np.abs(grid.lon) < 5.0)
    assert int(source.sum()) >= 12
    world.objects["tectonics.breakup_seaways"] = [{
        "id": "breakup_seaway:residence:test",
        "kind": "breakup_seaway",
        "basis": "p110b_historical_supercontinent_breakup_pressure",
        "cells": np.where(source)[0].astype(int).tolist(),
        "cell_count": int(np.count_nonzero(source)),
        "quality_score": 1.2,
        "component_area_fraction": 0.22,
        "area_fraction": float(grid.cell_area[source].sum() / grid.cell_area.sum()),
        "topology_score": 0.45,
        "historical_breakup_pressure": 1.0,
        "historical_residence_debt_myr": 520.0,
        "historical_residence_controller_active": True,
    }]

    shaped = terrain._open_breakup_seaway_objects(
        world,
        surface,
        0.0,
        crust,
        domain,
        stability,
    )

    area = grid.cell_area
    before_land = surface >= 0.0
    after_land = shaped >= 0.0
    before_components = terrain._components(grid, before_land)
    after_components = terrain._components(grid, after_land)
    before_largest = max(float(area[c].sum()) for c in before_components) / float(
        area[before_land].sum())
    after_largest = max(float(area[c].sum()) for c in after_components) / float(
        area[after_land].sum())

    assert before_largest > 0.95
    assert world.g("terrain.last_breakup_seaway_openings") >= 1.0
    assert world.g("terrain.last_breakup_seaway_area_fraction") > 0.0
    assert len(after_components) > len(before_components)
    assert after_largest < 0.75
    opened = world.objects["terrain.breakup_seaway_opened_corridors"][0]
    assert opened["historical_residence_controller_active"]
    assert opened["historical_residence_debt_myr"] == 520.0


def test_p109_final_extreme_peak_floor_polish_restores_parented_summit_tail():
    spec = get_preset("earthlike")
    spec.grid_cells = 1000
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    land = grid.lat > -28.0
    surface = np.where(land, 720.0, -3600.0).astype(np.float64)
    highland = land & (grid.lat > 18.0) & (grid.lon > -80.0) & (grid.lon < 80.0)
    surface[highland] = 3300.0
    crust = np.where(land, CONT, OCEAN).astype(np.float64)
    detail = np.zeros(grid.n, dtype=np.float64)
    detail[highland] = CONT_DETAIL_OROGEN
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[highland] = 1.0

    area = grid.cell_area
    before_land = surface >= 0.0
    before_fraction = float(
        area[before_land & (surface >= 4500.0)].sum()
        / max(float(area[before_land].sum()), 1.0e-12)
    )

    shaped = terrain._p109_final_extreme_peak_floor_polish(
        world,
        surface,
        0.0,
        crust,
        detail,
        inventory,
    )

    after_land = shaped >= 0.0
    after_fraction = float(
        area[after_land & (shaped >= 4500.0)].sum()
        / max(float(area[after_land].sum()), 1.0e-12)
    )
    assert before_fraction < 0.0032
    assert after_fraction >= 0.0032
    assert np.array_equal(before_land, after_land)
    assert world.g("terrain.last_p109_final_extreme_peak_floor_area_fraction") > 0.0


def test_p110b_historical_breakup_seaway_objects_persist_for_retention_window():
    spec = get_preset("earthlike")
    spec.grid_cells = 120
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    tectonics = TectonicsModule()

    def retained_at(time_myr: float, last_active_myr: float) -> list[dict]:
        world = WorldState(grid=grid, spec=spec, time_myr=time_myr)
        world.objects["tectonics.breakup_seaways"] = [{
            "id": "breakup_seaway:test",
            "kind": "breakup_seaway",
            "basis": "p110b_historical_supercontinent_breakup_pressure",
            "lineage_key": "continent:0:history:test",
            "birth_myr": 900.0,
            "last_active_myr": last_active_myr,
            "quality_score": 0.4,
            "component_area_fraction": 0.20,
            "cells": [0, 1, 2, 3],
        }]
        ocean = np.zeros(grid.n, dtype=np.float64)
        return tectonics._breakup_seaway_objects(
            world,
            grid,
            ocean,
            np.zeros(grid.n, dtype=np.float64),
            np.full(grid.n, 0.2, dtype=np.float64),
            np.zeros(grid.n, dtype=np.float64),
            np.full(grid.n, -1.0, dtype=np.float64),
            (np.arange(grid.n) % 6).astype(np.int64),
            {},
            [],
            [],
            time_myr,
        )

    retained = retained_at(1500.0, 1000.0)
    assert len(retained) == 1
    assert retained[0]["retained_from_previous_step"]
    assert retained[0]["stage"] == "inherited_rifted_seaway"
    assert retained[0]["retention_age_myr"] == 500.0

    assert retained_at(1700.0, 1000.0) == []


def test_p110b_historical_breakup_seaway_crustalizes_weak_corridor():
    spec = get_preset("earthlike")
    spec.grid_cells = 120
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=1500.0)
    tectonics = TectonicsModule()

    ctype = np.full(grid.n, CONT, dtype=np.float64)
    thick = np.full(grid.n, CONT_THICK, dtype=np.float64)
    age = np.full(grid.n, 1600.0, dtype=np.float64)
    origin = np.full(grid.n, ORIGIN_PRIMORDIAL, dtype=np.float64)
    reworked = np.zeros(grid.n, dtype=np.float64)
    stability = np.full(grid.n, 0.42, dtype=np.float64)
    domain = np.full(grid.n, DOMAIN_CONTINENTAL_INTERIOR, dtype=np.float64)
    continent_id = np.full(grid.n, 7.0, dtype=np.float64)
    terrane_id = np.full(grid.n, -1.0, dtype=np.float64)
    province = np.full(grid.n, PROVINCE_PLATFORM, dtype=np.float64)
    block_id = np.ones(grid.n, dtype=np.float64)
    block_code = np.full(
        grid.n, INTERNAL_BLOCK_RIFTED_MARGIN, dtype=np.float64)
    basement_age = np.full(grid.n, 2000.0, dtype=np.float64)
    basement_stability = np.full(grid.n, 0.55, dtype=np.float64)
    basement_thick = np.full(grid.n, CONT_THICK, dtype=np.float64)
    orog_age = np.full(grid.n, 900.0, dtype=np.float64)
    volc_age = np.full(grid.n, 700.0, dtype=np.float64)
    domain[0] = DOMAIN_CRATON
    stability[0] = 0.95

    tectonics._apply_breakup_seaway_crustal_opening(
        world,
        grid,
        [{
            "basis": "p110b_historical_supercontinent_breakup_pressure",
            "birth_myr": 1200.0,
            "cells": [0, 1, 2, 3],
        }],
        ctype,
        thick,
        age,
        origin,
        reworked,
        stability,
        domain,
        continent_id,
        terrane_id,
        province,
        block_id,
        block_code,
        basement_age,
        basement_stability,
        basement_thick,
        orog_age,
        volc_age,
        1500.0,
    )

    assert ctype[0] == CONT
    assert domain[0] == DOMAIN_CRATON
    assert np.all(ctype[1:4] == OCEAN)
    assert np.all(domain[1:4] == DOMAIN_OCEANIC)
    assert np.all(thick[1:4] <= OCEAN_THICK + 1200.0)
    assert np.all(age[1:4] <= 300.0)
    assert np.all(origin[1:4] == ORIGIN_RIDGE)
    assert np.all(continent_id[1:4] == -1.0)
    assert np.all(block_code[1:4] == INTERNAL_BLOCK_NONE)
    assert world.g("tectonics.last_p110b_breakup_crustal_opening_cells") == 3.0


def test_p1134_lowres_false_oceanic_island_cleanup_preserves_supported_islands():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.target_land_fraction = 0.30
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=spec.t_end_myr)
    terrain = TerrainModule()
    sea_level = 0.0

    base_land = grid.lat > 8.0
    far_from_land = ~terrain._dilate_mask(grid, base_land, passes=2)
    ocean_pool = np.where(far_from_land & (grid.lat < -25.0))[0]
    unsupported = int(ocean_pool[0])
    unsupported_mask = np.zeros(grid.n, dtype=bool)
    unsupported_mask[unsupported] = True
    far_from_unsupported = ~terrain._dilate_mask(
        grid,
        unsupported_mask,
        passes=4,
    )
    supported = int(ocean_pool[far_from_unsupported[ocean_pool]][0])

    surface = np.full(grid.n, sea_level - 4200.0, dtype=np.float64)
    surface[base_land] = sea_level + 430.0
    surface[[unsupported, supported]] = sea_level + 70.0
    crust_type = np.full(grid.n, OCEAN, dtype=np.float64)
    crust_type[base_land] = CONT
    crust_domain = np.full(grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    crust_domain[base_land] = DOMAIN_CONTINENTAL_INTERIOR
    crust_origin = np.full(grid.n, ORIGIN_RIDGE, dtype=np.float64)
    crust_origin[base_land] = ORIGIN_PRIMORDIAL
    world.objects["tectonics.plumes"] = [{
        "cell": supported,
        "cells": [supported],
    }]
    world.networks["tectonics.boundaries"] = {
        "active_margin": [unsupported],
        "trench": [unsupported],
    }

    adjusted, telemetry = terrain._p1134_lowres_false_oceanic_island_cleanup(
        world,
        surface,
        sea_level,
        crust_type,
        crust_domain,
        crust_origin,
    )

    assert telemetry["accepted"] is True
    assert telemetry["adjusted_component_count"] == 1.0
    assert adjusted[unsupported] < sea_level
    assert adjusted[supported] >= sea_level
    assert telemetry["component_count_after"] < telemetry["component_count_before"]


def test_p1134_lowres_false_oceanic_island_cleanup_respects_land_floor():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.target_land_fraction = 0.30
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=spec.t_end_myr)
    terrain = TerrainModule()
    sea_level = 0.0

    island = int(np.argmin(grid.lat))
    land_cells = np.argsort(-grid.lat)[:229]
    surface = np.full(grid.n, sea_level - 4200.0, dtype=np.float64)
    surface[land_cells] = sea_level + 430.0
    surface[island] = sea_level + 70.0
    crust_type = np.full(grid.n, OCEAN, dtype=np.float64)
    crust_type[land_cells] = CONT
    crust_domain = np.full(grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    crust_domain[land_cells] = DOMAIN_CONTINENTAL_INTERIOR
    crust_origin = np.full(grid.n, ORIGIN_RIDGE, dtype=np.float64)
    crust_origin[land_cells] = ORIGIN_PRIMORDIAL

    adjusted, telemetry = terrain._p1134_lowres_false_oceanic_island_cleanup(
        world,
        surface,
        sea_level,
        crust_type,
        crust_domain,
        crust_origin,
    )

    assert telemetry["accepted"] is False
    assert telemetry["rejected_by_land_floor"] is True
    assert adjusted[island] == surface[island]
    assert telemetry["land_fraction_after"] == telemetry["land_fraction_before"]


def test_p110b_lineage_survival_separates_exposed_oceanic_arcs_from_continents():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    land_a = (
        (grid.lat > 20.0)
        & (grid.lat < 66.0)
        & (grid.lon > -165.0)
        & (grid.lon < -90.0)
    )
    land_b = (
        (grid.lat > -60.0)
        & (grid.lat < -18.0)
        & (grid.lon > 75.0)
        & (grid.lon < 160.0)
    )
    island_arc = (
        (grid.lat > -8.0)
        & (grid.lat < 26.0)
        & (grid.lon > -25.0)
        & (grid.lon < 35.0)
    )
    land = land_a | land_b | island_arc
    assert int(land_a.sum()) > 20
    assert int(land_b.sum()) > 20
    assert int(island_arc.sum()) > 20

    sea_level = -2100.0
    elevation = np.full(grid.n, -4200.0, dtype=np.float64)
    elevation[land_a | land_b] = -1000.0
    elevation[island_arc] = -1600.0
    crust_type = np.zeros(grid.n, dtype=np.float64)
    crust_type[land_a | land_b] = 1.0
    domain = np.zeros(grid.n, dtype=np.float64)
    domain[land_a | land_b] = DOMAIN_CRATON
    stability = np.full(grid.n, 0.15, dtype=np.float64)
    stability[land_a | land_b] = 0.86
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[land_a] = 101.0
    continent_id[land_b] = 202.0

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 10).astype(np.float64),
        "terrain.elevation_m": elevation,
        "crust.type": crust_type,
        "crust.domain": domain,
        "crust.stability": stability,
        "crust.age_myr": np.where(crust_type >= 0.5, 1800.0, 80.0),
        "ocean.depth_province": np.where(land, 0.0, 4.0),
        "terrain.province": np.where(land, 2.0, 0.0),
        "terrain.continental_detail": np.where(land_a | land_b, 2.0, 7.0),
        "terrain.mountain_inventory": np.zeros(grid.n, dtype=np.float64),
        "tectonics.continent_id": continent_id,
        "tectonics.terrane_id": np.full(grid.n, -1.0, dtype=np.float64),
    })
    world.objects["terrain.arc_plume_landforms"] = [
        {"kind": "island_arc", "cells": np.where(island_arc)[0].astype(int).tolist()}
    ]
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(island_arc)[0].astype(np.int64)
    }
    world.set_g("ocean.sea_level_m", sea_level)
    world.set_g("terrain.last_p110b_final_state_archetype_code", 1.0)
    world.set_g("terrain.last_p110b_largest_share_preferred_ceiling", 0.54)
    world.set_g("terrain.last_p110b_largest_share_soft_ceiling", 0.59)
    world.set_g("terrain.last_p110b_min_nonlargest_large_components", 2.0)

    metrics = terminal_audit_metrics(world)
    plan = metrics["p110a_modern_planform"]
    lineage = plan["p110b_lineage_survival"]

    assert lineage["major_component_count"] == 3
    assert lineage["continental_major_component_count"] == 2
    assert lineage["continental_lineage_supported_major_component_count"] == 2
    assert lineage["oceanic_landform_major_component_count"] == 1
    assert lineage["provenance_supported_major_component_count"] == 3
    assert lineage["unsupported_major_component_indices"] == []
    assert lineage["unsupported_continental_major_component_indices"] == []
    assert any(
        row["lineage_support_basis"] == "oceanic_landform_provenance"
        and row["component_class"] == "oceanic_exposed_landform"
        for row in lineage["component_rows"]
    )
    assert "p110b_lineage_support_incomplete" not in plan["warning_flags"]
    assert "p110b_landform_provenance_incomplete" not in plan["warning_flags"]


def test_p110b_terminal_continent_lineage_fission_labels_separated_children():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    parent_main = (
        (grid.lat > 20.0)
        & (grid.lat < 66.0)
        & (grid.lon > -165.0)
        & (grid.lon < -70.0)
    )
    parent_child = (
        (grid.lat > 10.0)
        & (grid.lat < 56.0)
        & (grid.lon > 10.0)
        & (grid.lon < 80.0)
    )
    other_parent = (
        (grid.lat > -66.0)
        & (grid.lat < -20.0)
        & (grid.lon > 95.0)
        & (grid.lon < 165.0)
    )
    land = parent_main | parent_child | other_parent
    assert parent_main.sum() > parent_child.sum() > 10
    surface = np.where(land, 600.0, -3600.0).astype(np.float64)
    crust = land.astype(np.float64)
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[parent_main] = 5.0
    continent_id[parent_child] = 5.0
    continent_id[other_parent] = 6.0

    new_ids, next_id = terrain._p110b_terminal_continent_lineage_fission(
        world,
        surface,
        0.0,
        crust,
        continent_id,
        10,
    )

    child_values = set(int(x) for x in np.unique(new_ids[parent_child]))
    assert set(int(x) for x in np.unique(new_ids[parent_main])) == {5}
    assert set(int(x) for x in np.unique(new_ids[other_parent])) == {6}
    assert child_values == {10}
    assert next_id == 11
    splits = world.objects["tectonics.continent_lineage_splits"]
    assert len(splits) == 1
    assert splits[0]["parent_continent_id"] == 5
    assert splits[0]["child_continent_id"] == 10
    assert world.g("terrain.last_p110b_lineage_fission_new_continent_count") == 1.0


def test_p110b_terminal_continent_lineage_fission_splits_monolithic_landmass():
    spec = get_preset("earthlike")
    spec.grid_cells = 1200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    land = (
        (np.abs(grid.lat) < 50.0)
        & (grid.lon > -145.0)
        & (grid.lon < 145.0)
    )
    assert int(land.sum()) > 200
    surface = np.where(land, 650.0, -3600.0).astype(np.float64)
    crust = land.astype(np.float64)
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[land] = 5.0

    new_ids, next_id = terrain._p110b_terminal_continent_lineage_fission(
        world,
        surface,
        0.0,
        crust,
        continent_id,
        20,
    )

    ids = sorted(int(x) for x in np.unique(new_ids[land]) if int(x) >= 0)
    splits = world.objects["tectonics.continent_lineage_splits"]

    assert 5 in ids
    assert len(ids) >= 2
    assert next_id > 20
    assert any(
        split["parent_process"] == "p110b_terminal_internal_domain_fission"
        for split in splits
    )
    assert world.g(
        "terrain.last_p110b_internal_lineage_fission_new_continent_count"
    ) >= 1.0


def test_p110b_internal_domain_seaway_polish_opens_protected_corridor():
    spec = get_preset("earthlike")
    spec.grid_cells = 1200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    world.set_g("terrain.enable_p110b_internal_domain_seaway_opening", 1.0)
    terrain = TerrainModule()

    land = (
        (np.abs(grid.lat) < 34.0)
        & (grid.lon > -142.0)
        & (grid.lon < 142.0)
    )
    assert int(land.sum()) > 200
    surface = np.where(land, 520.0, -3600.0).astype(np.float64)
    crust = land.astype(np.float64)
    domain = np.where(land, DOMAIN_CONTINENTAL_INTERIOR, 0.0).astype(np.float64)
    stability = np.where(land, 0.45, 0.1).astype(np.float64)
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[land & (grid.lon < -25.0)] = 5.0
    continent_id[land & (grid.lon >= -25.0) & (grid.lon <= 25.0)] = 6.0
    continent_id[land & (grid.lon > 25.0)] = 7.0

    shaped = terrain._p110b_internal_domain_seaway_polish(
        world,
        surface,
        0.0,
        crust,
        domain,
        stability,
        continent_id,
    )

    before_components = len(terrain._components(grid, surface >= 0.0))
    after_land = shaped >= 0.0
    after_components = len(terrain._components(grid, after_land))

    assert world.g("terrain.last_p110b_internal_domain_seaway_openings") >= 1.0
    assert after_components > before_components
    assert np.count_nonzero((surface >= 0.0) & ~after_land) > 0
    assert len(world.objects["terrain.p110b_internal_domain_seaways"]) >= 1
    assert len(world.objects["terrain.breakup_seaway_opened_corridors"]) >= 1


def test_p110b_internal_domain_polish_defaults_to_boundary_not_seaway():
    spec = get_preset("earthlike")
    spec.grid_cells = 1200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    land = (
        (np.abs(grid.lat) < 34.0)
        & (grid.lon > -142.0)
        & (grid.lon < 142.0)
    )
    surface = np.where(land, 520.0, -3600.0).astype(np.float64)
    crust = land.astype(np.float64)
    domain = np.where(land, DOMAIN_CONTINENTAL_INTERIOR, 0.0).astype(np.float64)
    stability = np.where(land, 0.45, 0.1).astype(np.float64)
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[land & (grid.lon < -25.0)] = 5.0
    continent_id[land & (grid.lon >= -25.0) & (grid.lon <= 25.0)] = 6.0
    continent_id[land & (grid.lon > 25.0)] = 7.0

    shaped = terrain._p110b_internal_domain_seaway_polish(
        world,
        surface,
        0.0,
        crust,
        domain,
        stability,
        continent_id,
    )

    assert world.g("terrain.last_p110b_internal_domain_seaway_openings") == 0.0
    assert np.array_equal(shaped >= 0.0, surface >= 0.0)
    assert len(world.objects["terrain.p110b_internal_domain_seaways"]) == 0
    assert world.g("terrain.last_p110b_internal_domain_boundary_count") >= 1.0
    assert len(world.objects["terrain.p110b_internal_domain_boundaries"]) >= 1


def test_p110a_planform_metrics_are_spherical_and_area_weighted():
    spec = get_preset("earthlike")
    spec.grid_cells = 1200
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    def center_xyz(lat_deg: float, lon_deg: float) -> np.ndarray:
        lat = np.deg2rad(lat_deg)
        lon = np.deg2rad(lon_deg)
        return np.array([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ])

    land = np.zeros(grid.n, dtype=bool)
    for lat, lon in (
        (0.0, 179.0),
        (35.0, -80.0),
        (-35.0, -70.0),
        (35.0, 35.0),
        (-35.0, 55.0),
        (0.0, -5.0),
    ):
        land |= (grid.xyz @ center_xyz(lat, lon)) >= np.cos(np.deg2rad(23.0))
    elev = np.where(land, 600.0, -4200.0).astype(np.float64)
    basin_id = np.full(grid.n, -1.0, dtype=np.float64)
    ocean = ~land
    basin_id[ocean & (grid.lon < -90.0)] = 0.0
    basin_id[ocean & (grid.lon >= -90.0) & (grid.lon < 0.0)] = 1.0
    basin_id[ocean & (grid.lon >= 0.0) & (grid.lon < 90.0)] = 2.0
    basin_id[ocean & (grid.lon >= 90.0)] = 3.0

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 12).astype(float),
        "terrain.elevation_m": elev,
        "crust.type": land.astype(float),
        "crust.age_myr": np.where(land, 1600.0, 80.0),
        "ocean.depth_province": np.where(ocean, 4.0, 0.0),
        "ocean.basin_id": basin_id,
        "terrain.province": np.where(land, 2.0, 0.0),
        "terrain.continental_detail": np.where(land, 2.0, 0.0),
    })

    metrics = terminal_audit_metrics(world, [])
    planform = metrics["p110a_modern_planform"]

    assert planform["land"]["component_count"] == 6
    assert planform["land"]["major_component_count"] == 6
    assert planform["summary"]["largest_land_component_share"] < 0.30
    assert planform["ocean_basins"]["source"] == "ocean.basin_id"
    assert planform["summary"]["major_ocean_basin_count"] >= 3
    assert planform["ocean_connectivity"]["closed_ocean_ring_score"] < 0.05


def test_p110a_modern_planform_target_includes_production_plate_counts():
    spec = get_preset("earthlike")
    spec.grid_cells = 120
    spec.n_plates = 36
    world = WorldState(
        grid=SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS),
        spec=spec,
        time_myr=4500.0,
    )
    terrain = TerrainModule()

    target = terrain._p110a_modern_planform_target(world)

    assert target["modern_earthlike"]
    assert target["production_plate_count"]
    assert target["n_plates"] == 36
    assert target["p110b_archetype"] in {
        "modern_multipolar",
        "afro_eurasia_analog",
        "post_supercontinent_breakup",
        "archipelago_active_margin",
    }
    assert target["p110b_largest_share_soft_ceiling"] < 0.65
    assert world.g("terrain.last_p110b_final_state_archetype_code") >= 1.0

    spec.n_plates = 60
    target = terrain._p110a_modern_planform_target(world)
    assert target["modern_earthlike"]
    assert target["production_plate_count"]

    world.set_g("terrain.p110b_final_state_archetype_code", 3.0)
    target = terrain._p110a_modern_planform_target(world)
    assert target["modern_earthlike"]
    assert target["p110b_archetype"] == "post_supercontinent_breakup"
    assert target["p110b_archetype_override"]
    assert target["p110b_min_nonlargest_large_components"] == 3

    world.set_g("terrain.allow_p110a_supercontinent_final_state", 1.0)
    target = terrain._p110a_modern_planform_target(world)
    assert not target["modern_earthlike"]
    assert target["supercontinent_allowed"]
    assert target["p110b_archetype"] == "supercontinent_final"


def test_p110a_ocean_basin_metrics_prefer_lifecycle_objects():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    land = np.abs(grid.lat) > 55.0
    elev = np.where(land, 500.0, -4200.0)
    ocean = ~land
    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 12).astype(float),
        "terrain.elevation_m": elev.astype(float),
        "crust.type": land.astype(float),
        "crust.age_myr": np.where(land, 1600.0, 90.0),
        "ocean.depth_province": np.where(ocean, 4.0, 0.0),
        "ocean.basin_id": np.where(ocean, 0.0, -1.0),
        "terrain.province": np.where(land, 2.0, 0.0),
        "terrain.continental_detail": np.where(land, 2.0, 0.0),
    })
    world.objects["tectonics.ocean_basins"] = [
        {
            "id": f"basin:{idx}",
            "kind": "ocean_basin_lifecycle",
            "area_fraction": frac,
            "cell_count": 100,
        }
        for idx, frac in enumerate((0.30, 0.24, 0.16, 0.05))
    ]
    world.objects["ocean.basins"] = [
        {
            "id": idx,
            "kind": "ocean_basin",
            "area_fraction": frac,
            "fraction_of_ocean": frac / 0.80,
            "cell_count": 90,
        }
        for idx, frac in enumerate((0.32, 0.22, 0.18))
    ]

    metrics = terminal_audit_metrics(world, [])
    basins = metrics["p110a_modern_planform"]["ocean_basins"]

    assert basins["source"] == "ocean_basin_objects"
    assert basins["basin_count"] == 3
    assert {row["source"] for row in basins["basin_rows_top10"]} == {"ocean.basins"}
    assert basins["major_basin_count"] >= 3
    assert basins["largest_basin_share_of_ocean"] < 0.50


def test_p110b_ocean_gateway_topology_metrics_distinguish_lifecycle_gateways():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    land = np.abs(grid.lat) > 58.0
    ocean = ~land
    gateway = ocean & (np.abs(grid.lon) < 7.5) & (np.abs(grid.lat) < 25.0)
    assert int(gateway.sum()) >= 3
    basin_id = np.full(grid.n, -1.0, dtype=np.float64)
    basin_id[ocean & (grid.lon < 0.0)] = 0.0
    basin_id[ocean & (grid.lon >= 0.0)] = 1.0
    gateway_id = np.full(grid.n, -1.0, dtype=np.float64)
    gateway_id[gateway] = 0.0
    depth_province = np.where(ocean, 4.0, 0.0)
    depth_province[gateway] = 7.0
    elev = np.where(land, 600.0, -4200.0)

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 12).astype(float),
        "terrain.elevation_m": elev.astype(float),
        "crust.type": land.astype(float),
        "crust.age_myr": np.where(land, 1600.0, 90.0),
        "ocean.depth_province": depth_province.astype(float),
        "ocean.basin_id": basin_id,
        "ocean.gateway_id": gateway_id,
        "ocean.gateway_system_id": gateway_id,
        "ocean.margin_type": np.where(gateway, 4.0, np.where(ocean, 5.0, 0.0)),
        "ocean.shelf_width": np.where(gateway, 1.0, 0.0),
        "terrain.province": np.where(land, 2.0, 0.0),
        "terrain.continental_detail": np.where(land, 2.0, 0.0),
    })
    ocean_area = float(grid.cell_area[ocean].sum())
    west = ocean & (basin_id == 0.0)
    east = ocean & (basin_id == 1.0)
    world.objects["ocean.basins"] = [
        {
            "id": 0,
            "type": "ocean_basin",
            "cell_count": int(west.sum()),
            "area_fraction": float(grid.cell_area[west].sum() / grid.cell_area.sum()),
            "fraction_of_ocean": float(grid.cell_area[west].sum() / ocean_area),
            "restricted_fraction": 0.04,
        },
        {
            "id": 1,
            "type": "ocean_basin",
            "cell_count": int(east.sum()),
            "area_fraction": float(grid.cell_area[east].sum() / grid.cell_area.sum()),
            "fraction_of_ocean": float(grid.cell_area[east].sum() / ocean_area),
            "restricted_fraction": 0.02,
        },
    ]
    world.objects["ocean.gateways"] = [
        {
            "id": 0,
            "type": "ocean_gateway",
            "cell_count": int(gateway.sum()),
            "area_fraction_of_ocean": float(grid.cell_area[gateway].sum() / ocean_area),
            "basin_ids": [0, 1],
            "wilson_phase_codes": [2, 4],
            "cells": np.where(gateway)[0].astype(int).tolist(),
        }
    ]
    world.objects["ocean.gateway_systems"] = [
        {
            "id": 0,
            "type": "ocean_gateway_system",
            "kind": "restricted_remnant_gateway",
            "fragment_count": 1,
            "cell_count": int(gateway.sum()),
            "area_fraction_of_ocean": float(grid.cell_area[gateway].sum() / ocean_area),
            "basin_ids": [0, 1],
            "interbasin": True,
            "wilson_phase_codes": [2, 4],
            "phase_backed": True,
            "restricted_fraction": 1.0,
            "cells": np.where(gateway)[0].astype(int).tolist(),
        }
    ]
    world.objects["tectonics.ocean_gateways"] = [
        {
            "id": "gateway:ridge:test",
            "kind": "ocean_gateway",
            "status": "opening",
            "phase_code": 2.0,
            "cells": np.where(gateway)[0].astype(int).tolist(),
        },
        {
            "id": "gateway:trench:test",
            "kind": "ocean_gateway",
            "status": "closing",
            "phase_code": 4.0,
            "cells": np.where(gateway)[0].astype(int).tolist(),
        },
    ]

    metrics = terminal_audit_metrics(world)
    gateway_metrics = metrics["p110a_modern_planform"]["ocean_gateway_topology"]

    assert metrics["object_count_by_set"]["ocean.basins"] == 2
    assert metrics["object_count_by_set"]["ocean.gateways"] == 1
    assert metrics["object_count_by_set"]["ocean.gateway_systems"] == 1
    assert gateway_metrics["terminal_gateway_count"] == 1
    assert gateway_metrics["gateway_field_component_count"] == 1
    assert gateway_metrics["terminal_interbasin_gateway_count"] == 1
    assert gateway_metrics["terminal_phase_backed_gateway_count"] == 1
    assert gateway_metrics["terminal_gateway_system_count"] == 1
    assert gateway_metrics["gateway_system_field_component_count"] == 1
    assert gateway_metrics["terminal_interbasin_gateway_system_count"] == 1
    assert gateway_metrics["terminal_phase_backed_gateway_system_count"] == 1
    assert gateway_metrics["gateway_fragment_to_system_ratio"] == 1.0
    assert gateway_metrics["gateway_system_kind_counts"] == {
        "restricted_remnant_gateway": 1
    }
    assert gateway_metrics["tectonic_gateway_count"] == 2
    assert gateway_metrics["tectonic_gateway_status_counts"] == {
        "closing": 1,
        "opening": 1,
    }
    assert gateway_metrics["restricted_ocean_fraction"] > 0.0
    assert gateway_metrics["unbacked_major_disconnected_ocean_component_count"] == 0


def test_p110b_ocean_gateway_systems_group_fragments_by_basin_and_phase():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    terrain = TerrainModule()

    def nearest(lat_deg: float, lon_deg: float, used: set[int]) -> int:
        lat = np.deg2rad(lat_deg)
        lon = np.deg2rad(lon_deg)
        xyz = np.array([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ])
        order = np.argsort(-(grid.xyz @ xyz))
        for cell in order:
            c = int(cell)
            if c not in used:
                used.add(c)
                return c
        raise AssertionError("no unique grid cell")

    used: set[int] = set()
    gateway_cells = [
        nearest(2.0, 0.0, used),
        nearest(8.0, 8.0, used),
        nearest(-42.0, 120.0, used),
        nearest(-49.0, 128.0, used),
    ]
    gateway_labels = np.full(grid.n, -1, dtype=int)
    for gid, cell in enumerate(gateway_cells):
        gateway_labels[cell] = gid
    gateway_area = np.asarray(
        [float(grid.cell_area[cell]) for cell in gateway_cells],
        dtype=np.float64,
    )
    basin_id = np.zeros(grid.n, dtype=int)
    tect_phase = np.full(grid.n, -1.0, dtype=np.float64)
    depth = np.full(grid.n, 4.0, dtype=np.float64)
    margin = np.full(grid.n, 5.0, dtype=np.float64)

    for cell in gateway_cells[:2]:
        basin_id[cell] = 0
        basin_id[int(grid.neighbors[cell][0])] = 1
        tect_phase[cell] = 2.0
        margin[cell] = 1.0
    for cell in gateway_cells[2:]:
        local = [cell, *[int(nb) for nb in grid.neighbors[cell]]]
        basin_id[local] = 2
        tect_phase[cell] = 5.0
        depth[cell] = 7.0
        margin[cell] = 4.0

    system_id, systems = terrain._ocean_gateway_system_objects(
        grid,
        gateway_labels,
        gateway_area,
        basin_id,
        tect_phase,
        depth,
        margin,
        float(grid.cell_area.sum()),
    )

    assert int(np.count_nonzero(system_id >= 0)) == 4
    assert len(systems) == 2
    assert sorted(system["fragment_count"] for system in systems) == [2, 2]
    kinds = {system["kind"] for system in systems}
    assert "open_seaway" in kinds
    assert "inland_relict_sea" in kinds


def test_p110a_open_ocean_partition_splits_dominant_ocean_with_tiny_seas():
    spec = get_preset("earthlike")
    spec.grid_cells = 1200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    def cap(lat_deg: float, lon_deg: float, radius_deg: float) -> np.ndarray:
        lat = np.deg2rad(lat_deg)
        lon = np.deg2rad(lon_deg)
        xyz = np.array([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ])
        return (grid.xyz @ xyz) >= np.cos(np.deg2rad(radius_deg))

    ocean = np.abs(grid.lat) < 48.0
    ocean |= cap(68.0, -35.0, 5.5)
    ocean |= cap(-68.0, 110.0, 5.5)
    land = ~ocean
    basin_id, basin_area = terrain._component_labels(grid, ocean)
    assert basin_area.size >= 3
    ocean_area = float(grid.cell_area[ocean].sum())
    assert float(np.max(basin_area) / ocean_area) > 0.90

    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=8)
    split_id, split_area = terrain._partition_open_ocean_basins(
        world, basin_id, basin_area, ocean, coast_pass)
    split_shares = split_area / max(float(split_area.sum()), 1.0e-12)

    assert split_area.size > basin_area.size
    assert int(np.count_nonzero(split_shares >= 0.06)) >= 3
    assert float(np.max(split_shares)) < 0.75
    assert np.all(split_id[land] < 0)


def test_p110a_component_balanced_payback_grows_secondary_continents():
    spec = get_preset("earthlike")
    spec.grid_cells = 2400
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    def cap(lat_deg: float, lon_deg: float, radius_deg: float) -> np.ndarray:
        lat = np.deg2rad(lat_deg)
        lon = np.deg2rad(lon_deg)
        xyz = np.array([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ])
        return (grid.xyz @ xyz) >= np.cos(np.deg2rad(radius_deg))

    land_a = cap(18.0, -125.0, 34.0)
    land_b = cap(-18.0, 45.0, 34.0)
    land_c = cap(35.0, 132.0, 12.0)
    land_d = cap(-48.0, -20.0, 8.0)
    shelf_c = cap(35.0, 132.0, 22.0)
    shelf_d = cap(-48.0, -20.0, 16.0)
    land = land_a | land_b | land_c | land_d
    shelf = (shelf_c | shelf_d) & ~land

    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[land] = 520.0
    surface[shelf] = -260.0
    crust = np.full(grid.n, OCEAN, dtype=np.float64)
    crust[land | shelf] = CONT
    domain = np.full(grid.n, DOMAIN_CONTINENTAL_MARGIN, dtype=np.float64)
    domain[land_a | land_b] = DOMAIN_CRATON
    domain[land_c | shelf_c] = DOMAIN_CONTINENTAL_INTERIOR
    domain[land_d | shelf_d] = DOMAIN_CONTINENTAL_MARGIN
    stability = np.full(grid.n, 0.45, dtype=np.float64)
    stability[land_a | land_b] = 0.90
    stability[land_c | shelf_c] = 0.72
    stability[land_d | shelf_d] = 0.62

    def shares(mask: np.ndarray) -> list[float]:
        comps = terrain._components(grid, mask)
        comps.sort(key=lambda c: float(grid.cell_area[c].sum()), reverse=True)
        total_land = max(float(grid.cell_area[mask].sum()), 1.0)
        return [float(grid.cell_area[c].sum()) / total_land for c in comps]

    before = shares(surface >= 0.0)
    assert len(before) >= 4
    assert int(sum(v >= 0.030 for v in before)) == 3
    assert int(sum(v >= 0.080 for v in before[1:])) == 1

    shaped = terrain._p110a_component_balanced_land_payback(
        world, surface, 0.0, crust, domain, stability)
    after_land = shaped >= 0.0
    after = shares(after_land)

    assert float(grid.cell_area[after_land].sum()) > float(grid.cell_area[land].sum())
    assert int(sum(v >= 0.030 for v in after)) >= 4
    assert int(sum(v >= 0.080 for v in after[1:])) >= 2
    assert after[0] <= before[0]
    assert world.g("terrain.last_p110a_component_payback_added_fraction") > 0.0
    assert world.g("terrain.last_p110a_component_payback_major_count_after") >= 4.0
    assert world.g("terrain.last_p110a_component_payback_nonlargest_large_after") >= 2.0


def test_p110a_final_modern_seaway_polish_splits_weak_supercontinent():
    spec = get_preset("earthlike")
    spec.grid_cells = 1800
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    def center_xyz(lat_deg: float, lon_deg: float) -> np.ndarray:
        lat = np.deg2rad(lat_deg)
        lon = np.deg2rad(lon_deg)
        return np.array([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ])

    land = np.zeros(grid.n, dtype=bool)
    for lon in (-105.0, -35.0, 35.0, 105.0):
        land |= (grid.xyz @ center_xyz(0.0, lon)) >= np.cos(np.deg2rad(28.0))
    weak_bridge = np.zeros(grid.n, dtype=bool)
    for lo, hi in ((-82.0, -58.0), (-12.0, 12.0), (58.0, 82.0)):
        bridge = (grid.lon > lo) & (grid.lon < hi) & (np.abs(grid.lat) < 11.0)
        weak_bridge |= bridge
        land |= bridge

    surface = np.where(land, 650.0, -4300.0).astype(np.float64)
    surface[weak_bridge] = 420.0
    crust_type = np.where(land, CONT, OCEAN).astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, 0.0).astype(np.float64)
    crust_stability = np.where(land, 0.60, 0.0).astype(np.float64)
    crust_domain[land & ~weak_bridge] = DOMAIN_CRATON
    crust_stability[land & ~weak_bridge] = 0.86
    crust_domain[weak_bridge] = DOMAIN_SUTURE
    crust_stability[weak_bridge] = 0.20

    before_comps = terrain._components(grid, surface >= 0.0)
    before_largest = max(
        float(grid.cell_area[c].sum()) for c in before_comps
    ) / float(grid.cell_area[surface >= 0.0].sum())

    shaped = terrain._p110a_final_modern_seaway_polish(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        max_passes=4,
    )
    after_land = shaped >= 0.0
    after_comps = terrain._components(grid, after_land)
    after_areas = sorted(
        (float(grid.cell_area[c].sum()) for c in after_comps),
        reverse=True,
    )
    after_largest = after_areas[0] / float(grid.cell_area[after_land].sum())

    assert before_largest > 0.99
    assert world.g("terrain.last_p110a_final_seaway_polish_passes") >= 1.0
    assert after_largest < 0.65
    assert len(after_comps) >= 2
    assert np.mean(shaped[weak_bridge] < 0.0) > 0.05
    endpoint_seaways = world.objects["terrain.p1114_modern_endpoint_seaways"]
    assert endpoint_seaways
    assert endpoint_seaways[0]["kind"] == "p1114_modern_endpoint_seaway"
    assert endpoint_seaways[0]["area_fraction"] > 0.0
    protected = terrain._breakup_seaway_protection_mask(
        world, grid, shaped < 0.0, passes=0)
    assert int(np.count_nonzero(protected)) >= int(endpoint_seaways[0]["cell_count"])


def test_p110a_final_modern_seaway_polish_uses_continent_domain_boundary():
    spec = get_preset("earthlike")
    spec.grid_cells = 1800
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    land = (
        (np.abs(grid.lat) < 46.0)
        & (grid.lon > -150.0)
        & (grid.lon < 150.0)
    )
    surface = np.where(land, 2200.0, -4300.0).astype(np.float64)
    crust_type = np.where(land, CONT, OCEAN).astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, 0.0).astype(np.float64)
    crust_stability = np.where(land, 0.70, 0.0).astype(np.float64)
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id[land & (grid.lon < 0.0)] = 101.0
    continent_id[land & (grid.lon >= 0.0)] = 202.0
    world.fields["tectonics.continent_id"] = continent_id

    before_comps = terrain._components(grid, surface >= 0.0)
    before_largest = max(
        float(grid.cell_area[c].sum()) for c in before_comps
    ) / float(grid.cell_area[surface >= 0.0].sum())

    shaped = terrain._p110a_final_modern_seaway_polish(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        max_passes=1,
    )
    after_land = shaped >= 0.0
    after_comps = terrain._components(grid, after_land)
    after_areas = sorted(
        (float(grid.cell_area[c].sum()) for c in after_comps),
        reverse=True,
    )
    after_largest = after_areas[0] / float(grid.cell_area[after_land].sum())

    assert before_largest > 0.99
    assert world.g("terrain.last_p110a_final_seaway_polish_passes") == 1.0
    assert len(after_comps) >= 2
    assert after_largest < 0.62
    assert np.count_nonzero((surface >= 0.0) & ~after_land) > 0
    endpoint_seaways = world.objects["terrain.p1114_modern_endpoint_seaways"]
    assert endpoint_seaways
    assert endpoint_seaways[0]["basis"] == "continent_domain_boundary"


def test_p1114_strict_endpoint_polish_continues_past_p110b_soft_stop():
    spec = get_preset("earthlike")
    spec.grid_cells = 2400
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    def center_xyz(lat_deg: float, lon_deg: float) -> np.ndarray:
        lat = np.deg2rad(lat_deg)
        lon = np.deg2rad(lon_deg)
        return np.array([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ])

    land = np.zeros(grid.n, dtype=bool)
    land |= (grid.xyz @ center_xyz(5.0, -95.0)) >= np.cos(np.deg2rad(32.0))
    land |= (grid.xyz @ center_xyz(5.0, -15.0)) >= np.cos(np.deg2rad(22.0))
    weak_bridge = (
        (grid.lon > -75.0)
        & (grid.lon < -35.0)
        & (np.abs(grid.lat - 5.0) < 6.0)
    )
    land |= weak_bridge
    for lat, lon, radius in (
        (15.0, 65.0, 25.0),
        (-30.0, 125.0, 20.0),
        (-35.0, -145.0, 18.0),
        (45.0, 155.0, 15.0),
    ):
        land |= (grid.xyz @ center_xyz(lat, lon)) >= np.cos(np.deg2rad(radius))

    surface = np.where(land, 700.0, -4300.0).astype(np.float64)
    surface[weak_bridge] = 250.0
    crust_type = np.where(land, CONT, OCEAN).astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, 0.0).astype(np.float64)
    crust_stability = np.where(land, 0.72, 0.0).astype(np.float64)
    crust_domain[land & ~weak_bridge] = DOMAIN_CRATON
    crust_stability[land & ~weak_bridge] = 0.88
    crust_domain[weak_bridge] = DOMAIN_SUTURE
    crust_stability[weak_bridge] = 0.20

    def shares(mask: np.ndarray) -> list[float]:
        comps = terrain._components(grid, mask)
        comps.sort(key=lambda c: float(grid.cell_area[c].sum()), reverse=True)
        total_land = max(float(grid.cell_area[mask].sum()), 1.0)
        return [float(grid.cell_area[c].sum()) / total_land for c in comps]

    before = shares(surface >= 0.0)
    planform = terrain._p110a_modern_planform_target(
        world, float(spec.target_land_fraction))
    assert 0.45 < before[0] < planform["p110b_largest_share_soft_ceiling"]
    assert int(sum(v >= 0.030 for v in before)) >= 4
    assert int(sum(v >= 0.080 for v in before[1:])) >= 2

    shaped = terrain._p110a_final_modern_seaway_polish(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        max_passes=4,
    )
    after_land = shaped >= 0.0
    after = shares(after_land)

    assert world.g("terrain.last_p1114_strict_endpoint_polish_attempted") == 1.0
    assert world.g("terrain.last_p1114_strict_endpoint_polish_applied") == 1.0
    assert world.g("terrain.last_p110a_final_seaway_polish_passes") >= 1.0
    assert after[0] >= 0.30
    assert after[0] < 0.45
    endpoint_seaways = world.objects["terrain.p1114_modern_endpoint_seaways"]
    assert endpoint_seaways
    assert endpoint_seaways[-1]["p111_strict_preferred"] is True
    assert endpoint_seaways[-1]["basis"] in {
        "continent_domain_boundary",
        "rift_suture_or_margin_weak_zone",
        "weak_lowland_bridge",
    }
    assert np.count_nonzero((surface >= 0.0) & ~after_land) > 0


def test_p110a_final_modern_seaway_polish_default_allows_six_passes():
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    surface = np.full(grid.n, 500.0, dtype=np.float64)
    crust_type = np.full(grid.n, CONT, dtype=np.float64)
    crust_domain = np.full(grid.n, DOMAIN_SUTURE, dtype=np.float64)
    crust_stability = np.full(grid.n, 0.20, dtype=np.float64)
    calls: list[int] = []

    def fake_open(world_arg, surface_arg, sea_level_arg, *_args, **_kwargs):
        out = np.asarray(surface_arg, dtype=np.float64).copy()
        idx = len(calls)
        out[idx] = float(sea_level_arg) - 100.0
        calls.append(idx)
        world_arg.set_g("terrain.last_modern_seaway_applied", 1.0)
        return out

    terrain._open_modern_earthlike_seaways = fake_open  # type: ignore[method-assign]

    shaped = terrain._p110a_final_modern_seaway_polish(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
    )

    assert world.g("terrain.last_p110a_final_seaway_polish_passes") == 6.0
    assert len(calls) == 6
    assert np.count_nonzero(shaped < 0.0) == 6


def test_p1115_open_ocean_shoal_deepen_preserves_object_backed_relief():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 55.0
    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    surface[land] = 800.0
    ocean = ~land
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=8)
    far_open = np.where(
        ocean
        & ((coast_pass < 0) | (coast_pass >= 5))
        & (grid.lat < -25.0)
    )[0]
    assert far_open.size >= 3
    object_backed = int(far_open[0])
    volcanic_only = int(far_open[1])
    unsupported = int(far_open[-1])
    surface[[object_backed, volcanic_only, unsupported]] = -1800.0

    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.continent_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.fields["crust.domain"][object_backed] = DOMAIN_LIP
    world.fields["crust.origin"][object_backed] = ORIGIN_PLUME_IMPACT
    world.fields["tectonics.volcanism_age_myr"][volcanic_only] = world.time_myr
    world.set_g("terrain.enable_p1115_object_backed_shoal_preservation", 1.0)

    shaped = terrain._deepen_modern_earthlike_open_ocean_shoals(
        world, surface, sea_level)

    assert shaped[object_backed] == surface[object_backed]
    assert shaped[volcanic_only] <= sea_level - 3200.0
    assert shaped[unsupported] <= sea_level - 3200.0
    assert (
        world.g("terrain.last_p1115_object_backed_open_ocean_shoal_area_fraction")
        > 0.0
    )
    assert (
        world.g("terrain.last_p1115_unsupported_open_ocean_shoal_deepened_fraction")
        > 0.0
    )


def test_p111625_abyssal_depth_calibration_deepens_slope_rise_and_far_ocean():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 58.0
    ocean = ~land
    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    surface[land] = 700.0
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=5)
    slope_cells = np.where(ocean & (coast_pass == 2) & (grid.lon < -30.0))[0]
    rise_cells = np.where(ocean & (coast_pass == 3) & (grid.lon < 30.0))[0]
    far_cells = np.where(
        ocean
        & ((coast_pass < 0) | (coast_pass >= 4))
        & (grid.lat < -25.0)
    )[0]
    assert slope_cells.size >= 1
    assert rise_cells.size >= 1
    assert far_cells.size >= 2
    slope_cell = int(slope_cells[0])
    rise_cell = int(rise_cells[0])
    far_cell = int(far_cells[-1])
    backed_cell = int(far_cells[0])
    surface[[slope_cell, rise_cell, far_cell, backed_cell]] = -180.0

    world.fields["crust.age_myr"] = np.full(grid.n, 100.0, dtype=np.float64)
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = land.astype(np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.continent_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.fields["crust.domain"][backed_cell] = DOMAIN_LIP
    world.fields["crust.origin"][backed_cell] = ORIGIN_PLUME_IMPACT
    world.set_g("terrain.enable_p1115_object_backed_shoal_preservation", 1.0)
    world.set_g("terrain._p1115_final_ocean_floor_expression_active", 1.0)

    shaped = terrain._deepen_modern_earthlike_open_ocean_shoals(
        world, surface, sea_level)
    world.set_g("terrain._p1115_final_ocean_floor_expression_active", 0.0)

    assert shaped[slope_cell] <= sea_level - 900.0
    assert shaped[rise_cell] <= sea_level - 2400.0
    assert shaped[far_cell] <= sea_level - 4900.0
    assert shaped[backed_cell] == surface[backed_cell]
    assert np.array_equal(shaped >= sea_level, surface >= sea_level)
    assert world.g("terrain.last_p111625_abyssal_depth_calibration_used") == 1.0
    assert (
        world.g("terrain.last_p111625_abyssal_depth_calibration_adjusted_area_fraction")
        > 0.0
    )
    assert world.g("terrain.last_p111625_ocean_mean_depth_after_m") > world.g(
        "terrain.last_p111625_ocean_mean_depth_before_m"
    )
    assert world.g("terrain.last_p111625_ocean_p50_depth_after_m") > world.g(
        "terrain.last_p111625_ocean_p50_depth_before_m"
    )


def test_p111626_ridge_flank_narrowing_preserves_axis_and_deepens_flanks():
    spec = get_preset("earthlike")
    spec.grid_cells = 1200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 62.0
    ocean = ~land
    ridge = ocean & (np.abs(grid.lat) < 9.0) & (np.abs(grid.lon) < 145.0)
    assert int(ridge.sum()) >= 20
    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    surface[land] = 800.0
    surface[ridge] = -2200.0
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(ridge)[0].astype(np.int64),
    }
    world.fields["crust.age_myr"] = np.full(grid.n, 90.0, dtype=np.float64)
    world.fields["crust.age_myr"][ridge] = 12.0
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = land.astype(np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.continent_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.set_g("terrain._p1115_final_ocean_floor_expression_active", 1.0)
    world.set_g("terrain.enable_p111626_ridge_flank_narrowing", 1.0)

    shaped = terrain._deepen_modern_earthlike_open_ocean_shoals(
        world, surface, sea_level)
    world.set_g("terrain._p1115_final_ocean_floor_expression_active", 0.0)
    depth = sea_level - shaped
    ridge_depth = depth[ridge]

    assert np.array_equal(shaped >= sea_level, surface >= sea_level)
    assert world.g("terrain.last_p111626_ridge_flank_narrowing_used") == 1.0
    assert world.g("terrain.last_p111626_ridge_flank_candidate_area_fraction") > 0.0
    assert world.g("terrain.last_p111626_ridge_flank_adjusted_area_fraction") > 0.0
    assert world.g("terrain.last_p111626_ridge_open_shoal_after_fraction") < world.g(
        "terrain.last_p111626_ridge_open_shoal_before_fraction"
    )
    assert float(np.min(ridge_depth)) < 2600.0
    assert float(np.percentile(ridge_depth, 80)) > 3000.0


def test_p111628_final_floor_deepens_ridge_shoulders_not_axis():
    spec = get_preset("earthlike")
    spec.grid_cells = 1200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 64.0
    ocean = ~land
    ridge = ocean & (np.abs(grid.lat) < 11.0) & (np.abs(grid.lon) < 150.0)
    axis_hint = ridge & (np.abs(grid.lat) < 2.0)
    shoulder = ridge & (np.abs(grid.lat) > 7.0)
    assert int(axis_hint.sum()) >= 4
    assert int(shoulder.sum()) >= 8

    surface = np.full(grid.n, -4300.0, dtype=np.float64)
    surface[land] = 700.0
    surface[ridge] = -2200.0
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(ridge)[0].astype(np.int64),
    }
    age = np.full(grid.n, 90.0, dtype=np.float64)
    age[axis_hint] = 2.0
    age[shoulder] = 50.0
    world.fields["crust.age_myr"] = age
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = land.astype(np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.continent_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.set_g("terrain._p1115_final_ocean_floor_expression_active", 1.0)

    shaped = terrain._deepen_modern_earthlike_open_ocean_shoals(
        world, surface, sea_level)
    world.set_g("terrain._p1115_final_ocean_floor_expression_active", 0.0)
    depth = sea_level - shaped

    assert np.array_equal(shaped >= sea_level, surface >= sea_level)
    assert world.g("terrain.last_p111628_final_shoal_floor_used") == 1.0
    assert world.g("terrain.last_p111628_ridge_shoulder_adjusted_area_fraction") > 0.0
    assert float(np.percentile(depth[axis_hint], 35)) < 2700.0
    assert float(np.percentile(depth[shoulder], 55)) > 3600.0


def test_p111628_final_semantic_cleanup_deepens_only_unsupported_shoals():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 62.0
    ocean = ~land
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=8)
    far_open = np.where(
        ocean
        & ((coast_pass < 0) | (coast_pass >= 5))
        & (np.abs(grid.lat) < 35.0)
    )[0]
    assert far_open.size >= 3

    supported_ridge = int(far_open[0])
    supported_chain = int(far_open[len(far_open) // 3])
    support = np.zeros(grid.n, dtype=bool)
    support[[supported_ridge, supported_chain]] = True
    support_halo = terrain._dilate_mask(grid, support, passes=1)
    far_open_mask = np.zeros(grid.n, dtype=bool)
    far_open_mask[far_open] = True
    halo_candidates = np.where(far_open_mask & support_halo & ~support)[0]
    assert halo_candidates.size > 0
    supported_halo = int(halo_candidates[0])
    unsupported_candidates = [int(c) for c in far_open if not support_halo[int(c)]]
    assert unsupported_candidates
    unsupported = unsupported_candidates[0]

    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    surface[land] = 800.0
    surface[[supported_ridge, supported_chain, supported_halo, unsupported]] = -2200.0
    world.fields["crust.age_myr"] = np.full(grid.n, 90.0, dtype=np.float64)
    world.fields["crust.age_myr"][unsupported] = 0.0
    world.fields["crust.age_myr"][supported_halo] = 0.0
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = land.astype(np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.set_g("terrain._p1115_final_ocean_floor_expression_active", 1.0)

    shaped = terrain._p111628_final_semantic_open_ocean_shoal_cleanup(
        world,
        surface,
        sea_level,
        ocean_fabric=[
            {"kind": "spreading_center", "cells": [supported_ridge]},
        ],
        arc_plume_landforms=[
            {"kind": "seamount_chain", "cells": [supported_chain]},
        ],
        margin_landforms=[],
    )
    world.set_g("terrain._p1115_final_ocean_floor_expression_active", 0.0)

    assert np.array_equal(shaped >= sea_level, surface >= sea_level)
    assert shaped[supported_ridge] == surface[supported_ridge]
    assert shaped[supported_chain] == surface[supported_chain]
    assert shaped[supported_halo] <= sea_level - 2650.0
    assert shaped[unsupported] <= sea_level - 2850.0
    assert world.g("terrain.last_p111628_final_semantic_cleanup_used") == 1.0
    assert (
        world.g(
            "terrain.last_p111628_final_semantic_cleanup_adjusted_area_fraction"
        )
        > 0.0
    )
    assert world.g("terrain.last_p111630_object_halo_shoal_cleanup_used") == 1.0
    assert (
        world.g(
            "terrain.last_p111630_object_halo_shoal_cleanup_adjusted_area_fraction"
        )
        > 0.0
    )


def test_p173_terminal_process_island_promotion_uses_object_backed_ocean_highs():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 66.0
    ocean = ~land
    surface = np.full(grid.n, -4300.0, dtype=np.float64)
    surface[land] = 700.0
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=5)
    far_open = np.where(
        ocean
        & ((coast_pass < 0) | (coast_pass >= 3))
        & (np.abs(grid.lat) < 50.0)
    )[0]
    assert far_open.size >= 8
    far_mask = np.zeros(grid.n, dtype=bool)
    far_mask[far_open] = True
    backed = None
    for candidate in far_open:
        nbs = [
            int(nb)
            for nb in grid.neighbors[int(candidate)]
            if far_mask[int(nb)]
        ]
        if len(nbs) >= 2:
            backed = np.asarray([int(candidate), *nbs[:4]], dtype=np.int64)
            break
    assert backed is not None
    backed_set = set(int(c) for c in backed.tolist())
    shoulder = []
    for c in backed:
        for nb in grid.neighbors[int(c)]:
            nb = int(nb)
            if far_mask[nb] and nb not in backed_set:
                shoulder.append(nb)
    shoulder = np.asarray(sorted(set(shoulder))[:8], dtype=np.int64)
    assert shoulder.size >= 1
    support_set = backed_set | set(int(c) for c in shoulder.tolist())
    outer = []
    for c in shoulder:
        for nb in grid.neighbors[int(c)]:
            nb = int(nb)
            if far_mask[nb] and nb not in support_set:
                outer.append(nb)
    outer = np.asarray(sorted(set(outer))[:8], dtype=np.int64)
    support_set |= set(int(c) for c in outer.tolist())
    unsupported = next(int(c) for c in far_open if int(c) not in support_set)
    surface[backed] = -760.0
    surface[shoulder] = -1700.0
    surface[outer] = -2400.0
    surface[unsupported] = -700.0

    shaped = terrain._p173_terminal_process_island_promotion(
        world,
        surface,
        sea_level,
        arc_plume_landforms=[
            {
                "kind": "seamount_chain",
                "cells": np.concatenate([backed, shoulder, outer]).astype(int).tolist(),
            },
        ],
    )

    promoted = (shaped >= sea_level) & (surface < sea_level)
    assert np.count_nonzero(promoted & np.isin(np.arange(grid.n), backed)) >= 1
    assert shaped[unsupported] == surface[unsupported]
    assert world.g("terrain.last_p173_terminal_process_island_added_cells") >= 1.0
    assert world.g("terrain.last_p173_terminal_process_island_shelf_cells") >= 1.0
    assert world.g("terrain.last_p173_terminal_process_island_reef_cells") >= 0.0
    assert world.objects["terrain.p173_terminal_process_islands"]


def test_p11153_object_backed_shoal_preservation_does_not_keep_halo():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 55.0
    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    surface[land] = 800.0
    ocean = ~land
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=8)
    far_open = ocean & ((coast_pass < 0) | (coast_pass >= 5)) & (grid.lat < -25.0)
    core = None
    halo = None
    for candidate in np.where(far_open)[0]:
        nbs = [int(nb) for nb in grid.neighbors[int(candidate)] if far_open[int(nb)]]
        if nbs:
            core = int(candidate)
            halo = nbs[0]
            break
    assert core is not None
    assert halo is not None
    surface[[core, halo]] = -1800.0

    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.fields["crust.domain"][core] = DOMAIN_LIP
    world.fields["crust.origin"][core] = ORIGIN_PLUME_IMPACT
    world.set_g("terrain.enable_p1115_object_backed_shoal_preservation", 1.0)

    shaped = terrain._deepen_modern_earthlike_open_ocean_shoals(
        world, surface, sea_level)

    assert shaped[core] == surface[core]
    assert shaped[halo] <= sea_level - 3200.0


def test_p11153_broad_volcanism_only_affects_preterminal_ocean_fabric():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 65.0
    ocean = ~land
    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    surface[land] = 700.0
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=8)
    candidates = np.where(
        ocean
        & ((coast_pass < 0) | (coast_pass >= 5))
        & (grid.lat < -20.0)
    )[0]
    assert candidates.size
    volcanic = int(candidates[0])

    world.fields["crust.age_myr"] = np.full(grid.n, 80.0, dtype=np.float64)
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"][volcanic] = world.time_myr - 500.0
    empty = np.zeros(grid.n, dtype=bool)

    preterminal = terrain._apply_coherent_ocean_floor_fabric(
        world,
        surface,
        sea_level,
        ocean,
        coast_pass,
        empty,
        empty,
        empty,
        p1115_object_relief=False,
    )
    final = terrain._apply_coherent_ocean_floor_fabric(
        world,
        surface,
        sea_level,
        ocean,
        coast_pass,
        empty,
        empty,
        empty,
        p1115_object_relief=True,
    )

    assert preterminal[volcanic] > surface[volcanic]
    assert final[volcanic] == surface[volcanic]


def test_p1116_final_open_ridge_cleanup_keeps_legible_mid_ocean_high():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 60.0
    ocean = ~land
    ridge = ocean & (np.abs(grid.lat) < 4.0)
    assert ridge.any()
    ridge_cell = int(np.where(ridge)[0][0])
    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    surface[land] = 800.0
    surface[ridge] = -2100.0
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(ridge)[0].astype(np.int64),
    }
    world.fields["crust.age_myr"] = np.where(ridge, 2.0, 80.0)
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = land.astype(float)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)

    preterminal = terrain._deepen_modern_earthlike_open_ocean_shoals(
        world, surface, sea_level)
    world.set_g("terrain._p1115_final_ocean_floor_expression_active", 1.0)
    final = terrain._deepen_modern_earthlike_open_ocean_shoals(
        world, surface, sea_level)
    world.set_g("terrain._p1115_final_ocean_floor_expression_active", 0.0)

    assert preterminal[ridge_cell] <= sea_level - 2600.0
    assert final[ridge_cell] == sea_level - 2200.0


def test_p1116_final_ocean_fabric_deepens_old_oceanic_crust():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 70.0
    ocean = ~land
    surface = np.full(grid.n, -2500.0, dtype=np.float64)
    surface[land] = 800.0
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=8)
    young = ocean & (grid.lat > 20.0) & (grid.lat < 45.0)
    old = ocean & (grid.lat < -35.0)
    assert young.any()
    assert old.any()
    age = np.full(grid.n, 90.0, dtype=np.float64)
    age[young] = 8.0
    age[old] = 180.0
    world.fields["crust.age_myr"] = age
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = land.astype(float)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.set_g("terrain.enable_p1116_age_depth_final_expression", 1.0)
    empty = np.zeros(grid.n, dtype=bool)

    shaped = terrain._apply_coherent_ocean_floor_fabric(
        world,
        surface,
        sea_level,
        ocean,
        coast_pass,
        empty,
        empty,
        empty,
        p1115_object_relief=True,
    )

    assert np.median(sea_level - shaped[old]) > np.median(sea_level - shaped[young])
    assert np.median(sea_level - shaped[old]) >= 3300.0


def test_p1115_final_ocean_floor_expression_defaults_to_endpoint_guardrails():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 55.0
    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    surface[land] = 900.0
    ocean = ~land
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=8)
    far_open = np.where(
        ocean
        & ((coast_pass < 0) | (coast_pass >= 5))
        & (grid.lat < -25.0)
    )[0]
    assert far_open.size >= 3
    plateau_cell = int(far_open[0])
    unsupported_shoal = int(far_open[1])
    deep_cell = int(far_open[-1])
    surface[unsupported_shoal] = -1800.0
    surface[deep_cell] = -6500.0

    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.continent_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.fields["crust.age_myr"] = np.full(grid.n, 80.0, dtype=np.float64)
    world.fields["crust.domain"][plateau_cell] = DOMAIN_LIP
    world.fields["crust.origin"][plateau_cell] = ORIGIN_PLUME_IMPACT

    p1115_profile = {}
    shaped = terrain._p1115_final_ocean_floor_hierarchy_expression(
        world, surface, sea_level, profile_stage_seconds=p1115_profile)

    assert np.array_equal(shaped >= sea_level, surface >= sea_level)
    assert shaped[plateau_cell] == surface[plateau_cell]
    assert shaped[unsupported_shoal] <= sea_level - 2850.0
    assert shaped[deep_cell] == surface[deep_cell]
    assert world.g("terrain.last_p1115_ocean_floor_land_mask_preserved") == 1.0
    assert world.g("terrain.last_p1115_ocean_floor_adjusted_area_fraction") > 0.0
    assert world.g("terrain.last_p1115_endpoint_readability_guardrail_mode") == 1.0
    assert world.g("terrain.last_p1115_process_relief_enabled") == 0.0
    assert world.g("terrain.last_p1115_process_relief_adjusted_area_fraction") == 0.0
    assert p1115_profile["boundary_zone_derivation"] >= 0.0
    assert p1115_profile["endpoint_semantic_shoal_cleanup"] >= 0.0
    assert p1115_profile["land_mask_guard"] >= 0.0


def test_p1115_legacy_ocean_floor_expression_requires_explicit_opt_out():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 55.0
    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    surface[land] = 900.0
    ocean = ~land
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=8)
    far_open = np.where(
        ocean
        & ((coast_pass < 0) | (coast_pass >= 5))
        & (grid.lat < -25.0)
    )[0]
    assert far_open.size >= 2
    plateau_cell = int(far_open[0])
    deep_cell = int(far_open[-1])
    surface[deep_cell] = -6500.0

    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.continent_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.fields["crust.age_myr"] = np.full(grid.n, 80.0, dtype=np.float64)
    world.fields["crust.domain"][plateau_cell] = DOMAIN_LIP
    world.fields["crust.origin"][plateau_cell] = ORIGIN_PLUME_IMPACT
    world.set_g("terrain.enable_p1115_endpoint_readability_guardrails", 0.0)

    p1115_profile = {}
    shaped = terrain._p1115_final_ocean_floor_hierarchy_expression(
        world, surface, sea_level, profile_stage_seconds=p1115_profile)

    assert np.array_equal(shaped >= sea_level, surface >= sea_level)
    assert shaped[plateau_cell] > surface[plateau_cell]
    assert shaped[deep_cell] == surface[deep_cell]
    assert world.g("terrain.last_p1115_endpoint_readability_guardrail_mode") == 0.0
    assert world.g("terrain.last_p1115_process_relief_enabled") == 1.0
    assert world.g("terrain.last_p1115_process_relief_adjusted_area_fraction") > 0.0
    assert p1115_profile["coherent_ocean_floor_fabric"] >= 0.0
    assert p1115_profile["land_mask_guard"] >= 0.0


def test_p11163_final_trench_objects_drive_narrow_bathymetry_and_depth_province():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 58.0
    ocean = ~land
    surface = np.full(grid.n, -3200.0, dtype=np.float64)
    surface[land] = 900.0
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=8)
    trench_line = (
        ocean
        & (np.abs(grid.lat + 28.0) < 4.0)
        & (np.abs(grid.lon) < 105.0)
        & ((coast_pass < 0) | (coast_pass >= 4))
    )
    trench_line = terrain._generalized_boundary_source(
        grid, trench_line, spacing=2) & trench_line
    assert int(trench_line.sum()) >= 4
    margin_landforms = [
        {"kind": "trench", "cells": np.where(trench_line)[0].astype(int).tolist()}
    ]

    world.fields["crust.domain"] = np.full(grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.continent_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.fields["crust.age_myr"] = np.full(grid.n, 90.0, dtype=np.float64)

    shaped = terrain._p1115_final_ocean_floor_hierarchy_expression(
        world, surface, sea_level, margin_landforms=margin_landforms)
    final_axis = terrain._p11163_final_trench_axis(
        world, shaped < sea_level, coast_pass, margin_landforms)
    final_depth = sea_level - shaped

    assert np.array_equal(shaped >= sea_level, surface >= sea_level)
    assert final_axis.any()
    assert np.median(final_depth[final_axis]) >= 5200.0
    assert (
        world.g("terrain.last_p11163_trench_axis_area_fraction_ocean")
        < 0.095
    )

    ocean_geo_fields = {
        "ocean.depth_province": np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(float),
        "ocean.margin_type": np.zeros(grid.n, dtype=float),
    }
    patched = terrain._p11163_mark_final_trench_depth_province(
        world, shaped, sea_level, ocean_geo_fields, margin_landforms)

    assert np.all(
        patched["ocean.depth_province"][final_axis] == float(OCEAN_DEPTH_TRENCH)
    )


def test_p11152_ocean_feature_metrics_do_not_treat_ridge_province_as_object_backed():
    spec = get_preset("earthlike")
    spec.grid_cells = 360
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    ocean = grid.lat < 65.0
    depth = np.where(ocean, 2000.0, 0.0)
    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.int32)
    broad_ridge_province = ocean & (np.abs(grid.lat) < 26.0)
    depth_province[broad_ridge_province] = OCEAN_DEPTH_RIDGE
    world.fields["ocean.shelf_width"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)

    metrics = _ocean_feature_metrics(world, depth, ocean, depth_province)

    assert metrics["p1115_ridge_province_area_fraction_ocean"] > 0.25
    assert metrics["p1115_ridge_area_fraction_ocean"] == 0.0
    assert metrics["p1115_object_backed_open_ocean_shoal_fraction"] == 0.0
    assert metrics["p1115_unsupported_open_ocean_shoal_fraction"] > 0.90


def test_p11152_ocean_fabric_objects_thin_broad_ridge_province_to_axis():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 70.0
    ocean = ~land
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[land] = 600.0
    ridge_province = ocean & (np.abs(grid.lat) < 24.0)
    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.float64)
    depth_province[ridge_province] = OCEAN_DEPTH_RIDGE
    age = np.full(grid.n, 80.0, dtype=np.float64)
    age[ocean & (np.abs(grid.lat) < 4.0)] = 2.0
    world.fields["crust.age_myr"] = age
    world.networks["tectonics.boundaries"] = {}
    ocean_geo_fields = {
        "ocean.basin_id": np.zeros(grid.n, dtype=np.float64),
        "ocean.depth_province": depth_province,
        "ocean.shelf_width": np.where(ocean, 5.0, 0.0),
    }

    objects = terrain._ocean_fabric_objects(
        world,
        surface,
        sea_level,
        np.zeros(grid.n, dtype=np.float64),
        np.full(grid.n, 600.0, dtype=np.float64),
        ocean_geo_fields,
    )
    spreading = np.zeros(grid.n, dtype=bool)
    for obj in objects:
        if obj.get("kind") == "spreading_center":
            cells = np.asarray(obj.get("cells", []), dtype=np.int64)
            spreading[cells] = True

    assert spreading.any()
    ridge_area = float(grid.cell_area[ridge_province].sum())
    spreading_area = float(grid.cell_area[spreading].sum())
    assert spreading_area < 0.55 * ridge_area


def test_p111627_ocean_geography_narrows_broad_ridge_depth_province():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 72.0
    ocean = ~land
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[land] = 650.0
    broad_ridge = ocean & (np.abs(grid.lat) < 24.0)
    age = np.full(grid.n, 90.0, dtype=np.float64)
    age[ocean & (np.abs(grid.lat) < 3.5)] = 2.0
    world.fields["crust.age_myr"] = age
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(broad_ridge)[0].astype(np.int64),
    }

    fields, _objects, _diag = terrain._ocean_geography(world, surface, sea_level)
    ridge_depth = fields["ocean.depth_province"].astype(int) == OCEAN_DEPTH_RIDGE
    ocean_area = float(grid.cell_area[ocean].sum())
    raw_area = float(grid.cell_area[broad_ridge].sum() / ocean_area)
    ridge_area = float(grid.cell_area[ridge_depth].sum() / ocean_area)

    assert ridge_depth.any()
    assert ridge_area < 0.72 * raw_area
    assert world.g("terrain.last_p111627_depth_ridge_axis_area_fraction") > 0.0
    assert (
        world.g("terrain.last_p111627_depth_ridge_province_area_fraction_after")
        < world.g("terrain.last_p111627_depth_ridge_province_area_fraction_before")
    )


def test_p111627_ocean_fabric_uses_axis_instead_of_broad_raw_ridge_source():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 72.0
    ocean = ~land
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[land] = 650.0
    broad_ridge = ocean & (np.abs(grid.lat) < 24.0)
    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.float64)
    depth_province[broad_ridge] = OCEAN_DEPTH_RIDGE
    age = np.full(grid.n, 90.0, dtype=np.float64)
    age[ocean & (np.abs(grid.lat) < 3.5)] = 2.0
    world.fields["crust.age_myr"] = age
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(broad_ridge)[0].astype(np.int64),
    }
    ocean_geo_fields = {
        "ocean.basin_id": np.zeros(grid.n, dtype=np.float64),
        "ocean.depth_province": depth_province,
        "ocean.shelf_width": np.where(ocean, 5.0, 0.0),
    }

    objects = terrain._ocean_fabric_objects(
        world,
        surface,
        sea_level,
        np.zeros(grid.n, dtype=np.float64),
        np.full(grid.n, 600.0, dtype=np.float64),
        ocean_geo_fields,
    )
    spreading = np.zeros(grid.n, dtype=bool)
    for obj in objects:
        if obj.get("kind") == "spreading_center":
            cells = np.asarray(obj.get("cells", []), dtype=np.int64)
            spreading[cells] = True

    ocean_area = float(grid.cell_area[ocean].sum())
    raw_area = float(grid.cell_area[broad_ridge].sum() / ocean_area)
    spreading_area = float(grid.cell_area[spreading].sum() / ocean_area)
    assert spreading.any()
    assert spreading_area < 0.50 * raw_area
    assert (
        world.g("terrain.last_p111627_spreading_center_axis_area_fraction")
        < world.g("terrain.last_p111627_spreading_center_source_area_fraction")
    )


def test_p111627_ocean_feature_metrics_prefer_semantic_ridge_over_broad_raw_boundary():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    ocean = grid.lat < 72.0
    raw_ridge = ocean & (np.abs(grid.lat) < 24.0)
    semantic_ridge = ocean & (np.abs(grid.lat) < 4.0)
    depth = np.where(ocean, 3200.0, 0.0)
    depth[semantic_ridge] = 2300.0
    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.int32)
    depth_province[semantic_ridge] = OCEAN_DEPTH_RIDGE
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(raw_ridge)[0].astype(np.int64),
    }
    world.objects["terrain.ocean_fabric"] = [
        {
            "kind": "spreading_center",
            "cells": np.where(semantic_ridge)[0].astype(int).tolist(),
        }
    ]
    world.fields["ocean.shelf_width"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)

    metrics = _ocean_feature_metrics(world, depth, ocean, depth_province)
    raw_area = float(grid.cell_area[raw_ridge].sum() / grid.cell_area[ocean].sum())
    semantic_area = float(
        grid.cell_area[semantic_ridge].sum() / grid.cell_area[ocean].sum()
    )

    assert metrics["p1115_ridge_area_fraction_ocean"] < 0.45 * raw_area
    assert abs(metrics["p1115_ridge_area_fraction_ocean"] - semantic_area) < 0.02


def test_p111628_ocean_feature_metrics_attribute_shoal_sources():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    ocean = grid.lat < 72.0
    depth = np.where(ocean, 4200.0, 0.0)
    back_arc = ocean & (np.abs(grid.lat - 10.0) < 5.0) & (np.abs(grid.lon) < 40.0)
    plain = ocean & (np.abs(grid.lat + 20.0) < 5.0) & (np.abs(grid.lon - 80.0) < 35.0)
    hill = ocean & (np.abs(grid.lat + 36.0) < 5.0) & (np.abs(grid.lon + 80.0) < 35.0)
    microcontinent = (
        ocean
        & (np.abs(grid.lat - 36.0) < 6.0)
        & (np.abs(grid.lon + 140.0) < 35.0)
    )
    oceanic_plateau = (
        ocean
        & (np.abs(grid.lat + 2.0) < 6.0)
        & (np.abs(grid.lon - 140.0) < 35.0)
    )
    seamount_chain = (
        ocean
        & (np.abs(grid.lat - 50.0) < 6.0)
        & (np.abs(grid.lon - 80.0) < 35.0)
    )
    island_arc = (
        ocean
        & (np.abs(grid.lat + 52.0) < 6.0)
        & (np.abs(grid.lon) < 35.0)
    )
    assert back_arc.any() and plain.any() and hill.any()
    assert (
        microcontinent.any()
        and oceanic_plateau.any()
        and seamount_chain.any()
        and island_arc.any()
    )
    depth[
        back_arc
        | plain
        | hill
        | microcontinent
        | oceanic_plateau
        | seamount_chain
        | island_arc
    ] = 1800.0
    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.int32)
    world.fields["ocean.shelf_width"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.objects["terrain.arc_plume_landforms"] = [
        {"kind": "back_arc_basin", "cells": np.where(back_arc)[0].astype(int).tolist()},
        {
            "kind": "microcontinent",
            "cells": np.where(microcontinent)[0].astype(int).tolist(),
        },
        {
            "kind": "oceanic_plateau",
            "cells": np.where(oceanic_plateau)[0].astype(int).tolist(),
        },
        {
            "kind": "seamount_chain",
            "cells": np.where(seamount_chain)[0].astype(int).tolist(),
        },
        {"kind": "island_arc", "cells": np.where(island_arc)[0].astype(int).tolist()},
    ]
    world.objects["terrain.ocean_fabric"] = [
        {"kind": "abyssal_plain", "cells": np.where(plain)[0].astype(int).tolist()},
        {"kind": "abyssal_hill", "cells": np.where(hill)[0].astype(int).tolist()},
    ]
    world.set_g("terrain.last_p111628_final_semantic_cleanup_used", 1.0)
    world.set_g(
        "terrain.last_p111628_final_semantic_cleanup_candidate_area_fraction",
        0.012,
    )
    world.set_g(
        "terrain.last_p111628_final_semantic_cleanup_adjusted_area_fraction",
        0.011,
    )
    world.set_g(
        "terrain.last_p111628_final_semantic_cleanup_ridge_area_fraction",
        0.004,
    )
    world.set_g(
        "terrain.last_p111628_final_semantic_cleanup_abyss_area_fraction",
        0.007,
    )
    world.set_g("terrain.last_p111630_object_halo_shoal_cleanup_used", 1.0)
    world.set_g(
        "terrain.last_p111630_object_halo_shoal_cleanup_candidate_area_fraction",
        0.009,
    )
    world.set_g(
        "terrain.last_p111630_object_halo_shoal_cleanup_adjusted_area_fraction",
        0.008,
    )

    metrics = _ocean_feature_metrics(world, depth, ocean, depth_province)

    assert metrics["p111628_final_semantic_cleanup_used"]
    assert metrics["p111628_final_semantic_cleanup_candidate_area_fraction"] == 0.012
    assert metrics["p111628_final_semantic_cleanup_adjusted_area_fraction"] == 0.011
    assert metrics["p111628_final_semantic_cleanup_ridge_area_fraction"] == 0.004
    assert metrics["p111628_final_semantic_cleanup_abyss_area_fraction"] == 0.007
    assert metrics["p111630_object_halo_shoal_cleanup_used"]
    assert metrics["p111630_object_halo_shoal_cleanup_candidate_area_fraction"] == 0.009
    assert metrics["p111630_object_halo_shoal_cleanup_adjusted_area_fraction"] == 0.008
    assert metrics["p111630_object_backed_shoal_narrow_core_fraction"] > 0.0
    assert metrics["p111630_object_backed_shoal_broad_core_fraction"] > 0.0
    assert metrics["p111630_shoal_microcontinent_fraction"] > 0.0
    assert metrics["p111630_shoal_oceanic_plateau_fraction"] > 0.0
    assert metrics["p111630_shoal_seamount_chain_fraction"] > 0.0
    assert metrics["p111630_shoal_island_arc_fraction"] > 0.0
    assert metrics["p111628_shoal_back_arc_basin_fraction"] > 0.0
    assert metrics["p111628_unsupported_shoal_back_arc_basin_fraction"] > 0.0
    assert metrics["p111628_shoal_abyssal_plain_fraction"] > 0.0
    assert metrics["p111628_unsupported_shoal_abyssal_plain_fraction"] > 0.0
    assert metrics["p111628_shoal_abyssal_hill_fraction"] > 0.0
    assert metrics["p111628_unsupported_shoal_abyssal_hill_fraction"] > 0.0


def test_p11152_ocean_fabric_objects_fallback_old_abyssal_plain():
    spec = get_preset("earthlike")
    spec.grid_cells = 420
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 72.0
    ocean = ~land
    surface = np.full(grid.n, -4300.0, dtype=np.float64)
    surface[land] = 500.0
    world.fields["crust.age_myr"] = np.where(ocean, 105.0, 1600.0)
    world.networks["tectonics.boundaries"] = {}
    ocean_geo_fields = {
        "ocean.basin_id": np.zeros(grid.n, dtype=np.float64),
        "ocean.depth_province": np.where(
            ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.float64),
        "ocean.shelf_width": np.where(ocean, 5.0, 0.0),
    }

    objects = terrain._ocean_fabric_objects(
        world,
        surface,
        sea_level,
        np.zeros(grid.n, dtype=np.float64),
        np.full(grid.n, 300.0, dtype=np.float64),
        ocean_geo_fields,
    )

    assert any(obj.get("kind") == "abyssal_plain" for obj in objects)


def test_p11152_microcontinent_requires_more_than_continent_id():
    spec = get_preset("earthlike")
    spec.grid_cells = 360
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    surface = np.full(grid.n, -1200.0, dtype=np.float64)
    crust_type = np.ones(grid.n, dtype=np.float64)
    crust_domain = np.full(grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    crust_origin = np.zeros(grid.n, dtype=np.float64)
    crust_thick = np.full(grid.n, CONT_THICK, dtype=np.float64)
    sediment = np.full(grid.n, 500.0, dtype=np.float64)
    world.fields["crust.age_myr"] = np.full(grid.n, 90.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.continent_id"] = np.ones(grid.n, dtype=np.float64)
    world.fields["tectonics.rift_potential"] = np.zeros(grid.n, dtype=np.float64)
    ocean_geo_fields = {
        "ocean.depth_province": np.full(
            grid.n, OCEAN_DEPTH_ABYSS, dtype=np.float64),
    }

    objects = terrain._arc_plume_landform_objects(
        world,
        surface,
        sea_level,
        crust_type,
        crust_domain,
        crust_origin,
        crust_thick,
        sediment,
        ocean_geo_fields,
    )
    assert not any(obj.get("kind") == "microcontinent" for obj in objects)

    terrane_cells = np.where((np.abs(grid.lat) < 12.0) & (np.abs(grid.lon) < 24.0))[0]
    assert terrane_cells.size >= 2
    world.fields["tectonics.terrane_id"][terrane_cells] = 3.0
    objects = terrain._arc_plume_landform_objects(
        world,
        surface,
        sea_level,
        crust_type,
        crust_domain,
        crust_origin,
        crust_thick,
        sediment,
        ocean_geo_fields,
    )
    assert any(obj.get("kind") == "microcontinent" for obj in objects)


def test_p111631_microcontinent_footprint_cap_limits_broad_candidates():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    surface = np.full(grid.n, -1200.0, dtype=np.float64)
    broad = (np.abs(grid.lat) < 45.0) & (np.abs(grid.lon) < 120.0)
    crust_type = np.zeros(grid.n, dtype=np.float64)
    crust_type[broad] = 1.0
    crust_domain = np.full(grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    crust_domain[broad] = DOMAIN_ACCRETED_TERRANE
    crust_origin = np.zeros(grid.n, dtype=np.float64)
    crust_thick = np.full(grid.n, CONT_THICK, dtype=np.float64)
    sediment = np.full(grid.n, 500.0, dtype=np.float64)
    world.fields["crust.age_myr"] = np.full(grid.n, 90.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.terrane_id"][broad] = 7.0
    world.fields["tectonics.continent_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.rift_potential"] = np.zeros(grid.n, dtype=np.float64)
    ocean_geo_fields = {
        "ocean.depth_province": np.full(
            grid.n, OCEAN_DEPTH_ABYSS, dtype=np.float64),
    }

    objects = terrain._arc_plume_landform_objects(
        world,
        surface,
        sea_level,
        crust_type,
        crust_domain,
        crust_origin,
        crust_thick,
        sediment,
        ocean_geo_fields,
    )
    micro_objects = [obj for obj in objects if obj.get("kind") == "microcontinent"]
    ocean_area = float(grid.cell_area.sum())
    object_area = sum(
        float(grid.cell_area[np.asarray(obj["cells"], dtype=int)].sum())
        for obj in micro_objects
    ) / ocean_area

    assert micro_objects
    assert world.g("terrain.last_p111631_microcontinent_raw_area_fraction_ocean") > 0.10
    assert world.g("terrain.last_p111631_microcontinent_area_fraction_ocean") <= 0.021
    assert world.g("terrain.last_p111631_microcontinent_removed_area_fraction_ocean") > 0.08
    assert world.g("terrain.last_p111631_microcontinent_component_count") <= 18.0
    assert object_area <= 0.021


def test_p111631_microcontinent_final_relief_does_not_shoal_entire_candidate():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 72.0
    ocean = ~land
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[land] = 700.0
    broad = (
        ocean
        & (np.abs(grid.lat) < 35.0)
        & (np.abs(grid.lon) < 120.0)
    )
    surface[broad] = -1700.0
    coast_pass = terrain._ocean_pass_distance(grid, land, ocean, max_passes=5)

    world.fields["crust.age_myr"] = np.full(grid.n, 90.0, dtype=np.float64)
    world.fields["crust.type"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["crust.type"][broad] = 1.0
    world.fields["crust.domain"] = np.full(
        grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    world.fields["crust.domain"][broad] = DOMAIN_ACCRETED_TERRANE
    world.fields["crust.origin"] = np.zeros(grid.n, dtype=np.float64)
    world.fields["tectonics.terrane_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.terrane_id"][broad] = 4.0
    world.fields["tectonics.continent_id"] = np.full(grid.n, -1.0, dtype=np.float64)
    world.fields["tectonics.volcanism_age_myr"] = np.full(
        grid.n, -1.0, dtype=np.float64)

    shaped = terrain._apply_coherent_ocean_floor_fabric(
        world,
        surface,
        sea_level,
        ocean,
        coast_pass,
        np.zeros(grid.n, dtype=bool),
        np.zeros(grid.n, dtype=bool),
        np.zeros(grid.n, dtype=bool),
        p1115_object_relief=True,
    )
    depth = sea_level - shaped
    ocean_area = max(float(grid.cell_area[ocean].sum()), 1.0e-12)
    shallow_raw_fraction = float(
        grid.cell_area[broad & (depth > 0.0) & (depth < 2600.0)].sum()
        / ocean_area
    )

    assert world.g("terrain.last_p111631_microcontinent_relief_raw_area_fraction_ocean") > 0.10
    assert world.g("terrain.last_p111631_microcontinent_relief_area_fraction_ocean") <= 0.021
    assert world.g("terrain.last_p111631_microcontinent_relief_removed_area_fraction_ocean") > 0.08
    assert 0.0 < shallow_raw_fraction <= 0.03


def test_p108_audit_reports_boundary_width_and_high_mountain_coherence():
    spec = get_preset("earthlike")
    spec.grid_cells = 420
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    broad_ridge = np.abs(grid.lat) < 11.0
    land = grid.lat > 28.0
    elev = np.full(grid.n, -3800.0, dtype=np.float64)
    elev[land] = 800.0
    high_belt = land & (np.abs(grid.lon) < 30.0)
    elev[high_belt] = 3400.0
    isolated = int(np.argmax((grid.lat < 65.0) * (grid.lon > 120.0) * grid.lat))
    land[isolated] = True
    elev[isolated] = 5200.0

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 8).astype(float),
        "tectonics.plate_rank": np.zeros(grid.n, dtype=float),
        "tectonics.protected_plate_id": np.full(grid.n, -1.0),
        "tectonics.boundary_province_kind": np.where(high_belt, 8.0, -1.0),
        "terrain.elevation_m": elev,
        "crust.type": land.astype(float),
        "crust.age_myr": np.where(land, 1500.0, 50.0),
        "ocean.depth_province": np.where(elev < 0.0, 4.0, 0.0),
        "terrain.province": np.where(land, 5.0, 0.0),
        "terrain.continental_detail": np.where(land, 5.0, 0.0),
        "terrain.orogenic_load": np.where(high_belt, 0.8, 0.0),
        "tectonics.orogeny_age_myr": np.where(high_belt, 4450.0, -1.0),
    })
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(broad_ridge)[0].astype(np.int64),
        "collision": np.where(high_belt)[0].astype(np.int64),
        "transform": np.where((np.abs(grid.lat - 16.0) < 4.0) & (np.abs(grid.lon) < 60.0))[0].astype(np.int64),
    }
    world.set_g("terrain.last_p111611_bridge_candidate_cell_count", 7.0)
    world.set_g("terrain.last_p111611_bridge_candidate_area_fraction", 0.014)
    world.set_g("terrain.last_p111611_peak_hierarchy_shoulder_cell_count", 5.0)
    world.set_g("terrain.last_p111611_high_pair_count", 3.0)
    world.set_g("terrain.last_p111611_safe_path_count", 2.0)
    world.set_g("terrain.last_p111611_blocked_high_pair_count", 1.0)
    world.set_g("terrain.last_p111611_diagnostic_system_count", 2.0)
    world.set_g("terrain.last_p111611_bridge_deferred_no_safe_path", 1.0)
    world.set_g("terrain.last_p111622_spine_guided_repair_used", 1.0)
    world.set_g("terrain.last_p111622_spine_guided_candidate_cell_count", 11.0)
    world.set_g("terrain.last_p111622_spine_guided_candidate_area_fraction", 0.021)
    world.set_g("terrain.last_p111622_spine_guided_bridge_cell_count", 4.0)
    world.set_g("terrain.last_p111622_spine_guided_bridge_area_fraction", 0.008)
    world.set_g("terrain.last_p111622_high_component_delta", 3.0)
    world.set_g("terrain.last_p111622_top3_high_share_delta", 0.12)
    world.set_g("terrain.last_p111622_fragmentation_delta", 0.18)
    world.set_g("terrain.last_p111624_spine_to_terrain_response_used", 1.0)
    world.set_g("terrain.last_p111624_spine_response_candidate_cell_count", 9.0)
    world.set_g("terrain.last_p111624_spine_response_candidate_area_fraction", 0.017)
    world.set_g("terrain.last_p111624_spine_response_bridge_cell_count", 5.0)
    world.set_g("terrain.last_p111624_spine_response_bridge_area_fraction", 0.009)
    world.set_g("terrain.last_p111624_spine_saddle_response_cell_count", 6.0)
    world.set_g("terrain.last_p111624_spine_saddle_response_area_fraction", 0.011)
    world.set_g("terrain.last_p111624_spine_saddle_guard_reverted", 0.0)
    world.set_g("terrain.last_p111624_high_component_delta", 2.0)
    world.set_g("terrain.last_p111624_top3_high_share_delta", 0.10)
    world.set_g("terrain.last_p111624_fragmentation_delta", 0.16)
    world.set_g("terrain.last_p111624_guard_accepted_small_spine_response", 1.0)

    metrics = terminal_audit_metrics(world, [])
    width = metrics["boundary_width_diagnostics"]["boundary_ridge"]
    mountains = metrics["high_mountain_coherence"]

    assert width["cell_count"] > 0
    assert width["fraction_width_gt1"] > 0.0
    assert metrics["boundary_width_diagnostics"]["summary"]["max_p90_width_steps"] >= 1.0
    assert mountains["high_mountain_component_count"] >= 1
    assert mountains["high_mountain_parent_overlap_fraction"] > 0.0
    assert 0.0 <= mountains["top3_high_mountain_component_share"] <= 1.0
    assert (
        0.0
        <= mountains["small_high_mountain_component_area_fraction_of_high"]
        <= 1.0
    )
    assert 0.0 <= mountains["high_mountain_fragmentation_index"] <= 1.0
    assert mountains["isolated_extreme_peak_area_fraction_of_extreme"] >= 0.0
    assert mountains["p111611_bridge_candidate_cell_count"] == 7.0
    assert mountains["p111611_bridge_candidate_area_fraction"] == 0.014
    assert mountains["p111611_peak_hierarchy_shoulder_cell_count"] == 5.0
    assert mountains["p111611_high_pair_count"] == 3.0
    assert mountains["p111611_safe_path_count"] == 2.0
    assert mountains["p111611_blocked_high_pair_count"] == 1.0
    assert mountains["p111611_diagnostic_system_count"] == 2.0
    assert mountains["p111611_bridge_deferred_no_safe_path"] is True
    assert mountains["p111622_spine_guided_repair_used"] is True
    assert mountains["p111622_spine_guided_candidate_cell_count"] == 11.0
    assert mountains["p111622_spine_guided_candidate_area_fraction"] == 0.021
    assert mountains["p111622_spine_guided_bridge_cell_count"] == 4.0
    assert mountains["p111622_spine_guided_bridge_area_fraction"] == 0.008
    assert mountains["p111622_high_component_delta"] == 3.0
    assert mountains["p111622_top3_high_share_delta"] == 0.12
    assert mountains["p111622_fragmentation_delta"] == 0.18
    assert mountains["p111624_spine_to_terrain_response_used"] is True
    assert mountains["p111624_spine_response_candidate_cell_count"] == 9.0
    assert mountains["p111624_spine_response_candidate_area_fraction"] == 0.017
    assert mountains["p111624_spine_response_bridge_cell_count"] == 5.0
    assert mountains["p111624_spine_response_bridge_area_fraction"] == 0.009
    assert mountains["p111624_spine_saddle_response_cell_count"] == 6.0
    assert mountains["p111624_spine_saddle_response_area_fraction"] == 0.011
    assert mountains["p111624_spine_saddle_guard_reverted"] is False
    assert mountains["p111624_high_component_delta"] == 2.0
    assert mountains["p111624_top3_high_share_delta"] == 0.10
    assert mountains["p111624_fragmentation_delta"] == 0.16
    assert mountains["p111624_guard_accepted_small_spine_response"] is True


def test_p1116_audit_reports_boundary_linework_and_oceanic_age_structure():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TectonicsModule()
    terrain = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 72.0
    ocean = ~land
    ridge = ocean & (np.abs(grid.lat) < 3.5)
    transform = ocean & (np.abs(grid.lon - 45.0) < 4.0) & (np.abs(grid.lat) < 42.0)
    trench = ocean & (np.abs(grid.lat + 42.0) < 3.5) & (np.abs(grid.lon) < 120.0)
    subduction = ocean & (np.abs(grid.lat + 39.0) < 4.5) & (np.abs(grid.lon) < 125.0)

    age = np.where(ocean, 6.0 + np.abs(grid.lat) * 4.0, 1600.0)
    depth = np.where(ocean, 1700.0 + age * 8.0, 0.0)
    depth[ridge] = 2100.0
    depth[trench] = 6400.0
    elev = np.where(ocean, sea_level - depth, 900.0)
    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.float64)
    depth_province[ridge] = OCEAN_DEPTH_RIDGE
    depth_province[trench] = 6.0

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 10).astype(float),
        "tectonics.plate_rank": np.zeros(grid.n, dtype=float),
        "tectonics.protected_plate_id": np.full(grid.n, -1.0),
        "tectonics.boundary_province_kind": np.full(grid.n, -1.0),
        "terrain.elevation_m": elev,
        "crust.type": land.astype(float),
        "crust.age_myr": age,
        "ocean.depth_province": depth_province,
        "terrain.province": np.zeros(grid.n, dtype=float),
        "terrain.continental_detail": np.zeros(grid.n, dtype=float),
    })
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(ridge)[0].astype(np.int64),
        "transform": np.where(transform)[0].astype(np.int64),
        "trench": np.where(trench)[0].astype(np.int64),
        "subduction": np.where(subduction)[0].astype(np.int64),
    }
    boundary_polylines = module._boundary_polyline_objects(
        grid,
        world.networks["tectonics.boundaries"],
        [],
        world.fields["tectonics.plate_id"],
        world.time_myr,
    )
    world.objects["tectonics.boundary_polylines"] = boundary_polylines
    module._record_p129_boundary_polyline_metrics(
        world,
        grid,
        boundary_polylines,
    )
    assert terrain._tectonic_process_mask(world, "ridge").any()
    assert terrain._tectonic_process_mask(world, "trench").any()

    metrics = terminal_audit_metrics(world, [])
    linework = metrics["p1116_boundary_linework"]
    crust = metrics["p1116_oceanic_crust_structure"]

    assert metrics["acceptance"]["has_p1116_boundary_linework_metrics"]
    assert metrics["acceptance"]["has_p1116_oceanic_crust_structure_metrics"]
    assert linework["ridge"]["cell_count"] > 0
    assert linework["ridge"]["line_coherence_score"] > 0.0
    assert linework["trench_subduction_attachment_fraction"] > 0.75
    assert linework["p129_boundary_polyline_ready"] is True
    assert linework["p129_boundary_polyline_object_count"] > 0.0
    assert linework["p129_ridge_polyline_count"] > 0.0
    assert linework["p129_trench_polyline_count"] > 0.0
    assert linework["p129_boundary_polyline_all"]["cell_count"] > 0
    assert linework["p131_gap_bridge_cell_count"] >= 0.0
    assert linework["p131_gap_bridge_object_count"] >= 0.0
    assert isinstance(linework["p131_gap_bridge_used"], bool)
    assert "ridge" in linework["p131_gap_bridge_cell_count_by_kind"]
    assert metrics["object_count_by_set"]["tectonics.boundary_polylines"] == len(
        boundary_polylines)
    consumption = linework["p130_terrain_polyline_consumption"]
    assert consumption["ridge"]["polyline_used"] is True
    assert consumption["ridge"]["polyline_cell_count"] > 0.0
    assert consumption["trench"]["polyline_used"] is True
    assert consumption["trench"]["polyline_cell_count"] > 0.0
    assert crust["age_depth_correlation"] > 0.25
    assert crust["age_distance_from_ridge_correlation"] > 0.25
    assert crust["old_minus_young_median_depth_m"] > 500.0
    assert crust["far_minus_near_ridge_median_age_myr"] > 30.0
    assert crust["nonridge_minus_ridge_median_depth_m"] > 300.0
    assert crust["trench_minus_nontrench_median_depth_m"] > 2500.0


def test_p129_boundary_polyline_objects_order_boundary_components():
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TectonicsModule()

    center = next(
        idx for idx, neighbors in enumerate(grid.neighbors)
        if len(neighbors) >= 3
    )
    neighbors = [int(nb) for nb in grid.neighbors[center]]
    ridge_cells = [neighbors[0], center, neighbors[1]]
    trench_cells = [center, neighbors[2]]
    plate = (np.arange(grid.n) % 4).astype(np.float64)
    boundary_objects = [
        {
            "id": "ridge:test",
            "kind": "ridge",
            "stage": "young_ocean",
            "cells": ridge_cells,
            "parent_plate_ids": [1, 2],
            "birth_myr": 4400.0,
        },
        {
            "id": "trench:test",
            "kind": "trench",
            "stage": "subduction",
            "cells": trench_cells,
            "parent_plate_ids": [2, 3],
            "birth_myr": 4300.0,
        },
    ]

    objects = module._boundary_polyline_objects(
        grid,
        {},
        boundary_objects,
        plate,
        world.time_myr,
    )
    module._record_p129_boundary_polyline_metrics(world, grid, objects)

    ridge_axes = [
        obj for obj in objects
        if obj.get("kind") == "ridge_polyline"
    ]
    trench_axes = [
        obj for obj in objects
        if obj.get("kind") == "trench_polyline"
    ]

    assert len(ridge_axes) == 1
    assert len(trench_axes) == 1
    assert ridge_axes[0]["boundary_kind"] == "ridge"
    assert ridge_axes[0]["source_boundary_object_id"] == "ridge:test"
    assert ridge_axes[0]["source_component_cell_count"] == 3
    assert ridge_axes[0]["axis_cell_count"] >= 2
    assert len(ridge_axes[0]["endpoint_cells"]) == 2
    assert 0.0 <= ridge_axes[0]["directness"] <= 1.0
    assert ridge_axes[0]["sinuosity"] >= 1.0
    assert trench_axes[0]["role"] == "subduction_trench_axis"
    assert world.g("tectonics.last_p129_object_count") == 2.0
    assert world.g("tectonics.last_p129_ridge_polyline_count") == 1.0
    assert world.g("tectonics.last_p129_trench_polyline_count") == 1.0
    assert world.g("tectonics.last_p129_polyline_ready") == 1.0


def test_p131_boundary_polyline_bridges_guarded_short_gap():
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TectonicsModule()

    center = next(
        idx for idx, neighbors in enumerate(grid.neighbors)
        if len(neighbors) >= 4
    )
    neighbors = [int(nb) for nb in grid.neighbors[center]]
    source_pair = None
    for left_index, left in enumerate(neighbors):
        adjacent = {int(nb) for nb in grid.neighbors[left]}
        for right in neighbors[left_index + 1:]:
            if int(right) not in adjacent:
                source_pair = (int(left), int(right))
                break
        if source_pair is not None:
            break
    assert source_pair is not None

    left, right = source_pair
    plate = (np.arange(grid.n) % 4).astype(np.float64)
    boundary_objects = [
        {
            "id": "ridge:gapped",
            "kind": "ridge",
            "stage": "young_ocean",
            "cells": [left, right],
            "parent_plate_ids": [1, 2],
            "birth_myr": 4400.0,
        }
    ]

    objects = module._boundary_polyline_objects(
        grid,
        {},
        boundary_objects,
        plate,
        world.time_myr,
    )
    module._record_p129_boundary_polyline_metrics(world, grid, objects)

    ridge_axes = [
        obj for obj in objects
        if obj.get("kind") == "ridge_polyline"
    ]

    assert len(ridge_axes) == 1
    assert center in ridge_axes[0]["cells"]
    assert ridge_axes[0]["raw_source_component_cell_count"] == 2
    assert ridge_axes[0]["source_component_cell_count"] == 3
    assert ridge_axes[0]["p131_gap_bridge_cell_count"] == 1
    assert ridge_axes[0]["p131_gap_bridge_used"] is True
    assert world.g("tectonics.last_p131_gap_bridge_cell_count") == 1.0
    assert world.g("tectonics.last_p131_gap_bridge_object_count") == 1.0
    assert world.g("tectonics.last_p131_gap_bridge_used") == 1.0
    assert world.g("tectonics.last_p131_ridge_gap_bridge_cell_count") == 1.0


def test_p11162_coherent_trench_axis_bridges_sampled_subduction_line():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    land = grid.lat > 70.0
    ocean = ~land
    coast_pass = module._ocean_pass_distance(grid, land, ocean, max_passes=8)
    full_trench = ocean & (np.abs(grid.lat + 28.0) < 4.0) & (np.abs(grid.lon) < 130.0)
    sampled = module._generalized_boundary_source(grid, full_trench, spacing=2) & full_trench
    assert int(sampled.sum()) >= 4

    world.networks["tectonics.boundaries"] = {
        "trench": np.where(sampled)[0].astype(np.int64),
        "subduction": np.where(sampled)[0].astype(np.int64),
    }

    raw_components = module._components(grid, sampled)
    axis = module._coherent_subduction_trench_axis(world, ocean, coast_pass)
    axis_components = module._components(grid, axis)
    raw_same = module._same_neighbor_count(grid, sampled)
    axis_same = module._same_neighbor_count(grid, axis)

    assert int(axis.sum()) > int(sampled.sum())
    assert len(axis_components) < len(raw_components)
    assert float(np.mean(axis_same[axis] <= 0)) < float(np.mean(raw_same[sampled] <= 0))


def test_p11162_audit_prefers_object_trench_linework_over_raw_boundary_scatter():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    sea_level = 0.0
    land = grid.lat > 72.0
    ocean = ~land
    object_trench = ocean & (np.abs(grid.lat + 34.0) < 3.5) & (np.abs(grid.lon) < 135.0)
    raw_trench = object_trench & ((np.arange(grid.n) % 5) == 0)
    subduction = ocean & (np.abs(grid.lat + 32.0) < 5.0) & (np.abs(grid.lon) < 140.0)
    ridge = ocean & (np.abs(grid.lat) < 3.5)

    age = np.where(ocean, 20.0 + np.abs(grid.lat) * 2.5, 1600.0)
    depth = np.where(ocean, 3000.0 + 4.0 * age, 0.0)
    depth[ridge] = 2200.0
    depth[object_trench] = 6200.0
    elev = np.where(ocean, sea_level - depth, 850.0)
    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.float64)
    depth_province[ridge] = OCEAN_DEPTH_RIDGE
    depth_province[object_trench] = OCEAN_DEPTH_TRENCH

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 10).astype(float),
        "tectonics.plate_rank": np.zeros(grid.n, dtype=float),
        "tectonics.protected_plate_id": np.full(grid.n, -1.0),
        "tectonics.boundary_province_kind": np.full(grid.n, -1.0),
        "terrain.elevation_m": elev,
        "crust.type": land.astype(float),
        "crust.age_myr": age,
        "ocean.depth_province": depth_province,
        "terrain.province": np.zeros(grid.n, dtype=float),
        "terrain.continental_detail": np.zeros(grid.n, dtype=float),
    })
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(ridge)[0].astype(np.int64),
        "trench": np.where(raw_trench)[0].astype(np.int64),
        "subduction": np.where(subduction)[0].astype(np.int64),
    }
    world.objects["terrain.margin_landforms"] = [
        {"kind": "trench", "cells": np.where(object_trench)[0].astype(int).tolist()}
    ]

    linework = terminal_audit_metrics(world, [])["p1116_boundary_linework"]

    assert linework["schema"] == "aevum.p1116_boundary_linework.v3"
    assert linework["trench_source"] == "terrain_object"
    assert linework["trench"]["cell_count"] >= int(object_trench.sum())
    assert (
        linework["raw_trench_boundary"]["isolated_cell_fraction"]
        > linework["trench"]["isolated_cell_fraction"]
    )
    assert (
        linework["raw_trench_boundary"]["component_count"]
        > linework["trench"]["component_count"]
    )


def test_p11164_line_endpoint_bridge_connects_short_semantic_gaps():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()

    corridor = (
        (np.abs(grid.lat + 18.0) < 4.0)
        & (grid.lon > -145.0)
        & (grid.lon < 145.0)
    )
    anchors = module._generalized_boundary_source(grid, corridor, spacing=3) & corridor
    anchors &= (np.arange(grid.n) % 11) != 0
    raw_components = module._components(grid, anchors)
    assert len(raw_components) >= 3

    bridged = module._p11164_bridge_line_endpoints(
        grid,
        anchors,
        corridor,
        max_steps=5,
        max_added_area_fraction=0.80,
    )
    bridged_components = module._components(grid, bridged)
    raw_same = module._same_neighbor_count(grid, anchors)
    bridged_same = module._same_neighbor_count(grid, bridged)
    new_cells = bridged & ~anchors

    assert int(bridged.sum()) > int(anchors.sum())
    assert len(bridged_components) < len(raw_components)
    assert float(np.mean(bridged_same[bridged] <= 0)) < float(
        np.mean(raw_same[anchors] <= 0)
    )
    assert not np.any(new_cells & (bridged_same >= 4))


def test_p11164_audit_prefers_object_ridge_linework_over_raw_scatter():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    sea_level = 0.0
    land = grid.lat > 74.0
    ocean = ~land
    object_ridge = ocean & (np.abs(grid.lat - 2.0) < 3.5) & (np.abs(grid.lon) < 140.0)
    raw_ridge = object_ridge & ((np.arange(grid.n) % 2) == 0)
    transform = ocean & (np.abs(grid.lat - 12.0) < 3.0) & (np.abs(grid.lon) < 90.0)
    depth = np.where(ocean, 3700.0 + np.abs(grid.lat) * 12.0, 0.0)
    depth[object_ridge] = 2100.0
    elev = np.where(ocean, sea_level - depth, 700.0)
    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.float64)
    depth_province[object_ridge] = OCEAN_DEPTH_RIDGE

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 9).astype(float),
        "tectonics.plate_rank": np.zeros(grid.n, dtype=float),
        "tectonics.protected_plate_id": np.full(grid.n, -1.0),
        "tectonics.boundary_province_kind": np.full(grid.n, -1.0),
        "terrain.elevation_m": elev,
        "crust.type": land.astype(float),
        "crust.age_myr": np.where(ocean, 5.0 + np.abs(grid.lat) * 3.0, 1800.0),
        "ocean.depth_province": depth_province,
        "terrain.province": np.zeros(grid.n, dtype=float),
        "terrain.continental_detail": np.zeros(grid.n, dtype=float),
    })
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(raw_ridge)[0].astype(np.int64),
        "transform": np.where(transform)[0].astype(np.int64),
    }
    world.objects["terrain.ocean_fabric"] = [
        {"kind": "spreading_center", "cells": np.where(object_ridge)[0].astype(int).tolist()},
        {"kind": "transform_fault", "cells": np.where(transform)[0].astype(int).tolist()},
    ]

    linework = terminal_audit_metrics(world, [])["p1116_boundary_linework"]

    assert linework["schema"] == "aevum.p1116_boundary_linework.v3"
    assert linework["ridge_source"] == "terrain_object"
    assert linework["ridge"]["cell_count"] >= int(object_ridge.sum())
    assert (
        linework["raw_ridge_boundary"]["isolated_cell_fraction"]
        > linework["ridge"]["isolated_cell_fraction"]
    )
    assert (
        linework["raw_ridge_boundary"]["component_count"]
        > linework["ridge"]["component_count"]
    )


def test_p11165_ridge_adjacent_transform_survives_old_crust_age_window():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    sea_level = 0.0
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    crust_type = np.zeros(grid.n, dtype=np.float64)
    sediment = np.full(grid.n, 180.0, dtype=np.float64)
    ridge = (np.abs(grid.lat) < 4.0) & (np.abs(grid.lon) < 140.0)
    transform = (
        (np.abs(grid.lon - 24.0) < 12.0)
        & (np.abs(grid.lat) < 22.0)
    )
    ridge_transform_junction = transform & ridge
    near_ridge_transform = (
        transform & module._dilate_mask(grid, ridge, passes=1)
    )
    assert int(near_ridge_transform.sum()) >= 4
    assert int(ridge_transform_junction.sum()) >= 2

    age = np.full(grid.n, 140.0, dtype=np.float64)
    age[ridge] = 6.0
    depth_province = np.full(grid.n, OCEAN_DEPTH_ABYSS, dtype=np.float64)
    depth_province[ridge] = OCEAN_DEPTH_RIDGE

    world.fields.update({
        "crust.age_myr": age,
        "crust.type": crust_type,
    })
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(ridge)[0].astype(np.int64),
        "transform": np.where(transform)[0].astype(np.int64),
    }
    ocean_geo_fields = {
        "ocean.basin_id": np.zeros(grid.n, dtype=np.float64),
        "ocean.depth_province": depth_province,
        "ocean.shelf_width": np.zeros(grid.n, dtype=np.float64),
    }

    objects = module._ocean_fabric_objects(
        world,
        surface,
        sea_level,
        crust_type,
        sediment,
        ocean_geo_fields,
    )
    transform_cells = np.zeros(grid.n, dtype=bool)
    for obj in objects:
        if obj.get("kind") == "transform_fault":
            transform_cells |= module._cell_mask(grid, obj.get("cells", []))

    assert np.any(transform_cells & near_ridge_transform)
    assert np.any(transform_cells & ridge_transform_junction)
    assert (
        transform_cells & near_ridge_transform
    ).sum() >= max(2, int(0.50 * near_ridge_transform.sum()))


def test_p11166_trench_arc_refinement_bridges_short_gap_without_broadening():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    start = int(np.argmin(np.abs(grid.lat) + np.abs(grid.lon + 40.0)))
    path = [start]
    seen = {start}
    for _ in range(18):
        c = path[-1]
        candidates = [int(nb) for nb in grid.neighbors[c] if int(nb) not in seen]
        if not candidates:
            break
        candidates.sort(key=lambda nb: (abs(grid.lat[nb] - grid.lat[start]), -grid.lon[nb]))
        nxt = candidates[0]
        path.append(nxt)
        seen.add(nxt)
    assert len(path) >= 16

    support = np.zeros(grid.n, dtype=bool)
    support[np.asarray(path, dtype=np.int64)] = True
    axis = support.copy()
    for idx in (5, 10, 14):
        axis[path[idx]] = False
    ocean = np.ones(grid.n, dtype=bool)
    coast_pass = np.full(grid.n, 2, dtype=np.int64)

    refined = module._p11166_refine_subduction_trench_arc(
        world,
        axis,
        support,
        ocean,
        coast_pass,
    )

    assert world.g("terrain.last_p11166_trench_arc_refinement_accepted") == 1.0
    assert np.all(refined[axis])
    assert int(refined.sum()) > int(axis.sum())
    assert int(refined.sum()) <= int(np.ceil(axis.sum() * 1.16))
    assert len(module._components(grid, refined)) < len(module._components(grid, axis))


def test_p111632_parent_guided_trench_axis_bridges_terminal_subduction_arc():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    ocean = np.ones(grid.n, dtype=bool)
    coast_pass = np.full(grid.n, 4, dtype=np.int64)
    allowed = (
        (np.abs(grid.lat + 10.0) < 25.0)
        & (grid.lon > -150.0)
        & (grid.lon < 150.0)
    )
    start = int(np.argmin(
        np.where(allowed, (grid.lon + 130.0) ** 2 + (grid.lat + 12.0) ** 2, np.inf)
    ))
    target = int(np.argmin(
        np.where(allowed, (grid.lon - 130.0) ** 2 + (grid.lat + 8.0) ** 2, np.inf)
    ))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        neighbors = sorted(
            (int(nb) for nb in grid.neighbors[c]),
            key=lambda nb: (abs(grid.lat[nb] + 10.0), grid.lon[nb]),
        )
        for nb in neighbors:
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    path = np.asarray(path[::-1], dtype=np.int64)
    assert path.size >= 24

    parent = np.zeros(grid.n, dtype=bool)
    parent[path] = True
    axis = parent.copy()
    for idx in range(2, path.size - 2, 5):
        axis[int(path[idx])] = False
    assert len(module._components(grid, axis)) >= 4

    world.objects["tectonics.boundary_objects"] = [{
        "kind": "subduction_parent",
        "cells": path.astype(int).tolist(),
    }]

    refined = module._p111632_parent_guided_trench_axis(
        world,
        axis,
        parent | axis,
        ocean,
        coast_pass,
    )

    assert world.g("terrain.last_p111632_parent_trench_refinement_accepted") == 1.0
    assert int(refined.sum()) > int(axis.sum())
    assert len(module._components(grid, refined)) < len(module._components(grid, axis))
    assert world.g("terrain.last_p111632_parent_trench_components_after") < world.g(
        "terrain.last_p111632_parent_trench_components_before"
    )
    assert world.g(
        "terrain.last_p111632_parent_trench_parent_overlap_after"
    ) >= world.g("terrain.last_p111632_parent_trench_parent_overlap_before")
    assert world.g("terrain.last_p111632_parent_trench_area_fraction_ocean_after") <= (
        1.45
        * world.g("terrain.last_p111632_parent_trench_area_fraction_ocean_before")
    )


def test_p111632_parent_guided_trench_axis_requires_trench_anchor():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    ocean = np.ones(grid.n, dtype=bool)
    coast_pass = np.full(grid.n, 4, dtype=np.int64)
    parent = np.zeros(grid.n, dtype=bool)
    parent[np.where(
        (np.abs(grid.lat + 18.0) < 6.0)
        & (grid.lon > -120.0)
        & (grid.lon < 120.0)
    )[0][:24]] = True
    assert parent.any()
    world.objects["tectonics.boundary_objects"] = [{
        "kind": "subduction_parent",
        "cells": np.where(parent)[0].astype(int).tolist(),
    }]

    refined = module._p111632_parent_guided_trench_axis(
        world,
        np.zeros(grid.n, dtype=bool),
        parent,
        ocean,
        coast_pass,
    )

    assert not refined.any()
    assert world.g("terrain.last_p111632_parent_trench_refinement_accepted") == 0.0


def test_p11167_ordered_convergent_parent_axis_bridges_local_arc_gaps():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    module = TectonicsModule()

    corridor = (
        (np.abs(grid.lat + 16.0) < 4.5)
        & (grid.lon > -145.0)
        & (grid.lon < 145.0)
    )
    source = module._thin_boundary_mask(grid, corridor, spacing=3) & corridor
    source &= (np.arange(grid.n) % 13) != 0
    assert int(source.sum()) >= 8
    assert len(module._connected_components(grid, source)) >= 3

    parent = module._p11167_ordered_convergent_parent_axis(
        grid,
        source,
        support=corridor,
        excluded=np.zeros(grid.n, dtype=bool),
        sample_spacing=2,
        max_axis_seeds=96,
    )

    assert parent.any()
    assert np.all(module._dilate_mask(grid, corridor, passes=2)[parent])
    assert len(module._connected_components(grid, parent)) < len(
        module._connected_components(grid, source)
    )
    assert int(parent.sum()) <= max(int(source.sum() * 4.0), int(source.sum()) + 12)


def test_p11167_audit_prefers_convergent_parent_linework_over_raw_scatter():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    sea_level = 0.0
    land = grid.lat > 72.0
    ocean = ~land
    parent = ocean & (np.abs(grid.lat + 24.0) < 4.0) & (np.abs(grid.lon) < 135.0)
    raw = parent & ((np.arange(grid.n) % 5) == 0)
    trench = ocean & (np.abs(grid.lat + 27.0) < 3.5) & (np.abs(grid.lon) < 132.0)
    ridge = ocean & (np.abs(grid.lat) < 3.5)
    assert len(module._components(grid, raw)) > len(module._components(grid, parent))

    age = np.where(ocean, 18.0 + np.abs(grid.lat) * 2.0, 1600.0)
    depth = np.where(ocean, 3000.0 + age * 4.0, 0.0)
    depth[ridge] = 2200.0
    depth[trench] = 6200.0
    elev = np.where(ocean, sea_level - depth, 800.0)
    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.float64)
    depth_province[ridge] = OCEAN_DEPTH_RIDGE
    depth_province[trench] = OCEAN_DEPTH_TRENCH

    world.fields.update({
        "tectonics.plate_id": (np.arange(grid.n) % 8).astype(float),
        "tectonics.plate_rank": np.zeros(grid.n, dtype=float),
        "tectonics.protected_plate_id": np.full(grid.n, -1.0),
        "tectonics.boundary_province_kind": np.full(grid.n, -1.0),
        "terrain.elevation_m": elev,
        "crust.type": land.astype(float),
        "crust.age_myr": age,
        "ocean.depth_province": depth_province,
        "terrain.province": np.zeros(grid.n, dtype=float),
        "terrain.continental_detail": np.zeros(grid.n, dtype=float),
    })
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(ridge)[0].astype(np.int64),
        "convergent": np.where(raw)[0].astype(np.int64),
        "subduction": np.where(raw)[0].astype(np.int64),
        "trench": np.where(trench)[0].astype(np.int64),
    }
    world.objects["tectonics.boundary_objects"] = [
        {
            "kind": "convergent_parent",
            "cells": np.where(parent)[0].astype(int).tolist(),
        }
    ]

    linework = terminal_audit_metrics(world, [])["p1116_boundary_linework"]

    assert linework["convergent_source"] == "tectonic_parent_object"
    assert linework["convergent"]["cell_count"] >= int(parent.sum())
    assert (
        linework["raw_convergent_boundary"]["component_count"]
        > linework["convergent"]["component_count"]
    )
    assert (
        linework["raw_convergent_boundary"]["isolated_cell_fraction"]
        > linework["convergent"]["isolated_cell_fraction"]
    )


def test_p11168_parent_line_hierarchy_builds_guarded_crest_branch_foreland():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    parent_cells = np.where(
        (np.abs(grid.lat) < 6.0)
        & (grid.lon > -95.0)
        & (grid.lon < 95.0)
    )[0]
    parent_cells = parent_cells[:34]
    assert parent_cells.size >= 12
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=2)

    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:convergent_parent",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 820.0, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 950.0, 0.0) + np.where(parent, 300.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    assert hierarchy["crest"].any()
    assert hierarchy["branch"].any()
    assert hierarchy["foreland"].any()
    codes = hierarchy["hierarchy"]
    assert np.all(codes[hierarchy["crest"]] == 3.0)
    assert np.all(codes[hierarchy["branch"]] == 2.0)
    assert np.all(codes[hierarchy["foreland"]] == 1.0)
    assert world.g("terrain.last_p11168_parent_orogen_hierarchy_accepted") == 1.0
    assert world.g("terrain.last_p11168_crest_area_fraction") <= 0.036
    assert world.g("terrain.last_p11168_branch_area_fraction") <= 0.042
    assert world.g("terrain.last_p11168_foreland_area_fraction") <= 0.052


def test_p130_parent_hierarchy_consumes_boundary_polyline_parent_axis():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    parent_cells = np.where(
        (np.abs(grid.lat) < 6.0)
        & (grid.lon > -95.0)
        & (grid.lon < 95.0)
    )[0]
    parent_cells = parent_cells[:34]
    assert parent_cells.size >= 12
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=2)

    world.objects["tectonics.boundary_polylines"] = [{
        "id": "test:convergent_parent_polyline",
        "kind": "convergent_parent_polyline",
        "boundary_kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 820.0, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 950.0, 0.0) + np.where(parent, 300.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    assert hierarchy["crest"].any()
    assert hierarchy["branch"].any()
    assert hierarchy["foreland"].any()
    assert world.g("terrain.last_p11168_parent_orogen_hierarchy_accepted") == 1.0
    assert not world.objects.get("tectonics.boundary_objects")


def test_p146_parent_orogen_corridor_initializes_continuous_zones():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()
    world.set_g("terrain.enable_p147_boundary_group_parent_orogen_corridors", 1.0)
    world.set_g("terrain.enable_p148_p147_trial_promotion_guard", 0.0)
    world.set_g("terrain.p147_parent_orogen_corridor_max_side_to_trunk_ratio", 3.0)

    allowed = (
        (np.abs(grid.lat) < 18.0)
        & (grid.lon > -120.0)
        & (grid.lon < 90.0)
    )
    start = int(np.argmin(np.where(
        allowed,
        (grid.lat + 8.0) ** 2 + (grid.lon + 110.0) ** 2,
        np.inf,
    )))
    target = int(np.argmin(np.where(
        allowed,
        (grid.lat - 12.0) ** 2 + (grid.lon - 80.0) ** 2,
        np.inf,
    )))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        cell = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[cell]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = cell
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    parent_cells = np.asarray(path[::-1], dtype=np.int64)
    assert parent_cells.size >= 16
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=2)

    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p146_curved_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[active] = 1500.0
    surface[parent] = 2500.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1150.0, 0.0) + np.where(parent, 650.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    corridors = world.objects.get("terrain.parent_orogen_corridors", [])
    assert world.g("terrain.last_p146_parent_orogen_corridor_accepted") == 1.0
    assert world.g("terrain.last_p146_parent_orogen_corridor_object_count") >= 1.0
    assert world.g("terrain.last_p146_candidate_corridor_component_count") >= 1.0
    assert world.g("terrain.last_p146_retained_corridor_component_count") >= 1.0
    assert world.g("terrain.last_p146_retained_corridor_cell_count") > 0.0
    assert world.g("terrain.last_p146_trunk_parent_overlap_fraction") >= 0.15
    assert world.g("terrain.last_p146_reject_code") == 0.0
    assert world.g("terrain.last_p147_boundary_group_corridor_accepted") == 1.0
    assert world.g("terrain.last_p147_object_group_count") >= 1.0
    assert world.g("terrain.last_p147_accepted_group_count") >= 1.0
    assert world.g("terrain.last_p146_corridor_trunk_cell_count") > 0.0
    assert world.g("terrain.last_p146_corridor_branch_cell_count") > 0.0
    assert world.g("terrain.last_p146_corridor_foreland_cell_count") > 0.0
    assert corridors
    assert any(obj.get("trunk_cells") for obj in corridors)
    assert any(obj.get("branch_cells") for obj in corridors)
    assert any(obj.get("foreland_cells") for obj in corridors)
    assert hierarchy["crest"].any()
    assert hierarchy["branch"].any()
    assert hierarchy["foreland"].any()


def test_p147_boundary_group_corridors_keep_two_parent_belts_separate():
    spec = get_preset("earthlike")
    spec.grid_cells = 1200
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()
    world.set_g("terrain.enable_p147_boundary_group_parent_orogen_corridors", 1.0)
    world.set_g("terrain.enable_p148_p147_trial_promotion_guard", 0.0)
    world.set_g("terrain.p147_parent_orogen_corridor_max_side_to_trunk_ratio", 3.0)

    def path_between(lat0, lon0, lat1, lon1, allowed):
        start = int(np.argmin(np.where(
            allowed,
            (grid.lat - float(lat0)) ** 2 + (grid.lon - float(lon0)) ** 2,
            np.inf,
        )))
        target = int(np.argmin(np.where(
            allowed,
            (grid.lat - float(lat1)) ** 2 + (grid.lon - float(lon1)) ** 2,
            np.inf,
        )))
        prev = {start: -1}
        queue = [start]
        head = 0
        while head < len(queue) and target not in prev:
            cell = queue[head]
            head += 1
            for nb in sorted(int(x) for x in grid.neighbors[cell]):
                if nb in prev or not allowed[nb]:
                    continue
                prev[nb] = cell
                queue.append(nb)
        assert target in prev
        out = []
        cur = target
        while cur >= 0:
            out.append(cur)
            cur = prev[cur]
        return np.asarray(out[::-1], dtype=np.int64)

    belt_a = path_between(
        -10.0,
        -135.0,
        -6.0,
        -35.0,
        (grid.lat > -22.0) & (grid.lat < 4.0)
        & (grid.lon > -150.0) & (grid.lon < -20.0),
    )
    belt_b = path_between(
        10.0,
        25.0,
        16.0,
        135.0,
        (grid.lat > -3.0) & (grid.lat < 26.0)
        & (grid.lon > 10.0) & (grid.lon < 150.0),
    )
    assert belt_a.size >= 12
    assert belt_b.size >= 12
    parent = np.zeros(grid.n, dtype=bool)
    parent[belt_a] = True
    parent[belt_b] = True
    active = module._dilate_mask(grid, parent, passes=2)

    world.objects["tectonics.boundary_objects"] = [
        {
            "id": "test:p147_parent_belt_a",
            "kind": "convergent_parent",
            "cells": belt_a.astype(int).tolist(),
        },
        {
            "id": "test:p147_parent_belt_b",
            "kind": "convergent_parent",
            "cells": belt_b.astype(int).tolist(),
        },
    ]
    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[active] = 1500.0
    surface[parent] = 2500.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1150.0, 0.0) + np.where(parent, 650.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    corridors = world.objects.get("terrain.parent_orogen_corridors", [])
    source_ids = {
        str(obj.get("source_parent_group_id", ""))
        for obj in corridors
        if obj.get("source_parent_group_id")
    }
    assert world.g("terrain.last_p147_boundary_group_corridor_accepted") == 1.0
    assert world.g("terrain.last_p147_object_group_count") == 2.0
    assert world.g("terrain.last_p147_attempted_group_count") == 2.0
    assert world.g("terrain.last_p147_accepted_group_count") == 2.0
    assert world.g("terrain.last_p147_rejected_group_count") == 0.0
    assert world.g("terrain.last_p147_aggregate_guard_rejected") == 0.0
    assert world.g("terrain.last_p147_side_to_trunk_ratio") <= 3.0
    assert source_ids == {
        "test:p147_parent_belt_a",
        "test:p147_parent_belt_b",
    }
    assert all(str(obj.get("schema")) == "aevum.p147_parent_orogen_corridor.v1"
               for obj in corridors)
    assert hierarchy["crest"].any()
    assert hierarchy["branch"].any()
    assert hierarchy["foreland"].any()


def test_p148_trial_guard_rejects_fragmenting_p147_promotion():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()
    world.set_g("terrain.enable_p147_boundary_group_parent_orogen_corridors", 1.0)
    world.set_g("terrain.enable_p149_p147_staged_promotion", 0.0)
    world.set_g("terrain.enable_p150_group_trunk_quality_repair", 0.0)
    world.set_g("terrain.enable_p151_guarded_object_aware_cap", 0.0)
    world.set_g("terrain.p147_parent_orogen_corridor_max_side_to_trunk_ratio", 3.0)

    allowed = (
        (np.abs(grid.lat) < 18.0)
        & (grid.lon > -120.0)
        & (grid.lon < 90.0)
    )
    start = int(np.argmin(np.where(
        allowed,
        (grid.lat + 8.0) ** 2 + (grid.lon + 110.0) ** 2,
        np.inf,
    )))
    target = int(np.argmin(np.where(
        allowed,
        (grid.lat - 12.0) ** 2 + (grid.lon - 80.0) ** 2,
        np.inf,
    )))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        cell = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[cell]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = cell
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    parent_cells = np.asarray(path[::-1], dtype=np.int64)
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=2)

    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p148_trial_guard_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[active] = 1500.0
    surface[parent] = 2500.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1150.0, 0.0) + np.where(parent, 650.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    assert world.g("terrain.last_p147_attempted_group_count") >= 1.0
    assert world.g("terrain.last_p148_p147_trial_guard_enabled") == 1.0
    assert world.g("terrain.last_p148_p147_trial_guard_rejected") == 1.0
    assert world.g("terrain.last_p146_parent_orogen_corridor_accepted") == 0.0
    assert world.g("terrain.last_p146_reject_code") == 9.0
    assert world.g("terrain.last_p148_trial_added_corridor_cell_count") > 0.0
    assert hierarchy["crest"].any()
    assert hierarchy["branch"].any()
    assert hierarchy["foreland"].any()


def test_p149_staged_promotion_commits_trunk_without_side_bands():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()
    world.set_g("terrain.enable_p147_boundary_group_parent_orogen_corridors", 1.0)
    world.set_g("terrain.enable_p150_group_trunk_quality_repair", 0.0)
    world.set_g("terrain.enable_p151_guarded_object_aware_cap", 0.0)
    world.set_g("terrain.p147_parent_orogen_corridor_max_side_to_trunk_ratio", 3.0)

    allowed = (
        (np.abs(grid.lat) < 18.0)
        & (grid.lon > -120.0)
        & (grid.lon < 90.0)
    )
    start = int(np.argmin(np.where(
        allowed,
        (grid.lat + 8.0) ** 2 + (grid.lon + 110.0) ** 2,
        np.inf,
    )))
    target = int(np.argmin(np.where(
        allowed,
        (grid.lat - 12.0) ** 2 + (grid.lon - 80.0) ** 2,
        np.inf,
    )))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        cell = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[cell]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = cell
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    parent_cells = np.asarray(path[::-1], dtype=np.int64)
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=2)

    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p149_staged_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[active] = 1500.0
    surface[parent] = 2500.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1150.0, 0.0) + np.where(parent, 650.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    corridors = world.objects.get("terrain.parent_orogen_corridors", [])
    assert world.g("terrain.last_p148_p147_trial_guard_rejected") == 1.0
    assert world.g("terrain.last_p149_staged_promotion_accepted") == 1.0
    assert world.g("terrain.last_p149_trunk_trial_accepted") == 1.0
    assert world.g("terrain.last_p149_side_trial_rejected") == 1.0
    assert world.g("terrain.last_p146_parent_orogen_corridor_accepted") == 1.0
    assert world.g("terrain.last_p149_committed_trunk_cell_count") > 0.0
    assert world.g("terrain.last_p149_committed_branch_cell_count") == 0.0
    assert world.g("terrain.last_p149_committed_foreland_cell_count") == 0.0
    assert corridors
    assert all(
        str(obj.get("schema")) == "aevum.p149_staged_parent_orogen_corridor.v1"
        for obj in corridors
    )
    assert hierarchy["crest"].any()
    assert hierarchy["branch"].any()
    assert hierarchy["foreland"].any()


def test_p150_group_trunk_repair_bridges_fragmented_parent_group():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()
    world.set_g("terrain.enable_p147_boundary_group_parent_orogen_corridors", 1.0)
    world.set_g("terrain.enable_p148_p147_trial_promotion_guard", 0.0)
    world.set_g("terrain.enable_p150_group_trunk_quality_repair", 1.0)
    world.set_g("terrain.p147_parent_orogen_corridor_max_side_to_trunk_ratio", 4.0)
    world.set_g("terrain.p150_group_trunk_bridge_max_fraction", 0.020)
    world.set_g("terrain.p150_group_trunk_bridge_max_steps", 48.0)

    allowed = (
        (np.abs(grid.lat) < 15.0)
        & (grid.lon > -130.0)
        & (grid.lon < 80.0)
    )
    start = int(np.argmin(np.where(
        allowed,
        (grid.lat + 4.0) ** 2 + (grid.lon + 120.0) ** 2,
        np.inf,
    )))
    target = int(np.argmin(np.where(
        allowed,
        (grid.lat - 8.0) ** 2 + (grid.lon - 70.0) ** 2,
        np.inf,
    )))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        cell = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[cell]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = cell
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    parent_cells = np.asarray(path[::-1], dtype=np.int64)
    assert parent_cells.size >= 24

    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = np.zeros(grid.n, dtype=bool)
    segment_size = 3
    active[parent_cells[:segment_size]] = True
    mid = parent_cells.size // 2
    active[parent_cells[mid:mid + segment_size]] = True
    active[parent_cells[-segment_size:]] = True

    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p150_fragmented_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 40.0, dtype=np.float64)
    surface[active] = 1800.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1600.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    assert world.g("terrain.last_p150_group_trunk_repair_enabled") == 1.0
    assert world.g("terrain.last_p150_group_trunk_trial_count") >= 1.0
    assert world.g("terrain.last_p150_group_trunk_bridge_accepted_count") >= 1.0
    assert world.g("terrain.last_p150_group_trunk_bridge_cell_count") > 0.0
    assert world.g("terrain.last_p150_group_trunk_bridge_path_count") > 0.0
    assert world.g("terrain.last_p150_group_trunk_class_small_before_max") >= 0.999
    assert world.g("terrain.last_p150_group_trunk_class_small_after_max") < world.g(
        "terrain.last_p150_group_trunk_class_small_before_max"
    )
    assert world.g("terrain.last_p147_boundary_group_corridor_accepted") == 1.0
    assert world.g("terrain.last_p146_parent_orogen_corridor_accepted") == 1.0
    assert hierarchy["crest"].any()


def test_p132_parent_anchor_promotes_exact_boundary_spine_gap():
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    center, left, right = 9, 14, 22
    cells = np.asarray([center, left, right], dtype=np.int64)
    for cell in cells:
        assert cell in range(grid.n)
    assert right in {int(nb) for nb in grid.neighbors[center]}
    assert right in {int(nb) for nb in grid.neighbors[left]}

    surface = np.full(grid.n, -1000.0, dtype=np.float64)
    crust = np.zeros(grid.n, dtype=np.float64)
    active = np.zeros(grid.n, dtype=bool)
    score = np.zeros(grid.n, dtype=np.float64)
    surface[cells] = 2500.0
    crust[cells] = CONT
    active[cells] = True
    score[cells] = 1000.0
    for cell in cells:
        for neighbor in grid.neighbors[int(cell)]:
            nb = int(neighbor)
            surface[nb] = max(float(surface[nb]), 500.0)
            crust[nb] = CONT

    world.objects["tectonics.boundary_polylines"] = [{
        "id": "test:p132_exact_parent_anchor",
        "kind": "convergent_parent_polyline",
        "boundary_kind": "convergent_parent",
        "cells": [int(right)],
    }]

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    assert world.g("terrain.last_p132_parent_anchor_spine_promotion_accepted") == 1.0
    assert world.g("terrain.last_p132_candidate_cell_count") == 1.0
    assert world.g("terrain.last_p132_promoted_spine_cell_count") == 1.0
    assert world.g("terrain.last_p132_parent_aligned_spine_fraction_after") > world.g(
        "terrain.last_p132_parent_aligned_spine_fraction_before")
    assert world.g("terrain.last_p132_linework_score_after") >= world.g(
        "terrain.last_p132_linework_score_before")
    assert world.g("terrain.last_p132_spine_components_after") <= world.g(
        "terrain.last_p132_spine_components_before")
    assert world.g("terrain.last_p132_short_spine_components_after") <= world.g(
        "terrain.last_p132_short_spine_components_before")
    assert hierarchy["spine"][right] == 3.0
    assert np.count_nonzero(hierarchy["crest"] | hierarchy["branch"]) == 3


def test_p111610_parent_hierarchy_refinement_bridges_short_axis_gap():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 12.0)
        & (grid.lon > -90.0)
        & (grid.lon < 90.0)
    )
    start = int(np.argmin(np.where(allowed, grid.lat ** 2 + (grid.lon + 70.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(allowed, grid.lat ** 2 + (grid.lon - 70.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    path = np.asarray(path[::-1], dtype=np.int64)
    assert path.size >= 16

    active = np.zeros(grid.n, dtype=bool)
    active[path] = True
    gap = path[path.size // 2]
    active[int(gap)] = False
    parent_cells = path.astype(np.int64)
    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:convergent_parent_axis_gap",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 920.0, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1200.0, 0.0)
    score[parent_cells] += 400.0

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    gap_mask = np.zeros(grid.n, dtype=bool)
    gap_mask[int(gap)] = True
    gap_near = module._dilate_mask(grid, gap_mask, passes=1)
    assert np.any(hierarchy["hierarchy"][gap_near] >= 2.0)
    assert world.g("terrain.last_p111610_hierarchy_refinement_accepted") == 1.0
    assert world.g("terrain.last_p111610_bridge_cell_count") >= 1.0
    assert world.g("terrain.last_p111610_crest_components_after_refine") < world.g(
        "terrain.last_p111610_crest_components_before_refine"
    )


def test_p111613_hierarchy_geometry_refinement_prunes_unsupported_foreland_and_extends_branch():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    active_axis = np.where(
        (np.abs(grid.lat) < 8.0)
        & (grid.lon > -95.0)
        & (grid.lon < -25.0)
    )[0][:18]
    inactive_axis = np.where(
        (grid.lat > 42.0)
        & (grid.lat < 58.0)
        & (grid.lon > 80.0)
        & (grid.lon < 135.0)
    )[0][:14]
    assert active_axis.size >= 8
    assert inactive_axis.size >= 5
    parent_cells = np.unique(np.r_[active_axis, inactive_axis]).astype(np.int64)
    active = np.zeros(grid.n, dtype=bool)
    active[active_axis] = True
    active_near = module._dilate_mask(grid, active, passes=1)
    shoulder = active_near & ~active

    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p111613_mixed_parent_line",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 720.0, dtype=np.float64)
    surface[active] = 2600.0
    surface[shoulder] = 2150.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.full(grid.n, 120.0, dtype=np.float64)
    score[active] = 1650.0
    score[shoulder] = 1050.0
    score[parent_cells] += 500.0

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    assert world.g("terrain.last_p111613_hierarchy_geometry_refinement_accepted") == 1.0
    assert world.g("terrain.last_p111613_branch_extension_cell_count") > 0.0
    assert world.g("terrain.last_p111613_foreland_removed_cell_count") > 0.0
    assert world.g("terrain.last_p111613_foreland_removed_component_count") > 0.0
    assert world.g("terrain.last_p111613_foreland_components_after_refine") < world.g(
        "terrain.last_p111613_foreland_components_before_refine"
    )
    assert np.any(hierarchy["branch"] & shoulder)
    inactive_mask = np.zeros(grid.n, dtype=bool)
    inactive_mask[inactive_axis] = True
    inactive_near = module._dilate_mask(grid, inactive_mask, passes=3)
    assert not np.any(hierarchy["foreland"] & inactive_near)


def test_p111614_hierarchy_spines_are_thin_subsets_of_refined_hierarchy():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    parent_cells = np.where(
        (np.abs(grid.lat) < 10.0)
        & (grid.lon > -120.0)
        & (grid.lon < 120.0)
    )[0][:44]
    assert parent_cells.size >= 16
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=2)

    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p111614_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 820.0, dtype=np.float64)
    surface[active] = 2350.0
    surface[parent] = 2750.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1300.0, 0.0) + np.where(parent, 500.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    spine = hierarchy["spine"]
    crest_spine = hierarchy["crest_spine"]
    branch_spine = hierarchy["branch_spine"]
    assert world.g("terrain.last_p111614_spine_refinement_accepted") == 1.0
    assert np.any(crest_spine)
    assert np.all(hierarchy["crest"][crest_spine])
    assert np.all(hierarchy["branch"][branch_spine])
    assert not np.any(crest_spine & branch_spine)
    assert np.all(spine[crest_spine] == 3.0)
    assert np.all(spine[branch_spine] == 2.0)
    assert int(np.count_nonzero(spine)) < int(
        np.count_nonzero(hierarchy["crest"] | hierarchy["branch"])
    )
    assert world.g("terrain.last_p111614_crest_width_ratio") >= 1.0
    assert world.g("terrain.last_p111614_spine_overlap_valid") == 1.0


def test_p111623_branch_spine_continuity_reduces_peak_system_fragmentation():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    parent_cells = np.where(
        (np.abs(grid.lat) < 10.0)
        & (grid.lon > -120.0)
        & (grid.lon < 120.0)
    )[0][:44]
    assert parent_cells.size >= 16
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=2)

    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p111623_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 820.0, dtype=np.float64)
    surface[active] = 2350.0
    surface[parent] = 2750.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1300.0, 0.0) + np.where(parent, 500.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    assert world.g("terrain.last_p111623_branch_spine_continuity_accepted") == 1.0
    assert world.g("terrain.last_p111623_bridge_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111623_bridge_cell_count") > 0.0
    assert world.g("terrain.last_p111623_path_count") > 0.0
    assert world.g("terrain.last_p111623_peak_components_after_refine") <= world.g(
        "terrain.last_p111623_peak_components_before_refine"
    )
    assert world.g("terrain.last_p111623_branch_components_after_refine") <= (
        world.g("terrain.last_p111623_branch_components_before_refine") + 1.0
    )
    assert world.g("terrain.last_p111623_spine_components_after_refine") <= world.g(
        "terrain.last_p111623_spine_components_before_refine"
    )
    assert not np.any(hierarchy["branch"] & hierarchy["foreland"])


def test_p111635_spine_linework_smoothing_connects_parent_supported_fragments():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    parent_cells = np.where(
        (np.abs(grid.lat) < 10.0)
        & (grid.lon > -120.0)
        & (grid.lon < 120.0)
    )[0][:32]
    assert parent_cells.size >= 16
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=2)

    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p111635_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 820.0, dtype=np.float64)
    surface[active] = 2350.0
    surface[parent] = 2750.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1300.0, 0.0) + np.where(parent, 500.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    assert world.g("terrain.last_p111635_spine_linework_smoothing_accepted") == 1.0
    assert world.g("terrain.last_p111635_bridge_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111635_bridge_cell_count") > 0.0
    assert world.g("terrain.last_p111635_path_count") > 0.0
    assert world.g("terrain.last_p111635_spine_components_after_refine") < world.g(
        "terrain.last_p111635_spine_components_before_refine"
    )
    assert world.g(
        "terrain.last_p111635_short_spine_components_after_refine"
    ) <= world.g("terrain.last_p111635_short_spine_components_before_refine")
    assert world.g(
        "terrain.last_p111635_branch_attachment_fraction_after"
    ) >= world.g("terrain.last_p111635_branch_attachment_fraction_before")
    assert world.g("terrain.last_p11168_crest_area_fraction") <= 0.036
    assert world.g("terrain.last_p11168_branch_area_fraction") <= 0.042
    assert not np.any(hierarchy["crest_spine"] & hierarchy["branch_spine"])
    assert np.all(hierarchy["crest"][hierarchy["crest_spine"]])
    assert np.all(hierarchy["branch"][hierarchy["branch_spine"]])


def _find_disjoint_three_step_paths(grid: SphereGrid) -> tuple[list[int], list[int]]:
    def topology(path: list[int]) -> tuple[int, int]:
        path_set = set(path)
        endpoints = 0
        junctions = 0
        for c in path:
            degree = sum(1 for nb in grid.neighbors[int(c)] if int(nb) in path_set)
            if degree <= 1:
                endpoints += 1
            elif degree > 2:
                junctions += 1
        return endpoints, junctions

    for start in range(grid.n):
        start_neighbors = set(int(nb) for nb in grid.neighbors[start])
        by_target: dict[int, list[tuple[int, int]]] = {}
        for mid_a in sorted(start_neighbors):
            for mid_b in sorted(int(nb) for nb in grid.neighbors[mid_a]):
                if mid_b == start or mid_b in start_neighbors:
                    continue
                for target in sorted(int(nb) for nb in grid.neighbors[mid_b]):
                    if (
                        target == start
                        or target in start_neighbors
                        or target == mid_a
                    ):
                        continue
                    by_target.setdefault(target, []).append((mid_a, mid_b))
        for target, paths in by_target.items():
            for base_mid in paths:
                for ranked_mid in paths:
                    if base_mid == ranked_mid or set(base_mid) & set(ranked_mid):
                        continue
                    base = [
                        int(start),
                        int(base_mid[0]),
                        int(base_mid[1]),
                        int(target),
                    ]
                    ranked = [
                        int(start),
                        int(ranked_mid[0]),
                        int(ranked_mid[1]),
                        int(target),
                    ]
                    if topology(base) != topology(ranked):
                        continue
                    return (
                        base,
                        ranked,
                    )
    raise AssertionError("could not find disjoint three-step paths")


def test_p1141_ranked_shortest_path_prefers_ranked_equal_length_corridor():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()

    selected: tuple[int, tuple[int, int], tuple[int, int], int] | None = None
    for start in range(grid.n):
        start_neighbors = set(int(nb) for nb in grid.neighbors[start])
        by_target: dict[int, list[tuple[int, int]]] = {}
        for mid_a in sorted(start_neighbors):
            for mid_b in sorted(int(nb) for nb in grid.neighbors[mid_a]):
                if mid_b == start or mid_b in start_neighbors:
                    continue
                for target in sorted(int(nb) for nb in grid.neighbors[mid_b]):
                    if (
                        target == start
                        or target in start_neighbors
                        or target == mid_a
                    ):
                        continue
                    by_target.setdefault(target, []).append((mid_a, mid_b))
        for target, paths in by_target.items():
            for low_path in paths:
                for high_path in paths:
                    if low_path == high_path or set(low_path) & set(high_path):
                        continue
                    allowed = np.zeros(grid.n, dtype=bool)
                    allowed[[start, *low_path, *high_path, target]] = True
                    rank = np.zeros(grid.n, dtype=np.float64)
                    rank[list(high_path)] = 10.0
                    path = module._ranked_shortest_path_within_mask(
                        grid,
                        start,
                        target,
                        allowed,
                        rank,
                        max_steps=3,
                    )
                    if path.tolist() == [start, *high_path, target]:
                        selected = (start, low_path, high_path, target)
                        break
                if selected is not None:
                    break
            if selected is not None:
                break
        if selected is not None:
            break
    assert selected is not None
    start, low_path, high_path, target = selected
    allowed = np.zeros(grid.n, dtype=bool)
    allowed[[start, *low_path, *high_path, target]] = True
    rank = np.zeros(grid.n, dtype=np.float64)
    rank[list(high_path)] = 10.0

    path = module._ranked_shortest_path_within_mask(
        grid,
        start,
        target,
        allowed,
        rank,
        max_steps=3,
    )

    assert path.tolist() == [start, *high_path, target]


def test_p1142_ranked_orogenic_path_guard_rejects_less_direct_equal_length_path():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()

    selected: tuple[list[int], list[int]] | None = None
    for start in range(grid.n):
        start_neighbors = set(int(nb) for nb in grid.neighbors[start])
        by_target: dict[int, list[list[int]]] = {}
        for mid_a in sorted(start_neighbors):
            for mid_b in sorted(int(nb) for nb in grid.neighbors[mid_a]):
                if mid_b == start or mid_b in start_neighbors:
                    continue
                for target in sorted(int(nb) for nb in grid.neighbors[mid_b]):
                    if (
                        target == start
                        or target in start_neighbors
                        or target == mid_a
                    ):
                        continue
                    by_target.setdefault(target, []).append(
                        [start, mid_a, mid_b, target]
                    )
        for paths in by_target.values():
            for base in paths:
                for ranked in paths:
                    if base == ranked:
                        continue
                    base_penalty = module._path_directness_penalty(
                        grid, np.asarray(base, dtype=np.int64))
                    ranked_penalty = module._path_directness_penalty(
                        grid, np.asarray(ranked, dtype=np.int64))
                    if ranked_penalty > base_penalty + 0.12:
                        selected = (base, ranked)
                        break
                if selected is not None:
                    break
            if selected is not None:
                break
        if selected is not None:
            break
    assert selected is not None
    base, ranked = selected

    accepted, directness_rejected = module._ranked_orogenic_path_guard(
        grid,
        np.asarray(base, dtype=np.int64),
        np.asarray(ranked, dtype=np.int64),
        directness_slack=0.080,
    )

    assert accepted is False
    assert directness_rejected is True


def test_p1144_component_guard_rejects_ranked_path_that_loses_support():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    base, ranked = _find_disjoint_three_step_paths(grid)

    rank = np.zeros(grid.n, dtype=np.float64)
    rank[ranked] = 10.0
    relief = np.full(grid.n, 1200.0, dtype=np.float64)
    support = np.zeros(grid.n, dtype=bool)
    support[base] = True

    accepted, reason = module._ranked_orogenic_component_guard(
        grid,
        np.asarray(base, dtype=np.int64),
        np.asarray(ranked, dtype=np.int64),
        rank,
        relief,
        support,
        base_near_fraction_floor=0.0,
    )

    assert accepted is False
    assert reason == "support"


def test_p1144_component_guard_accepts_supported_ranked_relief_gain():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    base, ranked = _find_disjoint_three_step_paths(grid)

    rank = np.zeros(grid.n, dtype=np.float64)
    rank[ranked] = 10.0
    relief = np.full(grid.n, 1000.0, dtype=np.float64)
    relief[ranked] = 1400.0
    support = np.zeros(grid.n, dtype=bool)
    support[base] = True
    support[ranked] = True

    accepted, reason = module._ranked_orogenic_component_guard(
        grid,
        np.asarray(base, dtype=np.int64),
        np.asarray(ranked, dtype=np.int64),
        rank,
        relief,
        support,
    )

    assert accepted is True
    assert reason == "accepted"


def test_p111637a_spine_object_promotion_attaches_branch_paths_without_new_peaks():
    spec = get_preset("earthlike")
    spec.grid_cells = 1200
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    parent_cells = np.where(
        (np.abs(grid.lat) < 8.0)
        & (grid.lon > -140.0)
        & (grid.lon < 140.0)
    )[0][:70]
    assert parent_cells.size >= 20
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=2)

    world.set_g("terrain.enable_p111635_spine_linework_smoothing", 0.0)
    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p111637a_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 820.0, dtype=np.float64)
    surface[active] = 2100.0
    surface[parent] = 2900.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1300.0, 0.0) + np.where(parent, 500.0, 0.0)

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    assert world.g("terrain.last_p1141_ranked_orogenic_spine_paths_accepted") == 0.0
    assert world.g("terrain.last_p1141_ranked_path_used_count") == 0.0
    assert world.g("terrain.last_p111637a_spine_object_promotion_accepted") == 1.0
    assert world.g("terrain.last_p111637a_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111637a_bridge_cell_count") > 0.0
    assert world.g("terrain.last_p111637a_path_count") > 0.0
    assert world.g("terrain.last_p111637a_spine_components_after") <= world.g(
        "terrain.last_p111637a_spine_components_before"
    )
    assert world.g("terrain.last_p111637a_short_spine_components_after") <= world.g(
        "terrain.last_p111637a_short_spine_components_before"
    )
    assert world.g("terrain.last_p111637a_branch_attachment_fraction_after") > world.g(
        "terrain.last_p111637a_branch_attachment_fraction_before"
    )
    assert world.g("terrain.last_p111637a_spine_top3_share_after") >= world.g(
        "terrain.last_p111637a_spine_top3_share_before"
    )
    assert world.g("terrain.last_p111637a_linework_score_after") >= world.g(
        "terrain.last_p111637a_linework_score_before"
    )
    assert np.all(hierarchy["crest"][hierarchy["crest_spine"]])
    assert np.all(hierarchy["branch"][hierarchy["branch_spine"]])
    assert not np.any(hierarchy["crest_spine"] & hierarchy["branch_spine"])


def test_p1138_terminal_spine_stub_cleanup_demotes_unsupported_short_stubs():
    spec = get_preset("earthlike")
    spec.grid_cells = 2500
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 10.0)
        & (grid.lon > -80.0)
        & (grid.lon < 80.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 62.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 62.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        current = queue[head]
        head += 1
        for neighbor in sorted(int(x) for x in grid.neighbors[current]):
            if neighbor in prev or not allowed[neighbor]:
                continue
            prev[neighbor] = current
            queue.append(neighbor)
    assert target in prev
    path = []
    current = target
    while current >= 0:
        path.append(current)
        current = prev[current]
    path = np.asarray(path[::-1], dtype=np.int64)[:18]
    assert path.size >= 12
    path_mask = np.zeros(grid.n, dtype=bool)
    path_mask[path] = True

    blocked = module._dilate_mask(grid, path_mask, passes=3)
    stub_seed = int(np.argmax(np.where(~blocked, grid.lat, -np.inf)))
    stub_neighbor = next(
        int(nb) for nb in grid.neighbors[stub_seed] if not blocked[int(nb)])
    blocked[stub_seed] = True
    blocked[stub_neighbor] = True
    high_stub = int(np.argmax(np.where(~module._dilate_mask(
        grid,
        np.isin(np.arange(grid.n), [stub_seed, stub_neighbor]) | blocked,
        passes=2,
    ), -grid.lat, -np.inf)))

    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    hierarchy[path] = 2.0
    spine[path] = 2.0
    crest_cells = path[::5]
    hierarchy[crest_cells] = 3.0
    spine[crest_cells] = 3.0
    stub_cells = np.asarray([stub_seed, stub_neighbor, high_stub], dtype=np.int64)
    hierarchy[stub_cells] = [2.0, 3.0, 3.0]
    spine[stub_cells] = [2.0, 3.0, 3.0]
    crust = np.full(grid.n, CONT, dtype=np.float64)
    surface = np.full(grid.n, 640.0, dtype=np.float64)
    surface[path] = 2300.0
    surface[crest_cells] = 2900.0
    surface[[stub_seed, stub_neighbor]] = 900.0
    surface[high_stub] = 3300.0

    cleaned, telemetry = module._p1138_terminal_spine_stub_cleanup(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["demoted_cell_count"] == 3.0
    assert telemetry["component_count_after"] < telemetry["component_count_before"]
    assert telemetry["short_component_count_after"] < telemetry[
        "short_component_count_before"]
    assert telemetry["linework_score_after"] > telemetry["linework_score_before"]
    assert np.all(cleaned[path] == spine[path])
    assert np.all(cleaned[stub_cells] == 0.0)
    assert np.all(hierarchy[stub_cells] >= 2.0)


def test_p1139_terminal_spine_leaf_tip_cleanup_demotes_low_short_spurs_only():
    spec = get_preset("earthlike")
    spec.grid_cells = 2500
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 12.0)
        & (grid.lon > -95.0)
        & (grid.lon < 95.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 76.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 76.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        current = queue[head]
        head += 1
        for neighbor in sorted(int(x) for x in grid.neighbors[current]):
            if neighbor in prev or not allowed[neighbor]:
                continue
            prev[neighbor] = current
            queue.append(neighbor)
    assert target in prev
    path = []
    current = target
    while current >= 0:
        path.append(current)
        current = prev[current]
    path = np.asarray(path[::-1], dtype=np.int64)[:24]
    assert path.size >= 16
    path_set = set(int(x) for x in path)

    def find_spur(used: set[int]) -> tuple[int, int, int]:
        for junction in path[3:-3]:
            junction = int(junction)
            if junction in used:
                continue
            for leaf_inner in sorted(int(x) for x in grid.neighbors[junction]):
                if leaf_inner in path_set or leaf_inner in used:
                    continue
                inner_path_neighbors = [
                    int(x) for x in grid.neighbors[leaf_inner]
                    if int(x) in path_set
                ]
                if inner_path_neighbors != [junction]:
                    continue
                for leaf_tip in sorted(int(x) for x in grid.neighbors[leaf_inner]):
                    if (
                        leaf_tip == junction
                        or leaf_tip in path_set
                        or leaf_tip in used
                    ):
                        continue
                    tip_path_neighbors = [
                        int(x) for x in grid.neighbors[leaf_tip]
                        if int(x) in path_set
                    ]
                    tip_used_neighbors = [
                        int(x) for x in grid.neighbors[leaf_tip]
                        if int(x) in used
                    ]
                    if tip_path_neighbors or tip_used_neighbors:
                        continue
                    used.update({junction, leaf_inner, leaf_tip})
                    used.update(int(x) for x in grid.neighbors[leaf_inner])
                    used.update(int(x) for x in grid.neighbors[leaf_tip])
                    return junction, leaf_inner, leaf_tip
        raise AssertionError("no spur found")

    used: set[int] = set()
    low_junction, low_inner, low_tip = find_spur(used)
    high_junction, high_inner, high_tip = find_spur(used)

    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    hierarchy[path] = 2.0
    spine[path] = 2.0
    crest_cells = path[::6]
    hierarchy[crest_cells] = 3.0
    spine[crest_cells] = 3.0
    low_spur = np.asarray([low_inner, low_tip], dtype=np.int64)
    high_spur = np.asarray([high_inner, high_tip], dtype=np.int64)
    hierarchy[low_spur] = 2.0
    spine[low_spur] = 2.0
    hierarchy[high_spur] = 2.0
    spine[high_spur] = 2.0
    hierarchy[[low_junction, high_junction]] = 3.0
    spine[[low_junction, high_junction]] = 3.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[path] = 2500.0
    surface[crest_cells] = 3100.0
    surface[low_spur] = [1550.0, 1450.0]
    surface[high_spur] = [2550.0, 2450.0]

    cleaned, telemetry = module._p1139_terminal_spine_leaf_tip_cleanup(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["demoted_cell_count"] == 2.0
    assert telemetry["endpoint_count_after"] < telemetry["endpoint_count_before"]
    assert telemetry["linework_score_after"] >= telemetry["linework_score_before"]
    assert np.all(cleaned[path] == spine[path])
    assert np.all(cleaned[low_spur] == 0.0)
    assert np.all(cleaned[high_spur] == spine[high_spur])
    assert np.all(hierarchy[low_spur] == 2.0)


def test_p1140_terminal_spine_triangular_kink_cleanup_demotes_low_redundant_branch_only():
    spec = get_preset("earthlike")
    spec.grid_cells = 2500
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 14.0)
        & (grid.lon > -95.0)
        & (grid.lon < 95.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 72.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 72.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        current = queue[head]
        head += 1
        for neighbor in sorted(int(x) for x in grid.neighbors[current]):
            if neighbor in prev or not allowed[neighbor]:
                continue
            prev[neighbor] = current
            queue.append(neighbor)
    assert target in prev
    path = []
    current = target
    while current >= 0:
        path.append(current)
        current = prev[current]
    path = np.asarray(path[::-1], dtype=np.int64)[:24]
    assert path.size >= 16
    path_set = set(int(x) for x in path)

    def find_triangular_branch(used: set[int]) -> tuple[int, int, int]:
        for left, right in zip(path[3:-4], path[4:-3]):
            left = int(left)
            right = int(right)
            if left in used or right in used:
                continue
            mutual = sorted(
                int(nb)
                for nb in grid.neighbors[left]
                if int(nb) in set(int(x) for x in grid.neighbors[right])
            )
            for candidate in mutual:
                if candidate in path_set or candidate in used:
                    continue
                path_neighbors = {
                    int(nb)
                    for nb in grid.neighbors[candidate]
                    if int(nb) in path_set
                }
                if path_neighbors != {left, right}:
                    continue
                used.update({left, right, candidate})
                used.update(int(nb) for nb in grid.neighbors[candidate])
                return left, right, candidate
        raise AssertionError("no triangular branch found")

    used: set[int] = set()
    low_left, low_right, low_kink = find_triangular_branch(used)
    high_left, high_right, high_kink = find_triangular_branch(used)

    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    hierarchy[path] = 3.0
    spine[path] = 3.0
    hierarchy[[low_kink, high_kink]] = 2.0
    spine[[low_kink, high_kink]] = 2.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    surface = np.full(grid.n, 720.0, dtype=np.float64)
    surface[path] = 2650.0
    surface[[low_left, low_right, high_left, high_right]] = 2900.0
    surface[low_kink] = 1250.0
    surface[high_kink] = 2750.0

    cleaned, telemetry = module._p1140_terminal_spine_triangular_kink_cleanup(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["demoted_cell_count"] == 1.0
    assert telemetry["component_count_after"] == telemetry["component_count_before"]
    assert telemetry["endpoint_count_after"] <= telemetry["endpoint_count_before"]
    assert telemetry["redundant_kink_count_after"] < telemetry[
        "redundant_kink_count_before"]
    assert telemetry["linework_score_after"] >= telemetry["linework_score_before"]
    assert np.all(cleaned[path] == spine[path])
    assert cleaned[low_kink] == 0.0
    assert cleaned[high_kink] == spine[high_kink]
    assert hierarchy[low_kink] == 2.0


def test_p1149_orogenic_spine_object_inventory_builds_trunk_branch_topology():
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    center = next(
        idx for idx, neighbors in enumerate(grid.neighbors)
        if len(neighbors) >= 3
    )
    neighbors = [int(nb) for nb in grid.neighbors[center]]
    crest_cells = np.asarray([center, neighbors[0]], dtype=np.int64)
    branch_cells = np.asarray([neighbors[1]], dtype=np.int64)
    peak_cells = np.r_[crest_cells, branch_cells]

    surface = np.full(grid.n, -1200.0, dtype=np.float64)
    surface[peak_cells] = 1800.0
    crust = np.zeros(grid.n, dtype=np.float64)
    crust[peak_cells] = CONT
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[crest_cells] = 3.0
    hierarchy[branch_cells] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[crest_cells] = 3.0
    spine[branch_cells] = 2.0

    inventory = module._orogenic_spine_object_inventory(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    objects = inventory["objects"]
    telemetry = inventory["telemetry"]
    roles = {str(obj["role"]) for obj in objects}
    branch_objects = [
        obj for obj in objects
        if obj.get("role") == "branch"
    ]

    assert roles == {"system", "trunk", "branch"}
    assert telemetry["system_count"] == 1.0
    assert telemetry["trunk_count"] == 1.0
    assert telemetry["branch_count"] == 1.0
    assert telemetry["attached_branch_fraction"] == 1.0
    assert telemetry["spine_cell_count"] == 3.0
    assert branch_objects[0]["attached_to_trunk"] is True
    assert branch_objects[0]["parent_trunk_object_ids"]


def test_p128_orogenic_axis_polyline_inventory_orders_main_and_branch_axes():
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    center = next(
        idx for idx, neighbors in enumerate(grid.neighbors)
        if len(neighbors) >= 3
    )
    neighbors = [int(nb) for nb in grid.neighbors[center]]
    crest_cells = np.asarray([neighbors[0], center, neighbors[1]], dtype=np.int64)
    branch_cells = np.asarray([neighbors[2]], dtype=np.int64)
    peak_cells = np.r_[crest_cells, branch_cells]

    surface = np.full(grid.n, -1200.0, dtype=np.float64)
    surface[peak_cells] = 2200.0
    crust = np.zeros(grid.n, dtype=np.float64)
    crust[peak_cells] = CONT
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[crest_cells] = 3.0
    hierarchy[branch_cells] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[crest_cells] = 3.0
    spine[branch_cells] = 2.0

    spine_inventory = module._orogenic_spine_object_inventory(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )
    axis_inventory = module._orogenic_axis_polyline_inventory(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
        list(spine_inventory["objects"]),
    )

    objects = axis_inventory["objects"]
    telemetry = axis_inventory["telemetry"]
    main_axes = [
        obj for obj in objects
        if obj.get("role") == "main_crest_axis"
    ]
    branch_axes = [
        obj for obj in objects
        if obj.get("role") == "branch_axis"
    ]

    assert telemetry["object_count"] == 2.0
    assert telemetry["main_axis_count"] == 1.0
    assert telemetry["branch_axis_count"] == 1.0
    assert telemetry["attached_branch_axis_fraction"] == 1.0
    assert telemetry["polyline_ready"] == 1.0
    assert main_axes[0]["source_component_cell_count"] == 3
    assert main_axes[0]["axis_cell_count"] >= 2
    assert len(main_axes[0]["endpoint_cells"]) == 2
    assert 0.0 <= main_axes[0]["directness"] <= 1.0
    assert main_axes[0]["sinuosity"] >= 1.0
    assert branch_axes[0]["attached_to_main_axis"] is True
    assert branch_axes[0]["parent_axis_id"] == main_axes[0]["id"]


def test_p115_orogenic_spine_repair_candidates_bridge_disconnected_trunks():
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    mid = -1
    left = -1
    right = -1
    for candidate_mid, neighbors in enumerate(grid.neighbors):
        neighbors = [int(nb) for nb in neighbors]
        for i, a in enumerate(neighbors):
            for b in neighbors[i + 1:]:
                if b not in set(int(x) for x in grid.neighbors[a]):
                    mid = int(candidate_mid)
                    left = int(a)
                    right = int(b)
                    break
            if mid >= 0:
                break
        if mid >= 0:
            break
    assert mid >= 0

    peak_cells = np.asarray([left, mid, right], dtype=np.int64)
    surface = np.full(grid.n, -1200.0, dtype=np.float64)
    surface[peak_cells] = 1900.0
    crust = np.zeros(grid.n, dtype=np.float64)
    crust[peak_cells] = CONT
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[[left, right]] = 3.0
    hierarchy[mid] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[[left, right]] = 3.0

    inventory = module._orogenic_spine_repair_candidate_inventory(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    objects = inventory["objects"]
    telemetry = inventory["telemetry"]
    trunk_candidates = [
        obj for obj in objects
        if obj.get("role") == "trunk_bridge"
    ]

    assert telemetry["candidate_count"] >= 1.0
    assert telemetry["trunk_bridge_candidate_count"] >= 1.0
    assert telemetry["viable_candidate_count"] >= 1.0
    assert telemetry["candidate_cell_count"] >= 1.0
    assert trunk_candidates
    assert trunk_candidates[0]["proxy_viable"] is True
    assert mid in trunk_candidates[0]["cells"]


def test_p116_orogenic_spine_candidate_promotion_applies_guarded_bridge():
    spec = get_preset("earthlike")
    spec.grid_cells = 240
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    mid = -1
    left = -1
    right = -1
    for candidate_mid, neighbors in enumerate(grid.neighbors):
        neighbors = [int(nb) for nb in neighbors]
        for i, a in enumerate(neighbors):
            for b in neighbors[i + 1:]:
                if b not in set(int(x) for x in grid.neighbors[a]):
                    mid = int(candidate_mid)
                    left = int(a)
                    right = int(b)
                    break
            if mid >= 0:
                break
        if mid >= 0:
            break
    assert mid >= 0

    peak_cells = np.asarray([left, mid, right], dtype=np.int64)
    surface = np.full(grid.n, -1200.0, dtype=np.float64)
    surface[peak_cells] = 1900.0
    crust = np.zeros(grid.n, dtype=np.float64)
    crust[peak_cells] = CONT
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[[left, right]] = 3.0
    hierarchy[mid] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[[left, right]] = 3.0

    inventory = module._orogenic_spine_repair_candidate_inventory(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )
    promoted, telemetry = module._p116_orogenic_spine_candidate_promotion(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
        list(inventory["objects"]),
    )

    assert telemetry["accepted"] == 1.0
    assert telemetry["guard_reverted"] == 0.0
    assert telemetry["input_candidate_count"] >= 1.0
    assert telemetry["applied_candidate_count"] >= 1.0
    assert telemetry["applied_cell_count"] >= 1.0
    assert telemetry["component_count_after"] <= telemetry["component_count_before"]
    assert promoted[mid] == 2.0
    assert promoted[left] == 3.0
    assert promoted[right] == 3.0


def test_p118_non_orogenic_midland_suppression_preserves_supported_highlands():
    spec = get_preset("earthlike")
    spec.grid_cells = 2200
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    surface = np.full(grid.n, 120.0, dtype=np.float64)
    unsupported = np.arange(10, 36, dtype=np.int64)
    supported_orogen = np.arange(60, 66, dtype=np.int64)
    supported_mountain = np.arange(90, 96, dtype=np.int64)
    protected_plateau = np.arange(120, 126, dtype=np.int64)
    surface[unsupported] = 860.0
    surface[supported_orogen] = 880.0
    surface[supported_mountain] = 890.0
    surface[protected_plateau] = 910.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    detail = np.full(grid.n, CONT_DETAIL_PLATFORM, dtype=np.float64)
    detail[protected_plateau] = CONT_DETAIL_PLATEAU
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[supported_mountain] = 1.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[supported_orogen] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    shoulder = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)

    adjusted, telemetry = module._p118_non_orogenic_midland_stripe_suppression(
        world,
        surface,
        0.0,
        crust,
        detail,
        inventory,
        hierarchy,
        spine,
        shoulder,
        apron,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["candidate_cell_count"] == float(unsupported.size)
    assert telemetry["adjusted_cell_count"] == float(unsupported.size)
    assert telemetry["unsupported_midland_area_fraction_after"] < (
        telemetry["unsupported_midland_area_fraction_before"])
    assert np.all(adjusted[unsupported] < surface[unsupported])
    np.testing.assert_allclose(adjusted[supported_orogen], surface[supported_orogen])
    np.testing.assert_allclose(adjusted[supported_mountain], surface[supported_mountain])
    np.testing.assert_allclose(adjusted[protected_plateau], surface[protected_plateau])
    assert np.array_equal(adjusted >= 0.0, surface >= 0.0)


def test_p119_non_orogenic_interior_anti_raster_smoothing_is_guarded():
    spec = get_preset("earthlike")
    spec.grid_cells = 2200
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    neighbors: dict[int, list[int]] = {}
    for i, j in grid.edges:
        neighbors.setdefault(int(i), []).append(int(j))
        neighbors.setdefault(int(j), []).append(int(i))

    path = [10]
    used = {10}
    while len(path) < 72:
        current = path[-1]
        options = [nb for nb in neighbors.get(current, []) if nb not in used]
        if not options:
            break
        nxt = min(options)
        path.append(nxt)
        used.add(nxt)
    assert len(path) >= 36
    rough_path = np.asarray(path, dtype=np.int64)

    surface = np.full(grid.n, 280.0, dtype=np.float64)
    surface[rough_path[::2]] = 940.0
    surface[rough_path[1::2]] = 120.0
    supported_orogen = np.arange(300, 306, dtype=np.int64)
    supported_mountain = np.arange(360, 366, dtype=np.int64)
    protected_plateau = np.arange(420, 426, dtype=np.int64)
    surface[supported_orogen] = 880.0
    surface[supported_mountain] = 890.0
    surface[protected_plateau] = 910.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    detail = np.full(grid.n, CONT_DETAIL_PLATFORM, dtype=np.float64)
    detail[protected_plateau] = CONT_DETAIL_PLATEAU
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[supported_mountain] = 1.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[supported_orogen] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    shoulder = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)

    adjusted, telemetry = (
        module._p119_non_orogenic_interior_anti_raster_smoothing(
            world,
            surface,
            0.0,
            crust,
            detail,
            inventory,
            hierarchy,
            spine,
            shoulder,
            apron,
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["candidate_cell_count"] > 0.0
    assert telemetry["selected_cell_count"] > 0.0
    assert telemetry["adjusted_cell_count"] > 0.0
    assert telemetry["roughness_after_m"] < telemetry["roughness_before_m"]
    assert telemetry[
        "non_orogenic_highland_500_area_fraction_after"] <= telemetry[
            "non_orogenic_highland_500_area_fraction_before"]
    assert telemetry["p90_land_relief_after_m"] <= (
        telemetry["p90_land_relief_before_m"] + 1.0e-9)
    np.testing.assert_allclose(adjusted[supported_orogen], surface[supported_orogen])
    np.testing.assert_allclose(adjusted[supported_mountain], surface[supported_mountain])
    np.testing.assert_allclose(adjusted[protected_plateau], surface[protected_plateau])
    assert np.array_equal(adjusted >= 0.0, surface >= 0.0)


def test_p120_audit_reports_continental_semantic_region_geometry():
    spec = get_preset("earthlike")
    spec.grid_cells = 420
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    world.set_field("terrain.elevation_m", np.full(grid.n, 250.0))
    world.set_field("crust.type", np.full(grid.n, CONT, dtype=np.float64))
    raw = np.full(grid.n, CONT_DETAIL_PLATFORM, dtype=np.float64)
    raw[np.arange(0, grid.n, 11)] = CONT_DETAIL_SHIELD
    world.set_field("terrain.province", raw.copy())
    world.set_field("terrain.continental_detail", raw.copy())
    world.set_field("terrain.continental_detail_region_code",
                    np.full(grid.n, CONT_DETAIL_PLATFORM, dtype=np.float64))
    world.set_field("tectonics.internal_geographic_block_code", raw.copy())
    world.set_field("terrain.internal_geographic_block_region_code",
                    np.full(grid.n, 2.0, dtype=np.float64))
    world.set_field("terrain.inland_geomorphology_region_code",
                    np.full(grid.n, 2.0, dtype=np.float64))
    world.set_field("terrain.continental_province_code",
                    np.full(grid.n, 2.0, dtype=np.float64))
    world.set_field("terrain.mountain_inventory", np.zeros(grid.n))
    world.set_field("terrain.orogenic_parent_hierarchy", np.zeros(grid.n))
    world.set_field("terrain.orogenic_hierarchy_spine", np.zeros(grid.n))
    world.set_field("terrain.orogenic_shoulder_halo", np.zeros(grid.n))
    world.set_field("terrain.orogenic_highland_apron", np.zeros(grid.n))

    metrics = terminal_audit_metrics(world)
    p120 = metrics["p120_continental_semantic_geometry"]

    assert p120["schema"] == "aevum.p120_continental_semantic_geometry.v1"
    assert p120["region_fields_available"] is True
    assert p120["unsupported_interior_cell_count"] == grid.n
    assert p120["unsupported_detail_region_tiny_delta"] < 0.0
    assert p120["unsupported_internal_region_tiny_delta"] < 0.0
    raw_detail = p120["field_metrics"]["continental_detail_raw"][
        "unsupported_interior"]
    region_detail = p120["field_metrics"]["continental_detail_region"][
        "unsupported_interior"]
    raw_internal = p120["field_metrics"]["internal_block_raw"][
        "unsupported_interior"]
    region_internal = p120["field_metrics"]["internal_block_region"][
        "unsupported_interior"]
    assert region_detail["same_neighbor_fraction"] > (
        raw_detail["same_neighbor_fraction"])
    assert region_internal["same_neighbor_fraction"] > (
        raw_internal["same_neighbor_fraction"])
    assert metrics["acceptance"][
        "has_p120_continental_semantic_geometry_metrics"] is True


def test_p120_continental_semantic_geometry_repair_breaks_unsupported_stripe():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            xyz = []
            for y in range(self.height):
                for x in range(self.width):
                    xyz.append([float(x), float(y), 0.0])
            self.xyz = np.asarray(xyz, dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

    grid = RectGrid(90, 5)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    stripe = np.asarray([2 * grid.width + x for x in range(6, 84)],
                        dtype=np.int64)
    protected = int(stripe[20])

    surface = np.full(grid.n, -1000.0, dtype=np.float64)
    surface[stripe] = 220.0
    crust = np.zeros(grid.n, dtype=np.float64)
    crust[stripe] = CONT
    terrain_province = np.zeros(grid.n, dtype=np.float64)
    terrain_province[stripe] = 2.0
    detail_raw = np.zeros(grid.n, dtype=np.float64)
    detail_raw[stripe] = CONT_DETAIL_PLATFORM
    detail_region = detail_raw.copy()
    internal_region = np.zeros(grid.n, dtype=np.float64)
    internal_region[stripe] = INTERNAL_BLOCK_STABLE_PLATFORM
    inland_region = np.zeros(grid.n, dtype=np.float64)
    inland_region[stripe] = INLAND_PROVINCE_PLATFORM
    province = np.zeros(grid.n, dtype=np.float64)
    province[stripe] = CONT_PROVINCE_PLATFORM
    for start in range(8, stripe.size - 4, 10):
        cut = stripe[start:start + 3]
        province[cut] = CONT_PROVINCE_INTRACRATONIC_BASIN
        inland_region[cut] = INLAND_PROVINCE_SAG_BASIN

    inventory = np.zeros(grid.n, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[protected] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    halo = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)

    (
        terrain_out,
        internal_out,
        detail_out,
        telemetry,
    ) = module._p120_continental_semantic_geometry_repair(
        world,
        surface,
        0.0,
        crust,
        terrain_province,
        detail_raw,
        internal_region,
        detail_region,
        inland_region,
        province,
        inventory,
        hierarchy,
        spine,
        halo,
        apron,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["detail_elongated_after"] < (
        telemetry["detail_elongated_before"])
    assert telemetry["internal_elongated_after"] < (
        telemetry["internal_elongated_before"])
    assert telemetry["detail_tiny_after"] <= telemetry["detail_tiny_before"]
    assert telemetry["internal_tiny_after"] <= telemetry["internal_tiny_before"]
    assert detail_out[protected] == detail_region[protected]
    assert internal_out[protected] == internal_region[protected]
    assert terrain_out[protected] == terrain_province[protected]
    assert np.count_nonzero(detail_out[stripe] == CONT_DETAIL_BASIN) > 0
    assert np.count_nonzero(
        internal_out[stripe] == INTERNAL_BLOCK_INTRACRATONIC_BASIN) > 0


def test_p122_orogenic_detail_footprint_compaction_preserves_spine_core():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            xyz = []
            for y in range(self.height):
                for x in range(self.width):
                    xyz.append([float(x), float(y), 0.0])
            self.xyz = np.asarray(xyz, dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

    grid = RectGrid(72, 8)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    stripe = np.asarray([4 * grid.width + x for x in range(6, 66)],
                        dtype=np.int64)
    core = stripe[28:35]
    rift_flank = stripe[42:58]

    surface = np.full(grid.n, -1000.0, dtype=np.float64)
    surface[stripe] = 620.0
    surface[core] = 820.0
    crust = np.zeros(grid.n, dtype=np.float64)
    crust[stripe] = CONT
    terrain_province = np.zeros(grid.n, dtype=np.float64)
    terrain_province[stripe] = 5.0
    detail_raw = np.zeros(grid.n, dtype=np.float64)
    detail_raw[stripe] = CONT_DETAIL_PLATFORM
    detail_raw[core] = CONT_DETAIL_OROGEN
    detail_region = np.zeros(grid.n, dtype=np.float64)
    detail_region[stripe] = CONT_DETAIL_OROGEN
    internal_region = np.zeros(grid.n, dtype=np.float64)
    internal_region[stripe] = INTERNAL_BLOCK_MOBILE_BELT
    inland_region = np.zeros(grid.n, dtype=np.float64)
    inland_region[stripe] = INLAND_PROVINCE_PLATFORM
    inland_region[rift_flank] = INLAND_PROVINCE_RIFT
    province = np.zeros(grid.n, dtype=np.float64)
    province[stripe] = CONT_PROVINCE_PLATFORM
    province[rift_flank] = CONT_PROVINCE_RIFT_SYSTEM
    province[core] = CONT_PROVINCE_ACTIVE_OROGEN

    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[core] = 1.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[core] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[core] = 3.0

    (
        terrain_out,
        internal_out,
        detail_out,
        telemetry,
    ) = module._p122_orogenic_detail_footprint_compaction(
        world,
        surface,
        0.0,
        crust,
        terrain_province,
        detail_raw,
        internal_region,
        detail_region,
        inland_region,
        province,
        inventory,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["low_relief_orogen_area_fraction_after"] < (
        telemetry["low_relief_orogen_area_fraction_before"])
    assert telemetry["noncore_orogen_area_fraction_after"] < (
        telemetry["noncore_orogen_area_fraction_before"])
    assert np.all(detail_out[core] == CONT_DETAIL_OROGEN)
    assert np.count_nonzero(
        detail_out[rift_flank] == CONT_DETAIL_RIFT_BASIN) > 0
    assert np.count_nonzero(
        internal_out[rift_flank] == INTERNAL_BLOCK_RIFTED_MARGIN) > 0
    assert np.count_nonzero(detail_out[stripe] == CONT_DETAIL_PLATFORM) > 0
    assert np.count_nonzero(terrain_out[stripe] == 2.0) > 0
    assert np.all(detail_out[surface < 0.0] == detail_region[surface < 0.0])


def test_p123_spine_high_relief_continuity_expression_fills_spine_saddle():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

    grid = RectGrid(80, 30)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    row = 15
    spine_cells = np.asarray(
        [row * grid.width + x for x in range(10, 21)], dtype=np.int64)
    saddle = row * grid.width + 15
    high_left = np.asarray(
        [row * grid.width + x for x in range(10, 15)], dtype=np.int64)
    high_right = np.asarray(
        [row * grid.width + x for x in range(16, 21)], dtype=np.int64)

    surface = np.full(grid.n, 200.0, dtype=np.float64)
    surface[high_left] = 3100.0
    surface[high_right] = 3100.0
    surface[saddle] = 2200.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[spine_cells] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[spine_cells] = 3.0

    adjusted, telemetry = (
        module._p123_spine_high_relief_continuity_expression(
            world,
            surface,
            0.0,
            crust,
            hierarchy,
            spine,
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["spine_2400_component_count_after"] < (
        telemetry["spine_2400_component_count_before"])
    assert telemetry["spine_3000_component_count_after"] < (
        telemetry["spine_3000_component_count_before"])
    assert telemetry["high_component_count_after"] < (
        telemetry["high_component_count_before"])
    assert adjusted[saddle] >= 3000.0
    assert np.array_equal(adjusted >= 0.0, surface >= 0.0)
    off_spine = np.ones(grid.n, dtype=bool)
    off_spine[spine_cells] = False
    assert np.all(adjusted[off_spine] == surface[off_spine])


def test_p134_terminal_spine_relief_gap_bridging_merges_short_relief_gaps():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

    grid = RectGrid(80, 30)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    row = 15
    crest_spine = np.asarray(
        [row * grid.width + x for x in range(10, 39)], dtype=np.int64)
    gap_3000 = np.asarray(
        [row * grid.width + x for x in range(13, 16)], dtype=np.int64)
    gap_2400 = np.asarray(
        [row * grid.width + x for x in range(33, 36)], dtype=np.int64)
    high_3000 = np.asarray(
        [row * grid.width + x for x in [10, 11, 12, 16, 17, 18]],
        dtype=np.int64,
    )
    high_2400 = np.asarray(
        [row * grid.width + x for x in [30, 31, 32, 36, 37, 38]],
        dtype=np.int64,
    )

    surface = np.full(grid.n, 200.0, dtype=np.float64)
    surface[crest_spine] = 2100.0
    surface[high_3000] = 3150.0
    surface[gap_3000] = 2520.0
    surface[high_2400] = 2550.0
    surface[gap_2400] = 1850.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[crest_spine] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[crest_spine] = 3.0

    adjusted, telemetry = module._p134_terminal_spine_relief_gap_bridging(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["candidate_cell_count"] == 6.0
    assert telemetry["adjusted_cell_count"] == 6.0
    assert telemetry["bridge_gap_component_count_2400"] == 1.0
    assert telemetry["bridge_gap_component_count_3000"] == 1.0
    assert telemetry["spine_2400_component_count_after"] < (
        telemetry["spine_2400_component_count_before"])
    assert telemetry["spine_3000_component_count_after"] < (
        telemetry["spine_3000_component_count_before"])
    assert telemetry["high_component_count_after"] < (
        telemetry["high_component_count_before"])
    assert np.all(adjusted[gap_3000] >= 3000.0)
    assert np.all(adjusted[gap_2400] >= 2400.0)
    assert np.array_equal(adjusted >= 0.0, surface >= 0.0)
    off_spine = np.ones(grid.n, dtype=bool)
    off_spine[crest_spine] = False
    assert np.all(adjusted[off_spine] == surface[off_spine])


def test_p135_terminal_crest_core_relief_gap_bridging_merges_3000m_core_gap():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

    grid = RectGrid(80, 30)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    row = 15
    crest_spine = np.asarray(
        [row * grid.width + x for x in range(10, 21)], dtype=np.int64)
    high_left = np.asarray(
        [row * grid.width + x for x in range(10, 13)], dtype=np.int64)
    core_gap = np.asarray(
        [row * grid.width + x for x in range(13, 16)], dtype=np.int64)
    high_right = np.asarray(
        [row * grid.width + x for x in range(16, 21)], dtype=np.int64)

    surface = np.full(grid.n, 200.0, dtype=np.float64)
    surface[crest_spine] = 2500.0
    surface[high_left] = 3150.0
    surface[high_right] = 3100.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[crest_spine] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[crest_spine] = 3.0

    adjusted, telemetry = (
        module._p135_terminal_crest_core_relief_gap_bridging(
            world,
            surface,
            0.0,
            crust,
            hierarchy,
            spine,
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["gap_component_count"] == 1.0
    assert telemetry["candidate_cell_count"] == 3.0
    assert telemetry["adjusted_cell_count"] == 3.0
    assert telemetry["crest_3000_component_count_after"] < (
        telemetry["crest_3000_component_count_before"])
    assert telemetry["spine_3000_component_count_after"] < (
        telemetry["spine_3000_component_count_before"])
    assert telemetry["high_component_count_after"] < (
        telemetry["high_component_count_before"])
    assert np.all(adjusted[core_gap] >= 3000.0)
    assert np.array_equal(adjusted >= 0.0, surface >= 0.0)
    off_spine = np.ones(grid.n, dtype=bool)
    off_spine[crest_spine] = False
    assert np.all(adjusted[off_spine] == surface[off_spine])


def test_p136_terminal_high_peak_speckle_rebalance_lowers_tiny_high_components():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

    grid = RectGrid(80, 30)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    row = 15
    main_core = np.asarray(
        [row * grid.width + x for x in range(10, 16)], dtype=np.int64)
    tiny_pair = np.asarray(
        [row * grid.width + x for x in range(24, 26)], dtype=np.int64)
    tiny_single = np.asarray([row * grid.width + 34], dtype=np.int64)
    extreme_core = np.asarray(
        [row * grid.width + x for x in range(44, 47)], dtype=np.int64)
    all_peak = np.concatenate([main_core, tiny_pair, tiny_single, extreme_core])

    surface = np.full(grid.n, 200.0, dtype=np.float64)
    surface[main_core] = 3150.0
    surface[tiny_pair] = 3020.0
    surface[tiny_single] = 3005.0
    surface[extreme_core] = 4680.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[all_peak] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[all_peak] = 3.0

    adjusted, telemetry = module._p136_terminal_high_peak_speckle_rebalance(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["candidate_component_count"] == 2.0
    assert telemetry["lowered_component_count"] == 2.0
    assert telemetry["candidate_cell_count"] == 3.0
    assert telemetry["adjusted_cell_count"] == 3.0
    assert telemetry["high_component_count_after"] < (
        telemetry["high_component_count_before"])
    assert telemetry["high_cell_count_after"] < telemetry["high_cell_count_before"]
    assert np.all(adjusted[tiny_pair] == 2850.0)
    assert np.all(adjusted[tiny_single] == 2850.0)
    assert np.all(adjusted[main_core] == surface[main_core])
    assert np.all(adjusted[extreme_core] == surface[extreme_core])
    assert np.array_equal(adjusted >= 0.0, surface >= 0.0)


def test_p124_orogenic_spine_geometry_regularization_fills_supported_gap():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

    grid = RectGrid(80, 30)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    row = 15
    crest_cells = np.asarray(
        [row * grid.width + x for x in range(10, 21)], dtype=np.int64)
    gap = row * grid.width + 15
    left = np.asarray(
        [row * grid.width + x for x in range(10, 15)], dtype=np.int64)
    right = np.asarray(
        [row * grid.width + x for x in range(16, 21)], dtype=np.int64)

    surface = np.full(grid.n, 200.0, dtype=np.float64)
    surface[left] = 2600.0
    surface[right] = 2600.0
    surface[gap] = 2500.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[crest_cells] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[left] = 3.0
    spine[right] = 3.0

    regularized, telemetry = (
        module._p124_orogenic_spine_geometry_regularization(
            world,
            surface,
            0.0,
            crust,
            hierarchy,
            spine,
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] >= 1.0
    assert telemetry["added_cell_count"] == 1.0
    assert telemetry["component_count_after"] < telemetry["component_count_before"]
    assert telemetry["endpoint_count_after"] < telemetry["endpoint_count_before"]
    assert telemetry["linework_score_after"] > telemetry["linework_score_before"]
    assert regularized[gap] == 3.0
    unchanged = np.ones(grid.n, dtype=bool)
    unchanged[gap] = False
    assert np.all(regularized[unchanged] == spine[unchanged])


def test_p125_terminal_orogenic_semantic_land_consistency_clears_submerged_only():
    grid = SimpleNamespace(
        n=6,
        cell_area=np.ones(6, dtype=np.float64),
        lat=np.asarray([0.0, 12.0, 70.0, -20.0, 5.0, -64.0], dtype=np.float64),
        lon=np.asarray([0.0, 20.0, 175.0, 40.0, -174.0, 10.0], dtype=np.float64),
    )
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    surface = np.asarray([500.0, 1200.0, -20.0, -400.0, 150.0, -50.0])
    hierarchy = np.asarray([3.0, 2.0, 3.0, 0.0, 1.0, 0.0])
    spine = np.asarray([3.0, 2.0, 3.0, 0.0, 0.0, 0.0])
    halo = np.asarray([1.0, 0.0, 1.0, 1.0, 0.0, 0.0])
    apron = np.asarray([0.0, 1.0, 1.0, 1.0, 0.0, 0.0])

    fields, telemetry = module._p125_terminal_orogenic_semantic_land_consistency(
        world,
        surface,
        0.0,
        hierarchy,
        spine,
        halo,
        apron,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] == 2.0
    assert telemetry["cleared_hierarchy_cell_count"] == 1.0
    assert telemetry["cleared_spine_cell_count"] == 1.0
    assert telemetry["cleared_halo_cell_count"] == 2.0
    assert telemetry["cleared_apron_cell_count"] == 2.0
    assert telemetry["submerged_hierarchy_cell_count_after"] == 0.0
    assert telemetry["submerged_spine_cell_count_after"] == 0.0
    assert telemetry["submerged_halo_cell_count_after"] == 0.0
    assert telemetry["submerged_apron_cell_count_after"] == 0.0
    assert telemetry["polar_cleared_cell_count"] == 1.0
    assert telemetry["edge_cleared_cell_count"] == 1.0
    assert telemetry["land_semantic_cleared_cell_count"] == 0.0
    land = surface >= 0.0
    assert np.array_equal(
        fields["terrain.orogenic_parent_hierarchy"][land], hierarchy[land])
    assert np.array_equal(fields["terrain.orogenic_hierarchy_spine"][land], spine[land])
    assert np.array_equal(fields["terrain.orogenic_shoulder_halo"][land], halo[land])
    assert np.array_equal(fields["terrain.orogenic_highland_apron"][land], apron[land])
    submerged = ~land
    assert np.all(fields["terrain.orogenic_parent_hierarchy"][submerged] == 0.0)
    assert np.all(fields["terrain.orogenic_hierarchy_spine"][submerged] == 0.0)
    assert np.all(fields["terrain.orogenic_shoulder_halo"][submerged] == 0.0)
    assert np.all(fields["terrain.orogenic_highland_apron"][submerged] == 0.0)


def test_p126_terminal_orogenic_fringe_support_regularization_clears_fragments():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

        def cell(self, x: int, y: int) -> int:
            return int(y) * self.width + int(x)

    grid = RectGrid(30, 20)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    surface = np.full(grid.n, 500.0, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    halo = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)

    crest = [grid.cell(x, 3) for x in range(3, 7)]
    hierarchy[crest] = 3.0
    spine[crest] = 3.0
    attached_foreland = [grid.cell(x, 4) for x in range(3, 8)]
    detached_foreland = [grid.cell(10, 1), grid.cell(10, 2)]
    hierarchy[attached_foreland] = 1.0
    hierarchy[detached_foreland] = 1.0
    halo_keep = [grid.cell(2, y) for y in range(2, 7)]
    halo_clear = [grid.cell(8, 4), grid.cell(8, 5)]
    halo[halo_keep] = 1.0
    halo[halo_clear] = 1.0
    apron_keep = [grid.cell(x, 5) for x in range(4, 8)]
    apron_clear = [grid.cell(9, 6)]
    apron[apron_keep] = 1.0
    apron[apron_clear] = 1.0

    fields, telemetry = (
        module._p126_terminal_orogenic_fringe_support_regularization(
            world,
            surface,
            0.0,
            hierarchy,
            spine,
            halo,
            apron,
        )
    )

    out_hierarchy = fields["terrain.orogenic_parent_hierarchy"]
    out_halo = fields["terrain.orogenic_shoulder_halo"]
    out_apron = fields["terrain.orogenic_highland_apron"]
    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] == 5.0
    assert telemetry["cleared_foreland_cell_count"] == 2.0
    assert telemetry["cleared_halo_cell_count"] == 2.0
    assert telemetry["cleared_apron_cell_count"] == 1.0
    assert telemetry["peak_hierarchy_changed_cell_count"] == 0.0
    assert telemetry["spine_changed_cell_count"] == 0.0
    assert telemetry["foreland_component_count_after"] < (
        telemetry["foreland_component_count_before"])
    assert telemetry["halo_component_count_after"] < (
        telemetry["halo_component_count_before"])
    assert telemetry["apron_component_count_after"] < (
        telemetry["apron_component_count_before"])
    assert np.all(out_hierarchy[crest] == hierarchy[crest])
    assert np.all(spine[crest] == 3.0)
    assert np.all(out_hierarchy[attached_foreland] == 1.0)
    assert np.all(out_hierarchy[detached_foreland] == 0.0)
    assert np.all(out_halo[halo_keep] == 1.0)
    assert np.all(out_halo[halo_clear] == 0.0)
    assert np.all(out_apron[apron_keep] == 1.0)
    assert np.all(out_apron[apron_clear] == 0.0)


def test_p133_terminal_orogenic_fringe_band_ordering_clears_misordered_bands():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

        def cell(self, x: int, y: int) -> int:
            return int(y) * self.width + int(x)

    grid = RectGrid(30, 20)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    surface = np.full(grid.n, 500.0, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    halo = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)

    crest = [grid.cell(x, 5) for x in range(3, 8)]
    hierarchy[crest] = 3.0
    spine[crest] = 3.0
    foreland_keep = [grid.cell(x, 7) for x in range(3, 8)]
    foreland_clear = [grid.cell(18, 14)]
    hierarchy[foreland_keep] = 1.0
    hierarchy[foreland_clear] = 1.0
    surface[foreland_keep] = 850.0
    surface[foreland_clear] = 2850.0

    halo_keep = [grid.cell(2, y) for y in range(4, 8)]
    halo_clear = [grid.cell(19, 14)]
    halo[halo_keep] = 1.0
    halo[halo_clear] = 1.0
    surface[halo_keep] = 700.0
    surface[halo_clear] = 180.0

    apron_keep = [grid.cell(x, 8) for x in range(4, 8)]
    apron_clear = [grid.cell(20, 14)]
    apron[apron_keep] = 1.0
    apron[apron_clear] = 1.0
    surface[apron_keep] = 2200.0
    surface[apron_clear] = 1200.0

    fields, telemetry = module._p133_terminal_orogenic_fringe_band_ordering(
        world,
        surface,
        0.0,
        hierarchy,
        spine,
        halo,
        apron,
    )

    out_hierarchy = fields["terrain.orogenic_parent_hierarchy"]
    out_halo = fields["terrain.orogenic_shoulder_halo"]
    out_apron = fields["terrain.orogenic_highland_apron"]
    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] == 3.0
    assert telemetry["band_violation_cell_count"] == 3.0
    assert telemetry["cleared_foreland_cell_count"] == 1.0
    assert telemetry["cleared_halo_cell_count"] == 1.0
    assert telemetry["cleared_apron_cell_count"] == 1.0
    assert telemetry["peak_hierarchy_changed_cell_count"] == 0.0
    assert telemetry["spine_changed_cell_count"] == 0.0
    assert telemetry["foreland_component_count_after"] < (
        telemetry["foreland_component_count_before"])
    assert telemetry["halo_component_count_after"] < (
        telemetry["halo_component_count_before"])
    assert telemetry["apron_component_count_after"] < (
        telemetry["apron_component_count_before"])
    assert np.all(out_hierarchy[crest] == hierarchy[crest])
    assert np.all(spine[crest] == 3.0)
    assert np.all(out_hierarchy[foreland_keep] == 1.0)
    assert np.all(out_hierarchy[foreland_clear] == 0.0)
    assert np.all(out_halo[halo_keep] == 1.0)
    assert np.all(out_halo[halo_clear] == 0.0)
    assert np.all(out_apron[apron_keep] == 1.0)
    assert np.all(out_apron[apron_clear] == 0.0)


def test_p137_terminal_orogenic_fringe_gap_consolidation_bridges_short_gaps():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

        def cell(self, x: int, y: int) -> int:
            return int(y) * self.width + int(x)

    grid = RectGrid(32, 18)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    surface = np.full(grid.n, 700.0, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    halo = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)

    crest = [grid.cell(x, 5) for x in range(4, 26)]
    hierarchy[crest] = 3.0
    spine[crest] = 3.0
    surface[crest] = 2600.0

    fore_left = [grid.cell(6, 8), grid.cell(7, 8)]
    fore_gap = [grid.cell(8, 8)]
    fore_right = [grid.cell(9, 8), grid.cell(10, 8)]
    hierarchy[fore_left + fore_right] = 1.0
    surface[fore_left + fore_gap + fore_right] = 900.0

    halo_left = [grid.cell(14, 7), grid.cell(15, 7)]
    halo_gap = [grid.cell(16, 7)]
    halo_right = [grid.cell(17, 7), grid.cell(18, 7)]
    halo[halo_left + halo_right] = 1.0
    surface[halo_left + halo_gap + halo_right] = 1100.0

    apron_left = [grid.cell(20, 7), grid.cell(21, 7)]
    apron_gap = [grid.cell(22, 7)]
    apron_right = [grid.cell(23, 7), grid.cell(24, 7)]
    apron[apron_left + apron_right] = 1.0
    surface[apron_left + apron_gap + apron_right] = 1900.0

    fields, telemetry = module._p137_terminal_orogenic_fringe_gap_consolidation(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
        halo,
        apron,
    )

    out_hierarchy = fields["terrain.orogenic_parent_hierarchy"]
    out_halo = fields["terrain.orogenic_shoulder_halo"]
    out_apron = fields["terrain.orogenic_highland_apron"]
    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] == 3.0
    assert telemetry["added_foreland_cell_count"] == 1.0
    assert telemetry["added_halo_cell_count"] == 1.0
    assert telemetry["added_apron_cell_count"] == 1.0
    assert telemetry["foreland_gap_component_count"] == 1.0
    assert telemetry["halo_gap_component_count"] == 1.0
    assert telemetry["apron_gap_component_count"] == 1.0
    assert telemetry["foreland_component_count_after"] < (
        telemetry["foreland_component_count_before"])
    assert telemetry["halo_component_count_after"] < (
        telemetry["halo_component_count_before"])
    assert telemetry["apron_component_count_after"] < (
        telemetry["apron_component_count_before"])
    assert telemetry["peak_hierarchy_changed_cell_count"] == 0.0
    assert telemetry["spine_changed_cell_count"] == 0.0
    assert np.all(out_hierarchy[fore_gap] == 1.0)
    assert np.all(out_halo[halo_gap] == 1.0)
    assert np.all(out_apron[apron_gap] == 1.0)
    assert np.all(out_hierarchy[crest] == 3.0)
    assert np.all(spine[crest] == 3.0)


def test_p138_terminal_orogenic_fringe_component_thickening_grows_tiny_foreland():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

        def cell(self, x: int, y: int) -> int:
            return int(y) * self.width + int(x)

    grid = RectGrid(30, 20)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    surface = np.full(grid.n, 700.0, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    halo = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)

    crest = [grid.cell(x, 3) for x in range(4, 13)]
    hierarchy[crest] = 3.0
    spine[crest] = 3.0
    surface[crest] = 2600.0
    foreland = [grid.cell(x, 6) for x in range(4, 10)]
    growth_cells = [grid.cell(x, 7) for x in range(4, 7)]
    hierarchy[foreland] = 1.0
    surface[foreland + growth_cells] = 850.0

    fields, telemetry = (
        module._p138_terminal_orogenic_fringe_component_thickening(
            world,
            surface,
            0.0,
            crust,
            hierarchy,
            spine,
            halo,
            apron,
        )
    )

    out_hierarchy = fields["terrain.orogenic_parent_hierarchy"]
    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] == 3.0
    assert telemetry["added_foreland_cell_count"] == 3.0
    assert telemetry["grown_foreland_component_count"] == 1.0
    assert telemetry["foreland_component_count_before"] == 1.0
    assert telemetry["foreland_component_count_after"] == 1.0
    assert telemetry["foreland_small_area_fraction_before"] == 1.0
    assert telemetry["foreland_small_area_fraction_after"] == 0.0
    assert telemetry["peak_hierarchy_changed_cell_count"] == 0.0
    assert telemetry["spine_changed_cell_count"] == 0.0
    assert np.count_nonzero(out_hierarchy == 1.0) == len(foreland) + 3
    assert np.all(out_hierarchy[crest] == 3.0)
    assert np.all(spine[crest] == 3.0)


def test_p139_terminal_branch_range_component_thickening_promotes_supported_halo():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

        def cell(self, x: int, y: int) -> int:
            return int(y) * self.width + int(x)

    grid = RectGrid(30, 20)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    surface = np.full(grid.n, 700.0, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    halo = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)

    crest = [grid.cell(x, 3) for x in range(4, 14)]
    hierarchy[crest] = 3.0
    spine[crest] = 3.0
    surface[crest] = 2600.0
    branch = [grid.cell(x, 6) for x in range(4, 12)]
    hierarchy[branch] = 2.0
    spine[branch[:4]] = 2.0
    surface[branch] = 1800.0
    shoulder = [grid.cell(6, 7)]
    halo[shoulder] = 1.0
    surface[shoulder] = 1850.0

    fields, telemetry = (
        module._p139_terminal_branch_range_component_thickening(
            world,
            surface,
            0.0,
            crust,
            hierarchy,
            spine,
            halo,
            apron,
        )
    )

    out_hierarchy = fields["terrain.orogenic_parent_hierarchy"]
    out_halo = fields["terrain.orogenic_shoulder_halo"]
    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] == 1.0
    assert telemetry["added_branch_cell_count"] == 1.0
    assert telemetry["cleared_halo_cell_count"] == 1.0
    assert telemetry["grown_branch_component_count"] == 1.0
    assert telemetry["branch_cell_count_before"] == 8.0
    assert telemetry["branch_cell_count_after"] == 9.0
    assert telemetry["branch_component_count_before"] == 1.0
    assert telemetry["branch_component_count_after"] == 1.0
    assert telemetry["branch_small_area_fraction_before"] == 1.0
    assert telemetry["branch_small_area_fraction_after"] == 0.0
    assert telemetry["crest_changed_cell_count"] == 0.0
    assert telemetry["spine_changed_cell_count"] == 0.0
    assert np.all(out_hierarchy[shoulder] == 2.0)
    assert np.all(out_halo[shoulder] == 0.0)
    assert np.all(out_hierarchy[crest] == 3.0)
    assert np.all(spine[branch[:4]] == 2.0)


def test_p127_terminal_orogenic_spine_node_thinning_demotes_overwide_node():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

        def cell(self, x: int, y: int) -> int:
            return int(y) * self.width + int(x)

    grid = RectGrid(50, 50)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    surface = np.full(grid.n, 700.0, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    block = [
        grid.cell(x, y)
        for y in range(20, 23)
        for x in range(20, 23)
    ]
    center = grid.cell(21, 21)
    hierarchy[block] = 3.0
    spine[block] = 3.0

    thinned, telemetry = module._p127_terminal_orogenic_spine_node_thinning(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] >= 1.0
    assert telemetry["demoted_cell_count"] == 1.0
    assert telemetry["component_count_after"] == telemetry["component_count_before"]
    assert telemetry["endpoint_count_after"] == telemetry["endpoint_count_before"]
    assert telemetry["high_degree_count_after"] < (
        telemetry["high_degree_count_before"])
    assert telemetry["junction_pressure_after"] < (
        telemetry["junction_pressure_before"])
    assert telemetry["linework_score_after"] > telemetry["linework_score_before"]
    assert thinned[center] == 0.0
    unchanged = np.ones(grid.n, dtype=bool)
    unchanged[center] = False
    assert np.all(thinned[unchanged] == spine[unchanged])


def test_p140_terminal_branch_spine_component_thickening_completes_branch_axis():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

        def cell(self, x: int, y: int) -> int:
            return int(y) * self.width + int(x)

    grid = RectGrid(50, 50)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    surface = np.full(grid.n, 700.0, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)

    crest = [grid.cell(x, 20) for x in range(10, 19)]
    branch = [grid.cell(x, 23) for x in range(10, 19)]
    hierarchy[crest] = 3.0
    hierarchy[branch] = 2.0
    spine[crest] = 3.0
    spine[branch[:3]] = 2.0
    surface[crest] = 2600.0
    surface[branch] = 1800.0

    thickened, telemetry = (
        module._p140_terminal_branch_spine_component_thickening(
            world,
            surface,
            0.0,
            crust,
            hierarchy,
            spine,
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] == 6.0
    assert telemetry["promoted_branch_spine_cell_count"] == 6.0
    assert telemetry["grown_branch_spine_component_count"] == 1.0
    assert telemetry["branch_spine_cell_count_before"] == 3.0
    assert telemetry["branch_spine_cell_count_after"] == 9.0
    assert telemetry["branch_spine_component_count_before"] == 1.0
    assert telemetry["branch_spine_component_count_after"] == 1.0
    assert telemetry["branch_spine_small_area_fraction_before"] == 1.0
    assert telemetry["branch_spine_small_area_fraction_after"] == 0.0
    assert telemetry["all_spine_component_count_after"] <= (
        telemetry["all_spine_component_count_before"])
    assert telemetry["hierarchy_changed_cell_count"] == 0.0
    assert telemetry["crest_spine_changed_cell_count"] == 0.0
    assert np.all(thickened[branch] == 2.0)
    assert np.all(thickened[crest] == 3.0)
    non_branch = np.ones(grid.n, dtype=bool)
    non_branch[branch] = False
    assert np.all(thickened[non_branch] == spine[non_branch])


def test_p141_terminal_crest_spine_component_thickening_bridges_short_main_gap():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

        def cell(self, x: int, y: int) -> int:
            return int(y) * self.width + int(x)

    grid = RectGrid(50, 50)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    surface = np.full(grid.n, 700.0, dtype=np.float64)
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)

    crest = [grid.cell(x, 20) for x in range(4, 13)]
    left_spine = [grid.cell(x, 20) for x in range(4, 7)]
    right_spine = [grid.cell(x, 20) for x in range(10, 13)]
    gap = [grid.cell(x, 20) for x in range(7, 10)]
    hierarchy[crest] = 3.0
    spine[left_spine + right_spine] = 3.0
    surface[crest] = 2400.0

    thickened, telemetry = (
        module._p141_terminal_crest_spine_component_thickening(
            world,
            surface,
            0.0,
            crust,
            hierarchy,
            spine,
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] == 3.0
    assert telemetry["promoted_crest_spine_cell_count"] == 3.0
    assert telemetry["connected_crest_spine_component_count"] == 1.0
    assert telemetry["grown_crest_spine_component_count"] == 0.0
    assert telemetry["crest_spine_cell_count_before"] == 6.0
    assert telemetry["crest_spine_cell_count_after"] == 9.0
    assert telemetry["crest_spine_component_count_before"] == 2.0
    assert telemetry["crest_spine_component_count_after"] == 1.0
    assert telemetry["crest_spine_small_area_fraction_before"] == 1.0
    assert telemetry["crest_spine_small_area_fraction_after"] == 0.0
    assert telemetry["all_spine_component_count_after"] < (
        telemetry["all_spine_component_count_before"])
    assert telemetry["hierarchy_changed_cell_count"] == 0.0
    assert telemetry["branch_spine_changed_cell_count"] == 0.0
    assert np.all(thickened[crest] == 3.0)
    assert np.all(thickened[gap] == 3.0)
    non_crest = np.ones(grid.n, dtype=bool)
    non_crest[crest] = False
    assert np.all(thickened[non_crest] == spine[non_crest])


def test_p121_semantic_region_elevation_response_is_guarded():
    class RectGrid:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self.n = self.width * self.height
            self.cell_area = np.ones(self.n, dtype=np.float64)
            self.xyz = np.zeros((self.n, 3), dtype=np.float64)
            edges: list[tuple[int, int]] = []
            neighbors: list[list[int]] = [[] for _ in range(self.n)]
            for y in range(self.height):
                for x in range(self.width):
                    c = y * self.width + x
                    if x + 1 < self.width:
                        nb = y * self.width + x + 1
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
                    if y + 1 < self.height:
                        nb = (y + 1) * self.width + x
                        edges.append((c, nb))
                        neighbors[c].append(nb)
                        neighbors[nb].append(c)
            self.edges = np.asarray(edges, dtype=np.int64)
            self.neighbors = [
                np.asarray(items, dtype=np.int64) for items in neighbors
            ]

    grid = RectGrid(80, 30)
    world = SimpleNamespace(grid=grid)
    module = TerrainModule()
    surface = np.full(grid.n, 520.0, dtype=np.float64)
    for y in range(5, 25):
        for x in range(8, 72):
            c = y * grid.width + x
            surface[c] = 900.0 if (x + y) % 2 == 0 else 420.0
    protected = 12 * grid.width + 30
    surface[protected] = 940.0

    crust = np.full(grid.n, CONT, dtype=np.float64)
    raw_detail = np.full(grid.n, CONT_DETAIL_PLATFORM, dtype=np.float64)
    terrain_province = np.full(grid.n, 2.0, dtype=np.float64)
    internal_region = np.full(
        grid.n, INTERNAL_BLOCK_STABLE_PLATFORM, dtype=np.float64)
    detail_region = np.full(grid.n, CONT_DETAIL_PLATFORM, dtype=np.float64)
    inventory = np.zeros(grid.n, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[protected] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    halo = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)

    adjusted, telemetry = module._p121_semantic_region_elevation_response(
        world,
        surface,
        0.0,
        crust,
        raw_detail,
        terrain_province,
        internal_region,
        detail_region,
        inventory,
        hierarchy,
        spine,
        halo,
        apron,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["candidate_cell_count"] > 0.0
    assert telemetry["selected_cell_count"] > 0.0
    assert telemetry["adjusted_cell_count"] > 0.0
    assert telemetry["roughness_after_m"] < telemetry["roughness_before_m"]
    assert telemetry[
        "non_orogenic_highland_500_area_fraction_after"] <= telemetry[
            "non_orogenic_highland_500_area_fraction_before"]
    assert telemetry[
        "non_orogenic_highland_1000_area_fraction_after"] <= telemetry[
            "non_orogenic_highland_1000_area_fraction_before"]
    assert telemetry["p90_land_relief_after_m"] <= (
        telemetry["p90_land_relief_before_m"] + 1.0e-9)
    assert adjusted[protected] == surface[protected]
    assert np.array_equal(adjusted >= 0.0, surface >= 0.0)


def test_p111637b_final_spine_elevation_consistency_lifts_supported_spine_only():
    spec = get_preset("earthlike")
    spec.grid_cells = 4000
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    def walk(seed: int, length: int, blocked: set[int]) -> np.ndarray:
        path = [int(seed)]
        while len(path) < length:
            current = path[-1]
            next_cells = [
                int(nb)
                for nb in grid.neighbors[current]
                if int(nb) not in blocked and int(nb) not in path
            ]
            assert next_cells
            path.append(next_cells[0])
        return np.asarray(path, dtype=np.int64)

    safe_path = walk(300, 5, set())
    blocked_cells = set(safe_path.tolist())
    bridge_seed = 2500
    while bridge_seed in blocked_cells:
        bridge_seed += 17
    bridge_path = walk(bridge_seed, 5, blocked_cells)
    blocked_cells.update(bridge_path.tolist())
    unsupported = int(np.argmax(grid.lat))
    while (
        unsupported in blocked_cells
        or np.any(np.isin(grid.neighbors[unsupported], list(blocked_cells)))
    ):
        unsupported = int((unsupported + 137) % grid.n)

    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    hierarchy[safe_path] = 2.0
    spine[safe_path] = 2.0
    hierarchy[bridge_path] = 2.0
    spine[bridge_path] = 2.0
    hierarchy[unsupported] = 2.0
    spine[unsupported] = 2.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    surface = np.full(grid.n, -1800.0, dtype=np.float64)
    surface[safe_path[0]] = 760.0
    surface[safe_path[1]] = 120.0
    surface[safe_path[2:]] = -260.0
    surface[bridge_path[0]] = 760.0
    surface[bridge_path[-1]] = 820.0
    surface[bridge_path[1:-1]] = -260.0
    surface[unsupported] = -260.0

    adjusted, telemetry = module._p111637b_final_spine_elevation_consistency(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
        previous_land=surface >= 0.0,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] >= 4.0
    assert telemetry["adjusted_cell_count"] >= 4.0
    assert telemetry["bridge_blocked_submerged_spine_cell_count"] >= 3.0
    assert telemetry["submerged_spine_cell_count_after"] < telemetry[
        "submerged_spine_cell_count_before"]
    assert telemetry["underexpressed_spine_cell_count_after"] < telemetry[
        "underexpressed_spine_cell_count_before"]
    assert telemetry["land_delta_area_fraction"] <= 0.0025
    assert np.all(adjusted[safe_path[2:]] >= 520.0)
    assert adjusted[safe_path[1]] >= 520.0
    assert np.all(adjusted[bridge_path[1:-1]] == surface[bridge_path[1:-1]])
    assert adjusted[unsupported] == surface[unsupported]


def test_p111637b_final_spine_elevation_consistency_expresses_midres_ridge():
    spec = get_preset("earthlike")
    spec.grid_cells = 2500
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 16.0)
        & (grid.lon > -70.0)
        & (grid.lon < 70.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 55.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 55.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        current = queue[head]
        head += 1
        for neighbor in sorted(int(x) for x in grid.neighbors[current]):
            if neighbor in prev or not allowed[neighbor]:
                continue
            prev[neighbor] = current
            queue.append(neighbor)
    assert target in prev
    path = []
    current = target
    while current >= 0:
        path.append(current)
        current = prev[current]
    path = np.asarray(path[::-1], dtype=np.int64)[:18]
    assert path.size >= 12

    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    hierarchy[path] = 2.0
    spine[path] = 2.0
    crest_cells = path[::6]
    hierarchy[crest_cells] = 3.0
    spine[crest_cells] = 3.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    surface = np.full(grid.n, 640.0, dtype=np.float64)
    surface[path] = np.linspace(1380.0, 1880.0, path.size)
    surface[crest_cells] = 2220.0
    old_land = surface >= 0.0
    before_2400 = int(np.count_nonzero(
        (spine[path] >= 2.0) & (surface[path] >= 2400.0)))

    adjusted, telemetry = module._p111637b_final_spine_elevation_consistency(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
        previous_land=old_land,
    )

    after_2400 = int(np.count_nonzero(adjusted[path] >= 2400.0))
    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] >= path.size
    assert telemetry["adjusted_cell_count"] > 0.0
    assert telemetry["land_delta_area_fraction"] == 0.0
    assert telemetry["underexpressed_spine_cell_count_after"] < telemetry[
        "underexpressed_spine_cell_count_before"]
    assert after_2400 > before_2400
    assert np.array_equal(adjusted >= 0.0, old_land)
    path_mask = np.zeros(grid.n, dtype=bool)
    path_mask[path] = True
    assert np.all(adjusted[~path_mask] == surface[~path_mask])


def test_p111637c_anti_raster_orogenic_expression_smooths_existing_land_only():
    spec = get_preset("earthlike")
    spec.grid_cells = 2500
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    path = [500]
    while len(path) < 18:
        current = path[-1]
        next_cells = [int(nb) for nb in grid.neighbors[current] if int(nb) not in path]
        assert next_cells
        path.append(next_cells[0])
    path = np.asarray(path, dtype=np.int64)

    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    halo = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)
    hierarchy[path] = 2.0
    spine[path[::5]] = 2.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    surface = np.full(grid.n, -1800.0, dtype=np.float64)
    surface[path] = np.where(np.arange(path.size) % 2 == 0, 760.0, 1660.0)
    extreme = int(path[-1])
    surface[extreme] = 4700.0
    old_land = surface >= 0.0

    adjusted, telemetry = module._p111637c_anti_raster_orogenic_expression(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
        halo,
        apron,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["candidate_cell_count"] > 0.0
    assert telemetry["adjusted_cell_count"] > 0.0
    assert telemetry["roughness_after_m"] < telemetry["roughness_before_m"]
    assert np.array_equal(adjusted >= 0.0, old_land)
    assert np.all(adjusted[~old_land] == surface[~old_land])
    assert adjusted[extreme] == surface[extreme]


def test_p111615_spine_width_morphology_adds_halo_without_pruning_high_peaks():
    spec = get_preset("earthlike")
    spec.grid_cells = 4000
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 12.0)
        & (grid.lon > -140.0)
        & (grid.lon < 140.0)
    )
    start = int(np.argmin(np.where(
        allowed,
        grid.lat ** 2 + (grid.lon + 125.0) ** 2,
        np.inf,
    )))
    target = int(np.argmin(np.where(
        allowed,
        grid.lat ** 2 + (grid.lon - 125.0) ** 2,
        np.inf,
    )))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    parent_cells = np.asarray(path[::-1], dtype=np.int64)
    assert parent_cells.size >= 18
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=4)

    parent_core = module._dilate_mask(grid, parent, passes=1)
    parent_near = module._dilate_mask(grid, parent, passes=3)
    high_cells = np.where(active & parent_near & ~parent_core)[0]
    assert high_cells.size > 0
    high_cell = int(high_cells[0])
    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p111615_broad_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 820.0, dtype=np.float64)
    surface[active] = 1850.0
    surface[parent] = 2700.0
    surface[high_cell] = 3300.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1200.0, 0.0) + np.where(parent, 500.0, 0.0)
    score[high_cell] = 1900.0

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    halo = hierarchy["shoulder_halo"] > 0.0
    hierarchy_mask = hierarchy["crest"] | hierarchy["branch"] | hierarchy["foreland"]
    assert world.g("terrain.last_p111615_belt_morphology_refinement_accepted") == 1.0
    assert world.g("terrain.last_p111615_shoulder_halo_cell_count") > 0.0
    assert world.g("terrain.last_p111615_peak_hierarchy_cell_count_after") < world.g(
        "terrain.last_p111615_peak_hierarchy_cell_count_before"
    )
    assert world.g("terrain.last_p111615_high_peak_removed_cell_count") == 0.0
    assert world.g("terrain.last_p111615_halo_hierarchy_overlap_valid") == 1.0
    assert np.any(halo)
    assert not np.any(halo & hierarchy_mask)
    assert hierarchy["hierarchy"][high_cell] >= 2.0
    assert world.g("terrain.last_p111614_branch_width_ratio") <= world.g(
        "terrain.last_p111615_branch_width_ratio_before"
    )


def test_p111616_high_latitude_apron_splits_broad_orogen_without_removing_extreme_peaks():
    spec = get_preset("earthlike")
    spec.grid_cells = 5000
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    allowed = (
        (grid.lat > 58.0)
        & (grid.lat < 74.0)
        & (grid.lon > -145.0)
        & (grid.lon < 95.0)
    )
    start = int(np.argmin(np.where(
        allowed,
        (grid.lat - 64.0) ** 2 + (grid.lon + 130.0) ** 2,
        np.inf,
    )))
    target = int(np.argmin(np.where(
        allowed,
        (grid.lat - 66.0) ** 2 + (grid.lon - 80.0) ** 2,
        np.inf,
    )))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    parent_cells = np.asarray(path[::-1], dtype=np.int64)
    assert parent_cells.size >= 14
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=4)

    parent_core = module._dilate_mask(grid, parent, passes=1)
    parent_near = module._dilate_mask(grid, parent, passes=3)
    candidate_cells = np.where(active & parent_near & ~parent_core)[0]
    assert candidate_cells.size >= 3
    apron_seed = int(candidate_cells[0])
    halo_seed = int(candidate_cells[1])
    extreme_seed = int(candidate_cells[2])

    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p111616_highlat_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 720.0, dtype=np.float64)
    surface[active] = 2750.0
    surface[parent] = 3050.0
    surface[halo_seed] = 1900.0
    surface[apron_seed] = 2850.0
    surface[extreme_seed] = 3600.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1450.0, 0.0) + np.where(parent, 550.0, 0.0)
    score[halo_seed] = 1500.0
    score[apron_seed] = 1600.0
    score[extreme_seed] = 1900.0

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    apron = hierarchy["highland_apron"] > 0.0
    halo = hierarchy["shoulder_halo"] > 0.0
    hierarchy_mask = hierarchy["crest"] | hierarchy["branch"] | hierarchy["foreland"]
    assert world.g("terrain.last_p111616_highlat_morphology_refinement_accepted") == 1.0
    assert world.g("terrain.last_p111616_highlat_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111616_highland_apron_cell_count") > 0.0
    assert world.g("terrain.last_p111616_highlat_halo_extension_cell_count") >= 0.0
    assert world.g("terrain.last_p111616_peak_hierarchy_cell_count_after") < world.g(
        "terrain.last_p111616_peak_hierarchy_cell_count_before"
    )
    assert world.g("terrain.last_p111616_extreme_peak_reclassified_cell_count") == 0.0
    assert world.g("terrain.last_p111616_highland_apron_overlap_valid") == 1.0
    assert np.any(apron)
    assert not np.any(apron & halo)
    assert not np.any(apron & hierarchy_mask)
    assert hierarchy["hierarchy"][extreme_seed] >= 2.0


def test_p111636_polar_edge_repair_reclassifies_residual_polar_peak_width():
    spec = get_preset("earthlike")
    spec.grid_cells = 5000
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    module = TerrainModule()

    allowed = (
        (grid.lat > 58.0)
        & (grid.lat < 74.0)
        & (grid.lon > -145.0)
        & (grid.lon < 95.0)
    )
    start = int(np.argmin(np.where(
        allowed,
        (grid.lat - 64.0) ** 2 + (grid.lon + 130.0) ** 2,
        np.inf,
    )))
    target = int(np.argmin(np.where(
        allowed,
        (grid.lat - 66.0) ** 2 + (grid.lon - 80.0) ** 2,
        np.inf,
    )))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    parent_cells = np.asarray(path[::-1], dtype=np.int64)
    assert parent_cells.size >= 14
    parent = np.zeros(grid.n, dtype=bool)
    parent[parent_cells] = True
    active = module._dilate_mask(grid, parent, passes=4)
    extreme_seed = int(parent_cells[parent_cells.size // 2])

    world.set_g("terrain.enable_p111616_highlat_morphology_refinement", 0.0)
    world.set_g("terrain.enable_p111635_spine_linework_smoothing", 0.0)
    world.objects["tectonics.boundary_objects"] = [{
        "id": "test:p111636_polar_residual_parent_axis",
        "kind": "convergent_parent",
        "cells": parent_cells.astype(int).tolist(),
    }]
    surface = np.full(grid.n, 720.0, dtype=np.float64)
    surface[active] = 2650.0
    surface[parent] = 3050.0
    surface[extreme_seed] = 3600.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    score = np.where(active, 1450.0, 0.0) + np.where(parent, 550.0, 0.0)
    score[extreme_seed] = 2100.0

    hierarchy = module._p11168_parent_orogen_hierarchy_masks(
        world,
        surface,
        0.0,
        crust,
        active,
        score,
    )

    apron = hierarchy["highland_apron"] > 0.0
    hierarchy_mask = hierarchy["crest"] | hierarchy["branch"] | hierarchy["foreland"]
    assert world.g("terrain.last_p111616_highlat_morphology_refinement_accepted") == 0.0
    assert world.g("terrain.last_p111636_polar_edge_refinement_accepted") == 1.0
    assert world.g("terrain.last_p111636_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111636_reclassified_cell_count") > 0.0
    assert world.g("terrain.last_p111636_polar_reclassified_cell_count") > 0.0
    assert world.g("terrain.last_p111636_extreme_reclassified_cell_count") == 0.0
    assert world.g("terrain.last_p111636_polar_peak_area_fraction_after") < world.g(
        "terrain.last_p111636_polar_peak_area_fraction_before"
    )
    assert world.g("terrain.last_p160_polar_edge_spine_thinning_accepted") == 1.0
    assert world.g("terrain.last_p160_guard_reverted") == 0.0
    assert world.g("terrain.last_p160_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p160_demoted_cell_count") > 0.0
    assert world.g("terrain.last_p160_polar_demoted_cell_count") > 0.0
    assert world.g("terrain.last_p160_extreme_demoted_cell_count") == 0.0
    assert world.g("terrain.last_p160_polar_spine_area_fraction_after") < world.g(
        "terrain.last_p160_polar_spine_area_fraction_before"
    )
    assert world.g("terrain.last_p160_spine_component_count_after") <= world.g(
        "terrain.last_p160_spine_component_count_before"
    )
    assert world.g("terrain.last_p160_short_spine_component_count_after") <= world.g(
        "terrain.last_p160_short_spine_component_count_before"
    )
    assert np.any(apron)
    assert not np.any(apron & hierarchy_mask)
    assert hierarchy["spine"][extreme_seed] >= 2.0
    assert hierarchy["hierarchy"][extreme_seed] >= 2.0


def test_p11168_audit_reports_parent_orogen_hierarchy_metrics():
    spec = get_preset("earthlike")
    spec.grid_cells = 300
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[:5] = 3.0
    hierarchy[5:11] = 2.0
    hierarchy[11:19] = 1.0
    halo = np.zeros(grid.n, dtype=np.float64)
    halo[19:24] = 1.0
    apron = np.zeros(grid.n, dtype=np.float64)
    apron[24:29] = 1.0
    world.fields.update({
        "terrain.elevation_m": np.full(grid.n, 600.0, dtype=np.float64),
        "crust.type": np.full(grid.n, CONT, dtype=np.float64),
        "terrain.orogenic_parent_hierarchy": hierarchy,
        "terrain.orogenic_shoulder_halo": halo,
        "terrain.orogenic_highland_apron": apron,
    })
    world.set_g("terrain.last_p11168_parent_orogen_hierarchy_accepted", 1.0)
    world.set_g("terrain.last_p11168_parent_line_cell_count", 14.0)
    world.set_g("terrain.last_p11168_parent_line_area_fraction", 0.015)
    world.set_g("terrain.last_p11168_active_orogen_candidate_cell_count", 11.0)
    world.set_g("terrain.last_p11168_parent_overlap_fraction", 0.75)
    world.set_g("terrain.last_p146_parent_orogen_corridor_accepted", 1.0)
    world.set_g("terrain.last_p146_parent_orogen_corridor_object_count", 2.0)
    world.set_g("terrain.last_p146_parent_orogen_corridor_component_count", 2.0)
    world.set_g("terrain.last_p146_candidate_corridor_component_count", 5.0)
    world.set_g("terrain.last_p146_candidate_corridor_cell_count", 31.0)
    world.set_g("terrain.last_p146_retained_corridor_component_count", 2.0)
    world.set_g("terrain.last_p146_retained_corridor_cell_count", 21.0)
    world.set_g("terrain.last_p146_pruned_corridor_component_count", 3.0)
    world.set_g("terrain.last_p146_trunk_parent_overlap_fraction", 0.42)
    world.set_g("terrain.last_p146_reject_code", 0.0)
    world.set_g("terrain.last_p146_corridor_trunk_cell_count", 7.0)
    world.set_g("terrain.last_p146_corridor_branch_cell_count", 5.0)
    world.set_g("terrain.last_p146_corridor_foreland_cell_count", 9.0)
    world.set_g("terrain.last_p146_corridor_bridge_cell_count", 3.0)
    world.set_g("terrain.last_p146_corridor_bridge_path_count", 1.0)
    world.set_g("terrain.last_p147_boundary_group_corridor_accepted", 1.0)
    world.set_g("terrain.last_p147_parent_group_count", 4.0)
    world.set_g("terrain.last_p147_polyline_group_count", 3.0)
    world.set_g("terrain.last_p147_object_group_count", 1.0)
    world.set_g("terrain.last_p147_attempted_group_count", 4.0)
    world.set_g("terrain.last_p147_accepted_group_count", 2.0)
    world.set_g("terrain.last_p147_rejected_group_count", 2.0)
    world.set_g("terrain.last_p147_low_coverage_rejected_group_count", 1.0)
    world.set_g("terrain.last_p147_candidate_corridor_component_count", 8.0)
    world.set_g("terrain.last_p147_candidate_corridor_cell_count", 44.0)
    world.set_g("terrain.last_p147_retained_corridor_component_count", 3.0)
    world.set_g("terrain.last_p147_retained_corridor_cell_count", 30.0)
    world.set_g("terrain.last_p147_trunk_parent_overlap_mean", 0.35)
    world.set_g("terrain.last_p147_trunk_parent_overlap_min", 0.22)
    world.set_g("terrain.last_p147_side_to_trunk_ratio", 0.8)
    world.set_g("terrain.last_p147_aggregate_guard_rejected", 0.0)
    world.set_g("terrain.last_p147_corridor_object_count", 2.0)
    world.set_g("terrain.last_p147_corridor_cell_count", 30.0)
    world.set_g("terrain.last_p148_p147_trial_guard_enabled", 1.0)
    world.set_g("terrain.last_p148_p147_trial_guard_accepted", 0.0)
    world.set_g("terrain.last_p148_p147_trial_guard_rejected", 1.0)
    world.set_g("terrain.last_p148_p147_trial_reject_code", 3.0)
    world.set_g("terrain.last_p148_baseline_class_score", 0.55)
    world.set_g("terrain.last_p148_trial_class_score", 0.52)
    world.set_g("terrain.last_p148_baseline_class_small_area_fraction", 0.40)
    world.set_g("terrain.last_p148_trial_class_small_area_fraction", 0.45)
    world.set_g("terrain.last_p148_baseline_crest_component_count", 3.0)
    world.set_g("terrain.last_p148_trial_crest_component_count", 4.0)
    world.set_g("terrain.last_p148_baseline_branch_component_count", 2.0)
    world.set_g("terrain.last_p148_trial_branch_component_count", 3.0)
    world.set_g("terrain.last_p148_trial_added_corridor_cell_count", 11.0)
    world.set_g("terrain.last_p149_staged_promotion_enabled", 1.0)
    world.set_g("terrain.last_p149_staged_promotion_accepted", 1.0)
    world.set_g("terrain.last_p149_trunk_trial_accepted", 1.0)
    world.set_g("terrain.last_p149_trunk_trial_rejected", 0.0)
    world.set_g("terrain.last_p149_side_trial_accepted", 0.0)
    world.set_g("terrain.last_p149_side_trial_rejected", 1.0)
    world.set_g("terrain.last_p149_reject_code", 0.0)
    world.set_g("terrain.last_p149_trunk_trial_class_score", 0.61)
    world.set_g(
        "terrain.last_p149_trunk_trial_class_small_area_fraction",
        0.24,
    )
    world.set_g("terrain.last_p149_side_trial_class_score", 0.40)
    world.set_g(
        "terrain.last_p149_side_trial_class_small_area_fraction",
        0.55,
    )
    world.set_g("terrain.last_p149_committed_trunk_cell_count", 12.0)
    world.set_g("terrain.last_p149_committed_branch_cell_count", 0.0)
    world.set_g("terrain.last_p149_committed_foreland_cell_count", 0.0)
    world.set_g("terrain.last_p150_group_trunk_repair_enabled", 1.0)
    world.set_g("terrain.last_p150_group_trunk_trial_count", 4.0)
    world.set_g("terrain.last_p150_group_trunk_bridge_accepted_count", 2.0)
    world.set_g("terrain.last_p150_group_trunk_small_rejected_count", 1.0)
    world.set_g("terrain.last_p150_group_trunk_bridge_cell_count", 9.0)
    world.set_g("terrain.last_p150_group_trunk_bridge_path_count", 3.0)
    world.set_g("terrain.last_p150_group_trunk_class_score_before_mean", 0.22)
    world.set_g("terrain.last_p150_group_trunk_class_score_after_mean", 0.46)
    world.set_g("terrain.last_p150_group_trunk_class_small_before_max", 1.0)
    world.set_g("terrain.last_p150_group_trunk_class_small_after_max", 0.62)
    world.set_g("terrain.last_p150_component_preserving_cap_enabled", 1.0)
    world.set_g("terrain.last_p150_aggregate_trunk_bridge_accepted", 1.0)
    world.set_g("terrain.last_p150_aggregate_trunk_bridge_cell_count", 7.0)
    world.set_g("terrain.last_p150_aggregate_trunk_bridge_path_count", 2.0)
    world.set_g("terrain.last_p150_aggregate_trunk_class_score_before", 0.30)
    world.set_g("terrain.last_p150_aggregate_trunk_class_score_after", 0.57)
    world.set_g("terrain.last_p150_aggregate_trunk_class_small_before", 1.0)
    world.set_g("terrain.last_p150_aggregate_trunk_class_small_after", 0.48)
    world.set_g("terrain.last_p151_guarded_object_cap_enabled", 1.0)
    world.set_g("terrain.last_p151_trial_component_cap_evaluated", 1.0)
    world.set_g("terrain.last_p151_trial_component_cap_accepted", 1.0)
    world.set_g("terrain.last_p151_trial_component_cap_rejected", 0.0)
    world.set_g("terrain.last_p151_trial_selected_component_cap", 1.0)
    world.set_g("terrain.last_p151_trial_old_class_score", 0.21)
    world.set_g("terrain.last_p151_trial_component_class_score", 0.43)
    world.set_g(
        "terrain.last_p151_trial_old_class_small_area_fraction",
        1.0,
    )
    world.set_g(
        "terrain.last_p151_trial_component_class_small_area_fraction",
        0.27,
    )
    world.set_g("terrain.last_p151_trial_reject_code", 0.0)
    world.set_g("terrain.last_p151_final_component_cap_evaluated", 1.0)
    world.set_g("terrain.last_p151_final_component_cap_accepted", 0.0)
    world.set_g("terrain.last_p151_final_component_cap_rejected", 1.0)
    world.set_g("terrain.last_p151_final_old_class_score", 0.52)
    world.set_g("terrain.last_p151_final_component_class_score", 0.50)
    world.set_g(
        "terrain.last_p151_final_old_class_small_area_fraction",
        0.32,
    )
    world.set_g(
        "terrain.last_p151_final_component_class_small_area_fraction",
        0.51,
    )
    world.set_g("terrain.last_p151_final_reject_code", 2.0)
    world.set_g("terrain.last_p11168_crest_cell_count", 5.0)
    world.set_g("terrain.last_p11168_branch_cell_count", 6.0)
    world.set_g("terrain.last_p11168_foreland_cell_count", 8.0)
    world.set_g("terrain.last_p11168_crest_area_fraction", 0.01)
    world.set_g("terrain.last_p11168_branch_area_fraction", 0.012)
    world.set_g("terrain.last_p11168_foreland_area_fraction", 0.016)
    world.set_g("terrain.last_p11168_crest_component_count", 1.0)
    world.set_g("terrain.last_p11168_branch_component_count", 1.0)
    world.set_g("terrain.last_p11168_foreland_component_count", 2.0)
    world.set_g("terrain.last_p111610_hierarchy_refinement_accepted", 1.0)
    world.set_g("terrain.last_p111610_crest_components_before_refine", 3.0)
    world.set_g("terrain.last_p111610_crest_components_after_refine", 1.0)
    world.set_g("terrain.last_p111610_branch_components_before_refine", 2.0)
    world.set_g("terrain.last_p111610_branch_components_after_refine", 1.0)
    world.set_g("terrain.last_p111610_bridge_cell_count", 4.0)
    world.set_g("terrain.last_p111610_bridge_area_fraction", 0.008)
    world.set_g("terrain.last_p111613_hierarchy_geometry_refinement_accepted", 1.0)
    world.set_g("terrain.last_p111613_branch_components_before_refine", 4.0)
    world.set_g("terrain.last_p111613_branch_components_after_refine", 2.0)
    world.set_g("terrain.last_p111613_foreland_components_before_refine", 9.0)
    world.set_g("terrain.last_p111613_foreland_components_after_refine", 3.0)
    world.set_g("terrain.last_p111613_branch_extension_cell_count", 6.0)
    world.set_g("terrain.last_p111613_branch_extension_area_fraction", 0.012)
    world.set_g("terrain.last_p111613_foreland_removed_cell_count", 7.0)
    world.set_g("terrain.last_p111613_foreland_removed_area_fraction", 0.014)
    world.set_g("terrain.last_p111613_foreland_removed_component_count", 5.0)
    world.set_g("terrain.last_p111623_branch_spine_continuity_accepted", 1.0)
    world.set_g("terrain.last_p111623_branch_components_before_refine", 5.0)
    world.set_g("terrain.last_p111623_branch_components_after_refine", 4.0)
    world.set_g("terrain.last_p111623_peak_components_before_refine", 6.0)
    world.set_g("terrain.last_p111623_peak_components_after_refine", 3.0)
    world.set_g("terrain.last_p111623_bridge_cell_count", 8.0)
    world.set_g("terrain.last_p111623_bridge_area_fraction", 0.01)
    world.set_g("terrain.last_p111623_bridge_candidate_cell_count", 13.0)
    world.set_g("terrain.last_p111623_bridge_candidate_area_fraction", 0.018)
    world.set_g("terrain.last_p111623_path_count", 4.0)
    world.set_g("terrain.last_p111623_spine_components_before_refine", 7.0)
    world.set_g("terrain.last_p111623_spine_components_after_refine", 3.0)
    world.set_g("terrain.last_p111623_integrated_spine_cell_count", 9.0)
    world.set_g("terrain.last_p111623_integrated_spine_area_fraction", 0.011)
    world.set_g("terrain.last_p111623_integrated_spine_accepted", 1.0)
    world.set_g("terrain.last_p111614_spine_refinement_accepted", 1.0)
    world.set_g("terrain.last_p111614_crest_spine_cell_count", 4.0)
    world.set_g("terrain.last_p111614_branch_spine_cell_count", 3.0)
    world.set_g("terrain.last_p111614_crest_spine_area_fraction", 0.006)
    world.set_g("terrain.last_p111614_branch_spine_area_fraction", 0.004)
    world.set_g("terrain.last_p111614_crest_spine_component_count", 1.0)
    world.set_g("terrain.last_p111614_branch_spine_component_count", 2.0)
    world.set_g("terrain.last_p111614_crest_width_ratio", 2.5)
    world.set_g("terrain.last_p111614_branch_width_ratio", 3.0)
    world.set_g("terrain.last_p111614_spine_overlap_valid", 1.0)
    world.set_g("terrain.last_p111633c_final_spine_refinement_accepted", 1.0)
    world.set_g("terrain.last_p111633c_peak_components_before_refine", 9.0)
    world.set_g("terrain.last_p111633c_peak_components_after_refine", 6.0)
    world.set_g("terrain.last_p111633c_spine_less_components_before_refine", 5.0)
    world.set_g("terrain.last_p111633c_spine_less_components_after_refine", 1.0)
    world.set_g("terrain.last_p111633c_orphan_peak_reclassified_cell_count", 4.0)
    world.set_g("terrain.last_p111633c_orphan_high_peak_reclassified_cell_count", 1.0)
    world.set_g("terrain.last_p111633c_orphan_extreme_peak_reclassified_cell_count", 0.0)
    world.set_g("terrain.last_p111633c_final_spine_cell_count", 12.0)
    world.set_g("terrain.last_p111633c_final_spine_component_count", 4.0)
    world.set_g("terrain.last_p111635_spine_linework_smoothing_accepted", 1.0)
    world.set_g("terrain.last_p111635_spine_components_before_refine", 8.0)
    world.set_g("terrain.last_p111635_spine_components_after_refine", 5.0)
    world.set_g("terrain.last_p111635_short_spine_components_before_refine", 4.0)
    world.set_g("terrain.last_p111635_short_spine_components_after_refine", 2.0)
    world.set_g("terrain.last_p111635_branch_attachment_fraction_before", 0.25)
    world.set_g("terrain.last_p111635_branch_attachment_fraction_after", 0.75)
    world.set_g("terrain.last_p111635_bridge_cell_count", 6.0)
    world.set_g("terrain.last_p111635_bridge_area_fraction", 0.008)
    world.set_g("terrain.last_p111635_bridge_candidate_cell_count", 14.0)
    world.set_g("terrain.last_p111635_bridge_candidate_area_fraction", 0.021)
    world.set_g("terrain.last_p111635_path_count", 3.0)
    world.set_g("terrain.last_p111636_polar_edge_refinement_accepted", 1.0)
    world.set_g("terrain.last_p111636_polar_peak_area_fraction_before", 0.014)
    world.set_g("terrain.last_p111636_polar_peak_area_fraction_after", 0.009)
    world.set_g("terrain.last_p111636_edge_peak_component_count_before", 7.0)
    world.set_g("terrain.last_p111636_edge_peak_component_count_after", 4.0)
    world.set_g("terrain.last_p111636_candidate_cell_count", 12.0)
    world.set_g("terrain.last_p111636_candidate_area_fraction", 0.018)
    world.set_g("terrain.last_p111636_reclassified_cell_count", 9.0)
    world.set_g("terrain.last_p111636_reclassified_area_fraction", 0.013)
    world.set_g("terrain.last_p111636_polar_reclassified_cell_count", 6.0)
    world.set_g("terrain.last_p111636_edge_reclassified_cell_count", 3.0)
    world.set_g("terrain.last_p111636_extreme_reclassified_cell_count", 0.0)
    world.set_g(
        "terrain.last_p125_terminal_orogenic_semantic_land_consistency_accepted",
        1.0,
    )
    world.set_g("terrain.last_p125_guard_reverted", 0.0)
    world.set_g("terrain.last_p125_candidate_cell_count", 11.0)
    world.set_g("terrain.last_p125_cleared_hierarchy_cell_count", 4.0)
    world.set_g("terrain.last_p125_cleared_spine_cell_count", 2.0)
    world.set_g(
        "terrain.last_p126_terminal_orogenic_fringe_support_regularization_accepted",
        1.0,
    )
    world.set_g("terrain.last_p126_guard_reverted", 0.0)
    world.set_g("terrain.last_p126_candidate_cell_count", 17.0)
    world.set_g("terrain.last_p126_candidate_area_fraction", 0.021)
    world.set_g("terrain.last_p126_cleared_foreland_cell_count", 5.0)
    world.set_g("terrain.last_p126_cleared_halo_cell_count", 9.0)
    world.set_g("terrain.last_p126_cleared_apron_cell_count", 3.0)
    world.set_g("terrain.last_p126_foreland_component_count_before", 12.0)
    world.set_g("terrain.last_p126_foreland_component_count_after", 7.0)
    world.set_g("terrain.last_p126_halo_component_count_before", 20.0)
    world.set_g("terrain.last_p126_halo_component_count_after", 11.0)
    world.set_g("terrain.last_p126_apron_component_count_before", 6.0)
    world.set_g("terrain.last_p126_apron_component_count_after", 3.0)
    world.set_g("terrain.last_p126_tiny_foreland_component_count_before", 8.0)
    world.set_g("terrain.last_p126_tiny_foreland_component_count_after", 3.0)
    world.set_g("terrain.last_p126_tiny_halo_component_count_before", 13.0)
    world.set_g("terrain.last_p126_tiny_halo_component_count_after", 4.0)
    world.set_g("terrain.last_p126_tiny_apron_component_count_before", 4.0)
    world.set_g("terrain.last_p126_tiny_apron_component_count_after", 1.0)
    world.set_g("terrain.last_p126_peak_hierarchy_changed_cell_count", 0.0)
    world.set_g("terrain.last_p126_spine_changed_cell_count", 0.0)
    world.set_g(
        "terrain.last_p133_terminal_orogenic_fringe_band_ordering_accepted",
        1.0,
    )
    world.set_g("terrain.last_p133_guard_reverted", 0.0)
    world.set_g("terrain.last_p133_candidate_cell_count", 14.0)
    world.set_g("terrain.last_p133_candidate_area_fraction", 0.017)
    world.set_g("terrain.last_p133_cleared_foreland_cell_count", 4.0)
    world.set_g("terrain.last_p133_cleared_halo_cell_count", 6.0)
    world.set_g("terrain.last_p133_cleared_apron_cell_count", 4.0)
    world.set_g("terrain.last_p133_foreland_component_count_before", 7.0)
    world.set_g("terrain.last_p133_foreland_component_count_after", 5.0)
    world.set_g("terrain.last_p133_halo_component_count_before", 11.0)
    world.set_g("terrain.last_p133_halo_component_count_after", 8.0)
    world.set_g("terrain.last_p133_apron_component_count_before", 3.0)
    world.set_g("terrain.last_p133_apron_component_count_after", 2.0)
    world.set_g("terrain.last_p133_tiny_foreland_component_count_before", 3.0)
    world.set_g("terrain.last_p133_tiny_foreland_component_count_after", 1.0)
    world.set_g("terrain.last_p133_tiny_halo_component_count_before", 4.0)
    world.set_g("terrain.last_p133_tiny_halo_component_count_after", 2.0)
    world.set_g("terrain.last_p133_tiny_apron_component_count_before", 1.0)
    world.set_g("terrain.last_p133_tiny_apron_component_count_after", 0.0)
    world.set_g("terrain.last_p133_band_violation_cell_count", 9.0)
    world.set_g("terrain.last_p133_peak_hierarchy_changed_cell_count", 0.0)
    world.set_g("terrain.last_p133_spine_changed_cell_count", 0.0)
    world.set_g(
        "terrain.last_p137_terminal_orogenic_fringe_gap_consolidation_accepted",
        1.0,
    )
    world.set_g("terrain.last_p137_guard_reverted", 0.0)
    world.set_g("terrain.last_p137_candidate_cell_count", 6.0)
    world.set_g("terrain.last_p137_candidate_area_fraction", 0.009)
    world.set_g("terrain.last_p137_added_foreland_cell_count", 3.0)
    world.set_g("terrain.last_p137_added_foreland_area_fraction", 0.004)
    world.set_g("terrain.last_p137_added_halo_cell_count", 2.0)
    world.set_g("terrain.last_p137_added_halo_area_fraction", 0.003)
    world.set_g("terrain.last_p137_added_apron_cell_count", 1.0)
    world.set_g("terrain.last_p137_added_apron_area_fraction", 0.002)
    world.set_g("terrain.last_p137_foreland_gap_component_count", 2.0)
    world.set_g("terrain.last_p137_halo_gap_component_count", 1.0)
    world.set_g("terrain.last_p137_apron_gap_component_count", 1.0)
    world.set_g("terrain.last_p137_foreland_component_count_before", 5.0)
    world.set_g("terrain.last_p137_foreland_component_count_after", 3.0)
    world.set_g("terrain.last_p137_halo_component_count_before", 8.0)
    world.set_g("terrain.last_p137_halo_component_count_after", 7.0)
    world.set_g("terrain.last_p137_apron_component_count_before", 2.0)
    world.set_g("terrain.last_p137_apron_component_count_after", 1.0)
    world.set_g("terrain.last_p137_tiny_foreland_component_count_before", 1.0)
    world.set_g("terrain.last_p137_tiny_foreland_component_count_after", 0.0)
    world.set_g("terrain.last_p137_tiny_halo_component_count_before", 2.0)
    world.set_g("terrain.last_p137_tiny_halo_component_count_after", 1.0)
    world.set_g("terrain.last_p137_tiny_apron_component_count_before", 1.0)
    world.set_g("terrain.last_p137_tiny_apron_component_count_after", 0.0)
    world.set_g("terrain.last_p137_peak_hierarchy_changed_cell_count", 0.0)
    world.set_g("terrain.last_p137_spine_changed_cell_count", 0.0)
    world.set_g(
        "terrain.last_p138_terminal_orogenic_fringe_component_thickening_accepted",
        1.0,
    )
    world.set_g("terrain.last_p138_guard_reverted", 0.0)
    world.set_g("terrain.last_p138_candidate_cell_count", 9.0)
    world.set_g("terrain.last_p138_candidate_area_fraction", 0.014)
    world.set_g("terrain.last_p138_added_foreland_cell_count", 6.0)
    world.set_g("terrain.last_p138_added_halo_cell_count", 3.0)
    world.set_g("terrain.last_p138_added_apron_cell_count", 0.0)
    world.set_g("terrain.last_p138_grown_foreland_component_count", 2.0)
    world.set_g("terrain.last_p138_grown_halo_component_count", 1.0)
    world.set_g("terrain.last_p138_grown_apron_component_count", 0.0)
    world.set_g("terrain.last_p138_foreland_component_count_before", 3.0)
    world.set_g("terrain.last_p138_foreland_component_count_after", 3.0)
    world.set_g("terrain.last_p138_halo_component_count_before", 7.0)
    world.set_g("terrain.last_p138_halo_component_count_after", 7.0)
    world.set_g("terrain.last_p138_apron_component_count_before", 1.0)
    world.set_g("terrain.last_p138_apron_component_count_after", 1.0)
    world.set_g("terrain.last_p138_foreland_small_area_fraction_before", 0.64)
    world.set_g("terrain.last_p138_foreland_small_area_fraction_after", 0.31)
    world.set_g("terrain.last_p138_halo_small_area_fraction_before", 0.42)
    world.set_g("terrain.last_p138_halo_small_area_fraction_after", 0.22)
    world.set_g("terrain.last_p138_apron_small_area_fraction_before", 0.0)
    world.set_g("terrain.last_p138_apron_small_area_fraction_after", 0.0)
    world.set_g("terrain.last_p138_peak_hierarchy_changed_cell_count", 0.0)
    world.set_g("terrain.last_p138_spine_changed_cell_count", 0.0)
    world.set_g(
        "terrain.last_p139_terminal_branch_range_component_thickening_accepted",
        1.0,
    )
    world.set_g("terrain.last_p139_guard_reverted", 0.0)
    world.set_g("terrain.last_p139_candidate_cell_count", 1.0)
    world.set_g("terrain.last_p139_candidate_area_fraction", 0.002)
    world.set_g("terrain.last_p139_added_branch_cell_count", 1.0)
    world.set_g("terrain.last_p139_added_branch_area_fraction", 0.002)
    world.set_g("terrain.last_p139_cleared_halo_cell_count", 1.0)
    world.set_g("terrain.last_p139_cleared_apron_cell_count", 0.0)
    world.set_g("terrain.last_p139_grown_branch_component_count", 1.0)
    world.set_g("terrain.last_p139_branch_cell_count_before", 44.0)
    world.set_g("terrain.last_p139_branch_cell_count_after", 45.0)
    world.set_g("terrain.last_p139_branch_component_count_before", 7.0)
    world.set_g("terrain.last_p139_branch_component_count_after", 7.0)
    world.set_g("terrain.last_p139_branch_small_area_fraction_before", 0.36)
    world.set_g("terrain.last_p139_branch_small_area_fraction_after", 0.18)
    world.set_g("terrain.last_p139_halo_component_count_before", 5.0)
    world.set_g("terrain.last_p139_halo_component_count_after", 5.0)
    world.set_g("terrain.last_p139_apron_component_count_before", 2.0)
    world.set_g("terrain.last_p139_apron_component_count_after", 2.0)
    world.set_g("terrain.last_p139_crest_changed_cell_count", 0.0)
    world.set_g("terrain.last_p139_spine_changed_cell_count", 0.0)
    world.set_g("terrain.last_p111637a_spine_object_promotion_accepted", 1.0)
    world.set_g("terrain.last_p111637a_spine_components_before", 9.0)
    world.set_g("terrain.last_p111637a_spine_components_after", 6.0)
    world.set_g("terrain.last_p111637a_short_spine_components_before", 2.0)
    world.set_g("terrain.last_p111637a_short_spine_components_after", 0.0)
    world.set_g("terrain.last_p111637a_branch_attachment_fraction_before", 0.50)
    world.set_g("terrain.last_p111637a_branch_attachment_fraction_after", 0.90)
    world.set_g("terrain.last_p111637a_spine_top3_share_before", 0.62)
    world.set_g("terrain.last_p111637a_spine_top3_share_after", 0.84)
    world.set_g("terrain.last_p111637a_linework_score_before", 0.58)
    world.set_g("terrain.last_p111637a_linework_score_after", 0.81)
    world.set_g("terrain.last_p111637a_bridge_cell_count", 7.0)
    world.set_g("terrain.last_p111637a_bridge_area_fraction", 0.010)
    world.set_g("terrain.last_p111637a_candidate_cell_count", 21.0)
    world.set_g("terrain.last_p111637a_candidate_area_fraction", 0.031)
    world.set_g("terrain.last_p111637a_path_count", 4.0)
    world.set_g("terrain.last_p1145_whole_mask_spine_planner_accepted", 1.0)
    world.set_g("terrain.last_p1145_bridge_cell_count", 5.0)
    world.set_g("terrain.last_p1145_bridge_area_fraction", 0.007)
    world.set_g("terrain.last_p1145_candidate_cell_count", 18.0)
    world.set_g("terrain.last_p1145_candidate_area_fraction", 0.026)
    world.set_g("terrain.last_p1145_path_count", 2.0)
    world.set_g("terrain.last_p1145_linework_score_before", 0.81)
    world.set_g("terrain.last_p1145_linework_score_after", 0.91)
    world.set_g("terrain.last_p1145_spine_components_before", 6.0)
    world.set_g("terrain.last_p1145_spine_components_after", 4.0)
    world.set_g("terrain.last_p1145_short_spine_components_before", 1.0)
    world.set_g("terrain.last_p1145_short_spine_components_after", 0.0)
    world.set_g("terrain.last_p1145_branch_attachment_fraction_before", 0.90)
    world.set_g("terrain.last_p1145_branch_attachment_fraction_after", 1.0)
    world.set_g("terrain.last_p1145_spine_top3_share_before", 0.84)
    world.set_g("terrain.last_p1145_spine_top3_share_after", 0.92)
    world.set_g("terrain.last_p1146_reject_code", 15.0)
    world.set_g("terrain.last_p1146_support_component_count", 6.0)
    world.set_g("terrain.last_p1146_multi_spine_support_component_count", 2.0)
    world.set_g("terrain.last_p1146_attempted_path_count", 5.0)
    world.set_g("terrain.last_p1146_found_path_count", 4.0)
    world.set_g("terrain.last_p1146_bridge_path_count", 2.0)
    world.set_g("terrain.last_p1146_attempted_bridge_cell_count", 5.0)
    world.set_g("terrain.last_p1146_trial_linework_score", 0.91)
    world.set_g("terrain.last_p1146_trial_spine_components", 4.0)
    world.set_g("terrain.last_p1146_trial_short_spine_components", 0.0)
    world.set_g("terrain.last_p1146_trial_branch_attachment_fraction", 1.0)
    world.set_g("terrain.last_p1146_trial_spine_top3_share", 0.92)
    world.set_g("terrain.last_p1147_pair_option_count", 9.0)
    world.set_g("terrain.last_p1147_pair_selected_count", 2.0)
    world.set_g("terrain.last_p1147_pair_balance_rejected_count", 3.0)
    world.set_g("terrain.last_p1147_pair_profile_rejected_count", 4.0)
    world.set_g("terrain.last_p1148_terminal_proxy_enabled", 1.0)
    world.set_g("terrain.last_p1148_terminal_proxy_score_before", 0.88)
    world.set_g("terrain.last_p1148_terminal_proxy_score_after", 0.93)
    world.set_g("terrain.last_p1148_terminal_proxy_component_count_before", 5.0)
    world.set_g("terrain.last_p1148_terminal_proxy_component_count_after", 4.0)
    world.set_g("terrain.last_p1148_terminal_proxy_short_count_before", 1.0)
    world.set_g("terrain.last_p1148_terminal_proxy_short_count_after", 0.0)
    world.set_g("terrain.last_p111615_belt_morphology_refinement_accepted", 1.0)
    world.set_g("terrain.last_p111615_shoulder_halo_cell_count", 5.0)
    world.set_g("terrain.last_p111615_shoulder_halo_area_fraction", 0.007)
    world.set_g("terrain.last_p111615_shoulder_halo_component_count", 1.0)
    world.set_g("terrain.last_p111615_crest_pruned_cell_count", 2.0)
    world.set_g("terrain.last_p111615_branch_pruned_cell_count", 3.0)
    world.set_g("terrain.last_p111615_crest_width_ratio_before", 3.5)
    world.set_g("terrain.last_p111615_branch_width_ratio_before", 4.0)
    world.set_g("terrain.last_p111615_peak_hierarchy_cell_count_before", 16.0)
    world.set_g("terrain.last_p111615_peak_hierarchy_cell_count_after", 11.0)
    world.set_g("terrain.last_p111615_high_peak_removed_cell_count", 0.0)
    world.set_g("terrain.last_p111615_halo_hierarchy_overlap_valid", 1.0)
    world.set_g("terrain.last_p111616_highlat_morphology_refinement_accepted", 1.0)
    world.set_g("terrain.last_p111616_highlat_candidate_cell_count", 9.0)
    world.set_g("terrain.last_p111616_highland_apron_cell_count", 5.0)
    world.set_g("terrain.last_p111616_highland_apron_area_fraction", 0.009)
    world.set_g("terrain.last_p111616_highland_apron_component_count", 1.0)
    world.set_g("terrain.last_p111616_highlat_halo_extension_cell_count", 4.0)
    world.set_g("terrain.last_p111616_highlat_halo_extension_area_fraction", 0.006)
    world.set_g("terrain.last_p111616_crest_reclassified_cell_count", 3.0)
    world.set_g("terrain.last_p111616_branch_reclassified_cell_count", 6.0)
    world.set_g("terrain.last_p111616_peak_hierarchy_cell_count_before", 20.0)
    world.set_g("terrain.last_p111616_peak_hierarchy_cell_count_after", 11.0)
    world.set_g("terrain.last_p111616_extreme_peak_reclassified_cell_count", 0.0)
    world.set_g("terrain.last_p111616_highland_apron_overlap_valid", 1.0)

    metrics = terminal_audit_metrics(world)
    p11168 = metrics["p11168_orogenic_parent_hierarchy"]

    assert p11168["accepted"] is True
    assert p11168["parent_line_cell_count"] == 14.0
    assert p11168["parent_overlap_fraction"] == 0.75
    assert p11168["p146_parent_orogen_corridor_accepted"] is True
    assert p11168["p146_parent_orogen_corridor_object_count"] == 2.0
    assert p11168["p146_parent_orogen_corridor_component_count"] == 2.0
    assert p11168["p146_candidate_corridor_component_count"] == 5.0
    assert p11168["p146_candidate_corridor_cell_count"] == 31.0
    assert p11168["p146_retained_corridor_component_count"] == 2.0
    assert p11168["p146_retained_corridor_cell_count"] == 21.0
    assert p11168["p146_pruned_corridor_component_count"] == 3.0
    assert p11168["p146_trunk_parent_overlap_fraction"] == 0.42
    assert p11168["p146_reject_code"] == 0.0
    assert p11168["p146_corridor_trunk_cell_count"] == 7.0
    assert p11168["p146_corridor_branch_cell_count"] == 5.0
    assert p11168["p146_corridor_foreland_cell_count"] == 9.0
    assert p11168["p146_corridor_bridge_cell_count"] == 3.0
    assert p11168["p146_corridor_bridge_path_count"] == 1.0
    assert p11168["p147_boundary_group_corridor_accepted"] is True
    assert p11168["p147_parent_group_count"] == 4.0
    assert p11168["p147_polyline_group_count"] == 3.0
    assert p11168["p147_object_group_count"] == 1.0
    assert p11168["p147_attempted_group_count"] == 4.0
    assert p11168["p147_accepted_group_count"] == 2.0
    assert p11168["p147_rejected_group_count"] == 2.0
    assert p11168["p147_low_coverage_rejected_group_count"] == 1.0
    assert p11168["p147_candidate_corridor_component_count"] == 8.0
    assert p11168["p147_candidate_corridor_cell_count"] == 44.0
    assert p11168["p147_retained_corridor_component_count"] == 3.0
    assert p11168["p147_retained_corridor_cell_count"] == 30.0
    assert p11168["p147_trunk_parent_overlap_mean"] == 0.35
    assert p11168["p147_trunk_parent_overlap_min"] == 0.22
    assert p11168["p147_side_to_trunk_ratio"] == 0.8
    assert p11168["p147_aggregate_guard_rejected"] is False
    assert p11168["p147_corridor_object_count"] == 2.0
    assert p11168["p147_corridor_cell_count"] == 30.0
    assert p11168["p148_p147_trial_guard_enabled"] is True
    assert p11168["p148_p147_trial_guard_accepted"] is False
    assert p11168["p148_p147_trial_guard_rejected"] is True
    assert p11168["p148_p147_trial_reject_code"] == 3.0
    assert p11168["p148_baseline_class_score"] == 0.55
    assert p11168["p148_trial_class_score"] == 0.52
    assert p11168["p148_baseline_class_small_area_fraction"] == 0.40
    assert p11168["p148_trial_class_small_area_fraction"] == 0.45
    assert p11168["p148_baseline_crest_component_count"] == 3.0
    assert p11168["p148_trial_crest_component_count"] == 4.0
    assert p11168["p148_baseline_branch_component_count"] == 2.0
    assert p11168["p148_trial_branch_component_count"] == 3.0
    assert p11168["p148_trial_added_corridor_cell_count"] == 11.0
    assert p11168["p149_staged_promotion_enabled"] is True
    assert p11168["p149_staged_promotion_accepted"] is True
    assert p11168["p149_trunk_trial_accepted"] is True
    assert p11168["p149_trunk_trial_rejected"] is False
    assert p11168["p149_side_trial_accepted"] is False
    assert p11168["p149_side_trial_rejected"] is True
    assert p11168["p149_reject_code"] == 0.0
    assert p11168["p149_trunk_trial_class_score"] == 0.61
    assert p11168["p149_trunk_trial_class_small_area_fraction"] == 0.24
    assert p11168["p149_side_trial_class_score"] == 0.40
    assert p11168["p149_side_trial_class_small_area_fraction"] == 0.55
    assert p11168["p149_committed_trunk_cell_count"] == 12.0
    assert p11168["p149_committed_branch_cell_count"] == 0.0
    assert p11168["p149_committed_foreland_cell_count"] == 0.0
    assert p11168["p150_group_trunk_repair_enabled"] is True
    assert p11168["p150_group_trunk_trial_count"] == 4.0
    assert p11168["p150_group_trunk_bridge_accepted_count"] == 2.0
    assert p11168["p150_group_trunk_small_rejected_count"] == 1.0
    assert p11168["p150_group_trunk_bridge_cell_count"] == 9.0
    assert p11168["p150_group_trunk_bridge_path_count"] == 3.0
    assert p11168["p150_group_trunk_class_score_before_mean"] == 0.22
    assert p11168["p150_group_trunk_class_score_after_mean"] == 0.46
    assert p11168["p150_group_trunk_class_small_before_max"] == 1.0
    assert p11168["p150_group_trunk_class_small_after_max"] == 0.62
    assert p11168["p150_component_preserving_cap_enabled"] is True
    assert p11168["p150_aggregate_trunk_bridge_accepted"] is True
    assert p11168["p150_aggregate_trunk_bridge_cell_count"] == 7.0
    assert p11168["p150_aggregate_trunk_bridge_path_count"] == 2.0
    assert p11168["p150_aggregate_trunk_class_score_before"] == 0.30
    assert p11168["p150_aggregate_trunk_class_score_after"] == 0.57
    assert p11168["p150_aggregate_trunk_class_small_before"] == 1.0
    assert p11168["p150_aggregate_trunk_class_small_after"] == 0.48
    assert p11168["p151_guarded_object_cap_enabled"] is True
    assert p11168["p151_trial_component_cap_evaluated"] is True
    assert p11168["p151_trial_component_cap_accepted"] is True
    assert p11168["p151_trial_component_cap_rejected"] is False
    assert p11168["p151_trial_selected_component_cap"] is True
    assert p11168["p151_trial_old_class_score"] == 0.21
    assert p11168["p151_trial_component_class_score"] == 0.43
    assert p11168["p151_trial_old_class_small_area_fraction"] == 1.0
    assert p11168["p151_trial_component_class_small_area_fraction"] == 0.27
    assert p11168["p151_trial_reject_code"] == 0.0
    assert p11168["p151_final_component_cap_evaluated"] is True
    assert p11168["p151_final_component_cap_accepted"] is False
    assert p11168["p151_final_component_cap_rejected"] is True
    assert p11168["p151_final_old_class_score"] == 0.52
    assert p11168["p151_final_component_class_score"] == 0.50
    assert p11168["p151_final_old_class_small_area_fraction"] == 0.32
    assert p11168["p151_final_component_class_small_area_fraction"] == 0.51
    assert p11168["p151_final_reject_code"] == 2.0
    assert p11168["crest_component_count"] == 1.0
    assert p11168["field_crest_area_fraction"] > 0.0
    assert p11168["field_branch_area_fraction"] > 0.0
    assert p11168["field_foreland_area_fraction"] > 0.0
    assert p11168["field_crest_cell_count"] == 5
    assert p11168["field_branch_cell_count"] == 6
    assert p11168["field_foreland_cell_count"] == 8
    assert p11168["field_crest_component_count"] > 0
    assert p11168["field_branch_component_count"] > 0
    assert p11168["field_foreland_component_count"] > 0
    assert p11168["p111610_refinement_accepted"] is True
    assert p11168["p111610_crest_components_before_refine"] == 3.0
    assert p11168["p111610_crest_components_after_refine"] == 1.0
    assert p11168["p111610_bridge_cell_count"] == 4.0
    assert p11168["p111613_refinement_accepted"] is True
    assert p11168["p111613_branch_components_before_refine"] == 4.0
    assert p11168["p111613_branch_components_after_refine"] == 2.0
    assert p11168["p111613_foreland_components_before_refine"] == 9.0
    assert p11168["p111613_foreland_components_after_refine"] == 3.0
    assert p11168["p111613_branch_extension_cell_count"] == 6.0
    assert p11168["p111613_foreland_removed_cell_count"] == 7.0
    assert p11168["p111613_foreland_removed_component_count"] == 5.0
    assert p11168["p111623_refinement_accepted"] is True
    assert p11168["p111623_branch_components_before_refine"] == 5.0
    assert p11168["p111623_branch_components_after_refine"] == 4.0
    assert p11168["p111623_peak_components_before_refine"] == 6.0
    assert p11168["p111623_peak_components_after_refine"] == 3.0
    assert p11168["p111623_bridge_cell_count"] == 8.0
    assert p11168["p111623_bridge_area_fraction"] == 0.01
    assert p11168["p111623_bridge_candidate_cell_count"] == 13.0
    assert p11168["p111623_bridge_candidate_area_fraction"] == 0.018
    assert p11168["p111623_path_count"] == 4.0
    assert p11168["p111623_spine_components_before_refine"] == 7.0
    assert p11168["p111623_spine_components_after_refine"] == 3.0
    assert p11168["p111623_integrated_spine_cell_count"] == 9.0
    assert p11168["p111623_integrated_spine_area_fraction"] == 0.011
    assert p11168["p111623_integrated_spine_accepted"] is True
    assert p11168["p111614_spine_refinement_accepted"] is True
    assert p11168["p111614_crest_spine_cell_count"] == 4.0
    assert p11168["p111614_branch_spine_cell_count"] == 3.0
    assert p11168["p111614_crest_spine_component_count"] == 1.0
    assert p11168["p111614_branch_spine_component_count"] == 2.0
    assert p11168["p111614_crest_width_ratio"] == 2.5
    assert p11168["p111614_branch_width_ratio"] == 3.0
    assert p11168["p111614_spine_overlap_valid"] is True
    assert p11168["p111633c_final_spine_refinement_accepted"] is True
    assert p11168["p111633c_peak_components_before_refine"] == 9.0
    assert p11168["p111633c_peak_components_after_refine"] == 6.0
    assert p11168["p111633c_spine_less_components_before_refine"] == 5.0
    assert p11168["p111633c_spine_less_components_after_refine"] == 1.0
    assert p11168["p111633c_orphan_peak_reclassified_cell_count"] == 4.0
    assert p11168["p111633c_orphan_high_peak_reclassified_cell_count"] == 1.0
    assert p11168["p111633c_orphan_extreme_peak_reclassified_cell_count"] == 0.0
    assert p11168["p111633c_final_spine_cell_count"] == 12.0
    assert p11168["p111633c_final_spine_component_count"] == 4.0
    assert p11168["p111635_spine_linework_smoothing_accepted"] is True
    assert p11168["p111635_spine_components_before_refine"] == 8.0
    assert p11168["p111635_spine_components_after_refine"] == 5.0
    assert p11168["p111635_short_spine_components_before_refine"] == 4.0
    assert p11168["p111635_short_spine_components_after_refine"] == 2.0
    assert p11168["p111635_branch_attachment_fraction_before"] == 0.25
    assert p11168["p111635_branch_attachment_fraction_after"] == 0.75
    assert p11168["p111635_bridge_cell_count"] == 6.0
    assert p11168["p111635_bridge_candidate_cell_count"] == 14.0
    assert p11168["p111635_path_count"] == 3.0
    assert p11168["p111636_polar_edge_refinement_accepted"] is True
    assert p11168["p111636_polar_peak_area_fraction_before"] == 0.014
    assert p11168["p111636_polar_peak_area_fraction_after"] == 0.009
    assert p11168["p111636_edge_peak_component_count_before"] == 7.0
    assert p11168["p111636_edge_peak_component_count_after"] == 4.0
    assert p11168["p111636_candidate_cell_count"] == 12.0
    assert p11168["p111636_reclassified_cell_count"] == 9.0
    assert p11168["p111636_polar_reclassified_cell_count"] == 6.0
    assert p11168["p111636_edge_reclassified_cell_count"] == 3.0
    assert p11168["p111636_extreme_reclassified_cell_count"] == 0.0
    assert (
        p11168[
            "p125_terminal_orogenic_semantic_land_consistency_accepted"
        ] is True
    )
    assert p11168["p125_guard_reverted"] is False
    assert p11168["p125_candidate_cell_count"] == 11.0
    assert p11168["p125_cleared_hierarchy_cell_count"] == 4.0
    assert p11168["p125_cleared_spine_cell_count"] == 2.0
    assert (
        p11168[
            "p126_terminal_orogenic_fringe_support_regularization_accepted"
        ] is True
    )
    assert p11168["p126_guard_reverted"] is False
    assert p11168["p126_candidate_cell_count"] == 17.0
    assert p11168["p126_cleared_foreland_cell_count"] == 5.0
    assert p11168["p126_cleared_halo_cell_count"] == 9.0
    assert p11168["p126_cleared_apron_cell_count"] == 3.0
    assert p11168["p126_foreland_component_count_before"] == 12.0
    assert p11168["p126_foreland_component_count_after"] == 7.0
    assert p11168["p126_halo_component_count_before"] == 20.0
    assert p11168["p126_halo_component_count_after"] == 11.0
    assert p11168["p126_apron_component_count_before"] == 6.0
    assert p11168["p126_apron_component_count_after"] == 3.0
    assert p11168["p126_tiny_foreland_component_count_before"] == 8.0
    assert p11168["p126_tiny_foreland_component_count_after"] == 3.0
    assert p11168["p126_tiny_halo_component_count_before"] == 13.0
    assert p11168["p126_tiny_halo_component_count_after"] == 4.0
    assert p11168["p126_tiny_apron_component_count_before"] == 4.0
    assert p11168["p126_tiny_apron_component_count_after"] == 1.0
    assert p11168["p126_peak_hierarchy_changed_cell_count"] == 0.0
    assert p11168["p126_spine_changed_cell_count"] == 0.0
    assert (
        p11168[
            "p133_terminal_orogenic_fringe_band_ordering_accepted"
        ] is True
    )
    assert p11168["p133_guard_reverted"] is False
    assert p11168["p133_candidate_cell_count"] == 14.0
    assert p11168["p133_cleared_foreland_cell_count"] == 4.0
    assert p11168["p133_cleared_halo_cell_count"] == 6.0
    assert p11168["p133_cleared_apron_cell_count"] == 4.0
    assert p11168["p133_foreland_component_count_before"] == 7.0
    assert p11168["p133_foreland_component_count_after"] == 5.0
    assert p11168["p133_halo_component_count_before"] == 11.0
    assert p11168["p133_halo_component_count_after"] == 8.0
    assert p11168["p133_apron_component_count_before"] == 3.0
    assert p11168["p133_apron_component_count_after"] == 2.0
    assert p11168["p133_tiny_foreland_component_count_before"] == 3.0
    assert p11168["p133_tiny_foreland_component_count_after"] == 1.0
    assert p11168["p133_tiny_halo_component_count_before"] == 4.0
    assert p11168["p133_tiny_halo_component_count_after"] == 2.0
    assert p11168["p133_tiny_apron_component_count_before"] == 1.0
    assert p11168["p133_tiny_apron_component_count_after"] == 0.0
    assert p11168["p133_band_violation_cell_count"] == 9.0
    assert p11168["p133_peak_hierarchy_changed_cell_count"] == 0.0
    assert p11168["p133_spine_changed_cell_count"] == 0.0
    assert (
        p11168[
            "p137_terminal_orogenic_fringe_gap_consolidation_accepted"
        ] is True
    )
    assert p11168["p137_guard_reverted"] is False
    assert p11168["p137_candidate_cell_count"] == 6.0
    assert p11168["p137_added_foreland_cell_count"] == 3.0
    assert p11168["p137_added_halo_cell_count"] == 2.0
    assert p11168["p137_added_apron_cell_count"] == 1.0
    assert p11168["p137_foreland_gap_component_count"] == 2.0
    assert p11168["p137_halo_gap_component_count"] == 1.0
    assert p11168["p137_apron_gap_component_count"] == 1.0
    assert p11168["p137_foreland_component_count_before"] == 5.0
    assert p11168["p137_foreland_component_count_after"] == 3.0
    assert p11168["p137_halo_component_count_before"] == 8.0
    assert p11168["p137_halo_component_count_after"] == 7.0
    assert p11168["p137_apron_component_count_before"] == 2.0
    assert p11168["p137_apron_component_count_after"] == 1.0
    assert p11168["p137_tiny_foreland_component_count_before"] == 1.0
    assert p11168["p137_tiny_foreland_component_count_after"] == 0.0
    assert p11168["p137_tiny_halo_component_count_before"] == 2.0
    assert p11168["p137_tiny_halo_component_count_after"] == 1.0
    assert p11168["p137_tiny_apron_component_count_before"] == 1.0
    assert p11168["p137_tiny_apron_component_count_after"] == 0.0
    assert p11168["p137_peak_hierarchy_changed_cell_count"] == 0.0
    assert p11168["p137_spine_changed_cell_count"] == 0.0
    assert (
        p11168[
            "p138_terminal_orogenic_fringe_component_thickening_accepted"
        ] is True
    )
    assert p11168["p138_guard_reverted"] is False
    assert p11168["p138_candidate_cell_count"] == 9.0
    assert p11168["p138_added_foreland_cell_count"] == 6.0
    assert p11168["p138_added_halo_cell_count"] == 3.0
    assert p11168["p138_added_apron_cell_count"] == 0.0
    assert p11168["p138_grown_foreland_component_count"] == 2.0
    assert p11168["p138_grown_halo_component_count"] == 1.0
    assert p11168["p138_grown_apron_component_count"] == 0.0
    assert p11168["p138_foreland_component_count_before"] == 3.0
    assert p11168["p138_foreland_component_count_after"] == 3.0
    assert p11168["p138_halo_component_count_before"] == 7.0
    assert p11168["p138_halo_component_count_after"] == 7.0
    assert p11168["p138_apron_component_count_before"] == 1.0
    assert p11168["p138_apron_component_count_after"] == 1.0
    assert p11168["p138_foreland_small_area_fraction_before"] == 0.64
    assert p11168["p138_foreland_small_area_fraction_after"] == 0.31
    assert p11168["p138_halo_small_area_fraction_before"] == 0.42
    assert p11168["p138_halo_small_area_fraction_after"] == 0.22
    assert p11168["p138_apron_small_area_fraction_before"] == 0.0
    assert p11168["p138_apron_small_area_fraction_after"] == 0.0
    assert p11168["p138_peak_hierarchy_changed_cell_count"] == 0.0
    assert p11168["p138_spine_changed_cell_count"] == 0.0
    assert (
        p11168[
            "p139_terminal_branch_range_component_thickening_accepted"
        ] is True
    )
    assert p11168["p139_guard_reverted"] is False
    assert p11168["p139_candidate_cell_count"] == 1.0
    assert p11168["p139_added_branch_cell_count"] == 1.0
    assert p11168["p139_cleared_halo_cell_count"] == 1.0
    assert p11168["p139_cleared_apron_cell_count"] == 0.0
    assert p11168["p139_grown_branch_component_count"] == 1.0
    assert p11168["p139_branch_cell_count_before"] == 44.0
    assert p11168["p139_branch_cell_count_after"] == 45.0
    assert p11168["p139_branch_component_count_before"] == 7.0
    assert p11168["p139_branch_component_count_after"] == 7.0
    assert p11168["p139_branch_small_area_fraction_before"] == 0.36
    assert p11168["p139_branch_small_area_fraction_after"] == 0.18
    assert p11168["p139_halo_component_count_before"] == 5.0
    assert p11168["p139_halo_component_count_after"] == 5.0
    assert p11168["p139_apron_component_count_before"] == 2.0
    assert p11168["p139_apron_component_count_after"] == 2.0
    assert p11168["p139_crest_changed_cell_count"] == 0.0
    assert p11168["p139_spine_changed_cell_count"] == 0.0
    assert p11168["p111637a_spine_object_promotion_accepted"] is True
    assert p11168["p111637a_spine_components_before"] == 9.0
    assert p11168["p111637a_spine_components_after"] == 6.0
    assert p11168["p111637a_short_spine_components_before"] == 2.0
    assert p11168["p111637a_short_spine_components_after"] == 0.0
    assert p11168["p111637a_branch_attachment_fraction_before"] == 0.50
    assert p11168["p111637a_branch_attachment_fraction_after"] == 0.90
    assert p11168["p111637a_spine_top3_share_before"] == 0.62
    assert p11168["p111637a_spine_top3_share_after"] == 0.84
    assert p11168["p111637a_linework_score_before"] == 0.58
    assert p11168["p111637a_linework_score_after"] == 0.81
    assert p11168["p111637a_bridge_cell_count"] == 7.0
    assert p11168["p111637a_candidate_cell_count"] == 21.0
    assert p11168["p111637a_path_count"] == 4.0
    assert p11168["p1145_whole_mask_spine_planner_accepted"] is True
    assert p11168["p1145_bridge_cell_count"] == 5.0
    assert p11168["p1145_candidate_cell_count"] == 18.0
    assert p11168["p1145_path_count"] == 2.0
    assert p11168["p1145_linework_score_before"] == 0.81
    assert p11168["p1145_linework_score_after"] == 0.91
    assert p11168["p1145_spine_components_before"] == 6.0
    assert p11168["p1145_spine_components_after"] == 4.0
    assert p11168["p1145_short_spine_components_before"] == 1.0
    assert p11168["p1145_short_spine_components_after"] == 0.0
    assert p11168["p1145_branch_attachment_fraction_before"] == 0.90
    assert p11168["p1145_branch_attachment_fraction_after"] == 1.0
    assert p11168["p1145_spine_top3_share_before"] == 0.84
    assert p11168["p1145_spine_top3_share_after"] == 0.92
    assert p11168["p1146_reject_code"] == 15
    assert p11168["p1146_reject_reason"] == "no_material_improvement"
    assert p11168["p1146_support_component_count"] == 6.0
    assert p11168["p1146_multi_spine_support_component_count"] == 2.0
    assert p11168["p1146_attempted_path_count"] == 5.0
    assert p11168["p1146_found_path_count"] == 4.0
    assert p11168["p1146_bridge_path_count"] == 2.0
    assert p11168["p1146_attempted_bridge_cell_count"] == 5.0
    assert p11168["p1146_trial_linework_score"] == 0.91
    assert p11168["p1146_trial_spine_components"] == 4.0
    assert p11168["p1146_trial_short_spine_components"] == 0.0
    assert p11168["p1146_trial_branch_attachment_fraction"] == 1.0
    assert p11168["p1146_trial_spine_top3_share"] == 0.92
    assert p11168["p1147_pair_option_count"] == 9.0
    assert p11168["p1147_pair_selected_count"] == 2.0
    assert p11168["p1147_pair_balance_rejected_count"] == 3.0
    assert p11168["p1147_pair_profile_rejected_count"] == 4.0
    assert p11168["p1148_terminal_proxy_enabled"] is True
    assert p11168["p1148_terminal_proxy_score_before"] == 0.88
    assert p11168["p1148_terminal_proxy_score_after"] == 0.93
    assert p11168["p1148_terminal_proxy_component_count_before"] == 5.0
    assert p11168["p1148_terminal_proxy_component_count_after"] == 4.0
    assert p11168["p1148_terminal_proxy_short_count_before"] == 1.0
    assert p11168["p1148_terminal_proxy_short_count_after"] == 0.0
    assert p11168["p111615_refinement_accepted"] is True
    assert p11168["p111615_shoulder_halo_cell_count"] == 5.0
    assert p11168["p111615_field_shoulder_halo_area_fraction"] > 0.0
    assert p11168["p111615_crest_pruned_cell_count"] == 2.0
    assert p11168["p111615_branch_pruned_cell_count"] == 3.0
    assert p11168["p111615_crest_width_ratio_before"] == 3.5
    assert p11168["p111615_branch_width_ratio_before"] == 4.0
    assert p11168["p111615_peak_hierarchy_cell_count_before"] == 16.0
    assert p11168["p111615_peak_hierarchy_cell_count_after"] == 11.0
    assert p11168["p111615_high_peak_removed_cell_count"] == 0.0
    assert p11168["p111615_halo_hierarchy_overlap_valid"] is True
    assert p11168["p111616_refinement_accepted"] is True
    assert p11168["p111616_highlat_candidate_cell_count"] == 9.0
    assert p11168["p111616_highland_apron_cell_count"] == 5.0
    assert p11168["p111616_highland_apron_area_fraction"] == 0.009
    assert p11168["p111616_highland_apron_component_count"] == 1.0
    assert p11168["p111616_field_highland_apron_area_fraction"] > 0.0
    assert p11168["p111616_highlat_halo_extension_cell_count"] == 4.0
    assert p11168["p111616_crest_reclassified_cell_count"] == 3.0
    assert p11168["p111616_branch_reclassified_cell_count"] == 6.0
    assert p11168["p111616_peak_hierarchy_cell_count_before"] == 20.0
    assert p11168["p111616_peak_hierarchy_cell_count_after"] == 11.0
    assert p11168["p111616_extreme_peak_reclassified_cell_count"] == 0.0
    assert p11168["p111616_highland_apron_overlap_valid"] is True


def test_p111633_audit_reports_orogenic_belt_morphology_metrics():
    spec = get_preset("earthlike")
    spec.grid_cells = 300
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[:5] = 3.0
    hierarchy[5:11] = 2.0
    hierarchy[11:17] = 1.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[:3] = 3.0
    spine[5:8] = 2.0
    halo = np.zeros(grid.n, dtype=np.float64)
    halo[17:22] = 1.0
    apron = np.zeros(grid.n, dtype=np.float64)
    apron[22:25] = 1.0
    elev = np.full(grid.n, 600.0, dtype=np.float64)
    elev[:5] = 3200.0
    elev[5:8] = 2600.0
    elev[8:11] = 2100.0
    elev[11:17] = 900.0
    elev[17:22] = 1400.0
    elev[22:25] = 2700.0
    world.fields.update({
        "terrain.elevation_m": elev,
        "crust.type": np.full(grid.n, CONT, dtype=np.float64),
        "terrain.orogenic_parent_hierarchy": hierarchy,
        "terrain.orogenic_hierarchy_spine": spine,
        "terrain.orogenic_shoulder_halo": halo,
        "terrain.orogenic_highland_apron": apron,
    })
    world.set_g("terrain.last_p111633_belt_relief_response_used", 1.0)
    world.set_g("terrain.last_p111633_belt_relief_candidate_cell_count", 7.0)
    world.set_g("terrain.last_p111633_belt_relief_candidate_area_fraction", 0.021)
    world.set_g("terrain.last_p111633_belt_relief_bridge_cell_count", 3.0)
    world.set_g("terrain.last_p111633_belt_relief_bridge_area_fraction", 0.009)
    world.set_g("terrain.last_p111633_belt_relief_1800_component_count_before", 4.0)
    world.set_g("terrain.last_p111633_belt_relief_1800_component_count_after", 2.0)
    world.set_g("terrain.last_p111633_belt_relief_2400_component_count_before", 3.0)
    world.set_g("terrain.last_p111633_belt_relief_2400_component_count_after", 2.0)
    world.set_g("terrain.last_p111633_belt_relief_guard_reverted", 0.0)

    metrics = terminal_audit_metrics(world)
    belt = metrics["p111633_orogenic_belt_morphology"]

    assert belt["schema"] == "aevum.p111633_orogenic_belt_morphology.v1"
    assert belt["hierarchy"]["cell_count"] == 11
    assert belt["crest"]["cell_count"] == 5
    assert belt["branch"]["cell_count"] == 6
    assert belt["spine"]["cell_count"] == 6
    assert belt["relief_1800"]["cell_count"] == 11
    assert belt["relief_2400"]["cell_count"] == 8
    assert belt["relief_3000"]["cell_count"] == 5
    assert belt["relief_3000"]["coverage_fraction_of_peak_hierarchy"] > 0.0
    assert belt["relief_2400"]["spine_overlap_fraction_of_tier"] > 0.0
    assert belt["near_high_saddle_cell_count"] == 3
    assert belt["high_peak_fragmentation_pressure"] >= 0.0
    assert 0.0 <= belt["graded_belt_continuity_score"] <= 1.0
    assert belt["p111633_belt_relief_response_used"] is True
    assert belt["p111633_belt_relief_candidate_cell_count"] == 7.0
    assert belt["p111633_belt_relief_bridge_cell_count"] == 3.0
    assert belt["p111633_belt_relief_1800_component_count_before"] == 4.0
    assert belt["p111633_belt_relief_1800_component_count_after"] == 2.0
    assert belt["p111633_belt_relief_2400_component_count_before"] == 3.0
    assert belt["p111633_belt_relief_2400_component_count_after"] == 2.0
    assert belt["p111633_belt_relief_guard_reverted"] is False
    assert metrics["acceptance"]["has_p111633_orogenic_belt_morphology_metrics"]


def test_p111634_audit_reports_spine_aligned_elevation_morphology_metrics():
    spec = get_preset("earthlike")
    spec.grid_cells = 300
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[:7] = 3.0
    hierarchy[7:14] = 2.0
    hierarchy[14:20] = 1.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[:4] = 3.0
    spine[7:11] = 2.0
    halo = np.zeros(grid.n, dtype=np.float64)
    halo[20:26] = 1.0
    apron = np.zeros(grid.n, dtype=np.float64)
    apron[26:30] = 1.0
    elev = np.full(grid.n, 540.0, dtype=np.float64)
    elev[:4] = 3300.0
    elev[4:7] = 2700.0
    elev[7:11] = 2550.0
    elev[11:14] = 2100.0
    elev[14:20] = 820.0
    elev[20:26] = 1350.0
    elev[26:30] = 1800.0
    elev[40:42] = 3200.0
    hierarchy[40:42] = 2.0

    world.fields.update({
        "terrain.elevation_m": elev,
        "crust.type": np.full(grid.n, CONT, dtype=np.float64),
        "terrain.orogenic_parent_hierarchy": hierarchy,
        "terrain.orogenic_hierarchy_spine": spine,
        "terrain.orogenic_shoulder_halo": halo,
        "terrain.orogenic_highland_apron": apron,
    })
    world.set_g("terrain.last_p111634_spine_aligned_response_used", 1.0)
    world.set_g("terrain.last_p111634_axis_raise_candidate_cell_count", 9.0)
    world.set_g("terrain.last_p111634_axis_raise_candidate_area_fraction", 0.03)
    world.set_g("terrain.last_p111634_axis_raised_cell_count", 4.0)
    world.set_g("terrain.last_p111634_axis_raised_area_fraction", 0.012)
    world.set_g("terrain.last_p111634_offaxis_high_candidate_cell_count", 2.0)
    world.set_g("terrain.last_p111634_offaxis_high_candidate_area_fraction", 0.006)
    world.set_g("terrain.last_p111634_offaxis_high_softened_cell_count", 2.0)
    world.set_g("terrain.last_p111634_offaxis_high_softened_area_fraction", 0.006)
    world.set_g("terrain.last_p111634_high_component_count_before", 3.0)
    world.set_g("terrain.last_p111634_high_component_count_after", 2.0)
    world.set_g("terrain.last_p111634_high_near_spine_fraction_before", 0.55)
    world.set_g("terrain.last_p111634_high_near_spine_fraction_after", 0.82)
    world.set_g("terrain.last_p111634_guard_reverted", 0.0)
    world.set_g("terrain.last_p111637c_anti_raster_orogenic_expression_accepted", 1.0)
    world.set_g("terrain.last_p111637c_guard_reverted", 0.0)
    world.set_g("terrain.last_p111637c_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p111637c_candidate_cell_count", 14.0)
    world.set_g("terrain.last_p111637c_candidate_area_fraction", 0.041)
    world.set_g("terrain.last_p111637c_adjusted_cell_count", 9.0)
    world.set_g("terrain.last_p111637c_adjusted_area_fraction", 0.027)
    world.set_g("terrain.last_p111637c_roughness_before_m", 640.0)
    world.set_g("terrain.last_p111637c_roughness_after_m", 470.0)
    world.set_g("terrain.last_p111637c_mean_abs_delta_m", 85.0)
    world.set_g("terrain.last_p111637c_max_abs_delta_m", 160.0)
    world.set_g("terrain.last_p111637c_high_component_count_before", 5.0)
    world.set_g("terrain.last_p111637c_high_component_count_after", 4.0)
    world.set_g(
        "terrain.last_p123_spine_high_relief_continuity_expression_accepted",
        1.0,
    )
    world.set_g("terrain.last_p123_guard_reverted", 0.0)
    world.set_g("terrain.last_p123_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p123_candidate_cell_count", 7.0)
    world.set_g("terrain.last_p123_candidate_area_fraction", 0.014)
    world.set_g("terrain.last_p123_adjusted_cell_count", 5.0)
    world.set_g("terrain.last_p123_adjusted_area_fraction", 0.010)
    world.set_g("terrain.last_p123_spine_2400_component_count_before", 6.0)
    world.set_g("terrain.last_p123_spine_2400_component_count_after", 4.0)
    world.set_g("terrain.last_p123_spine_3000_component_count_before", 5.0)
    world.set_g("terrain.last_p123_spine_3000_component_count_after", 3.0)
    world.set_g("terrain.last_p123_spine_2400_coverage_before", 0.62)
    world.set_g("terrain.last_p123_spine_2400_coverage_after", 0.74)
    world.set_g("terrain.last_p123_high_component_count_before", 9.0)
    world.set_g("terrain.last_p123_high_component_count_after", 6.0)
    world.set_g("terrain.last_p123_mean_lift_m", 220.0)
    world.set_g("terrain.last_p123_max_lift_m", 520.0)
    world.set_g("terrain.last_p123_p90_land_relief_before_m", 2100.0)
    world.set_g("terrain.last_p123_p90_land_relief_after_m", 2140.0)
    world.set_g("terrain.last_p123_p98_land_relief_before_m", 3800.0)
    world.set_g("terrain.last_p123_p98_land_relief_after_m", 3880.0)
    world.set_g(
        "terrain.last_p134_terminal_spine_relief_gap_bridging_accepted",
        1.0,
    )
    world.set_g("terrain.last_p134_guard_reverted", 0.0)
    world.set_g("terrain.last_p134_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p134_candidate_cell_count", 8.0)
    world.set_g("terrain.last_p134_candidate_area_fraction", 0.016)
    world.set_g("terrain.last_p134_adjusted_cell_count", 6.0)
    world.set_g("terrain.last_p134_adjusted_area_fraction", 0.012)
    world.set_g("terrain.last_p134_bridge_gap_component_count_2400", 2.0)
    world.set_g("terrain.last_p134_bridge_gap_component_count_3000", 1.0)
    world.set_g("terrain.last_p134_bridge_gap_cell_count_2400", 5.0)
    world.set_g("terrain.last_p134_bridge_gap_cell_count_3000", 3.0)
    world.set_g("terrain.last_p134_spine_2400_component_count_before", 7.0)
    world.set_g("terrain.last_p134_spine_2400_component_count_after", 4.0)
    world.set_g("terrain.last_p134_spine_3000_component_count_before", 5.0)
    world.set_g("terrain.last_p134_spine_3000_component_count_after", 3.0)
    world.set_g("terrain.last_p134_spine_2400_coverage_before", 0.68)
    world.set_g("terrain.last_p134_spine_2400_coverage_after", 0.78)
    world.set_g("terrain.last_p134_spine_3000_coverage_before", 0.36)
    world.set_g("terrain.last_p134_spine_3000_coverage_after", 0.44)
    world.set_g("terrain.last_p134_high_component_count_before", 8.0)
    world.set_g("terrain.last_p134_high_component_count_after", 5.0)
    world.set_g("terrain.last_p134_mean_lift_m", 310.0)
    world.set_g("terrain.last_p134_max_lift_m", 720.0)
    world.set_g("terrain.last_p134_p90_land_relief_before_m", 2140.0)
    world.set_g("terrain.last_p134_p90_land_relief_after_m", 2190.0)
    world.set_g("terrain.last_p134_p98_land_relief_before_m", 3880.0)
    world.set_g("terrain.last_p134_p98_land_relief_after_m", 3940.0)
    world.set_g(
        "terrain.last_p135_terminal_crest_core_relief_gap_bridging_accepted",
        1.0,
    )
    world.set_g("terrain.last_p135_guard_reverted", 0.0)
    world.set_g("terrain.last_p135_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p135_candidate_cell_count", 4.0)
    world.set_g("terrain.last_p135_candidate_area_fraction", 0.008)
    world.set_g("terrain.last_p135_adjusted_cell_count", 4.0)
    world.set_g("terrain.last_p135_adjusted_area_fraction", 0.008)
    world.set_g("terrain.last_p135_gap_component_count", 1.0)
    world.set_g("terrain.last_p135_gap_cell_count", 4.0)
    world.set_g("terrain.last_p135_crest_3000_component_count_before", 7.0)
    world.set_g("terrain.last_p135_crest_3000_component_count_after", 6.0)
    world.set_g("terrain.last_p135_spine_3000_component_count_before", 8.0)
    world.set_g("terrain.last_p135_spine_3000_component_count_after", 7.0)
    world.set_g("terrain.last_p135_high_component_count_before", 8.0)
    world.set_g("terrain.last_p135_high_component_count_after", 7.0)
    world.set_g("terrain.last_p135_crest_3000_coverage_before", 0.42)
    world.set_g("terrain.last_p135_crest_3000_coverage_after", 0.51)
    world.set_g("terrain.last_p135_spine_3000_coverage_before", 0.36)
    world.set_g("terrain.last_p135_spine_3000_coverage_after", 0.44)
    world.set_g("terrain.last_p135_mean_lift_m", 420.0)
    world.set_g("terrain.last_p135_max_lift_m", 560.0)
    world.set_g("terrain.last_p135_p90_land_relief_before_m", 2190.0)
    world.set_g("terrain.last_p135_p90_land_relief_after_m", 2225.0)
    world.set_g("terrain.last_p135_p98_land_relief_before_m", 3940.0)
    world.set_g("terrain.last_p135_p98_land_relief_after_m", 3980.0)
    world.set_g(
        "terrain.last_p136_terminal_high_peak_speckle_rebalance_accepted",
        1.0,
    )
    world.set_g("terrain.last_p136_guard_reverted", 0.0)
    world.set_g("terrain.last_p136_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p136_candidate_cell_count", 6.0)
    world.set_g("terrain.last_p136_candidate_area_fraction", 0.012)
    world.set_g("terrain.last_p136_adjusted_cell_count", 6.0)
    world.set_g("terrain.last_p136_adjusted_area_fraction", 0.012)
    world.set_g("terrain.last_p136_candidate_component_count", 3.0)
    world.set_g("terrain.last_p136_lowered_component_count", 3.0)
    world.set_g("terrain.last_p136_high_component_count_before", 8.0)
    world.set_g("terrain.last_p136_high_component_count_after", 4.0)
    world.set_g("terrain.last_p136_spine_3000_component_count_before", 8.0)
    world.set_g("terrain.last_p136_spine_3000_component_count_after", 4.0)
    world.set_g("terrain.last_p136_crest_3000_component_count_before", 7.0)
    world.set_g("terrain.last_p136_crest_3000_component_count_after", 4.0)
    world.set_g("terrain.last_p136_spine_2400_component_count_before", 4.0)
    world.set_g("terrain.last_p136_spine_2400_component_count_after", 4.0)
    world.set_g("terrain.last_p136_high_cell_count_before", 36.0)
    world.set_g("terrain.last_p136_high_cell_count_after", 30.0)
    world.set_g("terrain.last_p136_mean_lower_m", 190.0)
    world.set_g("terrain.last_p136_max_lower_m", 260.0)
    world.set_g("terrain.last_p136_p90_land_relief_before_m", 2600.0)
    world.set_g("terrain.last_p136_p90_land_relief_after_m", 2600.0)
    world.set_g("terrain.last_p136_p98_land_relief_before_m", 3230.0)
    world.set_g("terrain.last_p136_p98_land_relief_after_m", 3230.0)

    metrics = terminal_audit_metrics(world)
    aligned = metrics["p111634_spine_aligned_elevation_morphology"]

    assert aligned["schema"] == (
        "aevum.p111634_spine_aligned_elevation_morphology.v1")
    assert aligned["peak_hierarchy_cell_count"] > 0
    assert aligned["spine_cell_count"] == 8
    assert aligned["high_cell_count"] > 0
    assert aligned["high_near_spine_fraction_d2"] >= 0.0
    assert aligned["high_far_from_spine_fraction"] >= 0.0
    assert 0.0 <= aligned["ridge_axis_relief_continuity_score"] <= 1.0
    assert 0.0 <= aligned["gradient_order_score"] <= 1.0
    assert aligned["crest_spine_median_relief_m"] >= (
        aligned["branch_spine_median_relief_m"])
    assert aligned["p111634_spine_aligned_response_used"] is True
    assert aligned["p111634_axis_raise_candidate_cell_count"] == 9.0
    assert aligned["p111634_axis_raised_cell_count"] == 4.0
    assert aligned["p111634_offaxis_high_candidate_cell_count"] == 2.0
    assert aligned["p111634_offaxis_high_softened_cell_count"] == 2.0
    assert aligned["p111634_high_component_count_before"] == 3.0
    assert aligned["p111634_high_component_count_after"] == 2.0
    assert aligned["p111634_high_near_spine_fraction_after"] == 0.82
    assert aligned["p111634_guard_reverted"] is False
    assert aligned["p111637c_anti_raster_orogenic_expression_accepted"] is True
    assert aligned["p111637c_guard_reverted"] is False
    assert aligned["p111637c_land_mask_preserved"] is True
    assert aligned["p111637c_candidate_cell_count"] == 14.0
    assert aligned["p111637c_candidate_area_fraction"] == 0.041
    assert aligned["p111637c_adjusted_cell_count"] == 9.0
    assert aligned["p111637c_adjusted_area_fraction"] == 0.027
    assert aligned["p111637c_roughness_before_m"] == 640.0
    assert aligned["p111637c_roughness_after_m"] == 470.0
    assert aligned["p111637c_mean_abs_delta_m"] == 85.0
    assert aligned["p111637c_max_abs_delta_m"] == 160.0
    assert aligned["p111637c_high_component_count_before"] == 5.0
    assert aligned["p111637c_high_component_count_after"] == 4.0
    assert aligned[
        "p123_spine_high_relief_continuity_expression_accepted"] is True
    assert aligned["p123_guard_reverted"] is False
    assert aligned["p123_land_mask_preserved"] is True
    assert aligned["p123_candidate_cell_count"] == 7.0
    assert aligned["p123_candidate_area_fraction"] == 0.014
    assert aligned["p123_adjusted_cell_count"] == 5.0
    assert aligned["p123_adjusted_area_fraction"] == 0.010
    assert aligned["p123_spine_2400_component_count_before"] == 6.0
    assert aligned["p123_spine_2400_component_count_after"] == 4.0
    assert aligned["p123_spine_3000_component_count_before"] == 5.0
    assert aligned["p123_spine_3000_component_count_after"] == 3.0
    assert aligned["p123_spine_2400_coverage_before"] == 0.62
    assert aligned["p123_spine_2400_coverage_after"] == 0.74
    assert aligned["p123_high_component_count_before"] == 9.0
    assert aligned["p123_high_component_count_after"] == 6.0
    assert aligned["p123_mean_lift_m"] == 220.0
    assert aligned["p123_max_lift_m"] == 520.0
    assert aligned["p123_p90_land_relief_before_m"] == 2100.0
    assert aligned["p123_p90_land_relief_after_m"] == 2140.0
    assert aligned["p123_p98_land_relief_before_m"] == 3800.0
    assert aligned["p123_p98_land_relief_after_m"] == 3880.0
    assert aligned[
        "p134_terminal_spine_relief_gap_bridging_accepted"] is True
    assert aligned["p134_guard_reverted"] is False
    assert aligned["p134_land_mask_preserved"] is True
    assert aligned["p134_candidate_cell_count"] == 8.0
    assert aligned["p134_adjusted_cell_count"] == 6.0
    assert aligned["p134_bridge_gap_component_count_2400"] == 2.0
    assert aligned["p134_bridge_gap_component_count_3000"] == 1.0
    assert aligned["p134_bridge_gap_cell_count_2400"] == 5.0
    assert aligned["p134_bridge_gap_cell_count_3000"] == 3.0
    assert aligned["p134_spine_2400_component_count_before"] == 7.0
    assert aligned["p134_spine_2400_component_count_after"] == 4.0
    assert aligned["p134_spine_3000_component_count_before"] == 5.0
    assert aligned["p134_spine_3000_component_count_after"] == 3.0
    assert aligned["p134_spine_2400_coverage_before"] == 0.68
    assert aligned["p134_spine_2400_coverage_after"] == 0.78
    assert aligned["p134_spine_3000_coverage_before"] == 0.36
    assert aligned["p134_spine_3000_coverage_after"] == 0.44
    assert aligned["p134_high_component_count_before"] == 8.0
    assert aligned["p134_high_component_count_after"] == 5.0
    assert aligned["p134_mean_lift_m"] == 310.0
    assert aligned["p134_max_lift_m"] == 720.0
    assert aligned["p134_p90_land_relief_before_m"] == 2140.0
    assert aligned["p134_p90_land_relief_after_m"] == 2190.0
    assert aligned["p134_p98_land_relief_before_m"] == 3880.0
    assert aligned["p134_p98_land_relief_after_m"] == 3940.0
    assert aligned[
        "p135_terminal_crest_core_relief_gap_bridging_accepted"] is True
    assert aligned["p135_guard_reverted"] is False
    assert aligned["p135_land_mask_preserved"] is True
    assert aligned["p135_candidate_cell_count"] == 4.0
    assert aligned["p135_adjusted_cell_count"] == 4.0
    assert aligned["p135_gap_component_count"] == 1.0
    assert aligned["p135_gap_cell_count"] == 4.0
    assert aligned["p135_crest_3000_component_count_before"] == 7.0
    assert aligned["p135_crest_3000_component_count_after"] == 6.0
    assert aligned["p135_spine_3000_component_count_before"] == 8.0
    assert aligned["p135_spine_3000_component_count_after"] == 7.0
    assert aligned["p135_high_component_count_before"] == 8.0
    assert aligned["p135_high_component_count_after"] == 7.0
    assert aligned["p135_crest_3000_coverage_before"] == 0.42
    assert aligned["p135_crest_3000_coverage_after"] == 0.51
    assert aligned["p135_spine_3000_coverage_before"] == 0.36
    assert aligned["p135_spine_3000_coverage_after"] == 0.44
    assert aligned["p135_mean_lift_m"] == 420.0
    assert aligned["p135_max_lift_m"] == 560.0
    assert aligned["p135_p90_land_relief_before_m"] == 2190.0
    assert aligned["p135_p90_land_relief_after_m"] == 2225.0
    assert aligned["p135_p98_land_relief_before_m"] == 3940.0
    assert aligned["p135_p98_land_relief_after_m"] == 3980.0
    assert aligned[
        "p136_terminal_high_peak_speckle_rebalance_accepted"] is True
    assert aligned["p136_guard_reverted"] is False
    assert aligned["p136_land_mask_preserved"] is True
    assert aligned["p136_candidate_cell_count"] == 6.0
    assert aligned["p136_adjusted_cell_count"] == 6.0
    assert aligned["p136_candidate_component_count"] == 3.0
    assert aligned["p136_lowered_component_count"] == 3.0
    assert aligned["p136_high_component_count_before"] == 8.0
    assert aligned["p136_high_component_count_after"] == 4.0
    assert aligned["p136_spine_3000_component_count_before"] == 8.0
    assert aligned["p136_spine_3000_component_count_after"] == 4.0
    assert aligned["p136_crest_3000_component_count_before"] == 7.0
    assert aligned["p136_crest_3000_component_count_after"] == 4.0
    assert aligned["p136_spine_2400_component_count_before"] == 4.0
    assert aligned["p136_spine_2400_component_count_after"] == 4.0
    assert aligned["p136_high_cell_count_before"] == 36.0
    assert aligned["p136_high_cell_count_after"] == 30.0
    assert aligned["p136_mean_lower_m"] == 190.0
    assert aligned["p136_max_lower_m"] == 260.0
    assert aligned["p136_p90_land_relief_before_m"] == 2600.0
    assert aligned["p136_p90_land_relief_after_m"] == 2600.0
    assert aligned["p136_p98_land_relief_before_m"] == 3230.0
    assert aligned["p136_p98_land_relief_after_m"] == 3230.0
    assert metrics["acceptance"][
        "has_p111634_spine_aligned_elevation_morphology_metrics"]


def test_p111635_audit_reports_orogenic_spine_linework_metrics():
    spec = get_preset("earthlike")
    spec.grid_cells = 300
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[:8] = 3.0
    hierarchy[8:18] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[:4] = 3.0
    spine[8:10] = 2.0
    spine[12:14] = 2.0
    elev = np.full(grid.n, 100.0, dtype=np.float64)
    elev[:8] = 3200.0
    elev[8:18] = 2300.0
    elev[30:33] = 1200.0
    mountain_inventory = np.zeros(grid.n, dtype=np.float64)
    mountain_inventory[:18] = 1.0
    world.fields.update({
        "terrain.elevation_m": elev,
        "crust.type": np.full(grid.n, CONT, dtype=np.float64),
        "terrain.orogenic_parent_hierarchy": hierarchy,
        "terrain.orogenic_hierarchy_spine": spine,
        "terrain.mountain_inventory": mountain_inventory,
    })
    world.set_g("terrain.last_p111635_spine_linework_smoothing_accepted", 1.0)
    world.set_g("terrain.last_p111635_spine_components_before_refine", 6.0)
    world.set_g("terrain.last_p111635_spine_components_after_refine", 3.0)
    world.set_g("terrain.last_p111635_short_spine_components_before_refine", 4.0)
    world.set_g("terrain.last_p111635_short_spine_components_after_refine", 1.0)
    world.set_g("terrain.last_p111635_branch_attachment_fraction_before", 0.2)
    world.set_g("terrain.last_p111635_branch_attachment_fraction_after", 0.8)
    world.set_g("terrain.last_p111635_bridge_cell_count", 5.0)
    world.set_g("terrain.last_p111635_bridge_area_fraction", 0.014)
    world.set_g("terrain.last_p111635_bridge_candidate_cell_count", 11.0)
    world.set_g("terrain.last_p111635_bridge_candidate_area_fraction", 0.032)
    world.set_g("terrain.last_p111635_path_count", 2.0)
    world.set_g("terrain.last_p111637a_spine_object_promotion_accepted", 1.0)
    world.set_g("terrain.last_p111637a_spine_components_before", 5.0)
    world.set_g("terrain.last_p111637a_spine_components_after", 3.0)
    world.set_g("terrain.last_p111637a_short_spine_components_before", 2.0)
    world.set_g("terrain.last_p111637a_short_spine_components_after", 0.0)
    world.set_g("terrain.last_p111637a_branch_attachment_fraction_before", 0.4)
    world.set_g("terrain.last_p111637a_branch_attachment_fraction_after", 0.9)
    world.set_g("terrain.last_p111637a_spine_top3_share_before", 0.55)
    world.set_g("terrain.last_p111637a_spine_top3_share_after", 0.82)
    world.set_g("terrain.last_p111637a_linework_score_before", 0.50)
    world.set_g("terrain.last_p111637a_linework_score_after", 0.78)
    world.set_g("terrain.last_p111637a_bridge_cell_count", 4.0)
    world.set_g("terrain.last_p111637a_bridge_area_fraction", 0.012)
    world.set_g("terrain.last_p111637a_candidate_cell_count", 10.0)
    world.set_g("terrain.last_p111637a_candidate_area_fraction", 0.030)
    world.set_g("terrain.last_p111637a_path_count", 2.0)
    world.set_g("terrain.last_p132_parent_anchor_spine_promotion_accepted", 1.0)
    world.set_g("terrain.last_p132_candidate_cell_count", 6.0)
    world.set_g("terrain.last_p132_candidate_area_fraction", 0.018)
    world.set_g("terrain.last_p132_promoted_spine_cell_count", 3.0)
    world.set_g("terrain.last_p132_promoted_spine_area_fraction", 0.009)
    world.set_g("terrain.last_p132_parent_aligned_spine_fraction_before", 0.42)
    world.set_g("terrain.last_p132_parent_aligned_spine_fraction_after", 0.66)
    world.set_g("terrain.last_p132_linework_score_before", 0.78)
    world.set_g("terrain.last_p132_linework_score_after", 0.81)
    world.set_g("terrain.last_p132_spine_components_before", 3.0)
    world.set_g("terrain.last_p132_spine_components_after", 3.0)
    world.set_g("terrain.last_p132_short_spine_components_before", 1.0)
    world.set_g("terrain.last_p132_short_spine_components_after", 0.0)
    world.set_g("terrain.last_p1148_terminal_proxy_enabled", 1.0)
    world.set_g("terrain.last_p1148_terminal_proxy_score_before", 0.86)
    world.set_g("terrain.last_p1148_terminal_proxy_score_after", 0.91)
    world.set_g("terrain.last_p1148_terminal_proxy_component_count_before", 4.0)
    world.set_g("terrain.last_p1148_terminal_proxy_component_count_after", 3.0)
    world.set_g("terrain.last_p1148_terminal_proxy_short_count_before", 1.0)
    world.set_g("terrain.last_p1148_terminal_proxy_short_count_after", 0.0)
    world.set_g("terrain.last_p1149_object_count", 3.0)
    world.set_g("terrain.last_p1149_system_count", 1.0)
    world.set_g("terrain.last_p1149_trunk_count", 1.0)
    world.set_g("terrain.last_p1149_branch_count", 1.0)
    world.set_g("terrain.last_p1149_fallback_trunk_count", 0.0)
    world.set_g("terrain.last_p1149_orphan_branch_count", 0.0)
    world.set_g("terrain.last_p1149_attached_branch_fraction", 1.0)
    world.set_g("terrain.last_p1149_mean_branch_count_per_system", 1.0)
    world.set_g("terrain.last_p1149_spine_cell_count", 8.0)
    world.set_g("terrain.last_p1149_trunk_cell_count", 4.0)
    world.set_g("terrain.last_p1149_branch_cell_count", 4.0)
    world.set_g("terrain.last_p1149_endpoint_count", 2.0)
    world.objects["terrain.orogenic_spine_objects"] = [
        {"kind": "orogen_spine_system", "cells": list(range(18))}
    ]
    world.set_g("terrain.last_p128_object_count", 2.0)
    world.set_g("terrain.last_p128_main_axis_count", 1.0)
    world.set_g("terrain.last_p128_branch_axis_count", 1.0)
    world.set_g("terrain.last_p128_fallback_main_axis_count", 0.0)
    world.set_g("terrain.last_p128_attached_branch_axis_fraction", 1.0)
    world.set_g("terrain.last_p128_source_spine_cell_count", 8.0)
    world.set_g("terrain.last_p128_axis_cell_count", 7.0)
    world.set_g("terrain.last_p128_axis_source_coverage_fraction", 0.875)
    world.set_g("terrain.last_p128_mean_path_coverage_fraction", 0.82)
    world.set_g("terrain.last_p128_mean_directness", 0.74)
    world.set_g("terrain.last_p128_mean_sinuosity", 1.35)
    world.set_g("terrain.last_p128_max_sinuosity", 1.62)
    world.set_g("terrain.last_p128_source_junction_cell_count", 3.0)
    world.set_g("terrain.last_p128_source_high_degree_cell_count", 1.0)
    world.set_g("terrain.last_p128_polyline_ready", 1.0)
    world.objects["terrain.orogenic_axis_polylines"] = [
        {"kind": "orogenic_main_crest_axis", "cells": [0, 1, 2, 3]},
        {"kind": "orogenic_branch_axis", "cells": [8, 9, 10]},
    ]
    world.set_g("terrain.last_p115_candidate_count", 2.0)
    world.set_g("terrain.last_p115_viable_candidate_count", 1.0)
    world.set_g("terrain.last_p115_trunk_bridge_attempt_count", 3.0)
    world.set_g("terrain.last_p115_branch_attachment_attempt_count", 2.0)
    world.set_g("terrain.last_p115_trunk_bridge_candidate_count", 1.0)
    world.set_g("terrain.last_p115_branch_attachment_candidate_count", 1.0)
    world.set_g("terrain.last_p115_rejected_proxy_count", 1.0)
    world.set_g("terrain.last_p115_multi_trunk_system_count", 1.0)
    world.set_g("terrain.last_p115_detached_branch_component_count", 1.0)
    world.set_g("terrain.last_p115_best_proxy_score_delta", 0.05)
    world.set_g("terrain.last_p115_best_component_delta", 1.0)
    world.set_g("terrain.last_p115_candidate_cell_count", 3.0)
    world.set_g("terrain.last_p115_viable_candidate_cell_count", 2.0)
    world.objects["terrain.orogenic_spine_repair_candidates"] = [
        {"kind": "orogenic_spine_repair_candidate", "cells": [1, 2]}
    ]
    world.set_g("terrain.last_p116_enabled", 1.0)
    world.set_g("terrain.last_p116_accepted", 1.0)
    world.set_g("terrain.last_p116_guard_reverted", 0.0)
    world.set_g("terrain.last_p116_input_candidate_count", 2.0)
    world.set_g("terrain.last_p116_input_viable_candidate_count", 1.0)
    world.set_g("terrain.last_p116_considered_candidate_count", 1.0)
    world.set_g("terrain.last_p116_selected_candidate_count", 1.0)
    world.set_g("terrain.last_p116_applied_candidate_count", 1.0)
    world.set_g("terrain.last_p116_rejected_candidate_count", 0.0)
    world.set_g("terrain.last_p116_rejected_overlap_count", 0.0)
    world.set_g("terrain.last_p116_rejected_support_count", 0.0)
    world.set_g("terrain.last_p116_rejected_profile_count", 0.0)
    world.set_g("terrain.last_p116_applied_cell_count", 2.0)
    world.set_g("terrain.last_p116_applied_area_fraction", 0.004)
    world.set_g("terrain.last_p116_area_budget_fraction", 0.004)
    world.set_g("terrain.last_p116_linework_score_before", 0.61)
    world.set_g("terrain.last_p116_linework_score_after", 0.70)
    world.set_g("terrain.last_p116_component_count_before", 5.0)
    world.set_g("terrain.last_p116_component_count_after", 4.0)
    world.set_g("terrain.last_p116_short_count_before", 1.0)
    world.set_g("terrain.last_p116_short_count_after", 0.0)
    world.set_g("terrain.last_p116_branch_attachment_fraction_before", 0.8)
    world.set_g("terrain.last_p116_branch_attachment_fraction_after", 1.0)
    world.set_g("terrain.last_p116_spine_top3_share_before", 0.7)
    world.set_g("terrain.last_p116_spine_top3_share_after", 0.82)
    world.set_g(
        "terrain.last_p124_orogenic_spine_geometry_regularization_accepted",
        1.0,
    )
    world.set_g("terrain.last_p124_guard_reverted", 0.0)
    world.set_g("terrain.last_p124_candidate_cell_count", 7.0)
    world.set_g("terrain.last_p124_candidate_area_fraction", 0.018)
    world.set_g("terrain.last_p124_added_cell_count", 3.0)
    world.set_g("terrain.last_p124_added_area_fraction", 0.007)
    world.set_g("terrain.last_p124_component_count_before", 5.0)
    world.set_g("terrain.last_p124_component_count_after", 4.0)
    world.set_g("terrain.last_p124_endpoint_count_before", 10.0)
    world.set_g("terrain.last_p124_endpoint_count_after", 8.0)
    world.set_g("terrain.last_p124_short_component_count_before", 2.0)
    world.set_g("terrain.last_p124_short_component_count_after", 1.0)
    world.set_g("terrain.last_p124_branch_attachment_fraction_before", 0.7)
    world.set_g("terrain.last_p124_branch_attachment_fraction_after", 0.85)
    world.set_g("terrain.last_p124_linework_score_before", 0.51)
    world.set_g("terrain.last_p124_linework_score_after", 0.62)
    world.set_g(
        "terrain.last_p127_terminal_orogenic_spine_node_thinning_accepted",
        1.0,
    )
    world.set_g("terrain.last_p127_guard_reverted", 0.0)
    world.set_g("terrain.last_p127_candidate_cell_count", 12.0)
    world.set_g("terrain.last_p127_candidate_area_fraction", 0.026)
    world.set_g("terrain.last_p127_demoted_cell_count", 4.0)
    world.set_g("terrain.last_p127_demoted_area_fraction", 0.009)
    world.set_g("terrain.last_p127_spine_cell_count_before", 31.0)
    world.set_g("terrain.last_p127_spine_cell_count_after", 27.0)
    world.set_g("terrain.last_p127_component_count_before", 4.0)
    world.set_g("terrain.last_p127_component_count_after", 4.0)
    world.set_g("terrain.last_p127_endpoint_count_before", 8.0)
    world.set_g("terrain.last_p127_endpoint_count_after", 8.0)
    world.set_g("terrain.last_p127_junction_count_before", 15.0)
    world.set_g("terrain.last_p127_junction_count_after", 10.0)
    world.set_g("terrain.last_p127_high_degree_count_before", 7.0)
    world.set_g("terrain.last_p127_high_degree_count_after", 2.0)
    world.set_g("terrain.last_p127_short_component_count_before", 1.0)
    world.set_g("terrain.last_p127_short_component_count_after", 1.0)
    world.set_g("terrain.last_p127_branch_attachment_fraction_before", 0.85)
    world.set_g("terrain.last_p127_branch_attachment_fraction_after", 0.85)
    world.set_g("terrain.last_p127_linework_score_before", 0.58)
    world.set_g("terrain.last_p127_linework_score_after", 0.64)
    world.set_g("terrain.last_p127_junction_pressure_before", 0.23)
    world.set_g("terrain.last_p127_junction_pressure_after", 0.07)
    world.set_g(
        "terrain.last_p140_terminal_branch_spine_component_thickening_accepted",
        1.0,
    )
    world.set_g("terrain.last_p140_guard_reverted", 0.0)
    world.set_g("terrain.last_p140_candidate_cell_count", 6.0)
    world.set_g("terrain.last_p140_candidate_area_fraction", 0.014)
    world.set_g("terrain.last_p140_promoted_branch_spine_cell_count", 6.0)
    world.set_g("terrain.last_p140_promoted_branch_spine_area_fraction", 0.014)
    world.set_g("terrain.last_p140_grown_branch_spine_component_count", 1.0)
    world.set_g("terrain.last_p140_branch_spine_cell_count_before", 30.0)
    world.set_g("terrain.last_p140_branch_spine_cell_count_after", 36.0)
    world.set_g("terrain.last_p140_branch_spine_component_count_before", 6.0)
    world.set_g("terrain.last_p140_branch_spine_component_count_after", 6.0)
    world.set_g(
        "terrain.last_p140_branch_spine_small_area_fraction_before",
        0.333,
    )
    world.set_g(
        "terrain.last_p140_branch_spine_small_area_fraction_after",
        0.178,
    )
    world.set_g("terrain.last_p140_all_spine_component_count_before", 2.0)
    world.set_g("terrain.last_p140_all_spine_component_count_after", 2.0)
    world.set_g("terrain.last_p140_junction_count_before", 32.0)
    world.set_g("terrain.last_p140_junction_count_after", 38.0)
    world.set_g("terrain.last_p140_high_degree_count_before", 7.0)
    world.set_g("terrain.last_p140_high_degree_count_after", 12.0)
    world.set_g("terrain.last_p140_junction_fraction_before", 0.421)
    world.set_g("terrain.last_p140_junction_fraction_after", 0.463)
    world.set_g("terrain.last_p140_high_degree_fraction_before", 0.092)
    world.set_g("terrain.last_p140_high_degree_fraction_after", 0.146)
    world.set_g("terrain.last_p140_branch_attachment_fraction_before", 1.0)
    world.set_g("terrain.last_p140_branch_attachment_fraction_after", 1.0)
    world.set_g("terrain.last_p140_hierarchy_changed_cell_count", 0.0)
    world.set_g("terrain.last_p140_crest_spine_changed_cell_count", 0.0)
    world.set_g(
        "terrain.last_p141_terminal_crest_spine_component_thickening_accepted",
        1.0,
    )
    world.set_g("terrain.last_p141_guard_reverted", 0.0)
    world.set_g("terrain.last_p141_candidate_cell_count", 3.0)
    world.set_g("terrain.last_p141_candidate_area_fraction", 0.007)
    world.set_g("terrain.last_p141_promoted_crest_spine_cell_count", 3.0)
    world.set_g("terrain.last_p141_promoted_crest_spine_area_fraction", 0.007)
    world.set_g("terrain.last_p141_connected_crest_spine_component_count", 1.0)
    world.set_g("terrain.last_p141_grown_crest_spine_component_count", 0.0)
    world.set_g("terrain.last_p141_crest_spine_cell_count_before", 49.0)
    world.set_g("terrain.last_p141_crest_spine_cell_count_after", 52.0)
    world.set_g("terrain.last_p141_crest_spine_component_count_before", 10.0)
    world.set_g("terrain.last_p141_crest_spine_component_count_after", 9.0)
    world.set_g(
        "terrain.last_p141_crest_spine_small_area_fraction_before",
        0.711,
    )
    world.set_g(
        "terrain.last_p141_crest_spine_small_area_fraction_after",
        0.558,
    )
    world.set_g("terrain.last_p141_all_spine_component_count_before", 5.0)
    world.set_g("terrain.last_p141_all_spine_component_count_after", 4.0)
    world.set_g("terrain.last_p141_junction_count_before", 38.0)
    world.set_g("terrain.last_p141_junction_count_after", 39.0)
    world.set_g("terrain.last_p141_high_degree_count_before", 12.0)
    world.set_g("terrain.last_p141_high_degree_count_after", 13.0)
    world.set_g("terrain.last_p141_junction_fraction_before", 0.463)
    world.set_g("terrain.last_p141_junction_fraction_after", 0.459)
    world.set_g("terrain.last_p141_high_degree_fraction_before", 0.146)
    world.set_g("terrain.last_p141_high_degree_fraction_after", 0.153)
    world.set_g("terrain.last_p141_branch_attachment_fraction_before", 1.0)
    world.set_g("terrain.last_p141_branch_attachment_fraction_after", 1.0)
    world.set_g("terrain.last_p141_hierarchy_changed_cell_count", 0.0)
    world.set_g("terrain.last_p141_branch_spine_changed_cell_count", 0.0)
    world.set_g(
        "terrain.last_p164_terminal_orogenic_axis_skeleton_simplification_accepted",
        1.0,
    )
    world.set_g("terrain.last_p164_guard_reverted", 0.0)
    world.set_g("terrain.last_p164_candidate_cell_count", 9.0)
    world.set_g("terrain.last_p164_candidate_area_fraction", 0.018)
    world.set_g("terrain.last_p164_demoted_cell_count", 4.0)
    world.set_g("terrain.last_p164_demoted_area_fraction", 0.008)
    world.set_g("terrain.last_p164_off_axis_candidate_cell_count", 12.0)
    world.set_g("terrain.last_p164_spine_cell_count_before", 36.0)
    world.set_g("terrain.last_p164_spine_cell_count_after", 32.0)
    world.set_g("terrain.last_p164_component_count_before", 3.0)
    world.set_g("terrain.last_p164_component_count_after", 3.0)
    world.set_g("terrain.last_p164_short_component_count_before", 1.0)
    world.set_g("terrain.last_p164_short_component_count_after", 1.0)
    world.set_g("terrain.last_p164_endpoint_count_before", 8.0)
    world.set_g("terrain.last_p164_endpoint_count_after", 8.0)
    world.set_g("terrain.last_p164_junction_count_before", 14.0)
    world.set_g("terrain.last_p164_junction_count_after", 9.0)
    world.set_g("terrain.last_p164_high_degree_count_before", 4.0)
    world.set_g("terrain.last_p164_high_degree_count_after", 2.0)
    world.set_g("terrain.last_p164_junction_fraction_before", 0.389)
    world.set_g("terrain.last_p164_junction_fraction_after", 0.281)
    world.set_g("terrain.last_p164_high_degree_fraction_before", 0.111)
    world.set_g("terrain.last_p164_high_degree_fraction_after", 0.062)
    world.set_g("terrain.last_p164_branch_attachment_fraction_before", 1.0)
    world.set_g("terrain.last_p164_branch_attachment_fraction_after", 1.0)
    world.set_g("terrain.last_p164_linework_score_before", 0.72)
    world.set_g("terrain.last_p164_linework_score_after", 0.73)
    world.set_g("terrain.last_p164_axis_score_before", 0.61)
    world.set_g("terrain.last_p164_axis_score_after", 0.68)
    world.set_g("terrain.last_p164_protected_extreme_demoted_cell_count", 0.0)
    world.set_g("terrain.last_p164_reject_code", 0.0)
    world.set_g("terrain.last_p142_class_aware_combined_score", 0.96)
    world.set_g("terrain.last_p142_class_aware_class_score", 0.58)
    world.set_g("terrain.last_p142_class_aware_crest_small_area_fraction", 0.71)
    world.set_g("terrain.last_p142_class_aware_branch_small_area_fraction", 0.65)
    world.set_g("terrain.last_p142_class_aware_class_small_area_fraction", 0.71)
    world.set_g("terrain.last_p142_class_aware_crest_component_count", 10.0)
    world.set_g("terrain.last_p142_class_aware_branch_component_count", 9.0)
    world.set_g("terrain.last_p142_class_aware_blind_spot", 1.0)
    world.set_g("terrain.last_p143_class_repair_needed", 1.0)
    world.set_g("terrain.last_p143_class_profile_improved", 1.0)
    world.set_g("terrain.last_p143_class_score_before", 0.58)
    world.set_g("terrain.last_p143_class_score_after", 0.66)
    world.set_g("terrain.last_p143_class_small_area_fraction_before", 0.71)
    world.set_g("terrain.last_p143_class_small_area_fraction_after", 0.55)
    world.set_g("terrain.last_p143_crest_small_area_fraction_before", 0.71)
    world.set_g("terrain.last_p143_crest_small_area_fraction_after", 0.55)
    world.set_g("terrain.last_p143_branch_small_area_fraction_before", 0.65)
    world.set_g("terrain.last_p143_branch_small_area_fraction_after", 0.50)
    world.set_g("terrain.last_p144_class_path_option_count", 7.0)
    world.set_g("terrain.last_p144_class_path_selected_count", 2.0)
    world.set_g("terrain.last_p144_crest_path_option_count", 4.0)
    world.set_g("terrain.last_p144_branch_path_option_count", 3.0)
    world.set_g("terrain.last_p144_class_attempted_path_count", 9.0)
    world.set_g("terrain.last_p144_class_found_path_count", 8.0)
    world.set_g("terrain.last_p144_crest_promoted_spine_cell_count", 5.0)
    world.set_g("terrain.last_p144_branch_promoted_spine_cell_count", 4.0)
    world.set_g("terrain.last_p144_class_promoted_spine_cell_count", 9.0)
    world.set_g("terrain.last_p144_class_path_profile_rejected_count", 2.0)
    world.set_g(
        "terrain.last_p145_class_hierarchy_component_consolidation_accepted",
        1.0,
    )
    world.set_g("terrain.last_p145_removed_crest_component_count", 2.0)
    world.set_g("terrain.last_p145_removed_branch_component_count", 3.0)
    world.set_g("terrain.last_p145_removed_crest_cell_count", 5.0)
    world.set_g("terrain.last_p145_removed_branch_cell_count", 7.0)
    world.set_g("terrain.last_p145_removed_crest_area_fraction", 0.011)
    world.set_g("terrain.last_p145_removed_branch_area_fraction", 0.014)
    world.set_g("terrain.last_p145_added_crest_cell_count", 6.0)
    world.set_g("terrain.last_p145_added_branch_cell_count", 8.0)
    world.set_g("terrain.last_p145_added_crest_area_fraction", 0.016)
    world.set_g("terrain.last_p145_added_branch_area_fraction", 0.019)
    world.set_g("terrain.last_p145_bridge_path_count", 3.0)
    world.set_g("terrain.last_p145_class_score_before", 0.45)
    world.set_g("terrain.last_p145_class_score_after", 0.62)
    world.set_g("terrain.last_p145_class_small_area_fraction_before", 1.0)
    world.set_g("terrain.last_p145_class_small_area_fraction_after", 0.69)
    world.set_g("terrain.last_p145_crest_small_area_fraction_before", 0.80)
    world.set_g("terrain.last_p145_crest_small_area_fraction_after", 0.52)
    world.set_g("terrain.last_p145_branch_small_area_fraction_before", 1.0)
    world.set_g("terrain.last_p145_branch_small_area_fraction_after", 0.69)
    world.set_g(
        "terrain.last_p118_non_orogenic_midland_stripe_suppression_accepted",
        1.0,
    )
    world.set_g("terrain.last_p118_guard_reverted", 0.0)
    world.set_g("terrain.last_p118_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p118_candidate_cell_count", 4.0)
    world.set_g("terrain.last_p118_candidate_area_fraction", 0.012)
    world.set_g("terrain.last_p118_adjusted_cell_count", 3.0)
    world.set_g("terrain.last_p118_adjusted_area_fraction", 0.009)
    world.set_g("terrain.last_p118_unsupported_midland_area_fraction_before", 0.018)
    world.set_g("terrain.last_p118_unsupported_midland_area_fraction_after", 0.006)
    world.set_g(
        "terrain.last_p118_non_orogenic_highland_500_area_fraction_before",
        0.024,
    )
    world.set_g(
        "terrain.last_p118_non_orogenic_highland_500_area_fraction_after",
        0.014,
    )
    world.set_g("terrain.last_p118_mean_lowering_m", 220.0)
    world.set_g("terrain.last_p118_max_lowering_m", 410.0)
    world.set_g("terrain.last_p118_mean_land_relief_before_m", 720.0)
    world.set_g("terrain.last_p118_mean_land_relief_after_m", 690.0)
    world.set_g("terrain.last_p118_p50_land_relief_before_m", 430.0)
    world.set_g("terrain.last_p118_p50_land_relief_after_m", 390.0)
    world.set_g("terrain.last_p118_p90_land_relief_before_m", 1850.0)
    world.set_g("terrain.last_p118_p90_land_relief_after_m", 1810.0)
    world.set_g(
        "terrain.last_p119_non_orogenic_interior_anti_raster_smoothing_accepted",
        1.0,
    )
    world.set_g("terrain.last_p119_guard_reverted", 0.0)
    world.set_g("terrain.last_p119_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p119_candidate_cell_count", 12.0)
    world.set_g("terrain.last_p119_candidate_area_fraction", 0.030)
    world.set_g("terrain.last_p119_selected_cell_count", 8.0)
    world.set_g("terrain.last_p119_selected_area_fraction", 0.020)
    world.set_g("terrain.last_p119_adjusted_cell_count", 7.0)
    world.set_g("terrain.last_p119_adjusted_area_fraction", 0.018)
    world.set_g("terrain.last_p119_compatible_edge_count", 19.0)
    world.set_g("terrain.last_p119_roughness_before_m", 310.0)
    world.set_g("terrain.last_p119_roughness_after_m", 240.0)
    world.set_g(
        "terrain.last_p119_non_orogenic_highland_500_area_fraction_before",
        0.064,
    )
    world.set_g(
        "terrain.last_p119_non_orogenic_highland_500_area_fraction_after",
        0.060,
    )
    world.set_g("terrain.last_p119_mean_abs_delta_m", 54.0)
    world.set_g("terrain.last_p119_max_abs_delta_m", 105.0)
    world.set_g("terrain.last_p119_p50_land_relief_before_m", 390.0)
    world.set_g("terrain.last_p119_p50_land_relief_after_m", 386.0)
    world.set_g("terrain.last_p119_p90_land_relief_before_m", 1810.0)
    world.set_g("terrain.last_p119_p90_land_relief_after_m", 1810.0)
    world.set_g(
        "terrain.last_p120_continental_semantic_geometry_repair_accepted",
        1.0,
    )
    world.set_g("terrain.last_p120_guard_reverted", 0.0)
    world.set_g("terrain.last_p120_candidate_cell_count", 16.0)
    world.set_g("terrain.last_p120_candidate_area_fraction", 0.042)
    world.set_g("terrain.last_p120_detail_changed_area_fraction", 0.012)
    world.set_g("terrain.last_p120_internal_changed_area_fraction", 0.014)
    world.set_g("terrain.last_p120_terrain_changed_area_fraction", 0.004)
    world.set_g("terrain.last_p120_detail_elongated_before", 0.18)
    world.set_g("terrain.last_p120_detail_elongated_after", 0.11)
    world.set_g("terrain.last_p120_internal_elongated_before", 0.24)
    world.set_g("terrain.last_p120_internal_elongated_after", 0.15)
    world.set_g("terrain.last_p120_terrain_elongated_before", 0.10)
    world.set_g("terrain.last_p120_terrain_elongated_after", 0.08)
    world.set_g("terrain.last_p120_detail_tiny_before", 0.020)
    world.set_g("terrain.last_p120_detail_tiny_after", 0.018)
    world.set_g("terrain.last_p120_internal_tiny_before", 0.030)
    world.set_g("terrain.last_p120_internal_tiny_after", 0.026)
    world.set_g("terrain.last_p120_terrain_tiny_before", 0.012)
    world.set_g("terrain.last_p120_terrain_tiny_after", 0.012)
    world.set_g(
        "terrain.last_p122_orogenic_detail_footprint_compaction_accepted",
        1.0,
    )
    world.set_g("terrain.last_p122_guard_reverted", 0.0)
    world.set_g("terrain.last_p122_candidate_cell_count", 9.0)
    world.set_g("terrain.last_p122_candidate_area_fraction", 0.021)
    world.set_g("terrain.last_p122_detail_changed_area_fraction", 0.017)
    world.set_g("terrain.last_p122_internal_changed_area_fraction", 0.015)
    world.set_g("terrain.last_p122_terrain_changed_area_fraction", 0.012)
    world.set_g("terrain.last_p122_detail_orogen_area_fraction_before", 0.110)
    world.set_g("terrain.last_p122_detail_orogen_area_fraction_after", 0.093)
    world.set_g(
        "terrain.last_p122_low_relief_orogen_area_fraction_before",
        0.078,
    )
    world.set_g(
        "terrain.last_p122_low_relief_orogen_area_fraction_after",
        0.055,
    )
    world.set_g("terrain.last_p122_noncore_orogen_area_fraction_before", 0.045)
    world.set_g("terrain.last_p122_noncore_orogen_area_fraction_after", 0.022)
    world.set_g("terrain.last_p122_preserved_core_orogen_cell_count", 11.0)
    world.set_g(
        "terrain.last_p121_semantic_region_elevation_response_accepted",
        1.0,
    )
    world.set_g("terrain.last_p121_guard_reverted", 0.0)
    world.set_g("terrain.last_p121_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p121_candidate_cell_count", 22.0)
    world.set_g("terrain.last_p121_candidate_area_fraction", 0.052)
    world.set_g("terrain.last_p121_selected_cell_count", 17.0)
    world.set_g("terrain.last_p121_selected_area_fraction", 0.041)
    world.set_g("terrain.last_p121_adjusted_cell_count", 15.0)
    world.set_g("terrain.last_p121_adjusted_area_fraction", 0.038)
    world.set_g("terrain.last_p121_compatible_edge_count", 44.0)
    world.set_g("terrain.last_p121_roughness_before_m", 180.0)
    world.set_g("terrain.last_p121_roughness_after_m", 120.0)
    world.set_g(
        "terrain.last_p121_non_orogenic_highland_500_area_fraction_before",
        0.070,
    )
    world.set_g(
        "terrain.last_p121_non_orogenic_highland_500_area_fraction_after",
        0.055,
    )
    world.set_g(
        "terrain.last_p121_non_orogenic_highland_1000_area_fraction_before",
        0.015,
    )
    world.set_g(
        "terrain.last_p121_non_orogenic_highland_1000_area_fraction_after",
        0.010,
    )
    world.set_g("terrain.last_p121_mean_abs_delta_m", 65.0)
    world.set_g("terrain.last_p121_max_abs_delta_m", 140.0)
    world.set_g("terrain.last_p121_p50_land_relief_before_m", 480.0)
    world.set_g("terrain.last_p121_p50_land_relief_after_m", 460.0)
    world.set_g("terrain.last_p121_p90_land_relief_before_m", 2200.0)
    world.set_g("terrain.last_p121_p90_land_relief_after_m", 2180.0)
    world.set_g("terrain.last_p111637a_terminal_spine_land_consistency_accepted", 1.0)
    world.set_g("terrain.last_p111637a_terminal_submerged_spine_cell_count_before", 5.0)
    world.set_g("terrain.last_p111637a_terminal_submerged_spine_cell_count_after", 1.0)
    world.set_g("terrain.last_p111637a_terminal_spine_lift_cell_count", 4.0)
    world.set_g("terrain.last_p111637a_terminal_spine_lift_area_fraction", 0.006)
    world.set_g("terrain.last_p111637b_final_spine_elevation_consistency_accepted", 1.0)
    world.set_g("terrain.last_p111637b_guard_reverted", 0.0)
    world.set_g("terrain.last_p111637b_candidate_cell_count", 12.0)
    world.set_g("terrain.last_p111637b_candidate_area_fraction", 0.018)
    world.set_g("terrain.last_p111637b_adjusted_cell_count", 9.0)
    world.set_g("terrain.last_p111637b_adjusted_area_fraction", 0.014)
    world.set_g("terrain.last_p111637b_submerged_spine_cell_count_before", 6.0)
    world.set_g("terrain.last_p111637b_submerged_spine_cell_count_after", 1.0)
    world.set_g("terrain.last_p111637b_bridge_blocked_submerged_spine_cell_count", 3.0)
    world.set_g("terrain.last_p111637b_bridge_blocked_submerged_spine_area_fraction", 0.005)
    world.set_g("terrain.last_p111637b_underexpressed_spine_cell_count_before", 5.0)
    world.set_g("terrain.last_p111637b_underexpressed_spine_cell_count_after", 2.0)
    world.set_g("terrain.last_p111637b_land_delta_area_fraction", 0.004)
    world.set_g("terrain.last_p111637b_mean_lift_m", 430.0)
    world.set_g("terrain.last_p111637b_max_lift_m", 780.0)
    world.set_g("terrain.last_p111637b_linework_score_before", 0.62)
    world.set_g("terrain.last_p111637b_linework_score_after", 0.68)

    metrics = terminal_audit_metrics(world)
    linework = metrics["p111635_orogenic_spine_linework"]

    assert linework["schema"] == "aevum.p111635_orogenic_spine_linework.v1"
    assert linework["peak_hierarchy_cell_count"] == 18
    assert linework["spine_cell_count"] == 8
    assert linework["spine_component_count"] >= 1
    assert linework["short_spine_component_count"] >= 0
    assert linework["spine_endpoint_count"] >= 0
    assert linework["spine_junction_count"] >= 0
    assert linework["spine_high_degree_count"] >= 0
    assert 0.0 <= linework["spine_high_degree_fraction"] <= 1.0
    assert 0.0 <= linework["linework_continuity_score"] <= 1.0
    assert 0.0 <= linework["branch_attachment_fraction"] <= 1.0
    assert linework["p117_highland_500_cell_count"] == 21
    assert linework["p117_highland_1000_cell_count"] == 21
    assert linework["p117_non_orogenic_highland_500_cell_count"] == 3
    assert linework["p117_non_orogenic_highland_1000_cell_count"] == 3
    assert 0.0 < linework["p117_highland_500_orogenic_support_fraction"] < 1.0
    assert 0.0 < linework["p117_highland_1000_orogenic_support_fraction"] < 1.0
    assert 0.0 < linework["p117_highland_500_mountain_inventory_fraction"] < 1.0
    assert linework["p111635_spine_linework_smoothing_accepted"] is True
    assert linework["p111635_spine_components_before_refine"] == 6.0
    assert linework["p111635_spine_components_after_refine"] == 3.0
    assert linework["p111635_short_spine_components_before_refine"] == 4.0
    assert linework["p111635_short_spine_components_after_refine"] == 1.0
    assert linework["p111635_branch_attachment_fraction_before"] == 0.2
    assert linework["p111635_branch_attachment_fraction_after"] == 0.8
    assert linework["p111635_bridge_cell_count"] == 5.0
    assert linework["p111635_bridge_candidate_cell_count"] == 11.0
    assert linework["p111635_path_count"] == 2.0
    assert linework["p111637a_spine_object_promotion_accepted"] is True
    assert linework["p111637a_spine_components_before"] == 5.0
    assert linework["p111637a_spine_components_after"] == 3.0
    assert linework["p111637a_short_spine_components_before"] == 2.0
    assert linework["p111637a_short_spine_components_after"] == 0.0
    assert linework["p111637a_branch_attachment_fraction_before"] == 0.4
    assert linework["p111637a_branch_attachment_fraction_after"] == 0.9
    assert linework["p111637a_spine_top3_share_before"] == 0.55
    assert linework["p111637a_spine_top3_share_after"] == 0.82
    assert linework["p111637a_linework_score_before"] == 0.50
    assert linework["p111637a_linework_score_after"] == 0.78
    assert linework["p111637a_bridge_cell_count"] == 4.0
    assert linework["p111637a_candidate_cell_count"] == 10.0
    assert linework["p111637a_path_count"] == 2.0
    assert linework["p132_parent_anchor_spine_promotion_accepted"] is True
    assert linework["p132_candidate_cell_count"] == 6.0
    assert linework["p132_promoted_spine_cell_count"] == 3.0
    assert linework["p132_parent_aligned_spine_fraction_before"] == 0.42
    assert linework["p132_parent_aligned_spine_fraction_after"] == 0.66
    assert linework["p132_linework_score_before"] == 0.78
    assert linework["p132_linework_score_after"] == 0.81
    assert linework["p132_spine_components_before"] == 3.0
    assert linework["p132_spine_components_after"] == 3.0
    assert linework["p132_short_spine_components_before"] == 1.0
    assert linework["p132_short_spine_components_after"] == 0.0
    assert linework["p1148_terminal_proxy_enabled"] is True
    assert linework["p1148_terminal_proxy_score_before"] == 0.86
    assert linework["p1148_terminal_proxy_score_after"] == 0.91
    assert linework["p1148_terminal_proxy_component_count_before"] == 4.0
    assert linework["p1148_terminal_proxy_component_count_after"] == 3.0
    assert linework["p1148_terminal_proxy_short_count_before"] == 1.0
    assert linework["p1148_terminal_proxy_short_count_after"] == 0.0
    assert linework["p1149_orogenic_spine_object_count"] == 3.0
    assert linework["p1149_orogenic_spine_system_count"] == 1.0
    assert linework["p1149_orogenic_spine_trunk_count"] == 1.0
    assert linework["p1149_orogenic_spine_branch_count"] == 1.0
    assert linework["p1149_fallback_trunk_count"] == 0.0
    assert linework["p1149_orphan_branch_count"] == 0.0
    assert linework["p1149_attached_branch_fraction"] == 1.0
    assert linework["p1149_mean_branch_count_per_system"] == 1.0
    assert linework["p1149_spine_cell_count"] == 8.0
    assert linework["p1149_trunk_cell_count"] == 4.0
    assert linework["p1149_branch_cell_count"] == 4.0
    assert linework["p1149_endpoint_count"] == 2.0
    assert metrics["object_count_by_set"]["terrain.orogenic_spine_objects"] == 1
    assert linework["p128_orogenic_axis_polyline_object_count"] == 2.0
    assert linework["p128_main_axis_count"] == 1.0
    assert linework["p128_branch_axis_count"] == 1.0
    assert linework["p128_fallback_main_axis_count"] == 0.0
    assert linework["p128_attached_branch_axis_fraction"] == 1.0
    assert linework["p128_source_spine_cell_count"] == 8.0
    assert linework["p128_axis_cell_count"] == 7.0
    assert linework["p128_axis_source_coverage_fraction"] == 0.875
    assert linework["p128_mean_path_coverage_fraction"] == 0.82
    assert linework["p128_mean_directness"] == 0.74
    assert linework["p128_mean_sinuosity"] == 1.35
    assert linework["p128_max_sinuosity"] == 1.62
    assert linework["p128_source_junction_cell_count"] == 3.0
    assert linework["p128_source_high_degree_cell_count"] == 1.0
    assert linework["p128_polyline_ready"] is True
    assert metrics["object_count_by_set"]["terrain.orogenic_axis_polylines"] == 2
    assert linework["p115_spine_repair_candidate_count"] == 2.0
    assert linework["p115_viable_spine_repair_candidate_count"] == 1.0
    assert linework["p115_trunk_bridge_attempt_count"] == 3.0
    assert linework["p115_branch_attachment_attempt_count"] == 2.0
    assert linework["p115_trunk_bridge_candidate_count"] == 1.0
    assert linework["p115_branch_attachment_candidate_count"] == 1.0
    assert linework["p115_rejected_proxy_count"] == 1.0
    assert linework["p115_multi_trunk_system_count"] == 1.0
    assert linework["p115_detached_branch_component_count"] == 1.0
    assert linework["p115_best_proxy_score_delta"] == 0.05
    assert linework["p115_best_component_delta"] == 1.0
    assert linework["p115_candidate_cell_count"] == 3.0
    assert linework["p115_viable_candidate_cell_count"] == 2.0
    assert metrics["object_count_by_set"][
        "terrain.orogenic_spine_repair_candidates"] == 1
    assert linework["p116_candidate_promotion_enabled"] is True
    assert linework["p116_candidate_promotion_accepted"] is True
    assert linework["p116_guard_reverted"] is False
    assert linework["p116_input_candidate_count"] == 2.0
    assert linework["p116_input_viable_candidate_count"] == 1.0
    assert linework["p116_considered_candidate_count"] == 1.0
    assert linework["p116_selected_candidate_count"] == 1.0
    assert linework["p116_applied_candidate_count"] == 1.0
    assert linework["p116_rejected_candidate_count"] == 0.0
    assert linework["p116_rejected_overlap_count"] == 0.0
    assert linework["p116_rejected_support_count"] == 0.0
    assert linework["p116_rejected_profile_count"] == 0.0
    assert linework["p116_applied_cell_count"] == 2.0
    assert linework["p116_applied_area_fraction"] == 0.004
    assert linework["p116_area_budget_fraction"] == 0.004
    assert linework["p116_linework_score_before"] == 0.61
    assert linework["p116_linework_score_after"] == 0.70
    assert linework["p116_component_count_before"] == 5.0
    assert linework["p116_component_count_after"] == 4.0
    assert linework["p116_short_count_before"] == 1.0
    assert linework["p116_short_count_after"] == 0.0
    assert linework["p116_branch_attachment_fraction_before"] == 0.8
    assert linework["p116_branch_attachment_fraction_after"] == 1.0
    assert linework["p116_spine_top3_share_before"] == 0.7
    assert linework["p116_spine_top3_share_after"] == 0.82
    assert linework[
        "p124_orogenic_spine_geometry_regularization_accepted"] is True
    assert linework["p124_guard_reverted"] is False
    assert linework["p124_candidate_cell_count"] == 7.0
    assert linework["p124_candidate_area_fraction"] == 0.018
    assert linework["p124_added_cell_count"] == 3.0
    assert linework["p124_added_area_fraction"] == 0.007
    assert linework["p124_component_count_before"] == 5.0
    assert linework["p124_component_count_after"] == 4.0
    assert linework["p124_endpoint_count_before"] == 10.0
    assert linework["p124_endpoint_count_after"] == 8.0
    assert linework["p124_short_component_count_before"] == 2.0
    assert linework["p124_short_component_count_after"] == 1.0
    assert linework["p124_branch_attachment_fraction_before"] == 0.7
    assert linework["p124_branch_attachment_fraction_after"] == 0.85
    assert linework["p124_linework_score_before"] == 0.51
    assert linework["p124_linework_score_after"] == 0.62
    assert (
        linework[
            "p127_terminal_orogenic_spine_node_thinning_accepted"
        ] is True
    )
    assert linework["p127_guard_reverted"] is False
    assert linework["p127_candidate_cell_count"] == 12.0
    assert linework["p127_candidate_area_fraction"] == 0.026
    assert linework["p127_demoted_cell_count"] == 4.0
    assert linework["p127_demoted_area_fraction"] == 0.009
    assert linework["p127_spine_cell_count_before"] == 31.0
    assert linework["p127_spine_cell_count_after"] == 27.0
    assert linework["p127_component_count_before"] == 4.0
    assert linework["p127_component_count_after"] == 4.0
    assert linework["p127_endpoint_count_before"] == 8.0
    assert linework["p127_endpoint_count_after"] == 8.0
    assert linework["p127_junction_count_before"] == 15.0
    assert linework["p127_junction_count_after"] == 10.0
    assert linework["p127_high_degree_count_before"] == 7.0
    assert linework["p127_high_degree_count_after"] == 2.0
    assert linework["p127_short_component_count_before"] == 1.0
    assert linework["p127_short_component_count_after"] == 1.0
    assert linework["p127_branch_attachment_fraction_before"] == 0.85
    assert linework["p127_branch_attachment_fraction_after"] == 0.85
    assert linework["p127_linework_score_before"] == 0.58
    assert linework["p127_linework_score_after"] == 0.64
    assert linework["p127_junction_pressure_before"] == 0.23
    assert linework["p127_junction_pressure_after"] == 0.07
    assert (
        linework[
            "p140_terminal_branch_spine_component_thickening_accepted"
        ] is True
    )
    assert linework["p140_guard_reverted"] is False
    assert linework["p140_candidate_cell_count"] == 6.0
    assert linework["p140_promoted_branch_spine_cell_count"] == 6.0
    assert linework["p140_grown_branch_spine_component_count"] == 1.0
    assert linework["p140_branch_spine_cell_count_before"] == 30.0
    assert linework["p140_branch_spine_cell_count_after"] == 36.0
    assert linework["p140_branch_spine_component_count_before"] == 6.0
    assert linework["p140_branch_spine_component_count_after"] == 6.0
    assert linework["p140_branch_spine_small_area_fraction_before"] == 0.333
    assert linework["p140_branch_spine_small_area_fraction_after"] == 0.178
    assert linework["p140_all_spine_component_count_before"] == 2.0
    assert linework["p140_all_spine_component_count_after"] == 2.0
    assert linework["p140_junction_count_before"] == 32.0
    assert linework["p140_junction_count_after"] == 38.0
    assert linework["p140_high_degree_count_before"] == 7.0
    assert linework["p140_high_degree_count_after"] == 12.0
    assert linework["p140_junction_fraction_before"] == 0.421
    assert linework["p140_junction_fraction_after"] == 0.463
    assert linework["p140_high_degree_fraction_before"] == 0.092
    assert linework["p140_high_degree_fraction_after"] == 0.146
    assert linework["p140_branch_attachment_fraction_before"] == 1.0
    assert linework["p140_branch_attachment_fraction_after"] == 1.0
    assert linework["p140_hierarchy_changed_cell_count"] == 0.0
    assert linework["p140_crest_spine_changed_cell_count"] == 0.0
    assert (
        linework[
            "p141_terminal_crest_spine_component_thickening_accepted"
        ] is True
    )
    assert linework["p141_guard_reverted"] is False
    assert linework["p141_candidate_cell_count"] == 3.0
    assert linework["p141_promoted_crest_spine_cell_count"] == 3.0
    assert linework["p141_connected_crest_spine_component_count"] == 1.0
    assert linework["p141_grown_crest_spine_component_count"] == 0.0
    assert linework["p141_crest_spine_cell_count_before"] == 49.0
    assert linework["p141_crest_spine_cell_count_after"] == 52.0
    assert linework["p141_crest_spine_component_count_before"] == 10.0
    assert linework["p141_crest_spine_component_count_after"] == 9.0
    assert linework["p141_crest_spine_small_area_fraction_before"] == 0.711
    assert linework["p141_crest_spine_small_area_fraction_after"] == 0.558
    assert linework["p141_all_spine_component_count_before"] == 5.0
    assert linework["p141_all_spine_component_count_after"] == 4.0
    assert linework["p141_junction_count_before"] == 38.0
    assert linework["p141_junction_count_after"] == 39.0
    assert linework["p141_high_degree_count_before"] == 12.0
    assert linework["p141_high_degree_count_after"] == 13.0
    assert linework["p141_junction_fraction_before"] == 0.463
    assert linework["p141_junction_fraction_after"] == 0.459
    assert linework["p141_high_degree_fraction_before"] == 0.146
    assert linework["p141_high_degree_fraction_after"] == 0.153
    assert linework["p141_branch_attachment_fraction_before"] == 1.0
    assert linework["p141_branch_attachment_fraction_after"] == 1.0
    assert linework["p141_hierarchy_changed_cell_count"] == 0.0
    assert linework["p141_branch_spine_changed_cell_count"] == 0.0
    assert (
        linework[
            "p164_terminal_orogenic_axis_skeleton_simplification_accepted"
        ]
        is True
    )
    assert linework["p164_guard_reverted"] is False
    assert linework["p164_candidate_cell_count"] == 9.0
    assert linework["p164_demoted_cell_count"] == 4.0
    assert linework["p164_off_axis_candidate_cell_count"] == 12.0
    assert linework["p164_spine_cell_count_before"] == 36.0
    assert linework["p164_spine_cell_count_after"] == 32.0
    assert linework["p164_component_count_before"] == 3.0
    assert linework["p164_component_count_after"] == 3.0
    assert linework["p164_endpoint_count_before"] == 8.0
    assert linework["p164_endpoint_count_after"] == 8.0
    assert linework["p164_junction_count_before"] == 14.0
    assert linework["p164_junction_count_after"] == 9.0
    assert linework["p164_high_degree_count_before"] == 4.0
    assert linework["p164_high_degree_count_after"] == 2.0
    assert linework["p164_axis_score_before"] == 0.61
    assert linework["p164_axis_score_after"] == 0.68
    assert linework["p164_protected_extreme_demoted_cell_count"] == 0.0
    assert linework["p142_class_aware_combined_score"] == 0.96
    assert linework["p142_class_aware_class_score"] == 0.58
    assert linework["p142_class_aware_crest_small_area_fraction"] == 0.71
    assert linework["p142_class_aware_branch_small_area_fraction"] == 0.65
    assert linework["p142_class_aware_class_small_area_fraction"] == 0.71
    assert linework["p142_class_aware_crest_component_count"] == 10.0
    assert linework["p142_class_aware_branch_component_count"] == 9.0
    assert linework["p142_class_aware_blind_spot"] is True
    assert linework["p143_class_repair_needed"] is True
    assert linework["p143_class_profile_improved"] is True
    assert linework["p143_class_score_before"] == 0.58
    assert linework["p143_class_score_after"] == 0.66
    assert linework["p143_class_small_area_fraction_before"] == 0.71
    assert linework["p143_class_small_area_fraction_after"] == 0.55
    assert linework["p143_crest_small_area_fraction_before"] == 0.71
    assert linework["p143_crest_small_area_fraction_after"] == 0.55
    assert linework["p143_branch_small_area_fraction_before"] == 0.65
    assert linework["p143_branch_small_area_fraction_after"] == 0.50
    assert linework["p144_class_path_option_count"] == 7.0
    assert linework["p144_class_path_selected_count"] == 2.0
    assert linework["p144_crest_path_option_count"] == 4.0
    assert linework["p144_branch_path_option_count"] == 3.0
    assert linework["p144_class_attempted_path_count"] == 9.0
    assert linework["p144_class_found_path_count"] == 8.0
    assert linework["p144_crest_promoted_spine_cell_count"] == 5.0
    assert linework["p144_branch_promoted_spine_cell_count"] == 4.0
    assert linework["p144_class_promoted_spine_cell_count"] == 9.0
    assert linework["p144_class_path_profile_rejected_count"] == 2.0
    assert linework[
        "p145_class_hierarchy_component_consolidation_accepted"] is True
    assert linework["p145_removed_crest_component_count"] == 2.0
    assert linework["p145_removed_branch_component_count"] == 3.0
    assert linework["p145_removed_crest_cell_count"] == 5.0
    assert linework["p145_removed_branch_cell_count"] == 7.0
    assert linework["p145_removed_crest_area_fraction"] == 0.011
    assert linework["p145_removed_branch_area_fraction"] == 0.014
    assert linework["p145_added_crest_cell_count"] == 6.0
    assert linework["p145_added_branch_cell_count"] == 8.0
    assert linework["p145_added_crest_area_fraction"] == 0.016
    assert linework["p145_added_branch_area_fraction"] == 0.019
    assert linework["p145_bridge_path_count"] == 3.0
    assert linework["p145_class_score_before"] == 0.45
    assert linework["p145_class_score_after"] == 0.62
    assert linework["p145_class_small_area_fraction_before"] == 1.0
    assert linework["p145_class_small_area_fraction_after"] == 0.69
    assert linework["p145_crest_small_area_fraction_before"] == 0.80
    assert linework["p145_crest_small_area_fraction_after"] == 0.52
    assert linework["p145_branch_small_area_fraction_before"] == 1.0
    assert linework["p145_branch_small_area_fraction_after"] == 0.69
    assert linework[
        "p118_non_orogenic_midland_stripe_suppression_accepted"] is True
    assert linework["p118_guard_reverted"] is False
    assert linework["p118_land_mask_preserved"] is True
    assert linework["p118_candidate_cell_count"] == 4.0
    assert linework["p118_candidate_area_fraction"] == 0.012
    assert linework["p118_adjusted_cell_count"] == 3.0
    assert linework["p118_adjusted_area_fraction"] == 0.009
    assert linework[
        "p118_unsupported_midland_area_fraction_before"] == 0.018
    assert linework[
        "p118_unsupported_midland_area_fraction_after"] == 0.006
    assert linework[
        "p118_non_orogenic_highland_500_area_fraction_before"] == 0.024
    assert linework[
        "p118_non_orogenic_highland_500_area_fraction_after"] == 0.014
    assert linework["p118_mean_lowering_m"] == 220.0
    assert linework["p118_max_lowering_m"] == 410.0
    assert linework["p118_mean_land_relief_before_m"] == 720.0
    assert linework["p118_mean_land_relief_after_m"] == 690.0
    assert linework["p118_p50_land_relief_before_m"] == 430.0
    assert linework["p118_p50_land_relief_after_m"] == 390.0
    assert linework["p118_p90_land_relief_before_m"] == 1850.0
    assert linework["p118_p90_land_relief_after_m"] == 1810.0
    assert linework[
        "p119_non_orogenic_interior_anti_raster_smoothing_accepted"] is True
    assert linework["p119_guard_reverted"] is False
    assert linework["p119_land_mask_preserved"] is True
    assert linework["p119_candidate_cell_count"] == 12.0
    assert linework["p119_candidate_area_fraction"] == 0.030
    assert linework["p119_selected_cell_count"] == 8.0
    assert linework["p119_selected_area_fraction"] == 0.020
    assert linework["p119_adjusted_cell_count"] == 7.0
    assert linework["p119_adjusted_area_fraction"] == 0.018
    assert linework["p119_compatible_edge_count"] == 19.0
    assert linework["p119_roughness_before_m"] == 310.0
    assert linework["p119_roughness_after_m"] == 240.0
    assert linework[
        "p119_non_orogenic_highland_500_area_fraction_before"] == 0.064
    assert linework[
        "p119_non_orogenic_highland_500_area_fraction_after"] == 0.060
    assert linework["p119_mean_abs_delta_m"] == 54.0
    assert linework["p119_max_abs_delta_m"] == 105.0
    assert linework["p119_p50_land_relief_before_m"] == 390.0
    assert linework["p119_p50_land_relief_after_m"] == 386.0
    assert linework["p119_p90_land_relief_before_m"] == 1810.0
    assert linework["p119_p90_land_relief_after_m"] == 1810.0
    assert linework[
        "p120_continental_semantic_geometry_repair_accepted"] is True
    assert linework["p120_guard_reverted"] is False
    assert linework["p120_candidate_cell_count"] == 16.0
    assert linework["p120_candidate_area_fraction"] == 0.042
    assert linework["p120_detail_changed_area_fraction"] == 0.012
    assert linework["p120_internal_changed_area_fraction"] == 0.014
    assert linework["p120_terrain_changed_area_fraction"] == 0.004
    assert linework["p120_detail_elongated_before"] == 0.18
    assert linework["p120_detail_elongated_after"] == 0.11
    assert linework["p120_internal_elongated_before"] == 0.24
    assert linework["p120_internal_elongated_after"] == 0.15
    assert linework["p120_terrain_elongated_before"] == 0.10
    assert linework["p120_terrain_elongated_after"] == 0.08
    assert linework["p120_detail_tiny_before"] == 0.020
    assert linework["p120_detail_tiny_after"] == 0.018
    assert linework["p120_internal_tiny_before"] == 0.030
    assert linework["p120_internal_tiny_after"] == 0.026
    assert linework["p120_terrain_tiny_before"] == 0.012
    assert linework["p120_terrain_tiny_after"] == 0.012
    assert linework[
        "p122_orogenic_detail_footprint_compaction_accepted"] is True
    assert linework["p122_guard_reverted"] is False
    assert linework["p122_candidate_cell_count"] == 9.0
    assert linework["p122_candidate_area_fraction"] == 0.021
    assert linework["p122_detail_changed_area_fraction"] == 0.017
    assert linework["p122_internal_changed_area_fraction"] == 0.015
    assert linework["p122_terrain_changed_area_fraction"] == 0.012
    assert linework["p122_detail_orogen_area_fraction_before"] == 0.110
    assert linework["p122_detail_orogen_area_fraction_after"] == 0.093
    assert linework[
        "p122_low_relief_orogen_area_fraction_before"] == 0.078
    assert linework[
        "p122_low_relief_orogen_area_fraction_after"] == 0.055
    assert linework["p122_noncore_orogen_area_fraction_before"] == 0.045
    assert linework["p122_noncore_orogen_area_fraction_after"] == 0.022
    assert linework["p122_preserved_core_orogen_cell_count"] == 11.0
    assert linework[
        "p121_semantic_region_elevation_response_accepted"] is True
    assert linework["p121_guard_reverted"] is False
    assert linework["p121_land_mask_preserved"] is True
    assert linework["p121_candidate_cell_count"] == 22.0
    assert linework["p121_candidate_area_fraction"] == 0.052
    assert linework["p121_selected_cell_count"] == 17.0
    assert linework["p121_selected_area_fraction"] == 0.041
    assert linework["p121_adjusted_cell_count"] == 15.0
    assert linework["p121_adjusted_area_fraction"] == 0.038
    assert linework["p121_compatible_edge_count"] == 44.0
    assert linework["p121_roughness_before_m"] == 180.0
    assert linework["p121_roughness_after_m"] == 120.0
    assert linework[
        "p121_non_orogenic_highland_500_area_fraction_before"] == 0.070
    assert linework[
        "p121_non_orogenic_highland_500_area_fraction_after"] == 0.055
    assert linework[
        "p121_non_orogenic_highland_1000_area_fraction_before"] == 0.015
    assert linework[
        "p121_non_orogenic_highland_1000_area_fraction_after"] == 0.010
    assert linework["p121_mean_abs_delta_m"] == 65.0
    assert linework["p121_max_abs_delta_m"] == 140.0
    assert linework["p121_p50_land_relief_before_m"] == 480.0
    assert linework["p121_p50_land_relief_after_m"] == 460.0
    assert linework["p121_p90_land_relief_before_m"] == 2200.0
    assert linework["p121_p90_land_relief_after_m"] == 2180.0
    assert linework["p111637a_terminal_spine_land_consistency_accepted"] is True
    assert linework["p111637a_terminal_submerged_spine_cell_count_before"] == 5.0
    assert linework["p111637a_terminal_submerged_spine_cell_count_after"] == 1.0
    assert linework["p111637a_terminal_spine_lift_cell_count"] == 4.0
    assert linework["p111637a_terminal_spine_lift_area_fraction"] == 0.006
    assert linework["p111637b_final_spine_elevation_consistency_accepted"] is True
    assert linework["p111637b_guard_reverted"] is False
    assert linework["p111637b_candidate_cell_count"] == 12.0
    assert linework["p111637b_candidate_area_fraction"] == 0.018
    assert linework["p111637b_adjusted_cell_count"] == 9.0
    assert linework["p111637b_adjusted_area_fraction"] == 0.014
    assert linework["p111637b_submerged_spine_cell_count_before"] == 6.0
    assert linework["p111637b_submerged_spine_cell_count_after"] == 1.0
    assert linework["p111637b_bridge_blocked_submerged_spine_cell_count"] == 3.0
    assert linework["p111637b_bridge_blocked_submerged_spine_area_fraction"] == 0.005
    assert linework["p111637b_underexpressed_spine_cell_count_before"] == 5.0
    assert linework["p111637b_underexpressed_spine_cell_count_after"] == 2.0
    assert linework["p111637b_land_delta_area_fraction"] == 0.004
    assert linework["p111637b_mean_lift_m"] == 430.0
    assert linework["p111637b_max_lift_m"] == 780.0
    assert linework["p111637b_linework_score_before"] == 0.62
    assert linework["p111637b_linework_score_after"] == 0.68
    assert metrics["acceptance"]["has_p111635_orogenic_spine_linework_metrics"]


def test_p111636_audit_reports_polar_edge_overclassification_metrics():
    spec = get_preset("earthlike")
    spec.grid_cells = 420
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    polar_peak = np.where((grid.lat > 62.0) & (grid.lon > -80.0) & (grid.lon < 80.0))[0]
    edge_peak = np.where((np.abs(grid.lat) < 18.0) & (grid.lon > 172.0))[0]
    assert polar_peak.size >= 3
    assert edge_peak.size >= 2
    hierarchy[polar_peak[:6]] = 2.0
    hierarchy[edge_peak[:4]] = 2.0
    spine[polar_peak[:2]] = 2.0
    elev = np.full(grid.n, 500.0, dtype=np.float64)
    elev[hierarchy > 0.0] = 1800.0
    world.fields.update({
        "terrain.elevation_m": elev,
        "crust.type": np.full(grid.n, CONT, dtype=np.float64),
        "terrain.orogenic_parent_hierarchy": hierarchy,
        "terrain.orogenic_hierarchy_spine": spine,
        "terrain.orogenic_shoulder_halo": np.zeros(grid.n, dtype=np.float64),
        "terrain.orogenic_highland_apron": np.zeros(grid.n, dtype=np.float64),
    })
    world.set_g("terrain.last_p111636_polar_edge_refinement_accepted", 1.0)
    world.set_g("terrain.last_p111636_polar_peak_area_fraction_before", 0.02)
    world.set_g("terrain.last_p111636_polar_peak_area_fraction_after", 0.012)
    world.set_g("terrain.last_p111636_edge_peak_component_count_before", 5.0)
    world.set_g("terrain.last_p111636_edge_peak_component_count_after", 3.0)
    world.set_g("terrain.last_p111636_candidate_cell_count", 10.0)
    world.set_g("terrain.last_p111636_candidate_area_fraction", 0.016)
    world.set_g("terrain.last_p111636_reclassified_cell_count", 8.0)
    world.set_g("terrain.last_p111636_reclassified_area_fraction", 0.011)
    world.set_g("terrain.last_p111636_polar_reclassified_cell_count", 5.0)
    world.set_g("terrain.last_p111636_edge_reclassified_cell_count", 3.0)
    world.set_g("terrain.last_p111636_extreme_reclassified_cell_count", 0.0)
    world.set_g("terrain.last_p160_polar_edge_spine_thinning_accepted", 1.0)
    world.set_g("terrain.last_p160_guard_reverted", 0.0)
    world.set_g("terrain.last_p160_candidate_cell_count", 7.0)
    world.set_g("terrain.last_p160_candidate_area_fraction", 0.009)
    world.set_g("terrain.last_p160_demoted_cell_count", 6.0)
    world.set_g("terrain.last_p160_demoted_area_fraction", 0.007)
    world.set_g("terrain.last_p160_polar_demoted_cell_count", 4.0)
    world.set_g("terrain.last_p160_edge_demoted_cell_count", 2.0)
    world.set_g("terrain.last_p160_extreme_demoted_cell_count", 0.0)
    world.set_g("terrain.last_p160_polar_spine_area_fraction_before", 0.011)
    world.set_g("terrain.last_p160_polar_spine_area_fraction_after", 0.008)
    world.set_g("terrain.last_p160_edge_spine_component_count_before", 3.0)
    world.set_g("terrain.last_p160_edge_spine_component_count_after", 2.0)
    world.set_g("terrain.last_p160_spine_component_count_before", 9.0)
    world.set_g("terrain.last_p160_spine_component_count_after", 8.0)
    world.set_g("terrain.last_p160_short_spine_component_count_before", 2.0)
    world.set_g("terrain.last_p160_short_spine_component_count_after", 1.0)
    world.set_g("terrain.last_p160_branch_attachment_fraction_before", 0.82)
    world.set_g("terrain.last_p160_branch_attachment_fraction_after", 0.84)
    world.set_g("terrain.last_p160_linework_score_before", 0.57)
    world.set_g("terrain.last_p160_linework_score_after", 0.59)
    world.set_g(
        "terrain.last_p165_terminal_polar_edge_orogenic_semantic_compaction_accepted",
        1.0,
    )
    world.set_g("terrain.last_p165_guard_reverted", 0.0)
    world.set_g("terrain.last_p165_candidate_cell_count", 11.0)
    world.set_g("terrain.last_p165_candidate_area_fraction", 0.013)
    world.set_g("terrain.last_p165_peak_reclassified_cell_count", 6.0)
    world.set_g("terrain.last_p165_peak_reclassified_area_fraction", 0.007)
    world.set_g("terrain.last_p165_fringe_cleared_cell_count", 5.0)
    world.set_g("terrain.last_p165_fringe_cleared_area_fraction", 0.006)
    world.set_g("terrain.last_p165_polar_semantic_area_fraction_before", 0.031)
    world.set_g("terrain.last_p165_polar_semantic_area_fraction_after", 0.025)
    world.set_g("terrain.last_p165_edge_semantic_area_fraction_before", 0.014)
    world.set_g("terrain.last_p165_edge_semantic_area_fraction_after", 0.012)
    world.set_g("terrain.last_p165_polar_peak_area_fraction_before", 0.019)
    world.set_g("terrain.last_p165_polar_peak_area_fraction_after", 0.015)
    world.set_g("terrain.last_p165_edge_peak_component_count_before", 6.0)
    world.set_g("terrain.last_p165_edge_peak_component_count_after", 5.0)
    world.set_g("terrain.last_p165_peak_component_count_before", 13.0)
    world.set_g("terrain.last_p165_peak_component_count_after", 12.0)
    world.set_g("terrain.last_p165_hierarchy_changed_cell_count", 6.0)
    world.set_g("terrain.last_p165_halo_changed_cell_count", 8.0)
    world.set_g("terrain.last_p165_apron_changed_cell_count", 3.0)
    world.set_g("terrain.last_p165_protected_spine_changed_cell_count", 0.0)
    world.set_g("terrain.last_p165_protected_extreme_reclassified_cell_count", 0.0)
    world.set_g("terrain.last_p165_reject_code", 0.0)
    world.set_g(
        "terrain.last_p125_terminal_orogenic_semantic_land_consistency_accepted",
        1.0,
    )
    world.set_g("terrain.last_p125_guard_reverted", 0.0)
    world.set_g("terrain.last_p125_candidate_cell_count", 9.0)
    world.set_g("terrain.last_p125_candidate_area_fraction", 0.014)
    world.set_g("terrain.last_p125_cleared_hierarchy_cell_count", 4.0)
    world.set_g("terrain.last_p125_cleared_hierarchy_area_fraction", 0.006)
    world.set_g("terrain.last_p125_cleared_spine_cell_count", 2.0)
    world.set_g("terrain.last_p125_cleared_spine_area_fraction", 0.003)
    world.set_g("terrain.last_p125_cleared_halo_cell_count", 3.0)
    world.set_g("terrain.last_p125_cleared_halo_area_fraction", 0.005)
    world.set_g("terrain.last_p125_cleared_apron_cell_count", 1.0)
    world.set_g("terrain.last_p125_cleared_apron_area_fraction", 0.002)
    world.set_g("terrain.last_p125_submerged_hierarchy_cell_count_before", 4.0)
    world.set_g("terrain.last_p125_submerged_hierarchy_cell_count_after", 0.0)
    world.set_g("terrain.last_p125_submerged_spine_cell_count_before", 2.0)
    world.set_g("terrain.last_p125_submerged_spine_cell_count_after", 0.0)
    world.set_g("terrain.last_p125_submerged_halo_cell_count_before", 3.0)
    world.set_g("terrain.last_p125_submerged_halo_cell_count_after", 0.0)
    world.set_g("terrain.last_p125_submerged_apron_cell_count_before", 1.0)
    world.set_g("terrain.last_p125_submerged_apron_cell_count_after", 0.0)
    world.set_g("terrain.last_p125_polar_cleared_cell_count", 2.0)
    world.set_g("terrain.last_p125_edge_cleared_cell_count", 1.0)
    world.set_g("terrain.last_p125_land_semantic_cleared_cell_count", 0.0)

    metrics = terminal_audit_metrics(world)
    polar_edge = metrics["p111636_polar_edge_orogen_overclassification"]

    assert polar_edge["schema"] == (
        "aevum.p111636_polar_edge_orogen_overclassification.v3"
    )
    assert polar_edge["polar_peak_cell_count"] > 0
    assert polar_edge["polar_spine_cell_count"] > 0
    assert polar_edge["edge_peak_component_count"] >= 1
    assert polar_edge["p111636_polar_edge_refinement_accepted"] is True
    assert polar_edge["p111636_polar_peak_area_fraction_delta"] == 0.008
    assert polar_edge["p111636_edge_peak_component_count_delta"] == 2.0
    assert polar_edge["p111636_candidate_cell_count"] == 10.0
    assert polar_edge["p111636_reclassified_cell_count"] == 8.0
    assert polar_edge["p111636_polar_reclassified_cell_count"] == 5.0
    assert polar_edge["p111636_edge_reclassified_cell_count"] == 3.0
    assert polar_edge["p111636_no_extreme_reclassification"] is True
    assert polar_edge["p111636_nonexpansive_polar_edge"] is True
    assert polar_edge["p160_polar_edge_spine_thinning_accepted"] is True
    assert polar_edge["p160_guard_reverted"] is False
    assert polar_edge["p160_candidate_cell_count"] == 7.0
    assert polar_edge["p160_demoted_cell_count"] == 6.0
    assert polar_edge["p160_polar_demoted_cell_count"] == 4.0
    assert polar_edge["p160_edge_demoted_cell_count"] == 2.0
    assert polar_edge["p160_extreme_demoted_cell_count"] == 0.0
    assert math.isclose(
        polar_edge["p160_polar_spine_area_fraction_delta"],
        0.003,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert polar_edge["p160_edge_spine_component_count_delta"] == 1.0
    assert polar_edge["p160_spine_component_count_after"] == 8.0
    assert polar_edge["p160_short_spine_component_count_after"] == 1.0
    assert polar_edge["p160_branch_attachment_fraction_after"] == 0.84
    assert polar_edge["p160_linework_score_after"] == 0.59
    assert polar_edge["p160_no_extreme_demotion"] is True
    assert polar_edge["p160_nonexpansive_polar_edge_spine"] is True
    assert (
        polar_edge[
            "p165_terminal_polar_edge_orogenic_semantic_compaction_accepted"
        ] is True
    )
    assert polar_edge["p165_guard_reverted"] is False
    assert polar_edge["p165_candidate_cell_count"] == 11.0
    assert polar_edge["p165_peak_reclassified_cell_count"] == 6.0
    assert polar_edge["p165_fringe_cleared_cell_count"] == 5.0
    assert math.isclose(
        polar_edge["p165_polar_semantic_area_fraction_delta"],
        0.006,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert math.isclose(
        polar_edge["p165_edge_semantic_area_fraction_delta"],
        0.002,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert math.isclose(
        polar_edge["p165_polar_peak_area_fraction_delta"],
        0.004,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert polar_edge["p165_edge_peak_component_count_delta"] == 1.0
    assert polar_edge["p165_peak_component_count_after"] == 12.0
    assert polar_edge["p165_hierarchy_changed_cell_count"] == 6.0
    assert polar_edge["p165_halo_changed_cell_count"] == 8.0
    assert polar_edge["p165_apron_changed_cell_count"] == 3.0
    assert polar_edge["p165_no_protected_spine_change"] is True
    assert polar_edge["p165_no_extreme_reclassification"] is True
    assert polar_edge["p165_nonexpansive_polar_edge_semantics"] is True
    assert polar_edge["submerged_orogenic_semantic_cell_count"] == 0
    assert (
        polar_edge[
            "p125_terminal_orogenic_semantic_land_consistency_accepted"
        ] is True
    )
    assert polar_edge["p125_guard_reverted"] is False
    assert polar_edge["p125_candidate_cell_count"] == 9.0
    assert polar_edge["p125_cleared_hierarchy_cell_count"] == 4.0
    assert polar_edge["p125_cleared_spine_cell_count"] == 2.0
    assert polar_edge["p125_cleared_halo_cell_count"] == 3.0
    assert polar_edge["p125_cleared_apron_cell_count"] == 1.0
    assert polar_edge["p125_submerged_hierarchy_cell_count_before"] == 4.0
    assert polar_edge["p125_submerged_hierarchy_cell_count_after"] == 0.0
    assert polar_edge["p125_submerged_spine_cell_count_before"] == 2.0
    assert polar_edge["p125_submerged_spine_cell_count_after"] == 0.0
    assert polar_edge["p125_polar_cleared_cell_count"] == 2.0
    assert polar_edge["p125_edge_cleared_cell_count"] == 1.0
    assert polar_edge["p125_no_submerged_terminal_orogenic_semantics"] is True
    assert polar_edge["p125_land_semantics_preserved"] is True
    assert metrics["acceptance"][
        "has_p111636_polar_edge_orogen_overclassification_metrics"]


def test_p107_ranked_policy_protects_process_parented_microplates():
    grid = SphereGrid.fibonacci(900, CONSTANTS.EARTH_RADIUS)
    module = TectonicsModule()
    seed = int(np.argmax(grid.xyz @ np.array([1.0, 0.0, 0.0])))
    micro_cells = [seed]
    seen = {seed}
    idx = 0
    while idx < len(micro_cells) and len(micro_cells) < 5:
        c = micro_cells[idx]
        idx += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in seen:
                continue
            seen.add(nb)
            micro_cells.append(nb)
            if len(micro_cells) >= 5:
                break
    micro = np.asarray(micro_cells, dtype=np.int64)
    plate = np.zeros(grid.n, dtype=np.int64)
    plate[micro] = 1
    plates = [
        {"id": 0, "pole": [0.0, 0.0, 1.0], "rate": 0.008},
        {"id": 1, "pole": [0.0, 1.0, 0.0], "rate": 0.004},
    ]
    origin = np.full(grid.n, ORIGIN_RIDGE, dtype=np.float64)
    terrane = np.full(grid.n, -1.0, dtype=np.float64)

    ocean_ctype = np.full(grid.n, OCEAN, dtype=np.float64)
    no_continent = np.full(grid.n, -1.0, dtype=np.float64)
    ocean_topology = module._plate_topology_objects(
        grid, plate, plates, ocean_ctype, no_continent, terrane, {}, 100.0)
    ocean_plate = plate.copy()
    ocean_context = {
        "ctype": ocean_ctype,
        "origin": origin,
        "age": np.full(grid.n, 80.0),
        "plate_topologies": ocean_topology,
        "enable_p107_ranked_plate_policy": True,
    }
    ocean_merges = module._merge_tiny_plates(grid, ocean_plate, ocean_context)

    cont_ctype = ocean_ctype.copy()
    cont_ctype[micro] = CONT
    continent = no_continent.copy()
    continent[micro] = 2.0
    cont_topology = module._plate_topology_objects(
        grid, plate, plates, cont_ctype, continent, terrane, {}, 100.0)
    protected_plate = plate.copy()
    protected_context = {
        "ctype": cont_ctype,
        "origin": origin,
        "age": np.full(grid.n, 800.0),
        "plate_topologies": cont_topology,
        "enable_p107_ranked_plate_policy": True,
    }
    protected_merges = module._merge_tiny_plates(
        grid, protected_plate, protected_context)

    assert len(ocean_merges) == 1
    assert np.all(ocean_plate[micro] == 0)
    assert protected_merges == []
    assert np.all(protected_plate[micro] == 1)
    protected = protected_context["protected_microplates"]
    assert protected[0]["plate_id"] == 1
    assert protected[0]["reason"] == "microcontinent_or_terrane_cargo"


def test_p107_boundary_province_skeleton_classifies_process_belts():
    grid = SphereGrid.fibonacci(360, CONSTANTS.EARTH_RADIUS)
    module = TectonicsModule()
    ridge = np.where(np.abs(grid.lat) < 5.0)[0][:12]
    margin = np.where((grid.lat > 10.0) & (grid.lat < 35.0))[0][:12]
    suture = np.where((grid.lat > 40.0) & (grid.lat < 65.0))[0][:12]
    ctype = np.full(grid.n, OCEAN, dtype=np.float64)
    ctype[margin] = CONT
    ctype[suture] = CONT
    age = np.full(grid.n, 100.0, dtype=np.float64)
    age[ridge] = 5.0
    boundary_objects = [
        {
            "id": "ridge:test",
            "kind": "ridge",
            "stage": "young_ocean",
            "cells": ridge.astype(int).tolist(),
            "boundary_continental_fraction": 0.0,
            "parent_plate_ids": [0, 1],
        },
        {
            "id": "active_margin:test",
            "kind": "active_margin",
            "stage": "active_margin",
            "cells": margin.astype(int).tolist(),
            "boundary_continental_fraction": 0.8,
            "parent_plate_ids": [1, 2],
        },
        {
            "id": "suture:test",
            "kind": "suture",
            "stage": "collision_suture",
            "cells": suture.astype(int).tolist(),
            "boundary_continental_fraction": 0.9,
            "parent_plate_ids": [2, 3],
        },
    ]

    provinces, field = module._boundary_province_objects(
        grid, boundary_objects, ctype, age, 4500.0)

    kinds = {obj["kind"] for obj in provinces}
    assert {"mid_ocean_ridge", "continental_arc_margin", "continent_continent_collision"} <= kinds
    assert set(np.unique(field[ridge])) == {1.0}
    assert set(np.unique(field[margin])) == {5.0}
    assert set(np.unique(field[suture])) == {8.0}
    assert all("length_estimate_m" in obj for obj in provinces)


def test_p107_continuous_boundary_skeleton_preserves_long_ridge_network():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    module = TectonicsModule()
    allowed = np.abs(grid.lat) < 18.0
    start = int(np.argmin((grid.lat - 0.0) ** 2 + (grid.lon + 130.0) ** 2))
    target = int(np.argmin((grid.lat - 0.0) ** 2 + (grid.lon - 120.0) ** 2))

    queue = [start]
    parent = {start: -1}
    while queue and target not in parent:
        c = queue.pop(0)
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in parent or not allowed[nb]:
                continue
            parent[nb] = c
            queue.append(nb)
    assert target in parent

    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = parent[cur]
    path = np.asarray(path, dtype=np.int64)
    raw = np.zeros(grid.n, dtype=bool)
    raw[path] = True

    thinned = module._thin_boundary_mask(grid, raw, spacing=3)
    continuous = module._continuous_boundary_skeleton(
        grid, raw, allowed=allowed, bridge_passes=1)

    assert len(module._connected_components(grid, raw)) == 1
    assert len(module._connected_components(grid, thinned)) > 1
    assert len(module._connected_components(grid, continuous)) == 1
    assert int(np.count_nonzero(continuous)) >= int(0.80 * path.size)


def test_p107_continuous_boundary_skeleton_bridges_one_cell_gap():
    grid = SphereGrid.fibonacci(360, CONSTANTS.EARTH_RADIUS)
    module = TectonicsModule()
    seed = int(np.argmin(grid.lat ** 2 + grid.lon ** 2))
    path = [seed]
    seen = {seed}
    while len(path) < 16:
        cur = path[-1]
        choices = [
            int(nb) for nb in grid.neighbors[cur]
            if int(nb) not in seen and abs(float(grid.lat[int(nb)])) < 35.0
        ]
        assert choices
        nxt = min(choices, key=lambda c: (abs(float(grid.lat[c])), abs(float(grid.lon[c]))))
        path.append(nxt)
        seen.add(nxt)

    raw = np.zeros(grid.n, dtype=bool)
    raw[np.asarray(path, dtype=np.int64)] = True
    raw[path[len(path) // 2]] = False

    bridged = module._continuous_boundary_skeleton(
        grid, raw, allowed=np.ones(grid.n, dtype=bool), bridge_passes=2)

    assert len(module._connected_components(grid, raw)) == 2
    assert len(module._connected_components(grid, bridged)) == 1


def test_p107_boundary_province_terrain_response_shapes_key_belts():
    spec = get_preset("earthlike")
    spec.grid_cells = 360
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    ridge = np.where(np.abs(grid.lat) < 4.0)[0][:10]
    trench = np.where((grid.lat < -15.0) & (grid.lat > -35.0))[0][:10]
    collision = np.where((grid.lat > 30.0) & (grid.lat < 55.0))[0][:10]
    province = np.full(grid.n, -1.0, dtype=np.float64)
    province[ridge] = 1.0
    province[trench] = 4.0
    province[collision] = 8.0
    world.fields["tectonics.boundary_province_kind"] = province

    surface = np.full(grid.n, -5000.0, dtype=np.float64)
    crust = np.full(grid.n, OCEAN, dtype=np.float64)
    surface[collision] = 900.0
    crust[collision] = CONT

    shaped = terrain._p107_boundary_province_terrain_response(
        world, surface, 0.0, crust)

    assert float(np.median(shaped[ridge])) > float(np.median(surface[ridge]))
    assert float(np.median(shaped[trench])) < -5600.0
    assert float(np.median(shaped[collision])) > 1500.0
    assert world.g("terrain.last_p107_boundary_province_response_fraction") > 0.0


def test_p108_high_mountain_coherence_raises_belts_and_suppresses_isolated_peaks():
    spec = get_preset("earthlike")
    spec.grid_cells = 500
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    axis = np.where((grid.lat > -8.0) & (grid.lat < 8.0) & (grid.lon > -70.0) & (grid.lon < 70.0))[0]
    axis = axis[:24]
    assert axis.size >= 8
    far_candidates = np.where((np.abs(grid.lat) > 45.0) & (np.abs(grid.lon) > 110.0))[0]
    assert far_candidates.size
    axis_centroid = grid.xyz[axis].mean(axis=0)
    axis_centroid /= max(float(np.linalg.norm(axis_centroid)), 1.0e-12)
    isolated = int(far_candidates[int(np.argmin(grid.xyz[far_candidates] @ axis_centroid))])

    surface = np.full(grid.n, 900.0, dtype=np.float64)
    surface[axis] = 2600.0
    surface[isolated] = 5400.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    domain = np.full(grid.n, 2.0, dtype=np.float64)
    origin = np.zeros(grid.n, dtype=np.float64)
    orog_age = np.full(grid.n, -1.0, dtype=np.float64)
    orog_age[axis] = 4450.0
    province = np.full(grid.n, -1.0, dtype=np.float64)
    province[axis] = 8.0
    world.fields["tectonics.boundary_province_kind"] = province
    world.fields["terrain.orogenic_load"] = np.where(province == 8.0, 0.9, 0.0)
    world.networks["tectonics.boundaries"] = {
        "collision": axis.astype(np.int64),
        "suture": axis.astype(np.int64),
    }

    shaped = terrain._p108_high_mountain_coherence_response(
        world, surface, 0.0, crust, domain, origin, orog_age)

    assert float(np.median(shaped[axis])) > float(np.median(surface[axis]))
    assert float(shaped[isolated]) < float(surface[isolated])
    assert world.g("terrain.last_p108_high_mountain_component_count") >= 1.0
    assert world.g("terrain.last_p108_high_mountain_parent_overlap_fraction") >= 0.0


def test_p108_post_inventory_cleanup_uses_mountain_range_parentage():
    spec = get_preset("earthlike")
    spec.grid_cells = 520
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    belt = np.where((grid.lat > -10.0) & (grid.lat < 10.0) & (grid.lon > -80.0) & (grid.lon < 80.0))[0]
    belt = belt[:28]
    assert belt.size >= 10
    orphan = int(np.where((grid.lat > 45.0) & (grid.lon > 120.0))[0][0])
    surface = np.full(grid.n, 900.0, dtype=np.float64)
    surface[belt] = 3180.0
    surface[orphan] = 4100.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    ranges = np.zeros(grid.n, dtype=np.float64)
    ranges[belt] = 1.0
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[belt] = 1.0

    shaped = terrain._p108_post_inventory_high_mountain_cleanup(
        world, surface, 0.0, crust, ranges, inventory)

    assert np.all(shaped[belt] >= 3000.0)
    assert float(shaped[orphan]) < 3000.0
    assert world.g("terrain.last_p108_high_mountain_component_count") >= 1.0
    assert world.g("terrain.last_p108_high_mountain_parent_overlap_fraction") > 0.0


def test_p1113_orogen_belt_repair_connects_parented_high_mountain_clusters():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 18.0)
        & (grid.lon > -115.0)
        & (grid.lon < 115.0)
    )
    start = int(np.argmin(np.where(allowed, grid.lat ** 2 + (grid.lon + 88.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(allowed, grid.lat ** 2 + (grid.lon - 88.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    path = np.asarray(path[::-1], dtype=np.int64)
    assert path.size >= 28

    corridor = np.zeros(grid.n, dtype=bool)
    corridor[path] = True
    corridor = terrain._dilate_mask(grid, corridor, passes=1) & allowed
    high_seed = np.zeros(grid.n, dtype=bool)
    for idx in np.linspace(4, path.size - 5, 5).astype(int):
        high_seed[int(path[idx])] = True
    high_clusters = terrain._dilate_mask(grid, high_seed, passes=1) & corridor

    surface = np.full(grid.n, 850.0, dtype=np.float64)
    surface[corridor] = 2420.0
    surface[high_clusters] = 3180.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    ranges = np.zeros(grid.n, dtype=np.float64)
    ranges[corridor] = 1.0
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[corridor] = 1.0

    before_high = surface >= 3000.0
    before_components = len(terrain._components(grid, before_high))
    shaped = terrain._p111_coherent_orogen_belt_repair(
        world,
        surface,
        0.0,
        crust,
        ranges,
        inventory,
    )
    after_high = shaped >= 3000.0
    after_components = len(terrain._components(grid, after_high))

    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert after_components < before_components
    assert int(after_high.sum()) > int(before_high.sum())
    assert world.g("terrain.last_p1113_bridge_area_fraction") > 0.0
    assert world.g("terrain.last_p1113_high_component_count_after") < world.g(
        "terrain.last_p1113_high_component_count_before"
    )
    assert world.g("terrain.last_p1113_fragmentation_index_after") < world.g(
        "terrain.last_p1113_fragmentation_index_before"
    )
    assert world.g("terrain.last_p1113_land_mask_preserved") == 1.0


def test_p11169_orogen_belt_repair_uses_parent_hierarchy_without_foreland_peaks():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 18.0)
        & (grid.lon > -115.0)
        & (grid.lon < 115.0)
    )
    start = int(np.argmin(np.where(allowed, grid.lat ** 2 + (grid.lon + 88.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(allowed, grid.lat ** 2 + (grid.lon - 88.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    path = np.asarray(path[::-1], dtype=np.int64)
    assert path.size >= 28

    crest = np.zeros(grid.n, dtype=bool)
    crest[path] = True
    branch = terrain._dilate_mask(grid, crest, passes=1) & allowed & ~crest
    foreland = terrain._dilate_mask(grid, crest, passes=2) & allowed & ~crest & ~branch
    corridor = crest | branch
    high_seed = np.zeros(grid.n, dtype=bool)
    for idx in np.linspace(4, path.size - 5, 5).astype(int):
        high_seed[int(path[idx])] = True
    high_clusters = terrain._dilate_mask(grid, high_seed, passes=1) & corridor

    surface = np.full(grid.n, 850.0, dtype=np.float64)
    surface[corridor] = 2420.0
    surface[foreland] = 2620.0
    surface[high_clusters] = 3180.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    ranges = np.zeros(grid.n, dtype=np.float64)
    ranges[corridor] = 1.0
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[corridor] = 1.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[foreland] = 1.0
    hierarchy[branch] = 2.0
    hierarchy[crest] = 3.0

    before_high = surface >= 3000.0
    shaped = terrain._p111_coherent_orogen_belt_repair(
        world,
        surface,
        0.0,
        crust,
        ranges,
        inventory,
        hierarchy,
    )
    after_high = shaped >= 3000.0

    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert len(terrain._components(grid, after_high)) < len(
        terrain._components(grid, before_high)
    )
    assert int(np.count_nonzero(after_high & crest)) > int(
        np.count_nonzero(before_high & crest)
    )
    assert np.all(shaped[foreland] < 3000.0)
    assert world.g("terrain.last_p11169_hierarchy_used") == 1.0
    assert world.g("terrain.last_p11169_hierarchy_bridge_area_fraction") > 0.0
    assert world.g("terrain.last_p11169_foreland_peak_cell_count") == 0.0
    assert world.g("terrain.last_p11169_bridge_restricted_to_peak_hierarchy") == 1.0
    assert world.g("terrain.last_p11169_fragmentation_delta") > 0.0
    assert world.g("terrain.last_p111611_bridge_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111611_peak_hierarchy_shoulder_cell_count") > 0.0
    assert world.g("terrain.last_p111611_safe_path_count") > 0.0
    assert world.g("terrain.last_p111611_bridge_deferred_no_safe_path") == 0.0


def test_p111622_orogen_repair_uses_spine_to_connect_low_crest_gaps():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 16.0)
        & (grid.lon > -95.0)
        & (grid.lon < 95.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 72.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 72.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    path = np.asarray(path[::-1], dtype=np.int64)
    assert path.size >= 18

    crest = np.zeros(grid.n, dtype=bool)
    crest[path] = True
    high_seed = np.zeros(grid.n, dtype=bool)
    for idx in np.linspace(4, path.size - 5, 3).astype(int):
        high_seed[int(path[idx])] = True
    high_clusters = terrain._dilate_mask(grid, high_seed, passes=1) & crest

    surface = np.full(grid.n, 820.0, dtype=np.float64)
    surface[crest] = 980.0
    surface[high_clusters] = 3180.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    ranges = np.zeros(grid.n, dtype=np.float64)
    ranges[crest] = 1.0
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[crest] = 1.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[crest] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[crest] = 3.0

    before_high = surface >= 3000.0
    shaped = terrain._p111_coherent_orogen_belt_repair(
        world,
        surface,
        0.0,
        crust,
        ranges,
        inventory,
        hierarchy,
        spine,
    )
    after_high = shaped >= 3000.0

    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert len(terrain._components(grid, after_high)) < len(
        terrain._components(grid, before_high)
    )
    assert int(np.count_nonzero(after_high & crest)) > int(
        np.count_nonzero(before_high & crest)
    )
    assert world.g("terrain.last_p111622_spine_guided_repair_used") == 1.0
    assert world.g("terrain.last_p111622_spine_guided_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111622_spine_guided_bridge_cell_count") > 0.0
    assert world.g("terrain.last_p111622_high_component_delta") > 0.0
    assert world.g("terrain.last_p111622_fragmentation_delta") > 0.0
    assert world.g("terrain.last_p11169_foreland_peak_cell_count") == 0.0
    assert world.g("terrain.last_p1113_land_mask_preserved") == 1.0


def test_p111624_spine_to_terrain_response_lifts_continuous_crest_gaps():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 15.0)
        & (grid.lon > -72.0)
        & (grid.lon < 72.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 54.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 54.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    path = np.asarray(path[::-1], dtype=np.int64)
    assert path.size >= 12

    crest = np.zeros(grid.n, dtype=bool)
    crest[path] = True
    high_seed = np.zeros(grid.n, dtype=bool)
    high_seed[int(path[2])] = True
    high_seed[int(path[path.size // 2])] = True
    high_seed[int(path[-3])] = True
    high_clusters = terrain._dilate_mask(grid, high_seed, passes=1) & crest

    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[crest] = 920.0
    surface[high_clusters] = 3180.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    ranges = np.zeros(grid.n, dtype=np.float64)
    ranges[high_clusters] = 1.0
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[high_clusters] = 1.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[crest] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[crest] = 3.0

    before_high = surface >= 3000.0
    shaped = terrain._p111_coherent_orogen_belt_repair(
        world,
        surface,
        0.0,
        crust,
        ranges,
        inventory,
        hierarchy,
        spine,
    )
    after_high = shaped >= 3000.0

    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert len(terrain._components(grid, after_high)) < len(
        terrain._components(grid, before_high)
    )
    assert int(np.count_nonzero(after_high & crest)) > int(
        np.count_nonzero(before_high & crest)
    )
    assert world.g("terrain.last_p111624_spine_to_terrain_response_used") == 1.0
    assert world.g("terrain.last_p111624_spine_response_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111624_spine_response_bridge_cell_count") > 0.0
    assert world.g("terrain.last_p111624_fragmentation_delta") > 0.0
    assert world.g("terrain.last_p11169_foreland_peak_cell_count") == 0.0
    assert world.g("terrain.last_p1113_land_mask_preserved") == 1.0


def test_p111633_belt_relief_response_connects_lower_orogen_tiers_without_new_peaks():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 12.0)
        & (grid.lon > -58.0)
        & (grid.lon < 58.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 42.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 42.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    path = np.asarray(path[::-1], dtype=np.int64)
    assert path.size >= 10

    short_path = path[:min(path.size, 11)]
    crest = np.zeros(grid.n, dtype=bool)
    crest[short_path] = True
    expressed = np.zeros(grid.n, dtype=bool)
    expressed[short_path[2:9:2]] = True

    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[crest] = 920.0
    surface[expressed] = 2140.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    ranges = np.zeros(grid.n, dtype=np.float64)
    ranges[crest] = 1.0
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[crest] = 1.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[crest] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[crest] = 3.0

    before_1800 = crest & (surface >= 1800.0)
    before_high = surface >= 3000.0
    shaped = terrain._p111_coherent_orogen_belt_repair(
        world,
        surface,
        0.0,
        crust,
        ranges,
        inventory,
        hierarchy,
        spine,
    )
    after_1800 = crest & (shaped >= 1800.0)
    after_high = shaped >= 3000.0
    lifted_gaps = crest & ~before_1800 & (shaped >= 1800.0)

    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert len(terrain._components(grid, after_1800)) < len(
        terrain._components(grid, before_1800)
    )
    assert int(np.count_nonzero(lifted_gaps)) > 0
    assert np.all(shaped[lifted_gaps] < 3000.0)
    assert int(np.count_nonzero(after_high)) == int(np.count_nonzero(before_high))
    assert world.g("terrain.last_p111633_belt_relief_response_used") == 1.0
    assert world.g("terrain.last_p111633_belt_relief_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111633_belt_relief_bridge_cell_count") > 0.0
    assert world.g(
        "terrain.last_p111633_belt_relief_1800_component_count_after"
    ) < world.g("terrain.last_p111633_belt_relief_1800_component_count_before")
    assert world.g("terrain.last_p111633_belt_relief_guard_reverted") == 0.0
    assert world.g("terrain.last_p1113_land_mask_preserved") == 1.0


def test_p111634_spine_aligned_response_lifts_branch_axis_and_softens_orphan_highs():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 12.0)
        & (grid.lon > -55.0)
        & (grid.lon < 55.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 38.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 38.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    branch_path = np.asarray(path[::-1], dtype=np.int64)[:12]
    assert branch_path.size >= 8

    orphan_seed = int(np.argmin(np.where(
        (grid.lat > 26.0)
        & (grid.lat < 42.0)
        & (grid.lon > 82.0)
        & (grid.lon < 112.0),
        (grid.lat - 34.0) ** 2 + (grid.lon - 96.0) ** 2,
        np.inf,
    )))
    orphan_high = np.zeros(grid.n, dtype=bool)
    orphan_high[orphan_seed] = True
    orphan_high = terrain._dilate_mask(grid, orphan_high, passes=1)

    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[branch_path] = 2200.0
    surface[orphan_high] = 3260.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    ranges = np.zeros(grid.n, dtype=np.float64)
    ranges[branch_path] = 1.0
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[branch_path] = 1.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[branch_path] = 2.0
    hierarchy[orphan_high] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[branch_path] = 2.0

    before_branch_max = float(np.max(surface[branch_path]))
    before_orphan_high = orphan_high & (surface >= 3000.0)
    shaped = terrain._p111_coherent_orogen_belt_repair(
        world,
        surface,
        0.0,
        crust,
        ranges,
        inventory,
        hierarchy,
        spine,
    )

    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert float(np.max(shaped[branch_path])) > before_branch_max
    assert np.all(shaped[branch_path] < 3000.0)
    assert int(np.count_nonzero(orphan_high & (shaped >= 3000.0))) < int(
        np.count_nonzero(before_orphan_high))
    assert world.g("terrain.last_p111634_spine_aligned_response_used") == 1.0
    assert world.g("terrain.last_p111634_axis_raise_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111634_axis_raised_cell_count") > 0.0
    assert world.g("terrain.last_p111634_offaxis_high_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111634_offaxis_high_softened_cell_count") > 0.0
    assert world.g("terrain.last_p111634_high_component_count_after") <= world.g(
        "terrain.last_p111634_high_component_count_before")
    assert world.g("terrain.last_p111634_guard_reverted") == 0.0
    assert world.g("terrain.last_p1113_land_mask_preserved") == 1.0


def test_p153_terminal_high_fleck_cleanup_softens_offaxis_spine_noise():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    world.set_g("terrain.enable_p111634_spine_aligned_elevation_response", 0.0)
    world.set_g("terrain.enable_p153_terminal_high_fleck_cleanup", 1.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 12.0)
        & (grid.lon > -55.0)
        & (grid.lon < 55.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 38.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 38.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    branch_path = np.asarray(path[::-1], dtype=np.int64)[:12]
    assert branch_path.size >= 8

    fleck = int(np.argmin(np.where(
        (grid.lat > 28.0)
        & (grid.lat < 46.0)
        & (grid.lon > 86.0)
        & (grid.lon < 118.0),
        (grid.lat - 36.0) ** 2 + (grid.lon - 102.0) ** 2,
        np.inf,
    )))

    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[branch_path] = 2200.0
    surface[fleck] = 3260.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    ranges = np.zeros(grid.n, dtype=np.float64)
    ranges[branch_path] = 1.0
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[branch_path] = 1.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[branch_path] = 2.0
    hierarchy[fleck] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[branch_path] = 2.0

    shaped = terrain._p111_coherent_orogen_belt_repair(
        world,
        surface,
        0.0,
        crust,
        ranges,
        inventory,
        hierarchy,
        spine,
    )

    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert float(shaped[fleck]) < 3000.0
    assert world.g("terrain.last_p153_terminal_high_fleck_cleanup_accepted") == 1.0
    assert world.g("terrain.last_p153_candidate_cell_count") >= 1.0
    assert world.g("terrain.last_p153_softened_cell_count") >= 1.0
    assert world.g("terrain.last_p153_high_component_count_after") < world.g(
        "terrain.last_p153_high_component_count_before")
    assert world.g("terrain.last_p153_fragmentation_index_after") < world.g(
        "terrain.last_p153_fragmentation_index_before")
    assert world.g("terrain.last_p153_guard_reverted") == 0.0
    assert world.g("terrain.last_p1113_land_mask_preserved") == 1.0


def test_p155_terminal_high_relief_gate_accepts_global_ab_improvement():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 10.0)
        & (grid.lon > -45.0)
        & (grid.lon < 45.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 30.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 30.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    branch_path = np.asarray(path[::-1], dtype=np.int64)[:8]
    assert branch_path.size >= 5

    fleck = int(np.argmin(np.where(
        (grid.lat > 30.0)
        & (grid.lat < 48.0)
        & (grid.lon > 86.0)
        & (grid.lon < 118.0),
        (grid.lat - 38.0) ** 2 + (grid.lon - 102.0) ** 2,
        np.inf,
    )))

    surface = np.full(grid.n, 760.0, dtype=np.float64)
    high_a = branch_path[:2]
    saddle = np.asarray([branch_path[2]], dtype=np.int64)
    high_b = branch_path[3:5]
    surface[high_a] = 3180.0
    surface[saddle] = 2700.0
    surface[high_b] = 3160.0
    surface[fleck] = 3260.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[branch_path] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[branch_path] = 2.0

    shaped, telemetry = terrain._p155_terminal_high_relief_consistency_gate(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert float(np.min(shaped[saddle])) >= 3000.0
    assert float(shaped[fleck]) < 3000.0
    assert telemetry["raised_cell_count"] >= 1.0
    assert telemetry["softened_cell_count"] >= 1.0
    assert (
        telemetry["high_component_count_after"]
        < telemetry["high_component_count_before"]
    )
    assert (
        telemetry["fragmentation_index_after"]
        <= telemetry["fragmentation_index_before"]
    )
    assert telemetry["parent_overlap_after"] >= telemetry["parent_overlap_before"]


def test_p162_terminal_high_mountain_fragment_cleanup_lowers_small_flecks():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    main_seed = int(np.argmin(grid.lat ** 2 + (grid.lon + 20.0) ** 2))
    main = terrain._dilate_mask(
        grid,
        np.eye(1, grid.n, main_seed, dtype=bool).reshape(grid.n),
        passes=2,
    )
    main_cells = np.where(main)[0][:15]
    assert main_cells.size >= 12

    used = np.zeros(grid.n, dtype=bool)
    used[main_cells] = True
    flecks: list[int] = []
    for cell in np.argsort(np.abs(grid.lat) + 0.002 * np.abs(grid.lon)):
        cell = int(cell)
        if used[cell]:
            continue
        if np.any(used[np.asarray(grid.neighbors[cell], dtype=np.int64)]):
            continue
        flecks.append(cell)
        used[cell] = True
        for neighbor in grid.neighbors[cell]:
            used[int(neighbor)] = True
        if len(flecks) >= 10:
            break
    assert len(flecks) == 10

    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[main_cells] = 3220.0
    surface[flecks] = 3260.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[main_cells] = 3.0
    hierarchy[flecks] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[main_cells[:10]] = 3.0

    shaped, telemetry = terrain._p162_terminal_high_mountain_fragment_cleanup(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["candidate_cell_count"] >= 1.0
    assert telemetry["softened_cell_count"] >= 1.0
    assert telemetry["extreme_softened_cell_count"] == 0.0
    assert (
        telemetry["high_component_count_after"]
        < telemetry["high_component_count_before"]
    )
    assert (
        telemetry["fragmentation_index_after"]
        < telemetry["fragmentation_index_before"]
    )
    assert (
        telemetry["spine_3000_component_count_after"]
        <= telemetry["spine_3000_component_count_before"]
    )
    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert np.count_nonzero(shaped[flecks] < 3000.0) >= 1
    assert np.all(shaped[main_cells] >= 3000.0)


def test_p164_terminal_orogenic_axis_skeleton_simplification_prunes_off_axis_nodes():
    spec = get_preset("earthlike")
    spec.grid_cells = 2200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    path: list[int] | None = None
    side_cells: list[int] = []
    starts = np.where(np.abs(grid.lat) < 20.0)[0]
    for start in starts[:500]:
        trial = [int(start)]
        for _ in range(11):
            candidates = [
                int(nb)
                for nb in grid.neighbors[trial[-1]]
                if int(nb) not in trial and abs(grid.lat[int(nb)]) < 35.0
            ]
            if not candidates:
                break
            candidates.sort(key=lambda nb: (abs(grid.lat[nb]), grid.lon[nb], nb))
            trial.append(candidates[0])
        if len(trial) < 12:
            continue
        trial_set = set(trial)
        sides: list[int] = []
        for left, right in zip(trial, trial[1:]):
            common = (
                set(int(nb) for nb in grid.neighbors[left])
                & set(int(nb) for nb in grid.neighbors[right])
            )
            common -= trial_set
            common -= set(sides)
            common = {cell for cell in common if abs(grid.lat[cell]) < 40.0}
            if common:
                sides.append(min(common, key=lambda cell: (abs(grid.lat[cell]), cell)))
            if len(sides) >= 3:
                break
        if len(sides) >= 3:
            path = trial
            side_cells = sides
            break

    assert path is not None
    axis = np.asarray(path, dtype=np.int64)
    sides = np.asarray(side_cells, dtype=np.int64)
    spine_cells = np.unique(np.concatenate([axis, sides]))

    surface = np.full(grid.n, 100.0, dtype=np.float64)
    surface[spine_cells] = 2400.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[spine_cells] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[spine_cells] = 3.0

    shaped, telemetry = (
        terrain._p164_terminal_orogenic_axis_skeleton_simplification(
            world,
            surface,
            0.0,
            crust,
            hierarchy,
            spine,
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["demoted_cell_count"] > 0.0
    assert telemetry["protected_extreme_demoted_cell_count"] == 0.0
    assert telemetry["component_count_after"] <= telemetry["component_count_before"]
    assert telemetry["short_component_count_after"] <= (
        telemetry["short_component_count_before"]
    )
    assert telemetry["endpoint_count_after"] <= telemetry["endpoint_count_before"]
    assert telemetry["junction_count_after"] < telemetry["junction_count_before"]
    assert telemetry["axis_score_after"] >= telemetry["axis_score_before"] * 0.985
    assert np.all(shaped[axis] >= 3.0)
    assert np.count_nonzero(shaped[sides] <= 0.0) > 0


def test_p165_terminal_polar_edge_orogenic_semantic_compaction_prunes_weak_band():
    spec = get_preset("earthlike")
    spec.grid_cells = 2200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (grid.lat > 60.0)
        & (grid.lat < 78.0)
        & (grid.lon > -150.0)
        & (grid.lon < 100.0)
    )
    start = int(np.argmin(np.where(
        allowed,
        (grid.lat - 65.0) ** 2 + (grid.lon + 125.0) ** 2,
        np.inf,
    )))
    target = int(np.argmin(np.where(
        allowed,
        (grid.lat - 69.0) ** 2 + (grid.lon - 65.0) ** 2,
        np.inf,
    )))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        cell = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[cell]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = cell
            queue.append(nb)
    assert target in prev
    path = []
    cursor = target
    while cursor >= 0:
        path.append(cursor)
        cursor = prev[cursor]
    axis = np.asarray(path[::-1], dtype=np.int64)
    assert axis.size >= 10

    axis_mask = np.zeros(grid.n, dtype=bool)
    axis_mask[axis] = True
    peak_band = terrain._dilate_mask(grid, axis_mask, passes=4) & allowed
    outer_peak = peak_band & ~terrain._dilate_mask(grid, axis_mask, passes=2)
    assert np.count_nonzero(outer_peak) >= 8
    fringe_pool = np.where(
        (grid.lat > 66.0)
        & (grid.lat < 78.0)
        & (grid.lon > 105.0)
        & (grid.lon < 165.0)
    )[0]
    assert fringe_pool.size >= 8
    fringe_cells = fringe_pool[:10]

    surface = np.full(grid.n, 350.0, dtype=np.float64)
    surface[peak_band] = 2300.0
    surface[axis] = 3600.0
    surface[fringe_cells[:5]] = 1350.0
    surface[fringe_cells[5:]] = 2350.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[peak_band] = 2.0
    hierarchy[axis] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[axis] = 3.0
    halo = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)
    halo[fringe_cells[:5]] = 1.0
    apron[fringe_cells[5:]] = 1.0

    compacted, telemetry = (
        terrain._p165_terminal_polar_edge_orogenic_semantic_compaction(
            world,
            surface,
            0.0,
            crust,
            hierarchy,
            spine,
            halo,
            apron,
        )
    )

    new_hierarchy = compacted["terrain.orogenic_parent_hierarchy"]
    new_halo = compacted["terrain.orogenic_shoulder_halo"]
    new_apron = compacted["terrain.orogenic_highland_apron"]
    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["candidate_cell_count"] > 0.0
    assert telemetry["peak_reclassified_cell_count"] > 0.0
    assert telemetry["fringe_cleared_cell_count"] > 0.0
    assert telemetry["protected_spine_changed_cell_count"] == 0.0
    assert telemetry["protected_extreme_reclassified_cell_count"] == 0.0
    assert telemetry["polar_semantic_area_fraction_after"] < (
        telemetry["polar_semantic_area_fraction_before"]
    )
    assert telemetry["polar_peak_area_fraction_after"] < (
        telemetry["polar_peak_area_fraction_before"]
    )
    assert np.all(new_hierarchy[axis] >= 3.0)
    assert np.array_equal(spine[axis], np.full(axis.size, 3.0))
    assert np.count_nonzero(new_hierarchy[outer_peak] < 2.0) >= 1
    assert np.count_nonzero(
        (new_halo[fringe_cells] <= 0.0) & (new_apron[fringe_cells] <= 0.0)
    ) >= 1


def test_p166_terminal_straight_highland_relief_softening_breaks_weak_line():
    spec = get_preset("earthlike")
    spec.grid_cells = 2200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 22.0)
        & (grid.lon > -130.0)
        & (grid.lon < 130.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 105.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 105.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        cell = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[cell]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = cell
            queue.append(nb)
    assert target in prev
    path = []
    cursor = target
    while cursor >= 0:
        path.append(cursor)
        cursor = prev[cursor]
    line = np.asarray(path[::-1], dtype=np.int64)
    assert line.size >= 18

    protected = line[:4]
    surface = np.full(grid.n, 260.0, dtype=np.float64)
    surface[line] = 2250.0
    surface[protected] = 3500.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[protected] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[protected] = 3.0
    halo = np.zeros(grid.n, dtype=np.float64)
    apron = np.zeros(grid.n, dtype=np.float64)

    shaped, telemetry = (
        terrain._p166_terminal_straight_highland_relief_softening(
            world,
            surface,
            0.0,
            crust,
            hierarchy,
            spine,
            halo,
            apron,
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["candidate_cell_count"] > 0.0
    assert telemetry["adjusted_cell_count"] > 0.0
    assert telemetry["straightness_score_after"] < (
        telemetry["straightness_score_before"]
    )
    assert telemetry["straight_highland_area_fraction_after"] < (
        telemetry["straight_highland_area_fraction_before"]
    )
    assert telemetry["protected_core_lowered_cell_count"] == 0.0
    assert telemetry["extreme_softened_cell_count"] == 0.0
    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert np.allclose(shaped[protected], surface[protected])
    assert np.count_nonzero(shaped[line[4:]] < 900.0) >= 1


def test_p167_terminal_inland_plateau_diversity_repair_cuts_basin_pockets():
    spec = get_preset("earthlike")
    spec.grid_cells = 2200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    center = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    angular = grid.xyz @ center
    plateau = angular > np.cos(np.deg2rad(33.0))
    assert np.count_nonzero(plateau) >= 60
    protected = np.where(angular > np.cos(np.deg2rad(8.0)))[0]
    assert protected.size >= 3

    surface = np.full(grid.n, 420.0, dtype=np.float64)
    surface[plateau] = 2250.0
    surface[protected] = 3550.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    stability = np.full(grid.n, 0.35, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[protected] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[protected] = 3.0

    shaped, telemetry = (
        terrain._p167_terminal_inland_plateau_diversity_repair(
            world,
            surface,
            0.0,
            crust,
            stability,
            hierarchy,
            spine,
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["candidate_cell_count"] > 0.0
    assert telemetry["adjusted_cell_count"] > 0.0
    assert telemetry["inland_plateau_area_fraction_after"] < (
        telemetry["inland_plateau_area_fraction_before"]
    )
    assert (
        telemetry["plateau_relief_iqr_after_m"]
        >= telemetry["plateau_relief_iqr_before_m"]
    )
    assert telemetry["protected_core_lowered_cell_count"] == 0.0
    assert telemetry["extreme_softened_cell_count"] == 0.0
    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert np.allclose(shaped[protected], surface[protected])
    assert np.count_nonzero(shaped[plateau] < 1400.0) >= 1


def test_p168_terminal_surface_derasterization_diversifies_weak_surface():
    spec = get_preset("earthlike")
    spec.grid_cells = 2200
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    land = grid.lat > -18.0
    ocean = ~land
    surface = np.where(land, 720.0, -3600.0).astype(np.float64)
    land_band = land & (np.abs(grid.lon + 45.0) < 5.0) & (grid.lat < 55.0)
    ocean_band = ocean & (np.abs(grid.lon - 55.0) < 6.0)
    surface[land_band] = 1500.0
    surface[ocean_band] = -2450.0

    protected_land = np.where(
        land & (np.abs(grid.lat - 35.0) < 8.0) & (np.abs(grid.lon + 120.0) < 8.0)
    )[0]
    assert protected_land.size >= 2
    surface[protected_land] = 4650.0
    ridge_cells = np.where(ocean & (np.abs(grid.lon - 120.0) < 6.0))[0]
    trench_cells = np.where(ocean & (np.abs(grid.lon + 150.0) < 6.0))[0]
    assert ridge_cells.size >= 2
    assert trench_cells.size >= 2

    crust = np.where(land, CONT, OCEAN).astype(np.float64)
    stability = np.full(grid.n, 0.45, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[protected_land] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[protected_land] = 3.0
    depth_province = np.full(grid.n, OCEAN_DEPTH_ABYSS, dtype=np.float64)
    depth_province[land] = 0.0
    depth_province[ridge_cells] = OCEAN_DEPTH_RIDGE
    depth_province[trench_cells] = OCEAN_DEPTH_TRENCH

    shaped, telemetry = (
        terrain._p168_terminal_surface_derasterization_and_inland_diversity(
            world,
            surface,
            0.0,
            crust,
            stability,
            hierarchy,
            spine,
            {"ocean.depth_province": depth_province},
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["adjusted_cell_count"] > 0.0
    assert telemetry["land_adjusted_cell_count"] > 0.0
    assert telemetry["ocean_adjusted_cell_count"] > 0.0
    assert telemetry["weak_land_iqr_after_m"] > telemetry["weak_land_iqr_before_m"]
    assert telemetry["protected_core_changed_cell_count"] == 0.0
    assert telemetry["ridge_trench_changed_cell_count"] == 0.0
    assert telemetry["extreme_changed_cell_count"] == 0.0
    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert np.allclose(shaped[protected_land], surface[protected_land])
    assert np.allclose(shaped[ridge_cells], surface[ridge_cells])
    assert np.allclose(shaped[trench_cells], surface[trench_cells])


def test_p168_terminal_surface_derasterization_breaks_lowland_stripe_and_plateau():
    spec = get_preset("earthlike")
    spec.grid_cells = 2600
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    land = grid.lat > -34.0
    ocean = ~land
    surface = np.where(land, 1450.0, -3600.0).astype(np.float64)
    low_stripe = (
        land
        & (np.abs(grid.lon + 42.0) < 5.5)
        & (grid.lat > -18.0)
        & (grid.lat < 62.0)
    )
    assert int(np.count_nonzero(low_stripe)) >= 8
    surface[low_stripe] = 420.0

    protected_land = np.where(
        land & (np.abs(grid.lat - 38.0) < 7.0) & (np.abs(grid.lon - 95.0) < 7.0)
    )[0]
    assert protected_land.size >= 2
    surface[protected_land] = 4650.0

    crust = np.where(land, CONT, OCEAN).astype(np.float64)
    stability = np.full(grid.n, 0.42, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[protected_land] = 3.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[protected_land] = 3.0
    depth_province = np.full(grid.n, OCEAN_DEPTH_ABYSS, dtype=np.float64)
    depth_province[land] = 0.0

    shaped, telemetry = (
        terrain._p168_terminal_surface_derasterization_and_inland_diversity(
            world,
            surface,
            0.0,
            crust,
            stability,
            hierarchy,
            spine,
            {"ocean.depth_province": depth_province},
        )
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["land_mask_preserved"] is True
    assert telemetry["broad_plateau_area_fraction_after"] < (
        telemetry["broad_plateau_area_fraction_before"]
    )
    assert float(np.mean(shaped[low_stripe] - surface[low_stripe])) > 45.0
    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert np.allclose(shaped[protected_land], surface[protected_land])


def test_p166167_audit_reports_terminal_relief_planform_metrics():
    spec = get_preset("earthlike")
    spec.grid_cells = 520
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    line = np.where((np.abs(grid.lat) < 8.0) & (np.abs(grid.lon) < 120.0))[0]
    plateau = np.where(
        (grid.lat > 20.0)
        & (grid.lat < 58.0)
        & (grid.lon > -90.0)
        & (grid.lon < 40.0)
    )[0]
    assert line.size >= 8
    assert plateau.size >= 12
    elevation = np.full(grid.n, 320.0, dtype=np.float64)
    elevation[line] = 1800.0
    elevation[plateau] = 2100.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    spine = np.zeros(grid.n, dtype=np.float64)
    world.fields.update({
        "terrain.elevation_m": elevation,
        "crust.type": np.full(grid.n, CONT, dtype=np.float64),
        "terrain.orogenic_parent_hierarchy": hierarchy,
        "terrain.orogenic_hierarchy_spine": spine,
        "terrain.orogenic_shoulder_halo": np.zeros(grid.n, dtype=np.float64),
        "terrain.orogenic_highland_apron": np.zeros(grid.n, dtype=np.float64),
    })
    world.set_g(
        "terrain.last_p166_terminal_straight_highland_relief_softening_accepted",
        1.0,
    )
    world.set_g("terrain.last_p166_guard_reverted", 0.0)
    world.set_g("terrain.last_p166_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p166_candidate_cell_count", 14.0)
    world.set_g("terrain.last_p166_adjusted_cell_count", 6.0)
    world.set_g("terrain.last_p166_straightness_score_before", 0.072)
    world.set_g("terrain.last_p166_straightness_score_after", 0.041)
    world.set_g(
        "terrain.last_p166_straight_highland_area_fraction_before",
        0.031,
    )
    world.set_g(
        "terrain.last_p166_straight_highland_area_fraction_after",
        0.020,
    )
    world.set_g("terrain.last_p166_protected_core_lowered_cell_count", 0.0)
    world.set_g("terrain.last_p166_extreme_softened_cell_count", 0.0)
    world.set_g(
        "terrain.last_p167_terminal_inland_plateau_diversity_repair_accepted",
        1.0,
    )
    world.set_g("terrain.last_p167_guard_reverted", 0.0)
    world.set_g("terrain.last_p167_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p167_candidate_cell_count", 22.0)
    world.set_g("terrain.last_p167_adjusted_cell_count", 9.0)
    world.set_g("terrain.last_p167_inland_plateau_area_fraction_before", 0.083)
    world.set_g("terrain.last_p167_inland_plateau_area_fraction_after", 0.061)
    world.set_g("terrain.last_p167_largest_plateau_component_share_before", 0.74)
    world.set_g("terrain.last_p167_largest_plateau_component_share_after", 0.52)
    world.set_g("terrain.last_p167_plateau_relief_iqr_before_m", 120.0)
    world.set_g("terrain.last_p167_plateau_relief_iqr_after_m", 330.0)
    world.set_g("terrain.last_p167_protected_core_lowered_cell_count", 0.0)
    world.set_g("terrain.last_p167_extreme_softened_cell_count", 0.0)
    world.set_g(
        "terrain.last_p168_terminal_surface_derasterization_and_inland_diversity_accepted",
        1.0,
    )
    world.set_g("terrain.last_p168_guard_reverted", 0.0)
    world.set_g("terrain.last_p168_land_mask_preserved", 1.0)
    world.set_g("terrain.last_p168_candidate_cell_count", 48.0)
    world.set_g("terrain.last_p168_adjusted_cell_count", 31.0)
    world.set_g("terrain.last_p168_land_adjusted_cell_count", 20.0)
    world.set_g("terrain.last_p168_ocean_adjusted_cell_count", 11.0)
    world.set_g("terrain.last_p168_weak_land_iqr_before_m", 80.0)
    world.set_g("terrain.last_p168_weak_land_iqr_after_m", 180.0)
    world.set_g("terrain.last_p168_weak_land_edge_delta_p95_before_m", 420.0)
    world.set_g("terrain.last_p168_weak_land_edge_delta_p95_after_m", 260.0)
    world.set_g("terrain.last_p168_weak_ocean_edge_delta_p95_before_m", 580.0)
    world.set_g("terrain.last_p168_weak_ocean_edge_delta_p95_after_m", 360.0)
    world.set_g("terrain.last_p168_broad_plateau_area_fraction_before", 0.044)
    world.set_g("terrain.last_p168_broad_plateau_area_fraction_after", 0.037)
    world.set_g("terrain.last_p168_protected_core_changed_cell_count", 0.0)
    world.set_g("terrain.last_p168_ridge_trench_changed_cell_count", 0.0)
    world.set_g("terrain.last_p168_extreme_changed_cell_count", 0.0)

    metrics = terminal_audit_metrics(world)
    planform = metrics["p166167_terminal_relief_planform"]

    assert planform["schema"] == "aevum.p166167_terminal_relief_planform.v1"
    assert planform["straight_highland_component_count"] >= 0
    assert planform["inland_plateau_area_fraction_world"] >= 0.0
    assert (
        planform[
            "p166_terminal_straight_highland_relief_softening_accepted"
        ] is True
    )
    assert planform["p166_straightness_score_after"] == 0.041
    assert (
        planform[
            "p167_terminal_inland_plateau_diversity_repair_accepted"
        ] is True
    )
    assert planform["p167_inland_plateau_area_fraction_after"] == 0.061
    assert planform["p167_plateau_relief_iqr_after_m"] == 330.0
    assert (
        planform[
            "p168_terminal_surface_derasterization_and_inland_diversity_accepted"
        ] is True
    )
    assert planform["p168_adjusted_cell_count"] == 31.0
    assert planform["p168_weak_land_iqr_after_m"] == 180.0
    assert planform["p168_weak_ocean_edge_delta_p95_after_m"] == 360.0
    assert planform["p168_ridge_trench_changed_cell_count"] == 0.0
    assert metrics["acceptance"][
        "has_p166167_terminal_relief_planform_metrics"]


def test_p156_high_saddle_path_connects_multi_cell_spine_gap():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 10.0)
        & (grid.lon > -50.0)
        & (grid.lon < 50.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 34.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 34.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    branch_path = np.asarray(path[::-1], dtype=np.int64)[:9]
    assert branch_path.size >= 6

    surface = np.full(grid.n, 760.0, dtype=np.float64)
    high_a = branch_path[:2]
    saddle = branch_path[2:4]
    high_b = branch_path[4:6]
    surface[high_a] = 3180.0
    surface[saddle] = 2700.0
    surface[high_b] = 3160.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[branch_path] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[branch_path] = 2.0

    shaped, telemetry = terrain._p155_terminal_high_relief_consistency_gate(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["p156_path_candidate_count"] >= 1.0
    assert telemetry["p156_path_candidate_cell_count"] >= 2.0
    assert telemetry["p156_selected_path_count"] >= 1.0
    assert telemetry["p156_selected_path_cell_count"] >= 2.0
    assert float(np.min(shaped[saddle])) >= 3000.0
    assert (
        telemetry["spine_3000_component_count_after"]
        < telemetry["spine_3000_component_count_before"]
    )
    assert (
        telemetry["high_component_count_after"]
        < telemetry["high_component_count_before"]
    )
    assert np.array_equal(surface >= 0.0, shaped >= 0.0)


def test_p157_seed_path_extends_process_backed_spine_high_coverage():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    allowed = (
        (np.abs(grid.lat) < 10.0)
        & (grid.lon > -50.0)
        & (grid.lon < 50.0)
    )
    start = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon + 34.0) ** 2, np.inf)))
    target = int(np.argmin(np.where(
        allowed, grid.lat ** 2 + (grid.lon - 34.0) ** 2, np.inf)))
    prev = {start: -1}
    queue = [start]
    head = 0
    while head < len(queue) and target not in prev:
        c = queue[head]
        head += 1
        for nb in sorted(int(x) for x in grid.neighbors[c]):
            if nb in prev or not allowed[nb]:
                continue
            prev[nb] = c
            queue.append(nb)
    assert target in prev
    path = []
    cur = target
    while cur >= 0:
        path.append(cur)
        cur = prev[cur]
    branch_path = np.asarray(path[::-1], dtype=np.int64)[:9]
    assert branch_path.size >= 6

    surface = np.full(grid.n, 760.0, dtype=np.float64)
    high_seed = branch_path[:2]
    seed_path = branch_path[2:6]
    surface[high_seed] = 3180.0
    surface[seed_path] = 2700.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[branch_path] = 2.0
    spine = np.zeros(grid.n, dtype=np.float64)
    spine[branch_path] = 2.0

    shaped, telemetry = terrain._p155_terminal_high_relief_consistency_gate(
        world,
        surface,
        0.0,
        crust,
        hierarchy,
        spine,
    )

    assert telemetry["accepted"] is True
    assert telemetry["guard_reverted"] is False
    assert telemetry["p156_path_candidate_count"] == 0.0
    assert telemetry["p157_seed_path_candidate_count"] >= 1.0
    assert telemetry["p157_seed_path_candidate_cell_count"] >= 2.0
    assert telemetry["p157_selected_seed_path_count"] >= 1.0
    assert telemetry["p157_selected_seed_path_cell_count"] >= 2.0
    assert float(np.min(shaped[seed_path])) >= 3000.0
    assert (
        telemetry["high_component_count_after"]
        == telemetry["high_component_count_before"]
    )
    assert (
        telemetry["spine_3000_component_count_after"]
        == telemetry["spine_3000_component_count_before"]
    )
    assert (
        telemetry["spine_3000_coverage_after"]
        > telemetry["spine_3000_coverage_before"]
    )
    assert np.array_equal(surface >= 0.0, shaped >= 0.0)


def test_p111611_orogen_repair_reports_blocked_hierarchy_path_without_forced_uplift():
    spec = get_preset("earthlike")
    spec.grid_cells = 900
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    west = (np.abs(grid.lat) < 14.0) & (grid.lon > -95.0) & (grid.lon < -45.0)
    east = (np.abs(grid.lat) < 14.0) & (grid.lon > 45.0) & (grid.lon < 95.0)
    assert int(west.sum()) >= 8
    assert int(east.sum()) >= 8
    west_center = int(np.argmin(np.where(west, grid.lat ** 2 + (grid.lon + 70.0) ** 2, np.inf)))
    east_center = int(np.argmin(np.where(east, grid.lat ** 2 + (grid.lon - 70.0) ** 2, np.inf)))
    high_seed = np.zeros(grid.n, dtype=bool)
    high_seed[[west_center, east_center]] = True
    high_clusters = terrain._dilate_mask(grid, high_seed, passes=1) & (west | east)
    corridor = terrain._dilate_mask(grid, high_clusters, passes=1) & (west | east)
    crest = high_seed | (
        terrain._dilate_mask(grid, high_seed, passes=1) & corridor
    )
    branch = corridor & ~crest
    foreland = terrain._dilate_mask(grid, corridor, passes=1) & (west | east) & ~corridor

    surface = np.full(grid.n, 850.0, dtype=np.float64)
    surface[corridor] = 2420.0
    surface[foreland] = 2620.0
    surface[high_clusters] = 3180.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    ranges = np.zeros(grid.n, dtype=np.float64)
    ranges[corridor] = 1.0
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[corridor] = 1.0
    hierarchy = np.zeros(grid.n, dtype=np.float64)
    hierarchy[foreland] = 1.0
    hierarchy[branch] = 2.0
    hierarchy[crest] = 3.0

    before_high = surface >= 3000.0
    shaped = terrain._p111_coherent_orogen_belt_repair(
        world,
        surface,
        0.0,
        crust,
        ranges,
        inventory,
        hierarchy,
    )
    after_high = shaped >= 3000.0

    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert int(np.count_nonzero(after_high)) == int(np.count_nonzero(before_high))
    assert world.g("terrain.last_p11169_hierarchy_used") == 1.0
    assert world.g("terrain.last_p11169_hierarchy_bridge_area_fraction") == 0.0
    assert world.g("terrain.last_p11169_foreland_peak_cell_count") == 0.0
    assert world.g("terrain.last_p111611_bridge_candidate_cell_count") > 0.0
    assert world.g("terrain.last_p111611_peak_hierarchy_shoulder_cell_count") > 0.0
    assert world.g("terrain.last_p111611_high_pair_count") > 0.0
    assert world.g("terrain.last_p111611_safe_path_count") == 0.0
    assert world.g("terrain.last_p111611_blocked_high_pair_count") > 0.0
    assert world.g("terrain.last_p111611_bridge_deferred_no_safe_path") == 1.0


def test_p109_continental_hypsometry_rebalance_lowers_overraised_interiors():
    spec = get_preset("earthlike")
    spec.grid_cells = 640
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    crust = np.full(grid.n, CONT, dtype=np.float64)
    surface = np.full(grid.n, 2600.0, dtype=np.float64)
    detail = np.full(grid.n, 1.0, dtype=np.float64)
    basin = grid.lat < -35.0
    platform = (grid.lat >= -35.0) & (grid.lat < 20.0)
    orogen = (grid.lat >= 20.0) & (grid.lat < 55.0) & (np.abs(grid.lon) < 80.0)
    plateau = (grid.lat >= 45.0) & (np.abs(grid.lon) > 90.0)
    detail[platform] = 2.0
    detail[basin] = 3.0
    detail[orogen] = 5.0
    detail[plateau] = 6.0
    surface[basin] = 2100.0
    surface[orogen] = 4200.0
    surface[plateau] = 5200.0
    mountain_inventory = np.zeros(grid.n, dtype=np.float64)
    mountain_inventory[orogen] = 1.0
    mountain_inventory[plateau] = 3.0

    shaped = terrain._p109_continental_hypsometry_rebalance(
        world,
        surface,
        0.0,
        crust,
        detail,
        mountain_inventory,
    )

    ordinary = mountain_inventory <= 0.0
    basin_mean = float(np.mean(shaped[basin & ordinary]))
    platform_mean = float(np.mean(shaped[platform & ordinary]))
    orogen_mean = float(np.mean(shaped[orogen]))
    plateau_mean = float(np.mean(shaped[plateau]))
    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert float(np.mean(shaped[ordinary])) < float(np.mean(surface[ordinary])) - 1200.0
    assert float(np.percentile(shaped, 50)) < 900.0
    assert basin_mean < platform_mean < orogen_mean < plateau_mean
    assert np.any(shaped[orogen] >= 1800.0)
    assert np.any(shaped[plateau] >= 3000.0)
    assert np.count_nonzero(shaped >= 4500.0) > 0
    assert np.count_nonzero(shaped >= 4500.0) < np.count_nonzero(surface >= 4500.0)
    assert world.g("terrain.last_p109_land_mask_preserved") == 1.0
    assert 0.0 < world.g("terrain.last_p109_extreme_retained_area_fraction") < 0.02
    assert world.g("terrain.last_p109_land_mean_after_m") < world.g(
        "terrain.last_p109_land_mean_before_m")


def test_p109_post_p110a_rebalance_restores_parented_highland_shoulder():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    crust = np.full(grid.n, CONT, dtype=np.float64)
    surface = np.full(grid.n, 620.0, dtype=np.float64)
    detail = np.full(grid.n, 2.0, dtype=np.float64)
    shoulder = (
        (grid.lat > -18.0)
        & (grid.lat < 42.0)
        & (np.abs(grid.lon) < 62.0)
    )
    plateau = (grid.lat > 52.0) & (np.abs(grid.lon) > 100.0)
    assert np.mean(shoulder) > 0.05
    assert np.mean(plateau) > 0.015

    surface[shoulder] = 1750.0
    surface[plateau] = 4300.0
    detail[plateau] = 6.0
    mountain_inventory = np.zeros(grid.n, dtype=np.float64)
    mountain_inventory[shoulder] = 1.0
    mountain_inventory[plateau] = 3.0

    shaped = terrain._p109_continental_hypsometry_rebalance(
        world,
        surface,
        0.0,
        crust,
        detail,
        mountain_inventory,
        post_p110a_planform=True,
    )

    area = grid.cell_area
    rel = shaped
    land = shaped >= 0.0
    land_area = float(area[land].sum())
    p90 = terrain._weighted_percentile_local(rel[land], area[land], 90.0)
    shoulder_band = (
        land
        & (rel >= 2000.0)
        & (rel < 3000.0)
    )
    ordinary = mountain_inventory <= 0.0
    assert np.array_equal(surface >= 0.0, shaped >= 0.0)
    assert p90 >= 1700.0
    assert float(area[shoulder_band].sum() / land_area) >= 0.01
    assert float(np.mean(shaped[ordinary])) < 850.0
    assert world.g("terrain.last_p109_post_p110a_tail_area_fraction") > 0.0
    assert world.g("terrain.last_p109_post_p110a_shoulder_area_fraction") > 0.0


def test_p109_final_p90_floor_polish_repairs_near_miss_without_planform_change():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    crust = np.full(grid.n, CONT, dtype=np.float64)
    surface = np.full(grid.n, 620.0, dtype=np.float64)
    detail = np.full(grid.n, 2.0, dtype=np.float64)
    inventory = np.zeros(grid.n, dtype=np.float64)

    order = np.argsort(-grid.lat)
    area = grid.cell_area
    total = float(area.sum())
    cumulative = np.cumsum(area[order])
    high = order[cumulative <= 0.085 * total]
    support = order[(cumulative > 0.085 * total) & (cumulative <= 0.140 * total)]
    assert high.size > 0
    assert support.size > 0

    surface[high] = 1760.0
    surface[support] = 1640.0
    detail[high] = 5.0
    detail[support] = 5.0
    inventory[high] = 1.0
    inventory[support] = 1.0
    p90_before = terrain._weighted_percentile_local(surface, area, 90.0)
    assert p90_before < 1700.0

    polished = terrain._p109_final_p90_floor_polish(
        world,
        surface,
        0.0,
        crust,
        detail,
        inventory,
    )

    p90_after = terrain._weighted_percentile_local(polished, area, 90.0)
    assert np.array_equal(surface >= 0.0, polished >= 0.0)
    assert p90_after >= 1700.0
    assert float(np.max(polished)) < 2100.0
    assert world.g("terrain.last_p109_final_p90_floor_area_fraction") > 0.0
    assert world.g("terrain.last_p109_final_p90_floor_before_m") == p90_before
    assert world.g("terrain.last_p109_final_p90_floor_after_m") == p90_after
    assert world.g("terrain.last_p109_land_p90_after_m") == p90_after


def test_p109_final_lowland_shoulder_polish_repairs_low_band_without_planform_change():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    area = grid.cell_area
    total = float(area.sum())
    order = np.argsort(grid.lat)
    cumulative = np.cumsum(area[order])
    low = order[cumulative <= 0.382 * total]
    mid = order[(cumulative > 0.382 * total) & (cumulative <= 0.620 * total)]
    high = order[cumulative > 0.860 * total]

    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[low] = 158.0
    surface[mid] = 330.0
    surface[high] = 1850.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    detail = np.full(grid.n, 2.0, dtype=np.float64)
    detail[low[: max(1, low.size // 5)]] = 4.0
    inventory = np.zeros(grid.n, dtype=np.float64)

    land = surface >= 0.0
    land_area = float(area[land].sum())
    low_before = float(area[land & (surface >= 0.0) & (surface < 200.0)].sum() / land_area)
    assert low_before > 0.35

    polished = terrain._p109_final_lowland_shoulder_polish(
        world,
        surface,
        0.0,
        crust,
        detail,
        inventory,
    )

    low_after = float(area[land & (polished >= 0.0) & (polished < 200.0)].sum() / land_area)
    mid_after = float(area[land & (polished >= 200.0) & (polished < 500.0)].sum() / land_area)
    assert np.array_equal(surface >= 0.0, polished >= 0.0)
    assert low_after <= 0.350
    assert mid_after > float(area[land & (surface >= 200.0) & (surface < 500.0)].sum() / land_area)
    assert world.g("terrain.last_p109_final_lowland_shoulder_area_fraction") > 0.0
    assert world.g("terrain.last_p109_final_lowland_shoulder_before_fraction") == low_before
    assert world.g("terrain.last_p109_final_lowland_shoulder_after_fraction") == low_after


def test_p109_final_lowland_shoulder_polish_uses_near_sea_level_fallback():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    area = grid.cell_area
    total = float(area.sum())
    order = np.argsort(grid.lat)
    cumulative = np.cumsum(area[order])
    low = order[cumulative <= 0.382 * total]

    surface = np.full(grid.n, 760.0, dtype=np.float64)
    surface[low] = 45.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    detail = np.full(grid.n, 2.0, dtype=np.float64)
    inventory = np.zeros(grid.n, dtype=np.float64)

    land = surface >= 0.0
    land_area = float(area[land].sum())
    low_before = float(
        area[land & (surface >= 0.0) & (surface < 200.0)].sum() / land_area)
    assert low_before > 0.35

    polished = terrain._p109_final_lowland_shoulder_polish(
        world,
        surface,
        0.0,
        crust,
        detail,
        inventory,
    )

    low_after = float(
        area[land & (polished >= 0.0) & (polished < 200.0)].sum() / land_area)
    assert np.array_equal(surface >= 0.0, polished >= 0.0)
    assert low_after < low_before
    assert world.g("terrain.last_p109_final_lowland_shoulder_area_fraction") > 0.0


def test_p109_final_p90_ceiling_polish_repairs_near_miss_without_planform_change():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.n_plates = 36
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    area = grid.cell_area
    total = float(area.sum())
    order = np.argsort(-grid.lat)
    cumulative = np.cumsum(area[order])
    high = order[cumulative <= 0.105 * total]
    summit = order[cumulative <= 0.012 * total]

    surface = np.full(grid.n, 620.0, dtype=np.float64)
    surface[high] = 2760.0
    surface[summit] = 4550.0
    crust = np.full(grid.n, CONT, dtype=np.float64)
    detail = np.full(grid.n, 2.0, dtype=np.float64)
    detail[summit] = 6.0
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[summit] = 3.0

    p90_before = terrain._weighted_percentile_local(surface, area, 90.0)
    assert p90_before > 2700.0

    polished = terrain._p109_final_p90_ceiling_polish(
        world,
        surface,
        0.0,
        crust,
        detail,
        inventory,
    )

    p90_after = terrain._weighted_percentile_local(polished, area, 90.0)
    assert np.array_equal(surface >= 0.0, polished >= 0.0)
    assert p90_after <= 2700.0
    assert float(np.max(polished[summit])) >= 4500.0
    assert world.g("terrain.last_p109_final_p90_ceiling_area_fraction") > 0.0
    assert world.g("terrain.last_p109_final_p90_ceiling_before_m") == p90_before
    assert world.g("terrain.last_p109_final_p90_ceiling_after_m") == p90_after


def test_p107_underresolved_plate_hierarchy_restoration_splits_large_plates():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    module = TectonicsModule()
    plate = np.zeros(grid.n, dtype=np.int64)
    plates = [
        {"id": pid, "pole": [0.0, 0.0, 1.0], "rate": 0.008}
        for pid in range(8)
    ]
    context = {
        "ctype": np.full(grid.n, OCEAN, dtype=np.float64),
        "origin": np.full(grid.n, ORIGIN_RIDGE, dtype=np.float64),
        "age": np.full(grid.n, 60.0, dtype=np.float64),
        "rift_potential": np.clip(1.0 - np.abs(grid.lat) / 90.0, 0.0, 1.0),
        "boundaries": {
            "ridge": np.where(np.abs(grid.lat) < 8.0)[0].astype(np.int64),
            "transform": np.where(np.abs(grid.lon) < 20.0)[0].astype(np.int64),
        },
        "plate_topologies": [],
        "enable_p107_ranked_plate_policy": True,
    }

    splits = module._split_underresolved_plate_hierarchy(
        grid,
        plate,
        plates,
        free_ids=list(range(1, 8)),
        t=4500.0,
        context=context,
        target_active=5,
    )

    assert len(splits) >= 4
    assert len(set(int(x) for x in np.unique(plate))) >= 5
    assert all(split["basis"] == "p107_ranked_hierarchy_restoration" for split in splits)
