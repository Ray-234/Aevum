import json
from types import SimpleNamespace

import numpy as np

from aevum.archive.world_archive import Frame, P171_REQUIRED_OBJECT_FIELDS, WorldArchive
from aevum.core.events import EventBus
from aevum.core.grid import SphereGrid
from aevum.core.state import WorldState
from aevum.core.units import CONSTANTS
from aevum.diagnostics.historical_geomorphology import (
    LAND_METRIC_KEYS,
    OCEAN_METRIC_KEYS,
    historical_geomorphology_summary,
    frame_geomorphology_metrics,
    write_historical_geomorphology_audit,
)
from aevum.diagnostics.historical_objects import (
    historical_object_persistence_summary,
    write_historical_object_audit,
)
from aevum.modules.tectonics import (
    DOMAIN_ACCRETED_TERRANE,
    DOMAIN_LIP,
    DOMAIN_OCEANIC,
    ORIGIN_PLUME_IMPACT,
)
from aevum.modules.terrain import (
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATFORM,
    CONT_DETAIL_RIFT_BASIN,
    CONT_DETAIL_SHIELD,
    CONT_PROVINCE_ACTIVE_OROGEN,
    CONT_PROVINCE_FORELAND_BASIN,
    CONT_PROVINCE_INTRACRATONIC_BASIN,
    CONT_PROVINCE_PASSIVE_MARGIN_LOWLAND,
    CONT_PROVINCE_PLATFORM,
    INLAND_PROVINCE_OLD_OROGEN,
    INLAND_PROVINCE_PLATFORM,
    INLAND_PROVINCE_RIFT,
    INLAND_PROVINCE_SAG_BASIN,
    INLAND_PROVINCE_SHIELD,
    OCEAN_DEPTH_ABYSS,
    OCEAN_DEPTH_RIDGE,
    OCEAN_DEPTH_TRENCH,
    RIFT_MARGIN_STAGE_ESCARPMENT,
    RIFT_MARGIN_STAGE_PASSIVE_LOWLAND,
    RIFT_MARGIN_STAGE_RIFT_BASIN,
    RIFT_MARGIN_STAGE_SHOULDER,
)
from aevum.spec.presets import get_preset


def _world_and_archive():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    spec.seed = 170
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    land = grid.lat > -8.0
    ocean = ~land
    area = grid.cell_area

    bland_elev = np.where(land, 1450.0, -4300.0).astype(np.float64)
    bland_fields = {
        "terrain.elevation_m": bland_elev,
        "crust.type": land.astype(np.float64),
        "crust.age_myr": np.where(land, 1800.0, 90.0),
        "crust.domain": np.where(land, 1.0, DOMAIN_OCEANIC).astype(np.float64),
        "crust.origin": np.zeros(grid.n, dtype=np.float64),
        "terrain.continental_detail": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.continental_detail_region_code": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.inland_geomorphology_region_code": np.where(
            land, INLAND_PROVINCE_PLATFORM, 0).astype(np.float64),
        "terrain.continental_province_code": np.where(land, 2.0, 0.0),
        "tectonics.orogeny_age_myr": np.full(grid.n, -1.0, dtype=np.float64),
        "terrain.rift_margin_stage": np.zeros(grid.n, dtype=np.float64),
        "ocean.depth_province": np.where(
            ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.float64),
    }

    diverse_elev = np.full(grid.n, -4400.0, dtype=np.float64)
    diverse_elev[land] = 650.0
    shield = land & (grid.lon < -70.0)
    basin = land & (grid.lon >= -70.0) & (grid.lon < -10.0)
    rift = land & (np.abs(grid.lon) < 10.0)
    orogen = land & (grid.lon > 30.0) & (grid.lon < 115.0)
    diverse_elev[shield] = 520.0
    diverse_elev[basin] = 180.0
    diverse_elev[rift] = 420.0
    diverse_elev[orogen] = 2350.0 + 4.0 * np.maximum(grid.lat[orogen], 0.0)

    ridge = ocean & (np.abs(grid.lon) < 8.0)
    trench = ocean & (grid.lon > 115.0) & (grid.lon < 145.0)
    plateau = ocean & (grid.lon < -135.0) & (grid.lat < -35.0)
    microcontinent = ocean & (grid.lon > 65.0) & (grid.lon < 95.0) & (grid.lat < -35.0)
    diverse_elev[ridge] = -2450.0
    diverse_elev[trench] = -6200.0
    diverse_elev[plateau] = -1150.0
    diverse_elev[microcontinent] = -900.0

    detail = np.zeros(grid.n, dtype=np.float64)
    detail[shield] = CONT_DETAIL_SHIELD
    detail[basin] = CONT_DETAIL_BASIN
    detail[rift] = CONT_DETAIL_RIFT_BASIN
    detail[orogen] = CONT_DETAIL_OROGEN
    detail[land & (detail == 0.0)] = CONT_DETAIL_PLATFORM

    inland_region = np.zeros(grid.n, dtype=np.float64)
    inland_region[shield] = INLAND_PROVINCE_SHIELD
    inland_region[basin] = INLAND_PROVINCE_SAG_BASIN
    inland_region[rift] = INLAND_PROVINCE_RIFT
    inland_region[orogen] = INLAND_PROVINCE_OLD_OROGEN
    inland_region[land & (inland_region == 0.0)] = INLAND_PROVINCE_PLATFORM

    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0).astype(np.float64)
    depth_province[ridge] = OCEAN_DEPTH_RIDGE
    depth_province[trench] = OCEAN_DEPTH_TRENCH

    crust_type = land.astype(np.float64)
    crust_type[microcontinent] = 1.0
    crust_domain = np.where(ocean, DOMAIN_OCEANIC, 1.0).astype(np.float64)
    crust_domain[plateau] = DOMAIN_LIP
    origin = np.zeros(grid.n, dtype=np.float64)
    origin[plateau] = ORIGIN_PLUME_IMPACT
    age = np.where(ocean, 120.0, 1800.0).astype(np.float64)
    age[ridge] = 0.0
    age[ocean & (grid.lon < 0.0)] = np.minimum(age[ocean & (grid.lon < 0.0)], 70.0)
    sediment = np.full(grid.n, 120.0, dtype=np.float64)
    sediment[basin | rift] = 900.0

    diverse_fields = {
        "terrain.elevation_m": diverse_elev,
        "crust.type": crust_type,
        "crust.age_myr": age,
        "crust.domain": crust_domain,
        "crust.origin": origin,
        "terrain.continental_detail": detail,
        "terrain.continental_detail_region_code": detail,
        "terrain.inland_geomorphology_region_code": inland_region,
        "terrain.continental_province_code": np.where(land, detail + 10.0, 0.0),
        "tectonics.orogeny_age_myr": np.where(orogen, 900.0, -1.0),
        "terrain.rift_margin_stage": np.where(rift, 2.0, 0.0),
        "sediment.thickness_m": sediment,
        "ocean.depth_province": depth_province,
    }

    frames = [
        Frame(
            time_myr=1500.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=bland_fields,
            diagnostics={
                "terrain": {
                    "ocean_fabric_kind_counts": {},
                    "arc_plume_landform_kind_counts": {},
                }
            },
        ),
        Frame(
            time_myr=4500.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=diverse_fields,
            diagnostics={
                "terrain": {
                    "ocean_fabric_kind_counts": {
                        "spreading_center": int(np.count_nonzero(ridge) > 0),
                        "fracture_zone": 2,
                    },
                    "arc_plume_landform_kind_counts": {
                        "hotspot_track": 1,
                        "seamount_chain": 1,
                        "oceanic_plateau": 1,
                    },
                }
            },
        ),
    ]
    assert float(area[land].sum()) > 0.0
    return world, SimpleNamespace(world=world, frames=frames)


def test_p170_historical_geomorphology_summary_is_deterministic_and_key_complete():
    world, archive = _world_and_archive()

    first = historical_geomorphology_summary(world, archive)
    second = historical_geomorphology_summary(world, archive)

    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
    assert first["schema"] == "aevum.p170_historical_geomorphology.v1"
    assert first["usable_frame_count"] == 2
    assert first["required_metric_keys_present"]
    assert first["acceptance"]["audit_completed"]
    assert "p174_continuity" in first
    assert first["p174_continuity"]["mature_support_frame_count"] >= 1
    assert "terminal_lowland_pop_candidate" in first["p174_continuity"]
    for row in first["frame_rows"]:
        assert set(LAND_METRIC_KEYS).issubset(row["land_metrics"])
        assert set(OCEAN_METRIC_KEYS).issubset(row["ocean_metrics"])

    bland, diverse = first["frame_rows"]
    assert bland["diagnostic_flags"]["ordinary_plateau_like"]
    assert bland["diagnostic_flags"]["ordinary_deep_ocean_like"]
    assert diverse["land_metrics"]["inland_detail_entropy"] > bland[
        "land_metrics"]["inland_detail_entropy"]
    assert diverse["land_metrics"]["old_orogen_expression_fraction"] > 0.0
    assert diverse["land_metrics"]["rift_basin_expression_fraction"] > 0.0
    assert 0.0 <= diverse["land_metrics"][
        "craton_shield_platform_split_fraction"] <= 1.0
    assert bland["land_metrics"]["lowland_plain_fraction"] == 0.0
    assert diverse["land_metrics"]["lowland_plain_fraction"] > 0.0
    assert diverse["land_metrics"]["broad_lowland_plain_component_count"] > 0
    assert diverse["land_metrics"]["lowland_plain_parented_fraction"] > 0.5
    assert bland["diagnostic_flags"]["lowland_plain_deficient"]
    assert not diverse["diagnostic_flags"]["lowland_plain_deficient"]
    assert diverse["ocean_metrics"]["ridge_visible_fraction"] > 0.0
    assert diverse["ocean_metrics"]["hotspot_track_count"] == 1
    assert diverse["ocean_metrics"]["seamount_chain_count"] == 2
    assert diverse["ocean_metrics"]["oceanic_plateau_fraction"] > 0.0


def test_p170_land_metrics_include_p104f_p174_frame_global_telemetry():
    world, archive = _world_and_archive()
    archive.frames[0].globals.update({
        "terrain.last_p104f_pre_p174_lowland_prep_area_fraction": 0.125,
        "terrain.last_p104f_pre_p174_lowland_prep_mean_lowering_m": 480.0,
        "terrain.last_p104f_pre_p174_lowland_support_area_fraction": 0.214,
        "terrain.last_p104f_pre_p174_lowland_support_largest_component_fraction": 0.119,
        "terrain.last_p104f_pre_p174_lowland_memory_seed_area_fraction": 0.151,
        "terrain.last_p174_lowland_plain_response_area_fraction": 0.091,
        "terrain.last_p174_lowland_plain_candidate_area_fraction": 0.144,
        "terrain.last_p174_lowland_plain_parent_area_fraction": 0.233,
        "terrain.last_p174_lowland_plain_continuity_memory_area_fraction": 0.087,
        "terrain.last_p174_lowland_plain_continuity_parent_area_fraction": 0.192,
        "terrain.last_p174_lowland_plain_response_mean_abs_delta_m": 310.0,
        "terrain.last_p174_lowland_plain_fraction_before": 0.018,
        "terrain.last_p174_lowland_plain_fraction_after": 0.107,
        "terrain.last_p174_lowland_plain_largest_component_fraction_after": 0.083,
        "terrain.last_p174_lowland_plain_parented_fraction_after": 0.94,
        "terrain.last_p174_lowland_plain_response_stage_code": 3.0,
        "terrain.last_p174_support_component_response_area_fraction": 0.052,
        "terrain.last_p174_support_component_response_largest_component_fraction": 0.031,
        "terrain.last_p174_component_gap_repair_infill_area_fraction": 0.006,
        "terrain.last_p174_component_gap_repair_diagnostic_connector_domain_area_fraction": 0.021,
        "terrain.last_p174_component_gap_repair_diagnostic_connector_response_area_fraction": 0.005,
        "terrain.last_p174_component_gap_repair_diagnostic_connector_path_count": 2.0,
        "terrain.last_p174_component_gap_repair_diagnostic_connector_target_component_area_fraction": 0.018,
        "terrain.last_p174_component_gap_repair_largest_growth_domain_area_fraction": 0.012,
        "terrain.last_p174_component_gap_repair_largest_growth_picked_area_fraction": 0.004,
    })

    summary = historical_geomorphology_summary(world, archive)
    metrics = summary["frame_rows"][0]["land_metrics"]

    assert metrics["p104f_pre_p174_lowland_prep_area_fraction"] == 0.125
    assert metrics["p104f_pre_p174_lowland_prep_mean_lowering_m"] == 480.0
    assert metrics["p104f_pre_p174_lowland_support_area_fraction"] == 0.214
    assert (
        metrics["p104f_pre_p174_lowland_support_largest_component_fraction"]
        == 0.119
    )
    assert metrics["p104f_pre_p174_lowland_memory_seed_area_fraction"] == 0.151
    assert metrics["p174_lowland_plain_response_area_fraction"] == 0.091
    assert metrics["p174_lowland_plain_candidate_area_fraction"] == 0.144
    assert metrics["p174_lowland_plain_parent_area_fraction"] == 0.233
    assert metrics["p174_lowland_plain_continuity_memory_area_fraction"] == 0.087
    assert metrics["p174_lowland_plain_continuity_parent_area_fraction"] == 0.192
    assert metrics["p174_lowland_plain_response_mean_abs_delta_m"] == 310.0
    assert metrics["p174_lowland_plain_fraction_before"] == 0.018
    assert metrics["p174_lowland_plain_fraction_after"] == 0.107
    assert metrics["p174_lowland_plain_largest_component_fraction_after"] == 0.083
    assert metrics["p174_lowland_plain_parented_fraction_after"] == 0.94
    assert metrics["p174_support_component_response_area_fraction"] == 0.052
    assert (
        metrics["p174_support_component_response_largest_component_fraction"]
        == 0.031
    )
    assert metrics["p174_lowland_plain_response_stage_code"] == 3.0
    assert metrics["p174_component_gap_repair_infill_area_fraction"] == 0.006
    assert (
        metrics[
            "p174_component_gap_repair_diagnostic_connector_domain_area_fraction"
        ]
        == 0.021
    )
    assert (
        metrics[
            "p174_component_gap_repair_diagnostic_connector_response_area_fraction"
        ]
        == 0.005
    )
    assert metrics[
        "p174_component_gap_repair_diagnostic_connector_path_count"
    ] == 2.0
    assert (
        metrics[
            "p174_component_gap_repair_diagnostic_connector_target_component_area_fraction"
        ]
        == 0.018
    )
    assert (
        metrics["p174_component_gap_repair_largest_growth_domain_area_fraction"]
        == 0.012
    )
    assert (
        metrics["p174_component_gap_repair_largest_growth_picked_area_fraction"]
        == 0.004
    )
    assert (
        summary["metric_extremes"]["land_metrics"][
            "p104f_pre_p174_lowland_prep_area_fraction"
        ]["max"]
        == 0.125
    )
    assert summary["required_metric_keys_present"]


def test_p170_parented_shoal_reads_margin_and_backarc_objects():
    spec = get_preset("earthlike")
    spec.grid_cells = 360
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    land = grid.lat > 72.0
    margin_core = np.where(
        (~land) & (np.abs(grid.lat) < 8.0) & (grid.lon < -10.0) & (grid.lon > -50.0)
    )[0][:3]
    backarc_core = np.where(
        (~land) & (np.abs(grid.lat) < 8.0) & (grid.lon > 10.0) & (grid.lon < 50.0)
    )[0][:3]
    rift_core = np.where(
        (~land) & (np.abs(grid.lat) < 8.0) & (grid.lon > 80.0) & (grid.lon < 120.0)
    )[0][:3]
    patch_margin = np.zeros(grid.n, dtype=bool)
    patch_backarc = np.zeros(grid.n, dtype=bool)
    patch_rift = np.zeros(grid.n, dtype=bool)
    for cell in margin_core:
        patch_margin[int(cell)] = True
        patch_margin[np.asarray(grid.neighbors[int(cell)], dtype=int)] = True
    for cell in backarc_core:
        patch_backarc[int(cell)] = True
        patch_backarc[np.asarray(grid.neighbors[int(cell)], dtype=int)] = True
    for cell in rift_core:
        patch_rift[int(cell)] = True
        patch_rift[np.asarray(grid.neighbors[int(cell)], dtype=int)] = True
    patch_margin &= ~land
    patch_backarc &= ~land
    patch_rift &= ~land
    shallow = patch_margin | patch_backarc | patch_rift
    assert margin_core.size > 0
    assert backarc_core.size > 0
    assert rift_core.size > 0
    assert np.count_nonzero(patch_margin) > margin_core.size
    assert np.count_nonzero(patch_backarc) > backarc_core.size
    assert np.count_nonzero(patch_rift) > rift_core.size

    elev = np.where(land, 650.0, -4200.0).astype(np.float64)
    elev[shallow] = -900.0
    fields = {
        "terrain.elevation_m": elev,
        "crust.type": land.astype(np.float64),
        "crust.age_myr": np.where(land, 1800.0, 80.0),
        "crust.domain": np.where(land, 1.0, DOMAIN_OCEANIC).astype(np.float64),
        "crust.origin": np.zeros(grid.n, dtype=np.float64),
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS).astype(np.float64),
    }
    unparented_frame = Frame(
        time_myr=1900.0,
        globals={"ocean.sea_level_m": 0.0},
        fields=fields,
    )
    parented_frame = Frame(
        time_myr=1900.0,
        globals={"ocean.sea_level_m": 0.0},
        fields=fields,
        objects={
            "terrain.margin_landforms": [
                {
                    "kind": "volcanic_arc",
                    "cells": margin_core.astype(int).tolist(),
                },
            ],
            "terrain.arc_plume_landforms": [
                {
                    "kind": "back_arc_basin",
                    "cells": backarc_core.astype(int).tolist(),
                },
            ],
            "terrain.rift_margin_sequences": [
                {
                    "kind": "rift_margin_sequence",
                    "cells": rift_core.astype(int).tolist(),
                },
            ],
        },
    )

    unparented = frame_geomorphology_metrics(grid, unparented_frame)
    parented = frame_geomorphology_metrics(grid, parented_frame)

    assert unparented["ocean_metrics"]["unparented_shoal_fraction"] > 0.0
    assert parented["ocean_metrics"]["unparented_shoal_fraction"] == 0.0


def test_p170_handles_missing_fields_and_writes_json(tmp_path):
    spec = get_preset("earthlike")
    spec.grid_cells = 120
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=40.0)
    elev = np.where(grid.lat > 0.0, 1200.0, -4100.0)
    archive = SimpleNamespace(
        world=world,
        frames=[
            Frame(
                time_myr=40.0,
                globals={"ocean.sea_level_m": 0.0},
                fields={"terrain.elevation_m": elev},
                diagnostics={},
            )
        ],
    )

    summary = write_historical_geomorphology_audit(world, archive, tmp_path)
    path = tmp_path / "p170_historical_geomorphology_audit.json"

    assert path.exists()
    loaded = json.loads(path.read_text())
    assert loaded == summary
    assert summary["usable_frame_count"] == 1
    assert summary["required_metric_keys_present"]
    assert "terrain.continental_detail" in summary["frame_rows"][0]["missing_fields"]
    assert "terrain.continental_detail" in WorldArchive.DEFAULT_KEYS
    assert "terrain.inland_geomorphology_region_code" in WorldArchive.DEFAULT_KEYS
    assert "terrain.rift_margin_stage" in WorldArchive.DEFAULT_KEYS
    assert "sediment.thickness_m" in WorldArchive.DEFAULT_KEYS
    assert "terrain.p104f_pre_p174_lowland_support_mask" in WorldArchive.DEFAULT_KEYS
    assert "terrain.p174_lowland_plain_continuity_memory" in WorldArchive.DEFAULT_KEYS


def test_p170_rift_metric_counts_only_rift_basin_stage():
    spec = get_preset("earthlike")
    spec.grid_cells = 160
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    land = grid.lat > 0.0
    fields = {
        "terrain.elevation_m": np.where(land, 700.0, -3600.0),
        "crust.type": land.astype(np.float64),
        "crust.domain": np.zeros(grid.n, dtype=np.float64),
        "terrain.continental_detail": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.continental_detail_region_code": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.inland_geomorphology_region_code": np.where(
            land, INLAND_PROVINCE_PLATFORM, 0).astype(np.float64),
        "terrain.continental_province_code": np.where(land, 2.0, 0.0),
        "tectonics.orogeny_age_myr": np.full(grid.n, -1.0, dtype=np.float64),
        "terrain.rift_margin_stage": np.where(
            land, RIFT_MARGIN_STAGE_PASSIVE_LOWLAND, 0).astype(np.float64),
    }

    for ignored_stage in (
        RIFT_MARGIN_STAGE_PASSIVE_LOWLAND,
        RIFT_MARGIN_STAGE_SHOULDER,
        RIFT_MARGIN_STAGE_ESCARPMENT,
    ):
        fields["terrain.rift_margin_stage"] = np.where(
            land, ignored_stage, 0).astype(np.float64)
        row = frame_geomorphology_metrics(
            grid,
            Frame(
                time_myr=1800.0,
                globals={"ocean.sea_level_m": 0.0},
                fields=fields,
                diagnostics={},
            ),
        )
        assert row["land_metrics"]["rift_basin_expression_fraction"] == 0.0

    fields["terrain.rift_margin_stage"] = np.where(
        land, RIFT_MARGIN_STAGE_RIFT_BASIN, 0).astype(np.float64)
    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1800.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=fields,
            diagnostics={},
        ),
    )
    assert row["land_metrics"]["rift_basin_expression_fraction"] > 0.0


def test_p174_lowland_plain_metrics_track_broad_parented_lowlands():
    spec = get_preset("earthlike")
    spec.grid_cells = 720
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    land = grid.lat > -25.0
    lowland = (
        land
        & (grid.lon > -145.0)
        & (grid.lon < -40.0)
        & (grid.lat > -5.0)
        & (grid.lat < 48.0)
    )
    assert np.count_nonzero(lowland) >= 8

    base_elev = np.where(land, 1450.0, -3900.0).astype(np.float64)
    parented_elev = base_elev.copy()
    parented_elev[lowland] = 160.0 + 0.8 * np.maximum(grid.lat[lowland], 0.0)
    unparented_elev = parented_elev.copy()

    fields = {
        "terrain.elevation_m": parented_elev,
        "crust.type": land.astype(np.float64),
        "crust.domain": np.where(land, 1.0, DOMAIN_OCEANIC).astype(np.float64),
        "terrain.continental_detail": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.continental_detail_region_code": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.inland_geomorphology_region_code": np.where(
            land, INLAND_PROVINCE_PLATFORM, 0).astype(np.float64),
        "terrain.continental_province_code": np.where(
            land, CONT_PROVINCE_PLATFORM, 0).astype(np.float64),
        "tectonics.orogeny_age_myr": np.full(grid.n, -1.0, dtype=np.float64),
        "terrain.rift_margin_stage": np.zeros(grid.n, dtype=np.float64),
        "sediment.thickness_m": np.where(lowland, 850.0, 80.0).astype(np.float64),
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS).astype(np.float64),
    }
    fields["terrain.continental_detail"][lowland] = CONT_DETAIL_BASIN
    fields["terrain.continental_detail_region_code"][lowland] = CONT_DETAIL_BASIN
    fields["terrain.inland_geomorphology_region_code"][lowland] = (
        INLAND_PROVINCE_SAG_BASIN)
    fields["terrain.continental_province_code"][lowland] = (
        CONT_PROVINCE_INTRACRATONIC_BASIN)

    parented = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=2400.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=fields,
            diagnostics={},
        ),
    )

    high_fields = {key: np.asarray(value).copy() for key, value in fields.items()}
    high_fields["terrain.elevation_m"] = base_elev
    high = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=2400.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=high_fields,
            diagnostics={},
        ),
    )

    unparented_fields = {key: np.asarray(value).copy() for key, value in fields.items()}
    unparented_fields["terrain.elevation_m"] = unparented_elev
    for key in (
        "terrain.continental_detail",
        "terrain.continental_detail_region_code",
        "terrain.inland_geomorphology_region_code",
        "terrain.continental_province_code",
        "sediment.thickness_m",
    ):
        unparented_fields[key][lowland] = 0.0
    unparented = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=2400.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=unparented_fields,
            diagnostics={},
        ),
    )
    memory_parented_fields = {
        key: np.asarray(value).copy() for key, value in unparented_fields.items()
    }
    memory_parented_fields["terrain.p174_lowland_plain_continuity_memory"] = (
        lowland.astype(np.float64)
    )
    memory_parented = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=2400.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=memory_parented_fields,
            diagnostics={},
        ),
    )
    rough_parented_fields = {
        key: np.asarray(value).copy() for key, value in fields.items()
    }
    rough_parented_fields["terrain.elevation_m"][lowland] = 180.0
    rough_stripes = lowland & (
        (np.floor((grid.lon + 180.0) / 18.0) % 2) == 0
    )
    rough_parented_fields["terrain.elevation_m"][rough_stripes] = 690.0
    rough_parented = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=2400.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=rough_parented_fields,
            diagnostics={},
        ),
    )

    assert high["land_metrics"]["lowland_plain_fraction"] == 0.0
    assert high["diagnostic_flags"]["lowland_plain_deficient"]
    assert parented["land_metrics"]["lowland_plain_fraction"] > 0.08
    assert parented["land_metrics"]["broad_lowland_plain_component_count"] >= 1
    assert parented["land_metrics"]["largest_lowland_plain_component_fraction"] > 0.04
    assert parented["land_metrics"]["lowland_plain_parented_fraction"] > 0.85
    assert not parented["diagnostic_flags"]["lowland_plain_deficient"]
    assert unparented["land_metrics"]["lowland_plain_fraction"] > 0.08
    assert unparented["land_metrics"]["lowland_plain_parented_fraction"] < 0.30
    assert unparented["diagnostic_flags"]["lowland_plain_deficient"]
    assert unparented["land_metrics"]["lowland_residual_parentage_limited"] == 1.0
    assert unparented["land_metrics"]["lowland_residual_dominant_code"] == 2.0
    assert memory_parented["land_metrics"]["lowland_plain_parented_fraction"] > 0.85
    assert not memory_parented["diagnostic_flags"]["lowland_plain_deficient"]
    assert rough_parented["land_metrics"][
        "lowland_elevation_parented_fraction"
    ] > rough_parented["land_metrics"]["lowland_plain_fraction"]
    assert rough_parented["land_metrics"][
        "largest_lowland_elevation_parented_component_fraction"
    ] >= rough_parented["land_metrics"][
        "largest_lowland_plain_component_fraction"
    ]
    assert rough_parented["land_metrics"][
        "lowland_local_relief_blocked_fraction"
    ] > 0.0


def test_p170_lowland_residual_attribution_marks_relief_boundary_split():
    spec = get_preset("earthlike")
    spec.grid_cells = 1600
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    land = (
        (grid.lon > -175.0)
        & (grid.lon < 140.0)
        & (grid.lat > -62.0)
        & (grid.lat < 66.0)
    )
    lowland = (
        land
        & (grid.lon > -150.0)
        & (grid.lon < 70.0)
        & (grid.lat > -38.0)
        & (grid.lat < 46.0)
    )
    ribs = lowland & ((np.floor((grid.lon + 180.0) / 9.0) % 3) == 0)
    assert np.count_nonzero(lowland) >= 450
    assert np.count_nonzero(ribs) >= 90

    elevation = np.where(land, 1400.0, -4100.0).astype(np.float64)
    elevation[lowland] = 180.0
    elevation[ribs] = 690.0
    fields = {
        "terrain.elevation_m": elevation,
        "crust.type": land.astype(np.float64),
        "crust.domain": np.where(land, 1.0, DOMAIN_OCEANIC).astype(np.float64),
        "terrain.continental_detail": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.continental_detail_region_code": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.inland_geomorphology_region_code": np.where(
            land, INLAND_PROVINCE_PLATFORM, 0).astype(np.float64),
        "terrain.continental_province_code": np.where(
            land, CONT_PROVINCE_PLATFORM, 0).astype(np.float64),
        "tectonics.orogeny_age_myr": np.full(grid.n, -1.0, dtype=np.float64),
        "terrain.rift_margin_stage": np.zeros(grid.n, dtype=np.float64),
        "sediment.thickness_m": np.where(lowland, 950.0, 80.0).astype(np.float64),
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS).astype(np.float64),
    }
    for key in (
        "terrain.continental_detail",
        "terrain.continental_detail_region_code",
    ):
        fields[key][lowland] = CONT_DETAIL_BASIN
    fields["terrain.inland_geomorphology_region_code"][lowland] = (
        INLAND_PROVINCE_SAG_BASIN)
    fields["terrain.continental_province_code"][lowland] = (
        CONT_PROVINCE_INTRACRATONIC_BASIN)
    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=2100.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=fields,
            diagnostics={},
        ),
    )
    land_metrics = row["land_metrics"]

    assert row["diagnostic_flags"]["lowland_plain_deficient"]
    assert land_metrics["lowland_plain_fraction"] >= 0.06
    assert land_metrics["largest_lowland_plain_component_fraction"] < 0.025
    assert (
        land_metrics["largest_lowland_elevation_parented_component_fraction"]
        >= 0.025
    )
    assert land_metrics["lowland_relief_boundary_gap_fraction"] > 0.003
    assert land_metrics["lowland_residual_relief_boundary_limited"] == 1.0
    assert land_metrics["lowland_residual_dominant_code"] == 4.0


def test_p170_lowland_classifier_keeps_active_province_edge_from_splitting_plain():
    spec = get_preset("earthlike")
    spec.grid_cells = 1100
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    land = (
        (grid.lon > -170.0)
        & (grid.lon < 110.0)
        & (grid.lat > -52.0)
        & (grid.lat < 58.0)
    )
    plain = (
        land
        & (grid.lon > -120.0)
        & (grid.lon < 35.0)
        & (grid.lat > -28.0)
        & (grid.lat < 34.0)
    )
    active_edge = plain & (
        (np.floor((grid.lon + 180.0) / 22.0) % 5) == 0
    )
    stable_plain = plain & ~active_edge
    assert np.count_nonzero(stable_plain) >= 120
    assert np.count_nonzero(active_edge) >= 12

    fields = {
        "terrain.elevation_m": np.where(land, 1400.0, -4100.0).astype(np.float64),
        "crust.type": land.astype(np.float64),
        "crust.domain": np.where(land, 1.0, DOMAIN_OCEANIC).astype(np.float64),
        "terrain.continental_detail": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.continental_detail_region_code": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.inland_geomorphology_region_code": np.where(
            land, INLAND_PROVINCE_PLATFORM, 0).astype(np.float64),
        "terrain.continental_province_code": np.where(
            land, CONT_PROVINCE_PLATFORM, 0).astype(np.float64),
        "tectonics.orogeny_age_myr": np.full(grid.n, -1.0, dtype=np.float64),
        "terrain.rift_margin_stage": np.zeros(grid.n, dtype=np.float64),
        "sediment.thickness_m": np.where(plain, 980.0, 80.0).astype(np.float64),
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS).astype(np.float64),
    }
    fields["terrain.elevation_m"][stable_plain] = 180.0
    fields["terrain.elevation_m"][active_edge] = 670.0
    for key in (
        "terrain.continental_detail",
        "terrain.continental_detail_region_code",
    ):
        fields[key][plain] = CONT_DETAIL_BASIN
    fields["terrain.inland_geomorphology_region_code"][plain] = (
        INLAND_PROVINCE_SAG_BASIN)
    fields["terrain.continental_province_code"][plain] = (
        CONT_PROVINCE_INTRACRATONIC_BASIN)
    fields["terrain.continental_province_code"][active_edge] = (
        CONT_PROVINCE_ACTIVE_OROGEN)

    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1800.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=fields,
            diagnostics={},
        ),
    )
    land_metrics = row["land_metrics"]

    assert land_metrics["lowland_active_exclusion_fraction"] == 0.0
    assert land_metrics["lowland_relief_boundary_gap_fraction"] < 0.003
    assert land_metrics["lowland_residual_relief_boundary_limited"] == 0.0
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.025
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]


def test_p170_active_lowland_context_accepts_relative_sediment_signal():
    spec = get_preset("earthlike")
    spec.grid_cells = 1100
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    land = (
        (grid.lon > -170.0)
        & (grid.lon < 115.0)
        & (grid.lat > -52.0)
        & (grid.lat < 58.0)
    )
    plain = (
        land
        & (grid.lon > -125.0)
        & (grid.lon < 25.0)
        & (grid.lat > -28.0)
        & (grid.lat < 34.0)
    )
    active_edge = plain & (
        (np.floor((grid.lon + 180.0) / 24.0) % 5) == 0
    )
    stable_plain = plain & ~active_edge
    high_sediment_background = land & (grid.lon > 35.0)
    assert np.count_nonzero(stable_plain) >= 100
    assert np.count_nonzero(active_edge) >= 10
    assert np.count_nonzero(high_sediment_background) >= np.count_nonzero(land) * 0.25

    sediment = np.where(high_sediment_background, 180.0, 80.0).astype(np.float64)
    sediment[active_edge] = 100.0
    fields = {
        "terrain.elevation_m": np.where(land, 1400.0, -4100.0).astype(np.float64),
        "crust.type": land.astype(np.float64),
        "crust.domain": np.where(land, 1.0, DOMAIN_OCEANIC).astype(np.float64),
        "terrain.continental_detail": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.continental_detail_region_code": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.inland_geomorphology_region_code": np.where(
            land, INLAND_PROVINCE_PLATFORM, 0).astype(np.float64),
        "terrain.continental_province_code": np.where(
            land, CONT_PROVINCE_PLATFORM, 0).astype(np.float64),
        "tectonics.orogeny_age_myr": np.full(grid.n, -1.0, dtype=np.float64),
        "terrain.rift_margin_stage": np.zeros(grid.n, dtype=np.float64),
        "sediment.thickness_m": sediment,
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS).astype(np.float64),
    }
    fields["terrain.elevation_m"][stable_plain] = 180.0
    fields["terrain.elevation_m"][active_edge] = 620.0
    for key in (
        "terrain.continental_detail",
        "terrain.continental_detail_region_code",
    ):
        fields[key][stable_plain] = CONT_DETAIL_BASIN
    fields["terrain.inland_geomorphology_region_code"][stable_plain] = (
        INLAND_PROVINCE_SAG_BASIN)
    fields["terrain.continental_province_code"][stable_plain] = (
        CONT_PROVINCE_INTRACRATONIC_BASIN)
    fields["terrain.continental_province_code"][active_edge] = (
        CONT_PROVINCE_ACTIVE_OROGEN)

    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1800.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=fields,
            diagnostics={},
        ),
    )
    land_metrics = row["land_metrics"]

    assert land_metrics["lowland_active_exclusion_fraction"] == 0.0
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.025
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]


def test_p170_accreted_active_margin_basin_context_counts_as_lowland_plain():
    spec = get_preset("earthlike")
    spec.grid_cells = 1200
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    land = (
        (grid.lon > -170.0)
        & (grid.lon < 120.0)
        & (grid.lat > -54.0)
        & (grid.lat < 60.0)
    )
    plain = (
        land
        & (grid.lon > -136.0)
        & (grid.lon < 38.0)
        & (grid.lat > -30.0)
        & (grid.lat < 34.0)
    )
    accreted_edge = plain & (
        (np.floor((grid.lon + 180.0) / 26.0) % 4) == 0
    )
    stable_plain = plain & ~accreted_edge
    true_highland = land & (grid.lon > 70.0) & (grid.lat > 10.0)
    assert np.count_nonzero(stable_plain) >= 120
    assert np.count_nonzero(accreted_edge) >= 14
    assert np.count_nonzero(true_highland) >= 12

    fields = {
        "terrain.elevation_m": np.where(land, 1450.0, -4100.0).astype(np.float64),
        "crust.type": land.astype(np.float64),
        "crust.domain": np.where(
            land, 1.0, DOMAIN_OCEANIC).astype(np.float64),
        "terrain.continental_detail": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.continental_detail_region_code": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.inland_geomorphology_region_code": np.where(
            land, INLAND_PROVINCE_PLATFORM, 0).astype(np.float64),
        "terrain.continental_province_code": np.where(
            land, CONT_PROVINCE_PLATFORM, 0).astype(np.float64),
        "tectonics.orogeny_age_myr": np.full(grid.n, -1.0, dtype=np.float64),
        "terrain.rift_margin_stage": np.zeros(grid.n, dtype=np.float64),
        "sediment.thickness_m": np.where(plain, 1150.0, 90.0).astype(np.float64),
        "terrain.p174_lowland_plain_continuity_memory": np.where(
            accreted_edge, 0.45, 0.0).astype(np.float64),
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS).astype(np.float64),
    }
    fields["terrain.elevation_m"][stable_plain] = 180.0
    fields["terrain.elevation_m"][accreted_edge] = 650.0
    fields["terrain.elevation_m"][true_highland] = 2350.0
    fields["crust.domain"][accreted_edge | true_highland] = (
        DOMAIN_ACCRETED_TERRANE)
    fields["terrain.continental_detail"][stable_plain] = CONT_DETAIL_BASIN
    fields["terrain.continental_detail_region_code"][stable_plain] = (
        CONT_DETAIL_BASIN)
    fields["terrain.inland_geomorphology_region_code"][plain] = (
        INLAND_PROVINCE_SAG_BASIN)
    fields["terrain.continental_province_code"][stable_plain] = (
        CONT_PROVINCE_INTRACRATONIC_BASIN)
    fields["terrain.continental_province_code"][accreted_edge | true_highland] = (
        CONT_PROVINCE_ACTIVE_OROGEN)
    fields["terrain.continental_detail"][accreted_edge | true_highland] = (
        CONT_DETAIL_OROGEN)
    fields["terrain.continental_detail_region_code"][
        accreted_edge | true_highland
    ] = CONT_DETAIL_OROGEN

    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1600.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=fields,
            diagnostics={},
        ),
    )
    land_metrics = row["land_metrics"]

    assert land_metrics["lowland_active_exclusion_fraction"] == 0.0
    assert land_metrics["lowland_plain_fraction"] >= 0.06
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.025
    assert land_metrics["lowland_plain_parented_fraction"] > 0.80
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]


def test_p174_lowland_parentage_accepts_foreland_and_passive_margin_plains():
    spec = get_preset("earthlike")
    spec.grid_cells = 480
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    land = grid.lat > -20.0
    foreland = land & (grid.lon > -130.0) & (grid.lon < -50.0) & (grid.lat > 5.0)
    passive = land & (grid.lon > 45.0) & (grid.lon < 145.0) & (grid.lat < 45.0)
    plain = foreland | passive
    assert np.count_nonzero(foreland) >= 4
    assert np.count_nonzero(passive) >= 4

    fields = {
        "terrain.elevation_m": np.where(land, 1200.0, -4100.0).astype(np.float64),
        "crust.type": land.astype(np.float64),
        "crust.domain": np.where(land, 1.0, DOMAIN_OCEANIC).astype(np.float64),
        "terrain.continental_detail": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.continental_detail_region_code": np.where(
            land, CONT_DETAIL_PLATFORM, 0).astype(np.float64),
        "terrain.inland_geomorphology_region_code": np.where(
            land, INLAND_PROVINCE_PLATFORM, 0).astype(np.float64),
        "terrain.continental_province_code": np.where(
            land, CONT_PROVINCE_PLATFORM, 0).astype(np.float64),
        "tectonics.orogeny_age_myr": np.full(grid.n, -1.0, dtype=np.float64),
        "terrain.rift_margin_stage": np.zeros(grid.n, dtype=np.float64),
        "sediment.thickness_m": np.full(grid.n, 120.0, dtype=np.float64),
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS).astype(np.float64),
    }
    fields["terrain.elevation_m"][plain] = 210.0
    fields["terrain.continental_detail"][foreland] = CONT_DETAIL_BASIN
    fields["terrain.continental_detail_region_code"][foreland] = CONT_DETAIL_BASIN
    fields["terrain.continental_province_code"][foreland] = (
        CONT_PROVINCE_FORELAND_BASIN)
    fields["terrain.continental_province_code"][passive] = (
        CONT_PROVINCE_PASSIVE_MARGIN_LOWLAND)
    fields["terrain.rift_margin_stage"][passive] = RIFT_MARGIN_STAGE_PASSIVE_LOWLAND

    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=3300.0,
            globals={"ocean.sea_level_m": 0.0},
            fields=fields,
            diagnostics={},
        ),
    )

    assert row["land_metrics"]["lowland_plain_fraction"] > 0.08
    assert row["land_metrics"]["lowland_plain_parented_fraction"] > 0.90
    assert row["land_metrics"]["broad_lowland_plain_component_count"] >= 1


def test_p171_archive_captures_process_object_snapshots_for_p170():
    spec = get_preset("earthlike")
    spec.grid_cells = 160
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=1000.0)
    land = grid.lat > 35.0
    ridge = (~land) & (np.abs(grid.lon) < 20.0)
    world.fields.update({
        "terrain.elevation_m": np.where(land, 600.0, -3600.0),
        "crust.type": land.astype(np.float64),
        "crust.age_myr": np.where(ridge, 0.0, np.where(land, 1200.0, 90.0)),
        "crust.domain": np.zeros(grid.n, dtype=np.float64),
        "crust.origin": np.zeros(grid.n, dtype=np.float64),
        "sediment.thickness_m": np.where(land, 520.0, 140.0),
        "terrain.p174_lowland_plain_continuity_memory": land.astype(np.float64),
        "ocean.depth_province": np.where(ridge, OCEAN_DEPTH_RIDGE, OCEAN_DEPTH_ABYSS),
    })
    fracture_cells = np.where((~land) & ~ridge)[0][:4]
    assert fracture_cells.size == 4
    world.objects["terrain.ocean_fabric"] = [
        {"id": "ridge:1", "kind": "spreading_center", "cells": np.where(ridge)[0].tolist()},
        {"id": "fracture:1", "kind": "fracture_zone", "cells": fracture_cells.tolist()},
    ]
    world.objects["terrain.arc_plume_landforms"] = [
        {
            "type": "hotspot_track",
            "cell": 12,
            "formation_myr": 900.0,
            "plate_id": 4,
            "parent_process": "plume:12",
        }
    ]
    archive = WorldArchive(world, EventBus())

    frame = archive.capture(diagnostics={})
    world.objects["terrain.ocean_fabric"][0]["kind"] = "mutated_after_capture"
    summary = historical_geomorphology_summary(world, archive)

    assert "sediment.thickness_m" in frame.fields
    assert "sediment.thickness_m" not in summary["frame_rows"][0]["missing_fields"]
    assert "terrain.p174_lowland_plain_continuity_memory" in frame.fields
    assert "terrain.ocean_fabric" in WorldArchive.DEFAULT_OBJECT_KEYS
    assert frame.objects["terrain.ocean_fabric"][0]["kind"] == "spreading_center"
    assert "birth_myr" not in world.objects["terrain.ocean_fabric"][0]
    for obj in frame.objects["terrain.ocean_fabric"]:
        assert set(P171_REQUIRED_OBJECT_FIELDS).issubset(obj)
        assert obj["p171_required_fields_present"]
    hotspot = frame.objects["terrain.arc_plume_landforms"][0]
    assert hotspot["kind"] == "hotspot_track"
    assert hotspot["cells"] == [12]
    assert hotspot["birth_myr"] == 900.0
    assert hotspot["age_myr"] == 100.0
    assert hotspot["parent_plate_id"] == 4
    assert hotspot["parent_process_id"] == "plume:12"
    assert set(P171_REQUIRED_OBJECT_FIELDS).issubset(hotspot)
    assert {"id", "kind", "cells", "age_myr"}.issubset(
        set(hotspot["p171_synthesized_fields"]))
    ocean = summary["frame_rows"][0]["ocean_metrics"]
    assert ocean["ridge_visible_fraction"] > 0.0
    assert ocean["fracture_zone_length_fraction"] > 0.0
    assert ocean["hotspot_track_count"] == 1
    assert summary["required_metric_keys_present"]


def test_p171_object_persistence_audit_counts_required_fields_and_recurring_ids(tmp_path):
    spec = get_preset("earthlike")
    spec.grid_cells = 120
    grid = SphereGrid.fibonacci(spec.grid_cells, CONSTANTS.EARTH_RADIUS)
    world = WorldState(grid=grid, spec=spec, time_myr=1000.0)
    ridge_cells = np.where(np.abs(grid.lon) < 12.0)[0].tolist()
    world.objects["terrain.ocean_fabric"] = [
        {
            "id": "ridge:stable",
            "kind": "spreading_center",
            "cells": ridge_cells,
            "birth_myr": 850.0,
            "parent_process_id": "spreading:center",
            "parent_plate_id": "plate:1+plate:2",
            "lineage_id": "ridge:lineage",
        }
    ]
    archive = WorldArchive(world, EventBus())
    first = archive.capture(diagnostics={})

    world.time_myr = 1100.0
    world.objects["terrain.ocean_fabric"][0]["last_active_myr"] = 1100.0
    world.objects["terrain.ocean_fabric"].append({
        "kind": "fracture_zone",
        "cells": [2, 4, 6],
        "mean_age_myr": 70.0,
    })
    second = archive.capture(diagnostics={})

    summary = write_historical_object_audit(world, archive, tmp_path)
    loaded = json.loads(
        (tmp_path / "p171_historical_object_persistence_audit.json").read_text())

    assert loaded == summary
    assert summary["schema"] == "aevum.p171_historical_object_persistence.v1"
    assert summary["frame_count"] == 2
    assert summary["total_object_observations"] == 3
    assert summary["required_fields_complete"]
    assert summary["missing_required_field_slot_count"] == 0
    assert summary["recurring_object_id_count"] == 1
    assert summary["acceptance"]["persistence_checked"]
    assert summary["acceptance"]["recurring_object_ids_present"]
    assert first.objects["terrain.ocean_fabric"][0]["age_myr"] == 150.0
    assert second.objects["terrain.ocean_fabric"][0]["age_myr"] == 250.0
    assert second.objects["terrain.ocean_fabric"][1]["p171_required_fields_present"]
    assert set(P171_REQUIRED_OBJECT_FIELDS).issubset(
        second.objects["terrain.ocean_fabric"][1])
    direct = historical_object_persistence_summary(world, archive)
    assert direct == summary
