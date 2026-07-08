from types import SimpleNamespace

import numpy as np

from aevum.archive.world_archive import Frame
from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.diagnostics.historical_geomorphology import frame_geomorphology_metrics
from aevum.modules.terrain import (
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_PROVINCE_ACTIVE_OROGEN,
    CONT_PROVINCE_FORELAND_BASIN,
    CONT_PROVINCE_INTRACRATONIC_BASIN,
    CONT_PROVINCE_OLD_OROGEN,
    CONT_PROVINCE_PASSIVE_MARGIN_LOWLAND,
    CONT_PROVINCE_PLATFORM,
    CONT_DETAIL_PLATFORM,
    DOMAIN_ACCRETED_TERRANE,
    DOMAIN_CONTINENTAL_INTERIOR,
    DOMAIN_LIP,
    DOMAIN_OCEANIC,
    INTERNAL_BLOCK_STABLE_PLATFORM,
    INLAND_PROVINCE_SAG_BASIN,
    INLAND_PROVINCE_PLATFORM,
    OCEAN_DEPTH_ABYSS,
    OCEAN_DEPTH_SHELF,
    OCEAN_MARGIN_NONE,
    RIFT_MARGIN_STAGE_ESCARPMENT,
    RIFT_MARGIN_STAGE_RIFT_BASIN,
    RIFT_MARGIN_STAGE_SHOULDER,
    TerrainModule,
)
from aevum.modules.tectonics import _same_neighbor_count as tectonics_same_neighbor_count


def test_same_neighbor_count_matches_bruteforce_for_sparse_masks():
    grid = SphereGrid.fibonacci(240, CONSTANTS.EARTH_RADIUS)
    mask = (
        ((grid.lon > -80.0) & (grid.lon < -20.0) & (grid.lat > -30.0))
        | ((grid.lon > 40.0) & (grid.lon < 85.0) & (grid.lat < 35.0))
    )
    expected = np.zeros(grid.n, dtype=np.int16)
    for cell in range(grid.n):
        expected[cell] = sum(1 for nb in grid.neighbors[cell] if mask[int(nb)])

    module = TerrainModule()
    np.testing.assert_array_equal(module._same_neighbor_count(grid, mask), expected)
    np.testing.assert_array_equal(tectonics_same_neighbor_count(grid, mask), expected)


def test_overlapping_object_ids_from_index_vectorized_view_matches_sparse_index():
    module = TerrainModule()
    index = {
        2: ("arc:1", "ridge:3"),
        5: ("arc:1",),
        8: ("trench:2",),
    }
    mask = np.zeros(12, dtype=bool)
    mask[[2, 5]] = True

    assert module._overlapping_object_ids_from_index(mask, index) == [
        "arc:1",
        "ridge:3",
    ]

    mask[:] = False
    assert module._overlapping_object_ids_from_index(mask, index) == []

    mask[[8, 10]] = True
    assert module._overlapping_object_ids_from_index(mask, index) == ["trench:2"]


def test_p172_process_object_lifecycle_reuses_previous_id_and_birth():
    grid = SphereGrid.fibonacci(80, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    previous = {
        "id": "landform:old_subdued_orogen:stable",
        "kind": "old_subdued_orogen",
        "cells": [10, 11, 12],
        "birth_myr": 900.0,
        "lineage_id": "landform:old_subdued_orogen:stable",
        "parent_process": "eroded_orogen",
        "parent_continent_ids": [3],
    }
    world = SimpleNamespace(
        grid=grid,
        time_myr=1100.0,
        objects={"terrain.continental_landforms": [previous]},
    )
    current = [{
        "id": "old_subdued_orogen:10:0",
        "kind": "old_subdued_orogen",
        "cells": [11, 12, 13],
        "parent_process": "eroded_orogen",
        "parent_continent_ids": [3],
    }]

    stabilized = module._stabilize_process_objects(
        world,
        "terrain.continental_landforms",
        current,
        id_prefix="landform",
    )

    assert stabilized[0]["id"] == "landform:old_subdued_orogen:stable"
    assert stabilized[0]["birth_myr"] == 900.0
    assert stabilized[0]["age_myr"] == 200.0
    assert stabilized[0]["lineage_id"] == "landform:old_subdued_orogen:stable"
    assert stabilized[0]["activity_state"] == "decaying"
    assert stabilized[0]["relief_stage"] == "eroded_orogen"


def test_p172_process_object_lifecycle_synthesizes_new_stable_metadata():
    grid = SphereGrid.fibonacci(80, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    world = SimpleNamespace(
        grid=grid,
        time_myr=1200.0,
        objects={"terrain.continental_landforms": []},
    )
    current = [{
        "id": "rift_basin:20:0",
        "kind": "rift_basin",
        "cells": [20, 21, 22],
        "mean_age_myr": 75.0,
        "parent_process_code": 7,
        "parent_continent_ids": [4, 5],
    }]

    first = module._stabilize_process_objects(
        world,
        "terrain.continental_landforms",
        current,
        id_prefix="landform",
    )[0]
    second = module._stabilize_process_objects(
        world,
        "terrain.continental_landforms",
        current,
        id_prefix="landform",
    )[0]

    assert first["id"].startswith("landform:rift_basin:")
    assert first["id"] == second["id"]
    assert first["birth_myr"] == 1125.0
    assert first["age_myr"] == 75.0
    assert first["parent_process_id"] == "7"
    assert first["parent_plate_id"] == "4+5"
    assert first["activity_state"] == "active"
    assert first["relief_stage"] == "rift_relief"


def test_p172_rift_margin_sequence_does_not_promote_all_inland_shoulder_candidates():
    grid = SphereGrid.fibonacci(220, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    land = grid.lat > 0.0
    surface = np.where(land, 520.0, -3600.0).astype(np.float64)
    crust_type = land.astype(np.float64)
    crust_domain = np.zeros(grid.n, dtype=np.float64)
    crust_stability = np.full(grid.n, 0.72, dtype=np.float64)
    continental_detail = np.where(
        land, CONT_DETAIL_PLATFORM, 0).astype(np.float64)
    fields = {
        "tectonics.orogeny_age_myr": np.full(grid.n, -1.0, dtype=np.float64),
        "tectonics.deformation_style": np.zeros(grid.n, dtype=np.float64),
        "tectonics.deformation_intensity": np.zeros(grid.n, dtype=np.float64),
    }
    world = SimpleNamespace(
        grid=grid,
        time_myr=1800.0,
        objects={},
        networks={},
        get_field=lambda name, default=0.0: fields.get(
            name,
            np.full(grid.n, default, dtype=np.float64),
        ),
        set_g=lambda *args, **kwargs: None,
    )

    def fake_inland_state(*args, **kwargs):
        empty = np.zeros(grid.n, dtype=bool)
        return {
            "province": np.zeros(grid.n, dtype=np.float64),
            "regional_relief_m": np.zeros(grid.n, dtype=np.float64),
            "shield": empty,
            "platform": empty,
            "sag_basin": empty,
            "old_orogen_root": empty,
            "rift_axis": empty,
            "rift_shoulder": land.copy(),
            "platform_swell": empty,
            "escarpment": empty,
            "plateau_margin": empty,
        }

    module._inland_geomorphology_state = fake_inland_state
    ocean_geo_fields = {
        "ocean.depth_province": np.where(
            land, 0, OCEAN_DEPTH_ABYSS).astype(np.float64),
        "ocean.margin_type": np.full(
            grid.n, OCEAN_MARGIN_NONE, dtype=np.float64),
        "ocean.shelf_width": np.where(
            land, 0, OCEAN_DEPTH_SHELF).astype(np.float64),
    }

    result = module._production_bathymetry_margin_sequence(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        np.zeros(grid.n, dtype=np.float64),
        crust_stability,
        continental_detail,
        np.zeros(grid.n, dtype=bool),
        [],
        [],
        ocean_geo_fields,
    )

    stage = result["fields"]["terrain.rift_margin_stage"].astype(int)
    rift_like = np.isin(stage, [
        RIFT_MARGIN_STAGE_SHOULDER,
        RIFT_MARGIN_STAGE_RIFT_BASIN,
        RIFT_MARGIN_STAGE_ESCARPMENT,
    ])
    land_area = max(float(grid.cell_area[land].sum()), 1.0e-12)
    assert float(grid.cell_area[land & rift_like].sum() / land_area) < 0.35


def test_p172_inland_state_does_not_turn_broad_rift_potential_into_all_rift():
    grid = SphereGrid.fibonacci(220, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    crust_type = np.ones(grid.n, dtype=np.float64)
    crust_domain = np.full(
        grid.n, DOMAIN_CONTINENTAL_INTERIOR, dtype=np.float64)
    crust_stability = np.full(grid.n, 0.78, dtype=np.float64)
    fields = {
        "crust.thickness_m": np.full(grid.n, 35000.0, dtype=np.float64),
        "crust.age_myr": np.full(grid.n, 1600.0, dtype=np.float64),
        "tectonics.rift_potential": np.full(grid.n, 0.65, dtype=np.float64),
        "tectonics.platform_subsidence": np.zeros(grid.n, dtype=np.float64),
        "tectonics.deformation_style": np.zeros(grid.n, dtype=np.float64),
        "tectonics.deformation_intensity": np.zeros(grid.n, dtype=np.float64),
    }
    world = SimpleNamespace(
        grid=grid,
        time_myr=2200.0,
        objects={},
        networks={},
        get_field=lambda name, default=0.0: fields.get(
            name,
            np.full(grid.n, default, dtype=np.float64),
        ),
        set_g=lambda *args, **kwargs: None,
    )

    state = module._inland_geomorphology_state(
        world,
        np.full(grid.n, 720.0, dtype=np.float64),
        0.0,
        crust_type,
        crust_domain,
        np.zeros(grid.n, dtype=np.float64),
        crust_stability,
        np.full(grid.n, -1.0, dtype=np.float64),
    )

    rift = np.asarray(state["rift_axis"], dtype=bool) | np.asarray(
        state["rift_shoulder"], dtype=bool)
    assert float(np.count_nonzero(rift) / grid.n) < 0.20


def test_p172_age_aware_lifecycle_response_uses_object_kind_and_age():
    grid = SphereGrid.fibonacci(180, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    metrics = {}
    world = SimpleNamespace(
        grid=grid,
        time_myr=2400.0,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    active_cells = np.arange(20, 36, dtype=np.int64)
    basin_cells = np.arange(70, 86, dtype=np.int64)
    shield_cells = np.arange(120, 136, dtype=np.int64)
    candidate = np.ones(grid.n, dtype=bool)
    province_code = np.zeros(grid.n, dtype=np.float64)
    sediment_signal = np.zeros(grid.n, dtype=np.float64)
    sediment_signal[basin_cells] = 1.0
    thick_signal = np.zeros(grid.n, dtype=np.float64)
    thick_signal[active_cells] = 1.0
    stability = np.full(grid.n, 0.5, dtype=np.float64)
    stability[shield_cells] = 0.9

    response = module._p172_age_aware_inland_lifecycle_relief_response(
        world,
        candidate,
        province_code,
        {"province_id": np.zeros(grid.n, dtype=np.float64), "objects": []},
        [
            {
                "kind": "active_orogen",
                "age_myr": 45.0,
                "relief_stage": "young_orogen",
                "cells": active_cells.tolist(),
            },
            {
                "kind": "intracratonic_basin",
                "age_myr": 900.0,
                "relief_stage": "basin_fill",
                "cells": basin_cells.tolist(),
            },
            {
                "kind": "shield",
                "age_myr": 1900.0,
                "relief_stage": "shield",
                "cells": shield_cells.tolist(),
            },
        ],
        sediment_signal,
        thick_signal,
        stability,
    )

    assert float(response[active_cells].mean()) > 300.0
    assert float(response[basin_cells].mean()) < -170.0
    assert float(response[shield_cells].mean()) > 100.0
    assert metrics["terrain.last_p172_age_aware_inland_response_object_count"] == 3.0
    assert metrics["terrain.last_p172_age_aware_inland_response_area_fraction"] > 0.0


def test_p104f_applies_p172_age_aware_object_response_upstream():
    grid = SphereGrid.fibonacci(240, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    active_cells = np.arange(30, 52, dtype=np.int64)
    basin_cells = np.arange(120, 142, dtype=np.int64)
    metrics = {}
    fields = {
        "crust.thickness_m": np.full(n, 35000.0, dtype=np.float64),
    }
    fields["crust.thickness_m"][active_cells] = 52000.0
    world = SimpleNamespace(
        grid=grid,
        time_myr=2600.0,
        get_field=lambda name, default=0.0: fields.get(
            name,
            np.full(n, default, dtype=np.float64),
        ),
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    surface = np.full(n, 820.0, dtype=np.float64)
    crust_type = np.ones(n, dtype=np.float64)
    crust_domain = np.full(n, DOMAIN_CONTINENTAL_INTERIOR, dtype=np.float64)
    crust_origin = np.zeros(n, dtype=np.float64)
    crust_stability = np.full(n, 0.62, dtype=np.float64)
    sediment = np.zeros(n, dtype=np.float64)
    sediment[basin_cells] = 1400.0
    province_code = np.full(n, 2.0, dtype=np.float64)
    province_code[active_cells] = float(CONT_PROVINCE_ACTIVE_OROGEN)
    province_code[basin_cells] = float(CONT_PROVINCE_INTRACRATONIC_BASIN)
    province_id = np.zeros(n, dtype=np.float64)
    province_id[active_cells] = 1.0
    province_id[basin_cells] = 2.0

    def fake_inland_state(*args, **kwargs):
        empty = np.zeros(n, dtype=bool)
        return {
            "province": np.full(n, INLAND_PROVINCE_PLATFORM, dtype=np.float64),
            "regional_relief_m": np.zeros(n, dtype=np.float64),
            "shield": empty,
            "platform": empty,
            "sag_basin": empty,
            "old_orogen_root": empty,
            "rift_axis": empty,
            "rift_shoulder": empty,
            "platform_swell": empty,
            "escarpment": empty,
            "plateau_margin": empty,
        }

    module._inland_geomorphology_state = fake_inland_state
    result = module._apply_inland_landform_region_elevation_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_origin,
        crust_stability,
        sediment,
        np.full(n, -1.0, dtype=np.float64),
        np.full(n, CONT_DETAIL_PLATFORM, dtype=np.float64),
        np.full(n, CONT_DETAIL_PLATFORM, dtype=np.float64),
        {
            "province_id": province_id,
            "province_code": province_code,
            "objects": [
                {
                    "kind": "active_orogen",
                    "province_id": 1,
                    "age_myr": 60.0,
                    "relief_stage": "young_orogen",
                },
                {
                    "kind": "intracratonic_basin",
                    "province_id": 2,
                    "age_myr": 1100.0,
                    "relief_stage": "basin_fill",
                },
            ],
        },
        [],
    )

    out = result["surface"]
    assert float(out[active_cells].mean()) > float(surface[active_cells].mean())
    assert float(out[basin_cells].mean()) < float(surface[basin_cells].mean())
    assert metrics["terrain.last_p172_age_aware_inland_response_object_count"] == 2.0
    assert metrics["terrain.last_p172_age_aware_inland_response_area_fraction"] > 0.0


def test_p174_lowland_plain_response_runs_inside_p104f_for_parented_basins():
    grid = SphereGrid.fibonacci(360, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    plain = (
        (grid.lon > -135.0)
        & (grid.lon < -20.0)
        & (grid.lat > -20.0)
        & (grid.lat < 42.0)
    )
    assert np.count_nonzero(plain) >= 8
    metrics = {}
    fields = {
        "crust.thickness_m": np.full(n, 35000.0, dtype=np.float64),
        "tectonics.platform_subsidence": np.zeros(n, dtype=np.float64),
    }
    fields["tectonics.platform_subsidence"][plain] = 0.55
    world = SimpleNamespace(
        grid=grid,
        time_myr=2800.0,
        get_field=lambda name, default=0.0: fields.get(
            name,
            np.full(n, default, dtype=np.float64),
        ),
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    surface = np.full(n, 1280.0, dtype=np.float64)
    surface[plain] = 1120.0
    crust_type = np.ones(n, dtype=np.float64)
    crust_domain = np.full(n, DOMAIN_CONTINENTAL_INTERIOR, dtype=np.float64)
    crust_origin = np.zeros(n, dtype=np.float64)
    crust_stability = np.full(n, 0.58, dtype=np.float64)
    sediment = np.full(n, 90.0, dtype=np.float64)
    sediment[plain] = 1450.0
    detail = np.full(n, CONT_DETAIL_PLATFORM, dtype=np.float64)
    detail[plain] = CONT_DETAIL_BASIN
    inland = np.full(n, INLAND_PROVINCE_PLATFORM, dtype=np.float64)
    inland[plain] = INLAND_PROVINCE_SAG_BASIN
    province_code = np.full(n, CONT_PROVINCE_PLATFORM, dtype=np.float64)
    province_code[plain] = CONT_PROVINCE_INTRACRATONIC_BASIN

    def fake_inland_state(*args, **kwargs):
        empty = np.zeros(n, dtype=bool)
        return {
            "province": inland.copy(),
            "regional_relief_m": np.zeros(n, dtype=np.float64),
            "shield": empty,
            "platform": empty,
            "sag_basin": plain.copy(),
            "old_orogen_root": empty,
            "rift_axis": empty,
            "rift_shoulder": empty,
            "platform_swell": empty,
            "escarpment": empty,
            "plateau_margin": empty,
        }

    module._inland_geomorphology_state = fake_inland_state
    result = module._apply_inland_landform_region_elevation_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_origin,
        crust_stability,
        sediment,
        np.full(n, -1.0, dtype=np.float64),
        detail,
        detail,
        {
            "province_id": np.zeros(n, dtype=np.float64),
            "province_code": province_code,
            "objects": [],
        },
        [
            {
                "kind": "foreland_basin",
                "age_myr": 420.0,
                "relief_stage": "basin_fill",
                "cells": np.where(plain)[0].astype(int).tolist(),
            },
        ],
    )

    out = result["surface"]
    assert float(out[plain].mean()) < 690.0
    assert float(out[~plain].mean()) > 930.0
    assert float(out[~plain].mean()) > float(out[plain].mean()) + 250.0
    assert metrics["terrain.last_p174_lowland_plain_response_area_fraction"] > 0.0
    assert metrics["terrain.last_p174_lowland_plain_fraction_after"] > metrics[
        "terrain.last_p174_lowland_plain_fraction_before"]
    assert metrics["terrain.last_p174_lowland_plain_parented_fraction_after"] > 0.90
    assert metrics["terrain.last_p174_lowland_plain_response_land_mask_preserved"] == 1.0


def test_p174_lowland_plain_response_requires_process_parentage():
    grid = SphereGrid.fibonacci(260, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    candidate = (
        (grid.lon > 30.0)
        & (grid.lon < 150.0)
        & (grid.lat > -25.0)
        & (grid.lat < 45.0)
    )
    assert np.count_nonzero(candidate) >= 6
    metrics = {}
    world = SimpleNamespace(
        grid=grid,
        time_myr=2200.0,
        get_field=lambda name, default=0.0: np.full(n, default, dtype=np.float64),
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    surface = np.full(n, 1120.0, dtype=np.float64)
    crust_type = np.ones(n, dtype=np.float64)
    crust_domain = np.full(n, DOMAIN_CONTINENTAL_INTERIOR, dtype=np.float64)
    crust_stability = np.full(n, 0.55, dtype=np.float64)
    sediment = np.full(n, 60.0, dtype=np.float64)
    detail = np.zeros(n, dtype=np.float64)
    province = np.zeros(n, dtype=np.float64)
    inland = np.zeros(n, dtype=np.float64)

    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        np.full(n, -1.0, dtype=np.float64),
        detail,
        province,
        inland,
        [],
    )

    assert np.allclose(out, surface)
    assert metrics["terrain.last_p174_lowland_plain_response_area_fraction"] == 0.0
    assert metrics["terrain.last_p174_lowland_plain_fraction_after"] == metrics[
        "terrain.last_p174_lowland_plain_fraction_before"]


def test_p174_lowland_response_is_p170_broad_parented_plain():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -170.0)
        & (grid.lon < 50.0)
        & (grid.lat > -50.0)
        & (grid.lat < 60.0)
    )
    plain = (
        land
        & (grid.lon > -100.0)
        & (grid.lon < -40.0)
        & (grid.lat > -10.0)
        & (grid.lat < 15.0)
    )
    assert np.count_nonzero(plain) >= 20
    metrics = {}
    fields = {
        "tectonics.platform_subsidence": np.where(plain, 0.60, 0.02).astype(
            np.float64
        ),
    }
    world = SimpleNamespace(
        grid=grid,
        time_myr=4500.0,
        get_field=lambda name, default=0.0: fields.get(
            name, np.full(n, default, dtype=np.float64)
        ),
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    surface = np.where(land, 1700.0, -4200.0).astype(np.float64)
    surface[plain] = 1200.0
    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.45, 0.0).astype(np.float64)
    sediment = np.where(land, 100.0, 40.0).astype(np.float64)
    sediment[plain] = 1600.0
    orog_age = np.full(n, -1.0, dtype=np.float64)
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[plain] = CONT_DETAIL_BASIN
    province = np.zeros(n, dtype=np.float64)
    province[land] = CONT_PROVINCE_PLATFORM
    province[plain] = CONT_PROVINCE_INTRACRATONIC_BASIN
    inland = np.zeros(n, dtype=np.float64)
    inland[land] = INLAND_PROVINCE_PLATFORM
    inland[plain] = INLAND_PROVINCE_SAG_BASIN

    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        orog_age,
        detail,
        province,
        inland,
        [],
        stage_code=2.0,
    )
    frame = Frame(
        time_myr=4500.0,
        globals={"ocean.sea_level_m": 0.0},
        fields={
            "terrain.elevation_m": out,
            "crust.type": crust_type,
            "crust.age_myr": np.where(land, 1800.0, 90.0),
            "crust.domain": crust_domain,
            "crust.origin": np.zeros(n, dtype=np.float64),
            "terrain.continental_detail": detail,
            "terrain.continental_detail_region_code": detail,
            "terrain.inland_geomorphology_region_code": inland,
            "terrain.continental_province_code": province,
            "tectonics.orogeny_age_myr": orog_age,
            "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
            "sediment.thickness_m": sediment,
            "ocean.depth_province": np.where(
                land, 0.0, OCEAN_DEPTH_ABYSS
            ).astype(np.float64),
        },
    )
    row = frame_geomorphology_metrics(grid, frame)
    land_metrics = row["land_metrics"]

    assert metrics["terrain.last_p174_lowland_plain_response_area_fraction"] > 0.0
    assert land_metrics["lowland_plain_fraction"] > 0.08
    assert land_metrics["broad_lowland_plain_component_count"] >= 1
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.025
    assert land_metrics["lowland_plain_parented_fraction"] > 0.80
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]


def test_p174_lowland_continuity_memory_preserves_supported_platform_plain():
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -160.0)
        & (grid.lon < 80.0)
        & (grid.lat > -55.0)
        & (grid.lat < 60.0)
    )
    basin = (
        land
        & (grid.lon > -95.0)
        & (grid.lon < -20.0)
        & (grid.lat > -18.0)
        & (grid.lat < 25.0)
    )
    assert np.count_nonzero(basin) >= 12
    metrics = {}
    fields = {
        "tectonics.platform_subsidence": np.where(basin, 0.55, 0.02).astype(
            np.float64
        ),
    }

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=2200.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    surface = np.where(land, 1600.0, -4200.0).astype(np.float64)
    surface[basin] = 1180.0
    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.48, 0.0).astype(np.float64)
    sediment = np.where(land, 80.0, 40.0).astype(np.float64)
    sediment[basin] = 1450.0
    orog_age = np.full(n, -1.0, dtype=np.float64)
    basin_detail = np.zeros(n, dtype=np.float64)
    basin_detail[land] = CONT_DETAIL_PLATFORM
    basin_detail[basin] = CONT_DETAIL_BASIN
    basin_province = np.zeros(n, dtype=np.float64)
    basin_province[land] = CONT_PROVINCE_PLATFORM
    basin_province[basin] = CONT_PROVINCE_INTRACRATONIC_BASIN
    basin_inland = np.zeros(n, dtype=np.float64)
    basin_inland[land] = INLAND_PROVINCE_PLATFORM
    basin_inland[basin] = INLAND_PROVINCE_SAG_BASIN

    module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        orog_age,
        basin_detail,
        basin_province,
        basin_inland,
        [],
        stage_code=1.0,
    )
    memory = fields["terrain.p174_lowland_plain_continuity_memory"]
    assert float(memory[basin].mean()) > 0.50

    metrics.clear()
    fields["tectonics.platform_subsidence"] = np.zeros(n, dtype=np.float64)
    surface_reworked = np.where(land, 1450.0, -4200.0).astype(np.float64)
    surface_reworked[basin] = 1320.0
    platform_detail = np.zeros(n, dtype=np.float64)
    platform_detail[land] = CONT_DETAIL_PLATFORM
    platform_province = np.zeros(n, dtype=np.float64)
    platform_province[land] = CONT_PROVINCE_PLATFORM
    platform_inland = np.zeros(n, dtype=np.float64)
    platform_inland[land] = INLAND_PROVINCE_PLATFORM
    weak_sediment = np.where(land, 90.0, 40.0).astype(np.float64)

    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface_reworked,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        weak_sediment,
        orog_age,
        platform_detail,
        platform_province,
        platform_inland,
        [],
        stage_code=1.0,
    )

    assert metrics[
        "terrain.last_p174_lowland_plain_continuity_parent_area_fraction"] > 0.0
    assert metrics[
        "terrain.last_p174_lowland_plain_continuity_memory_area_fraction"] > 0.0
    assert float(out[basin].mean()) < float(surface_reworked[basin].mean()) - 250.0


def test_p174_lowland_continuity_microbenchmark_bounds_component_drift():
    grid = SphereGrid.fibonacci(1600, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -170.0)
        & (grid.lon < 80.0)
        & (grid.lat > -58.0)
        & (grid.lat < 62.0)
    )
    basin = (
        land
        & (grid.lon > -120.0)
        & (grid.lon < -15.0)
        & (grid.lat > -28.0)
        & (grid.lat < 35.0)
    )
    assert np.count_nonzero(basin) >= 120
    metrics = {}
    fields = {
        "tectonics.platform_subsidence": np.where(basin, 0.55, 0.03).astype(
            np.float64
        ),
    }

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=0.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.50, 0.0).astype(np.float64)
    orog_age = np.full(n, -1.0, dtype=np.float64)
    rows = []

    for step, time_myr in enumerate((1600.0, 2000.0, 2400.0, 2800.0)):
        world.time_myr = time_myr
        if step == 0:
            surface = np.where(land, 1700.0, -4300.0).astype(np.float64)
            surface[basin] = 1250.0
            sediment = np.where(land, 90.0, 40.0).astype(np.float64)
            sediment[basin] = 1550.0
            detail = np.zeros(n, dtype=np.float64)
            detail[land] = CONT_DETAIL_PLATFORM
            detail[basin] = CONT_DETAIL_BASIN
            province = np.zeros(n, dtype=np.float64)
            province[land] = CONT_PROVINCE_PLATFORM
            province[basin] = CONT_PROVINCE_INTRACRATONIC_BASIN
            inland = np.zeros(n, dtype=np.float64)
            inland[land] = INLAND_PROVINCE_PLATFORM
            inland[basin] = INLAND_PROVINCE_SAG_BASIN
            fields["tectonics.platform_subsidence"] = np.where(
                basin, 0.55, 0.03
            ).astype(np.float64)
        else:
            surface = np.where(land, 1550.0, -4300.0).astype(np.float64)
            surface[basin] = 1300.0 + 80.0 * step
            ribs = basin & (
                ((np.floor((grid.lon + 180.0) / 10.0) + step) % 4) == 0
            )
            surface[ribs] = 1120.0 + 70.0 * step
            sediment = np.where(land, 100.0, 40.0).astype(np.float64)
            sediment[basin] = 220.0
            detail = np.zeros(n, dtype=np.float64)
            detail[land] = CONT_DETAIL_PLATFORM
            province = np.zeros(n, dtype=np.float64)
            province[land] = CONT_PROVINCE_PLATFORM
            inland = np.zeros(n, dtype=np.float64)
            inland[land] = INLAND_PROVINCE_PLATFORM
            fields["tectonics.platform_subsidence"] = np.zeros(
                n, dtype=np.float64
            )

        metrics.clear()
        out = module._p174_process_parented_lowland_plain_response(
            world,
            surface,
            0.0,
            crust_type,
            crust_domain,
            crust_stability,
            sediment,
            orog_age,
            detail,
            province,
            inland,
            [],
            stage_code=1.0,
        )
        frame = Frame(
            time_myr=time_myr,
            globals={"ocean.sea_level_m": 0.0},
            fields={
                "terrain.elevation_m": out,
                "crust.type": crust_type,
                "crust.age_myr": np.where(land, 1800.0, 100.0),
                "crust.domain": crust_domain,
                "crust.origin": np.zeros(n, dtype=np.float64),
                "terrain.continental_detail": detail,
                "terrain.continental_detail_region_code": detail,
                "terrain.inland_geomorphology_region_code": inland,
                "terrain.continental_province_code": province,
                "tectonics.orogeny_age_myr": orog_age,
                "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
                "sediment.thickness_m": sediment,
                "terrain.p174_lowland_plain_continuity_memory": fields.get(
                    "terrain.p174_lowland_plain_continuity_memory",
                    np.zeros(n, dtype=np.float64),
                ),
                "ocean.depth_province": np.where(
                    land, 0.0, OCEAN_DEPTH_ABYSS
                ).astype(np.float64),
            },
        )
        rows.append(frame_geomorphology_metrics(grid, frame))

    lowland = np.asarray([
        row["land_metrics"]["lowland_plain_fraction"] for row in rows
    ], dtype=np.float64)
    largest = np.asarray([
        row["land_metrics"]["largest_lowland_plain_component_fraction"]
        for row in rows
    ], dtype=np.float64)

    assert float(np.min(lowland)) > 0.12
    assert float(np.max(lowland) - np.min(lowland)) < 0.05
    assert float(np.min(largest)) > 0.10
    for row in rows:
        land_metrics = row["land_metrics"]
        assert land_metrics["broad_lowland_plain_component_count"] >= 1
        assert land_metrics["lowland_plain_parented_fraction"] > 0.95
        assert not row["diagnostic_flags"]["lowland_plain_deficient"]


def test_p174_lowland_fragment_bridge_microbenchmark_connects_adjacent_components():
    grid = SphereGrid.fibonacci(2500, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -175.0)
        & (grid.lon < 110.0)
        & (grid.lat > -62.0)
        & (grid.lat < 66.0)
    )
    belt = (
        land
        & (grid.lon > -145.0)
        & (grid.lon < 40.0)
        & (grid.lat > -28.0)
        & (grid.lat < 34.0)
    )
    low_fragments = np.zeros(n, dtype=bool)
    for start_lon in (-142.0, -110.0, -78.0, -46.0, -14.0, 18.0):
        low_fragments |= (
            belt
            & (grid.lon > start_lon)
            & (grid.lon < start_lon + 10.0)
        )
    high_separators = belt & ~low_fragments
    assert np.count_nonzero(low_fragments) >= 150
    assert np.count_nonzero(high_separators) >= 250

    fields = {
        "terrain.p174_lowland_plain_continuity_memory": low_fragments.astype(
            np.float64
        ),
        "tectonics.platform_subsidence": np.zeros(n, dtype=np.float64),
    }
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=2600.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    surface = np.where(land, 1500.0, -4300.0).astype(np.float64)
    surface[low_fragments] = 360.0
    surface[high_separators] = 2800.0
    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.55, 0.0).astype(np.float64)
    sediment = np.where(land, 100.0, 40.0).astype(np.float64)
    sediment[belt] = 180.0
    orog_age = np.full(n, -1.0, dtype=np.float64)
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    province = np.zeros(n, dtype=np.float64)
    province[land] = CONT_PROVINCE_PLATFORM
    inland = np.zeros(n, dtype=np.float64)
    inland[land] = INLAND_PROVINCE_PLATFORM

    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        orog_age,
        detail,
        province,
        inland,
        [],
        stage_code=1.0,
    )
    frame = Frame(
        time_myr=2600.0,
        globals={"ocean.sea_level_m": 0.0},
        fields={
            "terrain.elevation_m": out,
            "crust.type": crust_type,
            "crust.age_myr": np.where(land, 1800.0, 100.0),
            "crust.domain": crust_domain,
            "crust.origin": np.zeros(n, dtype=np.float64),
            "terrain.continental_detail": detail,
            "terrain.continental_detail_region_code": detail,
            "terrain.inland_geomorphology_region_code": inland,
            "terrain.continental_province_code": province,
            "tectonics.orogeny_age_myr": orog_age,
            "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
            "sediment.thickness_m": sediment,
            "terrain.p174_lowland_plain_continuity_memory": fields[
                "terrain.p174_lowland_plain_continuity_memory"
            ],
            "ocean.depth_province": np.where(
                land, 0.0, OCEAN_DEPTH_ABYSS
            ).astype(np.float64),
        },
    )
    row = frame_geomorphology_metrics(grid, frame)
    land_metrics = row["land_metrics"]

    assert metrics["terrain.last_p174_lowland_plain_candidate_area_fraction"] > 0.10
    assert metrics["terrain.last_p174_lowland_plain_response_area_fraction"] > 0.10
    assert land_metrics["lowland_plain_fraction"] > 0.09
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.08
    assert land_metrics["broad_lowland_plain_component_count"] >= 1
    assert land_metrics["lowland_plain_parented_fraction"] > 0.95
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]


def test_p174_stable_platform_support_microbenchmark_prevents_terminal_only_lowland():
    grid = SphereGrid.fibonacci(2200, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -170.0)
        & (grid.lon < 95.0)
        & (grid.lat > -60.0)
        & (grid.lat < 64.0)
    )
    plain = (
        land
        & (grid.lon > -125.0)
        & (grid.lon < 25.0)
        & (grid.lat > -34.0)
        & (grid.lat < 36.0)
    )
    assert np.count_nonzero(plain) >= 350

    metrics = {}
    fields = {
        "tectonics.platform_subsidence": np.zeros(n, dtype=np.float64),
    }

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=0.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.76, 0.0).astype(np.float64)
    orog_age = np.full(n, -1.0, dtype=np.float64)
    rows = []

    phases = (
        {
            "time_myr": 1800.0,
            "relief_m": 2150.0,
            "rib_relief_m": 2380.0,
            "sediment_m": 100.0,
            "subsidence": 0.0,
            "detail": CONT_DETAIL_PLATFORM,
            "province": CONT_PROVINCE_PLATFORM,
            "inland": INLAND_PROVINCE_PLATFORM,
        },
        {
            "time_myr": 2200.0,
            "relief_m": 2850.0,
            "rib_relief_m": 3120.0,
            "sediment_m": 1450.0,
            "subsidence": 0.52,
            "detail": CONT_DETAIL_BASIN,
            "province": CONT_PROVINCE_INTRACRATONIC_BASIN,
            "inland": INLAND_PROVINCE_SAG_BASIN,
        },
        {
            "time_myr": 2600.0,
            "relief_m": 1520.0,
            "rib_relief_m": 1740.0,
            "sediment_m": 130.0,
            "subsidence": 0.0,
            "detail": CONT_DETAIL_PLATFORM,
            "province": CONT_PROVINCE_PLATFORM,
            "inland": INLAND_PROVINCE_PLATFORM,
        },
        {
            "time_myr": 3000.0,
            "relief_m": 1380.0,
            "rib_relief_m": 1620.0,
            "sediment_m": 300.0,
            "subsidence": 0.10,
            "detail": CONT_DETAIL_PLATFORM,
            "province": CONT_PROVINCE_PLATFORM,
            "inland": INLAND_PROVINCE_PLATFORM,
        },
    )

    for index, phase in enumerate(phases):
        world.time_myr = phase["time_myr"]
        surface = np.where(land, 1850.0, -4300.0).astype(np.float64)
        surface[plain] = phase["relief_m"]
        ribs = plain & (
            (
                np.floor((grid.lon + 180.0) / 12.0)
                + np.floor((grid.lat + 90.0) / 18.0)
                + index
            )
            % 5
            == 0
        )
        surface[ribs] = phase["rib_relief_m"]
        sediment = np.where(land, 100.0, 40.0).astype(np.float64)
        sediment[plain] = phase["sediment_m"]
        detail = np.zeros(n, dtype=np.float64)
        detail[land] = CONT_DETAIL_PLATFORM
        detail[plain] = phase["detail"]
        province = np.zeros(n, dtype=np.float64)
        province[land] = CONT_PROVINCE_PLATFORM
        province[plain] = phase["province"]
        inland = np.zeros(n, dtype=np.float64)
        inland[land] = INLAND_PROVINCE_PLATFORM
        inland[plain] = phase["inland"]
        fields["tectonics.platform_subsidence"] = np.zeros(n, dtype=np.float64)
        fields["tectonics.platform_subsidence"][plain] = phase["subsidence"]

        metrics.clear()
        out = module._p174_process_parented_lowland_plain_response(
            world,
            surface,
            0.0,
            crust_type,
            crust_domain,
            crust_stability,
            sediment,
            orog_age,
            detail,
            province,
            inland,
            [],
            stage_code=1.0,
        )
        frame = Frame(
            time_myr=phase["time_myr"],
            globals={"ocean.sea_level_m": 0.0},
            fields={
                "terrain.elevation_m": out,
                "crust.type": crust_type,
                "crust.age_myr": np.where(land, 1800.0, 100.0),
                "crust.domain": crust_domain,
                "crust.origin": np.zeros(n, dtype=np.float64),
                "terrain.continental_detail": detail,
                "terrain.continental_detail_region_code": detail,
                "terrain.inland_geomorphology_region_code": inland,
                "terrain.continental_province_code": province,
                "tectonics.orogeny_age_myr": orog_age,
                "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
                "sediment.thickness_m": sediment,
                "terrain.p174_lowland_plain_continuity_memory": fields.get(
                    "terrain.p174_lowland_plain_continuity_memory",
                    np.zeros(n, dtype=np.float64),
                ),
                "ocean.depth_province": np.where(
                    land, 0.0, OCEAN_DEPTH_ABYSS
                ).astype(np.float64),
            },
        )
        rows.append((metrics.copy(), frame_geomorphology_metrics(grid, frame)))

    lowland = np.asarray([
        row["land_metrics"]["lowland_plain_fraction"] for _, row in rows
    ], dtype=np.float64)
    largest = np.asarray([
        row["land_metrics"]["largest_lowland_plain_component_fraction"]
        for _, row in rows
    ], dtype=np.float64)

    assert float(np.min(lowland)) > 0.10
    assert float(np.max(lowland) - np.min(lowland)) < 0.09
    assert float(np.min(largest)) > 0.085
    for phase_metrics, row in rows:
        land_metrics = row["land_metrics"]
        assert phase_metrics[
            "terrain.last_p174_lowland_plain_candidate_area_fraction"
        ] > 0.06
        assert phase_metrics[
            "terrain.last_p174_lowland_plain_continuity_memory_area_fraction"
        ] > 0.06
        assert land_metrics["broad_lowland_plain_component_count"] >= 1
        assert land_metrics["lowland_plain_parented_fraction"] > 0.95
        assert not row["diagnostic_flags"]["lowland_plain_deficient"]


def test_p174_parented_lowland_relief_microbenchmark_qualifies_rough_low_plain():
    grid = SphereGrid.fibonacci(2600, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -170.0)
        & (grid.lon < 120.0)
        & (grid.lat > -58.0)
        & (grid.lat < 62.0)
    )
    plain = (
        land
        & (grid.lon > -128.0)
        & (grid.lon < 28.0)
        & (grid.lat > -36.0)
        & (grid.lat < 38.0)
    )
    far_platform = land & ~module._dilate_mask(grid, plain, passes=2)
    assert np.count_nonzero(plain) >= 450
    assert np.count_nonzero(far_platform) >= 300

    surface = np.where(land, 1420.0, -4300.0).astype(np.float64)
    surface[plain] = 190.0
    ribs = plain & (
        (
            np.floor((grid.lon + 180.0) / 13.0)
            + np.floor((grid.lat + 90.0) / 17.0)
        )
        % 3
        == 0
    )
    surface[ribs] = 690.0
    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.66, 0.0).astype(np.float64)
    sediment = np.where(land, 120.0, 35.0).astype(np.float64)
    sediment[plain] = 1650.0
    orog_age = np.full(n, -1.0, dtype=np.float64)
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[plain] = CONT_DETAIL_BASIN
    province = np.zeros(n, dtype=np.float64)
    province[land] = CONT_PROVINCE_PLATFORM
    province[plain] = CONT_PROVINCE_INTRACRATONIC_BASIN
    inland = np.zeros(n, dtype=np.float64)
    inland[land] = INLAND_PROVINCE_PLATFORM
    inland[plain] = INLAND_PROVINCE_SAG_BASIN
    fields = {
        "tectonics.platform_subsidence": np.zeros(n, dtype=np.float64),
    }
    fields["tectonics.platform_subsidence"][plain] = 0.54
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=2600.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    frame_fields = {
        "crust.type": crust_type,
        "crust.age_myr": np.where(land, 2100.0, 120.0),
        "crust.domain": crust_domain,
        "crust.origin": np.zeros(n, dtype=np.float64),
        "terrain.continental_detail": detail,
        "terrain.continental_detail_region_code": detail,
        "terrain.inland_geomorphology_region_code": inland,
        "terrain.continental_province_code": province,
        "tectonics.orogeny_age_myr": orog_age,
        "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
        "sediment.thickness_m": sediment,
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS
        ).astype(np.float64),
    }
    before_row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=2600.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={"terrain.elevation_m": surface, **frame_fields},
        ),
    )

    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        orog_age,
        detail,
        province,
        inland,
        [],
        stage_code=1.0,
    )
    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=2600.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={
                "terrain.elevation_m": out,
                **frame_fields,
                "terrain.p174_lowland_plain_continuity_memory": fields.get(
                    "terrain.p174_lowland_plain_continuity_memory",
                    np.zeros(n, dtype=np.float64),
                ),
            },
        ),
    )
    before_metrics = before_row["land_metrics"]
    land_metrics = row["land_metrics"]

    assert before_metrics["lowland_local_relief_blocked_fraction"] > 0.20
    assert land_metrics["lowland_elevation_parented_fraction"] > 0.16
    assert land_metrics["lowland_local_relief_blocked_fraction"] < (
        0.45 * before_metrics["lowland_local_relief_blocked_fraction"]
    )
    assert land_metrics["lowland_plain_fraction"] > 0.12
    assert land_metrics["lowland_plain_fraction"] > (
        before_metrics["lowland_plain_fraction"] + 0.18
    )
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.08
    assert land_metrics["broad_lowland_plain_component_count"] >= 1
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]
    assert float(out[far_platform].mean()) > 900.0
    assert metrics["terrain.last_p174_lowland_plain_response_area_fraction"] > 0.08


def test_p174_lowland_component_stitching_expands_parented_elevation_component():
    grid = SphereGrid.fibonacci(2800, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -172.0)
        & (grid.lon < 118.0)
        & (grid.lat > -58.0)
        & (grid.lat < 62.0)
    )
    plain = (
        land
        & (grid.lon > -132.0)
        & (grid.lon < 38.0)
        & (grid.lat > -36.0)
        & (grid.lat < 38.0)
    )
    far_platform = land & ~module._dilate_mask(grid, plain, passes=2)
    assert np.count_nonzero(plain) >= 500
    assert np.count_nonzero(far_platform) >= 300

    surface = np.where(land, 1460.0, -4300.0).astype(np.float64)
    surface[plain] = 180.0
    internal_ribs = plain & (
        (
            np.floor((grid.lon + 180.0) / 11.0)
            + 2.0 * np.floor((grid.lat + 90.0) / 13.0)
        )
        % 4
        == 0
    )
    surface[internal_ribs] = 680.0
    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.64, 0.0).astype(np.float64)
    sediment = np.where(land, 120.0, 35.0).astype(np.float64)
    sediment[plain] = 1550.0
    orog_age = np.full(n, -1.0, dtype=np.float64)
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[plain] = CONT_DETAIL_BASIN
    province = np.zeros(n, dtype=np.float64)
    province[land] = CONT_PROVINCE_PLATFORM
    province[plain] = CONT_PROVINCE_INTRACRATONIC_BASIN
    inland = np.zeros(n, dtype=np.float64)
    inland[land] = INLAND_PROVINCE_PLATFORM
    inland[plain] = INLAND_PROVINCE_SAG_BASIN
    fields = {
        "tectonics.platform_subsidence": np.zeros(n, dtype=np.float64),
    }
    fields["tectonics.platform_subsidence"][plain] = 0.50
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=1900.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    frame_fields = {
        "crust.type": crust_type,
        "crust.age_myr": np.where(land, 1800.0, 120.0),
        "crust.domain": crust_domain,
        "crust.origin": np.zeros(n, dtype=np.float64),
        "terrain.continental_detail": detail,
        "terrain.continental_detail_region_code": detail,
        "terrain.inland_geomorphology_region_code": inland,
        "terrain.continental_province_code": province,
        "tectonics.orogeny_age_myr": orog_age,
        "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
        "sediment.thickness_m": sediment,
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS
        ).astype(np.float64),
    }
    before_row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1900.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={"terrain.elevation_m": surface, **frame_fields},
        ),
    )

    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        orog_age,
        detail,
        province,
        inland,
        [],
        stage_code=1.0,
    )
    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1900.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={
                "terrain.elevation_m": out,
                **frame_fields,
                "terrain.p174_lowland_plain_continuity_memory": fields.get(
                    "terrain.p174_lowland_plain_continuity_memory",
                    np.zeros(n, dtype=np.float64),
                ),
            },
        ),
    )
    before_metrics = before_row["land_metrics"]
    land_metrics = row["land_metrics"]

    assert before_metrics[
        "largest_lowland_elevation_parented_component_fraction"
    ] > 0.10
    assert before_metrics[
        "largest_lowland_elevation_parented_component_fraction"
    ] > 1.5 * before_metrics["largest_lowland_plain_component_fraction"]
    assert land_metrics["lowland_plain_fraction"] > 0.14
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.09
    assert land_metrics["largest_lowland_plain_component_fraction"] > (
        before_metrics["largest_lowland_plain_component_fraction"] + 0.05
    )
    assert land_metrics["broad_lowland_plain_component_count"] >= 1
    assert land_metrics["lowland_plain_parented_fraction"] > 0.90
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]
    assert float(out[far_platform].mean()) > 900.0
    assert metrics["terrain.last_p174_lowland_plain_response_area_fraction"] > 0.08


def test_p174_component_stitching_uses_multiring_same_parent_corridors():
    grid = SphereGrid.fibonacci(3600, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -176.0)
        & (grid.lon < 176.0)
        & (grid.lat > -62.0)
        & (grid.lat < 64.0)
    )
    low_fragments = np.zeros(n, dtype=bool)
    continuity_corridor = np.zeros(n, dtype=bool)
    fragment_starts: list[float] = []
    lon = -158.0
    while lon < 72.0:
        fragment_starts.append(lon)
        lon += 28.0
    for lat in (-40.0, -16.0, 8.0, 32.0):
        for lon in fragment_starts:
            low_fragments |= (
                land
                & (grid.lon > lon)
                & (grid.lon < lon + 16.0)
                & (grid.lat > lat)
                & (grid.lat < lat + 9.0)
            )
        for left, right in zip(fragment_starts[:-1], fragment_starts[1:]):
            continuity_corridor |= (
                land
                & (grid.lon > left + 16.0)
                & (grid.lon < right)
                & (grid.lat > lat + 2.0)
                & (grid.lat < lat + 7.0)
            )
    far_platform = land & (grid.lon > 118.0) & (grid.lat < -24.0)
    assert np.count_nonzero(low_fragments) >= 300
    assert np.count_nonzero(continuity_corridor) >= 90
    assert np.count_nonzero(far_platform) >= 80

    surface = np.where(land, 2800.0, -4300.0).astype(np.float64)
    surface[low_fragments] = 180.0
    surface[continuity_corridor] = 1250.0
    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.18, 0.0).astype(np.float64)
    sediment = np.where(land, 80.0, 35.0).astype(np.float64)
    sediment[low_fragments] = 1200.0
    sediment[continuity_corridor] = 180.0
    orog_age = np.full(n, -1.0, dtype=np.float64)
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[continuity_corridor] = 0.0
    detail[low_fragments] = CONT_DETAIL_BASIN
    province = np.zeros(n, dtype=np.float64)
    province[land] = CONT_PROVINCE_PLATFORM
    province[continuity_corridor] = 0.0
    province[low_fragments] = CONT_PROVINCE_INTRACRATONIC_BASIN
    inland = np.zeros(n, dtype=np.float64)
    inland[land] = INLAND_PROVINCE_PLATFORM
    inland[continuity_corridor] = 0.0
    inland[low_fragments] = INLAND_PROVINCE_SAG_BASIN
    fields = {
        "tectonics.platform_subsidence": np.zeros(n, dtype=np.float64),
        "terrain.p174_lowland_plain_continuity_memory": np.where(
            low_fragments,
            0.45,
            np.where(continuity_corridor, 0.12, 0.0),
        ).astype(np.float64),
    }
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=1500.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    frame_fields = {
        "crust.type": crust_type,
        "crust.age_myr": np.where(land, 1600.0, 110.0),
        "crust.domain": crust_domain,
        "crust.origin": np.zeros(n, dtype=np.float64),
        "terrain.continental_detail": detail,
        "terrain.continental_detail_region_code": detail,
        "terrain.inland_geomorphology_region_code": inland,
        "terrain.continental_province_code": province,
        "tectonics.orogeny_age_myr": orog_age,
        "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
        "sediment.thickness_m": sediment,
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS
        ).astype(np.float64),
    }
    before = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1500.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={"terrain.elevation_m": surface, **frame_fields},
        ),
    )["land_metrics"]
    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        orog_age,
        detail,
        province,
        inland,
        [],
        stage_code=2.0,
    )
    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1500.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={
                "terrain.elevation_m": out,
                **frame_fields,
                "terrain.p174_lowland_plain_continuity_memory": fields[
                    "terrain.p174_lowland_plain_continuity_memory"
                ],
            },
        ),
    )
    land_metrics = row["land_metrics"]

    assert before["lowland_plain_fraction"] > 0.10
    assert before["largest_lowland_plain_component_fraction"] < 0.01
    assert metrics[
        "terrain.last_p174_component_stitch_multiring_response_area_fraction"
    ] > 0.006
    assert metrics[
        "terrain.last_p174_component_stitch_response_area_fraction"
    ] >= metrics[
        "terrain.last_p174_component_stitch_multiring_response_area_fraction"
    ]
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.025
    assert float(out[far_platform].mean()) > 2600.0


def test_p174_upstream_p104f_prepares_broad_lowland_province_before_p174(monkeypatch):
    grid = SphereGrid.fibonacci(2600, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -165.0)
        & (grid.lon < 95.0)
        & (grid.lat > -58.0)
        & (grid.lat < 62.0)
    )
    basin = (
        land
        & (grid.lon > -110.0)
        & (grid.lon < 5.0)
        & (grid.lat > -30.0)
        & (grid.lat < 34.0)
    )
    platform_lowland = (
        land
        & (grid.lon > -140.0)
        & (grid.lon < 45.0)
        & (grid.lat > -42.0)
        & (grid.lat < 45.0)
    )
    assert np.count_nonzero(basin) >= 350
    assert np.count_nonzero(platform_lowland) >= 800

    surface = np.where(land, 2450.0, -4300.0).astype(np.float64)
    surface[basin] = 2820.0
    ribs = platform_lowland & (
        (np.floor((grid.lon + 180.0) / 14.0) % 5) == 0
    )
    surface[ribs] += 260.0

    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_origin = np.zeros(n, dtype=np.float64)
    crust_stability = np.where(land, 0.62, 0.0).astype(np.float64)
    crust_thick = np.where(land, 35000.0, 7000.0).astype(np.float64)
    crust_thick[basin] = 30000.0
    sediment = np.where(land, 180.0, 40.0).astype(np.float64)
    sediment[basin] = 1550.0
    platform_subsidence = np.zeros(n, dtype=np.float64)
    platform_subsidence[basin] = 0.58
    fields = {
        "crust.thickness_m": crust_thick,
        "crust.age_myr": np.where(land, 1900.0, 120.0).astype(np.float64),
        "tectonics.rift_potential": np.zeros(n, dtype=np.float64),
        "tectonics.platform_subsidence": platform_subsidence,
    }
    metrics = {}
    world = SimpleNamespace(
        grid=grid,
        time_myr=2600.0,
        objects={},
        networks={},
        get_field=lambda name, default=0.0: fields.get(
            name,
            np.full(n, default, dtype=np.float64),
        ),
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )

    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[basin] = CONT_DETAIL_BASIN
    province_code = np.zeros(n, dtype=np.float64)
    province_code[land] = CONT_PROVINCE_PLATFORM
    province_code[basin] = CONT_PROVINCE_INTRACRATONIC_BASIN
    province_graph = {
        "province_id": np.where(land, 1.0, 0.0).astype(np.float64),
        "province_code": province_code,
        "objects": [],
    }

    captured = {}

    def fake_p174(
        world_arg,
        p104f_surface,
        *args,
        **kwargs,
    ):
        captured["surface_before_p174"] = np.asarray(
            p104f_surface, dtype=np.float64
        ).copy()
        return p104f_surface

    monkeypatch.setattr(
        module,
        "_p174_process_parented_lowland_plain_response",
        fake_p174,
    )
    result = module._apply_inland_landform_region_elevation_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_origin,
        crust_stability,
        sediment,
        np.full(n, -1.0, dtype=np.float64),
        detail,
        detail,
        province_graph,
        [],
    )

    pre_p174 = captured["surface_before_p174"]
    assert np.array_equal(result["surface"] >= 0.0, surface >= 0.0)
    assert float(pre_p174[basin].mean()) < 1700.0
    assert float(pre_p174[platform_lowland & ~basin].mean()) < 2050.0
    assert (
        metrics["terrain.last_p104f_pre_p174_lowland_prep_area_fraction"]
        > 0.08
    )
    assert (
        metrics["terrain.last_p104f_pre_p174_lowland_prep_mean_lowering_m"]
        > 350.0
    )


def test_p174_generated_province_support_infers_broad_subsiding_basin_before_p174(monkeypatch):
    grid = SphereGrid.fibonacci(3000, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -170.0)
        & (grid.lon < 120.0)
        & (grid.lat > -56.0)
        & (grid.lat < 62.0)
    )
    subsiding_basin = (
        land
        & (grid.lon > -118.0)
        & (grid.lon < 8.0)
        & (grid.lat > -32.0)
        & (grid.lat < 34.0)
    )
    old_orogen = (
        land
        & (grid.lon > 35.0)
        & (grid.lon < 92.0)
        & (grid.lat > -44.0)
        & (grid.lat < 46.0)
    )
    passive_margin = (
        land
        & (grid.lon < -145.0)
        & (grid.lat > -45.0)
        & (grid.lat < 48.0)
    )
    assert np.count_nonzero(subsiding_basin) >= 350
    assert np.count_nonzero(old_orogen) >= 180
    assert np.count_nonzero(passive_margin) >= 60

    surface = np.where(land, 2380.0, -4300.0).astype(np.float64)
    surface[subsiding_basin] = 2860.0
    ribs = subsiding_basin & (
        (np.floor((grid.lon + 180.0) / 16.0) % 4) == 0
    )
    surface[ribs] += 220.0
    surface[old_orogen] = 1780.0
    surface[passive_margin] = 560.0

    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_origin = np.zeros(n, dtype=np.float64)
    crust_stability = np.where(land, 0.68, 0.0).astype(np.float64)
    sediment = np.where(land, 160.0, 40.0).astype(np.float64)
    sediment[subsiding_basin] = 1750.0
    sediment[passive_margin] = 950.0
    crust_thick = np.where(land, 35000.0, 7000.0).astype(np.float64)
    crust_thick[subsiding_basin] = 29500.0
    platform_subsidence = np.zeros(n, dtype=np.float64)
    platform_subsidence[subsiding_basin] = 0.64
    internal_blocks = np.zeros(n, dtype=np.float64)
    internal_blocks[land] = INTERNAL_BLOCK_STABLE_PLATFORM
    continent_id = np.where(land, 1.0, -1.0).astype(np.float64)
    fields = {
        "crust.thickness_m": crust_thick,
        "crust.age_myr": np.where(land, 2100.0, 90.0).astype(np.float64),
        "crust.stability": crust_stability,
        "crust.type": crust_type,
        "sediment.thickness_m": sediment,
        "tectonics.platform_subsidence": platform_subsidence,
        "tectonics.rift_potential": np.zeros(n, dtype=np.float64),
        "tectonics.internal_geographic_block_code": internal_blocks,
        "tectonics.continent_id": continent_id,
        "tectonics.deformation_style": np.zeros(n, dtype=np.float64),
        "tectonics.deformation_intensity": np.zeros(n, dtype=np.float64),
    }
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    world = SimpleNamespace(
        grid=grid,
        time_myr=2600.0,
        objects={},
        networks={},
        get_field=get_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    continental_landforms = [
        {
            "id": "landform:old_subdued_orogen:test",
            "kind": "old_subdued_orogen",
            "cells": np.where(old_orogen)[0].tolist(),
            "priority": 8,
        },
        {
            "id": "landform:passive_margin_lowland:test",
            "kind": "passive_margin_lowland",
            "cells": np.where(passive_margin)[0].tolist(),
            "priority": 9,
        },
    ]
    province_graph = module._production_continental_province_graph(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_origin,
        detail,
        passive_margin,
        continental_landforms,
        [],
    )
    province_code = province_graph["province_code"].astype(int)
    cont_area = max(float(grid.cell_area[land].sum()), 1.0e-12)

    def cont_fraction(mask):
        return float(grid.cell_area[mask].sum() / cont_area)

    basin_support = land & (
        province_code == CONT_PROVINCE_INTRACRATONIC_BASIN
    )
    platform_support = land & (province_code == CONT_PROVINCE_PLATFORM)
    passive_support = land & (
        province_code == CONT_PROVINCE_PASSIVE_MARGIN_LOWLAND
    )
    old_orogen_support = land & (province_code == CONT_PROVINCE_OLD_OROGEN)
    largest_basin = 0.0
    for comp in module._components(grid, basin_support):
        largest_basin = max(largest_basin, float(grid.cell_area[comp].sum()))

    assert cont_fraction(basin_support) > 0.08
    assert float(largest_basin / cont_area) > 0.06
    assert cont_fraction(platform_support) > 0.16
    assert cont_fraction(passive_support) > 0.015
    assert cont_fraction(old_orogen_support) > 0.04

    captured = {}

    def fake_p174(world_arg, p104f_surface, *args, **kwargs):
        captured["surface_before_p174"] = np.asarray(
            p104f_surface, dtype=np.float64
        ).copy()
        return p104f_surface

    monkeypatch.setattr(
        module,
        "_p174_process_parented_lowland_plain_response",
        fake_p174,
    )
    result = module._apply_inland_landform_region_elevation_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_origin,
        crust_stability,
        sediment,
        np.full(n, -1.0, dtype=np.float64),
        detail,
        detail,
        province_graph,
        continental_landforms,
    )

    pre_p174 = captured["surface_before_p174"]
    assert np.array_equal(result["surface"] >= 0.0, surface >= 0.0)
    assert float(pre_p174[subsiding_basin].mean()) < 1450.0
    assert (
        metrics["terrain.last_p104f_pre_p174_lowland_prep_area_fraction"]
        > 0.08
    )
    assert (
        metrics["terrain.last_p104f_pre_p174_lowland_prep_mean_lowering_m"]
        > 600.0
    )


def test_p174_upstream_lowland_support_forms_connected_candidate_before_p174(monkeypatch):
    grid = SphereGrid.fibonacci(3200, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -170.0)
        & (grid.lon < 120.0)
        & (grid.lat > -58.0)
        & (grid.lat < 62.0)
    )
    lowland_belt = (
        land
        & (grid.lon > -142.0)
        & (grid.lon < 45.0)
        & (grid.lat > -36.0)
        & (grid.lat < 38.0)
    )
    basin_west = (
        lowland_belt
        & (grid.lon > -132.0)
        & (grid.lon < -72.0)
        & (grid.lat > -28.0)
        & (grid.lat < 30.0)
    )
    basin_central = (
        lowland_belt
        & (grid.lon > -56.0)
        & (grid.lon < 0.0)
        & (grid.lat > -30.0)
        & (grid.lat < 32.0)
    )
    basin_east = (
        lowland_belt
        & (grid.lon > 14.0)
        & (grid.lon < 42.0)
        & (grid.lat > -24.0)
        & (grid.lat < 30.0)
    )
    basin = basin_west | basin_central | basin_east
    platform_corridor = lowland_belt & ~basin
    far_platform = land & ~module._dilate_mask(grid, lowland_belt, passes=2)
    assert np.count_nonzero(basin) >= 420
    assert np.count_nonzero(platform_corridor) >= 340
    assert np.count_nonzero(far_platform) >= 300

    surface = np.where(land, 2380.0, -4300.0).astype(np.float64)
    surface[basin] = 2680.0
    surface[platform_corridor] = 1980.0
    separator_ribs = platform_corridor & (
        (np.floor((grid.lon + 180.0) / 13.0) % 5) == 0
    )
    surface[separator_ribs] = 2220.0

    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_origin = np.zeros(n, dtype=np.float64)
    crust_stability = np.where(land, 0.66, 0.0).astype(np.float64)
    crust_thick = np.where(land, 35000.0, 7000.0).astype(np.float64)
    crust_thick[basin] = 30000.0
    sediment = np.where(land, 130.0, 40.0).astype(np.float64)
    sediment[basin] = 1750.0
    sediment[platform_corridor] = 520.0
    platform_subsidence = np.zeros(n, dtype=np.float64)
    platform_subsidence[basin] = 0.62
    platform_subsidence[platform_corridor] = 0.08
    fields = {
        "crust.thickness_m": crust_thick,
        "crust.age_myr": np.where(land, 1900.0, 120.0).astype(np.float64),
        "crust.stability": crust_stability,
        "crust.type": crust_type,
        "sediment.thickness_m": sediment,
        "tectonics.platform_subsidence": platform_subsidence,
        "tectonics.rift_potential": np.zeros(n, dtype=np.float64),
        "tectonics.internal_geographic_block_code": np.where(
            land, INTERNAL_BLOCK_STABLE_PLATFORM, 0.0
        ).astype(np.float64),
        "tectonics.continent_id": np.where(land, 1.0, -1.0).astype(np.float64),
        "tectonics.deformation_style": np.zeros(n, dtype=np.float64),
        "tectonics.deformation_intensity": np.zeros(n, dtype=np.float64),
    }
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=2100.0,
        objects={},
        networks={},
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[basin] = CONT_DETAIL_BASIN
    province_code = np.zeros(n, dtype=np.float64)
    province_code[land] = CONT_PROVINCE_PLATFORM
    province_code[basin] = CONT_PROVINCE_INTRACRATONIC_BASIN
    province_graph = {
        "province_id": np.where(land, 1.0, 0.0).astype(np.float64),
        "province_code": province_code,
        "objects": [],
    }
    captured = {}

    def fake_p174(world_arg, p104f_surface, *args, **kwargs):
        captured["surface_before_p174"] = np.asarray(
            p104f_surface, dtype=np.float64
        ).copy()
        return p104f_surface

    monkeypatch.setattr(
        module,
        "_p174_process_parented_lowland_plain_response",
        fake_p174,
    )
    result = module._apply_inland_landform_region_elevation_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_origin,
        crust_stability,
        sediment,
        np.full(n, -1.0, dtype=np.float64),
        detail,
        detail,
        province_graph,
        [],
    )

    pre_p174 = captured["surface_before_p174"]
    support = fields["terrain.p104f_pre_p174_lowland_support_mask"] > 0.5
    memory = fields["terrain.p174_lowland_plain_continuity_memory"] >= 0.25
    cont_area = max(float(grid.cell_area[land].sum()), 1.0e-12)
    largest_support = 0.0
    for comp in module._components(grid, support):
        largest_support = max(largest_support, float(grid.cell_area[comp].sum()))

    assert np.array_equal(result["surface"] >= 0.0, surface >= 0.0)
    assert float(pre_p174[basin].mean()) < 1550.0
    assert float(pre_p174[platform_corridor].mean()) < 1900.0
    assert float(pre_p174[far_platform].mean()) > 1900.0
    assert metrics[
        "terrain.last_p104f_pre_p174_lowland_support_area_fraction"
    ] > 0.14
    assert metrics[
        "terrain.last_p104f_pre_p174_lowland_support_largest_component_fraction"
    ] > 0.10
    assert float(largest_support / cont_area) > 0.10
    assert float(grid.cell_area[memory & support].sum() / cont_area) > 0.10
    assert metrics[
        "terrain.last_p104f_pre_p174_lowland_memory_seed_area_fraction"
    ] > 0.10


def test_p174_upstream_support_bridges_wide_same_basin_platform_apron(monkeypatch):
    grid = SphereGrid.fibonacci(3600, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -172.0)
        & (grid.lon < 128.0)
        & (grid.lat > -58.0)
        & (grid.lat < 64.0)
    )
    lowland_apron = (
        land
        & (grid.lon > -150.0)
        & (grid.lon < 52.0)
        & (grid.lat > -34.0)
        & (grid.lat < 36.0)
    )
    basin_west = (
        lowland_apron
        & (grid.lon > -142.0)
        & (grid.lon < -112.0)
        & (grid.lat > -24.0)
        & (grid.lat < 27.0)
    )
    basin_central = (
        lowland_apron
        & (grid.lon > -54.0)
        & (grid.lon < -24.0)
        & (grid.lat > -26.0)
        & (grid.lat < 29.0)
    )
    basin_east = (
        lowland_apron
        & (grid.lon > 26.0)
        & (grid.lon < 48.0)
        & (grid.lat > -20.0)
        & (grid.lat < 25.0)
    )
    basin = basin_west | basin_central | basin_east
    platform_apron = lowland_apron & ~module._dilate_mask(grid, basin, passes=3)
    apron_margin = lowland_apron & ~basin & ~platform_apron
    far_platform = land & ~module._dilate_mask(grid, lowland_apron, passes=3)
    assert np.count_nonzero(basin) >= 250
    assert np.count_nonzero(platform_apron) >= 340
    assert np.count_nonzero(apron_margin) >= 180
    assert np.count_nonzero(far_platform) >= 350

    surface = np.where(land, 2380.0, -4300.0).astype(np.float64)
    surface[basin] = 2500.0
    surface[platform_apron | apron_margin] = 1780.0
    surface[far_platform] = 2320.0

    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_origin = np.zeros(n, dtype=np.float64)
    crust_stability = np.where(land, 0.70, 0.0).astype(np.float64)
    crust_thick = np.where(land, 35000.0, 7000.0).astype(np.float64)
    crust_thick[basin] = 30500.0
    sediment = np.where(land, 120.0, 35.0).astype(np.float64)
    sediment[basin] = 1650.0
    sediment[apron_margin] = 470.0
    sediment[platform_apron] = 135.0
    platform_subsidence = np.zeros(n, dtype=np.float64)
    platform_subsidence[basin] = 0.52
    platform_subsidence[apron_margin] = 0.045
    platform_subsidence[platform_apron] = 0.012
    fields = {
        "crust.thickness_m": crust_thick,
        "crust.age_myr": np.where(land, 2100.0, 120.0).astype(np.float64),
        "crust.stability": crust_stability,
        "crust.type": crust_type,
        "sediment.thickness_m": sediment,
        "tectonics.platform_subsidence": platform_subsidence,
        "tectonics.rift_potential": np.zeros(n, dtype=np.float64),
        "tectonics.internal_geographic_block_code": np.where(
            land, INTERNAL_BLOCK_STABLE_PLATFORM, 0.0
        ).astype(np.float64),
        "tectonics.continent_id": np.where(land, 1.0, -1.0).astype(np.float64),
        "tectonics.deformation_style": np.zeros(n, dtype=np.float64),
        "tectonics.deformation_intensity": np.zeros(n, dtype=np.float64),
    }
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=1500.0,
        objects={},
        networks={},
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[basin] = CONT_DETAIL_BASIN
    province_code = np.zeros(n, dtype=np.float64)
    province_code[land] = CONT_PROVINCE_PLATFORM
    province_code[basin] = CONT_PROVINCE_INTRACRATONIC_BASIN
    province_graph = {
        "province_id": np.where(land, 1.0, 0.0).astype(np.float64),
        "province_code": province_code,
        "objects": [],
    }
    captured = {}

    def fake_p174(world_arg, p104f_surface, *args, **kwargs):
        captured["surface_before_p174"] = np.asarray(
            p104f_surface, dtype=np.float64
        ).copy()
        return p104f_surface

    monkeypatch.setattr(
        module,
        "_p174_process_parented_lowland_plain_response",
        fake_p174,
    )
    result = module._apply_inland_landform_region_elevation_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_origin,
        crust_stability,
        sediment,
        np.full(n, -1.0, dtype=np.float64),
        detail,
        detail,
        province_graph,
        [],
    )

    pre_p174 = captured["surface_before_p174"]
    support = fields["terrain.p104f_pre_p174_lowland_support_mask"] > 0.5
    memory = fields["terrain.p174_lowland_plain_continuity_memory"] >= 0.25
    cont_area = max(float(grid.cell_area[land].sum()), 1.0e-12)
    largest_support = 0.0
    for comp in module._components(grid, support):
        largest_support = max(largest_support, float(grid.cell_area[comp].sum()))

    assert np.array_equal(result["surface"] >= 0.0, surface >= 0.0)
    assert float(pre_p174[basin].mean()) < 1550.0
    assert float(pre_p174[platform_apron].mean()) < 1800.0
    assert float(pre_p174[far_platform].mean()) > 1850.0
    assert float(largest_support / cont_area) > 0.14
    assert float(grid.cell_area[support & platform_apron].sum() / cont_area) > 0.06
    assert float(grid.cell_area[support & far_platform].sum() / cont_area) < 0.01
    assert float(grid.cell_area[memory & support].sum() / cont_area) > 0.10
    assert metrics[
        "terrain.last_p104f_pre_p174_lowland_support_largest_component_fraction"
    ] > 0.14


def test_p174_support_to_plain_conversion_uses_connected_upstream_support():
    grid = SphereGrid.fibonacci(3200, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -170.0)
        & (grid.lon < 120.0)
        & (grid.lat > -58.0)
        & (grid.lat < 62.0)
    )
    support = (
        land
        & (grid.lon > -140.0)
        & (grid.lon < 48.0)
        & (grid.lat > -36.0)
        & (grid.lat < 38.0)
    )
    basin_seed = (
        support
        & (
            (
                (grid.lon > -124.0)
                & (grid.lon < -78.0)
                & (grid.lat > -26.0)
                & (grid.lat < 30.0)
            )
            | (
                (grid.lon > -38.0)
                & (grid.lon < 4.0)
                & (grid.lat > -30.0)
                & (grid.lat < 28.0)
            )
        )
    )
    support_corridor = support & ~basin_seed
    far_platform = land & ~module._dilate_mask(grid, support, passes=3)
    assert np.count_nonzero(support) >= 850
    assert np.count_nonzero(basin_seed) >= 260
    assert np.count_nonzero(support_corridor) >= 500
    assert np.count_nonzero(far_platform) >= 300

    surface = np.where(land, 2380.0, -4300.0).astype(np.float64)
    surface[support] = 2620.0
    surface[basin_seed] = 2840.0
    low_seed = basin_seed & (
        (np.floor((grid.lon + 180.0) / 17.0) % 3) == 0
    )
    surface[low_seed] = 980.0
    corridor_ribs = support_corridor & (
        (
            np.floor((grid.lon + 180.0) / 10.0)
            + np.floor((grid.lat + 90.0) / 12.0)
        )
        % 4
        == 0
    )
    surface[corridor_ribs] = 3050.0

    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.68, 0.0).astype(np.float64)
    sediment = np.where(land, 110.0, 35.0).astype(np.float64)
    sediment[basin_seed] = 1550.0
    sediment[support_corridor] = 360.0
    orog_age = np.full(n, -1.0, dtype=np.float64)
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[basin_seed] = CONT_DETAIL_BASIN
    province = np.zeros(n, dtype=np.float64)
    province[land] = CONT_PROVINCE_PLATFORM
    province[basin_seed] = CONT_PROVINCE_INTRACRATONIC_BASIN
    inland = np.zeros(n, dtype=np.float64)
    inland[land] = INLAND_PROVINCE_PLATFORM
    inland[basin_seed] = INLAND_PROVINCE_SAG_BASIN

    platform_subsidence = np.zeros(n, dtype=np.float64)
    platform_subsidence[basin_seed] = 0.10
    platform_subsidence[support_corridor] = 0.02
    fields = {
        "tectonics.platform_subsidence": platform_subsidence,
        "terrain.p104f_pre_p174_lowland_support_mask": support.astype(np.float64),
        "terrain.p174_lowland_plain_continuity_memory": np.where(
            basin_seed,
            0.45,
            0.0,
        ).astype(np.float64),
    }
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=2100.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    frame_fields = {
        "crust.type": crust_type,
        "crust.age_myr": np.where(land, 1900.0, 120.0),
        "crust.domain": crust_domain,
        "crust.origin": np.zeros(n, dtype=np.float64),
        "terrain.continental_detail": detail,
        "terrain.continental_detail_region_code": detail,
        "terrain.inland_geomorphology_region_code": inland,
        "terrain.continental_province_code": province,
        "tectonics.orogeny_age_myr": orog_age,
        "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
        "sediment.thickness_m": sediment,
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS
        ).astype(np.float64),
    }
    before_row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=2100.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={"terrain.elevation_m": surface, **frame_fields},
        ),
    )

    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        orog_age,
        detail,
        province,
        inland,
        [],
        stage_code=2.0,
    )
    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=2100.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={
                "terrain.elevation_m": out,
                **frame_fields,
                "terrain.p174_lowland_plain_continuity_memory": fields.get(
                    "terrain.p174_lowland_plain_continuity_memory",
                    np.zeros(n, dtype=np.float64),
                ),
                "terrain.p104f_pre_p174_lowland_support_mask": fields[
                    "terrain.p104f_pre_p174_lowland_support_mask"
                ],
            },
        ),
    )
    before_metrics = before_row["land_metrics"]
    land_metrics = row["land_metrics"]
    cont_area = max(float(grid.cell_area[land].sum()), 1.0e-12)
    support_area_fraction = float(grid.cell_area[support].sum() / cont_area)

    assert before_metrics["lowland_plain_fraction"] < 0.04
    assert land_metrics["lowland_plain_fraction"] > 0.10
    assert land_metrics["largest_lowland_plain_component_fraction"] > min(
        0.075,
        0.58 * support_area_fraction,
    )
    assert land_metrics["broad_lowland_plain_component_count"] >= 1
    assert land_metrics["lowland_plain_parented_fraction"] > 0.75
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]
    assert float(out[far_platform].mean()) > 1450.0
    assert metrics["terrain.last_p174_lowland_plain_response_area_fraction"] > 0.08


def test_p174_active_orogen_basin_support_converts_without_lowering_highland():
    grid = SphereGrid.fibonacci(2600, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -178.0)
        & (grid.lon < 176.0)
        & (grid.lat > -54.0)
        & (grid.lat < 60.0)
    )
    active_basin = (
        land
        & (grid.lon > -112.0)
        & (grid.lon < 42.0)
        & (grid.lat > -34.0)
        & (grid.lat < 36.0)
    )
    active_highland = (
        land
        & (grid.lon > 58.0)
        & (grid.lon < 112.0)
        & (grid.lat > -26.0)
        & (grid.lat < 42.0)
    )
    assert np.count_nonzero(active_basin) >= 520
    assert np.count_nonzero(active_highland) >= 90

    surface = np.where(land, 2100.0, -4300.0).astype(np.float64)
    surface[active_basin] = 1850.0
    surface[active_highland] = 2750.0
    low_seed = active_basin & (
        (
            np.floor((grid.lon + 180.0) / 16.0)
            + np.floor((grid.lat + 90.0) / 14.0)
        )
        % 5
        == 0
    )
    surface[low_seed] = 420.0

    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_domain[active_basin | active_highland] = DOMAIN_ACCRETED_TERRANE
    crust_stability = np.where(land, 0.62, 0.0).astype(np.float64)
    sediment = np.where(land, 120.0, 35.0).astype(np.float64)
    sediment[active_basin] = 1300.0
    orog_age = np.full(n, -1.0, dtype=np.float64)
    orog_age[active_highland] = 120.0
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[active_basin] = CONT_DETAIL_OROGEN
    detail[active_highland] = CONT_DETAIL_OROGEN
    province = np.zeros(n, dtype=np.float64)
    province[land] = CONT_PROVINCE_PLATFORM
    province[active_basin | active_highland] = CONT_PROVINCE_ACTIVE_OROGEN
    inland = np.zeros(n, dtype=np.float64)
    inland[land] = INLAND_PROVINCE_PLATFORM
    inland[active_basin] = INLAND_PROVINCE_SAG_BASIN

    platform_subsidence = np.zeros(n, dtype=np.float64)
    platform_subsidence[active_basin] = 0.06
    fields = {
        "tectonics.platform_subsidence": platform_subsidence,
        "terrain.p104f_pre_p174_lowland_support_mask": (
            active_basin.astype(np.float64)
        ),
        "terrain.p174_lowland_plain_continuity_memory": np.where(
            low_seed,
            0.45,
            0.0,
        ).astype(np.float64),
    }
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=1700.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    frame_fields = {
        "crust.type": crust_type,
        "crust.age_myr": np.where(land, 1600.0, 110.0),
        "crust.domain": crust_domain,
        "crust.origin": np.zeros(n, dtype=np.float64),
        "terrain.continental_detail": detail,
        "terrain.continental_detail_region_code": detail,
        "terrain.inland_geomorphology_region_code": inland,
        "terrain.continental_province_code": province,
        "tectonics.orogeny_age_myr": orog_age,
        "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
        "sediment.thickness_m": sediment,
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS
        ).astype(np.float64),
    }

    before = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1700.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={"terrain.elevation_m": surface, **frame_fields},
        ),
    )["land_metrics"]
    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        orog_age,
        detail,
        province,
        inland,
        [],
        stage_code=2.0,
    )
    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1700.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={
                "terrain.elevation_m": out,
                **frame_fields,
                "terrain.p174_lowland_plain_continuity_memory": fields[
                    "terrain.p174_lowland_plain_continuity_memory"
                ],
                "terrain.p104f_pre_p174_lowland_support_mask": fields[
                    "terrain.p104f_pre_p174_lowland_support_mask"
                ],
            },
        ),
    )
    land_metrics = row["land_metrics"]

    assert before["largest_lowland_plain_component_fraction"] < 0.025
    assert land_metrics["lowland_plain_fraction"] > 0.08
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.04
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]
    assert float(out[active_highland].mean()) > 2500.0
    assert metrics["terrain.last_p174_lowland_plain_response_area_fraction"] > 0.08


def test_p174_area_limited_lowland_topup_uses_adjacent_process_apron():
    grid = SphereGrid.fibonacci(2800, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -172.0)
        & (grid.lon < 132.0)
        & (grid.lat > -56.0)
        & (grid.lat < 62.0)
    )
    low_source = (
        land
        & (grid.lon > -118.0)
        & (grid.lon < -68.0)
        & (grid.lat > -12.0)
        & (grid.lat < 18.0)
    )
    process_apron = (
        land
        & (grid.lon > -142.0)
        & (grid.lon < 8.0)
        & (grid.lat > -34.0)
        & (grid.lat < 40.0)
    )
    process_apron &= ~low_source
    far_platform = land & (grid.lon > 78.0) & (grid.lat < -14.0)
    assert np.count_nonzero(low_source) >= 80
    assert np.count_nonzero(process_apron) >= 420
    assert np.count_nonzero(far_platform) >= 80

    surface = np.where(land, 2150.0, -4300.0).astype(np.float64)
    surface[low_source] = 260.0
    surface[process_apron] = 1320.0
    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.64, 0.0).astype(np.float64)
    sediment = np.where(land, 100.0, 35.0).astype(np.float64)
    sediment[low_source | process_apron] = 1150.0
    orog_age = np.full(n, -1.0, dtype=np.float64)
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[low_source | process_apron] = CONT_DETAIL_BASIN
    province = np.zeros(n, dtype=np.float64)
    province[land] = CONT_PROVINCE_PLATFORM
    province[low_source | process_apron] = CONT_PROVINCE_INTRACRATONIC_BASIN
    inland = np.zeros(n, dtype=np.float64)
    inland[land] = INLAND_PROVINCE_PLATFORM
    inland[low_source | process_apron] = INLAND_PROVINCE_SAG_BASIN
    platform_subsidence = np.zeros(n, dtype=np.float64)
    platform_subsidence[low_source | process_apron] = 0.04
    fields = {
        "tectonics.platform_subsidence": platform_subsidence,
        "terrain.p104f_pre_p174_lowland_support_mask": (
            low_source.astype(np.float64)
        ),
        "terrain.p174_lowland_plain_continuity_memory": np.where(
            low_source,
            0.45,
            0.0,
        ).astype(np.float64),
    }
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=1200.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    frame_fields = {
        "crust.type": crust_type,
        "crust.age_myr": np.where(land, 1200.0, 110.0),
        "crust.domain": crust_domain,
        "crust.origin": np.zeros(n, dtype=np.float64),
        "terrain.continental_detail": detail,
        "terrain.continental_detail_region_code": detail,
        "terrain.inland_geomorphology_region_code": inland,
        "terrain.continental_province_code": province,
        "tectonics.orogeny_age_myr": orog_age,
        "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
        "sediment.thickness_m": sediment,
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS
        ).astype(np.float64),
    }
    before = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1200.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={"terrain.elevation_m": surface, **frame_fields},
        ),
    )["land_metrics"]
    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        orog_age,
        detail,
        province,
        inland,
        [],
        stage_code=1.0,
    )
    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1200.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={
                "terrain.elevation_m": out,
                **frame_fields,
                "terrain.p174_lowland_plain_continuity_memory": fields[
                    "terrain.p174_lowland_plain_continuity_memory"
                ],
                "terrain.p104f_pre_p174_lowland_support_mask": fields[
                    "terrain.p104f_pre_p174_lowland_support_mask"
                ],
            },
        ),
    )
    land_metrics = row["land_metrics"]

    assert before["lowland_plain_fraction"] < 0.06
    assert land_metrics["lowland_plain_fraction"] >= 0.064
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.025
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]
    assert float(out[far_platform].mean()) > 1800.0


def test_p174_residual_relief_boundary_edge_cap_connects_near_floor_plain():
    grid = SphereGrid.fibonacci(3600, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = (
        (grid.lon > -166.0)
        & (grid.lon < 104.0)
        & (grid.lat > -52.0)
        & (grid.lat < 56.0)
    )
    problem_plain = (
        land
        & (grid.lon > -162.0)
        & (grid.lon < 80.0)
        & (grid.lat > -48.0)
        & (grid.lat < 54.0)
    )
    relief_edges = problem_plain & (
        (
            np.floor((grid.lon + 180.0) / 12.0)
            + np.floor((grid.lat + 90.0) / 14.0)
        )
        % 2
        == 0
    )
    support_source = problem_plain & ~relief_edges
    bonus_plain = (
        land
        & (
            (
                (grid.lon > 84.0)
                & (grid.lon < 102.0)
                & (grid.lat > -48.0)
                & (grid.lat < -30.0)
            )
            | (
                (grid.lon > 84.0)
                & (grid.lon < 102.0)
                & (grid.lat > -16.0)
                & (grid.lat < 2.0)
            )
            | (
                (grid.lon > -118.0)
                & (grid.lon < -96.0)
                & (grid.lat > -50.0)
                & (grid.lat < -34.0)
            )
        )
    )
    auxiliary_plain = (
        land
        & (
            (
                (grid.lon > 22.0)
                & (grid.lon < 76.0)
                & (grid.lat > -48.0)
                & (grid.lat < -12.0)
            )
            | (
                (grid.lon > 82.0)
                & (grid.lon < 118.0)
                & (grid.lat > 8.0)
                & (grid.lat < 44.0)
            )
            | (
                (grid.lon > -162.0)
                & (grid.lon < -124.0)
                & (grid.lat > 18.0)
                & (grid.lat < 54.0)
            )
        )
    )
    assert np.count_nonzero(problem_plain) >= 420
    assert np.count_nonzero(relief_edges) >= 120
    assert np.count_nonzero(support_source) >= 260
    assert np.count_nonzero(bonus_plain) >= 30
    assert np.count_nonzero(auxiliary_plain) >= 140

    surface = np.where(land, 1650.0, -4300.0).astype(np.float64)
    surface[problem_plain] = 140.0
    surface[relief_edges] = 660.0
    surface[auxiliary_plain] = 190.0
    surface[bonus_plain] = 190.0
    aux_ribs = auxiliary_plain & (
        (
            np.floor((grid.lon + 180.0) / 9.0)
            + np.floor((grid.lat + 90.0) / 9.0)
        )
        % 2
        == 0
    )
    surface[aux_ribs] = 650.0

    crust_type = land.astype(np.float64)
    crust_domain = np.where(
        land, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC
    ).astype(np.float64)
    crust_stability = np.where(land, 0.66, 0.0).astype(np.float64)
    sediment = np.where(land, 90.0, 35.0).astype(np.float64)
    sediment[problem_plain | auxiliary_plain | bonus_plain] = 980.0
    orog_age = np.full(n, -1.0, dtype=np.float64)
    detail = np.zeros(n, dtype=np.float64)
    detail[land] = CONT_DETAIL_PLATFORM
    detail[problem_plain | auxiliary_plain | bonus_plain] = CONT_DETAIL_BASIN
    province = np.zeros(n, dtype=np.float64)
    province[land] = CONT_PROVINCE_PLATFORM
    province[problem_plain | auxiliary_plain | bonus_plain] = (
        CONT_PROVINCE_INTRACRATONIC_BASIN
    )
    inland = np.zeros(n, dtype=np.float64)
    inland[land] = INLAND_PROVINCE_PLATFORM
    inland[problem_plain | auxiliary_plain | bonus_plain] = INLAND_PROVINCE_SAG_BASIN
    platform_subsidence = np.zeros(n, dtype=np.float64)
    platform_subsidence[support_source] = 0.06
    platform_subsidence[auxiliary_plain] = 0.04
    platform_subsidence[bonus_plain] = 0.04
    fields = {
        "tectonics.platform_subsidence": platform_subsidence,
        "terrain.p104f_pre_p174_lowland_support_mask": (
            support_source.astype(np.float64)
        ),
        "terrain.p174_lowland_plain_continuity_memory": np.where(
            support_source,
            0.45,
            0.0,
        ).astype(np.float64),
    }
    metrics = {}

    def get_field(name, default=0.0):
        return fields.get(name, np.full(n, default, dtype=np.float64))

    def set_field(name, values):
        fields[name] = np.asarray(values, dtype=np.float64)

    world = SimpleNamespace(
        grid=grid,
        time_myr=1600.0,
        get_field=get_field,
        set_field=set_field,
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    frame_fields = {
        "crust.type": crust_type,
        "crust.age_myr": np.where(land, 1700.0, 120.0),
        "crust.domain": crust_domain,
        "crust.origin": np.zeros(n, dtype=np.float64),
        "terrain.continental_detail": detail,
        "terrain.continental_detail_region_code": detail,
        "terrain.inland_geomorphology_region_code": inland,
        "terrain.continental_province_code": province,
        "tectonics.orogeny_age_myr": orog_age,
        "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
        "sediment.thickness_m": sediment,
        "ocean.depth_province": np.where(
            land, 0.0, OCEAN_DEPTH_ABYSS
        ).astype(np.float64),
    }
    before_row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1600.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={"terrain.elevation_m": surface, **frame_fields},
        ),
    )
    before = before_row["land_metrics"]

    out = module._p174_process_parented_lowland_plain_response(
        world,
        surface,
        0.0,
        crust_type,
        crust_domain,
        crust_stability,
        sediment,
        orog_age,
        detail,
        province,
        inland,
        [],
        stage_code=2.0,
    )
    row = frame_geomorphology_metrics(
        grid,
        Frame(
            time_myr=1600.0,
            globals={"ocean.sea_level_m": 0.0},
            fields={
                "terrain.elevation_m": out,
                **frame_fields,
                "terrain.p174_lowland_plain_continuity_memory": fields.get(
                    "terrain.p174_lowland_plain_continuity_memory",
                    np.zeros(n, dtype=np.float64),
                ),
                "terrain.p104f_pre_p174_lowland_support_mask": fields[
                    "terrain.p104f_pre_p174_lowland_support_mask"
                ],
            },
        ),
    )
    land_metrics = row["land_metrics"]

    assert before["lowland_plain_fraction"] >= 0.105
    assert before["largest_lowland_plain_component_fraction"] < 0.025
    assert before[
        "largest_lowland_elevation_parented_component_fraction"
    ] >= 0.025
    assert before["lowland_residual_relief_boundary_limited"] == 1.0
    assert land_metrics["lowland_plain_fraction"] >= before[
        "lowland_plain_fraction"
    ]
    assert land_metrics["largest_lowland_plain_component_fraction"] > 0.03
    assert not row["diagnostic_flags"]["lowland_plain_deficient"]
    assert float(out[relief_edges].mean()) < 420.0


def test_p173_ocean_object_lifecycle_stamps_ocean_stage_and_parent():
    grid = SphereGrid.fibonacci(90, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    world = SimpleNamespace(
        grid=grid,
        time_myr=1800.0,
        objects={"terrain.ocean_fabric": []},
    )
    current = [{
        "id": "fracture_zone:14:0",
        "kind": "fracture_zone",
        "cells": [14, 15, 16],
        "mean_age_myr": 85.0,
        "parent_basin_id": 3,
        "parent_process": "inactive_transform_trace_and_age_offset",
    }]

    stabilized = module._stabilize_process_objects(
        world,
        "terrain.ocean_fabric",
        current,
        id_prefix="ocean_fabric",
    )[0]

    assert stabilized["id"].startswith("ocean_fabric:fracture_zone:")
    assert stabilized["birth_myr"] == 1715.0
    assert stabilized["age_myr"] == 85.0
    assert stabilized["parent_plate_id"] == "basin:3"
    assert stabilized["activity_state"] == "decaying"
    assert stabilized["relief_stage"] == "inactive_fracture_zone"


def test_p173_age_aware_ocean_response_shapes_objects_and_preserves_mask():
    grid = SphereGrid.fibonacci(240, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = grid.lat > 65.0
    ocean = ~land
    ridge_cells = np.where(ocean & (grid.lon < -100.0))[0][:14]
    trench_cells = np.where(ocean & (grid.lon > -20.0) & (grid.lon < 40.0))[0][:14]
    plateau_cells = np.where(ocean & (grid.lon > 95.0))[0][:14]
    assert ridge_cells.size >= 6
    assert trench_cells.size >= 6
    assert plateau_cells.size >= 6
    surface = np.where(land, 540.0, -3900.0).astype(np.float64)
    surface[ridge_cells] = -4300.0
    surface[trench_cells] = -4200.0
    surface[plateau_cells] = -3600.0
    fields = {
        "crust.age_myr": np.where(land, 1200.0, 90.0).astype(np.float64),
        "sediment.thickness_m": np.zeros(n, dtype=np.float64),
    }
    fields["crust.age_myr"][ridge_cells] = 4.0
    fields["crust.age_myr"][trench_cells] = 80.0
    fields["crust.age_myr"][plateau_cells] = 350.0
    metrics = {}
    world = SimpleNamespace(
        grid=grid,
        time_myr=2300.0,
        objects={
            "terrain.ocean_fabric": [
                {
                    "kind": "spreading_center",
                    "age_myr": 5.0,
                    "cells": ridge_cells.tolist(),
                },
                {
                    "kind": "fracture_zone",
                    "age_myr": 95.0,
                    "cells": trench_cells[:6].tolist(),
                },
            ],
            "terrain.margin_landforms": [
                {
                    "kind": "trench",
                    "age_myr": 55.0,
                    "cells": trench_cells.tolist(),
                },
            ],
            "terrain.arc_plume_landforms": [
                {
                    "kind": "oceanic_plateau",
                    "age_myr": 420.0,
                    "cells": plateau_cells.tolist(),
                },
            ],
        },
        get_field=lambda name, default=0.0: fields.get(
            name,
            np.full(n, default, dtype=np.float64),
        ),
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    coast_pass = np.where(ocean, 5, 0).astype(np.int16)

    out = module._p173_age_aware_ocean_floor_lifecycle_response(
        world,
        surface,
        0.0,
        ocean,
        coast_pass,
    )

    assert float(out[ridge_cells].mean()) > float(surface[ridge_cells].mean()) + 250.0
    assert float(out[plateau_cells].mean()) > float(surface[plateau_cells].mean()) + 250.0
    trench_core = trench_cells[:6]
    assert float(out[trench_core].mean()) < float(surface[trench_core].mean()) - 250.0
    assert np.array_equal(out >= 0.0, surface >= 0.0)
    assert metrics["terrain.last_p173_age_aware_ocean_response_object_count"] == 4.0
    assert metrics["terrain.last_p173_age_aware_ocean_response_area_fraction"] > 0.0


def test_p173_unsupported_ocean_shoal_cleanup_preserves_object_backed_shoals():
    grid = SphereGrid.fibonacci(240, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = grid.lat > 70.0
    ocean = ~land
    backed = np.where(
        ocean & (np.abs(grid.lat) < 15.0) & (grid.lon < -15.0) & (grid.lon > -75.0)
    )[0][:12]
    backed_rift = np.where(
        ocean & (np.abs(grid.lat) < 15.0) & (grid.lon < -95.0) & (grid.lon > -155.0)
    )[0][:12]
    unsupported = np.where(
        ocean & (np.abs(grid.lat) < 15.0) & (grid.lon > 15.0) & (grid.lon < 75.0)
    )[0][:12]
    assert backed.size >= 6
    assert backed_rift.size >= 6
    assert unsupported.size >= 6

    surface = np.where(land, 640.0, -4200.0).astype(np.float64)
    surface[backed] = -900.0
    surface[backed_rift] = -900.0
    surface[unsupported] = -900.0
    fields = {
        "crust.age_myr": np.where(land, 1600.0, 90.0).astype(np.float64),
    }
    metrics = {}
    world = SimpleNamespace(
        grid=grid,
        objects={
            "terrain.ocean_fabric": [],
            "terrain.margin_landforms": [],
            "terrain.arc_plume_landforms": [
                {"kind": "seamount_chain", "cells": backed.astype(int).tolist()},
            ],
            "terrain.rift_margin_sequences": [
                {
                    "kind": "rift_margin_sequence",
                    "cells": backed_rift.astype(int).tolist(),
                },
            ],
        },
        get_field=lambda name, default=0.0: fields.get(
            name,
            np.full(n, default, dtype=np.float64),
        ),
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    coast_pass = np.where(ocean, 5, 0).astype(np.int16)

    out = module._p173_deepen_unsupported_ocean_shoals(
        world,
        surface,
        0.0,
        ocean,
        coast_pass,
        preserve=np.zeros(n, dtype=bool),
    )

    assert np.all(out[backed] == surface[backed])
    assert np.all(out[backed_rift] == surface[backed_rift])
    assert np.all(out[unsupported] <= -2600.0)
    assert np.all(out[land] == surface[land])
    assert metrics["terrain.last_p173_unsupported_ocean_shoal_candidate_fraction"] > 0.0
    assert metrics["terrain.last_p173_unsupported_ocean_shoal_deepened_fraction"] > 0.0


def test_p1732_young_open_ocean_age_depth_deepens_only_unsupported_residuals():
    grid = SphereGrid.fibonacci(320, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = grid.lat > 72.0
    ocean = ~land
    unsupported = np.where(
        ocean & (np.abs(grid.lat) < 12.0) & (grid.lon > 35.0) & (grid.lon < 95.0)
    )[0][:10]
    backed = np.where(
        ocean & (np.abs(grid.lat) < 12.0) & (grid.lon < -25.0) & (grid.lon > -85.0)
    )[0][:10]
    ridge_axis = np.where(
        ocean & (np.abs(grid.lat) < 12.0) & (grid.lon < -120.0)
    )[0][:8]
    assert unsupported.size >= 6
    assert backed.size >= 6
    assert ridge_axis.size >= 4

    surface = np.where(land, 620.0, -4300.0).astype(np.float64)
    surface[unsupported] = -180.0
    surface[backed] = -220.0
    surface[ridge_axis] = -900.0
    fields = {
        "crust.age_myr": np.where(land, 1600.0, 45.0).astype(np.float64),
    }
    fields["crust.age_myr"][ridge_axis] = 5.0
    metrics = {}
    world = SimpleNamespace(
        grid=grid,
        objects={
            "terrain.ocean_fabric": [],
            "terrain.margin_landforms": [],
            "terrain.arc_plume_landforms": [
                {"kind": "seamount_chain", "cells": backed.astype(int).tolist()},
            ],
            "terrain.rift_margin_sequences": [],
        },
        get_field=lambda name, default=0.0: fields.get(
            name,
            np.full(n, default, dtype=np.float64),
        ),
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )
    coast_pass = np.where(ocean, 5, 0).astype(np.int16)
    ridge_zone = np.zeros(n, dtype=bool)
    ridge_zone[ridge_axis] = True

    out = module._p173_enforce_young_open_ocean_age_depth(
        world,
        surface,
        0.0,
        ocean,
        coast_pass,
        preserve=np.zeros(n, dtype=bool),
        ridge_zone=ridge_zone,
    )

    assert np.all(out[unsupported] <= -2600.0)
    assert np.all(out[backed] == surface[backed])
    assert np.all(out[ridge_axis] == surface[ridge_axis])
    assert np.all(out[land] == surface[land])
    assert (
        metrics["terrain.last_p1732_young_open_ocean_age_depth_candidate_fraction"]
        > 0.0
    )
    assert (
        metrics["terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction"]
        > 0.0
    )
    assert (
        metrics["terrain.last_p1732_young_open_ocean_age_depth_mean_depth_after_m"]
        > metrics["terrain.last_p1732_young_open_ocean_age_depth_mean_depth_before_m"]
    )


def test_p1732_final_depth_floor_uses_current_objects_and_preserves_land_mask():
    grid = SphereGrid.fibonacci(320, CONSTANTS.EARTH_RADIUS)
    module = TerrainModule()
    n = grid.n
    land = grid.lat > 72.0
    ocean = ~land
    unsupported = np.where(
        ocean & (np.abs(grid.lat) < 12.0) & (grid.lon > 35.0) & (grid.lon < 95.0)
    )[0][:10]
    backed = np.where(
        ocean & (np.abs(grid.lat) < 12.0) & (grid.lon < -25.0) & (grid.lon > -85.0)
    )[0][:10]
    semantic_lip = np.where(
        ocean & (grid.lat < -15.0) & (grid.lat > -35.0) & (grid.lon > 95.0)
    )[0][:8]
    assert unsupported.size >= 6
    assert backed.size >= 6
    assert semantic_lip.size >= 4

    surface = np.where(land, 620.0, -4300.0).astype(np.float64)
    surface[unsupported] = -180.0
    surface[backed] = -220.0
    surface[semantic_lip] = -260.0
    fields = {
        "crust.age_myr": np.where(land, 1600.0, 45.0).astype(np.float64),
        "crust.domain": np.full(n, DOMAIN_OCEANIC, dtype=np.float64),
        "crust.origin": np.zeros(n, dtype=np.float64),
        "crust.type": np.zeros(n, dtype=np.float64),
        "tectonics.terrane_id": np.full(n, -1.0, dtype=np.float64),
    }
    fields["crust.domain"][semantic_lip] = float(DOMAIN_LIP)
    metrics = {}
    globals_ = {}
    world = SimpleNamespace(
        grid=grid,
        time_myr=1900.0,
        objects={
            "terrain.ocean_fabric": [],
            "terrain.margin_landforms": [],
            "terrain.arc_plume_landforms": [],
            "terrain.rift_margin_sequences": [],
            "tectonics.boundary_objects": [],
            "tectonics.boundary_polylines": [],
            "tectonics.spreading_centers": [],
            "tectonics.closing_margins": [],
        },
        networks={"tectonics.boundaries": {}},
        get_field=lambda name, default=0.0: fields.get(
            name,
            np.full(n, default, dtype=np.float64),
        ),
        g=lambda key, default=0.0: globals_.get(key, default),
        set_g=lambda key, value: metrics.__setitem__(key, value),
    )

    out = module._p1732_final_young_open_ocean_depth_floor(
        world,
        surface,
        0.0,
        arc_plume_landforms=[
            {"kind": "seamount_chain", "cells": backed.astype(int).tolist()},
        ],
    )

    assert np.array_equal(out >= 0.0, land)
    assert np.all(out[unsupported] <= -2600.0)
    assert np.all(out[backed] == surface[backed])
    assert np.all(out[semantic_lip] == surface[semantic_lip])
    assert metrics["terrain.last_p1732_young_open_ocean_age_depth_used"] == 1.0
    assert (
        metrics[
            "terrain.last_p1732_young_open_ocean_age_depth_land_mask_preserved"
        ]
        == 1.0
    )
