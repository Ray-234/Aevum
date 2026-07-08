import numpy as np
import pytest

from aevum import validation
from aevum.compiler.map_compiler import MapCompiler
from aevum.core.grid import SphereGrid
from aevum.core.state import WorldState
from aevum.core.units import CONSTANTS
from aevum.diagnostics.release_gate import P12RunConfig, run_p12_release_gate
from aevum.engine import Engine
from aevum.modules.biosphere import BiosphereModule
from aevum.modules.interior import InteriorModule
from aevum.modules.terrain import (
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATFORM,
    CONT_DETAIL_PLATEAU,
    TerrainModule,
)
from aevum.modules.tectonics import (
    DOMAIN_ACCRETED_TERRANE,
    DOMAIN_CONTINENTAL_INTERIOR,
    DOMAIN_CONTINENTAL_MARGIN,
    DOMAIN_CRATON,
    DOMAIN_LIP,
    DOMAIN_OCEANIC,
    ORIGIN_PLUME_IMPACT,
)
from aevum.render import to_raster
from aevum.spec.presets import get_preset


def _short(preset="earthlike", cells=2500, t_end=3000.0):
    spec = get_preset(preset)
    spec.grid_cells = cells
    spec.t_end_myr = t_end
    eng = Engine.build(spec)
    eng.run(n_frames=4)
    return eng


def test_earthlike_runs_and_conserves_carbon():
    eng = _short()
    assert validation.check_conservation(eng).passed
    assert validation.check_topology(eng).passed
    assert validation.check_causality(eng).passed


def test_scheduler_refreshes_state_at_final_time():
    eng = _short(cells=320, t_end=160.0)
    assert eng.world.time_myr == pytest.approx(160.0)
    assert eng.scheduler.history[-1].time_myr == pytest.approx(160.0)
    assert "tectonics" in eng.scheduler.history[-1].ran
    assert "terrain" in eng.scheduler.history[-1].ran
    assert eng.scheduler.history[-1].module_seconds["tectonics"] >= 0.0
    assert eng.scheduler.history[-1].module_seconds["terrain"] >= 0.0
    for scheduled in eng.scheduler.modules:
        assert scheduled.last_run_myr == pytest.approx(eng.world.time_myr)


def test_determinism_same_seed():
    a = _short().world.get_field("terrain.elevation_m")
    b = _short().world.get_field("terrain.elevation_m")
    assert np.allclose(a, b)


def test_compiles_a_map_with_layers():
    eng = _short()
    cm = MapCompiler(eng.world, eng.archive).compile(width=64, height=32, n_starts=4)
    assert cm.terrain.shape == (32, 64)
    assert cm.source_land_fraction.shape == (32, 64)
    assert cm.source_shelf_fraction.shape == (32, 64)
    assert cm.source_depth_province.shape == (32, 64)
    assert cm.source_terrain_province.shape == (32, 64)
    assert cm.source_continental_detail.shape == (32, 64)
    assert cm.food.max() > 0
    assert set(np.unique(cm.terrain)).issubset(set(range(6)))
    assert validation.check_compiler_consistency(eng).passed


def test_biosphere_biome_generalization_uses_neighbor_majority():
    grid = type("Grid", (), {})()
    grid.n = 6
    grid.edges = np.asarray([
        [0, 1],
        [0, 2],
        [0, 3],
        [0, 4],
        [4, 5],
    ], dtype=int)
    grid.neighbors = [
        np.asarray([1, 2, 3, 4], dtype=int),
        np.asarray([0], dtype=int),
        np.asarray([0], dtype=int),
        np.asarray([0], dtype=int),
        np.asarray([0, 5], dtype=int),
        np.asarray([4], dtype=int),
    ]
    biome = np.asarray([2.0, 4.0, 4.0, 4.0, 3.0, 1.0])
    ocean = np.asarray([False, False, False, False, False, True])
    fixed = np.asarray([False, False, False, False, False, True])

    out = BiosphereModule()._generalize_biomes(grid, biome, ocean, fixed)

    assert out[0] == 4.0
    assert out[5] == biome[5]


def test_terrain_mask_width_steps_uses_boundary_distance_and_cache():
    grid = type("Grid", (), {})()
    grid.n = 4
    grid.edges = np.asarray([[0, 1], [1, 2], [2, 3]], dtype=int)
    grid.neighbors = [
        np.asarray([1], dtype=int),
        np.asarray([0, 2], dtype=int),
        np.asarray([1, 3], dtype=int),
        np.asarray([2], dtype=int),
    ]
    module = TerrainModule()
    mask = np.asarray([True, True, True, False])

    width = module._mask_width_steps(grid, mask)
    width_again = module._mask_width_steps(grid, mask)

    assert width.tolist() == [3.0, 2.0, 1.0, 0.0]
    assert np.array_equal(width, width_again)
    assert len(module._mask_width_steps_cache) == 1


def test_truth_layers_are_cyclic_at_dateline():
    eng = _short(cells=1200, t_end=800.0)
    w = eng.world
    for name in ("tectonics.plate_id", "terrain.elevation_m", "crust.age_myr",
                 "ocean.basin_id", "ocean.depth_province"):
        raster = to_raster(w.grid, w.get_field(name), width=180, height=90)
        assert np.allclose(raster[:, 0], raster[:, -1])
    seam = validation._seam_continuity_metrics(w)
    assert seam["seam_edge_count"] > 0
    assert seam["render_duplicate_land_mismatch_fraction"] == pytest.approx(0.0)
    assert seam["render_duplicate_elevation_delta_p95_m"] == pytest.approx(0.0)


def test_artificial_antimeridian_mask_is_connected_on_sphere():
    grid = SphereGrid.fibonacci(1800, CONSTANTS.EARTH_RADIUS)
    seam_band = (np.abs(grid.lon) > 165.0) & (np.abs(grid.lat) < 28.0)
    assert len(validation._component_areas(grid, seam_band)) == 1


def test_ocean_geography_reports_large_closed_ocean_candidates():
    grid = SphereGrid.fibonacci(2400, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    land_band = np.abs(grid.lat) < 7.5
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[land_band] = 300.0
    world.set_field("terrain.elevation_m", surface)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("crust.type", land_band.astype(np.float64))
    world.set_field("ocean.basin_id", np.where(land_band, -1.0, 0.0))
    world.set_field("ocean.margin_type", np.where(land_band, 0.0, 1.0))
    world.set_field("ocean.depth_province", np.where(land_band, 0.0, 4.0))
    world.set_field("ocean.gateway_id", np.full(grid.n, -1.0))
    world.set_field("ocean.shelf_width", np.where(land_band, 0.0, 4.0))

    detail = validation._ocean_geography_metrics(world)

    assert detail["ocean_component_count"] == 2
    assert detail["closed_ocean_candidate_count_gt2pct_world"] == 1
    assert detail["max_closed_ocean_candidate_fraction_world"] > 0.20
    assert detail["closed_ocean_candidate_fraction_of_ocean"] > 0.30


def test_large_closed_ocean_repair_opens_seaway_to_main_ocean():
    grid = SphereGrid.fibonacci(2400, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    world.set_g("terrain.enable_closed_ocean_seaway_repair", 1.0)
    terrain = TerrainModule()

    land_band = np.abs(grid.lat) < 7.5
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[land_band] = 260.0
    crust_type = land_band.astype(np.float64)
    crust_domain = np.where(land_band, DOMAIN_CONTINENTAL_MARGIN, DOMAIN_OCEANIC)
    crust_stability = np.where(land_band, 0.35, 0.0)

    before = validation._component_areas(grid, surface < 0.0)
    repaired = terrain._connect_large_closed_oceans(
        world, surface, 0.0, crust_type, crust_domain, crust_stability)
    after = validation._component_areas(grid, repaired < 0.0)

    assert len(before) == 2
    assert len(after) == 1
    assert world.g("terrain.last_closed_ocean_seaway_openings") >= 1.0
    assert world.objects["terrain.closed_ocean_seaway_opened_corridors"]


def test_large_closed_ocean_repair_is_opt_in_by_default():
    grid = SphereGrid.fibonacci(1200, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    land_band = np.abs(grid.lat) < 7.5
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[land_band] = 260.0
    crust_type = land_band.astype(np.float64)
    crust_domain = np.where(land_band, DOMAIN_CONTINENTAL_MARGIN, DOMAIN_OCEANIC)
    crust_stability = np.where(land_band, 0.35, 0.0)

    repaired = terrain._connect_large_closed_oceans(
        world, surface, 0.0, crust_type, crust_domain, crust_stability)
    late = terrain._p111617_late_closed_ocean_repair(
        world, surface, 0.0, crust_type, crust_domain, crust_stability)

    assert np.array_equal(repaired, surface)
    assert np.array_equal(late, surface)
    assert world.g("terrain.last_closed_ocean_seaway_openings") == 0.0
    assert world.g("terrain.last_p111617_late_closed_ocean_attempted") == 0.0


def test_late_closed_ocean_repair_reports_improvement_and_preserves_corridors():
    grid = SphereGrid.fibonacci(2400, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    world.set_g("terrain.enable_closed_ocean_seaway_repair", 1.0)
    world.set_g("terrain.enable_p111617_late_closed_ocean_repair", 1.0)
    terrain = TerrainModule()

    land_band = np.abs(grid.lat) < 7.5
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[land_band] = 260.0
    crust_type = land_band.astype(np.float64)
    crust_domain = np.where(land_band, DOMAIN_CONTINENTAL_MARGIN, DOMAIN_OCEANIC)
    crust_stability = np.where(land_band, 0.35, 0.0)
    world.objects["terrain.closed_ocean_seaway_opened_corridors"] = [
        {"id": "prior_closed_ocean_seaway", "cells": [], "basis": "early_repair"}
    ]

    before = validation._component_areas(grid, surface < 0.0)
    repaired = terrain._p111617_late_closed_ocean_repair(
        world, surface, 0.0, crust_type, crust_domain, crust_stability)
    after = validation._component_areas(grid, repaired < 0.0)

    assert len(before) == 2
    assert len(after) == 1
    assert world.g("terrain.last_p111617_late_closed_ocean_attempted") == 1.0
    assert world.g("terrain.last_p111617_late_closed_ocean_applied") == 1.0
    assert world.g("terrain.last_p111617_late_closed_ocean_score_after") < (
        world.g("terrain.last_p111617_late_closed_ocean_score_before")
    )
    assert world.g("terrain.last_p111617_late_closed_ocean_openings") >= 1.0
    assert world.g("terrain.last_p111617_late_closed_ocean_area_fraction") > 0.0
    corridor_ids = [
        str(obj.get("id", ""))
        for obj in world.objects["terrain.closed_ocean_seaway_opened_corridors"]
    ]
    assert "prior_closed_ocean_seaway" in corridor_ids
    assert any(
        item.startswith("p111617_late_closed_ocean_seaway:")
        for item in corridor_ids
    )


def test_patchy_open_ocean_ridge_fragments_do_not_remain_shallow():
    grid = SphereGrid.fibonacci(1800, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    surface = np.full(grid.n, -900.0, dtype=np.float64)
    land = (grid.lon < -155.0) & (np.abs(grid.lat) < 35.0)
    surface[land] = 300.0
    ocean = surface < 0.0
    candidates = np.where(ocean & (np.abs(grid.lon) < 45.0) & (np.abs(grid.lat) < 55.0))[0]
    sparse_ridge = candidates[::70][:10]
    world.networks["tectonics.boundaries"] = {"ridge": sparse_ridge.astype(int).tolist()}

    repaired = terrain._deepen_modern_earthlike_open_ocean_shoals(world, surface, 0.0)

    assert sparse_ridge.size >= 4
    assert float(np.percentile(repaired[sparse_ridge], 50)) <= -3000.0


def test_patchy_open_ocean_trench_fragments_do_not_create_deep_pits():
    grid = SphereGrid.fibonacci(1800, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    surface = np.full(grid.n, -3200.0, dtype=np.float64)
    land = (grid.lon < -155.0) & (np.abs(grid.lat) < 35.0)
    surface[land] = 300.0
    ocean = surface < 0.0
    candidates = np.where(ocean & (np.abs(grid.lon) < 50.0) & (np.abs(grid.lat) < 55.0))[0]
    sparse_trench = candidates[::65][:10]
    world.networks["tectonics.boundaries"] = {"trench": sparse_trench.astype(int).tolist()}

    repaired = terrain._regionalize_ocean_floor(world, surface, 0.0)
    fields, _, _ = terrain._ocean_geography(world, repaired, 0.0)
    depth_province = fields["ocean.depth_province"].astype(int)

    assert sparse_trench.size >= 4
    assert float(np.percentile(repaired[sparse_trench], 50)) > -4100.0
    assert not np.isin(depth_province[sparse_trench], [6]).any()


def test_age_derived_open_ocean_ridge_restores_coherent_bathymetry():
    grid = SphereGrid.fibonacci(2400, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    surface = np.full(grid.n, -4300.0, dtype=np.float64)
    land = (grid.lon < -155.0) & (np.abs(grid.lat) < 40.0)
    surface[land] = 300.0
    ocean = surface < 0.0
    age = np.clip(np.abs(grid.lon) * 2.1, 0.0, 120.0)
    age[land] = 900.0
    world.set_field("crust.age_myr", age)

    repaired = terrain._regionalize_ocean_floor(world, surface, 0.0)
    fields, _, _ = terrain._ocean_geography(world, repaired, 0.0)
    depth_province = fields["ocean.depth_province"].astype(int)

    ridge_band = ocean & (np.abs(grid.lon) < 4.0) & (np.abs(grid.lat) < 55.0)
    old_flank = ocean & (np.abs(grid.lon) > 35.0) & (np.abs(grid.lon) < 70.0) & (
        np.abs(grid.lat) < 55.0
    )

    assert int(ridge_band.sum()) >= 8
    assert int(old_flank.sum()) >= 8
    assert float(np.median(repaired[ridge_band])) > float(np.median(repaired[old_flank])) + 700.0
    assert np.isin(depth_province[ridge_band], [5]).any()


def test_p1112_continental_hypsometry_retune_lowers_smooth_interiors_not_mountains():
    grid = SphereGrid.fibonacci(1800, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    spec.n_plates = 36
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    continent = (
        (np.abs(grid.lat) < 55.0)
        & (grid.lon > -130.0)
        & (grid.lon < 135.0)
    )
    basin = (
        continent
        & (grid.lon > -105.0)
        & (grid.lon < -35.0)
        & (np.abs(grid.lat) < 34.0)
    )
    active_range = continent & (grid.lon > 62.0) & (grid.lon < 92.0)
    plateau = (
        continent
        & (grid.lon > 96.0)
        & (grid.lon < 126.0)
        & (grid.lat > 4.0)
        & (grid.lat < 40.0)
    )
    platform = continent & ~(basin | active_range | plateau)

    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    surface[platform] = 920.0
    surface[basin] = 760.0
    surface[active_range] = 2850.0
    surface[plateau] = 3300.0
    crust_type = np.where(continent, 1.0, 0.0)
    detail = np.zeros(grid.n, dtype=np.float64)
    detail[platform] = CONT_DETAIL_PLATFORM
    detail[basin] = CONT_DETAIL_BASIN
    detail[active_range] = CONT_DETAIL_OROGEN
    detail[plateau] = CONT_DETAIL_PLATEAU
    inventory = np.zeros(grid.n, dtype=np.float64)
    inventory[active_range] = 1.0
    inventory[plateau] = 3.0
    world.set_field("crust.type", crust_type)
    world.set_field(
        "crust.domain",
        np.where(
            active_range,
            DOMAIN_CONTINENTAL_MARGIN,
            np.where(continent, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC),
        ),
    )
    world.set_field(
        "crust.stability",
        np.where(basin, 0.45, np.where(continent, 0.72, 0.0)),
    )
    world.set_field("terrain.continental_detail", detail)
    world.set_field("terrain.mountain_inventory", inventory)

    before_rel = surface.copy()
    before_land = surface >= 0.0
    before_low500 = float(
        grid.cell_area[before_land & (before_rel < 500.0)].sum()
        / grid.cell_area[before_land].sum()
    )
    before_mean = float(np.average(before_rel[before_land],
                                   weights=grid.cell_area[before_land]))
    before_high_median = float(np.median(before_rel[active_range | plateau]))

    repaired = terrain._p111_continental_hypsometry_finish_retune(
        world,
        surface,
        0.0,
        crust_type,
        detail,
        inventory,
    )
    after_land = repaired >= 0.0
    after_rel = repaired
    after_low500 = float(
        grid.cell_area[after_land & (after_rel < 500.0)].sum()
        / grid.cell_area[after_land].sum()
    )
    after_mean = float(np.average(after_rel[after_land],
                                  weights=grid.cell_area[after_land]))
    after_high_median = float(np.median(after_rel[active_range | plateau]))

    assert np.array_equal(before_land, after_land)
    assert after_low500 > before_low500 + 0.20
    assert after_mean < before_mean - 100.0
    assert after_high_median == pytest.approx(before_high_median)
    assert world.g("terrain.last_p1112_adjusted_area_fraction") > 0.04
    assert world.g("terrain.last_p1112_low500_fraction_after") > 0.35
    assert world.g("terrain.last_p1112_land_mask_preserved") == pytest.approx(1.0)
    assert abs(world.g("terrain.last_p1112_highland_preservation_delta_p50_m")) < 1.0


def test_p106_ocean_feature_metrics_cover_object_categories():
    grid = SphereGrid.fibonacci(1200, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    main_land = (grid.lon < -150.0) & (np.abs(grid.lat) < 35.0)
    ridge = np.where((np.abs(grid.lon) < 2.5) & (np.abs(grid.lat) < 45.0))[0]
    fracture = np.where((np.abs(grid.lat) < 2.5) & (np.abs(grid.lon) > 30.0))[0]
    seamount = np.where((np.abs(grid.lon - 80.0) < 4.0) & (np.abs(grid.lat) < 10.0))[0]
    plateau = np.where((np.abs(grid.lon + 85.0) < 7.0) & (np.abs(grid.lat + 10.0) < 8.0))[0]
    microcontinent = np.where((np.abs(grid.lon - 130.0) < 5.0) & (np.abs(grid.lat - 12.0) < 8.0))[0]
    parented_island = seamount[:1]
    surface[main_land] = 300.0
    surface[seamount] = -1200.0
    surface[parented_island] = 90.0
    surface[plateau] = -1700.0
    surface[microcontinent] = -900.0

    world.set_field("terrain.elevation_m", surface)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field(
        "crust.type",
        (main_land | np.isin(np.arange(grid.n), microcontinent)).astype(np.float64),
    )
    world.set_field("ocean.basin_id", np.zeros(grid.n, dtype=np.float64))
    world.set_field("ocean.margin_type", np.full(grid.n, 5.0))
    world.set_field("ocean.depth_province", np.full(grid.n, 4.0))
    world.set_field("ocean.gateway_id", np.full(grid.n, -1.0))
    world.set_field("ocean.shelf_width", np.full(grid.n, 4.0))
    world.objects["terrain.ocean_fabric"] = [
        {"kind": "spreading_center", "cells": ridge.astype(int).tolist(), "lon_span_deg": 5.0},
        {"kind": "transform_fault", "cells": fracture[:20].astype(int).tolist(), "lon_span_deg": 35.0},
        {"kind": "fracture_zone", "cells": fracture.astype(int).tolist(), "lon_span_deg": 120.0},
        {"kind": "abyssal_plain", "cells": np.where(surface < -3000.0)[0][:40].astype(int).tolist()},
        {"kind": "abyssal_hill", "cells": np.where(surface < -3000.0)[0][230:260].astype(int).tolist()},
        {"kind": "age_isochron", "cells": np.where(surface < -3000.0)[0][40:90].astype(int).tolist()},
    ]
    world.objects["terrain.arc_plume_landforms"] = [
        {"kind": "seamount_chain", "cells": seamount.astype(int).tolist(), "ocean_fraction": 1.0},
        {"kind": "hotspot_track", "cells": parented_island.astype(int).tolist(),
         "ocean_fraction": 0.0},
        {"kind": "oceanic_plateau", "cells": plateau.astype(int).tolist(), "ocean_fraction": 1.0},
        {"kind": "back_arc_basin", "cells": np.where(surface < -3000.0)[0][90:130].astype(int).tolist(),
         "ocean_fraction": 1.0},
        {"kind": "accreted_terrane", "cells": microcontinent.astype(int).tolist(),
         "ocean_fraction": 1.0},
    ]
    world.objects["terrain.margin_landforms"] = [
        {"kind": "passive_margin_wedge", "cells": np.where(surface < -3000.0)[0][130:170].astype(int).tolist()},
        {"kind": "delta_fan", "cells": np.where(surface < -3000.0)[0][170:200].astype(int).tolist()},
        {"kind": "trench", "cells": np.where(surface < -3000.0)[0][200:230].astype(int).tolist()},
    ]

    detail = validation._ocean_geography_metrics(world)

    assert detail["p106_core_ocean_fabric_complete"]
    assert detail["abyssal_hill_object_count"] >= 1
    assert detail["p106_object_category_count"] >= 10
    assert detail["seamount_chain_object_count"] >= 1
    assert detail["oceanic_plateau_object_count"] >= 1
    assert detail["back_arc_basin_object_count"] >= 1
    assert detail["passive_margin_wedge_object_count"] >= 1
    assert detail["parented_archipelago_component_count"] >= 1
    assert detail["unparented_far_ocean_shallow_component_count_lt1500m"] == 0
    assert detail["p106_full_ocean_feature_coverage"]


def test_p106_deep_pit_metric_ignores_parented_fracture_lineation():
    grid = SphereGrid.fibonacci(1200, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)

    parented = int(np.argmin((grid.lon + 40.0) ** 2 + grid.lat ** 2))
    unparented = int(np.argmin((grid.lon - 55.0) ** 2 + (grid.lat - 8.0) ** 2))
    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    surface[parented] = -4300.0
    surface[unparented] = -4300.0

    world.set_field("terrain.elevation_m", surface)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("crust.type", np.zeros(grid.n, dtype=np.float64))
    world.set_field("ocean.basin_id", np.zeros(grid.n, dtype=np.float64))
    world.set_field("ocean.margin_type", np.full(grid.n, 5.0))
    world.set_field("ocean.depth_province", np.full(grid.n, 4.0))
    world.set_field("ocean.gateway_id", np.full(grid.n, -1.0))
    world.set_field("ocean.shelf_width", np.full(grid.n, 4.0))
    world.objects["terrain.ocean_fabric"] = [
        {"kind": "fracture_zone", "cells": [parented]},
    ]
    world.objects["terrain.arc_plume_landforms"] = []
    world.objects["terrain.margin_landforms"] = []

    detail = validation._ocean_geography_metrics(world)

    assert parented != unparented
    assert detail["isolated_far_ocean_deep_pit_component_count_gt4100m"] == 1


def test_age_discontinuity_fallback_generates_fracture_zone_objects():
    grid = SphereGrid.fibonacci(2200, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    surface = np.full(grid.n, -4400.0, dtype=np.float64)
    ridge = (np.abs(grid.lon) < 3.0) & (np.abs(grid.lat) < 62.0)
    fracture_trace = (np.abs(grid.lat) < 3.5) & (np.abs(grid.lon) > 18.0) & (
        np.abs(grid.lon) < 105.0
    )
    surface[ridge] = -2600.0
    age = np.clip(np.abs(grid.lon) * 2.0, 0.0, 160.0)
    age[(grid.lat > 0.0) & (np.abs(grid.lon) > 18.0)] += 55.0
    age[ridge] = 0.0
    world.set_field("crust.age_myr", age)
    world.set_field("crust.type", np.zeros(grid.n, dtype=np.float64))
    world.objects["tectonics.spreading_centers"] = [
        {"id": "ridge:p106-age", "kind": "ridge", "cells": np.where(ridge)[0].astype(int).tolist()}
    ]

    fields, _, _ = terrain._ocean_geography(world, surface, 0.0)
    objects = terrain._ocean_fabric_objects(
        world, surface, 0.0, world.get_field("crust.type", 0.0),
        np.full(grid.n, 120.0), fields)
    fracture_objects = [obj for obj in objects if obj.get("kind") == "fracture_zone"]
    fracture_cells = np.asarray(
        [cell for obj in fracture_objects for cell in obj.get("cells", [])],
        dtype=int,
    )

    assert fracture_trace.sum() >= 10
    assert len(fracture_objects) >= 1
    assert fracture_cells.size >= 6
    assert float(np.ptp(grid.lon[fracture_cells])) > 45.0


def test_sparse_transform_tail_fallback_generates_fracture_zone_objects():
    grid = SphereGrid.fibonacci(1600, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    surface = np.full(grid.n, -4200.0, dtype=np.float64)
    ridge = (np.abs(grid.lon) < 3.0) & (np.abs(grid.lat) < 55.0)
    transform = (
        (np.abs(grid.lat) < 4.0)
        & (grid.lon > 5.0)
        & (grid.lon < 20.0)
    )
    surface[ridge] = -2600.0
    age = np.clip(np.abs(grid.lon) * 2.1, 0.0, 150.0)
    age[ridge] = 0.0
    world.set_field("crust.age_myr", age)
    world.set_field("crust.type", np.zeros(grid.n, dtype=np.float64))
    world.networks["tectonics.boundaries"] = {
        "ridge": np.where(ridge)[0].astype(int).tolist(),
        "transform": np.where(transform)[0].astype(int).tolist(),
    }

    fields, _, _ = terrain._ocean_geography(world, surface, 0.0)
    objects = terrain._ocean_fabric_objects(
        world, surface, 0.0, world.get_field("crust.type", 0.0),
        np.full(grid.n, 120.0), fields)
    fracture_objects = [obj for obj in objects if obj.get("kind") == "fracture_zone"]
    transform_objects = [obj for obj in objects if obj.get("kind") == "transform_fault"]

    assert len(transform_objects) >= 1
    assert len(fracture_objects) >= 1
    assert sum(int(obj.get("cell_count", 0)) for obj in fracture_objects) >= 2


def test_arc_plume_landforms_promote_p106_oceanic_objects():
    grid = SphereGrid.fibonacci(1800, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    surface = np.full(grid.n, -3400.0, dtype=np.float64)
    crust_type = np.zeros(grid.n, dtype=np.float64)
    crust_domain = np.full(grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    crust_origin = np.zeros(grid.n, dtype=np.float64)
    crust_thick = np.full(grid.n, 7000.0, dtype=np.float64)
    sediment = np.full(grid.n, 300.0, dtype=np.float64)
    volcanism_age = np.full(grid.n, -1.0, dtype=np.float64)
    terrane_id = np.full(grid.n, -1.0, dtype=np.float64)
    continent_id = np.full(grid.n, -1.0, dtype=np.float64)

    chain = (np.abs(grid.lon + 35.0) < 4.5) & (np.abs(grid.lat) < 35.0)
    plateau = (np.abs(grid.lon - 35.0) < 8.0) & (np.abs(grid.lat + 15.0) < 8.0)
    microcontinent = (np.abs(grid.lon - 95.0) < 6.0) & (np.abs(grid.lat - 15.0) < 8.0)
    surface[chain] = -2300.0
    surface[plateau] = -1800.0
    surface[microcontinent] = -1200.0
    volcanism_age[chain] = 4460.0
    crust_domain[plateau] = DOMAIN_LIP
    crust_origin[plateau] = ORIGIN_PLUME_IMPACT
    crust_type[microcontinent] = 1.0
    crust_domain[microcontinent] = DOMAIN_ACCRETED_TERRANE
    terrane_id[microcontinent] = 2.0
    continent_id[microcontinent] = 3.0
    world.set_field("tectonics.volcanism_age_myr", volcanism_age)
    world.set_field("tectonics.terrane_id", terrane_id)
    world.set_field("tectonics.continent_id", continent_id)
    world.objects["tectonics.plumes"] = [
        {"id": "plume:p106", "cell": int(np.where(chain)[0][0]), "cells": np.where(chain)[0].astype(int).tolist()}
    ]
    world.objects["tectonics.lips"] = [
        {"id": "lip:p106", "kind": "large_igneous_province",
         "cells": np.where(plateau)[0].astype(int).tolist()}
    ]
    world.objects["tectonics.volcanoes"] = [
        {"id": "volcano:p106", "cell": int(np.where(chain)[0][0]),
         "cells": np.where(chain)[0].astype(int).tolist()}
    ]
    ocean_geo_fields = {
        "ocean.basin_id": np.zeros(grid.n, dtype=np.float64),
        "ocean.margin_type": np.full(grid.n, 5.0),
        "ocean.depth_province": np.full(grid.n, 4.0),
        "ocean.shelf_width": np.full(grid.n, 4.0),
    }

    objects = terrain._arc_plume_landform_objects(
        world, surface, 0.0, crust_type, crust_domain, crust_origin,
        crust_thick, sediment, ocean_geo_fields)
    kinds = [obj.get("kind") for obj in objects]

    assert "seamount_chain" in kinds
    assert "oceanic_plateau" in kinds
    assert "microcontinent" in kinds


def test_parented_hotspot_chain_can_form_limited_archipelago():
    grid = SphereGrid.fibonacci(2200, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    spec.n_plates = 24
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    surface = np.full(grid.n, -3600.0, dtype=np.float64)
    land = (grid.lon < -165.0) & (np.abs(grid.lat) < 45.0)
    surface[land] = 300.0
    chain = (
        (np.abs(grid.lat) < 5.0)
        & (grid.lon > -60.0)
        & (grid.lon < 80.0)
        & ~land
    )
    surface[chain] = -1550.0
    volcanism_age = np.full(grid.n, -1.0, dtype=np.float64)
    volcanism_age[chain] = 4475.0
    world.set_field("tectonics.volcanism_age_myr", volcanism_age)
    world.set_field("crust.type", np.zeros(grid.n, dtype=np.float64))
    world.set_field("crust.domain", np.full(grid.n, DOMAIN_OCEANIC, dtype=np.float64))
    world.set_field("crust.origin", np.zeros(grid.n, dtype=np.float64))
    world.objects["tectonics.plumes"] = [
        {"id": "plume:archipelago", "cell": int(np.where(chain)[0][0]),
         "cells": np.where(chain)[0].astype(int).tolist()}
    ]
    world.objects["tectonics.volcanoes"] = [
        {"id": "volcano:archipelago", "cell": int(np.where(chain)[0][0]),
         "cells": np.where(chain)[0].astype(int).tolist()}
    ]

    repaired = terrain._regionalize_ocean_floor(world, surface, 0.0)
    emergent = chain & (repaired >= 0.0)
    world.set_field("terrain.elevation_m", repaired)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("crust.thickness_m", np.full(grid.n, 7000.0, dtype=np.float64))
    world.set_field("sediment.thickness_m", np.full(grid.n, 300.0, dtype=np.float64))
    ocean_geo_fields, _, _ = terrain._ocean_geography(world, repaired, 0.0)
    for name, values in ocean_geo_fields.items():
        world.set_field(name, values)
    world.objects["terrain.arc_plume_landforms"] = terrain._arc_plume_landform_objects(
        world,
        repaired,
        0.0,
        world.get_field("crust.type", 0.0),
        world.get_field("crust.domain", 0.0),
        world.get_field("crust.origin", 0.0),
        world.get_field("crust.thickness_m", 7000.0),
        world.get_field("sediment.thickness_m", 300.0),
        ocean_geo_fields,
    )
    world.objects["terrain.ocean_fabric"] = []
    world.objects["terrain.margin_landforms"] = []
    detail = validation._ocean_geography_metrics(world)

    assert int(chain.sum()) >= 12
    assert int(emergent.sum()) >= 1
    assert int(emergent.sum()) < int(chain.sum()) * 0.35
    assert float(np.percentile(repaired[chain & ~emergent], 50)) < -700.0
    assert detail["parented_archipelago_component_count"] >= 1
    assert detail["unparented_oceanic_island_component_count"] == 0


def test_final_trim_preserves_parented_archipelago_and_drowns_unparented_island():
    grid = SphereGrid.fibonacci(1600, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    terrain = TerrainModule()

    surface = np.full(grid.n, -2600.0, dtype=np.float64)
    main_land = (grid.lon < -150.0) & (np.abs(grid.lat) < 45.0)
    parented_island = np.where(
        (np.abs(grid.lon - 35.0) < 2.5) & (np.abs(grid.lat + 12.0) < 2.5)
    )[0][:1]
    unparented_island = np.where(
        (np.abs(grid.lon + 20.0) < 2.5) & (np.abs(grid.lat - 18.0) < 2.5)
    )[0][:1]
    surface[main_land] = 250.0
    surface[parented_island] = 80.0
    surface[unparented_island] = 80.0

    crust_type = np.zeros(grid.n, dtype=np.float64)
    crust_domain = np.full(grid.n, DOMAIN_OCEANIC, dtype=np.float64)
    crust_origin = np.zeros(grid.n, dtype=np.float64)
    crust_stability = np.zeros(grid.n, dtype=np.float64)
    volcanism_age = np.full(grid.n, -1.0, dtype=np.float64)
    volcanism_age[parented_island] = 4475.0
    world.set_field("tectonics.volcanism_age_myr", volcanism_age)

    trimmed = terrain._trim_final_unparented_oceanic_islands(
        world, surface, 0.0, crust_type, crust_domain, crust_origin, crust_stability)

    assert parented_island.size == 1
    assert unparented_island.size == 1
    assert trimmed[parented_island[0]] >= 0.0
    assert trimmed[unparented_island[0]] < 0.0


def test_stellar_flux_closes_orbital_energy_budget():
    eng = _short(cells=1200, t_end=500.0)
    w = eng.world
    flux = w.get_field("stellar.surface_flux")
    mean_flux = float(np.average(flux, weights=w.grid.cell_area))
    ecc_factor = 1.0 / np.sqrt(1.0 - w.g("orbit.eccentricity") ** 2)
    d = w.spec.orbit.semi_major_axis_au * CONSTANTS.AU
    expected = w.g("stellar.luminosity") / (16.0 * np.pi * d * d) * ecc_factor
    assert mean_flux == pytest.approx(expected, rel=1e-6)
    assert w.g("stellar.mean_flux") == pytest.approx(mean_flux, rel=1e-9)


def test_interior_heat_budget_terms_are_physical():
    eng = _short(cells=1200, t_end=500.0)
    w = eng.world
    module = InteriorModule()
    assert module._radiogenic(w, 0.0) > module._radiogenic(w, w.spec.t_end_myr)
    assert w.g("interior.radiogenic_heat") > 0.0
    assert w.g("interior.convective_heat_loss") > 0.0
    assert w.g("interior.surface_heat_flow") == pytest.approx(
        w.g("interior.convective_heat_loss") / w.spec.surface_area_m2,
        rel=1e-12,
    )
    assert w.g("interior.tectonic_vigor") >= 0.0


def test_interior_r1_potential_fields_are_bounded():
    eng = _short(cells=1200, t_end=800.0)
    w = eng.world

    bounded = {
        "mantle.heat_anomaly": (-1.0, 1.0),
        "mantle.upwelling_potential": (0.0, 1.0),
        "mantle.downwelling_potential": (0.0, 1.0),
        "tectonics.rift_potential": (0.0, 1.0),
        "tectonics.plume_potential": (0.0, 1.0),
    }
    for name, (lo, hi) in bounded.items():
        values = w.get_field(name)
        assert values.shape == (w.n_cells,)
        assert np.isfinite(values).all()
        assert float(values.min()) >= lo - 1e-9
        assert float(values.max()) <= hi + 1e-9

    lithosphere = w.get_field("lithosphere.thermal_thickness")
    assert lithosphere.shape == (w.n_cells,)
    assert np.isfinite(lithosphere).all()
    assert float(lithosphere.min()) >= 20000.0
    assert float(lithosphere.max()) <= 260000.0
    assert float(w.get_field("tectonics.plume_potential").max()) > 0.0


def test_r2_plate_motion_refresh_records_force_components():
    eng = _short(cells=1200, t_end=800.0)
    plates = eng.world.objects.get("tectonics.plates", [])
    active = set(eng.world.get_field("tectonics.plate_id").astype(int))
    refreshed = [
        plate for plate in plates
        if int(plate.get("id", -1)) in active and "motion_source" in plate
    ]

    assert refreshed
    sources = {plate["motion_source"] for plate in refreshed}
    assert sources <= {
        "r2_torque_proxy",
        "fallback_previous_motion",
        "r2_collision_locked",
        "p20_topology_split_inherited",
    }
    for plate in refreshed:
        components = plate.get("r2_force_components", {})
        assert set(components) >= {
            "slab_pull",
            "ridge_push",
            "collision_resistance",
            "basal_drag",
            "transform_friction",
            "torque_norm",
            "rate",
        }
        assert components["rate"] == pytest.approx(float(plate["rate"]), rel=1e-12)


def test_r3_boundary_objects_expose_persistence_and_polarity_metadata():
    eng = _short(cells=1600, t_end=1200.0)
    objects = eng.world.objects.get("tectonics.boundary_objects", [])

    assert objects
    assert any("age_myr" in obj for obj in objects)
    assert all("parent_plate_ids" in obj for obj in objects)
    assert all("boundary_continental_fraction" in obj for obj in objects)
    assert {obj.get("persistence") for obj in objects} <= {"new", "matched_overlap"}

    subduction = [
        obj for obj in objects
        if obj.get("kind") in {"trench", "active_margin"}
        and obj.get("polarity_basis") != "continental_collision_no_subduction_polarity"
    ]
    if subduction:
        assert any(obj.get("subducting_plate_id") is not None for obj in subduction)
        assert any(obj.get("overriding_plate_id") is not None for obj in subduction)


def test_p20_plate_topologies_expose_resolved_network_metadata():
    eng = _short(cells=1400, t_end=900.0)
    w = eng.world
    topologies = w.objects.get("tectonics.plate_topologies", [])
    active = set(w.get_field("tectonics.plate_id").astype(int))

    assert topologies
    assert {obj["numeric_id"] for obj in topologies} == active
    assert any(obj["neighbour_plate_ids"] for obj in topologies)
    assert any(
        max(obj["boundary_area_fractions"].values()) > 0.0
        for obj in topologies
    )
    for obj in topologies:
        assert obj["kind"] == "resolved_plate_topology"
        assert obj["component_count"] >= 1
        assert 0.0 < obj["area_fraction"] <= 1.0
        assert 0.0 <= obj["largest_component_fraction"] <= 1.0
        assert obj["continental_fraction"] + obj["oceanic_fraction"] == pytest.approx(1.0)
        assert set(obj["boundary_area_fractions"]) >= {
            "ridge",
            "trench",
            "collision",
            "suture",
            "active_margin",
            "passive_margin",
            "transform",
        }
        assert obj["motion_source"] in {
            "initial_or_unresolved",
            "r2_torque_proxy",
            "fallback_previous_motion",
            "r2_collision_locked",
            "p20_topology_split_inherited",
        }


def test_r4_wilson_lifecycle_objects_expose_parent_causality():
    eng = _short(cells=1800, t_end=2500.0)
    w = eng.world

    basins = w.objects.get("tectonics.ocean_basins", [])
    cycles = w.objects.get("tectonics.wilson_cycles", [])
    gateways = w.objects.get("tectonics.ocean_gateways", [])
    assert basins
    assert cycles
    assert gateways

    basin_ids = {obj["id"] for obj in basins}
    assert all(obj.get("lineage_key") for obj in basins)
    assert all(obj.get("boundary_object_ids") for obj in basins)
    assert all(obj.get("age_myr", -1.0) >= 0.0 for obj in basins)
    assert any(
        w.objects.get(name)
        for name in (
            "tectonics.rift_systems",
            "tectonics.passive_margins",
            "tectonics.spreading_centers",
            "tectonics.closing_margins",
            "tectonics.sutures",
        )
    )

    for cycle in cycles:
        assert cycle.get("ocean_basin_id") in basin_ids
        assert cycle.get("boundary_object_ids")
    for gateway in gateways:
        assert gateway.get("parent_basin_id") in basin_ids
        assert gateway.get("parent_boundary_object_id")
        assert (
            gateway.get("parent_rift_system_ids")
            or gateway.get("parent_margin_ids")
            or gateway.get("parent_closing_margin_ids")
            or gateway.get("parent_suture_ids")
        )


def test_tectonic_diagnostics_cover_current_world_and_archive():
    eng = _short(cells=1800, t_end=2500.0)
    result = validation.check_tectonic_diagnostics(eng)
    assert result.passed

    detail = result.detail
    assert set(detail) >= {
        "plate_fragmentation",
        "crust_distribution",
        "boundaries",
        "seafloor_age_geometry",
        "ocean_geography",
        "seam_continuity",
        "frame_continuity",
        "morphology",
        "warnings",
        "hard_failures",
    }
    assert detail["plate_fragmentation"]["n_active_plates"] >= 1
    assert detail["plate_fragmentation"]["total_plate_components"] >= 1
    assert 0.0 <= detail["crust_distribution"]["continental_area_fraction"] <= 1.0
    assert detail["crust_distribution"]["oceanic_thickness_p50_m"] >= 0.0
    assert detail["crust_distribution"]["oceanic_age_p95_myr"] < 320.0
    assert detail["crust_distribution"]["oceanic_old_fraction_gt300_myr"] < 0.05
    assert detail["crust_distribution"]["invalid_origin_cells"] == 0
    assert detail["crust_distribution"]["stability_out_of_range_cells"] == 0
    assert detail["crust_distribution"]["continental_ridge_origin_cells"] == 0
    assert 0.0 <= detail["crust_distribution"]["stable_craton_fraction_gt075"] <= 1.0
    assert set(detail["boundaries"]) >= {
        "ridge_cells",
        "trench_cells",
        "suture_cells",
        "active_margin_cells",
        "passive_margin_cells",
        "transform_cells",
    }
    assert set(detail["seafloor_age_geometry"]) >= {
        "ridge_ocean_cells",
        "ridge_ocean_age_p75_myr",
        "ocean_with_ridge_path_fraction",
        "age_distance_correlation",
        "old_ocean_near_trench_fraction_gt180_myr",
    }
    assert detail["ocean_geography"]["has_basin_id"]
    assert detail["ocean_geography"]["has_margin_type"]
    assert detail["ocean_geography"]["has_depth_province"]
    assert detail["ocean_geography"]["has_gateway_id"]
    assert detail["ocean_geography"]["has_shelf_width"]
    assert detail["ocean_geography"]["land_with_depth_province_cells"] == 0
    assert detail["ocean_geography"]["ocean_without_depth_province_cells"] == 0
    assert detail["ocean_geography"]["ocean_basin_object_count"] >= 1
    assert detail["seam_continuity"]["seam_edge_count"] > 0
    assert detail["seam_continuity"]["seam_exposed_land_component_mismatch_fraction"] == 0.0
    assert detail["seam_continuity"]["seam_ocean_basin_mismatch_fraction"] == 0.0
    assert set(detail["morphology"]) >= {
        "context",
        "exposed_land",
        "continental_crust",
        "exposed_continental_land",
        "crust_land_coupling",
        "warnings",
    }
    assert detail["morphology"]["exposed_land"]["component_count"] >= 0
    assert 0.0 <= detail["morphology"]["exposed_land"]["ribbon_area_fraction_gt_0_5"] <= 1.0
    assert 0.0 <= detail["morphology"]["continental_crust"][
        "largest_component_area_fraction_of_mask"
    ] <= 1.0
    assert detail["morphology"]["crust_land_coupling"]["high_oceanic_land_cells_gt1500m"] >= 0
    if detail["seafloor_age_geometry"]["ridge_ocean_cells"] > 0:
        assert detail["seafloor_age_geometry"]["ridge_ocean_age_p75_myr"] <= 30.0
    assert detail["frame_continuity"]["n_frame_pairs"] >= 1
    assert detail["frame_continuity"]["archive_max_components_per_plate"] >= 1
    assert detail["frame_continuity"]["late_archive_max_components_per_plate"] >= 1
    assert 0.0 <= detail["frame_continuity"]["archive_min_largest_component_share"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["late_archive_min_largest_component_share"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["max_crust_type_change_fraction"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["active_plate_id_jaccard_min"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["max_plate_area_distribution_delta"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["late_active_plate_id_jaccard_min"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["late_max_plate_area_distribution_delta"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["max_continental_fraction_delta"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["late_max_continental_fraction_delta"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["max_plate_crust_composition_delta"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["late_max_plate_crust_composition_delta"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["max_plate_crust_composition_delta_per_100myr"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["late_max_plate_crust_composition_delta_per_100myr"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["continent_id_jaccard_min"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["max_crust_domain_change_fraction"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["max_terrain_province_change_fraction"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["max_wilson_phase_change_fraction"] <= 1.0
    assert 0.0 <= detail["frame_continuity"]["max_exposed_land_change_fraction"] <= 1.0


def test_climate_diagnostics_cover_current_world_and_detect_artifacts():
    eng = _short(cells=1200, t_end=1200.0)
    result = validation.check_climate_diagnostics(eng)
    assert result.passed

    detail = result.detail
    assert set(detail) >= {
        "context",
        "temperature",
        "precipitation",
        "coasts",
        "seasonality",
        "ocean_currents",
        "geography",
        "warnings",
        "hard_failures",
    }
    assert detail["temperature"]["has_surface_temperature"]
    assert detail["temperature"]["nonfinite_temperature_cells"] == 0
    assert detail["temperature"]["max_adjacent_lat_band_delta_C"] >= 0.0
    assert detail["temperature"]["neighbor_temperature_delta_p99_C"] >= 0.0
    assert detail["precipitation"]["has_precipitation"]
    assert detail["precipitation"]["negative_precipitation_cells"] == 0
    assert detail["precipitation"]["precip_orographic_concentration"] >= 0.0
    assert detail["coasts"]["coastal_temperature_asymmetry_index"] >= 0.0
    assert detail["seasonality"]["has_seasonal_temperature"]
    assert not detail["seasonality"]["invalid_seasonal_temperature_shape"]
    assert detail["seasonality"]["has_seasonal_precipitation"]
    assert not detail["seasonality"]["invalid_seasonal_precipitation_shape"]
    assert detail["seasonality"]["annual_precip_matches_seasonal_aggregate"]
    assert detail["circulation"]["has_seasonal_wind"]
    assert not detail["circulation"]["invalid_seasonal_wind_shape"]
    assert detail["circulation"]["has_seasonal_pressure_proxy"]
    assert detail["circulation"]["has_moisture_access"]
    assert detail["circulation"]["has_monsoon_potential"]
    assert detail["circulation"]["itcz_DJF_lat_deg"] < 0.0
    assert detail["circulation"]["itcz_JJA_lat_deg"] > 0.0
    assert detail["ocean_currents"]["has_currents"]
    assert detail["ocean_currents"]["has_current_heat_transport"]
    assert detail["ocean_currents"]["current_over_solved_land_cells_gt_1e_8"] == 0
    assert detail["geography"]["has_continent_id"]
    assert detail["geography"]["has_basin_id"]
    assert detail["geography"]["has_coast_orientation"]
    assert "seasonal precipitation field is absent" not in detail["warnings"]

    w = eng.world
    w.fields["climate.surface_temperature"] = np.where(
        w.grid.lat > 0.0, 305.0, 245.0,
    )
    relief = np.zeros(w.grid.n)
    elev = w.get_field("terrain.elevation_m") - w.sea_level
    for c, nbs in enumerate(w.grid.neighbors):
        if nbs.size:
            relief[c] = float(np.max(np.abs(elev[nbs] - elev[c])))
    land = w.land_mask()
    rough_land = land & (relief >= np.percentile(relief[land], 85))
    w.fields["climate.precipitation"] = np.where(rough_land, 2500.0, 120.0)

    artifact_detail = validation.climate_diagnostics(eng)
    assert artifact_detail["temperature"]["max_adjacent_lat_band_delta_C"] > 30.0
    assert artifact_detail["precipitation"]["precip_orographic_concentration"] > 2.0
    assert any("latitude-band jump" in w for w in artifact_detail["warnings"])
    assert any("high-relief cells" in w for w in artifact_detail["warnings"])


def test_seasonal_temperature_tracks_obliquity_and_thermal_inertia():
    eng = _short(cells=1400, t_end=1200.0)
    w = eng.world
    seasonal = w.get_field("climate.seasonal_temperature")
    annual = w.get_field("climate.surface_temperature")
    seasonality = w.get_field("climate.temperature_seasonality")
    continentality = w.get_field("climate.continentality")
    land = w.land_mask()
    ocean = ~land
    nh_land = land & (w.grid.lat > 20.0)

    assert seasonal.shape == (4, w.grid.n)
    assert int(nh_land.sum()) > 5
    assert np.allclose(annual, seasonal.mean(axis=0))
    assert np.allclose(seasonality, seasonal.max(axis=0) - seasonal.min(axis=0))
    assert np.all((0.0 <= continentality) & (continentality <= 1.0))
    assert float(np.average(seasonal[2, nh_land], weights=w.grid.cell_area[nh_land])) > float(
        np.average(seasonal[0, nh_land], weights=w.grid.cell_area[nh_land])
    ) + 2.0
    assert float(np.percentile(seasonality[land], 50)) > float(
        np.percentile(seasonality[ocean], 50)
    )

    detail = validation.climate_diagnostics(eng)
    assert detail["seasonality"]["land_seasonal_temp_amplitude_p50_C"] > (
        detail["seasonality"]["ocean_seasonal_temp_amplitude_p50_C"]
    )


def test_tidally_locked_world_does_not_get_earthlike_seasons():
    eng = _short("tidally_locked", cells=900, t_end=700.0)
    w = eng.world
    seasonal = w.get_field("climate.seasonal_temperature")
    seasonality = w.get_field("climate.temperature_seasonality")
    seasonal_wind = w.fields["atmosphere.seasonal_wind"]
    itcz = w.fields["atmosphere.itcz_latitude"]

    assert seasonal.shape == (4, w.grid.n)
    assert np.allclose(seasonal, seasonal[0][None, :])
    assert float(seasonality.max()) == pytest.approx(0.0)
    assert seasonal_wind.shape == (4, w.grid.n, 3)
    assert np.allclose(seasonal_wind, seasonal_wind[0][None, :, :])
    assert np.allclose(itcz, 0.0)


def test_seasonal_winds_migrate_itcz_and_storm_tracks():
    eng = _short(cells=1400, t_end=1200.0)
    w = eng.world
    seasonal_wind = w.fields["atmosphere.seasonal_wind"]
    annual_wind = w.get_field("atmosphere.wind")
    itcz = w.fields["atmosphere.itcz_latitude"]
    itcz_intensity = w.fields["atmosphere.itcz_intensity"]
    storm = w.fields["atmosphere.storm_track_intensity"]

    assert seasonal_wind.shape == (4, w.grid.n, 3)
    assert itcz.shape == (4,)
    assert itcz_intensity.shape == (4, w.grid.n)
    assert storm.shape == (4, w.grid.n)
    assert np.allclose(annual_wind, seasonal_wind.mean(axis=0))
    normal_component = np.abs(np.einsum("snj,nj->sn", seasonal_wind, w.grid.xyz))
    assert float(normal_component.max()) < 1e-8

    assert itcz[0] < -10.0
    assert abs(float(itcz[1])) < 1.0
    assert itcz[2] > 10.0
    assert abs(float(itcz[3])) < 1.0

    def weighted_lat(field):
        weights = np.maximum(field, 0.0) * w.grid.cell_area
        return float(np.sum(w.grid.lat * weights) / np.sum(weights))

    assert weighted_lat(itcz_intensity[0]) < -8.0
    assert weighted_lat(itcz_intensity[2]) > 8.0

    nh_mid = (w.grid.lat >= 35.0) & (w.grid.lat <= 60.0)
    sh_mid = (w.grid.lat <= -35.0) & (w.grid.lat >= -60.0)
    nh_djf = float(np.average(storm[0, nh_mid], weights=w.grid.cell_area[nh_mid]))
    nh_jja = float(np.average(storm[2, nh_mid], weights=w.grid.cell_area[nh_mid]))
    sh_djf = float(np.average(storm[0, sh_mid], weights=w.grid.cell_area[sh_mid]))
    sh_jja = float(np.average(storm[2, sh_mid], weights=w.grid.cell_area[sh_mid]))
    assert nh_djf > nh_jja * 1.10
    assert sh_jja > sh_djf * 1.10

    detail = validation.climate_diagnostics(eng)
    assert detail["circulation"]["NH_winter_storm_track_ratio"] > 1.10
    assert detail["circulation"]["SH_winter_storm_track_ratio"] > 1.10


def _coastal_direction_to_ocean(world):
    grid = world.grid
    land = world.land_mask()
    ocean = ~land
    out = np.zeros((grid.n, 3))
    coastal = np.zeros(grid.n, dtype=bool)
    for c in np.where(land)[0]:
        nbs = grid.neighbors[int(c)]
        targets = nbs[ocean[nbs]]
        if targets.size == 0:
            continue
        vec = np.mean(grid.xyz[targets], axis=0)
        tangent = vec - float(vec @ grid.xyz[c]) * grid.xyz[c]
        norm = float(np.linalg.norm(tangent))
        if norm <= 1e-12:
            continue
        out[c] = tangent / norm
        coastal[c] = True
    return out, coastal


def test_geographic_circulation_anomalies_follow_land_sea_layout():
    eng = _short(cells=1400, t_end=1200.0)
    w = eng.world
    final_wind = w.fields["atmosphere.seasonal_wind"]
    background = w.fields["atmosphere.background_seasonal_wind"]
    pressure = w.fields["atmosphere.land_sea_pressure_proxy"]
    pressure_center_support = w.fields["atmosphere.pressure_center_support"]
    pressure_center_id = w.fields["atmosphere.pressure_center_id"]
    stationary_wave_support = w.fields["atmosphere.stationary_wave_pressure_support"]
    pressure_genesis_source = w.fields["atmosphere.pressure_genesis_source"]
    pressure_genesis_transfer = w.fields["atmosphere.pressure_genesis_wave_transfer"]
    ocean_low_source = w.fields["atmosphere.ocean_pressure_low_source_support"]
    thermal = w.fields["atmosphere.thermal_wind_anomaly"]
    orographic = w.fields["atmosphere.orographic_wind_anomaly"]
    geo_index = w.fields["atmosphere.geographic_circulation_index"]

    assert background.shape == (4, w.grid.n, 3)
    assert pressure.shape == (4, w.grid.n)
    assert pressure_center_support.shape == (4, w.grid.n)
    assert pressure_center_id.shape == (4, w.grid.n)
    assert stationary_wave_support.shape == (4, w.grid.n)
    assert pressure_genesis_source.shape == (4, w.grid.n)
    assert pressure_genesis_transfer.shape == (4, w.grid.n)
    assert ocean_low_source.shape == (4, w.grid.n)
    assert thermal.shape == (4, w.grid.n, 3)
    assert orographic.shape == (4, w.grid.n, 3)
    assert geo_index.shape == (w.grid.n,)
    assert np.allclose(final_wind.mean(axis=0), w.get_field("atmosphere.wind"))
    assert float(np.percentile(np.linalg.norm(final_wind - background, axis=2), 95)) > 0.5
    assert float(np.percentile(geo_index, 90)) > 0.10
    assert float(np.percentile(pressure_center_support, 95)) > 0.10
    assert float(np.percentile(stationary_wave_support, 95)) > 0.05
    assert np.all(np.isfinite(pressure_genesis_source))
    assert np.all(np.isfinite(pressure_genesis_transfer))
    assert np.all(np.isfinite(ocean_low_source))
    pressure_centers = w.objects.get("atmosphere.pressure_centers", [])
    assert isinstance(pressure_centers, list)
    assert pressure_centers
    pressure_center_kinds = {str(obj.get("kind")) for obj in pressure_centers}
    assert {"pressure_low", "pressure_high"} & pressure_center_kinds

    land = w.land_mask()
    temp_anom = (
        w.fields["climate.seasonal_temperature"]
        - w.fields["climate.seasonal_temperature"].mean(axis=0, keepdims=True)
    )
    corr = np.corrcoef(temp_anom[:, land].ravel(), pressure[:, land].ravel())[0, 1]
    assert corr < -0.55

    to_ocean, coastal_land = _coastal_direction_to_ocean(w)
    nh_coast = coastal_land & (w.grid.lat > 15.0)
    jja_low = nh_coast & (pressure[2] < -0.05)
    djf_high = nh_coast & (pressure[0] > 0.05)
    assert int(jja_low.sum()) > 5
    assert int(djf_high.sum()) > 5
    jja_onshore = np.einsum("ij,ij->i", thermal[2], -to_ocean)
    djf_offshore = np.einsum("ij,ij->i", thermal[0], to_ocean)
    assert float(np.mean(jja_onshore[jja_low])) > 0.15
    assert float(np.mean(djf_offshore[djf_high])) > 0.15


def test_geographic_circulation_is_weak_on_waterworld():
    water = _short("waterworld", cells=900, t_end=700.0)
    arid = _short("arid", cells=900, t_end=700.0)
    water_p90 = float(np.percentile(
        water.world.fields["atmosphere.geographic_circulation_index"], 90))
    arid_p90 = float(np.percentile(
        arid.world.fields["atmosphere.geographic_circulation_index"], 90))
    assert water_p90 < 0.12
    assert arid_p90 > water_p90 * 4.0


def test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields():
    eng = _short(cells=1400, t_end=1200.0)
    w = eng.world
    seasonal_precip = w.fields["climate.seasonal_precipitation"]
    monsoon_corridor = w.fields["climate.monsoon_rainfall_corridor"]
    storm_corridor = w.fields["climate.storm_track_rainfall_corridor"]
    rain_shadow = w.fields["climate.rain_shadow_index"]
    regional_response = w.fields["climate.regional_precipitation_response"]
    moisture_flow_source = w.fields["atmosphere.moisture_flow_source"]
    moisture_flow_pathway = w.fields["atmosphere.moisture_flow_pathway"]
    moisture_source_basin_id = w.fields["atmosphere.moisture_source_basin_id"]
    moisture_flow_network_id = w.fields["climate.moisture_flow_network_id"]
    moisture_flow_precip_response = w.fields[
        "climate.moisture_flow_precipitation_response"]
    moisture_budget_region_id = w.fields["climate.moisture_budget_region_id"]
    precipitation_response_region_id = w.fields[
        "climate.precipitation_response_region_id"]
    receiver_catchment_id = w.fields["climate.receiver_catchment_id"]
    source_basin_supply_index = w.fields["climate.source_basin_supply_index"]
    receiver_supply_balance = w.fields[
        "climate.receiver_catchment_supply_balance"]
    receiver_supply_feedback = w.fields[
        "climate.receiver_supply_precipitation_feedback"]
    moisture = w.fields["atmosphere.moisture_access"]
    monsoon = w.fields["atmosphere.monsoon_potential"]
    seasonal_pressure = w.fields["atmosphere.seasonal_pressure_proxy"]
    precip_pressure_feedback = w.fields["atmosphere.precipitation_pressure_feedback"]
    hydro_wind_anomaly = w.fields["atmosphere.hydro_coupled_wind_anomaly"]
    hydro_coupling_residual = w.fields["climate.hydro_coupling_residual"]
    hydro_feedback_iteration_delta = w.fields[
        "climate.hydro_feedback_iteration_delta"]
    source_warmth = w.fields["atmosphere.source_ocean_warmth"]
    dry_len = w.fields["climate.dry_season_length"]
    precip = w.fields["climate.precipitation"]
    hydro_regions = w.objects.get("climate.hydroclimate_regions", [])
    moisture_networks = w.objects.get("climate.moisture_flow_networks", [])
    precipitation_response_regions = w.objects.get(
        "climate.precipitation_response_regions", [])
    receiver_catchments = w.objects.get("climate.receiver_catchments", [])
    land = w.land_mask()
    ocean = ~land

    assert seasonal_precip.shape == (4, w.grid.n)
    assert monsoon_corridor.shape == (4, w.grid.n)
    assert storm_corridor.shape == (4, w.grid.n)
    assert rain_shadow.shape == (4, w.grid.n)
    assert regional_response.shape == (4, w.grid.n)
    assert moisture_flow_source.shape == (4, w.grid.n)
    assert moisture_flow_pathway.shape == (4, w.grid.n)
    assert moisture_source_basin_id.shape == (4, w.grid.n)
    assert moisture_flow_network_id.shape == (4, w.grid.n)
    assert moisture_flow_precip_response.shape == (4, w.grid.n)
    assert moisture_budget_region_id.shape == (4, w.grid.n)
    assert precipitation_response_region_id.shape == (4, w.grid.n)
    assert receiver_catchment_id.shape == (4, w.grid.n)
    assert source_basin_supply_index.shape == (4, w.grid.n)
    assert receiver_supply_balance.shape == (4, w.grid.n)
    assert receiver_supply_feedback.shape == (4, w.grid.n)
    assert moisture.shape == (4, w.grid.n)
    assert monsoon.shape == (4, w.grid.n)
    assert seasonal_pressure.shape == (4, w.grid.n)
    assert precip_pressure_feedback.shape == (4, w.grid.n)
    assert hydro_wind_anomaly.shape == (4, w.grid.n, 3)
    assert hydro_coupling_residual.shape == (w.grid.n,)
    assert hydro_feedback_iteration_delta.shape == (w.grid.n,)
    assert source_warmth.shape == (4, w.grid.n)
    assert np.isfinite(seasonal_precip).all()
    assert np.isfinite(monsoon_corridor).all()
    assert np.isfinite(storm_corridor).all()
    assert np.isfinite(rain_shadow).all()
    assert np.isfinite(regional_response).all()
    assert np.isfinite(moisture_flow_source).all()
    assert np.isfinite(moisture_flow_pathway).all()
    assert np.isfinite(moisture_source_basin_id).all()
    assert np.isfinite(moisture_flow_network_id).all()
    assert np.isfinite(moisture_flow_precip_response).all()
    assert np.isfinite(moisture_budget_region_id).all()
    assert np.isfinite(precipitation_response_region_id).all()
    assert np.isfinite(receiver_catchment_id).all()
    assert np.isfinite(source_basin_supply_index).all()
    assert np.isfinite(receiver_supply_balance).all()
    assert np.isfinite(receiver_supply_feedback).all()
    assert np.isfinite(moisture).all()
    assert np.isfinite(monsoon).all()
    assert np.isfinite(precip_pressure_feedback).all()
    assert np.isfinite(hydro_wind_anomaly).all()
    assert np.isfinite(hydro_coupling_residual).all()
    assert np.isfinite(hydro_feedback_iteration_delta).all()
    assert (seasonal_precip >= 0.0).all()
    assert float(np.max(np.abs(precip_pressure_feedback))) <= 0.075 + 1.0e-9
    assert float(np.max(hydro_coupling_residual)) <= 0.20 + 1.0e-9
    assert float(np.max(hydro_feedback_iteration_delta)) <= 0.20 + 1.0e-9
    assert float(np.percentile(np.linalg.norm(hydro_wind_anomaly, axis=2), 99)) <= 0.65 + 1.0e-9
    assert (monsoon_corridor >= 0.0).all()
    assert (storm_corridor >= 0.0).all()
    assert (rain_shadow >= 0.0).all()
    assert (moisture_flow_source >= 0.0).all()
    assert (moisture_flow_pathway >= 0.0).all()
    assert (moisture_flow_precip_response >= 0.0).all()
    assert (source_basin_supply_index >= 0.0).all()
    assert (receiver_supply_balance >= 0.0).all()
    assert (receiver_supply_feedback >= 0.0).all()
    assert float(np.max(source_basin_supply_index)) <= 1.5 + 1.0e-9
    assert float(np.max(receiver_supply_balance)) <= 1.0 + 1.0e-9
    assert float(np.percentile(receiver_supply_feedback[:, land], 5)) >= 0.88
    assert float(np.percentile(receiver_supply_feedback[:, land], 95)) <= 1.14
    assert np.allclose(precip, seasonal_precip.mean(axis=0))
    assert float(dry_len.min()) >= 0.0
    assert float(dry_len.max()) <= 4.0
    if land.any() and ocean.any():
        assert float(np.percentile(moisture[:, ocean], 50)) > float(
            np.percentile(moisture[:, land], 25)
        )
        assert float(np.percentile(monsoon[:, land], 90)) > 0.02
        assert float(np.percentile(monsoon_corridor[:, land], 90)) > 0.005
        assert float(np.percentile(storm_corridor[:, land], 90)) > 0.005
        assert float(np.percentile(rain_shadow[:, land], 90)) > 0.0
        assert float(np.percentile(moisture_flow_source[:, ocean], 90)) > 0.0
        assert float(np.percentile(moisture_flow_pathway[:, land], 90)) > 0.0
        assert np.any(moisture_source_basin_id[:, land] >= 0.0)
        assert float(np.percentile(source_basin_supply_index[:, land], 90)) > 0.0
        assert float(np.percentile(receiver_supply_balance[:, land], 50)) > 0.0
        assert float(np.percentile(receiver_supply_feedback[:, land], 95)) > 1.0
        assert float(np.percentile(receiver_supply_feedback[:, land], 5)) < 1.0
        assert np.allclose(receiver_supply_feedback[:, ocean], 1.0)
        assert float(np.percentile(moisture_flow_precip_response[:, land], 95)) > 1.0
        assert float(np.percentile(moisture_flow_precip_response[:, land], 5)) < 1.0
        assert float(np.percentile(regional_response[:, land], 95)) > float(
            np.percentile(regional_response[:, land], 5)
        )
    assert isinstance(hydro_regions, list)
    assert hydro_regions
    region_kinds = {str(obj.get("kind")) for obj in hydro_regions}
    assert "monsoon_rainfall_corridor" in region_kinds
    assert "storm_track_rainfall_corridor" in region_kinds
    assert "rain_shadow_region" in region_kinds
    assert {str(obj.get("season")) for obj in hydro_regions} <= {"DJF", "MAM", "JJA", "SON"}
    assert max(float(obj.get("area_fraction", 0.0)) for obj in hydro_regions) > 0.0
    assert all(int(obj.get("cell_count", 0)) >= 2 for obj in hydro_regions)
    assert all(str(obj.get("type")) == "hydroclimate_region" for obj in hydro_regions)
    assert isinstance(moisture_networks, list)
    assert moisture_networks
    assert {str(obj.get("season")) for obj in moisture_networks} <= {"DJF", "MAM", "JJA", "SON"}
    assert max(float(obj.get("area_fraction", 0.0)) for obj in moisture_networks) > 0.0
    assert all(str(obj.get("type")) == "moisture_flow_network" for obj in moisture_networks)
    assert isinstance(precipitation_response_regions, list)
    assert precipitation_response_regions
    assert {str(obj.get("season")) for obj in precipitation_response_regions} <= {
        "DJF", "MAM", "JJA", "SON"
    }
    assert max(float(obj.get("area_fraction", 0.0))
               for obj in precipitation_response_regions) > 0.0
    assert all(
        str(obj.get("type")) == "precipitation_response_region"
        for obj in precipitation_response_regions)
    assert isinstance(receiver_catchments, list)
    assert receiver_catchments
    assert {str(obj.get("season")) for obj in receiver_catchments} <= {
        "DJF", "MAM", "JJA", "SON"
    }
    assert max(float(obj.get("area_fraction", 0.0))
               for obj in receiver_catchments) > 0.0
    assert all(str(obj.get("type")) == "receiver_catchment"
               for obj in receiver_catchments)
    assert max(float(obj.get("mean_source_basin_supply_index", 0.0))
               for obj in receiver_catchments) > 0.0
    assert max(float(obj.get("precipitation_supply_balance", 0.0))
               for obj in receiver_catchments) > 0.0

    detail = validation.climate_diagnostics(eng)
    assert detail["seasonality"]["precip_seasonality_p75"] >= 1.0
    assert detail["seasonality"]["has_monsoon_rainfall_corridor"]
    assert detail["seasonality"]["has_storm_track_rainfall_corridor"]
    assert detail["seasonality"]["has_rain_shadow_index"]
    assert detail["seasonality"]["has_regional_precipitation_response"]
    assert detail["seasonality"]["has_moisture_flow_source"]
    assert detail["seasonality"]["has_moisture_flow_pathway"]
    assert detail["seasonality"]["has_moisture_source_basin_id"]
    assert detail["seasonality"]["has_moisture_flow_network_id"]
    assert detail["seasonality"]["has_moisture_flow_precipitation_response"]
    assert detail["seasonality"]["has_moisture_budget_region_id"]
    assert detail["seasonality"]["has_precipitation_response_region_id"]
    assert detail["seasonality"]["has_receiver_catchment_id"]
    assert detail["seasonality"]["has_source_basin_supply_index"]
    assert detail["seasonality"]["has_receiver_catchment_supply_balance"]
    assert detail["seasonality"]["has_receiver_supply_precipitation_feedback"]
    assert detail["seasonality"]["monsoon_rainfall_corridor_land_p90"] > 0.0
    assert detail["seasonality"]["storm_track_rainfall_corridor_land_p90"] > 0.0
    assert detail["seasonality"]["rain_shadow_index_land_p90"] > 0.0
    assert detail["seasonality"]["moisture_flow_source_ocean_p90"] > 0.0
    assert detail["seasonality"]["moisture_flow_pathway_land_p90"] > 0.0
    assert detail["seasonality"]["moisture_source_basin_attributed_land_fraction"] > 0.0
    assert detail["seasonality"]["source_basin_supply_index_land_p50"] >= 0.0
    assert detail["seasonality"]["source_basin_supply_attributed_land_fraction"] > 0.0
    assert detail["seasonality"]["receiver_catchment_supply_balance_land_p50"] > 0.0
    assert detail["seasonality"]["receiver_supply_precipitation_feedback_land_p95"] > 1.0
    assert detail["seasonality"]["receiver_supply_precipitation_feedback_land_p05"] < 1.0
    assert detail["seasonality"]["receiver_supply_precipitation_feedback_ocean_abs_p95"] < 1.0e-9
    assert detail["seasonality"]["moisture_flow_precipitation_response_land_p95"] > 1.0
    assert detail["seasonality"]["moisture_flow_precipitation_response_land_p05"] < 1.0
    assert detail["seasonality"]["moisture_budget_region_count_p50"] >= 1.0
    assert detail["seasonality"]["has_hydroclimate_region_objects"]
    assert detail["seasonality"]["hydroclimate_region_object_count"] == len(hydro_regions)
    assert detail["seasonality"]["hydroclimate_region_kind_count"] >= 3
    assert detail["seasonality"]["hydroclimate_region_season_count"] >= 2
    assert detail["seasonality"]["largest_hydroclimate_region_area_fraction"] > 0.0
    assert detail["seasonality"]["has_moisture_flow_network_objects"]
    assert detail["seasonality"]["moisture_flow_network_object_count"] == len(moisture_networks)
    assert detail["seasonality"]["moisture_flow_network_season_count"] >= 1
    assert detail["seasonality"]["largest_moisture_flow_network_area_fraction"] > 0.0
    assert detail["seasonality"]["has_precipitation_response_region_objects"]
    assert detail["seasonality"]["precipitation_response_region_object_count"] == len(
        precipitation_response_regions)
    assert detail["seasonality"]["precipitation_response_region_kind_count"] >= 2
    assert detail["seasonality"]["precipitation_response_region_season_count"] >= 1
    assert detail["seasonality"]["largest_precipitation_response_region_area_fraction"] > 0.0
    assert detail["seasonality"]["has_receiver_catchment_objects"]
    assert detail["seasonality"]["receiver_catchment_object_count"] == len(
        receiver_catchments)
    assert detail["seasonality"]["receiver_catchment_season_count"] >= 1
    assert detail["seasonality"]["largest_receiver_catchment_area_fraction"] > 0.0
    assert detail["seasonality"]["receiver_catchment_budget_attribution_p50"] > 0.0
    assert detail["seasonality"]["receiver_catchment_source_supply_p50"] >= 0.0
    assert detail["seasonality"]["receiver_catchment_supply_balance_p50"] > 0.0
    assert detail["circulation"]["moisture_access_land_p75"] > 0.0


def test_ocean_currents_are_basin_constrained_and_transport_heat():
    eng = _short(cells=1400, t_end=1200.0)
    w = eng.world
    ocean = w.fields["ocean.solved_mask"] > 0.5
    land = ~ocean
    currents = w.fields["ocean.currents"]
    heat = w.fields["ocean.current_heat_transport"]
    upwelling = w.fields["ocean.upwelling"]
    basin = w.fields["ocean.basin_id"]
    streamfunction = w.fields["ocean.current_streamfunction"]
    gyre = w.fields["ocean.gyre_id"]
    boundary_type = w.fields["ocean.boundary_current_type"]
    strait_exchange = w.fields["ocean.strait_exchange"]
    wind_response = w.fields["ocean.wind_stress_current_response"]
    sst_anomaly = w.fields["ocean.sst_anomaly"]
    seasonal_sst = w.fields["climate.seasonal_sst"]
    ocean_heat_flux = w.fields["climate.ocean_heat_flux"]
    evaporation_feedback = w.fields["climate.ocean_evaporation_feedback"]
    coupling_residual = w.fields["climate.coupling_residual"]

    assert currents.shape == (w.grid.n, 3)
    assert heat.shape == (w.grid.n,)
    assert upwelling.shape == (w.grid.n,)
    assert basin.shape == (w.grid.n,)
    assert streamfunction.shape == (w.grid.n,)
    assert gyre.shape == (w.grid.n,)
    assert boundary_type.shape == (w.grid.n,)
    assert strait_exchange.shape == (w.grid.n,)
    assert wind_response.shape == (w.grid.n, 3)
    assert sst_anomaly.shape == (w.grid.n,)
    assert seasonal_sst.shape == (4, w.grid.n)
    assert ocean_heat_flux.shape == (w.grid.n,)
    assert evaporation_feedback.shape == (w.grid.n,)
    assert coupling_residual.shape == (w.grid.n,)
    assert np.isfinite(currents).all()
    assert np.isfinite(heat).all()
    assert np.isfinite(upwelling).all()
    assert np.isfinite(streamfunction).all()
    assert np.isfinite(wind_response).all()
    assert np.isfinite(sst_anomaly).all()
    assert np.isfinite(seasonal_sst).all()
    assert np.isfinite(ocean_heat_flux).all()
    assert np.isfinite(evaporation_feedback).all()
    assert np.isfinite(coupling_residual).all()
    assert float(np.linalg.norm(currents[land], axis=1).max(initial=0.0)) == pytest.approx(0.0)
    assert float(np.linalg.norm(wind_response[land], axis=1).max(initial=0.0)) == pytest.approx(0.0)
    assert float(np.abs(streamfunction[land]).max(initial=0.0)) == pytest.approx(0.0)
    assert float(np.abs(ocean_heat_flux[land]).max(initial=0.0)) == pytest.approx(0.0)
    assert float(np.abs(evaporation_feedback[land]).max(initial=0.0)) == pytest.approx(0.0)
    assert float(coupling_residual[land].max(initial=0.0)) == pytest.approx(0.0)
    normal_component = np.abs(np.einsum("ij,ij->i", currents, w.grid.xyz))
    assert float(normal_component.max()) < 1e-8
    response_normal = np.abs(np.einsum("ij,ij->i", wind_response, w.grid.xyz))
    assert float(response_normal.max()) < 1e-8
    assert int(np.unique(basin[ocean]).size) >= 1
    assert int(np.unique(gyre[ocean][gyre[ocean] > 0]).size) >= 1
    assert float(np.percentile(np.abs(streamfunction[ocean]), 95)) > 0.35
    assert float(np.average(sst_anomaly[ocean], weights=w.grid.cell_area[ocean])) == pytest.approx(
        0.0, abs=1e-6)
    assert float(np.percentile(np.abs(sst_anomaly[ocean]), 95)) > 0.35
    assert float(np.average(ocean_heat_flux[ocean], weights=w.grid.cell_area[ocean])) == pytest.approx(
        0.0, abs=1e-6)
    assert float(np.average(evaporation_feedback[ocean], weights=w.grid.cell_area[ocean])) == pytest.approx(
        0.0, abs=1e-6)
    assert float(np.percentile(np.abs(ocean_heat_flux[ocean]), 95)) > 0.25
    assert float(np.percentile(np.abs(evaporation_feedback[ocean]), 95)) <= 1.0 + 1.0e-9
    assert float(np.percentile(np.abs(evaporation_feedback[ocean]), 75)) > 0.01
    assert float(seasonal_sst[:, ocean].min(initial=CONSTANTS.ZERO_C)) >= (
        CONSTANTS.ZERO_C - 1.8 - 1e-8)
    assert float(np.percentile(coupling_residual[ocean], 95)) < 1.5
    assert float((ocean & (np.abs(boundary_type) > 0.5)).sum()) > 0
    assert float(strait_exchange[land].max(initial=0.0)) == pytest.approx(0.0)
    assert float(np.percentile(np.linalg.norm(wind_response[ocean], axis=1), 50)) > 0.03
    assert float(np.percentile(np.linalg.norm(wind_response[ocean], axis=1), 95)) <= 0.28
    assert float(np.percentile(np.linalg.norm(currents[ocean], axis=1), 95)) > 0.15
    assert float(np.percentile(np.abs(heat[ocean]), 95)) > 0.35
    assert float((ocean & (heat > 0.5)).sum()) > 0
    assert float((ocean & (heat < -0.5)).sum()) > 0
    assert float(np.percentile(upwelling[ocean], 95)) > 0.02

    detail = validation.climate_diagnostics(eng)
    ocean_detail = detail["ocean_currents"]
    assert ocean_detail["current_over_solved_land_cells_gt_1e_8"] == 0
    assert ocean_detail["has_seasonal_sst"]
    assert ocean_detail["has_ocean_heat_flux"]
    assert ocean_detail["has_coupling_residual"]
    assert ocean_detail["has_current_streamfunction"]
    assert ocean_detail["has_wind_stress_current_response"]
    assert ocean_detail["current_streamfunction_ocean_abs_p95"] > 0.35
    assert ocean_detail["gyre_count"] >= 1
    assert abs(ocean_detail["ocean_heat_flux_ocean_mean_C"]) < 1e-6
    assert ocean_detail["ocean_heat_flux_ocean_abs_p95_C"] > 0.25
    assert ocean_detail["coupling_residual_p95"] < 1.5
    assert ocean_detail["wind_stress_response_ocean_p50_mps"] > 0.03
    assert ocean_detail["wind_stress_response_ocean_p95_mps"] <= 0.28
    assert ocean_detail["wind_stress_response_land_speed_max_mps"] == pytest.approx(0.0)
    assert ocean_detail["sst_anomaly_ocean_abs_p95_C"] > 0.35
    assert ocean_detail["ocean_heat_transport_abs_p95_C"] > 0.35
    assert detail["coasts"]["coastal_temperature_asymmetry_index"] > 0.55


def test_cryosphere_cloud_and_vegetation_feedback_fields_are_bounded():
    eng = _short(cells=1400, t_end=1200.0)
    w = eng.world
    ocean = w.fields["ocean.solved_mask"] > 0.5
    land = ~ocean
    seasonal_sea_ice = w.fields["cryosphere.seasonal_sea_ice"]
    sea_ice = w.fields["cryosphere.sea_ice"]
    seasonal_snow = w.fields["cryosphere.seasonal_snow"]
    snow = w.fields["cryosphere.snow_persistence"]
    seasonal_cloud = w.fields["climate.seasonal_cloud_albedo_proxy"]
    cloud = w.fields["climate.cloud_albedo_proxy"]
    vegetation = w.fields["biosphere.vegetation_climate_feedback"]

    for arr, shape in [
        (seasonal_sea_ice, (4, w.grid.n)),
        (seasonal_snow, (4, w.grid.n)),
        (seasonal_cloud, (4, w.grid.n)),
        (sea_ice, (w.grid.n,)),
        (snow, (w.grid.n,)),
        (cloud, (w.grid.n,)),
        (vegetation, (w.grid.n,)),
    ]:
        assert arr.shape == shape
        assert np.isfinite(arr).all()
        assert float(arr.min()) >= -1.0e-9
        assert float(arr.max()) <= 1.0 + 1.0e-9

    assert float(np.abs(seasonal_sea_ice[:, land]).max(initial=0.0)) == pytest.approx(0.0)
    assert float(np.abs(sea_ice[land]).max(initial=0.0)) == pytest.approx(0.0)
    assert float(np.abs(seasonal_snow[:, ocean]).max(initial=0.0)) == pytest.approx(0.0)
    assert float(np.abs(snow[ocean]).max(initial=0.0)) == pytest.approx(0.0)
    assert float(np.abs(vegetation[ocean]).max(initial=0.0)) == pytest.approx(0.0)
    if land.any():
        assert float(np.percentile(vegetation[land], 95)) > float(
            np.percentile(vegetation[land], 5)
        )
    assert 0.05 <= float(np.percentile(cloud, 50)) <= 0.95

    detail = validation.climate_diagnostics(eng)
    cryo = detail["cryosphere"]
    assert cryo["has_seasonal_sea_ice"]
    assert cryo["has_seasonal_snow"]
    assert cryo["has_snow_persistence"]
    assert cryo["has_seasonal_cloud_albedo_proxy"]
    assert cryo["has_cloud_albedo_proxy"]
    assert cryo["has_vegetation_climate_feedback"]
    assert cryo["seasonal_sea_ice_land_abs_max"] == pytest.approx(0.0)
    assert cryo["snow_persistence_ocean_abs_max"] == pytest.approx(0.0)
    assert cryo["vegetation_feedback_ocean_abs_max"] == pytest.approx(0.0)
    assert cryo["sea_ice_adjacent_lat_band_jump_max"] <= 0.75
    assert validation.check_climate_diagnostics(eng).passed


def test_frozen_world_can_form_broad_cryosphere_feedback_cover():
    eng = _short("frozen", cells=900, t_end=700.0)
    cryo = validation.climate_diagnostics(eng)["cryosphere"]
    assert validation.check_climate_diagnostics(eng).passed
    assert cryo["has_seasonal_sea_ice"]
    assert cryo["has_snow_persistence"]
    assert (
        cryo["seasonal_sea_ice_ocean_p95"] > 0.25
        or cryo["snow_persistence_land_p95"] > 0.25
    )


def test_geography_primitives_cover_sphere_and_follow_landforms():
    eng = _short(cells=1400, t_end=1200.0)
    w = eng.world
    ocean = w.fields["ocean.solved_mask"] > 0.5
    land = ~ocean
    continent = w.fields["climate.continent_id"].astype(int)
    basin = w.fields["ocean.basin_id"].astype(int)
    interiority = w.fields["climate.continent_interiority"]
    coast_orientation = w.fields["climate.coast_orientation"]
    coast_strength = w.fields["climate.coast_strength"]
    shelf = w.fields["ocean.shelf_index"]
    strait = w.fields["ocean.strait_index"]
    depth_province = w.fields["ocean.depth_province"].astype(int)
    margin_type = w.fields["ocean.margin_type"].astype(int)
    gateway = w.fields["ocean.gateway_id"]
    shelf_width = w.fields["ocean.shelf_width"]
    barrier = w.fields["terrain.barrier_index"]
    wind_gap = w.fields["terrain.wind_gap_index"]

    assert continent.shape == (w.grid.n,)
    assert basin.shape == (w.grid.n,)
    assert interiority.shape == (w.grid.n,)
    assert coast_orientation.shape == (w.grid.n, 3)
    assert shelf.shape == (w.grid.n,)
    assert strait.shape == (w.grid.n,)
    assert depth_province.shape == (w.grid.n,)
    assert margin_type.shape == (w.grid.n,)
    assert gateway.shape == (w.grid.n,)
    assert shelf_width.shape == (w.grid.n,)
    assert barrier.shape == (w.grid.n,)
    assert wind_gap.shape == (w.grid.n,)
    assert np.isfinite(coast_orientation).all()
    assert np.isfinite(interiority).all()
    assert np.isfinite(shelf).all()
    assert np.isfinite(strait).all()
    assert np.isfinite(depth_province).all()
    assert np.isfinite(margin_type).all()
    assert np.isfinite(gateway).all()
    assert np.isfinite(shelf_width).all()
    assert np.isfinite(barrier).all()

    assert (continent[land] >= 0).all()
    assert (continent[ocean] < 0).all()
    assert (basin[ocean] >= 0).all()
    assert (basin[land] < 0).all()
    assert (depth_province[land] == 0).all()
    assert (depth_province[ocean] > 0).all()
    assert (margin_type[land] == 0).all()
    assert (margin_type[ocean] > 0).all()
    assert (gateway[land] < 0.0).all()
    assert (shelf_width[land] == 0.0).all()
    assert int(np.unique(continent[land]).size) >= 1
    assert int(np.unique(basin[ocean]).size) >= 1
    assert (0.0 <= interiority).all() and (interiority <= 1.0).all()
    assert (0.0 <= shelf).all() and (shelf <= 1.0).all()
    assert (0.0 <= strait).all() and (strait <= 1.0).all()
    assert (0.0 <= barrier).all() and (barrier <= 1.0).all()
    assert (0.0 <= wind_gap).all() and (wind_gap <= 1.0).all()

    land_coast = np.zeros(w.grid.n, dtype=bool)
    ocean_coast = np.zeros(w.grid.n, dtype=bool)
    for c, nbs in enumerate(w.grid.neighbors):
        if land[c] and ocean[nbs].any():
            land_coast[c] = True
        if ocean[c] and land[nbs].any():
            ocean_coast[c] = True
    speed = np.linalg.norm(coast_orientation, axis=1)
    normal_component = np.abs(np.einsum("ij,ij->i", coast_orientation, w.grid.xyz))
    assert float(normal_component.max()) < 1e-8
    assert float(np.percentile(speed[land_coast], 50)) > 0.9
    assert float(np.percentile(speed[~land_coast], 95)) == pytest.approx(0.0)
    assert float(coast_strength[land & ~land_coast].max(initial=0.0)) == pytest.approx(0.0)

    coast_distance = w.fields["climate.coast_distance"]
    deep_ocean = ocean & (coast_distance > np.percentile(coast_distance[ocean], 65))
    assert float(np.percentile(shelf[ocean_coast], 75)) > float(
        np.percentile(shelf[deep_ocean], 75)
    ) + 0.20
    assert float(strait[land].max(initial=0.0)) == pytest.approx(0.0)

    elev = w.get_field("terrain.elevation_m") - w.sea_level
    high = land & (elev >= np.percentile(elev[land], 80))
    low = land & (elev <= np.percentile(elev[land], 35))
    assert float(np.percentile(barrier[high], 75)) > float(
        np.percentile(barrier[low], 75)
    ) + 0.10
    assert len(w.objects.get("climate.continents", [])) >= 1
    assert len(w.objects.get("ocean.basins", [])) >= 1
    assert "ocean.margins" in w.objects
    assert "ocean.gateways" in w.objects
    assert len(w.objects.get("climate.coastline_segments", [])) >= 1

    detail = validation.climate_diagnostics(eng)["geography"]
    assert detail["land_without_continent_id_cells"] == 0
    assert detail["ocean_without_basin_id_cells"] == 0
    assert detail["shelf_contrast_near_minus_deep"] > 0.20
    assert detail["barrier_contrast_high_minus_low"] > 0.10


def test_geography_primitives_handle_waterworld_and_arid_layouts():
    water = _short("waterworld", cells=900, t_end=700.0)
    arid = _short("arid", cells=900, t_end=700.0)
    water_geo = validation.climate_diagnostics(water)["geography"]
    arid_geo = validation.climate_diagnostics(arid)["geography"]

    assert validation.check_climate_diagnostics(water).passed
    assert validation.check_climate_diagnostics(arid).passed
    assert water_geo["ocean_basin_count"] >= 1
    assert water_geo["ocean_without_basin_id_cells"] == 0
    assert water_geo["land_without_continent_id_cells"] == 0
    assert arid_geo["continent_count"] >= 1
    assert arid_geo["land_without_continent_id_cells"] == 0
    assert arid_geo["largest_continent_fraction_of_land"] >= water_geo[
        "largest_continent_fraction_of_land"
    ]


def test_cold_boundary_currents_suppress_local_evaporation():
    eng = _short(cells=1400, t_end=1200.0)
    w = eng.world
    ocean = w.fields["ocean.solved_mask"] > 0.5
    heat = w.fields["ocean.current_heat_transport"]
    upwelling = w.fields["ocean.upwelling"]
    evap = w.fields["climate.evaporation"]
    subtropical = (np.abs(w.grid.lat) >= 10.0) & (np.abs(w.grid.lat) <= 45.0)
    cold = ocean & subtropical & (heat < -0.5) & (upwelling > 0.02)
    warm = ocean & subtropical & (heat > 0.5)

    assert int(cold.sum()) > 3
    assert int(warm.sum()) > 3
    cold_evap = float(np.average(evap[cold], weights=w.grid.cell_area[cold]))
    warm_evap = float(np.average(evap[warm], weights=w.grid.cell_area[warm]))
    assert cold_evap < warm_evap


def test_crust_origin_stability_and_explanations_are_available():
    eng = _short(cells=1800, t_end=2500.0)
    w = eng.world
    origin = w.get_field("crust.origin")
    reworked = w.get_field("crust.reworked_age_myr")
    stability = w.get_field("crust.stability")
    domain = w.get_field("crust.domain")
    continent_id = w.get_field("tectonics.continent_id")
    terrane_id = w.get_field("tectonics.terrane_id")
    terrain_province = w.get_field("terrain.province")
    continental_detail = w.get_field("terrain.continental_detail")
    wilson_phase = w.get_field("archive.wilson_cycle_phase")
    deformation_intensity = w.get_field("tectonics.deformation_intensity")
    deformation_style = w.get_field("tectonics.deformation_style")
    crust_type = w.get_field("crust.type").astype(int)
    assert origin.shape == (w.grid.n,)
    assert reworked.shape == (w.grid.n,)
    assert stability.shape == (w.grid.n,)
    assert domain.shape == (w.grid.n,)
    assert continent_id.shape == (w.grid.n,)
    assert terrane_id.shape == (w.grid.n,)
    assert terrain_province.shape == (w.grid.n,)
    assert continental_detail.shape == (w.grid.n,)
    assert wilson_phase.shape == (w.grid.n,)
    assert deformation_intensity.shape == (w.grid.n,)
    assert deformation_style.shape == (w.grid.n,)
    assert set(np.unique(origin.astype(int))).issubset({0, 1, 2, 3, 4, 5})
    assert np.all((0.0 <= stability) & (stability <= 1.0))
    assert set(np.unique(domain.astype(int))).issubset({0, 1, 2, 3, 4, 5, 6})
    assert np.all(continent_id[crust_type == 1] >= 0)
    assert np.all(continent_id[crust_type == 0] < 0)
    assert np.all(terrane_id[np.isfinite(terrane_id)] >= -1)
    assert set(np.unique(terrain_province.astype(int))).issubset(set(range(9)))
    assert set(np.unique(continental_detail.astype(int))).issubset(set(range(8)))
    assert set(np.unique(wilson_phase.astype(int))).issubset(set(range(7)))
    assert np.all((0.0 <= deformation_intensity) & (deformation_intensity <= 1.0))
    assert set(np.unique(deformation_style.astype(int))).issubset(set(range(7)))
    assert np.isin(terrain_province.astype(int), [1, 2, 3]).any()
    land = w.get_field("terrain.elevation_m") >= w.sea_level
    assert np.isin(continental_detail[land].astype(int), [1, 2, 3, 4, 5, 6, 7]).any()
    assert (origin == 0).any()      # ridge oceanic crust
    assert np.isin(origin, [1, 2, 3, 5]).any()

    boundaries = w.networks.get("tectonics.boundaries", {})
    assert "ridge" in boundaries and "trench" in boundaries and "suture" in boundaries
    assert len(boundaries["ridge"]) > 0
    assert len(boundaries["trench"]) > 0
    assert len(w.objects.get("tectonics.boundary_objects", [])) > 0
    assert len(w.objects.get("tectonics.wilson_cycles", [])) > 0
    assert isinstance(w.objects.get("tectonics.deforming_networks", []), list)
    assert len(w.objects.get("tectonics.deforming_networks", [])) > 0
    assert isinstance(w.objects.get("tectonics.ocean_gateways", []), list)
    assert len(w.objects.get("tectonics.ocean_gateways", [])) > 0
    assert isinstance(w.objects.get("tectonics.cratons", []), list)
    assert len(w.objects.get("tectonics.continents", [])) > 0
    assert isinstance(w.objects.get("tectonics.terranes", []), list)
    continent = w.objects["tectonics.continents"][0]
    assert {"id", "area_fraction", "core_fraction", "margin_fraction",
            "width_p50_steps"}.issubset(continent)
    assert isinstance(w.objects.get("tectonics.plumes", []), list)
    assert isinstance(w.objects.get("tectonics.lips", []), list)

    event_types = {e.type for e in eng.bus.events}
    assert event_types & {
        "ridge_birth",
        "ocean_basin_opening",
        "subduction_initiation",
        "suture_formation",
        "passive_margin_formation",
    }
    assert event_types & {"ocean_gateway_opened", "ocean_gateway_closed"}
    assert event_types & {"rift_birth", "trench_birth", "arc_birth", "orogen_built"}
    assert "plume_head" in event_types
    assert "large_igneous_province" in event_types

    story = eng.archive.explain_cell(0)
    assert "crust_origin" in story
    assert "crust_reworked_age_myr" in story
    assert "crust_stability" in story
    assert "tectonic_objects" in story


def test_local_plate_reorganization_preserves_deep_time_continuity():
    eng = _short(cells=2500, t_end=4500.0)
    detail = validation.tectonic_diagnostics(eng)
    plates = detail["plate_fragmentation"]
    crust = detail["crust_distribution"]
    frames = detail["frame_continuity"]

    assert plates["max_components_per_plate"] <= 2
    assert plates["total_plate_components"] <= plates["n_active_plates"] + 4
    assert plates["min_largest_component_share"] > 0.80
    assert frames["late_archive_max_components_per_plate"] <= 2
    assert frames["late_archive_min_largest_component_share"] > 0.80
    assert frames["late_max_plate_crust_composition_delta"] < 0.35
    assert frames["late_max_plate_crust_composition_delta_per_100myr"] < 0.09
    assert crust["stable_craton_fraction_gt075"] > 0.01
    assert len(eng.world.objects.get("tectonics.cratons", [])) > 0
    assert crust["oceanic_age_p95_myr"] < 320.0
    assert crust["oceanic_old_fraction_gt300_myr"] < 0.05
    assert detail["seafloor_age_geometry"]["ridge_ocean_age_p75_myr"] <= 30.0
    if detail["seafloor_age_geometry"]["ocean_with_ridge_path_fraction"] > 0.40:
        assert detail["seafloor_age_geometry"]["age_distance_correlation"] > 0.45
    assert frames["plate_label_persistence_mean"] > 0.25


def test_tectonic_diagnostics_respect_stagnant_lid_preset():
    eng = _short("stagnant_lid", cells=1200, t_end=1500.0)
    detail = validation.tectonic_diagnostics(eng)
    assert detail["plate_fragmentation"]["n_active_plates"] == 1
    assert detail["plate_fragmentation"]["max_components_per_plate"] >= 1
    assert detail["crust_distribution"]["continental_area_fraction"] > 0.0
    assert detail["boundaries"]["boundary_cell_fraction"] <= 1.0


@pytest.mark.parametrize("preset", ["waterworld", "arid", "stagnant_lid",
                                    "tidally_locked", "frozen"])
def test_divergent_worlds_are_valid(preset):
    eng = _short(preset, cells=2000, t_end=2500.0)
    assert validation.check_conservation(eng).passed
    assert validation.check_topology(eng).passed
    if preset == "waterworld":
        land_fraction = eng.world.land_fraction()
        assert 0.005 <= land_fraction <= 0.08


def test_explain_cell_has_causal_story():
    eng = _short()
    story = eng.archive.explain_cell(0)
    assert "elevation_m" in story and "crust_type" in story
    assert "events" in story


def test_earthlike_map_reasonableness_regression():
    eng = _short(cells=2500, t_end=4500.0)
    w = eng.world
    elev = w.get_field("terrain.elevation_m")
    rel = elev - w.sea_level
    land = rel >= 0
    crust_type = w.get_field("crust.type").astype(int)
    temp_c = w.get_field("climate.surface_temperature") - 273.15
    precip = w.get_field("climate.precipitation")
    biome = w.get_field("biosphere.biome").astype(int)

    assert 0.20 <= w.land_fraction() <= 0.38
    assert int(((crust_type == 0) & land & (rel > 1500)).sum()) == 0
    if land.any():
        oceanic_land_fraction = float(
            w.grid.cell_area[(crust_type == 0) & land].sum()
            / w.grid.cell_area[land].sum()
        )
        assert oceanic_land_fraction < 0.25
        land_rel = rel[land]
        assert float(np.percentile(land_rel, 95)) > 1500.0
        assert float(np.max(land_rel)) > 2200.0
        assert float(np.mean(land_rel)) < 1800.0
        continental_detail = w.get_field("terrain.continental_detail").astype(int)
        detail_land = continental_detail[land]
        assert set(np.unique(detail_land)).issubset(set(range(1, 8)))
        assert 0.02 <= float(np.isin(detail_land, [5, 6]).mean()) <= 0.35
        assert float((detail_land == 7).mean()) < 0.65
        morphology = validation.tectonic_diagnostics(eng)["morphology"]["exposed_land"]
        assert morphology["largest_component_area_fraction_of_mask"] < 0.75
        assert morphology["component_count"] <= 20
        assert morphology["narrow_neck_cells_per_1000_mask_cells"] < 25.0
    ocean = ~land
    frontier = land.copy()
    seen = frontier.copy()
    ocean_pass = np.full(w.grid.n, -1, dtype=int)
    for p in range(1, 6):
        nxt = np.zeros(w.grid.n, dtype=bool)
        for c in np.where(frontier)[0]:
            nxt[w.grid.neighbors[int(c)]] = True
        nxt &= ocean & ~seen
        ocean_pass[nxt] = p
        seen |= nxt
        frontier = nxt
    boundaries = w.networks.get("tectonics.boundaries", {})
    trench = np.zeros(w.grid.n, dtype=bool)
    trench[np.asarray(boundaries.get("trench", []), dtype=int)] = True
    near_shelf = ocean & (ocean_pass == 1) & ~trench
    far_ocean = ocean & ((ocean_pass < 0) | (ocean_pass >= 4)) & ~trench
    if int(near_shelf.sum()) > 20 and int(far_ocean.sum()) > 20:
        near_depth = -rel[near_shelf]
        far_depth = -rel[far_ocean]
        assert float(np.percentile(near_depth, 75)) < 800.0
        assert float(np.median(far_depth)) > float(np.percentile(near_depth, 75)) + 1200.0
        assert float((far_depth < 1500.0).mean()) < 0.12
    depth_province = w.fields["ocean.depth_province"].astype(int)
    margin_type = w.fields["ocean.margin_type"].astype(int)
    shelf_width = w.fields["ocean.shelf_width"]
    gateway = w.fields["ocean.gateway_id"]
    assert (depth_province[land] == 0).all()
    assert (depth_province[ocean] > 0).all()
    assert (margin_type[land] == 0).all()
    assert (margin_type[ocean] > 0).all()
    assert (gateway[land] < 0.0).all()
    province_shelf = ocean & (depth_province == 1)
    province_abyss = ocean & (depth_province == 4)
    province_trench = ocean & (depth_province == 6)
    ocean_depth = -rel
    if int(province_shelf.sum()) > 20 and int(province_abyss.sum()) > 20:
        shelf_p75 = float(np.percentile(ocean_depth[province_shelf], 75))
        abyss_p50 = float(np.median(ocean_depth[province_abyss]))
        assert shelf_p75 < 900.0
        assert abyss_p50 > shelf_p75 + 1200.0
    nearshore = ocean & (shelf_width == 1.0) & (depth_province != 6)
    if int(nearshore.sum()) > 20:
        assert float((nearshore & (ocean_depth > 2500.0)).sum() / nearshore.sum()) < 0.08
    if int(province_trench.sum()) > 10:
        assert float(np.median(ocean_depth[province_trench])) > 3000.0
    assert len(w.objects.get("ocean.basins", [])) >= 1
    assert "ocean.margins" in w.objects
    continental_kinds = [
        str(obj.get("kind", ""))
        for obj in w.objects.get("terrain.continental_landforms", [])
    ]
    ocean_fabric_kinds = [
        str(obj.get("kind", ""))
        for obj in w.objects.get("terrain.ocean_fabric", [])
    ]
    assert continental_kinds.count("interior_basin") >= 1
    assert continental_kinds.count("old_subdued_orogen") >= 1
    assert ocean_fabric_kinds.count("transform_fault") >= 1
    mean_temp_c = float(np.average(temp_c, weights=w.grid.cell_area))
    assert 8.0 <= mean_temp_c < 30.0
    assert float(temp_c.max()) < 60.0
    lat_bins = np.arange(-90, 91, 10)
    band_means = []
    for lo, hi in zip(lat_bins[:-1], lat_bins[1:]):
        mask = (w.grid.lat >= lo) & (w.grid.lat < hi)
        band_means.append(float(np.average(temp_c[mask], weights=w.grid.cell_area[mask])))
    assert max(abs(a - b) for a, b in zip(band_means[:-1], band_means[1:])) < 20.0
    assert float(np.percentile(precip[land], 75)) > 300.0
    assert np.isin(biome[land], [3, 4, 6]).any()

    cm = MapCompiler(w, eng.archive).compile(width=64, height=32, n_starts=4)
    compiler_detail = validation.compiler_consistency_metrics(eng, cm)
    assert compiler_detail["passed_envelope"]
    assert compiler_detail["broad_land_to_water_fraction"] < 0.08
    assert compiler_detail["broad_ocean_to_land_fraction"] < 0.06
    assert compiler_detail["shelf_as_deep_ocean_fraction"] < 0.20
    assert compiler_detail["terrain_elevation_sign_mismatch_fraction"] < 0.02
    assert float((cm.terrain == 5).mean()) < 0.15
    assert int((cm.resources != "").sum()) <= 24
    assert len(cm.starts) == 4
    assert max(cm.hazard[r, c] for r, c in cm.starts) <= 1.0
    assert cm.fairness["local_yield_cv"] < 0.45


def test_p12_release_gate_writes_summary_and_contact_sheet(tmp_path):
    config = P12RunConfig(
        presets=("earthlike",),
        cells=900,
        t_end_myr=1200.0,
        frames=3,
        hex_width=32,
        hex_height=16,
        starts=4,
        global_overrides={"diagnostics.test_global_override": 2.5},
    )
    summary = run_p12_release_gate(config, tmp_path)
    assert summary["schema"] == "aevum.p12_tectonics_release_summary.v1"
    assert summary["config"]["global_overrides"] == {
        "diagnostics.test_global_override": 2.5
    }
    assert summary["release_decision"]["status"] in {"pass", "warn"}
    assert len(summary["entries"]) == 1
    entry = summary["entries"][0]
    assert entry["preset"] == "earthlike"
    assert entry["release_gate"]["passed"]
    assert set(entry) >= {
        "morphology",
        "crust",
        "terrain_detail",
        "tectonic_object_telemetry",
        "p26_ribbon_drivers",
        "p26_rework_footprint",
        "p26_deforming_networks",
        "p29_inland_geomorphology",
        "earth_geomorphology_coverage",
        "earth_reference_distribution",
        "ocean_geography",
        "archive_continuity",
        "seam_continuity",
        "compiler",
        "climate_facing_prerequisites",
    }
    assert entry["terrain_detail"]["land_elevation_p95_m"] >= 0.0
    assert 0.0 <= entry["terrain_detail"]["orogen_or_plateau_fraction_of_land"] <= 1.0
    telemetry = entry["tectonic_object_telemetry"]
    assert set(telemetry) >= {
        "rift_system_object_count",
        "rift_system_total_area_fraction",
        "rift_system_max_cell_count",
        "breakup_seaway_object_count",
        "breakup_seaway_total_area_fraction",
        "breakup_seaway_max_topology_score",
        "tectonics_continent_shape_pressure",
        "tectonics_background_continent_shape_maintenance",
        "tectonics_passive_margin_progradation_cells",
        "tectonics_continent_gain_cells",
        "tectonics_continent_loss_cells",
        "tectonics_unforced_continent_gain_blocked",
        "tectonics_unforced_continent_loss_blocked",
        "terrain_breakup_seaway_openings",
        "terrain_breakup_seaway_area_fraction",
        "terrain_breakup_seaway_source_reuse",
        "terrain_breakup_seaway_opened_corridor_count",
        "terrain_breakup_seaway_opened_corridor_area_fraction",
        "terrain_breakup_seaway_attempt_count",
        "terrain_breakup_seaway_applied_attempt_count",
        "terrain_breakup_seaway_reject_reasons",
        "terrain_breakup_seaway_top_attempts",
        "terrain_land_component_stages",
        "terrain_largest_landmass_shave_area_fraction",
        "terrain_deformation_relief_area_fraction",
        "terrain_deformation_relief_mean_m",
        "breakup_component_evaluated_count",
        "breakup_component_eligible_count",
        "breakup_component_candidate_total",
        "breakup_component_accepted_total",
        "breakup_component_skip_reasons",
        "breakup_component_top",
    }
    assert telemetry["rift_system_object_count"] >= 0
    assert telemetry["rift_system_total_area_fraction"] >= 0.0
    assert telemetry["rift_system_max_cell_count"] >= 0
    assert telemetry["breakup_seaway_object_count"] >= 0
    assert telemetry["breakup_seaway_total_area_fraction"] >= 0.0
    assert telemetry["terrain_breakup_seaway_area_fraction"] >= 0.0
    assert telemetry["terrain_breakup_seaway_source_reuse"] >= 0.0
    assert telemetry["terrain_breakup_seaway_opened_corridor_count"] >= 0
    assert telemetry["terrain_breakup_seaway_opened_corridor_area_fraction"] >= 0.0
    assert telemetry["terrain_breakup_seaway_attempt_count"] >= 0
    assert telemetry["terrain_breakup_seaway_applied_attempt_count"] >= 0
    assert isinstance(telemetry["terrain_breakup_seaway_reject_reasons"], dict)
    assert isinstance(telemetry["terrain_breakup_seaway_top_attempts"], list)
    assert isinstance(telemetry["terrain_land_component_stages"], list)
    assert telemetry["breakup_component_evaluated_count"] >= 0
    assert telemetry["breakup_component_eligible_count"] >= 0
    assert isinstance(telemetry["breakup_component_skip_reasons"], dict)
    assert isinstance(telemetry["breakup_component_top"], list)
    ribbon_drivers = entry["p26_ribbon_drivers"]
    assert ribbon_drivers["schema"] == "aevum.p26_ribbon_drivers.v1"
    assert set(ribbon_drivers) >= {
        "summary",
        "top_exposed_land_components",
        "top_continental_crust_components",
        "top_exposed_continental_components",
    }
    assert ribbon_drivers["summary"]["exposed_land_component_count"] >= 0
    assert ribbon_drivers["summary"]["continental_component_count"] >= 0
    assert isinstance(ribbon_drivers["top_exposed_land_components"], list)
    rework = entry["p26_rework_footprint"]
    assert rework["schema"] == "aevum.p26_rework_footprint.v1"
    assert set(rework["metrics"]) >= {
        "active_rework_area_fraction_of_continental",
        "active_rework_outside_corridor_fraction_of_active",
        "active_rework_inside_corridor_fraction_of_active",
        "active_to_corridor_area_ratio",
        "overbroad_recent_rework",
    }
    deforming = entry["p26_deforming_networks"]
    assert deforming["schema"] == "aevum.p26_deforming_networks.v1"
    assert set(deforming["metrics"]) >= {
        "active_deformation_area_fraction_of_world",
        "active_deformation_area_fraction_of_continental",
        "core_area_fraction_of_active",
        "shoulder_area_fraction_of_active",
        "land_ribbon_deformation_coverage_fraction",
        "continental_ribbon_deformation_coverage_fraction",
        "deforming_network_object_count",
    }
    assert set(deforming["style_summaries"]) >= {
        "collision_core",
        "collision_shoulder",
        "subduction_core",
        "subduction_shoulder",
        "rift",
        "transform",
    }
    assert isinstance(deforming["top_deforming_network_objects"], list)
    assert isinstance(deforming["diagnostic_hints"], dict)
    inland = entry["p29_inland_geomorphology"]
    assert inland["schema"] == "aevum.p29_inland_geomorphology.v1"
    assert set(inland["metrics"]) >= {
        "inland_area_fraction_of_land",
        "inland_relief_p95_p05_m",
        "inland_flat_lowland_fraction",
        "inland_highland_fraction",
        "inland_detail_diversity_gt2pct",
        "inland_object_kind_count",
    }
    assert set(inland["diagnostic_hints"]) >= {
        "monotone_flat_inland",
        "low_inland_object_diversity",
        "broad_lowland_mode_present",
        "highland_tail_present",
    }
    assert isinstance(inland["top_inland_landform_objects"], list)
    earth_benches = summary["earth_geomorphology_benchmarks"]
    assert earth_benches["all_passed"]
    assert set(earth_benches["suites"]) == {"E1", "E2", "E3", "E4", "E5"}
    coverage = entry["earth_geomorphology_coverage"]
    assert coverage["schema"] == "aevum.earth_geomorphology_coverage.v1"
    assert coverage["feature_count"] >= 25
    assert coverage["missing_feature_count"] == 0
    assert coverage["acceptance"]["no_non_ice_feature_none"]
    assert coverage["acceptance"]["plate_boundary_and_ocean_features_partial"]
    assert coverage["acceptance"]["parentless_major_landform_fraction_ok"]
    assert coverage["parentless_major_landform_fraction"] < 0.10
    assert coverage["group_world_feature_counts"]["continental_interior"] >= 3
    assert coverage["group_world_feature_counts"]["margin_shelf"] >= 3
    assert coverage["group_world_feature_counts"]["ocean_basin_fabric"] >= 2
    assert coverage["group_world_feature_counts"]["arc_plume"] >= 3
    assert coverage["scores"]["ice_geomorphology_score"] >= 0.5
    earth_ref = entry["earth_reference_distribution"]
    assert earth_ref["schema"] == "aevum.earth_reference_distribution.v1"
    assert earth_ref["acceptance"]["screening_only"]
    assert not earth_ref["acceptance"]["source_is_direct_etopo_raster"]
    assert any(src["id"] == "NOAA_ETOPO_2022" for src in earth_ref["reference_sources"])
    assert set(earth_ref["metrics"]) >= {
        "land_fraction",
        "land_elevation_mean_m",
        "shelf_fraction_of_ocean",
        "abyss_fraction_of_ocean",
        "land_ribbon_fraction_gt_0_5",
    }
    assert isinstance(earth_ref["out_of_envelope"], list)
    assert not any(
        "ancient continental crust fraction" in warning
        for warning in entry["release_gate"]["warnings"]
    )
    assert (tmp_path / "p12_tectonics_release_summary.json").exists()
    assert (tmp_path / "p12_preset_matrix_contact_sheet.png").exists()
