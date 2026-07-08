import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.state import WorldState
from aevum.core.units import CONSTANTS
from aevum.modules.biosphere import BiosphereModule
from aevum.spec.presets import get_preset


def test_biome_generalization_preserves_climate_supported_forest_and_tropics():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("terrain.elevation_m", np.full(grid.n, 100.0))

    temp_c = np.full(grid.n, 18.0)
    precip = np.full(grid.n, 360.0)
    forest_patch = (
        (grid.lat > 32.0) & (grid.lat < 56.0)
        & (grid.lon > -60.0) & (grid.lon < 35.0)
    )
    tropical_patch = (
        (np.abs(grid.lat) < 12.0)
        & (grid.lon > 55.0) & (grid.lon < 145.0)
    )
    temp_c[forest_patch] = 15.0
    precip[forest_patch] = 760.0
    temp_c[tropical_patch] = 26.0
    precip[tropical_patch] = 980.0

    world.set_field("climate.surface_temperature", temp_c + 273.15)
    world.set_field(
        "climate.seasonal_temperature",
        np.repeat((temp_c + 273.15)[None, :], 4, axis=0),
    )
    world.set_field("climate.precipitation", precip)
    world.set_field(
        "climate.seasonal_precipitation",
        np.repeat(precip[None, :], 4, axis=0),
    )

    module = BiosphereModule()
    biome = module._biomes(world, np.zeros(grid.n))

    assert int(forest_patch.sum()) >= 8
    assert int(tropical_patch.sum()) >= 8
    assert int(((biome == 4) & forest_patch).sum()) >= 2
    assert int(((biome == 6) & tropical_patch).sum()) >= 2


def test_cold_dry_land_is_not_reclassified_as_desert():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("terrain.elevation_m", np.full(grid.n, 100.0))

    temp_c = np.full(grid.n, 12.0)
    precip = np.full(grid.n, 420.0)
    cold_dry = grid.lat > 62.0
    temp_c[cold_dry] = -9.0
    precip[cold_dry] = 120.0

    seasonal_temp_c = np.repeat(temp_c[None, :], 4, axis=0)
    seasonal_temp_c[0, cold_dry] = -18.0
    world.set_field("climate.surface_temperature", temp_c + 273.15)
    world.set_field("climate.seasonal_temperature", seasonal_temp_c + 273.15)
    world.set_field("climate.precipitation", precip)
    world.set_field(
        "climate.seasonal_precipitation",
        np.repeat(precip[None, :], 4, axis=0),
    )

    biome = BiosphereModule()._biomes(world, np.zeros(grid.n))

    assert int(cold_dry.sum()) >= 8
    assert int(((biome == 5) & cold_dry).sum()) >= int(0.6 * cold_dry.sum())
    assert int(((biome == 2) & cold_dry).sum()) == 0


def test_moist_temperate_lower_precip_band_can_be_forest():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("terrain.elevation_m", np.full(grid.n, 100.0))

    temp_c = np.full(grid.n, 12.0)
    precip = np.full(grid.n, 540.0)
    world.set_field("climate.surface_temperature", temp_c + 273.15)
    world.set_field(
        "climate.seasonal_temperature",
        np.repeat((temp_c + 273.15)[None, :], 4, axis=0),
    )
    world.set_field("climate.precipitation", precip)
    world.set_field(
        "climate.seasonal_precipitation",
        np.repeat(precip[None, :], 4, axis=0),
    )

    biome = BiosphereModule()._biomes(world, np.zeros(grid.n))

    assert int((biome == 4).sum()) >= int(0.95 * grid.n)
    assert int((biome == 3).sum()) == 0


def test_cool_moist_land_uses_lower_forest_precip_threshold():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("terrain.elevation_m", np.full(grid.n, 100.0))

    temp_c = np.full(grid.n, 6.0)
    precip = np.full(grid.n, 320.0)
    world.set_field("climate.surface_temperature", temp_c + 273.15)
    world.set_field(
        "climate.seasonal_temperature",
        np.repeat((temp_c + 273.15)[None, :], 4, axis=0),
    )
    world.set_field("climate.precipitation", precip)
    world.set_field(
        "climate.seasonal_precipitation",
        np.repeat(precip[None, :], 4, axis=0),
    )

    biome = BiosphereModule()._biomes(world, np.zeros(grid.n))

    assert int((biome == 4).sum()) >= int(0.95 * grid.n)


def test_high_latitude_cold_dry_land_prefers_tundra_over_desert():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("terrain.elevation_m", np.full(grid.n, 100.0))

    temp_c = np.full(grid.n, 7.0)
    precip = np.full(grid.n, 240.0)
    polar = np.abs(grid.lat) >= 60.0
    temp_c[polar] = -2.0
    precip[polar] = 80.0
    seasonal_temp_c = np.repeat(temp_c[None, :], 4, axis=0)
    seasonal_temp_c[0, polar & (grid.lat > 0.0)] = -10.0
    seasonal_temp_c[2, polar & (grid.lat < 0.0)] = -10.0
    world.set_field("climate.surface_temperature", temp_c + 273.15)
    world.set_field("climate.seasonal_temperature", seasonal_temp_c + 273.15)
    world.set_field("climate.precipitation", precip)
    world.set_field(
        "climate.seasonal_precipitation",
        np.repeat(precip[None, :], 4, axis=0),
    )

    biome = BiosphereModule()._biomes(world, np.zeros(grid.n))

    assert int(polar.sum()) >= 8
    assert int(((biome == 5) & polar).sum()) >= int(0.9 * polar.sum())
    assert int(((biome == 2) & polar).sum()) == 0


def test_cool_high_elevation_dry_land_prefers_alpine_over_desert():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("terrain.elevation_m", np.full(grid.n, 2600.0))

    temp_c = np.full(grid.n, 8.0)
    precip = np.full(grid.n, 130.0)
    seasonal_temp_c = np.repeat(temp_c[None, :], 4, axis=0)
    seasonal_temp_c[0] = -2.0
    world.set_field("climate.surface_temperature", temp_c + 273.15)
    world.set_field("climate.seasonal_temperature", seasonal_temp_c + 273.15)
    world.set_field("climate.precipitation", precip)
    world.set_field(
        "climate.seasonal_precipitation",
        np.repeat(precip[None, :], 4, axis=0),
    )

    biome = BiosphereModule()._biomes(world, np.zeros(grid.n))

    assert int((biome == 5).sum()) >= int(0.95 * grid.n)
    assert int((biome == 2).sum()) == 0


def test_seasonal_tropical_climate_remains_tropical_when_annual_rain_is_high():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("terrain.elevation_m", np.full(grid.n, 100.0))

    temp_c = np.full(grid.n, 25.0)
    precip = np.full(grid.n, 1120.0)
    seasonal_precip = np.repeat(precip[None, :], 4, axis=0)
    seasonal_precip[0] = 250.0
    seasonal_precip[1] = 850.0
    seasonal_precip[2] = 1530.0
    seasonal_precip[3] = 1850.0
    world.set_field("climate.surface_temperature", temp_c + 273.15)
    world.set_field(
        "climate.seasonal_temperature",
        np.repeat((temp_c + 273.15)[None, :], 4, axis=0),
    )
    world.set_field("climate.precipitation", seasonal_precip.mean(axis=0))
    world.set_field("climate.seasonal_precipitation", seasonal_precip)

    biome = BiosphereModule()._biomes(world, np.zeros(grid.n))

    assert int((biome == 6).sum()) >= int(0.95 * grid.n)


def test_drier_seasonal_tropical_edge_can_be_grassland():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec)
    world.set_g("ocean.sea_level_m", 0.0)
    world.set_field("terrain.elevation_m", np.full(grid.n, 100.0))

    temp_c = np.full(grid.n, 25.0)
    seasonal_precip = np.full((4, grid.n), 1120.0)
    seasonal_precip[0] = 120.0
    seasonal_precip[1] = 650.0
    seasonal_precip[2] = 1200.0
    seasonal_precip[3] = 1700.0
    world.set_field("climate.surface_temperature", temp_c + 273.15)
    world.set_field(
        "climate.seasonal_temperature",
        np.repeat((temp_c + 273.15)[None, :], 4, axis=0),
    )
    world.set_field("climate.precipitation", seasonal_precip.mean(axis=0))
    world.set_field("climate.seasonal_precipitation", seasonal_precip)

    biome = BiosphereModule()._biomes(world, np.zeros(grid.n))

    assert int((biome == 3).sum()) >= int(0.95 * grid.n)
