import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.modules.climate import ClimateModule


def test_low_latitude_seasonality_redistribution_preserves_annual_mean():
    grid = SphereGrid.fibonacci(360, CONSTANTS.EARTH_RADIUS)
    seasonal = np.repeat(np.full((1, grid.n), 800.0), 4, axis=0)
    low_lat_land = np.abs(grid.lat) < 18.0
    seasonal[:, low_lat_land] = np.array([[1200.0], [900.0], [700.0], [400.0]])
    ocean = grid.lat > 45.0
    land = ~ocean

    out = ClimateModule()._accentuate_low_latitude_seasonality(
        grid, seasonal, land)

    assert np.allclose(out.mean(axis=0), seasonal.mean(axis=0))
    assert np.allclose(out[:, ocean], seasonal[:, ocean])
    before_amp = seasonal[:, low_lat_land].max(axis=0) - seasonal[:, low_lat_land].min(axis=0)
    after_amp = out[:, low_lat_land].max(axis=0) - out[:, low_lat_land].min(axis=0)
    assert float(np.median(after_amp)) > float(np.median(before_amp))


def test_moisture_flow_network_extracts_downwind_land_corridors():
    grid = SphereGrid.fibonacci(520, CONSTANTS.EARTH_RADIUS)
    module = ClimateModule()
    module._edge_cache(grid)
    east, _ = module._tangent_basis(grid)
    seasonal_wind = np.repeat((8.0 * east)[None, :, :], 4, axis=0)
    land = (np.abs(grid.lat) < 36.0) & (grid.lon > -80.0) & (grid.lon < 125.0)
    ocean = ~land
    shape = (4, grid.n)
    moisture = np.zeros(shape, dtype=float)
    moisture[:, land] = 0.72
    monsoon = np.zeros(shape, dtype=float)
    monsoon[:, land & (np.abs(grid.lat) < 24.0)] = 0.35
    storm = np.zeros(shape, dtype=float)
    storm[:, land & (np.abs(grid.lat) > 22.0)] = 0.28
    shadow = np.zeros(shape, dtype=float)
    response = np.ones(shape, dtype=float)
    precip = np.full(shape, 500.0, dtype=float)
    source_warmth = np.zeros(shape, dtype=float)
    source_warmth[:, ocean & (grid.lon < -80.0)] = 1.0
    hydro = {
        "moisture_access": moisture,
        "monsoon_rainfall_corridor": monsoon,
        "storm_track_rainfall_corridor": storm,
        "rain_shadow_index": shadow,
        "regional_precipitation_response": response,
        "seasonal_precipitation": precip,
        "source_ocean_warmth": source_warmth,
    }
    geography = {
        "terrain.barrier_index": np.zeros(grid.n, dtype=float),
        "terrain.wind_gap_index": np.zeros(grid.n, dtype=float),
        "ocean.basin_id": np.where(ocean, 1.0, -1.0),
        "climate.continent_id": np.where(land, 2.0, -1.0),
        "climate.coast_distance": np.zeros(grid.n, dtype=float),
    }

    out = module._seasonal_moisture_flow_networks(
        grid,
        hydro,
        land,
        ocean,
        seasonal_wind,
        geography,
        np.zeros(grid.n, dtype=float),
        np.zeros(grid.n, dtype=float),
    )

    assert out["source"].shape == shape
    assert out["pathway"].shape == shape
    assert out["source_basin_id"].shape == shape
    assert out["network_id"].shape == shape
    assert float(np.percentile(out["source"][:, ocean], 90)) > 0.0
    assert float(np.percentile(out["pathway"][:, land], 90)) > 0.0
    active_land = land & (out["pathway"][0] > 0.05)
    assert np.any(out["source_basin_id"][:, active_land] == 1.0)
    assert out["objects"]
    assert {str(obj["type"]) for obj in out["objects"]} == {"moisture_flow_network"}
    assert {str(obj["season"]) for obj in out["objects"]} <= {"DJF", "MAM", "JJA", "SON"}
    assert max(float(obj["area_fraction"]) for obj in out["objects"]) > 0.0
    assert {
        int(obj.get("dominant_source_basin_id", -1))
        for obj in out["objects"]
    } == {1}
    assert np.any(out["network_id"] >= 0.0)


def test_moisture_flow_precipitation_response_preserves_land_budget():
    grid = SphereGrid.fibonacci(520, CONSTANTS.EARTH_RADIUS)
    module = ClimateModule()
    module._edge_cache(grid)
    land = (np.abs(grid.lat) < 52.0) & (grid.lon > -130.0)
    ocean = ~land
    shape = (4, grid.n)
    seasonal = np.full(shape, 640.0, dtype=float)
    seasonal[:, land & (grid.lat > 12.0)] = 720.0
    pathway = np.zeros(shape, dtype=float)
    corridor = land & (np.abs(grid.lat) < 25.0) & (grid.lon > -35.0) & (grid.lon < 95.0)
    pathway[:, corridor] = 1.0
    network_id = np.full(shape, -1.0, dtype=float)
    network_id[:, corridor] = 1.0
    monsoon = np.zeros(shape, dtype=float)
    monsoon[:, corridor] = 0.45
    storm = np.zeros(shape, dtype=float)
    shadow = np.zeros(shape, dtype=float)
    shadow[:, land & ~corridor] = 0.22
    hydro = {
        "seasonal_precipitation": seasonal.copy(),
        "precipitation": seasonal.mean(axis=0),
        "evaporation": np.full(grid.n, 180.0, dtype=float),
        "monsoon_rainfall_corridor": monsoon,
        "storm_track_rainfall_corridor": storm,
        "rain_shadow_index": shadow,
        "regional_precipitation_response": np.ones(shape, dtype=float),
    }
    moisture_flow = {
        "pathway": pathway,
        "network_id": network_id,
        "source_basin_id": np.where(np.repeat(corridor[None, :], 4, axis=0), 3.0, -1.0),
    }

    out, diag = module._apply_moisture_flow_precipitation_response(
        grid, hydro, moisture_flow, land, ocean)
    response = out["moisture_flow_precipitation_response"]
    shaped = out["seasonal_precipitation"]

    assert diag["enabled"]
    assert response.shape == shape
    assert np.isfinite(response).all()
    assert float(np.percentile(response[:, land], 95)) > 1.03
    assert float(np.percentile(response[:, land], 5)) < 0.97
    assert np.allclose(shaped[:, ocean], seasonal[:, ocean])
    for season in range(4):
        before = np.average(seasonal[season, land], weights=grid.cell_area[land])
        after = np.average(shaped[season, land], weights=grid.cell_area[land])
        assert np.isclose(after, before, atol=1e-6)
    assert float(np.mean(shaped[:, corridor])) > float(np.mean(seasonal[:, corridor]))
    assert np.allclose(out["precipitation"], shaped.mean(axis=0))
    assert out["moisture_budget_region_id"].shape == shape


def test_moisture_flow_precipitation_response_preserves_continent_budgets():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    module = ClimateModule()
    module._edge_cache(grid)
    land_west = (np.abs(grid.lat) < 48.0) & (grid.lon < -35.0)
    land_east = (np.abs(grid.lat) < 48.0) & (grid.lon > 35.0)
    land = land_west | land_east
    ocean = ~land
    shape = (4, grid.n)
    seasonal = np.full(shape, 780.0, dtype=float)
    seasonal[:, land_west] = 520.0
    seasonal[:, land_east] = 1320.0

    corridor = land_west & (np.abs(grid.lat) < 24.0) & (grid.lon > -120.0)
    pathway = np.zeros(shape, dtype=float)
    pathway[:, corridor] = 1.0
    network_id = np.full(shape, -1.0, dtype=float)
    network_id[:, corridor] = 1.0
    monsoon = np.zeros(shape, dtype=float)
    monsoon[:, corridor] = 0.55
    storm = np.zeros(shape, dtype=float)
    shadow = np.zeros(shape, dtype=float)
    shadow[:, land_west & ~corridor] = 0.18
    continent_id = np.full(grid.n, -1.0, dtype=float)
    continent_id[land_west] = 10.0
    continent_id[land_east] = 20.0
    hydro = {
        "seasonal_precipitation": seasonal.copy(),
        "precipitation": seasonal.mean(axis=0),
        "evaporation": np.full(grid.n, 220.0, dtype=float),
        "monsoon_rainfall_corridor": monsoon,
        "storm_track_rainfall_corridor": storm,
        "rain_shadow_index": shadow,
        "regional_precipitation_response": np.ones(shape, dtype=float),
    }
    moisture_flow = {
        "pathway": pathway,
        "network_id": network_id,
    }
    geography = {
        "climate.continent_id": continent_id,
    }

    out, diag = module._apply_moisture_flow_precipitation_response(
        grid, hydro, moisture_flow, land, ocean, geography)
    shaped = out["seasonal_precipitation"]
    budget_id = out["moisture_budget_region_id"]

    assert diag["enabled"]
    assert diag["budget_region_count_p50"] >= 2.0
    assert diag["budget_sector_split_count_p50"] >= 1.0
    assert diag["max_budget_region_mean_delta_mm_yr"] < 1e-6
    assert budget_id.shape == shape
    assert len({int(x) for x in np.unique(budget_id[0, land]) if int(x) > 0}) >= 2
    assert set(np.unique(budget_id[:, ocean])) == {-1.0}
    assert float(np.mean(shaped[:, corridor])) > float(np.mean(seasonal[:, corridor]))
    for season in range(4):
        west_before = np.average(
            seasonal[season, land_west], weights=grid.cell_area[land_west])
        west_after = np.average(
            shaped[season, land_west], weights=grid.cell_area[land_west])
        east_before = np.average(
            seasonal[season, land_east], weights=grid.cell_area[land_east])
        east_after = np.average(
            shaped[season, land_east], weights=grid.cell_area[land_east])
        assert np.isclose(west_after, west_before, atol=1e-6)
        assert np.isclose(east_after, east_before, atol=1e-6)


def test_moisture_budget_regions_split_residual_by_source_basin():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    module = ClimateModule()
    module._edge_cache(grid)
    land = (np.abs(grid.lat) < 50.0) & (np.abs(grid.lon) < 140.0)
    shape = (4, grid.n)
    continent_id = np.full(grid.n, -1.0, dtype=float)
    continent_id[land] = 4.0
    network_id = np.full(shape, -1.0, dtype=float)
    network_id[:, land] = 1.0
    pathway = np.zeros(shape, dtype=float)
    pathway[:, land] = 0.65
    source_basin_id = np.full(shape, -1.0, dtype=float)
    west_source = land & (grid.lon < 0.0)
    east_source = land & (grid.lon >= 0.0)
    source_basin_id[:, west_source] = 7.0
    source_basin_id[:, east_source] = 8.0

    budget_id, meta = module._seasonal_moisture_budget_regions(
        grid,
        land,
        {"climate.continent_id": continent_id},
        {
            "network_id": network_id,
            "pathway": pathway,
            "source_basin_id": source_basin_id,
        },
    )

    assert budget_id.shape == shape
    assert meta["budget_sector_split_count_p50"] >= 2.0
    for season in range(4):
        west_regions = {int(x) for x in np.unique(budget_id[season, west_source]) if x > 0}
        east_regions = {int(x) for x in np.unique(budget_id[season, east_source]) if x > 0}
        assert west_regions
        assert east_regions
        assert west_regions.isdisjoint(east_regions)


def test_precipitation_response_regions_bind_source_budget_and_flow_ids():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    module = ClimateModule()
    module._edge_cache(grid)
    land = np.abs(grid.lat) < 58.0
    shape = (4, grid.n)
    wet = land & (np.abs(grid.lat) < 24.0) & (grid.lon > -45.0) & (grid.lon < 70.0)
    dry = land & (grid.lat > 8.0) & (grid.lon < -80.0)

    response = np.ones(shape, dtype=float)
    response[:, wet] = 1.09
    response[:, dry] = 0.88
    seasonal = np.full(shape, 640.0, dtype=float)
    seasonal[:, wet] = 820.0
    seasonal[:, dry] = 280.0
    budget_id = np.full(shape, -1.0, dtype=float)
    budget_id[:, wet] = 21.0
    budget_id[:, dry] = 22.0
    pathway = np.zeros(shape, dtype=float)
    pathway[:, wet] = 0.95
    pathway[:, dry] = 0.35
    network_id = np.full(shape, -1.0, dtype=float)
    network_id[:, wet] = 5.0
    network_id[:, dry] = 6.0
    source_basin_id = np.full(shape, -1.0, dtype=float)
    source_basin_id[:, wet] = 7.0
    source_basin_id[:, dry] = 8.0

    out = module._precipitation_response_region_objects(
        grid,
        {
            "seasonal_precipitation": seasonal,
            "moisture_flow_precipitation_response": response,
            "moisture_budget_region_id": budget_id,
            "monsoon_rainfall_corridor": np.repeat(
                np.where(wet[None, :], 0.5, 0.0), 4, axis=0),
            "storm_track_rainfall_corridor": np.repeat(
                np.where(dry[None, :], 0.3, 0.0), 4, axis=0),
            "rain_shadow_index": np.repeat(
                np.where(dry[None, :], 0.2, 0.0), 4, axis=0),
        },
        {
            "pathway": pathway,
            "network_id": network_id,
            "source_basin_id": source_basin_id,
        },
        land,
    )

    region_id = out["region_id"]
    objects = out["objects"]
    assert region_id.shape == shape
    assert np.isfinite(region_id).all()
    assert np.any(region_id[:, wet] > 0.0)
    assert np.any(region_id[:, dry] > 0.0)
    assert objects
    kinds = {str(obj["kind"]) for obj in objects}
    assert "wet_precipitation_response_region" in kinds
    assert "dry_precipitation_response_region" in kinds
    assert {str(obj["type"]) for obj in objects} == {"precipitation_response_region"}
    assert {int(obj["dominant_source_basin_id"]) for obj in objects} >= {7, 8}
    assert {int(obj["dominant_budget_region_id"]) for obj in objects} >= {21, 22}
    assert {int(obj["dominant_flow_network_id"]) for obj in objects} >= {5, 6}
    assert min(float(obj["source_basin_attributed_fraction"]) for obj in objects) > 0.95
    assert min(float(obj["budget_region_attributed_fraction"]) for obj in objects) > 0.95


def test_receiver_catchments_merge_response_regions_by_budget_and_source():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    module = ClimateModule()
    module._edge_cache(grid)
    land = np.abs(grid.lat) < 58.0
    shape = (4, grid.n)
    west = land & (grid.lon < -20.0)
    east = land & (grid.lon > 20.0)

    seasonal = np.full(shape, 620.0, dtype=float)
    response = np.ones(shape, dtype=float)
    response[:, west] = 1.08
    response[:, east] = 0.92
    budget_id = np.full(shape, -1.0, dtype=float)
    budget_id[:, land] = 11.0
    pathway = np.zeros(shape, dtype=float)
    pathway[:, west | east] = 0.7
    network_id = np.full(shape, -1.0, dtype=float)
    network_id[:, west] = 3.0
    network_id[:, east] = 4.0
    source_basin_id = np.full(shape, -1.0, dtype=float)
    source_basin_id[:, west] = 7.0
    source_basin_id[:, east] = 8.0
    ocean = ~land
    source_basin_id[:, ocean & (grid.lon < 0.0)] = 7.0
    source_basin_id[:, ocean & (grid.lon >= 0.0)] = 8.0
    source = np.zeros(shape, dtype=float)
    source[:, ocean] = 0.8
    precip_region_id = np.full(shape, -1.0, dtype=float)
    precip_region_id[:, west & (grid.lat > 0.0)] = 101.0
    precip_region_id[:, west & (grid.lat <= 0.0)] = 102.0
    precip_region_id[:, east & (grid.lat > 0.0)] = 201.0
    precip_region_id[:, east & (grid.lat <= 0.0)] = 202.0

    out = module._receiver_catchment_objects(
        grid,
        {
            "seasonal_precipitation": seasonal,
            "moisture_flow_precipitation_response": response,
            "moisture_budget_region_id": budget_id,
        },
        {
            "source": source,
            "source_basin_id": source_basin_id,
            "pathway": pathway,
            "network_id": network_id,
        },
        {"region_id": precip_region_id},
        land,
    )

    catchment_id = out["catchment_id"]
    objects = out["objects"]
    assert catchment_id.shape == shape
    assert np.isfinite(catchment_id).all()
    assert np.any(catchment_id[:, west] > 0.0)
    assert np.any(catchment_id[:, east] > 0.0)
    assert objects
    assert {str(obj["type"]) for obj in objects} == {"receiver_catchment"}
    assert {int(obj["dominant_source_basin_id"]) for obj in objects} >= {7, 8}
    source_objects = [
        obj for obj in objects if int(obj["dominant_source_basin_id"]) >= 0
    ]
    assert min(float(obj["source_basin_purity"]) for obj in source_objects) >= 0.55
    assert max(
        float(obj["precipitation_response_attributed_fraction"]) for obj in objects
    ) > 0.0

    accounting = module._source_basin_receiver_accounting(
        grid,
        {
            "seasonal_precipitation": seasonal,
            "moisture_flow_precipitation_response": response,
        },
        {
            "source": source,
            "source_basin_id": source_basin_id,
            "pathway": pathway,
        },
        out,
        land,
        ocean,
    )
    supply = accounting["source_basin_supply_index"]
    balance = accounting["receiver_catchment_supply_balance"]
    accounting_objects = accounting["objects"]
    assert supply.shape == shape
    assert balance.shape == shape
    assert np.isfinite(supply).all()
    assert np.isfinite(balance).all()
    assert float(np.percentile(supply[:, west | east], 50)) > 0.0
    assert float(np.percentile(balance[:, land], 50)) > 0.0
    assert max(
        float(obj.get("mean_source_basin_supply_index", 0.0))
        for obj in accounting_objects
    ) > 0.0
    assert max(
        float(obj.get("precipitation_supply_balance", 0.0))
        for obj in accounting_objects
    ) > 0.0


def test_receiver_supply_feedback_preserves_local_budgets_and_ocean():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    module = ClimateModule()
    module._edge_cache(grid)
    land_west = (np.abs(grid.lat) < 55.0) & (grid.lon < -20.0)
    land_east = (np.abs(grid.lat) < 55.0) & (grid.lon > 20.0)
    land = land_west | land_east
    ocean = ~land
    shape = (4, grid.n)
    seasonal = np.full(shape, 500.0, dtype=float)
    seasonal[:, land_west] = 620.0
    seasonal[:, land_east] = 920.0
    seasonal[:, ocean] = 700.0
    budget_id = np.full(shape, -1.0, dtype=float)
    budget_id[:, land_west] = 11.0
    budget_id[:, land_east] = 22.0
    west_core = land_west & (np.abs(grid.lat) < 24.0)
    west_rest = land_west & ~west_core
    c4f_response = np.ones(shape, dtype=float)
    c4f_response[:, west_core] = 1.05
    source_supply = np.zeros(shape, dtype=float)
    source_supply[:, west_core] = 1.1
    source_supply[:, west_rest] = 0.25
    source_supply[:, land_east] = 0.25
    receiver_balance = np.zeros(shape, dtype=float)
    receiver_balance[:, west_core] = 0.85
    receiver_balance[:, west_rest] = 0.35
    receiver_balance[:, land_east] = 0.35

    out, diag = module._apply_receiver_supply_precipitation_feedback(
        grid,
        {
            "seasonal_precipitation": seasonal.copy(),
            "precipitation": seasonal.mean(axis=0),
            "evaporation": np.full(grid.n, 220.0),
            "moisture_budget_region_id": budget_id,
            "moisture_flow_precipitation_response": c4f_response,
        },
        {
            "source_basin_supply_index": source_supply,
            "receiver_catchment_supply_balance": receiver_balance,
        },
        land,
    )

    feedback = out["receiver_supply_precipitation_feedback"]
    shaped = out["seasonal_precipitation"]
    assert diag["enabled"]
    assert feedback.shape == shape
    assert np.isfinite(feedback).all()
    assert float(np.min(feedback[:, land])) >= 0.90
    assert float(np.max(feedback[:, land])) <= 1.12
    assert np.allclose(feedback[:, ocean], 1.0)
    assert np.allclose(shaped[:, ocean], seasonal[:, ocean])
    assert float(np.mean(shaped[:, west_core])) > float(np.mean(seasonal[:, west_core]))
    assert float(np.mean(shaped[:, west_rest])) < float(np.mean(seasonal[:, west_rest]))
    for season in range(4):
        for mask in (land_west, land_east):
            before = np.average(seasonal[season, mask], weights=grid.cell_area[mask])
            after = np.average(shaped[season, mask], weights=grid.cell_area[mask])
            assert np.isclose(after, before, atol=1e-6)
    assert np.allclose(out["precipitation"], shaped.mean(axis=0))
    assert diag["max_budget_region_mean_delta_mm_yr"] < 1e-6


def test_receiver_catchments_cover_land_without_budget_region():
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    module = ClimateModule()
    module._edge_cache(grid)
    land = np.abs(grid.lat) < 52.0
    shape = (4, grid.n)

    seasonal = np.full(shape, 520.0, dtype=float)
    response = np.ones(shape, dtype=float)
    budget_id = np.full(shape, -1.0, dtype=float)
    budgeted = land & (grid.lon < 0.0)
    residual = land & (grid.lon >= 0.0)
    budget_id[:, budgeted] = 11.0
    source_basin_id = np.full(shape, -1.0, dtype=float)
    source_basin_id[:, land] = 3.0

    out = module._receiver_catchment_objects(
        grid,
        {
            "seasonal_precipitation": seasonal,
            "moisture_flow_precipitation_response": response,
            "moisture_budget_region_id": budget_id,
        },
        {"source_basin_id": source_basin_id},
        {"region_id": np.full(shape, -1.0, dtype=float)},
        land,
    )

    catchment_id = out["catchment_id"]
    objects = out["objects"]
    assert np.all(catchment_id[:, budgeted] > 0.0)
    assert np.all(catchment_id[:, residual] > 0.0)
    residual_objects = [
        obj for obj in objects if bool(obj.get("residual_budget_region"))
    ]
    assert residual_objects
    assert {int(obj["dominant_budget_region_id"]) for obj in residual_objects} == {0}
    assert min(
        float(obj["budget_region_attributed_fraction"]) for obj in residual_objects
    ) == 1.0
