import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.diagnostics.earth_climate_reference import (
    SOURCE_MANIFEST_SCHEMA,
    build_source_manifest,
    earth_elevation_metrics,
    esa_cci_land_cover_available,
    esa_cci_land_cover_metrics,
    etopo2022_available,
    etopo2022_crosscheck_metrics,
    gloh2o_koppen_available,
    load_etopo2022_opendap_ascii,
    sample_koppen_to_grid,
    sample_etopo2022_to_grid,
    sample_etopo5_to_grid,
    sample_noaa_psl_to_grid,
    sample_resolve_ecoregions_to_grid,
    sample_worldclim_to_grid,
    noaa_psl_metrics,
    noaa_aoml_drifter_current_available,
    noaa_aoml_drifter_current_metrics,
    noaa_oisst_available,
    noaa_oisst_metrics,
    oscar_monthly_current_available,
    oscar_monthly_current_metrics,
    resolve_ecoregion_metrics,
    resolve_ecoregions_available,
    sample_oscar_monthly_current_to_grid,
    sample_esa_cci_land_cover_to_grid,
    sample_noaa_aoml_drifter_current_to_grid,
    sample_noaa_oisst_to_grid,
    worldclim_metrics,
)


def test_earth_climate_source_manifest_tracks_required_reference_layers():
    manifest = build_source_manifest()
    assert manifest["schema"] == SOURCE_MANIFEST_SCHEMA
    source_ids = {entry["source_id"] for entry in manifest["entries"]}
    assert {
        "NOAA_ETOPO5_LOCAL",
        "NOAA_ETOPO_2022",
        "WORLDCLIM_2_1",
        "NOAA_PSL_NCEP_NCAR_REANALYSIS_1",
        "NASA_JPL_OSCAR",
        "KOTTEK_RUBEL_KOPPEN_2006_ASCII",
        "GLOH2O_KOPPEN_GEIGER",
        "ESA_CCI_LAND_COVER",
        "RESOLVE_ECOREGIONS_2017",
        "NOAA_OISST_V2_LTM_1991_2020",
        "NOAA_AOML_DRIFTER_CURRENT_CLIMATOLOGY_V3",
    }.issubset(source_ids)
    local = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "NOAA_ETOPO5_LOCAL"
    )
    assert local["local_cache_exists"]
    assert local["parser_status"] == "implemented"
    assert local["checksum_status"] == "computed"
    etopo2022 = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "NOAA_ETOPO_2022"
    )
    assert etopo2022["parser_status"] == "implemented"
    assert etopo2022["local_cache_exists"]
    assert etopo2022["local_cache_file_count"] >= 1
    worldclim = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "WORLDCLIM_2_1"
    )
    assert worldclim["parser_status"] == "implemented"
    assert worldclim["local_cache_exists"]
    assert worldclim["checksum_status"] == "file_checksums_computed"
    koppen = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "KOTTEK_RUBEL_KOPPEN_2006_ASCII"
    )
    assert koppen["parser_status"] == "implemented"
    assert koppen["local_cache_exists"]
    gloh2o = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "GLOH2O_KOPPEN_GEIGER"
    )
    assert gloh2o["parser_status"] == "implemented"
    assert gloh2o["local_cache_exists"]
    resolve = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "RESOLVE_ECOREGIONS_2017"
    )
    assert resolve["parser_status"] == "implemented"
    assert resolve["local_cache_exists"]
    ncep = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "NOAA_PSL_NCEP_NCAR_REANALYSIS_1"
    )
    assert ncep["parser_status"] == "implemented"
    assert ncep["local_cache_exists"]
    oisst = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "NOAA_OISST_V2_LTM_1991_2020"
    )
    assert oisst["parser_status"] == "implemented"
    assert oisst["local_cache_exists"]
    currents = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "NOAA_AOML_DRIFTER_CURRENT_CLIMATOLOGY_V3"
    )
    assert currents["parser_status"] == "implemented"
    oscar = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "NASA_JPL_OSCAR"
    )
    assert oscar["parser_status"] == "implemented"
    assert not oscar["requires_account"]
    land_cover = next(
        entry for entry in manifest["entries"]
        if entry["source_id"] == "ESA_CCI_LAND_COVER"
    )
    assert land_cover["parser_status"] == "implemented"
    assert not land_cover["requires_account"]


def test_etopo5_samples_to_aevum_grid_with_plausible_hypsometry():
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    elevation = sample_etopo5_to_grid(grid)
    assert elevation.shape == (grid.n,)
    assert np.isfinite(elevation).all()
    assert float(elevation.min()) < -5000.0
    assert float(elevation.max()) > 1000.0

    metrics = earth_elevation_metrics(grid, elevation)
    assert 0.20 <= metrics["land_fraction"] <= 0.40
    assert metrics["ocean_depth_mean_m"] > 2500.0
    assert metrics["lowland_fraction_lt500m"] < metrics["lowland_fraction_lt1000m"]
    ocean_partition = (
        metrics["shelf_fraction_of_ocean"]
        + metrics["slope_rise_fraction_of_ocean"]
        + metrics["abyss_fraction_of_ocean"]
        + metrics["trench_and_hadal_fraction_of_ocean"]
    )
    np.testing.assert_allclose(ocean_partition, 1.0)


def test_etopo2022_opendap_ascii_parser_reads_small_fixture(tmp_path):
    path = tmp_path / "mini_etopo2022.asc"
    path.write_text(
        "\n".join([
            "Dataset {",
            "} test;",
            "---------------------------------------------",
            "lat[2]",
            "-0.5, 0.5",
            "",
            "lon[3]",
            "-1.0, 0.0, 1.0",
            "",
            "z.z[2][3]",
            "[0], -100.0, 5.0, 25.0",
            "[1], -200.0, 15.0, 35.0",
            "",
            "z.lat[2]",
            "-0.5, 0.5",
            "",
            "z.lon[3]",
            "-1.0, 0.0, 1.0",
            "",
        ])
    )
    raster, lat, lon = load_etopo2022_opendap_ascii(path)
    np.testing.assert_allclose(lat, [-0.5, 0.5])
    np.testing.assert_allclose(lon, [-1.0, 0.0, 1.0])
    np.testing.assert_allclose(
        raster,
        [[-100.0, 5.0, 25.0], [-200.0, 15.0, 35.0]],
    )


def test_etopo2022_crosscheck_samples_modern_relief_against_etopo5():
    assert etopo2022_available()
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    etopo5 = sample_etopo5_to_grid(grid)
    etopo2022 = sample_etopo2022_to_grid(grid)

    assert etopo2022.shape == (grid.n,)
    assert np.isfinite(etopo2022).all()
    assert float(etopo2022.min()) < -5000.0
    assert float(etopo2022.max()) > 1000.0

    metrics = earth_elevation_metrics(grid, etopo2022)
    assert 0.20 <= metrics["land_fraction"] <= 0.40
    assert metrics["ocean_depth_mean_m"] > 2500.0

    crosscheck = etopo2022_crosscheck_metrics(grid, etopo5, etopo2022)
    assert crosscheck["valid_area_fraction"] > 0.99
    assert abs(crosscheck["land_fraction_delta"]) < 0.08
    assert crosscheck["mean_abs_elevation_delta_m"] > 1.0
    assert crosscheck["land_ocean_class_mismatch_area_fraction"] < 0.08


def test_worldclim_samples_monthly_temperature_and_precipitation_to_grid():
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    fields = sample_worldclim_to_grid(grid)

    assert fields["monthly_temperature_C"].shape == (12, grid.n)
    assert fields["monthly_precip_mm"].shape == (12, grid.n)
    assert fields["seasonal_temperature_C"].shape == (4, grid.n)
    assert fields["seasonal_precip_mm_yr_equiv"].shape == (4, grid.n)

    land = np.isfinite(fields["annual_temperature_C"])
    precip_land = np.isfinite(fields["annual_precip_mm"])
    assert np.array_equal(np.isfinite(fields["seasonal_precip_mm_yr_equiv"]).any(axis=0), precip_land)
    assert np.array_equal(np.isfinite(fields["dry_month_count"]), precip_land)
    assert 0.20 <= float(land.mean()) <= 0.40
    assert -60.0 < float(np.nanmin(fields["annual_temperature_C"])) < 0.0
    assert 15.0 < float(np.nanmax(fields["annual_temperature_C"])) < 40.0
    assert float(np.nanmax(fields["annual_precip_mm"])) > 1500.0

    metrics = worldclim_metrics(grid, fields)
    assert 0.20 <= metrics["worldclim_land_area_fraction"] <= 0.40
    assert 0.0 <= metrics["land_annual_temperature_mean_C"] <= 20.0
    assert 400.0 <= metrics["land_annual_precip_mean_mm"] <= 1200.0
    assert metrics["land_dry_month_count_p90"] >= metrics["land_dry_month_count_p50"]


def test_koppen_reference_samples_major_classes_and_biome_proxy():
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    elevation = sample_etopo5_to_grid(grid)
    fields = sample_koppen_to_grid(grid, elevation_m=elevation)

    koppen = fields["koppen_class"]
    major = fields["koppen_major_class"]
    biome = fields["biome_class_proxy"]
    assert koppen.shape == (grid.n,)
    assert major.shape == (grid.n,)
    assert biome.shape == (grid.n,)
    assert fields["koppen_source_id"] in {
        "GLOH2O_KOPPEN_GEIGER",
        "KOTTEK_RUBEL_KOPPEN_2006_ASCII",
    }
    if gloh2o_koppen_available():
        assert fields["koppen_source_id"] == "GLOH2O_KOPPEN_GEIGER"
        assert fields["koppen_period"] == "1991_2020"
        assert int(np.nanmax(koppen)) <= 30
    assert set(np.unique(major)).issubset({0, 1, 2, 3, 4, 5})
    assert {1, 2, 3, 4, 5}.issubset(set(int(x) for x in np.unique(major)))
    assert set(np.unique(biome)).issubset({0, 1, 2, 3, 4, 5, 6})
    assert int((biome == 6).sum()) > 0
    assert int((biome == 2).sum()) > 0
    assert int((biome == 1).sum()) > 0


def test_resolve_ecoregions_sample_terrestrial_biomes_to_grid():
    assert resolve_ecoregions_available()
    grid = SphereGrid.fibonacci(520, CONSTANTS.EARTH_RADIUS)
    elevation = sample_etopo5_to_grid(grid)
    fields = sample_resolve_ecoregions_to_grid(grid, land_mask=elevation >= 0.0)

    biome = fields["resolve_biome_class"]
    ecoregion = fields["resolve_ecoregion_id"]
    assert biome.shape == (grid.n,)
    assert ecoregion.shape == (grid.n,)
    assert set(int(x) for x in np.unique(biome)).issubset(set(range(15)))
    assert int((biome > 0).sum()) > 40
    assert int((ecoregion > 0).sum()) == int((biome > 0).sum())
    metrics = resolve_ecoregion_metrics(
        grid,
        fields,
        land_mask=elevation >= 0.0,
    )
    assert metrics["resolve_assigned_land_area_fraction"] > 0.65
    assert metrics["resolve_observed_biome_count"] >= 6


def test_noaa_psl_samples_seasonal_wind_and_pressure_to_grid():
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    fields = sample_noaa_psl_to_grid(grid)

    assert fields["monthly_wind_u10_v10"].shape == (12, grid.n, 2)
    assert fields["seasonal_wind_u10_v10"].shape == (4, grid.n, 2)
    assert fields["monthly_slp_hPa"].shape == (12, grid.n)
    assert fields["seasonal_slp_hPa"].shape == (4, grid.n)
    assert fields["seasonal_slp_anomaly_hPa"].shape == (4, grid.n)
    assert np.isfinite(fields["seasonal_wind_u10_v10"]).all()
    assert np.isfinite(fields["seasonal_slp_anomaly_hPa"]).all()

    speed = np.linalg.norm(fields["seasonal_wind_u10_v10"], axis=2)
    assert 1.0 <= float(np.percentile(speed, 50)) <= 8.0
    assert 4.0 <= float(np.percentile(speed, 90)) <= 15.0
    metrics = noaa_psl_metrics(fields)
    assert 4.0 <= metrics["seasonal_wind_speed_p90_m_s"] <= 15.0
    assert 1.0 <= metrics["seasonal_slp_anomaly_abs_p90_hPa"] <= 10.0


def test_noaa_oisst_samples_sst_and_sea_ice_to_grid():
    assert noaa_oisst_available()
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    elevation = sample_etopo5_to_grid(grid)
    fields = sample_noaa_oisst_to_grid(grid, ocean_mask=elevation < 0.0)

    assert fields["monthly_sst_C"].shape == (12, grid.n)
    assert fields["seasonal_sst_C"].shape == (4, grid.n)
    assert fields["annual_sst_C"].shape == (grid.n,)
    assert fields["monthly_sea_ice_concentration_pct"].shape == (12, grid.n)
    valid = np.isfinite(fields["annual_sst_C"])
    assert 0.55 <= float(valid.mean()) <= 0.80
    assert -3.0 <= float(np.nanmin(fields["annual_sst_C"])) <= 5.0
    assert 25.0 <= float(np.nanmax(fields["annual_sst_C"])) <= 35.0
    metrics = noaa_oisst_metrics(grid, fields)
    assert 10.0 <= metrics["annual_sst_mean_C"] <= 25.0
    assert 20.0 <= metrics["tropical_sst_mean_C"] <= 32.0
    assert metrics["annual_sea_ice_area_fraction_gt15pct"] > 0.01


def test_noaa_aoml_drifter_current_samples_current_speed_to_grid():
    assert noaa_aoml_drifter_current_available()
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    fields = sample_noaa_aoml_drifter_current_to_grid(grid)

    assert fields["surface_current_u_v"].shape == (grid.n, 2)
    assert fields["annual_surface_current_speed_m_s"].shape == (grid.n,)
    assert fields["annual_surface_current_direction_deg"].shape == (grid.n,)
    speed = fields["annual_surface_current_speed_m_s"]
    valid = np.isfinite(speed)
    assert 0.35 <= float(valid.mean()) <= 0.75
    assert 0.01 <= float(np.nanpercentile(speed, 50)) <= 0.30
    assert 0.10 <= float(np.nanpercentile(speed, 90)) <= 0.80
    metrics = noaa_aoml_drifter_current_metrics(grid, fields)
    assert 0.10 <= metrics["current_speed_p90_m_s"] <= 0.80
    assert metrics["swift_current_area_fraction_gt0_3m_s"] > 0.0


def test_oscar_monthly_current_climatology_samples_if_cache_available():
    assert oscar_monthly_current_available()
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    fields = sample_oscar_monthly_current_to_grid(grid)

    assert fields["monthly_surface_current_u_v"].shape == (12, grid.n, 2)
    assert fields["seasonal_surface_current_u_v"].shape == (4, grid.n, 2)
    assert fields["annual_surface_current_u_v"].shape == (grid.n, 2)
    assert fields["monthly_surface_current_speed_m_s"].shape == (12, grid.n)
    assert fields["seasonal_surface_current_speed_m_s"].shape == (4, grid.n)
    assert fields["annual_surface_current_speed_m_s"].shape == (grid.n,)
    speed = fields["annual_surface_current_speed_m_s"]
    valid = np.isfinite(speed)
    assert 0.35 <= float(valid.mean()) <= 0.80
    assert 0.01 <= float(np.nanpercentile(speed, 50)) <= 0.35
    assert 0.10 <= float(np.nanpercentile(speed, 90)) <= 0.90
    metrics = oscar_monthly_current_metrics(grid, fields)
    assert 0.10 <= metrics["annual_current_speed_p90_m_s"] <= 0.90
    assert metrics["monthly_current_speed_p90_m_s"] >= metrics["annual_current_speed_p90_m_s"]


def test_esa_cci_land_cover_preview_samples_broad_classes_if_cache_available():
    assert esa_cci_land_cover_available()
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    elevation = sample_etopo5_to_grid(grid)
    fields = sample_esa_cci_land_cover_to_grid(grid)

    lccs = fields["esa_cci_lccs_class"]
    broad = fields["esa_cci_land_cover_broad_class"]
    assert lccs.shape == (grid.n,)
    assert broad.shape == (grid.n,)
    assert int((lccs > 0).sum()) > 300
    assert {1, 2, 3, 4, 7, 8}.issubset(set(int(x) for x in np.unique(broad)))
    metrics = esa_cci_land_cover_metrics(grid, fields, elevation_m=elevation)
    assert metrics["valid_area_fraction"] > 0.90
    assert metrics["observed_lccs_class_count"] >= 10
    assert 0.50 <= metrics["water_area_fraction"] <= 0.80
    assert 0.05 <= metrics["forest_area_fraction"] <= 0.25
    assert 0.02 <= metrics["cropland_area_fraction"] <= 0.25
