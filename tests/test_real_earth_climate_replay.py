import numpy as np

from aevum.core.grid import SphereGrid
from aevum.diagnostics.real_earth_climate_replay import (
    RealEarthClimateReplayConfig,
    run_real_earth_climate_replay,
)
from aevum.core.units import CONSTANTS


def _write_reference_npz(path, cells=96):
    grid = SphereGrid.fibonacci(cells, CONSTANTS.EARTH_RADIUS)
    lat = grid.lat
    lon = grid.lon
    land = (np.sin(np.radians(lon * 1.7)) + 0.55 * np.cos(np.radians(lat * 2.0))) > 0.15
    elevation = np.where(
        land,
        320.0 + 900.0 * np.maximum(np.cos(np.radians(lat)), 0.0),
        -3600.0 + 850.0 * np.cos(np.radians(lat)) ** 2,
    ).astype(np.float32)

    annual_land_temp = (25.0 - 0.45 * np.abs(lat) - 0.0045 * np.maximum(elevation, 0.0)).astype(np.float32)
    annual_sst = (27.0 - 0.34 * np.abs(lat)).astype(np.float32)
    annual_precip = np.where(
        land,
        350.0 + 1100.0 * np.exp(-(lat / 24.0) ** 2),
        900.0 + 800.0 * np.exp(-(lat / 30.0) ** 2),
    ).astype(np.float32)

    seasonal_temp = np.stack([
        annual_land_temp - 6.0 * np.sin(np.radians(lat)),
        annual_land_temp,
        annual_land_temp + 6.0 * np.sin(np.radians(lat)),
        annual_land_temp,
    ]).astype(np.float32)
    seasonal_sst = np.stack([
        annual_sst - 1.6 * np.sin(np.radians(lat)),
        annual_sst,
        annual_sst + 1.6 * np.sin(np.radians(lat)),
        annual_sst,
    ]).astype(np.float32)
    seasonal_precip = np.stack([
        annual_precip * np.where(lat < 0.0, 1.35, 0.75),
        annual_precip,
        annual_precip * np.where(lat > 0.0, 1.35, 0.75),
        annual_precip,
    ]).astype(np.float32)

    wind = np.zeros((4, cells, 2), dtype=np.float32)
    for season in range(4):
        wind[season, :, 0] = np.where(np.abs(lat) < 28.0, -4.0, 6.0)
        wind[season, :, 1] = (0.8 if season == 2 else -0.8) * np.sin(np.radians(lat))
    slp_anomaly = np.stack([
        0.08 * lat,
        0.03 * lat,
        -0.08 * lat,
        -0.03 * lat,
    ]).astype(np.float32)
    currents_uv = np.zeros((cells, 2), dtype=np.float32)
    currents_uv[:, 0] = np.where(~land, 0.12 * np.cos(np.radians(lat)), 0.0)
    current_speed = np.linalg.norm(currents_uv, axis=1).astype(np.float32)

    np.savez_compressed(
        path,
        lat=lat,
        lon=lon,
        cell_area=grid.cell_area,
        earth__elevation_m=elevation,
        earth__land_mask=land,
        earth__annual_temperature_C=annual_land_temp,
        earth__annual_sst_C=annual_sst,
        earth__annual_precip_mm=annual_precip,
        earth__seasonal_temperature_C=seasonal_temp,
        earth__seasonal_sst_C=seasonal_sst,
        earth__seasonal_precip_mm_yr_equiv=seasonal_precip,
        earth__annual_sea_ice_concentration_pct=np.where(np.abs(lat) > 65.0, 45.0, 0.0).astype(np.float32),
        earth__seasonal_wind_u10_v10=wind,
        earth__seasonal_slp_anomaly_hPa=slp_anomaly,
        earth__annual_surface_current_speed_m_s=current_speed,
        earth__biome_class_proxy=np.where(land, 4, 0).astype(np.int16),
    )


def test_real_earth_climate_replay_writes_comparable_summary(tmp_path):
    reference = tmp_path / "earth_reference_96cells.npz"
    outdir = tmp_path / "replay"
    _write_reference_npz(reference)

    summary = run_real_earth_climate_replay(
        RealEarthClimateReplayConfig(
            earth_reference_npz=reference,
            outdir=outdir,
            seed=11,
            render_assets=False,
        )
    )

    assert summary["schema"] == "aevum.real_earth_climate_replay.v1"
    assert summary["cells"] == 96
    assert summary["land_mask_mismatch_area_fraction"] == 0.0
    assert summary["validation"]["passed"]
    assert "geographic circulation index is absent" not in summary["validation"]["warnings"]
    assert (outdir / "terminal_climate_arrays.npz").exists()
    assert (outdir / "real_earth_climate_replay_summary.json").exists()
    with np.load(outdir / "terminal_climate_arrays.npz", allow_pickle=False) as arrays:
        assert arrays["atmosphere__pressure_center_support"].shape == (4, 96)
        assert arrays["atmosphere__pressure_center_id"].shape == (4, 96)
        assert arrays["atmosphere__stationary_wave_pressure_support"].shape == (4, 96)
        assert arrays["atmosphere__pressure_genesis_source"].shape == (4, 96)
        assert arrays["atmosphere__pressure_genesis_wave_transfer"].shape == (4, 96)
        assert arrays["atmosphere__ocean_pressure_low_source_support"].shape == (4, 96)
        assert arrays["atmosphere__ocean_pressure_high_source_support"].shape == (4, 96)
        assert arrays["atmosphere__land_pressure_source_support"].shape == (4, 96)
        assert arrays["atmosphere__terrain_pressure_wave_source_support"].shape == (4, 96)
        assert arrays["climate__seasonal_insolation_anomaly"].shape == (4, 96)
        assert arrays["climate__surface_heat_capacity_class"].shape == (96,)
        assert arrays["climate__land_thermal_anomaly"].shape == (4, 96)
        assert arrays["climate__ocean_mixed_layer_thermal_anomaly"].shape == (4, 96)
        assert arrays["climate__elevation_lapse_cooling"].shape == (96,)
        assert arrays["climate__snow_ice_albedo_support"].shape == (4, 96)
        assert arrays["climate__sst_gradient_support"].shape == (4, 96)
        assert arrays["climate__same_latitude_sst_anomaly"].shape == (4, 96)
        assert arrays["climate__land_sea_thermal_contrast"].shape == (4, 96)
        assert arrays["ocean__basin_id"].shape == (96,)
        assert np.all(np.isfinite(arrays["atmosphere__pressure_center_support"]))
        assert np.all(np.isfinite(
            arrays["atmosphere__stationary_wave_pressure_support"]))
        assert np.all(np.isfinite(arrays["atmosphere__pressure_genesis_source"]))
        assert np.all(np.isfinite(
            arrays["atmosphere__pressure_genesis_wave_transfer"]))
        assert np.all(np.isfinite(
            arrays["atmosphere__ocean_pressure_low_source_support"]))
        assert np.all(np.isfinite(arrays["climate__seasonal_insolation_anomaly"]))
        assert np.all(np.isfinite(arrays["climate__land_sea_thermal_contrast"]))
        ocean = arrays["terrain__elevation_m"] < 0.0
        basin = arrays["ocean__basin_id"].astype(int)
        assert np.all(basin[~ocean] < 0)
        assert np.all((basin[ocean] >= 0) & (basin[ocean] <= 4))
        assert np.unique(basin[ocean]).size >= 3
    assert np.isfinite(summary["residuals"]["surface_temperature_mae_C"])
    assert np.isfinite(
        summary["residuals"]["earth_surface_temperature_max_adjacent_lat_band_delta_C"]
    )
    assert np.isfinite(
        summary["residuals"]["replay_surface_temperature_max_adjacent_lat_band_delta_C"]
    )
    assert np.isfinite(summary["dynamics_residuals"]["annual_wind_speed_global_mae_m_s"])
