import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.diagnostics.real_earth_wind_replay import (
    RealEarthWindReplayConfig,
    _uv_to_tangent,
    run_real_earth_wind_replay,
)


def test_real_earth_wind_replay_writes_maps_and_metrics(tmp_path):
    grid = SphereGrid.fibonacci(96, CONSTANTS.EARTH_RADIUS)
    lat = grid.lat
    lon = grid.lon
    land = lat > 5.0
    earth_uv = np.zeros((4, grid.n, 2), dtype=np.float32)
    earth_slp = np.zeros((4, grid.n), dtype=np.float32)
    for season in range(4):
        earth_uv[season, :, 0] = np.where(np.abs(lat) < 30.0, -4.0, 7.0)
        earth_uv[season, :, 1] = (season - 1.5) * 0.4 * np.cos(np.radians(lat))
        earth_slp[season] = (season - 1.5) * 1.5 + 4.0 * np.sin(np.radians(lat))
    earth = tmp_path / "earth_reference.npz"
    np.savez(
        earth,
        lat=lat,
        lon=lon,
        cell_area=grid.cell_area,
        earth__elevation_m=np.where(land, 200.0, -3000.0),
        earth__seasonal_wind_u10_v10=earth_uv,
        earth__seasonal_slp_anomaly_hPa=earth_slp,
    )

    replay_wind = 1.08 * _uv_to_tangent(lat, lon, earth_uv)
    replay_pressure = earth_slp / 4.0
    replay_pressure_support = np.clip(np.abs(replay_pressure) / 8.0, 0.0, 1.0)
    replay_stationary_support = 0.75 * replay_pressure_support
    replay_pressure_source = 0.25 * np.sign(replay_pressure) * replay_pressure_support
    replay_pressure_transfer = 0.5 * replay_pressure_source
    seasonal = np.stack([
        -6.0 * np.sin(np.radians(lat)),
        np.zeros(grid.n),
        6.0 * np.sin(np.radians(lat)),
        np.zeros(grid.n),
    ]).astype(np.float32)
    ocean_anomaly = np.where(~land[None, :], 0.35 * seasonal, 0.0)
    land_anomaly = np.where(land[None, :], seasonal, 0.0)
    replay = tmp_path / "replay_arrays.npz"
    np.savez(
        replay,
        terrain__elevation_m=np.where(land, 200.0, -3000.0),
        climate__coast_distance=np.clip(np.abs(lat) / 90.0, 0.0, 1.0),
        terrain__barrier_index=np.where(land, np.clip((lat - 5.0) / 60.0, 0.0, 1.0), 0.0),
        terrain__wind_gap_index=np.where(land, 0.25, 0.0),
        ocean__basin_id=np.where(land, -1.0, (lon > 0.0).astype(float)),
        ocean__shelf_index=np.where(~land, 0.4, 0.0),
        ocean__strait_index=np.where(~land & (np.abs(lon) < 20.0), 0.8, 0.0),
        climate__seasonal_insolation_anomaly=2.0 * seasonal,
        climate__surface_heat_capacity_class=np.where(land, 0.35, 1.0),
        climate__land_thermal_anomaly=land_anomaly,
        climate__ocean_mixed_layer_thermal_anomaly=ocean_anomaly,
        climate__elevation_lapse_cooling=np.where(land, 1.4, 0.0),
        climate__snow_ice_albedo_support=np.repeat(
            np.clip(np.abs(lat)[None, :] / 90.0, 0.0, 1.0), 4, axis=0),
        climate__sst_gradient_support=np.where(~land[None, :], 0.5, 0.0),
        climate__same_latitude_sst_anomaly=ocean_anomaly,
        climate__land_sea_thermal_contrast=land_anomaly - ocean_anomaly,
        atmosphere__seasonal_wind=replay_wind,
        atmosphere__land_sea_pressure_proxy=replay_pressure,
        atmosphere__pressure_center_support=replay_pressure_support,
        atmosphere__stationary_wave_pressure_support=replay_stationary_support,
        atmosphere__pressure_genesis_source=replay_pressure_source,
        atmosphere__pressure_genesis_wave_transfer=replay_pressure_transfer,
        atmosphere__ocean_pressure_low_source_support=np.where(
            ~land[None, :], replay_pressure_support, 0.0),
        atmosphere__ocean_pressure_high_source_support=np.where(
            ~land[None, :], 0.5 * replay_pressure_support, 0.0),
        atmosphere__land_pressure_source_support=np.where(
            land[None, :], replay_pressure_support, 0.0),
        atmosphere__terrain_pressure_wave_source_support=np.where(
            land[None, :], 0.25 * replay_pressure_support, 0.0),
    )

    outdir = tmp_path / "out"
    summary = run_real_earth_wind_replay(
        RealEarthWindReplayConfig(
            earth_reference_npz=earth,
            replay_arrays_npz=replay,
            outdir=outdir,
        )
    )

    assert summary["schema"] == "aevum.real_earth_wind_replay.v1"
    assert summary["metrics"]["seasonal_speed_mae_m_s"] > 0.0
    assert summary["metrics"]["direction_cosine_p50"] > 0.99
    assert summary["metrics"]["speed_zonal_anomaly_corr_all"] > 0.99
    assert np.isfinite(summary["metrics"]["pressure_standardized_mae_all"])
    assert summary["metrics"]["pressure_zonal_anomaly_corr_all"] > 0.99
    assert "pressure_standardized_delta_seasons" in summary["assets"]
    assert "pressure_contact_sheet" in summary["assets"]
    assert "earth_pressure_zonal_anomaly_seasons" in summary["assets"]
    assert "replay_pressure_center_support_seasons" in summary["assets"]
    assert "replay_stationary_wave_pressure_support_seasons" in summary["assets"]
    assert "replay_pressure_genesis_source_seasons" in summary["assets"]
    assert "replay_pressure_genesis_wave_transfer_seasons" in summary["assets"]
    assert "replay_ocean_pressure_low_source_support_seasons" in summary["assets"]
    assert "replay_m0_boundary_support_contact_sheet" in summary["assets"]
    assert "replay_m1_energy_support_contact_sheet" in summary["assets"]
    assert (outdir / "real_earth_wind_replay_summary.json").exists()
    assert (outdir / "real_earth_wind_replay_seasonal_metrics.csv").exists()
    assert (outdir / "real_earth_wind_replay_contact_sheet.png").exists()
    assert (outdir / "real_earth_pressure_replay_contact_sheet.png").exists()
    assert (outdir / "replay_m0_boundary_support_contact_sheet.png").exists()
    assert (outdir / "replay_m1_energy_support_contact_sheet.png").exists()
    assert (outdir / "earth_wind_speed_zonal_anomaly_seasons.png").exists()
    assert (outdir / "pressure_standardized_delta_seasons.png").exists()
    assert (outdir / "pressure_zonal_anomaly_delta_seasons.png").exists()
    assert (outdir / "replay_pressure_center_support_seasons.png").exists()
    assert (outdir / "replay_stationary_wave_pressure_support_seasons.png").exists()
    assert (outdir / "replay_pressure_genesis_source_seasons.png").exists()
    assert (outdir / "replay_pressure_genesis_wave_transfer_seasons.png").exists()
    assert (outdir / "replay_ocean_pressure_low_source_support_seasons.png").exists()
