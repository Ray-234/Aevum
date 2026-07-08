import numpy as np
from PIL import Image

from aevum.diagnostics.earth_climate_comparison import (
    EarthClimateComparisonConfig,
    _earth_reference_row_label,
    earth_reference_metrics,
    run_earth_climate_comparison,
)


def test_earth_reference_metrics_reports_optional_esa_cci_land_cover(tmp_path):
    path = tmp_path / "out_earth_climate_reference_r6_landcover_20260706"
    path.mkdir()
    npz = path / "earth_reference_4cells.npz"
    n = 4
    np.savez(
        npz,
        cell_area=np.ones(n),
        lat=np.array([-45.0, -15.0, 15.0, 45.0]),
        earth__land_mask=np.array([True, True, False, False]),
        earth__annual_temperature_C=np.array([8.0, 18.0, np.nan, np.nan]),
        earth__annual_sst_C=np.array([np.nan, np.nan, 20.0, 22.0]),
        earth__annual_precip_mm=np.array([250.0, 750.0, np.nan, np.nan]),
        earth__seasonal_temperature_C=np.tile(
            np.array([7.0, 17.0, np.nan, np.nan]),
            (4, 1),
        ),
        earth__seasonal_sst_C=np.tile(
            np.array([np.nan, np.nan, 19.0, 21.0]),
            (4, 1),
        ),
        earth__seasonal_precip_mm_yr_equiv=np.tile(
            np.array([250.0, 750.0, np.nan, np.nan]),
            (4, 1),
        ),
        earth__biome_class_proxy=np.array([2, 4, 0, 0]),
        earth__annual_surface_current_speed_m_s=np.array([np.nan, np.nan, 0.1, 0.2]),
        earth__esa_cci_land_cover_broad_class=np.array([3, 2, 1, 1], dtype=np.uint8),
    )

    metrics = earth_reference_metrics(npz)

    assert metrics["land_cover_valid_area_fraction"] == 1.0
    assert metrics["land_cover_water_area_fraction"] == 0.5
    assert metrics["land_cover_cropland_area_fraction"] == 0.25
    assert metrics["land_cover_forest_area_fraction"] == 0.25


def test_earth_reference_row_label_uses_reference_generation(tmp_path):
    npz = (
        tmp_path
        / "out_earth_climate_reference_r6_landcover_20260706"
        / "earth_reference_8000cells.npz"
    )

    assert _earth_reference_row_label(npz) == "Earth R6"


def test_comparison_renders_generated_previews_from_archived_arrays(tmp_path):
    earth_dir = tmp_path / "out_earth_climate_reference_r6_landcover_20260706"
    earth_dir.mkdir()
    earth_npz = earth_dir / "earth_reference_8cells.npz"
    n = 8
    lat = np.linspace(-70.0, 70.0, n)
    lon = np.linspace(-160.0, 160.0, n)
    land = np.array([True, True, True, False, False, True, False, False])
    ocean = ~land
    np.savez(
        earth_npz,
        cell_area=np.ones(n),
        lat=lat,
        earth__land_mask=land,
        earth__annual_temperature_C=np.where(land, np.linspace(5.0, 25.0, n), np.nan),
        earth__annual_sst_C=np.where(ocean, np.linspace(12.0, 24.0, n), np.nan),
        earth__annual_precip_mm=np.where(land, np.linspace(100.0, 1100.0, n), np.nan),
        earth__seasonal_temperature_C=np.tile(
            np.where(land, np.linspace(4.0, 24.0, n), np.nan),
            (4, 1),
        ),
        earth__seasonal_sst_C=np.tile(
            np.where(ocean, np.linspace(11.0, 23.0, n), np.nan),
            (4, 1),
        ),
        earth__seasonal_precip_mm_yr_equiv=np.tile(
            np.where(land, np.linspace(120.0, 1000.0, n), np.nan),
            (4, 1),
        ),
        earth__biome_class_proxy=np.where(land, 4, 0),
        earth__annual_surface_current_speed_m_s=np.where(ocean, 0.16, np.nan),
    )

    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    arrays = assets / "terminal_climate_arrays.npz"
    elevation = np.where(land, 120.0, -2000.0)
    seasonal_temp = np.vstack([
        np.linspace(275.0 + season, 298.0 + season, n)
        for season in range(4)
    ])
    seasonal_precip = np.vstack([
        np.linspace(200.0 + 20.0 * season, 900.0 + 20.0 * season, n)
        for season in range(4)
    ])
    precip = seasonal_precip.mean(axis=0)
    currents = np.zeros((n, 3), dtype=np.float64)
    currents[ocean, 0] = 0.12
    np.savez(
        arrays,
        cell_area=np.ones(n),
        lat=lat,
        lon=lon,
        sea_level_m=np.asarray([0.0]),
        terrain__elevation_m=elevation,
        climate__surface_temperature=np.linspace(276.0, 299.0, n),
        climate__seasonal_temperature=seasonal_temp,
        climate__precipitation=precip,
        climate__seasonal_precipitation=seasonal_precip,
        climate__precipitation_seasonality=seasonal_precip.max(axis=0) / precip,
        biosphere__biome=np.where(land, 4, 0),
        ocean__currents=currents,
    )
    terminal = tmp_path / "terminal_summary.json"
    terminal.write_text(
        """
{
  "summaries": [
    {
      "preset": "earthlike_mobile_lid",
      "seed": 1,
      "arrays": "__ARRAYS__",
      "assets_dir": "__ASSETS__"
    }
  ]
}
""".replace("__ARRAYS__", str(arrays)).replace("__ASSETS__", str(assets))
    )

    report = run_earth_climate_comparison(
        EarthClimateComparisonConfig(
            earth_reference_npz=earth_npz,
            terminal_summary_json=terminal,
            outdir=tmp_path / "comparison",
            render_contact_sheet=True,
        )
    )

    assert report["run_count"] == 1
    assert (tmp_path / "comparison" / "earth_vs_generated_climate_contact_sheet.png").exists()
    for filename in ("temperature.png", "precip.png", "biomes.png", "currents.png"):
        path = assets / filename
        assert path.exists()
        with Image.open(path) as image:
            assert image.size[0] > 100
            assert image.size[1] > 50
