import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_fit_report
from aevum.diagnostics.earth_climate_fitting import (
    SCHEMA,
    EarthClimateFittingConfig,
    run_earth_climate_fitting_report,
)


def _comparison(generated, earth):
    return {
        "generated": generated,
        "earth": earth,
        "delta": generated - earth,
        "normalized_delta": abs(generated - earth),
    }


def test_earth_climate_fitting_report_identifies_dry_hydroclimate(tmp_path):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    (assets / "summary.json").write_text(json.dumps({
        "climate_step_diagnostics": {
            "moisture_access_land_p75": 0.48,
            "land_monsoon_index_p90": 1.8,
        },
        "climate_diagnostics": {
            "circulation": {
                "moisture_access_land_p75": 0.48,
                "monsoon_potential_land_p90": 0.07,
                "monsoon_potential_land_p99": 0.22,
                "source_ocean_warmth_ocean_p75": 0.84,
                "terrain_blocking_land_p75": 0.27,
            },
            "precipitation": {
                "land_wet_fraction_gt500mm": 0.01,
                "precip_orographic_concentration": 1.05,
            },
        },
    }))
    summary = {
        "earth_reference_npz": "earth_reference_8000cells.npz",
        "earth_metrics": {
            "land_precip_mean_mm_yr": 730.0,
            "land_precip_p50_mm_yr": 500.0,
            "land_precip_p90_mm_yr": 1700.0,
            "global_mean_temperature_C": 15.7,
            "land_mean_temperature_C": 9.2,
            "ocean_mean_temperature_C": 18.3,
            "current_speed_p90_m_s": 0.27,
        },
        "entries": [{
            "label": "earthlike_seed1",
            "preset": "earthlike_mobile_lid",
            "seed": 1,
            "mode": "earthlike_calibration",
            "earth_distance_score": 0.9,
            "flags": ["earthlike_land_precip_too_dry"],
            "metrics": {
                "seed": 1,
                "assets_dir": str(assets),
                "land_fraction": 0.27,
                "global_mean_temperature_C": 14.2,
                "land_mean_temperature_C": 11.0,
                "ocean_mean_temperature_C": 15.0,
                "land_precip_mean_mm_yr": 150.0,
                "land_precip_p50_mm_yr": 110.0,
                "land_precip_p90_mm_yr": 350.0,
                "precip_seasonality_land_p75": 2.0,
                "current_speed_p90_m_s": 0.45,
                "biome_desert_area_fraction": 0.22,
                "biome_forest_area_fraction": 0.0,
                "biome_tropical_area_fraction": 0.0,
            },
            "comparison": {
                "global_mean_temperature_C": _comparison(14.2, 15.7),
                "land_mean_temperature_C": _comparison(11.0, 9.2),
                "ocean_mean_temperature_C": _comparison(15.0, 18.3),
            },
        }],
    }
    comparison_path = tmp_path / "comparison.json"
    comparison_path.write_text(json.dumps(summary))

    report = run_earth_climate_fitting_report(
        EarthClimateFittingConfig(comparison_path, tmp_path / "out")
    )

    assert report["schema"] == SCHEMA
    assert report["earthlike_run_count"] == 1
    assert report["phase_assessment"]["F4_seasonal_hydroclimate"]["status"] == "dominant_blocker"
    assert report["phase_assessment"]["F5_koppen_biomes"]["status"] == "blocked_by_hydroclimate"
    assert report["overall_verdict"] == "fail"
    assert report["guardrail_assessment"]["failure_count"] >= 1
    assert report["generated_rows"][0]["land_precip_mean_ratio_to_earth"] < 0.45
    assert (tmp_path / "out" / "earth_climate_fitting_report.md").exists()
    assert (tmp_path / "out" / "earth_climate_fitting_levers.csv").exists()
    assert (tmp_path / "out" / "earth_climate_guardrails.csv").exists()
    assert report["guardrails_csv"].endswith("earth_climate_guardrails.csv")


def test_earth_climate_fitting_report_handles_no_earthlike_runs(tmp_path):
    summary = {
        "earth_reference_npz": "earth_reference_8000cells.npz",
        "earth_metrics": {},
        "entries": [{
            "label": "waterworld_seed1",
            "preset": "waterworld",
            "seed": 1,
            "mode": "diagnostic_only",
            "flags": [],
            "metrics": {},
            "comparison": {},
        }],
    }
    comparison_path = tmp_path / "comparison.json"
    comparison_path.write_text(json.dumps(summary))

    report = run_earth_climate_fitting_report(
        EarthClimateFittingConfig(comparison_path, tmp_path / "out")
    )

    assert report["earthlike_run_count"] == 0
    assert report["phase_assessment"]["F1_temperature_energy"]["status"] == "no_earthlike_runs"
    assert len(report["generated_rows"]) == 1
    assert report["guardrail_assessment"]["verdict"] == "pass"
    assert report["guardrail_assessment"]["skipped_count"] == 4
    assert report["guardrail_assessment"]["failure_count"] == 0


def test_earth_climate_fitting_report_scores_f3_with_summer_low_mid_arrays(tmp_path):
    assets = tmp_path / "earthlike_seed2"
    assets.mkdir()
    n = 10
    lat = np.array([-60.0, -25.0, -15.0, -8.0, 0.0, 8.0, 15.0, 25.0, 45.0, 60.0])
    land = np.ones(n, dtype=bool)
    moisture = np.full((4, n), 0.30, dtype=float)
    monsoon = np.zeros((4, n), dtype=float)
    low_mid_nh = land & (lat >= 5.0) & (lat <= 35.0)
    low_mid_sh = land & (lat <= -5.0) & (lat >= -35.0)
    moisture[2, low_mid_nh] = 0.82
    moisture[0, low_mid_sh] = 0.80
    monsoon[2, low_mid_nh] = 0.34
    monsoon[0, low_mid_sh] = 0.32
    monsoon[0, low_mid_nh] = 0.04
    monsoon[2, low_mid_sh] = 0.03
    np.savez(
        assets / "terminal_climate_arrays.npz",
        lat=lat,
        terrain__elevation_m=np.where(land, 100.0, -1000.0),
        sea_level_m=np.asarray([0.0]),
        atmosphere__moisture_access=moisture,
        atmosphere__monsoon_potential=monsoon,
    )
    (assets / "summary.json").write_text(json.dumps({
        "climate_step_diagnostics": {"moisture_access_land_p75": 0.60},
        "climate_diagnostics": {
            "circulation": {
                "moisture_access_land_p75": 0.60,
                "monsoon_potential_land_p90": 0.08,
                "monsoon_potential_land_p99": 0.20,
            },
            "precipitation": {
                "land_wet_fraction_gt500mm": 0.2,
                "precip_orographic_concentration": 1.1,
            },
        },
    }))
    summary = {
        "earth_reference_npz": "earth_reference_8000cells.npz",
        "earth_metrics": {
            "land_precip_mean_mm_yr": 730.0,
            "land_precip_p50_mm_yr": 500.0,
            "land_precip_p90_mm_yr": 1700.0,
            "global_mean_temperature_C": 15.7,
            "land_mean_temperature_C": 9.2,
            "ocean_mean_temperature_C": 18.3,
            "current_speed_p90_m_s": 0.27,
        },
        "entries": [{
            "label": "earthlike_seed2",
            "preset": "earthlike_mobile_lid",
            "seed": 2,
            "mode": "earthlike_calibration",
            "earth_distance_score": 0.2,
            "flags": [],
            "metrics": {
                "seed": 2,
                "assets_dir": str(assets),
                "arrays": str(assets / "terminal_climate_arrays.npz"),
                "land_fraction": 0.30,
                "global_mean_temperature_C": 15.5,
                "land_mean_temperature_C": 9.5,
                "ocean_mean_temperature_C": 18.0,
                "land_precip_mean_mm_yr": 560.0,
                "land_precip_p50_mm_yr": 470.0,
                "land_precip_p90_mm_yr": 1300.0,
                "precip_seasonality_land_p75": 2.0,
                "current_speed_p90_m_s": 0.27,
                "biome_desert_area_fraction": 0.08,
                "biome_forest_area_fraction": 0.04,
                "biome_tropical_area_fraction": 0.02,
            },
            "comparison": {
                "global_mean_temperature_C": _comparison(15.5, 15.7),
                "land_mean_temperature_C": _comparison(9.5, 9.2),
                "ocean_mean_temperature_C": _comparison(18.0, 18.3),
            },
        }],
    }
    comparison_path = tmp_path / "comparison.json"
    comparison_path.write_text(json.dumps(summary))

    report = run_earth_climate_fitting_report(
        EarthClimateFittingConfig(comparison_path, tmp_path / "out")
    )

    f3 = report["phase_assessment"]["F3_circulation_moisture_access"]
    row = report["generated_rows"][0]
    assert f3["status"] == "watch"
    assert f3["score"] == 0.0
    assert row["monsoon_potential_land_p90"] == 0.08
    assert row["monsoon_potential_summer_p90_low_mid"] > 0.30
    assert row["monsoon_potential_summer_minus_winter_p75_low_mid"] > 0.25


def test_earth_climate_fitting_guardrails_warn_on_thin_wet_tail(tmp_path):
    earth = {
        "land_precip_mean_mm_yr": 730.0,
        "land_precip_p50_mm_yr": 500.0,
        "land_precip_p90_mm_yr": 1700.0,
        "global_mean_temperature_C": 15.7,
        "land_mean_temperature_C": 9.2,
        "ocean_mean_temperature_C": 18.3,
        "current_speed_p90_m_s": 0.27,
    }
    summary = {
        "earth_reference_npz": "earth_reference_8000cells.npz",
        "earth_metrics": earth,
        "entries": [
            {
                "label": "earthlike_seed1",
                "preset": "earthlike_mobile_lid",
                "seed": 1,
                "mode": "earthlike_calibration",
                "earth_distance_score": 0.4,
                "flags": [],
                "metrics": {
                    "seed": 1,
                    "land_fraction": 0.28,
                    "global_mean_temperature_C": 15.0,
                    "land_mean_temperature_C": 10.0,
                    "ocean_mean_temperature_C": 17.0,
                    "land_precip_mean_mm_yr": 520.0,
                    "land_precip_p50_mm_yr": 500.0,
                    "land_precip_p90_mm_yr": 800.0,
                    "precip_seasonality_land_p75": 2.0,
                    "current_speed_p90_m_s": 0.29,
                    "biome_desert_area_fraction": 0.08,
                    "biome_forest_area_fraction": 0.04,
                    "biome_tropical_area_fraction": 0.02,
                },
                "comparison": {
                    "global_mean_temperature_C": _comparison(15.0, 15.7),
                    "land_mean_temperature_C": _comparison(10.0, 9.2),
                    "ocean_mean_temperature_C": _comparison(17.0, 18.3),
                },
            },
            {
                "label": "arid_seed1",
                "preset": "arid_world",
                "seed": 2,
                "mode": "diagnostic_only",
                "flags": [],
                "metrics": {
                    "seed": 2,
                    "land_fraction": 0.8,
                    "land_precip_mean_mm_yr": 120.0,
                    "land_precip_p50_mm_yr": 40.0,
                    "land_precip_p90_mm_yr": 350.0,
                    "current_speed_p90_m_s": 0.28,
                    "biome_desert_area_fraction": 0.7,
                    "biome_forest_area_fraction": 0.003,
                    "biome_tropical_area_fraction": 0.002,
                },
                "comparison": {},
            },
            {
                "label": "waterworld_seed1",
                "preset": "waterworld",
                "seed": 3,
                "mode": "diagnostic_only",
                "flags": [],
                "metrics": {
                    "seed": 3,
                    "land_fraction": 0.04,
                    "land_precip_mean_mm_yr": 600.0,
                    "land_precip_p50_mm_yr": 560.0,
                    "land_precip_p90_mm_yr": 900.0,
                    "current_speed_p90_m_s": 0.28,
                    "biome_desert_area_fraction": 0.0,
                    "biome_forest_area_fraction": 0.005,
                    "biome_tropical_area_fraction": 0.0,
                },
                "comparison": {},
            },
        ],
    }
    comparison_path = tmp_path / "comparison.json"
    comparison_path.write_text(json.dumps(summary))

    report = run_earth_climate_fitting_report(
        EarthClimateFittingConfig(comparison_path, tmp_path / "out")
    )

    assert report["overall_verdict"] == "pass_with_warnings"
    assert report["guardrail_assessment"]["failure_count"] == 0
    assert report["guardrail_assessment"]["warning_count"] >= 1
    assert any(
        item["metric"] == "land_precip_p90_ratio_to_earth"
        for item in report["guardrail_assessment"]["warnings"]
    )
    assert (tmp_path / "out" / "earth_climate_guardrails.csv").exists()


def test_earth_climate_fit_report_cli_can_fail_on_guardrail(tmp_path):
    summary = {
        "earth_reference_npz": "earth_reference_8000cells.npz",
        "earth_metrics": {
            "land_precip_mean_mm_yr": 730.0,
            "land_precip_p50_mm_yr": 500.0,
            "land_precip_p90_mm_yr": 1700.0,
            "global_mean_temperature_C": 15.7,
            "land_mean_temperature_C": 9.2,
            "ocean_mean_temperature_C": 18.3,
            "current_speed_p90_m_s": 0.27,
        },
        "entries": [{
            "label": "earthlike_seed1",
            "preset": "earthlike_mobile_lid",
            "seed": 1,
            "mode": "earthlike_calibration",
            "earth_distance_score": 0.9,
            "flags": [],
            "metrics": {
                "seed": 1,
                "land_fraction": 0.27,
                "global_mean_temperature_C": 14.0,
                "land_mean_temperature_C": 10.0,
                "ocean_mean_temperature_C": 16.0,
                "land_precip_mean_mm_yr": 200.0,
                "land_precip_p50_mm_yr": 120.0,
                "land_precip_p90_mm_yr": 350.0,
                "current_speed_p90_m_s": 0.3,
                "biome_desert_area_fraction": 0.25,
                "biome_forest_area_fraction": 0.0,
                "biome_tropical_area_fraction": 0.0,
            },
            "comparison": {
                "global_mean_temperature_C": _comparison(14.0, 15.7),
                "land_mean_temperature_C": _comparison(10.0, 9.2),
                "ocean_mean_temperature_C": _comparison(16.0, 18.3),
            },
        }],
    }
    comparison_path = tmp_path / "comparison.json"
    comparison_path.write_text(json.dumps(summary))

    args = SimpleNamespace(
        comparison_summary=comparison_path,
        out=tmp_path / "out",
        fail_on_guardrail=True,
    )
    with pytest.raises(SystemExit) as exc_info:
        cmd_earth_climate_fit_report(args)

    assert exc_info.value.code == 2
    assert (tmp_path / "out" / "earth_climate_guardrails.csv").exists()
