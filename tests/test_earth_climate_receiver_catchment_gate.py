import json
from types import SimpleNamespace

import numpy as np
import pytest

from aevum.cli import cmd_earth_climate_receiver_catchment_gate
from aevum.diagnostics.earth_climate_receiver_catchment_gate import (
    EarthClimateReceiverCatchmentGateConfig,
    run_earth_climate_receiver_catchment_gate,
)


def _write_receiver_case(tmp_path, *, include_archive=True):
    assets = tmp_path / "earthlike_seed1"
    assets.mkdir()
    n = 12
    area = np.ones(n, dtype=float)
    land = np.array([True] * 8 + [False] * 4)
    receiver = np.full((4, n), -1.0, dtype=float)
    receiver[:, :4] = 1.0
    receiver[:, 4:8] = 2.0
    source = np.full((4, n), -1.0, dtype=float)
    source[:, :4] = 7.0
    source[:, 4:8] = 8.0
    budget = np.full((4, n), -1.0, dtype=float)
    budget[:, :8] = 11.0
    precip_region = np.full((4, n), -1.0, dtype=float)
    precip_region[:, :4] = 101.0
    precip_region[:, 4:8] = 201.0
    response = np.ones((4, n), dtype=float)
    response[:, :4] = 1.08
    source_supply = np.zeros((4, n), dtype=float)
    source_supply[:, :4] = 0.75
    source_supply[:, 4:8] = 0.62
    source_supply[:, 8:] = 0.85
    receiver_balance = np.zeros((4, n), dtype=float)
    receiver_balance[:, :8] = 0.72
    receiver_feedback = np.ones((4, n), dtype=float)
    receiver_feedback[:, :4] = 1.03
    receiver_feedback[:, 4:8] = 0.97
    np.savez(
        assets / "terminal_climate_arrays.npz",
        cell_area=area,
        terrain__elevation_m=np.where(land, 100.0, -1000.0),
        sea_level_m=np.asarray([0.0]),
        climate__receiver_catchment_id=receiver,
        atmosphere__moisture_source_basin_id=source,
        climate__moisture_budget_region_id=budget,
        climate__precipitation_response_region_id=precip_region,
        climate__moisture_flow_precipitation_response=response,
        climate__source_basin_supply_index=source_supply,
        climate__receiver_catchment_supply_balance=receiver_balance,
        climate__receiver_supply_precipitation_feedback=receiver_feedback,
    )
    receiver_json = assets / "receiver_catchments.json"
    if include_archive:
        catchments = []
        for season_index, season in enumerate(("DJF", "MAM", "JJA", "SON")):
            for idx, source_id in [(1, 7), (2, 8)]:
                catchments.append({
                    "id": f"receiver_catchment:{season.lower()}:{season_index}-{idx}",
                    "type": "receiver_catchment",
                    "kind": "source_receiver_catchment",
                    "season": season,
                    "season_index": season_index,
                    "catchment_index": idx,
                    "area_fraction": 0.25,
                    "source_basin_attributed_fraction": 1.0,
                    "source_basin_purity": 1.0,
                    "dominant_source_basin_id": source_id,
                    "budget_region_attributed_fraction": 1.0,
                    "precipitation_response_attributed_fraction": 1.0,
                    "mean_source_basin_supply_index": 0.7,
                    "source_basin_supply_attributed_fraction": 1.0,
                    "source_basin_supply_mass_fraction": 0.5,
                    "supply_supported_precipitation_fraction": 0.65,
                    "precipitation_supply_balance": 0.72,
                })
        receiver_json.write_text(json.dumps({
            "schema": "aevum.receiver_catchments.v1",
            "catchments": catchments,
        }))
    summary = {
        "summaries": [{
            "label": "earthlike_seed1",
            "preset": "earthlike_mobile_lid",
            "seed": 1,
            "assets_dir": str(assets),
            "arrays": str(assets / "terminal_climate_arrays.npz"),
            "receiver_catchments_json": str(receiver_json) if include_archive else None,
        }],
    }
    path = tmp_path / "terminal_summary.json"
    path.write_text(json.dumps(summary))
    return path


def test_receiver_catchment_gate_passes_complete_archive(tmp_path):
    summary = _write_receiver_case(tmp_path)

    report = run_earth_climate_receiver_catchment_gate(
        EarthClimateReceiverCatchmentGateConfig(summary, tmp_path / "out")
    )

    assert report["verdict"] == "pass"
    assert report["failure_count"] == 0
    assert report["schema"] == "aevum.earth_climate_receiver_catchment_gate.v3"
    assert (tmp_path / "out" / "earth_climate_receiver_catchment_metrics.csv").exists()
    assert (tmp_path / "out" / "earth_climate_receiver_catchment_checks.csv").exists()


def test_receiver_catchment_gate_flags_missing_archive(tmp_path):
    summary = _write_receiver_case(tmp_path, include_archive=False)

    report = run_earth_climate_receiver_catchment_gate(
        EarthClimateReceiverCatchmentGateConfig(summary, tmp_path / "out")
    )

    assert report["verdict"] == "fail"
    assert any(
        check["metric"] == "receiver_archive_found"
        for check in report["failures"]
    )


def test_receiver_catchment_gate_cli_can_fail(tmp_path):
    summary = _write_receiver_case(tmp_path, include_archive=False)
    args = SimpleNamespace(
        terminal_summary=str(summary),
        out=str(tmp_path / "out"),
        fail_on_receiver_catchment=True,
    )

    with pytest.raises(SystemExit):
        cmd_earth_climate_receiver_catchment_gate(args)
