from types import SimpleNamespace

import numpy as np

from aevum.archive.world_archive import Frame
from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.diagnostics.p173_unsupported_shoal_attribution import (
    frame_unsupported_shoal_attribution,
    unsupported_shoal_attribution_summary,
    write_p173_unsupported_shoal_attribution,
)
from aevum.modules.terrain import (
    OCEAN_DEPTH_ABYSS,
    OCEAN_DEPTH_RIDGE,
    OCEAN_MARGIN_OPEN,
)
from aevum.spec.presets import get_preset


def _shoal_fixture():
    grid = SphereGrid.fibonacci(360, CONSTANTS.EARTH_RADIUS)
    n = grid.n
    land = grid.lat > 74.0
    ocean = ~land
    ridge = np.where(
        ocean & (np.abs(grid.lat) < 16.0) & (grid.lon < -100.0) & (grid.lon > -155.0)
    )[0][:10]
    backed = np.where(
        ocean & (np.abs(grid.lat) < 16.0) & (grid.lon < -15.0) & (grid.lon > -75.0)
    )[0][:10]
    unsupported = np.where(
        ocean & (np.abs(grid.lat) < 16.0) & (grid.lon > 45.0) & (grid.lon < 105.0)
    )[0][:10]
    assert ridge.size >= 6
    assert backed.size >= 6
    assert unsupported.size >= 6

    surface = np.where(land, 620.0, -4200.0).astype(np.float64)
    surface[ridge] = -1100.0
    surface[backed] = -900.0
    surface[unsupported] = -450.0
    age = np.where(land, 1600.0, 90.0).astype(np.float64)
    age[ridge] = 5.0
    depth_province = np.where(ocean, OCEAN_DEPTH_ABYSS, 0.0).astype(np.float64)
    depth_province[ridge] = float(OCEAN_DEPTH_RIDGE)
    fields = {
        "terrain.elevation_m": surface,
        "crust.age_myr": age,
        "crust.type": land.astype(np.float64),
        "crust.domain": np.zeros(n, dtype=np.float64),
        "crust.origin": np.zeros(n, dtype=np.float64),
        "ocean.depth_province": depth_province,
        "ocean.margin_type": np.where(ocean, OCEAN_MARGIN_OPEN, 0.0).astype(np.float64),
        "ocean.shelf_width": np.where(ocean, 5.0, 0.0).astype(np.float64),
        "terrain.rift_margin_stage": np.zeros(n, dtype=np.float64),
    }
    frame = Frame(
        time_myr=1900.0,
        globals={
            "ocean.sea_level_m": 0.0,
            "terrain.last_p1732_young_open_ocean_age_depth_used": 1.0,
            "terrain.last_p1732_young_open_ocean_age_depth_land_mask_preserved": 1.0,
            "terrain.last_p1732_young_open_ocean_age_depth_candidate_fraction": 0.07,
            "terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction": 0.06,
            "terrain.last_p1732_young_open_ocean_age_depth_mean_depth_before_m": 300.0,
            "terrain.last_p1732_young_open_ocean_age_depth_mean_depth_after_m": 3100.0,
        },
        fields=fields,
        objects={
            "terrain.arc_plume_landforms": [
                {"kind": "seamount_chain", "cells": backed.astype(int).tolist()},
            ],
        },
    )
    return grid, frame


def test_p1731_frame_attribution_splits_preserve_support_and_residual():
    grid, frame = _shoal_fixture()

    row = frame_unsupported_shoal_attribution(grid, frame)

    assert row["usable"]
    assert row["cleanup_candidate_fraction_of_ocean"] > 0.0
    assert row["structural_preserve_fraction_of_candidate"] > 0.0
    assert row["object_support_fraction_of_candidate"] > 0.0
    assert row["post_cleanup_residual_fraction_of_ocean"] > 0.0
    assert row["mask_counts"]["post_cleanup_residual"] > 0
    assert row["mask_fingerprints"]["post_cleanup_residual"] != "empty"
    assert row["residual_attribution"]["owner_hint"] in {
        "open_ocean_young_shallow_abyss_rise",
        "open_ocean_shallow_abyss_rise",
    }
    assert row["residual_attribution"]["by_depth_province"][0]["name"] == "abyss"
    assert row["residual_attribution"]["component_summary"]["component_count"] >= 1
    assert row["p1732_young_open_ocean_depth_floor"][
        "terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction"
    ] == 0.06


def test_p1731_summary_and_writer_report_peak_frame(tmp_path):
    grid, frame = _shoal_fixture()
    spec = get_preset("earthlike")
    world = SimpleNamespace(grid=grid, spec=spec)
    archive = SimpleNamespace(world=world, frames=[frame])

    summary = unsupported_shoal_attribution_summary(world, archive)
    written = write_p173_unsupported_shoal_attribution(world, archive, tmp_path)

    assert summary["schema"] == "aevum.p173_unsupported_shoal_attribution.v1"
    assert summary["usable_frame_count"] == 1
    assert summary["acceptance"]["attribution_completed"]
    assert summary["peak_residual_frame"]["time_myr"] == 1900.0
    assert summary["metric_extremes"][
        "post_cleanup_residual_fraction_of_ocean"]["max"] > 0.0
    assert summary["p1732_metric_extremes"][
        "terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction"
    ]["max"] == 0.06
    assert summary["peak_residual_frame"]["p1732_young_open_ocean_depth_floor"][
        "terrain.last_p1732_young_open_ocean_age_depth_candidate_fraction"
    ] == 0.07
    assert written["peak_residual_frame"]["owner_hint"] == summary[
        "peak_residual_frame"]["owner_hint"]
    assert (tmp_path / "p173_unsupported_shoal_attribution.json").exists()
