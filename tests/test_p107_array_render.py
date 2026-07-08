import json

import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.diagnostics.p107_array_render import render_p107_array_assets


def test_p107_array_renderer_writes_visual_qa_pack(tmp_path):
    grid = SphereGrid.fibonacci(180, CONSTANTS.EARTH_RADIUS)
    land = grid.lat > 0.0
    elevation = np.where(land, 400.0 + grid.lat * 10.0, -3500.0)
    elevation[(grid.lat < -20.0) & (grid.lon > 40.0)] = -6200.0
    plate = np.arange(grid.n) % 8
    crust_type = land.astype(np.int32)
    crust_age = np.where(land, 1600.0, 80.0)
    ocean_basin = np.where(land, -1, (grid.lon > 0.0).astype(np.int32))
    depth_province = np.where(land, 0, 4).astype(np.int32)
    depth_province[(grid.lat < -20.0) & (grid.lon > 40.0)] = 6
    terrain_province = np.where(land, 2, 0).astype(np.int32)
    detail = np.where(land, 3, 0).astype(np.int32)
    detail_region = np.where(land, 2, 0).astype(np.int32)
    internal_block_region = np.where(land, 2, 0).astype(np.int32)
    inland_region = np.where(land, 3, 0).astype(np.int32)
    terrain_cont_province = np.where(land, 5, 0).astype(np.int32)
    hierarchy = np.zeros(grid.n, dtype=np.int32)
    hierarchy[(grid.lat > 30.0) & land] = 1
    hierarchy[(grid.lat > 45.0) & land] = 2
    hierarchy[(grid.lat > 60.0) & land] = 3
    spine = np.zeros(grid.n, dtype=np.int32)
    spine[(grid.lat > 46.0) & (grid.lat < 50.0) & land] = 2
    spine[(grid.lat > 64.0) & (grid.lat < 68.0) & land] = 3
    halo = np.zeros(grid.n, dtype=np.int32)
    halo[(grid.lat > 25.0) & (grid.lat < 35.0) & land] = 1
    apron = np.zeros(grid.n, dtype=np.int32)
    apron[(grid.lat > 52.0) & (grid.lat < 58.0) & land] = 1
    ridge = ((np.abs(grid.lat) < 8.0) & ~land).astype(np.uint8)
    trench = ((grid.lat < -20.0) & (grid.lon > 40.0)).astype(np.uint8)

    arrays_path = tmp_path / "p107_terminal_arrays.npz"
    np.savez_compressed(
        arrays_path,
        grid_lat=grid.lat.astype(np.float32),
        grid_lon=grid.lon.astype(np.float32),
        grid_cell_area_m2=grid.cell_area.astype(np.float64),
        field__terrain_elevation_m=elevation.astype(np.float32),
        field__tectonics_plate_id=plate.astype(np.int32),
        field__crust_type=crust_type.astype(np.int32),
        field__crust_age_myr=crust_age.astype(np.float32),
        field__ocean_basin_id=ocean_basin.astype(np.int32),
        field__ocean_depth_province=depth_province.astype(np.int32),
        field__terrain_province=terrain_province.astype(np.int32),
        field__terrain_continental_detail=detail.astype(np.int32),
        field__terrain_continental_detail_region_code=(
            detail_region.astype(np.int32)
        ),
        field__terrain_internal_geographic_block_region_code=(
            internal_block_region.astype(np.int32)
        ),
        field__terrain_inland_geomorphology_region_code=(
            inland_region.astype(np.int32)
        ),
        field__terrain_continental_province_code=(
            terrain_cont_province.astype(np.int32)
        ),
        field__terrain_orogenic_parent_hierarchy=hierarchy.astype(np.int32),
        field__terrain_orogenic_hierarchy_spine=spine.astype(np.int32),
        field__terrain_orogenic_shoulder_halo=halo.astype(np.int32),
        field__terrain_orogenic_highland_apron=apron.astype(np.int32),
        boundary__ridge=ridge,
        boundary__trench=trench,
        object__province_mid_ocean_ridge=ridge,
        object__margin_trench=trench,
        object__island_arc=np.zeros(grid.n, dtype=np.uint8),
        object__parent_orogen_highland_apron=apron.astype(np.uint8),
    )
    metrics = {
        "schema": "aevum.p107_terminal_world_audit.v1",
        "cells": grid.n,
        "sea_level_m": 0.0,
        "array_archive": {
            "path": str(arrays_path),
            "manifest": {
                "fields": {
                    "terrain.elevation_m": "field__terrain_elevation_m",
                    "tectonics.plate_id": "field__tectonics_plate_id",
                    "crust.type": "field__crust_type",
                    "crust.age_myr": "field__crust_age_myr",
                    "ocean.basin_id": "field__ocean_basin_id",
                    "ocean.depth_province": "field__ocean_depth_province",
                    "terrain.province": "field__terrain_province",
                    "terrain.continental_detail": "field__terrain_continental_detail",
                    "terrain.continental_detail_region_code": (
                        "field__terrain_continental_detail_region_code"
                    ),
                    "terrain.internal_geographic_block_region_code": (
                        "field__terrain_internal_geographic_block_region_code"
                    ),
                    "terrain.inland_geomorphology_region_code": (
                        "field__terrain_inland_geomorphology_region_code"
                    ),
                    "terrain.continental_province_code": (
                        "field__terrain_continental_province_code"
                    ),
                    "terrain.orogenic_parent_hierarchy": (
                        "field__terrain_orogenic_parent_hierarchy"
                    ),
                    "terrain.orogenic_hierarchy_spine": (
                        "field__terrain_orogenic_hierarchy_spine"
                    ),
                    "terrain.orogenic_shoulder_halo": (
                        "field__terrain_orogenic_shoulder_halo"
                    ),
                    "terrain.orogenic_highland_apron": (
                        "field__terrain_orogenic_highland_apron"
                    ),
                }
            },
        },
        "p109_hypsometry_comparison": {
            "within_p109_envelope": True,
        },
        "p110a_modern_planform": {
            "summary": {
                "largest_land_component_share": 0.50,
                "second_land_component_share": 0.25,
                "third_land_component_share": 0.12,
            },
            "warning_flags": [],
        },
    }
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(metrics) + "\n")

    summary = render_p107_array_assets(
        metrics_path,
        tmp_path / "rendered",
        width=96,
        height=48,
    )

    assert summary["schema"] == "aevum.p107_array_render.v1"
    for name in [
        "elevation.png",
        "bathymetry_shelf_slope_abyss.png",
        "plates.png",
        "crust_age.png",
        "tectonic_boundaries.png",
        "object_masks.png",
        "continental_detail_raw_provinces.png",
        "continental_detail_region_provinces.png",
        "continental_detail_provinces.png",
        "internal_geographic_block_region_code.png",
        "inland_geomorphology_regions.png",
        "terrain_continental_province_code.png",
        "orogenic_parent_hierarchy.png",
        "orogenic_parent_hierarchy_overlay.png",
        "orogenic_hierarchy_spines.png",
        "orogenic_hierarchy_spine_overlay.png",
        "orogenic_shoulder_halo.png",
        "orogenic_highland_apron.png",
        "orogenic_belt_morphology_overlay.png",
        "p107_array_contact_sheet.png",
    ]:
        path = tmp_path / "rendered" / name
        assert path.exists()
        assert path.stat().st_size > 0
    saved = json.loads((tmp_path / "rendered" / "p107_array_render_summary.json").read_text())
    assert saved["cells"] == grid.n
