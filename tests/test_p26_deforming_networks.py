import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.state import WorldState
from aevum.diagnostics.p26_deforming_networks import deforming_network_summary
from aevum.modules.tectonics import (
    CONT,
    DEFORM_COLLISION_CORE,
    DEFORM_COLLISION_SHOULDER,
    DEFORM_RIFT,
    OCEAN,
)
from aevum.spec.presets import get_preset


def _unit(lat_deg: float, lon_deg: float) -> np.ndarray:
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return np.array([
        np.cos(lat) * np.cos(lon),
        np.cos(lat) * np.sin(lon),
        np.sin(lat),
    ])


def test_p26_deforming_network_summary_tracks_core_shoulder_and_ribbon_overlap():
    grid = SphereGrid.fibonacci(1800, 6.371e6)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    world.set_g("ocean.sea_level_m", 0.0)

    broad_land = (
        np.degrees(np.arccos(np.clip(grid.xyz @ _unit(5.0, -80.0), -1.0, 1.0)))
        < 26.0
    )
    ribbon = (np.abs(grid.lat) < 4.0) & (grid.lon > 35.0) & (grid.lon < 150.0)
    land = broad_land | ribbon
    core = ribbon & (grid.lon > 75.0) & (grid.lon < 110.0)
    shoulder = ribbon & ~core
    rift = broad_land & (grid.lon < -83.0) & (grid.lon > -98.0) & (grid.lat > -7.0)

    world.set_field("terrain.elevation_m", np.where(land, 600.0, -4200.0))
    world.set_field("crust.type", np.where(land, CONT, OCEAN))
    intensity = np.zeros(grid.n, dtype=np.float64)
    intensity[core] = 1.0
    intensity[shoulder] = 0.45
    intensity[rift] = 0.65
    style = np.zeros(grid.n, dtype=np.float64)
    style[core] = DEFORM_COLLISION_CORE
    style[shoulder] = DEFORM_COLLISION_SHOULDER
    style[rift] = DEFORM_RIFT
    world.set_field("tectonics.deformation_intensity", intensity)
    world.set_field("tectonics.deformation_style", style)
    world.objects["tectonics.deforming_networks"] = [
        {
            "id": "deforming_network:collision_core:0",
            "kind": "collision_core",
            "style_code": int(DEFORM_COLLISION_CORE),
            "cell_count": int(core.sum()),
            "area_fraction": float(grid.cell_area[core].sum() / grid.cell_area.sum()),
            "mean_intensity": 1.0,
            "continental_fraction": 1.0,
            "lat": 0.0,
            "lon": 92.0,
            "last_active_myr": 4500.0,
        },
        {
            "id": "deforming_network:collision_shoulder:0",
            "kind": "collision_shoulder",
            "style_code": int(DEFORM_COLLISION_SHOULDER),
            "cell_count": int(shoulder.sum()),
            "area_fraction": float(grid.cell_area[shoulder].sum() / grid.cell_area.sum()),
            "mean_intensity": 0.45,
            "continental_fraction": 1.0,
            "lat": 0.0,
            "lon": 65.0,
            "last_active_myr": 4500.0,
        },
    ]

    summary = deforming_network_summary(world, top_n=2)

    assert summary["schema"] == "aevum.p26_deforming_networks.v1"
    metrics = summary["metrics"]
    assert metrics["active_deformation_area_fraction_of_world"] > 0.0
    assert metrics["active_deformation_area_fraction_of_continental"] > 0.0
    assert metrics["core_area_fraction_of_active"] > 0.0
    assert metrics["shoulder_area_fraction_of_active"] > 0.0
    assert metrics["land_ribbon_deformation_coverage_fraction"] > 0.50
    assert metrics["deforming_network_object_count"] == 2
    assert summary["style_summaries"]["collision_core"]["mean_intensity"] == 1.0
    assert 0.40 <= summary["style_summaries"]["collision_shoulder"]["mean_intensity"] <= 0.50
    assert summary["style_summaries"]["rift"]["cell_count"] > 0
    assert summary["object_kind_counts"]["collision_core"] == 1
    assert summary["object_kind_counts"]["collision_shoulder"] == 1
    assert summary["diagnostic_hints"]["land_ribbon_is_deformation_coupled"]
