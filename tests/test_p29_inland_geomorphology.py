import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.state import WorldState
from aevum.diagnostics.p29_inland_geomorphology import inland_geomorphology_summary
from aevum.modules.tectonics import CONT, OCEAN, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_CRATON
from aevum.modules.terrain import (
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATFORM,
    CONT_DETAIL_PLATEAU,
    CONT_DETAIL_RIFT_BASIN,
    CONT_DETAIL_SHIELD,
)
from aevum.spec.presets import get_preset


def _world(n=900):
    spec = get_preset("earthlike")
    spec.grid_cells = n
    grid = SphereGrid.fibonacci(n, spec.radius_m)
    return WorldState(grid=grid, spec=spec, time_myr=4500.0)


def test_p29_inland_summary_distinguishes_mechanism_rich_interior():
    world = _world()
    grid = world.grid
    continent = (np.abs(grid.lat) < 44.0) & (np.abs(grid.lon) < 130.0)
    shield = continent & (grid.lon < -70.0) & (np.abs(grid.lat) < 25.0)
    basin = continent & (grid.lon > -35.0) & (grid.lon < -8.0) & (np.abs(grid.lat) < 25.0)
    orogen = continent & (grid.lon > 22.0) & (grid.lon < 42.0) & (np.abs(grid.lat) < 30.0)
    rift = continent & (grid.lon > 72.0) & (grid.lon < 88.0) & (np.abs(grid.lat) < 30.0)
    plateau = continent & (grid.lon > 48.0) & (grid.lon < 66.0) & (np.abs(grid.lat) < 24.0)

    elev = np.where(continent, 760.0, -4200.0)
    elev[shield] = 480.0
    elev[basin] = 220.0
    elev[orogen] = 1450.0
    elev[rift] = 360.0
    elev[plateau] = 2500.0
    detail = np.where(continent, CONT_DETAIL_PLATFORM, 0.0)
    detail[shield] = CONT_DETAIL_SHIELD
    detail[basin] = CONT_DETAIL_BASIN
    detail[orogen] = CONT_DETAIL_OROGEN
    detail[rift] = CONT_DETAIL_RIFT_BASIN
    detail[plateau] = CONT_DETAIL_PLATEAU
    world.set_field("terrain.elevation_m", elev)
    world.set_field("terrain.continental_detail", detail)
    world.set_field("crust.type", np.where(continent, CONT, OCEAN))
    world.set_field(
        "crust.domain",
        np.where(shield, DOMAIN_CRATON, np.where(continent, DOMAIN_CONTINENTAL_INTERIOR, 0.0)),
    )
    world.objects["terrain.continental_landforms"] = [
        {"id": "shield:test", "kind": "shield", "area_fraction": 0.02, "cell_count": 12},
        {"id": "platform:test", "kind": "platform", "area_fraction": 0.04, "cell_count": 30},
        {"id": "basin:test", "kind": "interior_basin", "area_fraction": 0.03, "cell_count": 20},
        {"id": "orogen:test", "kind": "old_subdued_orogen", "area_fraction": 0.02, "cell_count": 13},
        {"id": "rift:test", "kind": "rift_basin", "area_fraction": 0.02, "cell_count": 15},
        {"id": "plateau:test", "kind": "plateau", "area_fraction": 0.02, "cell_count": 12},
    ]

    summary = inland_geomorphology_summary(world)
    assert summary["schema"] == "aevum.p29_inland_geomorphology.v1"
    assert summary["metrics"]["inland_detail_diversity_gt2pct"] >= 5
    assert summary["metrics"]["inland_object_kind_count"] >= 6
    assert summary["metrics"]["inland_relief_p95_p05_m"] > 1200.0
    assert summary["diagnostic_hints"]["broad_lowland_mode_present"]
    assert summary["diagnostic_hints"]["highland_tail_present"]
    assert not summary["diagnostic_hints"]["monotone_flat_inland"]


def test_p29_inland_summary_flags_flat_single_mode_interior():
    world = _world()
    grid = world.grid
    continent = (np.abs(grid.lat) < 44.0) & (np.abs(grid.lon) < 130.0)
    world.set_field("terrain.elevation_m", np.where(continent, 720.0, -4200.0))
    world.set_field("terrain.continental_detail", np.where(continent, CONT_DETAIL_PLATFORM, 0.0))
    world.set_field("crust.type", np.where(continent, CONT, OCEAN))
    world.set_field("crust.domain", np.where(continent, DOMAIN_CONTINENTAL_INTERIOR, 0.0))
    world.objects["terrain.continental_landforms"] = [
        {"id": "platform:test", "kind": "platform", "area_fraction": 0.20, "cell_count": 200},
    ]

    summary = inland_geomorphology_summary(world)
    assert summary["diagnostic_hints"]["monotone_flat_inland"]
    assert summary["diagnostic_hints"]["low_inland_object_diversity"]
    assert summary["metrics"]["inland_relief_p95_p05_m"] == 0.0
    assert summary["metrics"]["inland_object_kind_count"] == 1
