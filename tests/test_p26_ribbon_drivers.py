import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.state import WorldState
from aevum.diagnostics.p26_ribbon_drivers import ribbon_driver_summary
from aevum.modules.tectonics import (
    CONT,
    DOMAIN_ACCRETED_TERRANE,
    DOMAIN_CONTINENTAL_INTERIOR,
    DOMAIN_OCEANIC,
    ORIGIN_ARC,
    ORIGIN_PRIMORDIAL,
    ORIGIN_RIDGE,
    ORIGIN_SUTURE,
)
from aevum.modules.terrain import CONT_DETAIL_ARC_MICROCONTINENT, CONT_DETAIL_PLATFORM
from aevum.spec.presets import get_preset


def _unit(lat_deg: float, lon_deg: float) -> np.ndarray:
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return np.array([
        np.cos(lat) * np.cos(lon),
        np.cos(lat) * np.sin(lon),
        np.sin(lat),
    ])


def test_p26_ribbon_drivers_identify_young_accretionary_ribbon():
    grid = SphereGrid.fibonacci(1800, 6.371e6)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    world.set_g("ocean.sea_level_m", 0.0)

    broad_center = _unit(0.0, -80.0)
    broad = np.degrees(np.arccos(np.clip(grid.xyz @ broad_center, -1.0, 1.0))) < 26.0
    ribbon = (np.abs(grid.lat) < 4.0) & (grid.lon > 40.0) & (grid.lon < 145.0)
    land = broad | ribbon

    world.set_field("terrain.elevation_m", np.where(land, 600.0, -4200.0))
    crust_type = np.where(land, CONT, 0.0)
    world.set_field("crust.type", crust_type)
    world.set_field("crust.age_myr", np.where(broad, 1800.0, np.where(ribbon, 180.0, 70.0)))
    world.set_field(
        "crust.domain",
        np.where(
            ribbon,
            DOMAIN_ACCRETED_TERRANE,
            np.where(broad, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC),
        ),
    )
    world.set_field(
        "crust.origin",
        np.where(ribbon, ORIGIN_ARC, np.where(broad, ORIGIN_PRIMORDIAL, ORIGIN_RIDGE)),
    )
    world.set_field("crust.stability", np.where(broad, 0.80, np.where(ribbon, 0.18, 0.0)))
    world.set_field("crust.reworked_age_myr", np.where(ribbon, 3700.0, -1.0))
    world.set_field("tectonics.orogeny_age_myr", np.full(grid.n, -1.0))
    world.set_field("tectonics.volcanism_age_myr", np.where(ribbon, 3820.0, -1.0))
    world.set_field(
        "terrain.continental_detail",
        np.where(
            ribbon,
            CONT_DETAIL_ARC_MICROCONTINENT,
            np.where(broad, CONT_DETAIL_PLATFORM, 0.0),
        ),
    )

    summary = ribbon_driver_summary(world, top_n=3)

    assert summary["schema"] == "aevum.p26_ribbon_drivers.v1"
    assert summary["summary"]["exposed_land_component_count"] >= 2
    assert (
        summary["summary"][
            "exposed_land_quiet_inherited_arc_suture_ribbon_share_gt650_myr"
        ]
        > 0.50
    )
    assert summary["summary"]["exposed_land_active_rework_ribbon_share_lt220_myr"] == 0.0
    assert summary["summary"]["primary_temporal_driver_hint"] == (
        "quiet inherited arc/suture provenance dominates ribbon area"
    )
    top = summary["top_exposed_land_components"][0]
    assert top["domain_shares"]["accreted_terrane"] > 0.90
    assert top["origin_shares"]["arc"] > 0.90
    assert top["continental_detail_shares"]["arc_microcontinent"] > 0.90
    assert top["mean_stability"] < 0.30
    assert top["ribbon_area_fraction_of_component"] > 0.25
    assert top["quiet_inherited_arc_suture_share_gt650_myr"] > 0.90
    assert top["active_rework_share_lt220_myr"] == 0.0
    assert top["mean_time_since_rework_myr"] > 600.0
    assert "accretionary" in top["driver_hint"]


def test_p26_ribbon_drivers_identify_recent_active_rework_ribbon():
    grid = SphereGrid.fibonacci(1800, 6.371e6)
    spec = get_preset("earthlike")
    spec.grid_cells = grid.n
    world = WorldState(grid=grid, spec=spec, time_myr=4500.0)
    world.set_g("ocean.sea_level_m", 0.0)

    ribbon = (np.abs(grid.lat) < 4.0) & (grid.lon > -125.0) & (grid.lon < -20.0)
    stable = np.degrees(np.arccos(np.clip(grid.xyz @ _unit(10.0, 80.0), -1.0, 1.0))) < 24.0
    land = ribbon | stable

    world.set_field("terrain.elevation_m", np.where(land, 700.0, -4300.0))
    world.set_field("crust.type", np.where(land, CONT, 0.0))
    world.set_field("crust.age_myr", np.where(stable, 2100.0, np.where(ribbon, 700.0, 80.0)))
    world.set_field(
        "crust.domain",
        np.where(
            ribbon,
            DOMAIN_ACCRETED_TERRANE,
            np.where(stable, DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_OCEANIC),
        ),
    )
    world.set_field(
        "crust.origin",
        np.where(ribbon, ORIGIN_SUTURE, np.where(stable, ORIGIN_PRIMORDIAL, ORIGIN_RIDGE)),
    )
    world.set_field("crust.stability", np.where(stable, 0.78, np.where(ribbon, 0.24, 0.0)))
    world.set_field("crust.reworked_age_myr", np.where(ribbon, 4460.0, -1.0))
    world.set_field("tectonics.orogeny_age_myr", np.where(ribbon, 4440.0, -1.0))
    world.set_field("tectonics.volcanism_age_myr", np.full(grid.n, -1.0))
    world.set_field(
        "terrain.continental_detail",
        np.where(ribbon, CONT_DETAIL_ARC_MICROCONTINENT, CONT_DETAIL_PLATFORM),
    )

    summary = ribbon_driver_summary(world, top_n=3)

    assert summary["summary"]["exposed_land_active_rework_ribbon_share_lt220_myr"] > 0.50
    assert (
        summary["summary"][
            "exposed_land_quiet_inherited_arc_suture_ribbon_share_gt650_myr"
        ]
        == 0.0
    )
    assert summary["summary"]["primary_temporal_driver_hint"] == (
        "recent collision/suture rework dominates ribbon area"
    )
