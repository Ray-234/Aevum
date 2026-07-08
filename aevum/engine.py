"""Engine assembly: build a world, wire modules into the scheduler, run deep time.

Correct boot order (project plan, section 11):
    PlanetSpec -> FeatureRegistry -> WorldState -> Scheduler
    -> full-history placeholder flow -> MapCompiler -> per-module replacement.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from aevum.archive.world_archive import WorldArchive
from aevum.core.events import EventBus
from aevum.core.grid import SphereGrid
from aevum.core.registry import FeatureRegistry
from aevum.core.scheduler import DeepTimeScheduler
from aevum.core.state import WorldState
from aevum.features import build_registry
from aevum.modules import (BiogeochemModule, BiosphereModule, ClimateModule,
                           ImpactModule, InteriorModule, ResourceModule,
                           StellarModule, TerrainModule, TectonicsModule)
from aevum.spec.planet_spec import PlanetSpec


@dataclass
class Engine:
    spec: PlanetSpec
    registry: FeatureRegistry
    world: WorldState
    bus: EventBus
    scheduler: DeepTimeScheduler
    archive: WorldArchive

    @classmethod
    def build(cls, spec: PlanetSpec) -> "Engine":
        registry = build_registry()
        problems = registry.validate()
        if problems:
            raise ValueError("feature registry invalid:\n" + "\n".join(problems))

        grid = SphereGrid.fibonacci(spec.grid_cells, spec.radius_m)
        world = WorldState(grid=grid, spec=spec, time_myr=0.0)
        bus = EventBus()
        archive = WorldArchive(world, bus)
        scheduler = DeepTimeScheduler(world, bus, t_end_myr=spec.t_end_myr,
                                      dt_start_myr=10.0, dt_end_myr=40.0)

        # modules run sequentially within each macro step (truth-layer order)
        scheduler.add(StellarModule(), interval_myr=25.0)
        scheduler.add(InteriorModule(), interval_myr=25.0)
        scheduler.add(ImpactModule(), interval_myr=20.0)
        scheduler.add(TectonicsModule(), interval_myr=20.0)
        scheduler.add(TerrainModule(), interval_myr=20.0)
        scheduler.add(ClimateModule(), interval_myr=40.0,
                      trigger=_climate_trigger)
        scheduler.add(BiogeochemModule(), interval_myr=25.0)
        scheduler.add(BiosphereModule(), interval_myr=25.0)
        scheduler.add(ResourceModule(), interval_myr=50.0)

        return cls(spec=spec, registry=registry, world=world, bus=bus,
                   scheduler=scheduler, archive=archive)

    def run(self, n_frames: int = 18, progress: bool = False) -> None:
        snap_interval = self.spec.t_end_myr / max(n_frames, 1)
        state = {"next_snap": 0.0}

        def on_step(rec):
            if rec.time_myr >= state["next_snap"]:
                self.archive.capture(diagnostics=rec.diagnostics)
                state["next_snap"] += snap_interval
            if progress and rec.ran:
                print(f"  t={rec.time_myr:7.0f} Myr  ran={','.join(rec.ran)}")

        self.scheduler.run(on_step=on_step)
        self.archive.capture(diagnostics={"final": True})


def _climate_trigger(world: WorldState) -> bool:
    """Re-solve climate when geography / forcing drifted past a threshold."""
    if "climate.surface_temperature" not in world.fields:
        return True
    land = world.land_fraction()
    co2 = world.g("atmosphere.co2", 280e-6)
    lum = world.g("stellar.luminosity", 3.8e26)
    last_land = world.g("climate.solved_land", -1.0)
    last_co2 = world.g("climate.solved_co2", co2)
    last_lum = world.g("climate.solved_lum", lum)
    solved_mask = world.fields.get("ocean.solved_mask")
    mask_drift = False
    if solved_mask is not None and np.asarray(solved_mask).shape == (world.grid.n,):
        solved_ocean = np.asarray(solved_mask, dtype=np.float64) > 0.5
        current_ocean = world.ocean_mask()
        mismatch = current_ocean != solved_ocean
        area = world.grid.cell_area
        mask_drift = float(area[mismatch].sum() / area.sum()) > 0.01
    trig = (abs(land - last_land) > 0.05
            or mask_drift
            or co2 / max(last_co2, 1e-12) > 1.6 or last_co2 / max(co2, 1e-12) > 1.6
            or abs(lum - last_lum) / max(last_lum, 1e-12) > 0.08)
    if trig:
        world.set_g("climate.solved_land", land)
        world.set_g("climate.solved_co2", co2)
        world.set_g("climate.solved_lum", lum)
    return trig
