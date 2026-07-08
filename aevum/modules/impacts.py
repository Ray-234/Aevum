"""Extraterrestrial impacts (time-varying stochastic process + scaling laws).

Bombardment decays from an early-heavy phase.  Each step draws a Poisson number
of impactors; diameters follow a power law.  Small impacts only roughen the
crater field; large ones emit ``impact`` events that downstream modules
(climate, biosphere) can read as forcing -- so an impact-driven extinction is
emergent, not hard-coded.
"""
from __future__ import annotations

import numpy as np

from aevum.core.events import Event
from aevum.core.module import Module, StepResult
from aevum.core.provenance import Provenance


class ImpactModule(Module):
    name = "impacts"
    produces = ["impacts.flux", "impacts.crater_field"]
    fidelity = "scaling_law"
    interval_myr = 20.0

    BASE_RATE = 0.02          # impactors >1 km per Myr at late times (per planet)
    EARLY_BOOST = 40.0
    EARLY_TAU = 250.0         # Myr decay of heavy bombardment
    MAJOR_DIAM_KM = 50.0      # crater diameter threshold to emit an event

    def init_state(self, world, rng_key) -> None:
        world.set_field("impacts.crater_field", np.zeros(world.grid.n))
        world.set_g("impacts.last_major_myr", -1e9)

    def step(self, world, t, dt, forcing, rng_key) -> StepResult:
        grid = world.grid
        scale = world.spec.impact_flux_scale
        rate = self.BASE_RATE * scale * (1.0 + self.EARLY_BOOST * np.exp(-t / self.EARLY_TAU))
        rng = rng_key.generator()
        n = rng.poisson(rate * dt)

        crater = world.field("impacts.crater_field").copy()
        events: list[Event] = []
        last_major = world.g("impacts.last_major_myr")
        max_diam = 0.0
        for _ in range(int(n)):
            # Diameter power law (cumulative N(>D) ~ D^-2); km.
            d_km = 1.0 * (rng.pareto(2.0) + 1.0)
            cell = int(rng.integers(grid.n))
            crater[cell] += d_km
            max_diam = max(max_diam, d_km)
            if d_km >= self.MAJOR_DIAM_KM:
                last_major = t
                events.append(Event(
                    "impact", t, self.name, location=cell, magnitude=d_km,
                    params={"crater_km": round(d_km, 1),
                            "lat": round(float(grid.lat[cell]), 1),
                            "lon": round(float(grid.lon[cell]), 1),
                            "ocean": bool(world.ocean_mask()[cell])}))

        world.provenance.record(Provenance(
            "impacts.crater_field", self.name, self.fidelity, "1",
            direct_cause=f"rate {rate:.3f}/Myr; {int(n)} impactors this step"))
        diag = {"flux_per_myr": float(rate), "n_impacts": int(n),
                "max_diam_km": round(max_diam, 1)}
        return StepResult(
            state_delta={"fields": {"impacts.crater_field": crater},
                         "globals": {"impacts.flux": float(rate),
                                     "impacts.last_major_myr": last_major}},
            events=events, diagnostics=diag)
