"""Biogeochemistry (COPSE-style long-term box model).

Carbon: volcanic/metamorphic degassing adds CO2; carbonate-silicate weathering
(temperature & runoff dependent) removes it.  This is the thermostat that keeps
liquid water under a faint young sun and drives icehouse/greenhouse swings.

Oxygen: organic-carbon burial is the source, oxidative weathering + volcanic
reductants the sink.  A persistent source>sink yields a Great-Oxidation-style
rise -- emergent, not scheduled.

Carbon is conserved across atmosphere + buried reservoirs (validation audits it).
"""
from __future__ import annotations

import numpy as np

from aevum.core.events import Event
from aevum.core.module import Module, StepResult
from aevum.core.provenance import Provenance

ATM_MOL = 1.8e20          # ~ moles in a 1-bar atmosphere column inventory proxy
CO2_REF_BAR = 280e-6
W0 = 6.0e18               # baseline weathering (mol/Myr), tuned to balance degassing
T_REF = 288.0
RUNOFF_REF = 400.0
SEAFLOOR_WEATHERING_SCALE = 0.45


class BiogeochemModule(Module):
    name = "biogeochem"
    produces = ["atmosphere.co2", "biogeochem.carbon", "biogeochem.oxygen",
                "biogeochem.weathering_flux"]
    fidelity = "copse_box"
    interval_myr = 25.0

    def init_state(self, world, rng_key) -> None:
        spec = world.spec
        co2_mol = spec.atmosphere.co2_bar * ATM_MOL
        world.set_g("bgc.carbon_atm_mol", co2_mol)
        world.set_g("bgc.carbon_buried_mol", 50.0 * ATM_MOL * CO2_REF_BAR * 1000.0)
        world.set_g("bgc.oxygen_mol", spec.atmosphere.o2_fraction * ATM_MOL)
        world.set_g("atmosphere.co2", spec.atmosphere.co2_bar)
        world.set_g("bgc.goe_done", 0.0)
        world.set_g("bgc.carbon_total_init",
                    world.g("bgc.carbon_atm_mol") + world.g("bgc.carbon_buried_mol"))

    def step(self, world, t, dt, forcing, rng_key) -> StepResult:
        spec = world.spec
        grid = world.grid
        area = grid.cell_area

        carbon_atm = world.g("bgc.carbon_atm_mol")
        carbon_buried = world.g("bgc.carbon_buried_mol")
        oxygen = world.g("bgc.oxygen_mol")
        co2_bar = carbon_atm / ATM_MOL

        T = world.get_field("climate.surface_temperature", 288.0)
        runoff = world.get_field("climate.runoff", RUNOFF_REF)
        meanT = float(np.average(T, weights=area))
        land = world.land_mask()
        mean_runoff = float(np.average(runoff[land], weights=area[land])) if land.any() else 0.0

        degass = world.g("interior.degassing_co2", W0)      # mol/Myr

        # carbonate-silicate weathering thermostat (bounded sensitivities)
        fT = np.exp(np.clip((meanT - T_REF) / 13.7, -4.0, 4.0))
        fR = max(mean_runoff / RUNOFF_REF, 0.0) ** 0.65
        fC = max(co2_bar / CO2_REF_BAR, 1e-6) ** 0.42
        land_factor = max(world.land_fraction(), 0.02) / 0.29
        litho = world.get_field("terrain.lithology", 0.0)
        reactive = float(np.average((litho != 1.0)[land], weights=area[land])) if land.any() else 0.0
        relief = world.get_field("terrain.elevation_m", 0.0) - world.sea_level
        relief_factor = 1.0
        if land.any():
            relief_factor += 0.35 * np.clip(np.percentile(relief[land], 90) / 2500.0, 0.0, 2.0)
        lithology_factor = 0.75 + 0.55 * reactive
        weathering = W0 * fT * fR * fC * land_factor * lithology_factor * relief_factor

        ocean_fraction = 1.0 - world.land_fraction()
        carbonate_saturation = max(co2_bar / max(0.0015 * world.g("atmosphere.pressure", 1.0), 1e-8) - 1.0, 0.0)
        temp_ocean = np.clip((meanT - 270.0) / 28.0, 0.0, 1.5)
        carbonate_burial = (
            1.10 * W0 * ocean_fraction
            * carbonate_saturation ** 0.75
            * temp_ocean
        )
        seafloor_weathering = (
            SEAFLOOR_WEATHERING_SCALE * W0 * ocean_fraction
            * max(co2_bar / CO2_REF_BAR, 1e-6) ** 0.32
            * np.clip(world.g("interior.tectonic_vigor", 1.0), 0.1, 2.0) ** 0.25
        )

        # Closed two-reservoir update (atmosphere/ocean <-> deep), so total carbon
        # is conserved exactly; the surface reservoir is flux-limited per step.
        removal = weathering + carbonate_burial + seafloor_weathering
        flux = (degass - removal) * dt                         # +ve: to atmosphere
        # Deep-time macro steps are tens of Myr; without damping the carbon box
        # can numerically overshoot into CO2 boom/bust cycles.  Limit each update
        # to a modest fraction of the surface reservoir while still conserving
        # carbon exactly between reservoirs.
        flux = float(np.clip(flux, -0.28 * carbon_atm, 0.28 * carbon_atm))
        flux = float(np.clip(flux, -carbon_atm + 1e-12, carbon_buried))
        carbon_atm = carbon_atm + flux
        carbon_buried = carbon_buried - flux
        co2_bar = carbon_atm / ATM_MOL

        # oxygen budget (relaxation toward a productivity-set equilibrium, with
        # reductant sinks); bounded so it cannot run away.
        npp = world.get_field("biosphere.npp", 0.0)         # kg C / m^2 / yr
        npp_total = float(np.sum(npp * area))               # kg C / yr
        biomass_norm = np.clip(npp_total / 6.0e13, 0.0, 1.5)
        oxygenic = (spec.biosphere.allow_oxygenic_photosynthesis
                    and world.g("bio.oxygenic_photosynthesis", 0.0) > 0.5)
        o2_frac = oxygen / ATM_MOL
        prod = 0.018 * biomass_norm if oxygenic else 0.0    # per Myr toward target
        sink = 0.05 * o2_frac + 0.02 * (degass / 6.0e18)    # oxidative + reductant
        o2_frac = float(np.clip(o2_frac + (prod - sink) * dt, 0.0, 0.30))
        oxygen = o2_frac * ATM_MOL

        events: list[Event] = []
        if world.g("bgc.goe_done") < 0.5 and o2_frac > 0.01:
            world.set_g("bgc.goe_done", 1.0)
            events.append(Event("great_oxidation", t, self.name, magnitude=o2_frac,
                                params={"o2_fraction": round(o2_frac, 4)}))

        world.provenance.record(Provenance(
            "atmosphere.co2", self.name, self.fidelity, "bar",
            direct_cause=f"degassing {degass:.2e} vs weathering {weathering:.2e} "
            f"+ carbonate burial {carbonate_burial:.2e} + seafloor weathering "
            f"{seafloor_weathering:.2e} mol/Myr "
            f"(thermostat at meanT {meanT-273.15:.1f} C)"))

        globals_ = {
            "atmosphere.co2": co2_bar,
            "bgc.carbon_atm_mol": carbon_atm,
            "bgc.carbon_buried_mol": carbon_buried,
            "bgc.oxygen_mol": oxygen,
            "biogeochem.oxygen_fraction": o2_frac,
            "biogeochem.weathering_flux": weathering,
        }
        diag = {"co2_ppm_eq": round(co2_bar / 1e-6, 1), "o2_fraction": round(o2_frac, 4),
                "weathering": weathering,
                "carbonate_burial": carbonate_burial,
                "seafloor_weathering": seafloor_weathering,
                "meanT_C": round(meanT - 273.15, 2)}
        return StepResult(state_delta={"globals": globals_}, events=events,
                          diagnostics=diag)
