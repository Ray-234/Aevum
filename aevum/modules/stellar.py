"""Stellar & orbital forcing.

Faint-young-sun brightening + annual-mean insolation as a function of latitude
and obliquity (or angle from the substellar point for tidally locked worlds).
"""
from __future__ import annotations

import numpy as np

from aevum.core.module import Module, StepResult
from aevum.core.provenance import Provenance
from aevum.core.units import CONSTANTS
from aevum.modules.geometry import substellar_cos


class StellarModule(Module):
    name = "stellar"
    produces = ["stellar.luminosity", "stellar.surface_flux",
                "stellar.mean_flux", "orbit.obliquity", "orbit.eccentricity",
                "orbit.eccentricity_factor"]
    fidelity = "analytic"
    interval_myr = 25.0

    def init_state(self, world, rng_key) -> None:
        spec = world.spec
        world.set_g("orbit.obliquity", spec.orbit.obliquity_deg)
        world.set_g("orbit.eccentricity", spec.orbit.eccentricity)
        world.set_g("atmosphere.pressure", spec.atmosphere.surface_pressure_bar)

    def _luminosity(self, world, t: float) -> float:
        spec = world.spec
        l_now = spec.stellar.luminosity_solar_now * CONSTANTS.SOLAR_LUMINOSITY
        # Faint-main-sequence brightening, normalized so luminosity_solar_now is
        # reached at the end of the requested planet history.  Lower-mass stars
        # brighten more slowly on the main sequence, so M-dwarf presets should
        # not inherit the full solar 40% correction.
        age0 = max(spec.stellar.age_at_start_gyr, 0.0)
        age_end = age0 + max(spec.t_end_myr, 1e-9) / 1000.0
        age = age0 + max(t, 0.0) / 1000.0
        frac = np.clip(age / max(age_end, 1e-9), 0.0, 1.0)
        brightening = 0.40 * np.clip(spec.stellar.mass_solar, 0.2, 2.0) ** 0.6
        return l_now / (1.0 + brightening * (1.0 - frac))

    def step(self, world, t, dt, forcing, rng_key) -> StepResult:
        spec = world.spec
        grid = world.grid
        lum = self._luminosity(world, t)
        d = spec.orbit.semi_major_axis_au * CONSTANTS.AU
        ecc = np.clip(world.g("orbit.eccentricity"), 0.0, 0.95)
        # Time-mean inverse-square forcing over an eccentric orbit scales as
        # 1/sqrt(1-e^2).  The result remains the sphere-averaged TOA flux.
        ecc_factor = 1.0 / np.sqrt(max(1.0 - ecc * ecc, 1e-6))
        s0 = lum / (16.0 * np.pi * d * d) * ecc_factor

        if spec.orbit.tidally_locked:
            mu = substellar_cos(grid.xyz, np.array([1.0, 0.0, 0.0]))
            # Day side scales with cos(zenith); night side gets ~0 direct flux.
            flux = 4.0 * s0 * np.clip(mu, 0.0, None)
        else:
            eps = np.radians(world.g("orbit.obliquity"))
            x = np.sin(np.radians(grid.lat))
            p2 = 0.5 * (3.0 * x * x - 1.0)
            # Second-Legendre annual-mean insolation; s2 flips sign past ~54 deg.
            # Low-order annual-mean insolation.  Keep the equator-to-pole
            # contrast moderate because unresolved atmospheric/ocean transport is
            # still a reduced proxy in v0.1.
            s2 = -0.30 * (1.0 - 1.5 * np.sin(eps) ** 2)
            flux = s0 * (1.0 + s2 * p2)
            flux = np.clip(flux, 0.0, None)

        # Keep the area-weighted global mean exactly tied to the orbital energy
        # budget.  This removes small grid and clipping biases before climate
        # sees the field.
        mean_flux = float(np.average(flux, weights=grid.cell_area))
        if mean_flux > 0.0:
            flux = flux * (s0 / mean_flux)
            mean_flux = float(np.average(flux, weights=grid.cell_area))

        fields = {"stellar.surface_flux": flux}
        globals_ = {"stellar.luminosity": lum,
                    "stellar.mean_flux": mean_flux,
                    "orbit.eccentricity_factor": float(ecc_factor)}
        world.provenance.record(Provenance(
            "stellar.surface_flux", self.name, self.fidelity, "W m^-2",
            direct_cause=f"L={lum:.3e} W at d={spec.orbit.semi_major_axis_au} AU"
            + (", tidally locked" if spec.orbit.tidally_locked else
               f", obliquity {world.g('orbit.obliquity'):.1f} deg")
            + f", eccentricity factor {ecc_factor:.3f}"))
        diag = {"luminosity_Lsun": lum / CONSTANTS.SOLAR_LUMINOSITY,
                "mean_flux_Wm2": mean_flux,
                "eccentricity_factor": float(ecc_factor)}
        return StepResult(state_delta={"fields": fields, "globals": globals_},
                          diagnostics=diag)
