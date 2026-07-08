"""Planetary interior thermal history (global box model).

dT/dt = (H_radiogenic(t) - Q_convective(T)) / heat_capacity

Tectonic vigour is a dimensionless function of mantle temperature (hotter ->
more vigorous convection -> faster plates -> more outgassing).  The regime
(mobile / sluggish / episodic / stagnant) is selected from vigour and the spec's
initial regime, so a cooling planet can transition between modes.
"""
from __future__ import annotations

import numpy as np

from aevum.core.module import Module, StepResult
from aevum.core.provenance import Provenance
from aevum.core.units import CONSTANTS, myr_to_seconds
from aevum.spec.planet_spec import TectonicRegime


class InteriorModule(Module):
    name = "interior"
    produces = ["interior.mantle_temperature", "interior.radiogenic_heat",
                "interior.convective_heat_loss", "interior.surface_heat_flow",
                "interior.tectonic_vigor", "interior.degassing_co2",
                "mantle.heat_anomaly", "mantle.upwelling_potential",
                "mantle.downwelling_potential", "lithosphere.thermal_thickness",
                "tectonics.rift_potential", "tectonics.plume_potential",
                "tectonics.regime"]
    fidelity = "box"
    interval_myr = 25.0

    # Reference scales.
    T_REF = 1300.0           # K, surface-side reference for convective loss
    DECAY_GYR = 3.8          # effective U/Th/K heat-production decay timescale
    PRESENT_RADIOGENIC_W = 2.0e13
    CONVECTIVE_REF_W = 4.7e13
    R1_PARAMETERS = {
        "heat_decay_per_25myr": 0.030,
        "heat_diffusion_mix_per_25myr": 0.110,
        "continental_insulation_source": 0.030,
        "ridge_heat_source": 0.020,
        "slab_cooling_source": 0.040,
        "old_ocean_cooling_source": 0.006,
        "thermal_lithosphere_min_m": 25000.0,
        "thermal_lithosphere_max_m": 240000.0,
    }

    def init_state(self, world, rng_key) -> None:
        spec = world.spec
        world.set_g("interior.mantle_temperature", spec.composition.initial_internal_temp_k)
        world.set_g("tectonics.regime_code", self._regime_code(spec.initial_tectonic_regime))
        n = world.n_cells
        world.set_field("mantle.heat_anomaly", np.zeros(n, dtype=np.float64))
        world.set_field("mantle.upwelling_potential", np.zeros(n, dtype=np.float64))
        world.set_field("mantle.downwelling_potential", np.zeros(n, dtype=np.float64))
        world.set_field("lithosphere.thermal_thickness",
                        np.full(n, 110000.0, dtype=np.float64))
        world.set_field("tectonics.rift_potential", np.zeros(n, dtype=np.float64))
        world.set_field("tectonics.plume_potential", np.zeros(n, dtype=np.float64))

    @staticmethod
    def _regime_code(regime: TectonicRegime) -> float:
        return {TectonicRegime.STAGNANT_LID: 0.0, TectonicRegime.SLUGGISH_LID: 1.0,
                TectonicRegime.EPISODIC_LID: 2.0, TectonicRegime.MOBILE_LID: 3.0}[regime]

    def _radiogenic(self, world, t: float) -> float:
        spec = world.spec
        # Treat radiogenic_abundance as a present-day Earth-relative inventory
        # and reconstruct higher early heat production from radioactive decay.
        present = (self.PRESENT_RADIOGENIC_W * spec.composition.radiogenic_abundance
                   * spec.composition.mass_earth)
        age_left_gyr = max(spec.t_end_myr - t, 0.0) / 1000.0
        return present * np.exp(age_left_gyr / self.DECAY_GYR)

    @staticmethod
    def _area_weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
        total = float(weights.sum())
        if total <= 0.0:
            return float(np.mean(values)) if values.size else 0.0
        return float(np.sum(values * weights) / total)

    @staticmethod
    def _neighbor_mean(grid, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype=np.float64)
        out = np.zeros(grid.n, dtype=np.float64)
        count = np.zeros(grid.n, dtype=np.float64)
        edges = grid.edges
        i, j = edges[:, 0], edges[:, 1]
        np.add.at(out, i, values[j])
        np.add.at(out, j, values[i])
        np.add.at(count, i, 1.0)
        np.add.at(count, j, 1.0)
        return np.divide(out, np.maximum(count, 1.0))

    @classmethod
    def _smooth(cls, grid, values: np.ndarray, passes: int = 1,
                mix: float = 0.35) -> np.ndarray:
        out = np.asarray(values, dtype=np.float64).copy()
        mix = float(np.clip(mix, 0.0, 1.0))
        for _ in range(max(0, int(passes))):
            out = (1.0 - mix) * out + mix * cls._neighbor_mean(grid, out)
        return out

    @staticmethod
    def _boundary_mask(world, *kinds: str) -> np.ndarray:
        mask = np.zeros(world.n_cells, dtype=bool)
        boundaries = world.networks.get("tectonics.boundaries", {})
        if not isinstance(boundaries, dict):
            return mask
        for kind in kinds:
            cells = np.asarray(boundaries.get(kind, []), dtype=int)
            cells = cells[(0 <= cells) & (cells < world.n_cells)]
            mask[cells] = True
        return mask

    def _tectonic_potential_fields(self, world, t: float, dt: float,
                                   vigor: float, regime_code: float
                                   ) -> tuple[dict[str, np.ndarray], dict[str, float]]:
        """Low-order deep-interior proxy fields for R1 tectonic refactor work.

        The field model is deliberately conservative: it is not a CFD mantle
        solver, but it preserves the main causal signs that later tectonics
        phases need to consume.  Continents insulate and accumulate positive
        anomaly, ridges focus upwelling, trenches/old slabs focus negative
        anomaly and downwelling, and lithospheric thickness responds to age and
        heat anomaly.
        """
        grid = world.grid
        area = grid.cell_area
        prev_heat = world.get_field("mantle.heat_anomaly", 0.0)
        ctype = world.get_field("crust.type", 0.0)
        age = np.maximum(world.get_field("crust.age_myr", 0.0), 0.0)
        stability = np.clip(world.get_field("crust.stability", 0.0), 0.0, 1.0)
        reworked = world.get_field("crust.reworked_age_myr", -1.0)

        continent = ctype >= 0.5
        ocean = ~continent
        mobile_factor = float(np.clip(regime_code / 3.0, 0.0, 1.0))
        vigor_factor = float(np.clip(vigor / 1.5, 0.15, 1.8))
        dt_factor = float(np.clip(dt / 25.0, 0.2, 2.0))

        ridge_mask = self._boundary_mask(world, "ridge", "divergent")
        trench_mask = self._boundary_mask(world, "trench", "subduction")
        passive_mask = self._boundary_mask(world, "passive_margin")
        suture_mask = self._boundary_mask(world, "suture", "collision")

        ridge_influence = np.clip(np.maximum(
            ridge_mask.astype(float),
            self._smooth(grid, ridge_mask.astype(float), passes=3, mix=0.55)),
            0.0, 1.0)
        trench_influence = np.clip(np.maximum(
            trench_mask.astype(float),
            self._smooth(grid, trench_mask.astype(float), passes=3, mix=0.55)),
            0.0, 1.0)
        passive_influence = np.clip(self._smooth(grid, passive_mask.astype(float),
                                                 passes=2, mix=0.50), 0.0, 1.0)
        suture_influence = np.clip(self._smooth(grid, suture_mask.astype(float),
                                                passes=2, mix=0.50), 0.0, 1.0)
        continent_interior = np.clip(self._smooth(grid, continent.astype(float),
                                                  passes=4, mix=0.50), 0.0, 1.0)

        old_ocean = np.zeros(grid.n, dtype=np.float64)
        old_ocean[ocean] = np.clip((age[ocean] - 55.0) / 130.0, 0.0, 1.0)
        continental_insulation = (
            continent.astype(float)
            * continent_interior
            * (0.35 + 0.65 * stability)
        )
        slab_cooling = trench_influence * (0.55 + 0.45 * old_ocean)

        source = (
            self.R1_PARAMETERS["continental_insulation_source"]
            * dt_factor * continental_insulation
            + self.R1_PARAMETERS["ridge_heat_source"]
            * dt_factor * mobile_factor * vigor_factor * ridge_influence
            - self.R1_PARAMETERS["slab_cooling_source"]
            * dt_factor * mobile_factor * slab_cooling
            - self.R1_PARAMETERS["old_ocean_cooling_source"]
            * dt_factor * old_ocean
        )
        source -= self._area_weighted_mean(source, area)

        diffusion_mix = float(np.clip(
            self.R1_PARAMETERS["heat_diffusion_mix_per_25myr"] * dt_factor
            * (0.65 + 0.20 * mobile_factor),
            0.02,
            0.24,
        ))
        heat = self._smooth(grid, prev_heat, passes=2, mix=diffusion_mix)
        decay = 1.0 - self.R1_PARAMETERS["heat_decay_per_25myr"] * dt_factor
        heat = np.clip(heat * max(decay, 0.0) + source, -1.0, 1.0)
        heat = self._smooth(grid, heat, passes=1, mix=0.08)

        positive_heat = np.clip(np.maximum(heat, 0.0) / 0.20, 0.0, 1.0)
        negative_heat = np.clip(np.maximum(-heat, 0.0) / 0.20, 0.0, 1.0)

        ocean_lithosphere = 34000.0 + 7100.0 * np.sqrt(np.clip(age, 0.0, 260.0))
        continental_lithosphere = 112000.0 + 98000.0 * stability
        lithosphere = np.where(continent, continental_lithosphere, ocean_lithosphere)
        lithosphere *= (1.0 - 0.23 * positive_heat + 0.12 * negative_heat)
        lithosphere *= (1.0 - 0.18 * ridge_influence)
        lithosphere = np.clip(
            lithosphere,
            self.R1_PARAMETERS["thermal_lithosphere_min_m"],
            self.R1_PARAMETERS["thermal_lithosphere_max_m"],
        )
        thin_lithosphere = np.clip((155000.0 - lithosphere) / 125000.0, 0.0, 1.0)

        upwelling = np.clip(
            0.60 * positive_heat
            + 0.25 * ridge_influence
            + 0.15 * thin_lithosphere,
            0.0,
            1.0,
        )
        upwelling = np.clip(self._smooth(grid, upwelling, passes=1, mix=0.18), 0.0, 1.0)

        downwelling = np.clip(
            0.74 * trench_influence
            + 0.18 * old_ocean
            + 0.26 * negative_heat,
            0.0,
            1.0,
        )
        downwelling = np.clip(self._smooth(grid, downwelling, passes=1, mix=0.12),
                              0.0, 1.0)

        recent_rework = np.zeros(grid.n, dtype=np.float64)
        reworked_mask = reworked >= 0.0
        if np.any(reworked_mask):
            recent_rework[reworked_mask] = np.clip(
                (260.0 - np.maximum(t - reworked[reworked_mask], 0.0)) / 260.0,
                0.0,
                1.0,
            )
        inherited_weakness = continent.astype(float) * np.clip(
            0.65 * (1.0 - stability) + 0.35 * recent_rework, 0.0, 1.0)
        margin_extension = np.clip(passive_influence + 0.55 * ridge_influence, 0.0, 1.0)
        rift_potential = np.clip(
            continent.astype(float)
            * (0.42 * upwelling
               + 0.35 * inherited_weakness
               + 0.20 * margin_extension
               + 0.10 * suture_influence
               - 0.22 * stability),
            0.0,
            1.0,
        )
        rift_potential = np.clip(self._smooth(grid, rift_potential, passes=1, mix=0.18),
                                 0.0, 1.0)

        plume_potential = np.clip(
            0.70 * positive_heat
            + 0.18 * upwelling
            + 0.12 * continental_insulation
            - 0.25 * downwelling,
            0.0,
            1.0,
        )
        plume_potential = np.clip(self._smooth(grid, plume_potential, passes=1, mix=0.10),
                                  0.0, 1.0)

        fields = {
            "mantle.heat_anomaly": heat,
            "mantle.upwelling_potential": upwelling,
            "mantle.downwelling_potential": downwelling,
            "lithosphere.thermal_thickness": lithosphere,
            "tectonics.rift_potential": rift_potential,
            "tectonics.plume_potential": plume_potential,
        }

        cont_heat = self._area_weighted_mean(heat[continent], area[continent]) if continent.any() else 0.0
        ocean_heat = self._area_weighted_mean(heat[ocean], area[ocean]) if ocean.any() else 0.0
        trench_cells = np.where(trench_mask)[0]
        far_cells = np.where(~self._smooth(grid, trench_mask.astype(float),
                                           passes=4, mix=0.60).astype(bool))[0]
        diag = {
            "heat_anomaly_mean": self._area_weighted_mean(heat, area),
            "heat_anomaly_min": float(np.min(heat)),
            "heat_anomaly_max": float(np.max(heat)),
            "continent_minus_ocean_heat": float(cont_heat - ocean_heat),
            "upwelling_p95": float(np.percentile(upwelling, 95)),
            "downwelling_p95": float(np.percentile(downwelling, 95)),
            "lithosphere_thickness_p50_km": float(np.percentile(lithosphere, 50) / 1000.0),
            "rift_potential_p95": float(np.percentile(rift_potential, 95)),
            "plume_potential_max": float(np.max(plume_potential)),
            "trench_downwelling_p90": (
                float(np.percentile(downwelling[trench_cells], 90))
                if trench_cells.size else 0.0
            ),
            "far_downwelling_p90": (
                float(np.percentile(downwelling[far_cells], 90))
                if far_cells.size else 0.0
            ),
        }
        return fields, diag

    def step(self, world, t, dt, forcing, rng_key) -> StepResult:
        spec = world.spec
        T = world.g("interior.mantle_temperature")
        H = self._radiogenic(world, t)

        # Convective heat loss scales super-linearly with super-adiabatic temp.
        # Larger planets have more area but proportionally more heat capacity;
        # water weakens the lithosphere and allows more efficient mobile-lid
        # cooling, while dry/strong lids retain heat.
        excess = max(T - self.T_REF, 0.0)
        strength = max(spec.crust_strength, 0.1)
        area_ratio = spec.composition.radius_earth ** 2
        water_lubrication = 0.75 + 0.25 * np.tanh(spec.composition.water_inventory_earth)
        cooling_efficiency = water_lubrication / strength
        Q = (self.CONVECTIVE_REF_W * area_ratio
             * (excess / 1700.0) ** 1.3 * cooling_efficiency)

        # Effective bulk heat capacity (J/K) of the mantle.
        C = 7.0e27 * spec.composition.mass_earth
        dT = (H - Q) * myr_to_seconds(dt) / C
        T = float(np.clip(T + dT, 800.0, 4500.0))

        heat_flow = Q / max(spec.surface_area_m2, 1.0)
        vigor = float(np.clip((excess / 1700.0) * cooling_efficiency, 0.0, 5.0))
        # Regime transition: cooling can drop a mobile lid to stagnant.
        code = world.g("tectonics.regime_code")
        new_code = code
        if spec.initial_tectonic_regime == TectonicRegime.MOBILE_LID:
            if vigor < 0.25:
                new_code = 1.0 if vigor > 0.1 else 0.0
            else:
                new_code = 3.0
        regime_changed = abs(new_code - code) > 0.5

        degassing = 6.0e18 * vigor * spec.composition.mass_earth
        potential_fields, potential_diag = self._tectonic_potential_fields(
            world, t, dt, vigor, new_code)

        globals_ = {
            "interior.mantle_temperature": T,
            "interior.radiogenic_heat": H,
            "interior.convective_heat_loss": Q,
            "interior.surface_heat_flow": heat_flow,
            "interior.tectonic_vigor": vigor,
            "interior.degassing_co2": degassing,
            "tectonics.regime_code": new_code,
        }
        world.provenance.record(Provenance(
            "interior.mantle_temperature", self.name, self.fidelity, "K",
            direct_cause=f"radiogenic {H:.2e} W vs convective loss {Q:.2e} W; "
            f"surface heat flow {heat_flow:.3f} W/m^2"))

        events = []
        if regime_changed:
            from aevum.core.events import Event
            events.append(Event(
                type="regime_transition", time_myr=t, producer=self.name,
                magnitude=new_code,
                params={"from": code, "to": new_code, "mantle_T": T},
            ))
        diag = {"mantle_T": T, "vigor": vigor, "regime_code": new_code,
                "radiogenic_W": H, "convective_W": Q,
                "heat_flow_Wm2": heat_flow}
        diag.update(potential_diag)
        return StepResult(state_delta={"globals": globals_, "fields": potential_fields},
                          events=events, diagnostics=diag)
