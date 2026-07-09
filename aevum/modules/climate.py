"""Atmosphere & ocean (energy-balance climate).

A diffusive energy-balance model gives surface temperature with an ice-albedo
feedback (so frozen worlds can runaway to snowball).  Banded winds drive a
moisture-advection precipitation model with orographic enhancement and leeward
rain shadows.  Sea ice and ice sheets follow from temperature; ocean currents
are a wind-driven proxy that is zeroed over land (never crosses continents).

This is the expensive module: the scheduler only re-solves it when geography,
orbit, atmosphere or topography have drifted past a threshold.
"""
from __future__ import annotations

import numpy as np

from aevum.core.module import Module, StepResult
from aevum.core.provenance import Provenance
from aevum.core.units import CONSTANTS


class ClimateModule(Module):
    name = "climate"
    produces = ["climate.surface_temperature", "climate.precipitation",
                "climate.evaporation", "climate.runoff", "atmosphere.wind",
                "ocean.currents", "cryosphere.sea_ice", "cryosphere.ice_sheet",
                "climate.seasonal_temperature", "climate.temperature_seasonality",
                "climate.seasonal_sst", "climate.ocean_heat_flux",
                "climate.ocean_evaporation_feedback",
                "climate.coupling_residual",
                "climate.seasonal_insolation_anomaly",
                "climate.surface_heat_capacity_class",
                "climate.land_thermal_anomaly",
                "climate.ocean_mixed_layer_thermal_anomaly",
                "climate.elevation_lapse_cooling",
                "climate.snow_ice_albedo_support",
                "climate.sst_gradient_support",
                "climate.same_latitude_sst_anomaly",
                "climate.land_sea_thermal_contrast",
                "climate.seasonal_precipitation",
                "climate.precipitation_seasonality", "climate.monsoon_index",
                "climate.dry_season_length", "climate.wet_season_peak",
                "climate.moisture_convergence",
                "climate.orographic_precipitation",
                "climate.monsoon_rainfall_corridor",
                "climate.storm_track_rainfall_corridor",
                "climate.rain_shadow_index",
                "climate.regional_precipitation_response",
                "atmosphere.moisture_flow_source",
                "atmosphere.moisture_flow_pathway",
                "atmosphere.moisture_source_basin_id",
                "climate.moisture_flow_network_id",
                "climate.moisture_flow_precipitation_response",
                "climate.moisture_budget_region_id",
                "climate.precipitation_response_region_id",
                "climate.receiver_catchment_id",
                "climate.source_basin_supply_index",
                "climate.receiver_catchment_supply_balance",
                "climate.receiver_supply_precipitation_feedback",
                "climate.hydro_coupling_residual",
                "climate.hydro_feedback_iteration_delta",
                "cryosphere.seasonal_sea_ice",
                "cryosphere.seasonal_snow",
                "cryosphere.snow_persistence",
                "climate.seasonal_cloud_albedo_proxy",
                "climate.cloud_albedo_proxy",
                "biosphere.vegetation_climate_feedback",
                "climate.continentality", "atmosphere.seasonal_wind",
                "atmosphere.itcz_latitude", "atmosphere.itcz_intensity",
                "atmosphere.storm_track_intensity",
                "atmosphere.background_seasonal_wind",
                "atmosphere.land_sea_pressure_proxy",
                "atmosphere.seasonal_pressure_proxy",
                "atmosphere.pressure_center_support",
                "atmosphere.pressure_center_id",
                "atmosphere.stationary_wave_pressure_support",
                "atmosphere.pressure_genesis_source",
                "atmosphere.pressure_genesis_wave_transfer",
                "atmosphere.ocean_pressure_low_source_support",
                "atmosphere.ocean_pressure_high_source_support",
                "atmosphere.land_pressure_source_support",
                "atmosphere.terrain_pressure_wave_source_support",
                "atmosphere.precipitation_pressure_feedback",
                "atmosphere.hydro_coupled_wind_anomaly",
                "atmosphere.thermal_wind_anomaly",
                "atmosphere.orographic_wind_anomaly",
                "atmosphere.geographic_circulation_index",
                "atmosphere.moisture_access",
                "atmosphere.monsoon_potential",
                "atmosphere.source_ocean_warmth",
                "atmosphere.terrain_blocking",
                "ocean.current_heat_transport", "ocean.upwelling",
                "ocean.basin_id", "ocean.solved_mask",
                "ocean.gyre_id", "ocean.current_streamfunction",
                "ocean.boundary_current_type", "ocean.strait_exchange",
                "ocean.wind_stress_current_response",
                "ocean.sst_anomaly",
                "climate.continent_id", "climate.continent_interiority",
                "climate.coast_orientation", "climate.coast_distance",
                "climate.coast_strength", "climate.coast_facing_east",
                "ocean.shelf_index", "ocean.strait_index",
                "terrain.barrier_index", "terrain.wind_gap_index"]
    fidelity = "ebm"
    interval_myr = 40.0

    A0 = 195.0          # OLR intercept (W/m^2) at 0 C
    B = 2.7             # OLR sensitivity (W/m^2/K)
    CO2_REF = 280e-6    # bar
    EBM_ITERS = 220
    HYDRO_FEEDBACK_ITERS = 3
    OCEAN_ATMOSPHERE_ITERS = 3
    SEASON_NAMES = ("DJF", "MAM", "JJA", "SON")

    def init_state(self, world, rng_key) -> None:
        self._edge_cache(world.grid)

    # ------------------------------------------------------------------
    def _edge_cache(self, grid):
        if not hasattr(self, "_edges") or self._edges_grid is not grid:
            self._edges = grid.edges
            self._edges_grid = grid
            deg = np.zeros(grid.n)
            np.add.at(deg, self._edges[:, 0], 1)
            np.add.at(deg, self._edges[:, 1], 1)
            self._deg = np.maximum(deg, 1)

    def _neighbor_mean(self, grid, f):
        i, j = self._edges[:, 0], self._edges[:, 1]
        acc = np.zeros_like(f)
        np.add.at(acc, i, f[j])
        np.add.at(acc, j, f[i])
        return acc / self._deg

    def _tangent_basis(self, grid):
        z = np.array([0.0, 0.0, 1.0])
        east = np.cross(z, grid.xyz)
        en = np.linalg.norm(east, axis=1, keepdims=True)
        east = np.where(en > 1e-6, east / np.maximum(en, 1e-9), np.array([1.0, 0, 0]))
        north = np.cross(grid.xyz, east)
        return east, north

    def _project_tangent(self, grid, vectors):
        vectors = np.asarray(vectors, dtype=np.float64)
        return vectors - np.sum(vectors * grid.xyz, axis=-1, keepdims=True) * grid.xyz

    def _tidally_locked_wind(self, world):
        grid = world.grid
        sub = np.array([1.0, 0.0, 0.0])
        mu = grid.xyz @ sub
        # surface flow from substellar high to antistellar low
        toward = sub[None, :] - mu[:, None] * grid.xyz
        tn = np.linalg.norm(toward, axis=1, keepdims=True)
        return -8.0 * np.where(tn > 1e-6, toward / np.maximum(tn, 1e-9), 0.0)

    def _seasonal_circulation(self, world):
        grid = world.grid
        if world.spec.orbit.tidally_locked:
            wind = np.repeat(self._tidally_locked_wind(world)[None, :, :], 4, axis=0)
            zeros = np.zeros((4, grid.n), dtype=np.float64)
            return wind, np.zeros(4, dtype=np.float64), zeros, zeros

        east, north = self._tangent_basis(grid)
        lat = grid.lat
        eps_deg = world.g("orbit.obliquity", world.spec.orbit.obliquity_deg)
        eps = np.radians(eps_deg)
        phases = np.array([-0.5 * np.pi, 0.0, 0.5 * np.pi, np.pi])
        decl = np.arcsin(np.clip(np.sin(eps) * np.sin(phases), -1.0, 1.0))
        decl_deg = np.degrees(decl)
        season_ratio = decl_deg / max(abs(eps_deg), 1.0)
        itcz_lat = np.clip(0.62 * decl_deg, -20.0, 20.0)

        rot = max(world.spec.orbit.rotation_period_hours, 1.0)
        rotation_factor = np.clip((24.0 / rot) ** 0.20, 0.65, 1.35)
        hadley_edge = np.clip(30.0 / rotation_factor, 22.0, 42.0)
        polar_edge = np.clip(62.0 / rotation_factor, 52.0, 72.0)

        winds: list[np.ndarray] = []
        itcz_intensity: list[np.ndarray] = []
        storm_tracks: list[np.ndarray] = []
        for shift, ratio in zip(itcz_lat, season_ratio):
            # Rain belts migrate seasonally, but the near-surface trade-wind
            # and midlatitude-westerly bands do not translate one-for-one with
            # the ITCZ.  Keeping the dynamical wind belts more anchored avoids
            # turning subtropical trades into broad seasonal westerlies.
            belt_shift = 0.36 * shift
            rel = lat - belt_shift
            convergence_rel = lat - shift
            abs_rel = np.abs(rel)

            # These exported winds are calibrated against near-surface Earth
            # reanalysis, not upper-troposphere jet speeds.  Keep the zonal
            # template continuous so land-sea and orographic anomalies can
            # shape the final wind field instead of being drowned by hard
            # latitude-band edges.
            trade_core = np.exp(-((abs_rel / max(0.72 * hadley_edge, 1.0)) ** 4))
            westerly_center = (
                0.5 * (hadley_edge + polar_edge)
                + np.where(lat < 0.0, 5.2, 0.0)
            )
            base_westerly_width = max(0.34 * (polar_edge - hadley_edge), 10.0)
            westerly_width = base_westerly_width * np.where(lat < 0.0, 1.18, 1.0)
            westerly_core = np.exp(-((abs_rel - westerly_center) / westerly_width) ** 2)
            polar_core = np.exp(-((abs_rel - (polar_edge + 8.0)) / 13.0) ** 2)
            winter_boost = np.clip(1.0 - 0.20 * np.sign(lat) * ratio, 0.80, 1.20)
            hemisphere_westerly = np.where(lat < 0.0, 1.22, 0.86)
            u = (
                -3.8 * trade_core
                + 6.5 * hemisphere_westerly * westerly_core * winter_boost
                - 2.4 * polar_core
            )

            convergence = np.exp(-(convergence_rel / 34.0) ** 2)
            v = -1.65 * np.tanh(convergence_rel / 13.0) * convergence
            winds.append(u[:, None] * east + v[:, None] * north)

            itcz = np.exp(-(convergence_rel / 13.0) ** 2)
            itcz_intensity.append(itcz)

            north_center = 45.0 + 4.0 * ratio
            south_center = -45.0 + 4.0 * ratio
            north_strength = 1.0 - 0.32 * ratio
            south_strength = 1.0 + 0.32 * ratio
            storm = (
                north_strength * np.exp(-((lat - north_center) / 16.0) ** 2)
                + south_strength * np.exp(-((lat - south_center) / 16.0) ** 2)
            )
            storm_tracks.append(storm)

        return (
            np.asarray(winds, dtype=np.float64),
            itcz_lat.astype(np.float64),
            np.asarray(itcz_intensity, dtype=np.float64),
            np.asarray(storm_tracks, dtype=np.float64),
        )

    def _winds(self, world):
        seasonal_wind, _, _, _ = self._seasonal_circulation(world)
        return seasonal_wind.mean(axis=0)

    def _solar_factor(self, world):
        # already encoded in stellar.surface_flux
        return world.get_field("stellar.surface_flux", 340.0)

    def _graph_gradient_vectors(self, grid, values):
        values = np.asarray(values, dtype=np.float64)
        i, j = self._edges[:, 0], self._edges[:, 1]
        d_ij = grid.xyz[j] - grid.xyz[i]
        dir_i = d_ij - np.sum(d_ij * grid.xyz[i], axis=1, keepdims=True) * grid.xyz[i]
        n_i = np.linalg.norm(dir_i, axis=1, keepdims=True)
        dir_i = np.where(n_i > 1e-12, dir_i / np.maximum(n_i, 1e-12), 0.0)

        d_ji = -d_ij
        dir_j = d_ji - np.sum(d_ji * grid.xyz[j], axis=1, keepdims=True) * grid.xyz[j]
        n_j = np.linalg.norm(dir_j, axis=1, keepdims=True)
        dir_j = np.where(n_j > 1e-12, dir_j / np.maximum(n_j, 1e-12), 0.0)

        grad = np.zeros((grid.n, 3), dtype=np.float64)
        delta_ij = values[j] - values[i]
        np.add.at(grad, i, delta_ij[:, None] * dir_i)
        np.add.at(grad, j, -delta_ij[:, None] * dir_j)
        grad = grad / self._deg[:, None]
        return self._project_tangent(grid, grad)

    def _latitude_band_anomaly(self, grid, values, *, bin_width_deg=5.0):
        values = np.asarray(values, dtype=np.float64)
        out = np.zeros_like(values, dtype=np.float64)
        band_id = np.floor((grid.lat + 90.0) / bin_width_deg).astype(int)
        for band in np.unique(band_id):
            cells = (band_id == int(band)) & np.isfinite(values)
            if not cells.any():
                continue
            mean = float(np.average(values[cells], weights=grid.cell_area[cells]))
            out[cells] = values[cells] - mean
        return out

    def _smooth_vector_field(self, grid, values, passes: int = 1, alpha: float = 0.15):
        out = np.asarray(values, dtype=np.float64).copy()
        for _ in range(passes):
            smoothed = np.column_stack([
                self._neighbor_mean(grid, out[:, k]) for k in range(3)
            ])
            out = (1.0 - alpha) * out + alpha * smoothed
            out = self._project_tangent(grid, out)
        return out

    def _masked_neighbor_mean(self, grid, values, mask):
        values = np.asarray(values, dtype=np.float64)
        mask = np.asarray(mask, dtype=bool)
        i, j = self._edges[:, 0], self._edges[:, 1]
        valid = mask[i] & mask[j]
        acc = np.zeros(grid.n, dtype=np.float64)
        deg = np.zeros(grid.n, dtype=np.float64)
        if valid.any():
            vi = i[valid]
            vj = j[valid]
            np.add.at(acc, vi, values[vj])
            np.add.at(acc, vj, values[vi])
            np.add.at(deg, vi, 1.0)
            np.add.at(deg, vj, 1.0)
        out = values.copy()
        has = deg > 0.0
        out[has] = acc[has] / deg[has]
        out[~mask] = 0.0
        return out

    def _smooth_field_masked(self, grid, values, mask, passes: int = 1,
                             alpha: float = 0.2):
        out = np.asarray(values, dtype=np.float64).copy()
        mask = np.asarray(mask, dtype=bool)
        out[~mask] = 0.0
        for _ in range(passes):
            mean = self._masked_neighbor_mean(grid, out, mask)
            out = np.where(mask, (1.0 - alpha) * out + alpha * mean, 0.0)
        return out

    def _weighted_neighbor_mean_masked(self, grid, values, mask, conductance):
        values = np.asarray(values, dtype=np.float64)
        mask = np.asarray(mask, dtype=bool)
        conductance = np.asarray(conductance, dtype=np.float64)
        if conductance.shape != (grid.n,):
            conductance = np.ones(grid.n, dtype=np.float64)
        i, j = self._edges[:, 0], self._edges[:, 1]
        weight = np.sqrt(
            np.clip(conductance[i], 0.0, None)
            * np.clip(conductance[j], 0.0, None)
        )
        valid = mask[i] & mask[j] & (weight > 1.0e-9)
        acc = np.zeros(grid.n, dtype=np.float64)
        wsum = np.zeros(grid.n, dtype=np.float64)
        if valid.any():
            vi = i[valid]
            vj = j[valid]
            w = weight[valid]
            np.add.at(acc, vi, w * values[vj])
            np.add.at(acc, vj, w * values[vi])
            np.add.at(wsum, vi, w)
            np.add.at(wsum, vj, w)
        out = values.copy()
        has = wsum > 1.0e-12
        out[has] = acc[has] / wsum[has]
        out[~mask] = 0.0
        return out

    def _diffuse_field_weighted(
        self,
        grid,
        values,
        mask,
        conductance,
        *,
        passes: int = 1,
        alpha: float = 0.2,
        seed_retention: float = 0.0,
        seed=None,
    ):
        mask = np.asarray(mask, dtype=bool)
        out = np.where(mask, np.asarray(values, dtype=np.float64), 0.0)
        seed_values = out if seed is None else np.asarray(seed, dtype=np.float64)
        seed_values = np.where(mask, seed_values, 0.0)
        conductance = np.asarray(conductance, dtype=np.float64)
        if conductance.shape != (grid.n,):
            conductance = np.ones(grid.n, dtype=np.float64)
        conductance = np.where(mask, np.clip(conductance, 0.0, 1.5), 0.0)
        conductance = np.maximum(conductance, 0.04 * mask.astype(np.float64))
        for _ in range(max(int(passes), 0)):
            mean = self._weighted_neighbor_mean_masked(
                grid, out, mask, conductance)
            local_alpha = alpha * np.clip(0.45 + 0.55 * conductance, 0.0, 1.0)
            out = np.where(
                mask,
                (1.0 - local_alpha) * out + local_alpha * mean,
                0.0,
            )
            if seed_retention > 0.0:
                out = np.where(
                    mask,
                    (1.0 - seed_retention) * out
                    + seed_retention * seed_values,
                    0.0,
                )
        return out

    def _directional_neighbor_mean_masked(
        self,
        grid,
        values,
        mask,
        conductance,
        axis,
        *,
        direction_strength: float = 1.25,
        direction_power: float = 2.0,
    ):
        values = np.asarray(values, dtype=np.float64)
        mask = np.asarray(mask, dtype=bool)
        conductance = np.asarray(conductance, dtype=np.float64)
        if conductance.shape != (grid.n,):
            conductance = np.ones(grid.n, dtype=np.float64)
        axis = np.asarray(axis, dtype=np.float64)
        if axis.shape != (grid.n, 3):
            east, _ = self._tangent_basis(grid)
            axis = east
        axis = self._project_tangent(grid, axis)
        axis_norm = np.linalg.norm(axis, axis=1, keepdims=True)
        east, _ = self._tangent_basis(grid)
        axis = np.where(axis_norm > 1.0e-9, axis / np.maximum(axis_norm, 1.0e-9), east)

        i, j = self._edges[:, 0], self._edges[:, 1]
        edge = grid.xyz[j] - grid.xyz[i]
        dir_i = edge - np.sum(edge * grid.xyz[i], axis=1, keepdims=True) * grid.xyz[i]
        norm_i = np.linalg.norm(dir_i, axis=1, keepdims=True)
        dir_i = np.where(norm_i > 1.0e-12, dir_i / np.maximum(norm_i, 1.0e-12), 0.0)
        dir_j = -edge - np.sum((-edge) * grid.xyz[j], axis=1, keepdims=True) * grid.xyz[j]
        norm_j = np.linalg.norm(dir_j, axis=1, keepdims=True)
        dir_j = np.where(norm_j > 1.0e-12, dir_j / np.maximum(norm_j, 1.0e-12), 0.0)

        align = 0.5 * (
            np.abs(np.sum(dir_i * axis[i], axis=1))
            + np.abs(np.sum(dir_j * axis[j], axis=1))
        )
        directional_weight = 0.25 + direction_strength * np.clip(
            align, 0.0, 1.0) ** direction_power
        weight = (
            np.sqrt(
                np.clip(conductance[i], 0.0, None)
                * np.clip(conductance[j], 0.0, None)
            )
            * directional_weight
        )
        valid = mask[i] & mask[j] & (weight > 1.0e-9)
        acc = np.zeros(grid.n, dtype=np.float64)
        wsum = np.zeros(grid.n, dtype=np.float64)
        if valid.any():
            vi = i[valid]
            vj = j[valid]
            w = weight[valid]
            np.add.at(acc, vi, w * values[vj])
            np.add.at(acc, vj, w * values[vi])
            np.add.at(wsum, vi, w)
            np.add.at(wsum, vj, w)
        out = values.copy()
        has = wsum > 1.0e-12
        out[has] = acc[has] / wsum[has]
        out[~mask] = 0.0
        return out

    def _diffuse_field_directional(
        self,
        grid,
        values,
        mask,
        conductance,
        axis,
        *,
        passes: int = 1,
        alpha: float = 0.2,
        seed_retention: float = 0.0,
        direction_strength: float = 1.25,
        direction_power: float = 2.0,
    ):
        mask = np.asarray(mask, dtype=bool)
        out = np.where(mask, np.asarray(values, dtype=np.float64), 0.0)
        seed_values = out.copy()
        conductance = np.asarray(conductance, dtype=np.float64)
        if conductance.shape != (grid.n,):
            conductance = np.ones(grid.n, dtype=np.float64)
        conductance = np.where(mask, np.clip(conductance, 0.0, 1.5), 0.0)
        conductance = np.maximum(conductance, 0.04 * mask.astype(np.float64))
        for _ in range(max(int(passes), 0)):
            mean = self._directional_neighbor_mean_masked(
                grid,
                out,
                mask,
                conductance,
                axis,
                direction_strength=direction_strength,
                direction_power=direction_power,
            )
            local_alpha = alpha * np.clip(0.45 + 0.55 * conductance, 0.0, 1.0)
            out = np.where(
                mask,
                (1.0 - local_alpha) * out + local_alpha * mean,
                0.0,
            )
            if seed_retention > 0.0:
                out = np.where(
                    mask,
                    (1.0 - seed_retention) * out
                    + seed_retention * seed_values,
                    0.0,
                )
        return out

    def _downwind_spread_field_directional(
        self,
        grid,
        seed,
        allowed,
        conductance,
        axis,
        *,
        passes: int = 1,
        alpha: float = 0.25,
        decay: float = 0.88,
        seed_retention: float = 0.08,
        direction_power: float = 2.0,
    ):
        """One-way reduced advection used for downwind source support."""
        allowed = np.asarray(allowed, dtype=bool)
        seed = np.where(allowed, np.asarray(seed, dtype=np.float64), 0.0)
        out = seed.copy()
        conductance = np.asarray(conductance, dtype=np.float64)
        if conductance.shape != (grid.n,):
            conductance = np.ones(grid.n, dtype=np.float64)
        conductance = np.where(allowed, np.clip(conductance, 0.0, 1.5), 0.0)
        axis = np.asarray(axis, dtype=np.float64)
        if axis.shape != (grid.n, 3):
            east, _ = self._tangent_basis(grid)
            axis = east
        axis = self._project_tangent(grid, axis)
        axis_norm = np.linalg.norm(axis, axis=1, keepdims=True)
        east, _ = self._tangent_basis(grid)
        axis = np.where(axis_norm > 1.0e-9, axis / np.maximum(axis_norm, 1.0e-9), east)

        i, j = self._edges[:, 0], self._edges[:, 1]
        edge = grid.xyz[j] - grid.xyz[i]
        dir_i = edge - np.sum(edge * grid.xyz[i], axis=1, keepdims=True) * grid.xyz[i]
        norm_i = np.linalg.norm(dir_i, axis=1, keepdims=True)
        dir_i = np.where(norm_i > 1.0e-12, dir_i / np.maximum(norm_i, 1.0e-12), 0.0)
        dir_j = -edge - np.sum((-edge) * grid.xyz[j], axis=1, keepdims=True) * grid.xyz[j]
        norm_j = np.linalg.norm(dir_j, axis=1, keepdims=True)
        dir_j = np.where(norm_j > 1.0e-12, dir_j / np.maximum(norm_j, 1.0e-12), 0.0)

        weight_ij = (
            np.sqrt(np.clip(conductance[i], 0.0, None) * np.clip(conductance[j], 0.0, None))
            * np.clip(np.sum(dir_i * axis[i], axis=1), 0.0, 1.0) ** direction_power
        )
        weight_ji = (
            np.sqrt(np.clip(conductance[i], 0.0, None) * np.clip(conductance[j], 0.0, None))
            * np.clip(np.sum(dir_j * axis[j], axis=1), 0.0, 1.0) ** direction_power
        )
        valid_ij = allowed[i] & allowed[j] & (weight_ij > 1.0e-9)
        valid_ji = allowed[i] & allowed[j] & (weight_ji > 1.0e-9)
        for _ in range(max(int(passes), 0)):
            acc = np.zeros(grid.n, dtype=np.float64)
            wsum = np.zeros(grid.n, dtype=np.float64)
            if valid_ij.any():
                vi = i[valid_ij]
                vj = j[valid_ij]
                w = weight_ij[valid_ij]
                np.add.at(acc, vj, w * out[vi])
                np.add.at(wsum, vj, w)
            if valid_ji.any():
                vi = i[valid_ji]
                vj = j[valid_ji]
                w = weight_ji[valid_ji]
                np.add.at(acc, vi, w * out[vj])
                np.add.at(wsum, vi, w)
            mean = np.zeros(grid.n, dtype=np.float64)
            has = wsum > 1.0e-12
            mean[has] = acc[has] / wsum[has]
            out = np.where(
                allowed,
                decay * ((1.0 - alpha) * out + alpha * mean)
                + seed_retention * seed,
                0.0,
            )
        return out

    def _expand_mask_within(self, grid, seed, allowed, passes: int = 1) -> np.ndarray:
        self._edge_cache(grid)
        out = np.asarray(seed, dtype=bool) & np.asarray(allowed, dtype=bool)
        allowed = np.asarray(allowed, dtype=bool)
        if not out.any() or passes <= 0:
            return out
        i = self._edges[:, 0]
        j = self._edges[:, 1]
        for _ in range(int(passes)):
            grow = out.copy()
            ij = out[i] & allowed[j]
            ji = out[j] & allowed[i]
            grow[j[ij]] = True
            grow[i[ji]] = True
            if np.array_equal(grow, out):
                break
            out = grow
        return out

    def _smooth_vector_field_masked(self, grid, values, mask, passes: int = 1,
                                    alpha: float = 0.15):
        out = np.asarray(values, dtype=np.float64).copy()
        mask = np.asarray(mask, dtype=bool)
        out[~mask] = 0.0
        for _ in range(passes):
            smoothed = np.column_stack([
                self._masked_neighbor_mean(grid, out[:, k], mask) for k in range(3)
            ])
            out = np.where(mask[:, None], (1.0 - alpha) * out + alpha * smoothed, 0.0)
            out = self._project_tangent(grid, out)
            out[~mask] = 0.0
        return out

    def _ocean_components(self, grid, ocean):
        return self._connected_components(grid, ocean)

    def _connected_components(self, grid, mask):
        mask = np.asarray(mask, dtype=bool)
        component_id = np.full(grid.n, -1, dtype=int)
        component_area: list[float] = []
        for start in np.where(mask)[0]:
            if component_id[start] >= 0:
                continue
            cid = len(component_area)
            stack = [int(start)]
            component_id[start] = cid
            acc = 0.0
            while stack:
                c = stack.pop()
                acc += float(grid.cell_area[c])
                for nb in grid.neighbors[c]:
                    nb = int(nb)
                    if mask[nb] and component_id[nb] < 0:
                        component_id[nb] = cid
                        stack.append(nb)
            component_area.append(acc)
        return component_id, np.asarray(component_area, dtype=np.float64)

    def _graph_hop_distance(self, grid, sources, allowed):
        sources = np.asarray(sources, dtype=bool)
        allowed = np.asarray(allowed, dtype=bool)
        dist = np.full(grid.n, -1, dtype=int)
        starts = np.where(sources & allowed)[0]
        dist[starts] = 0
        queue = [int(x) for x in starts]
        head = 0
        while head < len(queue):
            c = queue[head]
            head += 1
            nd = dist[c] + 1
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if allowed[nb] and dist[nb] < 0:
                    dist[nb] = nd
                    queue.append(nb)
        return dist

    def _component_objects(self, grid, component_id, component_area, mask,
                           adjacent_id=None, adjacent_mask=None,
                           boundary_mask=None, prefix="component"):
        objects: list[dict] = []
        total_area = max(float(grid.cell_area.sum()), 1e-12)
        ids = [int(x) for x in np.unique(component_id[mask]) if int(x) >= 0]
        for cid in ids:
            cells = mask & (component_id == cid)
            if not cells.any():
                continue
            weights = grid.cell_area[cells]
            centroid = np.average(grid.xyz[cells], axis=0, weights=weights)
            cn = float(np.linalg.norm(centroid))
            if cn > 1e-12:
                centroid = centroid / cn
            lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
            lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
            adjacent: set[int] = set()
            if adjacent_id is not None and adjacent_mask is not None:
                edge_cells = cells if boundary_mask is None else cells & boundary_mask
                for c in np.where(edge_cells)[0]:
                    nbs = grid.neighbors[int(c)]
                    for aid in adjacent_id[nbs[adjacent_mask[nbs]]].astype(int):
                        if aid >= 0:
                            adjacent.add(int(aid))
            objects.append({
                "id": cid,
                "type": prefix,
                "cell_count": int(cells.sum()),
                "area_fraction": float(component_area[cid] / total_area),
                "centroid_lat": lat,
                "centroid_lon": lon,
                "lat_min": float(np.min(grid.lat[cells])),
                "lat_max": float(np.max(grid.lat[cells])),
                "adjacent_ids": sorted(adjacent),
            })
        return objects

    def _real_earth_major_ocean_basins(self, grid, ocean):
        """Return semantic major-ocean ids for Earth replay diagnostics.

        Connected-component basins merge the modern ocean into one global water
        body.  That is topologically true but not useful for real-Earth pressure
        replay, where North Atlantic, North Pacific, Indian, Arctic, and
        Southern basin-scale forcing must be separable before R2a attribution.
        """
        ocean = np.asarray(ocean, dtype=bool)
        lon = ((np.asarray(grid.lon, dtype=np.float64) + 180.0) % 360.0) - 180.0
        lat = np.asarray(grid.lat, dtype=np.float64)
        basin_id = np.full(grid.n, -1, dtype=int)

        southern = ocean & (lat <= -55.0)
        arctic = ocean & (lat >= 66.0)
        active = ocean & ~southern & ~arctic

        indian = (
            active
            & (lon >= 20.0)
            & (lon < 150.0)
            & (lat > -55.0)
            & (lat < 32.0)
        )
        west_atlantic_limit = np.where(lat >= 0.0, -105.0, -70.0)
        atlantic = (
            active
            & ~indian
            & (lon >= west_atlantic_limit)
            & (lon < 25.0)
        )
        atlantic |= (
            active
            & (lat >= 30.0)
            & (lat <= 50.0)
            & (lon >= 20.0)
            & (lon <= 45.0)
        )
        pacific = active & ~indian & ~atlantic

        basin_id[atlantic] = 0
        basin_id[pacific] = 1
        basin_id[indian] = 2
        basin_id[arctic] = 3
        basin_id[southern] = 4
        basin_id[ocean & (basin_id < 0)] = 1

        names = {
            0: "Atlantic Ocean",
            1: "Pacific Ocean",
            2: "Indian Ocean",
            3: "Arctic Ocean",
            4: "Southern Ocean",
        }
        basin_area = np.zeros(max(names) + 1, dtype=np.float64)
        for bid in names:
            basin_area[bid] = float(np.sum(grid.cell_area[ocean & (basin_id == bid)]))
        return basin_id, basin_area, names

    @staticmethod
    def _dominant_nonnegative_id(values: np.ndarray) -> int | None:
        ids = np.asarray(values, dtype=int)
        ids = ids[ids >= 0]
        if ids.size == 0:
            return None
        counts = np.bincount(ids)
        return int(np.argmax(counts))

    def _pressure_center_diagnostics(self, grid, pressure_proxy, land, ocean,
                                     geography_fields):
        pressure_proxy = np.asarray(pressure_proxy, dtype=np.float64)
        land = np.asarray(land, dtype=bool)
        ocean = np.asarray(ocean, dtype=bool)
        center_support = np.zeros((4, grid.n), dtype=np.float64)
        center_id = np.zeros((4, grid.n), dtype=np.float64)
        stationary_support = np.zeros((4, grid.n), dtype=np.float64)
        if pressure_proxy.shape != (4, grid.n):
            return {
                "center_support": center_support,
                "center_id": center_id,
                "stationary_wave_support": stationary_support,
                "objects": [],
            }

        geography_fields = geography_fields or {}
        basin_id = np.asarray(
            geography_fields.get("ocean.basin_id", np.full(grid.n, -1.0)),
            dtype=np.float64,
        )
        if basin_id.shape != (grid.n,):
            basin_id = np.full(grid.n, -1.0, dtype=np.float64)
        continent_id = np.asarray(
            geography_fields.get("climate.continent_id", np.full(grid.n, -1.0)),
            dtype=np.float64,
        )
        if continent_id.shape != (grid.n,):
            continent_id = np.full(grid.n, -1.0, dtype=np.float64)
        coast_distance = np.asarray(
            geography_fields.get("climate.coast_distance", np.zeros(grid.n)),
            dtype=np.float64,
        )
        if coast_distance.shape != (grid.n,):
            coast_distance = np.zeros(grid.n, dtype=np.float64)
        coast_strength = np.asarray(
            geography_fields.get("climate.coast_strength", np.zeros(grid.n)),
            dtype=np.float64,
        )
        if coast_strength.shape != (grid.n,):
            coast_strength = np.zeros(grid.n, dtype=np.float64)
        barrier = np.asarray(
            geography_fields.get("terrain.barrier_index", np.zeros(grid.n)),
            dtype=np.float64,
        )
        if barrier.shape != (grid.n,):
            barrier = np.zeros(grid.n, dtype=np.float64)

        lat_gate = np.clip((88.0 - np.abs(grid.lat)) / 10.0, 0.0, 1.0)
        geography_anchor = np.clip(
            0.40
            + 0.30 * np.clip(coast_strength, 0.0, 1.0)
            + 0.22 * np.clip(barrier, 0.0, 1.0)
            + 0.18 * ocean.astype(float) * np.clip(coast_distance, 0.0, 1.0),
            0.0,
            1.25,
        )
        total_area = max(float(np.sum(grid.cell_area)), 1.0e-12)
        objects: list[dict] = []
        next_object_id = 1
        min_cells = max(3, int(np.ceil(grid.n * 3.5e-4)))
        min_area_fraction = 5.0e-4

        for season_index, season_name in enumerate(self.SEASON_NAMES):
            zonal = self._latitude_band_anomaly(grid, pressure_proxy[season_index])
            wave = self._smooth_field(grid, zonal, passes=2, alpha=0.10)
            finite = np.isfinite(wave)
            if not finite.any():
                continue
            abs_wave = np.abs(wave)
            lo, hi = np.percentile(abs_wave[finite], [50, 92])
            span = max(float(hi - lo), 1.0e-9)
            support = np.clip((abs_wave - lo) / span, 0.0, 1.35) * lat_gate
            support = np.where(finite, support, 0.0)
            center_support[season_index] = support
            stationary = np.clip(
                support * (0.55 + 0.45 * geography_anchor),
                0.0,
                1.35,
            )
            stationary_support[season_index] = self._smooth_field(
                grid, stationary, passes=1, alpha=0.08)

            active_threshold = max(0.24, float(np.percentile(support[finite], 76)))
            signed_masks = (
                ("pressure_low", wave < 0.0),
                ("pressure_high", wave > 0.0),
            )
            rows: list[tuple[float, int, np.ndarray, str]] = []
            for kind, sign_mask in signed_masks:
                active = finite & sign_mask & (support >= active_threshold)
                component_id, component_area = self._connected_components(grid, active)
                for cid in [int(x) for x in np.unique(component_id[active]) if int(x) >= 0]:
                    cells = active & (component_id == cid)
                    cell_count = int(np.count_nonzero(cells))
                    area = float(component_area[cid])
                    area_fraction = area / total_area
                    if cell_count < min_cells or area_fraction < min_area_fraction:
                        continue
                    weights = grid.cell_area[cells]
                    mean_support = float(np.average(support[cells], weights=weights))
                    score = area_fraction * (0.45 + mean_support)
                    rows.append((score, cid, cells, kind))

            rows.sort(key=lambda row: row[0], reverse=True)
            per_kind_counts = {"pressure_low": 0, "pressure_high": 0}
            for _, cid, cells, kind in rows:
                if per_kind_counts[kind] >= 12:
                    continue
                per_kind_counts[kind] += 1
                object_id = next_object_id
                next_object_id += 1
                weights = grid.cell_area[cells]
                support_weights = weights * np.maximum(support[cells], 1.0e-6)
                centroid = np.average(grid.xyz[cells], axis=0, weights=support_weights)
                cn = float(np.linalg.norm(centroid))
                if cn > 1.0e-12:
                    centroid = centroid / cn
                lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
                lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
                ocean_fraction = float(
                    np.average(ocean[cells].astype(float), weights=weights))
                land_fraction = float(
                    np.average(land[cells].astype(float), weights=weights))
                if ocean_fraction >= 0.67:
                    domain = "ocean"
                elif land_fraction >= 0.67:
                    domain = "land"
                else:
                    domain = "mixed"
                dominant_basin = self._dominant_nonnegative_id(basin_id[cells])
                dominant_continent = self._dominant_nonnegative_id(continent_id[cells])
                center_id[season_index, cells] = float(object_id)
                objects.append({
                    "id": int(object_id),
                    "component_id": int(cid),
                    "type": "atmosphere.pressure_center",
                    "kind": kind,
                    "season": season_name,
                    "season_index": int(season_index),
                    "domain": domain,
                    "cell_count": int(np.count_nonzero(cells)),
                    "area_fraction": float(np.sum(weights) / total_area),
                    "centroid_lat": lat,
                    "centroid_lon": lon,
                    "lat_min": float(np.min(grid.lat[cells])),
                    "lat_max": float(np.max(grid.lat[cells])),
                    "lon_min": float(np.min(grid.lon[cells])),
                    "lon_max": float(np.max(grid.lon[cells])),
                    "land_fraction": land_fraction,
                    "ocean_fraction": ocean_fraction,
                    "dominant_ocean_basin_id": (
                        int(dominant_basin) if dominant_basin is not None else -1),
                    "dominant_continent_id": (
                        int(dominant_continent)
                        if dominant_continent is not None else -1),
                    "mean_pressure_zonal_anomaly": float(
                        np.average(wave[cells], weights=weights)),
                    "min_pressure_zonal_anomaly": float(np.min(wave[cells])),
                    "max_pressure_zonal_anomaly": float(np.max(wave[cells])),
                    "mean_center_support": float(
                        np.average(support[cells], weights=weights)),
                    "p90_center_support": float(
                        np.percentile(support[cells], 90)),
                    "mean_stationary_wave_support": float(
                        np.average(stationary_support[season_index, cells],
                                   weights=weights)),
                    "mean_coast_distance": float(
                        np.average(np.clip(coast_distance[cells], 0.0, 1.0),
                                   weights=weights)),
                    "mean_coast_strength": float(
                        np.average(np.clip(coast_strength[cells], 0.0, 1.0),
                                   weights=weights)),
                    "mean_barrier_index": float(
                        np.average(np.clip(barrier[cells], 0.0, 1.0),
                                   weights=weights)),
                })

        return {
            "center_support": np.clip(center_support, 0.0, 1.35),
            "center_id": center_id,
            "stationary_wave_support": np.clip(stationary_support, 0.0, 1.35),
            "objects": objects,
        }

    def _pressure_object_support(
        self,
        grid,
        score,
        domain,
        labels,
        *,
        quantile: float,
        max_components: int,
        min_cells: int,
        smooth_passes: int,
        smooth_alpha: float,
        multi_large_labels: bool = False,
    ):
        """Extract broad pressure-source objects from a scalar support field."""
        score = np.asarray(score, dtype=np.float64)
        domain = np.asarray(domain, dtype=bool)
        labels = np.asarray(labels, dtype=np.int64)
        if score.shape != (grid.n,) or labels.shape != (grid.n,) or not domain.any():
            return np.zeros(grid.n, dtype=np.float64)

        support = np.zeros(grid.n, dtype=np.float64)
        total_area = max(float(np.sum(grid.cell_area)), 1.0e-12)
        ids = [int(x) for x in np.unique(labels[domain]) if int(x) >= 0]
        if not ids:
            ids = [0]
            labels = np.where(domain, 0, -1)

        for label_id in ids:
            cells = (
                domain
                & (labels == int(label_id))
                & np.isfinite(score)
                & (score > 0.0)
            )
            if int(np.count_nonzero(cells)) < min_cells:
                continue
            threshold = max(float(np.percentile(score[cells], quantile)), 1.0e-8)
            active = cells & (score >= threshold)
            component_id, component_area = self._connected_components(grid, active)
            rows: list[tuple[float, int, np.ndarray]] = []
            for cid in [int(x) for x in np.unique(component_id[active]) if int(x) >= 0]:
                component_cells = active & (component_id == cid)
                cell_count = int(np.count_nonzero(component_cells))
                area_fraction = float(component_area[cid] / total_area)
                if cell_count < min_cells or area_fraction < 2.0e-4:
                    continue
                rows.append((
                    float(np.sum(grid.cell_area[component_cells] * score[component_cells])),
                    cid,
                    component_cells,
                ))

            rows.sort(key=lambda row: row[0], reverse=True)
            label_area = float(np.sum(grid.cell_area[cells]))
            mean_abs_lat = float(np.average(
                np.abs(grid.lat[cells]), weights=grid.cell_area[cells]))
            broad_polar_basin = mean_abs_lat > 43.0 and label_area / total_area > 0.04
            broad_large_label = multi_large_labels and label_area / total_area > 0.035
            keep = int(
                max_components
                if (broad_polar_basin or broad_large_label)
                else min(max_components, 1))
            for _, _, component_cells in rows[:keep]:
                support[component_cells] = np.maximum(
                    support[component_cells], score[component_cells])

        support = self._smooth_field_masked(
            grid, support, domain, passes=smooth_passes, alpha=smooth_alpha)
        active_values = support[domain & np.isfinite(support)]
        if active_values.size:
            scale = float(np.percentile(active_values, 95))
            if scale > 1.0e-9:
                support = support / scale
        return np.where(domain, np.clip(support, 0.0, 1.25), 0.0)

    def _same_latitude_ocean_sst_anomaly(self, grid, seasonal_sst, ocean):
        seasonal_sst = np.asarray(seasonal_sst, dtype=np.float64)
        ocean = np.asarray(ocean, dtype=bool)
        out = np.zeros_like(seasonal_sst, dtype=np.float64)
        if seasonal_sst.shape != (4, grid.n) or not ocean.any():
            return out
        band_id = np.floor((grid.lat + 90.0) / 5.0).astype(int)
        for season in range(4):
            values = seasonal_sst[season]
            for band in np.unique(band_id):
                cells = (
                    ocean
                    & (band_id == int(band))
                    & np.isfinite(values)
                )
                if not cells.any():
                    continue
                mean = float(np.average(values[cells], weights=grid.cell_area[cells]))
                out[season, cells] = values[cells] - mean
        return out

    def _sst_front_support(self, grid, seasonal_sst, ocean):
        seasonal_sst = np.asarray(seasonal_sst, dtype=np.float64)
        ocean = np.asarray(ocean, dtype=bool)
        out = np.zeros_like(seasonal_sst, dtype=np.float64)
        if seasonal_sst.shape != (4, grid.n) or not ocean.any():
            return out
        for season in range(4):
            grad = self._graph_gradient_vectors(grid, seasonal_sst[season])
            mag = np.linalg.norm(grad, axis=1)
            cells = ocean & np.isfinite(mag)
            if not cells.any():
                continue
            lo, hi = np.percentile(mag[cells], [45, 95])
            out[season] = np.where(
                ocean,
                np.clip((mag - lo) / max(float(hi - lo), 1.0e-9), 0.0, 1.5),
                0.0,
            )
        return out

    def _southern_ocean_pressure_wave_gate(
        self,
        grid,
        basin_label,
        ocean,
        sst_front_season,
        shelf_index,
        same_lat_sst_anomaly,
        *,
        season_index: int,
    ):
        """Reduced wavenumber/front gate for Southern Ocean pressure sources."""
        basin_label = np.asarray(basin_label, dtype=np.int64)
        ocean = np.asarray(ocean, dtype=bool)
        sst_front_season = np.asarray(sst_front_season, dtype=np.float64)
        shelf_index = np.asarray(shelf_index, dtype=np.float64)
        same_lat_sst_anomaly = np.asarray(same_lat_sst_anomaly, dtype=np.float64)
        gate = np.ones(grid.n, dtype=np.float64)
        southern = ocean & (basin_label == 4) & (grid.lat < -45.0)
        if not southern.any() or sst_front_season.shape != (grid.n,):
            return gate, southern
        if shelf_index.shape != (grid.n,):
            shelf_index = np.zeros(grid.n, dtype=np.float64)
        if same_lat_sst_anomaly.shape != (grid.n,):
            same_lat_sst_anomaly = np.zeros(grid.n, dtype=np.float64)

        front_anom = self._latitude_band_anomaly(grid, sst_front_season)
        lo, hi = np.percentile(front_anom[southern], [35, 90])
        front_phase = np.clip(
            (front_anom - float(lo)) / max(float(hi - lo), 1.0e-9),
            0.0,
            1.0,
        )
        shelf_anom = self._latitude_band_anomaly(
            grid, np.where(ocean, np.clip(shelf_index, 0.0, 1.0), 0.0))
        slo, shi = np.percentile(shelf_anom[southern], [35, 90])
        shelf_phase = np.clip(
            (shelf_anom - float(slo)) / max(float(shi - slo), 1.0e-9),
            0.0,
            1.0,
        )
        warm_phase = np.clip(same_lat_sst_anomaly / 1.4, 0.0, 1.0)
        # A weak wavenumber term represents the reduced stationary-wave
        # tendency, but Southern Ocean source placement is primarily gated by
        # current front/shelf/SST support rather than fixed longitude phase.
        lon_wave = 0.5 + 0.5 * np.cos(
            np.radians(3.0 * grid.lon + 40.0 + 25.0 * float(season_index)))
        raw = (
            0.48 * front_phase
            + 0.32 * shelf_phase
            + 0.14 * warm_phase
            + 0.06 * lon_wave
        )
        wave_core = np.clip((raw - 0.11) / 0.23, 0.0, 1.0) ** 1.15
        gate[southern] = 0.05 + 0.95 * wave_core[southern]
        return gate, southern

    def _m2_pressure_genesis(
        self,
        world,
        pressure_proxy,
        seasonal_T,
        seasonal_sst,
        ocean,
        elev,
        geography_fields,
    ):
        """Object-based R2a pressure-source refinement.

        This pass keeps the existing thermal pressure proxy as the base, then
        adds only geography-supported pressure-center objects: winter subpolar
        ocean lows anchored to basin/front support, weak basin subtropical
        highs, and small terrain/continent-supported stationary refinements.
        It deliberately does not translate the adjusted pressure into wind;
        R2b remains a separate owner.
        """
        grid = world.grid
        pressure = np.asarray(pressure_proxy, dtype=np.float64)
        zeros = np.zeros((4, grid.n), dtype=np.float64)
        empty = {
            "pressure": pressure.copy() if pressure.shape == (4, grid.n) else zeros.copy(),
            "source": zeros.copy(),
            "wave_transfer": zeros.copy(),
            "ocean_low_support": zeros.copy(),
            "ocean_high_support": zeros.copy(),
            "land_support": zeros.copy(),
            "terrain_wave_support": zeros.copy(),
        }
        if pressure.shape != (4, grid.n):
            return empty
        seasonal_T = np.asarray(seasonal_T, dtype=np.float64)
        if seasonal_T.shape != (4, grid.n):
            return empty
        seasonal_sst = np.asarray(seasonal_sst, dtype=np.float64)
        if seasonal_sst.shape != (4, grid.n):
            seasonal_sst = seasonal_T

        ocean = np.asarray(ocean, dtype=bool)
        land = ~ocean
        land_extent = np.clip((world.land_fraction() - 0.035) / 0.22, 0.0, 1.0)
        if land_extent <= 1.0e-6:
            return empty

        geography_fields = geography_fields or {}
        coast_distance = np.asarray(
            geography_fields.get("climate.coast_distance", np.zeros(grid.n)),
            dtype=np.float64,
        )
        if coast_distance.shape != (grid.n,):
            coast_distance = np.zeros(grid.n, dtype=np.float64)
        continent_interiority = np.asarray(
            geography_fields.get("climate.continent_interiority", np.zeros(grid.n)),
            dtype=np.float64,
        )
        if continent_interiority.shape != (grid.n,):
            continent_interiority = np.zeros(grid.n, dtype=np.float64)
        continent_id = np.asarray(
            geography_fields.get("climate.continent_id", np.full(grid.n, -1.0)),
            dtype=np.float64,
        )
        if continent_id.shape != (grid.n,):
            continent_id = np.full(grid.n, -1.0, dtype=np.float64)
        continent_label = continent_id.astype(np.int64)
        basin_id = np.asarray(
            geography_fields.get("ocean.basin_id", np.full(grid.n, -1.0)),
            dtype=np.float64,
        )
        if basin_id.shape != (grid.n,):
            basin_id = np.full(grid.n, -1.0, dtype=np.float64)
        basin_label = basin_id.astype(np.int64)
        barrier = np.asarray(
            geography_fields.get("terrain.barrier_index", np.zeros(grid.n)),
            dtype=np.float64,
        )
        if barrier.shape != (grid.n,):
            barrier = np.zeros(grid.n, dtype=np.float64)
        coast_strength = np.asarray(
            geography_fields.get("climate.coast_strength", np.zeros(grid.n)),
            dtype=np.float64,
        )
        if coast_strength.shape != (grid.n,):
            coast_strength = np.zeros(grid.n, dtype=np.float64)
        shelf_index = np.asarray(
            geography_fields.get("ocean.shelf_index", np.zeros(grid.n)),
            dtype=np.float64,
        )
        if shelf_index.shape != (grid.n,):
            shelf_index = np.zeros(grid.n, dtype=np.float64)
        coast_orientation = np.asarray(
            geography_fields.get(
                "climate.coast_orientation",
                np.zeros((grid.n, 3), dtype=np.float64),
            ),
            dtype=np.float64,
        )
        if coast_orientation.shape != (grid.n, 3):
            coast_orientation = np.zeros((grid.n, 3), dtype=np.float64)

        basin_scale = np.zeros(grid.n, dtype=np.float64)
        basin_ids = [int(x) for x in np.unique(basin_label[ocean]) if int(x) >= 0]
        if basin_ids:
            basin_area = {
                bid: float(np.sum(grid.cell_area[ocean & (basin_label == bid)]))
                for bid in basin_ids
            }
            max_area = max(max(basin_area.values()), 1.0e-12)
            for bid, area in basin_area.items():
                basin_scale[basin_label == bid] = np.clip(
                    np.sqrt(area / max_area), 0.35, 1.0)

        lat_abs = np.abs(grid.lat)
        open_ocean = np.where(ocean, np.clip(coast_distance, 0.0, 1.0) ** 0.68, 0.0)
        subpolar = np.exp(-((lat_abs - 53.0) / 13.0) ** 2)
        subtropical = np.exp(-((lat_abs - 31.0) / 14.0) ** 2)
        midlatitude = np.exp(-((lat_abs - 48.0) / 16.0) ** 2)
        season_ratio = np.array([-1.0, 0.0, 1.0, 0.0], dtype=np.float64)

        same_lat_sst = self._same_latitude_ocean_sst_anomaly(grid, seasonal_sst, ocean)
        sst_front = self._sst_front_support(grid, seasonal_sst, ocean)
        sst_cell_anom = seasonal_sst - seasonal_sst.mean(axis=0, keepdims=True)
        temp_anom = seasonal_T - seasonal_T.mean(axis=0, keepdims=True)

        adjusted = pressure.copy()
        base_pressure = pressure.copy()
        source_field = np.zeros((4, grid.n), dtype=np.float64)
        low_support_field = np.zeros((4, grid.n), dtype=np.float64)
        high_support_field = np.zeros((4, grid.n), dtype=np.float64)
        land_support_field = np.zeros((4, grid.n), dtype=np.float64)
        terrain_support_field = np.zeros((4, grid.n), dtype=np.float64)
        southern_wave_gate_field = np.zeros((4, grid.n), dtype=np.float64)
        subpolar_high_transfer_field = np.zeros((4, grid.n), dtype=np.float64)
        mam_arctic_freeze_high_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        jja_north_pacific_high_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        mam_canadian_arctic_high_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        mam_central_arctic_high_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        djf_atlantic_gateway_low_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        djf_north_pacific_low_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        jja_north_pacific_central_high_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        jja_eurasian_thermal_low_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        jja_west_pacific_marginal_low_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        jja_west_asia_thermal_low_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        jja_east_china_sea_low_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        djf_north_america_winter_high_relief_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        mam_north_america_plains_high_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        mam_north_america_land_high_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        son_autumn_land_relief_transfer_field = np.zeros((4, grid.n), dtype=np.float64)
        son_north_atlantic_autumn_low_transfer_field = np.zeros(
            (4, grid.n), dtype=np.float64)
        land_shoulder_phase_transfer_field = np.zeros((4, grid.n), dtype=np.float64)
        min_cells = max(6, int(np.ceil(grid.n * 7.5e-4)))
        source_label = basin_label.copy()
        southern_ocean = ocean & (basin_label == 4) & (grid.lat < -45.0)
        if southern_ocean.any():
            sector = np.floor(((grid.lon + 180.0) % 360.0) / 60.0).astype(np.int64)
            source_label[southern_ocean] = 40 + sector[southern_ocean]
        east_axis, _ = self._tangent_basis(grid)
        downwind_allowed = ocean | land
        downwind_conductance = np.where(
            ocean,
            0.22 + 0.78 * open_ocean,
            0.12 + 0.42 * np.clip(coast_strength, 0.0, 1.0),
        )
        for season, ratio in enumerate(season_ratio):
            winter_hemi = np.clip(-np.sign(grid.lat) * ratio, 0.0, 1.0)
            summer_hemi = np.clip(np.sign(grid.lat) * ratio, 0.0, 1.0)
            shoulder_hemi = 1.0 - np.clip(winter_hemi + summer_hemi, 0.0, 1.0)
            southern_gate, southern_mask = self._southern_ocean_pressure_wave_gate(
                grid,
                basin_label,
                ocean,
                sst_front[season],
                shelf_index,
                same_lat_sst[season],
                season_index=season,
            )
            southern_wave_gate_field[season] = np.where(southern_mask, southern_gate, 0.0)
            front_support = np.clip(
                0.28
                + 0.55 * np.clip(sst_front[season], 0.0, 1.4)
                + 0.35 * np.clip(same_lat_sst[season] / 2.0, 0.0, 1.4),
                0.0,
                1.7,
            )
            land_thermal = (
                -np.clip(temp_anom[season] / 9.0, -1.0, 1.0)
                * land.astype(float)
                * (0.35 + 0.65 * np.clip(continent_interiority, 0.0, 1.0))
            )
            land_thermal = self._smooth_field_masked(
                grid, land_thermal, land, passes=3, alpha=0.10)
            cold_continent_seed = (
                np.clip(land_thermal, 0.0, 1.35)
                * land.astype(float)
                * midlatitude
                * winter_hemi
                * np.where(grid.lat >= 0.0, 1.0, 0.42)
                * (0.35 + 0.65 * np.clip(continent_interiority, 0.0, 1.0))
                * (1.0 - np.clip((lat_abs - 63.0) / 12.0, 0.0, 1.0))
            )
            downwind_cold_land = self._downwind_spread_field_directional(
                grid,
                cold_continent_seed,
                downwind_allowed,
                downwind_conductance,
                east_axis,
                passes=9,
                alpha=0.32,
                decay=0.88,
                seed_retention=0.10,
                direction_power=2.4,
            )
            lee_low_support = np.where(ocean, downwind_cold_land, 0.0)
            active_lee = lee_low_support[ocean & np.isfinite(lee_low_support)]
            if active_lee.size:
                lee_scale = float(np.percentile(active_lee, 95))
                if lee_scale > 1.0e-9:
                    lee_low_support = np.where(
                        ocean, np.clip(lee_low_support / lee_scale, 0.0, 1.35), 0.0)
            subpolar_coastal_front = (
                ocean.astype(float)
                * np.where(grid.lat >= 0.0, 1.0, 0.0)
                * subpolar
                * winter_hemi
                * np.clip((lat_abs - 38.0) / 10.0, 0.0, 1.0)
                * (1.0 - np.clip((lat_abs - 68.0) / 8.0, 0.0, 1.0))
                * np.clip(1.0 - open_ocean, 0.0, 1.0)
                * np.clip(sst_front[season] / 1.2, 0.0, 1.25)
                * (0.45 + 0.55 * basin_scale)
            )
            poleward_subpolar_front_low = (
                ocean.astype(float)
                * np.where(grid.lat >= 0.0, 1.0, 0.0)
                * winter_hemi
                * np.exp(-((lat_abs - 60.0) / 8.0) ** 2)
                * np.clip((lat_abs - 47.0) / 8.0, 0.0, 1.0)
                * (1.0 - np.clip((lat_abs - 67.0) / 4.0, 0.0, 1.0))
                * front_support
                * (0.60 + 0.40 * np.clip(1.0 - open_ocean, 0.0, 1.0))
                * (0.50 + 0.50 * basin_scale)
            )
            shoulder_warm_ocean_low = (
                ocean.astype(float)
                * np.where(grid.lat >= 0.0, 1.0, 0.0)
                * (1.0 - np.clip(winter_hemi + summer_hemi, 0.0, 1.0))
                * subpolar
                * np.clip(sst_cell_anom[season] / 1.15, 0.0, 1.25)
                * np.clip(sst_front[season] / 1.2, 0.0, 1.25)
                * (0.35 + 0.65 * open_ocean)
                * (0.45 + 0.55 * basin_scale)
            )
            low_score = (
                open_ocean
                * basin_scale
                * subpolar
                * winter_hemi
                * front_support
                * southern_gate
                * (0.88 + 0.26 * lee_low_support)
                + 0.10
                * lee_low_support
                * open_ocean
                * basin_scale
                * subpolar
                * winter_hemi
                * southern_gate
            )
            low_support = self._pressure_object_support(
                grid,
                low_score,
                ocean,
                source_label,
                quantile=72.0,
                max_components=3,
                min_cells=min_cells,
                smooth_passes=6,
                smooth_alpha=0.16,
            )
            shoulder_support = self._pressure_object_support(
                grid,
                shoulder_warm_ocean_low * front_support,
                ocean,
                basin_label,
                quantile=82.0,
                max_components=2,
                min_cells=min_cells,
                smooth_passes=4,
                smooth_alpha=0.14,
            )
            shoulder_amplitude = (
                0.22
                * np.clip(open_ocean / 0.45, 0.0, 1.0) ** 1.75
                * np.clip(sst_cell_anom[season] / 1.15, 0.0, 1.0)
                * np.clip(basin_scale, 0.35, 1.0) ** 3.0
            )
            low_support = np.maximum(
                low_support,
                np.clip(shoulder_support * shoulder_amplitude, 0.0, 0.42),
            )
            low_support = np.maximum(
                low_support,
                np.clip(0.30 * subpolar_coastal_front * front_support, 0.0, 0.85),
            )
            low_support = np.maximum(
                low_support,
                np.clip(1.05 * poleward_subpolar_front_low, 0.0, 0.82),
            )
            atlantic_gateway_sources = (
                ocean
                & (basin_label == 0)
                & (grid.lat >= 55.0)
                & (grid.lat <= 70.0)
                & (low_support > 0.45)
            )
            arctic_gateway_domain = (
                ocean
                & (grid.lat >= 48.0)
                & ((basin_label == 0) | (basin_label == 3))
            )
            arctic_gateway_target = (
                ocean
                & (basin_label == 3)
                & (grid.lat >= 62.0)
                & (grid.lat <= 80.0)
            )
            if atlantic_gateway_sources.any() and arctic_gateway_target.any():
                gateway_distance = self._graph_hop_distance(
                    grid, atlantic_gateway_sources, arctic_gateway_domain)
                gateway_distance_gate = np.where(
                    gateway_distance >= 0,
                    np.exp(-((np.maximum(gateway_distance, 0) / 6.0) ** 1.35)),
                    0.0,
                )
                atlantic_arctic_gateway_low = (
                    arctic_gateway_target.astype(float)
                    * winter_hemi
                    * gateway_distance_gate
                    * np.clip((lat_abs - 62.0) / 6.0, 0.0, 1.0)
                    * np.clip((82.0 - lat_abs) / 8.0, 0.0, 1.0)
                    * (0.45 + 0.55 * np.clip(shelf_index, 0.0, 1.0))
                    * (0.45 + 0.55 * np.clip(front_support / 0.8, 0.0, 1.4))
                )
                low_support = np.maximum(
                    low_support,
                    np.clip(0.75 * atlantic_arctic_gateway_low, 0.0, 0.58),
                )
            if southern_mask.any():
                low_support[southern_mask] *= southern_gate[southern_mask]
                active = low_support[ocean & np.isfinite(low_support)]
                if active.size:
                    scale = float(np.percentile(active, 95))
                    if scale > 1.0e-9:
                        low_support = np.where(
                            ocean, np.clip(low_support / scale, 0.0, 1.25), 0.0)

            southern_shoulder_low = (
                southern_mask.astype(float)
                * shoulder_hemi
                * (0.40 + 0.60 * np.clip(southern_gate, 0.0, 1.2))
                * (0.50 + 0.50 * np.clip(front_support / 0.75, 0.0, 1.3))
                * (0.35 + 0.65 * np.clip(shelf_index, 0.0, 1.0))
            )
            low_support = np.maximum(
                low_support,
                np.clip(0.42 * southern_shoulder_low, 0.0, 0.55),
            )
            son_north_atlantic_autumn_low_support = np.zeros(
                grid.n, dtype=np.float64)
            if season == 3:
                son_north_atlantic_lat_gate = (
                    np.exp(-((grid.lat - 59.0) / 10.0) ** 2)
                    * np.clip((grid.lat - 43.0) / 8.0, 0.0, 1.0)
                    * np.clip((74.0 - grid.lat) / 8.0, 0.0, 1.0)
                )
                son_north_atlantic_lon_gate = np.maximum(
                    np.exp(
                        -(
                            (
                                ((grid.lon + 25.0 + 180.0) % 360.0) - 180.0
                            )
                            / 34.0
                        ) ** 2
                    ),
                    0.62
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 10.0 + 180.0) % 360.0) - 180.0
                            )
                            / 38.0
                        ) ** 2
                    ),
                )
                son_north_atlantic_west_taper = np.clip(
                    (grid.lon + 58.0) / 20.0,
                    0.0,
                    1.0,
                )
                son_north_atlantic_east_taper = (
                    1.0 - np.clip((grid.lon - 28.0) / 18.0, 0.0, 1.0)
                )
                son_north_atlantic_autumn_low_support = np.where(
                    ocean
                    & (basin_label == 0)
                    & (grid.lat >= 43.0)
                    & (grid.lat <= 74.0)
                    & (grid.lon >= -70.0)
                    & (grid.lon <= 42.0),
                    son_north_atlantic_lat_gate
                    * son_north_atlantic_lon_gate
                    * son_north_atlantic_west_taper
                    * son_north_atlantic_east_taper
                    * (0.52 + 0.48 * np.clip(shelf_index, 0.0, 1.0))
                    * (0.58 + 0.42 * np.clip(sst_front[season] / 0.70, 0.0, 1.25))
                    * (
                        0.70
                        + 0.30 * np.clip(same_lat_sst[season] / 0.20, 0.0, 1.1)
                    ),
                    0.0,
                )
                low_support = np.maximum(
                    low_support,
                    np.clip(0.62 * son_north_atlantic_autumn_low_support, 0.0, 0.72),
                )
                son_north_atlantic_autumn_low_transfer_field[season] = np.clip(
                    son_north_atlantic_autumn_low_support,
                    0.0,
                    1.15,
                )

            high_score = (
                open_ocean
                * basin_scale
                * subtropical
                * (0.55 + 0.45 * summer_hemi)
                * (1.0 - 0.25 * np.clip(sst_front[season], 0.0, 1.0))
            )
            high_support = self._pressure_object_support(
                grid,
                high_score,
                ocean,
                basin_label,
                quantile=80.0,
                max_components=2,
                min_cells=min_cells,
                smooth_passes=6,
                smooth_alpha=0.16,
            )
            mam_subpolar_high_lat = (
                np.exp(-((grid.lat - 66.0) / 15.0) ** 2)
                * np.clip((grid.lat - 43.0) / 10.0, 0.0, 1.0)
                * np.clip((86.0 - grid.lat) / 7.0, 0.0, 1.0)
            )
            mam_subpolar_high_support = np.where(
                ocean
                & (season == 1)
                & (grid.lat >= 45.0)
                & (grid.lat <= 84.0),
                mam_subpolar_high_lat
                * (0.48 + 0.52 * np.clip(shelf_index, 0.0, 1.0))
                * (0.45 + 0.55 * np.clip(sst_front[season] / 0.75, 0.0, 1.3)),
                0.0,
            )
            pacific_subpolar_high_gate = np.maximum(
                np.exp(
                    -(
                        (
                            ((grid.lon + 160.0 + 180.0) % 360.0) - 180.0
                        )
                        / 58.0
                    ) ** 2
                ),
                0.70
                * np.exp(
                    -(
                        (
                            ((grid.lon - 165.0 + 180.0) % 360.0) - 180.0
                        )
                        / 48.0
                    ) ** 2
                ),
            )
            atlantic_subpolar_high_gate = np.maximum(
                np.exp(
                    -(
                        (
                            ((grid.lon + 35.0 + 180.0) % 360.0) - 180.0
                        )
                        / 56.0
                    ) ** 2
                ),
                0.62
                * np.exp(
                    -(
                        (
                            ((grid.lon - 20.0 + 180.0) % 360.0) - 180.0
                        )
                        / 48.0
                    ) ** 2
                ),
            )
            jja_subpolar_high_lat = (
                np.exp(-((grid.lat - 55.0) / 12.5) ** 2)
                * np.clip((grid.lat - 39.0) / 8.0, 0.0, 1.0)
                * np.clip((73.0 - grid.lat) / 8.0, 0.0, 1.0)
            )
            jja_subpolar_high_support = np.where(
                ocean
                & (season == 2)
                & (grid.lat >= 40.0)
                & (grid.lat <= 73.0),
                jja_subpolar_high_lat
                * np.maximum(
                    pacific_subpolar_high_gate,
                    atlantic_subpolar_high_gate,
                )
                * (0.42 + 0.58 * np.clip(sst_front[season] / 0.70, 0.0, 1.25))
                * (0.45 + 0.55 * np.clip(shelf_index, 0.0, 1.0)),
                0.0,
            )
            nh_subpolar_high_transfer_support = np.maximum(
                mam_subpolar_high_support,
                jja_subpolar_high_support,
            )
            high_support = np.maximum(
                high_support,
                np.clip(0.92 * nh_subpolar_high_transfer_support, 0.0, 0.95),
            )
            subpolar_high_transfer_field[season] = np.clip(
                nh_subpolar_high_transfer_support,
                0.0,
                1.25,
            )
            mam_arctic_freeze_high_support = np.zeros(grid.n, dtype=np.float64)
            if season == 1:
                mam_arctic_lat_gate = (
                    np.exp(-((grid.lat - 73.0) / 11.5) ** 2)
                    * np.clip((grid.lat - 58.0) / 9.0, 0.0, 1.0)
                    * np.clip((86.0 - grid.lat) / 6.0, 0.0, 1.0)
                )
                mam_arctic_freeze = np.clip(
                    (276.0 - seasonal_sst[season]) / 6.0,
                    0.0,
                    1.25,
                )
                beaufort_gate = np.exp(
                    -(
                        (
                            ((grid.lon + 145.0 + 180.0) % 360.0) - 180.0
                        )
                        / 45.0
                    ) ** 2
                )
                greenland_sea_gate = np.exp(
                    -(
                        (
                            ((grid.lon + 20.0 + 180.0) % 360.0) - 180.0
                        )
                        / 52.0
                    ) ** 2
                )
                baffin_gate = (
                    np.exp(
                        -(
                            (
                                ((grid.lon + 60.0 + 180.0) % 360.0) - 180.0
                            )
                            / 38.0
                        ) ** 2
                    )
                    * np.clip((76.0 - grid.lat) / 8.0, 0.0, 1.0)
                )
                arctic_lon_gate = np.maximum.reduce([
                    0.45 * (basin_label == 3).astype(float),
                    beaufort_gate,
                    0.88 * greenland_sea_gate,
                    0.95 * baffin_gate,
                ])
                mam_arctic_polar_high = np.where(
                    ocean
                    & (grid.lat >= 55.0)
                    & (grid.lat <= 86.0),
                    mam_arctic_lat_gate
                    * mam_arctic_freeze
                    * arctic_lon_gate
                    * (0.55 + 0.45 * np.clip(shelf_index, 0.0, 1.0))
                    * (0.65 + 0.35 * np.clip(sst_front[season] / 0.75, 0.0, 1.25)),
                    0.0,
                )
                baffin_labrador_lat_gate = (
                    np.exp(-((grid.lat - 64.0) / 8.5) ** 2)
                    * np.clip((grid.lat - 53.0) / 6.0, 0.0, 1.0)
                    * np.clip((77.0 - grid.lat) / 6.0, 0.0, 1.0)
                )
                baffin_labrador_lon_gate = np.maximum(
                    np.exp(
                        -(
                            (
                                ((grid.lon + 62.0 + 180.0) % 360.0) - 180.0
                            )
                            / 30.0
                        ) ** 2
                    ),
                    0.55
                    * np.exp(
                        -(
                            (
                                ((grid.lon + 45.0 + 180.0) % 360.0) - 180.0
                            )
                            / 34.0
                        ) ** 2
                    ),
                )
                baffin_labrador_freeze = np.clip(
                    (278.0 - seasonal_sst[season]) / 7.5,
                    0.0,
                    1.15,
                )
                baffin_labrador_high = np.where(
                    ocean
                    & (grid.lat >= 53.0)
                    & (grid.lat <= 77.0),
                    baffin_labrador_lat_gate
                    * baffin_labrador_lon_gate
                    * baffin_labrador_freeze
                    * (0.50 + 0.50 * np.clip(shelf_index, 0.0, 1.0))
                    * (0.45 + 0.55 * np.clip(sst_front[season] / 0.8, 0.0, 1.25)),
                    0.0,
                )
                mam_arctic_freeze_high_support = np.maximum(
                    mam_arctic_polar_high,
                    0.82 * baffin_labrador_high,
                )
                high_support = np.maximum(
                    high_support,
                    np.clip(mam_arctic_freeze_high_support, 0.0, 1.05),
                )
                mam_arctic_freeze_high_transfer_field[season] = np.clip(
                    mam_arctic_freeze_high_support,
                    0.0,
                    1.20,
                )
                mam_canadian_arctic_high_lat_gate = (
                    np.exp(-((grid.lat - 74.0) / 9.5) ** 2)
                    * np.clip((grid.lat - 62.0) / 7.0, 0.0, 1.0)
                    * np.clip((86.0 - grid.lat) / 6.0, 0.0, 1.0)
                )
                canadian_archipelago_gate = np.maximum.reduce([
                    np.exp(
                        -(
                            (
                                ((grid.lon + 145.0 + 180.0) % 360.0) - 180.0
                            )
                            / 44.0
                        ) ** 2
                    ),
                    0.82
                    * np.exp(
                        -(
                            (
                                ((grid.lon + 95.0 + 180.0) % 360.0) - 180.0
                            )
                            / 38.0
                        ) ** 2
                    ),
                    0.85
                    * np.exp(
                        -(
                            (
                                ((grid.lon + 60.0 + 180.0) % 360.0) - 180.0
                            )
                            / 34.0
                        ) ** 2
                    )
                    * np.clip((78.0 - grid.lat) / 9.0, 0.0, 1.0),
                    0.42
                    * np.exp(
                        -(
                            (
                                ((grid.lon + 25.0 + 180.0) % 360.0) - 180.0
                            )
                            / 44.0
                        ) ** 2
                    )
                    * np.clip((grid.lat - 72.0) / 5.0, 0.0, 1.0),
                ])
                canadian_arctic_east_taper = np.where(
                    grid.lon < -35.0,
                    1.0,
                    1.0 - np.clip((grid.lon + 35.0) / 70.0, 0.0, 1.0),
                )
                canadian_arctic_east_taper = np.clip(
                    canadian_arctic_east_taper,
                    0.0,
                    1.0,
                )
                mam_canadian_arctic_high_support = np.where(
                    ocean
                    & (grid.lat >= 58.0)
                    & (grid.lat <= 86.0)
                    & ((basin_label == 3) | (basin_label == 0)),
                    mam_canadian_arctic_high_lat_gate
                    * canadian_archipelago_gate
                    * canadian_arctic_east_taper
                    * baffin_labrador_freeze
                    * (0.55 + 0.45 * np.clip(shelf_index, 0.0, 1.0))
                    * (0.50 + 0.50 * np.clip(sst_front[season] / 0.8, 0.0, 1.25)),
                    0.0,
                )
                mam_canadian_arctic_high_support = np.clip(
                    mam_canadian_arctic_high_support,
                    0.0,
                    1.05,
                )
                mam_central_arctic_lat_gate = (
                    np.exp(-((grid.lat - 77.0) / 6.8) ** 2)
                    * np.clip((grid.lat - 69.0) / 5.0, 0.0, 1.0)
                    * np.clip((86.0 - grid.lat) / 5.0, 0.0, 1.0)
                )
                mam_central_arctic_lon_gate = np.maximum.reduce([
                    np.exp(
                        -(
                            (
                                ((grid.lon + 105.0 + 180.0) % 360.0) - 180.0
                            )
                            / 46.0
                        ) ** 2
                    ),
                    0.85
                    * np.exp(
                        -(
                            (
                                ((grid.lon + 55.0 + 180.0) % 360.0) - 180.0
                            )
                            / 42.0
                        ) ** 2
                    ),
                    0.42
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 170.0 + 180.0) % 360.0) - 180.0
                            )
                            / 34.0
                        ) ** 2
                    ),
                ])
                mam_central_arctic_barents_cut = np.where(
                    grid.lon < -25.0,
                    1.0,
                    1.0 - np.clip((grid.lon + 25.0) / 45.0, 0.0, 1.0),
                )
                mam_central_arctic_barents_cut = np.where(
                    grid.lon > 135.0,
                    1.0,
                    mam_central_arctic_barents_cut,
                )
                mam_central_arctic_barents_cut = np.clip(
                    mam_central_arctic_barents_cut,
                    0.0,
                    1.0,
                )
                mam_central_arctic_high_support = np.where(
                    ocean
                    & (basin_label == 3)
                    & (grid.lat >= 68.0)
                    & (grid.lat <= 86.0),
                    mam_central_arctic_lat_gate
                    * mam_central_arctic_lon_gate
                    * mam_central_arctic_barents_cut
                    * baffin_labrador_freeze
                    * (0.62 + 0.38 * np.clip(1.0 - shelf_index, 0.0, 1.0))
                    * (0.75 + 0.25 * np.clip(sst_front[season] / 0.8, 0.0, 1.2)),
                    0.0,
                )
                mam_central_arctic_high_support = np.clip(
                    mam_central_arctic_high_support,
                    0.0,
                    1.0,
                )
                high_support = np.maximum.reduce([
                    high_support,
                    np.clip(mam_canadian_arctic_high_support, 0.0, 1.05),
                    np.clip(mam_central_arctic_high_support, 0.0, 1.00),
                ])
                mam_canadian_arctic_high_transfer_field[season] = (
                    mam_canadian_arctic_high_support
                )
                mam_central_arctic_high_transfer_field[season] = (
                    mam_central_arctic_high_support
                )
            jja_north_pacific_high_support = np.zeros(grid.n, dtype=np.float64)
            if season == 2:
                jja_pacific_lat_gate = (
                    np.exp(-((grid.lat - 53.5) / 10.5) ** 2)
                    * np.clip((grid.lat - 41.0) / 7.0, 0.0, 1.0)
                    * np.clip((67.0 - grid.lat) / 8.0, 0.0, 1.0)
                )
                gulf_alaska_gate = np.exp(
                    -(
                        (
                            ((grid.lon + 145.0 + 180.0) % 360.0) - 180.0
                        )
                        / 42.0
                    ) ** 2
                )
                aleutian_gate = np.exp(
                    -(
                        (
                            ((grid.lon - 172.0 + 180.0) % 360.0) - 180.0
                        )
                        / 55.0
                    ) ** 2
                )
                northwest_pacific_gate = np.exp(
                    -(
                        (
                            ((grid.lon - 150.0 + 180.0) % 360.0) - 180.0
                        )
                        / 38.0
                    ) ** 2
                )
                jja_pacific_lon_gate = np.maximum.reduce([
                    gulf_alaska_gate,
                    0.88 * aleutian_gate,
                    0.35 * northwest_pacific_gate,
                ])
                jja_north_pacific_high_support = np.where(
                    ocean
                    & (basin_label == 1)
                    & (grid.lat >= 41.0)
                    & (grid.lat <= 67.0),
                    jja_pacific_lat_gate
                    * jja_pacific_lon_gate
                    * (0.50 + 0.50 * np.clip(shelf_index, 0.0, 1.0))
                    * (0.50 + 0.50 * np.clip(sst_front[season] / 0.65, 0.0, 1.25))
                    * (
                        0.72
                        + 0.28 * np.clip(-same_lat_sst[season] / 0.25, 0.0, 1.0)
                    ),
                    0.0,
                )
                high_support = np.maximum(
                    high_support,
                    np.clip(jja_north_pacific_high_support, 0.0, 1.05),
                )
                jja_north_pacific_high_transfer_field[season] = np.clip(
                    jja_north_pacific_high_support,
                    0.0,
                    1.20,
                )
                jja_central_pacific_lat_gate = (
                    np.exp(-((grid.lat - 49.0) / 11.5) ** 2)
                    * np.clip((grid.lat - 31.0) / 9.0, 0.0, 1.0)
                    * np.clip((66.0 - grid.lat) / 9.0, 0.0, 1.0)
                )
                jja_central_pacific_lon_gate = np.maximum.reduce([
                    np.exp(
                        -(
                            (
                                ((grid.lon + 150.0 + 180.0) % 360.0) - 180.0
                            )
                            / 38.0
                        ) ** 2
                    ),
                    0.82
                    * np.exp(
                        -(
                            (
                                ((grid.lon + 172.0 + 180.0) % 360.0) - 180.0
                            )
                            / 35.0
                        ) ** 2
                    ),
                    0.30
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 175.0 + 180.0) % 360.0) - 180.0
                            )
                            / 28.0
                        ) ** 2
                    ),
                ])
                jja_west_pacific_suppression = (
                    1.0
                    - 0.85
                    * np.clip((grid.lon - 130.0) / 32.0, 0.0, 1.0)
                    * np.clip((162.0 - grid.lon) / 32.0, 0.0, 1.0)
                )
                jja_west_pacific_suppression = np.where(
                    (grid.lon >= 130.0) & (grid.lon <= 162.0),
                    jja_west_pacific_suppression,
                    1.0,
                )
                jja_north_pacific_arctic_taper = (
                    1.0 - np.clip((grid.lat - 61.0) / 8.0, 0.0, 1.0)
                )
                jja_north_pacific_central_high_support = np.where(
                    ocean
                    & (basin_label == 1)
                    & (grid.lat >= 31.0)
                    & (grid.lat <= 66.0)
                    & ((grid.lon <= -110.0) | (grid.lon >= 162.0)),
                    jja_central_pacific_lat_gate
                    * jja_central_pacific_lon_gate
                    * jja_west_pacific_suppression
                    * jja_north_pacific_arctic_taper
                    * (0.68 + 0.32 * np.clip(shelf_index, 0.0, 1.0))
                    * (0.55 + 0.45 * np.clip(sst_front[season] / 0.65, 0.0, 1.3))
                    * (
                        0.72
                        + 0.28 * np.clip(-same_lat_sst[season] / 0.25, 0.0, 1.1)
                    ),
                    0.0,
                )
                jja_north_pacific_central_high_support = np.clip(
                    jja_north_pacific_central_high_support,
                    0.0,
                    1.0,
                )
                high_support = np.maximum(
                    high_support,
                    np.clip(jja_north_pacific_central_high_support, 0.0, 1.0),
                )
                jja_north_pacific_central_high_transfer_field[season] = (
                    jja_north_pacific_central_high_support
                )
                jja_eurasian_low_lat_gate = (
                    np.maximum(
                        np.exp(-((grid.lat - 31.0) / 16.0) ** 2),
                        0.72 * np.exp(-((grid.lat - 45.0) / 13.0) ** 2),
                    )
                    * np.clip((grid.lat - 8.0) / 7.0, 0.0, 1.0)
                    * np.clip((62.0 - grid.lat) / 10.0, 0.0, 1.0)
                )
                jja_eurasian_low_lon_gate = np.maximum.reduce([
                    0.95
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 48.0 + 180.0) % 360.0) - 180.0
                            )
                            / 25.0
                        ) ** 2
                    ),
                    0.95
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 78.0 + 180.0) % 360.0) - 180.0
                            )
                            / 24.0
                        ) ** 2
                    ),
                    np.exp(
                        -(
                            (
                                ((grid.lon - 116.0 + 180.0) % 360.0) - 180.0
                            )
                            / 30.0
                        ) ** 2
                    ),
                    0.68
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 138.0 + 180.0) % 360.0) - 180.0
                            )
                            / 25.0
                        ) ** 2
                    ),
                ])
                jja_eurasian_low_elevation_gate = (
                    1.0 - np.clip((elev - 1600.0) / 2400.0, 0.0, 1.0)
                )
                jja_eurasian_thermal_gate = np.clip(
                    (temp_anom[season] - 2.0) / 8.5,
                    0.0,
                    1.35,
                )
                jja_eurasian_thermal_low_support = np.where(
                    land
                    & (grid.lat >= 8.0)
                    & (grid.lat <= 62.0)
                    & (grid.lon >= 25.0)
                    & (grid.lon <= 155.0),
                    jja_eurasian_low_lat_gate
                    * jja_eurasian_low_lon_gate
                    * jja_eurasian_thermal_gate
                    * (0.55 + 0.45 * jja_eurasian_low_elevation_gate)
                    * (
                        0.62
                        + 0.38
                        * np.clip(continent_interiority + coast_strength, 0.0, 1.0)
                    )
                    * (1.0 - 0.18 * np.clip(barrier / 0.90, 0.0, 1.0)),
                    0.0,
                )
                jja_eurasian_thermal_low_transfer_field[season] = np.clip(
                    jja_eurasian_thermal_low_support,
                    0.0,
                    1.25,
                )
                jja_west_pacific_low_lat_gate = (
                    np.exp(-((grid.lat - 43.0) / 13.0) ** 2)
                    * np.clip((grid.lat - 24.0) / 8.0, 0.0, 1.0)
                    * np.clip((64.0 - grid.lat) / 8.0, 0.0, 1.0)
                )
                jja_west_pacific_low_lon_gate = np.maximum(
                    np.exp(
                        -(
                            (
                                ((grid.lon - 144.0 + 180.0) % 360.0) - 180.0
                            )
                            / 25.0
                        ) ** 2
                    ),
                    0.62
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 161.0 + 180.0) % 360.0) - 180.0
                            )
                            / 23.0
                        ) ** 2
                    ),
                )
                jja_west_pacific_low_west_taper = np.clip(
                    (grid.lon - 118.0) / 14.0,
                    0.0,
                    1.0,
                )
                jja_west_pacific_low_east_taper = (
                    1.0 - np.clip((grid.lon - 174.0) / 16.0, 0.0, 1.0)
                )
                jja_west_pacific_marginal_low_support = np.where(
                    ocean
                    & (grid.lat >= 24.0)
                    & (grid.lat <= 64.0)
                    & (grid.lon >= 118.0)
                    & (grid.lon <= 178.0)
                    & ((basin_label == 1) | (basin_label == 2)),
                    jja_west_pacific_low_lat_gate
                    * jja_west_pacific_low_lon_gate
                    * jja_west_pacific_low_west_taper
                    * jja_west_pacific_low_east_taper
                    * (0.48 + 0.52 * np.clip(shelf_index, 0.0, 1.0))
                    * (0.58 + 0.42 * np.clip(sst_front[season] / 0.65, 0.0, 1.25))
                    * (
                        0.76
                        + 0.24 * np.clip(same_lat_sst[season] / 0.20, 0.0, 1.1)
                    ),
                    0.0,
                )
                low_support = np.maximum(
                    low_support,
                    np.clip(0.68 * jja_west_pacific_marginal_low_support, 0.0, 0.75),
                )
                jja_west_pacific_marginal_low_transfer_field[season] = np.clip(
                    jja_west_pacific_marginal_low_support,
                    0.0,
                    1.15,
                )
                jja_west_asia_low_lat_gate = (
                    np.exp(-((grid.lat - 29.0) / 10.5) ** 2)
                    * np.clip((grid.lat - 13.0) / 7.0, 0.0, 1.0)
                    * np.clip((43.0 - grid.lat) / 7.0, 0.0, 1.0)
                )
                jja_west_asia_low_lon_gate = np.maximum.reduce([
                    np.exp(
                        -(
                            (
                                ((grid.lon - 47.0 + 180.0) % 360.0) - 180.0
                            )
                            / 19.0
                        ) ** 2
                    ),
                    0.90
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 60.0 + 180.0) % 360.0) - 180.0
                            )
                            / 18.0
                        ) ** 2
                    ),
                    0.58
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 73.0 + 180.0) % 360.0) - 180.0
                            )
                            / 16.0
                        ) ** 2
                    ),
                ])
                jja_west_asia_thermal_gate = np.clip(
                    (temp_anom[season] - 2.5) / 7.0,
                    0.0,
                    1.35,
                )
                jja_west_asia_low_elevation_gate = (
                    1.0 - np.clip((elev - 1800.0) / 2500.0, 0.0, 1.0)
                )
                jja_west_asia_dry_lowland_gate = (
                    (0.62 + 0.38 * np.clip(1.0 - coast_strength, 0.0, 1.0))
                    * (
                        0.72
                        + 0.28 * np.clip(1.0 - continent_interiority, 0.0, 1.0)
                    )
                )
                jja_west_asia_thermal_low_support = np.where(
                    land
                    & (grid.lat >= 13.0)
                    & (grid.lat <= 43.0)
                    & (grid.lon >= 28.0)
                    & (grid.lon <= 88.0),
                    jja_west_asia_low_lat_gate
                    * jja_west_asia_low_lon_gate
                    * jja_west_asia_thermal_gate
                    * (0.55 + 0.45 * jja_west_asia_low_elevation_gate)
                    * jja_west_asia_dry_lowland_gate
                    * (1.0 - 0.18 * np.clip(barrier / 0.90, 0.0, 1.0)),
                    0.0,
                )
                jja_west_asia_thermal_low_transfer_field[season] = np.clip(
                    jja_west_asia_thermal_low_support,
                    0.0,
                    1.15,
                )
                jja_east_china_sea_low_lat_gate = (
                    np.exp(-((grid.lat - 33.0) / 8.5) ** 2)
                    * np.clip((grid.lat - 21.0) / 6.0, 0.0, 1.0)
                    * np.clip((47.0 - grid.lat) / 7.0, 0.0, 1.0)
                )
                jja_east_china_sea_low_lon_gate = np.maximum(
                    np.exp(
                        -(
                            (
                                ((grid.lon - 128.0 + 180.0) % 360.0) - 180.0
                            )
                            / 17.0
                        ) ** 2
                    ),
                    0.74
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 140.0 + 180.0) % 360.0) - 180.0
                            )
                            / 16.0
                        ) ** 2
                    ),
                )
                jja_east_china_sea_low_support = np.where(
                    ocean
                    & (grid.lat >= 21.0)
                    & (grid.lat <= 47.0)
                    & (grid.lon >= 116.0)
                    & (grid.lon <= 154.0)
                    & ((basin_label == 1) | (basin_label == 2)),
                    jja_east_china_sea_low_lat_gate
                    * jja_east_china_sea_low_lon_gate
                    * (0.58 + 0.42 * np.clip(shelf_index / 0.45, 0.0, 1.25))
                    * (0.54 + 0.46 * np.clip(sst_front[season] / 0.55, 0.0, 1.25))
                    * (
                        0.70
                        + 0.30 * np.clip(same_lat_sst[season] / 0.25, 0.0, 1.1)
                    ),
                    0.0,
                )
                low_support = np.maximum(
                    low_support,
                    np.clip(0.75 * jja_east_china_sea_low_support, 0.0, 0.78),
                )
                jja_east_china_sea_low_transfer_field[season] = np.clip(
                    jja_east_china_sea_low_support,
                    0.0,
                    1.10,
                )
            djf_atlantic_gateway_low_support = np.zeros(
                grid.n, dtype=np.float64)
            djf_north_pacific_low_support = np.zeros(
                grid.n, dtype=np.float64)
            if season == 0:
                djf_atlantic_gateway_lat_gate = (
                    np.exp(-((grid.lat - 66.0) / 10.0) ** 2)
                    * np.clip((grid.lat - 55.0) / 7.0, 0.0, 1.0)
                    * np.clip((82.0 - grid.lat) / 8.0, 0.0, 1.0)
                )
                djf_atlantic_gateway_lon_gate = np.maximum.reduce([
                    np.exp(
                        -(
                            (
                                ((grid.lon + 20.0 + 180.0) % 360.0) - 180.0
                            )
                            / 42.0
                        ) ** 2
                    ),
                    0.78
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 45.0 + 180.0) % 360.0) - 180.0
                            )
                            / 45.0
                        ) ** 2
                    ),
                    0.55
                    * np.exp(
                        -(
                            (
                                ((grid.lon + 50.0 + 180.0) % 360.0) - 180.0
                            )
                            / 42.0
                        ) ** 2
                    ),
                ])
                djf_atlantic_gateway_base = np.where(
                    ocean
                    & (grid.lat >= 55.0)
                    & (grid.lat <= 82.0)
                    & (grid.lon >= -75.0)
                    & (grid.lon <= 95.0)
                    & ((basin_label == 0) | (basin_label == 3)),
                    djf_atlantic_gateway_lat_gate
                    * djf_atlantic_gateway_lon_gate
                    * (0.45 + 0.55 * np.clip(shelf_index, 0.0, 1.0))
                    * (0.45 + 0.55 * np.clip(sst_front[season] / 0.8, 0.0, 1.4)),
                    0.0,
                )
                djf_atlantic_arctic_lat_gate = (
                    np.exp(-((grid.lat - 72.0) / 8.5) ** 2)
                    * np.clip((grid.lat - 62.0) / 5.0, 0.0, 1.0)
                    * np.clip((82.0 - grid.lat) / 7.0, 0.0, 1.0)
                )
                djf_atlantic_arctic_lon_gate = np.maximum(
                    np.exp(
                        -(
                            (
                                ((grid.lon + 15.0 + 180.0) % 360.0) - 180.0
                            )
                            / 48.0
                        ) ** 2
                    ),
                    0.70
                    * np.exp(
                        -(
                            (
                                ((grid.lon - 45.0 + 180.0) % 360.0) - 180.0
                            )
                            / 50.0
                        ) ** 2
                    ),
                )
                djf_atlantic_arctic_branch = np.where(
                    ocean
                    & (basin_label == 3)
                    & (grid.lat >= 63.0)
                    & (grid.lat <= 82.0)
                    & (grid.lon >= -45.0)
                    & (grid.lon <= 90.0),
                    djf_atlantic_arctic_lat_gate
                    * djf_atlantic_arctic_lon_gate
                    * (0.50 + 0.50 * np.clip(shelf_index, 0.0, 1.0))
                    * (0.55 + 0.45 * np.clip(sst_front[season] / 0.8, 0.0, 1.4)),
                    0.0,
                )
                djf_atlantic_gateway_low_support = np.clip(
                    np.maximum(
                        djf_atlantic_gateway_base,
                        0.75 * djf_atlantic_arctic_branch,
                    ),
                    0.0,
                    1.15,
                )
                djf_north_pacific_lat_gate = (
                    np.exp(-((grid.lat - 53.5) / 8.5) ** 2)
                    * np.clip((grid.lat - 39.0) / 6.0, 0.0, 1.0)
                    * np.clip((66.0 - grid.lat) / 6.0, 0.0, 1.0)
                )
                djf_north_pacific_lon_gate = np.maximum(
                    np.exp(
                        -(
                            (
                                ((grid.lon - 175.0 + 180.0) % 360.0) - 180.0
                            )
                            / 46.0
                        ) ** 2
                    ),
                    0.62
                    * np.exp(
                        -(
                            (
                                ((grid.lon + 165.0 + 180.0) % 360.0) - 180.0
                            )
                            / 35.0
                        ) ** 2
                    ),
                )
                djf_bering_overlap = np.maximum(
                    np.exp(
                        -(
                            (
                                ((grid.lon - 175.0 + 180.0) % 360.0) - 180.0
                            )
                            / 32.0
                        ) ** 2
                    ),
                    np.exp(
                        -(
                            (
                                ((grid.lon + 170.0 + 180.0) % 360.0) - 180.0
                            )
                            / 30.0
                        ) ** 2
                    ),
                )
                djf_bering_suppression = (
                    1.0
                    - 0.55
                    * np.clip((grid.lat - 58.0) / 10.0, 0.0, 1.0)
                    * np.clip(djf_bering_overlap, 0.0, 1.0)
                )
                djf_north_pacific_low_support = np.where(
                    ocean
                    & (basin_label == 1)
                    & (grid.lat >= 39.0)
                    & (grid.lat <= 68.0),
                    djf_north_pacific_lat_gate
                    * djf_north_pacific_lon_gate
                    * djf_bering_suppression
                    * (0.55 + 0.45 * np.clip(sst_front[season] / 0.65, 0.0, 1.3))
                    * (0.48 + 0.52 * np.clip(shelf_index, 0.0, 1.0)),
                    0.0,
                )
                djf_north_pacific_low_support = np.clip(
                    djf_north_pacific_low_support,
                    0.0,
                    1.0,
                )
                djf_atlantic_gateway_low_transfer_field[season] = (
                    djf_atlantic_gateway_low_support
                )
                djf_north_pacific_low_transfer_field[season] = (
                    djf_north_pacific_low_support
                )

            land_center_score = (
                np.abs(land_thermal)
                * land.astype(float)
                * (0.45 + 0.55 * np.clip(continent_interiority, 0.0, 1.0))
                * (0.65 + 0.35 * np.clip(barrier, 0.0, 1.0))
            )
            land_center_support = self._pressure_object_support(
                grid,
                land_center_score,
                land,
                continent_label,
                quantile=68.0,
                max_components=2,
                min_cells=min_cells,
                smooth_passes=4,
                smooth_alpha=0.16,
                multi_large_labels=True,
            )
            land_center_signed = np.sign(land_thermal) * land_center_support

            spring_land_high_support = np.zeros(grid.n, dtype=np.float64)
            autumn_land_high_decay_support = np.zeros(grid.n, dtype=np.float64)
            if season == 1:
                spring_lat_gate = (
                    np.exp(-((grid.lat - 61.0) / 16.0) ** 2)
                    * np.clip((grid.lat - 42.0) / 9.0, 0.0, 1.0)
                    * np.clip((86.0 - grid.lat) / 8.0, 0.0, 1.0)
                )
                spring_cold_memory = np.clip(-temp_anom[season] / 4.0, 0.0, 1.2)
                spring_freeze_memory = np.clip(
                    (276.0 - seasonal_T[season]) / 14.0,
                    0.0,
                    1.25,
                )
                low_interior = np.clip(1.0 - continent_interiority, 0.0, 1.0)
                spring_broad_cold = (
                    spring_lat_gate
                    * (0.70 * spring_cold_memory + 0.30 * spring_freeze_memory)
                    * (0.42 + 0.58 * low_interior)
                )
                spring_snow_shield = (
                    spring_lat_gate
                    * (0.42 + 0.58 * spring_freeze_memory)
                    * (0.24 + 0.76 * low_interior ** 1.35)
                    * (0.72 + 0.28 * np.clip(coast_strength / 0.25, 0.0, 1.0))
                )
                spring_land_high_support = np.where(
                    land
                    & (grid.lat >= 42.0)
                    & (grid.lat <= 86.0),
                    0.50 * spring_broad_cold + 0.70 * spring_snow_shield,
                    0.0,
                )
            elif season == 3:
                autumn_lat_gate = (
                    np.exp(-((grid.lat - 57.0) / 17.0) ** 2)
                    * np.clip((grid.lat - 40.0) / 9.0, 0.0, 1.0)
                    * np.clip((82.0 - grid.lat) / 10.0, 0.0, 1.0)
                )
                summer_heat_memory = np.clip(temp_anom[2] / 10.0, 0.0, 1.3)
                autumn_cooling_memory = np.clip(
                    (temp_anom[2] - temp_anom[season]) / 11.0,
                    0.0,
                    1.3,
                )
                low_elevation_gate = 1.0 - np.clip(
                    (elev - 900.0) / 1500.0,
                    0.0,
                    1.0,
                )
                unfrozen_ground_memory = np.clip(
                    (seasonal_T[season] - 263.0) / 16.0,
                    0.0,
                    1.2,
                )
                autumn_land_high_decay_support = np.where(
                    land
                    & (grid.lat >= 40.0)
                    & (grid.lat <= 82.0),
                    autumn_lat_gate
                    * (0.52 * summer_heat_memory + 0.48 * autumn_cooling_memory)
                    * (0.72 + 0.28 * np.clip(continent_interiority, 0.0, 1.0))
                    * low_elevation_gate
                    * (0.58 + 0.42 * unfrozen_ground_memory),
                    0.0,
                )
                autumn_transition_maritime_gate = (
                    (0.55 + 0.45 * np.clip(1.0 - continent_interiority, 0.0, 1.0))
                    * (0.72 + 0.28 * np.clip(coast_strength / 0.35, 0.0, 1.0))
                )
                autumn_transition_barrier_escape = (
                    1.0 - 0.20 * np.clip(barrier / 0.90, 0.0, 1.0)
                )
                son_autumn_land_relief_support = np.where(
                    land
                    & (grid.lat >= 42.0)
                    & (grid.lat <= 78.0),
                    autumn_lat_gate
                    * (0.46 * summer_heat_memory + 0.54 * autumn_cooling_memory)
                    * low_elevation_gate ** 1.15
                    * (0.70 + 0.30 * unfrozen_ground_memory)
                    * autumn_transition_maritime_gate
                    * autumn_transition_barrier_escape,
                    0.0,
                )
                son_autumn_land_relief_transfer_field[season] = np.clip(
                    son_autumn_land_relief_support,
                    0.0,
                    1.20,
                )
            land_shoulder_phase_source = (
                0.38 * spring_land_high_support
                - 0.36 * autumn_land_high_decay_support
            )
            land_shoulder_phase_transfer_field[season] = np.clip(
                0.58 * spring_land_high_support
                - 0.55 * autumn_land_high_decay_support,
                -1.0,
                1.0,
            )
            djf_north_america_winter_high_relief_support = np.zeros(
                grid.n, dtype=np.float64)
            if season == 0:
                djf_north_america_relief_lat_gate = (
                    np.exp(-((grid.lat - 56.5) / 12.5) ** 2)
                    * np.clip((grid.lat - 42.0) / 7.0, 0.0, 1.0)
                    * np.clip((74.0 - grid.lat) / 8.0, 0.0, 1.0)
                )
                djf_north_america_relief_lon_gate = np.maximum(
                    np.exp(
                        -(
                            (
                                ((grid.lon + 103.0 + 180.0) % 360.0) - 180.0
                            )
                            / 45.0
                        ) ** 2
                    ),
                    0.55
                    * np.exp(
                        -(
                            (
                                ((grid.lon + 75.0 + 180.0) % 360.0) - 180.0
                            )
                            / 32.0
                        ) ** 2
                    ),
                )
                djf_north_america_low_elevation = (
                    1.0 - np.clip((elev - 1200.0) / 1600.0, 0.0, 1.0)
                ) ** 0.8
                djf_north_america_maritime_erosion = (
                    (0.45 + 0.55 * np.clip(coast_strength / 0.30, 0.0, 1.0))
                    * (
                        0.55
                        + 0.45
                        * np.clip((0.85 - continent_interiority) / 0.85, 0.0, 1.0)
                    )
                )
                djf_north_america_rocky_shelter = (
                    1.0 - 0.30 * np.clip(barrier / 0.8, 0.0, 1.0)
                )
                djf_north_america_winter_high_relief_support = np.where(
                    land
                    & (grid.lat >= 42.0)
                    & (grid.lat <= 74.0)
                    & (grid.lon >= -150.0)
                    & (grid.lon <= -50.0),
                    djf_north_america_relief_lat_gate
                    * djf_north_america_relief_lon_gate
                    * djf_north_america_low_elevation
                    * djf_north_america_maritime_erosion
                    * djf_north_america_rocky_shelter,
                    0.0,
                )
                djf_north_america_winter_high_relief_transfer_field[season] = (
                    np.clip(
                        djf_north_america_winter_high_relief_support,
                        0.0,
                        1.0,
                    )
                )
            mam_north_america_land_high_support = np.zeros(
                grid.n, dtype=np.float64)
            mam_north_america_plains_high_support = np.zeros(
                grid.n, dtype=np.float64)
            if season == 1:
                mam_north_america_lat_gate = (
                    np.exp(-((grid.lat - 57.0) / 12.5) ** 2)
                    * np.clip((grid.lat - 43.0) / 7.0, 0.0, 1.0)
                    * np.clip((73.0 - grid.lat) / 7.0, 0.0, 1.0)
                )
                mam_north_america_lon_gate = np.exp(
                    -(
                        (
                            ((grid.lon + 100.0 + 180.0) % 360.0) - 180.0
                        )
                        / 36.0
                    ) ** 2
                )
                mam_north_america_freeze_memory = np.clip(
                    (279.0 - seasonal_T[season]) / 13.0,
                    0.0,
                    1.25,
                )
                mam_north_america_cold_memory = np.clip(
                    -temp_anom[season] / 4.2,
                    0.0,
                    1.15,
                )
                mam_north_america_low_elevation = (
                    1.0 - np.clip((elev - 850.0) / 1250.0, 0.0, 1.0)
                ) ** 1.1
                mam_north_america_low_interior = (
                    0.55
                    + 0.45 * np.clip((0.65 - continent_interiority) / 0.65, 0.0, 1.0)
                )
                mam_north_america_coast_gate = (
                    0.70 + 0.30 * np.clip(coast_strength / 0.20, 0.0, 1.0)
                )
                mam_north_america_land_high_support = np.where(
                    land
                    & (grid.lat >= 43.0)
                    & (grid.lat <= 73.0)
                    & (grid.lon >= -165.0)
                    & (grid.lon <= -45.0),
                    mam_north_america_lat_gate
                    * mam_north_america_lon_gate
                    * (
                        0.58 * mam_north_america_freeze_memory
                        + 0.42 * mam_north_america_cold_memory
                    )
                    * mam_north_america_low_elevation
                    * mam_north_america_low_interior
                    * mam_north_america_coast_gate,
                    0.0,
                )
                mam_north_america_land_high_transfer_field[season] = np.clip(
                    mam_north_america_land_high_support,
                    0.0,
                    1.20,
                )
                mam_plains_lat_gate = (
                    np.exp(-((grid.lat - 51.5) / 13.0) ** 2)
                    * np.clip((grid.lat - 34.0) / 8.0, 0.0, 1.0)
                    * np.clip((70.0 - grid.lat) / 8.0, 0.0, 1.0)
                )
                mam_plains_lon_gate = np.maximum(
                    np.exp(
                        -(
                            (
                                ((grid.lon + 96.0 + 180.0) % 360.0) - 180.0
                            )
                            / 34.0
                        ) ** 2
                    ),
                    0.72
                    * np.exp(
                        -(
                            (
                                ((grid.lon + 78.0 + 180.0) % 360.0) - 180.0
                            )
                            / 30.0
                        ) ** 2
                    ),
                )
                mam_plains_low_elevation = (
                    1.0 - np.clip((elev - 1050.0) / 1500.0, 0.0, 1.0)
                ) ** 0.9
                mam_plains_low_interior = (
                    0.55
                    + 0.45 * np.clip((0.75 - continent_interiority) / 0.75, 0.0, 1.0)
                )
                mam_plains_coast_memory = (
                    0.78 + 0.22 * np.clip(coast_strength / 0.25, 0.0, 1.0)
                )
                mam_plains_terrain_shelter = (
                    1.0 - 0.25 * np.clip(barrier / 0.8, 0.0, 1.0)
                )
                mam_north_america_plains_high_support = np.where(
                    land
                    & (grid.lat >= 34.0)
                    & (grid.lat <= 70.0)
                    & (grid.lon >= -125.0)
                    & (grid.lon <= -55.0),
                    mam_plains_lat_gate
                    * mam_plains_lon_gate
                    * mam_plains_low_elevation
                    * mam_plains_low_interior
                    * mam_plains_coast_memory
                    * mam_plains_terrain_shelter,
                    0.0,
                )
                mam_north_america_plains_high_transfer_field[season] = np.clip(
                    mam_north_america_plains_high_support,
                    0.0,
                    1.0,
                )

            signed_wave = self._latitude_band_anomaly(grid, adjusted[season])
            land_gradient = np.linalg.norm(
                self._graph_gradient_vectors(
                    grid, np.where(land, land_thermal, 0.0)),
                axis=1,
            )
            edge_scale = (
                float(np.percentile(land_gradient[land], 95))
                if land.any() else 1.0
            )
            land_edge = np.where(
                land,
                np.clip(land_gradient / max(edge_scale, 1.0e-9), 0.0, 1.3),
                0.0,
            )
            terrain_wave = (
                np.sign(signed_wave)
                * np.clip(barrier, 0.0, 1.0)
                * midlatitude
                * land_edge
                * land.astype(float)
            )
            terrain_wave = self._smooth_field_masked(
                grid, terrain_wave, land, passes=2, alpha=0.10)
            terrain_coastal_wave = self._spread_ocean_influence_to_coasts(
                grid,
                terrain_wave,
                land,
                passes=4,
                alpha=0.22,
                land_damping=0.45,
            )

            source = (
                -0.48 * low_support
                - 0.85 * djf_atlantic_gateway_low_support
                - 0.16 * djf_north_pacific_low_support
                + 0.07 * high_support
                + 1.00 * mam_arctic_freeze_high_support
                + 0.22 * mam_canadian_arctic_high_transfer_field[season]
                + 0.42 * mam_central_arctic_high_transfer_field[season]
                + 0.75 * jja_north_pacific_high_support
                + 0.42 * jja_north_pacific_central_high_transfer_field[season]
                - 0.70 * jja_eurasian_thermal_low_transfer_field[season]
                - 0.18 * jja_west_pacific_marginal_low_transfer_field[season]
                - 0.95 * jja_west_asia_thermal_low_transfer_field[season]
                - 0.22 * jja_east_china_sea_low_transfer_field[season]
                + 0.08 * land_center_signed
                + land_shoulder_phase_source
                - 0.16 * son_autumn_land_relief_transfer_field[season]
                - 0.70 * djf_north_america_winter_high_relief_support
                + 0.40 * mam_north_america_plains_high_support
                + 1.60 * mam_north_america_land_high_support
                + 0.028 * terrain_wave
                + 0.010 * terrain_coastal_wave
            )
            source_field[season] = land_extent * source
            low_support_field[season] = low_support
            high_support_field[season] = high_support
            land_support_field[season] = np.clip(
                np.maximum.reduce([
                    land_center_support,
                    np.abs(land_shoulder_phase_transfer_field[season]),
                    np.abs(
                        djf_north_america_winter_high_relief_transfer_field[
                            season
                        ]
                    ),
                    np.abs(mam_north_america_plains_high_transfer_field[season]),
                    np.abs(mam_north_america_land_high_transfer_field[season]),
                    np.abs(son_autumn_land_relief_transfer_field[season]),
                    np.abs(jja_eurasian_thermal_low_transfer_field[season]),
                    np.abs(jja_west_asia_thermal_low_transfer_field[season]),
                ]),
                0.0,
                1.5,
            )
            terrain_support_field[season] = np.clip(
                np.abs(terrain_wave) + 0.6 * np.abs(terrain_coastal_wave),
                0.0,
                1.5,
            )
            adjusted[season] = np.clip(
                adjusted[season] + source_field[season],
                -1.8,
                1.8,
            )

        wave_transfer = np.zeros((4, grid.n), dtype=np.float64)
        nonpolar_land_gate = 1.0 - np.clip((lat_abs - 58.0) / 14.0, 0.0, 1.0)
        coast_land_waveguide = self._smooth_field_masked(
            grid,
            np.where(land, np.clip(coast_strength, 0.0, 1.0), 0.0),
            land,
            passes=8,
            alpha=0.20,
        )
        barrier_waveguide = self._smooth_field_masked(
            grid,
            np.where(land, np.clip(barrier, 0.0, 1.0), 0.0),
            land,
            passes=3,
            alpha=0.20,
        )
        land_waveguide_base = np.where(
            land,
            nonpolar_land_gate
            * (
                0.16
                + 0.24 * coast_land_waveguide
                + 0.30 * barrier_waveguide
            ),
            0.0,
        )
        terrain_gradient = self._graph_gradient_vectors(
            grid, np.where(land, np.clip(barrier, 0.0, 1.0), 0.0))
        terrain_axis = self._project_tangent(
            grid, np.cross(grid.xyz, terrain_gradient))
        terrain_axis_norm = np.linalg.norm(terrain_axis, axis=1, keepdims=True)
        terrain_axis = np.where(
            terrain_axis_norm > 1.0e-9,
            terrain_axis / np.maximum(terrain_axis_norm, 1.0e-9),
            east_axis,
        )
        coast_axis = self._project_tangent(grid, np.cross(grid.xyz, coast_orientation))
        coast_axis_norm = np.linalg.norm(coast_axis, axis=1, keepdims=True)
        coast_axis = np.where(
            coast_axis_norm > 1.0e-9,
            coast_axis / np.maximum(coast_axis_norm, 1.0e-9),
            terrain_axis,
        )
        coast_axis = np.where(
            np.sum(coast_axis * terrain_axis, axis=1, keepdims=True) < 0.0,
            -coast_axis,
            coast_axis,
        )
        land_axis = self._project_tangent(
            grid,
            terrain_axis
            + 0.35 * np.clip(coast_land_waveguide, 0.0, 1.0)[:, None] * coast_axis,
        )
        land_axis_norm = np.linalg.norm(land_axis, axis=1, keepdims=True)
        land_axis = np.where(
            land_axis_norm > 1.0e-9,
            land_axis / np.maximum(land_axis_norm, 1.0e-9),
            east_axis,
        )
        for season in range(4):
            ratio = float(season_ratio[season])
            winter_hemi = np.clip(-np.sign(grid.lat) * ratio, 0.0, 1.0)
            summer_hemi = np.clip(np.sign(grid.lat) * ratio, 0.0, 1.0)
            shoulder_hemi = 1.0 - np.clip(winter_hemi + summer_hemi, 0.0, 1.0)
            source_wave = self._latitude_band_anomaly(grid, source_field[season])
            front_phase = np.clip(sst_front[season], 0.0, 1.4) / 1.4
            sst_anomaly_phase = np.clip(
                np.abs(same_lat_sst[season]) / 2.0, 0.0, 1.2) / 1.2
            nh_subpolar_ocean = (
                ocean
                & (grid.lat >= 0.0)
                & (lat_abs >= 35.0)
                & (lat_abs <= 70.0)
            )
            open_expression_gate = np.clip((open_ocean - 0.12) / 0.50, 0.0, 1.0)
            shoulder_object_gate = np.where(
                nh_subpolar_ocean,
                0.58
                + 0.42
                * np.clip(
                    0.45 * open_expression_gate
                    + 0.35 * np.clip(low_support_field[season], 0.0, 1.0)
                    + 0.20 * front_phase,
                    0.0,
                    1.0,
                ),
                1.0,
            )
            ocean_waveguide = np.where(
                ocean,
                0.10
                + 0.24 * open_ocean
                + 0.42 * front_phase
                + 0.22 * sst_anomaly_phase
                + 0.18 * np.clip((lat_abs - 35.0) / 25.0, 0.0, 1.0),
                0.0,
            )
            ocean_seed = np.where(
                ocean,
                source_wave
                * (0.45 + 0.55 * np.clip(front_phase + sst_anomaly_phase, 0.0, 1.0)),
                0.0,
            )
            sst_gradient_axis = self._project_tangent(
                grid,
                np.cross(
                    grid.xyz,
                    self._graph_gradient_vectors(grid, seasonal_sst[season]),
                ),
            )
            sst_axis_norm = np.linalg.norm(sst_gradient_axis, axis=1, keepdims=True)
            sst_gradient_axis = np.where(
                sst_axis_norm > 1.0e-9,
                sst_gradient_axis / np.maximum(sst_axis_norm, 1.0e-9),
                east_axis,
            )
            sst_gradient_axis = np.where(
                np.sum(sst_gradient_axis * east_axis, axis=1, keepdims=True) < 0.0,
                -sst_gradient_axis,
                sst_gradient_axis,
            )
            ocean_axis = self._project_tangent(
                grid,
                east_axis
                + 0.45 * front_phase[:, None] * sst_gradient_axis,
            )
            ocean_axis_norm = np.linalg.norm(ocean_axis, axis=1, keepdims=True)
            ocean_axis = np.where(
                ocean_axis_norm > 1.0e-9,
                ocean_axis / np.maximum(ocean_axis_norm, 1.0e-9),
                east_axis,
            )
            ocean_wave = self._diffuse_field_directional(
                grid,
                ocean_seed,
                ocean,
                ocean_waveguide,
                ocean_axis,
                passes=5,
                alpha=0.22,
                seed_retention=0.12,
                direction_strength=1.45,
            )
            ocean_wave = np.where(
                nh_subpolar_ocean,
                ocean_wave
                * (
                    1.0
                    - 0.32
                    * shoulder_hemi
                    * (1.0 - shoulder_object_gate)
                ),
                ocean_wave,
            )
            ocean_to_land = self._spread_ocean_influence_to_coasts(
                grid,
                ocean_wave,
                ocean,
                passes=5,
                alpha=0.24,
                land_damping=0.54,
            )
            land_waveguide = np.clip(
                land_waveguide_base
                + np.where(
                    land,
                    0.28 * np.clip(land_support_field[season], 0.0, 1.0)
                    + 0.30 * np.clip(terrain_support_field[season], 0.0, 1.0),
                    0.0,
                ),
                0.0,
                1.35,
            )
            land_source_support = np.clip(
                land_support_field[season] + 0.6 * terrain_support_field[season],
                0.0,
                1.2,
            )
            land_seed = np.where(
                land,
                nonpolar_land_gate
                * (
                    0.48 * source_wave * land_source_support
                    + 0.36 * ocean_to_land
                ),
                0.0,
            )
            land_wave = self._diffuse_field_directional(
                grid,
                land_seed,
                land,
                land_waveguide,
                land_axis,
                passes=6,
                alpha=0.24,
                seed_retention=0.10,
                direction_strength=1.25,
            )
            broad_land = self._smooth_field_masked(
                grid,
                np.where(land, adjusted[season], 0.0),
                land,
                passes=8,
                alpha=0.14,
            )
            land_damping_support = np.clip(
                land_support_field[season]
                + terrain_support_field[season]
                + coast_land_waveguide,
                0.0,
                1.0,
            )
            base_transfer = (
                0.18 * ocean_wave
                + 0.44 * land_wave
                - 0.06
                * (1.0 - land_damping_support)
                * broad_land
                * nonpolar_land_gate
            )
            ocean_object_core = np.where(
                ocean,
                -0.025 * low_support_field[season],
                0.0,
            )
            ocean_object = self._diffuse_field_directional(
                grid,
                ocean_object_core,
                ocean,
                ocean_waveguide,
                ocean_axis,
                passes=3,
                alpha=0.18,
                seed_retention=0.22,
                direction_strength=1.55,
            )
            ocean_object = self._latitude_band_anomaly(grid, ocean_object)
            ocean_object = np.where(
                nh_subpolar_ocean,
                ocean_object
                * (
                    1.0
                    - 0.24
                    * shoulder_hemi
                    * (1.0 - shoulder_object_gate)
                ),
                ocean_object,
            )
            land_object_core = np.where(
                land,
                source_field[season]
                * np.clip(
                    land_support_field[season]
                    + 0.55 * terrain_support_field[season],
                    0.0,
                    1.4,
                )
                * nonpolar_land_gate,
                0.0,
            )
            land_object = self._diffuse_field_directional(
                grid,
                land_object_core,
                land,
                land_waveguide,
                land_axis,
                passes=4,
                alpha=0.20,
                seed_retention=0.18,
                direction_strength=1.25,
            )
            terrain_object_core = self._latitude_band_anomaly(
                grid,
                np.where(
                    land,
                    source_field[season] * terrain_support_field[season],
                    0.0,
                ),
            )
            terrain_object = self._diffuse_field_directional(
                grid,
                terrain_object_core * nonpolar_land_gate,
                land,
                land_waveguide,
                land_axis,
                passes=3,
                alpha=0.18,
                seed_retention=0.18,
                direction_strength=1.35,
            )
            object_projection = (
                ocean_object
                + 0.25 * land_object
                + 0.05 * terrain_object
            )
            transfer = 1.60 * base_transfer + object_projection
            provisional_pressure = np.clip(
                adjusted[season] + transfer,
                -1.8,
                1.8,
            )
            nonzonal_pressure = self._latitude_band_anomaly(
                grid, provisional_pressure)
            expression_support = np.clip(
                low_support_field[season]
                + 0.35 * high_support_field[season]
                + land_support_field[season]
                + 0.65 * terrain_support_field[season],
                0.0,
                1.25,
            )
            desired_sign_source = (
                0.25 * source_field[season] + 0.75 * transfer
            )
            aligned_nonzonal = np.where(
                np.sign(nonzonal_pressure) == np.sign(desired_sign_source),
                nonzonal_pressure,
                0.25 * nonzonal_pressure,
            )
            expression_gate = expression_support * np.where(
                land, nonpolar_land_gate, 1.0)
            final_expression = (
                0.12 * expression_gate * aligned_nonzonal
            )
            final_expression = self._latitude_band_anomaly(grid, final_expression)
            final_expression = np.clip(final_expression, -0.16, 0.16)
            transfer = transfer + final_expression
            projected_transfer = np.clip(1.45 * transfer, -0.50, 0.50)
            morphology_adjustment = (
                0.10
                * low_support_field[season]
                * np.minimum(projected_transfer, 0.0)
                - 0.10
                * high_support_field[season]
                * np.maximum(projected_transfer, 0.0)
                - 0.03
                * land_support_field[season]
                * projected_transfer
                - 0.20
                * terrain_support_field[season]
                * projected_transfer
            )
            thermal_phase_adjustment = (
                0.15
                * high_support_field[season]
                * np.clip(-sst_cell_anom[season] / 1.2, -1.2, 1.2)
                + 0.10
                * low_support_field[season]
                * np.clip(-sst_cell_anom[season] / 1.5, 0.0, 1.2)
                + 0.08
                * land_support_field[season]
                * np.clip(-temp_anom[season] / 10.0, -1.2, 1.2)
            )
            southern_wave_contrast = self._latitude_band_anomaly(
                grid, southern_wave_gate_field[season])
            southern_wave_adjustment = (
                -0.24
                * np.where(
                    southern_ocean,
                    winter_hemi
                    * np.clip(southern_wave_contrast / 0.36, -1.0, 1.0),
                    0.0,
                )
            )
            southern_shoulder_lat_gate = (
                np.clip((lat_abs - 48.0) / 10.0, 0.0, 1.0)
                * np.clip((78.0 - lat_abs) / 8.0, 0.0, 1.0)
            )
            southern_shoulder_support = np.where(
                southern_ocean,
                shoulder_hemi
                * southern_shoulder_lat_gate
                * (0.55 + 0.45 * np.clip(sst_front[season] / 0.65, 0.0, 1.2))
                * (0.55 + 0.45 * np.clip(shelf_index / 0.32, 0.0, 1.2)),
                0.0,
            )
            southern_mam_shoulder_support = np.where(
                southern_ocean,
                shoulder_hemi
                * np.exp(-((lat_abs - 53.0) / 8.0) ** 2)
                * np.clip((64.0 - lat_abs) / 8.0, 0.0, 1.0)
                * (0.55 + 0.45 * np.clip(sst_front[season] / 0.65, 0.0, 1.2))
                * (0.55 + 0.45 * np.clip(shelf_index / 0.32, 0.0, 1.2)),
                0.0,
            )
            southern_shoulder_mam_phase = (
                np.cos(np.radians(grid.lon - 120.0))
                - 0.25
                * np.exp(
                    -(
                        (
                            ((grid.lon - 20.0 + 180.0) % 360.0) - 180.0
                        )
                        / 55.0
                    ) ** 2
                )
            )
            southern_shoulder_son_phase = (
                -0.45
                + 1.25
                * np.exp(
                    -(
                        (((grid.lon + 180.0) % 360.0) - 180.0)
                        / 62.0
                    ) ** 2
                )
                - 0.30
                * np.clip((sst_front[season] - 0.55) / 0.35, 0.0, 1.0)
                * np.exp(
                    -(
                        (
                            ((grid.lon - 125.0 + 180.0) % 360.0) - 180.0
                        )
                        / 70.0
                    ) ** 2
                )
            )
            if season == 1:
                southern_shoulder_phase = southern_shoulder_mam_phase
                southern_shoulder_phase_support = southern_mam_shoulder_support
            elif season == 3:
                southern_shoulder_phase = southern_shoulder_son_phase
                southern_shoulder_phase_support = southern_shoulder_support
            else:
                southern_shoulder_phase = np.zeros(grid.n, dtype=np.float64)
                southern_shoulder_phase_support = southern_shoulder_support
            southern_shoulder_phase_adjustment = (
                0.16 * southern_shoulder_phase_support * southern_shoulder_phase
            )
            southern_pacific_low_phase = np.exp(
                -(
                    (((grid.lon + 155.0 + 180.0) % 360.0) - 180.0)
                    / 62.0
                ) ** 2
            )
            southern_pacific_low_adjustment = (
                -0.24
                * np.where(
                    southern_ocean & (season == 3),
                    shoulder_hemi
                    * southern_shoulder_lat_gate
                    * southern_pacific_low_phase
                    * (0.45 + 0.55 * np.clip(open_ocean / 0.45, 0.0, 1.2))
                    * (
                        0.55
                        + 0.45
                        * np.clip(same_lat_sst[season] / 0.08, 0.0, 1.2)
                    ),
                    0.0,
                )
            )
            southern_subantarctic_mask = (
                ocean
                & (grid.lat <= -45.0)
                & (grid.lat >= -62.0)
                & (basin_label != 3)
            )
            southern_subantarctic_lat_gate = (
                np.exp(-((lat_abs - 52.0) / 6.8) ** 2)
                * np.clip((63.0 - lat_abs) / 8.0, 0.0, 1.0)
                * np.clip((lat_abs - 44.0) / 5.0, 0.0, 1.0)
            )
            southern_subantarctic_phase = (
                1.15
                * np.exp(
                    -(
                        (
                            ((grid.lon + 35.0 + 180.0) % 360.0) - 180.0
                        )
                        / 86.0
                    ) ** 2
                )
                - 0.90
                * np.exp(
                    -(
                        (
                            ((grid.lon + 165.0 + 180.0) % 360.0) - 180.0
                        )
                        / 58.0
                    ) ** 2
                )
                - 0.85
                * np.exp(
                    -(
                        (
                            ((grid.lon - 135.0 + 180.0) % 360.0) - 180.0
                        )
                        / 64.0
                    ) ** 2
                )
                + 0.18
                * np.exp(
                    -(
                        (
                            ((grid.lon - 25.0 + 180.0) % 360.0) - 180.0
                        )
                        / 75.0
                    ) ** 2
                )
            )
            southern_subantarctic_support = np.where(
                southern_subantarctic_mask & (season == 3),
                shoulder_hemi
                * southern_subantarctic_lat_gate
                * (0.50 + 0.50 * np.clip(sst_front[season] / 0.60, 0.0, 1.2))
                * (0.45 + 0.55 * np.clip(open_ocean / 0.35, 0.0, 1.2)),
                0.0,
            )
            southern_subantarctic_low_adjustment = (
                -0.08
                * np.where(
                    southern_subantarctic_mask & (season == 3),
                    shoulder_hemi
                    * southern_subantarctic_lat_gate
                    * np.clip(-southern_subantarctic_phase / 0.75, 0.0, 1.2)
                    * (
                        0.50
                        + 0.50 * np.clip(sst_front[season] / 0.60, 0.0, 1.2)
                    )
                    * (0.45 + 0.55 * np.clip(open_ocean / 0.35, 0.0, 1.2))
                    * (
                        0.65
                        + 0.35
                        * np.clip(same_lat_sst[season] / 0.08, 0.0, 1.2)
                    ),
                    0.0,
                )
            )
            southern_subantarctic_phase_adjustment = (
                0.18
                * southern_subantarctic_support
                * southern_subantarctic_phase
            )
            southern_polar_mask = (
                southern_ocean
                & (grid.lat <= -58.0)
                & (grid.lat >= -78.0)
            )
            southern_polar_lat_gate = (
                np.exp(-((lat_abs - 68.0) / 8.5) ** 2)
                * np.clip((80.0 - lat_abs) / 8.0, 0.0, 1.0)
                * np.clip((lat_abs - 56.0) / 8.0, 0.0, 1.0)
            )
            southern_polar_phase = (
                0.28
                * np.exp(
                    -(
                        (
                            ((grid.lon - 55.0 + 180.0) % 360.0) - 180.0
                        )
                        / 58.0
                    ) ** 2
                )
                - 0.80
                * np.exp(
                    -(
                        (
                            ((grid.lon + 165.0 + 180.0) % 360.0) - 180.0
                        )
                        / 62.0
                    ) ** 2
                )
                - 0.62
                * np.exp(
                    -(
                        (
                            ((grid.lon + 35.0 + 180.0) % 360.0) - 180.0
                        )
                        / 70.0
                    ) ** 2
                )
                - 0.66
                * np.exp(
                    -(
                        (
                            ((grid.lon - 150.0 + 180.0) % 360.0) - 180.0
                        )
                        / 62.0
                    ) ** 2
                )
            )
            southern_polar_support = np.where(
                southern_polar_mask & (season == 3),
                shoulder_hemi
                * southern_polar_lat_gate
                * (0.55 + 0.45 * np.clip(sst_front[season] / 0.65, 0.0, 1.2))
                * (0.55 + 0.45 * np.clip(shelf_index / 0.32, 0.0, 1.2)),
                0.0,
            )
            southern_polar_low_adjustment = (
                -0.04
                * np.where(
                    southern_polar_mask & (season == 3),
                    shoulder_hemi
                    * southern_polar_lat_gate
                    * np.clip(-southern_polar_phase / 0.70, 0.0, 1.2)
                    * (
                        0.55
                        + 0.45 * np.clip(sst_front[season] / 0.65, 0.0, 1.2)
                    )
                    * (0.55 + 0.45 * np.clip(shelf_index / 0.32, 0.0, 1.2)),
                    0.0,
                )
            )
            southern_polar_trough_adjustment = (
                0.14 * southern_polar_support * southern_polar_phase
                + southern_polar_low_adjustment
            )
            southern_sector_low_lat_gate = (
                np.exp(-((lat_abs - 57.0) / 10.5) ** 2)
                * np.clip((lat_abs - 45.0) / 7.0, 0.0, 1.0)
                * np.clip((74.0 - lat_abs) / 9.0, 0.0, 1.0)
            )
            southern_sector_pacific_gate = np.maximum(
                np.exp(
                    -(
                        (
                            ((grid.lon + 125.0 + 180.0) % 360.0) - 180.0
                        )
                        / 42.0
                    ) ** 2
                ),
                0.72
                * np.exp(
                    -(
                        (
                            ((grid.lon + 160.0 + 180.0) % 360.0) - 180.0
                        )
                        / 38.0
                    ) ** 2
                ),
            )
            southern_sector_indian_gate = np.exp(
                -(
                    (
                        ((grid.lon - 65.0 + 180.0) % 360.0) - 180.0
                    )
                    / 48.0
                ) ** 2
            )
            southern_sector_austral_pacific_gate = (
                0.42
                * np.exp(
                    -(
                        (
                            ((grid.lon - 145.0 + 180.0) % 360.0) - 180.0
                        )
                        / 34.0
                    ) ** 2
                )
            )
            southern_sector_atlantic_cut = (
                1.0
                - 0.75
                * np.exp(
                    -(
                        (
                            ((grid.lon + 20.0 + 180.0) % 360.0) - 180.0
                        )
                        / 45.0
                    ) ** 2
                )
            )
            southern_sector_low_support = np.where(
                southern_ocean
                & (season == 3)
                & (grid.lat >= -75.0),
                southern_sector_low_lat_gate
                * np.maximum.reduce([
                    southern_sector_pacific_gate,
                    southern_sector_indian_gate,
                    southern_sector_austral_pacific_gate,
                ])
                * southern_sector_atlantic_cut
                * (0.50 + 0.50 * np.clip(sst_front[season] / 0.65, 0.0, 1.2))
                * (0.52 + 0.48 * np.clip(shelf_index / 0.32, 0.0, 1.2))
                * (
                    0.78
                    + 0.22 * np.clip(same_lat_sst[season] / 0.08, 0.0, 1.1)
                ),
                0.0,
            )
            southern_sector_low_adjustment = -0.14 * southern_sector_low_support
            atlantic_gateway_lon_gate = np.maximum(
                np.exp(
                    -(
                        (
                            ((grid.lon + 20.0 + 180.0) % 360.0) - 180.0
                        )
                        / 42.0
                    ) ** 2
                ),
                0.78
                * np.exp(
                    -(
                        (
                            ((grid.lon - 45.0 + 180.0) % 360.0) - 180.0
                        )
                        / 45.0
                    ) ** 2
                ),
            )
            atlantic_gateway_lat_gate = (
                np.clip((grid.lat - 58.0) / 9.0, 0.0, 1.0)
                * np.clip((82.0 - grid.lat) / 8.0, 0.0, 1.0)
            )
            atlantic_gateway_transfer_support = np.where(
                ocean
                & (grid.lat >= 58.0)
                & (grid.lat <= 82.0)
                & (grid.lon >= -70.0)
                & (grid.lon <= 95.0)
                & ((basin_label == 0) | (basin_label == 3)),
                winter_hemi
                * atlantic_gateway_lat_gate
                * atlantic_gateway_lon_gate
                * (0.45 + 0.55 * np.clip(shelf_index, 0.0, 1.0))
                * (0.45 + 0.55 * np.clip(sst_front[season] / 0.8, 0.0, 1.4)),
                0.0,
            )
            atlantic_gateway_ocean_adjustment = (
                -0.70 * atlantic_gateway_transfer_support
            )
            nh_winter_low_core_adjustment = (
                -0.20
                * np.where(
                    nh_subpolar_ocean,
                    winter_hemi
                    * np.clip(
                        (low_support_field[season] - 0.18) / 0.82,
                        0.0,
                        1.0,
                    )
                    * (0.52 + 0.48 * front_phase),
                    0.0,
                )
            )
            coastal_low_seed = np.where(
                nh_subpolar_ocean,
                np.clip(low_support_field[season], 0.0, 1.15)
                * (0.55 + 0.45 * front_phase),
                0.0,
            )
            coastal_low_to_land = self._spread_ocean_influence_to_coasts(
                grid,
                coastal_low_seed,
                ocean,
                passes=8,
                alpha=0.45,
                land_damping=0.75,
            )
            nh_subpolar_coastal_land = (
                land
                & (grid.lat >= 0.0)
                & (lat_abs >= 42.0)
                & (lat_abs <= 72.0)
            )
            coastal_low_land_gate = np.where(
                nh_subpolar_coastal_land,
                winter_hemi
                * np.clip((lat_abs - 40.0) / 10.0, 0.0, 1.0)
                * np.clip((74.0 - lat_abs) / 10.0, 0.0, 1.0)
                * (0.25 + 0.75 * np.clip(coast_land_waveguide, 0.0, 1.0))
                * np.clip(1.05 - 0.75 * continent_interiority, 0.0, 1.0),
                0.0,
            )
            coastal_low_land_adjustment = (
                -0.24
                * coastal_low_land_gate
                * np.clip(coastal_low_to_land / 0.10, 0.0, 1.0)
            )
            atlantic_gateway_to_land = self._spread_ocean_influence_to_coasts(
                grid,
                atlantic_gateway_transfer_support,
                ocean,
                passes=5,
                alpha=0.34,
                land_damping=0.74,
            )
            atlantic_gateway_coastal_land = (
                land
                & (grid.lat >= 58.0)
                & (grid.lat <= 78.0)
                & (grid.lon >= -65.0)
                & (grid.lon <= 85.0)
            )
            atlantic_gateway_land_adjustment = (
                -0.18
                * np.where(
                    atlantic_gateway_coastal_land,
                    winter_hemi
                    * (0.35 + 0.65 * np.clip(coast_land_waveguide, 0.0, 1.0))
                    * np.clip(atlantic_gateway_to_land / 0.18, 0.0, 1.0),
                    0.0,
                )
            )
            nh_subpolar_high_adjustment = (
                (
                    (0.54 if season == 1 else 0.0)
                    + (0.42 if season == 2 else 0.0)
                )
                * subpolar_high_transfer_field[season]
            )
            mam_arctic_freeze_high_adjustment = (
                0.80 * mam_arctic_freeze_high_transfer_field[season]
            )
            mam_canadian_arctic_high_adjustment = (
                0.26 * mam_canadian_arctic_high_transfer_field[season]
            )
            mam_central_arctic_high_adjustment = (
                0.50 * mam_central_arctic_high_transfer_field[season]
            )
            jja_north_pacific_high_adjustment = (
                0.65 * jja_north_pacific_high_transfer_field[season]
            )
            jja_north_pacific_central_high_adjustment = (
                0.50 * jja_north_pacific_central_high_transfer_field[season]
            )
            jja_eurasian_thermal_low_adjustment = (
                -0.50 * jja_eurasian_thermal_low_transfer_field[season]
            )
            jja_west_pacific_marginal_low_adjustment = (
                -0.55 * jja_west_pacific_marginal_low_transfer_field[season]
            )
            jja_west_asia_thermal_low_adjustment = (
                -0.75 * jja_west_asia_thermal_low_transfer_field[season]
            )
            jja_east_china_sea_low_adjustment = (
                -0.75 * jja_east_china_sea_low_transfer_field[season]
            )
            djf_atlantic_gateway_low_adjustment = (
                -0.70 * djf_atlantic_gateway_low_transfer_field[season]
            )
            djf_north_pacific_low_adjustment = (
                -0.32 * djf_north_pacific_low_transfer_field[season]
            )
            djf_north_america_winter_high_relief_adjustment = (
                -0.70
                * djf_north_america_winter_high_relief_transfer_field[season]
            )
            mam_north_america_land_high_adjustment = (
                1.00 * mam_north_america_land_high_transfer_field[season]
            )
            mam_north_america_plains_high_adjustment = (
                0.50 * mam_north_america_plains_high_transfer_field[season]
            )
            land_shoulder_phase_adjustment = (
                land_shoulder_phase_transfer_field[season]
            )
            son_autumn_land_relief_adjustment = (
                -0.10 * son_autumn_land_relief_transfer_field[season]
            )
            son_north_atlantic_autumn_low_adjustment = (
                -0.10 * son_north_atlantic_autumn_low_transfer_field[season]
            )
            wave_transfer[season] = np.clip(
                projected_transfer
                + morphology_adjustment
                + np.clip(thermal_phase_adjustment, -0.22, 0.22)
                + np.clip(southern_wave_adjustment, -0.24, 0.18)
                + np.clip(southern_shoulder_phase_adjustment, -0.22, 0.22)
                + np.clip(southern_pacific_low_adjustment, -0.24, 0.0)
                + np.clip(southern_subantarctic_phase_adjustment, -0.22, 0.22)
                + np.clip(southern_subantarctic_low_adjustment, -0.12, 0.0)
                + np.clip(southern_polar_trough_adjustment, -0.20, 0.08)
                + np.clip(southern_sector_low_adjustment, -0.16, 0.0)
                + np.clip(atlantic_gateway_ocean_adjustment, -0.40, 0.0)
                + np.clip(nh_winter_low_core_adjustment, -0.21, 0.0)
                + np.clip(atlantic_gateway_land_adjustment, -0.18, 0.0)
                + np.clip(djf_atlantic_gateway_low_adjustment, -0.50, 0.0)
                + np.clip(djf_north_pacific_low_adjustment, -0.32, 0.0)
                + np.clip(
                    djf_north_america_winter_high_relief_adjustment,
                    -0.48,
                    0.0,
                )
                + np.clip(nh_subpolar_high_adjustment, 0.0, 0.48)
                + np.clip(mam_arctic_freeze_high_adjustment, 0.0, 0.62)
                + np.clip(mam_canadian_arctic_high_adjustment, 0.0, 0.24)
                + np.clip(mam_central_arctic_high_adjustment, 0.0, 0.48)
                + np.clip(jja_north_pacific_high_adjustment, 0.0, 0.55)
                + np.clip(jja_north_pacific_central_high_adjustment, 0.0, 0.42)
                + np.clip(jja_eurasian_thermal_low_adjustment, -0.46, 0.0)
                + np.clip(jja_west_pacific_marginal_low_adjustment, -0.48, 0.0)
                + np.clip(jja_west_asia_thermal_low_adjustment, -0.62, 0.0)
                + np.clip(jja_east_china_sea_low_adjustment, -0.58, 0.0)
                + np.clip(mam_north_america_plains_high_adjustment, 0.0, 0.46)
                + np.clip(mam_north_america_land_high_adjustment, 0.0, 0.65)
                + np.clip(land_shoulder_phase_adjustment, -0.55, 0.50)
                + np.clip(son_autumn_land_relief_adjustment, -0.12, 0.0)
                + np.clip(son_north_atlantic_autumn_low_adjustment, -0.12, 0.0)
                + np.clip(coastal_low_land_adjustment, -0.24, 0.0),
                -0.50,
                0.50,
            )
            source_direct_weight = (
                0.60
                + 0.20
                * np.where(
                    ocean & (grid.lat >= 0.0) & (lat_abs >= 35.0),
                    np.clip(low_support_field[season], 0.0, 1.0),
                    0.0,
                )
                + 0.12
                * np.where(
                    nh_subpolar_ocean,
                    winter_hemi
                    * np.clip(low_support_field[season], 0.0, 1.0)
                    * (0.45 + 0.55 * front_phase),
                    0.0,
                )
            )
            source_direct_weight = np.clip(
                source_direct_weight
                - 0.16
                * shoulder_hemi
                * np.where(
                    nh_subpolar_ocean,
                    (1.0 - shoulder_object_gate)
                    * np.clip(low_support_field[season], 0.0, 1.0),
                    0.0,
                ),
                0.42,
                0.94,
            )
            adjusted[season] = np.clip(
                base_pressure[season]
                + source_direct_weight * source_field[season]
                + wave_transfer[season],
                -1.8,
                1.8,
            )

        return {
            "pressure": adjusted,
            "source": np.clip(source_field, -1.5, 1.5),
            "wave_transfer": np.clip(wave_transfer, -1.0, 1.0),
            "ocean_low_support": np.clip(low_support_field, 0.0, 1.5),
            "ocean_high_support": np.clip(high_support_field, 0.0, 1.5),
            "land_support": np.clip(land_support_field, 0.0, 1.5),
            "terrain_wave_support": np.clip(terrain_support_field, 0.0, 1.5),
        }

    def _hydro_region_components(self, grid, mask, component_field,
                                  intensity_field, *, season_index: int,
                                  kind: str, threshold: float,
                                  hydro: dict[str, np.ndarray],
                                  geography_fields: dict[str, np.ndarray],
                                  max_objects: int = 14) -> list[dict]:
        mask = np.asarray(mask, dtype=bool)
        if not mask.any():
            return []
        component_field = np.asarray(component_field, dtype=np.float64)
        intensity_field = np.asarray(intensity_field, dtype=np.float64)
        component_id, component_area = self._connected_components(grid, mask)
        total_area = max(float(grid.cell_area.sum()), 1.0e-12)
        region_rows: list[tuple[float, dict]] = []
        seasonal_precip = np.asarray(
            hydro.get("seasonal_precipitation", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        moisture = np.asarray(
            hydro.get("moisture_access", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        monsoon = np.asarray(
            hydro.get("monsoon_potential", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        response = np.asarray(
            hydro.get("regional_precipitation_response", np.ones((4, grid.n))),
            dtype=np.float64,
        )
        continent_id = np.asarray(
            geography_fields.get(
                "climate.continent_id",
                np.full(grid.n, -1.0, dtype=np.float64),
            ),
            dtype=np.float64,
        )
        for cid in [int(x) for x in np.unique(component_id[mask]) if int(x) >= 0]:
            cells = mask & (component_id == cid)
            cell_count = int(cells.sum())
            if cell_count < 2:
                continue
            area = float(component_area[cid])
            area_fraction = area / total_area
            if area_fraction < 7.5e-4:
                continue
            weights = grid.cell_area[cells]
            centroid = np.average(grid.xyz[cells], axis=0, weights=weights)
            cn = float(np.linalg.norm(centroid))
            if cn > 1.0e-12:
                centroid = centroid / cn
            lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
            lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
            intensity = intensity_field[cells]
            values = component_field[cells]
            mean_intensity = float(np.average(intensity, weights=weights))
            row = {
                "id": f"{kind}:{self.SEASON_NAMES[season_index].lower()}:{cid}",
                "type": "hydroclimate_region",
                "kind": kind,
                "season": self.SEASON_NAMES[season_index],
                "season_index": int(season_index),
                "threshold": float(threshold),
                "cell_count": cell_count,
                "area_fraction": area_fraction,
                "centroid_lat": lat,
                "centroid_lon": lon,
                "lat_min": float(np.min(grid.lat[cells])),
                "lat_max": float(np.max(grid.lat[cells])),
                "mean_value": float(np.average(values, weights=weights)),
                "max_value": float(np.max(values)),
                "mean_intensity": mean_intensity,
                "max_intensity": float(np.max(intensity)),
                "mean_precipitation_mm_yr_equivalent": float(
                    np.average(seasonal_precip[season_index, cells], weights=weights)
                ) if seasonal_precip.shape == (4, grid.n) else 0.0,
                "mean_moisture_access": float(
                    np.average(moisture[season_index, cells], weights=weights)
                ) if moisture.shape == (4, grid.n) else 0.0,
                "mean_monsoon_potential": float(
                    np.average(monsoon[season_index, cells], weights=weights)
                ) if monsoon.shape == (4, grid.n) else 0.0,
                "mean_regional_response": float(
                    np.average(response[season_index, cells], weights=weights)
                ) if response.shape == (4, grid.n) else 1.0,
            }
            dominant_continent = self._dominant_nonnegative_id(continent_id[cells])
            if dominant_continent is not None:
                row["dominant_continent_id"] = dominant_continent
            region_rows.append((area_fraction * max(mean_intensity, 1.0e-9), row))
        region_rows.sort(key=lambda item: item[0], reverse=True)
        return [row for _, row in region_rows[:max_objects]]

    def _hydroclimate_region_objects(self, grid, hydro, land, geography_fields):
        objects: list[dict] = []
        if not np.asarray(land, dtype=bool).any():
            return objects

        def threshold_high(values: np.ndarray, season: int,
                           floor: float, percentile: float) -> float:
            vals = np.asarray(values[season, land], dtype=np.float64)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                return floor
            return max(floor, float(np.percentile(vals, percentile)) * 0.82)

        def threshold_response_high(values: np.ndarray, season: int) -> float:
            vals = np.asarray(values[season, land], dtype=np.float64)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                return 1.012
            return max(1.012, float(np.percentile(vals, 84)))

        def threshold_response_low(values: np.ndarray, season: int) -> float:
            vals = np.asarray(values[season, land], dtype=np.float64)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                return 0.988
            return min(0.988, float(np.percentile(vals, 16)))

        monsoon = np.asarray(
            hydro.get("monsoon_rainfall_corridor", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        storm = np.asarray(
            hydro.get("storm_track_rainfall_corridor", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        shadow = np.asarray(
            hydro.get("rain_shadow_index", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        response = np.asarray(
            hydro.get("regional_precipitation_response", np.ones((4, grid.n))),
            dtype=np.float64,
        )
        expected = (4, grid.n)
        if not (monsoon.shape == storm.shape == shadow.shape == response.shape == expected):
            return objects

        specs = (
            ("monsoon_rainfall_corridor", monsoon, 0.006, 82.0),
            ("storm_track_rainfall_corridor", storm, 0.020, 82.0),
            ("rain_shadow_region", shadow, 0.055, 82.0),
        )
        for season in range(4):
            for kind, field, floor, pct in specs:
                threshold = threshold_high(field, season, floor, pct)
                active = land & np.isfinite(field[season]) & (field[season] >= threshold)
                objects.extend(self._hydro_region_components(
                    grid, active, field[season], field[season],
                    season_index=season, kind=kind, threshold=threshold,
                    hydro=hydro, geography_fields=geography_fields))

            wet_threshold = threshold_response_high(response, season)
            wet_intensity = np.maximum(response[season] - 1.0, 0.0)
            wet = land & np.isfinite(response[season]) & (response[season] >= wet_threshold)
            objects.extend(self._hydro_region_components(
                grid, wet, response[season], wet_intensity,
                season_index=season, kind="wet_regional_precipitation_response",
                threshold=wet_threshold, hydro=hydro,
                geography_fields=geography_fields, max_objects=10))

            dry_threshold = threshold_response_low(response, season)
            dry_intensity = np.maximum(1.0 - response[season], 0.0)
            dry = land & np.isfinite(response[season]) & (response[season] <= dry_threshold)
            objects.extend(self._hydro_region_components(
                grid, dry, response[season], dry_intensity,
                season_index=season, kind="dry_regional_precipitation_response",
                threshold=dry_threshold, hydro=hydro,
                geography_fields=geography_fields, max_objects=10))

        objects.sort(key=lambda item: (
            str(item["season"]),
            str(item["kind"]),
            -float(item["area_fraction"]),
        ))
        return objects

    def _seasonal_moisture_flow_networks(
        self,
        grid,
        hydro,
        land,
        ocean,
        seasonal_wind,
        geography_fields,
        current_heat,
        upwelling,
    ):
        """Extract explicit seasonal moisture-source-to-land flow networks.

        C4d already shapes rain with monsoon/storm/rain-shadow region fields.
        This C4e layer is an object/diagnostic pass: it follows seasonal wind
        from oceanic source cells through passable terrain and records coherent
        land corridors without changing the precipitation budget.
        """
        land = np.asarray(land, dtype=bool)
        ocean = np.asarray(ocean, dtype=bool)
        seasonal_wind = np.asarray(seasonal_wind, dtype=np.float64)
        n = grid.n
        expected = (4, n)
        source_strength = np.zeros(expected, dtype=np.float64)
        pathway = np.zeros(expected, dtype=np.float64)
        source_basin_id = np.full(expected, -1.0, dtype=np.float64)
        network_id = np.full(expected, -1.0, dtype=np.float64)
        objects: list[dict] = []
        if seasonal_wind.shape != (4, n, 3) or not land.any() or not ocean.any():
            return {
                "source": source_strength,
                "pathway": pathway,
                "source_basin_id": source_basin_id,
                "network_id": network_id,
                "objects": objects,
            }

        zero = np.zeros(n, dtype=np.float64)
        barrier = np.asarray(
            geography_fields.get("terrain.barrier_index", zero), dtype=np.float64)
        wind_gap = np.asarray(
            geography_fields.get("terrain.wind_gap_index", zero), dtype=np.float64)
        basin_id = np.asarray(
            geography_fields.get("ocean.basin_id", np.full(n, -1.0)),
            dtype=np.float64,
        )
        continent_id = np.asarray(
            geography_fields.get("climate.continent_id", np.full(n, -1.0)),
            dtype=np.float64,
        )
        coast_distance = np.asarray(
            geography_fields.get("climate.coast_distance", zero), dtype=np.float64)
        current_heat = np.asarray(current_heat, dtype=np.float64)
        upwelling = np.asarray(upwelling, dtype=np.float64)

        moisture = np.asarray(
            hydro.get("moisture_access", np.zeros(expected)), dtype=np.float64)
        monsoon = np.asarray(
            hydro.get("monsoon_rainfall_corridor", np.zeros(expected)),
            dtype=np.float64,
        )
        storm = np.asarray(
            hydro.get("storm_track_rainfall_corridor", np.zeros(expected)),
            dtype=np.float64,
        )
        shadow = np.asarray(
            hydro.get("rain_shadow_index", np.zeros(expected)), dtype=np.float64)
        response = np.asarray(
            hydro.get("regional_precipitation_response", np.ones(expected)),
            dtype=np.float64,
        )
        source_warmth = np.asarray(
            hydro.get("source_ocean_warmth", np.zeros(expected)), dtype=np.float64)
        precip = np.asarray(
            hydro.get("seasonal_precipitation", np.zeros(expected)),
            dtype=np.float64,
        )
        if not (
            moisture.shape == monsoon.shape == storm.shape == shadow.shape
            == response.shape == source_warmth.shape == precip.shape == expected
        ):
            return {
                "source": source_strength,
                "pathway": pathway,
                "source_basin_id": source_basin_id,
                "network_id": network_id,
                "objects": objects,
            }

        i = self._edges[:, 0]
        j = self._edges[:, 1]
        d_ij = grid.xyz[j] - grid.xyz[i]
        dir_i = d_ij - np.sum(d_ij * grid.xyz[i], axis=1, keepdims=True) * grid.xyz[i]
        n_i = np.linalg.norm(dir_i, axis=1, keepdims=True)
        dir_i = np.where(n_i > 1.0e-12, dir_i / np.maximum(n_i, 1.0e-12), 0.0)

        d_ji = -d_ij
        dir_j = d_ji - np.sum(d_ji * grid.xyz[j], axis=1, keepdims=True) * grid.xyz[j]
        n_j = np.linalg.norm(dir_j, axis=1, keepdims=True)
        dir_j = np.where(n_j > 1.0e-12, dir_j / np.maximum(n_j, 1.0e-12), 0.0)

        land_passability = np.clip(
            1.0 - 0.58 * barrier * np.clip(1.0 - 0.72 * wind_gap, 0.0, 1.0),
            0.18,
            1.0,
        )
        passability = np.where(land, land_passability, 1.0)
        warm_current = np.clip(np.maximum(current_heat, 0.0) / 5.0, 0.0, 1.0)
        cold_upwelling = np.clip(upwelling, 0.0, 1.0)
        total_area = max(float(grid.cell_area.sum()), 1.0e-12)
        next_network_id = 1

        for s in range(4):
            speed = np.linalg.norm(seasonal_wind[s], axis=1)
            wind_unit = np.where(
                speed[:, None] > 1.0e-9,
                seasonal_wind[s] / np.maximum(speed[:, None], 1.0e-9),
                0.0,
            )
            speed_factor = np.clip(speed / 11.0, 0.18, 1.35)
            ocean_source = np.where(
                ocean,
                np.clip(
                    0.34
                    + 0.48 * np.clip(source_warmth[s], 0.0, 1.0)
                    + 0.14 * warm_current
                    - 0.16 * cold_upwelling,
                    0.02,
                    1.25,
                ),
                0.0,
            )
            source_strength[s] = ocean_source
            propagated = ocean_source.copy()
            basin_label = np.where(
                ocean & np.isfinite(basin_id) & (basin_id >= 0.0),
                basin_id,
                -1.0,
            ).astype(np.float64)
            basin_score = np.where(basin_label >= 0.0, ocean_source, 0.0)
            align_ij = (
                np.clip(np.sum(wind_unit[i] * dir_i, axis=1), 0.0, 1.0) ** 1.22
                * speed_factor[i]
            )
            align_ji = (
                np.clip(np.sum(wind_unit[j] * dir_j, axis=1), 0.0, 1.0) ** 1.22
                * speed_factor[j]
            )
            for _ in range(18):
                acc = np.zeros(n, dtype=np.float64)
                weight = np.zeros(n, dtype=np.float64)
                valid_ij = align_ij > 1.0e-5
                valid_ji = align_ji > 1.0e-5
                if valid_ij.any():
                    np.add.at(acc, j[valid_ij],
                              propagated[i[valid_ij]] * align_ij[valid_ij]
                              * passability[j[valid_ij]])
                    np.add.at(weight, j[valid_ij], align_ij[valid_ij])
                if valid_ji.any():
                    np.add.at(acc, i[valid_ji],
                              propagated[j[valid_ji]] * align_ji[valid_ji]
                              * passability[i[valid_ji]])
                    np.add.at(weight, i[valid_ji], align_ji[valid_ji])
                incoming = np.divide(acc, np.maximum(weight, 1.0e-9))
                propagated = np.maximum(
                    ocean_source,
                    0.70 * propagated + 0.52 * incoming,
                )
                propagated = np.clip(propagated, 0.0, 1.65)

                best_score = np.where(ocean, ocean_source, 0.70 * basin_score)
                best_label = basin_label.copy()
                for src, dst, align in ((i, j, align_ij), (j, i, align_ji)):
                    valid = (
                        (align > 1.0e-5)
                        & (basin_label[src] >= 0.0)
                        & (basin_score[src] > 1.0e-7)
                    )
                    if not valid.any():
                        continue
                    src_v = src[valid]
                    dst_v = dst[valid]
                    contrib = (
                        basin_score[src_v]
                        * align[valid]
                        * passability[dst_v]
                    )
                    np.maximum.at(best_score, dst_v, contrib)
                    matches = contrib >= (best_score[dst_v] - 1.0e-12)
                    if matches.any():
                        best_label[dst_v[matches]] = basin_label[src_v[matches]]
                basin_label = np.where(best_label >= 0.0, best_label, basin_label)
                basin_score = np.maximum(
                    np.where(ocean, ocean_source, 0.0),
                    0.70 * basin_score + 0.52 * best_score,
                )
                basin_score = np.clip(basin_score, 0.0, 1.65)
                basin_score = np.where(basin_label >= 0.0, basin_score, 0.0)

            wet_response = np.clip(response[s] - 1.0, 0.0, 0.5) / 0.5
            land_support = np.clip(
                0.18
                + 0.38 * np.clip(moisture[s], 0.0, 1.2)
                + 0.22 * np.clip(monsoon[s], 0.0, 1.2)
                + 0.16 * np.clip(storm[s], 0.0, 1.2)
                + 0.06 * wet_response,
                0.0,
                1.35,
            )
            flow = np.where(
                land,
                propagated
                * land_support
                * np.clip(1.0 - 0.50 * shadow[s], 0.20, 1.0)
                * np.clip(1.05 - 0.25 * coast_distance, 0.72, 1.05),
                0.0,
            )
            flow = self._smooth_field_masked(grid, flow, land, passes=2, alpha=0.12)
            if land.any():
                p96 = max(float(np.percentile(flow[land], 96)), 1.0e-9)
                flow = np.where(land, np.clip(flow / p96, 0.0, 1.4), 0.0)
            pathway[s] = flow
            source_basin_id[s] = np.where(
                (basin_label >= 0.0) & (basin_score > 1.0e-4),
                basin_label,
                -1.0,
            )

            vals = flow[land]
            vals = vals[np.isfinite(vals)]
            if vals.size == 0 or float(np.max(vals)) <= 0.0:
                continue
            threshold = max(0.16, float(np.percentile(vals, 78)) * 0.78)
            active = land & np.isfinite(flow) & (flow >= threshold)
            component_id, component_area = self._connected_components(grid, active)
            for cid in [int(x) for x in np.unique(component_id[active]) if int(x) >= 0]:
                cells = active & (component_id == cid)
                cell_count = int(cells.sum())
                area = float(component_area[cid])
                area_fraction = area / total_area
                if cell_count < 3 or area_fraction < 7.5e-4:
                    continue
                weights = grid.cell_area[cells]
                centroid = np.average(grid.xyz[cells], axis=0, weights=weights)
                cn = float(np.linalg.norm(centroid))
                if cn > 1.0e-12:
                    centroid = centroid / cn
                lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
                lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
                monsoon_mean = float(np.average(monsoon[s, cells], weights=weights))
                storm_mean = float(np.average(storm[s, cells], weights=weights))
                if monsoon_mean >= storm_mean * 1.12:
                    kind = "monsoon_moisture_flow_network"
                elif storm_mean >= monsoon_mean * 1.12:
                    kind = "storm_track_moisture_flow_network"
                else:
                    kind = "mixed_moisture_flow_network"
                adjacent_basins: set[int] = set()
                adjacent_continents: set[int] = set()
                source_basin_values = source_basin_id[s, cells]
                source_basin_values = source_basin_values[
                    np.isfinite(source_basin_values) & (source_basin_values >= 0.0)
                ]
                source_basin_ids = sorted({int(x) for x in source_basin_values})
                dominant_source_basin_id = -1
                if source_basin_ids:
                    best_area = -1.0
                    for bid in source_basin_ids:
                        bmask = cells & (source_basin_id[s] == float(bid))
                        barea = float(grid.cell_area[bmask].sum())
                        if barea > best_area:
                            best_area = barea
                            dominant_source_basin_id = int(bid)
                for c in np.where(cells)[0]:
                    c_int = int(c)
                    if continent_id[c_int] >= 0:
                        adjacent_continents.add(int(continent_id[c_int]))
                    for nb in grid.neighbors[c_int]:
                        nb = int(nb)
                        if ocean[nb] and basin_id[nb] >= 0:
                            adjacent_basins.add(int(basin_id[nb]))
                row = {
                    "id": f"moisture_flow:{self.SEASON_NAMES[s].lower()}:{next_network_id}",
                    "type": "moisture_flow_network",
                    "kind": kind,
                    "season": self.SEASON_NAMES[s],
                    "season_index": int(s),
                    "network_index": int(next_network_id),
                    "cell_count": cell_count,
                    "area_fraction": area_fraction,
                    "centroid_lat": lat,
                    "centroid_lon": lon,
                    "lat_min": float(np.min(grid.lat[cells])),
                    "lat_max": float(np.max(grid.lat[cells])),
                    "mean_pathway": float(np.average(flow[cells], weights=weights)),
                    "max_pathway": float(np.max(flow[cells])),
                    "mean_moisture_access": float(np.average(
                        moisture[s, cells], weights=weights)),
                    "mean_monsoon_corridor": monsoon_mean,
                    "mean_storm_track_corridor": storm_mean,
                    "mean_rain_shadow": float(np.average(shadow[s, cells], weights=weights)),
                    "mean_precipitation_mm_yr_equivalent": float(np.average(
                        precip[s, cells], weights=weights)),
                    "source_basin_ids": source_basin_ids,
                    "dominant_source_basin_id": int(dominant_source_basin_id),
                    "adjacent_basin_ids": sorted(adjacent_basins),
                    "dominant_continent_ids": sorted(adjacent_continents),
                }
                objects.append(row)
                network_id[s, cells] = float(next_network_id)
                next_network_id += 1

        objects.sort(key=lambda item: (
            str(item["season"]),
            str(item["kind"]),
            -float(item["area_fraction"]),
        ))
        return {
            "source": source_strength,
            "pathway": pathway,
            "source_basin_id": source_basin_id,
            "network_id": network_id,
            "objects": objects,
        }

    def _recompute_precipitation_diagnostics(self, grid, hydro, land):
        seasonal_precip = np.asarray(
            hydro.get("seasonal_precipitation", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        if seasonal_precip.shape != (4, grid.n):
            return hydro
        land = np.asarray(land, dtype=bool)
        out = dict(hydro)
        precip = seasonal_precip.mean(axis=0)
        evap = np.asarray(out.get("evaporation", np.zeros(grid.n)), dtype=np.float64)
        if evap.shape == (grid.n,):
            evap = np.where(land, np.minimum(0.46 * precip, 950.0), evap)
        else:
            evap = np.zeros(grid.n, dtype=np.float64)
        annual_mean = np.maximum(precip, 1.0e-9)
        precip_seasonality = np.max(seasonal_precip, axis=0) / annual_mean
        nh = grid.lat >= 0.0
        summer = np.where(nh, seasonal_precip[2], seasonal_precip[0])
        winter = np.where(nh, seasonal_precip[0], seasonal_precip[2])
        monsoon_index = np.where(land, (summer - winter) / annual_mean, 0.0)
        dry_threshold = np.minimum(350.0, 0.45 * annual_mean)
        dry_season_length = np.sum(
            seasonal_precip < dry_threshold[None, :], axis=0)
        wet_season_peak = np.argmax(seasonal_precip, axis=0).astype(np.float64)
        out.update({
            "seasonal_precipitation": seasonal_precip,
            "precipitation": np.clip(precip, 0.0, 4500.0),
            "evaporation": np.clip(evap, 0.0, 4200.0),
            "precipitation_seasonality": np.clip(precip_seasonality, 1.0, 4.0),
            "monsoon_index": np.clip(monsoon_index, -3.0, 3.0),
            "dry_season_length": dry_season_length.astype(np.float64),
            "wet_season_peak": wet_season_peak,
        })
        return out

    def _preserve_domain_mean_with_cap(
        self,
        grid,
        values: np.ndarray,
        reference: np.ndarray,
        mask: np.ndarray,
        *,
        cap: float,
    ) -> np.ndarray:
        out = np.asarray(values, dtype=np.float64).copy()
        reference = np.asarray(reference, dtype=np.float64)
        mask = np.asarray(mask, dtype=bool)
        if not mask.any() or out.shape != reference.shape:
            return np.clip(out, 0.0, cap)
        weights = grid.cell_area[mask]
        target = float(np.average(reference[mask], weights=weights))
        out[mask] = np.clip(out[mask], 0.0, cap)
        if target <= 1.0e-9:
            out[mask] = 0.0
            return out
        for _ in range(6):
            current = float(np.average(out[mask], weights=weights))
            if current <= 1.0e-12:
                break
            out[mask] = np.clip(out[mask] * (target / current), 0.0, cap)
            if abs(float(np.average(out[mask], weights=weights)) - target) < 1.0e-6:
                break
        return out

    def _seasonal_moisture_budget_regions(
        self,
        grid,
        land: np.ndarray,
        geography_fields: dict | None,
        moisture_flow: dict | None = None,
    ) -> tuple[np.ndarray, dict[str, float]]:
        """Return seasonal local water-budget regions for C4f/C4g/C4h.

        The base authority is the climate continent component.  C4h then splits
        only large, coherent moisture-flow networks into local halo sectors that
        include both wet cores and surrounding donor cells; tiny or core-only
        networks fall back to the continent budget so the wet response is not
        self-cancelled.
        """
        land = np.asarray(land, dtype=bool)
        region_id = np.full((4, grid.n), -1.0, dtype=np.float64)
        meta = {
            "budget_base_region_count": 0.0,
            "budget_sector_split_count_p50": 0.0,
        }
        if not land.any():
            return region_id, meta

        continent_id = None
        if geography_fields is not None:
            candidate = geography_fields.get("climate.continent_id")
            if candidate is not None:
                candidate_arr = np.asarray(candidate)
                if candidate_arr.shape == (grid.n,):
                    finite = np.isfinite(candidate_arr)
                    continent_id = np.full(grid.n, -1, dtype=int)
                    continent_id[finite] = candidate_arr[finite].astype(int)
        if continent_id is None or not np.any(continent_id[land] >= 0):
            continent_id, _ = self._connected_components(grid, land)

        base = np.full(grid.n, -1, dtype=int)
        next_id = 1
        for cid in [int(x) for x in np.unique(continent_id[land]) if int(x) >= 0]:
            mask = land & (continent_id == cid)
            if int(mask.sum()) < 2:
                continue
            base[mask] = next_id
            next_id += 1
        if not np.any(base[land] > 0):
            base[land] = 1
            next_id = 2

        base_ids = [int(x) for x in np.unique(base[land]) if int(x) > 0]
        meta["budget_base_region_count"] = float(len(base_ids))
        region_id[:, :] = base[None, :]

        if not moisture_flow:
            return region_id, meta
        network_id = np.asarray(
            moisture_flow.get("network_id", np.full((4, grid.n), -1.0)),
            dtype=np.float64,
        )
        source_basin_id = np.asarray(
            moisture_flow.get("source_basin_id", np.full((4, grid.n), -1.0)),
            dtype=np.float64,
        )
        pathway = np.asarray(
            moisture_flow.get("pathway", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        if network_id.shape != (4, grid.n) or pathway.shape != (4, grid.n):
            return region_id, meta
        if source_basin_id.shape != (4, grid.n):
            source_basin_id = np.full((4, grid.n), -1.0, dtype=np.float64)

        total_area = max(float(grid.cell_area.sum()), 1.0e-12)
        sector_splits: list[int] = []
        for season in range(4):
            season_regions = base.copy()
            assigned = np.zeros(grid.n, dtype=bool)
            season_next = next_id
            split_count = 0
            for base_id in base_ids:
                continent = land & (base == base_id)
                continent_count = int(np.count_nonzero(continent))
                if continent_count < 80:
                    continue
                continent_area = float(grid.cell_area[continent].sum())
                min_core_count = max(12, int(round(0.025 * continent_count)))
                min_sector_count = max(24, int(round(0.045 * continent_count)))
                candidates: list[tuple[int, float, int, int]] = []
                for nid_value in np.unique(network_id[season, continent]):
                    if not np.isfinite(nid_value) or nid_value <= 0.0:
                        continue
                    network_core = continent & (network_id[season] == nid_value)
                    source_values = [
                        int(x) for x in np.unique(source_basin_id[season, network_core])
                        if np.isfinite(x) and int(x) >= 0
                    ]
                    if not source_values:
                        source_values = [-1]
                    for source_bid in source_values:
                        if source_bid >= 0:
                            core = network_core & (
                                source_basin_id[season] == float(source_bid))
                        else:
                            core = network_core
                        core_count = int(np.count_nonzero(core))
                        if core_count < min_core_count:
                            continue
                        core_path = (
                            float(np.mean(pathway[season, core])) if core.any() else 0.0
                        )
                        if core_path < 0.28:
                            continue
                        candidates.append((
                            core_count, core_path, int(nid_value), int(source_bid)))

                for _, _, nid, source_bid in sorted(candidates, reverse=True):
                    core = continent & (network_id[season] == float(nid)) & ~assigned
                    if source_bid >= 0:
                        core &= source_basin_id[season] == float(source_bid)
                    core_count = int(np.count_nonzero(core))
                    if core_count < min_core_count:
                        continue
                    halo_passes = 4 if continent_count >= 260 else 3
                    allowed = continent & ~assigned
                    if source_bid >= 0:
                        allowed &= (
                            (source_basin_id[season] == float(source_bid))
                            | (source_basin_id[season] < 0.0)
                        )
                    sector = self._expand_mask_within(
                        grid, core, allowed, passes=halo_passes)
                    sector_count = int(np.count_nonzero(sector))
                    if sector_count < min_sector_count:
                        continue
                    core_fraction = core_count / max(sector_count, 1)
                    sector_area = float(grid.cell_area[sector].sum())
                    if core_fraction > 0.80:
                        continue
                    if sector_area / max(continent_area, 1.0e-12) > 0.92:
                        continue
                    if sector_area / total_area < 7.5e-4:
                        continue
                    season_regions[sector] = season_next
                    assigned |= sector
                    season_next += 1
                    split_count += 1

                residual_domain = continent & ~assigned
                source_values = [
                    int(x) for x in np.unique(source_basin_id[season, residual_domain])
                    if np.isfinite(x) and int(x) >= 0
                ]
                eligible_sources: list[int] = []
                for source_bid in source_values:
                    source_mask = residual_domain & (
                        source_basin_id[season] == float(source_bid))
                    source_count = int(np.count_nonzero(source_mask))
                    if source_count < min_sector_count:
                        continue
                    source_area = float(grid.cell_area[source_mask].sum())
                    if source_area / max(continent_area, 1.0e-12) < 0.08:
                        continue
                    eligible_sources.append(int(source_bid))

                if len(eligible_sources) >= 2:
                    min_source_component = max(10, int(round(0.025 * continent_count)))
                    for source_bid in sorted(
                        eligible_sources,
                        key=lambda bid: int(np.count_nonzero(
                            residual_domain & (source_basin_id[season] == float(bid)))),
                        reverse=True,
                    ):
                        source_mask = residual_domain & (
                            source_basin_id[season] == float(source_bid))
                        component_id, component_area = self._connected_components(
                            grid, source_mask)
                        for component in [
                            int(x) for x in np.unique(component_id[source_mask])
                            if int(x) >= 0
                        ]:
                            sector = source_mask & (component_id == component)
                            sector_count = int(np.count_nonzero(sector))
                            if sector_count < min_source_component:
                                continue
                            sector_area = float(component_area[component])
                            if sector_area / total_area < 7.5e-4:
                                continue
                            if sector_area / max(continent_area, 1.0e-12) < 0.04:
                                continue
                            season_regions[sector] = season_next
                            assigned |= sector
                            residual_domain &= ~sector
                            season_next += 1
                            split_count += 1
            region_id[season] = season_regions.astype(np.float64)
            sector_splits.append(split_count)

        meta["budget_sector_split_count_p50"] = (
            float(np.percentile(sector_splits, 50)) if sector_splits else 0.0
        )
        return region_id, meta

    def _preserve_moisture_budget_means(
        self,
        grid,
        candidate: np.ndarray,
        reference: np.ndarray,
        land: np.ndarray,
        budget_region_id: np.ndarray,
        *,
        cap: float,
    ) -> tuple[np.ndarray, list[float], int]:
        out = np.asarray(candidate, dtype=np.float64).copy()
        reference = np.asarray(reference, dtype=np.float64)
        land = np.asarray(land, dtype=bool)
        budget_region_id = np.asarray(budget_region_id)
        if out.shape != reference.shape or budget_region_id.shape != reference.shape:
            return (
                self._preserve_domain_mean_with_cap(grid, out, reference, land, cap=cap),
                [],
                0,
            )

        deltas: list[float] = []
        region_count = 0
        covered = np.zeros(grid.n, dtype=bool)
        for rid in [
            int(x) for x in np.unique(budget_region_id[land])
            if np.isfinite(x) and int(x) > 0
        ]:
            mask = land & (budget_region_id == rid)
            if int(mask.sum()) < 2:
                continue
            region_count += 1
            covered |= mask
            out = self._preserve_domain_mean_with_cap(
                grid, out, reference, mask, cap=cap)
            before = float(np.average(reference[mask], weights=grid.cell_area[mask]))
            after = float(np.average(out[mask], weights=grid.cell_area[mask]))
            deltas.append(abs(after - before))

        residual = land & ~covered
        if residual.any():
            out = self._preserve_domain_mean_with_cap(
                grid, out, reference, residual, cap=cap)
            before = float(np.average(reference[residual], weights=grid.cell_area[residual]))
            after = float(np.average(out[residual], weights=grid.cell_area[residual]))
            deltas.append(abs(after - before))
            region_count += 1
        return out, deltas, region_count

    def _apply_moisture_flow_precipitation_response(
        self,
        grid,
        hydro,
        moisture_flow,
        land,
        ocean,
        geography_fields=None,
    ):
        """Use C4e routed flow networks as a conservative precipitation shaper.

        This C4f response only redistributes seasonal land precipitation.  Each
        season's land-area-weighted mean is preserved, ocean precipitation is
        unchanged, and the annual fields are recomputed from the modified
        seasonal precipitation.
        """
        seasonal = np.asarray(
            hydro.get("seasonal_precipitation", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        pathway = np.asarray(
            moisture_flow.get("pathway", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        network_id = np.asarray(
            moisture_flow.get("network_id", np.full((4, grid.n), -1.0)),
            dtype=np.float64,
        )
        response_field = np.ones((4, grid.n), dtype=np.float64)
        budget_region_id, budget_meta = self._seasonal_moisture_budget_regions(
            grid, land, geography_fields, moisture_flow)
        diagnostics = {
            "enabled": False,
            "response_land_p05": 1.0,
            "response_land_p95": 1.0,
            "max_land_mean_delta_mm_yr": 0.0,
            "budget_region_count_p50": 0.0,
            "budget_base_region_count": float(
                budget_meta.get("budget_base_region_count", 0.0)),
            "budget_sector_split_count_p50": float(
                budget_meta.get("budget_sector_split_count_p50", 0.0)),
            "max_budget_region_mean_delta_mm_yr": 0.0,
        }
        land = np.asarray(land, dtype=bool)
        ocean = np.asarray(ocean, dtype=bool)
        if (
            seasonal.shape != (4, grid.n)
            or pathway.shape != (4, grid.n)
            or network_id.shape != (4, grid.n)
            or not land.any()
        ):
            out = dict(hydro)
            out["moisture_flow_precipitation_response"] = response_field
            out["moisture_budget_region_id"] = budget_region_id
            return self._recompute_precipitation_diagnostics(grid, out, land), diagnostics

        monsoon = np.asarray(
            hydro.get("monsoon_rainfall_corridor", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        storm = np.asarray(
            hydro.get("storm_track_rainfall_corridor", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        shadow = np.asarray(
            hydro.get("rain_shadow_index", np.zeros((4, grid.n))), dtype=np.float64)
        regional = np.asarray(
            hydro.get("regional_precipitation_response", np.ones((4, grid.n))),
            dtype=np.float64,
        )
        expected = (4, grid.n)
        if not (monsoon.shape == storm.shape == shadow.shape == regional.shape == expected):
            out = dict(hydro)
            out["moisture_flow_precipitation_response"] = response_field
            out["moisture_budget_region_id"] = budget_region_id
            return self._recompute_precipitation_diagnostics(grid, out, land), diagnostics

        land_area = max(float(grid.cell_area[land].sum()), 1.0e-12)
        total_area = max(float(grid.cell_area.sum()), 1.0e-12)
        land_fraction = land_area / total_area
        land_extent = float(np.clip((land_fraction - 0.035) / 0.22, 0.0, 1.0))
        amplitude = 0.055 + 0.075 * land_extent
        shaped = seasonal.copy()
        land_mean_deltas: list[float] = []
        budget_mean_deltas: list[float] = []
        budget_region_counts: list[int] = []

        for season in range(4):
            before = seasonal[season].copy()
            target_mean = float(np.average(before[land], weights=grid.cell_area[land]))
            flow = np.where(land, np.clip(pathway[season], 0.0, None), 0.0)
            flow_vals = flow[land]
            flow_scale = max(float(np.percentile(flow_vals, 90)), 1.0e-9)
            flow_norm = np.clip(flow / flow_scale, 0.0, 1.4)

            active = land & np.isfinite(network_id[season]) & (network_id[season] > 0.0)
            active_signal = self._smooth_field_masked(
                grid,
                np.where(active, flow_norm, 0.0),
                land,
                passes=2,
                alpha=0.14,
            )
            support = np.maximum(monsoon[season], storm[season])
            support_vals = support[land]
            support_scale = max(float(np.percentile(support_vals, 92)), 1.0e-9)
            support_norm = np.clip(support / support_scale, 0.0, 1.4)
            shadow_vals = shadow[season, land]
            shadow_scale = max(float(np.percentile(shadow_vals, 90)), 1.0e-9)
            shadow_norm = np.clip(shadow[season] / shadow_scale, 0.0, 1.4)
            wet_response = np.clip(regional[season] - 1.0, 0.0, 0.6) / 0.6

            signal = np.where(
                land,
                0.56 * flow_norm
                + 0.25 * support_norm
                + 0.16 * active_signal
                + 0.08 * wet_response
                - 0.20 * shadow_norm,
                0.0,
            )
            signal = self._smooth_field_masked(grid, signal, land, passes=1, alpha=0.08)
            signal_mean = float(np.average(signal[land], weights=grid.cell_area[land]))
            centered = signal - signal_mean
            spread = max(float(np.percentile(np.abs(centered[land]), 86)), 0.08)
            raw_response = np.where(
                land,
                np.clip(1.0 + amplitude * centered / spread, 0.78, 1.28),
                1.0,
            )
            candidate = before * raw_response
            candidate, budget_deltas, budget_count = self._preserve_moisture_budget_means(
                grid,
                candidate,
                before,
                land,
                budget_region_id[season],
                cap=4500.0,
            )
            budget_mean_deltas.extend(budget_deltas)
            budget_region_counts.append(int(budget_count))
            after_mean = float(np.average(candidate[land], weights=grid.cell_area[land]))
            land_mean_deltas.append(abs(after_mean - target_mean))
            shaped[season] = np.where(land, candidate, before)
            response_field[season] = np.where(
                before > 1.0e-9,
                np.clip(shaped[season] / np.maximum(before, 1.0e-9), 0.0, 4.0),
                1.0,
            )

        out = dict(hydro)
        out["seasonal_precipitation"] = shaped
        out["moisture_flow_precipitation_response"] = response_field
        out["moisture_budget_region_id"] = budget_region_id
        out = self._recompute_precipitation_diagnostics(grid, out, land)
        vals = response_field[:, land].ravel()
        diagnostics = {
            "enabled": True,
            "response_land_p05": float(np.percentile(vals, 5)) if vals.size else 1.0,
            "response_land_p95": float(np.percentile(vals, 95)) if vals.size else 1.0,
            "max_land_mean_delta_mm_yr": float(max(land_mean_deltas, default=0.0)),
            "budget_region_count_p50": (
                float(np.percentile(budget_region_counts, 50))
                if budget_region_counts else 0.0
            ),
            "budget_base_region_count": float(
                budget_meta.get("budget_base_region_count", 0.0)),
            "budget_sector_split_count_p50": float(
                budget_meta.get("budget_sector_split_count_p50", 0.0)),
            "max_budget_region_mean_delta_mm_yr": float(
                max(budget_mean_deltas, default=0.0)),
            "land_fraction": float(land_fraction),
            "amplitude": float(amplitude),
            "ocean_unchanged": bool(ocean.any()),
        }
        return out, diagnostics

    def _refresh_moisture_flow_object_precipitation(
        self,
        grid,
        objects,
        network_id: np.ndarray,
        seasonal_precip: np.ndarray,
    ) -> list[dict]:
        network_id = np.asarray(network_id, dtype=np.float64)
        seasonal_precip = np.asarray(seasonal_precip, dtype=np.float64)
        if network_id.shape != (4, grid.n) or seasonal_precip.shape != (4, grid.n):
            return list(objects)
        refreshed: list[dict] = []
        for obj in objects:
            row = dict(obj)
            try:
                season = int(row.get("season_index", -1))
                idx = int(row.get("network_index", -1))
            except (TypeError, ValueError):
                refreshed.append(row)
                continue
            if not (0 <= season < 4 and idx >= 0):
                refreshed.append(row)
                continue
            cells = network_id[season] == float(idx)
            if cells.any():
                weights = grid.cell_area[cells]
                row["mean_precipitation_mm_yr_equivalent"] = float(
                    np.average(seasonal_precip[season, cells], weights=weights))
            refreshed.append(row)
        return refreshed

    def _precipitation_response_region_objects(
        self,
        grid,
        hydro,
        moisture_flow,
        land,
    ) -> dict[str, np.ndarray | list[dict]]:
        """Archive final C4f/C4i wet/dry precipitation-response regions.

        C4d hydroclimate regions describe the drivers.  This C4j object layer
        describes the final active precipitation response after local budget
        conservation, and binds each response patch to source basins, budget
        regions, and moisture-flow networks without changing precipitation.
        """
        land = np.asarray(land, dtype=bool)
        expected = (4, grid.n)
        region_id = np.full(expected, -1.0, dtype=np.float64)
        objects: list[dict] = []
        if not land.any():
            return {"region_id": region_id, "objects": objects}

        seasonal_precip = np.asarray(
            hydro.get("seasonal_precipitation", np.zeros(expected)), dtype=np.float64)
        response = np.asarray(
            hydro.get("moisture_flow_precipitation_response", np.ones(expected)),
            dtype=np.float64,
        )
        budget_id = np.asarray(
            hydro.get("moisture_budget_region_id", np.full(expected, -1.0)),
            dtype=np.float64,
        )
        pathway = np.asarray(
            moisture_flow.get("pathway", np.zeros(expected)), dtype=np.float64)
        source_basin_id = np.asarray(
            moisture_flow.get("source_basin_id", np.full(expected, -1.0)),
            dtype=np.float64,
        )
        network_id = np.asarray(
            moisture_flow.get("network_id", np.full(expected, -1.0)), dtype=np.float64)
        monsoon = np.asarray(
            hydro.get("monsoon_rainfall_corridor", np.zeros(expected)), dtype=np.float64)
        storm = np.asarray(
            hydro.get("storm_track_rainfall_corridor", np.zeros(expected)),
            dtype=np.float64,
        )
        shadow = np.asarray(
            hydro.get("rain_shadow_index", np.zeros(expected)), dtype=np.float64)
        if not (
            seasonal_precip.shape == response.shape == budget_id.shape
            == pathway.shape == source_basin_id.shape == network_id.shape
            == monsoon.shape == storm.shape == shadow.shape == expected
        ):
            return {"region_id": region_id, "objects": objects}

        total_area = max(float(grid.cell_area.sum()), 1.0e-12)
        next_region_id = 1

        def id_summary(values: np.ndarray, cells: np.ndarray,
                       min_valid: float) -> tuple[list[int], int, float, float]:
            vals = np.asarray(values[cells], dtype=np.float64)
            weights = grid.cell_area[cells]
            valid = np.isfinite(vals) & (vals >= min_valid)
            if not valid.any():
                return [], -1, 0.0, 0.0
            ids = sorted({int(x) for x in vals[valid]})
            total = max(float(weights.sum()), 1.0e-12)
            valid_area = float(weights[valid].sum())
            best_id = -1
            best_area = -1.0
            for idx in ids:
                area = float(weights[valid & (vals == float(idx))].sum())
                if area > best_area:
                    best_area = area
                    best_id = int(idx)
            return ids, int(best_id), valid_area / total, max(best_area, 0.0) / total

        def response_thresholds(season: int) -> tuple[float, float]:
            vals = response[season, land]
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                return 1.025, 0.975
            return max(1.025, float(np.percentile(vals, 84))), min(
                0.975, float(np.percentile(vals, 16)))

        for season in range(4):
            wet_threshold, dry_threshold = response_thresholds(season)
            specs = (
                ("wet_precipitation_response_region",
                 land & np.isfinite(response[season]) & (response[season] >= wet_threshold),
                 wet_threshold,
                 np.maximum(response[season] - 1.0, 0.0)),
                ("dry_precipitation_response_region",
                 land & np.isfinite(response[season]) & (response[season] <= dry_threshold),
                 dry_threshold,
                 np.maximum(1.0 - response[season], 0.0)),
            )
            for kind, active, threshold, intensity in specs:
                if not active.any():
                    continue
                component_id, component_area = self._connected_components(grid, active)
                component_rows: list[tuple[float, np.ndarray, dict]] = []
                for cid in [int(x) for x in np.unique(component_id[active]) if int(x) >= 0]:
                    cells = active & (component_id == cid)
                    cell_count = int(cells.sum())
                    area = float(component_area[cid])
                    area_fraction = area / total_area
                    if cell_count < 3 or area_fraction < 7.5e-4:
                        continue
                    weights = grid.cell_area[cells]
                    centroid = np.average(grid.xyz[cells], axis=0, weights=weights)
                    cn = float(np.linalg.norm(centroid))
                    if cn > 1.0e-12:
                        centroid = centroid / cn
                    lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
                    lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
                    values = response[season, cells]
                    anomaly = np.abs(values - 1.0)
                    source_ids, source_dom, source_attr, source_purity = id_summary(
                        source_basin_id[season], cells, 0.0)
                    budget_ids, budget_dom, budget_attr, budget_purity = id_summary(
                        budget_id[season], cells, 1.0)
                    flow_ids, flow_dom, flow_attr, flow_purity = id_summary(
                        network_id[season], cells, 1.0)
                    row = {
                        "id": (
                            "precipitation_response:"
                            f"{self.SEASON_NAMES[season].lower()}:{next_region_id}"
                        ),
                        "type": "precipitation_response_region",
                        "kind": kind,
                        "season": self.SEASON_NAMES[season],
                        "season_index": int(season),
                        "region_index": int(next_region_id),
                        "threshold": float(threshold),
                        "cell_count": cell_count,
                        "area_fraction": area_fraction,
                        "centroid_lat": lat,
                        "centroid_lon": lon,
                        "lat_min": float(np.min(grid.lat[cells])),
                        "lat_max": float(np.max(grid.lat[cells])),
                        "mean_response": float(np.average(values, weights=weights)),
                        "mean_abs_response_anomaly": float(
                            np.average(anomaly, weights=weights)),
                        "max_abs_response_anomaly": float(np.max(anomaly)),
                        "mean_precipitation_mm_yr_equivalent": float(
                            np.average(seasonal_precip[season, cells], weights=weights)),
                        "mean_pathway": float(np.average(
                            pathway[season, cells], weights=weights)),
                        "mean_monsoon_corridor": float(np.average(
                            monsoon[season, cells], weights=weights)),
                        "mean_storm_track_corridor": float(np.average(
                            storm[season, cells], weights=weights)),
                        "mean_rain_shadow": float(np.average(
                            shadow[season, cells], weights=weights)),
                        "source_basin_ids": source_ids,
                        "dominant_source_basin_id": int(source_dom),
                        "source_basin_attributed_fraction": float(source_attr),
                        "source_basin_purity": float(source_purity),
                        "budget_region_ids": budget_ids,
                        "dominant_budget_region_id": int(budget_dom),
                        "budget_region_attributed_fraction": float(budget_attr),
                        "budget_region_purity": float(budget_purity),
                        "flow_network_ids": flow_ids,
                        "dominant_flow_network_id": int(flow_dom),
                        "flow_network_attributed_fraction": float(flow_attr),
                        "flow_network_purity": float(flow_purity),
                    }
                    component_rows.append((
                        area_fraction * max(float(row["mean_abs_response_anomaly"]), 1.0e-9),
                        cells,
                        row,
                    ))
                    next_region_id += 1

                component_rows.sort(key=lambda item: item[0], reverse=True)
                for _, cells, row in component_rows[:12]:
                    objects.append(row)
                    region_id[season, cells] = float(row["region_index"])

        objects.sort(key=lambda item: (
            str(item["season"]),
            str(item["kind"]),
            -float(item["area_fraction"]),
        ))
        return {"region_id": region_id, "objects": objects}

    def _receiver_catchment_objects(
        self,
        grid,
        hydro,
        moisture_flow,
        precipitation_response,
        land,
    ) -> dict[str, np.ndarray | list[dict]]:
        """Merge response patches into stable seasonal receiver catchments.

        C4j regions deliberately describe wet/dry response patches.  C4k
        catchments are coarser receiving-side accounting objects: contiguous
        seasonal land domains inside local moisture-budget regions, split by
        source basin only when that source materially controls a sub-domain.
        This layer is diagnostic and does not change precipitation.
        """
        land = np.asarray(land, dtype=bool)
        expected = (4, grid.n)
        catchment_id = np.full(expected, -1.0, dtype=np.float64)
        objects: list[dict] = []
        if not land.any():
            return {"catchment_id": catchment_id, "objects": objects}

        seasonal_precip = np.asarray(
            hydro.get("seasonal_precipitation", np.zeros(expected)), dtype=np.float64)
        response = np.asarray(
            hydro.get("moisture_flow_precipitation_response", np.ones(expected)),
            dtype=np.float64,
        )
        budget_id = np.asarray(
            hydro.get("moisture_budget_region_id", np.full(expected, -1.0)),
            dtype=np.float64,
        )
        source_basin_id = np.asarray(
            moisture_flow.get("source_basin_id", np.full(expected, -1.0)),
            dtype=np.float64,
        )
        pathway = np.asarray(
            moisture_flow.get("pathway", np.zeros(expected)), dtype=np.float64)
        network_id = np.asarray(
            moisture_flow.get("network_id", np.full(expected, -1.0)), dtype=np.float64)
        precip_region_id = np.asarray(
            precipitation_response.get("region_id", np.full(expected, -1.0)),
            dtype=np.float64,
        )
        if not (
            seasonal_precip.shape == response.shape == budget_id.shape
            == source_basin_id.shape == pathway.shape == network_id.shape
            == precip_region_id.shape == expected
        ):
            return {"catchment_id": catchment_id, "objects": objects}

        total_area = max(float(grid.cell_area.sum()), 1.0e-12)
        next_catchment_id = 1

        def id_summary(values: np.ndarray, cells: np.ndarray,
                       min_valid: float) -> tuple[list[int], int, float, float]:
            vals = np.asarray(values[cells], dtype=np.float64)
            weights = grid.cell_area[cells]
            valid = np.isfinite(vals) & (vals >= min_valid)
            if not valid.any():
                return [], -1, 0.0, 0.0
            ids = sorted({int(x) for x in vals[valid]})
            total = max(float(weights.sum()), 1.0e-12)
            valid_area = float(weights[valid].sum())
            best_id = -1
            best_area = -1.0
            for idx in ids:
                area = float(weights[valid & (vals == float(idx))].sum())
                if area > best_area:
                    best_area = area
                    best_id = int(idx)
            return ids, int(best_id), valid_area / total, max(best_area, 0.0) / total

        def add_component(season: int, cells: np.ndarray,
                          forced_source_id: int | None,
                          forced_budget_id: int | None = None) -> None:
            nonlocal next_catchment_id
            cell_count = int(np.count_nonzero(cells))
            if cell_count < 1:
                return
            area = float(grid.cell_area[cells].sum())
            area_fraction = area / total_area
            weights = grid.cell_area[cells]
            centroid = np.average(grid.xyz[cells], axis=0, weights=weights)
            cn = float(np.linalg.norm(centroid))
            if cn > 1.0e-12:
                centroid = centroid / cn
            lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
            lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
            source_ids, source_dom, source_attr, source_purity = id_summary(
                source_basin_id[season], cells, 0.0)
            budget_ids, budget_dom, budget_attr, budget_purity = id_summary(
                budget_id[season], cells, 1.0)
            residual_budget_region = False
            if forced_budget_id is not None and budget_attr <= 0.0:
                budget_ids = [int(forced_budget_id)]
                budget_dom = int(forced_budget_id)
                budget_attr = 1.0
                budget_purity = 1.0
                residual_budget_region = True
            flow_ids, flow_dom, flow_attr, flow_purity = id_summary(
                network_id[season], cells, 1.0)
            precip_ids, precip_dom, precip_attr, precip_purity = id_summary(
                precip_region_id[season], cells, 1.0)
            values = response[season, cells]
            anomaly = np.abs(values - 1.0)
            kind = (
                "source_receiver_catchment"
                if source_dom >= 0 and source_purity >= 0.55
                else "mixed_receiver_catchment"
            )
            row = {
                "id": (
                    "receiver_catchment:"
                    f"{self.SEASON_NAMES[season].lower()}:{next_catchment_id}"
                ),
                "type": "receiver_catchment",
                "kind": kind,
                "season": self.SEASON_NAMES[season],
                "season_index": int(season),
                "catchment_index": int(next_catchment_id),
                "cell_count": cell_count,
                "area_fraction": area_fraction,
                "centroid_lat": lat,
                "centroid_lon": lon,
                "lat_min": float(np.min(grid.lat[cells])),
                "lat_max": float(np.max(grid.lat[cells])),
                "mean_precipitation_mm_yr_equivalent": float(
                    np.average(seasonal_precip[season, cells], weights=weights)),
                "mean_response": float(np.average(values, weights=weights)),
                "mean_abs_response_anomaly": float(
                    np.average(anomaly, weights=weights)),
                "mean_pathway": float(
                    np.average(pathway[season, cells], weights=weights)),
                "source_basin_ids": source_ids,
                "dominant_source_basin_id": int(
                    forced_source_id if forced_source_id is not None else source_dom),
                "source_basin_attributed_fraction": float(source_attr),
                "source_basin_purity": float(source_purity),
                "budget_region_ids": budget_ids,
                "dominant_budget_region_id": int(budget_dom),
                "budget_region_attributed_fraction": float(budget_attr),
                "budget_region_purity": float(budget_purity),
                "residual_budget_region": bool(residual_budget_region),
                "flow_network_ids": flow_ids,
                "dominant_flow_network_id": int(flow_dom),
                "flow_network_attributed_fraction": float(flow_attr),
                "flow_network_purity": float(flow_purity),
                "precipitation_response_region_ids": precip_ids,
                "dominant_precipitation_response_region_id": int(precip_dom),
                "precipitation_response_attributed_fraction": float(precip_attr),
                "precipitation_response_purity": float(precip_purity),
            }
            catchment_id[season, cells] = float(row["catchment_index"])
            objects.append(row)
            next_catchment_id += 1

        for season in range(4):
            season_budget = budget_id[season]
            budget_values = [
                int(x) for x in np.unique(season_budget[land])
                if np.isfinite(x) and int(x) > 0
            ]
            covered_by_budget = np.zeros(grid.n, dtype=bool)
            if not budget_values:
                budget_values = [1]
            for bid in budget_values:
                if bid == 1 and not np.any(season_budget[land] > 0):
                    base = land.copy()
                    forced_budget_id = 0
                else:
                    base = land & (season_budget == float(bid))
                    forced_budget_id = None
                if int(np.count_nonzero(base)) < 1:
                    continue
                covered_by_budget |= base
                base_area = max(float(grid.cell_area[base].sum()), 1.0e-12)
                source_candidates: list[tuple[float, int, np.ndarray]] = []
                for source_value in [
                    int(x) for x in np.unique(source_basin_id[season, base])
                    if np.isfinite(x) and int(x) >= 0
                ]:
                    mask = base & (source_basin_id[season] == float(source_value))
                    share = float(grid.cell_area[mask].sum()) / base_area
                    if share >= 0.10 and int(np.count_nonzero(mask)) >= 3:
                        source_candidates.append((share, int(source_value), mask))
                source_candidates.sort(reverse=True, key=lambda item: item[0])

                assigned = np.zeros(grid.n, dtype=bool)
                if len(source_candidates) >= 2:
                    for _, source_value, mask in source_candidates:
                        comp_id, _ = self._connected_components(grid, mask)
                        for cid in [
                            int(x) for x in np.unique(comp_id[mask])
                            if int(x) >= 0
                        ]:
                            cells = mask & (comp_id == cid) & ~assigned
                            add_component(season, cells, source_value, forced_budget_id)
                            assigned |= cells

                residual = base & ~assigned
                if residual.any():
                    comp_id, _ = self._connected_components(grid, residual)
                    for cid in [
                        int(x) for x in np.unique(comp_id[residual])
                        if int(x) >= 0
                    ]:
                        cells = residual & (comp_id == cid)
                        add_component(season, cells, None, forced_budget_id)

            residual_budget = land & ~covered_by_budget
            if residual_budget.any():
                comp_id, _ = self._connected_components(grid, residual_budget)
                for cid in [
                    int(x) for x in np.unique(comp_id[residual_budget])
                    if int(x) >= 0
                ]:
                    cells = residual_budget & (comp_id == cid)
                    add_component(season, cells, None, 0)

        objects.sort(key=lambda item: (
            str(item["season"]),
            -float(item["area_fraction"]),
            int(item["catchment_index"]),
        ))
        return {"catchment_id": catchment_id, "objects": objects}

    def _source_basin_receiver_accounting(
        self,
        grid,
        hydro,
        moisture_flow,
        receiver_catchments,
        land,
        ocean,
    ) -> dict[str, np.ndarray | list[dict] | dict[str, float]]:
        """Diagnose seasonal source-basin supply versus receiver catchments.

        This is an accounting layer, not a new precipitation solve.  It turns
        the C4i source-basin labels and C4k receiver catchments into bounded
        support and consistency fields that later gates can audit.
        """
        land = np.asarray(land, dtype=bool)
        ocean = np.asarray(ocean, dtype=bool)
        expected = (4, grid.n)
        source_basin_supply_index = np.zeros(expected, dtype=np.float64)
        receiver_supply_balance = np.zeros(expected, dtype=np.float64)
        objects = list(receiver_catchments.get("objects", []))
        diag = {
            "enabled": 0.0,
            "source_basin_supply_index_land_p50": 0.0,
            "source_basin_supply_attributed_land_fraction": 0.0,
            "receiver_catchment_supply_balance_land_p50": 0.0,
        }
        if not land.any() or not ocean.any():
            return {
                "source_basin_supply_index": source_basin_supply_index,
                "receiver_catchment_supply_balance": receiver_supply_balance,
                "objects": objects,
                "diagnostics": diag,
            }

        source_strength = np.asarray(
            moisture_flow.get("source", np.zeros(expected)), dtype=np.float64)
        pathway = np.asarray(
            moisture_flow.get("pathway", np.zeros(expected)), dtype=np.float64)
        source_basin_id = np.asarray(
            moisture_flow.get("source_basin_id", np.full(expected, -1.0)),
            dtype=np.float64,
        )
        catchment_id = np.asarray(
            receiver_catchments.get("catchment_id", np.full(expected, -1.0)),
            dtype=np.float64,
        )
        seasonal_precip = np.asarray(
            hydro.get("seasonal_precipitation", np.zeros(expected)), dtype=np.float64)
        response = np.asarray(
            hydro.get("moisture_flow_precipitation_response", np.ones(expected)),
            dtype=np.float64,
        )
        if not (
            source_strength.shape == pathway.shape == source_basin_id.shape
            == catchment_id.shape == seasonal_precip.shape == response.shape
            == expected
        ):
            return {
                "source_basin_supply_index": source_basin_supply_index,
                "receiver_catchment_supply_balance": receiver_supply_balance,
                "objects": objects,
                "diagnostics": diag,
            }

        basin_supply_share: dict[tuple[int, int], float] = {}
        for season in range(4):
            labels = source_basin_id[season]
            source = np.where(ocean, np.clip(source_strength[season], 0.0, 2.0), 0.0)
            labelled_ocean = (
                ocean
                & np.isfinite(labels)
                & (labels >= 0.0)
                & np.isfinite(source)
                & (source > 0.0)
            )
            if not labelled_ocean.any():
                continue

            source_vals = source[labelled_ocean]
            source_scale = max(float(np.percentile(source_vals, 90)), 1.0e-9)
            source_basin_supply_index[season, ocean] = np.clip(
                source[ocean] / source_scale, 0.0, 1.5)

            basin_mean: dict[int, float] = {}
            basin_mass: dict[int, float] = {}
            total_mass = 0.0
            for bid in [
                int(x) for x in np.unique(labels[labelled_ocean])
                if np.isfinite(x) and int(x) >= 0
            ]:
                basin_ocean = labelled_ocean & (labels == float(bid))
                weights = grid.cell_area[basin_ocean]
                mass = float(np.sum(weights * source[basin_ocean]))
                total_mass += mass
                basin_mass[int(bid)] = mass
                basin_mean[int(bid)] = (
                    float(np.average(source[basin_ocean], weights=weights))
                    if weights.size else 0.0
                )
            if total_mass > 1.0e-12:
                for bid, mass in basin_mass.items():
                    basin_supply_share[(season, bid)] = mass / total_mass

            raw_land = np.zeros(grid.n, dtype=np.float64)
            land_basin_ids = [
                int(x) for x in np.unique(labels[land])
                if np.isfinite(x) and int(x) >= 0
            ]
            for bid in land_basin_ids:
                cells = land & (labels == float(bid))
                if not cells.any():
                    continue
                mean_source = basin_mean.get(int(bid), source_scale)
                source_factor = np.clip(mean_source / source_scale, 0.25, 1.6)
                raw_land[cells] = (
                    np.clip(pathway[season, cells], 0.0, 1.5) * source_factor
                )
            positive = raw_land[land & np.isfinite(raw_land) & (raw_land > 0.0)]
            land_scale = max(
                float(np.percentile(positive, 90)) if positive.size else 0.0,
                1.0e-9,
            )
            source_basin_supply_index[season, land] = np.clip(
                raw_land[land] / land_scale, 0.0, 1.5)

            precip_land = seasonal_precip[season, land]
            precip_scale = max(
                float(np.percentile(precip_land[np.isfinite(precip_land)], 90))
                if np.isfinite(precip_land).any() else 0.0,
                1.0e-9,
            )
            supply_land = source_basin_supply_index[season, land]
            supply_scale = max(
                float(np.percentile(supply_land[supply_land > 0.0], 90))
                if np.any(supply_land > 0.0) else 0.0,
                1.0e-9,
            )
            precip_norm = np.clip(seasonal_precip[season] / precip_scale, 0.0, 2.0)
            supply_norm = np.clip(
                source_basin_supply_index[season] / supply_scale, 0.0, 2.0)

            for rid in [
                int(x) for x in np.unique(catchment_id[season, land])
                if np.isfinite(x) and int(x) > 0
            ]:
                cells = land & (catchment_id[season] == float(rid))
                if not cells.any():
                    continue
                weights = grid.cell_area[cells]
                mean_precip_norm = float(np.average(
                    precip_norm[cells], weights=weights))
                mean_supply_norm = float(np.average(
                    supply_norm[cells], weights=weights))
                balance = (
                    min(mean_precip_norm + 0.08, mean_supply_norm + 0.08)
                    / max(mean_precip_norm + 0.08, mean_supply_norm + 0.08)
                )
                response_anomaly = float(np.average(
                    np.abs(response[season, cells] - 1.0), weights=weights))
                if response_anomaly < 0.015:
                    support_area = float(np.sum(
                        weights[source_basin_supply_index[season, cells] > 0.05])
                    ) / max(float(np.sum(weights)), 1.0e-12)
                    balance = max(balance, 0.45 + 0.35 * support_area)
                receiver_supply_balance[season, cells] = np.clip(balance, 0.0, 1.0)

        updated_objects: list[dict] = []
        for obj in objects:
            if not isinstance(obj, dict):
                updated_objects.append(obj)
                continue
            row = dict(obj)
            season = int(row.get("season_index", -1))
            catchment_index = int(row.get("catchment_index", -1))
            if season < 0 or season >= 4 or catchment_index <= 0:
                updated_objects.append(row)
                continue
            cells = land & (catchment_id[season] == float(catchment_index))
            if not cells.any():
                updated_objects.append(row)
                continue
            weights = grid.cell_area[cells]
            supply = source_basin_supply_index[season, cells]
            precip_vals = seasonal_precip[season, cells]
            supported = np.clip(supply / (supply + 0.35), 0.0, 1.0)
            precip_mass = float(np.sum(weights * np.maximum(precip_vals, 0.0)))
            supported_precip_mass = float(np.sum(
                weights * np.maximum(precip_vals, 0.0) * supported))
            dominant_source = int(row.get("dominant_source_basin_id", -1))
            row.update({
                "mean_source_basin_supply_index": float(np.average(
                    supply, weights=weights)),
                "source_basin_supply_attributed_fraction": float(
                    np.sum(weights[supply > 0.05])
                    / max(float(np.sum(weights)), 1.0e-12)
                ),
                "source_basin_supply_mass_fraction": float(
                    basin_supply_share.get((season, dominant_source), 0.0)
                    if dominant_source >= 0 else 0.0
                ),
                "supply_supported_precipitation_fraction": float(
                    supported_precip_mass / max(precip_mass, 1.0e-12)
                ),
                "precipitation_supply_balance": float(np.average(
                    receiver_supply_balance[season, cells], weights=weights)),
            })
            updated_objects.append(row)

        land_supply = source_basin_supply_index[:, land]
        land_balance = receiver_supply_balance[:, land]
        attributed = (
            land[None, :]
            & np.isfinite(source_basin_id)
            & (source_basin_id >= 0.0)
            & (source_basin_supply_index > 0.05)
        )
        diag = {
            "enabled": 1.0,
            "source_basin_supply_index_land_p50": (
                float(np.percentile(land_supply[np.isfinite(land_supply)], 50))
                if np.isfinite(land_supply).any() else 0.0
            ),
            "source_basin_supply_attributed_land_fraction": float(
                np.count_nonzero(attributed)
                / max(4 * int(np.count_nonzero(land)), 1)
            ),
            "receiver_catchment_supply_balance_land_p50": (
                float(np.percentile(land_balance[np.isfinite(land_balance)], 50))
                if np.isfinite(land_balance).any() else 0.0
            ),
        }
        return {
            "source_basin_supply_index": source_basin_supply_index,
            "receiver_catchment_supply_balance": receiver_supply_balance,
            "objects": updated_objects,
            "diagnostics": diag,
        }

    def _apply_receiver_supply_precipitation_feedback(
        self,
        grid,
        hydro,
        source_receiver_accounting,
        land,
    ):
        """Feed C5e7 source/receiver accounting back into precipitation gently.

        This is a second, bounded land-only redistribution pass after C4f.  It
        preserves every seasonal moisture-budget-region mean, leaves ocean
        precipitation untouched, and only nudges precipitation toward cells with
        diagnosed source-basin supply and receiver-catchment support.
        """
        expected = (4, grid.n)
        seasonal = np.asarray(
            hydro.get("seasonal_precipitation", np.zeros(expected)), dtype=np.float64)
        budget_id = np.asarray(
            hydro.get("moisture_budget_region_id", np.full(expected, -1.0)),
            dtype=np.float64,
        )
        supply = np.asarray(
            source_receiver_accounting.get(
                "source_basin_supply_index", np.zeros(expected)),
            dtype=np.float64,
        )
        balance = np.asarray(
            source_receiver_accounting.get(
                "receiver_catchment_supply_balance", np.zeros(expected)),
            dtype=np.float64,
        )
        c4f_response = np.asarray(
            hydro.get("moisture_flow_precipitation_response", np.ones(expected)),
            dtype=np.float64,
        )
        land = np.asarray(land, dtype=bool)
        feedback = np.ones(expected, dtype=np.float64)
        diagnostics = {
            "enabled": False,
            "response_land_p05": 1.0,
            "response_land_p95": 1.0,
            "max_land_mean_delta_mm_yr": 0.0,
            "max_budget_region_mean_delta_mm_yr": 0.0,
            "budget_region_count_p50": 0.0,
            "amplitude": 0.0,
        }
        if (
            seasonal.shape != expected
            or budget_id.shape != expected
            or supply.shape != expected
            or balance.shape != expected
            or c4f_response.shape != expected
            or not land.any()
        ):
            out = dict(hydro)
            out["receiver_supply_precipitation_feedback"] = feedback
            return self._recompute_precipitation_diagnostics(grid, out, land), diagnostics

        land_area = max(float(grid.cell_area[land].sum()), 1.0e-12)
        total_area = max(float(grid.cell_area.sum()), 1.0e-12)
        land_fraction = land_area / total_area
        land_extent = float(np.clip((land_fraction - 0.035) / 0.22, 0.0, 1.0))
        amplitude = 0.018 + 0.027 * land_extent

        shaped = seasonal.copy()
        land_mean_deltas: list[float] = []
        budget_mean_deltas: list[float] = []
        budget_region_counts: list[int] = []
        for season in range(4):
            before = seasonal[season].copy()
            target_mean = float(np.average(before[land], weights=grid.cell_area[land]))
            supply_s = np.where(land, np.clip(supply[season], 0.0, 1.5), 0.0)
            supply_vals = supply_s[land & np.isfinite(supply_s)]
            positive_supply = supply_vals[supply_vals > 0.0]
            supply_scale = max(
                float(np.percentile(positive_supply, 90)) if positive_supply.size else 0.0,
                1.0e-9,
            )
            supply_norm = np.clip(supply_s / supply_scale, 0.0, 1.5)
            balance_norm = np.where(
                land,
                np.clip(balance[season], 0.0, 1.0),
                0.0,
            )
            wet_support = np.where(
                land,
                np.clip(c4f_response[season] - 1.0, 0.0, 0.35) / 0.35,
                0.0,
            )
            signal = np.where(
                land,
                0.58 * supply_norm + 0.30 * balance_norm + 0.12 * wet_support,
                0.0,
            )
            signal = self._smooth_field_masked(grid, signal, land, passes=1, alpha=0.07)
            signal_mean = float(np.average(signal[land], weights=grid.cell_area[land]))
            centered = signal - signal_mean
            spread = max(float(np.percentile(np.abs(centered[land]), 86)), 0.08)
            raw_feedback = np.where(
                land,
                np.clip(1.0 + amplitude * centered / spread, 0.92, 1.10),
                1.0,
            )
            candidate = before * raw_feedback
            candidate, budget_deltas, budget_count = self._preserve_moisture_budget_means(
                grid,
                candidate,
                before,
                land,
                budget_id[season],
                cap=4500.0,
            )
            budget_mean_deltas.extend(budget_deltas)
            budget_region_counts.append(int(budget_count))
            after_mean = float(np.average(candidate[land], weights=grid.cell_area[land]))
            land_mean_deltas.append(abs(after_mean - target_mean))
            shaped[season] = np.where(land, candidate, before)
            feedback[season] = np.where(
                before > 1.0e-9,
                np.clip(shaped[season] / np.maximum(before, 1.0e-9), 0.0, 4.0),
                1.0,
            )

        out = dict(hydro)
        out["seasonal_precipitation"] = shaped
        out["receiver_supply_precipitation_feedback"] = feedback
        out = self._recompute_precipitation_diagnostics(grid, out, land)
        vals = feedback[:, land].ravel()
        diagnostics = {
            "enabled": True,
            "response_land_p05": float(np.percentile(vals, 5)) if vals.size else 1.0,
            "response_land_p95": float(np.percentile(vals, 95)) if vals.size else 1.0,
            "max_land_mean_delta_mm_yr": float(max(land_mean_deltas, default=0.0)),
            "max_budget_region_mean_delta_mm_yr": float(
                max(budget_mean_deltas, default=0.0)),
            "budget_region_count_p50": (
                float(np.percentile(budget_region_counts, 50))
                if budget_region_counts else 0.0
            ),
            "amplitude": float(amplitude),
            "land_fraction": float(land_fraction),
        }
        return out, diagnostics

    def _geography_primitives(self, world, ocean, rel_elev, climate_elev):
        grid = world.grid
        land = ~ocean
        east, _ = self._tangent_basis(grid)
        continent_id, continent_area = self._connected_components(grid, land)
        basin_id, basin_area = self._connected_components(grid, ocean)
        basin_name_by_id: dict[int, str] = {}
        basin_source = "connected_components"
        total_area = max(float(grid.cell_area.sum()), 1.0e-12)
        terrain_basin = world.fields.get("ocean.basin_id")
        if terrain_basin is not None:
            candidate = np.asarray(terrain_basin)
            if candidate.shape == (grid.n,):
                candidate = candidate.astype(int)
                valid_ocean = bool((candidate[ocean] >= 0).all()) if ocean.any() else True
                valid_land = bool((candidate[land] < 0).all()) if land.any() else True
                if valid_ocean and valid_land:
                    basin_id = candidate
                    ids = [int(x) for x in np.unique(basin_id[ocean]) if int(x) >= 0]
                    basin_area = np.zeros(max(ids) + 1 if ids else 0, dtype=np.float64)
                    for bid in ids:
                        basin_area[bid] = float(grid.cell_area[ocean & (
                            basin_id == bid)].sum())
                    basin_source = "input_field"

        if world.g("diagnostics.real_earth_replay", 0.0) > 0.5:
            basin_id, basin_area, basin_name_by_id = (
                self._real_earth_major_ocean_basins(grid, ocean)
            )
            basin_source = "real_earth_major_oceans"

        i, j = self._edges[:, 0], self._edges[:, 1]
        land_coast = np.zeros(grid.n, dtype=bool)
        ocean_coast = np.zeros(grid.n, dtype=bool)
        land_to_ocean = land[i] & ocean[j]
        ocean_to_land = ocean[i] & land[j]
        land_coast[i[land_to_ocean]] = True
        ocean_coast[j[land_to_ocean]] = True
        ocean_coast[i[ocean_to_land]] = True
        land_coast[j[ocean_to_land]] = True
        coast_all = land_coast | ocean_coast

        coast_orientation, coast_strength = self._coastal_direction(grid, land, ocean)
        coast_orientation[~land_coast] = 0.0
        coast_strength = np.where(land_coast, coast_strength, 0.0)
        coast_facing_east = np.where(
            land_coast, np.sum(coast_orientation * east, axis=1), 0.0)

        land_dist = self._graph_hop_distance(grid, land_coast, land)
        ocean_dist = self._graph_hop_distance(grid, ocean_coast, ocean)
        all_dist = self._graph_hop_distance(grid, coast_all, np.ones(grid.n, dtype=bool))
        coast_distance = np.where(all_dist >= 0, all_dist.astype(float), 0.0)
        cd_scale = max(float(np.percentile(coast_distance[coast_distance > 0], 95))
                       if (coast_distance > 0).any() else 1.0, 1.0)
        coast_distance = np.clip(coast_distance / cd_scale, 0.0, 1.0)

        land_dist_f = np.where(land_dist >= 0, land_dist.astype(float), 0.0)
        ld_scale = max(float(np.percentile(land_dist_f[land], 95))
                       if land.any() else 1.0, 1.0)
        interior_base = np.clip(land_dist_f / ld_scale, 0.0, 1.0)
        continent_scale = np.zeros(grid.n, dtype=np.float64)
        if continent_area.size:
            max_continent = max(float(continent_area.max()), 1e-12)
            area_fraction = continent_area[continent_id[land]] / total_area
            absolute_scale = np.clip(np.sqrt(area_fraction / 0.10), 0.04, 1.0)
            relative_scale = np.sqrt(continent_area[continent_id[land]] / max_continent)
            continent_scale[land] = np.clip(
                relative_scale * (0.12 + 0.88 * absolute_scale), 0.03, 1.0)
        continent_interiority = np.where(
            land, np.clip(interior_base * (0.35 + 0.75 * continent_scale), 0.0, 1.0), 0.0)

        ocean_dist_f = np.where(ocean_dist >= 0, ocean_dist.astype(float), 99.0)
        depth = np.maximum(-rel_elev, 0.0)
        near_coast = np.exp(-ocean_dist_f / 4.0)
        shallow = np.clip(1.0 - depth / 3500.0, 0.0, 1.0)
        shelf_index = np.where(ocean, np.clip(0.68 * near_coast + 0.32 * shallow,
                                              0.0, 1.0), 0.0)
        terrain_shelf_width = world.fields.get("ocean.shelf_width")
        if terrain_shelf_width is not None:
            width = np.asarray(terrain_shelf_width, dtype=np.float64)
            if width.shape == (grid.n,):
                shelf_from_width = np.where(
                    ocean & (width > 0.0), np.clip(1.15 - 0.24 * (width - 1.0), 0.0, 1.0), 0.0)
                shelf_index = np.where(
                    ocean, np.maximum(shelf_index, 0.70 * shelf_from_width + 0.30 * shallow),
                    0.0)
        depth_province = world.fields.get("ocean.depth_province")
        if depth_province is not None:
            prov = np.asarray(depth_province, dtype=int)
            if prov.shape == (grid.n,):
                province_shelf = np.zeros(grid.n, dtype=np.float64)
                province_shelf[prov == 1] = 1.0
                province_shelf[prov == 2] = 0.62
                province_shelf[prov == 3] = 0.34
                province_shelf[prov == 7] = 0.72
                province_shelf[prov == 5] = 0.22
                province_shelf[prov == 6] = 0.0
                shelf_index = np.where(ocean, np.maximum(shelf_index, province_shelf), 0.0)
        shelf_index = self._smooth_field_masked(grid, shelf_index, ocean,
                                                passes=1, alpha=0.18)

        land_neighbor_frac = np.zeros(grid.n, dtype=np.float64)
        opposing_land = np.zeros(grid.n, dtype=np.float64)
        for c in np.where(ocean)[0]:
            nbs = grid.neighbors[int(c)]
            if nbs.size == 0:
                continue
            land_nbs = nbs[land[nbs]]
            land_neighbor_frac[c] = land_nbs.size / nbs.size
            if land_nbs.size >= 2:
                dirs = []
                for nb in land_nbs:
                    vec = grid.xyz[int(nb)] - float(grid.xyz[int(nb)] @ grid.xyz[c]) * grid.xyz[c]
                    norm = float(np.linalg.norm(vec))
                    if norm > 1e-12:
                        dirs.append(vec / norm)
                if dirs:
                    mean_dir = np.mean(np.asarray(dirs), axis=0)
                    opposing_land[c] = np.clip(1.0 - float(np.linalg.norm(mean_dir)),
                                               0.0, 1.0)
        narrow_sides = np.clip((land_neighbor_frac - 0.38) / 0.34, 0.0, 1.0)
        opposing_sides = np.clip((opposing_land - 0.34) / 0.42, 0.0, 1.0)
        candidate = (land_neighbor_frac > 0.30) | (opposing_land > 0.36)
        strait_index = np.where(
            ocean & candidate,
            np.clip(0.72 * narrow_sides + 0.55 * opposing_sides, 0.0, 1.0)
            * np.maximum(shelf_index, 0.0) ** 0.75,
            0.0,
        )
        strait_index = self._smooth_field_masked(grid, strait_index, ocean,
                                                 passes=1, alpha=0.10)
        strait_components, _ = self._connected_components(grid, strait_index > 0.22)
        for sid in [int(x) for x in np.unique(strait_components) if int(x) >= 0]:
            cells = strait_components == sid
            if int(cells.sum()) < 3:
                strait_index[cells] = 0.0

        topo = self._smooth_field(grid, np.maximum(climate_elev, 0.0),
                                  passes=6, alpha=0.24)
        topo_land = topo[land]
        topo_scale = max(float(np.percentile(topo_land, 95)) if topo_land.size else 1.0,
                         1.0)
        topo_norm = np.clip(topo / topo_scale, 0.0, 2.0)
        topo_grad = self._graph_gradient_vectors(grid, topo_norm)
        topo_mag = np.linalg.norm(topo_grad, axis=1)
        if land.any():
            g65 = float(np.percentile(topo_mag[land], 65))
            g95 = max(float(np.percentile(topo_mag[land], 95)), g65 + 1e-6)
            h60 = float(np.percentile(topo_norm[land], 60))
            h92 = max(float(np.percentile(topo_norm[land], 92)), h60 + 1e-6)
        else:
            g65, g95, h60, h92 = 0.0, 1.0, 0.0, 1.0
        grad_barrier = np.clip((topo_mag - g65) / (g95 - g65), 0.0, 1.0)
        high_barrier = np.clip((topo_norm - h60) / (h92 - h60), 0.0, 1.0)
        barrier_index = np.where(land, np.clip(0.62 * grad_barrier + 0.48 * high_barrier,
                                               0.0, 1.0), 0.0)
        barrier_index = self._smooth_field_masked(grid, barrier_index, land,
                                                  passes=2, alpha=0.22)
        neighbor_barrier = self._neighbor_mean(grid, barrier_index)
        low_pass = np.clip((1.35 * neighbor_barrier - barrier_index) / 0.34, 0.0, 1.0)
        low_topo = np.clip((0.90 - topo_norm) / 0.90, 0.0, 1.0)
        wind_gap_index = np.where(land, np.clip(1.55 * low_pass * low_topo, 0.0, 1.0),
                                  0.0)
        wind_gap_index = self._smooth_field_masked(grid, wind_gap_index, land,
                                                   passes=1, alpha=0.20)

        fields = {
            "climate.continent_id": continent_id.astype(float),
            "climate.continent_interiority": continent_interiority,
            "climate.coast_orientation": coast_orientation,
            "climate.coast_distance": coast_distance,
            "climate.coast_strength": coast_strength,
            "climate.coast_facing_east": coast_facing_east,
            "ocean.basin_id": basin_id.astype(float),
            "ocean.shelf_index": shelf_index,
            "ocean.strait_index": strait_index,
            "terrain.barrier_index": barrier_index,
            "terrain.wind_gap_index": wind_gap_index,
        }

        continents = self._component_objects(
            grid, continent_id, continent_area, land, adjacent_id=basin_id,
            adjacent_mask=ocean, boundary_mask=land_coast, prefix="continent")
        ocean_basins = self._component_objects(
            grid, basin_id, basin_area, ocean, adjacent_id=continent_id,
            adjacent_mask=land, boundary_mask=ocean_coast, prefix="ocean_basin")
        if basin_name_by_id:
            for obj in ocean_basins:
                bid = int(obj.get("id", -1))
                if bid in basin_name_by_id:
                    obj["name"] = basin_name_by_id[bid]
                    obj["semantic_basin"] = basin_name_by_id[bid]
                    obj["basin_source"] = basin_source
        terrain_ocean_basins = world.objects.get("ocean.basins", [])
        if (
            not basin_name_by_id
            and terrain_ocean_basins
            and any("mean_depth_m" in obj for obj in terrain_ocean_basins)
        ):
            ocean_basins = terrain_ocean_basins

        coastline_pairs: dict[tuple[int, int], int] = {}
        for c in np.where(land_coast)[0]:
            cid = int(continent_id[c])
            for nb in grid.neighbors[int(c)]:
                nb = int(nb)
                bid = int(basin_id[nb]) if ocean[nb] else -1
                if cid >= 0 and bid >= 0:
                    coastline_pairs[(cid, bid)] = coastline_pairs.get((cid, bid), 0) + 1
        coastline_segments = [
            {
                "continent_id": int(cid),
                "basin_id": int(bid),
                "boundary_edge_count": int(count),
            }
            for (cid, bid), count in sorted(coastline_pairs.items(),
                                            key=lambda item: item[1], reverse=True)
        ]

        strait_mask = ocean & (strait_index > 0.45)
        strait_id, strait_area = self._connected_components(grid, strait_mask)
        straits = self._component_objects(
            grid, strait_id, strait_area, strait_mask, adjacent_id=basin_id,
            adjacent_mask=ocean, prefix="strait")
        barrier_mask = land & (barrier_index > 0.55)
        barrier_id, barrier_area = self._connected_components(grid, barrier_mask)
        barrier_belts = self._component_objects(
            grid, barrier_id, barrier_area, barrier_mask, prefix="barrier_belt")

        objects = {
            "climate.continents": continents,
            "ocean.basins": ocean_basins,
            "climate.coastline_segments": coastline_segments,
            "ocean.straits": straits,
            "terrain.barrier_belts": barrier_belts,
        }
        diagnostics = {
            "continent_count": len(continents),
            "ocean_basin_count": len(ocean_basins),
            "coastal_land_cells": int(land_coast.sum()),
            "shelf_p75": float(np.percentile(shelf_index[ocean], 75)) if ocean.any() else 0.0,
            "strait_cells": int(strait_mask.sum()),
            "barrier_cells": int(barrier_mask.sum()),
            "wind_gap_cells": int((wind_gap_index > 0.25).sum()),
            "ocean_basin_source": basin_source,
        }
        return fields, objects, diagnostics

    def _spread_ocean_influence_to_coasts(self, grid, ocean_values, ocean,
                                          passes: int = 4, alpha: float = 0.28,
                                          land_damping: float = 0.62):
        ocean_values = np.asarray(ocean_values, dtype=np.float64)
        ocean = np.asarray(ocean, dtype=bool)
        influence = np.where(ocean, ocean_values, 0.0)
        for _ in range(passes):
            mean = self._neighbor_mean(grid, influence)
            influence = (1.0 - alpha) * influence + alpha * mean
            influence[ocean] = 0.82 * influence[ocean] + 0.18 * ocean_values[ocean]
            influence[~ocean] *= land_damping
        return influence

    def _basin_streamfunction(
        self,
        grid,
        ocean,
        basin_id,
        basin_area,
        annual_wind,
        east,
        north,
        shelf_index,
        basin_scale,
    ):
        wind_east = np.sum(annual_wind * east, axis=1)
        wind_north = np.sum(annual_wind * north, axis=1)
        grad_u = self._graph_gradient_vectors(grid, wind_east)
        grad_v = self._graph_gradient_vectors(grid, wind_north)
        curl = np.sum(grad_v * east, axis=1) - np.sum(grad_u * north, axis=1)
        if ocean.any():
            curl_scale = max(float(np.percentile(np.abs(curl[ocean]), 95)), 1e-9)
        else:
            curl_scale = 1.0
        curl = np.clip(curl / curl_scale, -1.8, 1.8)

        lat_abs = np.abs(grid.lat)
        subtropical = np.exp(-((lat_abs - 27.0) / 18.0) ** 2)
        subpolar = np.exp(-((lat_abs - 52.0) / 16.0) ** 2)
        hemisphere = np.where(grid.lat >= 0.0, 1.0, -1.0)
        planetary_curl = hemisphere * (0.46 * subtropical - 0.26 * subpolar)
        source = np.where(
            ocean,
            (0.70 * curl + planetary_curl)
            * basin_scale
            * np.clip(1.05 - 0.25 * shelf_index, 0.65, 1.10),
            0.0,
        )

        psi = np.zeros(grid.n, dtype=np.float64)
        total_area = max(float(grid.cell_area.sum()), 1.0e-12)
        ids = [int(x) for x in np.unique(basin_id[ocean].astype(int)) if int(x) >= 0]
        for bid in ids:
            cells = ocean & (basin_id.astype(int) == bid)
            if not cells.any():
                continue
            weights = grid.cell_area[cells]
            local_source = np.zeros(grid.n, dtype=np.float64)
            local_source[cells] = source[cells]
            src_mean = float(np.average(local_source[cells], weights=weights))
            local_source[cells] -= src_mean
            local = local_source.copy()
            area_fraction = (
                float(basin_area[bid] / total_area)
                if 0 <= bid < basin_area.size else 0.0
            )
            relax_weight = 0.070 + 0.030 * np.clip(area_fraction / 0.35, 0.0, 1.0)
            for _ in range(92):
                mean = self._masked_neighbor_mean(grid, local, cells)
                local = np.where(
                    cells,
                    (1.0 - relax_weight) * mean + relax_weight * local_source,
                    0.0,
                )
                local[cells] -= float(np.average(local[cells], weights=weights))
            p95 = max(float(np.percentile(np.abs(local[cells]), 95)), 1e-9)
            psi[cells] = np.clip(local[cells] / p95, -1.8, 1.8)

        psi = self._smooth_field_masked(grid, psi, ocean, passes=2, alpha=0.12)
        psi[~ocean] = 0.0

        gyre_id = np.full(grid.n, -1.0, dtype=np.float64)
        next_id = 1
        for bid in ids:
            cells = ocean & (basin_id.astype(int) == bid)
            if not cells.any():
                continue
            local_abs = np.abs(psi[cells])
            threshold = max(float(np.percentile(local_abs, 58)), 0.08)
            for sign in (-1.0, 1.0):
                mask = cells & (sign * psi >= threshold)
                if mask.any():
                    gyre_id[mask] = float(next_id)
                    next_id += 1
            weak = cells & (gyre_id < 0.0)
            gyre_id[weak] = 0.0
        return psi, gyre_id

    def _ocean_currents(self, world, ocean, seasonal_wind, geography_fields=None):
        grid = world.grid
        ocean = np.asarray(ocean, dtype=bool)
        east, north = self._tangent_basis(grid)
        basin_id, basin_area = self._ocean_components(grid, ocean)
        geography_fields = geography_fields or {}
        zero = np.zeros(grid.n, dtype=np.float64)
        shelf_index = np.asarray(
            geography_fields.get("ocean.shelf_index", zero), dtype=np.float64)
        strait_index = np.asarray(
            geography_fields.get("ocean.strait_index", zero), dtype=np.float64)
        basin_scale = np.zeros(grid.n, dtype=np.float64)
        if basin_area.size:
            max_area = max(float(basin_area.max()), 1e-12)
            basin_scale[ocean] = np.clip(
                np.sqrt(basin_area[basin_id[ocean]] / max_area), 0.18, 1.0)

        annual_wind = np.asarray(seasonal_wind, dtype=np.float64).mean(axis=0)
        wind_east = np.sum(annual_wind * east, axis=1)
        wind_north = np.sum(annual_wind * north, axis=1)
        lat_abs = np.abs(grid.lat)
        subtropics = np.exp(-((lat_abs - 27.0) / 18.0) ** 2)

        streamfunction, gyre_id = self._basin_streamfunction(
            grid, ocean, basin_id, basin_area, annual_wind, east, north,
            shelf_index, basin_scale)
        psi_grad = self._graph_gradient_vectors(grid, streamfunction)
        stream_flow = self._project_tangent(grid, np.cross(grid.xyz, psi_grad))
        if ocean.any():
            stream_scale = max(
                float(np.percentile(np.linalg.norm(stream_flow[ocean], axis=1), 95)),
                1e-9,
            )
        else:
            stream_scale = 1.0
        stream_flow = stream_flow / stream_scale

        wind_vector = self._project_tangent(
            grid, wind_east[:, None] * east + wind_north[:, None] * north)
        wind_speed = np.linalg.norm(wind_vector, axis=1)
        wind_p95 = max(
            float(np.percentile(wind_speed[ocean], 95)) if ocean.any() else 0.0,
            1.0e-9,
        )
        wind_unit = np.where(
            wind_speed[:, None] > 1.0e-9,
            wind_vector / np.maximum(wind_speed[:, None], 1.0e-9),
            0.0,
        )
        local_wind_drag = 0.018 * wind_vector
        direct_wind_response = (
            local_wind_drag * np.clip(1.0 - 0.45 * shelf_index, 0.45, 1.0)[:, None]
        )
        open_ocean_stress = (
            ocean.astype(float)
            * np.clip(1.0 - shelf_index, 0.0, 1.0)
            * np.clip(0.35 + 0.65 * basin_scale, 0.0, 1.0)
            * np.clip(wind_speed / wind_p95, 0.0, 1.4)
        )
        direct_wind_response += 0.006 * open_ocean_stress[:, None] * wind_unit
        direct_wind_response = self._project_tangent(grid, direct_wind_response)
        direct_wind_response[~ocean] = 0.0
        response_speed = np.linalg.norm(direct_wind_response, axis=1)
        direct_wind_response *= np.minimum(
            1.0, 0.34 / np.maximum(response_speed, 1.0e-9))[:, None]
        currents = (
            0.44 * basin_scale[:, None] * stream_flow
            + direct_wind_response
        )

        to_land, coast_strength = self._coastal_direction(grid, ocean, ~ocean)
        land_east = np.sum(to_land * east, axis=1)
        warm_edge = ocean & (land_east < -0.18)
        cold_edge = ocean & (land_east > 0.18)
        boundary_lat = np.clip((lat_abs - 4.0) / 24.0, 0.0, 1.0) * subtropics
        warm_strength = coast_strength * warm_edge.astype(float) * boundary_lat * basin_scale
        cold_strength = coast_strength * cold_edge.astype(float) * boundary_lat * basin_scale

        along = self._project_tangent(grid, np.cross(grid.xyz, to_land))
        along_norm = np.linalg.norm(along, axis=1, keepdims=True)
        along = np.where(along_norm > 1e-9, along / np.maximum(along_norm, 1e-9), 0.0)
        hemisphere = np.where(grid.lat >= 0.0, 1.0, -1.0)
        polarward = hemisphere[:, None] * north
        equatorward = -polarward

        warm_sign = np.sign(np.sum(along * polarward, axis=1))
        cold_sign = np.sign(np.sum(along * equatorward, axis=1))
        warm_sign = np.where(warm_sign == 0.0, 1.0, warm_sign)
        cold_sign = np.where(cold_sign == 0.0, 1.0, cold_sign)
        warm_along = warm_sign[:, None] * along
        cold_along = cold_sign[:, None] * along
        boundary_return = np.clip(np.linalg.norm(stream_flow, axis=1), 0.0, 1.4)
        currents += 0.46 * warm_strength[:, None] * (0.55 + 0.45 * boundary_return)[:, None] * warm_along
        currents += 0.38 * cold_strength[:, None] * (0.55 + 0.45 * boundary_return)[:, None] * cold_along
        boundary_influence = np.clip(
            1.25 * (warm_strength + cold_strength)
            + 0.22 * shelf_index
            + 0.20 * strait_index,
            0.0,
            1.0,
        )
        equatorial_current_core = np.exp(-(grid.lat / 10.0) ** 2)
        open_ocean_damping = (
            ocean.astype(float)
            * (1.0 - boundary_influence)
            * np.clip(1.0 - shelf_index, 0.0, 1.0)
            * (1.0 - 0.45 * equatorial_current_core)
        )
        currents *= np.clip(1.0 - 0.10 * open_ocean_damping, 0.88, 1.0)[:, None]
        currents = self._project_tangent(grid, currents)
        currents[~ocean] = 0.0
        currents = self._smooth_vector_field_masked(grid, currents, ocean,
                                                    passes=5, alpha=0.18)
        speed = np.linalg.norm(currents, axis=1)
        currents *= np.minimum(1.0, 1.25 / np.maximum(speed, 1e-9))[:, None]
        currents[~ocean] = 0.0

        cold_wind = np.maximum(np.sum(annual_wind * cold_along, axis=1) / 9.0, 0.0)
        curl_upwelling = np.clip(np.maximum(-streamfunction * np.sign(grid.lat), 0.0),
                                 0.0, 1.0)
        upwelling = np.clip(
            cold_strength * (0.62 + 0.30 * cold_wind + 0.25 * curl_upwelling),
            0.0,
            1.0,
        )
        upwelling = self._smooth_field_masked(grid, upwelling, ocean,
                                              passes=2, alpha=0.22)

        poleward_flow = np.sum(currents * polarward, axis=1)
        gyre_heat = 3.8 * np.tanh(poleward_flow / 0.28) * subtropics * ocean.astype(float)
        heat_ocean = (
            gyre_heat
            + 8.6 * warm_strength
            - 8.2 * cold_strength
            - 3.0 * upwelling
        )
        heat_ocean *= basin_scale
        heat_ocean[~ocean] = 0.0
        heat_ocean = self._smooth_field_masked(grid, heat_ocean, ocean,
                                               passes=4, alpha=0.24)
        if ocean.any():
            heat_ocean[ocean] -= float(np.average(
                heat_ocean[ocean], weights=grid.cell_area[ocean]))
        heat_ocean = np.where(ocean, np.clip(heat_ocean, -7.5, 7.5), 0.0)
        sst_anomaly = heat_ocean.copy()

        heat_effect = self._spread_ocean_influence_to_coasts(
            grid, heat_ocean, ocean, passes=3, alpha=0.42, land_damping=0.86)
        heat_effect -= float(np.average(heat_effect, weights=grid.cell_area))
        # Boundary currents should leave a visible same-latitude warm/cold
        # coastal imprint on adjacent land.  The ocean anomaly itself remains
        # controlled by the gyre/upwelling solve; this only strengthens the
        # exported coastal land expression after cross-shore spreading.
        heat_effect[~ocean] *= 1.50
        heat_effect = np.clip(heat_effect, -6.0, 6.0)
        # NOAA/AOML drifter calibration: the reduced gyre machinery above uses
        # a dynamically useful transport proxy, but the exported near-surface
        # current vectors should stay in observed drifter-speed ranges.  Keep
        # the diagnostic wind-stress response slightly stronger so the ocean
        # coupling layer remains visible without raising final current speeds.
        surface_currents = 0.64 * currents
        wind_stress_current_response = 0.80 * direct_wind_response
        boundary_current_type = np.zeros(grid.n, dtype=np.float64)
        boundary_current_type[warm_strength > 0.08] = 1.0
        boundary_current_type[cold_strength > 0.08] = -1.0
        speed = np.linalg.norm(surface_currents, axis=1)
        strait_exchange = np.where(ocean, np.clip(strait_index * speed / 0.35, 0.0, 1.0), 0.0)
        return {
            "currents": surface_currents,
            "current_heat_transport": heat_effect,
            "upwelling": upwelling,
            "basin_id": basin_id,
            "current_streamfunction": streamfunction,
            "gyre_id": gyre_id,
            "boundary_current_type": boundary_current_type,
            "strait_exchange": strait_exchange,
            "wind_stress_current_response": wind_stress_current_response,
            "sst_anomaly": sst_anomaly,
        }

    def _weak_ocean_atmosphere_coupling(
        self,
        world,
        seasonal_T,
        ocean,
        seasonal_wind,
        pressure_proxy,
        ocean_current,
        geography_fields,
    ):
        grid = world.grid
        ocean = np.asarray(ocean, dtype=bool)
        base_wind = np.asarray(seasonal_wind, dtype=np.float64)
        base_pressure = np.asarray(pressure_proxy, dtype=np.float64)
        wind = base_wind.copy()
        pressure = base_pressure.copy()
        current = dict(ocean_current)
        residual = np.zeros(grid.n, dtype=np.float64)
        ocean_heat_flux = np.zeros(grid.n, dtype=np.float64)
        evaporation_feedback = np.zeros(grid.n, dtype=np.float64)
        seasonal_sst = np.asarray(seasonal_T, dtype=np.float64).copy()

        for _ in range(self.OCEAN_ATMOSPHERE_ITERS):
            previous_sst = np.asarray(current["sst_anomaly"], dtype=np.float64)
            previous_flux = ocean_heat_flux.copy()
            upwelling = np.asarray(current["upwelling"], dtype=np.float64)
            ocean_heat_flux = np.where(
                ocean,
                np.clip(previous_sst - 1.15 * upwelling, -6.0, 6.0),
                0.0,
            )
            if ocean.any():
                ocean_heat_flux[ocean] -= float(np.average(
                    ocean_heat_flux[ocean], weights=grid.cell_area[ocean]))
            ocean_heat_flux = np.where(ocean, np.clip(ocean_heat_flux, -6.0, 6.0), 0.0)
            seasonal_sst = np.where(
                ocean[None, :],
                np.maximum(seasonal_T + ocean_heat_flux[None, :],
                           CONSTANTS.ZERO_C - 1.8),
                seasonal_T,
            )
            evaporation_feedback = self._ocean_evaporation_heat_feedback(
                grid, seasonal_sst, ocean, current["current_heat_transport"],
                upwelling, ocean_heat_flux)
            ocean_heat_flux = np.where(
                ocean,
                np.clip(ocean_heat_flux + evaporation_feedback, -6.0, 6.0),
                0.0,
            )
            if ocean.any():
                ocean_heat_flux[ocean] -= float(np.average(
                    ocean_heat_flux[ocean], weights=grid.cell_area[ocean]))
            ocean_heat_flux = np.where(ocean, np.clip(ocean_heat_flux, -6.0, 6.0), 0.0)
            seasonal_sst = np.where(
                ocean[None, :],
                np.maximum(seasonal_T + ocean_heat_flux[None, :],
                           CONSTANTS.ZERO_C - 1.8),
                seasonal_T,
            )

            warm = np.where(ocean, np.clip(ocean_heat_flux / 4.0, 0.0, 1.0), 0.0)
            cold = np.where(
                ocean,
                np.clip(np.maximum(-ocean_heat_flux / 3.2, 1.15 * upwelling), 0.0, 1.0),
                0.0,
            )
            evap_low = np.where(
                ocean, np.clip(-evaporation_feedback / 0.75, 0.0, 1.0), 0.0)
            pressure_source = np.where(
                ocean,
                -0.18 * warm - 0.035 * evap_low + 0.12 * cold,
                0.0,
            )
            pressure_source = self._spread_ocean_influence_to_coasts(
                grid, pressure_source, ocean, passes=3, alpha=0.30, land_damping=0.52)

            pressure_next = np.zeros_like(base_pressure)
            wind_next = np.zeros_like(base_wind)
            for s in range(4):
                pressure_next[s] = self._smooth_field(
                    grid, base_pressure[s] + pressure_source, passes=3, alpha=0.14)
                pressure_next[s] = np.clip(pressure_next[s], -1.8, 1.8)

                high_to_low = -self._graph_gradient_vectors(grid, pressure_next[s])
                mag = np.linalg.norm(high_to_low, axis=1)
                p95 = max(float(np.percentile(mag, 95)), 1.0e-9)
                anomaly = 0.42 * high_to_low / p95
                anomaly = self._smooth_vector_field(grid, anomaly, passes=2, alpha=0.16)
                speed = np.linalg.norm(anomaly, axis=1)
                anomaly *= np.minimum(1.0, 1.6 / np.maximum(speed, 1.0e-9))[:, None]
                wind_next[s] = self._project_tangent(grid, base_wind[s] + anomaly)
                wind_speed = np.linalg.norm(wind_next[s], axis=1)
                wind_next[s] *= np.minimum(
                    1.0, 24.0 / np.maximum(wind_speed, 1.0e-9))[:, None]

            next_current = self._ocean_currents(
                world, ocean, wind_next, geography_fields)
            residual = np.where(
                ocean,
                np.maximum(
                    np.abs(np.asarray(next_current["sst_anomaly"], dtype=np.float64)
                           - previous_sst),
                    0.65 * np.abs(ocean_heat_flux - previous_flux),
                ),
                0.0,
            )
            pressure = pressure_next
            wind = wind_next
            current = next_current

        final_sst = np.asarray(current["sst_anomaly"], dtype=np.float64)
        final_upwelling = np.asarray(current["upwelling"], dtype=np.float64)
        ocean_heat_flux = np.where(
            ocean,
            np.clip(final_sst - 1.15 * final_upwelling, -6.0, 6.0),
            0.0,
        )
        if ocean.any():
            ocean_heat_flux[ocean] -= float(np.average(
                ocean_heat_flux[ocean], weights=grid.cell_area[ocean]))
        ocean_heat_flux = np.where(ocean, np.clip(ocean_heat_flux, -6.0, 6.0), 0.0)
        seasonal_sst = np.where(
            ocean[None, :],
            np.maximum(seasonal_T + ocean_heat_flux[None, :],
                       CONSTANTS.ZERO_C - 1.8),
            seasonal_T,
        )
        evaporation_feedback = self._ocean_evaporation_heat_feedback(
            grid, seasonal_sst, ocean, current["current_heat_transport"],
            final_upwelling, ocean_heat_flux)
        ocean_heat_flux = np.where(
            ocean,
            np.clip(ocean_heat_flux + evaporation_feedback, -6.0, 6.0),
            0.0,
        )
        if ocean.any():
            ocean_heat_flux[ocean] -= float(np.average(
                ocean_heat_flux[ocean], weights=grid.cell_area[ocean]))
        ocean_heat_flux = np.where(ocean, np.clip(ocean_heat_flux, -6.0, 6.0), 0.0)
        seasonal_sst = np.where(
            ocean[None, :],
            np.maximum(seasonal_T + ocean_heat_flux[None, :],
                       CONSTANTS.ZERO_C - 1.8),
            seasonal_T,
        )
        return {
            "seasonal_wind": wind,
            "pressure_proxy": pressure,
            "ocean_current": current,
            "seasonal_sst": seasonal_sst,
            "ocean_heat_flux": ocean_heat_flux,
            "evaporation_feedback": evaporation_feedback,
            "coupling_residual": np.where(ocean, np.clip(residual, 0.0, 8.0), 0.0),
        }

    def _ocean_evaporation_heat_feedback(self, grid, seasonal_sst, ocean,
                                         current_heat, upwelling,
                                         ocean_heat_flux):
        ocean = np.asarray(ocean, dtype=bool)
        feedback = np.zeros(grid.n, dtype=np.float64)
        if not ocean.any():
            return feedback
        seasonal_sst = np.asarray(seasonal_sst, dtype=np.float64)
        heat = np.asarray(current_heat, dtype=np.float64)
        upwelling = np.asarray(upwelling, dtype=np.float64)
        ocean_heat_flux = np.asarray(ocean_heat_flux, dtype=np.float64)
        seasonal_evap = np.zeros((4, grid.n), dtype=np.float64)
        for s in range(4):
            seasonal_evap[s] = self._seasonal_evaporation(
                seasonal_sst[s], ocean, heat, upwelling, seasonal_sst[s],
                ocean_heat_flux)
        evap = seasonal_evap.mean(axis=0)
        e10, e50, e90 = np.percentile(evap[ocean], [10, 50, 90])
        scale = max(e90 - e10, 80.0)
        evap_anom = np.zeros(grid.n, dtype=np.float64)
        evap_anom[ocean] = (evap[ocean] - e50) / scale
        latent_cooling = -0.55 * evap_anom
        cold_suppression = np.where(
            ocean,
            0.18 * np.clip(np.maximum(-heat / 4.0, upwelling), 0.0, 1.0),
            0.0,
        )
        feedback = np.where(ocean, latent_cooling + cold_suppression, 0.0)
        feedback = self._smooth_field_masked(
            grid, feedback, ocean, passes=3, alpha=0.20)
        feedback = np.where(ocean, np.clip(feedback, -0.85, 0.65), 0.0)
        feedback[ocean] -= float(np.average(
            feedback[ocean], weights=grid.cell_area[ocean]))
        feedback = np.where(ocean, np.clip(feedback, -0.90, 0.90), 0.0)
        feedback[ocean] -= float(np.average(
            feedback[ocean], weights=grid.cell_area[ocean]))
        return np.where(ocean, feedback, 0.0)

    def _apply_ocean_current_hydro_adjustment(self, grid, evap, precip, ocean,
                                              heat_transport, upwelling):
        cold_ocean = np.where(
            ocean,
            np.clip(np.maximum(-heat_transport / 3.0, 1.6 * upwelling), 0.0, 1.0),
            0.0,
        )
        warm_ocean = np.where(ocean, np.clip(heat_transport / 5.0, 0.0, 1.0), 0.0)
        cold_influence = np.clip(self._spread_ocean_influence_to_coasts(
            grid, cold_ocean, ocean, passes=4, alpha=0.25, land_damping=0.68),
            0.0, 1.0)
        evap = np.where(ocean, evap * (1.0 + 0.08 * warm_ocean), evap)
        evap = np.where(ocean, evap * (1.0 - 0.55 * cold_ocean), evap)
        precip = precip * (1.0 - 0.16 * cold_influence)
        return np.clip(evap, 0.0, 4000.0), np.clip(precip, 0.0, 4000.0)

    def _seasonal_evaporation(self, T, ocean, current_heat, upwelling,
                              sst=None, ocean_heat_flux=None):
        source_T = np.asarray(T, dtype=np.float64)
        if sst is not None:
            source_T = np.where(ocean, np.asarray(sst, dtype=np.float64), source_T)
        heat = np.asarray(current_heat, dtype=np.float64)
        if ocean_heat_flux is not None:
            heat = np.where(
                ocean,
                0.60 * heat + 0.40 * np.asarray(ocean_heat_flux, dtype=np.float64),
                heat,
            )
        es = 6.11 * np.exp(17.27 * (source_T - 273.15)
                            / np.maximum(source_T - 35.85, 1.0))
        evap = np.where(ocean, np.clip(1.85 * es, 0.0, 4200.0), 0.0)
        warm = np.where(ocean, np.clip(heat / 5.0, 0.0, 1.0), 0.0)
        cold = np.where(
            ocean,
            np.clip(np.maximum(-heat / 3.5, 1.55 * upwelling), 0.0, 1.0),
            0.0,
        )
        ice = self._ice_fraction(source_T, ocean)
        evap *= 1.0 + 0.10 * warm
        evap *= 1.0 - 0.55 * cold
        evap *= 1.0 - 0.82 * ice
        return np.clip(evap, 0.0, 4200.0)

    def _advective_moisture_access(self, grid, source, wind, ocean,
                                   terrain_blocking):
        source = np.asarray(source, dtype=np.float64)
        wind = np.asarray(wind, dtype=np.float64)
        ocean = np.asarray(ocean, dtype=bool)
        terrain_blocking = np.asarray(terrain_blocking, dtype=np.float64)
        i, j = self._edges[:, 0], self._edges[:, 1]
        dpos = grid.xyz[j] - grid.xyz[i]
        flow = np.einsum("ij,ij->i", 0.5 * (wind[i] + wind[j]), dpos)
        pos = np.clip(flow, 0.0, None)
        neg = np.clip(-flow, 0.0, None)
        if float(np.max(pos, initial=0.0) + np.max(neg, initial=0.0)) <= 1e-12:
            pos = neg = np.ones_like(flow)

        access = np.clip(source, 0.0, 1.4)
        for _ in range(54):
            acc = np.zeros(grid.n, dtype=np.float64)
            weight = np.zeros(grid.n, dtype=np.float64)
            if pos.any():
                np.add.at(acc, j, access[i] * pos)
                np.add.at(weight, j, pos)
            if neg.any():
                np.add.at(acc, i, access[j] * neg)
                np.add.at(weight, i, neg)
            adv = acc / np.maximum(weight, 1e-9)
            mixed = 0.62 * adv + 0.38 * self._neighbor_mean(grid, access)
            access = 0.67 * access + 0.33 * mixed
            access[ocean] = np.maximum(access[ocean], source[ocean])
            # Climate-scale moisture is not just near-surface cell-to-cell
            # advection.  A weak free-atmosphere reservoir keeps Earthlike
            # interiors from becoming all-desert while mountain barriers still
            # remove moisture preferentially.
            access[~ocean] *= 1.0 - 0.072 * terrain_blocking[~ocean] ** 1.25
            access = np.clip(access, 0.0, 1.5)
        broad = self._smooth_field(grid, access, passes=7, alpha=0.18)
        land = ~ocean
        broad_tail = 0.34 * broad * (1.0 - 0.45 * terrain_blocking)
        access[land] = np.maximum(access[land], broad_tail[land])
        return np.clip(access, 0.0, 1.0)

    def _seasonal_pressure_moisture(self, world, seasonal_T, ocean, seasonal_wind,
                                    pressure_proxy, current_heat, upwelling,
                                    geography_fields, seasonal_sst=None,
                                    ocean_heat_flux=None):
        grid = world.grid
        land = ~ocean
        zero = np.zeros(grid.n, dtype=np.float64)
        barrier = np.asarray(
            geography_fields.get("terrain.barrier_index", zero), dtype=np.float64)
        wind_gap = np.asarray(
            geography_fields.get("terrain.wind_gap_index", zero), dtype=np.float64)
        interiority = np.asarray(
            geography_fields.get("climate.continent_interiority", zero),
            dtype=np.float64,
        )
        terrain_blocking = np.where(
            land, np.clip(barrier * (1.0 - 0.70 * wind_gap), 0.0, 1.0), 0.0)

        source_ocean_warmth = np.zeros((4, grid.n), dtype=np.float64)
        seasonal_pressure = np.asarray(pressure_proxy, dtype=np.float64).copy()
        moisture_access = np.zeros((4, grid.n), dtype=np.float64)
        monsoon_potential = np.zeros((4, grid.n), dtype=np.float64)
        temp_anom = seasonal_T - seasonal_T.mean(axis=0, keepdims=True)
        forward_cooling = temp_anom - np.roll(temp_anom, shift=-1, axis=0)
        phase_weight = np.where(land, 0.30 + 0.95 * interiority, 0.08)
        phase_pressure = forward_cooling * phase_weight[None, :]
        phase_scale = max(float(np.percentile(np.abs(phase_pressure), 95)), 1.0e-9)
        phase_pressure = 0.30 * phase_pressure / phase_scale
        for s in range(4):
            phase_pressure[s] = self._smooth_field(
                grid, phase_pressure[s], passes=4, alpha=0.20)

        heat_source = np.asarray(current_heat, dtype=np.float64)
        if ocean_heat_flux is not None:
            heat_source = np.where(
                ocean,
                0.58 * heat_source + 0.42 * np.asarray(ocean_heat_flux, dtype=np.float64),
                heat_source,
            )
        sst = np.asarray(seasonal_sst, dtype=np.float64) if seasonal_sst is not None else None

        warm_current = np.where(ocean, np.clip(heat_source / 5.0, 0.0, 1.0), 0.0)
        cold_current = np.where(
            ocean,
            np.clip(np.maximum(-heat_source / 3.5, 1.35 * upwelling), 0.0, 1.0),
            0.0,
        )

        for s in range(4):
            temp_c = seasonal_T[s] - 273.15
            source_T = seasonal_T[s] if sst is None else np.where(ocean, sst[s], seasonal_T[s])
            source_temp_c = source_T - 273.15
            ocean_source = np.where(
                ocean,
                np.clip((source_temp_c + 2.0) / 31.0, 0.0, 1.20)
                * (1.0 + 0.18 * warm_current)
                * (1.0 - 0.62 * cold_current)
                * (1.0 - 0.78 * self._ice_fraction(source_T, ocean)),
                0.0,
            )
            ocean_source = self._smooth_field_masked(
                grid, ocean_source, ocean, passes=2, alpha=0.20)
            source_ocean_warmth[s] = np.clip(ocean_source, 0.0, 1.0)

            ocean_pressure = np.where(
                ocean,
                -0.34 * source_ocean_warmth[s] + 0.16 * cold_current,
                0.0,
            )
            seasonal_pressure[s] = self._smooth_field(
                grid,
                seasonal_pressure[s] + phase_pressure[s] + ocean_pressure,
                passes=5,
                alpha=0.20)
            seasonal_pressure[s] = np.clip(seasonal_pressure[s], -1.8, 1.8)

            access = self._advective_moisture_access(
                grid, source_ocean_warmth[s], seasonal_wind[s], ocean,
                terrain_blocking)
            moisture_access[s] = access

            thermal_low = np.clip(-seasonal_pressure[s], 0.0, 1.6)
            cold_high = np.clip(seasonal_pressure[s], 0.0, 1.6)
            land_warmth = np.clip((temp_c + 4.0) / 28.0, 0.0, 1.0)
            moisture_gate = np.clip((access - 0.06) / 0.52, 0.0, 1.0)
            access_for_monsoon = np.clip(
                (0.72 * access + 0.28 * np.sqrt(access))
                * (0.20 + 0.80 * moisture_gate),
                0.0,
                1.2,
            )
            potential = (
                1.32 * thermal_low * access_for_monsoon * (0.30 + 0.70 * interiority)
                * (0.55 + 0.45 * land_warmth)
            )
            winter_export = 0.20 * cold_high * (0.35 + 0.65 * interiority)
            monsoon = np.where(land, potential - winter_export, 0.0)
            monsoon = self._smooth_field_masked(grid, monsoon, land,
                                                passes=2, alpha=0.22)
            monsoon_potential[s] = np.clip(monsoon, -0.8, 1.4)

        return {
            "seasonal_pressure": seasonal_pressure,
            "source_ocean_warmth": source_ocean_warmth,
            "terrain_blocking": terrain_blocking,
            "moisture_access": moisture_access,
            "monsoon_potential": monsoon_potential,
        }

    def _hydroclimate_pressure_wind_feedback(self, world, hydro, seasonal_T,
                                             ocean, seasonal_wind,
                                             pressure_proxy):
        grid = world.grid
        land = ~np.asarray(ocean, dtype=bool)
        seasonal_wind = np.asarray(seasonal_wind, dtype=np.float64)
        pressure_proxy = np.asarray(pressure_proxy, dtype=np.float64)
        feedback = np.zeros((4, grid.n), dtype=np.float64)
        wind_anomaly = np.zeros((4, grid.n, 3), dtype=np.float64)
        adjusted_pressure = pressure_proxy.copy()
        adjusted_wind = seasonal_wind.copy()
        if not land.any():
            return {
                "pressure_proxy": adjusted_pressure,
                "seasonal_wind": adjusted_wind,
                "pressure_feedback": feedback,
                "wind_anomaly": wind_anomaly,
                "residual": np.zeros(grid.n, dtype=np.float64),
            }
        # Tiny island worlds can have intense local rainfall, but that should
        # not drive a planet-scale pressure/wind feedback in this reduced
        # model.  Scale the hydro feedback by exposed-land extent; Earthlike
        # and arid continental cases are unaffected.
        land_fraction = float(np.mean(land))
        land_feedback_scale = np.clip((land_fraction - 0.025) / 0.18, 0.0, 1.0) ** 2

        seasonal_precip = np.asarray(
            hydro.get("seasonal_precipitation", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        moisture = np.asarray(
            hydro.get("moisture_access", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        monsoon = np.asarray(
            hydro.get("monsoon_potential", np.zeros((4, grid.n))),
            dtype=np.float64,
        )
        lat_abs = np.abs(grid.lat)
        subtropical_dry = np.exp(-((lat_abs - 28.0) / 13.0) ** 2)
        high_lat_stability = np.clip((lat_abs - 48.0) / 32.0, 0.0, 1.0)

        for s in range(4):
            land_precip = seasonal_precip[s, land]
            if land_precip.size < 4:
                continue
            p10, p35, p50, p85 = np.percentile(land_precip, [10, 35, 50, 85])
            wet = np.clip(
                (seasonal_precip[s] - p50) / max(p85 - p50, 1.0), 0.0, 1.0)
            dry = np.clip(
                (p35 - seasonal_precip[s]) / max(p35 - p10, 1.0), 0.0, 1.0)
            warm = np.clip((seasonal_T[s] - 276.0) / 22.0, 0.0, 1.0)
            low_moisture = np.clip(1.0 - moisture[s], 0.0, 1.0)
            wet_support = np.clip(
                0.62 * moisture[s] + 0.38 * np.clip(monsoon[s], 0.0, 1.2),
                0.0,
                1.25,
            )
            convective_low = wet * warm * wet_support
            subsidence_high = dry * low_moisture * (
                0.62 + 0.38 * np.maximum(subtropical_dry, high_lat_stability))
            raw = np.where(
                land,
                -0.055 * convective_low + 0.038 * subsidence_high,
                0.0,
            )
            raw *= land_feedback_scale
            raw = self._smooth_field_masked(grid, raw, land, passes=4, alpha=0.18)
            raw = np.where(land, np.clip(raw, -0.075, 0.055), 0.0)
            feedback[s] = raw
            adjusted_pressure[s] = np.clip(pressure_proxy[s] + raw, -1.8, 1.8)

            high_to_low = -self._graph_gradient_vectors(grid, raw)
            mag = np.linalg.norm(high_to_low, axis=1)
            p95 = float(np.percentile(mag, 95)) if mag.size else 0.0
            if p95 > 1.0e-12:
                anomaly = 0.30 * high_to_low / p95
                anomaly *= np.clip(np.abs(raw) / 0.060, 0.0, 1.0)[:, None]
                anomaly = self._smooth_vector_field(grid, anomaly, passes=2, alpha=0.16)
                speed = np.linalg.norm(anomaly, axis=1)
                anomaly *= np.minimum(
                    1.0, 0.65 / np.maximum(speed, 1.0e-9))[:, None]
            else:
                anomaly = np.zeros((grid.n, 3), dtype=np.float64)
            wind_anomaly[s] = anomaly
            adjusted_wind[s] = self._project_tangent(grid, seasonal_wind[s] + anomaly)
            wind_speed = np.linalg.norm(adjusted_wind[s], axis=1)
            adjusted_wind[s] *= np.minimum(
                1.0, 24.0 / np.maximum(wind_speed, 1.0e-9))[:, None]

        residual = np.mean(np.abs(feedback), axis=0)
        return {
            "pressure_proxy": adjusted_pressure,
            "seasonal_wind": adjusted_wind,
            "pressure_feedback": feedback,
            "wind_anomaly": wind_anomaly,
            "residual": np.where(land, np.clip(residual, 0.0, 0.20), 0.0),
        }

    def _seasonal_hydroclimate_feedback_loop(
        self,
        world,
        seasonal_T,
        ocean,
        seasonal_wind,
        elev,
        current_heat,
        upwelling,
        itcz_intensity,
        storm_tracks,
        pressure_proxy,
        geography_fields,
        seasonal_sst=None,
        ocean_heat_flux=None,
    ):
        grid = world.grid
        base_wind = np.asarray(seasonal_wind, dtype=np.float64)
        base_pressure = np.asarray(pressure_proxy, dtype=np.float64)
        current_wind = base_wind.copy()
        current_pressure = base_pressure.copy()
        prev_feedback = np.zeros((4, grid.n), dtype=np.float64)
        iteration_delta = np.zeros(grid.n, dtype=np.float64)
        feedback_state = {
            "pressure_proxy": current_pressure,
            "seasonal_wind": current_wind,
            "pressure_feedback": prev_feedback,
            "wind_anomaly": np.zeros((4, grid.n, 3), dtype=np.float64),
            "residual": np.zeros(grid.n, dtype=np.float64),
        }

        for iteration in range(self.HYDRO_FEEDBACK_ITERS):
            trial_hydro = self._seasonal_hydroclimate(
                world, seasonal_T, ocean, current_wind, elev, current_heat,
                upwelling, itcz_intensity, storm_tracks, current_pressure,
                geography_fields, seasonal_sst, ocean_heat_flux)
            feedback_state = self._hydroclimate_pressure_wind_feedback(
                world, trial_hydro, seasonal_T, ocean, base_wind, base_pressure)
            feedback = np.asarray(
                feedback_state["pressure_feedback"], dtype=np.float64)
            if iteration > 0:
                iteration_delta = np.mean(np.abs(feedback - prev_feedback), axis=0)
            prev_feedback = feedback

            blend = 1.0 if iteration == 0 else 0.62
            target_pressure = np.asarray(
                feedback_state["pressure_proxy"], dtype=np.float64)
            target_wind = np.asarray(
                feedback_state["seasonal_wind"], dtype=np.float64)
            current_pressure = np.clip(
                (1.0 - blend) * current_pressure + blend * target_pressure,
                -1.8,
                1.8,
            )
            current_wind = (
                (1.0 - blend) * current_wind + blend * target_wind
            )
            for s in range(4):
                current_wind[s] = self._project_tangent(grid, current_wind[s])
                speed = np.linalg.norm(current_wind[s], axis=1)
                current_wind[s] *= np.minimum(
                    1.0, 24.0 / np.maximum(speed, 1.0e-9))[:, None]

        final_hydro = self._seasonal_hydroclimate(
            world, seasonal_T, ocean, current_wind, elev, current_heat,
            upwelling, itcz_intensity, storm_tracks, current_pressure,
            geography_fields, seasonal_sst, ocean_heat_flux)
        final_feedback = self._hydroclimate_pressure_wind_feedback(
            world, final_hydro, seasonal_T, ocean, base_wind, base_pressure)
        final_feedback_arr = np.asarray(
            final_feedback["pressure_feedback"], dtype=np.float64)
        final_delta = np.mean(np.abs(final_feedback_arr - prev_feedback), axis=0)
        iteration_delta = np.maximum(iteration_delta, final_delta)
        # Use the feedback implied by the final hydro solve for diagnostics and
        # exported pressure/wind.  The extra solve above makes this a bounded
        # fixed-point step instead of a one-shot adjustment.
        final_pressure = np.asarray(final_feedback["pressure_proxy"], dtype=np.float64)
        final_wind = np.asarray(final_feedback["seasonal_wind"], dtype=np.float64)
        final_hydro = self._seasonal_hydroclimate(
            world, seasonal_T, ocean, final_wind, elev, current_heat,
            upwelling, itcz_intensity, storm_tracks, final_pressure,
            geography_fields, seasonal_sst, ocean_heat_flux)
        final_feedback["iteration_delta"] = np.where(
            ~np.asarray(ocean, dtype=bool),
            np.clip(iteration_delta, 0.0, 0.20),
            0.0,
        )
        final_feedback["iteration_count"] = int(self.HYDRO_FEEDBACK_ITERS)
        return final_hydro, final_feedback

    def _orographic_precip_terms(self, grid, wind, elev, land, barrier,
                                 wind_gap):
        i, j = self._edges[:, 0], self._edges[:, 1]
        dpos = grid.xyz[j] - grid.xyz[i]
        wmid = 0.5 * (wind[i] + wind[j])
        flow = np.einsum("ij,ij->i", wmid, dpos)
        dh = elev[j] - elev[i]
        up = flow * dh
        oro = np.zeros(grid.n, dtype=np.float64)
        lee = np.zeros(grid.n, dtype=np.float64)
        np.add.at(oro, j, np.maximum(up, 0.0))
        np.add.at(oro, i, np.maximum(-up, 0.0))
        np.add.at(lee, i, np.maximum(up, 0.0))
        np.add.at(lee, j, np.maximum(-up, 0.0))
        oro = self._smooth_field_masked(grid, oro, land, passes=2, alpha=0.22)
        lee = self._smooth_field_masked(grid, lee, land, passes=2, alpha=0.22)
        oro_scale = max(float(np.percentile(oro[land], 96)) if land.any() else 0.0,
                        1e-9)
        lee_scale = max(float(np.percentile(lee[land], 96)) if land.any() else 0.0,
                        1e-9)
        barrier_weight = np.clip(0.35 + 0.65 * barrier, 0.0, 1.0)
        pass_relief = np.clip(1.0 - 0.55 * wind_gap, 0.25, 1.0)
        oro_norm = np.where(
            land, np.clip(oro / oro_scale, 0.0, 1.0) * barrier_weight, 0.0)
        lee_norm = np.where(
            land, np.clip(lee / lee_scale, 0.0, 1.0) * barrier_weight * pass_relief,
            0.0,
        )
        return oro_norm, lee_norm

    def _slope_wind_exposure(self, grid, wind, elev, land, barrier, wind_gap):
        topo = self._smooth_field(grid, np.maximum(elev, 0.0), passes=2, alpha=0.18)
        topo_scale = max(float(np.percentile(topo[land], 95)) if land.any() else 1.0,
                         1.0)
        topo_norm = np.clip(topo / topo_scale, 0.0, 2.0)
        topo_grad = self._graph_gradient_vectors(grid, topo_norm)
        topo_mag = np.linalg.norm(topo_grad, axis=1)
        if land.any():
            g60 = float(np.percentile(topo_mag[land], 60))
            g95 = max(float(np.percentile(topo_mag[land], 95)), g60 + 1.0e-9)
        else:
            g60, g95 = 0.0, 1.0
        slope_strength = np.clip((topo_mag - g60) / (g95 - g60), 0.0, 1.0)
        slope_strength *= np.clip(0.38 + 0.72 * barrier, 0.0, 1.0)
        slope_strength = np.where(land, slope_strength, 0.0)
        wind_speed = np.linalg.norm(wind, axis=1)
        wind_unit = np.where(
            wind_speed[:, None] > 1.0e-9,
            wind / np.maximum(wind_speed[:, None], 1.0e-9),
            0.0,
        )
        grad_unit = np.where(
            topo_mag[:, None] > 1.0e-9,
            topo_grad / np.maximum(topo_mag[:, None], 1.0e-9),
            0.0,
        )
        wind_uphill = np.sum(wind_unit * grad_unit, axis=1)
        gap_relief = np.clip(1.0 - 0.45 * wind_gap, 0.25, 1.0)
        upslope = slope_strength * np.clip((wind_uphill - 0.06) / 0.54, 0.0, 1.0)
        downslope = (
            slope_strength
            * np.clip((-wind_uphill - 0.06) / 0.54, 0.0, 1.0)
            * gap_relief
        )
        upslope = self._smooth_field_masked(grid, upslope, land, passes=1, alpha=0.12)
        downslope = self._smooth_field_masked(grid, downslope, land,
                                              passes=1, alpha=0.12)
        return np.clip(upslope, 0.0, 1.0), np.clip(downslope, 0.0, 1.0)

    def _accentuate_low_latitude_seasonality(self, grid, seasonal_precip, land):
        seasonal_precip = np.asarray(seasonal_precip, dtype=np.float64)
        land = np.asarray(land, dtype=bool)
        annual = seasonal_precip.mean(axis=0)
        ratio = seasonal_precip / np.maximum(annual[None, :], 1.0e-9)
        lat_focus = np.exp(-((np.abs(grid.lat) - 12.0) / 14.0) ** 2)
        exponent = 1.0 + 0.35 * lat_focus * land.astype(np.float64)
        shaped = ratio ** exponent[None, :]
        shaped /= np.maximum(shaped.mean(axis=0, keepdims=True), 1.0e-9)
        return annual[None, :] * shaped

    def _seasonal_hydroclimate(self, world, seasonal_T, ocean, seasonal_wind,
                               elev, current_heat, upwelling, itcz_intensity,
                               storm_tracks, pressure_proxy, geography_fields,
                               seasonal_sst=None, ocean_heat_flux=None):
        grid = world.grid
        land = ~ocean
        zero = np.zeros(grid.n, dtype=np.float64)
        barrier = np.asarray(
            geography_fields.get("terrain.barrier_index", zero), dtype=np.float64)
        wind_gap = np.asarray(
            geography_fields.get("terrain.wind_gap_index", zero), dtype=np.float64)
        interiority = np.asarray(
            geography_fields.get("climate.continent_interiority", zero),
            dtype=np.float64,
        )
        pressure_moisture = self._seasonal_pressure_moisture(
            world, seasonal_T, ocean, seasonal_wind, pressure_proxy,
            current_heat, upwelling, geography_fields, seasonal_sst,
            ocean_heat_flux)
        moisture_access = pressure_moisture["moisture_access"]
        monsoon_potential = pressure_moisture["monsoon_potential"]
        source_ocean_warmth = pressure_moisture["source_ocean_warmth"]
        seasonal_pressure = pressure_moisture["seasonal_pressure"]
        coast_orientation = np.asarray(
            geography_fields.get(
                "climate.coast_orientation",
                np.zeros((grid.n, 3), dtype=np.float64),
            ),
            dtype=np.float64,
        )
        if coast_orientation.shape != (grid.n, 3):
            coast_orientation = np.zeros((grid.n, 3), dtype=np.float64)
        coast_strength = np.asarray(
            geography_fields.get("climate.coast_strength", zero), dtype=np.float64)
        if coast_strength.shape != (grid.n,):
            coast_strength = zero

        seasonal_precip = np.zeros((4, grid.n), dtype=np.float64)
        seasonal_evap = np.zeros((4, grid.n), dtype=np.float64)
        moisture_convergence = np.zeros((4, grid.n), dtype=np.float64)
        orographic_component = np.zeros((4, grid.n), dtype=np.float64)
        monsoon_corridor = np.zeros((4, grid.n), dtype=np.float64)
        storm_corridor = np.zeros((4, grid.n), dtype=np.float64)
        rain_shadow_index = np.zeros((4, grid.n), dtype=np.float64)
        regional_response = np.ones((4, grid.n), dtype=np.float64)
        lat_abs = np.abs(grid.lat)
        subtropical_dry = np.exp(-((lat_abs - 28.0) / 12.0) ** 2)

        for s in range(4):
            T = seasonal_T[s]
            sst_s = None if seasonal_sst is None else np.asarray(seasonal_sst)[s]
            access = moisture_access[s]
            monsoon = np.clip(monsoon_potential[s], 0.0, 1.2)
            access_eff = np.clip(
                0.58 * access + 0.42 * np.sqrt(np.maximum(access, 0.0))
                + 0.10 * monsoon,
                0.0,
                1.25,
            )
            tfac = np.clip((T - 252.0) / 42.0, 0.05, 1.08)
            evap_s = self._seasonal_evaporation(
                T, ocean, current_heat, upwelling, sst_s, ocean_heat_flux)
            seasonal_evap[s] = evap_s

            convergence = np.clip(
                0.74 * itcz_intensity[s] + 0.70 * storm_tracks[s]
                + 1.25 * monsoon,
                0.0,
                2.6,
            )
            dry_factor = np.clip(
                1.0 - 0.28 * subtropical_dry
                * (1.0 - 0.42 * itcz_intensity[s])
                * (1.0 - 0.60 * monsoon),
                0.58,
                1.05,
            )
            moisture_shadow = np.clip((0.74 - access_eff) / 0.74, 0.0, 1.0)
            continental_subsidence = (
                subtropical_dry
                * (0.35 + 0.65 * np.clip(interiority, 0.0, 1.0))
                * (0.45 + 0.55 * moisture_shadow)
                * (1.0 - 0.55 * itcz_intensity[s])
                * (1.0 - 0.70 * monsoon)
            )
            dry_factor = np.clip(dry_factor - 2.20 * continental_subsidence, 0.08, 1.05)

            land_precip = (
                (280.0 + 1180.0 * convergence) * access_eff
                + 210.0 * storm_tracks[s] * np.sqrt(np.maximum(access_eff, 0.0))
            ) * tfac * dry_factor
            wet_core = np.clip(convergence - 0.31, 0.0, None)
            wet_core_bonus = (
                780.0 * access_eff ** 1.22 * wet_core ** 1.30
                * tfac * (0.68 + 0.32 * itcz_intensity[s])
            )
            deep_convective_core = np.clip(convergence - 0.54, 0.0, None)
            deep_core_bonus = (
                520.0 * access_eff ** 1.55 * deep_convective_core ** 1.22
                * tfac * np.clip((T - 283.0) / 18.0, 0.0, 1.0)
            )
            warm_convective_tail = np.clip(convergence - 0.24, 0.0, None)
            warm_tail_bonus = (
                680.0 * access_eff ** 1.42 * warm_convective_tail ** 1.16
                * tfac * np.clip((T - 281.0) / 19.0, 0.0, 1.0)
            )
            tropical_rain_core = np.exp(-(lat_abs / 20.0) ** 2)
            tropical_tail_bonus = (
                1325.0 * tropical_rain_core
                * access_eff ** 1.36
                * np.clip(convergence - 0.18, 0.0, None) ** 1.18
                * tfac * np.clip((T - 282.0) / 18.0, 0.0, 1.0)
            )
            monsoon_bonus = (
                460.0 * access_eff * monsoon * (0.55 + 0.45 * convergence)
                * tfac
            )
            land_precip = (
                land_precip + wet_core_bonus + deep_core_bonus
                + warm_tail_bonus + tropical_tail_bonus + monsoon_bonus
            )
            land_precip *= np.clip(1.0 - 0.90 * continental_subsidence, 0.16, 1.0)
            ocean_precip = (
                390.0 + 0.34 * evap_s + 210.0 * itcz_intensity[s]
                + 120.0 * storm_tracks[s]
            ) * (1.0 - 0.18 * upwelling)

            oro, lee = self._orographic_precip_terms(
                grid, seasonal_wind[s], elev, land, barrier, wind_gap)
            upslope, downslope = self._slope_wind_exposure(
                grid, seasonal_wind[s], elev, land, barrier, wind_gap)
            oro_eff = np.clip(0.48 * oro + 0.82 * upslope, 0.0, 1.0)
            lee_eff = np.clip(0.48 * lee + 0.82 * downslope, 0.0, 1.0)
            wind_speed = np.linalg.norm(seasonal_wind[s], axis=1)
            wind_unit = np.where(
                wind_speed[:, None] > 1.0e-9,
                seasonal_wind[s] / np.maximum(wind_speed[:, None], 1.0e-9),
                0.0,
            )
            onshore = np.where(
                land,
                np.clip(np.sum(wind_unit * (-coast_orientation), axis=1), 0.0, 1.0)
                * np.clip(coast_strength, 0.0, 1.0),
                0.0,
            )
            onshore = self._smooth_field_masked(
                grid, onshore, land, passes=5, alpha=0.22)
            source_warm = self._smooth_field(
                grid, source_ocean_warmth[s], passes=3, alpha=0.18)
            thermal_low = np.clip(-seasonal_pressure[s], 0.0, 1.6)
            passability = np.clip(1.0 - 0.38 * barrier * (1.0 - wind_gap),
                                  0.35, 1.0)
            monsoon_seed = np.where(
                land,
                np.clip(monsoon, 0.0, 1.2)
                * np.clip(access_eff, 0.0, 1.2)
                * (
                    0.28
                    + 0.38 * onshore
                    + 0.20 * np.clip(source_warm, 0.0, 1.0)
                    + 0.14 * thermal_low
                )
                * passability,
                0.0,
            )
            monsoon_band = self._smooth_field_masked(
                grid, monsoon_seed, land, passes=7, alpha=0.20)
            monsoon_band = np.where(land, np.clip(monsoon_band, 0.0, 1.2), 0.0)

            storm_seed = np.where(
                land,
                storm_tracks[s] * np.clip(access_eff, 0.0, 1.2)
                * (0.45 + 0.35 * onshore + 0.20 * oro_eff)
                * (1.0 - 0.22 * lee_eff),
                0.0,
            )
            storm_band = self._smooth_field_masked(
                grid, storm_seed, land, passes=6, alpha=0.18)
            storm_band = np.where(land, np.clip(storm_band, 0.0, 1.3), 0.0)

            rain_shadow = np.where(
                land,
                lee_eff * (0.38 + 0.62 * barrier)
                * (1.0 - 0.24 * monsoon)
                * np.clip(1.0 - 0.42 * wind_gap, 0.25, 1.0),
                0.0,
            )
            rain_shadow = self._smooth_field_masked(
                grid, rain_shadow, land, passes=5, alpha=0.18)
            rain_shadow = np.where(land, np.clip(rain_shadow, 0.0, 1.2), 0.0)

            monsoon_corridor[s] = monsoon_band
            storm_corridor[s] = storm_band
            rain_shadow_index[s] = rain_shadow
            response = np.where(
                land,
                np.clip(
                    1.0 + 0.18 * monsoon_band + 0.12 * storm_band
                    - 0.20 * rain_shadow,
                    0.62,
                    1.42,
                ),
                1.0,
            )
            regional_response[s] = response
            oro_bonus = 620.0 * access_eff * oro_eff * (0.55 + 0.45 * convergence)
            land_precip = land_precip + oro_bonus
            land_precip *= 1.0 - 0.42 * lee_eff * (1.0 - 0.35 * monsoon)
            precip_s = np.where(land, land_precip, ocean_precip)
            precip_s = self._smooth_field(grid, precip_s, passes=3, alpha=0.15)
            if land.any():
                land_mean_before = float(np.average(precip_s[land], weights=grid.cell_area[land]))
                regional = np.where(land, precip_s * response, precip_s)
                regional_mean = float(np.average(regional[land], weights=grid.cell_area[land]))
                if regional_mean > 1.0e-9:
                    regional[land] *= land_mean_before / regional_mean
                precip_s = regional
                oro_shape = np.clip(
                    1.0 + 0.40 * oro_eff - 0.34 * lee_eff * (1.0 - 0.25 * monsoon),
                    0.58,
                    1.48,
                )
                shaped = np.where(land, precip_s * oro_shape, precip_s)
                land_mean_after = float(np.average(shaped[land], weights=grid.cell_area[land]))
                if land_mean_after > 1.0e-9:
                    shaped[land] *= land_mean_before / land_mean_after
                precip_s = shaped
            seasonal_precip[s] = np.clip(precip_s, 0.0, 4500.0)
            moisture_convergence[s] = np.where(
                land, np.clip(access * convergence, 0.0, 2.0), convergence)
            orographic_component[s] = np.where(
                land,
                np.clip(
                    oro_bonus + np.maximum(precip_s * (oro_eff - lee_eff), 0.0),
                    0.0,
                    1800.0,
                ),
                0.0,
            )

        seasonal_precip = self._accentuate_low_latitude_seasonality(
            grid, seasonal_precip, land)
        precip = seasonal_precip.mean(axis=0)
        evap = seasonal_evap.mean(axis=0)
        evap = np.where(land, np.minimum(0.46 * precip, 950.0), evap)
        annual_mean = np.maximum(precip, 1e-9)
        precip_seasonality = np.max(seasonal_precip, axis=0) / annual_mean

        nh = grid.lat >= 0.0
        summer = np.where(nh, seasonal_precip[2], seasonal_precip[0])
        winter = np.where(nh, seasonal_precip[0], seasonal_precip[2])
        monsoon_index = np.where(land, (summer - winter) / annual_mean, 0.0)
        dry_threshold = np.minimum(350.0, 0.45 * annual_mean)
        dry_season_length = np.sum(seasonal_precip < dry_threshold[None, :], axis=0)
        wet_season_peak = np.argmax(seasonal_precip, axis=0).astype(np.float64)

        return {
            **pressure_moisture,
            "seasonal_precipitation": seasonal_precip,
            "precipitation": np.clip(precip, 0.0, 4500.0),
            "evaporation": np.clip(evap, 0.0, 4200.0),
            "precipitation_seasonality": np.clip(precip_seasonality, 1.0, 4.0),
            "monsoon_index": np.clip(monsoon_index, -3.0, 3.0),
            "dry_season_length": dry_season_length.astype(np.float64),
            "wet_season_peak": wet_season_peak,
            "moisture_convergence": moisture_convergence.mean(axis=0),
            "orographic_precipitation": orographic_component.mean(axis=0),
            "monsoon_rainfall_corridor": monsoon_corridor,
            "storm_track_rainfall_corridor": storm_corridor,
            "rain_shadow_index": rain_shadow_index,
            "regional_precipitation_response": regional_response,
        }

    def _coastal_direction(self, grid, from_mask, to_mask):
        """Tangent unit vector from each from_mask cell toward adjacent to_mask cells."""
        out = np.zeros((grid.n, 3), dtype=np.float64)
        strength = np.zeros(grid.n, dtype=np.float64)
        for c in np.where(from_mask)[0]:
            nbs = grid.neighbors[int(c)]
            targets = nbs[to_mask[nbs]]
            if targets.size == 0:
                continue
            vec = np.mean(grid.xyz[targets], axis=0)
            tangent = vec - float(vec @ grid.xyz[c]) * grid.xyz[c]
            norm = float(np.linalg.norm(tangent))
            if norm <= 1e-12:
                continue
            out[c] = tangent / norm
            strength[c] = min(1.0, targets.size / max(1, nbs.size))
        return out, strength

    def _geographic_circulation_anomalies(self, world, seasonal_T, continentality,
                                          ocean, elev, background_wind,
                                          geography_fields=None):
        grid = world.grid
        land = ~ocean
        if world.spec.orbit.tidally_locked:
            zeros_v = np.zeros_like(background_wind)
            zeros_s = np.zeros((4, grid.n), dtype=np.float64)
            return background_wind, zeros_s, zeros_v, zeros_v, np.zeros(grid.n)

        annual_T = seasonal_T.mean(axis=0, keepdims=True)
        temp_anom = seasonal_T - annual_T
        previous_temp_anom = np.roll(temp_anom, shift=1, axis=0)
        cooling_tendency = previous_temp_anom - temp_anom
        # Do not let a tiny archipelago define a full continental or basin-scale
        # stationary wave just because it is the largest landmass in a
        # waterworld.  Earthlike and arid worlds are effectively unchanged;
        # island worlds keep only weak land-sea pressure anomalies.
        land_extent = np.clip((world.land_fraction() - 0.035) / 0.22, 0.0, 1.0)
        geography_strength = 0.20 + 0.80 * land_extent
        geography_fields = geography_fields or {}
        coast_distance = np.asarray(
            geography_fields.get("climate.coast_distance", np.zeros(grid.n)),
            dtype=np.float64,
        )
        if coast_distance.shape != (grid.n,):
            coast_distance = np.zeros(grid.n, dtype=np.float64)
        continent_interiority = np.asarray(
            geography_fields.get("climate.continent_interiority", continentality),
            dtype=np.float64,
        )
        if continent_interiority.shape != (grid.n,):
            continent_interiority = np.asarray(continentality, dtype=np.float64)
        basin_id = np.asarray(
            geography_fields.get("ocean.basin_id", np.full(grid.n, -1.0)),
            dtype=np.float64,
        )
        if basin_id.shape != (grid.n,):
            basin_id = np.full(grid.n, -1.0, dtype=np.float64)

        land_weight = np.where(land, 0.30 + 0.95 * continentality, 0.08)
        pressure = (-temp_anom + 0.70 * cooling_tendency) * land_weight[None, :]
        for s in range(4):
            pressure[s] = self._smooth_field(grid, pressure[s], passes=18, alpha=0.22)
        scale = max(float(np.percentile(np.abs(pressure), 95)), 2.0)
        pressure = np.clip(pressure / scale, -1.5, 1.5)
        pressure *= 0.18 + 0.82 * land_extent
        basin_scale = np.zeros(grid.n, dtype=np.float64)
        basin_ids = [int(x) for x in np.unique(basin_id[ocean]) if int(x) >= 0]
        if basin_ids:
            basin_areas = {
                bid: float(np.sum(grid.cell_area[ocean & (basin_id.astype(int) == bid)]))
                for bid in basin_ids
            }
            max_basin_area = max(max(basin_areas.values()), 1.0e-12)
            for bid, area in basin_areas.items():
                basin_scale[basin_id.astype(int) == bid] = np.clip(
                    np.sqrt(area / max_basin_area), 0.25, 1.0)
        lat_abs = np.abs(grid.lat)
        open_ocean_pressure = np.where(
            ocean, np.clip(coast_distance, 0.0, 1.0) ** 0.70, 0.0)
        subpolar_ocean = np.exp(-((lat_abs - 52.0) / 18.0) ** 2)
        subtropical_ocean = np.exp(-((lat_abs - 30.0) / 16.0) ** 2)
        season_ratio = np.array([-1.0, 0.0, 1.0, 0.0], dtype=np.float64)
        for s, ratio in enumerate(season_ratio):
            winter_hemi = np.clip(-np.sign(grid.lat) * ratio, 0.0, 1.0)
            summer_hemi = np.clip(np.sign(grid.lat) * ratio, 0.0, 1.0)
            ocean_pressure_source = (
                -0.50
                * open_ocean_pressure
                * basin_scale
                * subpolar_ocean
                * winter_hemi
                + 0.18
                * open_ocean_pressure
                * basin_scale
                * subtropical_ocean
                * (0.45 + 0.55 * summer_hemi)
            )
            pressure[s] += geography_strength * self._smooth_field_masked(
                grid, ocean_pressure_source, ocean, passes=8, alpha=0.18)
        pressure = np.clip(pressure, -1.5, 1.5)

        to_ocean, coast_land_strength = self._coastal_direction(grid, land, ocean)
        to_land, coast_ocean_strength = self._coastal_direction(grid, ocean, land)
        thermal = np.zeros_like(background_wind)
        for s in range(4):
            high_to_low = -self._graph_gradient_vectors(grid, pressure[s])
            mag = np.linalg.norm(high_to_low, axis=1)
            p95 = max(float(np.percentile(mag, 95)), 1e-9)
            high_to_low = high_to_low / p95

            warm_low = np.clip(-pressure[s], 0.0, 1.2)
            cold_high = np.clip(pressure[s], 0.0, 1.2)
            coastal = np.zeros((grid.n, 3), dtype=np.float64)
            # Summer heat lows: onshore flow near warm continental coasts.
            coastal += np.where(
                land[:, None],
                -to_ocean * (coast_land_strength * warm_low)[:, None],
                to_land * (coast_ocean_strength * warm_low)[:, None],
            )
            # Winter cold highs: offshore flow away from cold continental coasts.
            coastal += np.where(
                land[:, None],
                to_ocean * (coast_land_strength * cold_high)[:, None],
                -to_land * (coast_ocean_strength * cold_high)[:, None],
            )

            anomaly = geography_strength * (2.2 * high_to_low + 0.9 * coastal)
            bg_speed = np.linalg.norm(background_wind[s], axis=1)
            cap = 0.9 + 0.36 * bg_speed
            speed = np.linalg.norm(anomaly, axis=1)
            anomaly *= np.minimum(1.0, cap / np.maximum(speed, 1e-9))[:, None]
            thermal[s] = self._smooth_vector_field(grid, anomaly, passes=7, alpha=0.20)

        topo = self._smooth_field(grid, np.maximum(elev, 0.0), passes=4, alpha=0.25)
        topo_scale = max(float(np.percentile(topo[land], 95)) if land.any() else 1.0, 1.0)
        topo_norm = np.clip(topo / topo_scale, 0.0, 2.0)
        topo_grad = self._graph_gradient_vectors(grid, topo_norm)
        topo_mag = np.linalg.norm(topo_grad, axis=1)
        p70 = float(np.percentile(topo_mag[land], 70)) if land.any() else 0.0
        p95 = max(float(np.percentile(topo_mag[land], 95)) if land.any() else 1.0,
                  p70 + 1e-6)
        barrier = np.clip((topo_mag - p70) / (p95 - p70), 0.0, 1.0)
        barrier *= land.astype(float)
        normal = topo_grad / np.maximum(topo_mag[:, None], 1e-9)
        along = self._project_tangent(grid, np.cross(grid.xyz, normal))
        along_norm = np.linalg.norm(along, axis=1, keepdims=True)
        along = np.where(along_norm > 1e-9, along / np.maximum(along_norm, 1e-9), 0.0)

        orographic = np.zeros_like(background_wind)
        for s in range(4):
            pre = background_wind[s] + thermal[s]
            cross = np.sum(pre * normal, axis=1)
            damp_cross = -0.38 * barrier[:, None] * cross[:, None] * normal
            deflect = 0.16 * barrier[:, None] * cross[:, None] * along
            orographic[s] = self._project_tangent(grid, damp_cross + deflect)

        final = background_wind + thermal + orographic
        lat_abs = np.abs(grid.lat)
        east, _ = self._tangent_basis(grid)
        trade_band = np.exp(-((lat_abs - 16.0) / 20.0) ** 2)
        westerly_band = np.exp(-((lat_abs - 50.0) / 15.0) ** 2)
        open_ocean = np.where(ocean, np.clip(coast_distance, 0.0, 1.0) ** 0.70, 0.0)
        interior_land = np.where(land, np.clip(continent_interiority, 0.0, 1.0), 0.0)
        ocean_westerly_gain = (
            np.where(grid.lat < 0.0, 1.08, 0.48)
            * open_ocean
            * westerly_band
        )
        ocean_trade_gain = 0.44 * open_ocean * trade_band
        land_surface_drag = (
            land.astype(float)
            * (0.10 + 0.18 * interior_land)
            * (0.55 * trade_band + 0.85 * westerly_band)
        )
        terrain_roughness = np.where(
            land,
            0.35
            + 0.45 * interior_land
            + 0.25 * np.clip(
                barrier * (1.0 - 0.60 * np.clip(
                    geography_fields.get("terrain.wind_gap_index", np.zeros(grid.n)),
                    0.0,
                    1.0,
                )),
                0.0,
                1.0,
            ),
            0.0,
        )
        stationary_band_response = 0.45 + 0.55 * (0.45 * trade_band + 0.75 * westerly_band)
        open_ocean_storm_tail = (
            open_ocean
            * (0.85 * westerly_band + 0.10 * trade_band)
        )
        coriolis_support = np.sign(grid.lat) * np.clip((lat_abs - 8.0) / 28.0, 0.0, 1.0)
        pressure_steering_band = 0.45 + 0.55 * (
            0.35 * trade_band + 0.75 * westerly_band)
        pressure_steering_support = (
            geography_strength
            * pressure_steering_band
            * (0.40 + 0.45 * open_ocean + 0.15 * interior_land)
        )
        nonpolar_land = (
            land.astype(float)
            * (1.0 - np.clip((lat_abs - 55.0) / 18.0, 0.0, 1.0))
        )
        nonpolar_boundary_drag = (
            nonpolar_land
            * (
                (0.45 + 0.55 * interior_land) * stationary_band_response
                + 0.50 * terrain_roughness
            )
        )
        polar_katabatic_support = (
            land.astype(float)
            * np.clip((lat_abs - 58.0) / 22.0, 0.0, 1.0)
            * np.clip((topo_norm - 0.20) / 0.85, 0.0, 1.0)
            * np.clip(0.35 + 0.65 * interior_land, 0.0, 1.0)
        )
        pressure_stationary = np.zeros((4, grid.n), dtype=np.float64)
        broad_pressure_stationary = np.zeros((4, grid.n), dtype=np.float64)
        pressure_steering = np.zeros_like(final)
        katabatic_wind = np.zeros_like(final)
        for s in range(4):
            p_wave = np.abs(self._latitude_band_anomaly(grid, pressure[s]))
            p10, p90 = np.percentile(p_wave[np.isfinite(p_wave)], [10, 90])
            pressure_stationary[s] = np.clip(
                (p_wave - p10) / max(float(p90 - p10), 1.0e-9),
                0.0,
                1.0,
            )
            broad_wave = (
                self._smooth_field_masked(
                    grid, p_wave, land, passes=10, alpha=0.18)
                + self._smooth_field_masked(
                    grid, p_wave, ocean, passes=10, alpha=0.18)
            )
            b10, b90 = np.percentile(broad_wave[np.isfinite(broad_wave)], [10, 90])
            broad_pressure_stationary[s] = np.clip(
                (broad_wave - b10) / max(float(b90 - b10), 1.0e-9),
                0.0,
                1.0,
            )
            steering_pressure = self._smooth_field(
                grid,
                self._latitude_band_anomaly(grid, pressure[s]),
                passes=16,
                alpha=0.16,
            )
            pressure_grad = self._graph_gradient_vectors(grid, steering_pressure)
            pressure_grad_mag = np.linalg.norm(pressure_grad, axis=1)
            pressure_grad_scale = max(float(np.percentile(pressure_grad_mag, 95)), 1.0e-9)
            pressure_grad = pressure_grad / pressure_grad_scale
            geostrophic = (
                np.cross(grid.xyz, pressure_grad)
                * coriolis_support[:, None]
                * pressure_steering_support[:, None]
            )
            pressure_steering[s] = 1.5 * self._smooth_vector_field(
                grid, geostrophic, passes=4, alpha=0.10)
        stationary_wave = np.zeros_like(final)
        for s in range(4):
            seasonal_stationary = (
                (ocean_westerly_gain - ocean_trade_gain)[:, None] * east
                - land_surface_drag[:, None] * final[s]
            )
            seasonal_stationary = self._smooth_vector_field(
                grid, seasonal_stationary, passes=2, alpha=0.12)
            stationary_speed = np.linalg.norm(seasonal_stationary, axis=1)
            seasonal_stationary *= np.minimum(
                1.0, 1.35 / np.maximum(stationary_speed, 1.0e-9))[:, None]
            stationary_wave[s] = seasonal_stationary
            final[s] = self._smooth_vector_field(
                grid, final[s] + seasonal_stationary, passes=1, alpha=0.06)
            roughness_support = 0.65 * open_ocean + land.astype(float)
            speed_multiplier = (
                1.0
                + geography_strength
                * 0.18
                * (pressure_stationary[s] - 0.45)
                * roughness_support
                * stationary_band_response
                - geography_strength
                * 0.16
                * terrain_roughness
                * stationary_band_response
                + geography_strength
                * 0.06
                * open_ocean
                * stationary_band_response
                + geography_strength
                * 0.06
                * (broad_pressure_stationary[s] - 0.56)
                * (0.90 * land.astype(float) + 0.35 * open_ocean)
                * stationary_band_response
                - geography_strength
                * 0.16
                * terrain_roughness
                * stationary_band_response
                + geography_strength
                * 0.40
                * open_ocean_storm_tail
                * (0.70 + 0.30 * broad_pressure_stationary[s])
                - geography_strength
                * 0.33
                * land.astype(float)
                * (0.35 + 0.65 * interior_land)
                * (0.55 + 0.45 * pressure_stationary[s])
                * stationary_band_response
            )
            final[s] *= np.clip(speed_multiplier, 0.70, 1.38)[:, None]
            final[s] *= np.clip(
                1.0 - geography_strength * 0.24 * nonpolar_boundary_drag,
                0.70,
                1.0,
            )[:, None]
            final[s] = final[s] + pressure_steering[s]
            winter_hemi = np.clip(-np.sign(grid.lat) * season_ratio[s], 0.0, 1.0)
            katabatic_wind[s] = (
                8.0
                * polar_katabatic_support[:, None]
                * (0.85 + 0.25 * winter_hemi)[:, None]
                * (-normal)
            )
            final[s] = final[s] + katabatic_wind[s]
            speed = np.linalg.norm(final[s], axis=1)
            final[s] *= np.minimum(1.0, 24.0 / np.maximum(speed, 1e-9))[:, None]
            final[s] = self._project_tangent(grid, final[s])

        geo_speed = (
            np.linalg.norm(thermal, axis=2).mean(axis=0)
            + np.linalg.norm(orographic, axis=2).mean(axis=0)
            + np.linalg.norm(stationary_wave, axis=2).mean(axis=0)
            + np.linalg.norm(pressure_steering, axis=2).mean(axis=0)
            + np.linalg.norm(katabatic_wind, axis=2).mean(axis=0)
        )
        geo_index = geo_speed / 6.0
        return final, pressure, thermal, orographic, np.clip(geo_index, 0.0, 2.0)

    def _polar_continental_ice_cap_cooling(
        self,
        grid,
        land,
        lapse_elev,
        geography_fields,
    ):
        """Extra cooling for large polar continental ice-cap settings.

        The reduced EBM otherwise treats Antarctica too much like ordinary
        coastal highland because its cell-scale coastline keeps the generic
        continentality proxy low.  This term is component-based, so small polar
        islands and non-polar mountain belts are not pulled into the ice-cap
        regime.
        """
        land = np.asarray(land, dtype=bool)
        lapse_elev = np.asarray(lapse_elev, dtype=np.float64)
        continent_id = np.asarray(
            geography_fields.get("climate.continent_id", np.full(grid.n, -1.0)),
            dtype=np.int64,
        )
        if continent_id.shape != (grid.n,) or not land.any():
            return np.zeros(grid.n, dtype=np.float64)

        total_area = max(float(np.sum(grid.cell_area)), 1.0e-12)
        cooling = np.zeros(grid.n, dtype=np.float64)
        for cid in [int(x) for x in np.unique(continent_id[land]) if int(x) >= 0]:
            cells = land & (continent_id == cid)
            if int(np.count_nonzero(cells)) < 4:
                continue
            area_fraction = float(np.sum(grid.cell_area[cells]) / total_area)
            centroid_lat = float(np.average(
                grid.lat[cells], weights=grid.cell_area[cells]))
            mean_elev = float(np.average(
                lapse_elev[cells], weights=grid.cell_area[cells]))
            polar_score = np.clip((abs(centroid_lat) - 58.0) / 16.0, 0.0, 1.0)
            area_score = np.clip((area_fraction - 0.010) / 0.020, 0.0, 1.0)
            height_score = np.clip((mean_elev - 900.0) / 1300.0, 0.0, 1.0)
            cap_score = polar_score * area_score * height_score
            if cap_score <= 0.0:
                continue
            lat_score = np.clip((np.abs(grid.lat) - 55.0) / 25.0, 0.0, 1.0)
            local_height = np.clip((lapse_elev - 900.0) / 2600.0, 0.0, 1.0)
            component_cooling = cap_score * (
                5.0 + 15.0 * lat_score + 7.0 * local_height)
            cooling[cells] = np.maximum(cooling[cells], component_cooling[cells])

        cooling = self._smooth_field_masked(grid, cooling, land, passes=2, alpha=0.12)
        return np.where(land, np.clip(cooling, 0.0, 24.0), 0.0)

    def _energy_boundary_diagnostics(
        self,
        world,
        annual_flux,
        seasonal_T,
        seasonal_sst,
        ocean,
        continentality,
        lapse_cooling,
        c5_feedback,
    ):
        """Archive M1 energy-boundary support fields without changing physics."""
        grid = world.grid
        land = ~np.asarray(ocean, dtype=bool)
        ocean = np.asarray(ocean, dtype=bool)
        seasonal_T = np.asarray(seasonal_T, dtype=np.float64)
        expected = (4, grid.n)
        if seasonal_T.shape != expected:
            seasonal_T = np.zeros(expected, dtype=np.float64)
        seasonal_sst = np.asarray(seasonal_sst, dtype=np.float64)
        if seasonal_sst.shape != expected:
            seasonal_sst = np.where(ocean[None, :], seasonal_T, 0.0)
        seasonal_flux = self._seasonal_surface_flux(world, annual_flux)
        if seasonal_flux.shape != expected:
            seasonal_flux = np.zeros(expected, dtype=np.float64)
        insolation_anomaly = seasonal_flux - seasonal_flux.mean(axis=0, keepdims=True)

        continentality = np.asarray(continentality, dtype=np.float64)
        if continentality.shape != (grid.n,):
            continentality = np.zeros(grid.n, dtype=np.float64)
        heat_capacity = np.where(
            ocean,
            1.0,
            np.clip(0.58 - 0.38 * np.clip(continentality, 0.0, 1.0), 0.18, 0.58),
        )

        thermal_anomaly = seasonal_T - seasonal_T.mean(axis=0, keepdims=True)
        land_thermal_anomaly = np.where(land[None, :], thermal_anomaly, 0.0)
        sst_anomaly = seasonal_sst - seasonal_sst.mean(axis=0, keepdims=True)
        ocean_mixed_layer_anomaly = np.where(ocean[None, :], sst_anomaly, 0.0)

        seasonal_sea_ice = np.asarray(
            c5_feedback.get("seasonal_sea_ice", np.zeros(expected)),
            dtype=np.float64,
        )
        if seasonal_sea_ice.shape != expected:
            seasonal_sea_ice = np.zeros(expected, dtype=np.float64)
        seasonal_snow = np.asarray(
            c5_feedback.get("seasonal_snow", np.zeros(expected)),
            dtype=np.float64,
        )
        if seasonal_snow.shape != expected:
            seasonal_snow = np.zeros(expected, dtype=np.float64)
        snow_ice_support = np.clip(
            np.where(ocean[None, :], seasonal_sea_ice, seasonal_snow),
            0.0,
            1.0,
        )

        same_latitude_sst_anomaly = np.zeros(expected, dtype=np.float64)
        sst_gradient_support = np.zeros(expected, dtype=np.float64)
        band_id = np.floor((grid.lat + 90.0) / 5.0).astype(int)
        for season in range(4):
            for band in np.unique(band_id):
                cells = (
                    ocean
                    & (band_id == int(band))
                    & np.isfinite(seasonal_sst[season])
                )
                if not cells.any():
                    continue
                mean = float(np.average(
                    seasonal_sst[season, cells],
                    weights=grid.cell_area[cells],
                ))
                same_latitude_sst_anomaly[season, cells] = (
                    seasonal_sst[season, cells] - mean)

            grad = np.linalg.norm(
                self._graph_gradient_vectors(grid, seasonal_sst[season]),
                axis=1,
            )
            scale = max(float(np.percentile(grad[ocean], 95)) if ocean.any() else 1.0,
                        1.0e-9)
            sst_gradient_support[season] = np.where(
                ocean,
                np.clip(grad / scale, 0.0, 1.5),
                0.0,
            )

        land_sea_contrast = np.zeros(expected, dtype=np.float64)
        for season in range(4):
            land_signal = np.where(land, land_thermal_anomaly[season], 0.0)
            ocean_signal = np.where(ocean, ocean_mixed_layer_anomaly[season], 0.0)
            for _ in range(10):
                land_signal = np.where(
                    land,
                    land_thermal_anomaly[season],
                    0.76 * self._neighbor_mean(grid, land_signal),
                )
                ocean_signal = np.where(
                    ocean,
                    ocean_mixed_layer_anomaly[season],
                    0.76 * self._neighbor_mean(grid, ocean_signal),
                )
            land_sea_contrast[season] = np.clip(land_signal - ocean_signal, -45.0, 45.0)

        lapse_cooling = np.asarray(lapse_cooling, dtype=np.float64)
        if lapse_cooling.shape != (grid.n,):
            lapse_cooling = np.zeros(grid.n, dtype=np.float64)

        return {
            "climate.seasonal_insolation_anomaly": insolation_anomaly,
            "climate.surface_heat_capacity_class": heat_capacity,
            "climate.land_thermal_anomaly": land_thermal_anomaly,
            "climate.ocean_mixed_layer_thermal_anomaly": ocean_mixed_layer_anomaly,
            "climate.elevation_lapse_cooling": np.clip(lapse_cooling, 0.0, 80.0),
            "climate.snow_ice_albedo_support": snow_ice_support,
            "climate.sst_gradient_support": sst_gradient_support,
            "climate.same_latitude_sst_anomaly": same_latitude_sst_anomaly,
            "climate.land_sea_thermal_contrast": land_sea_contrast,
        }

    def step(self, world, t, dt, forcing, rng_key) -> StepResult:
        grid = world.grid
        self._edge_cache(grid)
        flux = self._solar_factor(world)
        elev_raw = world.get_field("terrain.elevation_m", 0.0)
        ocean = world.ocean_mask()
        rel_elev = elev_raw - world.sea_level
        elev = self._climate_elevation(grid, rel_elev)
        geography_fields, geography_objects, geography_diag = self._geography_primitives(
            world, ocean, rel_elev, elev)
        co2 = world.g("atmosphere.co2", self.CO2_REF)
        pressure = world.g("atmosphere.pressure", 1.0)

        # greenhouse forcing relative to reference
        f_gh = 5.35 * np.log(max(co2, 1e-9) / self.CO2_REF) + 8.0 * np.log(max(pressure, 1e-3))

        # heat capacity proxy: ocean mixes more -> stronger transport
        T = world.get_field("climate.surface_temperature", 288.0).copy()
        if T.mean() < 100 or T.mean() > 400:
            T = np.full(grid.n, 288.0)
        lapse = 7.2e-3                       # K/m
        elev_pos = np.maximum(elev, 0.0)
        raw_elev_pos = np.maximum(rel_elev, 0.0)
        lapse_elev = np.where(
            rel_elev > 0.0,
            0.65 * elev_pos + 0.35 * raw_elev_pos,
            elev_pos,
        )
        land_permanent_ice = np.clip((np.abs(grid.lat) - 62.0) / 20.0, 0.0, 1.0)
        land_permanent_ice *= np.clip((lapse_elev - 900.0) / 1800.0, 0.0, 1.0)
        land_ice_albedo = 0.50 + 0.08 * land_permanent_ice

        for _ in range(self.EBM_ITERS):
            ice = self._ice_fraction(T, ocean)
            # Include unresolved cloud/atmospheric reflection in the surface
            # classes.  A bare-ocean albedo of 0.10 made the EBM absorb far too
            # much shortwave energy and generated 90 C tropical oceans.
            base_albedo = np.where(ocean, 0.285, 0.26)
            ice_albedo = np.where(ocean, 0.52, land_ice_albedo)
            # Keep ice-albedo feedback continuous.  A hard jump at the first
            # trace of ice produced an unrealistic thermal wall at the ice edge.
            albedo = base_albedo * (1.0 - ice) + ice_albedo * ice
            absorbed = flux * (1.0 - albedo)
            olr = self.A0 - f_gh + self.B * (T - CONSTANTS.ZERO_C)
            net = absorbed - olr
            T = T + 0.25 * net / self.B
            # meridional/lateral heat transport
            transport = np.where(ocean, 0.62, 0.52)
            T = (1.0 - transport) * T + transport * self._neighbor_mean(grid, T)

        T_actual = T - lapse * lapse_elev        # elevation cooling for surface
        for _ in range(14):
            T_actual = 0.88 * T_actual + 0.12 * self._neighbor_mean(grid, T_actual)
        T_actual = self._smooth_zonal_background(grid, T_actual)
        # The reduced EBM is a climate-state proxy, not a local weather model.
        # Bound pathological numerical excursions while preserving snowball and
        # hothouse regimes for non-Earth presets.
        T_actual = np.clip(T_actual, 185.0, 330.0)
        seasonal_T, temp_seasonality, continentality = self._seasonal_temperature(
            world, T_actual, flux, ocean)
        T_actual = seasonal_T.mean(axis=0)

        sea_ice = self._ice_fraction(T_actual, ocean) * ocean
        # ice sheets: cold land with some moisture
        background_wind, itcz_lat, itcz_intensity, storm_tracks = self._seasonal_circulation(world)
        (seasonal_wind, pressure_proxy, thermal_anom, orographic_anom,
         geographic_index) = self._geographic_circulation_anomalies(
            world, seasonal_T, continentality, ocean, elev, background_wind,
            geography_fields)
        ocean_current = self._ocean_currents(
            world, ocean, seasonal_wind, geography_fields)
        coupling = self._weak_ocean_atmosphere_coupling(
            world, seasonal_T, ocean, seasonal_wind, pressure_proxy,
            ocean_current, geography_fields)
        seasonal_wind = coupling["seasonal_wind"]
        pressure_proxy = coupling["pressure_proxy"]
        ocean_current = coupling["ocean_current"]
        pressure_genesis = self._m2_pressure_genesis(
            world,
            pressure_proxy,
            seasonal_T,
            coupling["seasonal_sst"],
            ocean,
            elev,
            geography_fields,
        )
        pressure_proxy = pressure_genesis["pressure"]
        currents = ocean_current["currents"]
        current_heat = ocean_current["current_heat_transport"]
        upwelling = ocean_current["upwelling"]
        coupled_sst = np.asarray(coupling["seasonal_sst"], dtype=np.float64)
        ocean_heat_flux = np.asarray(coupling["ocean_heat_flux"], dtype=np.float64)
        raw_ocean_T = np.where(
            ocean[None, :],
            seasonal_T + ocean_heat_flux[None, :],
            seasonal_T + current_heat[None, :],
        )
        raw_T_actual = raw_ocean_T.mean(axis=0)
        sea_ice = self._ice_fraction(raw_T_actual, ocean) * ocean
        ocean_freezing_K = CONSTANTS.ZERO_C - 1.8
        # Ocean cells represent a mixed surface layer/SST for calibration
        # against OISST.  Sea-ice state still comes from the raw cold solution,
        # but exported ocean temperature should not fall far below seawater's
        # freezing point as if it were bare land or near-surface air.
        seasonal_T = np.where(
            ocean[None, :],
            np.maximum(coupled_sst, ocean_freezing_K),
            np.clip(seasonal_T + current_heat[None, :], 120.0, 360.0),
        )
        polar_land = (~ocean) & (np.abs(grid.lat) >= 50.0)
        polar_land_cooling = np.where(
            polar_land,
            5.2 * np.clip((np.abs(grid.lat) - 50.0) / 32.0, 0.0, 1.0)
            * (0.55 + 0.45 * np.clip(continentality, 0.0, 1.0)),
            0.0,
        )
        polar_land_cooling += self._polar_continental_ice_cap_cooling(
            grid, ~ocean, lapse_elev, geography_fields)
        seasonal_T = np.where(
            polar_land[None, :],
            seasonal_T - polar_land_cooling[None, :],
            seasonal_T,
        )
        temp_seasonality = seasonal_T.max(axis=0) - seasonal_T.min(axis=0)
        T_actual = seasonal_T.mean(axis=0)
        land_sea_pressure_proxy = pressure_proxy.copy()
        pressure_centers = self._pressure_center_diagnostics(
            grid, land_sea_pressure_proxy, ~ocean, ocean, geography_fields)
        hydro, hydro_feedback = self._seasonal_hydroclimate_feedback_loop(
            world, seasonal_T, ocean, seasonal_wind, elev, current_heat,
            upwelling, itcz_intensity, storm_tracks, pressure_proxy,
            geography_fields, coupling["seasonal_sst"],
            coupling["ocean_heat_flux"])
        seasonal_wind = hydro_feedback["seasonal_wind"]
        pressure_proxy = hydro_feedback["pressure_proxy"]
        if world.spec.orbit.tidally_locked:
            locked_wind = self._tidally_locked_wind(world)
            seasonal_wind = np.repeat(locked_wind[None, :, :], 4, axis=0)
        wind = seasonal_wind.mean(axis=0)
        moisture_flow = self._seasonal_moisture_flow_networks(
            grid, hydro, ~ocean, ocean, seasonal_wind, geography_fields,
            current_heat, upwelling)
        hydro, c4f_diag = self._apply_moisture_flow_precipitation_response(
            grid, hydro, moisture_flow, ~ocean, ocean, geography_fields)
        initial_precipitation_response = self._precipitation_response_region_objects(
            grid, hydro, moisture_flow, ~ocean)
        initial_receiver_catchments = self._receiver_catchment_objects(
            grid, hydro, moisture_flow, initial_precipitation_response, ~ocean)
        initial_source_receiver_accounting = self._source_basin_receiver_accounting(
            grid, hydro, moisture_flow, initial_receiver_catchments, ~ocean, ocean)
        hydro, receiver_feedback_diag = self._apply_receiver_supply_precipitation_feedback(
            grid, hydro, initial_source_receiver_accounting, ~ocean)
        hydroclimate_objects = self._hydroclimate_region_objects(
            grid, hydro, ~ocean, geography_fields)
        moisture_flow_objects = self._refresh_moisture_flow_object_precipitation(
            grid,
            moisture_flow["objects"],
            moisture_flow["network_id"],
            hydro["seasonal_precipitation"],
        )
        precipitation_response = self._precipitation_response_region_objects(
            grid, hydro, moisture_flow, ~ocean)
        receiver_catchments = self._receiver_catchment_objects(
            grid, hydro, moisture_flow, precipitation_response, ~ocean)
        source_receiver_accounting = self._source_basin_receiver_accounting(
            grid, hydro, moisture_flow, receiver_catchments, ~ocean, ocean)
        receiver_catchments = {
            **receiver_catchments,
            "objects": source_receiver_accounting["objects"],
        }
        evap = hydro["evaporation"]
        precip = hydro["precipitation"]
        runoff = np.where(~ocean, np.maximum(precip - evap, 0.0), 0.0)
        c5_feedback = self._cryosphere_cloud_vegetation_feedbacks(
            world,
            raw_ocean_T,
            seasonal_T,
            ocean,
            hydro,
            storm_tracks,
        )
        sea_ice = c5_feedback["sea_ice"]
        energy_boundary = self._energy_boundary_diagnostics(
            world,
            flux,
            seasonal_T,
            coupling["seasonal_sst"],
            ocean,
            continentality,
            lapse * lapse_elev,
            c5_feedback,
        )

        ice_sheet = np.where((~ocean) & (T_actual < 263.0),
                             np.clip((263.0 - T_actual) * 12.0, 0.0, 3000.0), 0.0)

        world.provenance.record(Provenance(
            "climate.surface_temperature", self.name, self.fidelity, "K",
            uncertainty=(float(T_actual.min()), float(T_actual.max())),
            direct_cause=f"EBM: CO2={co2*1e6:.0f} ppm-eq, greenhouse {f_gh:.1f} W/m^2, "
            f"ice-albedo feedback"))
        world.provenance.record(Provenance(
            "climate.precipitation", self.name, self.fidelity, "mm/yr",
            direct_cause="ocean evaporation advected by winds + orographic uplift / rain shadow"))

        mean_T = float(np.average(T_actual, weights=grid.cell_area))
        fields = {
            **geography_fields,
            **energy_boundary,
            "climate.surface_temperature": T_actual,
            "climate.precipitation": precip,
            "climate.evaporation": evap,
            "climate.runoff": runoff,
            "climate.seasonal_sst": coupling["seasonal_sst"],
            "climate.ocean_heat_flux": coupling["ocean_heat_flux"],
            "climate.ocean_evaporation_feedback": (
                coupling["evaporation_feedback"]),
            "climate.coupling_residual": coupling["coupling_residual"],
            "climate.seasonal_precipitation": hydro["seasonal_precipitation"],
            "climate.precipitation_seasonality": hydro["precipitation_seasonality"],
            "climate.monsoon_index": hydro["monsoon_index"],
            "climate.dry_season_length": hydro["dry_season_length"],
            "climate.wet_season_peak": hydro["wet_season_peak"],
            "climate.moisture_convergence": hydro["moisture_convergence"],
            "climate.orographic_precipitation": hydro["orographic_precipitation"],
            "climate.monsoon_rainfall_corridor": hydro["monsoon_rainfall_corridor"],
            "climate.storm_track_rainfall_corridor": (
                hydro["storm_track_rainfall_corridor"]),
            "climate.rain_shadow_index": hydro["rain_shadow_index"],
            "climate.regional_precipitation_response": (
                hydro["regional_precipitation_response"]),
            "atmosphere.moisture_flow_source": moisture_flow["source"],
            "atmosphere.moisture_flow_pathway": moisture_flow["pathway"],
            "atmosphere.moisture_source_basin_id": moisture_flow["source_basin_id"],
            "climate.moisture_flow_network_id": moisture_flow["network_id"],
            "climate.moisture_flow_precipitation_response": (
                hydro["moisture_flow_precipitation_response"]),
            "climate.moisture_budget_region_id": hydro["moisture_budget_region_id"],
            "climate.precipitation_response_region_id": (
                precipitation_response["region_id"]),
            "climate.receiver_catchment_id": receiver_catchments["catchment_id"],
            "climate.source_basin_supply_index": (
                source_receiver_accounting["source_basin_supply_index"]),
            "climate.receiver_catchment_supply_balance": (
                source_receiver_accounting["receiver_catchment_supply_balance"]),
            "climate.receiver_supply_precipitation_feedback": (
                hydro["receiver_supply_precipitation_feedback"]),
            "climate.hydro_coupling_residual": hydro_feedback["residual"],
            "climate.hydro_feedback_iteration_delta": (
                hydro_feedback["iteration_delta"]),
            "atmosphere.wind": wind,
            "ocean.currents": currents,
            "ocean.current_heat_transport": current_heat,
            "ocean.upwelling": upwelling,
            "ocean.gyre_id": ocean_current["gyre_id"],
            "ocean.current_streamfunction": ocean_current["current_streamfunction"],
            "ocean.boundary_current_type": ocean_current["boundary_current_type"],
            "ocean.strait_exchange": ocean_current["strait_exchange"],
            "ocean.wind_stress_current_response": (
                ocean_current["wind_stress_current_response"]),
            "ocean.sst_anomaly": ocean_current["sst_anomaly"],
            "ocean.solved_mask": ocean.astype(float),
            "cryosphere.sea_ice": sea_ice,
            "cryosphere.seasonal_sea_ice": c5_feedback["seasonal_sea_ice"],
            "cryosphere.seasonal_snow": c5_feedback["seasonal_snow"],
            "cryosphere.snow_persistence": c5_feedback["snow_persistence"],
            "cryosphere.ice_sheet": ice_sheet,
            "climate.seasonal_cloud_albedo_proxy": (
                c5_feedback["seasonal_cloud_albedo_proxy"]),
            "climate.cloud_albedo_proxy": c5_feedback["cloud_albedo_proxy"],
            "biosphere.vegetation_climate_feedback": (
                c5_feedback["vegetation_climate_feedback"]),
            "climate.seasonal_temperature": seasonal_T,
            "climate.temperature_seasonality": temp_seasonality,
            "climate.continentality": continentality,
            "atmosphere.seasonal_wind": seasonal_wind,
            "atmosphere.itcz_latitude": itcz_lat,
            "atmosphere.itcz_intensity": itcz_intensity,
            "atmosphere.storm_track_intensity": storm_tracks,
            "atmosphere.background_seasonal_wind": background_wind,
            "atmosphere.land_sea_pressure_proxy": land_sea_pressure_proxy,
            "atmosphere.seasonal_pressure_proxy": hydro["seasonal_pressure"],
            "atmosphere.pressure_center_support": (
                pressure_centers["center_support"]),
            "atmosphere.pressure_center_id": pressure_centers["center_id"],
            "atmosphere.stationary_wave_pressure_support": (
                pressure_centers["stationary_wave_support"]),
            "atmosphere.pressure_genesis_source": pressure_genesis["source"],
            "atmosphere.pressure_genesis_wave_transfer": (
                pressure_genesis["wave_transfer"]),
            "atmosphere.ocean_pressure_low_source_support": (
                pressure_genesis["ocean_low_support"]),
            "atmosphere.ocean_pressure_high_source_support": (
                pressure_genesis["ocean_high_support"]),
            "atmosphere.land_pressure_source_support": (
                pressure_genesis["land_support"]),
            "atmosphere.terrain_pressure_wave_source_support": (
                pressure_genesis["terrain_wave_support"]),
            "atmosphere.precipitation_pressure_feedback": (
                hydro_feedback["pressure_feedback"]),
            "atmosphere.hydro_coupled_wind_anomaly": (
                hydro_feedback["wind_anomaly"]),
            "atmosphere.thermal_wind_anomaly": thermal_anom,
            "atmosphere.orographic_wind_anomaly": orographic_anom,
            "atmosphere.geographic_circulation_index": geographic_index,
            "atmosphere.moisture_access": hydro["moisture_access"],
            "atmosphere.monsoon_potential": hydro["monsoon_potential"],
            "atmosphere.source_ocean_warmth": hydro["source_ocean_warmth"],
            "atmosphere.terrain_blocking": hydro["terrain_blocking"],
        }
        diag = {"mean_T_C": round(mean_T - 273.15, 2),
                "ice_cover": float(np.average((sea_ice + (ice_sheet > 0)) > 0.3,
                                              weights=grid.cell_area)),
                "mean_precip": round(float(np.average(precip, weights=grid.cell_area)), 1),
                "land_precip_seasonality_p75": round(
                    float(np.percentile(hydro["precipitation_seasonality"][~ocean], 75))
                    if (~ocean).any() else 0.0, 3),
                "land_monsoon_index_p90": round(
                    float(np.percentile(hydro["monsoon_index"][~ocean], 90))
                    if (~ocean).any() else 0.0, 3),
                "moisture_access_land_p75": round(
                    float(np.percentile(hydro["moisture_access"][:, ~ocean], 75))
                    if (~ocean).any() else 0.0, 3),
                "itcz_DJF_lat": round(float(itcz_lat[0]), 2),
                "itcz_JJA_lat": round(float(itcz_lat[2]), 2),
                "geographic_circulation_p90": round(
                    float(np.percentile(geographic_index, 90)), 3),
                "continent_count": geography_diag["continent_count"],
                "ocean_basin_count": geography_diag["ocean_basin_count"],
                "coastal_land_cells": geography_diag["coastal_land_cells"],
                "barrier_cells": geography_diag["barrier_cells"],
                "strait_cells": geography_diag["strait_cells"],
                "ocean_heat_transport_p95_C": round(
                    float(np.percentile(np.abs(current_heat), 95)), 2),
                "ocean_evaporation_feedback_p95_C": round(
                    float(np.percentile(
                        np.abs(coupling["evaporation_feedback"][ocean]), 95))
                    if ocean.any() else 0.0,
                    3,
                ),
                "upwelling_p95": round(float(np.percentile(upwelling, 95)), 3),
                "land_temp_seasonality_p50_C": round(
                    float(np.percentile(temp_seasonality[~ocean], 50)) if (~ocean).any()
                    else 0.0, 2),
                "ocean_temp_seasonality_p50_C": round(
                    float(np.percentile(temp_seasonality[ocean], 50)) if ocean.any()
                    else 0.0, 2)}
        diag["seasonal_sea_ice_ocean_p95"] = round(
            float(np.percentile(c5_feedback["seasonal_sea_ice"][:, ocean], 95))
            if ocean.any() else 0.0,
            3,
        )
        diag["snow_persistence_land_p95"] = round(
            float(np.percentile(c5_feedback["snow_persistence"][~ocean], 95))
            if (~ocean).any() else 0.0,
            3,
        )
        diag["cloud_albedo_proxy_p50"] = round(
            float(np.percentile(c5_feedback["cloud_albedo_proxy"], 50)),
            3,
        )
        diag["vegetation_climate_feedback_land_p50"] = round(
            float(np.percentile(c5_feedback["vegetation_climate_feedback"][~ocean], 50))
            if (~ocean).any() else 0.0,
            3,
        )
        diag["hydro_coupling_pressure_feedback_abs_p95"] = round(
            float(np.percentile(
                np.abs(hydro_feedback["pressure_feedback"][:, ~ocean]), 95))
            if (~ocean).any() else 0.0,
            4,
        )
        diag["hydro_coupling_wind_anomaly_p95_m_s"] = round(
            float(np.percentile(
                np.linalg.norm(hydro_feedback["wind_anomaly"][:, ~ocean, :], axis=2),
                95,
            )) if (~ocean).any() else 0.0,
            4,
        )
        diag["hydro_feedback_iteration_delta_p95"] = round(
            float(np.percentile(hydro_feedback["iteration_delta"][~ocean], 95))
            if (~ocean).any() else 0.0,
            5,
        )
        diag["hydro_feedback_iteration_count"] = int(
            hydro_feedback.get("iteration_count", self.HYDRO_FEEDBACK_ITERS))
        diag["moisture_flow_precipitation_response"] = c4f_diag
        diag["receiver_supply_precipitation_feedback"] = receiver_feedback_diag
        diag["source_basin_receiver_accounting"] = (
            source_receiver_accounting["diagnostics"])
        diag["hydroclimate_region_object_count"] = int(len(hydroclimate_objects))
        hydro_kind_counts: dict[str, int] = {}
        for obj in hydroclimate_objects:
            kind = str(obj.get("kind", "unknown"))
            hydro_kind_counts[kind] = hydro_kind_counts.get(kind, 0) + 1
        diag["hydroclimate_region_kind_counts"] = dict(sorted(hydro_kind_counts.items()))
        diag["moisture_flow_network_object_count"] = int(len(moisture_flow_objects))
        moisture_flow_kind_counts: dict[str, int] = {}
        for obj in moisture_flow_objects:
            kind = str(obj.get("kind", "unknown"))
            moisture_flow_kind_counts[kind] = (
                moisture_flow_kind_counts.get(kind, 0) + 1)
        diag["moisture_flow_network_kind_counts"] = dict(
            sorted(moisture_flow_kind_counts.items()))
        precipitation_response_objects = precipitation_response["objects"]
        diag["precipitation_response_region_object_count"] = int(
            len(precipitation_response_objects))
        precipitation_response_kind_counts: dict[str, int] = {}
        for obj in precipitation_response_objects:
            kind = str(obj.get("kind", "unknown"))
            precipitation_response_kind_counts[kind] = (
                precipitation_response_kind_counts.get(kind, 0) + 1)
        diag["precipitation_response_region_kind_counts"] = dict(
            sorted(precipitation_response_kind_counts.items()))
        receiver_catchment_objects = receiver_catchments["objects"]
        diag["receiver_catchment_object_count"] = int(
            len(receiver_catchment_objects))
        receiver_catchment_kind_counts: dict[str, int] = {}
        for obj in receiver_catchment_objects:
            kind = str(obj.get("kind", "unknown"))
            receiver_catchment_kind_counts[kind] = (
                receiver_catchment_kind_counts.get(kind, 0) + 1)
        diag["receiver_catchment_kind_counts"] = dict(
            sorted(receiver_catchment_kind_counts.items()))
        pressure_center_objects = pressure_centers["objects"]
        diag["pressure_center_object_count"] = int(len(pressure_center_objects))
        pressure_center_kind_counts: dict[str, int] = {}
        pressure_center_season_counts: dict[str, int] = {}
        for obj in pressure_center_objects:
            kind = str(obj.get("kind", "unknown"))
            season = str(obj.get("season", "unknown"))
            pressure_center_kind_counts[kind] = (
                pressure_center_kind_counts.get(kind, 0) + 1)
            pressure_center_season_counts[season] = (
                pressure_center_season_counts.get(season, 0) + 1)
        diag["pressure_center_kind_counts"] = dict(
            sorted(pressure_center_kind_counts.items()))
        diag["pressure_center_season_counts"] = dict(
            sorted(pressure_center_season_counts.items()))
        diag["pressure_center_support_p95"] = round(
            float(np.percentile(pressure_centers["center_support"], 95)), 3)
        diag["stationary_wave_pressure_support_p95"] = round(
            float(np.percentile(
                pressure_centers["stationary_wave_support"], 95)), 3)
        objects = dict(geography_objects)
        objects["atmosphere.pressure_centers"] = pressure_center_objects
        objects["climate.hydroclimate_regions"] = hydroclimate_objects
        objects["climate.moisture_flow_networks"] = moisture_flow_objects
        objects["climate.precipitation_response_regions"] = (
            precipitation_response_objects)
        objects["climate.receiver_catchments"] = receiver_catchment_objects
        return StepResult(
            state_delta={"fields": fields, "objects": objects},
            diagnostics=diag,
        )

    # ------------------------------------------------------------------
    def _ice_fraction(self, T, ocean):
        return np.clip((271.0 - T) / 12.0, 0.0, 1.0)

    def _cryosphere_cloud_vegetation_feedbacks(
        self,
        world,
        raw_ocean_T,
        seasonal_T,
        ocean,
        hydro,
        storm_tracks,
    ):
        grid = world.grid
        land = ~ocean
        expected = (4, grid.n)
        raw_ocean_T = np.asarray(raw_ocean_T, dtype=np.float64)
        seasonal_T = np.asarray(seasonal_T, dtype=np.float64)
        if raw_ocean_T.shape != expected:
            raw_ocean_T = np.repeat(seasonal_T.mean(axis=0, keepdims=True), 4, axis=0)
        if seasonal_T.shape != expected:
            seasonal_T = np.repeat(
                np.asarray(world.get_field("climate.surface_temperature", 288.0),
                           dtype=np.float64)[None, :],
                4,
                axis=0,
            )

        seasonal_precip = np.asarray(
            hydro.get("seasonal_precipitation", np.zeros(expected)),
            dtype=np.float64,
        )
        if seasonal_precip.shape != expected:
            seasonal_precip = np.zeros(expected, dtype=np.float64)
        moisture = np.asarray(
            hydro.get("moisture_access", np.zeros(expected)),
            dtype=np.float64,
        )
        if moisture.shape != expected:
            moisture = np.zeros(expected, dtype=np.float64)
        storms = np.asarray(storm_tracks, dtype=np.float64)
        if storms.shape != expected:
            storms = np.zeros(expected, dtype=np.float64)

        seasonal_sea_ice = self._ice_fraction(raw_ocean_T, ocean[None, :]) * ocean[None, :]
        for s in range(4):
            seasonal_sea_ice[s] = self._smooth_field_masked(
                grid, seasonal_sea_ice[s], ocean, passes=2, alpha=0.16)
        seasonal_sea_ice = np.clip(seasonal_sea_ice, 0.0, 1.0) * ocean[None, :]
        sea_ice = np.clip(
            0.58 * seasonal_sea_ice.mean(axis=0) + 0.42 * seasonal_sea_ice.max(axis=0),
            0.0,
            1.0,
        )
        sea_ice = self._smooth_field_masked(grid, sea_ice, ocean, passes=1, alpha=0.10)
        sea_ice = np.clip(sea_ice, 0.0, 1.0) * ocean

        snow_temp = np.clip((276.0 - seasonal_T) / 18.0, 0.0, 1.0)
        snow_precip = np.clip(seasonal_precip / 520.0, 0.0, 1.0)
        seasonal_snow = land[None, :] * snow_temp * (0.20 + 0.80 * snow_precip)
        for s in range(4):
            seasonal_snow[s] = self._smooth_field_masked(
                grid, seasonal_snow[s], land, passes=2, alpha=0.14)
        seasonal_snow = np.clip(seasonal_snow, 0.0, 1.0) * land[None, :]
        snow_persistence = np.clip(
            0.55 * seasonal_snow.mean(axis=0) + 0.45 * seasonal_snow.max(axis=0),
            0.0,
            1.0,
        )
        snow_persistence = self._smooth_field_masked(
            grid, snow_persistence, land, passes=1, alpha=0.12)
        snow_persistence = np.clip(snow_persistence, 0.0, 1.0) * land

        storm_norm = np.clip(storms / 1.35, 0.0, 1.0)
        precip_cloud = np.clip(seasonal_precip / 1200.0, 0.0, 1.0)
        dry_subtropics = np.exp(-((np.abs(grid.lat) - 27.0) / 11.0) ** 2)
        cold_stable = np.clip((268.0 - seasonal_T) / 28.0, 0.0, 1.0)
        seasonal_cloud = (
            0.08
            + 0.42 * np.clip(moisture, 0.0, 1.0)
            + 0.24 * storm_norm
            + 0.18 * precip_cloud
            + 0.10 * ocean[None, :]
            + 0.06 * cold_stable
            - 0.20 * dry_subtropics[None, :] * (1.0 - np.clip(moisture, 0.0, 1.0))
        )
        for s in range(4):
            seasonal_cloud[s] = self._smooth_field(
                grid, seasonal_cloud[s], passes=2, alpha=0.14)
        seasonal_cloud = np.clip(seasonal_cloud, 0.0, 1.0)
        cloud_albedo = np.clip(seasonal_cloud.mean(axis=0), 0.0, 1.0)

        annual_T = seasonal_T.mean(axis=0)
        annual_precip = np.asarray(
            hydro.get("precipitation", seasonal_precip.mean(axis=0)),
            dtype=np.float64,
        )
        if annual_precip.shape != (grid.n,):
            annual_precip = seasonal_precip.mean(axis=0)
        growing_temp = (
            np.clip((annual_T - 268.0) / 24.0, 0.0, 1.0)
            * np.clip((315.0 - annual_T) / 22.0, 0.0, 1.0)
        )
        moisture_support = np.clip(moisture.mean(axis=0), 0.0, 1.0)
        precip_support = np.clip(annual_precip / 1350.0, 0.0, 1.0)
        dry_length = np.asarray(
            hydro.get("dry_season_length", np.zeros(grid.n)), dtype=np.float64)
        if dry_length.shape != (grid.n,):
            dry_length = np.zeros(grid.n, dtype=np.float64)
        dry_penalty = np.clip((dry_length - 1.5) / 2.5, 0.0, 1.0)
        vegetation_feedback = (
            0.48 * growing_temp
            + 0.34 * moisture_support
            + 0.24 * precip_support
            - 0.22 * snow_persistence
            - 0.16 * dry_penalty
        )
        vegetation_feedback = np.clip(vegetation_feedback, 0.0, 1.0) * land
        vegetation_feedback = self._smooth_field_masked(
            grid, vegetation_feedback, land, passes=1, alpha=0.10)
        vegetation_feedback = np.clip(vegetation_feedback, 0.0, 1.0) * land

        return {
            "seasonal_sea_ice": seasonal_sea_ice,
            "sea_ice": sea_ice,
            "seasonal_snow": seasonal_snow,
            "snow_persistence": snow_persistence,
            "seasonal_cloud_albedo_proxy": seasonal_cloud,
            "cloud_albedo_proxy": cloud_albedo,
            "vegetation_climate_feedback": vegetation_feedback,
        }

    def _seasonal_surface_flux(self, world, annual_flux):
        """Four-season daily-mean insolation proxy.

        The annual EBM still uses `stellar.surface_flux`; this field supplies
        the seasonal anomaly pattern.  It is normalized so the four-season,
        area-weighted mean matches the orbital energy budget.
        """
        grid = world.grid
        if world.spec.orbit.tidally_locked:
            return np.repeat(annual_flux[None, :], 4, axis=0)

        mean_flux = world.g(
            "stellar.mean_flux",
            float(np.average(annual_flux, weights=grid.cell_area)),
        )
        solar_constant = 4.0 * max(mean_flux, 0.0)
        eps = np.radians(world.g("orbit.obliquity", world.spec.orbit.obliquity_deg))
        lat = np.radians(grid.lat)
        sin_lat = np.sin(lat)
        cos_lat = np.cos(lat)
        phases = np.array([-0.5 * np.pi, 0.0, 0.5 * np.pi, np.pi])
        decl = np.arcsin(np.clip(np.sin(eps) * np.sin(phases), -1.0, 1.0))
        ecc = np.clip(world.g("orbit.eccentricity", world.spec.orbit.eccentricity),
                      0.0, 0.4)
        # With no longitude of perihelion in the spec, place perihelion near
        # DJF for Earth-like orbits and keep the modifier zero-mean.
        ecc_mod = 1.0 + 2.0 * ecc * np.cos(phases + 0.5 * np.pi)

        seasons = []
        for dec, emod in zip(decl, ecc_mod):
            x = -np.tan(lat) * np.tan(dec)
            h0 = np.arccos(np.clip(x, -1.0, 1.0))
            h0 = np.where(x <= -1.0, np.pi, h0)  # polar day
            h0 = np.where(x >= 1.0, 0.0, h0)     # polar night
            q = (solar_constant * emod / np.pi) * (
                h0 * sin_lat * np.sin(dec) + cos_lat * np.cos(dec) * np.sin(h0)
            )
            seasons.append(np.clip(q, 0.0, None))

        seasonal_flux = np.asarray(seasons, dtype=np.float64)
        global_mean = float(np.average(seasonal_flux.mean(axis=0), weights=grid.cell_area))
        if global_mean > 0.0:
            seasonal_flux *= mean_flux / global_mean
        return seasonal_flux

    def _continentality(self, grid, ocean):
        maritime = ocean.astype(np.float64)
        for _ in range(70):
            maritime = np.maximum(ocean.astype(np.float64),
                                  0.93 * self._neighbor_mean(grid, maritime))
        return np.where(ocean, 0.0, np.clip(1.0 - maritime, 0.0, 1.0))

    def _seasonal_temperature(self, world, annual_T, annual_flux, ocean):
        grid = world.grid
        continentality = self._continentality(grid, ocean)
        if world.spec.orbit.tidally_locked:
            seasonal = np.repeat(annual_T[None, :], 4, axis=0)
            return seasonal, np.zeros(grid.n, dtype=np.float64), continentality

        seasonal_flux = self._seasonal_surface_flux(world, annual_flux)
        flux_anomaly = seasonal_flux - seasonal_flux.mean(axis=0, keepdims=True)
        lat_abs = np.abs(grid.lat) / 90.0

        land_response = 0.035 + 0.030 * continentality + 0.008 * lat_abs
        ocean_response = 0.010 + 0.006 * lat_abs
        response = np.where(ocean, ocean_response, land_response)
        anomaly = flux_anomaly * response[None, :]

        # Oceans and strongly maritime land respond with a seasonal lag.
        maritime_land = np.where(ocean, 0.0, 1.0 - continentality)
        lagged = anomaly.copy()
        previous = np.roll(anomaly, shift=1, axis=0)
        lagged[:, ocean] = 0.55 * anomaly[:, ocean] + 0.45 * previous[:, ocean]
        if maritime_land.any():
            lag = 0.20 * maritime_land
            lagged[:, ~ocean] = (
                (1.0 - lag[~ocean])[None, :] * anomaly[:, ~ocean]
                + lag[~ocean][None, :] * previous[:, ~ocean]
            )

        for s in range(lagged.shape[0]):
            lagged[s] = self._smooth_field(grid, lagged[s], passes=2, alpha=0.12)

        land_cap = 13.0 + 14.0 * continentality + 5.0 * lat_abs
        ocean_cap = 5.0 + 5.0 * lat_abs
        cap = np.where(ocean, ocean_cap, land_cap)
        lagged = np.clip(lagged, -cap[None, :], cap[None, :])
        lagged -= lagged.mean(axis=0, keepdims=True)

        seasonal = np.clip(annual_T[None, :] + lagged, 120.0, 360.0)
        temp_seasonality = seasonal.max(axis=0) - seasonal.min(axis=0)
        return seasonal, temp_seasonality, continentality

    def _smooth_zonal_background(self, grid, T):
        """Smooth only the latitudinal background temperature.

        The local graph diffusion handles cell-scale noise, but the reduced EBM
        can still create a sharp zonal temperature wall where insolation,
        elevation and ice feedback line up.  We smooth the zonal mean and add
        back only that broad correction, preserving local terrain anomalies.
        """
        bins = np.linspace(-90.0, 90.0, 37)
        band = np.clip(np.digitize(grid.lat, bins) - 1, 0, len(bins) - 2)
        means = np.zeros(len(bins) - 1)
        for b in range(means.size):
            mask = band == b
            if mask.any():
                means[b] = np.average(T[mask], weights=grid.cell_area[mask])
            else:
                means[b] = means[b - 1] if b else float(np.average(T, weights=grid.cell_area))
        smooth = means.copy()
        for _ in range(5):
            prev = np.r_[smooth[0], smooth[:-1]]
            nxt = np.r_[smooth[1:], smooth[-1]]
            smooth = 0.25 * prev + 0.5 * smooth + 0.25 * nxt
        correction = np.clip(smooth - means, -10.0, 10.0)
        return T + 0.8 * correction[band]

    def _smooth_field(self, grid, values, passes: int = 1, alpha: float = 0.2):
        out = np.asarray(values, dtype=np.float64).copy()
        for _ in range(passes):
            out = (1.0 - alpha) * out + alpha * self._neighbor_mean(grid, out)
        return out

    def _climate_elevation(self, grid, elev):
        """Low-pass topography for climate-scale lapse/orographic effects."""
        out = self._smooth_field(grid, elev, passes=4, alpha=0.26)
        # Preserve the sign of broad ocean basins while removing one-cell peaks.
        return np.where(elev > 0.0, np.maximum(out, 0.0), np.minimum(out, 0.0))

    def _hydroclimate(self, world, T, ocean, wind, elev):
        grid = world.grid
        # Saturation-limited evaporation (Clausius-Clapeyron-ish), mm/yr.
        es = 6.11 * np.exp(17.27 * (T - 273.15) / np.maximum(T - 35.85, 1.0))
        evap = np.where(ocean, np.clip(2.0 * es, 0.0, 4000.0), 0.0)

        # Moisture field: ocean evaporation is the source, but land needs a
        # slowly decaying atmospheric reservoir plus recycling.  The previous
        # implementation pinned M to ocean values and left land with only a
        # weak diffusion tail, creating an all-desert Earth analogue.
        M = evap.copy()
        ocean_source = evap.copy()
        for _ in range(80):
            M = 0.72 * M + 0.28 * self._neighbor_mean(grid, M)
            M[ocean] = 0.65 * M[ocean] + 0.35 * ocean_source[ocean]

        # Graph-distance proxy for maritime influence.  It is deliberately
        # broad: this is a reduced climate model, so storm tracks/monsoons must
        # stand in for unresolved synoptic transport.
        maritime = ocean.astype(float)
        for _ in range(70):
            maritime = np.maximum(ocean.astype(float),
                                  0.93 * self._neighbor_mean(grid, maritime))

        lat_abs = np.abs(grid.lat)
        itcz = np.exp(-(grid.lat / 30.0) ** 2)
        storm_tracks = np.exp(-((lat_abs - 45.0) / 18.0) ** 2)
        monsoon = np.exp(-(grid.lat / 36.0) ** 2) * (1.0 - maritime)
        land_warmth = np.clip((T - 273.15) / 30.0, 0.0, 1.0)

        land_recycled = np.where(
            ~ocean,
            (260.0 + 920.0 * itcz + 600.0 * storm_tracks + 520.0 * monsoon)
            * (0.45 + 0.55 * maritime)
            * (0.45 + 0.55 * land_warmth),
            0.0,
        )

        # temperature controls how much falls out
        tfac = np.clip((T - 250.0) / 40.0, 0.05, 1.0)
        precip = (0.70 * M + land_recycled) * tfac

        # orographic enhancement / rain shadow along the wind.  Use the
        # climate-scale elevation field, not single-cell relief, so mountain
        # belts affect storm tracks as regions rather than as pixel stripes.
        i, j = self._edges[:, 0], self._edges[:, 1]
        dpos = grid.xyz[j] - grid.xyz[i]
        wmid = 0.5 * (wind[i] + wind[j])
        flow = np.einsum("ij,ij->i", wmid, dpos)
        dh = elev[j] - elev[i]
        oro = np.zeros(grid.n)
        # air moving uphill (flow>0 toward j higher) rains on the upwind/high side
        up = flow * dh
        np.add.at(oro, j, np.maximum(up, 0.0))
        np.add.at(oro, i, np.maximum(-up, 0.0))
        oro_norm = oro / np.maximum(oro.max(), 1e-9)
        precip = precip * (1.0 + 0.9 * oro_norm)
        # rain shadow: leeward descending air dries
        lee = np.zeros(grid.n)
        np.add.at(lee, i, np.maximum(up, 0.0))
        np.add.at(lee, j, np.maximum(-up, 0.0))
        lee_norm = lee / np.maximum(lee.max(), 1e-9)
        precip = precip * (1.0 - 0.35 * lee_norm)

        precip = np.where(ocean, 0.45 * precip + 520.0, precip)
        precip = self._smooth_field(grid, precip, passes=4, alpha=0.16)
        evap = np.where(~ocean, np.minimum(0.45 * precip, 900.0), evap)
        return evap, np.clip(precip, 0.0, 4000.0)
