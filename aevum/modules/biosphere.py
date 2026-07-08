"""Biosphere: functional groups, niches, dispersal, innovation & extinction.

We model functional groups (traits), not millions of species.  A group's niche
is set by climate/substrate/redox; its *distribution* is additionally limited by
dispersal and isolation (continent/ocean connectivity), so the same niche on two
separated continents need not be filled.

Evolutionary innovations are condition-gated stochastic events:

    P(innovation) ~ available_niche * total_biomass * env_duration * mutation
                    / complexity_cost

Mass extinctions emerge from environmental change rate, anoxia and impacts --
they are not pinned to Earth dates.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from aevum.core.events import Event
from aevum.core.module import Module, StepResult
from aevum.core.provenance import Provenance


class BiosphereModule(Module):
    name = "biosphere"
    produces = ["biosphere.functional_groups", "biosphere.lineages",
                "biosphere.npp", "biosphere.biome", "biosphere.distribution"]
    fidelity = "niche"
    interval_myr = 25.0

    INNOVATIONS = ["oxygenic_photosynthesis", "eukaryotic_complexity",
                   "multicellularity", "predation", "land_colonization",
                   "vascular_forests", "skeletons_reefs"]

    def init_state(self, world, rng_key) -> None:
        world.objects["biosphere.functional_groups"] = []
        world.objects["biosphere.lineages"] = []
        world.set_field("biosphere.npp", np.zeros(world.grid.n))
        world.set_field("biosphere.biome", np.zeros(world.grid.n))
        world.set_field("biosphere.richness", np.zeros(world.grid.n))
        world.set_g("bio.prev_meanT", 288.0)
        self._occ: dict[int, np.ndarray] = {}
        self._next_id = 0

    # ------------------------------------------------------------------
    def _components(self, grid, mask):
        edges = grid.edges
        sel = mask[edges[:, 0]] & mask[edges[:, 1]]
        e = edges[sel]
        if e.size == 0:
            data = np.ones(0)
            mat = coo_matrix((data, (np.zeros(0, int), np.zeros(0, int))),
                             shape=(grid.n, grid.n))
        else:
            data = np.ones(e.shape[0] * 2)
            rows = np.concatenate([e[:, 0], e[:, 1]])
            cols = np.concatenate([e[:, 1], e[:, 0]])
            mat = coo_matrix((data, (rows, cols)), shape=(grid.n, grid.n))
        ncomp, labels = connected_components(mat, directed=False)
        labels = labels.copy()
        labels[~mask] = -1
        return labels

    def _neighbor_mean(self, grid, f):
        edges = grid.edges
        i, j = edges[:, 0], edges[:, 1]
        acc = np.zeros_like(f, dtype=np.float64)
        deg = np.zeros(grid.n, dtype=np.float64)
        np.add.at(acc, i, f[j])
        np.add.at(acc, j, f[i])
        np.add.at(deg, i, 1.0)
        np.add.at(deg, j, 1.0)
        return acc / np.maximum(deg, 1.0)

    def _smooth_field(self, grid, values, passes: int = 1, alpha: float = 0.25):
        out = np.asarray(values, dtype=np.float64).copy()
        for _ in range(passes):
            out = (1.0 - alpha) * out + alpha * self._neighbor_mean(grid, out)
        return out

    def _new_group(self, world, t, name, parent, traits) -> dict:
        g = {"id": self._next_id, "name": name, "parent": parent,
             "birth_myr": round(t, 1), "death_myr": None, "alive": True,
             "complexity": traits.get("complexity", 1), **traits}
        self._next_id += 1
        world.objects["biosphere.functional_groups"].append(g)
        world.objects["biosphere.lineages"].append(dict(g))
        return g

    # ------------------------------------------------------------------
    def step(self, world, t, dt, forcing, rng_key) -> StepResult:
        spec = world.spec
        grid = world.grid
        groups = world.objects["biosphere.functional_groups"]
        origin = spec.biosphere.life_origin_myr

        if origin is None or t < origin:
            return StepResult(state_delta={"fields": {
                "biosphere.npp": np.zeros(grid.n),
                "biosphere.biome": self._biomes(world, np.zeros(grid.n)),
                "biosphere.richness": np.zeros(grid.n)}})

        T = self._smooth_field(
            grid, world.get_field("climate.surface_temperature", 288.0),
            passes=3, alpha=0.22,
        )
        precip = self._smooth_field(
            grid, world.get_field("climate.precipitation", 500.0),
            passes=4, alpha=0.24,
        )
        ocean = world.ocean_mask()
        land = ~ocean
        o2 = world.g("biogeochem.oxygen_fraction", 0.0)
        flags = world.globals

        # seed primordial life; re-seed if every group has gone extinct
        if not any(g["alive"] for g in groups):
            self._new_group(world, t, "primordial_microbes", None,
                            {"t_opt": 325.0, "t_tol": 70.0, "marine": True,
                             "terrestrial": False, "requires_o2": False,
                             "water_need": 0.0, "complexity": 1})

        ocean_lab = self._components(grid, ocean)
        land_lab = self._components(grid, land)

        # --- potential primary productivity (Miami-type) ----------------
        tc = T - 273.15
        npp_t = 1.0 / (1.0 + np.exp(1.315 - 0.119 * tc))
        npp_w = 1.0 - np.exp(-0.000664 * np.clip(precip, 0, None))
        dev_marine = 0.15 + 0.5 * world.g("bio.oxygenic_photosynthesis", 0.0) \
            + 0.35 * world.g("bio.multicellularity", 0.0)
        dev_land = world.g("bio.land_colonization", 0.0) * (
            0.4 + 0.6 * world.g("bio.vascular_forests", 0.0))
        npp = np.where(ocean, npp_t * dev_marine,
                       np.minimum(npp_t, npp_w) * dev_land) * 1.2
        npp = np.clip(npp, 0.0, 5.0)

        # --- per-group distribution (niche + dispersal) -----------------
        richness = np.zeros(grid.n)
        alive_groups = [g for g in groups if g["alive"]]
        rng = rng_key.child("dist").generator()
        for g in alive_groups:
            suit = self._suitability(g, T, precip, ocean, o2)
            domain_lab = ocean_lab if g["marine"] else land_lab
            mask_suit = suit > 0.25
            occ_prev = self._occ.get(g["id"])
            if occ_prev is None or not occ_prev.any():
                cand = np.where(mask_suit & (domain_lab >= 0))[0]
                if cand.size == 0:
                    self._occ[g["id"]] = np.zeros(grid.n, bool)
                    continue
                seed = int(rng.choice(cand))
                comps = {domain_lab[seed]}
            else:
                comps = set(np.unique(domain_lab[occ_prev & (domain_lab >= 0)]))
            in_comp = np.isin(domain_lab, list(comps)) if comps else np.zeros(grid.n, bool)
            occ = mask_suit & in_comp & (domain_lab >= 0)
            self._occ[g["id"]] = occ
            richness += occ

        # --- innovations -------------------------------------------------
        events: list[Event] = []
        biomass_total = float(np.sum(npp * grid.cell_area))
        events += self._maybe_innovate(world, t, dt, rng_key, biomass_total, o2, land)

        # --- extinctions -------------------------------------------------
        events += self._maybe_extinct(world, t, dt, rng_key, T)

        world.provenance.record(Provenance(
            "biosphere.npp", self.name, self.fidelity, "kg/m^2/yr",
            direct_cause=f"climate-limited productivity; ecosystem dev "
            f"marine={dev_marine:.2f} land={dev_land:.2f}; O2={o2:.3f}"))
        world.provenance.record(Provenance(
            "biosphere.distribution", self.name, self.fidelity, "1",
            direct_cause="niche suitability constrained by continent/ocean connectivity"))

        world.set_g("bio.prev_meanT", float(np.average(T, weights=grid.cell_area)))
        n_alive = sum(1 for g in groups if g["alive"])
        fields = {"biosphere.npp": npp, "biosphere.richness": richness,
                  "biosphere.biome": self._biomes(world, npp)}
        diag = {"n_groups_alive": n_alive, "max_richness": int(richness.max()),
                "biomass_index": round(biomass_total / 1e13, 3)}
        return StepResult(state_delta={"fields": fields}, events=events,
                          diagnostics=diag)

    # ------------------------------------------------------------------
    def _suitability(self, g, T, precip, ocean, o2):
        ts = np.exp(-((T - g["t_opt"]) / g["t_tol"]) ** 2)
        if g["marine"]:
            ws = np.where(ocean, 1.0, 0.0)
        else:
            need = g.get("water_need", 200.0)
            ws = np.where(~ocean, np.clip(precip / max(need, 1.0), 0.0, 1.0), 0.0)
        ok = 1.0 if (not g["requires_o2"] or o2 > 0.01) else 0.0
        return ts * ws * ok

    def _maybe_innovate(self, world, t, dt, rng_key, biomass, o2, land):
        rng = rng_key.child("innovate").generator()
        events: list[Event] = []
        done = {g["name"] for g in world.objects["biosphere.functional_groups"]}
        groups = world.objects["biosphere.functional_groups"]
        complexity_cost = 1.0 + 0.5 * max((g["complexity"] for g in groups), default=1)
        niche_space = 1.0  # simplified available ecological space
        base = niche_space * (biomass / 1e14) * (dt / 25.0) / complexity_cost
        for innov in self.INNOVATIONS:
            if innov in done or not self._innov_allowed(world, innov, o2, land):
                continue
            p = np.clip(base * 0.5, 0.0, 0.9)
            if rng.random() < p:
                traits, flag = self._innov_traits(world, innov, t)
                parent = groups[rng.integers(len(groups))]["id"] if groups else None
                self._new_group(world, t, innov, parent, traits)
                if flag:
                    world.set_g(flag, 1.0)
                events.append(Event("innovation", t, self.name, magnitude=1.0,
                                    params={"innovation": innov}))
                break   # at most one major innovation per step
        return events

    def _innov_allowed(self, world, innov, o2, land):
        spec = world.spec
        if innov == "oxygenic_photosynthesis":
            return spec.biosphere.allow_oxygenic_photosynthesis
        if innov in ("eukaryotic_complexity", "multicellularity", "predation"):
            return o2 > 0.005
        if innov in ("land_colonization", "vascular_forests"):
            return (spec.biosphere.allow_land_colonization and o2 > 0.01
                    and land.any())
        if innov == "skeletons_reefs":
            return o2 > 0.02
        return True

    def _innov_traits(self, world, innov, t):
        table = {
            "oxygenic_photosynthesis": ({"t_opt": 300.0, "t_tol": 40.0, "marine": True,
                                         "terrestrial": False, "requires_o2": False,
                                         "water_need": 0.0, "complexity": 1},
                                        "bio.oxygenic_photosynthesis"),
            "eukaryotic_complexity": ({"t_opt": 295.0, "t_tol": 30.0, "marine": True,
                                       "terrestrial": False, "requires_o2": True,
                                       "water_need": 0.0, "complexity": 2}, None),
            "multicellularity": ({"t_opt": 293.0, "t_tol": 25.0, "marine": True,
                                  "terrestrial": False, "requires_o2": True,
                                  "water_need": 0.0, "complexity": 3},
                                 "bio.multicellularity"),
            "predation": ({"t_opt": 293.0, "t_tol": 28.0, "marine": True,
                           "terrestrial": False, "requires_o2": True,
                           "water_need": 0.0, "complexity": 3}, None),
            "land_colonization": ({"t_opt": 295.0, "t_tol": 30.0, "marine": False,
                                   "terrestrial": True, "requires_o2": True,
                                   "water_need": 150.0, "complexity": 3},
                                  "bio.land_colonization"),
            "vascular_forests": ({"t_opt": 297.0, "t_tol": 22.0, "marine": False,
                                  "terrestrial": True, "requires_o2": True,
                                  "water_need": 600.0, "complexity": 4},
                                 "bio.vascular_forests"),
            "skeletons_reefs": ({"t_opt": 299.0, "t_tol": 12.0, "marine": True,
                                 "terrestrial": False, "requires_o2": True,
                                 "water_need": 0.0, "complexity": 4}, None),
        }
        return table[innov]

    def _maybe_extinct(self, world, t, dt, rng_key, T):
        rng = rng_key.child("extinct").generator()
        events: list[Event] = []
        grid = world.grid
        meanT = float(np.average(T, weights=grid.cell_area))
        dTdt = abs(meanT - world.g("bio.prev_meanT", meanT)) / max(dt, 1.0)
        impact_recent = (t - world.g("impacts.last_major_myr", -1e9)) < dt * 1.5
        stress = np.clip(dTdt * 3.0, 0.0, 1.0)
        if impact_recent:
            stress = max(stress, 0.6)
        casualties = 0
        for g in world.objects["biosphere.functional_groups"]:
            if not g["alive"]:
                continue
            occ = self._occ.get(g["id"])
            area_frac = (float((occ * grid.cell_area).sum() / grid.cell_area.sum())
                         if occ is not None else 0.0)
            p_ext = stress * (0.5 + 0.5 * (1.0 - min(area_frac * 20, 1.0)))
            if area_frac < 1e-4:
                p_ext = max(p_ext, 0.3)
            if rng.random() < p_ext * 0.5:
                g["alive"] = False
                g["death_myr"] = round(t, 1)
                for ln in world.objects["biosphere.lineages"]:
                    if ln["id"] == g["id"]:
                        ln["death_myr"] = round(t, 1)
                casualties += 1
        if casualties >= 3:
            events.append(Event("mass_extinction", t, self.name, magnitude=casualties,
                                params={"casualties": casualties,
                                        "dT_dt_K_per_Myr": round(dTdt, 3),
                                        "impact_driven": bool(impact_recent)}))
        elif casualties > 0:
            events.append(Event("extinction", t, self.name, magnitude=casualties,
                                params={"casualties": casualties}))
        return events

    def _biomes(self, world, npp):
        """0 ocean, 1 ice, 2 desert, 3 grassland, 4 forest, 5 tundra, 6 tropical."""
        grid = world.grid
        T = self._smooth_field(
            grid, world.get_field("climate.surface_temperature", 288.0),
            passes=3, alpha=0.22,
        )
        precip = self._smooth_field(
            grid, world.get_field("climate.precipitation", 500.0),
            passes=4, alpha=0.24,
        )
        seasonal_T = world.fields.get("climate.seasonal_temperature")
        seasonal_precip = world.fields.get("climate.seasonal_precipitation")
        if seasonal_T is not None and np.asarray(seasonal_T).shape == (4, grid.n):
            seasonal_T = np.asarray(seasonal_T, dtype=np.float64)
            coldest_c = np.min(seasonal_T, axis=0) - 273.15
            warmest_c = np.max(seasonal_T, axis=0) - 273.15
        else:
            coldest_c = T - 273.15
            warmest_c = T - 273.15
        if seasonal_precip is not None and np.asarray(seasonal_precip).shape == (4, grid.n):
            seasonal_precip = np.asarray(seasonal_precip, dtype=np.float64)
            driest = np.min(seasonal_precip, axis=0)
            wettest = np.max(seasonal_precip, axis=0)
            dry_seasons = np.sum(
                seasonal_precip < np.minimum(350.0, 0.45 * np.maximum(precip, 1e-9))[None, :],
                axis=0,
            )
        else:
            driest = precip
            wettest = precip
            dry_seasons = np.zeros(grid.n, dtype=int)
        ocean = world.ocean_mask()
        ice_sheet = world.get_field("cryosphere.ice_sheet", 0.0)
        ice = ice_sheet > 1200.0
        sea_ice = world.get_field("cryosphere.sea_ice", 0.0) > 0.5
        biome = np.zeros(grid.n)
        biome[ocean] = 0.0
        biome[ocean & sea_ice] = 1.0
        land = ~ocean
        tc = T - 273.15
        elevation = world.get_field("terrain.elevation_m", 0.0)
        sea_level = world.g(
            "ocean.sea_level_m",
            world.g("ocean.sea_level", 0.0),
        )
        relative_elev = np.asarray(elevation, dtype=np.float64) - float(sea_level)
        high_lat_cold = land & (np.abs(grid.lat) >= 60.0) & (
            (tc < 2.0) | (coldest_c < -8.0)
        )
        alpine_stress = land & (
            ((relative_elev >= 2000.0) & ((tc < 12.0) | (coldest_c < 3.0)))
            | ((relative_elev >= 3200.0) & (tc < 17.0))
        )
        cold_stress = land & (
            (tc < -5.0) | (coldest_c < -14.0) | high_lat_cold | alpine_stress
        )
        desert_precip_threshold = np.where(
            tc < 8.0,
            170.0,
            np.where(tc < 14.0, 210.0, 250.0),
        )
        arid = land & ~cold_stress & (
            (precip < desert_precip_threshold)
            | ((precip < 520.0) & (dry_seasons >= 3) & (warmest_c > 8.0))
            | ((precip < 420.0) & (tc > 18.0) & (driest < 180.0))
        )
        tropical = land & (tc >= 20.0) & (precip >= 840.0) & (wettest >= 850.0)
        seasonal_tropical_grass = tropical & (
            ((dry_seasons >= 3) & (precip < 1050.0))
            | ((driest < 180.0) & (precip < 1050.0))
        )
        forest_precip_threshold = np.where(
            tc < 8.0,
            300.0,
            np.where(tc < 14.0, 420.0, 520.0),
        )
        moist_temperate = (
            land
            & (precip >= forest_precip_threshold)
            & (tc >= -5.0)
            & (tc < 24.0)
        )
        warm_seasonal = (
            land
            & (precip >= 520.0)
            & (precip < 900.0)
            & (tc >= 24.0)
        )
        grass = (
            land
            & (precip >= desert_precip_threshold)
            & (precip < forest_precip_threshold)
            & (tc >= -5.0)
        )

        biome[cold_stress] = 5.0
        biome[arid] = 2.0
        biome[grass & ~arid] = 3.0
        biome[warm_seasonal & ~arid] = 3.0
        biome[moist_temperate & ~arid] = 4.0
        biome[tropical & ~seasonal_tropical_grass & ~arid] = 6.0
        biome[seasonal_tropical_grass & ~arid] = 3.0
        biome[ice] = 1.0
        return self._generalize_biomes(grid, biome, ocean, ice | (ocean & sea_ice))

    def _generalize_biomes(self, grid, biome, ocean, fixed):
        out = biome.astype(np.float64).copy()
        land_domain = ~np.asarray(ocean, dtype=bool)
        fixed = np.asarray(fixed, dtype=bool)
        candidate = land_domain & ~fixed
        neighbor_domain = candidate.copy()
        edges = np.asarray(grid.edges, dtype=int)
        edge_i = edges[:, 0]
        edge_j = edges[:, 1]
        for _ in range(2):
            codes = np.nan_to_num(out, nan=0.0).astype(int)
            max_code = max(7, int(codes.max(initial=0)) + 1)
            counts = np.zeros((grid.n, max_code), dtype=np.int16)
            for code in range(max_code):
                j_matches = neighbor_domain[edge_j] & (codes[edge_j] == code)
                if j_matches.any():
                    np.add.at(counts[:, code], edge_i[j_matches], 1)
                i_matches = neighbor_domain[edge_i] & (codes[edge_i] == code)
                if i_matches.any():
                    np.add.at(counts[:, code], edge_j[i_matches], 1)
            neighbor_counts = counts.sum(axis=1)
            best = np.argmax(counts, axis=1)
            best_counts = counts[np.arange(grid.n), best]
            update = candidate & (neighbor_counts >= 3) & (best_counts >= 3)
            semantic_keep = candidate & np.isin(codes, [4, 5, 6])
            update &= ~semantic_keep
            nxt = out.copy()
            nxt[update] = best[update].astype(np.float64)
            out = nxt
        out[ocean] = biome[ocean]
        out[fixed] = biome[fixed]
        return out
