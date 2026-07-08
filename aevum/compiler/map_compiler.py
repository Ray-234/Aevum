"""MapCompiler: compile the truth world into a strategy hex map.

Pipeline (project plan, section 8):
  1. resample physical fields to the hex grid;
  2. classify terrain by *relative* scale (not fixed Earth altitudes);
  3. extract rivers from the real drainage network;
  4. identify harbours from coast/depth/exposure;
  5. place resource tiles from deposit objects;
  6. compute food/production/trade yields, movement cost, hazards;
  7. choose fair starting positions (select, do not move mountains).

Balancing lives ONLY here; the truth layer is never mutated.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

from aevum.archive.world_archive import WorldArchive
from aevum.compiler.hexgrid import HexGrid
from aevum.core.state import WorldState

TERRAIN_NAMES = ["ocean", "coast", "plains", "hills", "mountains", "ice"]


@dataclass
class CompiledMap:
    hexgrid: HexGrid
    cell_index: np.ndarray            # (H, W) nearest sphere cell
    terrain: np.ndarray               # (H, W) terrain class
    elevation: np.ndarray
    temperature_C: np.ndarray
    precip: np.ndarray
    biome: np.ndarray
    river: np.ndarray                 # bool
    food: np.ndarray
    production: np.ndarray
    trade: np.ndarray
    hazard: np.ndarray
    resources: np.ndarray             # object dtype: per-hex commodity string or ""
    harbors: np.ndarray               # bool
    source_land_fraction: Optional[np.ndarray] = None
    source_shelf_fraction: Optional[np.ndarray] = None
    source_depth_province: Optional[np.ndarray] = None
    source_terrain_province: Optional[np.ndarray] = None
    source_continental_detail: Optional[np.ndarray] = None
    starts: list[tuple[int, int]] = field(default_factory=list)
    fairness: dict[str, Any] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)


class MapCompiler:
    def __init__(self, world: WorldState, archive: WorldArchive) -> None:
        self.world = world
        self.archive = archive

    def compile(self, width: int = 96, height: int = 48,
                n_starts: int = 6, tech_level: float = 1.0) -> CompiledMap:
        w = self.world
        grid = w.grid
        hg = HexGrid(width, height)
        LAT, LON = hg.latlon()
        sample_idx = self._hex_sample_indices(hg, grid, LAT, LON)
        idx = sample_idx[:, :, 0]

        def rs(name, default=0.0):
            return w.get_field(name, default)[idx]

        def rs_mean(name, default=0.0):
            return w.get_field(name, default)[sample_idx].mean(axis=2)

        def rs_mode(name, default=0.0):
            return self._sample_mode(w.get_field(name, default)[sample_idx].astype(int))

        sea = w.sea_level
        elev_samples = w.get_field("terrain.elevation_m")[sample_idx]
        elev_raw = elev_samples.mean(axis=2)
        center_rel = elev_samples[:, :, 0] - sea
        source_land_fraction = (elev_samples >= sea).mean(axis=2)
        depth_samples = w.get_field("ocean.depth_province", 0.0)[sample_idx].astype(int)
        terrain_samples = w.get_field("terrain.province", 0.0)[sample_idx].astype(int)
        detail_field = (
            "terrain.continental_detail_region_code"
            if "terrain.continental_detail_region_code" in w.fields
            else "terrain.continental_detail"
        )
        detail_samples = w.get_field(detail_field, 0.0)[sample_idx].astype(int)
        source_depth_province = self._sample_mode(depth_samples)
        source_terrain_province = self._sample_mode(terrain_samples)
        source_continental_detail = self._sample_mode(detail_samples)
        source_shelf_fraction = np.isin(depth_samples, [1, 2, 7]).mean(axis=2)

        elev = self._hex_smooth(hg, elev_raw, passes=2, alpha=0.28)
        temp = self._hex_smooth(
            hg, rs_mean("climate.surface_temperature", 288.0) - 273.15,
            passes=1, alpha=0.20,
        )
        precip = self._hex_smooth(
            hg, rs_mean("climate.precipitation", 0.0),
            passes=1, alpha=0.22,
        )
        biome = self._hex_majority(hg, rs_mode("biosphere.biome", 0.0), passes=1)
        npp = self._hex_smooth(hg, rs_mean("biosphere.npp", 0.0), passes=1, alpha=0.20)
        litho = rs_mode("terrain.lithology", 0.0)
        ice_sheet = rs_mean("cryosphere.ice_sheet", 0.0)
        sea_ice = rs_mean("cryosphere.sea_ice", 0.0)

        center_land = center_rel >= 0.0
        land = (source_land_fraction >= 0.45) | (
            center_land & (source_land_fraction >= 0.25)
        )
        ocean = ~land

        # --- terrain classification (relative) ---------------------------
        terrain = np.full((height, width), 2, dtype=int)        # plains default
        terrain[ocean] = 0
        coast = ocean & (
            self._adjacent_mask(hg, land)
            | (source_shelf_fraction >= 0.25)
            | ((source_land_fraction > 0.05) & (source_land_fraction < 0.45))
        )
        terrain[coast] = 1
        if land.any():
            land_elev = elev[land] - sea
            topo_score = self._hex_smooth(hg, elev - sea, passes=1, alpha=0.36)
            p65 = np.percentile(topo_score[land], 65)
            p92 = np.percentile(topo_score[land], 92)
            le = elev - sea
            hills_province = (
                np.isin(source_terrain_province, [3, 4, 5, 6, 8])
                | np.isin(source_continental_detail, [4, 5, 6, 7])
            )
            mountain_province = (
                np.isin(source_terrain_province, [5, 7])
                | (np.isin(source_continental_detail, [5, 6]) & (topo_score > p65))
            )
            terrain[land & ((topo_score > p65) | hills_province) & (le > 120.0)] = 3
            terrain[land & ((topo_score > p92) | mountain_province)
                    & ((le > 1100.0) | (land_elev.max() < 1200.0))] = 4
        terrain = self._generalize_terrain(hg, terrain)
        if land.any():
            le = elev - sea
            lowland_source = (
                land
                & np.isin(source_terrain_province, [1, 2, 3])
                & (le < 700.0)
            )
            terrain[lowland_source & (terrain == 4)] = 3
            terrain[lowland_source & (le < 300.0) & (terrain >= 3)] = 2
        terrain = self._prune_tiny_land_components(hg, terrain)
        shallow_water = (terrain <= 1) & (source_shelf_fraction >= 0.35)
        terrain[shallow_water] = 1
        land_ice = land & (ice_sheet > 1200.0) & (temp < -20.0)
        pack_ice = ocean & (sea_ice > 0.95) & (temp < -35.0)
        terrain[land_ice | pack_ice] = 5
        final_land = ((terrain >= 2) & (terrain <= 4)) | land_ice
        final_water = (terrain <= 1) | pack_ice
        elev = np.where(final_land, np.maximum(elev, sea + 20.0), elev)
        elev = np.where(final_water, np.minimum(elev, sea - 20.0), elev)

        # --- rivers ------------------------------------------------------
        river = self._rivers_on_hex(hg, idx, (terrain >= 2) & (terrain <= 4))

        # --- harbours: coast adjacent to land, sheltered ----------------
        harbors = self._harbors(hg, terrain)

        # --- yields ------------------------------------------------------
        food = self._food(terrain, npp, precip, temp)
        production = self._production(terrain, biome, litho, npp)
        trade = self._trade(terrain, river, harbors)
        hazard = self._hazard(idx)

        # --- resources ---------------------------------------------------
        resources = self._resources(idx, terrain, tech_level)

        cm = CompiledMap(
            hexgrid=hg, cell_index=idx, terrain=terrain, elevation=elev,
            temperature_C=temp, precip=precip, biome=biome, river=river,
            food=food, production=production, trade=trade, hazard=hazard,
            resources=resources, harbors=harbors,
            source_land_fraction=source_land_fraction,
            source_shelf_fraction=source_shelf_fraction,
            source_depth_province=source_depth_province,
            source_terrain_province=source_terrain_province,
            source_continental_detail=source_continental_detail,
            meta={"sea_level_m": sea,
                  "land_fraction": float(final_land.mean()),
                  "land_or_coast_fraction": float((terrain >= 1).mean()),
                  "source_land_fraction_mean": float(source_land_fraction.mean()),
                  "coast_fraction": float((terrain == 1).mean()),
                  "source_continental_detail_field": detail_field,
                  "width": width, "height": height})

        # --- start positions & fairness ---------------------------------
        cm.starts, cm.fairness = self._starts(cm, n_starts)
        return cm

    # ------------------------------------------------------------------
    def _hex_sample_indices(self, hg, grid, lat, lon) -> np.ndarray:
        """Sample center plus four sub-hex offsets for province-aware aggregation."""
        H, W = lat.shape
        dlat = 180.0 / max(H, 1)
        dlon = 360.0 / max(W, 1)
        offsets = [
            (0.0, 0.0),
            (-0.33 * dlat, 0.0),
            (0.33 * dlat, 0.0),
            (0.0, -0.33 * dlon),
            (0.0, 0.33 * dlon),
        ]
        samples = []
        for dla, dlo in offsets:
            la = np.clip(lat + dla, -89.999, 89.999)
            lo = ((lon + dlo + 180.0) % 360.0) - 180.0
            samples.append(grid.nearest_latlon(la.ravel(), lo.ravel()).reshape(H, W))
        return np.stack(samples, axis=2)

    def _sample_mode(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values)
        H, W, _ = values.shape
        out = np.zeros((H, W), dtype=values.dtype)
        for r in range(H):
            for c in range(W):
                labels, counts = np.unique(values[r, c], return_counts=True)
                out[r, c] = labels[int(np.argmax(counts))]
        return out

    def _land_component_sizes(self, hg, land):
        land = np.asarray(land, dtype=bool)
        sizes = np.zeros(land.shape, dtype=np.int32)
        H, W = land.shape
        seen = np.zeros_like(land, dtype=bool)
        for r in range(H):
            for c in range(W):
                if seen[r, c] or not land[r, c]:
                    continue
                stack = [(r, c)]
                seen[r, c] = True
                comp = []
                while stack:
                    rr, cc = stack.pop()
                    comp.append((rr, cc))
                    for nr, nc in hg.neighbors(rr, cc):
                        if land[nr, nc] and not seen[nr, nc]:
                            seen[nr, nc] = True
                            stack.append((nr, nc))
                size = len(comp)
                for rr, cc in comp:
                    sizes[rr, cc] = size
        return sizes

    def _prune_tiny_land_components(self, hg, terrain):
        """Drop isolated one- to few-hex islands from the playable abstraction."""
        out = terrain.copy()
        land = (out >= 2) & (out <= 4)
        if float(land.mean()) < 0.12:
            return out
        sizes = self._land_component_sizes(hg, land)
        component_sizes = sorted({int(s) for s in sizes[land]}, reverse=True)
        if not component_sizes:
            return out
        protected_cutoff = component_sizes[min(len(component_sizes), 6) - 1]
        min_size = max(4, int(0.0025 * terrain.size))
        tiny = land & (sizes < min_size) & (sizes < protected_cutoff)
        out[tiny] = 0
        coast = (out <= 1) & self._adjacent_mask(hg, (out >= 2) & (out <= 4))
        out[out <= 1] = 0
        out[coast] = 1
        return out

    def _hex_smooth(self, hg, values, passes: int = 1, alpha: float = 0.25):
        out = np.asarray(values, dtype=np.float64).copy()
        H, W = out.shape
        for _ in range(passes):
            mean = np.zeros_like(out)
            for r in range(H):
                for c in range(W):
                    nbs = hg.neighbors(r, c)
                    vals = [out[rr, cc] for rr, cc in nbs]
                    mean[r, c] = float(np.mean(vals)) if vals else out[r, c]
            out = (1.0 - alpha) * out + alpha * mean
        return out

    def _hex_majority(self, hg, values, passes: int = 1):
        out = np.asarray(values).copy()
        H, W = out.shape
        for _ in range(passes):
            nxt = out.copy()
            for r in range(H):
                for c in range(W):
                    vals = [out[r, c]] + [out[rr, cc] for rr, cc in hg.neighbors(r, c)]
                    labels, counts = np.unique(vals, return_counts=True)
                    nxt[r, c] = labels[int(np.argmax(counts))]
            out = nxt
        return out

    def _hex_majority_bool(self, hg, values, passes: int = 1):
        out = np.asarray(values, dtype=bool).copy()
        H, W = out.shape
        for _ in range(passes):
            nxt = out.copy()
            for r in range(H):
                for c in range(W):
                    vals = [out[r, c]] + [out[rr, cc] for rr, cc in hg.neighbors(r, c)]
                    nxt[r, c] = sum(bool(v) for v in vals) >= (len(vals) / 2.0)
            out = nxt
        return out

    def _adjacent_mask(self, hg, mask):
        mask = np.asarray(mask, dtype=bool)
        out = np.zeros_like(mask)
        H, W = mask.shape
        for r in range(H):
            for c in range(W):
                out[r, c] = any(mask[rr, cc] for rr, cc in hg.neighbors(r, c))
        return out

    def _generalize_terrain(self, hg, terrain):
        out = terrain.copy()
        H, W = out.shape
        for _ in range(2):
            nxt = out.copy()
            for r in range(H):
                for c in range(W):
                    cur = int(out[r, c])
                    if cur == 5:
                        continue
                    nvals = np.array([out[rr, cc] for rr, cc in hg.neighbors(r, c)], dtype=int)
                    if cur <= 1:
                        cand = nvals[nvals <= 1]
                    else:
                        cand = nvals[(nvals >= 2) & (nvals <= 4)]
                    if cand.size < 3:
                        continue
                    labels, counts = np.unique(cand, return_counts=True)
                    k = int(np.argmax(counts))
                    if counts[k] >= 3:
                        nxt[r, c] = int(labels[k])
            out = nxt
        land = (out >= 2) & (out <= 4)
        coast = (out <= 1) & self._adjacent_mask(hg, land)
        out[out <= 1] = 0
        out[coast] = 1
        return out

    def _rivers_on_hex(self, hg, idx, land):
        w = self.world
        rivers = w.networks.get("hydrology.rivers")
        river = np.zeros(idx.shape, dtype=bool)
        if rivers is None:
            return river
        rc = set(int(x) for x in rivers.get("river_cells", []))
        rank_field = rivers.get("river_rank")
        flat = idx.ravel()
        mask = np.array([c in rc for c in flat]).reshape(idx.shape)
        if rank_field is not None:
            rank = rank_field[idx]
            mask &= (rank > 0.18)
        river = mask & land
        # Bridge one-hex gaps introduced by nearest-cell resampling so rivers
        # read as short connected drainage traces rather than isolated pixels.
        H, W = river.shape
        bridged = river.copy()
        for r in range(H):
            for c in range(W):
                if river[r, c] or not land[r, c]:
                    continue
                n = sum(1 for rr, cc in hg.neighbors(r, c) if river[rr, cc])
                if n >= 2:
                    bridged[r, c] = True
        return bridged

    def _harbors(self, hg, terrain):
        H, W = terrain.shape
        harbors = np.zeros_like(terrain, dtype=bool)
        land = terrain >= 2
        coast = terrain == 1
        for r in range(H):
            for c in range(W):
                if not coast[r, c]:
                    continue
                # Harbour if coast touches land, including across the longitude
                # seam.  The game map is a wrapped projection of a sphere; the
                # left and right edges are neighbours, not map borders.
                if any(land[rr, cc] for rr, cc in hg.neighbors(r, c)):
                    harbors[r, c] = True
        return harbors

    def _food(self, terrain, npp, precip, temp):
        # Climate-based fertility baseline (so playability does not depend on the
        # biosphere model being complete) plus a productivity bonus.
        tfit = np.exp(-((temp - 24.0) / 30.0) ** 2)            # broad warm-temperate optimum
        wfit = np.clip(precip / 650.0, 0.0, 1.0)
        base = 2.4 * tfit * wfit
        food = np.clip(base + 0.8 * npp, 0.0, 3.0)
        food[terrain == 0] = 0.8                               # ocean fishing
        food[terrain == 1] = np.maximum(food[terrain == 1], 1.2)
        food[terrain == 4] = np.minimum(food[terrain == 4], 0.5)
        food[terrain == 5] = 0.0
        return np.round(food, 2)

    def _production(self, terrain, biome, litho, npp):
        prod = np.zeros_like(npp)
        prod[terrain == 3] = 1.5                            # hills mining
        prod[terrain == 4] = 1.8
        prod[biome == 4] += 1.2                             # forests
        prod[biome == 6] += 1.0
        prod[litho == 0] += 0.4                             # igneous (ores)
        prod[terrain == 0] = 0.2
        return np.round(np.clip(prod, 0, 4), 2)

    def _trade(self, terrain, river, harbors):
        trade = np.zeros(terrain.shape)
        trade[terrain == 1] = 1.0
        trade[river] += 1.5
        trade[harbors] += 1.5
        trade[terrain == 0] += 0.5
        return np.round(np.clip(trade, 0, 4), 2)

    def _hazard(self, idx):
        w = self.world
        volc = w.get_field("tectonics.volcanism_age_myr", -1.0)[idx]
        t = w.time_myr
        hazard = np.zeros(idx.shape)
        hazard[volc >= (t - 100)] += 1.0
        boundaries = w.networks.get("tectonics.boundaries", {})
        conv = set(int(x) for x in boundaries.get("convergent", []))
        flat = idx.ravel()
        quake = np.array([c in conv for c in flat]).reshape(idx.shape)
        hazard[quake] += 1.0
        return np.round(np.clip(hazard, 0, 3), 2)

    def _resources(self, idx, terrain, tech_level):
        w = self.world
        H, Wd = idx.shape
        res = np.empty((H, Wd), dtype=object)
        res[:] = ""
        # map each deposit's sphere cell to the nearest hex
        cell_to_hex: dict[int, tuple[int, int]] = {}
        for r in range(H):
            for c in range(Wd):
                cell_to_hex.setdefault(int(idx[r, c]), (r, c))
        candidates: dict[tuple[int, int], tuple[float, str]] = {}
        for dep in w.objects.get("resources.deposits", []):
            hx = cell_to_hex.get(dep["cell"])
            if hx is None:
                continue
            r, c = hx
            if terrain[r, c] not in (1, 2, 3, 4):
                continue
            # accessibility: preservation x exposure(shallow burial) x tech
            burial = dep["burial_depth_m"]
            exposure = np.clip(1.0 - burial / (2000.0 * tech_level), 0.0, 1.0)
            access = dep["preservation"] * exposure
            if access < 0.35:
                continue
            commodity = max(dep["commodity_vector"].items(), key=lambda kv: kv[1])[0]
            score = access * max(dep.get("grade", 0.1), 0.1) * max(dep.get("tonnage_index", 0.1), 0.1)
            if score > candidates.get(hx, (-1.0, ""))[0]:
                candidates[hx] = (float(score), commodity)
        landish = int(np.isin(terrain, [2, 3, 4]).sum())
        cap = max(8, min(len(candidates), int(0.02 * max(landish, 1))))
        for (r, c), (_, commodity) in sorted(candidates.items(), key=lambda kv: kv[1][0], reverse=True)[:cap]:
            res[r, c] = commodity
        return res

    # ------------------------------------------------------------------
    def _starts(self, cm: CompiledMap, n_starts: int):
        """Farthest-point selection over habitable hexes (select, don't sculpt)."""
        if n_starts <= 0:
            return [], {"note": "starts disabled", "n_starts": 0}
        H, W = cm.terrain.shape
        land = (cm.terrain >= 2) & (cm.terrain <= 3)
        supported_land = land.copy()
        if float(land.mean()) >= 0.12:
            sizes = self._land_component_sizes(cm.hexgrid, land)
            supported_land &= sizes >= max(5, int(0.004 * H * W))
        # Relax the fertility bar on harsh worlds until enough sites exist.
        coords = np.empty((0, 2), int)
        for thr in (0.8, 0.6, 0.4, 0.2, 0.0):
            habitable = supported_land & (cm.food > thr) & (cm.hazard <= 1.0)
            coords = np.argwhere(habitable)
            if coords.shape[0] >= n_starts:
                break
        if coords.shape[0] < n_starts:
            habitable = land & (cm.food > 0.0)
            coords = np.argwhere(habitable)
        if coords.shape[0] < n_starts:
            return [], {"note": "too few habitable hexes", "n_habitable": int(coords.shape[0])}

        # seed at the most fertile hex, then greedily add farthest points
        base_scores = (cm.food + 0.5 * cm.production + 0.3 * cm.trade
                       - 0.9 * cm.hazard - 0.25 * (cm.terrain == 3))
        scores = self._local_start_scores(cm, base_scores)
        coord_scores = np.array([scores[r, c] for r, c in coords])
        quality_floor = np.percentile(coord_scores, 50)
        quality = np.array([scores[r, c] >= quality_floor for r, c in coords])
        if int(quality.sum()) >= n_starts:
            coords = coords[quality]
            coord_scores = coord_scores[quality]

        # recover actual coordinate of best
        best = coords[int(np.argmax(coord_scores))]
        start_list = [(int(best[0]), int(best[1]))]

        def torus_dist(a, b):
            dr = abs(a[0] - b[0])
            dc = abs(a[1] - b[1])
            dc = min(dc, W - dc)
            return (dr ** 2 + dc ** 2) ** 0.5

        while len(start_list) < n_starts:
            best_c, best_d = None, -1
            for rc in coords:
                rc = (int(rc[0]), int(rc[1]))
                d = min(torus_dist(rc, s) for s in start_list)
                d *= (0.5 + max(scores[rc[0], rc[1]], 0.0))
                if d > best_d:
                    best_d, best_c = d, rc
            start_list.append(best_c)

        yields = np.array([cm.food[r, c] + cm.production[r, c] + cm.trade[r, c]
                           for r, c in start_list])
        local_yields = np.array([scores[r, c] for r, c in start_list])
        if len(start_list) > 1:
            dists = [
                min(torus_dist(a, b) for j, b in enumerate(start_list) if j != i)
                for i, a in enumerate(start_list)
            ]
        else:
            dists = []
        min_separation = float(min(dists)) if dists else 0.0
        fairness = {
            "n_starts": len(start_list),
            "yield_min": round(float(yields.min()), 2),
            "yield_max": round(float(yields.max()), 2),
            "yield_cv": round(float(yields.std() / max(yields.mean(), 1e-6)), 3),
            "local_yield_min": round(float(local_yields.min()), 2),
            "local_yield_cv": round(float(local_yields.std() / max(local_yields.mean(), 1e-6)), 3),
            "min_separation": round(min_separation, 2),
        }
        return start_list, fairness

    def _local_start_scores(self, cm: CompiledMap, base_scores: np.ndarray) -> np.ndarray:
        """Blend tile score with one- and two-ring neighbourhood support."""
        H, W = cm.terrain.shape
        out = np.zeros_like(base_scores, dtype=float)
        for r in range(H):
            for c in range(W):
                vals = [base_scores[r, c]]
                ring1 = cm.hexgrid.neighbors(r, c)
                vals.extend(base_scores[rr, cc] for rr, cc in ring1)
                for rr, cc in ring1:
                    vals.extend(0.5 * base_scores[rr2, cc2]
                                for rr2, cc2 in cm.hexgrid.neighbors(rr, cc))
                out[r, c] = float(np.mean(vals))
        return out

    # ------------------------------------------------------------------
    def explain(self, cm: CompiledMap, r: int, c: int) -> dict:
        cell = int(cm.cell_index[r, c])
        story = self.archive.explain_cell(cell)
        story["hex"] = {"row": r, "col": c,
                        "terrain": TERRAIN_NAMES[int(cm.terrain[r, c])],
                        "food": float(cm.food[r, c]),
                        "production": float(cm.production[r, c]),
                        "trade": float(cm.trade[r, c]),
                        "hazard": float(cm.hazard[r, c]),
                        "resource": cm.resources[r, c]}
        return story
