"""Tectonics & crust (spherical plates + event rules).

v0.1 scheme -- a Lagrangian advection rasterised onto the fixed grid:

  1. every plate carries an Euler pole + rotation rate (scaled by mantle vigour);
  2. each step we rotate every crust parcel by its plate's increment;
  3. we rasterise the moved parcels back onto the fixed grid:
       * a fixed cell with no parcel nearby  -> new oceanic crust (sea-floor
         spreading at a divergent ridge, age = 0);
       * a cell reached by parcels of >=2 plates -> a convergent boundary:
            continental+continental -> collision / orogeny (crust thickens),
            continental+oceanic      -> subduction (volcanic arc),
            oceanic+oceanic          -> island arc.

This produces drifting continents, opening/closing oceans, mountain belts and
arcs -- and emits rift / subduction / collision / reorganization events so every
landform has a cause.  It is NOT GPlates; it is a tractable causal proxy.
"""
from __future__ import annotations

import heapq

import numpy as np
from scipy.spatial import cKDTree

from aevum.core.events import Event
from aevum.core.module import Module, StepResult
from aevum.core.provenance import Provenance
from aevum.modules.geometry import random_unit_vectors, rotation_matrix

CONT, OCEAN = 1.0, 0.0
ORIGIN_RIDGE = 0.0
ORIGIN_PRIMORDIAL = 1.0
ORIGIN_ARC = 2.0
ORIGIN_SUTURE = 3.0
ORIGIN_PLUME_IMPACT = 4.0
ORIGIN_CRATON = 5.0
PROVINCE_NONE = 0.0
PROVINCE_SHIELD = 1.0
PROVINCE_PLATFORM = 2.0
PROVINCE_INTERIOR_BASIN = 3.0
INTERNAL_BLOCK_NONE = 0.0
INTERNAL_BLOCK_CRATON_CORE = 1.0
INTERNAL_BLOCK_STABLE_PLATFORM = 2.0
INTERNAL_BLOCK_INTRACRATONIC_BASIN = 3.0
INTERNAL_BLOCK_MOBILE_BELT = 4.0
INTERNAL_BLOCK_RIFTED_MARGIN = 5.0
INTERNAL_BLOCK_ACCRETED_TERRANE = 6.0
INTERNAL_BLOCK_NAMES = {
    int(INTERNAL_BLOCK_NONE): "none",
    int(INTERNAL_BLOCK_CRATON_CORE): "craton_core",
    int(INTERNAL_BLOCK_STABLE_PLATFORM): "stable_platform",
    int(INTERNAL_BLOCK_INTRACRATONIC_BASIN): "intracratonic_basin",
    int(INTERNAL_BLOCK_MOBILE_BELT): "mobile_belt",
    int(INTERNAL_BLOCK_RIFTED_MARGIN): "rifted_margin",
    int(INTERNAL_BLOCK_ACCRETED_TERRANE): "accreted_terrane",
}
CONT_THICK = 36000.0
OCEAN_THICK = 7000.0
COLLIDE_UPLIFT = 900.0      # m of crustal thickening per Myr at collisions
ARC_UPLIFT = 220.0          # m / Myr at subduction arcs
MAX_CONT_THICK = 70000.0
DOMAIN_OCEANIC = 0.0
DOMAIN_CRATON = 1.0
DOMAIN_CONTINENTAL_INTERIOR = 2.0
DOMAIN_CONTINENTAL_MARGIN = 3.0
DOMAIN_ACCRETED_TERRANE = 4.0
DOMAIN_SUTURE = 5.0
DOMAIN_LIP = 6.0
WILSON_INACTIVE = 0.0
WILSON_RIFT = 1.0
WILSON_OPENING = 2.0
WILSON_MATURE = 3.0
WILSON_CLOSING = 4.0
WILSON_COLLISION = 5.0
WILSON_SUTURE = 6.0
DEFORM_NONE = 0.0
DEFORM_COLLISION_CORE = 1.0
DEFORM_COLLISION_SHOULDER = 2.0
DEFORM_SUBDUCTION_CORE = 3.0
DEFORM_SUBDUCTION_SHOULDER = 4.0
DEFORM_RIFT = 5.0
DEFORM_TRANSFORM = 6.0
P110B_HISTORICAL_BREAKUP_SEAWAY_RETENTION_MYR = 650.0
P110B_SUPERCONTINENT_RESIDENCE_DEBT_START_MYR = 900.0
P110B_SUPERCONTINENT_RESIDENCE_FULL_PRESSURE_MYR = 500.0

P39_OLD_PLATFORM_FUNNEL_METRICS = (
    "continental_cells",
    "legacy_anchor_cells",
    "source_ok_cells",
    "source_blocked_cells",
    "thick_ok_cells",
    "thick_blocked_cells",
    "quiet_cells",
    "recently_reworked_cells",
    "age_mature_cells",
    "age_old_cells",
    "age_mature_blocked_cells",
    "age_old_blocked_cells",
    "width_ge3_cells",
    "width_ge4_cells",
    "width_ge3_blocked_cells",
    "width_ge4_blocked_cells",
    "stability_ge030_cells",
    "stability_ge036_cells",
    "stability_ge030_blocked_cells",
    "stability_ge036_blocked_cells",
    "old_platform_core_cells",
    "old_platform_margin_cells",
    "mature_platform_cells",
    "process_anchor_cells",
    "lineage_broad_cells",
    "process_anchor_lineaged_cells",
    "process_anchor_unlineaged_cells",
)


def _continent_nucleus_count(spec) -> int:
    """Choose continent nuclei count from land budget and plate count."""
    target = float(spec.target_land_fraction)
    if target <= 0.08:
        return max(1, min(3, spec.n_plates // 4))
    if target >= 0.62:
        return max(3, min(7, spec.n_plates // 2))
    return max(3, min(4, spec.n_plates // 4))


def _unit_vector_from_lat_lon(lat_deg: float, lon_deg: float) -> np.ndarray:
    lat = np.radians(float(lat_deg))
    lon = np.radians(float(lon_deg))
    return np.array([
        np.cos(lat) * np.cos(lon),
        np.cos(lat) * np.sin(lon),
        np.sin(lat),
    ], dtype=np.float64)


P52_PROTO_CRUST_AXES = np.asarray([
    _unit_vector_from_lat_lon(-18.0, -118.0),
    _unit_vector_from_lat_lon(24.0, 24.0),
    _unit_vector_from_lat_lon(-38.0, 112.0),
    _unit_vector_from_lat_lon(56.0, -44.0),
    _unit_vector_from_lat_lon(8.0, 164.0),
    _unit_vector_from_lat_lon(46.0, 82.0),
    _unit_vector_from_lat_lon(-6.0, -8.0),
], dtype=np.float64)

P52_PROTO_CRUST_WEIGHTS = np.asarray(
    [1.00, 0.96, 0.91, 0.78, 0.66, 0.58, 0.52], dtype=np.float64)


def _normalize01(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    lo = float(np.min(values)) if values.size else 0.0
    hi = float(np.max(values)) if values.size else 0.0
    return (values - lo) / max(hi - lo, 1.0e-12)


def _stable_tie_breaker(cells: np.ndarray | int) -> np.ndarray | float:
    arr = np.asarray(cells, dtype=np.int64)
    value = ((arr * 1103515245 + 12345) & 0x7FFFFFFF).astype(np.float64)
    value = value / float(0x7FFFFFFF)
    return float(value) if np.isscalar(cells) else value


def _proto_crust_potential(grid, spec=None) -> np.ndarray:
    """Deterministic low-degree proxy for early buoyant continental nuclei.

    This is not a mantle-convection solver.  It is the initial-condition
    counterpart to the later R1 potential fields: broad low-degree mantle/
    lithosphere structure creates a few separated proto-crust highs, with water
    inventory and land budget only modulating amplitude.  No random stream is
    consumed, so initial continental cargo is reproducible from physical
    parameters and grid geometry rather than seed jitter.
    """
    dots = grid.xyz @ P52_PROTO_CRUST_AXES.T
    domes = np.maximum(dots, 0.0) ** 2.6
    field = domes @ P52_PROTO_CRUST_WEIGHTS
    harmonic = (
        0.18 * np.sin(2.4 * dots[:, 1] - 0.3)
        + 0.13 * np.cos(2.1 * dots[:, 3] + 0.8)
        + 0.08 * np.sin(3.0 * dots[:, 5] + 1.7)
    )
    mid_lat = np.clip(1.0 - (np.abs(grid.lat) / 82.0) ** 2, 0.0, 1.0)
    target = float(getattr(spec, "target_land_fraction", 0.29)) if spec is not None else 0.29
    water = 1.0
    if spec is not None and getattr(spec, "composition", None) is not None:
        water = float(getattr(spec.composition, "water_inventory_earth", 1.0))
    dry_bias = np.clip((0.8 - water) / 1.6, -0.25, 0.35)
    field = field + harmonic + (0.10 + 0.08 * dry_bias) * mid_lat
    field += 0.06 * np.cos(np.radians(grid.lon * (1.0 + 0.3 * target)))
    return _normalize01(field)


def _deterministic_margin_noise(grid, spec=None) -> np.ndarray:
    dots = grid.xyz @ P52_PROTO_CRUST_AXES.T
    target = float(getattr(spec, "target_land_fraction", 0.29)) if spec is not None else 0.29
    noise = (
        0.44 * np.sin(3.1 * dots[:, 0] + 0.30)
        + 0.31 * np.cos(2.7 * dots[:, 2] - 0.85)
        + 0.23 * np.sin((2.2 + target) * dots[:, 4] + 1.40)
        + 0.14 * np.cos(np.radians(2.0 * grid.lat - 0.7 * grid.lon))
    )
    return _normalize01(noise)


def _separated_seed_cells(grid, n_seeds: int, rng=None, *, avoid_poles: bool = True,
                          potential: np.ndarray | None = None) -> np.ndarray:
    """Deterministic potential maxima with farthest-point spacing."""
    del rng
    n_seeds = max(1, min(int(n_seeds), grid.n))
    candidates = np.arange(grid.n)
    if avoid_poles and grid.n > n_seeds * 8:
        midlat = np.where(np.abs(grid.lat) < 72.0)[0]
        if midlat.size >= n_seeds:
            candidates = midlat
    score0 = (
        _proto_crust_potential(grid) if potential is None
        else _normalize01(np.asarray(potential, dtype=np.float64))
    )
    score0 = score0 + 1.0e-9 * _stable_tie_breaker(np.arange(grid.n))
    first = int(candidates[int(np.argmax(score0[candidates]))])
    seeds = [first]
    max_dot = grid.xyz[candidates] @ grid.xyz[first]
    for _ in range(1, n_seeds):
        spacing = np.arccos(np.clip(max_dot, -1.0, 1.0)) / np.pi
        score = (
            0.58 * score0[candidates]
            + 0.42 * spacing
            + 1.0e-9 * _stable_tie_breaker(candidates)
        )
        too_close = spacing < max(0.18, 0.46 / max(np.sqrt(n_seeds), 1.0))
        score[too_close] -= 0.75
        nxt = int(candidates[int(np.argmax(score))])
        seeds.append(nxt)
        max_dot = np.maximum(max_dot, grid.xyz[candidates] @ grid.xyz[nxt])
    return np.asarray(seeds, dtype=np.int64)


def _plate_seed_cells(grid, n_plates: int, continent_seeds: np.ndarray,
                      rng=None) -> np.ndarray:
    """Choose plate seeds while keeping continent nuclei inside plate interiors."""
    del rng
    n_plates = max(1, min(int(n_plates), grid.n))
    selected: list[int] = []
    for seed in np.asarray(continent_seeds, dtype=np.int64):
        seed = int(seed)
        if seed not in selected:
            selected.append(seed)
        if len(selected) >= n_plates:
            return np.asarray(selected[:n_plates], dtype=np.int64)
    candidates = np.arange(grid.n)
    if selected:
        max_dot = np.max(grid.xyz[candidates] @ grid.xyz[np.asarray(selected)].T, axis=1)
    else:
        potential = _proto_crust_potential(grid)
        first = int(candidates[int(np.argmax(potential[candidates]))])
        selected.append(first)
        max_dot = grid.xyz[candidates] @ grid.xyz[first]
    plate_potential = _normalize01(
        0.55 * (1.0 - _proto_crust_potential(grid))
        + 0.25 * np.abs(np.sin(np.radians(grid.lon * 1.7)))
        + 0.20 * np.clip(1.0 - (np.abs(grid.lat) / 86.0) ** 2, 0.0, 1.0)
    )
    while len(selected) < n_plates:
        spacing = np.arccos(np.clip(max_dot, -1.0, 1.0)) / np.pi
        score = (
            0.70 * spacing
            + 0.30 * plate_potential[candidates]
            + 1.0e-9 * _stable_tie_breaker(candidates)
        )
        score[np.asarray(selected, dtype=np.int64)] = -np.inf
        nxt = int(candidates[int(np.argmax(score))])
        selected.append(nxt)
        max_dot = np.maximum(max_dot, grid.xyz[candidates] @ grid.xyz[nxt])
    return np.asarray(selected, dtype=np.int64)


def _legacy_random_separated_seed_cells(grid, n_seeds: int, rng,
                                        *, avoid_poles: bool = True) -> np.ndarray:
    """Legacy farthest-point sample with a random first seed."""
    n_seeds = max(1, min(int(n_seeds), grid.n))
    candidates = np.arange(grid.n)
    if avoid_poles and grid.n > n_seeds * 8:
        midlat = np.where(np.abs(grid.lat) < 72.0)[0]
        if midlat.size >= n_seeds:
            candidates = midlat
    first = int(rng.choice(candidates))
    seeds = [first]
    min_dot = grid.xyz[candidates] @ grid.xyz[first]
    for _ in range(1, n_seeds):
        score = -min_dot + 0.015 * rng.random(candidates.size)
        nxt = int(candidates[int(np.argmax(score))])
        seeds.append(nxt)
        min_dot = np.maximum(min_dot, grid.xyz[candidates] @ grid.xyz[nxt])
    return np.asarray(seeds, dtype=np.int64)


def _legacy_random_plate_seed_cells(grid, n_plates: int,
                                    continent_seeds: np.ndarray, rng) -> np.ndarray:
    """Legacy plate seeds that preserve continent nuclei as first plate seeds."""
    n_plates = max(1, min(int(n_plates), grid.n))
    selected: list[int] = []
    for seed in np.asarray(continent_seeds, dtype=np.int64):
        seed = int(seed)
        if seed not in selected:
            selected.append(seed)
        if len(selected) >= n_plates:
            return np.asarray(selected[:n_plates], dtype=np.int64)
    candidates = np.arange(grid.n)
    if selected:
        min_dot = np.max(grid.xyz[candidates] @ grid.xyz[np.asarray(selected)].T, axis=1)
    else:
        first = int(rng.choice(candidates))
        selected.append(first)
        min_dot = grid.xyz[candidates] @ grid.xyz[first]
    while len(selected) < n_plates:
        score = -min_dot + 0.01 * rng.random(candidates.size)
        score[np.asarray(selected, dtype=np.int64)] = -np.inf
        nxt = int(candidates[int(np.argmax(score))])
        selected.append(nxt)
        min_dot = np.maximum(min_dot, grid.xyz[candidates] @ grid.xyz[nxt])
    return np.asarray(selected, dtype=np.int64)


def _graph_width_steps(grid, mask: np.ndarray) -> np.ndarray:
    """Graph-distance width proxy inside a boolean mask."""
    mask = np.asarray(mask, dtype=bool)
    width = np.zeros(grid.n, dtype=np.float64)
    if not mask.any():
        return width
    boundary = np.zeros(grid.n, dtype=bool)
    for c in np.where(mask)[0]:
        if (~mask[grid.neighbors[int(c)]]).any():
            boundary[int(c)] = True
    if not boundary.any():
        width[mask] = np.sqrt(float(mask.sum()))
        return width
    queue = [int(c) for c in np.where(boundary)[0]]
    seen = np.zeros(grid.n, dtype=bool)
    seen[queue] = True
    width[queue] = 1.0
    head = 0
    while head < len(queue):
        c = queue[head]
        head += 1
        for nb in grid.neighbors[c]:
            nb = int(nb)
            if not mask[nb] or seen[nb]:
                continue
            width[nb] = width[c] + 1.0
            seen[nb] = True
            queue.append(nb)
    return width


def _same_neighbor_count(grid, mask: np.ndarray) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    if not grid.edges.size:
        return np.zeros(grid.n, dtype=np.int16)
    i, j = grid.edges[:, 0], grid.edges[:, 1]
    hits_i = i[mask[j]]
    hits_j = j[mask[i]]
    if hits_i.size and hits_j.size:
        hits = np.concatenate((hits_i, hits_j))
    elif hits_i.size:
        hits = hits_i
    elif hits_j.size:
        hits = hits_j
    else:
        return np.zeros(grid.n, dtype=np.int16)
    return np.bincount(hits, minlength=int(grid.n)).astype(np.int16)


def _component_labels(grid, mask: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    mask = np.asarray(mask, dtype=bool)
    labels = np.full(grid.n, -1, dtype=np.int64)
    comps: list[np.ndarray] = []
    for start in np.where(mask)[0]:
        if labels[start] >= 0:
            continue
        cid = len(comps)
        stack = [int(start)]
        labels[start] = cid
        cells: list[int] = []
        while stack:
            c = stack.pop()
            cells.append(c)
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if mask[nb] and labels[nb] < 0:
                    labels[nb] = cid
                    stack.append(nb)
        comps.append(np.asarray(cells, dtype=np.int64))
    return labels, comps


def _stable_component_ids(grid, mask: np.ndarray, previous: np.ndarray,
                          next_id: int) -> tuple[np.ndarray, list[np.ndarray], int]:
    """Assign stable numeric ids to current components by maximum area overlap."""
    _, comps = _component_labels(grid, mask)
    labels = np.full(grid.n, -1.0, dtype=np.float64)
    used: set[int] = set()
    previous = np.asarray(previous, dtype=np.float64)
    for comp in sorted(comps, key=lambda c: float(grid.cell_area[c].sum()), reverse=True):
        old = previous[comp].astype(int)
        old = old[old >= 0]
        assigned: int | None = None
        if old.size:
            values, counts = np.unique(old, return_counts=True)
            for idx in np.argsort(counts)[::-1]:
                candidate = int(values[idx])
                if candidate not in used:
                    assigned = candidate
                    break
        if assigned is None:
            assigned = int(next_id)
            next_id += 1
        used.add(assigned)
        labels[comp] = float(assigned)
    return labels, comps, next_id


def _continent_area_targets(grid, n_seeds: int, target_fraction: float, rng=None,
                            *, seeds=None, potential=None) -> np.ndarray:
    """Deterministic uneven but bounded continent-size targets."""
    del rng
    total_target = float(grid.cell_area.sum()) * float(target_fraction)
    if n_seeds <= 1:
        return np.asarray([max(total_target, 0.0)], dtype=np.float64)
    if potential is not None and seeds is not None:
        potential_arr = _normalize01(np.asarray(potential, dtype=np.float64))
        seed_arr = np.asarray(seeds, dtype=np.int64)
        strength = potential_arr[seed_arr]
        raw = 0.80 + 0.52 * _normalize01(strength)
    else:
        idx = np.arange(n_seeds, dtype=np.float64)
        raw = 1.0 + 0.22 * np.cos(1.73 * (idx + 1.0) + 4.1 * target_fraction)
    raw = np.clip(raw, 0.68, 1.45)
    weights = raw / raw.sum()
    targets = weights * total_target
    min_area = total_target * min(0.18, 0.62 / max(float(n_seeds), 1.0))
    max_area = total_target * (0.32 if n_seeds >= 4 else 0.38)
    targets = np.clip(targets, min_area, max_area)
    targets *= total_target / max(float(targets.sum()), 1e-12)
    return targets


def _compact_growth_priority(grid, cell: int, seed: int, owner: np.ndarray,
                             cid: int, radius_scale: float,
                             noise: np.ndarray, rng=None,
                             potential: np.ndarray | None = None) -> float:
    del rng
    dot = float(np.clip(grid.xyz[cell] @ grid.xyz[seed], -1.0, 1.0))
    angular = np.arccos(dot)
    same_neighbors = int((owner[grid.neighbors[cell]] == cid).sum())
    attachment_penalty = 0.22 * max(0, 2 - same_neighbors)
    potential_bonus = 0.0
    if potential is not None:
        potential_bonus = -0.10 * float(potential[cell])
    return (
        angular / max(radius_scale, 1e-6)
        + attachment_penalty
        + 0.12 * float(noise[cell])
        + potential_bonus
        + 0.003 * float(_stable_tie_breaker(int(cell)))
    )


def _grow_compact_continents(grid, seeds, target_fraction, rng=None,
                             *, spec=None) -> tuple[np.ndarray, np.ndarray]:
    """Grow compact continent objects to a target area.

    This is deliberately still a reduced model.  It avoids the old pure random
    flood by keeping growth centered on separated nuclei, favouring candidates
    with multiple same-continent neighbours, and adding only low-frequency
    margin irregularity.
    """
    del rng
    seeds = np.asarray(seeds, dtype=np.int64)
    n = int(seeds.size)
    is_cont = np.zeros(grid.n, dtype=bool)
    labels = np.full(grid.n, -1, dtype=np.int64)
    if n == 0 or target_fraction <= 0.0:
        return is_cont, labels

    area = grid.cell_area
    total = float(area.sum())
    potential = _proto_crust_potential(grid, spec)
    targets = _continent_area_targets(
        grid, n, target_fraction, seeds=seeds, potential=potential)
    acc = np.zeros(n, dtype=np.float64)
    noise = _deterministic_margin_noise(grid, spec)

    mean_cell_area = float(np.mean(area))
    radius_scale = np.sqrt(np.maximum(targets, mean_cell_area) / total) * np.pi
    radius_scale = np.maximum(radius_scale, 2.0 * np.median(grid.edge_lengths) / grid.radius_m)

    heap: list[tuple[float, int, int]] = []
    queued = np.zeros((n, grid.n), dtype=bool)
    for cid, seed in enumerate(seeds):
        seed = int(seed)
        heapq.heappush(heap, (0.0, cid, seed))
        queued[cid, seed] = True

    target_total = float(targets.sum())
    current_total = 0.0
    while heap and current_total < target_total:
        _, cid, c = heapq.heappop(heap)
        c = int(c)
        if labels[c] >= 0 or acc[cid] >= targets[cid]:
            continue
        foreign_neighbors = labels[grid.neighbors[c]]
        if np.any((foreign_neighbors >= 0) & (foreign_neighbors != cid)):
            continue
        labels[c] = cid
        is_cont[c] = True
        acc[cid] += float(area[c])
        current_total += float(area[c])
        for nb in grid.neighbors[c]:
            nb = int(nb)
            if labels[nb] >= 0 or queued[cid, nb] or acc[cid] >= targets[cid]:
                continue
            priority = _compact_growth_priority(
                grid, nb, int(seeds[cid]), labels, cid, float(radius_scale[cid]),
                noise, potential=potential,
            )
            heapq.heappush(heap, (priority, cid, nb))
            queued[cid, nb] = True

    # If competition between nearby seeds leaves a small area deficit, fill it
    # with the best perimeter candidates while keeping the same compact score.
    if current_total < target_total:
        candidates: list[tuple[float, int, int]] = []
        for c in np.where(is_cont)[0]:
            cid = int(labels[c])
            for nb in grid.neighbors[int(c)]:
                nb = int(nb)
                if labels[nb] >= 0:
                    continue
                priority = _compact_growth_priority(
                    grid, nb, int(seeds[cid]), labels, cid,
                    float(radius_scale[cid]), noise, potential=potential,
                )
                candidates.append((priority, cid, nb))
        heapq.heapify(candidates)
        while candidates and current_total < target_total:
            _, cid, c = heapq.heappop(candidates)
            c = int(c)
            if labels[c] >= 0:
                continue
            foreign_neighbors = labels[grid.neighbors[c]]
            if np.any((foreign_neighbors >= 0) & (foreign_neighbors != cid)):
                continue
            labels[c] = cid
            is_cont[c] = True
            acc[cid] += float(area[c])
            current_total += float(area[c])
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if labels[nb] < 0:
                    priority = _compact_growth_priority(
                        grid, nb, int(seeds[cid]), labels, cid,
                        float(radius_scale[cid]), noise, potential=potential,
                    )
                    heapq.heappush(candidates, (priority, cid, nb))
    return is_cont, labels


def _legacy_random_continent_area_targets(grid, n_seeds: int,
                                          target_fraction: float, rng) -> np.ndarray:
    total_target = float(grid.cell_area.sum()) * float(target_fraction)
    if n_seeds <= 1:
        return np.asarray([max(total_target, 0.0)], dtype=np.float64)
    raw = rng.lognormal(mean=0.0, sigma=0.38, size=n_seeds)
    raw = np.clip(raw, 0.45, 2.15)
    weights = raw / raw.sum()
    targets = weights * total_target
    min_area = total_target * 0.045
    max_area = total_target * 0.36
    targets = np.clip(targets, min_area, max_area)
    targets *= total_target / max(float(targets.sum()), 1e-12)
    return targets


def _legacy_random_compact_growth_priority(
    grid, cell: int, seed: int, owner: np.ndarray, cid: int,
    radius_scale: float, noise: np.ndarray, rng,
) -> float:
    dot = float(np.clip(grid.xyz[cell] @ grid.xyz[seed], -1.0, 1.0))
    angular = np.arccos(dot)
    same_neighbors = int((owner[grid.neighbors[cell]] == cid).sum())
    attachment_penalty = 0.18 * max(0, 2 - same_neighbors)
    return (
        angular / max(radius_scale, 1e-6)
        + attachment_penalty
        + 0.12 * float(noise[cell])
        + 0.035 * float(rng.random())
    )


def _legacy_random_grow_compact_continents(
    grid, seeds, target_fraction, rng
) -> tuple[np.ndarray, np.ndarray]:
    seeds = np.asarray(seeds, dtype=np.int64)
    n = int(seeds.size)
    is_cont = np.zeros(grid.n, dtype=bool)
    labels = np.full(grid.n, -1, dtype=np.int64)
    if n == 0 or target_fraction <= 0.0:
        return is_cont, labels

    area = grid.cell_area
    total = float(area.sum())
    targets = _legacy_random_continent_area_targets(grid, n, target_fraction, rng)
    acc = np.zeros(n, dtype=np.float64)
    noise_axes = random_unit_vectors(rng, 5)
    noise_weights = rng.uniform(0.35, 1.0, size=5)
    phase = rng.uniform(0.0, 2.0 * np.pi, size=5)
    dots = grid.xyz @ noise_axes.T
    noise = np.sum(noise_weights * np.sin(2.2 * dots + phase), axis=1)
    noise = (noise - noise.min()) / max(float(noise.max() - noise.min()), 1e-12)

    mean_cell_area = float(np.mean(area))
    radius_scale = np.sqrt(np.maximum(targets, mean_cell_area) / total) * np.pi
    radius_scale = np.maximum(radius_scale, 2.0 * np.median(grid.edge_lengths) / grid.radius_m)

    heap: list[tuple[float, int, int]] = []
    queued = np.zeros((n, grid.n), dtype=bool)
    for cid, seed in enumerate(seeds):
        seed = int(seed)
        heapq.heappush(heap, (0.0, cid, seed))
        queued[cid, seed] = True

    target_total = float(targets.sum())
    current_total = 0.0
    while heap and current_total < target_total:
        _, cid, c = heapq.heappop(heap)
        c = int(c)
        if labels[c] >= 0 or acc[cid] >= targets[cid]:
            continue
        labels[c] = cid
        is_cont[c] = True
        acc[cid] += float(area[c])
        current_total += float(area[c])
        for nb in grid.neighbors[c]:
            nb = int(nb)
            if labels[nb] >= 0 or queued[cid, nb] or acc[cid] >= targets[cid]:
                continue
            priority = _legacy_random_compact_growth_priority(
                grid, nb, int(seeds[cid]), labels, cid,
                float(radius_scale[cid]), noise, rng,
            )
            heapq.heappush(heap, (priority, cid, nb))
            queued[cid, nb] = True

    if current_total < target_total:
        candidates: list[tuple[float, int, int]] = []
        for c in np.where(is_cont)[0]:
            cid = int(labels[c])
            for nb in grid.neighbors[int(c)]:
                nb = int(nb)
                if labels[nb] >= 0:
                    continue
                priority = _legacy_random_compact_growth_priority(
                    grid, nb, int(seeds[cid]), labels, cid,
                    float(radius_scale[cid]), noise, rng,
                )
                candidates.append((priority, cid, nb))
        heapq.heapify(candidates)
        while candidates and current_total < target_total:
            _, cid, c = heapq.heappop(candidates)
            c = int(c)
            if labels[c] >= 0:
                continue
            labels[c] = cid
            is_cont[c] = True
            acc[cid] += float(area[c])
            current_total += float(area[c])
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if labels[nb] < 0:
                    priority = _legacy_random_compact_growth_priority(
                        grid, nb, int(seeds[cid]), labels, cid,
                        float(radius_scale[cid]), noise, rng,
                    )
                    heapq.heappush(candidates, (priority, cid, nb))
    return is_cont, labels


def _grow_continents(grid, seeds, target_fraction, rng) -> np.ndarray:
    """Grow compact contiguous continents to a target area."""
    return _grow_compact_continents(grid, seeds, target_fraction, rng)[0]


def _grow_within_mask(grid, allowed, seeds, target_area, rng) -> np.ndarray:
    """Grow connected regions from seeds, constrained to an existing mask."""
    allowed = np.asarray(allowed, dtype=bool)
    out = np.zeros(grid.n, dtype=bool)
    if not allowed.any() or target_area <= 0.0:
        return out
    visited = np.zeros(grid.n, dtype=bool)
    heap: list[tuple[float, int]] = []
    for s in seeds:
        s = int(s)
        if not allowed[s]:
            continue
        heapq.heappush(heap, (rng.random(), s))
        visited[s] = True
    acc = 0.0
    while heap and acc < target_area:
        _, c = heapq.heappop(heap)
        if not allowed[c]:
            continue
        out[c] = True
        acc += float(grid.cell_area[c])
        for nb in grid.neighbors[c]:
            nb = int(nb)
            if allowed[nb] and not visited[nb]:
                visited[nb] = True
                heapq.heappush(heap, (rng.random(), nb))
    return out


def _grow_compact_within_mask(grid, allowed, seeds, target_area, rng=None,
                              *, spec=None) -> np.ndarray:
    """Compact growth constrained to an allowed mask."""
    del rng
    allowed = np.asarray(allowed, dtype=bool)
    out = np.zeros(grid.n, dtype=bool)
    if not allowed.any() or target_area <= 0.0:
        return out
    seeds = [int(s) for s in np.asarray(seeds, dtype=np.int64) if allowed[int(s)]]
    if not seeds:
        potential = _proto_crust_potential(grid, spec)
        allowed_cells = np.where(allowed)[0]
        seeds = [int(allowed_cells[int(np.argmax(potential[allowed_cells]))])]

    area = grid.cell_area
    acc = 0.0
    owner = np.full(grid.n, -1, dtype=np.int64)
    seed = int(seeds[0])
    radius_scale = np.sqrt(max(float(target_area), float(np.mean(area))) / float(area.sum())) * np.pi
    radius_scale = max(radius_scale, 2.0 * np.median(grid.edge_lengths) / grid.radius_m)
    potential = _proto_crust_potential(grid, spec)
    noise = _deterministic_margin_noise(grid, spec)

    heap: list[tuple[float, int]] = []
    queued = np.zeros(grid.n, dtype=bool)
    for s in seeds:
        heapq.heappush(heap, (0.0, int(s)))
        queued[int(s)] = True
    while heap and acc < target_area:
        _, c = heapq.heappop(heap)
        if out[c] or not allowed[c]:
            continue
        out[c] = True
        owner[c] = 0
        acc += float(area[c])
        for nb in grid.neighbors[c]:
            nb = int(nb)
            if not allowed[nb] or out[nb] or queued[nb]:
                continue
            priority = _compact_growth_priority(
                grid, nb, seed, owner, 0, radius_scale, noise,
                potential=potential)
            heapq.heappush(heap, (priority, nb))
            queued[nb] = True
    return out


def _legacy_random_grow_compact_within_mask(grid, allowed, seeds, target_area, rng) -> np.ndarray:
    allowed = np.asarray(allowed, dtype=bool)
    out = np.zeros(grid.n, dtype=bool)
    if not allowed.any() or target_area <= 0.0:
        return out
    seeds = [int(s) for s in np.asarray(seeds, dtype=np.int64) if allowed[int(s)]]
    if not seeds:
        seeds = [int(rng.choice(np.where(allowed)[0]))]

    area = grid.cell_area
    acc = 0.0
    owner = np.full(grid.n, -1, dtype=np.int64)
    seed = int(seeds[0])
    radius_scale = np.sqrt(max(float(target_area), float(np.mean(area))) / float(area.sum())) * np.pi
    radius_scale = max(radius_scale, 2.0 * np.median(grid.edge_lengths) / grid.radius_m)
    noise_axes = random_unit_vectors(rng, 3)
    noise_weights = rng.uniform(0.4, 1.0, size=3)
    phase = rng.uniform(0.0, 2.0 * np.pi, size=3)
    dots = grid.xyz @ noise_axes.T
    noise = np.sum(noise_weights * np.sin(2.0 * dots + phase), axis=1)
    noise = (noise - noise.min()) / max(float(noise.max() - noise.min()), 1e-12)

    heap: list[tuple[float, int]] = []
    queued = np.zeros(grid.n, dtype=bool)
    for s in seeds:
        heapq.heappush(heap, (0.0, int(s)))
        queued[int(s)] = True
    while heap and acc < target_area:
        _, c = heapq.heappop(heap)
        if out[c] or not allowed[c]:
            continue
        out[c] = True
        owner[c] = 0
        acc += float(area[c])
        for nb in grid.neighbors[c]:
            nb = int(nb)
            if not allowed[nb] or out[nb] or queued[nb]:
                continue
            priority = _legacy_random_compact_growth_priority(
                grid, nb, seed, owner, 0, radius_scale, noise, rng)
            heapq.heappush(heap, (priority, nb))
            queued[nb] = True
    return out


class TectonicsModule(Module):
    name = "tectonics"
    produces = ["tectonics.plate_id", "tectonics.plate_velocity", "crust.type",
                "crust.thickness_m", "crust.age_myr", "crust.domain",
                "tectonics.continent_id", "tectonics.terrane_id",
                "tectonics.cratonic_province_code",
                "tectonics.internal_geographic_block_id",
                "tectonics.internal_geographic_block_code",
                "tectonics.basement_age_floor_myr",
                "tectonics.basement_stability_floor",
                "tectonics.basement_thickness_floor_m",
                "tectonics.stable_cratonic_lithosphere_support",
                "tectonics.plate_rank", "tectonics.protected_plate_id",
                "tectonics.boundaries",
                "tectonics.boundary_province_kind",
                "tectonics.deformation_intensity", "tectonics.deformation_style",
                "tectonics.deforming_networks",
                "tectonics.plates", "tectonics.boundary_objects",
                "tectonics.boundary_polylines",
                "tectonics.boundary_provinces",
                "tectonics.wilson_cycles", "tectonics.ocean_gateways",
                "tectonics.ocean_basins", "tectonics.rift_systems",
                "tectonics.breakup_seaways",
                "tectonics.passive_margins", "tectonics.spreading_centers",
                "tectonics.closing_margins", "tectonics.sutures",
                "tectonics.cratons", "tectonics.shields",
                "tectonics.platforms", "tectonics.interior_basins",
                "tectonics.internal_geographic_blocks",
                "tectonics.platform_subsidence",
                "archive.wilson_cycle_phase",
                "tectonics.continents", "tectonics.terranes",
                "tectonics.plate_topologies", "tectonics.microplates",
                "tectonics.plate_accretion_events",
                "tectonics.plumes", "tectonics.lips",
                "tectonics.orogens", "tectonics.volcanoes"]
    fidelity = "voronoi_euler"
    interval_myr = 20.0

    REORG_INTERVAL = 300.0
    MIN_PLATE_AREA_FRAC = 0.008
    MICRO_PLATE_AREA_FRAC = 0.0005
    MAJOR_PLATE_AREA_FRAC = 0.05
    MAX_PLATE_AREA_FRAC = 0.24
    R2_PARAMETERS = {
        "slab_pull_weight": 1.00,
        "ridge_push_weight": 0.62,
        "collision_resistance_weight": 1.20,
        "basal_drag_weight": 0.26,
        "transform_friction_weight": 0.22,
        "torque_rate_scale_rad_per_myr": 0.040,
        "min_rotation_rate_rad_per_myr": 0.0030,
        "max_rotation_rate_rad_per_myr": 0.0140,
        "fallback_rate_decay": 0.985,
        "motion_memory_weight": 0.59,
        "rate_memory_weight": 0.55,
    }

    def init_state(self, world, rng_key) -> None:
        grid = world.grid
        spec = world.spec
        rng = rng_key.child("init").generator()

        # Continental crust starts as compact proto-continents with broad
        # cratonic kernels, not as unconstrained random graph tendrils.
        n_nuclei = _continent_nucleus_count(spec)
        deterministic_initial_cargo = (
            world.g("tectonics.enable_deterministic_initial_cargo", 0.0) > 0.0
        )
        if deterministic_initial_cargo and 0.08 < float(spec.target_land_fraction) < 0.62:
            n_nuclei = max(n_nuclei, min(4, int(spec.n_plates)))
        proto_potential = _proto_crust_potential(grid, spec)
        margin_noise = _deterministic_margin_noise(grid, spec)
        if deterministic_initial_cargo:
            seeds = _separated_seed_cells(
                grid, n_nuclei, rng, potential=proto_potential)
            is_cont, continent_label = _grow_compact_continents(
                grid, seeds, spec.target_land_fraction, rng, spec=spec)
        else:
            seeds = _legacy_random_separated_seed_cells(grid, n_nuclei, rng)
            is_cont, continent_label = _legacy_random_grow_compact_continents(
                grid, seeds, spec.target_land_fraction, rng)
        ctype = np.where(is_cont, CONT, OCEAN)

        # Plates: Voronoi over seed cells, with continental nuclei starting
        # inside plate interiors rather than on arbitrary plate seams.
        if deterministic_initial_cargo:
            pseeds = _plate_seed_cells(grid, spec.n_plates, seeds, rng)
        else:
            pseeds = _legacy_random_plate_seed_cells(grid, spec.n_plates, seeds, rng)
        ptree = cKDTree(grid.xyz[pseeds])
        _, plate_id = ptree.query(grid.xyz)
        plate_id = plate_id.astype(np.int64)
        deterministic_initial_plate_motions = (
            world.g("tectonics.enable_deterministic_initial_plate_motions", 0.0) > 0.0
        )
        if not deterministic_initial_plate_motions:
            # Legacy production path, kept until P53+P54 preserve mature
            # cratonic province coverage through full-history audits.
            poles = random_unit_vectors(rng, spec.n_plates)
            base = 0.006 * rng.uniform(0.5, 2.0, size=spec.n_plates)
            plates = [{
                "pole": poles[i].tolist(),
                "rate": float(base[i]),
                "id": int(i),
                "motion_source": "legacy_random_initial_euler",
            } for i in range(spec.n_plates)]
            initial_motion = {
                "schema": "aevum.tectonics.initial_plate_motion.legacy.v1",
                "boundaries": {},
                "diagnostics": {},
                "telemetry": {
                    "motion_active_plates": float(spec.n_plates),
                    "motion_torque_driven_fraction": 0.0,
                    "motion_fallback_fraction": 0.0,
                    "motion_mean_rate_rad_per_myr": float(np.mean(base)),
                    "motion_rate_cv": float(np.std(base) / max(float(np.mean(base)), 1e-12)),
                    "motion_net_torque_ratio": 0.0,
                    "boundary_proxy_ridge_cells": 0.0,
                    "boundary_proxy_trench_cells": 0.0,
                    "boundary_proxy_collision_cells": 0.0,
                    "boundary_proxy_transform_cells": 0.0,
                },
            }

        interiority = _graph_width_steps(grid, is_cont)
        if deterministic_initial_cargo:
            thickness = np.where(is_cont, CONT_THICK, OCEAN_THICK).astype(np.float64)
            thickness += is_cont * (
                420.0 * (proto_potential - 0.5)
                + 520.0 * (interiority / max(float(interiority.max()), 1.0))
                + 260.0 * (margin_noise - 0.5)
            )
            ocean_age_pattern = _normalize01(
                0.55 * (1.0 - proto_potential)
                + 0.45 * np.abs(np.sin(np.radians(grid.lon * 1.3 + grid.lat * 0.7)))
            )
            age = np.where(
                is_cont,
                35.0 + 85.0 * proto_potential + 18.0 * np.clip(interiority, 0.0, 4.0),
                4.0 + 26.0 * ocean_age_pattern,
            ).astype(np.float64)
        else:
            thickness = np.where(is_cont, CONT_THICK, OCEAN_THICK).astype(np.float64)
            thickness += rng.normal(0, 650, size=grid.n) * is_cont
            age = np.where(is_cont, rng.uniform(0, 110, grid.n),
                           rng.uniform(0, 30, grid.n)).astype(np.float64)

        world.set_field("crust.type", ctype)
        world.set_field("crust.thickness_m", thickness)
        world.set_field("crust.age_myr", age)
        origin = np.where(is_cont, ORIGIN_PRIMORDIAL, ORIGIN_RIDGE).astype(np.float64)
        proto_craton = np.zeros(grid.n, dtype=bool)
        for cid, seed in enumerate(seeds):
            component = is_cont & (continent_label == cid)
            if not component.any():
                continue
            comp_area = float(grid.cell_area[component].sum())
            if deterministic_initial_cargo:
                seed_strength = float(proto_potential[int(seed)])
                core_fraction = float(np.clip(0.50 + 0.16 * seed_strength, 0.50, 0.64))
            else:
                core_fraction = float(rng.uniform(0.50, 0.64))
            core_target = comp_area * core_fraction
            if deterministic_initial_cargo:
                core = _grow_compact_within_mask(
                    grid, component, [int(seed)], core_target, rng, spec=spec)
            else:
                core = _legacy_random_grow_compact_within_mask(
                    grid, component, [int(seed)], core_target, rng)
            proto_craton |= core
        origin[proto_craton] = ORIGIN_CRATON
        if deterministic_initial_cargo:
            core_strength = np.clip(
                0.60 * proto_potential[proto_craton]
                + 0.40 * _normalize01(interiority)[proto_craton],
                0.0,
                1.0,
            )
            thickness[proto_craton] = np.maximum(
                thickness[proto_craton] + 3600.0 + 4700.0 * core_strength,
                CONT_THICK + 3500.0,
            )
            age[proto_craton] = np.maximum(
                age[proto_craton],
                120.0 + 130.0 * core_strength,
            )
        else:
            thickness[proto_craton] = np.maximum(
                thickness[proto_craton] + rng.uniform(
                    3500.0, 8500.0, int(proto_craton.sum())),
                CONT_THICK + 3500.0,
            )
            age[proto_craton] = np.maximum(
                age[proto_craton],
                rng.uniform(90.0, 240.0, int(proto_craton.sum())),
            )
        reworked = np.where(is_cont, -1.0, 0.0).astype(np.float64)
        stability = self._crust_stability(ctype, age, origin)
        stability[proto_craton] = np.maximum(stability[proto_craton], 0.82)
        stability = self._shape_guard_cratonic_stability(
            grid, ctype, origin, stability)
        domain = self._crust_domain_field(grid, ctype, origin, stability)
        continent_id = continent_label.astype(np.float64)
        continent_id[~is_cont] = -1.0
        terrane_id = np.full(grid.n, -1.0, dtype=np.float64)
        initial_width = _graph_width_steps(grid, is_cont)
        province_code = np.full(grid.n, PROVINCE_NONE, dtype=np.float64)
        province_code[
            is_cont
            & (proto_craton | (origin == ORIGIN_CRATON))
            & (initial_width >= 2.0)
        ] = PROVINCE_SHIELD
        initial_platform = (
            is_cont
            & (province_code == PROVINCE_NONE)
            & (origin == ORIGIN_PRIMORDIAL)
            & (initial_width >= 3.0)
            & (stability >= 0.30)
        )
        province_code[initial_platform] = PROVINCE_PLATFORM
        internal_block_code = self._initial_internal_geographic_block_code(
            grid,
            is_cont,
            continent_label,
            proto_craton,
            initial_width,
            proto_potential,
            margin_noise,
            thickness,
            stability,
        )
        province_code[
            internal_block_code == INTERNAL_BLOCK_CRATON_CORE
        ] = PROVINCE_SHIELD
        province_code[
            internal_block_code == INTERNAL_BLOCK_INTRACRATONIC_BASIN
        ] = PROVINCE_INTERIOR_BASIN
        province_code[
            np.isin(internal_block_code, [
                INTERNAL_BLOCK_STABLE_PLATFORM,
                INTERNAL_BLOCK_MOBILE_BELT,
            ])
            & (province_code == PROVINCE_NONE)
        ] = PROVINCE_PLATFORM
        basement_age_floor = np.where(is_cont, age, 0.0).astype(np.float64)
        basement_stability_floor = np.where(is_cont, stability, 0.0).astype(np.float64)
        basement_thickness_floor = np.where(is_cont, thickness, 0.0).astype(np.float64)
        basement_stability_floor[proto_craton] = np.maximum(
            basement_stability_floor[proto_craton], 0.82)
        basement_thickness_floor[proto_craton] = np.maximum(
            basement_thickness_floor[proto_craton], CONT_THICK + 3500.0)
        basement_stability_floor[initial_platform] = np.maximum(
            basement_stability_floor[initial_platform], 0.32)
        basement_thickness_floor[initial_platform] = np.maximum(
            basement_thickness_floor[initial_platform], CONT_THICK + 1200.0)
        if deterministic_initial_plate_motions:
            plates, initial_motion = self._initial_plate_motions_from_torque(
                world, grid, plate_id, ctype, age, thickness, proto_potential)
        internal_block_id, internal_block_objects, next_internal_block_id = (
            self._internal_geographic_block_ids_and_objects(
                world,
                grid,
                ctype,
                internal_block_code,
                np.zeros(grid.n, dtype=np.float64),
                1,
                age,
                thickness,
                stability,
                continent_id,
                0.0,
            )
        )
        initial_plate_topologies = self._plate_topology_objects(
            grid, plate_id, plates, ctype, continent_id, terrane_id, {}, 0.0)
        initial_plate_rank, initial_protected_plate_id, initial_microplates = (
            self._plate_rank_fields_and_objects(
                grid,
                plate_id,
                initial_plate_topologies,
                ctype,
                origin,
                age,
                {},
                0.0,
                protected_records=[],
                enable_protection=False,
            )
        )

        world.set_field("tectonics.plate_id", plate_id.astype(np.float64))
        world.set_field("tectonics.plate_rank", initial_plate_rank)
        world.set_field("tectonics.protected_plate_id", initial_protected_plate_id)
        world.set_field("tectonics.boundary_province_kind", np.full(grid.n, -1.0))
        world.set_field("tectonics.orogeny_age_myr", np.full(grid.n, -1.0))
        world.set_field("tectonics.volcanism_age_myr", np.full(grid.n, -1.0))
        world.set_field("crust.origin", origin)
        world.set_field("crust.reworked_age_myr", reworked)
        world.set_field("crust.stability", stability)
        world.set_field("crust.domain", domain)
        world.set_field("tectonics.platform_subsidence", np.zeros(grid.n, dtype=np.float64))
        world.set_field("tectonics.continent_id", continent_id)
        world.set_field("tectonics.terrane_id", terrane_id)
        world.set_field("tectonics.cratonic_province_code", province_code)
        world.set_field("tectonics.internal_geographic_block_id", internal_block_id)
        world.set_field("tectonics.internal_geographic_block_code", internal_block_code)
        world.set_field("tectonics.basement_age_floor_myr", basement_age_floor)
        world.set_field("tectonics.basement_stability_floor", basement_stability_floor)
        world.set_field("tectonics.basement_thickness_floor_m", basement_thickness_floor)
        initial_cratonic_support = np.zeros(grid.n, dtype=np.float64)
        initial_cratonic_support[proto_craton] = 1.0
        initial_cratonic_support[initial_platform] = 0.74
        world.set_field(
            "tectonics.stable_cratonic_lithosphere_support",
            initial_cratonic_support,
        )
        world.objects["tectonics.plates"] = plates
        world.objects["tectonics.continents"] = self._continent_objects(
            grid, continent_id, ctype, age, thickness, origin, stability, 0.0)
        world.objects["tectonics.plate_topologies"] = initial_plate_topologies
        world.objects["tectonics.microplates"] = initial_microplates
        world.objects["tectonics.plate_accretion_events"] = []
        world.objects["tectonics.boundary_provinces"] = []
        world.objects["tectonics.initial_plate_motion"] = initial_motion
        world.objects["tectonics.internal_geographic_blocks"] = internal_block_objects
        world.objects["tectonics.terranes"] = []
        world.objects["tectonics.continent_lifecycle_events"] = []
        world.objects["tectonics.breakup_seaways"] = []
        cont_frac = grid.cell_area[ctype == CONT].sum() / grid.cell_area.sum()
        cont_area = max(float(grid.cell_area[is_cont].sum()), 1.0e-12)
        world.set_g("tectonics.last_p52_initial_nuclei_count", float(n_nuclei))
        world.set_g(
            "tectonics.last_p52_initial_continent_width_p50_steps",
            float(np.percentile(initial_width[is_cont], 50.0)) if is_cont.any() else 0.0,
        )
        world.set_g(
            "tectonics.last_p52_initial_continent_width_p90_steps",
            float(np.percentile(initial_width[is_cont], 90.0)) if is_cont.any() else 0.0,
        )
        world.set_g(
            "tectonics.last_p52_initial_narrow_fraction_le2",
            float(grid.cell_area[is_cont & (initial_width <= 2.0)].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p52_initial_craton_fraction",
            float(grid.cell_area[proto_craton].sum() / cont_area),
        )
        for key, value in initial_motion.get("telemetry", {}).items():
            world.set_g(f"tectonics.last_p53_initial_{key}", float(value))
        world.set_g("crust.cont_fraction_init", float(cont_frac))
        world.set_g("tectonics.last_reorg", 0.0)
        world.set_g("tectonics.next_continent_id", float(max(n_nuclei, 1)))
        world.set_g("tectonics.next_terrane_id", 0.0)
        world.set_g("tectonics.next_internal_geographic_block_id",
                    float(next_internal_block_id))
        for key in (
            "tectonics.cumulative_p35_continental_debt_candidate_cells",
            "tectonics.cumulative_p35_continental_debt_restored_cells",
            "tectonics.cumulative_p35_continental_debt_area_fraction",
            "tectonics.cumulative_p37_shape_guard_input_cells",
            "tectonics.cumulative_p37_shape_guard_accepted_cells",
            "tectonics.cumulative_p37_shape_guard_rejected_cells",
            "tectonics.cumulative_p38_lineage_guard_input_cells",
            "tectonics.cumulative_p38_lineage_guard_accepted_cells",
            "tectonics.cumulative_p38_lineage_guard_rejected_cells",
            "tectonics.cumulative_p50_planform_recycled_cells",
            "tectonics.cumulative_p50_planform_filled_cells",
            "tectonics.cumulative_p50_planform_area_fraction",
        ):
            world.set_g(key, 0.0)

        # neighbour radius (chord) used for rasterisation.
        nn = np.median(grid.edge_lengths) / grid.radius_m
        world.set_g("tectonics.raster_radius", float(1.5 * nn))

    # ------------------------------------------------------------------
    def _plate_rotations(self, plates, vigor, dt, regime_code):
        regime_factor = {0.0: 0.02, 1.0: 0.4, 2.0: 0.7, 3.0: 1.0}.get(regime_code, 1.0)
        Rs = np.zeros((len(plates), 3, 3))
        for i, p in enumerate(plates):
            angle = p["rate"] * max(vigor, 0.05) * regime_factor * dt
            Rs[i] = rotation_matrix(np.asarray(p["pole"]), angle)
        return Rs

    def step(self, world, t, dt, forcing, rng_key) -> StepResult:
        grid = world.grid
        spec = world.spec
        xyz = grid.xyz
        plate = world.field("tectonics.plate_id").astype(int)
        ctype = world.field("crust.type").copy()
        thick = world.field("crust.thickness_m").copy()
        age = world.field("crust.age_myr").copy()
        origin = world.get_field("crust.origin", ORIGIN_RIDGE).copy()
        reworked = world.get_field("crust.reworked_age_myr", -1.0).copy()
        stability = world.get_field("crust.stability", 0.0).copy()
        province_code = world.get_field("tectonics.cratonic_province_code", 0.0).copy()
        internal_block_id = world.get_field(
            "tectonics.internal_geographic_block_id", 0.0).copy()
        internal_block_code = world.get_field(
            "tectonics.internal_geographic_block_code", 0.0).copy()
        if internal_block_id.shape != (grid.n,):
            internal_block_id = np.zeros(grid.n, dtype=np.float64)
        if internal_block_code.shape != (grid.n,):
            internal_block_code = np.zeros(grid.n, dtype=np.float64)
        if "tectonics.basement_age_floor_myr" in world.fields:
            basement_age = world.field("tectonics.basement_age_floor_myr").copy()
        else:
            basement_age = np.where(ctype == CONT, age, 0.0).astype(np.float64)
        if "tectonics.basement_stability_floor" in world.fields:
            basement_stability = world.field("tectonics.basement_stability_floor").copy()
        else:
            basement_stability = np.where(ctype == CONT, stability, 0.0).astype(np.float64)
        if "tectonics.basement_thickness_floor_m" in world.fields:
            basement_thick = world.field("tectonics.basement_thickness_floor_m").copy()
        else:
            basement_thick = np.where(ctype == CONT, thick, 0.0).astype(np.float64)
        plates = world.objects["tectonics.plates"]
        vigor = world.g("interior.tectonic_vigor", 1.0)
        regime_code = world.g("tectonics.regime_code", 3.0)

        events: list[Event] = []

        # Periodic plate reorganization should be a local tectonic adjustment,
        # not a global redraw of every plate.  The earlier Voronoi reset made
        # deep-time frames jump discontinuously.  Here we preserve labels and
        # geometry except for small merge/split repairs and a mild Euler-pole
        # refresh, which is closer to ridge jumps, plate capture and local
        # plate breakup.
        reorg_detail = None
        plate_before_reorg = plate.copy()
        if (t - world.g("tectonics.last_reorg") >= self.REORG_INTERVAL
                and regime_code >= 1.0):
            plate, plates, reorg_detail = self._local_reorganize_plates(
                world, grid, plate, plates, t
            )
            world.set_g("tectonics.last_reorg", t)
            events.append(Event("plate_reorganization", t, self.name,
                                magnitude=float(len(np.unique(plate))),
                                params=reorg_detail))
        self._record_p58_reorg_basement_ledger(
            world, grid, plate_before_reorg, plate, ctype, basement_age,
            basement_stability, basement_thick, reorg_detail, t)

        Rs = self._plate_rotations(plates, vigor, dt, regime_code)
        adv = np.einsum("nij,nj->ni", Rs[plate], xyz)        # moved crust parcels

        r = world.g("tectonics.raster_radius", 0.04)
        K = min(10, grid.n)
        tree = cKDTree(adv)
        dist, idx = tree.query(xyz, k=K)

        primary = idx[:, 0]
        primary, basement_primary = self._apply_p59_old_basement_raster_coverage(
            world, grid, primary, idx, dist, plate, ctype, origin, stability,
            basement_age, basement_stability, basement_thick, t)
        within = dist < r
        plate0 = plate[primary]
        type0 = ctype[primary]
        plate_nbr = plate[idx]
        type_nbr = ctype[idx]
        diff_plate = (plate_nbr != plate0[:, None]) & within
        boundary = diff_plate.any(axis=1)
        diff_rank = np.argmax(diff_plate, axis=1)
        diff_src = idx[np.arange(grid.n), diff_rank]
        before = np.linalg.norm(xyz[diff_src] - xyz[primary], axis=1)
        after = np.linalg.norm(adv[diff_src] - adv[primary], axis=1)
        rel_delta = after - before
        # Reserve a neutral band for transform-like shear.  Treating every
        # non-divergent contact as convergence made trenches too dense and
        # removed transform boundaries entirely.
        motion_eps = np.maximum(before * 0.0015, 1e-6)
        separating = boundary & (rel_delta > motion_eps)
        converging = boundary & (rel_delta < -motion_eps)
        seafloor_new = dist[:, 0] > r
        cont_within = ((type_nbr == CONT) & within).any(axis=1)
        cont_diffplate = ((type_nbr == CONT) & diff_plate).any(axis=1)
        collision = converging & cont_diffplate & (type0 == CONT)
        subduction = converging & ~collision
        ridge_seed = seafloor_new | (separating & (type0 == OCEAN))
        ridge = self._coherent_ridge_mask(grid, ridge_seed, type0 == OCEAN)
        protected_continent = ((type0 == CONT)
                               & ((origin[primary] == ORIGIN_CRATON)
                                  | (stability[primary] > 0.75)))
        ridge &= ~protected_continent
        trench = subduction
        suture = collision
        active_margin = trench & cont_within
        passive_margin = separating & cont_within & ~ridge
        transform = boundary & ~(separating | converging)

        # ---- new crust state -------------------------------------------------
        ctype_new = type0.astype(np.float64).copy()
        thick_new = thick[primary].copy()
        age_new = np.minimum(age[primary] + dt, t)
        origin_new = origin[primary].copy()
        reworked_new = reworked[primary].copy()
        stability_new = stability[primary].copy()
        province_code_new = province_code[primary].copy()
        internal_block_id_new = internal_block_id[primary].copy()
        internal_block_code_new = internal_block_code[primary].copy()
        basement_age_new = np.minimum(basement_age[basement_primary] + dt, t)
        basement_stability_new = basement_stability[basement_primary].copy()
        basement_thick_new = basement_thick[basement_primary].copy()
        inherited_craton = (origin_new == ORIGIN_CRATON) | (stability_new > 0.75)

        plate_new = self._repair_plate_speckle(grid, plate0.copy())
        plate_new, component_merges = self._merge_detached_plate_components(grid, plate_new)
        plate_contact = self._plate_contact_mask(grid, plate_new)

        arc_process = subduction & cont_within
        if world.g("tectonics.enable_process_provenance_localization", 0.0) > 0.0:
            collision_provenance = self._process_provenance_mask(
                grid, collision, plate_contact, spacing=2, shoulder_passes=1)
            arc_provenance = self._process_provenance_mask(
                grid, arc_process, plate_contact, spacing=2, shoulder_passes=1)
        else:
            collision_provenance = collision
            arc_provenance = arc_process
        world.set_g("tectonics.last_p45_collision_process_cells", float(np.count_nonzero(collision)))
        world.set_g(
            "tectonics.last_p45_collision_provenance_cells",
            float(np.count_nonzero(collision_provenance)),
        )
        world.set_g("tectonics.last_p45_arc_process_cells", float(np.count_nonzero(arc_process)))
        world.set_g(
            "tectonics.last_p45_arc_provenance_cells",
            float(np.count_nonzero(arc_provenance)),
        )

        ctype_new[ridge] = OCEAN
        thick_new[ridge] = OCEAN_THICK
        age_new[ridge] = 0.0
        origin_new[ridge] = ORIGIN_RIDGE
        reworked_new[ridge] = t
        stability_new[ridge] = 0.0
        province_code_new[ridge] = PROVINCE_NONE
        internal_block_id_new[ridge] = 0.0
        internal_block_code_new[ridge] = INTERNAL_BLOCK_NONE
        basement_age_new[ridge] = 0.0
        basement_stability_new[ridge] = 0.0
        basement_thick_new[ridge] = 0.0

        cont_width0 = _graph_width_steps(grid, type0 == CONT)
        broad_parent_continent = (type0 == CONT) & (cont_width0 >= 3.0)
        parent_support = _same_neighbor_count(grid, broad_parent_continent)
        attached_arc = subduction & cont_within & (parent_support >= 1)
        cont_wins = collision | (subduction & (type0 == CONT)) | attached_arc
        newly_cont = cont_wins & (type0 == OCEAN)
        oceanic_arc = subduction & ~cont_wins & (type0 == OCEAN)
        ctype_new[cont_wins] = CONT
        thick_new[newly_cont] = np.maximum(thick_new[newly_cont], CONT_THICK)
        origin_new[newly_cont] = ORIGIN_ARC
        internal_block_id_new[newly_cont] = 0.0
        internal_block_code_new[newly_cont] = INTERNAL_BLOCK_ACCRETED_TERRANE
        provenance_rework = newly_cont | collision_provenance | arc_provenance
        reworked_new[provenance_rework & ~inherited_craton] = t
        thick_new[oceanic_arc] = np.maximum(thick_new[oceanic_arc], OCEAN_THICK + 3500.0)
        origin_new[oceanic_arc] = ORIGIN_ARC
        reworked_new[oceanic_arc] = t
        province_code_new[oceanic_arc] = PROVINCE_NONE
        internal_block_id_new[oceanic_arc] = 0.0
        internal_block_code_new[oceanic_arc] = INTERNAL_BLOCK_NONE
        basement_age_new[oceanic_arc] = 0.0
        basement_stability_new[oceanic_arc] = 0.0
        basement_thick_new[oceanic_arc] = 0.0
        p58_post_boundary_ctype = ctype_new.copy()
        p58_post_boundary_basement_age = basement_age_new.copy()
        p58_post_boundary_basement_stability = basement_stability_new.copy()
        p58_post_boundary_basement_thick = basement_thick_new.copy()

        thick_new[collision] = np.minimum(thick_new[collision] + COLLIDE_UPLIFT * dt,
                                          MAX_CONT_THICK)
        origin_new[collision_provenance & ~inherited_craton] = ORIGIN_SUTURE
        mobile_rework = collision_provenance & (ctype_new == CONT) & ~inherited_craton
        internal_block_id_new[mobile_rework] = 0.0
        internal_block_code_new[mobile_rework] = INTERNAL_BLOCK_MOBILE_BELT
        arc = arc_process
        thick_new[arc] = np.minimum(thick_new[arc] + ARC_UPLIFT * dt, 60000.0)
        origin_new[arc_provenance & ~inherited_craton] = ORIGIN_ARC
        reworked_new[arc_provenance & ~inherited_craton] = t
        accreted_rework = arc_provenance & (ctype_new == CONT) & ~inherited_craton
        internal_block_id_new[accreted_rework] = 0.0
        internal_block_code_new[accreted_rework] = INTERNAL_BLOCK_ACCRETED_TERRANE

        # ---- orogeny / volcanism age tracking --------------------------------
        orog_age = world.field("tectonics.orogeny_age_myr").copy()
        volc_age = world.field("tectonics.volcanism_age_myr").copy()
        orog_age[collision_provenance & ~inherited_craton] = t
        volc_age[(arc_provenance | oceanic_arc) & ~inherited_craton] = t

        # ---- conserve continental crust area ---------------------------------
        # This must update the full crust state, not only the type bit.  Earlier
        # versions could flip thick continental/orogenic cells back to OCEAN
        # without resetting thickness, producing 3 km "oceanic crust" highlands.
        rift_potential = world.get_field("tectonics.rift_potential", 0.0)
        rifted_margin_sources = self._rifted_margin_candidates(
            grid, ctype_new, stability_new, rift_potential,
            margin_sources=passive_margin | separating,
        )
        rifted_block = (
            (passive_margin | separating | rifted_margin_sources)
            & (ctype_new == CONT)
            & (internal_block_code_new != INTERNAL_BLOCK_CRATON_CORE)
            & ~inherited_craton
        )
        internal_block_id_new[rifted_block] = 0.0
        internal_block_code_new[rifted_block] = INTERNAL_BLOCK_RIFTED_MARGIN
        world.set_g("tectonics.last_preconserve_craton_promoted_cells", 0.0)
        world.set_g("tectonics.last_preconserve_craton_anchor_cells", 0.0)
        ctype_new, thick_new, age_new, origin_new, reworked_new, stability_new, orog_age, volc_age, frac = (
            self._conserve_continental(world, grid, ctype_new, thick_new, age_new,
                                       origin_new, reworked_new, stability_new,
                                       orog_age, volc_age, t, dt, rng_key,
                                       accretion_sources=subduction | collision,
                                       erosion_sources=ridge | separating | rifted_margin_sources,
                                       basement_age=basement_age_new,
                                       basement_stability=basement_stability_new,
                                       basement_thick=basement_thick_new)
        )
        ridge_axis = self._coherent_ridge_mask(
            grid, ridge & (plate_contact | (age_new <= dt)), ctype_new == OCEAN
        )
        if not ridge_axis.any() and (ridge & (ctype_new == OCEAN)).any():
            ridge_axis = ridge & (ctype_new == OCEAN)
        trench_axis = trench & plate_contact
        if not trench_axis.any() and trench.any():
            trench_axis = trench

        age_new, reworked_new = self._impose_seafloor_age_from_ridges(
            grid, ctype_new, age_new, reworked_new, ridge_axis, regime_code, vigor, t
        )
        age_new, reworked_new = self._recycle_old_oceanic_crust(
            grid, ctype_new, age_new, reworked_new, trench_axis, regime_code, vigor,
            t, rng_key.child("ocean_recycle").generator()
        )
        oceanic = ctype_new == OCEAN
        too_thick_ocean = oceanic & (thick_new > OCEAN_THICK + 3000.0)
        thick_new[too_thick_ocean] = OCEAN_THICK + 1500.0
        orog_age[too_thick_ocean] = -1.0
        origin_new[too_thick_ocean] = ORIGIN_RIDGE
        internal_block_id_new[too_thick_ocean] = 0.0
        internal_block_code_new[too_thick_ocean] = INTERNAL_BLOCK_NONE
        p58_pre_repair_ctype = ctype_new.copy()
        p58_pre_repair_basement_age = basement_age_new.copy()
        p58_pre_repair_basement_stability = basement_stability_new.copy()
        p58_pre_repair_basement_thick = basement_thick_new.copy()
        basement_age_new, basement_stability_new, basement_thick_new = (
            self._repair_basement_cargo_after_conservation(
                world, grid, ctype_new, age_new, thick_new, stability_new,
                basement_age_new, basement_stability_new, basement_thick_new, t)
        )
        self._record_p58_raster_basement_ledger(
            world, grid, plate, ctype, basement_age, basement_stability,
            basement_thick, basement_primary, ridge, oceanic_arc,
            p58_post_boundary_ctype, p58_post_boundary_basement_age,
            p58_post_boundary_basement_stability, p58_post_boundary_basement_thick,
            p58_pre_repair_ctype, p58_pre_repair_basement_age,
            p58_pre_repair_basement_stability, p58_pre_repair_basement_thick,
            ctype_new, basement_age_new, basement_stability_new, basement_thick_new,
            t)
        stability_new = self._shape_guard_cratonic_stability(
            grid, ctype_new, origin_new,
            self._crust_stability(ctype_new, age_new, origin_new))
        origin_new, thick_new = self._mature_continental_provinces(
            grid, ctype_new, thick_new, age_new, origin_new, reworked_new,
            stability_new, orog_age, volc_age, t, dt)
        stability_new = self._shape_guard_cratonic_stability(
            grid, ctype_new, origin_new,
            self._crust_stability(ctype_new, age_new, origin_new))
        origin_new, thick_new, age_new, stability_new = (
            self._consolidate_parented_margin_slivers(
                world, grid, ctype_new, thick_new, age_new, origin_new,
                reworked_new, stability_new, orog_age, volc_age, t, dt,
                active_boundary_mask=collision | subduction | separating | ridge)
        )
        origin_new, thick_new, age_new, stability_new = (
            self._mature_parent_supported_collage_belts(
                world, grid, ctype_new, thick_new, age_new, origin_new,
                reworked_new, stability_new, orog_age, volc_age, t, dt,
                active_boundary_mask=collision | subduction | separating | ridge)
        )
        stability_new = self._shape_guard_cratonic_stability(
            grid, ctype_new, origin_new,
            np.maximum(
                stability_new, self._crust_stability(ctype_new, age_new, origin_new))
        )
        cont_width_new = _graph_width_steps(grid, ctype_new == CONT)
        craton_threshold = min(2200.0, max(650.0, 0.45 * t))
        craton_quiet_myr = max(260.0, 0.22 * t)
        craton_sources = np.isin(origin_new, [ORIGIN_PRIMORDIAL, ORIGIN_SUTURE, ORIGIN_CRATON])
        craton_quiet = (reworked_new < 0.0) | (reworked_new <= t - craton_quiet_myr)
        craton_interior = cont_width_new >= 3.0
        cratonized = ((ctype_new == CONT) & craton_sources & craton_quiet
                      & craton_interior & (age_new > craton_threshold)
                      & (origin_new != ORIGIN_CRATON))
        origin_new[cratonized] = ORIGIN_CRATON
        stability_new[cratonized] = np.maximum(stability_new[cratonized], 0.82)
        stability_new = self._shape_guard_cratonic_stability(
            grid, ctype_new, origin_new, stability_new)

        plume_potential = world.get_field("tectonics.plume_potential", 0.0)
        plume_mask, lip_mask, plume_objects, lip_objects, plume_events = self._plume_activity(
            grid, ctype_new, thick_new, origin_new, reworked_new, volc_age,
            plume_potential, t, dt, regime_code, vigor
        )
        events.extend(plume_events)
        stability_new = self._shape_guard_cratonic_stability(
            grid, ctype_new, origin_new,
            self._crust_stability(ctype_new, age_new, origin_new))
        stability_new[origin_new == ORIGIN_CRATON] = np.maximum(
            stability_new[origin_new == ORIGIN_CRATON], 0.82
        )
        stability_new = self._shape_guard_cratonic_stability(
            grid, ctype_new, origin_new, stability_new)
        origin_new, age_new, thick_new, stability_new = self._maintain_cratonic_kernels(
            grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
            stability_new, t, dt, world=world)
        stability_new = self._shape_guard_cratonic_stability(
            grid, ctype_new, origin_new, stability_new)
        origin_new, age_new, thick_new, reworked_new, stability_new = (
            self._preserve_inherited_province_physical_cargo(
                world, grid, ctype_new, origin_new, age_new, thick_new,
                reworked_new, stability_new, province_code_new, t,
                basement_age_floor=basement_age_new,
                basement_stability_floor=basement_stability_new,
                basement_thickness_floor=basement_thick_new,
                active_boundary_mask=collision | subduction | separating | ridge)
        )
        stability_new = self._shape_guard_cratonic_stability(
            grid, ctype_new, origin_new, stability_new)

        # ---- plate velocity vector field -------------------------------------
        vel = self._plate_velocity(grid, plate_new, plates, vigor, regime_code)

        # ---- boundary network -------------------------------------------------
        separating_axis = separating & plate_contact
        converging_axis = converging & plate_contact
        collision_axis = collision & plate_contact
        subduction_axis = subduction & plate_contact
        suture_axis = suture & plate_contact
        active_margin_axis = active_margin & plate_contact
        passive_margin_axis = passive_margin & plate_contact
        transform_axis = transform & plate_contact
        if not transform_axis.any():
            excluded_transform = (
                ridge_axis
                | trench_axis
                | collision_axis
                | subduction_axis
                | suture_axis
                | active_margin_axis
                | passive_margin_axis
            )
            transform_axis = self._ridge_offset_transform_axis(
                grid, ridge_axis, plate_contact, ctype_new, excluded_transform)
        p107_boundary_network = world.g(
            "tectonics.enable_p107_boundary_network_continuity",
            world.g("tectonics.enable_p107_ranked_plate_policy", 0.0),
        ) > 0.0
        p108_width_guard = world.g("tectonics.enable_p108_boundary_width_guard", 0.0) > 0.0
        if p107_boundary_network:
            ridge_out = self._continuous_boundary_skeleton(
                grid, ridge_axis, allowed=(ctype_new == OCEAN), bridge_passes=2,
                sample_spacing=(5 if p108_width_guard else 3),
                max_axis_seeds=(360 if p108_width_guard else 520))
        else:
            ridge_out = self._thin_boundary_mask(grid, ridge_axis, spacing=3)
        separating_out = self._thin_boundary_mask(grid, separating_axis, spacing=3)
        converging_out = self._thin_boundary_mask(grid, converging_axis, spacing=3)
        collision_out = self._thin_boundary_mask(grid, collision_axis, spacing=3)
        subduction_out = self._thin_boundary_mask(grid, subduction_axis, spacing=3)
        trench_spacing = 2 if p107_boundary_network else 3
        trench_out = self._thin_boundary_mask(grid, trench_axis, spacing=trench_spacing)
        suture_out = self._thin_boundary_mask(grid, suture_axis, spacing=3)
        active_margin_out = self._thin_boundary_mask(grid, active_margin_axis, spacing=2)
        passive_margin_out = self._thin_boundary_mask(grid, passive_margin_axis, spacing=3)
        if p107_boundary_network:
            transform_out = self._continuous_boundary_skeleton(
                grid, transform_axis, allowed=(plate_contact | (ctype_new == OCEAN)),
                bridge_passes=1,
                sample_spacing=(4 if p108_width_guard else 3),
                max_axis_seeds=(260 if p108_width_guard else 520))
        else:
            transform_out = self._thin_boundary_mask(grid, transform_axis, spacing=2)
        min_transform_cells = max(
            4,
            int(np.ceil(0.45 * float(np.count_nonzero(ridge_out)))),
        )
        if int(np.count_nonzero(transform_out)) < min_transform_cells:
            excluded_transform_out = (
                ridge_out
                | trench_out
                | collision_out
                | subduction_out
                | suture_out
                | active_margin_out
                | passive_margin_out
                | transform_out
            )
            transform_out |= self._ridge_offset_transform_axis(
                grid, ridge_out, plate_contact, ctype_new, excluded_transform_out)
            if p107_boundary_network:
                transform_out = self._continuous_boundary_skeleton(
                    grid, transform_out,
                    allowed=(plate_contact | (ctype_new == OCEAN)),
                    bridge_passes=1,
                    sample_spacing=(4 if p108_width_guard else 3),
                    max_axis_seeds=(260 if p108_width_guard else 520))
        convergent_parent_source = (
            converging_axis
            | collision_axis
            | subduction_axis
            | trench_axis
            | suture_axis
            | active_margin_axis
        )
        convergent_parent_out = self._p11167_ordered_convergent_parent_axis(
            grid,
            convergent_parent_source,
            support=convergent_parent_source | converging_out | subduction_out
            | trench_out | active_margin_out | collision_out | suture_out,
            excluded=ridge_out | transform_out | passive_margin_out,
            sample_spacing=2,
            max_axis_seeds=(360 if p108_width_guard else 520),
        )
        subduction_parent_source = subduction_axis | trench_axis | active_margin_axis
        subduction_parent_out = self._p11167_ordered_convergent_parent_axis(
            grid,
            subduction_parent_source,
            support=subduction_parent_source | subduction_out | trench_out
            | active_margin_out,
            excluded=ridge_out | transform_out | collision_out | suture_out
            | passive_margin_out,
            sample_spacing=2,
            max_axis_seeds=(300 if p108_width_guard else 420),
        )
        self._record_p11167_parent_line_metrics(
            world,
            grid,
            convergent_parent_source,
            convergent_parent_out,
            subduction_parent_source,
            subduction_parent_out,
        )
        boundaries = {
            "divergent": np.where(ridge_out | separating_out)[0].astype(np.int64),
            "ridge": np.where(ridge_out)[0].astype(np.int64),
            "convergent": np.where(converging_out)[0].astype(np.int64),
            "collision": np.where(collision_out)[0].astype(np.int64),
            "subduction": np.where(subduction_out)[0].astype(np.int64),
            "trench": np.where(trench_out)[0].astype(np.int64),
            "suture": np.where(suture_out)[0].astype(np.int64),
            "active_margin": np.where(active_margin_out)[0].astype(np.int64),
            "passive_margin": np.where(passive_margin_out)[0].astype(np.int64),
            "transform": np.where(transform_out)[0].astype(np.int64),
        }
        prev_boundary_objects = world.objects.get("tectonics.boundary_objects", [])
        boundary_objects = self._boundary_objects(
            grid, boundaries, plate_new, t,
            ctype=ctype_new,
            age=age_new,
            thick=thick_new,
            velocity=vel,
            previous_objects=prev_boundary_objects,
        )
        boundary_objects.extend(self._boundary_objects(
            grid,
            {
                "convergent_parent": np.where(convergent_parent_out)[0].astype(np.int64),
                "subduction_parent": np.where(subduction_parent_out)[0].astype(np.int64),
            },
            plate_new,
            t,
            ctype=ctype_new,
            age=age_new,
            thick=thick_new,
            velocity=vel,
                previous_objects=prev_boundary_objects,
            ))
        boundary_polylines = self._boundary_polyline_objects(
            grid,
            boundaries,
            boundary_objects,
            plate_new,
            t,
        )
        self._record_p129_boundary_polyline_metrics(
            world,
            grid,
            boundary_polylines,
        )
        boundary_provinces, boundary_province_kind = self._boundary_province_objects(
            grid, boundary_objects, ctype_new, age_new, t)
        prev_continent_id = world.get_field("tectonics.continent_id", -1.0)
        next_continent_id = int(world.g("tectonics.next_continent_id", 0.0))
        continent_id, _, next_continent_id = _stable_component_ids(
            grid, ctype_new == CONT, prev_continent_id, next_continent_id)
        world.set_g("tectonics.next_continent_id", float(next_continent_id))
        domain_new = self._crust_domain_field(grid, ctype_new, origin_new, stability_new)
        craton_objects = self._craton_objects(grid, origin_new == ORIGIN_CRATON, age_new, t)
        platform_subsidence = self._platform_subsidence_field(
            world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
            stability_new, domain_new, t)
        origin_new, age_new, thick_new, reworked_new, stability_new = (
            self._stabilize_basement_cratonic_current_state(
                world, grid, ctype_new, origin_new, age_new, thick_new,
                reworked_new, stability_new, domain_new, basement_age_new,
                basement_stability_new, basement_thick_new, t,
                active_boundary_mask=collision | subduction | separating | ridge)
        )
        domain_new = self._crust_domain_field(grid, ctype_new, origin_new, stability_new)
        craton_objects = self._craton_objects(grid, origin_new == ORIGIN_CRATON, age_new, t)
        platform_subsidence = self._platform_subsidence_field(
            world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
            stability_new, domain_new, t)
        cratonic_province_objects = self._cratonic_province_objects(
            world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
            stability_new, domain_new, continent_id, t,
            platform_subsidence=platform_subsidence,
            province_code=province_code_new,
            basement_age_floor=basement_age_new,
            basement_stability_floor=basement_stability_new,
            basement_thickness_floor=basement_thick_new)
        province_code_new = self._update_cratonic_province_memory(
            world, grid, ctype_new, origin_new, age_new, reworked_new,
            stability_new, domain_new, province_code_new, cratonic_province_objects,
            t)
        basement_age_new, basement_stability_new, basement_thick_new = (
            self._update_basement_cargo_archive(
                world, grid, ctype_new, origin_new, age_new, thick_new,
                reworked_new, stability_new, domain_new, province_code_new,
                cratonic_province_objects, basement_age_new,
                basement_stability_new, basement_thick_new, t,
                active_boundary_mask=collision | subduction | separating | ridge)
        )
        stable_cratonic_support = self._stable_cratonic_lithosphere_support(
            world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
            stability_new, domain_new, province_code_new, cratonic_province_objects,
            basement_age_new, basement_stability_new, basement_thick_new, t,
            active_boundary_mask=collision | subduction | separating | ridge)
        origin_new, age_new, thick_new, reworked_new, stability_new = (
            self._project_supported_cratonic_current_state(
                world, grid, ctype_new, origin_new, age_new, thick_new,
                reworked_new, stability_new, domain_new, province_code_new,
                stable_cratonic_support, basement_age_new, basement_stability_new,
                basement_thick_new, t, dt,
                active_boundary_mask=collision | subduction | separating | ridge)
        )
        domain_new = self._crust_domain_field(grid, ctype_new, origin_new, stability_new)
        craton_objects = self._craton_objects(grid, origin_new == ORIGIN_CRATON, age_new, t)
        age_new, reworked_new = self._refresh_supported_ancient_age_continuity(
            world, grid, ctype_new, age_new, reworked_new, stability_new, domain_new,
            province_code_new, stable_cratonic_support, basement_age_new,
            basement_stability_new, basement_thick_new, t,
            active_boundary_mask=collision | subduction | separating | ridge)
        if world.g("tectonics.enable_p65_post_projection_cratonic_refresh", 0.0) > 0.0:
            domain_new = self._crust_domain_field(grid, ctype_new, origin_new, stability_new)
            craton_objects = self._craton_objects(grid, origin_new == ORIGIN_CRATON, age_new, t)
            platform_subsidence = self._platform_subsidence_field(
                world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
                stability_new, domain_new, t)
            cratonic_province_objects = self._cratonic_province_objects(
                world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
                stability_new, domain_new, continent_id, t,
                platform_subsidence=platform_subsidence,
                province_code=province_code_new,
                basement_age_floor=basement_age_new,
                basement_stability_floor=basement_stability_new,
                basement_thickness_floor=basement_thick_new)
            province_before = np.asarray(province_code_new, dtype=np.float64).copy()
            province_code_new = self._update_cratonic_province_memory(
                world, grid, ctype_new, origin_new, age_new, reworked_new,
                stability_new, domain_new, province_code_new, cratonic_province_objects,
                t)
            cratonic_province_objects = self._cratonic_province_objects(
                world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
                stability_new, domain_new, continent_id, t,
                platform_subsidence=platform_subsidence,
                province_code=province_code_new,
                basement_age_floor=basement_age_new,
                basement_stability_floor=basement_stability_new,
                basement_thickness_floor=basement_thick_new)
            cratonic_province_objects, province_code_new = (
                self._augment_post_projection_supported_platform_objects(
                    world, grid, ctype_new, age_new, thick_new, reworked_new,
                    stability_new, domain_new, continent_id, province_code_new,
                    cratonic_province_objects, stable_cratonic_support,
                    basement_age_new, basement_stability_new, basement_thick_new,
                    t, active_boundary_mask=collision | subduction | separating | ridge)
            )
            cont_mask_p65 = ctype_new == CONT
            cont_area_p65 = max(float(grid.cell_area[cont_mask_p65].sum()), 1.0e-12)
            province_changed_p65 = cont_mask_p65 & (province_code_new != province_before)
            world.set_g(
                "tectonics.last_p65_province_code_changed_share_of_continental_crust",
                float(grid.cell_area[province_changed_p65].sum() / cont_area_p65),
            )
            world.set_g(
                "tectonics.last_p65_refreshed_p48_shield_share_of_continental_crust",
                float(world.g("tectonics.last_p48_shield_share_of_continental_crust", 0.0)),
            )
            world.set_g(
                "tectonics.last_p65_refreshed_p48_platform_share_of_continental_crust",
                float(world.g("tectonics.last_p48_platform_share_of_continental_crust", 0.0)),
            )
            world.set_g(
                "tectonics.last_p65_refreshed_p48_mature_share_of_continental_crust",
                float(world.g("tectonics.last_p48_shield_share_of_continental_crust", 0.0))
                + float(world.g("tectonics.last_p48_platform_share_of_continental_crust", 0.0))
                + float(world.g("tectonics.last_p48_interior_basin_share_of_continental_crust", 0.0)),
            )
        if world.g("tectonics.enable_p68_cratonic_balance_physical_ensemble", 0.0) > 0.0:
            (
                origin_new, age_new, thick_new, reworked_new, stability_new,
                province_code_new, basement_age_new, basement_stability_new,
                basement_thick_new,
            ) = self._rebalance_physical_ensemble_cratonic_state(
                world, grid, ctype_new, origin_new, age_new, thick_new,
                reworked_new, stability_new, domain_new, continent_id,
                province_code_new, cratonic_province_objects,
                stable_cratonic_support, basement_age_new, basement_stability_new,
                basement_thick_new, t,
                active_boundary_mask=collision | subduction | separating | ridge)
            domain_new = self._crust_domain_field(grid, ctype_new, origin_new, stability_new)
            craton_objects = self._craton_objects(grid, origin_new == ORIGIN_CRATON, age_new, t)
            platform_subsidence = self._platform_subsidence_field(
                world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
                stability_new, domain_new, t)
            cratonic_province_objects = self._cratonic_province_objects(
                world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
                stability_new, domain_new, continent_id, t,
                platform_subsidence=platform_subsidence,
                province_code=province_code_new,
                basement_age_floor=basement_age_new,
                basement_stability_floor=basement_stability_new,
                basement_thickness_floor=basement_thick_new)
            province_code_new = self._update_cratonic_province_memory(
                world, grid, ctype_new, origin_new, age_new, reworked_new,
                stability_new, domain_new, province_code_new, cratonic_province_objects,
                t)
            cratonic_province_objects = self._cratonic_province_objects(
                world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
                stability_new, domain_new, continent_id, t,
                platform_subsidence=platform_subsidence,
                province_code=province_code_new,
                basement_age_floor=basement_age_new,
                basement_stability_floor=basement_stability_new,
                basement_thickness_floor=basement_thick_new)
            stable_cratonic_support = self._stable_cratonic_lithosphere_support(
                world, grid, ctype_new, age_new, thick_new, origin_new, reworked_new,
                stability_new, domain_new, province_code_new, cratonic_province_objects,
                basement_age_new, basement_stability_new, basement_thick_new, t,
                active_boundary_mask=collision | subduction | separating | ridge)
            cont_mask_p68 = ctype_new == CONT
            cont_area_p68 = max(float(grid.cell_area[cont_mask_p68].sum()), 1.0e-12)
            world.set_g(
                "tectonics.last_p68_final_p48_shield_share_of_continental_crust",
                float(world.g("tectonics.last_p48_shield_share_of_continental_crust", 0.0)),
            )
            world.set_g(
                "tectonics.last_p68_final_p48_platform_share_of_continental_crust",
                float(world.g("tectonics.last_p48_platform_share_of_continental_crust", 0.0)),
            )
            world.set_g(
                "tectonics.last_p68_final_p48_mature_share_of_continental_crust",
                float(world.g("tectonics.last_p48_shield_share_of_continental_crust", 0.0))
                + float(world.g("tectonics.last_p48_platform_share_of_continental_crust", 0.0))
                + float(world.g("tectonics.last_p48_interior_basin_share_of_continental_crust", 0.0)),
            )
            world.set_g(
                "tectonics.last_p68_final_ancient_fraction_of_continental_crust",
                float(grid.cell_area[cont_mask_p68 & (age_new >= 2500.0)].sum() / cont_area_p68),
            )
            world.set_g(
                "tectonics.last_p68_final_old_archive_fraction_of_continental_crust",
                float(grid.cell_area[
                    cont_mask_p68
                    & (basement_age_new >= 2500.0)
                    & (basement_stability_new >= 0.70)
                ].sum() / cont_area_p68),
            )
        internal_block_code_new = self._update_internal_geographic_block_code(
            world,
            grid,
            ctype_new,
            origin_new,
            age_new,
            thick_new,
            reworked_new,
            stability_new,
            domain_new,
            province_code_new,
            internal_block_code_new,
            platform_subsidence,
            t,
            active_boundary_mask=collision | subduction | separating | ridge,
        )
        next_internal_block_id = int(
            world.g("tectonics.next_internal_geographic_block_id", 1.0))
        (
            internal_block_id_new,
            internal_block_objects,
            next_internal_block_id,
        ) = self._internal_geographic_block_ids_and_objects(
            world,
            grid,
            ctype_new,
            internal_block_code_new,
            internal_block_id_new,
            next_internal_block_id,
            age_new,
            thick_new,
            stability_new,
            continent_id,
            t,
        )
        world.set_g(
            "tectonics.next_internal_geographic_block_id",
            float(next_internal_block_id),
        )
        terrane_mask = self._terrane_mask(grid, ctype_new, origin_new, stability_new)
        prev_terrane_id = world.get_field("tectonics.terrane_id", -1.0)
        next_terrane_id = int(world.g("tectonics.next_terrane_id", 0.0))
        terrane_id, _, next_terrane_id = _stable_component_ids(
            grid, terrane_mask, prev_terrane_id, next_terrane_id)
        world.set_g("tectonics.next_terrane_id", float(next_terrane_id))
        continent_objects = self._continent_objects(
            grid, continent_id, ctype_new, age_new, thick_new, origin_new, stability_new, t)
        terrane_objects = self._terrane_objects(
            grid, terrane_id, continent_id, ctype_new, origin_new, age_new, thick_new, t)
        plate_topologies = self._plate_topology_objects(
            grid, plate_new, plates, ctype_new, continent_id, terrane_id, boundaries, t)
        (
            plate_rank,
            protected_plate_id,
            microplate_objects,
        ) = self._plate_rank_fields_and_objects(
            grid,
            plate_new,
            plate_topologies,
            ctype_new,
            origin_new,
            age_new,
            boundaries,
            t,
            protected_records=(
                reorg_detail.get("protected_microplates", [])
                if isinstance(reorg_detail, dict) else []
            ),
            enable_protection=bool(
                world.g("tectonics.enable_p107_ranked_plate_policy", 0.0) > 0.0),
        )
        plate_accretion_events = (
            world.objects.get("tectonics.plate_accretion_events", [])
            + self._plate_accretion_event_objects(reorg_detail, t)
        )[-128:]
        new_lineage_objects, lineage_events = self._continent_lifecycle_events(
            grid,
            prev_continent_id,
            continent_id,
            prev_terrane_id,
            terrane_id,
            world.objects.get("tectonics.terranes", []),
            terrane_objects,
            reorg_detail,
            t,
        )
        continent_lifecycle_events = (
            world.objects.get("tectonics.continent_lifecycle_events", [])
            + new_lineage_objects
        )[-96:]
        lifecycle_objects = self._wilson_lifecycle_objects(
            boundary_objects,
            previous={
                "ocean_basins": world.objects.get("tectonics.ocean_basins", []),
                "rift_systems": world.objects.get("tectonics.rift_systems", []),
                "passive_margins": world.objects.get("tectonics.passive_margins", []),
                "spreading_centers": world.objects.get("tectonics.spreading_centers", []),
                "closing_margins": world.objects.get("tectonics.closing_margins", []),
                "sutures": world.objects.get("tectonics.sutures", []),
            },
            t=t,
        )
        breakup_seaways = self._breakup_seaway_objects(
            world,
            grid,
            ctype_new,
            domain_new,
            stability_new,
            rift_potential,
            continent_id,
            plate_new,
            boundaries,
            plate_topologies,
            lifecycle_objects["tectonics.rift_systems"],
            t,
        )
        self._apply_breakup_seaway_crustal_opening(
            world,
            grid,
            breakup_seaways,
            ctype_new,
            thick_new,
            age_new,
            origin_new,
            reworked_new,
            stability_new,
            domain_new,
            continent_id,
            terrane_id,
            province_code_new,
            internal_block_id_new,
            internal_block_code_new,
            basement_age_new,
            basement_stability_new,
            basement_thick_new,
            orog_age,
            volc_age,
            t,
        )
        frac = float(np.mean(ctype_new == CONT))
        wilson_cycles = lifecycle_objects["tectonics.wilson_cycles"]
        ocean_gateways = lifecycle_objects["tectonics.ocean_gateways"]
        wilson_phase = self._wilson_phase_field(grid, wilson_cycles)
        if self._should_emit_wilson_events(t, dt):
            events += self._wilson_events_from_objects(
                boundary_objects, wilson_cycles, ocean_gateways, t)
        events += lineage_events
        events += self._cratonization_events(grid, cratonized, t)
        deformation_intensity, deformation_style, deforming_networks = (
            self._deforming_network_state(
                grid,
                collision,
                subduction,
                separating,
                transform,
                cont_within,
                plate_contact,
                ctype_new,
                t,
            )
        )

        # ---- sampled events ---------------------------------------------------
        rng = rng_key.child("events").generator()
        events += self._sample_events(grid, t, rng, collision, "collision",
                                      "orogeny from continental collision", orog_age)
        events += self._sample_events(grid, t, rng, arc, "subduction",
                                      "volcanic arc above subduction", volc_age)
        rift_cells = (ridge_axis | separating_axis) & self._adjacent_to_continent(grid, ctype_new)
        events += self._sample_events(grid, t, rng, rift_cells, "rift",
                                      "continental rifting / new ocean basin", None)

        # ---- notable objects --------------------------------------------------
        volcanoes = world.object_set("tectonics.volcanoes")
        for c in np.where(arc)[0][::max(1, int(arc.sum()) // 20 or 1)][:20]:
            volcanoes.append({"cell": int(c), "age_myr": round(t, 1), "kind": "arc"})

        world.provenance.record(Provenance(
            "crust.thickness_m", self.name, self.fidelity, "m",
            direct_cause=f"vigor={vigor:.2f}, regime={regime_code:.0f}; "
            f"{int(collision.sum())} collision / {int(subduction.sum())} subduction cells",
            upstream_events=[e.id for e in events[:4]]))

        fields = {
            "crust.type": ctype_new,
            "crust.thickness_m": thick_new,
            "crust.age_myr": age_new,
            "crust.origin": origin_new,
            "crust.reworked_age_myr": reworked_new,
            "crust.stability": stability_new,
            "crust.domain": domain_new,
            "tectonics.continent_id": continent_id,
            "tectonics.terrane_id": terrane_id,
            "tectonics.cratonic_province_code": province_code_new,
            "tectonics.internal_geographic_block_id": internal_block_id_new,
            "tectonics.internal_geographic_block_code": internal_block_code_new,
            "tectonics.basement_age_floor_myr": basement_age_new,
            "tectonics.basement_stability_floor": basement_stability_new,
            "tectonics.basement_thickness_floor_m": basement_thick_new,
            "tectonics.stable_cratonic_lithosphere_support": stable_cratonic_support,
            "tectonics.plate_id": plate_new.astype(np.float64),
            "tectonics.plate_rank": plate_rank,
            "tectonics.protected_plate_id": protected_plate_id,
            "tectonics.boundary_province_kind": boundary_province_kind,
            "tectonics.plate_velocity": vel,
            "tectonics.orogeny_age_myr": orog_age,
            "tectonics.volcanism_age_myr": volc_age,
            "tectonics.deformation_intensity": deformation_intensity,
            "tectonics.deformation_style": deformation_style,
            "tectonics.platform_subsidence": platform_subsidence,
            "archive.wilson_cycle_phase": wilson_phase,
        }
        plume_objects = world.objects.get("tectonics.plumes", [])[-24:] + plume_objects
        lip_objects = world.objects.get("tectonics.lips", [])[-24:] + lip_objects
        objects = {
            "tectonics.boundary_objects": boundary_objects,
            "tectonics.boundary_polylines": boundary_polylines,
            "tectonics.boundary_provinces": boundary_provinces,
            "tectonics.wilson_cycles": wilson_cycles,
            "tectonics.ocean_gateways": ocean_gateways,
            "tectonics.ocean_basins": lifecycle_objects["tectonics.ocean_basins"],
            "tectonics.rift_systems": lifecycle_objects["tectonics.rift_systems"],
            "tectonics.breakup_seaways": breakup_seaways,
            "tectonics.passive_margins": lifecycle_objects["tectonics.passive_margins"],
            "tectonics.spreading_centers": lifecycle_objects["tectonics.spreading_centers"],
            "tectonics.closing_margins": lifecycle_objects["tectonics.closing_margins"],
            "tectonics.sutures": lifecycle_objects["tectonics.sutures"],
            "tectonics.deforming_networks": deforming_networks,
            "tectonics.cratons": craton_objects,
            "tectonics.shields": cratonic_province_objects["tectonics.shields"],
            "tectonics.platforms": cratonic_province_objects["tectonics.platforms"],
            "tectonics.interior_basins": cratonic_province_objects["tectonics.interior_basins"],
            "tectonics.internal_geographic_blocks": internal_block_objects,
            "tectonics.continents": continent_objects,
            "tectonics.terranes": terrane_objects,
            "tectonics.continent_lifecycle_events": continent_lifecycle_events,
            "tectonics.plate_topologies": plate_topologies,
            "tectonics.microplates": microplate_objects,
            "tectonics.plate_accretion_events": plate_accretion_events,
            "tectonics.plumes": plume_objects,
            "tectonics.lips": lip_objects,
        }
        networks = {"tectonics.boundaries": boundaries}
        diag = {"cont_fraction": frac, "n_collision": int(collision.sum()),
                "n_subduction": int(subduction.sum()), "n_ridge": int(ridge.sum()),
                "n_trench": int(trench.sum()), "n_suture": int(suture.sum()),
                "n_plate_component_merges": component_merges,
                "n_boundary_objects": len(boundary_objects),
                "n_boundary_polylines": len(boundary_polylines),
                "n_boundary_provinces": len(boundary_provinces),
                "n_wilson_cycles": len(wilson_cycles),
                "n_lifecycle_ocean_basins": len(lifecycle_objects["tectonics.ocean_basins"]),
                "n_breakup_seaways": len(breakup_seaways),
                "n_ocean_gateways": len(ocean_gateways),
                "n_cratons": len(craton_objects),
                "n_shields": len(cratonic_province_objects["tectonics.shields"]),
                "n_platforms": len(cratonic_province_objects["tectonics.platforms"]),
                "n_interior_basins": len(cratonic_province_objects["tectonics.interior_basins"]),
                "n_internal_geographic_blocks": len(internal_block_objects),
                "n_continents": len(continent_objects),
                "n_terranes": len(terrane_objects),
                "n_continent_lifecycle_events": len(continent_lifecycle_events),
                "n_plate_topologies": len(plate_topologies),
                "n_lips": len(lip_objects),
                "n_rifted_margin_candidates": int(rifted_margin_sources.sum()),
                "raw_boundary_cells": int(np.unique(np.concatenate([
                    np.where(ridge_axis | separating_axis)[0],
                    np.where(converging_axis)[0],
                    np.where(collision_axis)[0],
                    np.where(subduction_axis)[0],
                ])).size),
                "generalized_boundary_cells": int(np.unique(np.concatenate([
                    boundaries["divergent"],
                    boundaries["convergent"],
                    boundaries["collision"],
                    boundaries["subduction"],
                ])).size)}
        return StepResult(state_delta={"fields": fields, "networks": networks,
                                       "objects": objects},
                          events=events, diagnostics=diag)

    # ------------------------------------------------------------------
    def _plate_contact_mask(self, grid, plate):
        contact = np.zeros(grid.n, dtype=bool)
        i, j = grid.edges[:, 0], grid.edges[:, 1]
        diff = plate[i] != plate[j]
        contact[i[diff]] = True
        contact[j[diff]] = True
        return contact

    def _thin_boundary_mask(self, grid, mask, spacing=3):
        mask = np.asarray(mask, dtype=bool)
        if not mask.any():
            return mask.copy()
        out = np.zeros(grid.n, dtype=bool)
        for comp in self._connected_components(grid, mask):
            if comp.size <= 4:
                out[comp] = True
                continue
            blocked = np.zeros(grid.n, dtype=bool)
            scores = []
            for c in comp:
                c = int(c)
                scores.append((np.count_nonzero(mask[grid.neighbors[c]]), -c, c))
            for _, _, c in sorted(scores, reverse=True):
                if blocked[c]:
                    continue
                out[c] = True
                seed = np.zeros(grid.n, dtype=bool)
                seed[c] = True
                blocked |= self._dilate_mask(grid, seed, passes=spacing)
        return out

    def _continuous_boundary_skeleton(self, grid, mask, *, allowed=None,
                                      bridge_passes=1, min_component_cells=3,
                                      sample_spacing=3, max_axis_seeds=520):
        """Keep boundary axes continuous while removing isolated process noise.

        `_thin_boundary_mask` intentionally samples long masks at intervals.
        That is useful for sparse event markers, but it makes modern ridge and
        trench systems render and audit as dotted fragments.  P107 needs the
        opposite behavior: preserve a one-cell-to-few-cell continuous graph,
        bridge short gaps on the spherical neighbor graph, and drop only tiny
        isolated components that have no network continuity.
        """
        base = np.asarray(mask, dtype=bool)
        if not base.any():
            return base.copy()
        if allowed is None:
            allowed_mask = np.ones(grid.n, dtype=bool)
        else:
            allowed_mask = np.asarray(allowed, dtype=bool).copy()
            if allowed_mask.shape != (grid.n,):
                allowed_mask = np.ones(grid.n, dtype=bool)
        bridged = base & allowed_mask
        if not bridged.any():
            bridged = base.copy()
            allowed_mask |= bridged

        for _ in range(max(0, int(bridge_passes))):
            labels, _ = _component_labels(grid, bridged)
            add = np.zeros(grid.n, dtype=bool)
            for c in np.where(allowed_mask & ~bridged)[0]:
                nbs = np.asarray(grid.neighbors[int(c)], dtype=np.int64)
                near = sorted({int(labels[int(nb)]) for nb in nbs if labels[int(nb)] >= 0})
                if len(near) >= 2:
                    add[int(c)] = True
                    continue
                if len(near) == 1 and np.count_nonzero(bridged[nbs]) >= 2:
                    add[int(c)] = True
            if not add.any():
                break
            bridged |= add

        out = np.zeros(grid.n, dtype=bool)
        for comp in self._connected_components(grid, bridged):
            comp = np.asarray(comp, dtype=np.int64)
            if comp.size <= max(6, int(min_component_cells) * 2):
                out[comp] = True
                continue
            comp_mask = np.zeros(grid.n, dtype=bool)
            comp_mask[comp] = True
            sampled = self._thin_boundary_mask(
                grid, comp_mask, spacing=max(1, int(sample_spacing))) & comp_mask
            seeds = np.where(sampled)[0]
            if seeds.size < 2:
                out[comp] = True
                continue
            max_seeds = max(8, int(max_axis_seeds))
            if seeds.size > max_seeds:
                order = self._boundary_component_axis_order(grid, comp, seeds)
                step = int(np.ceil(seeds.size / float(max_seeds)))
                seeds = order[::step]
                if int(order[-1]) not in set(int(x) for x in seeds):
                    seeds = np.concatenate([seeds, order[-1:]])
            out |= self._connect_ordered_boundary_seeds(grid, comp_mask, seeds)

        comps = self._connected_components(grid, out)
        if not comps:
            return out
        kept = np.zeros(grid.n, dtype=bool)
        for comp in comps:
            comp = np.asarray(comp, dtype=np.int64)
            if comp.size >= int(min_component_cells):
                kept[comp] = True
                continue
            # Keep a small spur if it is attached to a larger component through
            # the original unbridged evidence; otherwise treat it as noise.
            touches = 0
            for c in comp:
                touches += int(np.count_nonzero(base[grid.neighbors[int(c)]]))
            if touches >= 2:
                kept[comp] = True
        if not kept.any():
            return out
        return kept

    def _p11167_ordered_convergent_parent_axis(self, grid, source, *,
                                               support=None,
                                               excluded=None,
                                               sample_spacing=2,
                                               max_axis_seeds=420):
        """Build a narrow parent line for convergent/subduction terrain response.

        Raw convergent masks are useful diagnostics, but they are too scattered
        to act as parent geometry for trenches, arcs, and orogens.  This helper
        turns only local process-supported cells into a one-cell graph and
        rejects broad area growth.
        """
        source = np.asarray(source, dtype=bool)
        if source.shape != (grid.n,) or not source.any():
            return np.zeros(grid.n, dtype=bool)
        if support is None:
            support_mask = source.copy()
        else:
            support_mask = np.asarray(support, dtype=bool)
            if support_mask.shape != (grid.n,):
                support_mask = source.copy()
        if excluded is None:
            excluded_mask = np.zeros(grid.n, dtype=bool)
        else:
            excluded_mask = np.asarray(excluded, dtype=bool)
            if excluded_mask.shape != (grid.n,):
                excluded_mask = np.zeros(grid.n, dtype=bool)

        base = source & ~excluded_mask
        if not base.any():
            return base
        support_mask = (support_mask | base) & ~excluded_mask
        allowed = self._dilate_mask(grid, support_mask, passes=2) & ~excluded_mask
        axis = self._continuous_boundary_skeleton(
            grid,
            base,
            allowed=allowed,
            bridge_passes=2,
            min_component_cells=3,
            sample_spacing=max(1, int(sample_spacing)),
            max_axis_seeds=max_axis_seeds,
        )
        axis &= allowed
        if not axis.any():
            return base

        support_axis = np.zeros(grid.n, dtype=bool)
        for comp in self._connected_components(grid, support_mask & allowed):
            comp = np.asarray(comp, dtype=np.int64)
            seeds = comp[base[comp]]
            if seeds.size < 2:
                continue
            if seeds.size > max_axis_seeds:
                order = self._boundary_component_axis_order(grid, comp, seeds)
                step = int(np.ceil(seeds.size / float(max(2, max_axis_seeds))))
                seeds = order[::step]
                if int(order[-1]) not in set(int(x) for x in seeds):
                    seeds = np.concatenate([seeds, order[-1:]])
            comp_allowed = np.zeros(grid.n, dtype=bool)
            comp_allowed[comp] = True
            support_axis |= self._connect_ordered_boundary_seeds(
                grid, comp_allowed & allowed, seeds)
        if support_axis.any():
            axis_components = len(self._connected_components(grid, axis))
            support_components = len(self._connected_components(grid, support_axis))
            base_components = len(self._connected_components(grid, base))
            if support_components < axis_components or axis_components >= base_components:
                axis = support_axis

        area = np.asarray(grid.cell_area, dtype=np.float64)
        base_area = max(float(area[base].sum()), 1.0e-12)
        axis_area = float(area[axis].sum())
        line_area_cap = max(base_area * 4.0, base_area + 12.0 * float(area.mean()))
        if axis_area > line_area_cap:
            tighter_allowed = self._dilate_mask(grid, base | support_mask, passes=1) & ~excluded_mask
            tighter = self._continuous_boundary_skeleton(
                grid,
                base,
                allowed=tighter_allowed,
                bridge_passes=1,
                min_component_cells=3,
                sample_spacing=max(2, int(sample_spacing)),
                max_axis_seeds=max(80, int(max_axis_seeds * 0.75)),
            )
            if tighter.any():
                tighter_area = float(area[tighter].sum())
                if tighter_area <= line_area_cap:
                    axis = tighter & tighter_allowed
        return axis & ~excluded_mask

    def _record_p11167_parent_line_metrics(self, world, grid,
                                           convergent_source,
                                           convergent_parent,
                                           subduction_source,
                                           subduction_parent):
        def record(prefix: str, before_mask, after_mask) -> None:
            before = np.asarray(before_mask, dtype=bool)
            after = np.asarray(after_mask, dtype=bool)
            before_components = len(self._connected_components(grid, before)) if before.any() else 0
            after_components = len(self._connected_components(grid, after)) if after.any() else 0
            area = np.asarray(grid.cell_area, dtype=np.float64)
            total = max(float(area.sum()), 1.0e-12)
            world.set_g(
                f"tectonics.last_p11167_{prefix}_parent_components_before",
                float(before_components),
            )
            world.set_g(
                f"tectonics.last_p11167_{prefix}_parent_components_after",
                float(after_components),
            )
            world.set_g(
                f"tectonics.last_p11167_{prefix}_parent_cells_before",
                float(np.count_nonzero(before)),
            )
            world.set_g(
                f"tectonics.last_p11167_{prefix}_parent_cells_after",
                float(np.count_nonzero(after)),
            )
            world.set_g(
                f"tectonics.last_p11167_{prefix}_parent_area_fraction_before",
                float(area[before].sum() / total) if before.any() else 0.0,
            )
            world.set_g(
                f"tectonics.last_p11167_{prefix}_parent_area_fraction_after",
                float(area[after].sum() / total) if after.any() else 0.0,
            )
            accepted = after.any() and (
                before_components <= 0 or after_components <= before_components
            )
            world.set_g(
                f"tectonics.last_p11167_{prefix}_parent_refinement_accepted",
                1.0 if accepted else 0.0,
            )

        record("convergent", convergent_source, convergent_parent)
        record("subduction", subduction_source, subduction_parent)

    def _boundary_component_axis_order(self, grid, comp, seeds):
        comp = np.asarray(comp, dtype=np.int64)
        seeds = np.asarray(seeds, dtype=np.int64)
        if seeds.size <= 1:
            return seeds
        xyz = grid.xyz[comp]
        centered = xyz - xyz.mean(axis=0)
        try:
            _, _, vh = np.linalg.svd(centered, full_matrices=False)
            axis = vh[0]
        except np.linalg.LinAlgError:
            axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        proj = grid.xyz[seeds] @ axis
        tie = _stable_tie_breaker(seeds)
        return seeds[np.lexsort((tie, proj))]

    def _connect_ordered_boundary_seeds(self, grid, allowed, seeds):
        allowed = np.asarray(allowed, dtype=bool)
        seeds = np.asarray(seeds, dtype=np.int64)
        seeds = seeds[(0 <= seeds) & (seeds < grid.n) & allowed[seeds]]
        if seeds.size == 0:
            return np.zeros(grid.n, dtype=bool)
        ordered = self._boundary_component_axis_order(
            grid, np.where(allowed)[0], seeds)
        out = np.zeros(grid.n, dtype=bool)
        out[ordered] = True
        for a, b in zip(ordered[:-1], ordered[1:]):
            path = self._shortest_path_within_mask(grid, int(a), int(b), allowed)
            if path.size:
                out[path] = True
        return out

    def _shortest_path_within_mask(self, grid, start, target, allowed):
        start = int(start)
        target = int(target)
        allowed = np.asarray(allowed, dtype=bool)
        if start == target:
            return np.asarray([start], dtype=np.int64)
        if not (0 <= start < grid.n and 0 <= target < grid.n):
            return np.asarray([], dtype=np.int64)
        if not (allowed[start] and allowed[target]):
            return np.asarray([], dtype=np.int64)
        parent = np.full(grid.n, -2, dtype=np.int64)
        parent[start] = -1
        queue = [start]
        head = 0
        while head < len(queue):
            c = queue[head]
            head += 1
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if not allowed[nb] or parent[nb] != -2:
                    continue
                parent[nb] = c
                if nb == target:
                    path = [target]
                    cur = c
                    while cur >= 0:
                        path.append(cur)
                        cur = int(parent[cur])
                    return np.asarray(path, dtype=np.int64)
                queue.append(nb)
        return np.asarray([], dtype=np.int64)

    def _crust_domain_field(self, grid, ctype, origin, stability):
        cont = ctype == CONT
        domain = np.full(grid.n, DOMAIN_OCEANIC, dtype=np.float64)
        if not cont.any():
            return domain
        width = _graph_width_steps(grid, cont)
        domain[cont] = DOMAIN_CONTINENTAL_MARGIN
        domain[cont & (width >= 3.0)] = DOMAIN_CONTINENTAL_INTERIOR
        domain[cont & (origin == ORIGIN_SUTURE)] = DOMAIN_SUTURE
        domain[cont & (origin == ORIGIN_ARC) & (width <= 3.0)] = DOMAIN_ACCRETED_TERRANE
        domain[cont & (origin == ORIGIN_PLUME_IMPACT)] = DOMAIN_LIP
        craton = cont & (width >= 3.0) & (
            (origin == ORIGIN_CRATON) | (stability > 0.78)
        )
        domain[craton] = DOMAIN_CRATON
        return domain

    def _mature_continental_provinces(self, grid, ctype, thick, age, origin,
                                      reworked, stability, orog_age, volc_age,
                                      t, dt):
        """Let quiet sutures and accreted arcs become continental interior.

        A collision belt should not remain an active suture domain forever.  If
        it is broad, old, and has not been reworked recently, keep the
        orogeny-age archive but let the crustal origin mature back into
        continental interior.  Thickened old belts also relax slowly toward a
        post-orogenic crustal thickness instead of staying pinned at the
        collision cap.
        """
        cont = ctype == CONT
        if not cont.any():
            return origin, thick
        width = _graph_width_steps(grid, cont)
        quiet_myr = min(650.0, max(280.0, 0.16 * float(t)))
        quiet_reworked = (reworked < 0.0) | (reworked <= float(t) - quiet_myr)
        quiet_volcanic = (volc_age < 0.0) | (volc_age <= float(t) - quiet_myr)
        old_orogen = (orog_age >= 0.0) & (orog_age <= float(t) - quiet_myr)
        broad_interior = cont & (width >= 4.0)

        mature_suture = (
            broad_interior
            & (origin == ORIGIN_SUTURE)
            & old_orogen
            & quiet_reworked
        )
        mature_arc = (
            broad_interior
            & (origin == ORIGIN_ARC)
            & quiet_reworked
            & quiet_volcanic
            & (stability >= 0.28)
        )
        matured = mature_suture | mature_arc
        if matured.any():
            origin[matured] = ORIGIN_PRIMORDIAL

        relax_mask = (
            broad_interior
            & quiet_reworked
            & (
                old_orogen
                | mature_suture
                | mature_arc
                | ((origin == ORIGIN_CRATON) & (thick > CONT_THICK + 16000.0))
            )
        )
        if relax_mask.any():
            target = np.full(grid.n, CONT_THICK + 9000.0, dtype=np.float64)
            target[mature_arc] = CONT_THICK + 2500.0
            target[origin == ORIGIN_CRATON] = CONT_THICK + 6500.0
            rate = float(np.clip(dt / 420.0, 0.02, 0.10))
            excess = np.maximum(thick - target, 0.0)
            thick[relax_mask] -= excess[relax_mask] * rate
        return origin, thick

    def _consolidate_parented_margin_slivers(self, world, grid, ctype, thick, age,
                                             origin, reworked, stability,
                                             orog_age, volc_age, t, dt,
                                             active_boundary_mask=None):
        """Mature quiet parented accretion/suture strips into continental margin.

        Real continents retain the archive of accreted terranes and sutures, but
        old strips welded to a broad parent should not remain forever classified
        as active accretionary ribbons.  This pass changes crustal state only:
        it does not add or remove land/ocean cells, and it avoids current
        boundary belts plus recently reworked volcanic or collisional crust.
        """
        del dt
        world.set_g("tectonics.last_p33_parented_sliver_candidate_cells", 0.0)
        world.set_g("tectonics.last_p33_parented_sliver_consolidated_cells", 0.0)
        world.set_g("tectonics.last_p33_parented_sliver_area_fraction", 0.0)
        world.set_g("tectonics.last_p33_parented_sliver_active_preserved_cells", 0.0)
        world.set_g("tectonics.last_p33_parented_sliver_detached_preserved_cells", 0.0)

        mature_start_myr = max(3000.0, 0.70 * float(world.spec.t_end_myr))
        if float(t) < mature_start_myr:
            return origin, thick, age, stability

        cont = ctype == CONT
        if not cont.any():
            return origin, thick, age, stability

        width = _graph_width_steps(grid, cont)
        quiet_myr = min(560.0, max(260.0, 0.15 * float(t)))
        old_rework = (reworked < 0.0) | (reworked <= float(t) - quiet_myr)
        old_orogen = (orog_age < 0.0) | (orog_age <= float(t) - quiet_myr)
        old_volcanic = (volc_age < 0.0) | (volc_age <= float(t) - quiet_myr)
        quiet = old_rework & old_orogen & old_volcanic

        active = np.zeros(grid.n, dtype=bool)
        if active_boundary_mask is not None:
            mask = np.asarray(active_boundary_mask, dtype=bool)
            if mask.shape == (grid.n,) and mask.any():
                active = self._dilate_mask(
                    grid, mask, allowed=np.ones(grid.n, dtype=bool), passes=2)

        parent_anchor = (
            cont
            & (
                (origin == ORIGIN_CRATON)
                | ((origin == ORIGIN_PRIMORDIAL) & (stability >= 0.42))
                | (stability >= 0.66)
            )
        )
        if not parent_anchor.any():
            return origin, thick, age, stability

        parent_neighbors = _same_neighbor_count(grid, parent_anchor)
        narrow_accretionary = (
            cont
            & (width <= 2.0)
            & np.isin(origin, [ORIGIN_ARC, ORIGIN_SUTURE])
        )
        parented = narrow_accretionary & (parent_neighbors >= 1)
        candidate = parented & quiet
        world.set_g(
            "tectonics.last_p33_parented_sliver_candidate_cells",
            float(np.count_nonzero(candidate)),
        )
        world.set_g(
            "tectonics.last_p33_parented_sliver_active_preserved_cells",
            float(np.count_nonzero(parented & ~quiet)),
        )
        world.set_g(
            "tectonics.last_p33_parented_sliver_detached_preserved_cells",
            float(np.count_nonzero(narrow_accretionary & (parent_neighbors < 1))),
        )
        candidate &= ~active
        if not candidate.any():
            return origin, thick, age, stability

        max_area = (
            float(grid.cell_area.sum())
            * min(0.0025, 0.00075 * max(float(t), 1.0) / 4500.0)
        )
        cells = np.where(candidate)[0]
        cells = sorted((int(c) for c in cells), key=lambda c: (
            0 if origin[c] == ORIGIN_SUTURE else 1,
            -int(parent_neighbors[c]),
            float(width[c]),
            float(stability[c]),
            c,
        ))
        chosen: list[int] = []
        acc_area = 0.0
        for c in cells:
            if acc_area >= max_area:
                break
            chosen.append(c)
            acc_area += float(grid.cell_area[c])
        if not chosen:
            return origin, thick, age, stability

        chosen_arr = np.asarray(chosen, dtype=np.int64)
        for c in chosen:
            parents = [
                int(nb) for nb in grid.neighbors[c]
                if parent_anchor[int(nb)]
            ]
            if not parents:
                continue
            parents.sort(key=lambda p: (
                0 if origin[p] == ORIGIN_CRATON else 1,
                -float(stability[p]),
                -float(age[p]),
                -float(thick[p]),
                p,
            ))
            p = parents[0]
            origin[c] = ORIGIN_PRIMORDIAL
            thick[c] = max(thick[c], CONT_THICK, min(float(thick[p]), 52000.0))
            age[c] = max(age[c], min(float(age[p]), float(t)))
            stability[c] = max(stability[c], min(float(stability[p]) * 0.90, 0.72), 0.48)

        world.set_g(
            "tectonics.last_p33_parented_sliver_consolidated_cells",
            float(chosen_arr.size),
        )
        world.set_g(
            "tectonics.last_p33_parented_sliver_area_fraction",
            float(grid.cell_area[chosen_arr].sum() / max(float(grid.cell_area.sum()), 1.0)),
        )
        return origin, thick, age, stability

    def _mature_parent_supported_collage_belts(
        self, world, grid, ctype, thick, age, origin, reworked, stability,
        orog_age, volc_age, t, dt, active_boundary_mask=None,
    ):
        """Normalize quiet accreted/suture collage belts into continental platform.

        This is a state maturation pass, not area conservation.  It preserves the
        land/ocean mask and only retags old, quiet, parent-supported collage so
        welded continental margins do not remain active terrane/suture ribbons
        indefinitely.
        """
        world.set_g("tectonics.last_p46_mature_collage_candidate_cells", 0.0)
        world.set_g("tectonics.last_p46_mature_collage_consolidated_cells", 0.0)
        world.set_g("tectonics.last_p46_mature_collage_area_fraction", 0.0)
        world.set_g("tectonics.last_p46_mature_collage_active_preserved_cells", 0.0)
        world.set_g("tectonics.last_p46_mature_collage_detached_preserved_cells", 0.0)
        # This is now part of the default mature-continent path.  P45 localizes
        # active process provenance; without this complementary maturation pass,
        # old welded accretionary collage can stay semantically active and starve
        # later parented continental recovery at production resolution.
        if world.g("tectonics.enable_mature_parent_supported_collage", 1.0) <= 0.0:
            return origin, thick, age, stability

        mature_start_myr = max(3000.0, 0.70 * float(world.spec.t_end_myr))
        if float(t) < mature_start_myr:
            return origin, thick, age, stability

        cont = ctype == CONT
        if not cont.any():
            return origin, thick, age, stability

        width = _graph_width_steps(grid, cont)
        accretionary = cont & np.isin(origin, [ORIGIN_ARC, ORIGIN_SUTURE])
        if not accretionary.any():
            return origin, thick, age, stability

        quiet_myr = min(760.0, max(320.0, 0.16 * float(t)))
        old_rework = (reworked < 0.0) | (reworked <= float(t) - quiet_myr)
        old_orogen = (orog_age < 0.0) | (orog_age <= float(t) - quiet_myr)
        old_volcanic = (volc_age < 0.0) | (volc_age <= float(t) - quiet_myr)
        quiet = old_rework & old_orogen & old_volcanic

        active = np.zeros(grid.n, dtype=bool)
        if active_boundary_mask is not None:
            mask = np.asarray(active_boundary_mask, dtype=bool)
            if mask.shape == (grid.n,) and mask.any():
                active = self._dilate_mask(
                    grid, mask, allowed=np.ones(grid.n, dtype=bool), passes=2)

        parent_anchor = (
            cont
            & (width >= 3.0)
            & (
                (origin == ORIGIN_CRATON)
                | ((origin == ORIGIN_PRIMORDIAL) & (stability >= 0.40))
                | (stability >= 0.66)
            )
        )
        if not parent_anchor.any():
            world.set_g(
                "tectonics.last_p46_mature_collage_detached_preserved_cells",
                float(np.count_nonzero(accretionary & quiet)),
            )
            return origin, thick, age, stability

        parent_reach = self._dilate_mask(
            grid, parent_anchor, allowed=cont, passes=4)
        parent_neighbors = _same_neighbor_count(grid, parent_anchor)
        cont_neighbors = _same_neighbor_count(grid, cont)
        supported = (
            accretionary
            & parent_reach
            & (width <= 3.0)
            & (cont_neighbors >= 2)
            & (age >= min(1200.0, max(520.0, 0.22 * float(t))))
            & (thick >= CONT_THICK)
        )
        candidate = supported & quiet & ~active
        world.set_g(
            "tectonics.last_p46_mature_collage_candidate_cells",
            float(np.count_nonzero(candidate)),
        )
        world.set_g(
            "tectonics.last_p46_mature_collage_active_preserved_cells",
            float(np.count_nonzero(supported & ~quiet | supported & active)),
        )
        world.set_g(
            "tectonics.last_p46_mature_collage_detached_preserved_cells",
            float(np.count_nonzero(accretionary & quiet & ~parent_reach)),
        )
        if not candidate.any():
            return origin, thick, age, stability

        total = float(grid.cell_area.sum())
        step_budget = total * min(
            0.010,
            max(0.0025, 0.0045 * max(float(dt), 1.0) / 20.0),
        )
        cells = sorted((int(c) for c in np.where(candidate)[0]), key=lambda c: (
            0 if origin[c] == ORIGIN_SUTURE else 1,
            0 if width[c] >= 2.0 else 1,
            -int(parent_neighbors[c]),
            -int(cont_neighbors[c]),
            -float(age[c]),
            -float(stability[c]),
            c,
        ))
        chosen: list[int] = []
        acc_area = 0.0
        for c in cells:
            if acc_area >= step_budget:
                break
            chosen.append(c)
            acc_area += float(grid.cell_area[c])
        if not chosen:
            return origin, thick, age, stability

        chosen_arr = np.asarray(chosen, dtype=np.int64)
        was_suture = origin[chosen_arr] == ORIGIN_SUTURE
        origin[chosen_arr] = ORIGIN_PRIMORDIAL
        stability[chosen_arr] = np.maximum(stability[chosen_arr], 0.46)
        target = np.where(
            was_suture,
            CONT_THICK + 7200.0,
            CONT_THICK + 3200.0,
        )
        excess = np.maximum(thick[chosen_arr] - target, 0.0)
        thick[chosen_arr] -= excess * float(np.clip(dt / 360.0, 0.03, 0.12))
        age[chosen_arr] = np.maximum(age[chosen_arr], min(float(t), 1200.0))

        world.set_g(
            "tectonics.last_p46_mature_collage_consolidated_cells",
            float(chosen_arr.size),
        )
        world.set_g(
            "tectonics.last_p46_mature_collage_area_fraction",
            float(grid.cell_area[chosen_arr].sum() / max(total, 1.0)),
        )
        return origin, thick, age, stability

    def _maintain_cratonic_kernels(self, grid, ctype, age, thick, origin,
                                   reworked, stability, t, dt, blocked=None,
                                   world=None):
        """Preserve deterministic old stable continental kernels.

        Continental lithosphere should retain at least a small ancient stable
        core through repeated plate reorganization.  When rasterized advection
        and collision bookkeeping erase all high-stability cells, promote the
        oldest broad interior cells back to cratonic origin gradually, using
        only age, width, reworking age, and thickness as process proxies.
        """
        if world is not None:
            world.set_g("tectonics.last_p47_craton_target_fraction", 0.0)
            world.set_g("tectonics.last_p47_pre_stable_craton_fraction", 0.0)
            world.set_g("tectonics.last_p47_promoted_craton_cells", 0.0)
            world.set_g("tectonics.last_p47_promoted_craton_area_fraction", 0.0)
        cont = ctype == CONT
        if not cont.any() or t < 900.0:
            return origin, age, thick, stability
        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1e-12)
        width = _graph_width_steps(grid, cont)
        stable = cont & (width >= 3.0) & (stability > 0.75)
        blocked_mask = np.zeros(grid.n, dtype=bool)
        if blocked is not None:
            candidate_block = np.asarray(blocked, dtype=bool)
            if candidate_block.shape == (grid.n,) and candidate_block.any():
                blocked_mask = self._dilate_mask(
                    grid, candidate_block, allowed=np.ones(grid.n, dtype=bool),
                    passes=1)
        target_fraction = self._cratonic_kernel_target_fraction(t)
        target_area = target_fraction * cont_area
        pre_stable_fraction = float(area[stable].sum() / cont_area)
        if world is not None:
            world.set_g("tectonics.last_p47_craton_target_fraction", target_fraction)
            world.set_g(
                "tectonics.last_p47_pre_stable_craton_fraction",
                pre_stable_fraction,
            )
        deficit = target_area - float(area[stable].sum())
        if deficit <= 0.0:
            return origin, age, thick, stability

        quiet_bonus = np.where(reworked < 0.0, float(t), np.maximum(float(t) - reworked, 0.0))
        candidates = np.where(cont & (width >= 3.0) & ~stable & ~blocked_mask)[0]
        if candidates.size == 0:
            return origin, age, thick, stability

        step_budget = min(deficit, 0.012 * cont_area * max(dt, 1.0) / 20.0)
        ranked = sorted((int(c) for c in candidates), key=lambda c: (
            -float(age[c]),
            -float(width[c]),
            -float(quiet_bonus[c]),
            -float(thick[c]),
            c,
        ))
        promoted: list[int] = []
        acc = 0.0
        min_age = min(float(t), max(1800.0, 0.62 * float(t)))
        for c in ranked:
            if acc >= step_budget:
                break
            promoted.append(c)
            acc += float(area[c])
        if not promoted:
            return origin, age, thick, stability
        cells = np.asarray(promoted, dtype=int)
        origin[cells] = ORIGIN_CRATON
        age[cells] = np.maximum(age[cells], min_age)
        thick[cells] = np.maximum(thick[cells], CONT_THICK + 4500.0)
        stability[cells] = np.maximum(stability[cells], 0.82)
        if world is not None:
            world.set_g("tectonics.last_p47_promoted_craton_cells", float(cells.size))
            world.set_g(
                "tectonics.last_p47_promoted_craton_area_fraction",
                float(area[cells].sum() / max(cont_area, 1e-12)),
            )
        return origin, age, thick, stability

    def _cratonic_kernel_target_fraction(self, t):
        tt = max(float(t), 0.0)
        if tt < 900.0:
            return 0.0
        maturity = float(np.clip((tt - 900.0) / 2400.0, 0.0, 1.0))
        # Mature Earth-like continents should retain broad shield cores without
        # turning every quiet platform into high-stability craton.  The target is
        # a fraction of continental crust area, reached gradually by the
        # per-step budget.  Older platform aprons are handled separately from
        # this hard craton-core target.
        return float(0.035 + 0.045 * maturity)

    def _continent_objects(self, grid, continent_id, ctype, age, thick,
                           origin, stability, t):
        cont = ctype == CONT
        ids = [int(x) for x in np.unique(continent_id[cont].astype(int)) if x >= 0]
        area = grid.cell_area
        total = float(area.sum())
        width = _graph_width_steps(grid, cont)
        out = []
        for cid in ids:
            mask = cont & (continent_id.astype(int) == cid)
            cells = np.where(mask)[0].astype(np.int64)
            if cells.size == 0:
                continue
            cell_area = area[cells]
            comp_area = float(cell_area.sum())
            centroid = np.average(grid.xyz[cells], axis=0, weights=cell_area)
            norm = float(np.linalg.norm(centroid))
            if norm > 1e-12:
                centroid = centroid / norm
            lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
            lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
            craton = mask & ((origin == ORIGIN_CRATON) | (stability > 0.78))
            margin = mask & (width <= 2.0)
            out.append({
                "id": f"continent:{cid}",
                "numeric_id": int(cid),
                "kind": "continent",
                "cell_count": int(cells.size),
                "area_fraction": comp_area / total if total else 0.0,
                "lat": round(lat, 2),
                "lon": round(lon, 2),
                "width_p50_steps": float(np.percentile(width[cells], 50)),
                "width_p90_steps": float(np.percentile(width[cells], 90)),
                "core_fraction": (
                    float(area[craton].sum() / comp_area) if comp_area > 0.0 else 0.0
                ),
                "margin_fraction": (
                    float(area[margin].sum() / comp_area) if comp_area > 0.0 else 0.0
                ),
                "mean_age_myr": float(np.average(age[cells], weights=cell_area)),
                "mean_thickness_m": float(np.average(thick[cells], weights=cell_area)),
                "mean_stability": float(np.average(stability[cells], weights=cell_area)),
                "last_active_myr": round(float(t), 1),
                "cells": cells.astype(int).tolist() if cells.size <= 600 else [],
            })
        out.sort(key=lambda obj: obj["area_fraction"], reverse=True)
        return out

    def _terrane_mask(self, grid, ctype, origin, stability):
        cont = ctype == CONT
        if not cont.any():
            return np.zeros(grid.n, dtype=bool)
        width = _graph_width_steps(grid, cont)
        young_arc = cont & (origin == ORIGIN_ARC) & (stability < 0.55)
        narrow_suture = cont & (origin == ORIGIN_SUTURE) & (width <= 2.0) & (stability < 0.65)
        return young_arc | narrow_suture

    def _terrane_objects(self, grid, terrane_id, continent_id, ctype, origin,
                         age, thick, t):
        terrane = terrane_id >= 0
        ids = [int(x) for x in np.unique(terrane_id[terrane].astype(int)) if x >= 0]
        area = grid.cell_area
        total = float(area.sum())
        out = []
        for tid in ids:
            mask = terrane & (terrane_id.astype(int) == tid)
            cells = np.where(mask)[0].astype(np.int64)
            if cells.size == 0:
                continue
            cell_area = area[cells]
            comp_area = float(cell_area.sum())
            adjacent_ids: list[int] = []
            for c in cells:
                adjacent_ids.extend(
                    int(continent_id[nb])
                    for nb in grid.neighbors[int(c)]
                    if continent_id[nb] >= 0 and not mask[nb]
                )
            parent = -1
            if adjacent_ids:
                vals, counts = np.unique(np.asarray(adjacent_ids, dtype=int), return_counts=True)
                parent = int(vals[int(np.argmax(counts))])
            kind = "accreted_terrane"
            if np.any(origin[cells] == ORIGIN_ARC):
                kind = "island_arc_or_accreted_arc"
            if np.any(origin[cells] == ORIGIN_SUTURE):
                kind = "suture_terrane"
            centroid = np.average(grid.xyz[cells], axis=0, weights=cell_area)
            norm = float(np.linalg.norm(centroid))
            if norm > 1e-12:
                centroid = centroid / norm
            lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
            lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
            out.append({
                "id": f"terrane:{tid}",
                "numeric_id": int(tid),
                "kind": kind,
                "parent_continent_id": int(parent),
                "cell_count": int(cells.size),
                "area_fraction": comp_area / total if total else 0.0,
                "lat": round(lat, 2),
                "lon": round(lon, 2),
                "mean_age_myr": float(np.average(age[cells], weights=cell_area)),
                "mean_thickness_m": float(np.average(thick[cells], weights=cell_area)),
                "last_active_myr": round(float(t), 1),
                "cells": cells.astype(int).tolist() if cells.size <= 400 else [],
            })
        out.sort(key=lambda obj: obj["area_fraction"], reverse=True)
        return out

    def _continent_lifecycle_events(self, grid, prev_continent_id, continent_id,
                                    prev_terrane_id, terrane_id,
                                    previous_terranes, terrane_objects,
                                    reorg_detail, t):
        """Create object/event lineage from id overlap and plate-capture cargo.

        The raster id fields carry current membership; this object stream records
        why identities changed so archives can distinguish genuine continent
        split/merge/capture from plate-label cleanup.
        """
        area = grid.cell_area
        total = float(area.sum())
        prev_continent = np.asarray(prev_continent_id, dtype=np.float64).astype(int)
        curr_continent = np.asarray(continent_id, dtype=np.float64).astype(int)
        objects: list[dict] = []

        def add(event_type, entity_kind, parent_ids, child_ids, cells=None, **extra):
            parent_ids_i = [int(x) for x in parent_ids if int(x) >= 0]
            child_ids_i = [int(x) for x in child_ids if int(x) >= 0]
            cells_arr = np.asarray([] if cells is None else cells, dtype=np.int64)
            area_fraction = (
                float(area[cells_arr].sum() / total)
                if cells_arr.size and total > 0.0 else 0.0
            )
            obj = {
                "id": (
                    f"{event_type}:{entity_kind}:"
                    f"{'-'.join(map(str, parent_ids_i)) or 'none'}>"
                    f"{'-'.join(map(str, child_ids_i)) or 'none'}:"
                    f"{round(float(t), 1)}"
                ),
                "kind": event_type,
                "entity_kind": entity_kind,
                "parent_ids": parent_ids_i,
                "child_ids": child_ids_i,
                "area_fraction": area_fraction,
                "time_myr": round(float(t), 1),
                "cells": cells_arr.astype(int).tolist() if cells_arr.size <= 600 else [],
                **extra,
            }
            objects.append(obj)
            location = int(cells_arr[0]) if cells_arr.size else None
            return Event(
                event_type, t, self.name, location=location,
                magnitude=area_fraction, params={
                    "entity_kind": entity_kind,
                    "parent_ids": parent_ids_i,
                    "child_ids": child_ids_i,
                    **extra,
                },
            )

        events: list[Event] = []
        prev_active = sorted(int(x) for x in np.unique(prev_continent) if x >= 0)
        curr_active = sorted(int(x) for x in np.unique(curr_continent) if x >= 0)
        min_abs = 0.0006 * total
        min_rel = 0.08

        prev_to_curr: dict[int, list[tuple[int, float]]] = {pid: [] for pid in prev_active}
        curr_to_prev: dict[int, list[tuple[int, float]]] = {cid: [] for cid in curr_active}
        for pid in prev_active:
            pmask = prev_continent == pid
            parea = float(area[pmask].sum())
            if parea <= 0.0:
                continue
            for cid in curr_active:
                overlap = pmask & (curr_continent == cid)
                oarea = float(area[overlap].sum())
                if oarea >= min_abs and oarea / parea >= min_rel:
                    prev_to_curr[pid].append((cid, oarea))
        for cid in curr_active:
            cmask = curr_continent == cid
            carea = float(area[cmask].sum())
            if carea <= 0.0:
                continue
            for pid in prev_active:
                overlap = cmask & (prev_continent == pid)
                oarea = float(area[overlap].sum())
                if oarea >= min_abs and oarea / carea >= min_rel:
                    curr_to_prev[cid].append((pid, oarea))

        for pid, children in prev_to_curr.items():
            if len(children) >= 2:
                child_ids = [cid for cid, _ in sorted(children, key=lambda x: (-x[1], x[0]))]
                cells = np.where(prev_continent == pid)[0]
                events.append(add("continent_split", "continent", [pid], child_ids, cells))
        for cid, parents in curr_to_prev.items():
            if len(parents) >= 2:
                parent_ids = [pid for pid, _ in sorted(parents, key=lambda x: (-x[1], x[0]))]
                cells = np.where(curr_continent == cid)[0]
                events.append(add("continent_merge", "continent", parent_ids, [cid], cells))
        for cid in curr_active:
            if not curr_to_prev.get(cid):
                cells = np.where(curr_continent == cid)[0]
                if float(area[cells].sum()) >= min_abs:
                    events.append(add("continent_birth", "continent", [], [cid], cells))
        for pid in prev_active:
            if not prev_to_curr.get(pid):
                cells = np.where(prev_continent == pid)[0]
                if float(area[cells].sum()) >= min_abs:
                    events.append(add("continent_loss", "continent", [pid], [], cells))

        prev_terrane_by_id = {
            int(obj.get("numeric_id", -1)): obj
            for obj in previous_terranes or []
        }
        for obj in terrane_objects or []:
            tid = int(obj.get("numeric_id", -1))
            if tid < 0:
                continue
            prev = prev_terrane_by_id.get(tid)
            parent = int(obj.get("parent_continent_id", -1))
            prev_parent = int(prev.get("parent_continent_id", -1)) if prev else -1
            if parent >= 0 and prev_parent >= 0 and parent != prev_parent:
                cells = np.asarray(obj.get("cells", []), dtype=np.int64)
                events.append(add(
                    "terrane_capture", "terrane", [prev_parent], [parent], cells,
                    terrane_id=tid,
                    previous_parent_continent_id=prev_parent,
                    parent_continent_id=parent,
                ))

        for merge in (reorg_detail or {}).get("merged", []):
            if merge.get("cargo_policy") != "capture_plate_label_preserve_crust_cargo":
                continue
            if float(merge.get("source_continental_fraction", 0.0)) < 0.18:
                continue
            events.append(add(
                "microcontinent_plate_capture",
                "plate_label",
                [int(merge.get("from", -1))],
                [int(merge.get("to", -1))],
                [],
                source_plate_id=int(merge.get("from", -1)),
                target_plate_id=int(merge.get("to", -1)),
                source_continental_fraction=float(merge.get("source_continental_fraction", 0.0)),
                cargo_policy=merge.get("cargo_policy"),
            ))

        return objects, events

    def _plate_topology_objects(self, grid, plate, plates, ctype, continent_id,
                                terrane_id, boundaries, t):
        """Summarize resolved plate topology for downstream object logic.

        This mirrors the role of a resolved plate topology object in a
        reconstruction tool: the raster remains the render target, but plate
        decisions can inspect persistent plate areas, neighbours, boundary
        composition, and continental cargo from one object record.
        """
        plate = np.asarray(plate, dtype=int)
        ctype = np.asarray(ctype, dtype=np.float64)
        continent_id = np.asarray(continent_id, dtype=np.float64)
        terrane_id = np.asarray(terrane_id, dtype=np.float64)
        total = float(grid.cell_area.sum())
        ids = sorted(int(x) for x in np.unique(plate) if x >= 0)

        boundary_masks: dict[str, np.ndarray] = {}
        for kind in (
            "ridge",
            "trench",
            "collision",
            "suture",
            "active_margin",
            "passive_margin",
            "transform",
            "divergent",
            "convergent",
        ):
            mask = np.zeros(grid.n, dtype=bool)
            if isinstance(boundaries, dict):
                cells = np.asarray(boundaries.get(kind, []), dtype=int)
                cells = cells[(0 <= cells) & (cells < grid.n)]
                mask[cells] = True
            boundary_masks[kind] = mask

        neighbours: dict[int, set[int]] = {pid: set() for pid in ids}
        shared_edge_counts: dict[int, dict[int, int]] = {pid: {} for pid in ids}
        for i, j in grid.edges:
            pi = int(plate[int(i)])
            pj = int(plate[int(j)])
            if pi == pj or pi < 0 or pj < 0:
                continue
            neighbours.setdefault(pi, set()).add(pj)
            neighbours.setdefault(pj, set()).add(pi)
            shared_edge_counts.setdefault(pi, {})[pj] = (
                shared_edge_counts.setdefault(pi, {}).get(pj, 0) + 1
            )
            shared_edge_counts.setdefault(pj, {})[pi] = (
                shared_edge_counts.setdefault(pj, {}).get(pi, 0) + 1
            )

        plate_by_id = {int(obj.get("id", idx)): obj for idx, obj in enumerate(plates)}
        out = []
        for pid in ids:
            mask = plate == pid
            cells = np.where(mask)[0].astype(np.int64)
            if cells.size == 0:
                continue
            cell_area = grid.cell_area[cells]
            plate_area = float(cell_area.sum())
            comps = self._connected_components(grid, mask)
            largest = max((float(grid.cell_area[comp].sum()) for comp in comps), default=0.0)
            centroid = np.average(grid.xyz[cells], axis=0, weights=cell_area)
            centroid = centroid / max(float(np.linalg.norm(centroid)), 1e-12)
            lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
            lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
            cont = mask & (ctype == CONT)
            ocean = mask & (ctype == OCEAN)
            cargo_continent_ids = sorted(
                int(x) for x in np.unique(continent_id[mask].astype(int)) if x >= 0
            )
            cargo_terrane_ids = sorted(
                int(x) for x in np.unique(terrane_id[mask].astype(int)) if x >= 0
            )
            boundary_area_fractions = {
                kind: (
                    float(grid.cell_area[mask & bmask].sum() / plate_area)
                    if plate_area > 0.0 else 0.0
                )
                for kind, bmask in boundary_masks.items()
            }
            shared = {
                int(other): int(count)
                for other, count in sorted(shared_edge_counts.get(pid, {}).items())
            }
            motion = plate_by_id.get(pid, {})
            out.append({
                "id": f"plate_topology:{pid}",
                "numeric_id": int(pid),
                "kind": "resolved_plate_topology",
                "cell_count": int(cells.size),
                "area_fraction": plate_area / total if total else 0.0,
                "lat": round(lat, 2),
                "lon": round(lon, 2),
                "component_count": int(len(comps)),
                "largest_component_fraction": (
                    float(np.clip(largest / plate_area, 0.0, 1.0))
                    if plate_area > 0.0 else 0.0
                ),
                "neighbour_plate_ids": sorted(int(x) for x in neighbours.get(pid, set())),
                "shared_edge_counts": shared,
                "continental_fraction": (
                    float(grid.cell_area[cont].sum() / plate_area) if plate_area > 0.0 else 0.0
                ),
                "oceanic_fraction": (
                    float(grid.cell_area[ocean].sum() / plate_area) if plate_area > 0.0 else 0.0
                ),
                "continent_ids": cargo_continent_ids,
                "terrane_ids": cargo_terrane_ids,
                "boundary_area_fractions": boundary_area_fractions,
                "dominant_boundary_kind": max(
                    boundary_area_fractions,
                    key=lambda kind: boundary_area_fractions[kind],
                ),
                "pole": motion.get("pole"),
                "rate": float(motion.get("rate", 0.0)),
                "motion_source": motion.get("motion_source", "initial_or_unresolved"),
                "r2_force_components": motion.get("r2_force_components", {}),
                "last_resolved_myr": round(float(t), 1),
                "cells": cells.astype(int).tolist() if cells.size <= 500 else [],
            })
        out.sort(key=lambda obj: obj["area_fraction"], reverse=True)
        return out

    def _plate_rank_fields_and_objects(self, grid, plate, plate_topologies,
                                       ctype, origin, age, boundaries, t,
                                       *, protected_records=None,
                                       enable_protection=False):
        plate = np.asarray(plate, dtype=int)
        total = max(float(grid.cell_area.sum()), 1.0e-12)
        rank = np.full(grid.n, -1.0, dtype=np.float64)
        protected_field = np.full(grid.n, -1.0, dtype=np.float64)
        topology_by_id = {
            int(obj.get("numeric_id", -1)): obj
            for obj in plate_topologies or []
        }
        protected_by_id = {
            int(record.get("plate_id", -1)): str(record.get("reason", "protected"))
            for record in protected_records or []
        }
        context = {
            "ctype": ctype,
            "origin": origin,
            "age": age,
            "boundaries": boundaries,
            "plate_topologies": plate_topologies or [],
            "enable_p107_ranked_plate_policy": bool(enable_protection),
        }
        objects: list[dict] = []
        for pid in sorted(int(x) for x in np.unique(plate) if int(x) >= 0):
            mask = plate == pid
            area_fraction = float(grid.cell_area[mask].sum() / total)
            if area_fraction >= self.MAJOR_PLATE_AREA_FRAC:
                code = 3.0
                rank_name = "major"
            elif area_fraction >= self.MIN_PLATE_AREA_FRAC:
                code = 2.0
                rank_name = "minor"
            elif area_fraction >= self.MICRO_PLATE_AREA_FRAC:
                code = 1.0
                rank_name = "micro"
            else:
                code = 0.0
                rank_name = "tiny"
            rank[mask] = code
            reason = protected_by_id.get(pid, "")
            if not reason and enable_protection:
                reason = self._p107_protected_microplate_reason(
                    grid, plate, pid, area_fraction, context, topology_by_id.get(pid))
            if reason:
                protected_field[mask] = float(pid)
            if rank_name == "micro" or reason:
                cells = np.where(mask)[0].astype(np.int64)
                centroid = np.average(grid.xyz[cells], axis=0, weights=grid.cell_area[cells])
                centroid = centroid / max(float(np.linalg.norm(centroid)), 1.0e-12)
                lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
                lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
                topology = topology_by_id.get(pid, {})
                objects.append({
                    "id": f"microplate:{pid}:{round(float(t), 1)}",
                    "kind": "protected_microplate" if reason else "microplate",
                    "plate_id": int(pid),
                    "rank": rank_name,
                    "area_fraction": area_fraction,
                    "cell_count": int(mask.sum()),
                    "protected": bool(reason),
                    "protection_reason": reason,
                    "lat": round(lat, 2),
                    "lon": round(lon, 2),
                    "continent_ids": list(topology.get("continent_ids", []))
                    if isinstance(topology, dict) else [],
                    "terrane_ids": list(topology.get("terrane_ids", []))
                    if isinstance(topology, dict) else [],
                    "boundary_area_fractions": dict(topology.get("boundary_area_fractions", {}))
                    if isinstance(topology, dict) else {},
                    "last_resolved_myr": round(float(t), 1),
                    "cells": cells.astype(int).tolist() if cells.size <= 500 else [],
                })
        objects.sort(key=lambda obj: (-float(obj["area_fraction"]), int(obj["plate_id"])))
        return rank, protected_field, objects

    def _plate_accretion_event_objects(self, reorg_detail, t):
        detail = reorg_detail if isinstance(reorg_detail, dict) else {}
        out: list[dict] = []
        for idx, merge in enumerate(detail.get("merged", []) or []):
            source = int(merge.get("from", -1))
            target = int(merge.get("to", -1))
            out.append({
                "id": f"plate_accretion:{source}:{target}:{round(float(t), 1)}:{idx}",
                "kind": "plate_label_accretion",
                "source_plate_id": source,
                "target_plate_id": target,
                "area_fraction": float(merge.get("area_fraction", 0.0)),
                "basis": str(merge.get("basis", "")),
                "cargo_policy": str(merge.get("cargo_policy", "")),
                "source_continental_fraction": float(
                    merge.get("source_continental_fraction", 0.0)),
                "time_myr": round(float(t), 1),
            })
        return out[-96:]

    def _boundary_objects(self, grid, boundaries, plate, t, ctype=None, age=None,
                          thick=None, velocity=None, previous_objects=None):
        objects = []
        kind_stage = {
            "divergent": "continental_rift",
            "ridge": "young_ocean",
            "trench": "subduction",
            "convergent_parent": "ordered_convergent_parent",
            "subduction_parent": "ordered_subduction_parent",
            "suture": "collision_suture",
            "passive_margin": "passive_margin",
            "active_margin": "active_margin",
            "transform": "transform",
        }
        previous_by_kind: dict[str, list[dict]] = {}
        for obj in previous_objects or []:
            previous_by_kind.setdefault(str(obj.get("kind", "")), []).append(obj)
        used_previous: set[str] = set()
        for kind, stage in kind_stage.items():
            cells = np.asarray(boundaries.get(kind, []), dtype=np.int64)
            if cells.size == 0:
                continue
            if kind == "divergent":
                ridge_cells = np.asarray(boundaries.get("ridge", []), dtype=np.int64)
                if ridge_cells.size:
                    ridge_mask = np.zeros(grid.n, dtype=bool)
                    ridge_mask[ridge_cells] = True
                    cells = cells[~ridge_mask[cells]]
                    if cells.size == 0:
                        continue
            mask = np.zeros(grid.n, dtype=bool)
            mask[cells] = True
            comps = self._connected_components(grid, mask)
            comps.sort(key=lambda c: float(grid.cell_area[c].sum()), reverse=True)
            if kind in {"ridge", "transform"}:
                comps = self._boundary_lifecycle_components(
                    grid, comps, plate, kind=kind)
            for n, comp in enumerate(comps[:12]):
                if comp.size < 3:
                    continue
                pids = self._boundary_adjacent_plate_ids(grid, plate, comp)
                centroid = grid.xyz[comp].mean(axis=0)
                centroid = centroid / max(np.linalg.norm(centroid), 1e-12)
                lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
                lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
                matched = self._match_previous_boundary_object(
                    comp, previous_by_kind.get(kind, []), used_previous)
                if matched is not None:
                    obj_id = str(matched.get("id"))
                    birth = float(matched.get("birth_myr", t))
                    used_previous.add(obj_id)
                    persistence = "matched_overlap"
                else:
                    obj_id = f"{kind}:{int(comp.min())}:{n}:{round(float(t), 1)}"
                    birth = float(t)
                    persistence = "new"
                context = self._boundary_context(
                    grid, comp, plate, pids, ctype=ctype, age=age, thick=thick,
                    velocity=velocity, kind=kind)
                objects.append({
                    "id": obj_id,
                    "kind": kind,
                    "stage": stage,
                    "cells": comp.astype(int).tolist(),
                    "cell_count": int(comp.size),
                    "area_fraction": float(grid.cell_area[comp].sum() / grid.cell_area.sum()),
                    "birth_myr": round(float(birth), 1),
                    "last_active_myr": round(float(t), 1),
                    "age_myr": round(float(max(t - birth, 0.0)), 1),
                    "persistence": persistence,
                    "parent_plate_ids": pids,
                    "lat": round(lat, 2),
                    "lon": round(lon, 2),
                    **context,
                })
        return objects

    def _boundary_polyline_objects(self, grid, boundaries, boundary_objects,
                                   plate, t):
        """Build ordered polyline objects for major tectonic boundary classes."""
        kind_order = (
            "ridge",
            "transform",
            "trench",
            "suture",
            "active_margin",
            "passive_margin",
            "convergent_parent",
            "subduction_parent",
        )
        kind_to_polyline = {
            "ridge": "ridge_polyline",
            "transform": "transform_polyline",
            "trench": "trench_polyline",
            "suture": "suture_polyline",
            "active_margin": "active_margin_polyline",
            "passive_margin": "passive_margin_polyline",
            "convergent_parent": "convergent_parent_polyline",
            "subduction_parent": "subduction_parent_polyline",
        }
        kind_to_role = {
            "ridge": "divergent_ridge_axis",
            "transform": "transform_fault_axis",
            "trench": "subduction_trench_axis",
            "suture": "collision_suture_axis",
            "active_margin": "active_margin_axis",
            "passive_margin": "passive_margin_axis",
            "convergent_parent": "convergent_parent_axis",
            "subduction_parent": "subduction_parent_axis",
        }
        area = np.asarray(grid.cell_area, dtype=np.float64)
        total_area = max(float(area.sum()), 1.0e-12)
        xyz = np.asarray(getattr(grid, "xyz", np.zeros((0, 3))), dtype=np.float64)
        xyz_valid = xyz.shape == (grid.n, 3)
        source_by_kind: dict[str, list[dict]] = {kind: [] for kind in kind_order}
        for obj in boundary_objects or []:
            kind = str(obj.get("kind", ""))
            if kind not in source_by_kind:
                continue
            source_by_kind[kind].append(obj)
        for kind in kind_order:
            if source_by_kind[kind]:
                continue
            cells = np.asarray(boundaries.get(kind, []), dtype=np.int64)
            cells = cells[(0 <= cells) & (cells < grid.n)]
            if cells.size == 0:
                continue
            source_by_kind[kind].append({
                "id": f"{kind}:network:{round(float(t), 1)}",
                "kind": kind,
                "cells": cells.astype(int).tolist(),
                "parent_plate_ids": [],
                "stage": "",
                "birth_myr": round(float(t), 1),
            })

        def clean_cells(cells) -> np.ndarray:
            out = np.asarray(cells, dtype=np.int64)
            out = out[(0 <= out) & (out < grid.n)]
            return np.unique(out)

        def step_distance(left: int, right: int) -> float:
            if xyz_valid:
                value = float(np.linalg.norm(xyz[int(left)] - xyz[int(right)]))
                if np.isfinite(value) and value > 0.0:
                    return value
            return 1.0

        def bfs(start: int, allowed: set[int]) -> tuple[dict[int, int], dict[int, int]]:
            parent: dict[int, int] = {int(start): -1}
            dist: dict[int, int] = {int(start): 0}
            queue = [int(start)]
            head = 0
            while head < len(queue):
                cell = queue[head]
                head += 1
                for nb_raw in grid.neighbors[cell]:
                    nb = int(nb_raw)
                    if nb not in allowed or nb in dist:
                        continue
                    dist[nb] = dist[cell] + 1
                    parent[nb] = cell
                    queue.append(nb)
            return dist, parent

        def farthest(dist: dict[int, int], candidates: list[int]) -> int:
            reachable = [int(cell) for cell in candidates if int(cell) in dist]
            if not reachable:
                reachable = list(dist.keys())
            return max(reachable, key=lambda cell: (dist.get(cell, -1), -cell))

        def ordered_component_path(cells: np.ndarray) -> tuple[list[int], list[int], int, int, int]:
            cells = clean_cells(cells)
            if cells.size == 0:
                return [], [], 0, 0, 0
            if cells.size == 1:
                only = int(cells[0])
                return [only], [only, only], 0, 0, 0
            mask = np.zeros(grid.n, dtype=bool)
            mask[cells] = True
            degree = np.zeros(grid.n, dtype=np.int16)
            for cell in cells:
                degree[int(cell)] = int(np.count_nonzero(mask[grid.neighbors[int(cell)]]))
            endpoints = [int(cell) for cell in cells if int(degree[int(cell)]) <= 1]
            ordered_cells = [int(cell) for cell in cells]
            allowed = set(ordered_cells)
            seed = endpoints[0] if endpoints else ordered_cells[0]
            dist_seed, _ = bfs(seed, allowed)
            first_candidates = endpoints if len(endpoints) >= 2 else list(dist_seed)
            first = farthest(dist_seed, first_candidates)
            dist_first, parent = bfs(first, allowed)
            second_candidates = endpoints if len(endpoints) >= 2 else list(dist_first)
            second = farthest(dist_first, second_candidates)
            path = [second]
            cursor = second
            while cursor != first and cursor in parent:
                cursor = int(parent[cursor])
                if cursor < 0:
                    break
                path.append(cursor)
            if path[-1] != first:
                path = sorted(ordered_cells)
                first = path[0]
                second = path[-1]
            else:
                path.reverse()
            junction_count = int(np.count_nonzero(mask & (degree >= 3)))
            high_degree_count = int(np.count_nonzero(mask & (degree >= 4)))
            return path, [int(first), int(second)], int(
                len(endpoints)), junction_count, high_degree_count

        def path_metrics(source_cells: np.ndarray, path_cells: list[int]) -> dict[str, float]:
            path_length = 0.0
            for left, right in zip(path_cells, path_cells[1:]):
                path_length += step_distance(left, right)
            direct = (
                step_distance(path_cells[0], path_cells[-1])
                if len(path_cells) >= 2 else 0.0
            )
            if path_length <= 0.0:
                directness = 1.0
                sinuosity = 1.0
            else:
                directness = float(np.clip(direct / path_length, 0.0, 1.0))
                sinuosity = float(path_length / max(direct, 1.0e-9))
            return {
                "path_length_m": float(path_length),
                "direct_distance_m": float(direct),
                "directness": float(directness),
                "sinuosity": float(max(sinuosity, 1.0)),
                "path_coverage_fraction": float(
                    len(path_cells) / max(float(source_cells.size), 1.0)
                ),
            }

        def bridge_short_gaps(kind: str, source_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            source_mask = np.asarray(source_mask, dtype=bool)
            if not source_mask.any():
                return source_mask.copy(), np.zeros(grid.n, dtype=bool)
            policy = {
                "ridge": (2, 0.18, 36),
                "trench": (2, 0.22, 24),
                "convergent_parent": (2, 0.20, 48),
                "subduction_parent": (2, 0.20, 48),
                "transform": (1, 0.12, 24),
                "suture": (1, 0.12, 16),
            }
            if kind not in policy:
                return source_mask.copy(), np.zeros(grid.n, dtype=bool)
            labels, comps = _component_labels(grid, source_mask)
            if len(comps) <= 1:
                return source_mask.copy(), np.zeros(grid.n, dtype=bool)
            max_gap_cells, max_fraction, max_total = policy[kind]
            source_count = int(np.count_nonzero(source_mask))
            max_added = min(int(max_total), max(1, int(np.ceil(
                float(max_fraction) * max(float(source_count), 1.0)))))
            allowed = self._dilate_mask(
                grid,
                source_mask,
                passes=max(1, int(max_gap_cells)),
            )
            component_degree = _same_neighbor_count(grid, source_mask)
            endpoints: list[tuple[int, int]] = []
            for component_index, comp in enumerate(comps):
                comp = np.asarray(comp, dtype=np.int64)
                comp_endpoints = [
                    int(cell) for cell in comp
                    if int(component_degree[int(cell)]) <= 1
                ]
                if not comp_endpoints:
                    comp_endpoints = sorted(int(cell) for cell in comp)[:2]
                for cell in sorted(comp_endpoints)[:4]:
                    endpoints.append((int(cell), int(component_index)))
            if len(endpoints) <= 1:
                return source_mask.copy(), np.zeros(grid.n, dtype=bool)

            parent_root = list(range(len(comps)))

            def find_root(value: int) -> int:
                while parent_root[value] != value:
                    parent_root[value] = parent_root[parent_root[value]]
                    value = parent_root[value]
                return value

            def union(left: int, right: int) -> None:
                root_left = find_root(left)
                root_right = find_root(right)
                if root_left != root_right:
                    parent_root[root_right] = root_left

            def candidate_path(start: int, start_label: int) -> tuple[int, int, list[int], list[int]] | None:
                parent = np.full(grid.n, -2, dtype=np.int64)
                depth = np.zeros(grid.n, dtype=np.int16)
                parent[start] = -1
                queue = [int(start)]
                head = 0
                max_depth = int(max_gap_cells) + 1
                while head < len(queue):
                    cell = queue[head]
                    head += 1
                    if int(depth[cell]) >= max_depth:
                        continue
                    for nb_raw in grid.neighbors[cell]:
                        nb = int(nb_raw)
                        if not allowed[nb] or parent[nb] != -2:
                            continue
                        parent[nb] = cell
                        depth[nb] = depth[cell] + 1
                        nb_label = int(labels[nb])
                        if nb_label >= 0 and nb_label != start_label:
                            path = [nb]
                            cursor = cell
                            while cursor >= 0:
                                path.append(int(cursor))
                                cursor = int(parent[cursor])
                            path.reverse()
                            bridge = [
                                int(path_cell) for path_cell in path
                                if not source_mask[int(path_cell)]
                            ]
                            if 0 < len(bridge) <= int(max_gap_cells):
                                return int(start_label), nb_label, path, bridge
                            return None
                        queue.append(nb)
                return None

            candidates: list[tuple[int, int, int, int, list[int], list[int]]] = []
            for start, label in endpoints:
                found = candidate_path(start, label)
                if found is None:
                    continue
                left_label, right_label, path, bridge = found
                candidates.append((
                    len(bridge),
                    len(path),
                    min(int(path[0]), int(path[-1])),
                    max(int(path[0]), int(path[-1])),
                    [int(cell) for cell in path],
                    [int(cell) for cell in bridge],
                ))
            if not candidates:
                return source_mask.copy(), np.zeros(grid.n, dtype=bool)

            bridged = source_mask.copy()
            bridge_mask = np.zeros(grid.n, dtype=bool)
            added_count = 0
            for _bridge_len, _path_len, _start_key, _end_key, path, bridge in sorted(candidates):
                source_labels = sorted({
                    int(labels[int(cell)]) for cell in path
                    if int(labels[int(cell)]) >= 0
                })
                if len(source_labels) < 2:
                    continue
                roots = {find_root(label) for label in source_labels}
                if len(roots) < 2:
                    continue
                bridge_arr = np.asarray(bridge, dtype=np.int64)
                new_bridge = bridge_arr[~bridged[bridge_arr]]
                if new_bridge.size == 0:
                    continue
                if added_count + int(new_bridge.size) > max_added:
                    continue
                bridged[np.asarray(path, dtype=np.int64)] = True
                bridge_mask[new_bridge] = True
                added_count += int(new_bridge.size)
                first = source_labels[0]
                for label in source_labels[1:]:
                    union(first, label)
            if not bridge_mask.any():
                return source_mask.copy(), bridge_mask
            return bridged, bridge_mask

        objects: list[dict] = []
        sequence = 0
        for kind in kind_order:
            for source in source_by_kind[kind]:
                source_cells = clean_cells(source.get("cells", []))
                if source_cells.size == 0:
                    continue
                mask = np.zeros(grid.n, dtype=bool)
                mask[source_cells] = True
                raw_mask = mask.copy()
                mask, bridge_mask = bridge_short_gaps(kind, mask)
                components = self._connected_components(grid, mask)
                components.sort(
                    key=lambda comp: float(area[np.asarray(comp, dtype=np.int64)].sum()),
                    reverse=True,
                )
                for component_index, comp in enumerate(components):
                    comp = clean_cells(comp)
                    if comp.size == 0:
                        continue
                    if comp.size < 2 and kind not in {"transform", "suture"}:
                        continue
                    path_cells, endpoint_cells, endpoint_count, junction_count, high_degree_count = (
                        ordered_component_path(comp)
                    )
                    if not path_cells:
                        continue
                    metrics = path_metrics(comp, path_cells)
                    path_arr = np.asarray(path_cells, dtype=np.int64)
                    lat, lon = self._cell_centroid_lat_lon(grid, path_arr)
                    pids = list(source.get("parent_plate_ids", []))
                    if not pids:
                        pids = self._boundary_adjacent_plate_ids(grid, plate, comp)
                    objects.append({
                        "id": (
                            f"p129:{kind}:{int(comp.min())}:"
                            f"{component_index}:{sequence}:{round(float(t), 1)}"
                        ),
                        "type": "tectonic_boundary_polyline",
                        "kind": kind_to_polyline[kind],
                        "role": kind_to_role[kind],
                        "boundary_kind": kind,
                        "source_boundary_object_id": source.get("id"),
                        "source_stage": source.get("stage", ""),
                        "raw_source_component_cell_count": int(np.count_nonzero(raw_mask[comp])),
                        "source_component_cell_count": int(comp.size),
                        "axis_cell_count": int(len(path_cells)),
                        "cell_count": int(len(path_cells)),
                        "p131_gap_bridge_cell_count": int(np.count_nonzero(bridge_mask[comp])),
                        "p131_gap_bridge_used": bool(np.any(bridge_mask[comp])),
                        "area_fraction": float(area[comp].sum() / total_area),
                        "path_length_m": metrics["path_length_m"],
                        "direct_distance_m": metrics["direct_distance_m"],
                        "directness": metrics["directness"],
                        "sinuosity": metrics["sinuosity"],
                        "path_coverage_fraction": metrics["path_coverage_fraction"],
                        "endpoint_cells": endpoint_cells,
                        "source_endpoint_count": int(endpoint_count),
                        "source_junction_cell_count": int(junction_count),
                        "source_high_degree_cell_count": int(high_degree_count),
                        "parent_plate_ids": pids,
                        "birth_myr": source.get("birth_myr", round(float(t), 1)),
                        "last_active_myr": round(float(t), 1),
                        "lat": lat,
                        "lon": lon,
                        "cells": [int(cell) for cell in path_cells],
                    })
                    sequence += 1
        objects.sort(key=lambda obj: (
            str(obj.get("boundary_kind", "")),
            int(obj.get("source_component_cell_count", 0)) * -1,
            str(obj.get("id", "")),
        ))
        return objects

    def _record_p129_boundary_polyline_metrics(self, world, grid, objects):
        area = np.asarray(grid.cell_area, dtype=np.float64)
        total = max(float(area.sum()), 1.0e-12)
        source_cells = 0
        axis_cells = 0
        path_coverages: list[float] = []
        directness_values: list[float] = []
        sinuosity_values: list[float] = []
        junction_cells = 0
        high_degree_cells = 0
        p131_bridge_cells = 0
        p131_bridge_objects = 0
        by_boundary_kind: dict[str, int] = {}
        p131_bridge_cells_by_kind: dict[str, int] = {}
        for obj in objects or []:
            boundary_kind = str(obj.get("boundary_kind", ""))
            by_boundary_kind[boundary_kind] = by_boundary_kind.get(boundary_kind, 0) + 1
            source_cells += int(obj.get("source_component_cell_count", 0))
            axis_cells += int(obj.get("axis_cell_count", 0))
            path_coverages.append(float(obj.get("path_coverage_fraction", 0.0)))
            directness_values.append(float(obj.get("directness", 0.0)))
            sinuosity_values.append(float(obj.get("sinuosity", 1.0)))
            junction_cells += int(obj.get("source_junction_cell_count", 0))
            high_degree_cells += int(obj.get("source_high_degree_cell_count", 0))
            bridge_cells = int(obj.get("p131_gap_bridge_cell_count", 0))
            p131_bridge_cells += bridge_cells
            if bridge_cells > 0:
                p131_bridge_objects += 1
                p131_bridge_cells_by_kind[boundary_kind] = (
                    p131_bridge_cells_by_kind.get(boundary_kind, 0) + bridge_cells
                )
        world.set_g("tectonics.last_p129_object_count", float(len(objects or [])))
        world.set_g("tectonics.last_p129_source_cell_count", float(source_cells))
        world.set_g("tectonics.last_p129_axis_cell_count", float(axis_cells))
        world.set_g(
            "tectonics.last_p129_axis_source_coverage_fraction",
            float(axis_cells / max(float(source_cells), 1.0)),
        )
        world.set_g(
            "tectonics.last_p129_mean_path_coverage_fraction",
            float(np.mean(path_coverages) if path_coverages else 0.0),
        )
        world.set_g(
            "tectonics.last_p129_mean_directness",
            float(np.mean(directness_values) if directness_values else 0.0),
        )
        world.set_g(
            "tectonics.last_p129_mean_sinuosity",
            float(np.mean(sinuosity_values) if sinuosity_values else 1.0),
        )
        world.set_g(
            "tectonics.last_p129_max_sinuosity",
            float(np.max(sinuosity_values) if sinuosity_values else 1.0),
        )
        world.set_g(
            "tectonics.last_p129_source_junction_cell_count",
            float(junction_cells),
        )
        world.set_g(
            "tectonics.last_p129_source_high_degree_cell_count",
            float(high_degree_cells),
        )
        world.set_g(
            "tectonics.last_p129_polyline_ready",
            1.0 if objects else 0.0,
        )
        world.set_g(
            "tectonics.last_p131_gap_bridge_cell_count",
            float(p131_bridge_cells),
        )
        world.set_g(
            "tectonics.last_p131_gap_bridge_object_count",
            float(p131_bridge_objects),
        )
        world.set_g(
            "tectonics.last_p131_gap_bridge_used",
            1.0 if p131_bridge_cells > 0 else 0.0,
        )
        for kind in (
            "ridge",
            "transform",
            "trench",
            "suture",
            "active_margin",
            "passive_margin",
            "convergent_parent",
            "subduction_parent",
        ):
            world.set_g(
                f"tectonics.last_p129_{kind}_polyline_count",
                float(by_boundary_kind.get(kind, 0)),
            )
            world.set_g(
                f"tectonics.last_p131_{kind}_gap_bridge_cell_count",
                float(p131_bridge_cells_by_kind.get(kind, 0)),
            )

    def _boundary_province_objects(self, grid, boundary_objects, ctype, age, t):
        ctype = np.asarray(ctype, dtype=np.float64)
        age = np.asarray(age, dtype=np.float64)
        province_codes = {
            "mid_ocean_ridge": 1,
            "ridge_transform_fault": 2,
            "inactive_fracture_zone": 3,
            "ocean_ocean_subduction_trench": 4,
            "continental_arc_margin": 5,
            "island_arc_trench": 6,
            "backarc_spreading_center": 7,
            "continent_continent_collision": 8,
            "suture_zone": 9,
            "strike_slip_boundary": 10,
            "passive_margin": 11,
            "continental_rift": 12,
        }
        priority = {
            "continent_continent_collision": 90,
            "suture_zone": 80,
            "ocean_ocean_subduction_trench": 75,
            "island_arc_trench": 72,
            "continental_arc_margin": 70,
            "mid_ocean_ridge": 65,
            "backarc_spreading_center": 62,
            "ridge_transform_fault": 58,
            "strike_slip_boundary": 55,
            "passive_margin": 40,
            "continental_rift": 38,
            "inactive_fracture_zone": 20,
        }
        field = np.full(grid.n, -1.0, dtype=np.float64)
        field_priority = np.full(grid.n, -1, dtype=np.int16)
        median_edge = float(np.median(grid.edge_lengths)) if grid.edges.size else 0.0
        out: list[dict] = []
        for obj in boundary_objects or []:
            cells = np.asarray(obj.get("cells", []), dtype=np.int64)
            cells = cells[(0 <= cells) & (cells < grid.n)]
            if cells.size == 0:
                continue
            province_kind = self._classify_boundary_province(obj, ctype, age, cells)
            code = int(province_codes.get(province_kind, -1))
            if code < 0:
                continue
            rank = int(priority.get(province_kind, 0))
            update = rank >= field_priority[cells]
            if np.any(update):
                field[cells[update]] = float(code)
                field_priority[cells[update]] = rank
            bcf = float(obj.get("boundary_continental_fraction", 0.0))
            province_age = float(obj.get("age_myr", max(float(t) - float(obj.get("birth_myr", t)), 0.0)))
            out.append({
                "id": f"boundary_province:{province_kind}:{obj.get('id')}",
                "kind": province_kind,
                "source_boundary_kind": str(obj.get("kind", "")),
                "source_boundary_object_id": obj.get("id"),
                "adjacent_plate_ids": list(obj.get("parent_plate_ids", [])),
                "cell_count": int(cells.size),
                "area_fraction": float(grid.cell_area[cells].sum() / max(float(grid.cell_area.sum()), 1.0e-12)),
                "length_estimate_m": float(max(cells.size - 1, 1) * median_edge),
                "width_estimate_cells": self._boundary_province_width_cells(province_kind),
                "age_myr": round(province_age, 1),
                "birth_myr": obj.get("birth_myr"),
                "last_active_myr": obj.get("last_active_myr", round(float(t), 1)),
                "lifecycle_stage": obj.get("stage", ""),
                "boundary_continental_fraction": bcf,
                "mean_crust_age_myr": float(np.average(age[cells], weights=grid.cell_area[cells])),
                "kinematics": self._boundary_province_kinematics(province_kind),
                "lat": obj.get("lat"),
                "lon": obj.get("lon"),
                "cells": cells.astype(int).tolist() if cells.size <= 700 else [],
            })
        out.sort(key=lambda item: (-float(item["area_fraction"]), str(item["id"])))
        return out, field

    def _classify_boundary_province(self, obj, ctype, age, cells):
        source = str(obj.get("kind", ""))
        bcf = float(obj.get("boundary_continental_fraction", 0.0))
        ocean_fraction = float(np.mean(ctype[cells] < 0.5)) if cells.size else 0.0
        young_fraction = float(np.mean(age[cells] <= 45.0)) if cells.size else 0.0
        if source == "ridge":
            return "mid_ocean_ridge"
        if source == "transform":
            return "ridge_transform_fault" if young_fraction >= 0.20 else "strike_slip_boundary"
        if source == "trench":
            if bcf >= 0.35:
                return "continental_arc_margin"
            return "ocean_ocean_subduction_trench" if ocean_fraction >= 0.75 else "island_arc_trench"
        if source == "active_margin":
            return "continental_arc_margin" if bcf >= 0.25 else "island_arc_trench"
        if source == "suture":
            return "continent_continent_collision" if bcf >= 0.45 else "suture_zone"
        if source == "passive_margin":
            return "passive_margin"
        if source == "divergent":
            return "backarc_spreading_center" if ocean_fraction >= 0.65 else "continental_rift"
        return source or "unknown_boundary_province"

    @staticmethod
    def _boundary_province_width_cells(province_kind):
        return {
            "mid_ocean_ridge": 2,
            "ridge_transform_fault": 1,
            "inactive_fracture_zone": 1,
            "ocean_ocean_subduction_trench": 1,
            "continental_arc_margin": 4,
            "island_arc_trench": 3,
            "backarc_spreading_center": 2,
            "continent_continent_collision": 6,
            "suture_zone": 3,
            "strike_slip_boundary": 1,
            "passive_margin": 4,
            "continental_rift": 4,
        }.get(str(province_kind), 1)

    @staticmethod
    def _boundary_province_kinematics(province_kind):
        if province_kind in {"mid_ocean_ridge", "backarc_spreading_center", "continental_rift"}:
            return "extensional"
        if province_kind in {
            "ocean_ocean_subduction_trench",
            "continental_arc_margin",
            "island_arc_trench",
            "continent_continent_collision",
            "suture_zone",
        }:
            return "convergent"
        if province_kind in {"ridge_transform_fault", "strike_slip_boundary", "inactive_fracture_zone"}:
            return "strike_slip"
        if province_kind == "passive_margin":
            return "passive"
        return "mixed"

    def _ridge_offset_transform_axis(self, grid, ridge_axis, plate_contact,
                                     ctype, excluded):
        """Create a narrow transform proxy near ridge offsets when shear is unresolved."""
        ridge_axis = np.asarray(ridge_axis, dtype=bool)
        plate_contact = np.asarray(plate_contact, dtype=bool)
        ctype = np.asarray(ctype, dtype=np.float64)
        excluded = np.asarray(excluded, dtype=bool)
        out = np.zeros(grid.n, dtype=bool)
        if not ridge_axis.any():
            return out
        allowed = plate_contact & (ctype == OCEAN) & ~ridge_axis & ~excluded
        if not allowed.any():
            return out
        near_ridge = self._dilate_mask(grid, ridge_axis, allowed=allowed, passes=2)
        candidates = np.where(near_ridge & allowed)[0]
        if candidates.size == 0:
            candidates = np.where(allowed)[0]
        if candidates.size == 0:
            return out
        target = min(
            int(candidates.size),
            max(4, int(np.ceil(0.65 * float(np.count_nonzero(ridge_axis))))),
        )
        scores = []
        for c in candidates:
            c = int(c)
            ridge_neighbors = int(np.count_nonzero(ridge_axis[grid.neighbors[c]]))
            scores.append((ridge_neighbors, -c, c))
        for _, _, c in sorted(scores, reverse=True)[:target]:
            out[int(c)] = True
        return out

    def _boundary_lifecycle_components(self, grid, comps, plate, *, kind: str):
        """Keep sparse ridge/transform evidence as lifecycle objects.

        High-resolution boundary masks can thin ridge and transform axes into
        isolated cells on the spherical graph.  The lifecycle layer still needs
        those cells grouped by plate pair, while larger connected components
        should pass through unchanged.
        """
        large = [np.asarray(comp, dtype=np.int64) for comp in comps if comp.size >= 3]
        small = [np.asarray(comp, dtype=np.int64) for comp in comps if comp.size < 3]
        if not small:
            return large
        grouped: dict[tuple[int, ...], list[int]] = {}
        for comp in small:
            pids = tuple(self._boundary_adjacent_plate_ids(grid, plate, comp))
            if len(pids) < 2:
                pids = tuple(sorted(int(x) for x in np.unique(plate[comp])))
            if len(pids) < 2:
                pids = (int(plate[int(comp[0])]),)
            grouped.setdefault(pids, []).extend(int(c) for c in comp)
        aggregates = []
        for pids, cells in sorted(grouped.items(), key=lambda item: (item[0], item[1])):
            unique = np.asarray(sorted(set(cells)), dtype=np.int64)
            if unique.size < 3:
                continue
            aggregates.append(unique)
        if not aggregates:
            all_small = np.asarray(
                sorted({int(c) for comp in small for c in comp}),
                dtype=np.int64,
            )
            if all_small.size >= 3:
                aggregates.append(all_small)
        return large + aggregates

    def _boundary_adjacent_plate_ids(self, grid, plate, comp):
        ids: set[int] = set(int(x) for x in np.unique(plate[comp]))
        comp_mask = np.zeros(grid.n, dtype=bool)
        comp_mask[comp] = True
        for c in comp:
            for nb in grid.neighbors[int(c)]:
                ids.add(int(plate[int(nb)]))
        return sorted(ids)

    def _match_previous_boundary_object(self, comp, previous, used_previous):
        current = set(int(c) for c in comp)
        best = None
        best_score = 0.0
        for obj in previous:
            obj_id = str(obj.get("id"))
            if obj_id in used_previous:
                continue
            prev_cells = set(int(c) for c in obj.get("cells", []))
            if not prev_cells:
                continue
            overlap = len(current & prev_cells)
            if overlap == 0:
                continue
            score = overlap / max(min(len(current), len(prev_cells)), 1)
            if score > best_score:
                best_score = score
                best = obj
        return best if best_score >= 0.25 else None

    def _boundary_context(self, grid, comp, plate, pids, ctype=None, age=None,
                          thick=None, velocity=None, kind=""):
        zone = np.zeros(grid.n, dtype=bool)
        zone[comp] = True
        for c in comp:
            zone[grid.neighbors[int(c)]] = True
        out: dict[str, object] = {
            "subducting_plate_id": None,
            "overriding_plate_id": None,
            "subduction_polarity": None,
            "polarity_basis": None,
            "boundary_continental_fraction": 0.0,
            "mean_oceanic_age_by_plate": {},
            "continental_fraction_by_plate": {},
            "relative_motion": {},
        }
        if ctype is not None:
            ctype_comp = np.asarray(ctype)[comp]
            out["boundary_continental_fraction"] = float(np.mean(ctype_comp >= 0.5))
        plate_stats: dict[int, dict[str, float]] = {}
        for pid in pids:
            local = zone & (plate == pid)
            if not local.any():
                continue
            area = grid.cell_area[local]
            area_sum = float(area.sum())
            ctype_local = np.asarray(ctype)[local] if ctype is not None else np.zeros(local.sum())
            age_local = np.asarray(age)[local] if age is not None else np.zeros(local.sum())
            thick_local = np.asarray(thick)[local] if thick is not None else np.zeros(local.sum())
            ocean = ctype_local < 0.5
            cont_frac = float(np.average((ctype_local >= 0.5).astype(float), weights=area))
            if ocean.any():
                ocean_weights = area[ocean]
                ocean_age = float(np.average(age_local[ocean], weights=ocean_weights))
                ocean_thick = float(np.average(thick_local[ocean], weights=ocean_weights))
            else:
                ocean_age = 0.0
                ocean_thick = float(np.average(thick_local, weights=area)) if area_sum else 0.0
            plate_stats[int(pid)] = {
                "area": area_sum,
                "continental_fraction": cont_frac,
                "oceanic_fraction": 1.0 - cont_frac,
                "mean_oceanic_age": ocean_age,
                "mean_oceanic_thickness": ocean_thick,
            }
        out["mean_oceanic_age_by_plate"] = {
            str(pid): round(stats["mean_oceanic_age"], 2)
            for pid, stats in plate_stats.items()
        }
        out["continental_fraction_by_plate"] = {
            str(pid): round(stats["continental_fraction"], 3)
            for pid, stats in plate_stats.items()
        }
        if velocity is not None and len(pids) >= 2:
            vel = np.asarray(velocity)
            mean_vel = {}
            for pid in pids:
                local = zone & (plate == pid)
                if local.any():
                    mean_vel[int(pid)] = np.average(vel[local], axis=0,
                                                    weights=grid.cell_area[local])
            if len(mean_vel) >= 2:
                ids = sorted(mean_vel)[:2]
                rel = mean_vel[ids[0]] - mean_vel[ids[1]]
                out["relative_motion"] = {
                    "plate_pair": ids,
                    "speed_m_per_s": float(np.linalg.norm(rel)),
                }
        if (
            kind in {
                "trench",
                "subduction",
                "active_margin",
                "convergent_parent",
                "subduction_parent",
            }
            and len(plate_stats) >= 2
        ):
            subducting, overriding, basis = self._resolve_subduction_polarity(plate_stats)
            out["subducting_plate_id"] = subducting
            out["overriding_plate_id"] = overriding
            out["subduction_polarity"] = (
                f"plate_{subducting}_under_plate_{overriding}"
                if subducting is not None and overriding is not None else None
            )
            out["polarity_basis"] = basis
        if kind == "transform":
            out["shear_dominance"] = 1.0
        return out

    def _resolve_subduction_polarity(self, plate_stats):
        ids = sorted(plate_stats)
        if len(ids) < 2:
            return None, None, "insufficient_adjacent_plates"
        # Consider the two largest adjacent plates for the local boundary.
        ids = sorted(ids, key=lambda pid: plate_stats[pid]["area"], reverse=True)[:2]
        a, b = ids[0], ids[1]
        sa, sb = plate_stats[a], plate_stats[b]
        a_cont, b_cont = sa["continental_fraction"], sb["continental_fraction"]
        if a_cont < 0.35 <= b_cont:
            return a, b, "oceanic_plate_subducts_beneath_continent"
        if b_cont < 0.35 <= a_cont:
            return b, a, "oceanic_plate_subducts_beneath_continent"
        if a_cont < 0.35 and b_cont < 0.35:
            if sa["mean_oceanic_age"] > sb["mean_oceanic_age"]:
                return a, b, "older_oceanic_lithosphere_subducts"
            if sb["mean_oceanic_age"] > sa["mean_oceanic_age"]:
                return b, a, "older_oceanic_lithosphere_subducts"
            if sa["mean_oceanic_thickness"] >= sb["mean_oceanic_thickness"]:
                return a, b, "thicker_oceanic_lithosphere_subducts_tiebreak"
            return b, a, "thicker_oceanic_lithosphere_subducts_tiebreak"
        return None, None, "continental_collision_no_subduction_polarity"

    def _wilson_lineage_key(self, obj):
        pids = sorted(int(x) for x in obj.get("parent_plate_ids", []) if int(x) >= 0)
        if len(pids) >= 2:
            return f"plates:{pids[0]}-{pids[1]}"
        if pids:
            return f"plate:{pids[0]}"
        return f"boundary:{obj.get('id')}"

    def _previous_lifecycle_by_lineage(self, previous, key):
        out: dict[str, dict] = {}
        for obj in previous.get(key, []) if isinstance(previous, dict) else []:
            lineage = obj.get("lineage_key")
            if lineage is not None and lineage not in out:
                out[str(lineage)] = obj
        return out

    def _lifecycle_birth(self, previous_by_lineage, lineage, fallback):
        prev = previous_by_lineage.get(lineage)
        if prev is not None:
            return float(prev.get("birth_myr", fallback))
        return float(fallback)

    def _wilson_lifecycle_objects(self, boundary_objects, previous, t):
        phase_by_kind = {
            "divergent": ("continental_rift", WILSON_RIFT, "rifting"),
            "ridge": ("spreading_ocean", WILSON_OPENING, "opening"),
            "passive_margin": ("mature_ocean", WILSON_MATURE, "open"),
            "trench": ("closing_ocean", WILSON_CLOSING, "closing"),
            "active_margin": ("closing_arc_margin", WILSON_COLLISION, "closing"),
            "suture": ("suture_relict", WILSON_SUTURE, "closed"),
        }
        priority = {
            WILSON_RIFT: 1,
            WILSON_OPENING: 2,
            WILSON_MATURE: 3,
            WILSON_CLOSING: 4,
            WILSON_COLLISION: 5,
            WILSON_SUTURE: 6,
        }
        previous = previous or {}
        prev_basins = self._previous_lifecycle_by_lineage(previous, "ocean_basins")
        prev_rifts = self._previous_lifecycle_by_lineage(previous, "rift_systems")
        prev_passive = self._previous_lifecycle_by_lineage(previous, "passive_margins")
        prev_spreading = self._previous_lifecycle_by_lineage(previous, "spreading_centers")
        prev_closing = self._previous_lifecycle_by_lineage(previous, "closing_margins")
        prev_sutures = self._previous_lifecycle_by_lineage(previous, "sutures")

        basin_by_lineage: dict[str, dict] = {}
        rift_systems: list[dict] = []
        passive_margins: list[dict] = []
        spreading_centers: list[dict] = []
        closing_margins: list[dict] = []
        sutures: list[dict] = []

        def lifecycle_obj(prefix, lineage, boundary, stage, phase_code, prev_map,
                          extra=None):
            birth = self._lifecycle_birth(
                prev_map, lineage, boundary.get("birth_myr", t))
            obj = {
                "id": (
                    prev_map.get(lineage, {}).get("id")
                    if lineage in prev_map else f"{prefix}:{lineage}"
                ),
                "kind": prefix,
                "lineage_key": lineage,
                "stage": stage,
                "phase_code": float(phase_code),
                "boundary_object_id": boundary.get("id"),
                "parent_boundary_ids": [boundary.get("id")],
                "parent_plate_ids": boundary.get("parent_plate_ids", []),
                "cells": boundary.get("cells", []),
                "area_fraction": float(boundary.get("area_fraction", 0.0)),
                "lat": boundary.get("lat"),
                "lon": boundary.get("lon"),
                "birth_myr": round(float(birth), 1),
                "last_active_myr": round(float(t), 1),
                "age_myr": round(float(max(t - birth, 0.0)), 1),
            }
            if extra:
                obj.update(extra)
            return obj

        for boundary in boundary_objects:
            kind = boundary.get("kind")
            if kind not in phase_by_kind:
                continue
            stage, phase_code, status = phase_by_kind[kind]
            lineage = self._wilson_lineage_key(boundary)
            basin = basin_by_lineage.get(lineage)
            if basin is None:
                birth = self._lifecycle_birth(
                    prev_basins, lineage, boundary.get("birth_myr", t))
                basin = {
                    "id": (
                        prev_basins.get(lineage, {}).get("id")
                        if lineage in prev_basins else f"ocean_basin:{lineage}"
                    ),
                    "kind": "ocean_basin_lifecycle",
                    "lineage_key": lineage,
                    "stage": stage,
                    "status": status,
                    "phase_code": float(phase_code),
                    "boundary_object_ids": [],
                    "parent_plate_ids": boundary.get("parent_plate_ids", []),
                    "cells": [],
                    "area_fraction": 0.0,
                    "birth_myr": round(float(birth), 1),
                    "last_active_myr": round(float(t), 1),
                    "age_myr": round(float(max(t - birth, 0.0)), 1),
                    "lat": boundary.get("lat"),
                    "lon": boundary.get("lon"),
                    "parent_rift_system_ids": [],
                    "parent_margin_ids": [],
                    "parent_spreading_center_ids": [],
                    "parent_closing_margin_ids": [],
                    "parent_suture_ids": [],
                }
                basin_by_lineage[lineage] = basin
            if priority[phase_code] >= priority.get(float(basin.get("phase_code", 0.0)), 0):
                basin["stage"] = stage
                basin["status"] = status
                basin["phase_code"] = float(phase_code)
                basin["lat"] = boundary.get("lat")
                basin["lon"] = boundary.get("lon")
            basin["boundary_object_ids"].append(boundary.get("id"))
            basin["area_fraction"] += float(boundary.get("area_fraction", 0.0))
            basin["cells"].extend(boundary.get("cells", []))

            if kind == "divergent":
                rift = lifecycle_obj("rift_system", lineage, boundary,
                                     "continental_rift", WILSON_RIFT,
                                     prev_rifts)
                rift["parent_basin_id"] = basin["id"]
                rift_systems.append(rift)
                basin["parent_rift_system_ids"].append(rift["id"])
            elif kind == "ridge":
                rift_stage = (
                    "narrow_sea"
                    if float(boundary.get("boundary_continental_fraction", 0.0)) > 0.15
                    else "spreading_ocean"
                )
                rift = lifecycle_obj("rift_system", lineage, boundary, rift_stage,
                                     WILSON_OPENING, prev_rifts)
                spread = lifecycle_obj("spreading_center", lineage, boundary,
                                       "active_spreading_center", WILSON_OPENING,
                                       prev_spreading)
                rift["parent_basin_id"] = basin["id"]
                spread["parent_basin_id"] = basin["id"]
                rift_systems.append(rift)
                spreading_centers.append(spread)
                basin["parent_rift_system_ids"].append(rift["id"])
                basin["parent_spreading_center_ids"].append(spread["id"])
            elif kind == "passive_margin":
                margin = lifecycle_obj("passive_margin", lineage, boundary,
                                       "mature_passive_margin", WILSON_MATURE,
                                       prev_passive,
                                       extra={"parent_basin_id": basin["id"]})
                passive_margins.append(margin)
                basin["parent_margin_ids"].append(margin["id"])
            elif kind in {"trench", "active_margin"}:
                closing = lifecycle_obj("closing_margin", lineage, boundary,
                                        "subduction_closing_margin", phase_code,
                                        prev_closing,
                                        extra={
                                            "parent_basin_id": basin["id"],
                                            "subducting_plate_id": boundary.get("subducting_plate_id"),
                                            "overriding_plate_id": boundary.get("overriding_plate_id"),
                                            "polarity_basis": boundary.get("polarity_basis"),
                                        })
                closing_margins.append(closing)
                basin["parent_closing_margin_ids"].append(closing["id"])
            elif kind == "suture":
                suture = lifecycle_obj("suture", lineage, boundary,
                                       "closed_suture", WILSON_SUTURE,
                                       prev_sutures,
                                       extra={"parent_basin_id": basin["id"]})
                sutures.append(suture)
                basin["parent_suture_ids"].append(suture["id"])

        ocean_basins = list(basin_by_lineage.values())
        for basin in ocean_basins:
            basin["cells"] = sorted(set(int(c) for c in basin.get("cells", [])))
            basin["boundary_object_ids"] = [x for x in dict.fromkeys(basin["boundary_object_ids"]) if x is not None]
            for key in (
                "parent_rift_system_ids", "parent_margin_ids",
                "parent_spreading_center_ids", "parent_closing_margin_ids",
                "parent_suture_ids",
            ):
                basin[key] = [x for x in dict.fromkeys(basin[key]) if x is not None]
        ocean_basins.sort(key=lambda obj: (priority.get(float(obj.get("phase_code", 0.0)), 0),
                                           obj.get("area_fraction", 0.0)),
                          reverse=True)
        wilson_cycles = self._wilson_cycle_objects_from_basins(ocean_basins)
        ocean_gateways = self._ocean_gateway_objects_from_lifecycle(
            boundary_objects, wilson_cycles, ocean_basins, t)
        return {
            "tectonics.ocean_basins": ocean_basins[:32],
            "tectonics.rift_systems": rift_systems[:32],
            "tectonics.passive_margins": passive_margins[:32],
            "tectonics.spreading_centers": spreading_centers[:32],
            "tectonics.closing_margins": closing_margins[:32],
            "tectonics.sutures": sutures[:32],
            "tectonics.wilson_cycles": wilson_cycles[:32],
            "tectonics.ocean_gateways": ocean_gateways[:32],
        }

    def _breakup_seaway_objects(self, world, grid, ctype, domain, stability,
                                rift_potential, continent_id, plate, boundaries,
                                plate_topologies, rift_systems, t):
        """Resolve supercontinent-breakup corridors from weak rifted lithosphere.

        P21 keeps the terrain pass from inventing arbitrary cuts.  Tectonics
        identifies where a large continent has a plausible multi-rift corridor;
        terrain may later drown that corridor only if it actually improves
        exposed-land topology.
        """
        world.objects["tectonics.breakup_component_telemetry"] = []
        target = float(np.clip(world.spec.target_land_fraction, 0.0, 1.0))
        ctype = np.asarray(ctype, dtype=np.float64)
        (
            historical_breakup_pressure,
            historical_breakup_active,
            historical_supercontinent_share,
        ) = self._p110b_historical_supercontinent_breakup_pressure(
            world, grid, ctype, target, t)
        historical_residence_debt_myr = float(world.g(
            "tectonics.last_p110b_historical_supercontinent_residence_debt_myr",
            0.0,
        ))
        historical_residence_controller_active = bool(world.g(
            "tectonics.last_p110b_historical_residence_controller_active",
            0.0,
        ) >= 0.5)
        previous_breakup_seaways = list(
            world.objects.get("tectonics.breakup_seaways", [])
        )
        retained_historical = (
            self._retained_historical_breakup_seaway_objects(
                previous_breakup_seaways,
                t,
                historical_breakup_active=bool(historical_breakup_active),
            )
        )
        if (
            not str(world.spec.name).startswith("earthlike")
            or target < 0.20
            or target > 0.38
        ):
            return []
        if (
            float(t) < 2400.0
            and not bool(historical_breakup_active)
            and not retained_historical
        ):
            return []

        domain = np.asarray(domain, dtype=np.float64)
        stability = np.asarray(stability, dtype=np.float64)
        rift = np.nan_to_num(
            np.asarray(rift_potential, dtype=np.float64),
            nan=0.0, posinf=1.0, neginf=0.0)
        continent_id = np.asarray(continent_id, dtype=np.float64).astype(int)
        plate = np.asarray(plate, dtype=int)
        cont = ctype == CONT
        if not cont.any():
            return retained_historical

        total = float(grid.cell_area.sum())
        cont_area = float(grid.cell_area[cont].sum())
        if cont_area <= 0.0:
            return retained_historical

        previous = {
            str(obj.get("lineage_key")): obj
            for obj in previous_breakup_seaways
            if obj.get("lineage_key") is not None
        }
        rift_by_cell: dict[int, list[str]] = {}
        for obj in rift_systems or []:
            rid = str(obj.get("id", ""))
            for c in obj.get("cells", []):
                ci = int(c)
                if 0 <= ci < grid.n:
                    rift_by_cell.setdefault(ci, []).append(rid)

        boundary_source = self._cells_for_boundary(
            grid.n, boundaries, "ridge", "divergent", "passive_margin",
            "transform")
        suture_source = self._cells_for_boundary(
            grid.n, boundaries, "suture", "collision")
        plate_by_id = {
            int(obj.get("numeric_id", -1)): obj
            for obj in plate_topologies or []
        }

        objects: list[dict] = []
        component_telemetry: list[dict] = []
        components = self._connected_components(grid, cont)
        components.sort(key=lambda comp: float(grid.cell_area[comp].sum()), reverse=True)
        for comp_index, comp in enumerate(components[:8]):
            comp_mask = np.zeros(grid.n, dtype=bool)
            comp_mask[comp] = True
            comp_area = float(grid.cell_area[comp].sum())
            comp_share_of_cont = comp_area / max(cont_area, 1.0)
            comp_share_total = comp_area / max(total, 1.0)
            boundary_hit_count = int(np.count_nonzero(boundary_source[comp]))
            boundary_hit_area_fraction = float(
                grid.cell_area[comp[boundary_source[comp]]].sum() / max(total, 1.0)
            ) if boundary_hit_count else 0.0
            medium_rifted_component = (
                comp_share_of_cont >= 0.16
                and comp_share_total >= 0.055
                and (
                    float(np.percentile(rift[comp], 88)) >= 0.48
                    or float(np.mean(stability[comp] < 0.52)) >= 0.10
                    or (
                        boundary_hit_count >= 3
                        and boundary_hit_area_fraction >= 0.00010
                    )
                )
            )
            telemetry = {
                "component_index": int(comp_index),
                "cell_count": int(comp.size),
                "component_share_of_continental_crust": float(comp_share_of_cont),
                "component_area_fraction": float(comp_share_total),
                "historical_breakup_pressure": float(historical_breakup_pressure),
                "historical_breakup_active": bool(historical_breakup_active),
                "historical_largest_continent_share": float(
                    historical_supercontinent_share),
                "historical_residence_debt_myr": float(
                    historical_residence_debt_myr),
                "historical_residence_controller_active": bool(
                    historical_residence_controller_active),
                "boundary_hit_count": int(boundary_hit_count),
                "boundary_hit_area_fraction": float(boundary_hit_area_fraction),
                "rift_p88": float(np.percentile(rift[comp], 88)),
                "weak_stability_fraction": float(np.mean(stability[comp] < 0.52)),
                "medium_rifted_component": bool(medium_rifted_component),
                "eligible": True,
                "skip_reason": "",
                "weak_area_fraction": 0.0,
                "candidate_count": 0,
                "accepted_object_count": 0,
                "best_topology_score": 0.0,
                "best_quality_score": 0.0,
            }
            if (
                comp_share_of_cont < 0.46
                and comp_share_total < 0.115
                and not medium_rifted_component
                and not (
                    historical_breakup_active
                    and comp_share_of_cont >= 0.54
                    and comp_share_total >= 0.12
                )
            ):
                telemetry["eligible"] = False
                telemetry["skip_reason"] = "below_large_or_medium_rifted_component_threshold"
                component_telemetry.append(telemetry)
                continue

            boundary = comp_mask & self._component_boundary_mask(grid, comp_mask)
            if int(boundary.sum()) < 4:
                telemetry["eligible"] = False
                telemetry["skip_reason"] = "too_few_component_boundary_cells"
                component_telemetry.append(telemetry)
                continue
            core = comp_mask & (
                (domain == DOMAIN_CRATON)
                | ((stability > 0.80) & (rift < 0.35))
            )
            interior = comp_mask & (
                (domain == DOMAIN_CONTINENTAL_INTERIOR)
                | ((stability > 0.62) & (rift < 0.45))
            )
            weak = comp_mask & (
                (rift >= max(0.42, float(np.percentile(rift[comp], 64))))
                | np.isin(domain, [DOMAIN_CONTINENTAL_MARGIN,
                                   DOMAIN_ACCRETED_TERRANE,
                                   DOMAIN_SUTURE])
                | (stability < 0.52)
            )
            if historical_residence_controller_active:
                weak |= comp_mask & ~core & (
                    (stability < 0.78)
                    | (rift >= max(0.16, float(np.percentile(rift[comp], 48))))
                )
            forced_rift_belt = np.zeros(grid.n, dtype=bool)
            if historical_breakup_active and comp_share_of_cont >= 0.54:
                xyz = np.asarray(grid.xyz[comp], dtype=np.float64)
                weights = np.asarray(grid.cell_area[comp], dtype=np.float64)
                centroid = np.average(xyz, axis=0, weights=weights)
                centered = xyz - centroid
                cov = (centered * weights[:, None]).T @ centered / max(
                    float(weights.sum()), 1.0e-12)
                eigvals, eigvecs = np.linalg.eigh(cov)
                axis = eigvecs[:, int(np.argmax(eigvals))]
                projection = xyz @ axis
                if float(np.nanmax(projection) - np.nanmin(projection)) > 1.0e-9:
                    median = self._weighted_percentile_local(
                        projection, weights, 50.0)
                    half_width = max(
                        0.5 * (
                            self._weighted_percentile_local(projection, weights, 58.0)
                            - self._weighted_percentile_local(projection, weights, 42.0)
                        ),
                        1.0e-6,
                    )
                    center_values = [median]
                    if historical_residence_controller_active:
                        center_values.extend([
                            self._weighted_percentile_local(projection, weights, 38.0),
                            self._weighted_percentile_local(projection, weights, 62.0),
                        ])
                        if historical_residence_debt_myr >= 450.0:
                            center_values.extend([
                                self._weighted_percentile_local(projection, weights, 28.0),
                                self._weighted_percentile_local(projection, weights, 72.0),
                            ])
                    width_scale = 1.65 if historical_residence_controller_active else 1.0
                    for center_value in center_values:
                        local = (
                            np.abs(projection - center_value)
                            <= half_width * width_scale
                        )
                        forced_rift_belt[comp[local]] = True
                    forced_rift_belt &= comp_mask & ~core & (stability < 0.84)
                    if forced_rift_belt.any():
                        weak |= forced_rift_belt
            boundary_zone = self._dilate_mask(
                grid, boundary_source | suture_source, allowed=comp_mask, passes=2)
            if forced_rift_belt.any():
                boundary_zone |= self._dilate_mask(
                    grid, forced_rift_belt, allowed=comp_mask, passes=1)
            weak |= boundary_zone
            weak_area_fraction = float(grid.cell_area[weak].sum() / max(comp_area, 1.0))
            telemetry["weak_area_fraction"] = weak_area_fraction
            telemetry["forced_rift_belt_area_fraction"] = float(
                grid.cell_area[forced_rift_belt].sum() / max(comp_area, 1.0))
            min_weak_area_fraction = 0.018 if historical_breakup_active else 0.035
            if weak_area_fraction < min_weak_area_fraction:
                telemetry["eligible"] = False
                telemetry["skip_reason"] = "weak_area_fraction_below_threshold"
                component_telemetry.append(telemetry)
                continue

            cost = (
                1.45
                + 1.15 * interior.astype(float)
                + 5.50 * core.astype(float)
                + 1.10 * np.clip(stability, 0.0, 1.0)
                - 1.55 * np.clip(rift, 0.0, 1.0)
                - 1.15 * (domain == DOMAIN_CONTINENTAL_MARGIN).astype(float)
                - 1.10 * (domain == DOMAIN_SUTURE).astype(float)
                - 0.90 * (domain == DOMAIN_ACCRETED_TERRANE).astype(float)
                - 0.75 * boundary_zone.astype(float)
            )
            cost = np.where(comp_mask, np.maximum(cost, 0.20), np.inf)

            parent_cids = sorted(
                int(x) for x in np.unique(continent_id[comp]) if int(x) >= 0)
            parent_pid_counts = {
                int(pid): int(np.count_nonzero(plate[comp] == pid))
                for pid in np.unique(plate[comp])
                if int(pid) >= 0
            }
            parent_pids = [
                int(pid)
                for pid, _ in sorted(
                    parent_pid_counts.items(), key=lambda item: (-item[1], item[0]))
            ][:4]

            candidates: list[dict[str, object]] = []
            axes = self._breakup_candidate_axes(grid, comp, parent_pids, plate_by_id)
            max_path_cells = max(8, min(int(comp.size * 0.36), max(24, grid.n // 16)))

            def partition_metrics(cells: np.ndarray) -> tuple[float, float, float, int]:
                cut = np.zeros(grid.n, dtype=bool)
                cut[np.asarray(cells, dtype=np.int64)] = True
                cut &= comp_mask
                parts = self._connected_components(grid, comp_mask & ~cut)
                if len(parts) < 2:
                    return 1.0, 0.0, 0.0, len(parts)
                part_areas = np.array(
                    [float(grid.cell_area[p].sum()) for p in parts],
                    dtype=np.float64,
                )
                part_areas.sort()
                part_areas = part_areas[::-1]
                remaining_area = float(part_areas.sum())
                largest_after = float(part_areas[0] / max(remaining_area, 1.0))
                split_balance = float(part_areas[1] / max(part_areas[0], 1.0))
                partition_gain = max(0.0, 1.0 - largest_after)
                topology_score = partition_gain + 0.70 * min(split_balance, 1.0)
                if split_balance < 0.08:
                    topology_score *= 0.35
                return largest_after, split_balance, topology_score, len(parts)

            def add_candidate(axis_index, lineage_suffix: str, cells, *,
                              length_penalty: float) -> None:
                arr = np.unique(np.asarray(cells, dtype=np.int64))
                arr = arr[comp_mask[arr]]
                if arr.size < 4:
                    return
                path_area = float(grid.cell_area[arr].sum())
                weak_fraction = float(grid.cell_area[arr[weak[arr]]].sum()
                                      / max(path_area, 1.0))
                core_fraction = float(grid.cell_area[arr[core[arr]]].sum()
                                      / max(path_area, 1.0))
                mean_rift = float(np.average(rift[arr], weights=grid.cell_area[arr]))
                boundary_fraction = float(grid.cell_area[arr[boundary_zone[arr]]].sum()
                                          / max(path_area, 1.0))
                largest_after, split_balance, topology_score, part_count = (
                    partition_metrics(arr)
                )
                quality = (
                    1.9 * weak_fraction
                    + 1.1 * mean_rift
                    + 0.7 * boundary_fraction
                    + 1.65 * topology_score
                    + 0.35 * min(split_balance, 1.0)
                    - 2.8 * core_fraction
                    - length_penalty * float(arr.size)
                )
                relaxed_history_candidate = bool(
                    historical_breakup_active
                    and (
                        (
                            topology_score >= 0.18
                            and weak_fraction >= 0.12
                            and core_fraction <= 0.34
                        )
                        or (
                            historical_residence_controller_active
                            and topology_score >= 0.13
                            and weak_fraction >= 0.08
                            and core_fraction <= 0.42
                        )
                    )
                )
                if not relaxed_history_candidate and (
                    core_fraction > 0.42
                    or weak_fraction < 0.18
                    or (mean_rift < 0.30 and boundary_fraction < 0.22)
                ):
                    return
                candidates.append({
                    "axis_index": axis_index,
                    "lineage_suffix": lineage_suffix,
                    "cells": arr,
                    "path_area": path_area,
                    "weak_fraction": weak_fraction,
                    "core_fraction": core_fraction,
                    "mean_rift": mean_rift,
                    "boundary_fraction": boundary_fraction,
                    "largest_after_partition": largest_after,
                    "split_balance": split_balance,
                    "topology_score": topology_score,
                    "partition_count": part_count,
                    "quality": quality,
                })

            for axis_index, axis in enumerate(axes):
                values = grid.xyz @ axis
                local = values[comp]
                if local.size < 8:
                    continue
                lo = float(np.percentile(local, 14))
                hi = float(np.percentile(local, 86))
                low = boundary & (values <= lo)
                high = boundary & (values >= hi)
                if int(low.sum()) < 2 or int(high.sum()) < 2:
                    continue
                for starts, goals in ((low, high), (high, low)):
                    path = self._breakup_path_between_sets(
                        grid, starts, goals, comp_mask, cost,
                        max_path_cells=max_path_cells)
                    if len(path) < 4:
                        continue
                    add_candidate(
                        axis_index,
                        f"axis:{axis_index}",
                        path,
                        length_penalty=0.018,
                    )

            rift_floor = max(0.50, float(np.percentile(rift[comp], 72)))
            rift_axis = (
                comp_mask
                & ~core
                & (rift >= rift_floor)
                & (stability < 0.72)
            )
            rift_axis_comps = self._connected_components(grid, rift_axis)
            rift_axis_comps.sort(
                key=lambda cells: float(grid.cell_area[cells].sum()),
                reverse=True,
            )
            apron_allowed = (
                comp_mask
                & ~core
                & (
                    (rift >= 0.42)
                    | ((stability < 0.62) & ((rift >= 0.24) | boundary_zone))
                    | np.isin(domain, [DOMAIN_CONTINENTAL_MARGIN,
                                       DOMAIN_ACCRETED_TERRANE,
                                       DOMAIN_SUTURE])
                )
            )
            for rift_index, rift_comp in enumerate(rift_axis_comps[:8]):
                rift_mask = np.zeros(grid.n, dtype=bool)
                rift_mask[rift_comp] = True
                rift_mask = self._dilate_mask(
                    grid, rift_mask, allowed=apron_allowed, passes=1)
                cells = np.where(rift_mask)[0].astype(np.int64)
                path_area = float(grid.cell_area[cells].sum())
                area_share = path_area / max(comp_area, 1.0)
                if (
                    area_share < 0.006
                    or area_share > 0.180
                    or path_area / max(total, 1.0) > 0.075
                ):
                    continue
                largest_after, split_balance, topology_score, _ = partition_metrics(cells)
                if (
                    topology_score < 0.22
                    or split_balance < 0.08
                    or largest_after > 0.86
                ):
                    continue
                add_candidate(
                    f"rift:{rift_index}",
                    f"rift:{rift_index}",
                    cells,
                    length_penalty=0.006,
                )

            if (
                historical_breakup_active
                and comp_share_of_cont >= 0.62
                and comp_share_total >= 0.14
            ):
                historical_apron = (
                    comp_mask
                    & ~core
                    & (
                        forced_rift_belt
                        | boundary_zone
                        | (stability < 0.72)
                        | (
                            historical_residence_controller_active
                            & (stability < 0.80)
                        )
                        | (rift >= max(0.18, float(np.percentile(rift[comp], 55))))
                        | np.isin(domain, [
                            DOMAIN_CONTINENTAL_MARGIN,
                            DOMAIN_ACCRETED_TERRANE,
                            DOMAIN_SUTURE,
                        ])
                    )
                )
                if historical_apron.any():
                    fallback_axes: list[np.ndarray] = list(axes)
                    xyz = np.asarray(grid.xyz[comp], dtype=np.float64)
                    weights = np.asarray(grid.cell_area[comp], dtype=np.float64)
                    centroid = np.average(xyz, axis=0, weights=weights)
                    centered = xyz - centroid
                    cov = (centered * weights[:, None]).T @ centered / max(
                        float(weights.sum()), 1.0e-12)
                    eigvals, eigvecs = np.linalg.eigh(cov)
                    for eig_index in np.argsort(eigvals)[::-1]:
                        axis = np.asarray(eigvecs[:, int(eig_index)], dtype=np.float64)
                        axis = self._unit(axis)
                        if not axis.any():
                            continue
                        if any(abs(float(axis @ old)) > 0.92 for old in fallback_axes):
                            continue
                        fallback_axes.append(axis)
                    for axis_index, axis in enumerate(fallback_axes[:8]):
                        values = grid.xyz @ axis
                        local = values[comp]
                        local_weights = grid.cell_area[comp]
                        if local.size < 8:
                            continue
                        median = self._weighted_percentile_local(
                            local, local_weights, 50.0)
                        belt_specs = [
                            (5.0, 0.0),
                            (7.5, -8.0),
                            (7.5, 8.0),
                            (10.0, 0.0),
                        ]
                        if historical_residence_controller_active:
                            belt_specs.extend([
                                (11.0, -16.0),
                                (11.0, 16.0),
                                (13.0, 0.0),
                            ])
                            if historical_residence_debt_myr >= 450.0:
                                belt_specs.extend([
                                    (12.0, -26.0),
                                    (12.0, 26.0),
                                ])
                        for width_pct, offset_pct in belt_specs:
                            center_pct = float(np.clip(50.0 + offset_pct, 18.0, 82.0))
                            center = (
                                median if offset_pct == 0.0
                                else self._weighted_percentile_local(
                                    local, local_weights, center_pct)
                            )
                            lo = self._weighted_percentile_local(
                                local, local_weights,
                                max(1.0, center_pct - width_pct))
                            hi = self._weighted_percentile_local(
                                local, local_weights,
                                min(99.0, center_pct + width_pct))
                            if hi <= lo:
                                span = max(float(np.nanmax(local) - np.nanmin(local)),
                                           1.0e-9)
                                lo = center - 0.035 * span
                                hi = center + 0.035 * span
                            belt = (
                                historical_apron
                                & (values >= lo)
                                & (values <= hi)
                            )
                            belt = self._dilate_mask(
                                grid,
                                belt,
                                allowed=historical_apron,
                                passes=1,
                            )
                            cells = np.where(belt)[0].astype(np.int64)
                            if cells.size < 4:
                                continue
                            path_area = float(grid.cell_area[cells].sum())
                            area_share = path_area / max(comp_area, 1.0)
                            max_area_share = (
                                0.235
                                if historical_residence_controller_active
                                else 0.155
                            )
                            if area_share < 0.006 or area_share > max_area_share:
                                continue
                            largest_after, split_balance, topology_score, _ = (
                                partition_metrics(cells)
                            )
                            min_topology = (
                                0.11
                                if historical_residence_controller_active
                                else 0.16
                            )
                            min_balance = (
                                0.030
                                if historical_residence_controller_active
                                else 0.045
                            )
                            if topology_score < min_topology or split_balance < min_balance:
                                continue
                            add_candidate(
                                f"history:{axis_index}",
                                f"history:{axis_index}:{center_pct:.1f}",
                                cells,
                                length_penalty=0.004,
                            )

            candidates.sort(
                key=lambda item: (
                    -float(item.get("topology_score", 0.0)),
                    -float(item["quality"]),
                    float(item["path_area"]),
                    int(np.min(item["cells"])),
                )
            )
            telemetry["candidate_count"] = int(len(candidates))
            if candidates:
                telemetry["best_topology_score"] = float(
                    max(float(c.get("topology_score", 0.0)) for c in candidates)
                )
                telemetry["best_quality_score"] = float(
                    max(float(c.get("quality", 0.0)) for c in candidates)
                )
            used = np.zeros(grid.n, dtype=bool)
            accepted_here = 0
            max_accept_index = (
                3
                if (
                    historical_residence_controller_active
                    and historical_residence_debt_myr >= 450.0
                    and grid.n >= 5000
                )
                else 2
            )
            for local_index, cand in enumerate(candidates[:8]):
                cells = np.asarray(cand["cells"], dtype=np.int64)
                if used[cells].mean() > 0.30:
                    continue
                lineage_suffix = str(
                    cand.get("lineage_suffix", f"axis:{cand['axis_index']}")
                )
                lineage = (
                    f"continent:{parent_cids[0] if parent_cids else comp_index}:"
                    f"{lineage_suffix}"
                )
                prev = previous.get(lineage, {})
                birth = float(prev.get("birth_myr", t))
                lat, lon = self._cell_centroid_lat_lon(grid, cells)
                parent_rifts = sorted({
                    rid
                    for c in cells
                    for rid in rift_by_cell.get(int(c), [])
                })
                objects.append({
                    "id": prev.get("id", f"breakup_seaway:{lineage}"),
                    "kind": "breakup_seaway",
                    "stage": "incipient_multi_rift_seaway",
                    "lineage_key": lineage,
                    "basis": (
                        "p110b_historical_supercontinent_breakup_pressure"
                        if historical_breakup_active
                        else "p21_supercontinent_breakup_rift_weakness_path"
                    ),
                    "historical_breakup_pressure": float(
                        historical_breakup_pressure),
                    "historical_largest_continent_share": float(
                        historical_supercontinent_share),
                    "historical_residence_debt_myr": float(
                        historical_residence_debt_myr),
                    "historical_residence_controller_active": bool(
                        historical_residence_controller_active),
                    "parent_continent_ids": parent_cids,
                    "parent_plate_ids": parent_pids,
                    "parent_rift_system_ids": parent_rifts,
                    "component_share_of_continental_crust": float(comp_share_of_cont),
                    "component_area_fraction": float(comp_share_total),
                    "cell_count": int(cells.size),
                    "area_fraction": float(cand["path_area"] / max(total, 1.0)),
                    "mean_rift_potential": float(cand["mean_rift"]),
                    "weak_fraction": float(cand["weak_fraction"]),
                    "core_fraction": float(cand["core_fraction"]),
                    "boundary_fraction": float(cand["boundary_fraction"]),
                    "topology_score": float(cand.get("topology_score", 0.0)),
                    "split_balance": float(cand.get("split_balance", 0.0)),
                    "largest_after_partition": float(
                        cand.get("largest_after_partition", 1.0)),
                    "partition_count": int(cand.get("partition_count", 1)),
                    "quality_score": float(cand["quality"]),
                    "lat": lat,
                    "lon": lon,
                    "birth_myr": round(float(birth), 1),
                    "last_active_myr": round(float(t), 1),
                    "age_myr": round(float(max(t - birth, 0.0)), 1),
                    "cells": cells.astype(int).tolist(),
                })
                used[cells] = True
                accepted_here += 1
                if local_index >= max_accept_index:
                    break
            telemetry["accepted_object_count"] = int(accepted_here)
            if not candidates:
                telemetry["skip_reason"] = "no_candidate_paths"
            elif accepted_here == 0:
                telemetry["skip_reason"] = "all_candidates_overlapped_or_rejected"
            component_telemetry.append(telemetry)

        active_lineages = {
            str(obj.get("lineage_key"))
            for obj in objects
            if obj.get("lineage_key") is not None
        }
        for retained in retained_historical:
            lineage = str(retained.get("lineage_key"))
            if lineage and lineage in active_lineages:
                continue
            objects.append(retained)
            if lineage:
                active_lineages.add(lineage)
        objects.sort(
            key=lambda obj: (
                -float(obj.get("quality_score", 0.0)),
                -float(obj.get("component_area_fraction", 0.0)),
                str(obj.get("id")),
            )
        )
        world.objects["tectonics.breakup_component_telemetry"] = component_telemetry[:16]
        return objects[:12]

    def _retained_historical_breakup_seaway_objects(
        self,
        previous: list[dict],
        t: float,
        *,
        historical_breakup_active: bool,
    ) -> list[dict]:
        """Keep young historical breakup seaways available across timesteps."""
        retained: list[dict] = []
        time_now = float(t)
        if time_now < 700.0:
            return retained
        for obj in previous or []:
            if not isinstance(obj, dict):
                continue
            if not str(obj.get("basis", "")).startswith(
                "p110b_historical_supercontinent_breakup_pressure"
            ):
                continue
            cells = [int(c) for c in obj.get("cells", []) if int(c) >= 0]
            if not cells:
                continue
            last_active = float(obj.get(
                "last_active_myr",
                obj.get("birth_myr", time_now),
            ))
            retention_age = max(time_now - last_active, 0.0)
            if (
                retention_age
                > P110B_HISTORICAL_BREAKUP_SEAWAY_RETENTION_MYR
            ):
                continue
            birth = float(obj.get("birth_myr", last_active))
            carried = dict(obj)
            carried["cells"] = cells
            carried["retained_from_previous_step"] = True
            carried["retention_age_myr"] = round(float(retention_age), 1)
            carried["stage"] = (
                "active_historical_breakup_seaway"
                if historical_breakup_active
                else "inherited_rifted_seaway"
            )
            if historical_breakup_active:
                carried["last_active_myr"] = round(time_now, 1)
                carried["retention_age_myr"] = 0.0
            carried["age_myr"] = round(float(max(time_now - birth, 0.0)), 1)
            retained.append(carried)
        return retained

    def _apply_breakup_seaway_crustal_opening(
        self,
        world,
        grid,
        breakup_seaways: list[dict],
        ctype: np.ndarray,
        thick: np.ndarray,
        age: np.ndarray,
        origin: np.ndarray,
        reworked: np.ndarray,
        stability: np.ndarray,
        domain: np.ndarray,
        continent_id: np.ndarray,
        terrane_id: np.ndarray,
        province_code: np.ndarray,
        internal_block_id: np.ndarray,
        internal_block_code: np.ndarray,
        basement_age: np.ndarray,
        basement_stability: np.ndarray,
        basement_thick: np.ndarray,
        orog_age: np.ndarray,
        volc_age: np.ndarray,
        t: float,
    ) -> None:
        """Convert accepted historical breakup corridors into young oceanic crust."""
        opening = np.zeros(grid.n, dtype=bool)
        ocean_age = np.full(grid.n, np.inf, dtype=np.float64)
        time_now = float(t)
        for obj in breakup_seaways or []:
            if not isinstance(obj, dict):
                continue
            if not str(obj.get("basis", "")).startswith(
                "p110b_historical_supercontinent_breakup_pressure"
            ):
                continue
            cells = np.asarray([
                int(c) for c in obj.get("cells", [])
                if 0 <= int(c) < int(grid.n)
            ], dtype=np.int64)
            if cells.size == 0:
                continue
            local = np.zeros(grid.n, dtype=bool)
            local[cells] = True
            local &= (
                (ctype == CONT)
                & (domain != DOMAIN_CRATON)
                & (stability < 0.90)
            )
            if not local.any():
                continue
            birth = float(obj.get("birth_myr", time_now))
            local_age = max(time_now - birth, 0.0)
            opening |= local
            ocean_age[local] = np.minimum(ocean_age[local], local_age)

        if not opening.any():
            world.set_g("tectonics.last_p110b_breakup_crustal_opening_cells", 0.0)
            world.set_g(
                "tectonics.last_p110b_breakup_crustal_opening_area_fraction",
                0.0,
            )
            return

        ctype[opening] = OCEAN
        thick[opening] = np.minimum(thick[opening], OCEAN_THICK + 1200.0)
        finite_age = np.where(np.isfinite(ocean_age), ocean_age, 0.0)
        age[opening] = np.minimum(age[opening], finite_age[opening])
        origin[opening] = ORIGIN_RIDGE
        reworked[opening] = time_now
        stability[opening] = np.minimum(stability[opening], 0.18)
        domain[opening] = DOMAIN_OCEANIC
        continent_id[opening] = -1
        terrane_id[opening] = -1
        province_code[opening] = PROVINCE_NONE
        internal_block_id[opening] = 0.0
        internal_block_code[opening] = INTERNAL_BLOCK_NONE
        basement_age[opening] = 0.0
        basement_stability[opening] = 0.0
        basement_thick[opening] = 0.0
        orog_age[opening] = -1.0
        volc_age[opening] = -1.0
        world.set_g(
            "tectonics.last_p110b_breakup_crustal_opening_cells",
            float(np.count_nonzero(opening)),
        )
        world.set_g(
            "tectonics.last_p110b_breakup_crustal_opening_area_fraction",
            float(grid.cell_area[opening].sum() / max(float(grid.cell_area.sum()), 1.0)),
        )

    def _p110b_historical_supercontinent_breakup_pressure(
        self,
        world,
        grid,
        ctype,
        target: float,
        t: float,
    ) -> tuple[float, bool, float]:
        """Track prolonged over-connected continental crust through history."""
        enabled = (
            str(world.spec.name).startswith("earthlike")
            and 0.20 <= float(target) <= 0.38
            and world.g("terrain.allow_p110a_supercontinent_final_state", 0.0) < 0.5
        )
        time_now = float(t)
        previous_time = float(world.g(
            "tectonics.last_p110b_historical_supercontinent_pressure_time_myr",
            time_now,
        ))
        dt = max(time_now - previous_time, 0.0)
        cont = np.asarray(ctype, dtype=np.float64) == CONT
        cont_area = float(grid.cell_area[cont].sum())
        total = max(float(grid.cell_area.sum()), 1.0e-12)
        largest_share = 0.0
        component_count = 0
        if enabled and cont_area > 0.0:
            comps = self._connected_components(grid, cont)
            component_count = int(len(comps))
            if comps:
                largest_share = max(
                    float(grid.cell_area[comp].sum()) for comp in comps
                ) / max(cont_area, 1.0e-12)

        supercontinent_like = bool(
            enabled
            and cont_area / total >= 0.10
            and largest_share >= 0.64
            and time_now >= 350.0
        )
        previous_consecutive = float(world.g(
            "tectonics.last_p110b_historical_supercontinent_consecutive_myr", 0.0))
        previous_cumulative = float(world.g(
            "tectonics.last_p110b_historical_supercontinent_cumulative_myr", 0.0))
        if supercontinent_like:
            consecutive = previous_consecutive + dt
            cumulative = previous_cumulative + dt
        else:
            consecutive = max(0.0, previous_consecutive - 0.75 * dt)
            cumulative = previous_cumulative
        residence_debt = (
            max(0.0, consecutive - P110B_SUPERCONTINENT_RESIDENCE_DEBT_START_MYR)
            if supercontinent_like else 0.0
        )
        residence_controller_active = bool(
            enabled
            and time_now >= 700.0
            and time_now < 4100.0
            and residence_debt > 0.0
        )

        pressure = 0.0
        if enabled and time_now >= 700.0 and supercontinent_like:
            pressure = max(
                (largest_share - 0.64) / 0.24,
                (consecutive - 260.0) / 650.0,
                (cumulative - 620.0) / 1200.0,
                residence_debt / P110B_SUPERCONTINENT_RESIDENCE_FULL_PRESSURE_MYR,
            )
            pressure = float(np.clip(pressure, 0.0, 1.0))
        active = bool(pressure >= 0.30 and time_now < 3600.0)

        world.set_g(
            "tectonics.last_p110b_historical_supercontinent_pressure_time_myr",
            time_now,
        )
        world.set_g(
            "tectonics.last_p110b_historical_largest_continent_share",
            float(largest_share),
        )
        world.set_g(
            "tectonics.last_p110b_historical_continent_component_count",
            float(component_count),
        )
        world.set_g(
            "tectonics.last_p110b_historical_supercontinent_consecutive_myr",
            float(consecutive),
        )
        world.set_g(
            "tectonics.last_p110b_historical_supercontinent_cumulative_myr",
            float(cumulative),
        )
        world.set_g(
            "tectonics.last_p110b_historical_breakup_pressure",
            float(pressure),
        )
        world.set_g(
            "tectonics.last_p110b_historical_supercontinent_residence_debt_myr",
            float(residence_debt),
        )
        world.set_g(
            "tectonics.last_p110b_historical_residence_controller_active",
            1.0 if residence_controller_active else 0.0,
        )
        world.set_g(
            "tectonics.last_p110b_historical_breakup_active",
            1.0 if active else 0.0,
        )
        return float(pressure), active, float(largest_share)

    @staticmethod
    def _weighted_percentile_local(
        values: np.ndarray,
        weights: np.ndarray,
        percentile: float,
    ) -> float:
        vals = np.asarray(values, dtype=np.float64)
        w = np.asarray(weights, dtype=np.float64)
        finite = np.isfinite(vals) & np.isfinite(w) & (w > 0.0)
        if not finite.any():
            return 0.0
        vals = vals[finite]
        w = w[finite]
        order = np.argsort(vals)
        vals = vals[order]
        w = w[order]
        cdf = np.cumsum(w)
        cutoff = float(np.clip(percentile, 0.0, 100.0)) / 100.0 * float(cdf[-1])
        idx = int(np.searchsorted(cdf, cutoff, side="left"))
        return float(vals[np.clip(idx, 0, vals.size - 1)])

    def _component_boundary_mask(self, grid, mask):
        out = np.zeros(grid.n, dtype=bool)
        for c in np.where(mask)[0]:
            if np.any(~mask[grid.neighbors[int(c)]]):
                out[int(c)] = True
        return out

    def _breakup_candidate_axes(self, grid, comp, parent_plate_ids, plate_by_id):
        cells = np.asarray(comp, dtype=np.int64)
        weights = grid.cell_area[cells]
        centroid = np.average(grid.xyz[cells], axis=0, weights=weights)
        centroid = self._unit(centroid)
        axes: list[np.ndarray] = []

        def add_axis(vec):
            axis = self._unit(vec)
            if not axis.any():
                return
            axis = axis - float(axis @ centroid) * centroid
            axis = self._unit(axis)
            if not axis.any():
                return
            for old in axes:
                if abs(float(old @ axis)) > 0.92:
                    return
            axes.append(axis)

        for pid in parent_plate_ids:
            topo = plate_by_id.get(int(pid), {})
            pole = np.asarray(topo.get("pole", [0.0, 0.0, 0.0]), dtype=np.float64)
            add_axis(np.cross(centroid, pole))
        for base in (
            np.asarray([1.0, 0.0, 0.0]),
            np.asarray([0.0, 1.0, 0.0]),
            np.asarray([0.0, 0.0, 1.0]),
        ):
            add_axis(np.cross(centroid, base))
            add_axis(base)
        if not axes:
            axes.append(self._unit(np.asarray([1.0, 0.0, 0.0])))
        return axes[:6]

    def _breakup_path_between_sets(self, grid, starts, goals, allowed,
                                   cell_cost, *, max_path_cells):
        start_mask = np.asarray(starts, dtype=bool)
        goal_mask = np.asarray(goals, dtype=bool)
        allowed = np.asarray(allowed, dtype=bool)
        start_cells = np.where(start_mask & allowed)[0]
        if start_cells.size == 0 or not np.any(goal_mask & allowed):
            return []
        n = grid.n
        dist = np.full(n, np.inf, dtype=np.float64)
        prev = np.full(n, -1, dtype=np.int64)
        hops = np.zeros(n, dtype=np.int32)
        heap: list[tuple[float, int]] = []
        for s in start_cells:
            s = int(s)
            dist[s] = float(cell_cost[s])
            heapq.heappush(heap, (float(dist[s]), s))
        expanded = 0
        while heap:
            d, c = heapq.heappop(heap)
            c = int(c)
            if d != dist[c]:
                continue
            expanded += 1
            if expanded > 7000:
                break
            if goal_mask[c] and not start_mask[c] and hops[c] >= 3:
                path: list[int] = []
                cur = c
                while cur >= 0:
                    path.append(int(cur))
                    cur = int(prev[cur])
                path.reverse()
                return path if len(path) <= max_path_cells else []
            if hops[c] >= max_path_cells:
                continue
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if not allowed[nb]:
                    continue
                nd = d + float(cell_cost[nb]) + 1e-8 * nb
                if nd < dist[nb]:
                    dist[nb] = nd
                    prev[nb] = c
                    hops[nb] = hops[c] + 1
                    heapq.heappush(heap, (nd, nb))
        return []

    @staticmethod
    def _cell_centroid_lat_lon(grid, cells):
        cells = np.asarray(cells, dtype=np.int64)
        if cells.size == 0:
            return 0.0, 0.0
        weights = grid.cell_area[cells]
        centroid = np.average(grid.xyz[cells], axis=0, weights=weights)
        centroid = centroid / max(float(np.linalg.norm(centroid)), 1e-12)
        lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
        lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
        return round(lat, 2), round(lon, 2)

    def _wilson_cycle_objects_from_basins(self, ocean_basins):
        cycles = []
        for basin in ocean_basins:
            cycles.append({
                "id": f"wilson:{basin['id']}",
                "stage": basin["stage"],
                "phase_code": float(basin["phase_code"]),
                "boundary_object_id": (
                    basin["boundary_object_ids"][0]
                    if basin.get("boundary_object_ids") else None
                ),
                "boundary_object_ids": basin.get("boundary_object_ids", []),
                "ocean_basin_id": basin["id"],
                "lineage_key": basin["lineage_key"],
                "kind": "ocean_basin_lifecycle",
                "cells": basin.get("cells", []),
                "plate_ids": basin.get("parent_plate_ids", []),
                "last_transition_myr": basin.get("last_active_myr"),
                "age_myr": basin.get("age_myr", 0.0),
            })
        return cycles

    def _ocean_gateway_objects_from_lifecycle(self, boundary_objects, wilson_cycles,
                                              ocean_basins, t):
        cycle_by_boundary: dict[str, dict] = {}
        for cycle in wilson_cycles:
            for bid in cycle.get("boundary_object_ids", []):
                cycle_by_boundary[str(bid)] = cycle
        basin_by_id = {b["id"]: b for b in ocean_basins}
        status_by_phase = {
            WILSON_OPENING: "opening",
            WILSON_MATURE: "open",
            WILSON_CLOSING: "closing",
            WILSON_COLLISION: "restricted",
            WILSON_SUTURE: "closed",
        }
        out = []
        for obj in boundary_objects:
            cycle = cycle_by_boundary.get(str(obj.get("id")))
            if cycle is None:
                continue
            phase_code = float(cycle.get("phase_code", WILSON_INACTIVE))
            if phase_code == WILSON_INACTIVE:
                continue
            cells = obj.get("cells", [])
            if len(cells) < 3:
                continue
            basin = basin_by_id.get(cycle.get("ocean_basin_id"), {})
            gateway_id = f"gateway:{obj['id']}"
            out.append({
                "id": gateway_id,
                "kind": "ocean_gateway",
                "boundary_object_id": obj["id"],
                "parent_boundary_object_id": obj["id"],
                "wilson_cycle_id": cycle["id"],
                "parent_basin_id": cycle.get("ocean_basin_id"),
                "parent_rift_system_ids": basin.get("parent_rift_system_ids", []),
                "parent_margin_ids": basin.get("parent_margin_ids", []),
                "parent_closing_margin_ids": basin.get("parent_closing_margin_ids", []),
                "parent_suture_ids": basin.get("parent_suture_ids", []),
                "phase": cycle["stage"],
                "phase_code": phase_code,
                "status": status_by_phase.get(phase_code, "unknown"),
                "cells": cells,
                "plate_ids": obj.get("parent_plate_ids", []),
                "lat": obj.get("lat"),
                "lon": obj.get("lon"),
                "area_fraction": obj.get("area_fraction", 0.0),
                "birth_myr": obj.get("birth_myr", round(float(t), 1)),
                "last_active_myr": round(float(t), 1),
            })
        out.sort(key=lambda obj: float(obj.get("area_fraction", 0.0)), reverse=True)
        return out[:24]

    def _wilson_cycle_objects(self, boundary_objects, t):
        stage_by_kind = {
            "ridge": ("widening_ocean", WILSON_OPENING),
            "passive_margin": ("mature_ocean", WILSON_MATURE),
            "trench": ("closing_ocean", WILSON_CLOSING),
            "active_margin": ("arc_continent_margin", WILSON_COLLISION),
            "suture": ("suture_post_orogenic", WILSON_SUTURE),
        }
        cycles = []
        for obj in boundary_objects:
            kind = obj["kind"]
            if kind not in stage_by_kind:
                continue
            stage, phase_code = stage_by_kind[kind]
            cycles.append({
                "id": f"wilson:{obj['id']}",
                "stage": stage,
                "phase_code": float(phase_code),
                "boundary_object_id": obj["id"],
                "kind": kind,
                "cells": obj["cells"],
                "plate_ids": obj["parent_plate_ids"],
                "last_transition_myr": round(float(t), 1),
            })
        return cycles

    def _wilson_phase_field(self, grid, wilson_cycles):
        phase = np.full(grid.n, WILSON_INACTIVE, dtype=np.float64)
        priority = {
            WILSON_RIFT: 1,
            WILSON_OPENING: 2,
            WILSON_MATURE: 3,
            WILSON_CLOSING: 4,
            WILSON_COLLISION: 5,
            WILSON_SUTURE: 6,
        }
        current = np.zeros(grid.n, dtype=np.int16)
        for cycle in wilson_cycles:
            cells = np.asarray(cycle.get("cells", []), dtype=int)
            cells = cells[(0 <= cells) & (cells < grid.n)]
            if cells.size == 0:
                continue
            code = float(cycle.get("phase_code", WILSON_INACTIVE))
            rank = priority.get(code, 0)
            replace = rank >= current[cells]
            if replace.any():
                target = cells[replace]
                phase[target] = code
                current[target] = rank
        return phase

    def _ocean_gateway_objects(self, boundary_objects, wilson_cycles, t):
        cycle_by_boundary = {w["boundary_object_id"]: w for w in wilson_cycles}
        status_by_phase = {
            WILSON_OPENING: "opening",
            WILSON_MATURE: "open",
            WILSON_CLOSING: "closing",
            WILSON_COLLISION: "restricted",
            WILSON_SUTURE: "closed",
        }
        out = []
        for obj in boundary_objects:
            cycle = cycle_by_boundary.get(obj["id"])
            if cycle is None:
                continue
            phase_code = float(cycle.get("phase_code", WILSON_INACTIVE))
            if phase_code == WILSON_INACTIVE:
                continue
            cells = obj.get("cells", [])
            if len(cells) < 3:
                continue
            gateway_id = f"gateway:{obj['id']}"
            out.append({
                "id": gateway_id,
                "kind": "ocean_gateway",
                "boundary_object_id": obj["id"],
                "wilson_cycle_id": cycle["id"],
                "phase": cycle["stage"],
                "phase_code": phase_code,
                "status": status_by_phase.get(phase_code, "unknown"),
                "cells": cells,
                "plate_ids": obj.get("parent_plate_ids", []),
                "lat": obj.get("lat"),
                "lon": obj.get("lon"),
                "area_fraction": obj.get("area_fraction", 0.0),
                "birth_myr": obj.get("birth_myr", round(float(t), 1)),
                "last_active_myr": round(float(t), 1),
            })
        out.sort(key=lambda obj: float(obj.get("area_fraction", 0.0)), reverse=True)
        return out[:24]

    def _should_emit_wilson_events(self, t, dt):
        if t < 100.0:
            return False
        return int((t - dt) // 250.0) != int(t // 250.0)

    def _wilson_events_from_objects(self, boundary_objects, wilson_cycles,
                                    ocean_gateways, t):
        event_by_kind = {
            "ridge": ("ridge_birth", "rift_birth", "ocean_basin_opening",
                      "ocean_gateway_opened"),
            "passive_margin": ("passive_margin_formation", "ocean_basin_maturation",
                               "ocean_gateway_opened"),
            "trench": ("trench_birth", "subduction_initiation",
                       "ocean_basin_closure", "ocean_gateway_closed"),
            "active_margin": ("arc_birth", "subduction_initiation",
                              "terrane_accretion"),
            "suture": ("suture_formation", "suture_formed", "continent_collision",
                       "orogen_built", "ocean_gateway_closed"),
        }
        gateway_by_boundary = {
            g["boundary_object_id"]: g for g in ocean_gateways
        }
        out = []
        seen = set()
        for obj in boundary_objects:
            kinds = event_by_kind.get(obj["kind"], ())
            if not kinds or obj["kind"] in seen:
                continue
            seen.add(obj["kind"])
            location = obj["cells"][0] if obj["cells"] else None
            for type_ in kinds:
                out.append(Event(
                    type_, t, self.name, location=location,
                    magnitude=float(obj["area_fraction"]),
                    params={
                        "boundary_object_id": obj["id"],
                        "ocean_gateway_id": gateway_by_boundary.get(obj["id"], {}).get("id"),
                        "stage": obj["stage"],
                        "plate_ids": obj["parent_plate_ids"],
                        "lat": obj["lat"],
                        "lon": obj["lon"],
                    },
                ))
        if wilson_cycles:
            mature = next((w for w in wilson_cycles if w["stage"] == "mature_ocean"), None)
            if mature is not None:
                out.append(Event(
                    "ocean_basin_maturation", t, self.name,
                    location=mature["cells"][0] if mature["cells"] else None,
                    magnitude=float(len(mature["cells"])),
                    params={"wilson_cycle_id": mature["id"], "stage": mature["stage"]},
                ))
        return out

    def _craton_objects(self, grid, craton_mask, age, t):
        comps = self._connected_components(grid, np.asarray(craton_mask, dtype=bool))
        comps.sort(key=lambda c: float(grid.cell_area[c].sum()), reverse=True)
        out = []
        for n, comp in enumerate(comps[:16]):
            if comp.size < 2:
                continue
            out.append({
                "id": f"craton:{int(comp.min())}:{n}",
                "kind": "craton",
                "cells": comp.astype(int).tolist(),
                "area_fraction": float(grid.cell_area[comp].sum() / grid.cell_area.sum()),
                "oldest_age_myr": round(float(age[comp].max()), 1),
                "last_active_myr": round(float(t), 1),
            })
        return out

    def _platform_subsidence_field(self, world, grid, ctype, age, thick, origin,
                                   reworked, stability, domain, t):
        """Long-lived mature-platform accommodation proxy.

        The field is a deterministic tectonic state, not terrain paint.  It
        favors broad, old, quiet continental interiors with relatively thin crust
        and moderate stability, persists weakly from the previous step, and is
        suppressed on shields, sutures, margins, LIPs and accreted terranes.
        """
        cont = np.asarray(ctype, dtype=np.float64) == CONT
        zeros = np.zeros(grid.n, dtype=np.float64)
        world.set_g("tectonics.last_p49_platform_subsidence_area_fraction", 0.0)
        world.set_g("tectonics.last_p49_platform_subsidence_mean", 0.0)
        world.set_g("tectonics.last_p49_platform_subsidence_max", 0.0)
        if not cont.any():
            return zeros

        area = grid.cell_area
        width = _graph_width_steps(grid, cont)
        tt = max(float(t), 0.0)
        quiet_age = np.where(
            reworked < 0.0,
            tt,
            np.maximum(tt - np.asarray(reworked, dtype=np.float64), 0.0),
        )
        quiet_cut = min(1200.0, max(450.0, 0.18 * tt))
        mature_cut = min(2200.0, max(900.0, 0.36 * tt))
        active_or_accreted = np.isin(
            domain,
            [DOMAIN_CONTINENTAL_MARGIN, DOMAIN_ACCRETED_TERRANE,
             DOMAIN_SUTURE, DOMAIN_LIP],
        )
        shield = (
            (origin == ORIGIN_CRATON)
            | (domain == DOMAIN_CRATON)
            | ((stability >= 0.84) & (age >= 1800.0))
        )
        platform = (
            cont
            & (width >= 4.0)
            & (age >= mature_cut)
            & (quiet_age >= quiet_cut)
            & (stability >= 0.34)
            & ~active_or_accreted
            & ~shield
        )
        if not platform.any():
            return zeros

        platform_thick = thick[platform]
        thick_low = float(np.percentile(platform_thick, 25.0))
        thick_high = max(thick_low + 1800.0, float(np.percentile(platform_thick, 75.0)))
        thin = np.clip((thick_high - thick) / max(thick_high - thick_low, 1.0), 0.0, 1.0)
        stability_window = np.clip((0.88 - stability) / 0.48, 0.0, 1.0)
        width_maturity = np.clip((width - 3.0) / 5.0, 0.0, 1.0)
        age_maturity = np.clip((age - mature_cut) / max(3200.0 - mature_cut, 1.0), 0.0, 1.0)
        quiet_maturity = np.clip(quiet_age / max(quiet_cut * 1.8, 1.0), 0.0, 1.0)

        candidate = np.zeros(grid.n, dtype=np.float64)
        candidate[platform] = (
            0.48 * thin[platform]
            + 0.27 * stability_window[platform]
            + 0.10 * width_maturity[platform]
            + 0.08 * age_maturity[platform]
            + 0.07 * quiet_maturity[platform]
        )
        smooth = self._smooth_field(
            grid, candidate, allowed=platform, passes=1, alpha=0.22)
        candidate[platform] = 0.76 * candidate[platform] + 0.24 * smooth[platform]

        prev = np.asarray(world.get_field("tectonics.platform_subsidence", 0.0),
                          dtype=np.float64)
        if prev.shape == (grid.n,):
            field = np.maximum(candidate, 0.62 * prev)
        else:
            field = candidate
        field[~platform] = 0.0
        field = np.clip(field, 0.0, 1.0)

        active = cont & (field >= 0.38)
        cont_area = max(float(area[cont].sum()), 1e-12)
        world.set_g(
            "tectonics.last_p49_platform_subsidence_area_fraction",
            float(area[active].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p49_platform_subsidence_mean",
            float(np.average(field[active], weights=area[active])) if active.any() else 0.0,
        )
        world.set_g("tectonics.last_p49_platform_subsidence_max", float(field.max()))
        return field

    def _update_cratonic_province_memory(self, world, grid, ctype, origin, age,
                                         reworked, stability, domain, previous_code,
                                         province_objects, t):
        """Persist shield/platform province cargo across rasterized motion.

        The object layer may temporarily miss a mature province when an old
        continent is reworked at its edge.  This field carries the province as
        crustal cargo and lets later steps re-identify it from inherited
        context, while clearing memory from oceanic and active accreted cells.
        """
        cont = np.asarray(ctype, dtype=np.float64) == CONT
        code = np.asarray(previous_code, dtype=np.float64).copy()
        if code.shape != (grid.n,):
            code = np.zeros(grid.n, dtype=np.float64)
        code[~cont] = PROVINCE_NONE
        width = _graph_width_steps(grid, cont)
        quiet_age = np.where(
            np.asarray(reworked, dtype=np.float64) < 0.0,
            float(t),
            np.maximum(float(t) - np.asarray(reworked, dtype=np.float64), 0.0),
        )
        active_or_accreted = np.isin(
            domain,
            [DOMAIN_CONTINENTAL_MARGIN, DOMAIN_ACCRETED_TERRANE,
             DOMAIN_SUTURE, DOMAIN_LIP],
        )
        stale = (
            cont
            & (code > PROVINCE_NONE)
            & (
                (width < 2.0)
                | ((quiet_age < 180.0) & active_or_accreted & (stability < 0.68))
            )
        )
        code[stale] = PROVINCE_NONE
        craton = cont & (origin == ORIGIN_CRATON) & (width >= 2.0)
        code[craton] = PROVINCE_SHIELD

        object_sets = (
            ("tectonics.shields", PROVINCE_SHIELD),
            ("tectonics.platforms", PROVINCE_PLATFORM),
            ("tectonics.interior_basins", PROVINCE_INTERIOR_BASIN),
        )
        for object_set, value in object_sets:
            for obj in province_objects.get(object_set, []):
                cells = np.asarray(obj.get("cells", []), dtype=int)
                cells = cells[(0 <= cells) & (cells < grid.n)]
                if cells.size:
                    target = cells[cont[cells]]
                    if value != PROVINCE_SHIELD:
                        target = target[code[target] != PROVINCE_SHIELD]
                    code[target] = value

        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1e-12)
        world.set_g(
            "tectonics.last_p54_memory_shield_share_of_continental_crust",
            float(area[cont & (code == PROVINCE_SHIELD)].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p54_memory_platform_share_of_continental_crust",
            float(area[cont & (code == PROVINCE_PLATFORM)].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p54_memory_basin_share_of_continental_crust",
            float(area[cont & (code == PROVINCE_INTERIOR_BASIN)].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p54_memory_total_share_of_continental_crust",
            float(area[cont & (code > PROVINCE_NONE)].sum() / cont_area),
        )
        return code

    def _p58_old_basement_mask(self, ctype, basement_age, basement_stability,
                               basement_thick):
        ctype_arr = np.asarray(ctype, dtype=np.float64)
        age_floor = np.asarray(basement_age, dtype=np.float64)
        stability_floor = np.asarray(basement_stability, dtype=np.float64)
        thick_floor = np.asarray(basement_thick, dtype=np.float64)
        return (
            (ctype_arr == CONT)
            & (age_floor >= 2500.0)
            & (stability_floor >= 0.70)
            & (thick_floor > 0.0)
        )

    def _apply_p59_old_basement_raster_coverage(
        self, world, grid, primary, idx, dist, plate, ctype, origin, stability,
        basement_age, basement_stability, basement_thick, t,
    ):
        """Keep broad old-basement parcels represented during raster remap.

        The fixed-grid remap asks each output cell for its nearest advected
        parcel.  That can leave some old stable basement source parcels with no
        output representative.  For broad old basement, this is a numerical
        remap artifact rather than geological recycling, so one nearby output
        cell is reassigned to the unrepresented source unless doing so would
        orphan another old-basement source.
        """
        telemetry = (
            "tectonics.last_p59_candidate_old_source_cells",
            "tectonics.last_p59_unrepresented_old_source_cells_before",
            "tectonics.last_p59_unrepresented_old_source_cells_after",
            "tectonics.last_p59_unrepresented_old_source_area_fraction_before",
            "tectonics.last_p59_unrepresented_old_source_area_fraction_after",
            "tectonics.last_p59_raster_reassigned_output_cells",
            "tectonics.last_p59_same_plate_reassigned_output_cells",
            "tectonics.last_p59_duplicate_old_source_stolen_cells",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        if world.g("tectonics.enable_p59_old_basement_raster_coverage", 0.0) <= 0.0:
            arr = np.asarray(primary, dtype=int)
            return arr, arr.copy()

        primary_arr = np.asarray(primary, dtype=int).copy()
        basement_primary = primary_arr.copy()
        idx_arr = np.asarray(idx, dtype=int)
        dist_arr = np.asarray(dist, dtype=np.float64)
        plate_arr = np.asarray(plate, dtype=int)
        if (
            primary_arr.shape != (grid.n,)
            or idx_arr.ndim != 2
            or idx_arr.shape[0] != grid.n
            or dist_arr.shape != idx_arr.shape
            or plate_arr.shape != (grid.n,)
        ):
            return primary_arr, basement_primary

        source_old = self._p58_old_basement_mask(
            ctype, basement_age, basement_stability, basement_thick)
        if not source_old.any():
            return primary_arr, basement_primary

        width = _graph_width_steps(grid, source_old)
        origin_arr = np.asarray(origin, dtype=np.float64)
        stability_arr = np.asarray(stability, dtype=np.float64)
        basement_age_arr = np.asarray(basement_age, dtype=np.float64)
        basement_stability_arr = np.asarray(basement_stability, dtype=np.float64)
        broad_old = (
            source_old
            & (
                (width >= 2.0)
                | (origin_arr == ORIGIN_CRATON)
                | (stability_arr >= 0.78)
                | (basement_stability_arr >= 0.80)
            )
        )
        if not broad_old.any():
            return primary_arr, basement_primary

        counts = np.bincount(basement_primary, minlength=grid.n).astype(np.int64)
        missing_before = broad_old & (counts <= 0)
        area = grid.cell_area
        total = max(float(area.sum()), 1.0e-12)
        world.set_g(
            "tectonics.last_p59_candidate_old_source_cells",
            float(np.count_nonzero(broad_old)),
        )
        world.set_g(
            "tectonics.last_p59_unrepresented_old_source_cells_before",
            float(np.count_nonzero(missing_before)),
        )
        world.set_g(
            "tectonics.last_p59_unrepresented_old_source_area_fraction_before",
            float(area[missing_before].sum() / total),
        )
        if not missing_before.any():
            return primary_arr, basement_primary

        median_chord = float(np.median(grid.edge_lengths) / grid.radius_m)
        radius = max(
            float(world.g("tectonics.raster_radius", 0.04)) * 1.35,
            2.35 * median_chord,
        )
        output_taken = np.zeros(grid.n, dtype=bool)
        reassigned = 0
        same_plate = 0
        duplicate_old_stolen = 0
        sources = sorted(
            (int(c) for c in np.where(missing_before)[0]),
            key=lambda c: (
                -float(width[c]),
                -float(basement_stability_arr[c]),
                -float(basement_age_arr[c]),
                c,
            ),
        )
        for src in sources:
            if counts[src] > 0:
                continue
            rows, ranks = np.where(idx_arr == int(src))
            if rows.size == 0:
                continue
            candidates: list[tuple[float, int, int]] = []
            for row, rank in zip(rows.astype(int), ranks.astype(int)):
                if output_taken[row]:
                    continue
                d = float(dist_arr[row, rank])
                if not np.isfinite(d) or d > radius:
                    continue
                current = int(basement_primary[row])
                if current == src:
                    continue
                current_old = bool(source_old[current])
                if current_old and counts[current] <= 1:
                    continue
                same = int(plate_arr[current] == plate_arr[src])
                candidates.append((
                    0.0 if same else 1.0,
                    0.0 if not current_old else 1.0,
                    d,
                    row,
                ))
            if not candidates:
                continue
            candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
            _, _, _, row = candidates[0]
            current = int(basement_primary[row])
            if source_old[current]:
                duplicate_old_stolen += 1
            if plate_arr[current] == plate_arr[src]:
                same_plate += 1
            counts[current] = max(int(counts[current]) - 1, 0)
            basement_primary[row] = int(src)
            counts[src] += 1
            output_taken[row] = True
            reassigned += 1

        missing_after = broad_old & (counts <= 0)
        world.set_g(
            "tectonics.last_p59_unrepresented_old_source_cells_after",
            float(np.count_nonzero(missing_after)),
        )
        world.set_g(
            "tectonics.last_p59_unrepresented_old_source_area_fraction_after",
            float(area[missing_after].sum() / total),
        )
        world.set_g(
            "tectonics.last_p59_raster_reassigned_output_cells",
            float(reassigned),
        )
        world.set_g(
            "tectonics.last_p59_same_plate_reassigned_output_cells",
            float(same_plate),
        )
        world.set_g(
            "tectonics.last_p59_duplicate_old_source_stolen_cells",
            float(duplicate_old_stolen),
        )
        return primary_arr, basement_primary

    def _record_p58_reorg_basement_ledger(
        self, world, grid, plate_before, plate_after, ctype,
        basement_age, basement_stability, basement_thick, reorg_detail, t,
    ):
        """Record whether plate reorganization relabels old basement cargo."""
        telemetry = (
            "tectonics.last_p58_reorg_old_basement_area_fraction_before",
            "tectonics.last_p58_reorg_old_basement_area_fraction_after",
            "tectonics.last_p58_reorg_old_basement_plate_changed_area_fraction",
            "tectonics.last_p58_reorg_merge_count",
            "tectonics.last_p58_reorg_split_count",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        before = np.asarray(plate_before, dtype=int)
        after = np.asarray(plate_after, dtype=int)
        if before.shape != (grid.n,) or after.shape != (grid.n,):
            world.objects["tectonics.p58_last_reorg_basement_ledger"] = {
                "schema": "aevum.tectonics.p58_reorg_basement_ledger.v1",
                "time_myr": round(float(t), 3),
                "valid": False,
            }
            return

        old = self._p58_old_basement_mask(
            ctype, basement_age, basement_stability, basement_thick)
        area = grid.cell_area
        total = max(float(area.sum()), 1.0e-12)
        changed = old & (before != after)
        detail = reorg_detail if isinstance(reorg_detail, dict) else {}
        merges = detail.get("merged", []) or []
        splits = detail.get("split", []) or []
        rows: list[dict] = []
        for merge in merges:
            source = int(merge.get("from", -1))
            target = int(merge.get("to", -1))
            src_old = old & (before == source)
            rows.append({
                "kind": "merge",
                "from": source,
                "to": target,
                "source_old_area_fraction": float(area[src_old].sum() / total),
                "source_old_cells": int(np.count_nonzero(src_old)),
            })
        for split in splits:
            parent = int(split.get("parent", -1))
            child = int(split.get("new_id", -1))
            child_old = old & (after == child)
            parent_old_after = old & (after == parent)
            rows.append({
                "kind": "split",
                "parent": parent,
                "new_id": child,
                "child_old_area_fraction": float(area[child_old].sum() / total),
                "parent_old_area_fraction_after": float(
                    area[parent_old_after].sum() / total),
                "child_old_cells": int(np.count_nonzero(child_old)),
            })

        old_area = float(area[old].sum() / total)
        changed_area = float(area[changed].sum() / total)
        world.set_g(
            "tectonics.last_p58_reorg_old_basement_area_fraction_before",
            old_area,
        )
        world.set_g(
            "tectonics.last_p58_reorg_old_basement_area_fraction_after",
            old_area,
        )
        world.set_g(
            "tectonics.last_p58_reorg_old_basement_plate_changed_area_fraction",
            changed_area,
        )
        world.set_g("tectonics.last_p58_reorg_merge_count", float(len(merges)))
        world.set_g("tectonics.last_p58_reorg_split_count", float(len(splits)))
        world.objects["tectonics.p58_last_reorg_basement_ledger"] = {
            "schema": "aevum.tectonics.p58_reorg_basement_ledger.v1",
            "time_myr": round(float(t), 3),
            "valid": True,
            "old_basement_area_fraction": old_area,
            "old_basement_plate_changed_area_fraction": changed_area,
            "merge_count": int(len(merges)),
            "split_count": int(len(splits)),
            "rows": rows,
        }

    def _record_p58_raster_basement_ledger(
        self, world, grid, plate, ctype, basement_age, basement_stability,
        basement_thick, primary, ridge, oceanic_arc,
        post_boundary_ctype, post_boundary_basement_age,
        post_boundary_basement_stability, post_boundary_basement_thick,
        pre_repair_ctype, pre_repair_basement_age,
        pre_repair_basement_stability, pre_repair_basement_thick,
        post_repair_ctype, post_repair_basement_age,
        post_repair_basement_stability, post_repair_basement_thick, t,
    ):
        """Record old-basement survival through raster remap and repair stages."""
        telemetry = (
            "tectonics.last_p58_pre_step_old_basement_area_fraction",
            "tectonics.last_p58_raster_output_old_basement_area_fraction",
            "tectonics.last_p58_raster_unique_source_old_basement_area_fraction",
            "tectonics.last_p58_raster_unrepresented_old_basement_area_fraction",
            "tectonics.last_p58_raster_duplicate_old_basement_area_fraction",
            "tectonics.last_p58_ridge_old_basement_cleared_area_fraction",
            "tectonics.last_p58_arc_old_basement_cleared_area_fraction",
            "tectonics.last_p58_post_boundary_old_basement_area_fraction",
            "tectonics.last_p58_pre_repair_old_basement_area_fraction",
            "tectonics.last_p58_post_repair_old_basement_area_fraction",
            "tectonics.last_p58_conservation_delta_old_basement_area_fraction",
            "tectonics.last_p58_repair_gain_old_basement_area_fraction",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        primary_arr = np.asarray(primary, dtype=int)
        plate_arr = np.asarray(plate, dtype=int)
        if primary_arr.shape != (grid.n,) or plate_arr.shape != (grid.n,):
            world.objects["tectonics.p58_last_raster_basement_ledger"] = {
                "schema": "aevum.tectonics.p58_raster_basement_ledger.v1",
                "time_myr": round(float(t), 3),
                "valid": False,
            }
            return

        area = grid.cell_area
        total = max(float(area.sum()), 1.0e-12)
        source_old = self._p58_old_basement_mask(
            ctype, basement_age, basement_stability, basement_thick)
        copied_old = source_old[primary_arr]
        unique_sources = np.unique(primary_arr[copied_old])
        unique_sources = unique_sources[source_old[unique_sources]]
        represented_source = np.zeros(grid.n, dtype=bool)
        represented_source[unique_sources] = True
        source_area = float(area[source_old].sum())
        raster_output_area = float(area[copied_old].sum())
        represented_area = float(area[represented_source].sum())
        unrepresented_area = max(source_area - represented_area, 0.0)
        duplicate_area = max(raster_output_area - represented_area, 0.0)

        ridge_mask = np.asarray(ridge, dtype=bool)
        if ridge_mask.shape != (grid.n,):
            ridge_mask = np.zeros(grid.n, dtype=bool)
        arc_mask = np.asarray(oceanic_arc, dtype=bool)
        if arc_mask.shape != (grid.n,):
            arc_mask = np.zeros(grid.n, dtype=bool)
        ridge_cleared = copied_old & ridge_mask
        arc_cleared = copied_old & arc_mask
        post_boundary_old = self._p58_old_basement_mask(
            post_boundary_ctype, post_boundary_basement_age,
            post_boundary_basement_stability, post_boundary_basement_thick)
        pre_repair_old = self._p58_old_basement_mask(
            pre_repair_ctype, pre_repair_basement_age,
            pre_repair_basement_stability, pre_repair_basement_thick)
        post_repair_old = self._p58_old_basement_mask(
            post_repair_ctype, post_repair_basement_age,
            post_repair_basement_stability, post_repair_basement_thick)
        post_boundary_area = float(area[post_boundary_old].sum())
        pre_repair_area = float(area[pre_repair_old].sum())
        post_repair_area = float(area[post_repair_old].sum())

        source_rows: list[dict] = []
        for pid in sorted(int(x) for x in np.unique(plate_arr[source_old])):
            src_mask = source_old & (plate_arr == pid)
            out_mask = copied_old & (plate_arr[primary_arr] == pid)
            src_unique = np.unique(primary_arr[out_mask])
            src_unique = src_unique[source_old[src_unique]]
            represented = np.zeros(grid.n, dtype=bool)
            represented[src_unique] = True
            src_area = float(area[src_mask].sum())
            src_represented = float(area[represented & src_mask].sum())
            source_rows.append({
                "plate_id": int(pid),
                "source_old_area_fraction": float(src_area / total),
                "raster_output_old_area_fraction": float(
                    area[out_mask].sum() / total),
                "represented_source_old_area_fraction": float(src_represented / total),
                "unrepresented_source_old_area_fraction": float(
                    max(src_area - src_represented, 0.0) / total),
            })

        world.set_g(
            "tectonics.last_p58_pre_step_old_basement_area_fraction",
            float(source_area / total),
        )
        world.set_g(
            "tectonics.last_p58_raster_output_old_basement_area_fraction",
            float(raster_output_area / total),
        )
        world.set_g(
            "tectonics.last_p58_raster_unique_source_old_basement_area_fraction",
            float(represented_area / total),
        )
        world.set_g(
            "tectonics.last_p58_raster_unrepresented_old_basement_area_fraction",
            float(unrepresented_area / total),
        )
        world.set_g(
            "tectonics.last_p58_raster_duplicate_old_basement_area_fraction",
            float(duplicate_area / total),
        )
        world.set_g(
            "tectonics.last_p58_ridge_old_basement_cleared_area_fraction",
            float(area[ridge_cleared].sum() / total),
        )
        world.set_g(
            "tectonics.last_p58_arc_old_basement_cleared_area_fraction",
            float(area[arc_cleared].sum() / total),
        )
        world.set_g(
            "tectonics.last_p58_post_boundary_old_basement_area_fraction",
            float(post_boundary_area / total),
        )
        world.set_g(
            "tectonics.last_p58_pre_repair_old_basement_area_fraction",
            float(pre_repair_area / total),
        )
        world.set_g(
            "tectonics.last_p58_post_repair_old_basement_area_fraction",
            float(post_repair_area / total),
        )
        world.set_g(
            "tectonics.last_p58_conservation_delta_old_basement_area_fraction",
            float((pre_repair_area - post_boundary_area) / total),
        )
        world.set_g(
            "tectonics.last_p58_repair_gain_old_basement_area_fraction",
            float((post_repair_area - pre_repair_area) / total),
        )
        world.objects["tectonics.p58_last_raster_basement_ledger"] = {
            "schema": "aevum.tectonics.p58_raster_basement_ledger.v1",
            "time_myr": round(float(t), 3),
            "valid": True,
            "pre_step_old_basement_area_fraction": float(source_area / total),
            "raster_output_old_basement_area_fraction": float(
                raster_output_area / total),
            "represented_source_old_basement_area_fraction": float(
                represented_area / total),
            "unrepresented_source_old_basement_area_fraction": float(
                unrepresented_area / total),
            "duplicate_old_basement_area_fraction": float(duplicate_area / total),
            "ridge_cleared_old_basement_area_fraction": float(
                area[ridge_cleared].sum() / total),
            "arc_cleared_old_basement_area_fraction": float(
                area[arc_cleared].sum() / total),
            "post_boundary_old_basement_area_fraction": float(
                post_boundary_area / total),
            "pre_repair_old_basement_area_fraction": float(pre_repair_area / total),
            "post_repair_old_basement_area_fraction": float(post_repair_area / total),
            "conservation_delta_old_basement_area_fraction": float(
                (pre_repair_area - post_boundary_area) / total),
            "repair_gain_old_basement_area_fraction": float(
                (post_repair_area - pre_repair_area) / total),
            "source_plates": source_rows,
        }

    def _repair_basement_cargo_after_conservation(
        self, world, grid, ctype, age, thick, stability,
        basement_age, basement_stability, basement_thick, t,
    ):
        """Repair basement archive after conservation and recycling decisions.

        P56 added persistent basement floors.  P57 makes the conservation path
        explicit: true oceanic cells lose basement archive, while restored
        continental cells first inherit adjacent continental archive before
        falling back to their current, young crust state.
        """
        telemetry = (
            "tectonics.last_p57_true_ocean_basement_cleared_cells",
            "tectonics.last_p57_missing_basement_candidate_cells",
            "tectonics.last_p57_neighbor_basement_cargo_inherited_cells",
            "tectonics.last_p57_neighbor_old_basement_cargo_inherited_cells",
            "tectonics.last_p57_initialized_basement_cells",
            "tectonics.last_p57_repaired_old_basement_share_of_continental_crust",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        cont = np.asarray(ctype, dtype=np.float64) == CONT
        age_floor = np.asarray(basement_age, dtype=np.float64).copy()
        stability_floor = np.asarray(basement_stability, dtype=np.float64).copy()
        thick_floor = np.asarray(basement_thick, dtype=np.float64).copy()
        if age_floor.shape != (grid.n,):
            age_floor = np.where(cont, age, 0.0).astype(np.float64)
        if stability_floor.shape != (grid.n,):
            stability_floor = np.where(cont, stability, 0.0).astype(np.float64)
        if thick_floor.shape != (grid.n,):
            thick_floor = np.where(cont, thick, 0.0).astype(np.float64)

        non_cont = ~cont
        stale_ocean = non_cont & (
            (age_floor > 0.0) | (stability_floor > 0.0) | (thick_floor > 0.0)
        )
        world.set_g(
            "tectonics.last_p57_true_ocean_basement_cleared_cells",
            float(np.count_nonzero(stale_ocean)),
        )
        age_floor[non_cont] = 0.0
        stability_floor[non_cont] = 0.0
        thick_floor[non_cont] = 0.0

        missing = cont & (thick_floor <= 0.0)
        world.set_g(
            "tectonics.last_p57_missing_basement_candidate_cells",
            float(np.count_nonzero(missing)),
        )
        copied = np.zeros(grid.n, dtype=bool)
        copied_old = np.zeros(grid.n, dtype=bool)
        for _ in range(2):
            pending = np.where(cont & (thick_floor <= 0.0))[0]
            if pending.size == 0:
                break
            copied_this_pass = 0
            for c in pending:
                nbrs = np.asarray(grid.neighbors[int(c)], dtype=int)
                donors = nbrs[(cont[nbrs]) & (thick_floor[nbrs] > 0.0)]
                if donors.size == 0:
                    continue
                best = int(donors[np.argmax(
                    age_floor[donors]
                    + 900.0 * stability_floor[donors]
                    + 0.00002 * thick_floor[donors]
                )])
                age_floor[int(c)] = max(age_floor[int(c)], min(age_floor[best], float(t)))
                stability_floor[int(c)] = max(
                    stability_floor[int(c)],
                    min(0.96 * stability_floor[best], 0.88),
                )
                thick_floor[int(c)] = max(
                    thick_floor[int(c)],
                    min(thick_floor[best], thick[int(c)]),
                )
                copied[int(c)] = True
                copied_old[int(c)] = (
                    age_floor[int(c)] >= 2500.0
                    and stability_floor[int(c)] >= 0.70
                )
                copied_this_pass += 1
            if copied_this_pass == 0:
                break

        still_missing = cont & (thick_floor <= 0.0)
        age_floor[still_missing] = np.minimum(age[still_missing], float(t))
        stability_floor[still_missing] = stability[still_missing]
        thick_floor[still_missing] = thick[still_missing]

        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1.0e-12)
        repaired_old = cont & (age_floor >= 2500.0) & (stability_floor >= 0.70)
        world.set_g(
            "tectonics.last_p57_neighbor_basement_cargo_inherited_cells",
            float(np.count_nonzero(copied)),
        )
        world.set_g(
            "tectonics.last_p57_neighbor_old_basement_cargo_inherited_cells",
            float(np.count_nonzero(copied_old)),
        )
        world.set_g(
            "tectonics.last_p57_initialized_basement_cells",
            float(np.count_nonzero(still_missing)),
        )
        world.set_g(
            "tectonics.last_p57_repaired_old_basement_share_of_continental_crust",
            float(area[repaired_old].sum() / cont_area),
        )
        return age_floor, stability_floor, thick_floor

    def _update_basement_cargo_archive(
        self, world, grid, ctype, origin, age, thick, reworked, stability, domain,
        province_code, province_objects, basement_age, basement_stability,
        basement_thick, t, *, active_boundary_mask=None,
    ):
        """Maintain persistent cratonic/platform basement cargo floors.

        The current crust state can be reworked by collision, arcs, or local
        rasterization, but broad old continental interiors should keep a basement
        archive that later steps can use to recover shield/platform identity.
        """
        telemetry = (
            "tectonics.last_p56_basement_archive_old_share_of_continental_crust",
            "tectonics.last_p56_basement_archive_shield_share_of_continental_crust",
            "tectonics.last_p56_basement_archive_platform_share_of_continental_crust",
            "tectonics.last_p56_basement_archive_active_blocked_share_of_continental_crust",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        cont = np.asarray(ctype, dtype=np.float64) == CONT
        age_floor = np.asarray(basement_age, dtype=np.float64).copy()
        stability_floor = np.asarray(basement_stability, dtype=np.float64).copy()
        thick_floor = np.asarray(basement_thick, dtype=np.float64).copy()
        if age_floor.shape != (grid.n,):
            age_floor = np.where(cont, age, 0.0).astype(np.float64)
        if stability_floor.shape != (grid.n,):
            stability_floor = np.where(cont, stability, 0.0).astype(np.float64)
        if thick_floor.shape != (grid.n,):
            thick_floor = np.where(cont, thick, 0.0).astype(np.float64)

        age_floor[~cont] = 0.0
        stability_floor[~cont] = 0.0
        thick_floor[~cont] = 0.0
        if not cont.any():
            return age_floor, stability_floor, thick_floor

        tt = max(float(t), 0.0)
        age_floor[cont] = np.clip(np.maximum(age_floor[cont], age[cont]), 0.0, tt)
        stability_floor[cont] = np.clip(
            np.maximum(stability_floor[cont], np.minimum(stability[cont], 0.86)),
            0.0,
            1.0,
        )
        thick_floor[cont] = np.maximum(thick_floor[cont], np.minimum(thick[cont], MAX_CONT_THICK))

        width = _graph_width_steps(grid, cont)
        quiet_age = np.where(
            np.asarray(reworked, dtype=np.float64) < 0.0,
            tt,
            np.maximum(tt - np.asarray(reworked, dtype=np.float64), 0.0),
        )
        code = np.asarray(province_code, dtype=np.float64)
        if code.shape != (grid.n,):
            code = np.zeros(grid.n, dtype=np.float64)

        shield_obj = np.zeros(grid.n, dtype=bool)
        platform_obj = np.zeros(grid.n, dtype=bool)
        basin_obj = np.zeros(grid.n, dtype=bool)
        for obj in province_objects.get("tectonics.shields", []):
            cells = np.asarray(obj.get("cells", []), dtype=int)
            cells = cells[(0 <= cells) & (cells < grid.n)]
            shield_obj[cells] = True
        for obj in province_objects.get("tectonics.platforms", []):
            cells = np.asarray(obj.get("cells", []), dtype=int)
            cells = cells[(0 <= cells) & (cells < grid.n)]
            platform_obj[cells] = True
        for obj in province_objects.get("tectonics.interior_basins", []):
            cells = np.asarray(obj.get("cells", []), dtype=int)
            cells = cells[(0 <= cells) & (cells < grid.n)]
            basin_obj[cells] = True
        platform_obj &= ~shield_obj
        basin_obj &= ~shield_obj

        active_core = np.zeros(grid.n, dtype=bool)
        if active_boundary_mask is not None:
            active = np.asarray(active_boundary_mask, dtype=bool)
            if active.shape == (grid.n,) and active.any():
                active_core = active & cont
        active_zone = self._dilate_mask(grid, active_core, allowed=cont, passes=1)
        active_or_accreted = np.isin(
            domain,
            [DOMAIN_CONTINENTAL_MARGIN, DOMAIN_ACCRETED_TERRANE,
             DOMAIN_SUTURE, DOMAIN_LIP],
        )
        shield_age_floor = min(tt, max(2500.0, 0.58 * tt))
        platform_age_floor = min(tt, max(2500.0, 0.52 * tt))
        old_basement = cont & (age_floor >= platform_age_floor)
        active_blocked = (
            cont
            & (active_zone | active_or_accreted)
            & (width < 4.0)
            & (quiet_age < 360.0)
            & (stability < 0.62)
        )

        shield_candidate = (
            cont
            & (width >= 2.0)
            & ~active_blocked
            & (
                (code == PROVINCE_SHIELD)
                | shield_obj
                | ((origin == ORIGIN_CRATON) & (width >= 2.0))
                | (
                    (age_floor >= shield_age_floor)
                    & (stability_floor >= 0.78)
                    & (thick_floor >= CONT_THICK + 3800.0)
                )
            )
        )
        platform_candidate = (
            cont
            & (width >= 3.0)
            & ~active_blocked
            & ~shield_candidate
            & (
                (code == PROVINCE_PLATFORM)
                | (code == PROVINCE_INTERIOR_BASIN)
                | platform_obj
                | basin_obj
                | (
                    old_basement
                    & (stability_floor >= 0.30)
                    & (thick_floor >= CONT_THICK + 900.0)
                    & ((quiet_age >= 500.0) | (width >= 4.0))
                )
            )
        )

        if shield_candidate.any():
            age_floor[shield_candidate] = np.maximum(age_floor[shield_candidate], shield_age_floor)
            stability_floor[shield_candidate] = np.maximum(stability_floor[shield_candidate], 0.84)
            thick_floor[shield_candidate] = np.maximum(
                thick_floor[shield_candidate], CONT_THICK + 5200.0)
        if platform_candidate.any():
            age_floor[platform_candidate] = np.maximum(
                age_floor[platform_candidate], platform_age_floor)
            stability_floor[platform_candidate] = np.maximum(
                stability_floor[platform_candidate], 0.76)
            thick_floor[platform_candidate] = np.maximum(
                thick_floor[platform_candidate], CONT_THICK + 2600.0)

        narrow_stale = active_blocked & (width < 2.0)
        age_floor[narrow_stale] = np.minimum(age_floor[narrow_stale], age[narrow_stale])
        stability_floor[narrow_stale] = np.minimum(
            stability_floor[narrow_stale], np.maximum(stability[narrow_stale], 0.0))
        thick_floor[narrow_stale] = np.minimum(
            thick_floor[narrow_stale], np.maximum(thick[narrow_stale], CONT_THICK))

        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1e-12)
        old_archived = cont & (age_floor >= 2500.0) & (stability_floor >= 0.70)
        shield_archived = cont & (age_floor >= 2500.0) & (stability_floor >= 0.80)
        platform_archived = (
            cont
            & ~shield_archived
            & (age_floor >= 2500.0)
            & (stability_floor >= 0.70)
        )
        world.set_g(
            "tectonics.last_p56_basement_archive_old_share_of_continental_crust",
            float(area[old_archived].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p56_basement_archive_shield_share_of_continental_crust",
            float(area[shield_archived].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p56_basement_archive_platform_share_of_continental_crust",
            float(area[platform_archived].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p56_basement_archive_active_blocked_share_of_continental_crust",
            float(area[active_blocked].sum() / cont_area),
        )
        return age_floor, stability_floor, thick_floor

    def _preserve_inherited_province_physical_cargo(
        self, world, grid, ctype, origin, age, thick, reworked, stability,
        province_code, t, *, basement_age_floor=None, basement_stability_floor=None,
        basement_thickness_floor=None, active_boundary_mask=None,
    ):
        """Carry inherited mature-province basement properties with crust parcels.

        P54 preserves the shield/platform label.  This step preserves the
        matching basement cargo: ancient age, thick stable lithosphere, and a
        quiet rework archive for broad interiors.  Active narrow margins can
        still be reset; broad shield/platform interiors should not be
        youngened merely because rasterized plate motion crossed an active
        boundary proxy nearby.
        """
        telemetry = (
            "tectonics.last_p55_shield_cargo_preserved_share_of_continental_crust",
            "tectonics.last_p55_platform_cargo_preserved_share_of_continental_crust",
            "tectonics.last_p55_total_cargo_preserved_share_of_continental_crust",
            "tectonics.last_p55_active_blocked_share_of_continental_crust",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        if world.g("tectonics.enable_inherited_province_physical_cargo", 0.0) <= 0.0:
            return origin, age, thick, reworked, stability

        cont = np.asarray(ctype, dtype=np.float64) == CONT
        if not cont.any():
            return origin, age, thick, reworked, stability

        tt = max(float(t), 0.0)
        mature_start = max(2500.0, 0.58 * float(world.spec.t_end_myr))
        if tt < mature_start:
            return origin, age, thick, reworked, stability

        code = np.asarray(province_code, dtype=np.float64)
        if code.shape != (grid.n,):
            return origin, age, thick, reworked, stability

        width = _graph_width_steps(grid, cont)
        active_core = np.zeros(grid.n, dtype=bool)
        if active_boundary_mask is not None:
            active = np.asarray(active_boundary_mask, dtype=bool)
            if active.shape == (grid.n,) and active.any():
                active_core = active & cont
        active_zone = self._dilate_mask(grid, active_core, allowed=cont, passes=1)
        quiet_age = np.where(
            np.asarray(reworked, dtype=np.float64) < 0.0,
            tt,
            np.maximum(tt - np.asarray(reworked, dtype=np.float64), 0.0),
        )
        basement_age = np.asarray(
            basement_age_floor if basement_age_floor is not None else np.zeros(grid.n),
            dtype=np.float64,
        )
        basement_stability = np.asarray(
            basement_stability_floor if basement_stability_floor is not None else np.zeros(grid.n),
            dtype=np.float64,
        )
        basement_thick = np.asarray(
            basement_thickness_floor if basement_thickness_floor is not None else np.zeros(grid.n),
            dtype=np.float64,
        )
        if basement_age.shape != (grid.n,):
            basement_age = np.zeros(grid.n, dtype=np.float64)
        if basement_stability.shape != (grid.n,):
            basement_stability = np.zeros(grid.n, dtype=np.float64)
        if basement_thick.shape != (grid.n,):
            basement_thick = np.zeros(grid.n, dtype=np.float64)

        shield_age_floor = min(tt, max(2500.0, 0.58 * tt))
        platform_age_floor = min(tt, max(2500.0, 0.52 * tt))
        shield_quiet_floor = min(1200.0, max(650.0, 0.22 * tt))
        platform_quiet_floor = min(1000.0, max(500.0, 0.18 * tt))
        basement_shield = (
            cont
            & (width >= 2.0)
            & (basement_age >= shield_age_floor)
            & (basement_stability >= 0.80)
            & (basement_thick >= CONT_THICK + 4200.0)
        )
        basement_platform = (
            cont
            & (width >= 3.0)
            & (basement_age >= platform_age_floor)
            & (basement_stability >= 0.70)
            & (basement_thick >= CONT_THICK + 1800.0)
            & ~basement_shield
        )

        memory_shield = cont & ((code == PROVINCE_SHIELD) | basement_shield)
        memory_platform = cont & (
            (code == PROVINCE_PLATFORM) | (code == PROVINCE_INTERIOR_BASIN)
            | basement_platform
        )
        shield_blocked = (
            memory_shield
            & active_core
            & (width < 3.0)
            & (stability < 0.70)
        )
        platform_blocked = (
            memory_platform
            & active_core
            & (width < 4.0)
            & (stability < 0.54)
        )

        shield_preserve = (
            memory_shield
            & (width >= 2.0)
            & ~shield_blocked
            & (
                (quiet_age >= 180.0)
                | (width >= 4.0)
                | (stability >= 0.58)
                | ~active_zone
            )
        )
        platform_preserve = (
            memory_platform
            & (width >= 3.0)
            & ~platform_blocked
            & (
                (quiet_age >= 240.0)
                | (width >= 4.0)
                | (stability >= 0.34)
                | ~active_zone
            )
        )

        if not (shield_preserve.any() or platform_preserve.any()):
            area = grid.cell_area
            cont_area = max(float(area[cont].sum()), 1e-12)
            blocked = shield_blocked | platform_blocked
            world.set_g(
                "tectonics.last_p55_active_blocked_share_of_continental_crust",
                float(area[blocked].sum() / cont_area),
            )
            return origin, age, thick, reworked, stability

        if shield_preserve.any():
            origin[shield_preserve] = ORIGIN_CRATON
            age[shield_preserve] = np.maximum.reduce([
                age[shield_preserve],
                np.full(np.count_nonzero(shield_preserve), shield_age_floor),
                basement_age[shield_preserve],
            ])
            thick[shield_preserve] = np.maximum.reduce([
                thick[shield_preserve],
                np.full(np.count_nonzero(shield_preserve), CONT_THICK + 5200.0),
                basement_thick[shield_preserve],
            ])
            stability[shield_preserve] = np.maximum.reduce([
                stability[shield_preserve],
                np.full(np.count_nonzero(shield_preserve), 0.84),
                basement_stability[shield_preserve],
            ])
            recent = shield_preserve & (
                (reworked >= 0.0) & ((tt - np.asarray(reworked, dtype=np.float64)) < shield_quiet_floor)
            )
            reworked[recent] = max(tt - shield_quiet_floor, -1.0)

        if platform_preserve.any():
            age[platform_preserve] = np.maximum.reduce([
                age[platform_preserve],
                np.full(np.count_nonzero(platform_preserve), platform_age_floor),
                basement_age[platform_preserve],
            ])
            thick[platform_preserve] = np.maximum.reduce([
                thick[platform_preserve],
                np.full(np.count_nonzero(platform_preserve), CONT_THICK + 2600.0),
                basement_thick[platform_preserve],
            ])
            stability[platform_preserve] = np.maximum.reduce([
                stability[platform_preserve],
                np.full(np.count_nonzero(platform_preserve), 0.76),
                basement_stability[platform_preserve],
            ])
            platform_origin = platform_preserve & (origin != ORIGIN_CRATON)
            origin[platform_origin] = ORIGIN_PRIMORDIAL
            recent = platform_preserve & (
                (reworked >= 0.0) & ((tt - np.asarray(reworked, dtype=np.float64)) < platform_quiet_floor)
            )
            reworked[recent] = max(tt - platform_quiet_floor, -1.0)

        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1e-12)
        blocked = shield_blocked | platform_blocked
        world.set_g(
            "tectonics.last_p55_shield_cargo_preserved_share_of_continental_crust",
            float(area[shield_preserve].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p55_platform_cargo_preserved_share_of_continental_crust",
            float(area[platform_preserve].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p55_total_cargo_preserved_share_of_continental_crust",
            float(area[shield_preserve | platform_preserve].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p55_active_blocked_share_of_continental_crust",
            float(area[blocked].sum() / cont_area),
        )
        return origin, age, thick, reworked, stability

    def _stable_cratonic_lithosphere_support(
        self, world, grid, ctype, age, thick, origin, reworked, stability, domain,
        province_code, cratonic_province_objects, basement_age, basement_stability,
        basement_thick, t, *, active_boundary_mask=None,
    ):
        """Non-mutating stable-craton continuity support field.

        P61 proves that old basement cargo can recover current stable crust, but
        applying that mutation inside the full-history sequence currently
        perturbs the archive/province trajectory.  P62 keeps the signals
        separated: this field records where current crust, basement cargo, and
        shield/platform/basin objects agree on stable cratonic lithosphere
        without changing `crust.*` arrays.
        """
        telemetry = (
            "tectonics.last_p62_current_stable_craton_fraction_of_continental_crust",
            "tectonics.last_p62_supported_cratonic_lithosphere_fraction_of_continental_crust",
            "tectonics.last_p62_object_supported_fraction_of_continental_crust",
            "tectonics.last_p62_basement_supported_fraction_of_continental_crust",
            "tectonics.last_p62_proxy_only_fraction_of_continental_crust",
            "tectonics.last_p62_active_blocked_share_of_continental_crust",
            "tectonics.last_p62_shield_support_share_of_continental_crust",
            "tectonics.last_p62_platform_support_share_of_continental_crust",
            "tectonics.last_p62_basin_support_share_of_continental_crust",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        support = np.zeros(grid.n, dtype=np.float64)
        cont = np.asarray(ctype, dtype=np.float64) == CONT
        if not cont.any():
            return support

        def optional_field(values) -> np.ndarray:
            arr = np.asarray(values, dtype=np.float64)
            if arr.shape != (grid.n,):
                return np.zeros(grid.n, dtype=np.float64)
            return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)

        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1.0e-12)
        width = _graph_width_steps(grid, cont)
        tt = max(float(t), 0.0)
        def optional_nullable_field(values) -> np.ndarray:
            if values is None:
                return np.zeros(grid.n, dtype=np.float64)
            return optional_field(values)

        domain_arr = optional_nullable_field(domain)
        province_arr = optional_nullable_field(province_code)
        basement_age_arr = optional_nullable_field(basement_age)
        basement_stability_arr = np.clip(optional_nullable_field(basement_stability), 0.0, 1.0)
        basement_thick_arr = optional_nullable_field(basement_thick)

        active_core = np.zeros(grid.n, dtype=bool)
        if active_boundary_mask is not None:
            active = np.asarray(active_boundary_mask, dtype=bool)
            if active.shape == (grid.n,) and active.any():
                active_core = active & cont
        active_zone = self._dilate_mask(grid, active_core, allowed=cont, passes=1)
        active_mobile_domain = np.isin(
            domain_arr,
            [DOMAIN_ACCRETED_TERRANE, DOMAIN_SUTURE, DOMAIN_LIP],
        )
        passive_margin_domain = domain_arr == DOMAIN_CONTINENTAL_MARGIN
        core_domain = ~np.isin(
            domain_arr,
            [DOMAIN_CONTINENTAL_MARGIN, DOMAIN_ACCRETED_TERRANE,
             DOMAIN_SUTURE, DOMAIN_LIP],
        )
        stable_passive_margin = (
            passive_margin_domain
            & (width >= 2.0)
            & (basement_stability_arr >= 0.78)
        )

        platform_floor = min(tt, max(2500.0, 0.52 * tt)) if tt > 0.0 else 2500.0
        old_stable_basement_any_domain = (
            cont
            & (width >= 2.0)
            & (basement_age_arr >= platform_floor)
            & (basement_stability_arr >= 0.70)
            & (basement_thick_arr >= CONT_THICK + 1200.0)
        )
        active_blocked = (
            old_stable_basement_any_domain
            & (
                active_mobile_domain
                | (
                    active_zone
                    & (width < 4.0)
                    & (np.asarray(stability, dtype=np.float64) < 0.66)
                )
            )
        )
        old_basement = (
            old_stable_basement_any_domain
            & (core_domain | stable_passive_margin)
            & ~active_blocked
        )

        current_stable = (
            cont
            & (width >= 2.0)
            & (np.asarray(age, dtype=np.float64) >= 2500.0)
            & (np.asarray(stability, dtype=np.float64) >= 0.75)
            & ~active_blocked
        )
        current_stable &= (
            (np.asarray(origin, dtype=np.float64) == ORIGIN_CRATON)
            | (domain_arr == DOMAIN_CRATON)
            | (width >= 3.0)
        )

        def object_cells(set_name: str) -> np.ndarray:
            mask = np.zeros(grid.n, dtype=bool)
            for obj in (cratonic_province_objects or {}).get(set_name, []) or []:
                cells = np.asarray(obj.get("cells", []), dtype=int)
                if cells.size == 0:
                    continue
                cells = cells[(0 <= cells) & (cells < grid.n)]
                mask[cells] = True
            return mask & cont & ~active_blocked

        shield_object = object_cells("tectonics.shields")
        platform_object = object_cells("tectonics.platforms")
        basin_object = object_cells("tectonics.interior_basins")
        memory = (
            cont
            & np.isin(
                province_arr,
                [PROVINCE_SHIELD, PROVINCE_PLATFORM, PROVINCE_INTERIOR_BASIN],
            )
            & (old_basement | current_stable)
            & ~active_blocked
        )

        support[current_stable] = np.maximum(support[current_stable], 0.82)
        support[old_basement] = np.maximum(support[old_basement], 0.72)
        support[memory] = np.maximum(support[memory], 0.74)
        support[basin_object] = np.maximum(support[basin_object], 0.76)
        support[platform_object] = np.maximum(support[platform_object], 0.84)
        support[shield_object] = np.maximum(support[shield_object], 0.96)
        support[active_blocked] = 0.0
        support = np.clip(support, 0.0, 1.0)

        supported = support >= 0.70
        object_supported = supported & (shield_object | platform_object | basin_object)
        basement_supported = supported & old_basement
        proxy_only = supported & ~current_stable
        world.set_g(
            "tectonics.last_p62_current_stable_craton_fraction_of_continental_crust",
            float(area[current_stable].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p62_supported_cratonic_lithosphere_fraction_of_continental_crust",
            float(area[supported].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p62_object_supported_fraction_of_continental_crust",
            float(area[object_supported].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p62_basement_supported_fraction_of_continental_crust",
            float(area[basement_supported].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p62_proxy_only_fraction_of_continental_crust",
            float(area[proxy_only].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p62_active_blocked_share_of_continental_crust",
            float(area[active_blocked].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p62_shield_support_share_of_continental_crust",
            float(area[supported & shield_object].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p62_platform_support_share_of_continental_crust",
            float(area[supported & platform_object].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p62_basin_support_share_of_continental_crust",
            float(area[supported & basin_object].sum() / cont_area),
        )
        return support

    def _project_supported_cratonic_current_state(
        self, world, grid, ctype, origin, age, thick, reworked, stability, domain,
        province_code, stable_cratonic_support, basement_age, basement_stability,
        basement_thick, t, dt, *, active_boundary_mask=None,
    ):
        """Area-limited current-state projection from P62 support.

        This is deliberately gated and runs after basement/province archive
        writeback.  It converts only mature, supported, inactive cratonic
        lithosphere into current stable crust.  Shield-like support can regain
        cratonic origin and thickness; covered-platform support receives old
        age/stability floors but is not relabelled as a shield.
        """
        telemetry = (
            "tectonics.last_p64_candidate_share_of_continental_crust",
            "tectonics.last_p64_projected_share_of_continental_crust",
            "tectonics.last_p64_shield_projected_share_of_continental_crust",
            "tectonics.last_p64_platform_projected_share_of_continental_crust",
            "tectonics.last_p64_active_blocked_share_of_continental_crust",
            "tectonics.last_p64_stable_craton_fraction_before",
            "tectonics.last_p64_stable_craton_fraction_after",
            "tectonics.last_p64_projected_cells",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        if world.g("tectonics.enable_p64_support_guided_craton_projection", 0.0) <= 0.0:
            return origin, age, thick, reworked, stability

        cont = np.asarray(ctype, dtype=np.float64) == CONT
        if not cont.any():
            return origin, age, thick, reworked, stability

        mature_time = max(3000.0, 0.86 * float(world.spec.t_end_myr))
        if float(t) < mature_time:
            return origin, age, thick, reworked, stability

        def optional_field(values) -> np.ndarray:
            if values is None:
                return np.zeros(grid.n, dtype=np.float64)
            arr = np.asarray(values, dtype=np.float64)
            if arr.shape != (grid.n,):
                return np.zeros(grid.n, dtype=np.float64)
            return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)

        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1.0e-12)
        width = _graph_width_steps(grid, cont)
        support = np.clip(optional_field(stable_cratonic_support), 0.0, 1.0)
        domain_arr = optional_field(domain)
        province_arr = optional_field(province_code)
        basement_age_arr = optional_field(basement_age)
        basement_stability_arr = np.clip(optional_field(basement_stability), 0.0, 1.0)
        basement_thick_arr = optional_field(basement_thick)

        active_core = np.zeros(grid.n, dtype=bool)
        if active_boundary_mask is not None:
            active = np.asarray(active_boundary_mask, dtype=bool)
            if active.shape == (grid.n,) and active.any():
                active_core = active & cont
        active_zone = self._dilate_mask(grid, active_core, allowed=cont, passes=1)
        active_mobile_domain = np.isin(
            domain_arr,
            [DOMAIN_ACCRETED_TERRANE, DOMAIN_SUTURE, DOMAIN_LIP],
        )
        active_blocked = (
            cont
            & (support >= 0.70)
            & (
                active_mobile_domain
                | (
                    active_zone
                    & (width < 4.0)
                    & (np.asarray(stability, dtype=np.float64) < 0.66)
                )
            )
        )
        current_stable = (
            cont
            & (np.asarray(age, dtype=np.float64) >= 2500.0)
            & (np.asarray(stability, dtype=np.float64) >= 0.75)
        )
        before_stable_share = float(area[current_stable].sum() / cont_area)
        supported_proxy = (
            cont
            & ~current_stable
            & ~active_blocked
            & (support >= 0.70)
            & (width >= 2.0)
            & (basement_age_arr >= 2500.0)
            & (basement_stability_arr >= 0.70)
        )
        shield_like = (
            supported_proxy
            & (
                (support >= 0.95)
                | (province_arr == PROVINCE_SHIELD)
                | (
                    (basement_stability_arr >= 0.82)
                    & (basement_thick_arr >= CONT_THICK + 4200.0)
                    & (width >= 3.0)
                )
            )
        )
        platform_like = supported_proxy & ~shield_like
        candidate = shield_like | platform_like
        if not candidate.any():
            world.set_g(
                "tectonics.last_p64_active_blocked_share_of_continental_crust",
                float(area[active_blocked].sum() / cont_area),
            )
            world.set_g(
                "tectonics.last_p64_stable_craton_fraction_before",
                before_stable_share,
            )
            world.set_g(
                "tectonics.last_p64_stable_craton_fraction_after",
                before_stable_share,
            )
            return origin, age, thick, reworked, stability

        deficit = max(0.0, 0.120 * cont_area - float(area[current_stable].sum()))
        target = min(float(area[candidate].sum()), max(deficit, 0.0))
        if target <= 0.0:
            target = min(float(area[candidate].sum()), 0.004 * cont_area)
        # Mature cratonic stabilization is limited to the already-supported
        # proxy pool.  At the terminal mature stage the cap must be large enough
        # to project the remaining support signal, but still far below a broad
        # continent rewrite.
        step_cap = cont_area * min(0.060, max(0.045, 0.035 * max(float(dt), 1.0) / 20.0))
        target = min(target, step_cap)
        ranked = sorted(np.where(candidate)[0].astype(int), key=lambda c: (
            0 if shield_like[c] else 1,
            -float(support[c]),
            -float(basement_age_arr[c]),
            -float(basement_stability_arr[c]),
            -float(width[c]),
            c,
        ))
        chosen: list[int] = []
        acc = 0.0
        for c in ranked:
            if acc >= target:
                break
            chosen.append(int(c))
            acc += float(area[c])
        if not chosen:
            return origin, age, thick, reworked, stability

        cells = np.asarray(chosen, dtype=np.int64)
        shield_cells = cells[shield_like[cells]]
        platform_cells = cells[~shield_like[cells]]
        if shield_cells.size:
            origin[shield_cells] = ORIGIN_CRATON
            age[shield_cells] = np.maximum(age[shield_cells], basement_age_arr[shield_cells])
            thick[shield_cells] = np.maximum.reduce([
                thick[shield_cells],
                basement_thick_arr[shield_cells],
                np.full(shield_cells.size, CONT_THICK + 4800.0),
            ])
            stability[shield_cells] = np.maximum.reduce([
                stability[shield_cells],
                basement_stability_arr[shield_cells],
                np.full(shield_cells.size, 0.82),
            ])
        if platform_cells.size:
            age[platform_cells] = np.maximum(age[platform_cells], basement_age_arr[platform_cells])
            stability[platform_cells] = np.maximum.reduce([
                stability[platform_cells],
                basement_stability_arr[platform_cells],
                np.full(platform_cells.size, 0.755),
            ])
            non_craton_platform = platform_cells[origin[platform_cells] != ORIGIN_CRATON]
            origin[non_craton_platform] = ORIGIN_PRIMORDIAL

        quiet_floor = min(1100.0, max(560.0, 0.20 * float(t)))
        recent = cells[(reworked[cells] >= 0.0) & ((float(t) - reworked[cells]) < quiet_floor)]
        if recent.size:
            reworked[recent] = max(float(t) - quiet_floor, -1.0)

        after_stable = (
            cont
            & (np.asarray(age, dtype=np.float64) >= 2500.0)
            & (np.asarray(stability, dtype=np.float64) >= 0.75)
        )
        world.set_g(
            "tectonics.last_p64_candidate_share_of_continental_crust",
            float(area[candidate].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p64_projected_share_of_continental_crust",
            float(area[cells].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p64_shield_projected_share_of_continental_crust",
            float(area[shield_cells].sum() / cont_area) if shield_cells.size else 0.0,
        )
        world.set_g(
            "tectonics.last_p64_platform_projected_share_of_continental_crust",
            float(area[platform_cells].sum() / cont_area) if platform_cells.size else 0.0,
        )
        world.set_g(
            "tectonics.last_p64_active_blocked_share_of_continental_crust",
            float(area[active_blocked].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p64_stable_craton_fraction_before",
            before_stable_share,
        )
        world.set_g(
            "tectonics.last_p64_stable_craton_fraction_after",
            float(area[after_stable].sum() / cont_area),
        )
        world.set_g("tectonics.last_p64_projected_cells", float(cells.size))
        return origin, age, thick, reworked, stability

    def _refresh_supported_ancient_age_continuity(
        self, world, grid, ctype, age, reworked, stability, domain, province_code,
        stable_cratonic_support, basement_age, basement_stability, basement_thick, t,
        *, active_boundary_mask=None,
    ):
        """Restore narrow ancient-age continuity from supported old basement.

        P64 projects current stable cratonic state only where the support signal
        is strong enough.  P65 keeps that projection bounded but lets a small
        set of inactive, support-backed old-basement cells retain ancient age
        evidence.  It deliberately does not thicken crust, change origin, or
        raise stability; those stronger current-state mutations remain P64's
        responsibility.
        """
        telemetry = (
            "tectonics.last_p65_ancient_age_candidate_share_of_continental_crust",
            "tectonics.last_p65_ancient_age_refreshed_share_of_continental_crust",
            "tectonics.last_p65_active_blocked_share_of_continental_crust",
            "tectonics.last_p65_ancient_fraction_before",
            "tectonics.last_p65_ancient_fraction_after",
            "tectonics.last_p65_age_refreshed_cells",
            "tectonics.last_p65_province_code_changed_share_of_continental_crust",
            "tectonics.last_p65_refreshed_p48_shield_share_of_continental_crust",
            "tectonics.last_p65_refreshed_p48_platform_share_of_continental_crust",
            "tectonics.last_p65_refreshed_p48_mature_share_of_continental_crust",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        if world.g("tectonics.enable_p65_post_projection_cratonic_refresh", 0.0) <= 0.0:
            return age, reworked

        cont = np.asarray(ctype, dtype=np.float64) == CONT
        if not cont.any():
            return age, reworked

        mature_time = max(3000.0, 0.86 * float(world.spec.t_end_myr))
        if float(t) < mature_time:
            return age, reworked

        def optional_field(values) -> np.ndarray:
            if values is None:
                return np.zeros(grid.n, dtype=np.float64)
            arr = np.asarray(values, dtype=np.float64)
            if arr.shape != (grid.n,):
                return np.zeros(grid.n, dtype=np.float64)
            return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)

        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1.0e-12)
        width = _graph_width_steps(grid, cont)
        support = np.clip(optional_field(stable_cratonic_support), 0.0, 1.0)
        domain_arr = optional_field(domain)
        province_arr = optional_field(province_code)
        basement_age_arr = optional_field(basement_age)
        basement_stability_arr = np.clip(optional_field(basement_stability), 0.0, 1.0)
        basement_thick_arr = optional_field(basement_thick)
        age_arr = np.asarray(age, dtype=np.float64)
        stability_arr = np.asarray(stability, dtype=np.float64)

        active_core = np.zeros(grid.n, dtype=bool)
        if active_boundary_mask is not None:
            active = np.asarray(active_boundary_mask, dtype=bool)
            if active.shape == (grid.n,) and active.any():
                active_core = active & cont
        active_zone = self._dilate_mask(grid, active_core, allowed=cont, passes=1)
        active_mobile_domain = np.isin(
            domain_arr,
            [DOMAIN_ACCRETED_TERRANE, DOMAIN_SUTURE, DOMAIN_LIP],
        )
        core_domain = ~np.isin(
            domain_arr,
            [DOMAIN_CONTINENTAL_MARGIN, DOMAIN_ACCRETED_TERRANE,
             DOMAIN_SUTURE, DOMAIN_LIP],
        )
        stable_passive_margin = (
            (domain_arr == DOMAIN_CONTINENTAL_MARGIN)
            & (width >= 2.0)
            & (basement_stability_arr >= 0.80)
        )
        memory_province = np.isin(
            province_arr,
            [PROVINCE_SHIELD, PROVINCE_PLATFORM, PROVINCE_INTERIOR_BASIN],
        )
        old_basement_archive = (
            cont
            & (basement_age_arr >= 2500.0)
            & (basement_stability_arr >= 0.70)
        )
        archive_context = (
            old_basement_archive
            & (width >= 2.0)
            & (
                (support >= 0.70)
                | memory_province
                | (
                    (basement_stability_arr >= 0.72)
                    & (
                        (width >= 3.0)
                        | (basement_thick_arr >= CONT_THICK + 600.0)
                    )
                )
            )
        )
        old_supported = (
            archive_context
            & (core_domain | stable_passive_margin)
        )
        active_blocked = (
            old_supported
            & (
                active_mobile_domain
                | (
                    active_zone
                    & (width < 4.0)
                    & (stability_arr < 0.66)
                )
            )
        )
        candidate = (
            old_supported
            & ~active_blocked
            & (age_arr < 2500.0)
            & (
                (stability_arr < 0.75)
                | np.isin(province_arr, [PROVINCE_PLATFORM, PROVINCE_INTERIOR_BASIN])
            )
        )
        before_ancient = cont & (age_arr >= 2500.0)
        before_share = float(area[before_ancient].sum() / cont_area)
        world.set_g(
            "tectonics.last_p65_ancient_age_candidate_share_of_continental_crust",
            float(area[candidate].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p65_active_blocked_share_of_continental_crust",
            float(area[active_blocked].sum() / cont_area),
        )
        world.set_g("tectonics.last_p65_ancient_fraction_before", before_share)
        if not candidate.any():
            world.set_g("tectonics.last_p65_ancient_fraction_after", before_share)
            return age, reworked

        deficit = max(0.0, 0.202 * cont_area - float(area[before_ancient].sum()))
        target = min(float(area[candidate].sum()), deficit + 0.002 * cont_area)
        target = min(target, 0.010 * cont_area)
        if target <= 0.0:
            world.set_g("tectonics.last_p65_ancient_fraction_after", before_share)
            return age, reworked

        ranked = sorted(np.where(candidate)[0].astype(int), key=lambda c: (
            0 if stability_arr[c] < 0.75 else 1,
            -float(support[c]),
            -float(basement_age_arr[c]),
            -float(basement_stability_arr[c]),
            -float(width[c]),
            c,
        ))
        chosen: list[int] = []
        acc = 0.0
        for c in ranked:
            if acc >= target:
                break
            chosen.append(int(c))
            acc += float(area[c])
        if not chosen:
            world.set_g("tectonics.last_p65_ancient_fraction_after", before_share)
            return age, reworked

        cells = np.asarray(chosen, dtype=np.int64)
        age[cells] = np.maximum(age[cells], np.minimum(float(t), basement_age_arr[cells]))
        quiet_floor = min(900.0, max(520.0, 0.16 * float(t)))
        recent = cells[(reworked[cells] >= 0.0) & ((float(t) - reworked[cells]) < quiet_floor)]
        if recent.size:
            reworked[recent] = max(float(t) - quiet_floor, -1.0)

        after_ancient = cont & (np.asarray(age, dtype=np.float64) >= 2500.0)
        world.set_g(
            "tectonics.last_p65_ancient_age_refreshed_share_of_continental_crust",
            float(area[cells].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p65_ancient_fraction_after",
            float(area[after_ancient].sum() / cont_area),
        )
        world.set_g("tectonics.last_p65_age_refreshed_cells", float(cells.size))
        return age, reworked

    def _augment_post_projection_supported_platform_objects(
        self, world, grid, ctype, age, thick, reworked, stability, domain,
        continent_id, province_code, cratonic_province_objects,
        stable_cratonic_support, basement_age, basement_stability, basement_thick, t,
        *, active_boundary_mask=None,
    ):
        """Add small post-projection mature platform objects from archive evidence."""
        for key in (
            "tectonics.last_p65_added_platform_object_share_of_continental_crust",
            "tectonics.last_p65_added_platform_object_count",
        ):
            world.set_g(key, 0.0)

        cont = np.asarray(ctype, dtype=np.float64) == CONT
        if not cont.any():
            return cratonic_province_objects, province_code

        def optional_field(values) -> np.ndarray:
            if values is None:
                return np.zeros(grid.n, dtype=np.float64)
            arr = np.asarray(values, dtype=np.float64)
            if arr.shape != (grid.n,):
                return np.zeros(grid.n, dtype=np.float64)
            return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)

        area = grid.cell_area
        total_area = max(float(area.sum()), 1.0e-12)
        cont_area = max(float(area[cont].sum()), 1.0e-12)
        width = _graph_width_steps(grid, cont)
        support = np.clip(optional_field(stable_cratonic_support), 0.0, 1.0)
        province_arr = optional_field(province_code)
        basement_age_arr = optional_field(basement_age)
        basement_stability_arr = np.clip(optional_field(basement_stability), 0.0, 1.0)
        basement_thick_arr = optional_field(basement_thick)
        domain_arr = optional_field(domain)

        active_core = np.zeros(grid.n, dtype=bool)
        if active_boundary_mask is not None:
            active = np.asarray(active_boundary_mask, dtype=bool)
            if active.shape == (grid.n,) and active.any():
                active_core = active & cont
        active_zone = self._dilate_mask(grid, active_core, allowed=cont, passes=1)
        active_mobile_domain = np.isin(
            domain_arr,
            [DOMAIN_ACCRETED_TERRANE, DOMAIN_SUTURE, DOMAIN_LIP],
        )
        core_or_passive = ~np.isin(
            domain_arr,
            [DOMAIN_ACCRETED_TERRANE, DOMAIN_SUTURE, DOMAIN_LIP],
        )
        active_blocked = (
            cont
            & (
                active_mobile_domain
                | (
                    active_zone
                    & (width < 4.0)
                    & (np.asarray(stability, dtype=np.float64) < 0.66)
                )
            )
        )

        existing = np.zeros(grid.n, dtype=bool)
        for set_name in (
            "tectonics.shields",
            "tectonics.platforms",
            "tectonics.interior_basins",
        ):
            for obj in (cratonic_province_objects or {}).get(set_name, []) or []:
                cells = np.asarray(obj.get("cells", []), dtype=int)
                cells = cells[(0 <= cells) & (cells < grid.n)]
                existing[cells] = True
        existing &= cont
        existing_share = float(area[existing].sum() / cont_area)
        target_area = max(0.0, 0.162 * cont_area - float(area[existing].sum()))
        target_area = min(target_area, 0.010 * cont_area)
        if target_area <= 0.0:
            return cratonic_province_objects, province_code

        quiet_age = np.where(
            np.asarray(reworked, dtype=np.float64) < 0.0,
            max(float(t), 0.0),
            np.maximum(float(t) - np.asarray(reworked, dtype=np.float64), 0.0),
        )
        quiet_floor = min(900.0, max(420.0, 0.15 * float(t)))
        old_platform_evidence = (
            cont
            & ~existing
            & ~active_blocked
            & core_or_passive
            & (width >= 2.0)
            & (basement_age_arr >= 2500.0)
            & (basement_stability_arr >= 0.70)
            & (
                (support >= 0.70)
                | np.isin(province_arr, [PROVINCE_PLATFORM, PROVINCE_INTERIOR_BASIN])
                | (
                    (np.asarray(age, dtype=np.float64) >= 2500.0)
                    & (np.asarray(stability, dtype=np.float64) >= 0.45)
                )
            )
            & (
                (quiet_age >= quiet_floor)
                | (support >= 0.84)
                | np.isin(province_arr, [PROVINCE_PLATFORM, PROVINCE_INTERIOR_BASIN])
            )
        )
        if not old_platform_evidence.any():
            return cratonic_province_objects, province_code

        comps = self._connected_components(grid, old_platform_evidence)
        comps = [comp for comp in comps if comp.size >= 2]
        comps.sort(key=lambda comp: (
            -float(np.average(support[comp], weights=area[comp])),
            -float(np.average(basement_age_arr[comp], weights=area[comp])),
            -float(area[comp].sum()),
            int(comp.min()),
        ))
        selected = np.zeros(grid.n, dtype=bool)
        acc = 0.0
        for comp in comps:
            comp_area = float(area[comp].sum())
            if acc > 0.0 and acc + comp_area > target_area + 0.004 * cont_area:
                continue
            selected[comp] = True
            acc += comp_area
            if acc >= target_area:
                break
        if not selected.any():
            return cratonic_province_objects, province_code

        object_age = np.maximum(np.asarray(age, dtype=np.float64), basement_age_arr)
        object_thick = np.maximum(np.asarray(thick, dtype=np.float64), basement_thick_arr)
        object_stability = np.maximum(
            np.asarray(stability, dtype=np.float64), basement_stability_arr)
        new_objects = self._cratonic_province_component_objects(
            grid,
            "platform",
            selected,
            "support_guided_post_projection_mature_platform",
            object_age,
            object_thick,
            object_stability,
            width,
            np.asarray(continent_id, dtype=np.float64),
            t,
            total_area,
            max_objects=8,
        )
        if not new_objects:
            return cratonic_province_objects, province_code

        out = {
            "tectonics.shields": list((cratonic_province_objects or {}).get("tectonics.shields", [])),
            "tectonics.platforms": list((cratonic_province_objects or {}).get("tectonics.platforms", [])),
            "tectonics.interior_basins": list((cratonic_province_objects or {}).get("tectonics.interior_basins", [])),
        }
        out["tectonics.platforms"].extend(new_objects)
        province_code[selected] = np.where(
            province_code[selected] == PROVINCE_SHIELD,
            province_code[selected],
            PROVINCE_PLATFORM,
        )

        added_share = float(area[selected].sum() / cont_area)
        platform_cells = np.zeros(grid.n, dtype=bool)
        for obj in out["tectonics.platforms"]:
            cells = np.asarray(obj.get("cells", []), dtype=int)
            cells = cells[(0 <= cells) & (cells < grid.n)]
            platform_cells[cells] = True
        world.set_g(
            "tectonics.last_p48_platform_object_count",
            float(len(out["tectonics.platforms"])),
        )
        world.set_g(
            "tectonics.last_p48_platform_share_of_continental_crust",
            float(area[platform_cells & cont].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p65_added_platform_object_share_of_continental_crust",
            added_share,
        )
        world.set_g(
            "tectonics.last_p65_added_platform_object_count",
            float(len(new_objects)),
        )
        world.set_g(
            "tectonics.last_p65_existing_mature_share_before_object_augmentation",
            existing_share,
        )
        return out, province_code

    def _rebalance_physical_ensemble_cratonic_state(
        self, world, grid, ctype, origin, age, thick, reworked, stability, domain,
        continent_id, province_code, cratonic_province_objects,
        stable_cratonic_support, basement_age, basement_stability, basement_thick,
        t, *, active_boundary_mask=None,
    ):
        """Balance cratonic continuity against over-preserved ensemble members.

        P67 showed that a reference Earth-like member can pass continuity floors
        while dry/cool/fewer-plate members become nearly all ancient basement,
        and wet/high-plate members under-express platforms.  P68 keeps the same
        evidence chain but adds a gated balance pass: hard craton cores remain
        protected; weak/mobile, marginal, or non-core mature cells can be
        represented as reworked provinces or covered platforms.
        """
        telemetry = (
            "tectonics.last_p68_old_archive_fraction_before",
            "tectonics.last_p68_old_archive_fraction_after",
            "tectonics.last_p68_stable_fraction_before",
            "tectonics.last_p68_stable_fraction_after",
            "tectonics.last_p68_ancient_fraction_before",
            "tectonics.last_p68_ancient_fraction_after",
            "tectonics.last_p68_mature_object_fraction_before",
            "tectonics.last_p68_target_stable_fraction",
            "tectonics.last_p68_reworked_balance_share_of_continental_crust",
            "tectonics.last_p68_platform_reclassified_share_of_continental_crust",
            "tectonics.last_p68_protected_core_share_of_continental_crust",
            "tectonics.last_p68_active_blocked_share_of_continental_crust",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        cont = np.asarray(ctype, dtype=np.float64) == CONT
        if not cont.any():
            return (
                origin, age, thick, reworked, stability, province_code,
                basement_age, basement_stability, basement_thick,
            )

        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1.0e-12)
        width = _graph_width_steps(grid, cont)
        age_arr = np.asarray(age, dtype=np.float64)
        thick_arr = np.asarray(thick, dtype=np.float64)
        stability_arr = np.asarray(stability, dtype=np.float64)
        origin_arr = np.asarray(origin, dtype=np.float64)
        domain_arr = np.asarray(domain, dtype=np.float64)
        province_arr = np.asarray(province_code, dtype=np.float64)
        basement_age_arr = np.asarray(basement_age, dtype=np.float64)
        basement_stability_arr = np.clip(np.asarray(basement_stability, dtype=np.float64), 0.0, 1.0)
        basement_thick_arr = np.asarray(basement_thick, dtype=np.float64)
        support = np.clip(np.asarray(stable_cratonic_support, dtype=np.float64), 0.0, 1.0)
        if support.shape != (grid.n,):
            support = np.zeros(grid.n, dtype=np.float64)

        mature_object = np.zeros(grid.n, dtype=bool)
        shield_object = np.zeros(grid.n, dtype=bool)
        platform_object = np.zeros(grid.n, dtype=bool)
        for set_name, target in (
            ("tectonics.shields", shield_object),
            ("tectonics.platforms", platform_object),
            ("tectonics.interior_basins", mature_object),
        ):
            for obj in (cratonic_province_objects or {}).get(set_name, []) or []:
                cells = np.asarray(obj.get("cells", []), dtype=int)
                cells = cells[(0 <= cells) & (cells < grid.n)]
                target[cells] = True
                mature_object[cells] = True
        mature_object |= shield_object | platform_object

        active_core = np.zeros(grid.n, dtype=bool)
        if active_boundary_mask is not None:
            active = np.asarray(active_boundary_mask, dtype=bool)
            if active.shape == (grid.n,) and active.any():
                active_core = active & cont
        active_zone = self._dilate_mask(grid, active_core, allowed=cont, passes=1)
        active_mobile_domain = np.isin(
            domain_arr,
            [DOMAIN_CONTINENTAL_MARGIN, DOMAIN_ACCRETED_TERRANE,
             DOMAIN_SUTURE, DOMAIN_LIP],
        )
        active_blocked = active_mobile_domain | (active_zone & (width < 3.0))

        old_archive = cont & (basement_age_arr >= 2500.0) & (basement_stability_arr >= 0.70)
        ancient = cont & (age_arr >= 2500.0)
        current_stable = cont & (age_arr >= 2500.0) & (stability_arr >= 0.75)
        hard_core = (
            cont
            & (width >= 3.0)
            & current_stable
            & (
                (origin_arr == ORIGIN_CRATON)
                | (domain_arr == DOMAIN_CRATON)
            )
            & (basement_stability_arr >= 0.82)
            & (basement_thick_arr >= CONT_THICK + 4200.0)
        )
        protected_core = hard_core | (
            cont
            & (support >= 0.98)
            & (width >= 4.0)
            & (origin_arr == ORIGIN_CRATON)
            & (basement_stability_arr >= 0.86)
            & (basement_thick_arr >= CONT_THICK + 5000.0)
        )

        old_share_before = float(area[old_archive].sum() / cont_area)
        stable_share_before = float(area[current_stable].sum() / cont_area)
        ancient_share_before = float(area[ancient].sum() / cont_area)
        mature_share_before = float(area[mature_object & cont].sum() / cont_area)
        world.set_g("tectonics.last_p68_old_archive_fraction_before", old_share_before)
        world.set_g("tectonics.last_p68_stable_fraction_before", stable_share_before)
        world.set_g("tectonics.last_p68_ancient_fraction_before", ancient_share_before)
        world.set_g(
            "tectonics.last_p68_mature_object_fraction_before",
            mature_share_before,
        )
        world.set_g(
            "tectonics.last_p68_protected_core_share_of_continental_crust",
            float(area[protected_core & cont].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p68_active_blocked_share_of_continental_crust",
            float(area[active_blocked & cont].sum() / cont_area),
        )

        member_target = float(np.clip(getattr(world.spec, "target_land_fraction", 0.29), 0.0, 1.0))
        water = float(getattr(world.spec.composition, "water_inventory_earth", 1.0))
        plate_count = max(int(getattr(world.spec, "n_plates", 12)), 1)
        crust_strength = float(getattr(world.spec, "crust_strength", 1.0))
        dry_cool_bias = np.clip(
            0.45 * (member_target - 0.29) / 0.08
            + 0.35 * (1.0 - water) / 0.30
            + 0.20 * (crust_strength - 1.0) / 0.18,
            0.0,
            1.0,
        )
        wet_mobile_bias = np.clip(
            0.45 * (water - 1.0) / 0.30
            + 0.35 * (plate_count - 12.0) / 4.0
            + 0.20 * (1.0 - crust_strength) / 0.14,
            0.0,
            1.0,
        )

        target_old_share = 0.70 - 0.02 * dry_cool_bias
        target_ancient_share = 0.66 - 0.02 * dry_cool_bias
        target_mature_share = 0.60 - 0.01 * dry_cool_bias
        protected_stable_share = float(area[current_stable & protected_core].sum() / cont_area)
        target_stable_share = min(
            0.445,
            max(0.420 + 0.020 * dry_cool_bias, protected_stable_share + 0.012),
        )
        floor_share = max(
            0.22,
            protected_stable_share + 0.12,
        )
        target_old_share = max(target_old_share, floor_share)
        target_ancient_share = max(target_ancient_share, floor_share)
        world.set_g("tectonics.last_p68_target_stable_fraction", target_stable_share)

        rework_candidates = (
            cont
            & ~protected_core
            & (old_archive | ancient | mature_object)
            & (
                active_mobile_domain
                |
                (width <= 4.0)
                | (support < 0.90)
                | (province_arr != PROVINCE_SHIELD)
                | (origin_arr != ORIGIN_CRATON)
            )
        )
        stable_excess_area = max(
            0.0,
            float(area[current_stable].sum()) - target_stable_share * cont_area,
        )
        excess_area = max(0.0, float(area[old_archive].sum()) - target_old_share * cont_area)
        excess_area = max(
            excess_area,
            stable_excess_area,
        )
        excess_area = max(
            excess_area,
            max(0.0, float(area[ancient].sum()) - target_ancient_share * cont_area),
        )
        excess_area = max(
            excess_area,
            max(0.0, float(area[mature_object & cont].sum()) - target_mature_share * cont_area),
        )
        reworked_balance = np.zeros(grid.n, dtype=bool)
        if excess_area > 0.0 and rework_candidates.any():
            ranked = sorted(np.where(rework_candidates)[0].astype(int), key=lambda c: (
                0 if (stable_excess_area > 0.0 and current_stable[c]) else 1,
                0 if mature_object[c] else 1,
                0 if province_arr[c] != PROVINCE_SHIELD else 1,
                0 if origin_arr[c] != ORIGIN_CRATON else 1,
                0 if support[c] < 0.84 else 1,
                float(width[c]),
                float(stability_arr[c]),
                float(basement_stability_arr[c]),
                c,
            ))
            acc = 0.0
            chosen: list[int] = []
            max_rework = min(excess_area + 0.015 * cont_area, 0.48 * cont_area)
            for c in ranked:
                if acc >= max_rework:
                    break
                chosen.append(int(c))
                acc += float(area[c])
            if chosen:
                cells = np.asarray(chosen, dtype=np.int64)
                reworked_balance[cells] = True
                age[cells] = np.minimum(
                    age[cells],
                    1750.0 + 140.0 * np.clip(width[cells], 0.0, 4.0),
                )
                stability[cells] = np.minimum(stability[cells], 0.66)
                thick[cells] = np.minimum(thick[cells], CONT_THICK + 1800.0)
                origin[cells] = np.where(
                    origin[cells] == ORIGIN_CRATON,
                    ORIGIN_PRIMORDIAL,
                    origin[cells],
                )
                basement_age[cells] = np.minimum(basement_age[cells], 2250.0)
                basement_stability[cells] = np.minimum(basement_stability[cells], 0.66)
                basement_thick[cells] = np.minimum(basement_thick[cells], CONT_THICK + 1600.0)
                province_code[cells] = PROVINCE_NONE
                reworked[cells] = max(float(t) - 240.0, -1.0)

        platform_share_before = float(area[platform_object & cont].sum() / cont_area)
        platform_target = 0.055 + 0.020 * wet_mobile_bias
        platform_deficit = max(0.0, platform_target * cont_area - float(area[platform_object & cont].sum()))
        platform_candidates = (
            cont
            & ~reworked_balance
            & ~active_blocked
            & old_archive
            & (mature_object | shield_object | (support >= 0.72))
            & (width >= 2.0)
            & (
                (wet_mobile_bias > 0.15)
                | (platform_share_before < 0.035)
            )
        )
        platform_reclassified = np.zeros(grid.n, dtype=bool)
        if platform_deficit > 0.0 and platform_candidates.any():
            ranked = sorted(np.where(platform_candidates)[0].astype(int), key=lambda c: (
                0 if shield_object[c] else 1,
                0 if support[c] < 0.92 else 1,
                -float(width[c]),
                float(basement_thick_arr[c]),
                c,
            ))
            acc = 0.0
            chosen: list[int] = []
            max_platform = min(platform_deficit + 0.006 * cont_area, 0.080 * cont_area)
            for c in ranked:
                if acc >= max_platform:
                    break
                chosen.append(int(c))
                acc += float(area[c])
            if chosen:
                cells = np.asarray(chosen, dtype=np.int64)
                platform_reclassified[cells] = True
                origin[cells] = np.where(
                    origin[cells] == ORIGIN_CRATON,
                    ORIGIN_PRIMORDIAL,
                    origin[cells],
                )
                age[cells] = np.maximum(age[cells], 2500.0)
                stability[cells] = np.clip(stability[cells], 0.72, 0.78)
                thick[cells] = np.clip(thick[cells], CONT_THICK + 1800.0, CONT_THICK + 3200.0)
                basement_age[cells] = np.maximum(basement_age[cells], 2500.0)
                basement_stability[cells] = np.clip(basement_stability[cells], 0.72, 0.78)
                basement_thick[cells] = np.clip(
                    basement_thick[cells], CONT_THICK + 2000.0, CONT_THICK + 3200.0)
                province_code[cells] = PROVINCE_PLATFORM
                quiet_floor = min(900.0, max(560.0, 0.16 * float(t)))
                recent = cells[(reworked[cells] >= 0.0) & ((float(t) - reworked[cells]) < quiet_floor)]
                if recent.size:
                    reworked[recent] = max(float(t) - quiet_floor, -1.0)

        final_old_archive = (
            cont
            & (np.asarray(basement_age, dtype=np.float64) >= 2500.0)
            & (np.asarray(basement_stability, dtype=np.float64) >= 0.70)
        )
        final_stable = (
            cont
            & (np.asarray(age, dtype=np.float64) >= 2500.0)
            & (np.asarray(stability, dtype=np.float64) >= 0.75)
        )
        final_ancient = cont & (np.asarray(age, dtype=np.float64) >= 2500.0)
        world.set_g(
            "tectonics.last_p68_old_archive_fraction_after",
            float(area[final_old_archive].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p68_stable_fraction_after",
            float(area[final_stable].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p68_ancient_fraction_after",
            float(area[final_ancient].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p68_reworked_balance_share_of_continental_crust",
            float(area[reworked_balance].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p68_platform_reclassified_share_of_continental_crust",
            float(area[platform_reclassified].sum() / cont_area),
        )
        return (
            origin, age, thick, reworked, stability, province_code,
            basement_age, basement_stability, basement_thick,
        )

    def _stabilize_basement_cratonic_current_state(
        self, world, grid, ctype, origin, age, thick, reworked, stability, domain,
        basement_age, basement_stability, basement_thick, t,
        *, active_boundary_mask=None,
    ):
        """Recover current crust state from persistent old basement archive.

        P60 can name covered platforms from basement cargo without changing
        current crust.  This pass is narrower: if a cell carries old stable
        continental basement and is not in an active mobile belt, the current
        crust should retain old/stable lithosphere state.  Exposed strong cores
        regain cratonic origin; covered/passive-margin platforms regain old
        age and stability but remain platforms rather than shields.
        """
        telemetry = (
            "tectonics.last_p61_old_basement_current_state_candidate_share_of_continental_crust",
            "tectonics.last_p61_shield_current_state_stabilized_share_of_continental_crust",
            "tectonics.last_p61_platform_current_state_stabilized_share_of_continental_crust",
            "tectonics.last_p61_active_blocked_share_of_continental_crust",
            "tectonics.last_p61_stable_craton_fraction_before",
            "tectonics.last_p61_stable_craton_fraction_after",
            "tectonics.last_p61_current_state_stabilized_cells",
        )
        for key in telemetry:
            world.set_g(key, 0.0)

        if world.g("tectonics.enable_p61_basement_current_state_continuity", 0.0) <= 0.0:
            return origin, age, thick, reworked, stability

        cont = np.asarray(ctype, dtype=np.float64) == CONT
        if not cont.any():
            return origin, age, thick, reworked, stability

        age_floor = np.asarray(basement_age, dtype=np.float64)
        stability_floor = np.asarray(basement_stability, dtype=np.float64)
        thick_floor = np.asarray(basement_thick, dtype=np.float64)
        if (
            age_floor.shape != (grid.n,)
            or stability_floor.shape != (grid.n,)
            or thick_floor.shape != (grid.n,)
        ):
            return origin, age, thick, reworked, stability

        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1.0e-12)
        before_stable = cont & (age >= 2500.0) & (stability >= 0.75)
        width = _graph_width_steps(grid, cont)
        tt = max(float(t), 0.0)
        platform_age_floor = min(tt, max(2500.0, 0.52 * tt)) if tt > 0.0 else 2500.0
        shield_age_floor = min(tt, max(2500.0, 0.58 * tt)) if tt > 0.0 else 2500.0
        active_core = np.zeros(grid.n, dtype=bool)
        if active_boundary_mask is not None:
            active = np.asarray(active_boundary_mask, dtype=bool)
            if active.shape == (grid.n,) and active.any():
                active_core = active & cont
        active_zone = self._dilate_mask(grid, active_core, allowed=cont, passes=1)
        active_mobile_domain = np.isin(domain, [DOMAIN_ACCRETED_TERRANE, DOMAIN_SUTURE, DOMAIN_LIP])
        passive_margin_domain = domain == DOMAIN_CONTINENTAL_MARGIN
        interior_or_craton_domain = np.isin(domain, [DOMAIN_CONTINENTAL_INTERIOR, DOMAIN_CRATON])
        old_stable_basement = (
            cont
            & (age_floor >= platform_age_floor)
            & (stability_floor >= 0.70)
            & (thick_floor >= CONT_THICK + 1200.0)
            & (width >= 2.0)
        )
        active_blocked = (
            old_stable_basement
            & (
                active_mobile_domain
                | (
                    active_zone
                    & (width < 4.0)
                    & (stability < 0.66)
                )
            )
        )
        stable_passive_margin = (
            old_stable_basement
            & passive_margin_domain
            & (stability_floor >= 0.78)
            & ~active_blocked
        )
        interior_old_platform = (
            old_stable_basement
            & interior_or_craton_domain
            & ~active_blocked
        )
        shield_candidate = (
            interior_old_platform
            & (age_floor >= shield_age_floor)
            & (stability_floor >= 0.80)
            & (thick_floor >= CONT_THICK + 4200.0)
            & (width >= 2.0)
        )
        platform_candidate = (
            (stable_passive_margin | interior_old_platform)
            & ~shield_candidate
            & ~active_blocked
        )
        current_unstable = (age < 2500.0) | (stability < 0.75) | (thick < CONT_THICK + 1800.0)
        shield_stabilize = shield_candidate & current_unstable
        platform_stabilize = platform_candidate & current_unstable
        stabilized = shield_stabilize | platform_stabilize
        if not stabilized.any():
            after_stable = before_stable
            world.set_g(
                "tectonics.last_p61_old_basement_current_state_candidate_share_of_continental_crust",
                float(area[old_stable_basement & ~active_blocked].sum() / cont_area),
            )
            world.set_g(
                "tectonics.last_p61_active_blocked_share_of_continental_crust",
                float(area[active_blocked].sum() / cont_area),
            )
            world.set_g(
                "tectonics.last_p61_stable_craton_fraction_before",
                float(area[before_stable].sum() / cont_area),
            )
            world.set_g(
                "tectonics.last_p61_stable_craton_fraction_after",
                float(area[after_stable].sum() / cont_area),
            )
            return origin, age, thick, reworked, stability

        if shield_stabilize.any():
            origin[shield_stabilize] = ORIGIN_CRATON
            age[shield_stabilize] = np.maximum(age[shield_stabilize], age_floor[shield_stabilize])
            thick[shield_stabilize] = np.maximum.reduce([
                thick[shield_stabilize],
                thick_floor[shield_stabilize],
                np.full(np.count_nonzero(shield_stabilize), CONT_THICK + 5200.0),
            ])
            stability[shield_stabilize] = np.maximum.reduce([
                stability[shield_stabilize],
                stability_floor[shield_stabilize],
                np.full(np.count_nonzero(shield_stabilize), 0.84),
            ])
        if platform_stabilize.any():
            age[platform_stabilize] = np.maximum(
                age[platform_stabilize], age_floor[platform_stabilize])
            stability[platform_stabilize] = np.maximum.reduce([
                stability[platform_stabilize],
                stability_floor[platform_stabilize],
                np.full(np.count_nonzero(platform_stabilize), 0.76),
            ])
            platform_origin = platform_stabilize & (origin != ORIGIN_CRATON)
            origin[platform_origin] = ORIGIN_PRIMORDIAL

        quiet_floor = min(1000.0, max(520.0, 0.18 * tt))
        recent = stabilized & (
            (reworked >= 0.0) & ((tt - np.asarray(reworked, dtype=np.float64)) < quiet_floor)
        )
        reworked[recent] = max(tt - quiet_floor, -1.0)
        after_stable = cont & (age >= 2500.0) & (stability >= 0.75)
        world.set_g(
            "tectonics.last_p61_old_basement_current_state_candidate_share_of_continental_crust",
            float(area[old_stable_basement & ~active_blocked].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p61_shield_current_state_stabilized_share_of_continental_crust",
            float(area[shield_stabilize].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p61_platform_current_state_stabilized_share_of_continental_crust",
            float(area[platform_stabilize].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p61_active_blocked_share_of_continental_crust",
            float(area[active_blocked].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p61_stable_craton_fraction_before",
            float(area[before_stable].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p61_stable_craton_fraction_after",
            float(area[after_stable].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p61_current_state_stabilized_cells",
            float(np.count_nonzero(stabilized)),
        )
        return origin, age, thick, reworked, stability

    def _choose_internal_block_area(self, grid, mask: np.ndarray,
                                    score: np.ndarray,
                                    target_area: float,
                                    *,
                                    highest: bool = True) -> np.ndarray:
        chosen = np.zeros(grid.n, dtype=bool)
        cells = np.where(np.asarray(mask, dtype=bool))[0]
        if cells.size == 0 or target_area <= 0.0:
            return chosen
        score = np.asarray(score, dtype=np.float64)
        if score.shape != (grid.n,):
            score = np.zeros(grid.n, dtype=np.float64)
        order = cells[np.argsort(score[cells])]
        if highest:
            order = order[::-1]
        picked: list[int] = []
        acc = 0.0
        for c in order:
            picked.append(int(c))
            acc += float(grid.cell_area[int(c)])
            if acc >= target_area:
                break
        if picked:
            chosen[np.asarray(picked, dtype=np.int64)] = True
        return chosen

    def _initial_internal_geographic_block_code(
        self,
        grid,
        is_cont: np.ndarray,
        continent_label: np.ndarray,
        proto_craton: np.ndarray,
        width: np.ndarray,
        proto_potential: np.ndarray,
        margin_noise: np.ndarray,
        thick: np.ndarray,
        stability: np.ndarray,
    ) -> np.ndarray:
        """Seed inherited continental internal blocks from deterministic cargo.

        P104C keeps random seeds out of this layer.  Broad proto-crust potential
        supplies old cores, width marks continental interior versus margins,
        and thickness/stability separate covered platforms from sag basins and
        mobile belts.  Later plate motion inherits and locally rewrites this
        cargo.
        """
        is_cont = np.asarray(is_cont, dtype=bool)
        code = np.zeros(grid.n, dtype=np.float64)
        if not is_cont.any():
            return code
        labels = np.asarray(continent_label, dtype=np.int64)
        proto_craton = np.asarray(proto_craton, dtype=bool)
        width = np.asarray(width, dtype=np.float64)
        potential = _normalize01(proto_potential)
        margin = _normalize01(margin_noise)
        thick = np.asarray(thick, dtype=np.float64)
        stability = np.clip(np.asarray(stability, dtype=np.float64), 0.0, 1.0)
        area = grid.cell_area

        for cid in sorted(int(x) for x in np.unique(labels[is_cont]) if int(x) >= 0):
            component = is_cont & (labels == cid)
            if not component.any():
                continue
            comp = np.where(component)[0]
            comp_area = max(float(area[comp].sum()), 1.0e-12)
            assigned = np.zeros(grid.n, dtype=bool)

            core = component & proto_craton
            if not core.any():
                core_score = (
                    0.58 * potential
                    + 0.28 * np.clip(width / max(float(width[component].max()), 1.0),
                                     0.0, 1.0)
                    + 0.14 * stability
                )
                core = self._choose_internal_block_area(
                    grid,
                    component & (width >= 2.0),
                    core_score,
                    max(float(area[comp[0]]), 0.38 * comp_area),
                    highest=True,
                )
            code[core] = INTERNAL_BLOCK_CRATON_CORE
            assigned |= core

            if comp.size >= 6:
                margin_cut = float(np.percentile(margin[comp], 58.0))
                rift_candidate = (
                    component
                    & ~assigned
                    & (width <= 2.0)
                    & (margin >= margin_cut)
                )
                rift_score = (
                    0.70 * margin
                    + 0.22 * (1.0 - stability)
                    + 0.08 * _stable_tie_breaker(np.arange(grid.n))
                )
                rift = self._choose_internal_block_area(
                    grid,
                    rift_candidate,
                    rift_score,
                    max(float(area[comp[0]]), 0.055 * comp_area),
                    highest=True,
                )
                code[rift] = INTERNAL_BLOCK_RIFTED_MARGIN
                assigned |= rift

            if comp.size >= 7:
                mobile_cut = float(np.percentile(margin[comp], 50.0))
                mobile_candidate = (
                    component
                    & ~assigned
                    & (
                        (width <= 3.0)
                        | (margin >= mobile_cut)
                        | (stability <= float(np.percentile(stability[comp], 38.0)))
                    )
                )
                mobile_score = (
                    0.42 * margin
                    + 0.34 * (1.0 - stability)
                    + 0.18 * (1.0 - potential)
                    + 0.06 * _stable_tie_breaker(np.arange(grid.n))
                )
                mobile = self._choose_internal_block_area(
                    grid,
                    mobile_candidate,
                    mobile_score,
                    max(float(area[comp[0]]), 0.105 * comp_area),
                    highest=True,
                )
                code[mobile] = INTERNAL_BLOCK_MOBILE_BELT
                assigned |= mobile

            if comp.size >= 8:
                thick_p45 = float(np.percentile(thick[comp], 45.0))
                basin_candidate = (
                    component
                    & ~assigned
                    & (width >= 2.0)
                    & (
                        (thick <= thick_p45)
                        | (stability <= float(np.percentile(stability[comp], 48.0)))
                    )
                )
                basin_score = (
                    0.42 * (1.0 - _normalize01(thick))
                    + 0.30 * (1.0 - stability)
                    + 0.18 * margin
                    + 0.10 * np.clip(width / max(float(width[component].max()), 1.0),
                                     0.0, 1.0)
                )
                basin = self._choose_internal_block_area(
                    grid,
                    basin_candidate,
                    basin_score,
                    max(float(area[comp[0]]), 0.115 * comp_area),
                    highest=True,
                )
                code[basin] = INTERNAL_BLOCK_INTRACRATONIC_BASIN
                assigned |= basin

            platform = component & ~assigned
            code[platform] = INTERNAL_BLOCK_STABLE_PLATFORM

        code[~is_cont] = INTERNAL_BLOCK_NONE
        return code

    def _update_internal_geographic_block_code(
        self,
        world,
        grid,
        ctype: np.ndarray,
        origin: np.ndarray,
        age: np.ndarray,
        thick: np.ndarray,
        reworked: np.ndarray,
        stability: np.ndarray,
        domain: np.ndarray,
        province_code: np.ndarray,
        internal_block_code: np.ndarray,
        platform_subsidence: np.ndarray,
        current_time_myr: float,
        *,
        active_boundary_mask: np.ndarray,
    ) -> np.ndarray:
        del world, thick
        cont = np.asarray(ctype, dtype=np.float64) == CONT
        code = np.asarray(internal_block_code, dtype=np.float64).copy()
        if code.shape != (grid.n,):
            code = np.zeros(grid.n, dtype=np.float64)
        code[~cont] = INTERNAL_BLOCK_NONE
        if not cont.any():
            return code

        origin = np.asarray(origin, dtype=np.float64)
        stability = np.clip(np.asarray(stability, dtype=np.float64), 0.0, 1.0)
        age = np.asarray(age, dtype=np.float64)
        reworked = np.asarray(reworked, dtype=np.float64)
        domain = np.asarray(domain, dtype=np.float64)
        province = np.asarray(province_code, dtype=np.float64)
        if province.shape != (grid.n,):
            province = np.zeros(grid.n, dtype=np.float64)
        subsidence = np.clip(
            np.nan_to_num(
                np.asarray(platform_subsidence, dtype=np.float64),
                nan=0.0,
                posinf=1.0,
                neginf=0.0,
            ),
            0.0,
            1.0,
        )
        if subsidence.shape != (grid.n,):
            subsidence = np.zeros(grid.n, dtype=np.float64)
        active = np.asarray(active_boundary_mask, dtype=bool)
        if active.shape != (grid.n,):
            active = np.zeros(grid.n, dtype=bool)
        width = _graph_width_steps(grid, cont)
        tt = float(current_time_myr)
        quiet_age = np.where(
            reworked < 0.0,
            tt,
            np.maximum(tt - reworked, 0.0),
        )
        recent_rework = cont & (quiet_age < 520.0)
        stable_core = (
            cont
            & (code == INTERNAL_BLOCK_CRATON_CORE)
            & ((origin == ORIGIN_CRATON) | (stability >= 0.76))
        )

        code[
            cont
            & (province == PROVINCE_SHIELD)
            & (width >= 2.0)
            & (stability >= 0.68)
        ] = INTERNAL_BLOCK_CRATON_CORE
        code[
            cont
            & ~stable_core
            & (province == PROVINCE_INTERIOR_BASIN)
        ] = INTERNAL_BLOCK_INTRACRATONIC_BASIN
        code[
            cont
            & ~stable_core
            & (subsidence >= 0.52)
            & (width >= 3.0)
        ] = INTERNAL_BLOCK_INTRACRATONIC_BASIN
        code[
            cont
            & ~stable_core
            & np.isin(domain, [DOMAIN_SUTURE])
            & (active | recent_rework | (stability < 0.34))
        ] = INTERNAL_BLOCK_MOBILE_BELT
        code[
            cont
            & ~stable_core
            & np.isin(domain, [DOMAIN_ACCRETED_TERRANE])
            & (active | recent_rework | ((origin == ORIGIN_ARC) & (stability < 0.48)))
        ] = INTERNAL_BLOCK_ACCRETED_TERRANE
        code[
            cont
            & ~stable_core
            & (origin == ORIGIN_ARC)
            & (active | recent_rework | (stability < 0.48))
        ] = INTERNAL_BLOCK_ACCRETED_TERRANE
        code[
            cont
            & ~stable_core
            & (domain == DOMAIN_CONTINENTAL_MARGIN)
            & (width <= 3.0)
            & (active | recent_rework | (code == INTERNAL_BLOCK_RIFTED_MARGIN))
            & (subsidence < 0.42)
        ] = INTERNAL_BLOCK_RIFTED_MARGIN
        code[
            cont
            & ~stable_core
            & active
            & (domain != DOMAIN_CRATON)
            & (
                recent_rework
                | (stability < 0.38)
                | (origin == ORIGIN_SUTURE)
                | (code == INTERNAL_BLOCK_MOBILE_BELT)
            )
        ] = INTERNAL_BLOCK_MOBILE_BELT
        code[
            cont
            & (code == INTERNAL_BLOCK_NONE)
            & (province == PROVINCE_PLATFORM)
        ] = INTERNAL_BLOCK_STABLE_PLATFORM
        code[
            cont
            & (code == INTERNAL_BLOCK_NONE)
            & (width <= 2.0)
        ] = INTERNAL_BLOCK_RIFTED_MARGIN
        code[
            cont
            & (code == INTERNAL_BLOCK_NONE)
        ] = INTERNAL_BLOCK_STABLE_PLATFORM
        area = grid.cell_area
        cont_area = max(float(area[cont].sum()), 1.0e-12)
        for comp in self._connected_components(grid, cont):
            comp_area = max(float(area[comp].sum()), 1.0e-12)
            if comp_area / cont_area < 0.08 or comp.size < 8:
                continue
            comp_mask = np.zeros(grid.n, dtype=bool)
            comp_mask[comp] = True
            mobile_like = (
                comp_mask
                & np.isin(code, [
                    int(INTERNAL_BLOCK_MOBILE_BELT),
                    int(INTERNAL_BLOCK_ACCRETED_TERRANE),
                ])
            )
            mobile_area = float(area[mobile_like].sum())
            if mobile_area / comp_area <= 0.58:
                continue
            excess_area = mobile_area - 0.52 * comp_area
            wide_mobile = mobile_like & (width >= 2.0)
            if not wide_mobile.any() or excess_area <= 0.0:
                continue
            basin_score = (
                0.46 * subsidence
                + 0.24 * (1.0 - stability)
                + 0.18 * np.clip(width / max(float(width[comp].max()), 1.0),
                                 0.0, 1.0)
                + 0.12 * _stable_tie_breaker(np.arange(grid.n))
            )
            basin = self._choose_internal_block_area(
                grid,
                wide_mobile,
                basin_score,
                min(0.13 * comp_area, 0.46 * excess_area),
                highest=True,
            )
            code[basin] = INTERNAL_BLOCK_INTRACRATONIC_BASIN
            wide_mobile &= ~basin
            platform_score = (
                0.58 * np.clip(width / max(float(width[comp].max()), 1.0),
                               0.0, 1.0)
                + 0.24 * stability
                + 0.18 * (1.0 - subsidence)
            )
            platform = self._choose_internal_block_area(
                grid,
                wide_mobile,
                platform_score,
                min(0.20 * comp_area, 0.72 * excess_area),
                highest=True,
            )
            code[platform] = INTERNAL_BLOCK_STABLE_PLATFORM
        return code

    def _internal_geographic_block_ids_and_objects(
        self,
        world,
        grid,
        ctype: np.ndarray,
        internal_block_code: np.ndarray,
        previous_id: np.ndarray,
        next_id: int,
        age: np.ndarray,
        thick: np.ndarray,
        stability: np.ndarray,
        continent_id: np.ndarray,
        t: float,
    ) -> tuple[np.ndarray, list[dict], int]:
        cont = np.asarray(ctype, dtype=np.float64) == CONT
        code = np.asarray(internal_block_code, dtype=np.float64).astype(int)
        if code.shape != (grid.n,):
            code = np.zeros(grid.n, dtype=int)
        code[~cont] = int(INTERNAL_BLOCK_NONE)
        previous = np.asarray(previous_id, dtype=np.float64)
        if previous.shape != (grid.n,):
            previous = np.zeros(grid.n, dtype=np.float64)
        labels = np.zeros(grid.n, dtype=np.float64)
        objects: list[dict] = []
        used_ids: set[int] = set()
        total_area = max(float(grid.cell_area.sum()), 1.0e-12)
        cont_area = max(float(grid.cell_area[cont].sum()), 1.0e-12)
        width = _graph_width_steps(grid, cont)

        for block_code in sorted(
            int(x) for x in np.unique(code[cont]) if int(x) > 0
        ):
            mask = cont & (code == block_code)
            comps = self._connected_components(grid, mask)
            comps.sort(key=lambda c: float(grid.cell_area[c].sum()), reverse=True)
            for index, comp in enumerate(comps):
                old = previous[comp].astype(int)
                old = old[old > 0]
                assigned: int | None = None
                if old.size:
                    values, counts = np.unique(old, return_counts=True)
                    for idx in np.argsort(counts)[::-1]:
                        candidate = int(values[idx])
                        if candidate not in used_ids:
                            assigned = candidate
                            break
                if assigned is None:
                    assigned = int(next_id)
                    next_id += 1
                used_ids.add(assigned)
                labels[comp] = float(assigned)
                comp_area = max(float(grid.cell_area[comp].sum()), 1.0e-12)
                parent_continents = sorted(
                    int(x) for x in np.unique(continent_id[comp].astype(int))
                    if int(x) >= 0
                )
                lat, lon = self._cell_centroid_lat_lon(grid, comp)
                kind = INTERNAL_BLOCK_NAMES.get(block_code, "unknown")
                objects.append({
                    "id": f"internal_geographic_block:{assigned}",
                    "type": "internal_geographic_block",
                    "block_id": int(assigned),
                    "block_code": int(block_code),
                    "kind": kind,
                    "parent_process": {
                        int(INTERNAL_BLOCK_CRATON_CORE): "early_cratonization",
                        int(INTERNAL_BLOCK_STABLE_PLATFORM): "stable_platform_cover",
                        int(INTERNAL_BLOCK_INTRACRATONIC_BASIN): "intracratonic_sag",
                        int(INTERNAL_BLOCK_MOBILE_BELT): "orogenic_reworking",
                        int(INTERNAL_BLOCK_RIFTED_MARGIN): "continental_extension",
                        int(INTERNAL_BLOCK_ACCRETED_TERRANE): "arc_or_terrane_accretion",
                    }.get(block_code, "unknown"),
                    "parent_continent_ids": parent_continents,
                    "cell_count": int(comp.size),
                    "area_fraction": float(comp_area / total_area),
                    "share_of_continental_crust": float(comp_area / cont_area),
                    "centroid_lat": lat,
                    "centroid_lon": lon,
                    "mean_age_myr": float(np.average(age[comp], weights=grid.cell_area[comp])),
                    "mean_thickness_m": float(np.average(thick[comp], weights=grid.cell_area[comp])),
                    "mean_stability": float(np.average(stability[comp], weights=grid.cell_area[comp])),
                    "width_p50_steps": float(np.percentile(width[comp], 50.0)),
                    "last_active_myr": round(float(t), 1),
                    "cells": comp.astype(int).tolist() if comp.size <= 900 else [],
                    "component_index_in_class": int(index),
                })

        class_count = len({int(x) for x in np.unique(code[cont]) if int(x) > 0})
        coverage = float(grid.cell_area[cont & (code > 0)].sum() / cont_area)
        major_class_counts: list[int] = []
        largest_fractions: list[float] = []
        cid_arr = np.asarray(continent_id, dtype=np.float64)
        for cid in sorted(int(x) for x in np.unique(cid_arr[cont].astype(int)) if int(x) >= 0):
            mask = cont & (cid_arr.astype(int) == cid)
            comp_area = float(grid.cell_area[mask].sum())
            if comp_area / cont_area < 0.08:
                continue
            values = [int(x) for x in np.unique(code[mask]) if int(x) > 0]
            major_class_counts.append(len(values))
            if values:
                largest_fractions.append(max(
                    float(grid.cell_area[mask & (code == value)].sum() / comp_area)
                    for value in values
                ))
        world.set_g("tectonics.last_p104c_internal_block_object_count",
                    float(len(objects)))
        world.set_g("tectonics.last_p104c_internal_block_class_count",
                    float(class_count))
        world.set_g("tectonics.last_p104c_internal_block_coverage_fraction",
                    coverage)
        world.set_g(
            "tectonics.last_p104c_min_internal_block_class_count_per_major_continent",
            float(min(major_class_counts) if major_class_counts else 0),
        )
        world.set_g(
            "tectonics.last_p104c_max_largest_internal_block_class_fraction",
            float(max(largest_fractions) if largest_fractions else 0.0),
        )
        objects.sort(key=lambda obj: (-float(obj["area_fraction"]), int(obj["block_id"])))
        return labels, objects, int(next_id)

    def _cratonic_province_objects(self, world, grid, ctype, age, thick, origin,
                                   reworked, stability, domain, continent_id, t,
                                   platform_subsidence=None,
                                   province_code=None,
                                   basement_age_floor=None,
                                   basement_stability_floor=None,
                                   basement_thickness_floor=None):
        """Split mature continental interiors into shield, platform and basin objects.

        This is a tectonic provenance layer for terrain.  It does not mutate crust
        or land state; it names broad, quiet continental provinces from process
        proxies that already exist in the plate model.
        """
        out = {
            "tectonics.shields": [],
            "tectonics.platforms": [],
            "tectonics.interior_basins": [],
        }
        telemetry_keys = {
            "shield": (
                "tectonics.last_p48_shield_object_count",
                "tectonics.last_p48_shield_share_of_continental_crust",
            ),
            "platform": (
                "tectonics.last_p48_platform_object_count",
                "tectonics.last_p48_platform_share_of_continental_crust",
            ),
            "interior_basin": (
                "tectonics.last_p48_interior_basin_object_count",
                "tectonics.last_p48_interior_basin_share_of_continental_crust",
            ),
        }
        for count_key, share_key in telemetry_keys.values():
            world.set_g(count_key, 0.0)
            world.set_g(share_key, 0.0)
        for key in (
            "tectonics.last_p60_basement_old_candidate_share_of_continental_crust",
            "tectonics.last_p60_basement_shield_candidate_share_of_continental_crust",
            "tectonics.last_p60_basement_platform_candidate_share_of_continental_crust",
            "tectonics.last_p60_basement_sag_candidate_share_of_continental_crust",
        ):
            world.set_g(key, 0.0)

        cont = np.asarray(ctype, dtype=np.float64) == CONT
        if not cont.any():
            return out
        area = grid.cell_area
        total_area = max(float(area.sum()), 1e-12)
        cont_area = max(float(area[cont].sum()), 1e-12)
        width = _graph_width_steps(grid, cont)
        tt = max(float(t), 0.0)
        quiet_age = np.where(
            reworked < 0.0,
            tt,
            np.maximum(tt - np.asarray(reworked, dtype=np.float64), 0.0),
        )
        quiet_cut = min(1100.0, max(450.0, 0.18 * tt))
        mature_cut = min(2100.0, max(900.0, 0.36 * tt))
        shield_cut = min(2600.0, max(1700.0, 0.52 * tt))
        quiet = quiet_age >= quiet_cut
        province_memory = np.asarray(
            province_code if province_code is not None else np.zeros(grid.n),
            dtype=np.float64,
        )
        if province_memory.shape != (grid.n,):
            province_memory = np.zeros(grid.n, dtype=np.float64)
        memory_shield = province_memory == PROVINCE_SHIELD
        memory_platform = (
            (province_memory == PROVINCE_PLATFORM)
            | (province_memory == PROVINCE_INTERIOR_BASIN)
        )
        active_or_accreted = np.isin(
            domain,
            [DOMAIN_CONTINENTAL_MARGIN, DOMAIN_ACCRETED_TERRANE,
             DOMAIN_SUTURE, DOMAIN_LIP],
        )
        def optional_field(values) -> np.ndarray:
            arr = np.asarray(
                values if values is not None else np.zeros(grid.n),
                dtype=np.float64,
            )
            if arr.shape != (grid.n,):
                return np.zeros(grid.n, dtype=np.float64)
            return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)

        basement_age = optional_field(basement_age_floor)
        basement_stability = np.clip(optional_field(basement_stability_floor), 0.0, 1.0)
        basement_thick = optional_field(basement_thickness_floor)
        subsidence = np.zeros(grid.n, dtype=np.float64)
        if platform_subsidence is not None:
            subsidence = np.asarray(platform_subsidence, dtype=np.float64)
            if subsidence.shape != (grid.n,):
                subsidence = np.zeros(grid.n, dtype=np.float64)
            subsidence = np.clip(
                np.nan_to_num(subsidence, nan=0.0, posinf=1.0, neginf=0.0),
                0.0,
                1.0,
            )
        basement_platform_floor = min(tt, max(2500.0, 0.52 * tt)) if tt > 0.0 else 2500.0
        basement_shield_floor = min(tt, max(2500.0, 0.58 * tt)) if tt > 0.0 else 2500.0
        core_domain = ~active_or_accreted
        stable_passive_margin = (
            (domain == DOMAIN_CONTINENTAL_MARGIN)
            & (width >= 2.0)
            & (basement_stability >= 0.80)
        )
        basement_domain = core_domain | stable_passive_margin
        basement_old = (
            cont
            & basement_domain
            & (
                ((width >= 2.0) & core_domain)
                | ((width >= 2.0) & stable_passive_margin)
            )
            & (basement_age >= basement_platform_floor)
            & (basement_stability >= 0.70)
            & (basement_thick >= CONT_THICK + 1800.0)
        )
        basement_shield = (
            cont
            & core_domain
            & (width >= 2.0)
            & (basement_age >= basement_shield_floor)
            & (basement_stability >= 0.80)
            & (basement_thick >= CONT_THICK + 4200.0)
        )
        covered_strong_basement = basement_shield & (width >= 4.0) & (subsidence >= 0.35)
        basement_platform = (
            (basement_old & ~basement_shield)
            | covered_strong_basement
            | (basement_old & stable_passive_margin)
        )
        basement_shield &= ~covered_strong_basement
        basement_sag = (
            basement_platform
            & (width >= 4.0)
            & (
                (subsidence >= 0.45)
                | (basement_thick <= CONT_THICK + 2600.0)
            )
        )
        world.set_g(
            "tectonics.last_p60_basement_old_candidate_share_of_continental_crust",
            float(area[basement_old].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p60_basement_shield_candidate_share_of_continental_crust",
            float(area[basement_shield].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p60_basement_platform_candidate_share_of_continental_crust",
            float(area[basement_platform].sum() / cont_area),
        )
        world.set_g(
            "tectonics.last_p60_basement_sag_candidate_share_of_continental_crust",
            float(area[basement_sag].sum() / cont_area),
        )
        object_age = np.maximum(age, np.where(basement_old, basement_age, 0.0))
        object_thick = np.maximum(thick, np.where(basement_old, basement_thick, 0.0))
        object_stability = np.maximum(
            stability, np.where(basement_old, basement_stability, 0.0))
        shield_core_support = (
            (
                memory_shield
                & (
                    (width >= 4.0)
                    | (stability >= 0.80)
                    | (thick >= CONT_THICK + 4200.0)
                )
            )
            | (
                ((origin == ORIGIN_CRATON) | (domain == DOMAIN_CRATON))
                & (stability >= 0.82)
                & (thick >= CONT_THICK + 3600.0)
            )
            | (
                (stability >= 0.88)
                & (age >= shield_cut)
                & (thick >= CONT_THICK + 3200.0)
            )
            | basement_shield
        )
        shield = (
            cont
            & (width >= 3.0)
            & quiet
            & ((age >= shield_cut) | memory_shield | basement_shield)
            & shield_core_support
        )
        platform_context = (
            cont
            & (
                (width >= 4.0)
                | (basement_platform & (width >= 2.0))
                | (
                    (width >= 3.0)
                    & (
                        memory_platform
                        | basement_platform
                        | ((age >= mature_cut) & (stability >= 0.50))
                    )
                )
            )
            & quiet
            & ((age >= mature_cut) | memory_platform | basement_platform)
            & (
                (stability >= 0.34)
                | (memory_platform & (stability >= 0.26))
                | (basement_platform & (basement_stability >= 0.70))
            )
            & (~active_or_accreted | basement_platform)
        )
        shield_core_for_platform = shield & (width >= 4.0)
        platform = platform_context & ~shield_core_for_platform
        interior_basin = np.zeros(grid.n, dtype=bool)
        if platform.any():
            platform_thick = thick[platform]
            platform_stability = stability[platform]
            basin_cut = float(np.percentile(platform_thick, 45.0))
            stability_cut = max(
                0.70,
                min(0.94, float(np.percentile(platform_stability, 70.0))),
            )
            interior_basin = (
                platform
                & (thick <= basin_cut)
                & (stability <= stability_cut)
                & (width >= 4.0)
            )
            if platform_subsidence is not None:
                interior_basin |= (
                    platform
                    & (subsidence >= 0.46)
                    & (thick <= float(np.percentile(platform_thick, 68.0)))
                )
            interior_basin |= basement_sag & platform

        province_defs = (
            ("shield", "tectonics.shields", shield,
             "ancient_stable_cratonic_shield_core"),
            ("platform", "tectonics.platforms", platform,
             "mature_quiet_continental_platform"),
            ("interior_basin", "tectonics.interior_basins", interior_basin,
             "intracontinental_subsidence_on_mature_platform"),
        )
        for kind, set_name, mask, process in province_defs:
            objects = self._cratonic_province_component_objects(
                grid, kind, np.asarray(mask, dtype=bool), process,
                object_age, object_thick, object_stability, width, continent_id,
                t, total_area,
                max_objects=24 if kind == "platform" else 16,
            )
            out[set_name] = objects
            cells = np.zeros(grid.n, dtype=bool)
            for obj in objects:
                obj_cells = np.asarray(obj.get("cells", []), dtype=int)
                obj_cells = obj_cells[(0 <= obj_cells) & (obj_cells < grid.n)]
                cells[obj_cells] = True
            count_key, share_key = telemetry_keys[kind]
            world.set_g(count_key, float(len(objects)))
            world.set_g(share_key, float(area[cells].sum() / cont_area))
        return out

    def _cratonic_province_component_objects(self, grid, kind, mask, process,
                                             age, thick, stability, width,
                                             continent_id, t, total_area,
                                             *, max_objects=16):
        comps = self._connected_components(grid, np.asarray(mask, dtype=bool))
        comps.sort(key=lambda c: float(grid.cell_area[c].sum()), reverse=True)
        out = []
        min_cells = 2 if grid.n <= 3000 else 3
        for n, comp in enumerate(comps[:max_objects]):
            if comp.size < min_cells:
                continue
            cell_area = grid.cell_area[comp]
            comp_area = max(float(cell_area.sum()), 1e-12)
            parent_continents = sorted(
                int(x) for x in np.unique(continent_id[comp].astype(int))
                if int(x) >= 0
            )
            lat, lon = self._cell_centroid_lat_lon(grid, comp)
            out.append({
                "id": f"{kind}:{int(comp.min())}:{n}",
                "type": "cratonic_province",
                "kind": kind,
                "parent_process": process,
                "parent_continent_ids": parent_continents,
                "cell_count": int(comp.size),
                "area_fraction": float(comp_area / total_area),
                "centroid_lat": lat,
                "centroid_lon": lon,
                "mean_age_myr": float(np.average(age[comp], weights=cell_area)),
                "mean_thickness_m": float(np.average(thick[comp], weights=cell_area)),
                "mean_stability": float(np.average(stability[comp], weights=cell_area)),
                "width_p50_steps": float(np.percentile(width[comp], 50.0)),
                "last_active_myr": round(float(t), 1),
                "cells": comp.astype(int).tolist(),
            })
        return out

    def _cratonization_events(self, grid, cratonized, t):
        comps = self._connected_components(grid, np.asarray(cratonized, dtype=bool))
        comps.sort(key=lambda c: float(grid.cell_area[c].sum()), reverse=True)
        out = []
        for comp in comps[:3]:
            if comp.size < 4:
                continue
            out.append(Event(
                "cratonization", t, self.name, location=int(comp[0]),
                magnitude=float(grid.cell_area[comp].sum() / grid.cell_area.sum()),
                params={"cells": int(comp.size)},
            ))
        return out

    def _plume_activity(self, grid, ctype, thick, origin, reworked, volc_age,
                        plume_potential, t, dt, regime_code, vigor):
        plume_objects = []
        lip_objects = []
        events = []
        plume_mask = np.zeros(grid.n, dtype=bool)
        lip_mask = np.zeros(grid.n, dtype=bool)
        if regime_code < 1.0 or t < 250.0:
            return plume_mask, lip_mask, plume_objects, lip_objects, events
        if int((t - dt) // 550.0) == int(t // 550.0):
            return plume_mask, lip_mask, plume_objects, lip_objects, events

        score = np.asarray(plume_potential, dtype=np.float64).copy()
        if score.shape != (grid.n,):
            score = np.zeros(grid.n, dtype=np.float64)
        score = np.clip(score, 0.0, 1.0)
        eligible = (origin != ORIGIN_CRATON) | (score >= 0.72)
        local_max = np.zeros(grid.n, dtype=bool)
        for c in range(grid.n):
            if not eligible[c]:
                continue
            nbrs = grid.neighbors[c]
            if nbrs.size == 0 or score[c] >= float(np.max(score[nbrs])) - 1e-12:
                local_max[c] = True
        threshold = max(0.20, float(np.percentile(score[eligible], 98.0))
                        if np.any(eligible) else 1.0)
        candidates = np.where(local_max & eligible & (score >= threshold))[0]
        if candidates.size == 0:
            ranked = np.where(eligible & (score >= 0.20))[0]
            candidates = ranked
        if candidates.size == 0:
            return plume_mask, lip_mask, plume_objects, lip_objects, events
        n_plumes = 1 if vigor < 1.35 else 2
        order = sorted(((-float(score[c]), int(c)) for c in candidates))
        picks: list[int] = []
        min_spacing = np.cos(np.radians(18.0))
        for _, c in order:
            if any(float(grid.xyz[c] @ grid.xyz[p]) > min_spacing for p in picks):
                continue
            picks.append(int(c))
            if len(picks) >= n_plumes:
                break
        for k, seed in enumerate(picks):
            seed = int(seed)
            seed_mask = np.zeros(grid.n, dtype=bool)
            seed_mask[seed] = True
            region = self._dilate_mask(grid, seed_mask, passes=2)
            plume_mask |= region
            lip = region & (ctype == CONT)
            lip_mask |= lip
            thick[lip] = np.minimum(thick[lip] + 4500.0, MAX_CONT_THICK)
            origin[lip] = ORIGIN_PLUME_IMPACT
            reworked[lip] = t
            volc_age[region] = t
            plume_id = f"plume:{round(float(t), 1)}:{seed}"
            plume_objects.append({
                "id": plume_id,
                "kind": "mantle_plume",
                "cell": seed,
                "cells": np.where(region)[0].astype(int).tolist(),
                "birth_myr": round(float(t), 1),
                "last_active_myr": round(float(t), 1),
                "vigor": round(float(vigor), 3),
                "potential": round(float(score[seed]), 4),
            })
            lip_objects.append({
                "id": f"lip:{round(float(t), 1)}:{seed}",
                "kind": "large_igneous_province",
                "plume_id": plume_id,
                "cells": np.where(lip)[0].astype(int).tolist(),
                "area_fraction": float(grid.cell_area[lip].sum() / grid.cell_area.sum())
                if lip.any() else 0.0,
                "birth_myr": round(float(t), 1),
            })
            events.append(Event(
                "plume_head", t, self.name, location=seed,
                magnitude=float(region.sum()),
                params={"plume_id": plume_id, "vigor": round(float(vigor), 3),
                        "potential": round(float(score[seed]), 4)},
            ))
            if lip.any():
                events.append(Event(
                    "large_igneous_province", t, self.name, location=seed,
                    magnitude=float(grid.cell_area[lip].sum() / grid.cell_area.sum()),
                    params={"plume_id": plume_id, "cells": int(lip.sum())},
                ))
        return plume_mask, lip_mask, plume_objects, lip_objects, events

    def _merge_detached_plate_components(self, grid, plate, max_passes=2):
        """Attach disconnected plate-label fragments to adjacent plate bodies.

        The Eulerian raster step can leave one plate id in several disconnected
        islands.  Plate ids are identity fields for coherent rigid lids, so
        detached islands should be captured by a neighbouring plate instead of
        surviving as remote fragments with the same label.
        """
        out = plate.astype(int).copy()
        changed = 0
        for _ in range(max_passes):
            labels = out.copy()
            largest_body = np.zeros(grid.n, dtype=bool)
            detached: list[tuple[int, np.ndarray]] = []

            for pid in np.unique(labels):
                comps = self._connected_components(grid, labels == pid)
                if len(comps) <= 1:
                    if comps:
                        largest_body[comps[0]] = True
                    continue
                areas = [float(grid.cell_area[c].sum()) for c in comps]
                keep = int(np.argmax(areas))
                largest_body[comps[keep]] = True
                for k, comp in enumerate(comps):
                    if k != keep:
                        detached.append((int(pid), comp))

            if not detached:
                break

            pass_changed = 0
            for pid, comp in detached:
                target = self._best_adjacent_plate_for_component(
                    grid, labels, largest_body, comp, pid
                )
                if target is None:
                    continue
                out[comp] = target
                pass_changed += int(comp.size)
            changed += pass_changed
            if pass_changed == 0:
                break

        if changed:
            out = self._repair_plate_speckle(grid, out)
        return out, changed

    def _connected_components(self, grid, mask):
        nodes = np.where(mask)[0]
        if nodes.size == 0:
            return []
        seen = np.zeros(grid.n, dtype=bool)
        comps = []
        for start in nodes:
            if seen[start]:
                continue
            stack = [int(start)]
            seen[start] = True
            comp = []
            while stack:
                c = stack.pop()
                comp.append(c)
                for nb in grid.neighbors[c]:
                    nb = int(nb)
                    if mask[nb] and not seen[nb]:
                        seen[nb] = True
                        stack.append(nb)
            comps.append(np.asarray(comp, dtype=np.int64))
        return comps

    def _best_adjacent_plate_for_component(self, grid, labels, largest_body, comp, pid):
        scores: dict[int, int] = {}
        fallback: dict[int, int] = {}
        comp_mask = np.zeros(grid.n, dtype=bool)
        comp_mask[comp] = True
        for c in comp:
            for nb in grid.neighbors[int(c)]:
                nb = int(nb)
                target = int(labels[nb])
                if target == pid or comp_mask[nb]:
                    continue
                fallback[target] = fallback.get(target, 0) + 1
                if largest_body[nb]:
                    scores[target] = scores.get(target, 0) + 1
        source = scores if scores else fallback
        if not source:
            return None
        return max(source.items(), key=lambda item: (item[1], -item[0]))[0]

    def _local_reorganize_plates(self, world, grid, plate, plates, t):
        plate = plate.copy()
        plates = [dict(p) for p in plates]
        plate = self._repair_plate_speckle(grid, plate, max_passes=2)

        reorg_context = {
            "ctype": world.get_field("crust.type", OCEAN),
            "origin": world.get_field("crust.origin", ORIGIN_RIDGE),
            "age": world.get_field("crust.age_myr", 0.0),
            "stability": world.get_field("crust.stability", 0.0),
            "rift_potential": world.get_field("tectonics.rift_potential", 0.0),
            "boundaries": world.networks.get("tectonics.boundaries", {}),
            "plate_topologies": world.objects.get("tectonics.plate_topologies", []),
            "enable_p107_ranked_plate_policy": bool(
                world.g("tectonics.enable_p107_ranked_plate_policy", 0.0) > 0.0),
            "time_myr": float(t),
        }
        merges = self._merge_tiny_plates(grid, plate, reorg_context)
        active = set(int(x) for x in np.unique(plate))
        free_ids = [i for i in range(len(plates)) if i not in active]
        splits = []
        if free_ids:
            splits = self._split_large_plates(
                grid, plate, plates, free_ids, t, reorg_context)
        if reorg_context.get("enable_p107_ranked_plate_policy", False):
            active_now = set(int(x) for x in np.unique(plate))
            free_ids = [i for i in range(len(plates)) if i not in active_now]
            if free_ids:
                hierarchy_splits = self._split_underresolved_plate_hierarchy(
                    grid,
                    plate,
                    plates,
                    free_ids,
                    t,
                    reorg_context,
                    target_active=self._p107_target_active_plate_count(world),
                )
                splits.extend(hierarchy_splits)

        torque_refresh = self._refresh_plate_motions_from_torque(
            world, grid, plate, plates, active | {s["new_id"] for s in splits}, t)
        world.objects["tectonics.plates"] = plates
        active_after = sorted(int(x) for x in np.unique(plate))
        return plate, plates, {
            "mode": "local",
            "merged": merges,
            "protected_microplates": reorg_context.get("protected_microplates", []),
            "split": splits,
            "torque_refresh": torque_refresh,
            "active_plates": active_after,
        }

    @staticmethod
    def _unit(vec):
        vec = np.asarray(vec, dtype=np.float64)
        norm = float(np.linalg.norm(vec))
        if norm <= 1e-12:
            return np.zeros(3, dtype=np.float64)
        return vec / norm

    @classmethod
    def _tangent_toward(cls, at, target):
        at = cls._unit(at)
        target = cls._unit(target)
        tangent = target - float(target @ at) * at
        return cls._unit(tangent)

    def _plate_centroid(self, grid, cells):
        cells = np.asarray(cells, dtype=int)
        if cells.size == 0:
            return np.array([1.0, 0.0, 0.0], dtype=np.float64)
        weights = grid.cell_area[cells]
        centroid = np.average(grid.xyz[cells], axis=0, weights=weights)
        norm = float(np.linalg.norm(centroid))
        if norm <= 1e-12:
            centroid = grid.xyz[int(cells[0])]
        return self._unit(centroid)

    @staticmethod
    def _cells_for_boundary(world_n, boundaries, *kinds):
        mask = np.zeros(world_n, dtype=bool)
        if not isinstance(boundaries, dict):
            return mask
        for kind in kinds:
            cells = np.asarray(boundaries.get(kind, []), dtype=int)
            cells = cells[(0 <= cells) & (cells < world_n)]
            mask[cells] = True
        return mask

    def _initial_plate_motion_boundaries(self, grid, plate, ctype, age):
        """Infer a deterministic boundary proxy before real Wilson-cycle objects exist."""
        plate = np.asarray(plate, dtype=int)
        ctype = np.asarray(ctype, dtype=np.float64)
        age = np.asarray(age, dtype=np.float64)
        ocean_age = age[ctype == OCEAN]
        if ocean_age.size:
            young_cut = float(np.percentile(ocean_age, 34.0))
            old_cut = float(np.percentile(ocean_age, 66.0))
        else:
            young_cut = 10.0
            old_cut = 60.0
        age_span = max(old_cut - young_cut, 1.0)
        buckets: dict[str, list[int]] = {
            "ridge": [],
            "divergent": [],
            "trench": [],
            "subduction": [],
            "collision": [],
            "suture": [],
            "active_margin": [],
            "transform": [],
        }

        def add(kind: str, *cells: int) -> None:
            buckets[kind].extend(int(c) for c in cells)

        i = grid.edges[:, 0].astype(int)
        j = grid.edges[:, 1].astype(int)
        contacts = plate[i] != plate[j]
        for a, b in zip(i[contacts], j[contacts]):
            a = int(a)
            b = int(b)
            a_cont = bool(ctype[a] == CONT)
            b_cont = bool(ctype[b] == CONT)
            if a_cont and b_cont:
                add("collision", a, b)
                add("suture", a, b)
                continue
            if a_cont != b_cont:
                ocean_cell = b if a_cont else a
                cont_cell = a if a_cont else b
                add("trench", ocean_cell)
                add("subduction", ocean_cell)
                add("active_margin", cont_cell)
                continue

            a_age = float(age[a])
            b_age = float(age[b])
            if min(a_age, b_age) <= young_cut:
                if a_age <= b_age:
                    add("ridge", a)
                    add("divergent", a)
                    if b_age <= old_cut:
                        add("ridge", b)
                        add("divergent", b)
                else:
                    add("ridge", b)
                    add("divergent", b)
                    if a_age <= old_cut:
                        add("ridge", a)
                        add("divergent", a)
            elif abs(a_age - b_age) >= 0.35 * age_span:
                older = a if a_age >= b_age else b
                add("trench", older)
                add("subduction", older)
            else:
                add("transform", a, b)

        out: dict[str, list[int]] = {}
        for kind, cells in buckets.items():
            if cells:
                out[kind] = sorted(set(cells))
            else:
                out[kind] = []
        return out

    def _initial_plate_motions_from_torque(self, world, grid, plate, ctype, age,
                                           thick, proto_potential):
        """Create deterministic initial Euler poles/rates from process proxies."""
        plate = np.asarray(plate, dtype=int)
        active = sorted(int(x) for x in np.unique(plate))
        n_plates = max(int(world.spec.n_plates), max(active, default=-1) + 1)
        potential = _normalize01(np.asarray(proto_potential, dtype=np.float64))
        total_area = float(grid.cell_area.sum())
        proto_axis = self._unit(np.average(
            grid.xyz, axis=0, weights=np.maximum(potential, 0.02)))
        if not proto_axis.any():
            proto_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
        spin_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)

        previous: list[dict] = []
        for pid in range(n_plates):
            cells = np.where(plate == pid)[0]
            if cells.size:
                centroid = self._plate_centroid(grid, cells)
                pole = self._unit(
                    0.62 * np.cross(centroid, spin_axis)
                    + 0.38 * np.cross(centroid, proto_axis)
                )
                if not pole.any():
                    pole = self._unit(np.cross(centroid, np.array([1.0, 0.0, 0.0])))
                if not pole.any():
                    pole = spin_axis.copy()
                weights = grid.cell_area[cells]
                area_frac = float(weights.sum() / max(total_area, 1.0))
                local_ctype = np.asarray(ctype)[cells]
                cont_frac = float(np.average((local_ctype == CONT).astype(float),
                                             weights=weights))
                ocean = local_ctype == OCEAN
                if ocean.any():
                    ocean_age = float(np.average(np.asarray(age)[cells][ocean],
                                                 weights=weights[ocean]))
                    old_ocean = float(np.clip(ocean_age / 120.0, 0.0, 1.0))
                else:
                    old_ocean = 0.0
                mobility = (1.0 - 0.55 * cont_frac) * (0.55 + 0.45 * old_ocean)
                size_drag = 1.0 / np.sqrt(max(area_frac, 0.025))
                rate = (
                    self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"]
                    + 0.0022 * mobility * size_drag
                )
            else:
                pole = spin_axis.copy()
                rate = self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"]
            previous.append({
                "id": int(pid),
                "pole": pole.tolist(),
                "rate": float(np.clip(
                    rate,
                    self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"],
                    self.R2_PARAMETERS["max_rotation_rate_rad_per_myr"],
                )),
                "motion_source": "p53_initial_process_fallback",
            })

        boundaries = self._initial_plate_motion_boundaries(grid, plate, ctype, age)
        solution = self._torque_proxy_plate_motions(
            grid,
            plate,
            ctype,
            age,
            thick,
            boundaries,
            previous_plates=previous,
            active_ids=active,
            vigor=world.g("interior.tectonic_vigor", 1.0),
            regime_code=world.g("tectonics.regime_code", 3.0),
        )
        source_map = {
            "r2_torque_proxy": "p53_initial_torque_proxy",
            "r2_collision_locked": "p53_initial_collision_locked",
            "fallback_previous_motion": "p53_initial_process_fallback",
        }
        plates = [dict(p) for p in previous]
        force_rows = solution["diagnostics"].get("force_components", {})
        for pid, update in solution["updates"].items():
            if pid < 0 or pid >= len(plates):
                continue
            cells = np.where(plate == int(pid))[0]
            if cells.size:
                weights = grid.cell_area[cells]
                cont_frac = float(np.average(
                    (np.asarray(ctype)[cells] == CONT).astype(float),
                    weights=weights,
                ))
            else:
                cont_frac = 0.0
            min_rate = self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"]
            cargo_drag = 1.0 + 2.35 * (max(cont_frac, 0.0) ** 1.35)
            rate = min_rate + (float(update["rate"]) - min_rate) / cargo_drag
            source = source_map.get(str(update["source"]), f"p53_initial_{update['source']}")
            plates[pid].update({
                "pole": update["pole"],
                "rate": float(np.clip(
                    rate,
                    self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"],
                    self.R2_PARAMETERS["max_rotation_rate_rad_per_myr"],
                )),
                "motion_source": source,
                "p53_force_components": {
                    **force_rows.get(pid, {}),
                    "continental_cargo_fraction": cont_frac,
                    "continental_cargo_drag": cargo_drag,
                },
            })

        active_plates = [plates[pid] for pid in active if 0 <= pid < len(plates)]
        rates = np.asarray([float(p.get("rate", 0.0)) for p in active_plates],
                           dtype=np.float64)
        torque_driven = sum(
            1 for p in active_plates
            if str(p.get("motion_source", "")) == "p53_initial_torque_proxy"
        )
        fallback = sum(
            1 for p in active_plates
            if str(p.get("motion_source", "")) == "p53_initial_process_fallback"
        )
        boundary_counts = {
            f"boundary_proxy_{kind}_cells": float(len(cells))
            for kind, cells in boundaries.items()
        }
        telemetry = {
            "motion_active_plates": float(len(active_plates)),
            "motion_torque_driven_fraction": (
                float(torque_driven / len(active_plates)) if active_plates else 0.0
            ),
            "motion_fallback_fraction": (
                float(fallback / len(active_plates)) if active_plates else 0.0
            ),
            "motion_mean_rate_rad_per_myr": float(rates.mean()) if rates.size else 0.0,
            "motion_rate_cv": (
                float(rates.std() / max(float(rates.mean()), 1e-12))
                if rates.size else 0.0
            ),
            "motion_net_torque_ratio": float(
                solution["diagnostics"].get("net_torque_ratio", 0.0)),
            **boundary_counts,
        }
        return plates, {
            "schema": "aevum.tectonics.initial_plate_motion.p53.v1",
            "boundaries": boundaries,
            "diagnostics": solution["diagnostics"],
            "telemetry": telemetry,
        }

    def _torque_proxy_plate_motions(self, grid, plate, ctype, age, thick,
                                    boundaries, previous_plates=None,
                                    active_ids=None, vigor=1.0,
                                    regime_code=3.0):
        plate = np.asarray(plate, dtype=int)
        ctype = np.asarray(ctype, dtype=np.float64)
        age = np.asarray(age, dtype=np.float64)
        thick = np.asarray(thick, dtype=np.float64)
        if active_ids is None:
            active = sorted(int(x) for x in np.unique(plate))
        else:
            active = sorted(int(x) for x in active_ids)
        total_area = float(grid.cell_area.sum())
        regime_factor = {0.0: 0.02, 1.0: 0.4, 2.0: 0.7, 3.0: 1.0}.get(float(regime_code), 1.0)

        ridge_mask = self._cells_for_boundary(grid.n, boundaries, "ridge", "divergent")
        trench_mask = self._cells_for_boundary(grid.n, boundaries, "trench", "subduction")
        collision_mask = self._cells_for_boundary(grid.n, boundaries, "collision", "suture")
        transform_mask = self._cells_for_boundary(grid.n, boundaries, "transform")

        updates: dict[int, dict] = {}
        force_rows: dict[int, dict[str, float]] = {}
        net_torque = np.zeros(3, dtype=np.float64)
        total_component_norm = 0.0

        for pid in active:
            cells = np.where(plate == pid)[0]
            if cells.size == 0:
                continue
            centroid = self._plate_centroid(grid, cells)
            area_frac = float(grid.cell_area[cells].sum() / max(total_area, 1.0))
            torque = np.zeros(3, dtype=np.float64)
            net_force = np.zeros(3, dtype=np.float64)
            components = {
                "slab_pull": 0.0,
                "ridge_push": 0.0,
                "collision_resistance": 0.0,
                "basal_drag": 0.0,
                "transform_friction": 0.0,
            }

            def add_force(mask, strength_values, direction_sign, component):
                nonlocal torque, net_force
                bcells = np.where(mask & (plate == pid))[0]
                if bcells.size == 0:
                    return
                for c in bcells:
                    c = int(c)
                    toward_centroid = self._tangent_toward(grid.xyz[c], centroid)
                    direction = direction_sign * toward_centroid
                    strength = float(strength_values[c])
                    weight = float(grid.cell_area[c] / max(total_area, 1.0))
                    force = strength * weight * direction
                    net_force += force
                    torque += np.cross(grid.xyz[c], force)
                    components[component] += abs(strength) * weight

            oceanic = ctype == OCEAN
            continental = ctype == CONT
            old_ocean = np.clip((age - 45.0) / 155.0, 0.0, 1.0)
            slab_strength = (
                self.R2_PARAMETERS["slab_pull_weight"]
                * (0.25 + 0.75 * old_ocean)
                * np.where(oceanic, 1.0, 0.20)
            )
            add_force(trench_mask, slab_strength, -1.0, "slab_pull")

            young_ocean = np.clip((85.0 - age) / 85.0, 0.0, 1.0)
            ridge_strength = (
                self.R2_PARAMETERS["ridge_push_weight"]
                * (0.55 + 0.45 * young_ocean)
                * np.where(oceanic, 1.0, 0.45)
            )
            add_force(ridge_mask, ridge_strength, 1.0, "ridge_push")

            collision_strength = (
                self.R2_PARAMETERS["collision_resistance_weight"]
                * np.where(continental, 1.0, 0.25)
                * np.clip(thick / 62000.0, 0.30, 1.25)
            )
            collision_cells = np.where(collision_mask & (plate == pid))[0]
            if collision_cells.size:
                components["collision_resistance"] = float(np.sum(
                    collision_strength[collision_cells]
                    * grid.cell_area[collision_cells]
                    / max(total_area, 1.0)
                ))

            transform_cells = np.where(transform_mask & (plate == pid))[0]
            transform_load = float(grid.cell_area[transform_cells].sum() / max(total_area, 1.0))
            components["transform_friction"] = (
                self.R2_PARAMETERS["transform_friction_weight"] * transform_load
            )
            components["basal_drag"] = self.R2_PARAMETERS["basal_drag_weight"] * np.sqrt(
                max(area_frac, 1e-6))
            drag = 1.0 + components["basal_drag"] + components["transform_friction"]
            drag += 0.80 * components["collision_resistance"]

            torque_eff = torque / max(drag, 1e-9)
            torque_norm = float(np.linalg.norm(torque_eff))
            net_torque += torque_eff
            total_component_norm += torque_norm
            if torque_norm <= 1e-12:
                prev = previous_plates[pid] if previous_plates and pid < len(previous_plates) else {}
                old_pole = self._unit(prev.get("pole", [0.0, 0.0, 1.0]))
                if not old_pole.any():
                    old_pole = np.array([0.0, 0.0, 1.0], dtype=np.float64)
                old_rate = float(prev.get("rate", self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"]))
                pole = old_pole
                if components["collision_resistance"] > 0.0:
                    rate = self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"]
                    source = "r2_collision_locked"
                else:
                    rate = float(np.clip(
                        old_rate * self.R2_PARAMETERS["fallback_rate_decay"],
                        self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"],
                        self.R2_PARAMETERS["max_rotation_rate_rad_per_myr"],
                    ))
                    source = "fallback_previous_motion"
            else:
                pole = torque_eff / torque_norm
                mobility = max(float(vigor), 0.05) * regime_factor
                rate = (
                    self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"]
                    +
                    self.R2_PARAMETERS["torque_rate_scale_rad_per_myr"]
                    * mobility
                    * torque_norm
                    / max(area_frac, 0.025)
                )
                rate = float(np.clip(
                    rate,
                    self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"],
                    self.R2_PARAMETERS["max_rotation_rate_rad_per_myr"],
                ))
                if previous_plates and pid < len(previous_plates):
                    prev = previous_plates[pid]
                    old_pole = self._unit(prev.get("pole", pole))
                    if old_pole.any() and float(old_pole @ pole) > -0.65:
                        memory = self.R2_PARAMETERS["motion_memory_weight"]
                        pole = self._unit(memory * old_pole + (1.0 - memory) * pole)
                    old_rate = float(prev.get("rate", rate))
                    rate_memory = self.R2_PARAMETERS["rate_memory_weight"]
                    rate = float(np.clip(
                        rate_memory * old_rate + (1.0 - rate_memory) * rate,
                        self.R2_PARAMETERS["min_rotation_rate_rad_per_myr"],
                        self.R2_PARAMETERS["max_rotation_rate_rad_per_myr"],
                    ))
                source = "r2_torque_proxy"

            updates[pid] = {
                "pole": pole.tolist(),
                "rate": rate,
                "torque_norm": torque_norm,
                "net_force": net_force.tolist(),
                "source": source,
            }
            force_rows[pid] = {k: float(v) for k, v in components.items()}
            force_rows[pid]["area_fraction"] = area_frac
            force_rows[pid]["torque_norm"] = torque_norm
            force_rows[pid]["rate"] = rate
            force_rows[pid]["source"] = source

        diagnostics = {
            "active_plates": len(updates),
            "net_torque_norm": float(np.linalg.norm(net_torque)),
            "net_torque_ratio": float(np.linalg.norm(net_torque) / max(total_component_norm, 1e-12)),
            "force_components": force_rows,
            "parameters": dict(self.R2_PARAMETERS),
        }
        return {"updates": updates, "diagnostics": diagnostics}

    def _refresh_plate_motions_from_torque(self, world, grid, plate, plates, active_ids, t):
        solution = self._torque_proxy_plate_motions(
            grid,
            plate,
            world.get_field("crust.type", OCEAN),
            world.get_field("crust.age_myr", 0.0),
            world.get_field("crust.thickness_m", OCEAN_THICK),
            world.networks.get("tectonics.boundaries", {}),
            previous_plates=plates,
            active_ids=active_ids,
            vigor=world.g("interior.tectonic_vigor", 1.0),
            regime_code=world.g("tectonics.regime_code", 3.0),
        )
        for pid, update in solution["updates"].items():
            if pid < 0 or pid >= len(plates):
                continue
            plates[pid].update({
                "pole": update["pole"],
                "rate": update["rate"],
                "last_reorg_myr": round(t, 1),
                "motion_source": update["source"],
                "r2_force_components": solution["diagnostics"]["force_components"][pid],
            })
        return solution["diagnostics"]

    def _merge_tiny_plates(self, grid, plate, context=None):
        area = grid.cell_area
        total = float(area.sum())
        merges = []
        context = context or {}
        active = [int(x) for x in np.unique(plate)]
        if len(active) <= 1:
            return merges
        topologies = {
            int(obj.get("numeric_id", -1)): obj
            for obj in context.get("plate_topologies", []) or []
        }
        for pid in active:
            mask = plate == pid
            frac = float(area[mask].sum() / total)
            if frac >= self.MIN_PLATE_AREA_FRAC:
                continue
            topology = topologies.get(pid)
            protect_reason = self._p107_protected_microplate_reason(
                grid, plate, pid, frac, context, topology)
            if protect_reason:
                context.setdefault("protected_microplates", []).append({
                    "plate_id": int(pid),
                    "area_fraction": frac,
                    "reason": protect_reason,
                })
                continue
            target, edge_fraction, edge_count = self._dominant_capture_target(
                grid, plate, pid, context, topology)
            if target is None:
                continue
            source_cont = self._plate_continental_fraction(
                grid, plate, pid, context, topology)
            plate[mask] = target
            merges.append({
                "from": int(pid),
                "to": int(target),
                "area_fraction": frac,
                "basis": "p20_topology_shared_edge_capture",
                "shared_edge_fraction": round(float(edge_fraction), 4),
                "shared_edge_count": int(edge_count),
                "source_continental_fraction": round(float(source_cont), 4),
                "cargo_policy": (
                    "capture_plate_label_preserve_crust_cargo"
                    if source_cont >= 0.18 else "oceanic_microplate_capture"
                ),
            })
        return merges

    def _p107_protected_microplate_reason(self, grid, plate, pid, area_fraction,
                                          context, topology=None):
        if not bool(context.get("enable_p107_ranked_plate_policy", False)):
            return ""
        if float(area_fraction) < self.MICRO_PLATE_AREA_FRAC:
            return ""
        mask = plate == int(pid)
        if not np.any(mask):
            return ""
        n = grid.n
        ctype = np.asarray(context.get("ctype", np.full(n, OCEAN)), dtype=np.float64)
        origin = np.asarray(context.get("origin", np.full(n, ORIGIN_RIDGE)), dtype=np.float64)
        age = np.asarray(context.get("age", np.zeros(n)), dtype=np.float64)
        rift = np.asarray(context.get("rift_potential", np.zeros(n)), dtype=np.float64)
        if ctype.shape != (n,):
            ctype = np.full(n, OCEAN, dtype=np.float64)
        if origin.shape != (n,):
            origin = np.full(n, ORIGIN_RIDGE, dtype=np.float64)
        if age.shape != (n,):
            age = np.zeros(n, dtype=np.float64)
        if rift.shape != (n,):
            rift = np.zeros(n, dtype=np.float64)

        source_cont = self._plate_continental_fraction(grid, plate, pid, context, topology)
        boundary = topology.get("boundary_area_fractions", {}) if isinstance(topology, dict) else {}
        transform_load = float(boundary.get("transform", 0.0))
        arc_load = float(boundary.get("active_margin", 0.0) + boundary.get("trench", 0.0))
        collision_load = float(boundary.get("collision", 0.0) + boundary.get("suture", 0.0))
        rift_load = float(
            boundary.get("ridge", 0.0)
            + boundary.get("divergent", 0.0)
            + boundary.get("passive_margin", 0.0)
        )
        topo_continents = topology.get("continent_ids", []) if isinstance(topology, dict) else []
        topo_terranes = topology.get("terrane_ids", []) if isinstance(topology, dict) else []
        plate_area = max(float(grid.cell_area[mask].sum()), 1.0e-12)
        arc_origin_fraction = float(
            grid.cell_area[mask & (origin == ORIGIN_ARC)].sum() / plate_area)
        plume_origin_fraction = float(
            grid.cell_area[mask & (origin == ORIGIN_PLUME_IMPACT)].sum() / plate_area)
        ocean_fraction = float(grid.cell_area[mask & (ctype == OCEAN)].sum() / plate_area)
        mean_age = float(np.average(age[mask], weights=grid.cell_area[mask]))
        mean_rift = float(np.average(np.clip(rift[mask], 0.0, 1.0), weights=grid.cell_area[mask]))

        if source_cont >= 0.10 and (topo_continents or topo_terranes):
            return "microcontinent_or_terrane_cargo"
        if arc_origin_fraction >= 0.18 or arc_load >= 0.015:
            return "island_arc_or_active_margin_microplate"
        if collision_load >= 0.015 and source_cont >= 0.04:
            return "collision_wedge_microplate"
        if transform_load >= 0.010:
            return "transform_bounded_microplate"
        if rift_load >= 0.015 and mean_rift >= 0.18:
            return "rift_fragment_or_backarc_microplate"
        if ocean_fraction >= 0.70 and mean_age <= 70.0 and (rift_load + transform_load) >= 0.015:
            return "young_ocean_microplate"
        if plume_origin_fraction >= 0.25 and ocean_fraction >= 0.50:
            return "plume_or_plateau_microplate"
        return ""

    def _dominant_capture_target(self, grid, plate, pid, context, topology=None):
        current_counts: dict[int, float] = {}
        i, j = grid.edges[:, 0], grid.edges[:, 1]
        left = (plate[i] == pid) & (plate[j] != pid)
        right = (plate[j] == pid) & (plate[i] != pid)
        for target, count in zip(*np.unique(plate[j[left]], return_counts=True)):
            current_counts[int(target)] = current_counts.get(int(target), 0.0) + float(count)
        for target, count in zip(*np.unique(plate[i[right]], return_counts=True)):
            current_counts[int(target)] = current_counts.get(int(target), 0.0) + float(count)

        current_counts = {
            int(target): float(count)
            for target, count in current_counts.items()
            if int(target) != int(pid) and float(count) > 0.0
        }
        if not current_counts:
            return None, 0.0, 0
        total_edges = float(sum(current_counts.values()))
        target, edge_count = max(
            current_counts.items(),
            key=lambda item: (item[1], -item[0]),
        )
        edge_fraction = float(edge_count / max(total_edges, 1e-9))
        return int(target), edge_fraction, int(round(edge_count))

    def _plate_continental_fraction(self, grid, plate, pid, context, topology=None):
        mask = plate == int(pid)
        plate_area = float(grid.cell_area[mask].sum())
        if plate_area <= 0.0:
            return 0.0
        ctype = np.asarray(context.get("ctype", np.full(grid.n, OCEAN)), dtype=np.float64)
        if ctype.shape == (grid.n,):
            return float(grid.cell_area[mask & (ctype == CONT)].sum() / plate_area)
        if isinstance(topology, dict):
            return float(topology.get("continental_fraction", 0.0))
        return 0.0

    def _split_large_plates(self, grid, plate, plates, free_ids, t, context=None):
        area = grid.cell_area
        total = float(area.sum())
        context = context or {}
        splits = []
        for new_id in free_ids:
            active = np.unique(plate)
            candidates = []
            for pid in active:
                pid = int(pid)
                frac = float(area[plate == pid].sum() / total)
                if frac < self.MAX_PLATE_AREA_FRAC:
                    continue
                pressure = self._plate_split_pressure(pid, frac, context)
                candidates.append((pressure, frac, -pid, pid))
            if not candidates:
                break
            _, frac, _, parent = max(candidates)
            if parent < 0 or frac < self.MAX_PLATE_AREA_FRAC:
                break
            parent_mask = plate == parent
            cells = np.where(parent_mask)[0]
            if cells.size < 8:
                break
            score = self._plate_split_cell_score(grid, plate, parent, context)
            seed = self._deterministic_plate_split_seed(grid, parent_mask, score)
            target_area = float(area[parent_mask].sum() * 0.35)
            region = self._grow_region_within_plate(
                grid, plate, parent, seed, target_area, score=score)
            if region.sum() < 4:
                continue
            plate[region] = new_id
            parent_plate = plates[parent]
            parent_cells = np.where(plate == parent)[0]
            child_cells = np.where(plate == new_id)[0]
            parent_centroid = self._plate_centroid(grid, parent_cells)
            child_centroid = self._plate_centroid(grid, child_cells)
            split_axis = self._unit(np.cross(parent_centroid, child_centroid))
            parent_pole = self._unit(parent_plate.get("pole", [0.0, 0.0, 1.0]))
            pole = self._unit(0.82 * parent_pole + 0.18 * split_axis)
            if not pole.any():
                pole = parent_pole if parent_pole.any() else np.array([0.0, 0.0, 1.0])
            parent_area = float(area[parent_mask].sum())
            child_area = float(area[region].sum())
            child_share = child_area / max(parent_area, 1.0)
            rate_scale = float(np.clip(0.92 + 0.28 * (child_share - 0.35), 0.82, 1.10))
            plates[new_id].update({
                "pole": pole.tolist(),
                "rate": float(np.clip(parent_plate["rate"] * rate_scale,
                                      0.0025, 0.014)),
                "parent_id": int(parent),
                "birth_myr": round(t, 1),
                "last_reorg_myr": round(t, 1),
                "motion_source": "p20_topology_split_inherited",
                "split_basis": "topology_rift_boundary_score",
            })
            splits.append({
                "parent": int(parent),
                "new_id": int(new_id),
                "area_fraction": float(child_area / total),
                "parent_area_fraction_before": float(parent_area / total),
                "seed_cell": int(seed),
                "seed_score": round(float(score[seed]), 6),
                "basis": "p20_topology_rift_boundary_score",
            })
        return splits

    def _p107_target_active_plate_count(self, world):
        requested = max(int(getattr(world.spec, "n_plates", 12)), 1)
        if requested >= 48:
            return min(requested, max(28, int(np.ceil(0.50 * requested))))
        if requested >= 24:
            return min(requested, max(20, int(np.ceil(0.55 * requested))))
        return min(requested, max(10, int(np.ceil(0.60 * requested))))

    def _split_underresolved_plate_hierarchy(self, grid, plate, plates, free_ids, t,
                                             context=None, target_active=0):
        area = grid.cell_area
        total = max(float(area.sum()), 1.0e-12)
        context = context or {}
        splits = []
        target_active = int(target_active)
        for new_id in list(free_ids):
            active = set(int(x) for x in np.unique(plate) if int(x) >= 0)
            if len(active) >= target_active:
                break
            candidates = []
            for pid in active:
                mask = plate == pid
                frac = float(area[mask].sum() / total)
                if frac < max(0.055, 2.5 * self.MIN_PLATE_AREA_FRAC):
                    continue
                pressure = self._plate_split_pressure(pid, frac, context)
                pressure += 0.45 * np.sqrt(max(frac, 0.0))
                candidates.append((pressure, frac, -pid, pid))
            if not candidates:
                break
            _, frac, _, parent = max(candidates)
            parent_mask = plate == parent
            parent_area = float(area[parent_mask].sum())
            if parent_area <= 0.0:
                continue
            min_child_area = total * max(0.020, 2.4 * self.MIN_PLATE_AREA_FRAC)
            if parent_area <= min_child_area * 2.2:
                continue
            score = self._plate_split_cell_score(grid, plate, parent, context)
            seed = self._deterministic_plate_split_seed(grid, parent_mask, score)
            target_fraction = float(np.clip(0.20 + 0.10 * (frac / 0.18), 0.20, 0.32))
            target_area = max(
                parent_area * target_fraction,
                min_child_area,
            )
            target_area = min(target_area, parent_area * 0.36)
            region = self._grow_region_within_plate(
                grid, plate, parent, seed, target_area, score=score)
            if int(region.sum()) < 4:
                continue
            child_area = float(area[region].sum())
            child_frac = child_area / total
            if child_frac < max(0.018, 2.0 * self.MIN_PLATE_AREA_FRAC):
                continue
            plate[region] = new_id
            parent_plate = plates[parent]
            parent_cells = np.where(plate == parent)[0]
            child_cells = np.where(plate == new_id)[0]
            parent_centroid = self._plate_centroid(grid, parent_cells)
            child_centroid = self._plate_centroid(grid, child_cells)
            split_axis = self._unit(np.cross(parent_centroid, child_centroid))
            parent_pole = self._unit(parent_plate.get("pole", [0.0, 0.0, 1.0]))
            pole = self._unit(0.76 * parent_pole + 0.24 * split_axis)
            if not pole.any():
                pole = parent_pole if parent_pole.any() else np.array([0.0, 0.0, 1.0])
            rate_scale = float(np.clip(0.88 + 0.55 * child_frac, 0.84, 1.12))
            plates[new_id].update({
                "pole": pole.tolist(),
                "rate": float(np.clip(parent_plate["rate"] * rate_scale,
                                      0.0025, 0.014)),
                "parent_id": int(parent),
                "birth_myr": round(t, 1),
                "last_reorg_myr": round(t, 1),
                "motion_source": "p107_ranked_hierarchy_restoration",
                "split_basis": "p107_ranked_hierarchy_restoration",
            })
            splits.append({
                "parent": int(parent),
                "new_id": int(new_id),
                "area_fraction": float(child_frac),
                "parent_area_fraction_before": float(parent_area / total),
                "seed_cell": int(seed),
                "seed_score": round(float(score[seed]), 6),
                "basis": "p107_ranked_hierarchy_restoration",
                "target_active_plate_count": int(target_active),
            })
        return splits

    def _plate_split_pressure(self, pid, area_fraction, context):
        topologies = {
            int(obj.get("numeric_id", -1)): obj
            for obj in context.get("plate_topologies", []) or []
        }
        topo = topologies.get(int(pid), {})
        boundary = topo.get("boundary_area_fractions", {}) if isinstance(topo, dict) else {}
        rift_load = float(boundary.get("ridge", 0.0) + boundary.get("divergent", 0.0)
                          + boundary.get("passive_margin", 0.0)
                          + 0.5 * boundary.get("transform", 0.0))
        continent_load = float(topo.get("continental_fraction", 0.0)) if topo else 0.0
        excess = max(float(area_fraction) - self.MAX_PLATE_AREA_FRAC, 0.0)
        return float(excess + 0.35 * rift_load + 0.08 * continent_load)

    def _plate_split_cell_score(self, grid, plate, parent, context):
        parent_mask = plate == parent
        n = grid.n
        score = np.zeros(n, dtype=np.float64)

        rift = np.asarray(context.get("rift_potential", np.zeros(n)), dtype=np.float64)
        if rift.shape == (n,):
            rift = np.nan_to_num(rift, nan=0.0, posinf=1.0, neginf=0.0)
            score += 5.0 * np.clip(rift, 0.0, 1.0)

        boundaries = context.get("boundaries", {})
        boundary_source = self._cells_for_boundary(
            n, boundaries, "ridge", "divergent", "passive_margin", "transform")
        boundary_zone = self._dilate_mask(
            grid, boundary_source, allowed=parent_mask, passes=2)
        score += 2.7 * boundary_zone.astype(float)

        ctype = np.asarray(context.get("ctype", np.full(n, OCEAN)), dtype=np.float64)
        stability = np.asarray(context.get("stability", np.zeros(n)), dtype=np.float64)
        if ctype.shape == (n,):
            oceanic = ctype == OCEAN
            continental = ctype == CONT
            score += 0.9 * oceanic.astype(float)
            if stability.shape == (n,):
                weak_cont = continental & (stability < 0.55)
                stable_core = continental & (stability > 0.78)
                score += 1.3 * weak_cont.astype(float)
                score -= 5.5 * stable_core.astype(float)

        score[~parent_mask] = -np.inf
        if not np.isfinite(score[parent_mask]).any():
            score[parent_mask] = 0.0
        return score

    def _deterministic_plate_split_seed(self, grid, parent_mask, score):
        cells = np.where(parent_mask)[0]
        if cells.size == 0:
            return 0
        parent_centroid = self._plate_centroid(grid, cells)
        angular = np.degrees(np.arccos(np.clip(grid.xyz[cells] @ parent_centroid,
                                               -1.0, 1.0)))
        edge_bonus = np.clip(angular / 90.0, 0.0, 1.0)
        ranked = sorted(
            (float(score[c]) + 0.08 * float(edge_bonus[i]), -int(c), int(c))
            for i, c in enumerate(cells)
        )
        return int(ranked[-1][2])

    def _grow_region_within_plate(self, grid, plate, parent, seed, target_area,
                                  score=None):
        region = np.zeros(grid.n, dtype=bool)
        parent_mask = np.asarray(plate == parent, dtype=bool)
        seed = int(seed)
        if seed < 0 or seed >= grid.n or not parent_mask[seed]:
            cells = np.where(parent_mask)[0]
            if cells.size == 0:
                return region
            seed = int(cells[0])
        region[seed] = True
        acc = 0.0
        score_arr = (
            np.zeros(grid.n, dtype=np.float64)
            if score is None else np.nan_to_num(np.asarray(score, dtype=np.float64),
                                                nan=0.0, posinf=0.0, neginf=0.0)
        )
        seed_xyz = grid.xyz[int(seed)]
        acc += float(grid.cell_area[seed])
        mean_cell_area = max(float(np.mean(grid.cell_area[parent_mask])), 1.0)
        candidate_set: set[int] = set()
        for nb in grid.neighbors[seed]:
            nb = int(nb)
            if parent_mask[nb] and not region[nb]:
                candidate_set.add(nb)

        while candidate_set and acc < target_area:
            candidates = np.fromiter(candidate_set, dtype=np.int64)
            if candidates.size == 0:
                break
            remaining_cells = max(1, int(np.ceil((target_area - acc) / mean_cell_area)))
            batch_size = min(
                int(candidates.size),
                max(1, min(remaining_cells, max(8, int(np.ceil(remaining_cells * 0.16))))),
            )
            candidate_xyz = grid.xyz[candidates]
            angle = np.arccos(np.clip(candidate_xyz @ seed_xyz, -1.0, 1.0))
            priority = angle - 0.18 * score_arr[candidates] + 1.0e-8 * candidates
            if batch_size < candidates.size:
                chosen_idx = np.argpartition(priority, batch_size - 1)[:batch_size]
                chosen = candidates[chosen_idx]
            else:
                chosen = candidates
            region[chosen] = True
            acc += float(grid.cell_area[chosen].sum())
            for c in chosen:
                candidate_set.discard(int(c))
            for c in chosen:
                for nb in grid.neighbors[int(c)]:
                    nb = int(nb)
                    if parent_mask[nb] and not region[nb]:
                        candidate_set.add(nb)
        return region

    def _repair_plate_speckle(self, grid, plate, max_passes=1):
        out = plate.astype(int).copy()
        edges = grid.edges
        i, j = edges[:, 0], edges[:, 1]
        deg = np.zeros(grid.n, dtype=int)
        np.add.at(deg, i, 1)
        np.add.at(deg, j, 1)
        for _ in range(max_passes):
            ids = np.unique(out)
            counts = np.zeros((ids.size, grid.n), dtype=np.int16)
            for row, pid in enumerate(ids):
                mask = out == pid
                np.add.at(counts[row], i, mask[j])
                np.add.at(counts[row], j, mask[i])
            best_row = np.argmax(counts, axis=0)
            best = ids[best_row]
            best_count = counts[best_row, np.arange(grid.n)]
            threshold = np.maximum(3, deg - 1)
            nxt = np.where((best != out) & (best_count >= threshold), best, out)
            if np.array_equal(nxt, out):
                break
            out = nxt
        return out

    def _coherent_ridge_mask(self, grid, ridge_seed, ocean_seed):
        """Bridge one-cell gaps in oceanic ridge axes without flooding oceans."""
        ridge = np.asarray(ridge_seed, dtype=bool).copy()
        allowed = np.asarray(ocean_seed, dtype=bool) | ridge
        if not ridge.any():
            return ridge

        for _ in range(2):
            add = np.zeros(grid.n, dtype=bool)
            for c in np.where(allowed & ~ridge)[0]:
                nbs = grid.neighbors[c]
                if np.count_nonzero(ridge[nbs]) >= 2:
                    add[c] = True
            if not add.any():
                break
            ridge |= add

        # Keep ridge axes from becoming wide boundary swaths.  Cells with a
        # neighbouring ridge survive; isolated seed pixels are dropped unless
        # they are the only ridge evidence for this step.
        keep = ridge.copy()
        if int(ridge.sum()) > 8:
            keep[:] = False
            for c in np.where(ridge)[0]:
                keep[c] = np.count_nonzero(ridge[grid.neighbors[c]]) > 0
            ridge &= keep
        return ridge

    def _graph_distance_from_sources(self, grid, sources, allowed):
        """Shortest great-circle graph distance from sources through allowed cells."""
        sources = np.asarray(sources, dtype=bool)
        allowed = np.asarray(allowed, dtype=bool)
        dist = np.full(grid.n, np.inf, dtype=np.float64)
        heap: list[tuple[float, int]] = []
        for c in np.where(sources & allowed)[0]:
            dist[c] = 0.0
            heapq.heappush(heap, (0.0, int(c)))

        while heap:
            d, c = heapq.heappop(heap)
            if d != dist[c]:
                continue
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if not allowed[nb]:
                    continue
                nd = d + grid.great_circle_distance(c, nb)
                if nd < dist[nb]:
                    dist[nb] = nd
                    heapq.heappush(heap, (nd, nb))
        return dist

    def _mobile_lid_ocean_age_cap(self, regime_code, vigor):
        cap_by_regime = {1.0: 720.0, 2.0: 420.0, 3.0: 260.0}
        base_cap = cap_by_regime.get(float(regime_code), 420.0)
        return base_cap / np.clip(vigor, 0.6, 1.6) ** 0.35

    def _impose_seafloor_age_from_ridges(self, grid, ctype, age, reworked,
                                         ridge, regime_code, vigor, t):
        """Make oceanic age a field radiating away from current ridge axes."""
        if regime_code < 1.0:
            return age, reworked
        oceanic = ctype == OCEAN
        ridge_ocean = np.asarray(ridge, dtype=bool) & oceanic
        if not oceanic.any() or not ridge_ocean.any():
            return age, reworked

        dist = self._graph_distance_from_sources(grid, ridge_ocean, oceanic)
        reached = oceanic & np.isfinite(dist)
        if not reached.any():
            return age, reworked

        out = age.copy()
        rw = reworked.copy()
        regime_factor = {1.0: 0.55, 2.0: 0.78, 3.0: 1.0}.get(float(regime_code), 0.78)
        half_spread_m_per_myr = 36000.0 * regime_factor * np.clip(vigor, 0.55, 1.7) ** 0.55
        target = dist[reached] / max(half_spread_m_per_myr, 1.0)
        cap = self._mobile_lid_ocean_age_cap(regime_code, vigor)
        out[reached] = np.minimum(np.minimum(target, cap), t)
        out[ridge_ocean] = 0.0
        rw[ridge_ocean] = t
        return out, rw

    def _deforming_network_state(self, grid, collision, subduction, separating,
                                 transform, cont_within, plate_contact, ctype,
                                 t):
        """Represent active deformation separately from crust provenance.

        `collision` and `subduction` are broad rasterized process swaths.  This
        state keeps the shared plate-contact axis as the high-confidence core,
        while preserving surrounding broad process cells as lower-intensity
        shoulders without changing crust origin.
        """
        collision = np.asarray(collision, dtype=bool)
        subduction = np.asarray(subduction, dtype=bool)
        separating = np.asarray(separating, dtype=bool)
        transform = np.asarray(transform, dtype=bool)
        cont_within = np.asarray(cont_within, dtype=bool)
        plate_contact = np.asarray(plate_contact, dtype=bool)
        ctype = np.asarray(ctype, dtype=np.float64)

        active = collision | subduction | separating | transform
        intensity = np.zeros(grid.n, dtype=np.float64)
        style = np.full(grid.n, DEFORM_NONE, dtype=np.float64)
        if not active.any():
            return intensity, style, []

        collision_core = self._process_axis(grid, collision, plate_contact, spacing=2)
        collision_shoulder = collision & ~collision_core
        subduction_core = self._process_axis(grid, subduction, plate_contact, spacing=2)
        subduction_shoulder = subduction & ~subduction_core
        rift = self._process_axis(
            grid, separating & ~collision & ~subduction, plate_contact, spacing=3)
        transform_zone = self._process_axis(
            grid, transform & ~collision & ~subduction, plate_contact, spacing=2)

        intensity[collision_shoulder] = np.maximum(intensity[collision_shoulder], 0.45)
        style[collision_shoulder] = DEFORM_COLLISION_SHOULDER
        intensity[subduction_shoulder] = np.maximum(intensity[subduction_shoulder], 0.42)
        style[subduction_shoulder] = DEFORM_SUBDUCTION_SHOULDER
        intensity[rift] = np.maximum(intensity[rift], 0.65)
        style[rift] = DEFORM_RIFT
        intensity[transform_zone] = np.maximum(intensity[transform_zone], 0.55)
        style[transform_zone] = DEFORM_TRANSFORM
        intensity[subduction_core] = np.maximum(intensity[subduction_core], 0.90)
        style[subduction_core] = DEFORM_SUBDUCTION_CORE
        intensity[collision_core] = np.maximum(intensity[collision_core], 1.0)
        style[collision_core] = DEFORM_COLLISION_CORE

        objects = self._deforming_network_objects(
            grid, style, intensity, ctype, t)
        return intensity, style, objects

    def _process_axis(self, grid, process, plate_contact, *, spacing=2):
        process = np.asarray(process, dtype=bool)
        axis = process & np.asarray(plate_contact, dtype=bool)
        if axis.any():
            return self._thin_boundary_mask(grid, axis, spacing=spacing)
        return self._thin_boundary_mask(grid, process, spacing=spacing)

    def _process_provenance_mask(self, grid, process, plate_contact, *,
                                 spacing=2, shoulder_passes=1):
        """Localize crust-provenance rewrites to process axes.

        Broad deformation shoulders can affect relief and thickness, but they
        should not make an entire mature continent become a young suture or
        accreted terrane every timestep.  This returns a thin boundary axis plus
        a small shoulder, clipped to the process swath.
        """
        process = np.asarray(process, dtype=bool)
        if not process.any():
            return np.zeros(grid.n, dtype=bool)
        axis = self._process_axis(grid, process, plate_contact, spacing=spacing)
        axis &= process
        if shoulder_passes > 0 and axis.any():
            axis = self._dilate_mask(
                grid, axis, allowed=process, passes=shoulder_passes)
        return axis & process

    def _deforming_network_objects(self, grid, style, intensity, ctype, t):
        objects: list[dict] = []
        area = grid.cell_area
        total_area = float(area.sum())
        style_names = {
            int(DEFORM_COLLISION_CORE): "collision_core",
            int(DEFORM_COLLISION_SHOULDER): "collision_shoulder",
            int(DEFORM_SUBDUCTION_CORE): "subduction_core",
            int(DEFORM_SUBDUCTION_SHOULDER): "subduction_shoulder",
            int(DEFORM_RIFT): "rift",
            int(DEFORM_TRANSFORM): "transform",
        }
        for code, name in style_names.items():
            mask = style.astype(int) == int(code)
            if not mask.any():
                continue
            _, comps = _component_labels(grid, mask)
            for idx, cells in enumerate(comps):
                if cells.size == 0:
                    continue
                weights = area[cells]
                comp_area = float(weights.sum())
                centroid = np.average(grid.xyz[cells], axis=0, weights=weights)
                norm = float(np.linalg.norm(centroid))
                if norm > 1e-12:
                    centroid = centroid / norm
                lat = float(np.degrees(np.arcsin(np.clip(centroid[2], -1.0, 1.0))))
                lon = float(np.degrees(np.arctan2(centroid[1], centroid[0])))
                objects.append({
                    "id": f"deforming_network:{name}:{idx}",
                    "kind": name,
                    "style_code": int(code),
                    "cell_count": int(cells.size),
                    "area_fraction": comp_area / total_area if total_area else 0.0,
                    "mean_intensity": float(np.average(intensity[cells], weights=weights)),
                    "continental_fraction": float(np.mean(ctype[cells] == CONT)),
                    "lat": round(lat, 2),
                    "lon": round(lon, 2),
                    "last_active_myr": round(float(t), 1),
                    "cells": cells.astype(int).tolist() if cells.size <= 600 else [],
                })
        objects.sort(key=lambda obj: (-float(obj["area_fraction"]), str(obj["id"])))
        return objects[:64]

    def _dilate_mask(self, grid, mask, allowed=None, passes=1):
        out = np.asarray(mask, dtype=bool).copy()
        if allowed is None:
            allowed = np.ones(grid.n, dtype=bool)
        allowed = np.asarray(allowed, dtype=bool)
        for _ in range(passes):
            add = np.zeros(grid.n, dtype=bool)
            for c in np.where(out)[0]:
                add[grid.neighbors[c]] = True
            out |= add & allowed
        return out

    def _smooth_field(self, grid, values, *, allowed=None, passes=1, alpha=0.25):
        out = np.asarray(values, dtype=np.float64).copy()
        if allowed is None:
            allowed = np.ones(grid.n, dtype=bool)
        else:
            allowed = np.asarray(allowed, dtype=bool)
        for _ in range(int(passes)):
            nxt = out.copy()
            for c in np.where(allowed)[0]:
                nbrs = np.asarray(grid.neighbors[int(c)], dtype=int)
                nbrs = nbrs[allowed[nbrs]]
                if nbrs.size == 0:
                    continue
                nxt[int(c)] = (
                    (1.0 - float(alpha)) * out[int(c)]
                    + float(alpha) * float(np.mean(out[nbrs]))
                )
            out = nxt
            out[~allowed] = 0.0
        return out

    def _recycle_old_oceanic_crust(self, grid, ctype, age, reworked, trench,
                                   regime_code, vigor, t, rng):
        """Proxy for mobile-lid oceanic lithosphere turnover.

        Once the global redraw is removed, old oceanic crust needs an explicit
        recycling timescale.  Mobile lids should rarely preserve very old ocean
        floor; sluggish and episodic lids can keep it longer; stagnant lids skip
        this entirely.
        """
        if regime_code < 1.0:
            return age, reworked
        cap = self._mobile_lid_ocean_age_cap(regime_code, vigor)
        oceanic = ctype == OCEAN
        old = oceanic & (age > cap)
        if not old.any():
            return age, reworked
        out = age.copy()
        rw = reworked.copy()
        trench_zone = self._dilate_mask(grid, trench, allowed=oceanic, passes=3)
        recycled = old & trench_zone
        if recycled.any():
            # Slabs are consumed at trench margins; the fixed-grid proxy records
            # the rollover by replacing only trench-neighbour old oceanic cells.
            out[recycled] = rng.uniform(0.45 * cap, 0.95 * cap, size=int(recycled.sum()))
            rw[recycled] = t
        # Isolated stale ages can remain in disconnected basins with no current
        # ridge path.  Cap them deterministically without marking a false trench
        # event, so diagnostics stay physical while avoiding immortal ocean.
        stale = old & ~recycled
        out[stale] = np.minimum(out[stale], 0.98 * cap)
        return out, rw

    def _rifted_margin_candidates(self, grid, ctype, stability, rift_potential,
                                  margin_sources=None):
        """Deterministic breakup candidates from rift potential and weakness."""
        cont = ctype == CONT
        out = np.zeros(grid.n, dtype=bool)
        if not cont.any():
            return out
        potential = np.asarray(rift_potential, dtype=np.float64)
        if potential.shape != (grid.n,):
            potential = np.zeros(grid.n, dtype=np.float64)
        potential = np.nan_to_num(potential, nan=0.0, posinf=1.0, neginf=0.0)
        potential = np.clip(potential, 0.0, 1.0)
        cont_potential = potential[cont]
        if cont_potential.size == 0:
            return out
        threshold = max(0.48, float(np.percentile(cont_potential, 82.0)))
        weak = cont & (stability < 0.72)
        out = weak & (potential >= threshold)

        if margin_sources is not None:
            source = np.asarray(margin_sources, dtype=bool)
            if source.shape == (grid.n,) and source.any():
                margin_zone = self._dilate_mask(grid, source, allowed=cont, passes=2)
                out |= margin_zone & weak & (potential >= max(0.34, threshold * 0.72))
        return out

    def _continental_shape_pressure(self, grid, ctype, origin, stability) -> float:
        """Return a deterministic pressure to heal implausibly ribboned continents.

        Mature continents should not remain dominated by one- and two-cell
        accreted/suture strips.  This pressure is not a random reshaper; it is a
        trigger for area-neutral continental welding when unstable narrow
        margins outweigh broad interior area.
        """
        cont = ctype == CONT
        if not cont.any():
            return 0.0
        area = grid.cell_area
        cont_area = float(area[cont].sum())
        if cont_area <= 0.0:
            return 0.0
        width = _graph_width_steps(grid, cont)
        narrow = cont & (width <= 2.0)
        unsupported = narrow & (
            np.isin(origin, [ORIGIN_ARC, ORIGIN_SUTURE])
            | (stability < 0.55)
        )
        broad = cont & (width >= 4.0)
        narrow_fraction = float(area[narrow].sum() / cont_area)
        unsupported_fraction = float(area[unsupported].sum() / cont_area)
        broad_fraction = float(area[broad].sum() / cont_area)
        pressure = max(
            (narrow_fraction - 0.35) / 0.35,
            (unsupported_fraction - 0.14) / 0.30,
            (0.30 - broad_fraction) / 0.24,
            0.0,
        )
        return float(np.clip(pressure, 0.0, 1.0))

    def _conserve_continental(self, world, grid, ctype, thick, age,
                              origin, reworked, stability, orog_age, volc_age,
                              t, dt, rng_key, accretion_sources=None,
                              erosion_sources=None, basement_age=None,
                              basement_stability=None, basement_thick=None):
        area = grid.cell_area
        target = world.g("crust.cont_fraction_init", world.spec.target_land_fraction)
        total = area.sum()
        frac = area[ctype == CONT].sum() / total
        accretion_seed = (
            np.zeros(grid.n, dtype=bool) if accretion_sources is None
            else np.asarray(accretion_sources, dtype=bool)
        )
        erosion_seed = (
            np.zeros(grid.n, dtype=bool) if erosion_sources is None
            else np.asarray(erosion_sources, dtype=bool)
        )
        if accretion_seed.shape != (grid.n,):
            accretion_seed = np.zeros(grid.n, dtype=bool)
        if erosion_seed.shape != (grid.n,):
            erosion_seed = np.zeros(grid.n, dtype=bool)
        has_accretion_process = bool(accretion_seed.any())
        has_erosion_process = bool(erosion_seed.any())
        world.set_g("tectonics.last_continent_gain_cells", 0.0)
        world.set_g("tectonics.last_continent_loss_cells", 0.0)
        world.set_g("tectonics.last_unforced_continent_gain_blocked", 0.0)
        world.set_g("tectonics.last_unforced_continent_loss_blocked", 0.0)
        world.set_g("tectonics.last_background_continent_shape_maintenance", 0.0)
        world.set_g("tectonics.last_passive_margin_progradation_cells", 0.0)
        world.set_g("tectonics.last_p34_parented_platform_candidate_cells", 0.0)
        world.set_g("tectonics.last_p34_parented_platform_restored_cells", 0.0)
        world.set_g("tectonics.last_p34_parented_platform_area_fraction", 0.0)
        world.set_g("tectonics.last_p34_parented_platform_protected_cells", 0.0)
        world.set_g("tectonics.last_p35_continental_debt_candidate_cells", 0.0)
        world.set_g("tectonics.last_p35_continental_debt_restored_cells", 0.0)
        world.set_g("tectonics.last_p35_continental_debt_area_fraction", 0.0)
        world.set_g("tectonics.last_p35_stable_parent_cells", 0.0)
        world.set_g("tectonics.last_p35_old_platform_anchor_cells", 0.0)
        world.set_g("tectonics.last_p35_continental_gap_fraction", 0.0)
        world.set_g("tectonics.last_p37_shape_guard_input_cells", 0.0)
        world.set_g("tectonics.last_p37_shape_guard_accepted_cells", 0.0)
        world.set_g("tectonics.last_p37_shape_guard_rejected_cells", 0.0)
        world.set_g("tectonics.last_p38_lineage_guard_input_cells", 0.0)
        world.set_g("tectonics.last_p38_lineage_guard_accepted_cells", 0.0)
        world.set_g("tectonics.last_p38_lineage_guard_rejected_cells", 0.0)
        for key in P39_OLD_PLATFORM_FUNNEL_METRICS:
            world.set_g(f"tectonics.last_p39_old_platform_funnel_{key}", 0.0)
        world.set_g("tectonics.last_p40_parented_restoration_inherited_rework_cells", 0.0)
        world.set_g("tectonics.last_p40_parented_restoration_reactivated_cells", 0.0)
        world.set_g("tectonics.last_p40_shape_maintenance_inherited_rework_cells", 0.0)
        world.set_g("tectonics.last_p40_shape_maintenance_reactivated_cells", 0.0)
        world.set_g("tectonics.last_p41_parented_restoration_priority_cells", 0.0)
        world.set_g("tectonics.last_p41_active_only_accretion_deferred_cells", 0.0)
        world.set_g("tectonics.last_p43_mature_active_only_candidate_cells", 0.0)
        world.set_g("tectonics.last_p43_mature_active_only_accepted_cells", 0.0)
        world.set_g("tectonics.last_p43_mature_active_only_rejected_cells", 0.0)
        world.set_g("tectonics.last_p43_mature_active_only_accepted_area_fraction", 0.0)
        world.set_g("tectonics.last_p43_mature_active_only_cap_area_fraction", 0.0)
        world.set_g("tectonics.last_p104e_earthlike_floor_target_fraction", 0.0)
        world.set_g("tectonics.last_p104e_earthlike_floor_gap_fraction", 0.0)
        world.set_g("tectonics.last_p104e_earthlike_floor_candidate_cells", 0.0)
        world.set_g("tectonics.last_p104e_earthlike_floor_restored_cells", 0.0)
        world.set_g(
            "tectonics.last_p104e_earthlike_floor_restored_area_fraction", 0.0)
        world.set_g("tectonics.last_p44_parented_recovery_candidate_cells", 0.0)
        world.set_g("tectonics.last_p44_parented_recovery_restored_cells", 0.0)
        world.set_g("tectonics.last_p44_parented_recovery_area_fraction", 0.0)
        world.set_g("tectonics.last_p44_parented_recovery_protected_cells", 0.0)
        world.set_g("tectonics.last_p50_planform_recycle_candidate_cells", 0.0)
        world.set_g("tectonics.last_p50_planform_fill_candidate_cells", 0.0)
        world.set_g("tectonics.last_p50_planform_recycled_cells", 0.0)
        world.set_g("tectonics.last_p50_planform_filled_cells", 0.0)
        world.set_g("tectonics.last_p50_planform_area_fraction", 0.0)
        world.set_g("tectonics.last_p50_planform_before_narrow_fraction", 0.0)
        world.set_g("tectonics.last_p50_planform_after_narrow_fraction", 0.0)
        world.set_g("tectonics.last_p50_planform_narrow_fraction_delta", 0.0)
        world.set_g("tectonics.last_p57_parented_basement_cargo_inherited_cells", 0.0)
        world.set_g("tectonics.last_p57_parented_old_basement_cargo_inherited_cells", 0.0)
        world.set_g("tectonics.last_p57_erosion_basement_cargo_cleared_cells", 0.0)
        basement_enabled = (
            isinstance(basement_age, np.ndarray)
            and isinstance(basement_stability, np.ndarray)
            and isinstance(basement_thick, np.ndarray)
            and basement_age.shape == (grid.n,)
            and basement_stability.shape == (grid.n,)
            and basement_thick.shape == (grid.n,)
        )
        shape_pressure = self._continental_shape_pressure(grid, ctype, origin, stability)
        world.set_g("tectonics.last_continent_shape_pressure", shape_pressure)
        shape_budget = total * min(
            0.00075,
            (0.00028 + 0.00045 * shape_pressure) * max(dt, 1.0) / 20.0,
        )
        accretion_zone = None
        erosion_zone = None
        if has_accretion_process:
            accretion_zone = self._dilate_mask(
                grid, accretion_seed, allowed=np.ones(grid.n, dtype=bool), passes=2,
            )
        if has_erosion_process:
            erosion_zone = self._dilate_mask(
                grid, erosion_seed, allowed=np.ones(grid.n, dtype=bool), passes=2,
            )
        restoration_zone = self._passive_margin_progradation_zone(
            grid, ctype, origin, stability, erosion_zone, shape_pressure)
        p34_restoration_zone = self._parented_platform_progradation_zone(
            world, grid, ctype, age, origin, reworked, stability, orog_age, volc_age,
            erosion_zone, shape_pressure, target, frac, t)
        p35_restoration_zone = self._continental_debt_recovery_zone(
            world, grid, ctype, thick, age, origin, reworked, stability, erosion_zone,
            shape_pressure, target, frac, t)
        p44_restoration_zone = self._mature_parented_recovery_zone(
            world, grid, ctype, thick, age, origin, reworked, stability, orog_age, volc_age,
            erosion_zone, shape_pressure, target, frac, t)
        reference_earthlike = (
            str(world.spec.name).startswith("earthlike")
            and 0.25 <= float(world.spec.target_land_fraction) <= 0.32
            and 0.90 <= float(world.spec.composition.water_inventory_earth) <= 1.10
            and 11 <= int(world.spec.n_plates) <= 13
        )
        mature_stage = float(t) >= max(3000.0, 0.98 * float(world.spec.t_end_myr))
        earthlike_floor = 0.25
        below_earthlike_floor = bool(
            reference_earthlike and mature_stage and float(frac) < earthlike_floor
        )
        p104e_floor_restoration_zone = np.zeros(grid.n, dtype=bool)
        if below_earthlike_floor:
            cont = ctype == CONT
            width = _graph_width_steps(grid, cont)
            stable_parent = (
                cont
                & (
                    (origin == ORIGIN_CRATON)
                    | (stability >= 0.54)
                    | (width >= 2.0)
                )
            )
            cont_neighbors = _same_neighbor_count(grid, cont)
            stable_neighbors = _same_neighbor_count(grid, stable_parent)
            erosion_block = (
                np.zeros(grid.n, dtype=bool)
                if erosion_zone is None
                else np.asarray(erosion_zone, dtype=bool)
            )
            p104e_floor_restoration_zone = (
                ~cont
                & (cont_neighbors >= 2)
                & (stable_neighbors >= 1)
                & ~erosion_block
            )
            if not p104e_floor_restoration_zone.any():
                p104e_floor_restoration_zone = (
                    ~cont
                    & (cont_neighbors >= 1)
                    & ~erosion_block
                )
            world.set_g(
                "tectonics.last_p104e_earthlike_floor_target_fraction",
                earthlike_floor,
            )
            world.set_g(
                "tectonics.last_p104e_earthlike_floor_gap_fraction",
                max(earthlike_floor - float(frac), 0.0),
            )
            world.set_g(
                "tectonics.last_p104e_earthlike_floor_candidate_cells",
                float(np.count_nonzero(p104e_floor_restoration_zone)),
            )
        restoration_zone |= p34_restoration_zone
        restoration_zone |= p35_restoration_zone
        restoration_zone |= p44_restoration_zone
        restoration_zone |= p104e_floor_restoration_zone
        effective_target = (
            max(float(target), earthlike_floor)
            if below_earthlike_floor else float(target)
        )
        has_restoration_process = bool(
            restoration_zone.any()
            and (float(frac) < float(target) - 0.02 or below_earthlike_floor)
        )
        prefer_parented_restoration = bool(
            has_restoration_process
            and world.g("tectonics.enable_parented_restoration_priority", 0.0) > 0.0
        )
        if accretion_zone is None:
            accretion_zone = np.zeros(grid.n, dtype=bool)
        accretion_or_restoration_zone = accretion_zone | restoration_zone
        if abs(frac - effective_target) <= 0.02 and not below_earthlike_floor:
            if has_accretion_process or has_erosion_process or shape_pressure > 0.0:
                world.set_g(
                    "tectonics.last_background_continent_shape_maintenance",
                    float(not (has_accretion_process or has_erosion_process)),
                )
                ctype, thick, age, origin, reworked, stability, orog_age, volc_age = (
                    self._shape_aware_continental_maintenance(
                        grid, ctype, thick, age, origin, reworked, stability,
                        orog_age, volc_age, t, shape_budget,
                        add_allowed=accretion_or_restoration_zone,
                        active_accretion_allowed=accretion_zone,
                        erase_allowed=erosion_zone,
                        world=world)
                )
                frac = area[ctype == CONT].sum() / total
            return ctype, thick, age, origin, reworked, stability, orog_age, volc_age, float(frac)
        neigh = grid.neighbors
        is_cont = ctype == CONT
        budget = total * min(
            abs(float(frac) - effective_target),
            0.04 * max(dt, 1.0) / 20.0,
        )
        if budget <= 0.0:
            return ctype, thick, age, origin, reworked, stability, orog_age, volc_age, float(frac)
        if erosion_zone is None:
            erosion_zone = np.zeros(grid.n, dtype=bool)
        # border cells of each type
        if frac > effective_target:             # too much continent: erode borders
            if not has_erosion_process:
                world.set_g("tectonics.last_unforced_continent_loss_blocked", 1.0)
                return ctype, thick, age, origin, reworked, stability, orog_age, volc_age, float(frac)
            width = _graph_width_steps(grid, is_cont)
            cont_neighbors = _same_neighbor_count(grid, is_cont)
            protected_core = (origin == ORIGIN_CRATON) | (stability > 0.78)
            candidates = [int(c) for c in np.where(is_cont & ~protected_core)[0]
                          if erosion_zone[c] and np.any(~is_cont[neigh[c]])]
            if not candidates:
                candidates = [int(c) for c in np.where(is_cont & erosion_zone)[0]
                              if np.any(~is_cont[neigh[c]])]
            candidates.sort(key=lambda c: (
                0 if erosion_zone[c] else 1,
                0 if width[c] <= 2.0 else 1,
                0 if origin[c] in (ORIGIN_ARC, ORIGIN_SUTURE) else 1,
                float(width[c]),
                -int(cont_neighbors[c] <= 2),
                float(stability[c]),
                float(thick[c]),
                -float(age[c]),
                c,
            ))
            need = budget
            lost = 0
            p57_eroded_basement = 0
            for c in candidates:
                if need <= 0:
                    break
                if basement_enabled and (
                    basement_age[c] > 0.0
                    or basement_stability[c] > 0.0
                    or basement_thick[c] > 0.0
                ):
                    p57_eroded_basement += 1
                    basement_age[c] = 0.0
                    basement_stability[c] = 0.0
                    basement_thick[c] = 0.0
                ctype[c] = OCEAN
                thick[c] = OCEAN_THICK
                age[c] = min(age[c], 80.0)
                origin[c] = ORIGIN_RIDGE
                reworked[c] = t
                stability[c] = 0.0
                orog_age[c] = -1.0
                volc_age[c] = -1.0
                need -= area[c]
                lost += 1
            world.set_g("tectonics.last_continent_loss_cells", float(lost))
            world.set_g(
                "tectonics.last_p57_erosion_basement_cargo_cleared_cells",
                float(p57_eroded_basement),
            )
        else:                                   # too little: accrete at borders
            if not (has_accretion_process or has_restoration_process):
                world.set_g("tectonics.last_unforced_continent_gain_blocked", 1.0)
                return ctype, thick, age, origin, reworked, stability, orog_age, volc_age, float(frac)
            cont_neighbors = _same_neighbor_count(grid, is_cont)
            stable_cont = is_cont & ((origin == ORIGIN_CRATON) | (stability > 0.70))
            stable_neighbors = _same_neighbor_count(grid, stable_cont)
            mature_stage = float(t) >= max(3000.0, 0.98 * float(world.spec.t_end_myr))
            mature_active_guard = bool(
                mature_stage
                and world.g("tectonics.enable_mature_active_accretion_shape_guard", 0.0) > 0.0
            )
            broad_neighbors = np.zeros(grid.n, dtype=np.int64)
            active_only_cap_area = float(total)
            if mature_active_guard:
                width = _graph_width_steps(grid, is_cont)
                broad_cont = is_cont & (width >= 3.0)
                broad_neighbors = _same_neighbor_count(grid, broad_cont)
                gap = max(float(target) - float(frac), 0.0)
                cap_fraction = min(0.018, max(0.004, 0.20 * gap))
                active_only_cap_area = (
                    float(total) * cap_fraction * max(float(dt), 1.0) / 20.0
                )
                world.set_g(
                    "tectonics.last_p43_mature_active_only_cap_area_fraction",
                    float(active_only_cap_area / max(float(total), 1.0)),
                )
            candidates = [int(c) for c in np.where(~is_cont)[0]
                          if accretion_or_restoration_zone[c] and np.any(is_cont[neigh[c]])]
            candidates = [c for c in candidates if cont_neighbors[c] >= 1]
            if prefer_parented_restoration:
                candidates.sort(key=lambda c: (
                    0 if p104e_floor_restoration_zone[c] else 1,
                    0 if p44_restoration_zone[c] else 1,
                    0 if restoration_zone[c] else 1,
                    0 if cont_neighbors[c] >= 3 else 1,
                    -int(cont_neighbors[c]),
                    -int(stable_neighbors[c]),
                    0 if accretion_zone[c] else 1,
                    0 if volc_age[c] >= max(t - 250.0, 0.0) else 1,
                    -float(thick[c]),
                    -float(age[c]),
                    c,
                ))
            else:
                candidates.sort(key=lambda c: (
                    0 if p104e_floor_restoration_zone[c] else 1,
                    0 if p44_restoration_zone[c] else 1,
                    0 if accretion_zone[c] else 1,
                    0 if restoration_zone[c] else 1,
                    0 if cont_neighbors[c] >= 3 else 1,
                    -int(cont_neighbors[c]),
                    -int(stable_neighbors[c]),
                    0 if volc_age[c] >= max(t - 250.0, 0.0) else 1,
                    -float(thick[c]),
                    -float(age[c]),
                    c,
                ))
            need = budget
            gained = 0
            restored = 0
            p34_restored = 0
            p34_restored_area = 0.0
            p35_restored = 0
            p35_restored_area = 0.0
            p40_inherited_rework = 0
            p40_reactivated = 0
            p41_priority_restored = 0
            p41_deferred_active_only = 0
            p43_mature_active_candidates = 0
            p43_mature_active_accepted = 0
            p43_mature_active_rejected = 0
            p43_mature_active_area = 0.0
            p104e_floor_restored = 0
            p104e_floor_area = 0.0
            p44_parented_restored = 0
            p44_parented_area = 0.0
            p57_parented_basement = 0
            p57_parented_old_basement = 0
            inherit_rework_parent = self._quiet_mature_rework_inheritance_parent_mask(
                grid, ctype, thick, age, origin, reworked, stability, t)
            for c in candidates:
                if need <= 0:
                    break
                parent = None
                if restoration_zone[c] and (
                    p104e_floor_restoration_zone[c]
                    or p44_restoration_zone[c]
                    or prefer_parented_restoration
                    or not accretion_zone[c]
                ):
                    parent = self._continental_margin_parent(
                        grid, c, ctype, thick, age, origin, stability,
                        preferred=inherit_rework_parent)
                active_only_candidate = bool(accretion_zone[c] and parent is None)
                if prefer_parented_restoration and active_only_candidate:
                    p41_deferred_active_only += 1
                    continue
                if mature_active_guard and active_only_candidate:
                    p43_mature_active_candidates += 1
                    parent_supported = (
                        cont_neighbors[c] >= 2
                        or broad_neighbors[c] >= 1
                    )
                    if (
                        not parent_supported
                        or p43_mature_active_area + float(area[c]) > active_only_cap_area
                    ):
                        p43_mature_active_rejected += 1
                        continue
                ctype[c] = CONT
                if parent is not None:
                    p = int(parent)
                    parent_origin = origin[p]
                    origin[c] = (
                        ORIGIN_PRIMORDIAL if parent_origin == ORIGIN_CRATON
                        else parent_origin
                    )
                    thick[c] = max(thick[c], CONT_THICK, min(thick[p], 54000.0))
                    age[c] = max(age[c], min(age[p], t))
                    stability[c] = max(stability[c], min(stability[p] * 0.90, 0.74))
                    if orog_age[p] >= 0.0:
                        orog_age[c] = orog_age[p]
                    if parent_origin == ORIGIN_ARC and volc_age[p] >= 0.0:
                        volc_age[c] = volc_age[p]
                    if basement_enabled and basement_thick[p] > 0.0:
                        basement_age[c] = max(basement_age[c], min(basement_age[p], t))
                        basement_stability[c] = max(
                            basement_stability[c],
                            min(basement_stability[p] * 0.96, 0.88),
                        )
                        basement_thick[c] = max(
                            basement_thick[c],
                            min(basement_thick[p], thick[c]),
                        )
                        p57_parented_basement += 1
                        if basement_age[c] >= 2500.0 and basement_stability[c] >= 0.70:
                            p57_parented_old_basement += 1
                    if inherit_rework_parent[p]:
                        reworked[c] = reworked[p]
                        p40_inherited_rework += 1
                    else:
                        reworked[c] = t
                        p40_reactivated += 1
                    if prefer_parented_restoration and restoration_zone[c]:
                        p41_priority_restored += 1
                    if p44_restoration_zone[c]:
                        p44_parented_restored += 1
                        p44_parented_area += float(area[c])
                    if p104e_floor_restoration_zone[c]:
                        p104e_floor_restored += 1
                        p104e_floor_area += float(area[c])
                    restored += 1
                    if p34_restoration_zone[c]:
                        p34_restored += 1
                        p34_restored_area += float(area[c])
                    if p35_restoration_zone[c]:
                        p35_restored += 1
                        p35_restored_area += float(area[c])
                else:
                    thick[c] = max(thick[c], CONT_THICK)
                    age[c] = max(age[c], min(800.0, t))
                    origin[c] = ORIGIN_ARC
                    stability[c] = max(stability[c], 0.2)
                    reworked[c] = t
                    p40_reactivated += 1
                    if mature_active_guard and active_only_candidate:
                        p43_mature_active_accepted += 1
                        p43_mature_active_area += float(area[c])
                    if p104e_floor_restoration_zone[c]:
                        p104e_floor_restored += 1
                        p104e_floor_area += float(area[c])
                need -= area[c]
                gained += 1
            world.set_g("tectonics.last_continent_gain_cells", float(gained))
            world.set_g("tectonics.last_passive_margin_progradation_cells", float(restored))
            world.set_g(
                "tectonics.last_p40_parented_restoration_inherited_rework_cells",
                float(p40_inherited_rework),
            )
            world.set_g(
                "tectonics.last_p40_parented_restoration_reactivated_cells",
                float(p40_reactivated),
            )
            world.set_g(
                "tectonics.last_p41_parented_restoration_priority_cells",
                float(p41_priority_restored),
            )
            world.set_g(
                "tectonics.last_p41_active_only_accretion_deferred_cells",
                float(p41_deferred_active_only),
            )
            world.set_g(
                "tectonics.last_p43_mature_active_only_candidate_cells",
                float(p43_mature_active_candidates),
            )
            world.set_g(
                "tectonics.last_p43_mature_active_only_accepted_cells",
                float(p43_mature_active_accepted),
            )
            world.set_g(
                "tectonics.last_p43_mature_active_only_rejected_cells",
                float(p43_mature_active_rejected),
            )
            world.set_g(
                "tectonics.last_p43_mature_active_only_accepted_area_fraction",
                float(p43_mature_active_area / max(float(total), 1.0)),
            )
            world.set_g(
                "tectonics.last_p104e_earthlike_floor_restored_cells",
                float(p104e_floor_restored),
            )
            world.set_g(
                "tectonics.last_p104e_earthlike_floor_restored_area_fraction",
                float(p104e_floor_area / max(float(total), 1.0)),
            )
            world.set_g(
                "tectonics.last_p44_parented_recovery_restored_cells",
                float(p44_parented_restored),
            )
            world.set_g(
                "tectonics.last_p44_parented_recovery_area_fraction",
                float(p44_parented_area / max(float(total), 1.0)),
            )
            world.set_g(
                "tectonics.last_p57_parented_basement_cargo_inherited_cells",
                float(p57_parented_basement),
            )
            world.set_g(
                "tectonics.last_p57_parented_old_basement_cargo_inherited_cells",
                float(p57_parented_old_basement),
            )
            world.set_g("tectonics.last_p34_parented_platform_restored_cells", float(p34_restored))
            world.set_g(
                "tectonics.last_p34_parented_platform_area_fraction",
                float(p34_restored_area / max(float(total), 1.0)),
            )
            world.set_g("tectonics.last_p35_continental_debt_restored_cells", float(p35_restored))
            world.set_g(
                "tectonics.last_p35_continental_debt_area_fraction",
                float(p35_restored_area / max(float(total), 1.0)),
            )
            if p35_restored:
                world.set_g(
                    "tectonics.cumulative_p35_continental_debt_restored_cells",
                    world.g("tectonics.cumulative_p35_continental_debt_restored_cells", 0.0)
                    + float(p35_restored),
                )
                world.set_g(
                    "tectonics.cumulative_p35_continental_debt_area_fraction",
                    world.g("tectonics.cumulative_p35_continental_debt_area_fraction", 0.0)
                    + float(p35_restored_area / max(float(total), 1.0)),
                )
        ctype, thick, age, origin, reworked, stability, orog_age, volc_age = (
            self._shape_aware_continental_maintenance(
                grid, ctype, thick, age, origin, reworked, stability,
                orog_age, volc_age, t, shape_budget,
                add_allowed=accretion_or_restoration_zone,
                active_accretion_allowed=accretion_zone,
                erase_allowed=erosion_zone,
                world=world)
        )
        frac = area[ctype == CONT].sum() / total
        return ctype, thick, age, origin, reworked, stability, orog_age, volc_age, float(frac)

    def _stable_parent_anchor_mask(self, grid, ctype, origin, stability,
                                   age=None, thick=None, reworked=None, t=None):
        cont = ctype == CONT
        if not cont.any():
            return np.zeros(grid.n, dtype=bool)
        width = _graph_width_steps(grid, cont)
        anchor = (
            (origin == ORIGIN_CRATON)
            | ((origin == ORIGIN_PRIMORDIAL) & (stability >= 0.42))
            | ((width >= 3.0) & (stability >= 0.54))
            | (stability >= 0.72)
        )
        if age is not None and thick is not None:
            funnel = self._old_platform_anchor_funnel(
                None, grid, ctype, thick, age, origin, reworked, stability, t)
            masks = funnel["masks"]
            anchor |= (
                masks["old_platform_core"]
                | masks["old_platform_margin"]
                | masks["mature_platform"]
            )
        return cont & anchor

    def _old_platform_anchor_funnel(self, world, grid, ctype, thick, age, origin,
                                    reworked, stability, t):
        """Diagnose each gate used to recognize old quiet platform anchors."""
        cont = np.asarray(ctype, dtype=np.float64) == CONT
        zero = np.zeros(grid.n, dtype=bool)
        masks = {
            "cont": cont,
            "legacy_anchor": zero.copy(),
            "source_ok": zero.copy(),
            "thick_ok": zero.copy(),
            "quiet": zero.copy(),
            "age_mature": zero.copy(),
            "age_old": zero.copy(),
            "width_ge3": zero.copy(),
            "width_ge4": zero.copy(),
            "stability_ge030": zero.copy(),
            "stability_ge036": zero.copy(),
            "old_platform_core": zero.copy(),
            "old_platform_margin": zero.copy(),
            "mature_platform": zero.copy(),
            "process_anchor": zero.copy(),
            "lineage_broad": zero.copy(),
        }
        if not cont.any():
            metrics = {key: 0.0 for key in P39_OLD_PLATFORM_FUNNEL_METRICS}
            return {"metrics": metrics, "masks": masks}

        origin_arr = np.asarray(origin, dtype=np.float64)
        stability_arr = np.asarray(stability, dtype=np.float64)
        thick_arr = np.asarray(thick, dtype=np.float64)
        age_arr = np.asarray(age, dtype=np.float64)
        if (
            origin_arr.shape != (grid.n,)
            or stability_arr.shape != (grid.n,)
            or thick_arr.shape != (grid.n,)
            or age_arr.shape != (grid.n,)
        ):
            metrics = {key: 0.0 for key in P39_OLD_PLATFORM_FUNNEL_METRICS}
            metrics["continental_cells"] = float(np.count_nonzero(cont))
            return {"metrics": metrics, "masks": masks}

        width = _graph_width_steps(grid, cont)
        legacy_anchor = cont & (
            (origin_arr == ORIGIN_CRATON)
            | ((origin_arr == ORIGIN_PRIMORDIAL) & (stability_arr >= 0.42))
            | ((width >= 3.0) & (stability_arr >= 0.54))
            | (stability_arr >= 0.72)
        )

        if t is None:
            old_cut = 1600.0
            mature_cut = 1200.0
            quiet = np.ones(grid.n, dtype=bool)
        else:
            tt = float(t)
            old_cut = min(2200.0, max(900.0, 0.42 * tt))
            mature_cut = min(1800.0, max(700.0, 0.30 * tt))
            quiet_myr = min(900.0, max(360.0, 0.18 * tt))
            quiet = np.ones(grid.n, dtype=bool)
            if reworked is not None:
                rw = np.asarray(reworked, dtype=np.float64)
                if rw.shape == (grid.n,):
                    quiet = (rw < 0.0) | (rw <= tt - quiet_myr)

        source_ok = np.isin(
            origin_arr, [ORIGIN_PRIMORDIAL, ORIGIN_CRATON, ORIGIN_SUTURE])
        thick_ok = thick_arr >= CONT_THICK + 500.0
        age_mature = age_arr >= mature_cut
        age_old = age_arr >= old_cut
        width_ge3 = width >= 3.0
        width_ge4 = width >= 4.0
        stability_ge030 = stability_arr >= 0.30
        stability_ge036 = stability_arr >= 0.36

        old_platform_core = (
            cont
            & width_ge3
            & source_ok
            & quiet
            & thick_ok
            & age_old
            & stability_ge030
        )
        old_platform_margin = zero.copy()
        if old_platform_core.any():
            core_support = self._dilate_mask(
                grid, old_platform_core, allowed=cont, passes=2)
            old_platform_margin = (
                cont
                & source_ok
                & quiet
                & thick_ok
                & age_mature
                & stability_ge030
                & core_support
            )
        mature_platform = (
            cont
            & width_ge4
            & (origin_arr == ORIGIN_PRIMORDIAL)
            & quiet
            & age_mature
            & stability_ge036
        )
        process_anchor = cont & (
            legacy_anchor | old_platform_core | old_platform_margin | mature_platform
        )

        lineage_broad = zero.copy()
        if world is not None:
            previous = world.get_field("tectonics.continent_id", -1.0)
            previous = np.asarray(previous, dtype=np.float64)
            if previous.shape == (grid.n,):
                broad_ids = {
                    int(obj.get("numeric_id", -1))
                    for obj in world.objects.get("tectonics.continents", [])
                    if (
                        int(obj.get("numeric_id", -1)) >= 0
                        and float(obj.get("width_p50_steps", 0.0)) >= 2.75
                        and float(obj.get("width_p90_steps", 0.0)) >= 4.0
                    )
                }
                if broad_ids:
                    prev_ids = previous.astype(int)
                    lineage_broad = cont & np.isin(
                        prev_ids, np.asarray(sorted(broad_ids), dtype=int))

        masks = {
            "cont": cont,
            "legacy_anchor": legacy_anchor,
            "source_ok": cont & source_ok,
            "thick_ok": cont & thick_ok,
            "quiet": cont & quiet,
            "age_mature": cont & age_mature,
            "age_old": cont & age_old,
            "width_ge3": cont & width_ge3,
            "width_ge4": cont & width_ge4,
            "stability_ge030": cont & stability_ge030,
            "stability_ge036": cont & stability_ge036,
            "old_platform_core": old_platform_core,
            "old_platform_margin": old_platform_margin,
            "mature_platform": mature_platform,
            "process_anchor": process_anchor,
            "lineage_broad": lineage_broad,
        }
        metrics = {
            "continental_cells": float(np.count_nonzero(cont)),
            "legacy_anchor_cells": float(np.count_nonzero(legacy_anchor)),
            "source_ok_cells": float(np.count_nonzero(cont & source_ok)),
            "source_blocked_cells": float(np.count_nonzero(cont & ~source_ok)),
            "thick_ok_cells": float(np.count_nonzero(cont & thick_ok)),
            "thick_blocked_cells": float(np.count_nonzero(cont & ~thick_ok)),
            "quiet_cells": float(np.count_nonzero(cont & quiet)),
            "recently_reworked_cells": float(np.count_nonzero(cont & ~quiet)),
            "age_mature_cells": float(np.count_nonzero(cont & age_mature)),
            "age_old_cells": float(np.count_nonzero(cont & age_old)),
            "age_mature_blocked_cells": float(np.count_nonzero(cont & ~age_mature)),
            "age_old_blocked_cells": float(np.count_nonzero(cont & ~age_old)),
            "width_ge3_cells": float(np.count_nonzero(cont & width_ge3)),
            "width_ge4_cells": float(np.count_nonzero(cont & width_ge4)),
            "width_ge3_blocked_cells": float(np.count_nonzero(cont & ~width_ge3)),
            "width_ge4_blocked_cells": float(np.count_nonzero(cont & ~width_ge4)),
            "stability_ge030_cells": float(np.count_nonzero(cont & stability_ge030)),
            "stability_ge036_cells": float(np.count_nonzero(cont & stability_ge036)),
            "stability_ge030_blocked_cells": float(
                np.count_nonzero(cont & ~stability_ge030)),
            "stability_ge036_blocked_cells": float(
                np.count_nonzero(cont & ~stability_ge036)),
            "old_platform_core_cells": float(np.count_nonzero(old_platform_core)),
            "old_platform_margin_cells": float(np.count_nonzero(old_platform_margin)),
            "mature_platform_cells": float(np.count_nonzero(mature_platform)),
            "process_anchor_cells": float(np.count_nonzero(process_anchor)),
            "lineage_broad_cells": float(np.count_nonzero(lineage_broad)),
            "process_anchor_lineaged_cells": float(
                np.count_nonzero(process_anchor & lineage_broad)),
            "process_anchor_unlineaged_cells": float(
                np.count_nonzero(process_anchor & ~lineage_broad)),
        }
        return {
            "metrics": metrics,
            "masks": masks,
            "old_cut_myr": float(old_cut),
            "mature_cut_myr": float(mature_cut),
        }

    def _record_old_platform_anchor_funnel(self, world, funnel):
        for key in P39_OLD_PLATFORM_FUNNEL_METRICS:
            world.set_g(
                f"tectonics.last_p39_old_platform_funnel_{key}",
                float(funnel["metrics"].get(key, 0.0)),
            )

    def _quiet_mature_rework_inheritance_parent_mask(
        self, grid, ctype, thick, age, origin, reworked, stability, t
    ):
        """Parents whose quiet history can be inherited by maintenance fills.

        This targets the P36 gap: quiet mature platforms that are not already
        accepted by the legacy stable-anchor rule.  True cratons also inherit.
        Ordinary legacy-stable margins still count as reactivated when they
        gain new cells, avoiding widespread ancient fragments from small
        conservation repairs.
        """
        cont = np.asarray(ctype, dtype=np.float64) == CONT
        out = np.zeros(grid.n, dtype=bool)
        if not cont.any():
            return out
        thick_arr = np.asarray(thick, dtype=np.float64)
        age_arr = np.asarray(age, dtype=np.float64)
        origin_arr = np.asarray(origin, dtype=np.float64)
        reworked_arr = np.asarray(reworked, dtype=np.float64)
        stability_arr = np.asarray(stability, dtype=np.float64)
        if (
            thick_arr.shape != (grid.n,)
            or age_arr.shape != (grid.n,)
            or origin_arr.shape != (grid.n,)
            or reworked_arr.shape != (grid.n,)
            or stability_arr.shape != (grid.n,)
        ):
            return out
        del thick_arr, age_arr, reworked_arr
        legacy_anchor = self._stable_parent_anchor_mask(
            grid, ctype, origin, stability)
        funnel = self._old_platform_anchor_funnel(
            None, grid, ctype, thick, age, origin, reworked, stability, t)
        masks = funnel["masks"]
        quiet_mature = (
            masks["old_platform_core"]
            | masks["old_platform_margin"]
            | masks["mature_platform"]
        )
        return quiet_mature & (
            (origin_arr == ORIGIN_CRATON) | ~legacy_anchor
        )

    def _mature_parented_recovery_zone(self, world, grid, ctype, thick, age,
                                       origin, reworked, stability, orog_age,
                                       volc_age, erosion_zone, shape_pressure,
                                       target, frac, t):
        """Parent-first mature platform recovery before active-only accretion."""
        out = np.zeros(grid.n, dtype=bool)
        if world.g("tectonics.enable_mature_parented_recovery", 0.0) <= 0.0:
            return out
        mature_start_myr = max(3000.0, 0.98 * float(world.spec.t_end_myr))
        if float(t) < mature_start_myr:
            return out
        gap = float(target) - float(frac)
        del shape_pressure
        if gap < 0.025:
            return out

        cont = np.asarray(ctype, dtype=np.float64) == CONT
        if not cont.any():
            return out
        width = _graph_width_steps(grid, cont)
        funnel = self._old_platform_anchor_funnel(
            world, grid, ctype, thick, age, origin, reworked, stability, t)
        self._record_old_platform_anchor_funnel(world, funnel)
        masks = funnel["masks"]
        legacy_parent = masks["legacy_anchor"]
        process_parent = masks["process_anchor"]
        if not bool(world.g("tectonics.enable_old_platform_debt_recovery", 0.0)):
            process_parent = legacy_parent

        origin_arr = np.asarray(origin, dtype=np.float64)
        stability_arr = np.asarray(stability, dtype=np.float64)
        age_arr = np.asarray(age, dtype=np.float64)
        thick_arr = np.asarray(thick, dtype=np.float64)
        broad_mature_parent = (
            cont
            & (width >= 3.0)
            & np.isin(origin_arr, [ORIGIN_PRIMORDIAL, ORIGIN_CRATON, ORIGIN_SUTURE])
            & (stability_arr >= 0.30)
            & (age_arr >= min(1400.0, max(760.0, 0.22 * float(t))))
            & (thick_arr >= CONT_THICK)
        )
        broad_parent_margin = self._dilate_mask(
            grid, broad_mature_parent, allowed=cont, passes=3)
        parent = process_parent | broad_parent_margin
        if not parent.any():
            return out

        cont_neighbors = _same_neighbor_count(grid, cont)
        parent_neighbors = _same_neighbor_count(grid, parent)
        labels, _ = _component_labels(grid, cont)
        single_parent_component = np.zeros(grid.n, dtype=bool)
        for c in np.where(~cont & (cont_neighbors >= 1) & (parent_neighbors >= 1))[0]:
            parent_ids = {
                int(labels[int(nb)])
                for nb in grid.neighbors[int(c)]
                if cont[int(nb)] and labels[int(nb)] >= 0
            }
            single_parent_component[int(c)] = len(parent_ids) <= 1

        quiet_myr = min(900.0, max(360.0, 0.18 * float(t)))
        quiet_rework = (reworked < 0.0) | (reworked <= float(t) - quiet_myr)
        quiet_orogen = (orog_age < 0.0) | (orog_age <= float(t) - quiet_myr)
        quiet_volcanic = (volc_age < 0.0) | (volc_age <= float(t) - quiet_myr)
        quiet = quiet_rework & quiet_orogen & quiet_volcanic

        apron = self._dilate_mask(grid, parent, allowed=~cont, passes=2)
        embayment = (
            (~cont)
            & apron
            & single_parent_component
            & (
                (cont_neighbors >= 3)
                | ((cont_neighbors >= 2) & (parent_neighbors >= 1))
                | ((cont_neighbors >= 1) & (parent_neighbors >= 2))
            )
        )
        old_arc_or_shell = (
            embayment
            & (
                quiet
                | ((origin_arr == ORIGIN_ARC) & (age_arr >= 120.0))
                | (age_arr >= 180.0)
            )
        )
        protected = self._protected_breakup_corridor_mask(world, grid, passes=1)
        if erosion_zone is not None:
            erosion = np.asarray(erosion_zone, dtype=bool)
            if erosion.shape == (grid.n,) and erosion.any():
                protected |= self._dilate_mask(
                    grid, erosion, allowed=np.ones(grid.n, dtype=bool), passes=1)
        world.set_g(
            "tectonics.last_p44_parented_recovery_protected_cells",
            float(np.count_nonzero(old_arc_or_shell & protected)),
        )
        eligible = old_arc_or_shell & ~protected
        if not eligible.any():
            return out

        max_area = float(grid.cell_area.sum()) * min(
            0.018,
            max(0.0040, 0.0040 + 0.060 * max(gap, 0.0)),
        )
        cells = sorted((int(c) for c in np.where(eligible)[0]), key=lambda c: (
            0 if parent_neighbors[c] >= 2 else 1,
            0 if cont_neighbors[c] >= 3 else 1,
            -int(parent_neighbors[c]),
            -int(cont_neighbors[c]),
            -float(age_arr[c]),
            c,
        ))
        acc_area = 0.0
        for c in cells:
            if acc_area >= max_area:
                break
            out[c] = True
            acc_area += float(grid.cell_area[c])
        world.set_g(
            "tectonics.last_p44_parented_recovery_candidate_cells",
            float(np.count_nonzero(out)),
        )
        return out

    def _continental_debt_recovery_zone(self, world, grid, ctype, thick, age,
                                        origin, reworked, stability, erosion_zone,
                                        shape_pressure, target, frac, t):
        """Recover parented margin cells when fixed-grid advection lost continent.

        This is a conservation-debt repair, not a new tectonic accretion event:
        it only fills same-parent edge embayments around stable continental
        anchors when the mature Earth-like continental budget is far below its
        initial target.  Active breakup/ridge corridors and explicit seaways are
        excluded.
        """
        out = np.zeros(grid.n, dtype=bool)
        mature_start_myr = max(3000.0, 0.98 * float(world.spec.t_end_myr))
        if float(t) < mature_start_myr:
            return out
        gap = float(target) - float(frac)
        world.set_g("tectonics.last_p35_continental_gap_fraction", max(gap, 0.0))
        if gap < 0.085 or (shape_pressure < 0.55 and gap < 0.14):
            return out
        cont = ctype == CONT
        if not cont.any():
            return out
        funnel = self._old_platform_anchor_funnel(
            world, grid, ctype, thick, age, origin, reworked, stability, t)
        self._record_old_platform_anchor_funnel(world, funnel)
        legacy_parent = funnel["masks"]["legacy_anchor"]
        process_parent = funnel["masks"]["process_anchor"]
        allow_old_platform_anchors = bool(
            world.g("tectonics.enable_old_platform_debt_recovery", 0.0)
        )
        stable_parent = process_parent if allow_old_platform_anchors else legacy_parent
        world.set_g(
            "tectonics.last_p35_stable_parent_cells",
            float(np.count_nonzero(stable_parent)),
        )
        world.set_g(
            "tectonics.last_p35_old_platform_anchor_cells",
            float(np.count_nonzero(stable_parent & ~legacy_parent)),
        )
        if not stable_parent.any():
            return out

        cont_neighbors = _same_neighbor_count(grid, cont)
        stable_neighbors = _same_neighbor_count(grid, stable_parent)
        labels, _ = _component_labels(grid, cont)
        single_parent_component = np.zeros(grid.n, dtype=bool)
        for c in np.where(~cont & (cont_neighbors >= 2) & (stable_neighbors >= 1))[0]:
            parent_ids = {
                int(labels[int(nb)])
                for nb in grid.neighbors[int(c)]
                if cont[int(nb)] and labels[int(nb)] >= 0
            }
            single_parent_component[int(c)] = len(parent_ids) <= 1

        apron = self._dilate_mask(grid, stable_parent, allowed=~cont, passes=1)
        embayment = (
            (~cont)
            & apron
            & single_parent_component
            & (
                (cont_neighbors >= 3)
                | ((cont_neighbors >= 2) & (stable_neighbors >= 2))
            )
        )
        protected = self._protected_breakup_corridor_mask(world, grid, passes=1)
        if erosion_zone is not None:
            erosion = np.asarray(erosion_zone, dtype=bool)
            if erosion.shape == (grid.n,) and erosion.any():
                protected |= self._dilate_mask(
                    grid, erosion, allowed=np.ones(grid.n, dtype=bool), passes=1)
        eligible = embayment & ~protected
        if eligible.any():
            max_area = float(grid.cell_area.sum()) * min(0.0025, 0.0010 + 0.010 * gap)
            cells = sorted((int(c) for c in np.where(eligible)[0]), key=lambda c: (
                -int(cont_neighbors[c]),
                -int(stable_neighbors[c]),
                c,
            ))
            acc_area = 0.0
            for c in cells:
                if acc_area >= max_area:
                    break
                out[c] = True
                acc_area += float(grid.cell_area[c])
        world.set_g(
            "tectonics.last_p35_continental_debt_candidate_cells",
            float(np.count_nonzero(out)),
        )
        raw_candidate_count = float(np.count_nonzero(out))
        if raw_candidate_count > 0.0:
            world.set_g(
                "tectonics.cumulative_p35_continental_debt_candidate_cells",
                world.g("tectonics.cumulative_p35_continental_debt_candidate_cells", 0.0)
                + raw_candidate_count,
            )
        if bool(world.g("tectonics.enable_shape_guarded_debt_recovery", 0.0)):
            guarded = self._shape_guard_continental_restoration(
                world, grid, cont, out)
            guard_input = float(np.count_nonzero(out))
            guard_accepted = float(np.count_nonzero(guarded))
            guard_rejected = float(np.count_nonzero(out & ~guarded))
            world.set_g(
                "tectonics.last_p37_shape_guard_input_cells",
                guard_input,
            )
            world.set_g(
                "tectonics.last_p37_shape_guard_accepted_cells",
                guard_accepted,
            )
            world.set_g(
                "tectonics.last_p37_shape_guard_rejected_cells",
                guard_rejected,
            )
            if guard_input > 0.0:
                world.set_g(
                    "tectonics.cumulative_p37_shape_guard_input_cells",
                    world.g("tectonics.cumulative_p37_shape_guard_input_cells", 0.0)
                    + guard_input,
                )
                world.set_g(
                    "tectonics.cumulative_p37_shape_guard_accepted_cells",
                    world.g("tectonics.cumulative_p37_shape_guard_accepted_cells", 0.0)
                    + guard_accepted,
                )
                world.set_g(
                    "tectonics.cumulative_p37_shape_guard_rejected_cells",
                    world.g("tectonics.cumulative_p37_shape_guard_rejected_cells", 0.0)
                    + guard_rejected,
                )
                if bool(world.g("tectonics.enable_lineage_guarded_debt_recovery", 0.0)):
                    world.set_g(
                        "tectonics.cumulative_p38_lineage_guard_input_cells",
                        world.g("tectonics.cumulative_p38_lineage_guard_input_cells", 0.0)
                        + world.g("tectonics.last_p38_lineage_guard_input_cells", 0.0),
                    )
                    world.set_g(
                        "tectonics.cumulative_p38_lineage_guard_accepted_cells",
                        world.g("tectonics.cumulative_p38_lineage_guard_accepted_cells", 0.0)
                        + world.g("tectonics.last_p38_lineage_guard_accepted_cells", 0.0),
                    )
                    world.set_g(
                        "tectonics.cumulative_p38_lineage_guard_rejected_cells",
                        world.g("tectonics.cumulative_p38_lineage_guard_rejected_cells", 0.0)
                        + world.g("tectonics.last_p38_lineage_guard_rejected_cells", 0.0),
                    )
            out = guarded
        return out

    def _shape_guard_continental_restoration(self, world, grid, cont, candidate):
        """Accept restoration only when it preserves broad parent geometry.

        Continental-area debt repair must not turn a low-area world into many
        narrow preserved fragments.  This guard evaluates each candidate patch
        against its single parent continent component, using graph-width and
        narrow-area ratios as object-level shape proxies.
        """
        cont = np.asarray(cont, dtype=bool)
        candidate = np.asarray(candidate, dtype=bool) & ~cont
        out = np.zeros(grid.n, dtype=bool)
        if not candidate.any() or not cont.any():
            return out
        lineage_required = bool(
            world.g("tectonics.enable_lineage_guarded_debt_recovery", 0.0)
        )
        if lineage_required:
            world.set_g(
                "tectonics.last_p38_lineage_guard_input_cells",
                float(np.count_nonzero(candidate)),
            )
        labels, _ = _component_labels(grid, cont)
        parent_width = _graph_width_steps(grid, cont)
        area = grid.cell_area
        for patch in self._connected_components(grid, candidate):
            if patch.size == 0:
                continue
            patch_mask = np.zeros(grid.n, dtype=bool)
            patch_mask[patch] = True
            parent_ids: set[int] = set()
            for c in patch:
                for nb in grid.neighbors[int(c)]:
                    pid = int(labels[int(nb)])
                    if pid >= 0:
                        parent_ids.add(pid)
            if len(parent_ids) != 1:
                continue
            parent_id = next(iter(parent_ids))
            parent_mask = labels == parent_id
            parent_cells = np.where(parent_mask)[0]
            if parent_cells.size == 0:
                continue
            p50 = float(np.percentile(parent_width[parent_cells], 50))
            p90 = float(np.percentile(parent_width[parent_cells], 90))
            if p50 < 2.75 or p90 < 4.0:
                continue
            parent_area = max(float(area[parent_mask].sum()), 1e-12)
            parent_narrow = float(
                area[parent_mask & (parent_width <= 2.0)].sum() / parent_area)
            trial = parent_mask | patch_mask
            trial_width = _graph_width_steps(grid, trial)
            trial_cells = np.where(trial)[0]
            trial_area = max(float(area[trial].sum()), 1e-12)
            trial_narrow = float(
                area[trial & (trial_width <= 2.0)].sum() / trial_area)
            trial_p50 = float(np.percentile(trial_width[trial_cells], 50))
            if trial_narrow > max(parent_narrow + 0.08, 0.50):
                continue
            if trial_p50 < max(2.0, p50 - 1.0):
                continue
            if lineage_required and not self._lineage_guard_parent_component(
                world, grid, parent_mask, parent_width):
                continue
            out[patch] = True
        if lineage_required:
            world.set_g(
                "tectonics.last_p38_lineage_guard_accepted_cells",
                float(np.count_nonzero(out)),
            )
            world.set_g(
                "tectonics.last_p38_lineage_guard_rejected_cells",
                float(np.count_nonzero(candidate & ~out)),
            )
        return out

    def _lineage_guard_parent_component(self, world, grid, parent_mask,
                                        parent_width):
        """Require a restoration parent to map to a broad continent object."""
        previous = world.get_field("tectonics.continent_id", -1.0)
        previous = np.asarray(previous, dtype=np.float64)
        if previous.shape != (grid.n,):
            return False
        ids = previous[parent_mask].astype(int)
        valid = ids >= 0
        if not valid.any():
            return False
        area = grid.cell_area[parent_mask]
        parent_area = max(float(area.sum()), 1e-12)
        valid_ids = ids[valid]
        valid_area = area[valid]
        unique_ids = sorted(int(x) for x in np.unique(valid_ids) if int(x) >= 0)
        if not unique_ids:
            return False
        weighted = []
        for cid in unique_ids:
            cid_area = float(valid_area[valid_ids == cid].sum())
            weighted.append((cid_area, cid))
        weighted.sort(reverse=True)
        dominant_area, dominant_id = weighted[0]
        if dominant_area / parent_area < 0.62:
            return False

        objects = {
            int(obj.get("numeric_id", -1)): obj
            for obj in world.objects.get("tectonics.continents", [])
            if int(obj.get("numeric_id", -1)) >= 0
        }
        obj = objects.get(int(dominant_id))
        if obj is None:
            return False
        obj_p50 = float(obj.get("width_p50_steps", 0.0))
        obj_p90 = float(obj.get("width_p90_steps", 0.0))
        if obj_p50 < 2.75 or obj_p90 < 4.0:
            return False
        parent_cells = np.where(parent_mask)[0]
        if parent_cells.size == 0:
            return False
        local_p50 = float(np.percentile(parent_width[parent_cells], 50))
        local_p90 = float(np.percentile(parent_width[parent_cells], 90))
        return local_p50 >= 2.75 and local_p90 >= 4.0

    def _protected_breakup_corridor_mask(self, world, grid, passes=1):
        protected = np.zeros(grid.n, dtype=bool)
        for object_set in (
            "tectonics.breakup_seaways",
            "terrain.breakup_seaway_opened_corridors",
        ):
            for obj in world.objects.get(object_set, []):
                cells = obj.get("cells", [])
                if not cells:
                    continue
                arr = np.asarray(cells, dtype=np.int64)
                arr = arr[(arr >= 0) & (arr < grid.n)]
                protected[arr] = True
        if protected.any() and passes > 0:
            protected = self._dilate_mask(
                grid, protected, allowed=np.ones(grid.n, dtype=bool), passes=passes)
        return protected

    def _parented_platform_progradation_zone(self, world, grid, ctype, age,
                                             origin, reworked, stability,
                                             orog_age, volc_age, erosion_zone,
                                             shape_pressure, target, frac, t):
        """Parent-supported passive platform restoration.

        This is the upstream counterpart to terrain-level coastline payback:
        when mature Earth-like continents are underwide and below the target
        continental budget, quiet cells next to a stable parent margin can
        inherit continental cargo.  It excludes active rifts/ridges and explicit
        breakup corridors, so it represents shelf/platform progradation or old
        margin welding rather than arbitrary ocean-to-land painting.
        """
        out = np.zeros(grid.n, dtype=bool)
        mature_start_myr = max(3000.0, 0.98 * float(world.spec.t_end_myr))
        if float(t) < mature_start_myr:
            return out
        gap = float(target) - float(frac)
        if gap < 0.025 or (shape_pressure < 0.35 and gap < 0.06):
            return out
        cont = ctype == CONT
        if not cont.any():
            return out
        width = _graph_width_steps(grid, cont)
        stable_parent = cont & (
            (origin == ORIGIN_CRATON)
            | ((origin == ORIGIN_PRIMORDIAL) & (stability >= 0.45))
            | ((width >= 3.0) & (stability >= 0.56))
            | (stability >= 0.72)
        )
        if not stable_parent.any():
            return out

        cont_neighbors = _same_neighbor_count(grid, cont)
        stable_neighbors = _same_neighbor_count(grid, stable_parent)
        labels, _ = _component_labels(grid, cont)
        single_parent_component = np.zeros(grid.n, dtype=bool)
        for c in np.where(~cont & (cont_neighbors >= 1) & (stable_neighbors >= 1))[0]:
            parent_ids = {
                int(labels[int(nb)])
                for nb in grid.neighbors[int(c)]
                if cont[int(nb)] and labels[int(nb)] >= 0
            }
            single_parent_component[int(c)] = len(parent_ids) <= 1

        quiet_myr = min(700.0, max(320.0, 0.16 * float(t)))
        quiet_rework = (reworked < 0.0) | (reworked <= float(t) - quiet_myr)
        quiet_orogen = (orog_age < 0.0) | (orog_age <= float(t) - quiet_myr)
        quiet_volcanic = (volc_age < 0.0) | (volc_age <= float(t) - quiet_myr)
        quiet = quiet_rework & quiet_orogen & quiet_volcanic

        old_attached_arc = (
            (~cont)
            & (origin == ORIGIN_ARC)
            & (age >= 140.0)
            & quiet
            & (stable_neighbors >= 1)
        )
        parented_platform_shell = (
            (~cont)
            & quiet
            & (age >= 60.0)
            & (
                (cont_neighbors >= 3)
                | ((cont_neighbors >= 2) & (stable_neighbors >= 2))
            )
        )
        apron = self._dilate_mask(grid, stable_parent, allowed=~cont, passes=2)
        eligible = (
            apron
            & single_parent_component
            & (old_attached_arc | parented_platform_shell)
        )

        protected = self._protected_breakup_corridor_mask(world, grid, passes=1)
        if erosion_zone is not None:
            erosion = np.asarray(erosion_zone, dtype=bool)
            if erosion.shape == (grid.n,) and erosion.any():
                protected |= self._dilate_mask(
                    grid, erosion, allowed=np.ones(grid.n, dtype=bool), passes=1)
        world.set_g(
            "tectonics.last_p34_parented_platform_protected_cells",
            float(np.count_nonzero(eligible & protected)),
        )
        eligible &= ~protected
        if eligible.any():
            max_area = float(grid.cell_area.sum()) * min(
                0.0040,
                max(0.0012, 0.0012 + 0.0060 * max(gap, 0.0)),
            )
            cells = sorted((int(c) for c in np.where(eligible)[0]), key=lambda c: (
                0 if parented_platform_shell[c] else 1,
                -int(cont_neighbors[c]),
                -int(stable_neighbors[c]),
                0 if old_attached_arc[c] else 1,
                -float(age[c]),
                c,
            ))
            acc_area = 0.0
            for c in cells:
                if acc_area >= max_area:
                    break
                out[c] = True
                acc_area += float(grid.cell_area[c])
        world.set_g(
            "tectonics.last_p34_parented_platform_candidate_cells",
            float(np.count_nonzero(out)),
        )
        return out

    def _passive_margin_progradation_zone(self, grid, ctype, origin, stability,
                                          erosion_zone, shape_pressure):
        """Cells where continental cargo can inherit from a broad parent margin.

        This is a deterministic process proxy for passive-margin/platform
        progradation and fixed-grid cargo restoration.  It only activates when
        existing continents are too narrow or unstable, and it excludes current
        rift/ridge erosion zones so breakup is not immediately undone.
        """
        cont = ctype == CONT
        out = np.zeros(grid.n, dtype=bool)
        if shape_pressure < 0.35 or not cont.any():
            return out
        width = _graph_width_steps(grid, cont)
        area = grid.cell_area
        cont_area = float(area[cont].sum())
        if cont_area <= 0.0:
            return out
        unstable_accretionary = (
            cont
            & (width <= 2.0)
            & (
                np.isin(origin, [ORIGIN_ARC, ORIGIN_SUTURE])
                | (stability < 0.30)
            )
        )
        unstable_fraction = float(area[unstable_accretionary].sum() / cont_area)
        if unstable_fraction < 0.04:
            return out
        stable_parent = cont & (
            ((width >= 3.0) & (stability > 0.50))
            | (stability > 0.62)
            | (origin == ORIGIN_CRATON)
        )
        if not stable_parent.any():
            return out
        cont_neighbors = _same_neighbor_count(grid, cont)
        stable_neighbors = _same_neighbor_count(grid, stable_parent)
        labels, _ = _component_labels(grid, cont)
        single_parent_component = np.zeros(grid.n, dtype=bool)
        for c in np.where(~cont & (cont_neighbors >= 2) & (stable_neighbors >= 1))[0]:
            parent_ids = {
                int(labels[int(nb)])
                for nb in grid.neighbors[int(c)]
                if cont[int(nb)] and labels[int(nb)] >= 0
            }
            single_parent_component[int(c)] = len(parent_ids) <= 1
        passes = 2 if shape_pressure > 0.70 else 1
        apron = self._dilate_mask(grid, stable_parent, allowed=~cont, passes=passes)
        out = (
            apron
            & ~cont
            & single_parent_component
            & (
                (cont_neighbors >= 3)
                | ((cont_neighbors >= 2) & (stable_neighbors >= 2))
            )
        )
        if erosion_zone is not None:
            erosion = np.asarray(erosion_zone, dtype=bool)
            if erosion.shape == (grid.n,) and erosion.any():
                protected_breakup = self._dilate_mask(
                    grid, erosion, allowed=np.ones(grid.n, dtype=bool), passes=1)
                out &= ~protected_breakup
        return out

    def _shape_aware_continental_maintenance(self, grid, ctype, thick, age, origin,
                                             reworked, stability, orog_age, volc_age,
                                             t, budget, add_allowed=None,
                                             erase_allowed=None,
                                             active_accretion_allowed=None,
                                             world=None):
        """Balanced reshape: erode unsupported ribbons and fill supported margins."""
        if budget <= 0.0:
            return ctype, thick, age, origin, reworked, stability, orog_age, volc_age
        area = grid.cell_area
        cont = ctype == CONT
        if not cont.any() or (~cont).sum() == 0:
            return ctype, thick, age, origin, reworked, stability, orog_age, volc_age
        width = _graph_width_steps(grid, cont)
        before_cont_area = max(float(area[cont].sum()), 1e-12)
        before_narrow_fraction = float(area[cont & (width <= 2.0)].sum() / before_cont_area)
        cont_neighbors = _same_neighbor_count(grid, cont)
        craton = cont & ((origin == ORIGIN_CRATON) | (stability > 0.78))
        neighbor_degree = np.asarray([len(n) for n in grid.neighbors])
        non_cont_neighbors = neighbor_degree - cont_neighbors
        border_cont = cont & (cont_neighbors < neighbor_degree)
        component_labels, _ = _component_labels(grid, cont)
        single_parent_component = np.zeros(grid.n, dtype=bool)
        for c in np.where(~cont & (cont_neighbors >= 3))[0]:
            parent_ids = {
                int(component_labels[int(nb)])
                for nb in grid.neighbors[int(c)]
                if cont[int(nb)] and component_labels[int(nb)] >= 0
            }
            single_parent_component[int(c)] = len(parent_ids) <= 1

        broad_parent = (
            cont
            & (width >= 3.0)
            & (
                craton
                | ((origin == ORIGIN_PRIMORDIAL) & (stability >= 0.38))
                | ((origin == ORIGIN_SUTURE) & (stability >= 0.34))
                | (stability >= 0.58)
            )
        )
        if world is not None:
            funnel = self._old_platform_anchor_funnel(
                world, grid, ctype, thick, age, origin, reworked, stability, t)
            broad_parent |= (
                funnel["masks"]["old_platform_core"]
                | funnel["masks"]["old_platform_margin"]
                | funnel["masks"]["mature_platform"]
            )
            broad_parent &= cont & (width >= 3.0)
        support_allowed = (
            cont
            & (
                craton
                | ((origin == ORIGIN_PRIMORDIAL) & (stability >= 0.34))
                | ((origin == ORIGIN_SUTURE) & (stability >= 0.46))
                | (stability >= 0.58)
            )
        )
        parent_support = self._dilate_mask(
            grid, broad_parent, allowed=support_allowed, passes=2) if broad_parent.any() else broad_parent
        parent_neighbors = _same_neighbor_count(grid, parent_support)

        base_erase_mask = (
            border_cont
            & ~craton
            & (width <= 2.0)
            & (
                (cont_neighbors <= 4)
                | np.isin(origin, [ORIGIN_ARC, ORIGIN_SUTURE])
                | (stability < 0.45)
            )
        )
        base_add_mask = (
            (~cont)
            & (cont_neighbors >= 3)
            & (non_cont_neighbors <= 2)
            & single_parent_component
        )
        if erase_allowed is not None:
            allowed = np.asarray(erase_allowed, dtype=bool)
            if allowed.shape == (grid.n,):
                base_erase_mask &= allowed
        if add_allowed is not None:
            allowed = np.asarray(add_allowed, dtype=bool)
            if allowed.shape == (grid.n,):
                base_add_mask &= allowed

        p50_erase_mask = np.zeros(grid.n, dtype=bool)
        p50_add_mask = np.zeros(grid.n, dtype=bool)
        p50_add_parent_id = np.full(grid.n, -1, dtype=np.int64)
        shape_pressure = 0.0
        if world is not None:
            shape_pressure = float(world.g("tectonics.last_continent_shape_pressure", 0.0))
        mature_start_myr = 0.0
        if world is not None:
            mature_start_myr = max(3000.0, 0.98 * float(world.spec.t_end_myr))
        p50_enabled = (
            world is not None
            and world.g("tectonics.enable_planform_rebalance", 1.0) > 0.0
            and float(t) >= mature_start_myr
            and shape_pressure >= 0.18
            and broad_parent.any()
        )
        if p50_enabled:
            single_broad_parent = np.zeros(grid.n, dtype=bool)
            broad_adjacent = np.where((~cont) & (cont_neighbors >= 1) & (parent_neighbors >= 1))[0]
            for c in broad_adjacent:
                parent_ids = {
                    int(component_labels[int(nb)])
                    for nb in grid.neighbors[int(c)]
                    if cont[int(nb)] and component_labels[int(nb)] >= 0
                }
                if len(parent_ids) == 1:
                    single_broad_parent[int(c)] = True
                    p50_add_parent_id[int(c)] = next(iter(parent_ids))

            protected = self._protected_breakup_corridor_mask(world, grid, passes=1)
            active_accretion_zone = np.zeros(grid.n, dtype=bool)
            if active_accretion_allowed is not None:
                active = np.asarray(active_accretion_allowed, dtype=bool)
                if active.shape == (grid.n,) and active.any():
                    active_accretion_zone = self._dilate_mask(
                        grid, active, allowed=np.ones(grid.n, dtype=bool), passes=1)
            active_erosion_zone = np.zeros(grid.n, dtype=bool)
            if erase_allowed is not None:
                erosion = np.asarray(erase_allowed, dtype=bool)
                if erosion.shape == (grid.n,) and erosion.any():
                    active_erosion_zone = self._dilate_mask(
                        grid, erosion, allowed=np.ones(grid.n, dtype=bool), passes=1)

            broad_apron = self._dilate_mask(grid, parent_support, allowed=~cont, passes=1)
            p50_add_mask = (
                (~cont)
                & broad_apron
                & single_broad_parent
                & ~protected
                & ~active_erosion_zone
                & (non_cont_neighbors <= 2)
                & (
                    ((cont_neighbors >= 3) & (parent_neighbors >= 1))
                    | ((cont_neighbors >= 2) & (parent_neighbors >= 2))
                )
            )
            p50_raw_erase = (
                border_cont
                & ~craton
                & (width <= 2.0)
                & ~protected
                & ~active_accretion_zone
                & (
                    np.isin(origin, [ORIGIN_ARC, ORIGIN_SUTURE])
                    | (stability < 0.52)
                    | (cont_neighbors <= 3)
                )
                & (
                    (parent_neighbors <= 1)
                    | (cont_neighbors <= 3)
                    | (stability < 0.36)
                )
            )
            safe_terminal = np.zeros(grid.n, dtype=bool)
            for c in np.where(p50_raw_erase)[0]:
                c = int(c)
                cont_nbrs = [
                    int(nb) for nb in grid.neighbors[c]
                    if cont[int(nb)]
                ]
                if len(cont_nbrs) <= 1:
                    safe_terminal[c] = True
                elif len(cont_nbrs) == 2:
                    a, b = cont_nbrs
                    if b in grid.neighbors[a]:
                        safe_terminal[c] = True
                        continue
                    common = set(int(x) for x in grid.neighbors[a]).intersection(
                        int(x) for x in grid.neighbors[b])
                    if any(x != c and cont[x] for x in common):
                        safe_terminal[c] = True
            p50_erase_mask = p50_raw_erase & safe_terminal
            if world is not None:
                world.set_g(
                    "tectonics.last_p50_planform_recycle_candidate_cells",
                    float(np.count_nonzero(p50_erase_mask)),
                )
                world.set_g(
                    "tectonics.last_p50_planform_fill_candidate_cells",
                    float(np.count_nonzero(p50_add_mask)),
                )

        erase_mask = base_erase_mask | p50_erase_mask
        add_mask = base_add_mask | p50_add_mask
        erase = np.where(erase_mask)[0]
        add = np.where(add_mask)[0]
        if erase.size == 0 or add.size == 0:
            if world is not None:
                world.set_g(
                    "tectonics.last_p50_planform_before_narrow_fraction",
                    before_narrow_fraction,
                )
                world.set_g(
                    "tectonics.last_p50_planform_after_narrow_fraction",
                    before_narrow_fraction,
                )
            return ctype, thick, age, origin, reworked, stability, orog_age, volc_age

        stable_neighbors = _same_neighbor_count(grid, craton | (cont & (stability > 0.70)))
        erase_order = sorted((int(c) for c in erase), key=lambda c: (
            0 if p50_erase_mask[c] else 1,
            float(width[c]),
            0 if origin[c] == ORIGIN_ARC else 1,
            int(cont_neighbors[c]),
            float(stability[c]),
            float(thick[c]),
            c,
        ))
        add_order = sorted((int(c) for c in add), key=lambda c: (
            0 if p50_add_mask[c] else 1,
            -int(parent_neighbors[c]) if p50_add_mask[c] else 0,
            -int(cont_neighbors[c]),
            -int(stable_neighbors[c]),
            -float(thick[c]),
            c,
        ))
        active_accretion_mask = None
        if active_accretion_allowed is not None:
            allowed = np.asarray(active_accretion_allowed, dtype=bool)
            if allowed.shape == (grid.n,):
                active_accretion_mask = allowed
        inherit_rework_parent = self._quiet_mature_rework_inheritance_parent_mask(
            grid, ctype, thick, age, origin, reworked, stability, t)
        preferred_parent = inherit_rework_parent.copy()
        if p50_enabled:
            preferred_parent |= parent_support

        changed_area = 0.0
        inherited_rework = 0
        reactivated = 0
        p50_recycled = 0
        p50_filled = 0
        p50_area = 0.0
        p50_added_parent_id = np.full(grid.n, -1, dtype=np.int64)
        n = min(len(erase_order), len(add_order))
        for idx in range(n):
            if changed_area >= budget:
                break
            e = erase_order[idx]
            a = add_order[idx]
            if ctype[e] != CONT or ctype[a] == CONT:
                continue
            if p50_add_mask[a]:
                pid = int(p50_add_parent_id[a])
                if pid >= 0:
                    bridge_conflict = False
                    for nb in grid.neighbors[a]:
                        nb = int(nb)
                        old_pid = int(component_labels[nb])
                        new_pid = int(p50_added_parent_id[nb])
                        if ctype[nb] == CONT and old_pid >= 0 and old_pid != pid:
                            bridge_conflict = True
                            break
                        if new_pid >= 0 and new_pid != pid:
                            bridge_conflict = True
                            break
                    if bridge_conflict:
                        continue
            parent = self._continental_margin_parent(
                grid, a, ctype, thick, age, origin, stability,
                preferred=preferred_parent)
            p50_pair = bool(p50_erase_mask[e] or p50_add_mask[a])
            ctype[e] = OCEAN
            thick[e] = OCEAN_THICK
            age[e] = min(age[e], 80.0)
            origin[e] = ORIGIN_RIDGE
            reworked[e] = t
            stability[e] = 0.0
            orog_age[e] = -1.0
            volc_age[e] = -1.0

            ctype[a] = CONT
            active_accretion = (
                active_accretion_mask is not None
                and bool(active_accretion_mask[a])
                and not p50_pair
            )
            if active_accretion:
                thick[a] = max(thick[a], CONT_THICK)
                age[a] = max(age[a], min(900.0, t))
                origin[a] = ORIGIN_ARC
                stability[a] = max(stability[a], 0.25)
                reworked[a] = t
                reactivated += 1
            elif parent is None:
                thick[a] = max(thick[a], CONT_THICK)
                age[a] = max(age[a], min(900.0, t))
                origin[a] = ORIGIN_PRIMORDIAL
                stability[a] = max(stability[a], 0.35)
                reworked[a] = t
                reactivated += 1
            else:
                p = int(parent)
                parent_origin = origin[p]
                inherited_origin = (
                    ORIGIN_PRIMORDIAL if parent_origin == ORIGIN_CRATON
                    else parent_origin
                )
                origin[a] = inherited_origin
                thick[a] = max(thick[a], CONT_THICK, min(thick[p], 56000.0))
                age[a] = max(age[a], min(age[p], t))
                inherited_stability = stability[p] * (
                    0.86 if parent_origin == ORIGIN_CRATON else 0.92
                )
                stability[a] = max(stability[a], min(inherited_stability, 0.76))
                if orog_age[p] >= 0.0:
                    orog_age[a] = orog_age[p]
                if parent_origin == ORIGIN_ARC and volc_age[p] >= 0.0:
                    volc_age[a] = volc_age[p]
                if inherit_rework_parent[p]:
                    reworked[a] = reworked[p]
                    inherited_rework += 1
                else:
                    reworked[a] = t
                    reactivated += 1
            pair_area = 0.5 * (float(area[e]) + float(area[a]))
            changed_area += pair_area
            if p50_pair:
                p50_recycled += 1
                p50_filled += 1
                p50_area += pair_area
                if p50_add_mask[a]:
                    p50_added_parent_id[a] = int(p50_add_parent_id[a])
        if world is not None:
            world.set_g(
                "tectonics.last_p40_shape_maintenance_inherited_rework_cells",
                world.g("tectonics.last_p40_shape_maintenance_inherited_rework_cells", 0.0)
                + float(inherited_rework),
            )
            world.set_g(
                "tectonics.last_p40_shape_maintenance_reactivated_cells",
                world.g("tectonics.last_p40_shape_maintenance_reactivated_cells", 0.0)
                + float(reactivated),
            )
            after_cont = ctype == CONT
            after_width = _graph_width_steps(grid, after_cont)
            after_cont_area = max(float(area[after_cont].sum()), 1e-12)
            after_narrow_fraction = float(
                area[after_cont & (after_width <= 2.0)].sum() / after_cont_area)
            world.set_g("tectonics.last_p50_planform_recycled_cells", float(p50_recycled))
            world.set_g("tectonics.last_p50_planform_filled_cells", float(p50_filled))
            world.set_g(
                "tectonics.last_p50_planform_area_fraction",
                float(p50_area / max(float(area.sum()), 1.0)),
            )
            world.set_g(
                "tectonics.last_p50_planform_before_narrow_fraction",
                before_narrow_fraction,
            )
            world.set_g(
                "tectonics.last_p50_planform_after_narrow_fraction",
                after_narrow_fraction,
            )
            world.set_g(
                "tectonics.last_p50_planform_narrow_fraction_delta",
                float(before_narrow_fraction - after_narrow_fraction),
            )
            if p50_recycled > 0:
                world.set_g(
                    "tectonics.cumulative_p50_planform_recycled_cells",
                    world.g("tectonics.cumulative_p50_planform_recycled_cells", 0.0)
                    + float(p50_recycled),
                )
                world.set_g(
                    "tectonics.cumulative_p50_planform_filled_cells",
                    world.g("tectonics.cumulative_p50_planform_filled_cells", 0.0)
                    + float(p50_filled),
                )
                world.set_g(
                    "tectonics.cumulative_p50_planform_area_fraction",
                    world.g("tectonics.cumulative_p50_planform_area_fraction", 0.0)
                    + float(p50_area / max(float(area.sum()), 1.0)),
                )
        return ctype, thick, age, origin, reworked, stability, orog_age, volc_age

    def _continental_margin_parent(self, grid, cell, ctype, thick, age,
                                   origin, stability, preferred=None):
        """Pick the deterministic parent crust for a shape-maintenance fill.

        Background shape maintenance moves a fixed-grid continental edge; it is
        not a new island-arc accretion event.  The added cell should therefore
        inherit the strongest adjacent continental margin instead of being
        labelled as newly generated arc crust.
        """
        candidates = [
            int(nb) for nb in grid.neighbors[int(cell)]
            if ctype[int(nb)] == CONT
        ]
        if not candidates:
            return None
        preferred_mask = None
        if preferred is not None:
            arr = np.asarray(preferred, dtype=bool)
            if arr.shape == (grid.n,):
                preferred_mask = arr
        candidates.sort(key=lambda c: (
            0 if preferred_mask is not None and preferred_mask[c] else 1,
            0 if origin[c] == ORIGIN_CRATON else 1,
            -float(stability[c]),
            -float(age[c]),
            -float(thick[c]),
            c,
        ))
        return candidates[0]

    def _crust_stability(self, ctype, age, origin):
        cont = ctype == CONT
        age_term = np.clip((age - 600.0) / 2200.0, 0.0, 1.0)
        origin_bonus = np.where(origin == ORIGIN_CRATON, 0.35,
                                np.where(origin == ORIGIN_PRIMORDIAL, 0.18, 0.0))
        stability = np.where(cont, 0.12 + 0.65 * age_term + origin_bonus, 0.0)
        stable_craton = cont & (age > 2500.0) & (stability > 0.70)
        out = np.clip(stability, 0.0, 1.0)
        out[stable_craton] = np.maximum(out[stable_craton], 0.85)
        return out

    def _shape_guard_cratonic_stability(self, grid, ctype, origin, stability):
        """Prevent narrow continental ribbons from masquerading as craton core."""
        del origin
        cont = np.asarray(ctype, dtype=np.float64) == CONT
        out = np.asarray(stability, dtype=np.float64).copy()
        if out.shape != (grid.n,) or not cont.any():
            return out
        width = _graph_width_steps(grid, cont)
        broad = cont & (width >= 3.0)
        broad_support = self._dilate_mask(grid, broad, allowed=cont, passes=1)
        narrow = cont & (width <= 2.0)
        unsupported_narrow = narrow & ~broad_support
        out[narrow] = np.minimum(out[narrow], 0.74)
        out[unsupported_narrow] = np.minimum(out[unsupported_narrow], 0.62)
        return out

    def _adjacent_to_continent(self, grid, ctype) -> np.ndarray:
        out = np.zeros(grid.n, dtype=bool)
        is_cont = ctype == CONT
        edges = grid.edges
        i, j = edges[:, 0], edges[:, 1]
        ci, cj = is_cont[i], is_cont[j]
        out[i[cj]] = True
        out[j[ci]] = True
        return out

    def _plate_velocity(self, grid, plate, plates, vigor, regime_code) -> np.ndarray:
        regime_factor = {0.0: 0.02, 1.0: 0.4, 2.0: 0.7, 3.0: 1.0}.get(regime_code, 1.0)
        vel = np.zeros((grid.n, 3))
        sec_per_myr = 3.15576e13
        for i, p in enumerate(plates):
            mask = plate == i
            if not mask.any():
                continue
            omega = np.asarray(p["pole"]) * p["rate"] * max(vigor, 0.05) * regime_factor
            omega_si = omega / sec_per_myr            # rad / s
            v = np.cross(omega_si, grid.xyz[mask] * grid.radius_m)
            vel[mask] = v
        return vel

    def _sample_events(self, grid, t, rng, mask, type_, cause, age_field, k=None):
        cells = np.where(mask)[0]
        if cells.size == 0:
            return []
        if k is None:
            # Emit representative events, not a fixed four per tectonic step.
            # Counts now scale with active boundary area and naturally vary
            # across time and tectonic regime.
            k = int(np.clip(cells.size // 450 + 1, 1, 8))
        pick = cells if cells.size <= k else rng.choice(cells, size=k, replace=False)
        out = []
        for c in pick:
            out.append(Event(type_, t, self.name, location=int(c),
                             magnitude=float(cells.size),
                             params={"cause": cause,
                                     "lat": round(float(grid.lat[c]), 1),
                                     "lon": round(float(grid.lon[c]), 1)}))
        return out
