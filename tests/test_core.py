from collections import deque

import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.rng import derive_seed, generator
from aevum.compiler.hexgrid import HexGrid
from aevum.diagnostics.morphology import analyze_mask_morphology
from aevum.features import build_registry
from aevum.modules.tectonics import (
    CONT,
    CONT_THICK,
    OCEAN,
    OCEAN_THICK,
    ORIGIN_ARC,
    ORIGIN_PRIMORDIAL,
    TectonicsModule,
    _grow_compact_continents,
    _separated_seed_cells,
)


def test_registry_resolves_dependencies():
    reg = build_registry()
    assert reg.validate() == []          # no unknown dependencies (fatal)
    assert len(reg) > 40


def test_registry_feedback_loops_exist():
    # feedback loops are expected (climate<->weathering, terrain<->precip, ...)
    reg = build_registry()
    assert len(reg.feedback_loops()) >= 1


def test_grid_area_closes():
    g = SphereGrid.fibonacci(4000, 6.371e6)
    sphere = 4 * np.pi * g.radius_m ** 2
    assert abs(g.total_area - sphere) / sphere < 1e-3


def test_grid_neighbors_symmetric():
    g = SphereGrid.fibonacci(2000, 6.371e6)
    nb = g.neighbors
    for i in range(0, g.n, 137):
        for j in nb[i]:
            assert i in nb[j]


def test_rng_namespacing_is_stable():
    a = derive_seed(42, "tectonics", 100.0, 1)
    b = derive_seed(42, "tectonics", 100.0, 1)
    c = derive_seed(42, "climate", 100.0, 1)
    assert a == b and a != c
    g1 = generator(7, "x", 0.0)
    g2 = generator(7, "x", 0.0)
    assert np.array_equal(g1.random(5), g2.random(5))


def test_hex_grid_wraps_longitude_edges():
    hg = HexGrid(width=12, height=6)
    for r in range(hg.height):
        assert (r, hg.width - 1) in hg.neighbors(r, 0)
        assert (r, 0) in hg.neighbors(r, hg.width - 1)


def test_morphology_distinguishes_compact_land_from_ribbon():
    g = SphereGrid.fibonacci(1200, 6.371e6)
    center = np.array([1.0, 0.0, 0.0])
    angular_distance = np.degrees(np.arccos(np.clip(g.xyz @ center, -1.0, 1.0)))
    compact = angular_distance < 26.0
    ribbon = (np.abs(g.lat) < 3.2) & (np.abs(g.lon) < 125.0)

    compact_diag = analyze_mask_morphology(g, compact, name="compact")
    ribbon_diag = analyze_mask_morphology(g, ribbon, name="ribbon")

    assert compact_diag.metrics["component_count"] == 1
    assert ribbon_diag.metrics["component_count"] >= 1
    assert ribbon_diag.metrics["narrow_fraction_le2_width"] > compact_diag.metrics[
        "narrow_fraction_le2_width"
    ]
    assert ribbon_diag.metrics["ribbon_index_area_mean"] > compact_diag.metrics[
        "ribbon_index_area_mean"
    ] + 0.25
    assert ribbon_diag.metrics["coastline_complexity_p95_component"] > compact_diag.metrics[
        "coastline_complexity_p95_component"
    ]


def test_morphology_detects_narrow_articulation_neck():
    g = SphereGrid.fibonacci(1600, 6.371e6)

    def unit(lat_deg, lon_deg):
        lat = np.radians(lat_deg)
        lon = np.radians(lon_deg)
        return np.array([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ])

    west = unit(0.0, -35.0)
    east = unit(0.0, 35.0)
    west_blob = np.degrees(np.arccos(np.clip(g.xyz @ west, -1.0, 1.0))) < 15.0
    east_blob = np.degrees(np.arccos(np.clip(g.xyz @ east, -1.0, 1.0))) < 15.0
    start = int(np.argmax(g.xyz @ west))
    end = int(np.argmax(g.xyz @ east))
    parent = np.full(g.n, -1, dtype=int)
    q = deque([start])
    parent[start] = start
    while q and parent[end] < 0:
        c = q.popleft()
        for nb in g.neighbors[c]:
            nb = int(nb)
            if parent[nb] >= 0:
                continue
            parent[nb] = c
            q.append(nb)
    bridge = np.zeros(g.n, dtype=bool)
    c = end
    while c != start:
        bridge[c] = True
        c = int(parent[c])
    bridge[start] = True
    mask = west_blob | east_blob | bridge

    diag = analyze_mask_morphology(g, mask, name="bridge")
    assert diag.metrics["component_count"] == 1
    assert diag.metrics["narrow_neck_cells"] > 0
    assert diag.metrics["narrow_neck_cells_per_1000_mask_cells"] > 0.0


def test_compact_continent_growth_starts_with_broad_interiors():
    g = SphereGrid.fibonacci(2500, 6.371e6)
    rng = np.random.default_rng(123)
    seeds = _separated_seed_cells(g, 4, rng)
    mask, labels = _grow_compact_continents(g, seeds, 0.29, rng)
    diag = analyze_mask_morphology(g, mask, name="compact_initial_continents")

    assert int(mask.sum()) > 0
    assert int(np.unique(labels[mask]).size) >= 3
    assert 0.26 <= diag.metrics["mask_area_fraction_of_total"] <= 0.32
    assert diag.metrics["width_p50_steps"] >= 2.0
    assert diag.metrics["width_p90_steps"] >= 5.0
    assert diag.metrics["narrow_fraction_le2_width"] < 0.55
    assert diag.metrics["ribbon_area_fraction_gt_0_5"] < 0.35


def test_continental_shape_maintenance_reduces_unstable_ribbon():
    g = SphereGrid.fibonacci(1800, 6.371e6)
    base = (np.abs(g.lat) < 24.0) & (np.abs(g.lon) < 46.0)
    ribbon = (np.abs(g.lat) < 3.5) & (g.lon >= 46.0) & (g.lon < 145.0)
    mask = base | ribbon
    before = analyze_mask_morphology(g, mask, name="unstable_continental_ribbon")

    ctype = np.full(g.n, OCEAN)
    ctype[mask] = CONT
    thick = np.full(g.n, OCEAN_THICK)
    thick[mask] = CONT_THICK
    age = np.zeros(g.n)
    age[mask] = 1200.0
    origin = np.full(g.n, 0.0)
    origin[base] = ORIGIN_PRIMORDIAL
    origin[ribbon] = ORIGIN_ARC
    reworked = np.full(g.n, -1.0)
    stability = np.zeros(g.n)
    stability[base] = 0.68
    stability[ribbon] = 0.18
    orog_age = np.full(g.n, -1.0)
    volc_age = np.full(g.n, -1.0)

    module = TectonicsModule()
    assert module._continental_shape_pressure(g, ctype, origin, stability) > 0.5
    result = module._shape_aware_continental_maintenance(
        g, ctype, thick, age, origin, reworked, stability, orog_age, volc_age,
        t=1200.0, budget=float(g.cell_area.sum()) * 0.08,
    )
    after_mask = result[0] == CONT
    added = after_mask & ~mask
    after = analyze_mask_morphology(g, after_mask, name="healed_continental_ribbon")

    before_area = float(g.cell_area[mask].sum())
    after_area = float(g.cell_area[after_mask].sum())
    assert added.any()
    assert not np.any(result[3][added] == ORIGIN_ARC)
    assert np.all(result[3][added] == ORIGIN_PRIMORDIAL)
    assert float(np.median(result[6][added])) < 0.0
    assert float(np.median(result[5][added])) >= 0.55
    assert abs(after_area - before_area) / before_area < 0.01
    assert after.metrics["ribbon_area_fraction_gt_0_5"] < before.metrics[
        "ribbon_area_fraction_gt_0_5"
    ]
    assert after.metrics["coastline_complexity_largest_component"] < before.metrics[
        "coastline_complexity_largest_component"
    ]


def test_continental_shape_maintenance_does_not_bridge_separate_continents():
    g = SphereGrid.fibonacci(1800, 6.371e6)
    left = (np.abs(g.lat) < 20.0) & (g.lon > -62.0) & (g.lon < -8.0)
    right = (np.abs(g.lat) < 20.0) & (g.lon > 8.0) & (g.lon < 62.0)
    ribbon = (np.abs(g.lat) < 3.5) & (g.lon >= 62.0) & (g.lon < 142.0)
    mask = left | right | ribbon
    before = analyze_mask_morphology(g, mask, name="two_continents_with_ribbon")

    ctype = np.full(g.n, OCEAN)
    ctype[mask] = CONT
    thick = np.full(g.n, OCEAN_THICK)
    thick[mask] = CONT_THICK
    age = np.zeros(g.n)
    age[mask] = 1200.0
    origin = np.full(g.n, 0.0)
    origin[left | right] = ORIGIN_PRIMORDIAL
    origin[ribbon] = ORIGIN_ARC
    reworked = np.full(g.n, -1.0)
    stability = np.zeros(g.n)
    stability[left | right] = 0.68
    stability[ribbon] = 0.18
    orog_age = np.full(g.n, -1.0)
    volc_age = np.full(g.n, -1.0)

    result = TectonicsModule()._shape_aware_continental_maintenance(
        g, ctype, thick, age, origin, reworked, stability, orog_age, volc_age,
        t=1200.0, budget=float(g.cell_area.sum()) * 0.08,
    )
    after = analyze_mask_morphology(g, result[0] == CONT, name="protected_gateway")

    assert before.metrics["component_count"] >= 2
    assert after.metrics["component_count"] >= before.metrics["component_count"]
    assert after.metrics["largest_component_area_fraction_of_mask"] <= before.metrics[
        "largest_component_area_fraction_of_mask"
    ]
