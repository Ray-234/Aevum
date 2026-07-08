"""Deterministic province-graph fixtures for continental physiography.

The fixtures are intentionally small raster graphs.  They are not production
terrain generation, but they make the next terrain rewrite testable: province
objects must be spatially coherent, process-parented, adjacent in plausible
ways, and tied to expected elevation/sediment ordering.
"""
from __future__ import annotations

from collections import Counter, defaultdict, deque
from typing import Any


SCHEMA = "aevum.province_graph_fixtures.v1"

REQUIRED_FIXTURES = (
    "craton_platform_basin_fixture",
    "active_orogen_foreland_fixture",
    "old_suture_orogen_fixture",
    "rift_axis_shoulder_fixture",
    "passive_margin_lowland_fixture",
    "volcanic_lip_plateau_fixture",
    "multi_province_continent_fixture",
)

REQUIRED_PROVINCE_CLASSES = {
    "shield",
    "platform",
    "intracratonic_basin",
    "active_orogen",
    "foreland_basin",
    "old_orogen",
    "old_suture",
    "rift_axis",
    "rift_shoulder",
    "rift_basin",
    "passive_margin_lowland",
    "continental_shelf",
    "volcanic_lip_plateau",
}

REQUIRED_PARENT_PROCESSES = {
    "cratonization",
    "platform_subsidence",
    "intracratonic_sag",
    "collision_orogeny",
    "flexural_loading",
    "suture_inheritance",
    "orogenic_decay",
    "continental_extension",
    "rift_shoulder_uplift",
    "rift_basin_subsidence",
    "passive_margin_subsidence",
    "shelf_sedimentation",
    "plume_lip_emplacement",
}


def _rows(*rows: str) -> tuple[tuple[str, ...], ...]:
    return tuple(tuple(row.split()) for row in rows)


def _province(
    province_class: str,
    parent_processes: tuple[str, ...],
    elevation_m: float,
    sediment_m: float,
    parent_objects: tuple[str, ...],
) -> dict[str, Any]:
    return {
        "class": province_class,
        "parent_processes": parent_processes,
        "parent_objects": parent_objects,
        "elevation_m": float(elevation_m),
        "sediment_m": float(sediment_m),
        "source": "deterministic_process_fixture",
        "expected_contiguous": True,
    }


FIXTURES: tuple[dict[str, Any], ...] = (
    {
        "name": "craton_platform_basin_fixture",
        "grid": _rows(
            "SH SH SH PL PL IB",
            "SH SH SH PL PL IB",
            "SH SH PL PL IB IB",
            "SH PL PL IB IB IB",
        ),
        "provinces": {
            "SH": _province("shield", ("cratonization",), 650.0, 120.0, ("archean_craton",)),
            "PL": _province("platform", ("platform_subsidence",), 340.0, 650.0, ("covered_craton",)),
            "IB": _province("intracratonic_basin", ("intracratonic_sag",), 120.0, 2200.0, ("sag_basin",)),
        },
        "expected_classes": ("shield", "platform", "intracratonic_basin"),
        "expected_class_edges": (
            ("shield", "platform"),
            ("platform", "intracratonic_basin"),
        ),
        "elevation_order": (("shield", "platform", "intracratonic_basin"),),
        "sediment_order": (("intracratonic_basin", "platform", "shield"),),
    },
    {
        "name": "active_orogen_foreland_fixture",
        "grid": _rows(
            "PL PL FB FB AO AO",
            "PL PL FB FB AO AO",
            "PL PL FB FB AO AO",
            "PL PL FB FB AO AO",
        ),
        "provinces": {
            "PL": _province("platform", ("platform_subsidence",), 260.0, 700.0, ("foreland_hinterland_platform",)),
            "FB": _province("foreland_basin", ("flexural_loading",), 90.0, 3200.0, ("active_orogen_load",)),
            "AO": _province("active_orogen", ("collision_orogeny",), 3100.0, 260.0, ("convergent_boundary",)),
        },
        "expected_classes": ("platform", "foreland_basin", "active_orogen"),
        "expected_class_edges": (
            ("platform", "foreland_basin"),
            ("foreland_basin", "active_orogen"),
        ),
        "elevation_order": (("active_orogen", "platform", "foreland_basin"),),
        "sediment_order": (("foreland_basin", "platform", "active_orogen"),),
    },
    {
        "name": "old_suture_orogen_fixture",
        "grid": _rows(
            "PF PF OO OO ST SH SH",
            "PF PF OO OO ST SH SH",
            "PF PF OO OO ST SH SH",
            "PF PF OO OO ST SH SH",
        ),
        "provinces": {
            "PF": _province("platform", ("platform_subsidence",), 300.0, 850.0, ("post_orogenic_platform",)),
            "OO": _province("old_orogen", ("orogenic_decay", "collision_orogeny"), 920.0, 420.0, ("paleo_orogen",)),
            "ST": _province("old_suture", ("suture_inheritance",), 220.0, 350.0, ("closed_ocean_suture",)),
            "SH": _province("shield", ("cratonization",), 560.0, 150.0, ("opposing_craton",)),
        },
        "expected_classes": ("platform", "old_orogen", "old_suture", "shield"),
        "expected_class_edges": (
            ("platform", "old_orogen"),
            ("old_orogen", "old_suture"),
            ("old_suture", "shield"),
        ),
        "elevation_order": (("old_orogen", "shield", "platform", "old_suture"),),
        "sediment_order": (("platform", "old_orogen", "old_suture", "shield"),),
    },
    {
        "name": "rift_axis_shoulder_fixture",
        "grid": _rows(
            "PLW PLW RS1 RB1 RA RB2 RS2 PLE PLE",
            "PLW PLW RS1 RB1 RA RB2 RS2 PLE PLE",
            "PLW PLW RS1 RB1 RA RB2 RS2 PLE PLE",
            "PLW PLW RS1 RB1 RA RB2 RS2 PLE PLE",
        ),
        "provinces": {
            "PLW": _province("platform", ("platform_subsidence",), 310.0, 650.0, ("western_rifted_platform",)),
            "PLE": _province("platform", ("platform_subsidence",), 310.0, 650.0, ("eastern_rifted_platform",)),
            "RS1": _province("rift_shoulder", ("rift_shoulder_uplift", "continental_extension"), 980.0, 230.0, ("western_rift_shoulder",)),
            "RB1": _province("rift_basin", ("rift_basin_subsidence", "continental_extension"), -80.0, 2700.0, ("western_half_graben",)),
            "RA": _province("rift_axis", ("continental_extension",), 60.0, 1900.0, ("rift_axis",)),
            "RB2": _province("rift_basin", ("rift_basin_subsidence", "continental_extension"), -70.0, 2600.0, ("eastern_half_graben",)),
            "RS2": _province("rift_shoulder", ("rift_shoulder_uplift", "continental_extension"), 940.0, 260.0, ("eastern_rift_shoulder",)),
        },
        "expected_classes": ("platform", "rift_shoulder", "rift_basin", "rift_axis"),
        "expected_class_edges": (
            ("platform", "rift_shoulder"),
            ("rift_shoulder", "rift_basin"),
            ("rift_basin", "rift_axis"),
        ),
        "elevation_order": (("rift_shoulder", "platform", "rift_axis", "rift_basin"),),
        "sediment_order": (("rift_basin", "rift_axis", "platform", "rift_shoulder"),),
    },
    {
        "name": "passive_margin_lowland_fixture",
        "grid": _rows(
            "SH PL PL PM PM CS CS",
            "SH PL PL PM PM CS CS",
            "SH PL PL PM PM CS CS",
            "SH PL PL PM PM CS CS",
        ),
        "provinces": {
            "SH": _province("shield", ("cratonization",), 520.0, 130.0, ("continental_core",)),
            "PL": _province("platform", ("platform_subsidence",), 260.0, 800.0, ("coastal_platform",)),
            "PM": _province("passive_margin_lowland", ("passive_margin_subsidence",), 35.0, 1800.0, ("passive_margin_prism",)),
            "CS": _province("continental_shelf", ("shelf_sedimentation", "passive_margin_subsidence"), -120.0, 3000.0, ("continental_shelf",)),
        },
        "expected_classes": ("shield", "platform", "passive_margin_lowland", "continental_shelf"),
        "expected_class_edges": (
            ("shield", "platform"),
            ("platform", "passive_margin_lowland"),
            ("passive_margin_lowland", "continental_shelf"),
        ),
        "elevation_order": (("shield", "platform", "passive_margin_lowland", "continental_shelf"),),
        "sediment_order": (("continental_shelf", "passive_margin_lowland", "platform", "shield"),),
    },
    {
        "name": "volcanic_lip_plateau_fixture",
        "grid": _rows(
            "PLW PLW VL VL VL PLE",
            "PLW VL VL VL VL PLE",
            "PLW VL VL VL PLE PLE",
            "PLW PLW VL PLE PLE PLE",
        ),
        "provinces": {
            "PLW": _province("platform", ("platform_subsidence",), 330.0, 760.0, ("western_pre_lip_platform",)),
            "PLE": _province("platform", ("platform_subsidence",), 330.0, 760.0, ("eastern_pre_lip_platform",)),
            "VL": _province("volcanic_lip_plateau", ("plume_lip_emplacement",), 1450.0, 210.0, ("large_igneous_province",)),
        },
        "expected_classes": ("platform", "volcanic_lip_plateau"),
        "expected_class_edges": (("platform", "volcanic_lip_plateau"),),
        "elevation_order": (("volcanic_lip_plateau", "platform"),),
        "sediment_order": (("platform", "volcanic_lip_plateau"),),
    },
    {
        "name": "multi_province_continent_fixture",
        "grid": _rows(
            "SH SH PL1 PL1 IB IB PM1 CS1",
            "SH SH PL1 PL1 IB IB PM1 CS1",
            "OO OO ST PL1 PL1 PM1 PM1 CS1",
            "OO OO ST PL1 FB FB AO AO",
            "PL2 RS1 RB1 RA RB2 RS2 PM2 CS2",
            "PL2 RS1 RB1 RA RB2 RS2 PM2 CS2",
            "PL2 PL2 VL VL VL PM2 PM2 CS2",
            "PL2 PL2 VL VL PM2 PM2 CS2 CS2",
        ),
        "provinces": {
            "SH": _province("shield", ("cratonization",), 620.0, 130.0, ("continental_core",)),
            "PL1": _province("platform", ("platform_subsidence",), 310.0, 780.0, ("northern_covered_craton",)),
            "PL2": _province("platform", ("platform_subsidence",), 310.0, 780.0, ("southern_covered_craton",)),
            "IB": _province("intracratonic_basin", ("intracratonic_sag",), 110.0, 2300.0, ("sag_basin",)),
            "PM1": _province("passive_margin_lowland", ("passive_margin_subsidence",), 40.0, 1750.0, ("northern_passive_margin_prism",)),
            "PM2": _province("passive_margin_lowland", ("passive_margin_subsidence",), 40.0, 1750.0, ("southern_passive_margin_prism",)),
            "CS1": _province("continental_shelf", ("shelf_sedimentation", "passive_margin_subsidence"), -130.0, 3100.0, ("northern_shelf_prism",)),
            "CS2": _province("continental_shelf", ("shelf_sedimentation", "passive_margin_subsidence"), -130.0, 3100.0, ("southern_shelf_prism",)),
            "OO": _province("old_orogen", ("orogenic_decay", "collision_orogeny"), 850.0, 420.0, ("paleo_orogen",)),
            "ST": _province("old_suture", ("suture_inheritance",), 210.0, 330.0, ("suture_boundary",)),
            "FB": _province("foreland_basin", ("flexural_loading",), 95.0, 3300.0, ("active_orogen_load",)),
            "AO": _province("active_orogen", ("collision_orogeny",), 3000.0, 270.0, ("active_collision",)),
            "RS1": _province("rift_shoulder", ("rift_shoulder_uplift", "continental_extension"), 960.0, 240.0, ("western_rift_shoulder",)),
            "RB1": _province("rift_basin", ("rift_basin_subsidence", "continental_extension"), -80.0, 2700.0, ("western_half_graben",)),
            "RA": _province("rift_axis", ("continental_extension",), 50.0, 1900.0, ("rift_axis",)),
            "RB2": _province("rift_basin", ("rift_basin_subsidence", "continental_extension"), -70.0, 2600.0, ("eastern_half_graben",)),
            "RS2": _province("rift_shoulder", ("rift_shoulder_uplift", "continental_extension"), 930.0, 250.0, ("eastern_rift_shoulder",)),
            "VL": _province("volcanic_lip_plateau", ("plume_lip_emplacement",), 1380.0, 220.0, ("large_igneous_province",)),
        },
        "expected_classes": tuple(sorted(REQUIRED_PROVINCE_CLASSES)),
        "expected_class_edges": (
            ("shield", "platform"),
            ("platform", "intracratonic_basin"),
            ("platform", "passive_margin_lowland"),
            ("passive_margin_lowland", "continental_shelf"),
            ("old_orogen", "old_suture"),
            ("foreland_basin", "active_orogen"),
            ("rift_shoulder", "rift_basin"),
            ("rift_basin", "rift_axis"),
            ("platform", "volcanic_lip_plateau"),
        ),
        "elevation_order": (
            ("active_orogen", "volcanic_lip_plateau", "rift_shoulder", "old_orogen", "platform", "intracratonic_basin", "passive_margin_lowland", "continental_shelf"),
        ),
        "sediment_order": (
            ("foreland_basin", "continental_shelf", "rift_basin", "intracratonic_basin", "rift_axis", "passive_margin_lowland", "platform", "old_orogen", "shield"),
        ),
    },
)


def _pair(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


def _component_count(grid: tuple[tuple[str, ...], ...], province_id: str) -> int:
    cells = {
        (r, c)
        for r, row in enumerate(grid)
        for c, value in enumerate(row)
        if value == province_id
    }
    seen: set[tuple[int, int]] = set()
    count = 0
    for start in sorted(cells):
        if start in seen:
            continue
        count += 1
        queue: deque[tuple[int, int]] = deque([start])
        seen.add(start)
        while queue:
            r, c = queue.popleft()
            for nr, nc in ((r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)):
                if (nr, nc) in cells and (nr, nc) not in seen:
                    seen.add((nr, nc))
                    queue.append((nr, nc))
    return count


def _class_means(
    grid: tuple[tuple[str, ...], ...],
    provinces: dict[str, dict[str, Any]],
    field: str,
) -> dict[str, float]:
    totals: defaultdict[str, float] = defaultdict(float)
    counts: Counter[str] = Counter()
    for row in grid:
        for province_id in row:
            spec = provinces[province_id]
            province_class = str(spec["class"])
            totals[province_class] += float(spec[field])
            counts[province_class] += 1
    return {
        province_class: totals[province_class] / counts[province_class]
        for province_class in sorted(counts)
    }


def _ordering_passes(
    class_means: dict[str, float],
    orderings: tuple[tuple[str, ...], ...],
) -> bool:
    for ordering in orderings:
        values = [class_means[province_class] for province_class in ordering]
        if any(left <= right for left, right in zip(values, values[1:])):
            return False
    return True


def _analyze_fixture(fixture: dict[str, Any]) -> dict[str, Any]:
    grid = fixture["grid"]
    provinces = fixture["provinces"]
    grid_ids = {value for row in grid for value in row}
    declared_ids = set(provinces)
    missing_declarations = sorted(grid_ids - declared_ids)
    unused_declarations = sorted(declared_ids - grid_ids)

    classes_present = sorted({str(provinces[province_id]["class"]) for province_id in grid_ids})
    parent_processes = sorted({
        str(process)
        for province_id in grid_ids
        for process in provinces[province_id]["parent_processes"]
    })
    unparented = sorted(
        province_id
        for province_id in grid_ids
        if not provinces[province_id]["parent_processes"])
    random_texture_only = sorted(
        province_id
        for province_id in grid_ids
        if str(provinces[province_id].get("source")) == "random_texture")

    component_counts = {
        province_id: _component_count(grid, province_id)
        for province_id in sorted(grid_ids)
        if provinces[province_id].get("expected_contiguous", True)
    }
    non_contiguous = sorted(
        province_id for province_id, count in component_counts.items() if count != 1)

    province_edges: set[tuple[str, str]] = set()
    class_edges: set[tuple[str, str]] = set()
    for r, row in enumerate(grid):
        for c, province_id in enumerate(row):
            for nr, nc in ((r + 1, c), (r, c + 1)):
                if nr >= len(grid) or nc >= len(row):
                    continue
                other_id = grid[nr][nc]
                if other_id == province_id:
                    continue
                province_edges.add(_pair(province_id, other_id))
                class_edges.add(_pair(
                    str(provinces[province_id]["class"]),
                    str(provinces[other_id]["class"]),
                ))

    expected_classes = set(fixture["expected_classes"])
    expected_class_edges = {_pair(a, b) for a, b in fixture["expected_class_edges"]}
    elevation_means = _class_means(grid, provinces, "elevation_m")
    sediment_means = _class_means(grid, provinces, "sediment_m")

    missing_classes = sorted(expected_classes - set(classes_present))
    missing_class_edges = sorted(expected_class_edges - class_edges)
    elevation_ordering_ok = _ordering_passes(elevation_means, fixture["elevation_order"])
    sediment_ordering_ok = _ordering_passes(sediment_means, fixture["sediment_order"])

    acceptance = {
        "expected_classes_present": not missing_classes,
        "expected_edges_present": not missing_class_edges,
        "parent_process_links_present": not unparented,
        "province_ids_contiguous_where_expected": not non_contiguous,
        "no_random_texture_only_provinces": not random_texture_only,
        "expected_elevation_ordering": elevation_ordering_ok,
        "expected_sediment_ordering": sediment_ordering_ok,
        "all_grid_ids_declared": not missing_declarations,
        "all_declared_ids_used": not unused_declarations,
    }
    cell_count = sum(len(row) for row in grid)
    province_cell_counts = Counter(value for row in grid for value in row)
    largest_province_fraction = max(province_cell_counts.values()) / cell_count
    return {
        "name": fixture["name"],
        "passed": all(acceptance.values()),
        "acceptance": acceptance,
        "metrics": {
            "cell_count": cell_count,
            "province_count": len(grid_ids),
            "province_class_count": len(classes_present),
            "class_edge_count": len(class_edges),
            "province_edge_count": len(province_edges),
            "parent_process_count": len(parent_processes),
            "largest_province_fraction": largest_province_fraction,
            "unparented_province_count": len(unparented),
            "random_texture_only_province_count": len(random_texture_only),
            "non_contiguous_province_count": len(non_contiguous),
            "missing_expected_class_count": len(missing_classes),
            "missing_expected_edge_count": len(missing_class_edges),
        },
        "classes_present": classes_present,
        "parent_processes": parent_processes,
        "province_edges": sorted(province_edges),
        "class_edges": sorted(class_edges),
        "elevation_class_means_m": elevation_means,
        "sediment_class_means_m": sediment_means,
        "missing_expected_classes": missing_classes,
        "missing_expected_class_edges": missing_class_edges,
        "non_contiguous_provinces": non_contiguous,
        "random_texture_only_provinces": random_texture_only,
        "unparented_provinces": unparented,
    }


def province_graph_fixture_summary() -> dict[str, Any]:
    fixtures = [_analyze_fixture(fixture) for fixture in FIXTURES]
    fixture_names = [fixture["name"] for fixture in fixtures]
    classes_covered = sorted({
        province_class
        for fixture in fixtures
        for province_class in fixture["classes_present"]
    })
    parent_processes_covered = sorted({
        process
        for fixture in fixtures
        for process in fixture["parent_processes"]
    })
    acceptance = {
        "expected_fixture_suite_complete": tuple(fixture_names) == REQUIRED_FIXTURES,
        "all_fixtures_pass": all(fixture["passed"] for fixture in fixtures),
        "required_province_classes_covered": REQUIRED_PROVINCE_CLASSES.issubset(classes_covered),
        "required_parent_processes_covered": REQUIRED_PARENT_PROCESSES.issubset(parent_processes_covered),
        "parent_process_links_present": all(
            fixture["acceptance"]["parent_process_links_present"] for fixture in fixtures),
        "province_ids_contiguous_where_expected": all(
            fixture["acceptance"]["province_ids_contiguous_where_expected"] for fixture in fixtures),
        "no_random_texture_only_provinces": all(
            fixture["acceptance"]["no_random_texture_only_provinces"] for fixture in fixtures),
        "expected_elevation_ordering": all(
            fixture["acceptance"]["expected_elevation_ordering"] for fixture in fixtures),
        "expected_sediment_ordering": all(
            fixture["acceptance"]["expected_sediment_ordering"] for fixture in fixtures),
    }
    return {
        "schema": SCHEMA,
        "status": "fixture_suite_ready" if all(acceptance.values()) else "fixture_suite_incomplete",
        "fixture_count": len(fixtures),
        "province_class_count": len(classes_covered),
        "parent_process_count": len(parent_processes_covered),
        "classes_covered": classes_covered,
        "parent_processes_covered": parent_processes_covered,
        "required_fixture_names": REQUIRED_FIXTURES,
        "required_province_classes": sorted(REQUIRED_PROVINCE_CLASSES),
        "required_parent_processes": sorted(REQUIRED_PARENT_PROCESSES),
        "fixtures": fixtures,
        "acceptance": acceptance,
        "next_gates": (
            "P72.generated_world_province_diversity_gate",
            "P73.real_earth_case_study_calibration",
            "P74.terrain_coupling_rewrite",
        ),
    }
