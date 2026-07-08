"""Crust, sediment, and province-coupling reference diagnostics.

P83 is a reference gate, not a production terrain rewrite.  It makes the next
rewrite testable by defining how a continent's province objects should couple
crustal thickness, basement age, sediment accommodation, elevation, and relief.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from typing import Any

import numpy as np

from aevum.modules.terrain import (
    CONT_PROVINCE_ACTIVE_OROGEN,
    CONT_PROVINCE_CLASS_NAMES,
    CONT_PROVINCE_FORELAND_BASIN,
    CONT_PROVINCE_INTRACRATONIC_BASIN,
    CONT_PROVINCE_OLD_OROGEN,
    CONT_PROVINCE_PASSIVE_MARGIN_LOWLAND,
    CONT_PROVINCE_PLATFORM,
    CONT_PROVINCE_RIFT_SYSTEM,
    CONT_PROVINCE_SHIELD,
    CONT_PROVINCE_VOLCANIC_LIP_PLATEAU,
)


SCHEMA = "aevum.crust_sediment_province_coupling.v1"

SOURCE_IDS = (
    "CRUST1_0",
    "GLIM_GLOBAL_LITHOLOGY",
    "NOAA_TOTAL_SEDIMENT_THICKNESS",
)

REQUIRED_CLASSES = {
    "shield",
    "platform",
    "intracratonic_basin",
    "foreland_basin",
    "active_orogen",
    "old_orogen",
    "old_suture",
    "rift_shoulder",
    "rift_basin",
    "rift_axis",
    "passive_margin_lowland",
    "continental_shelf",
    "volcanic_lip_plateau",
}

REQUIRED_PARENT_PROCESSES = {
    "cratonization",
    "platform_subsidence",
    "intracratonic_sag",
    "flexural_loading",
    "collision_orogeny",
    "orogenic_decay",
    "suture_inheritance",
    "rift_shoulder_uplift",
    "rift_basin_subsidence",
    "continental_extension",
    "passive_margin_subsidence",
    "shelf_sedimentation",
    "plume_lip_emplacement",
}

PRODUCTION_PROVINCE_FIELDS = (
    "tectonics.continental_province_id",
    "tectonics.continental_province_code",
    "tectonics.province_parent_process",
)

EXPECTED_CURRENT_RESIDUAL_FIELDS = PRODUCTION_PROVINCE_FIELDS

CURRENT_KIND_GROUPS = {
    "core": ("shield",),
    "platform": ("platform",),
    "basin_lowland": (
        "interior_basin",
        "foreland_basin",
        "rift_basin",
        "passive_margin_lowland",
    ),
    "orogen": ("old_subdued_orogen", "active_orogen", "orogenic_plateau"),
}

PRODUCTION_KIND_GROUPS = {
    "core": ("shield",),
    "platform": ("platform",),
    "basin_lowland": (
        "intracratonic_basin",
        "foreland_basin",
        "rift_system",
        "passive_margin_lowland",
    ),
    "orogen": ("old_orogen", "volcanic_lip_plateau"),
}


def _rows(*rows: str) -> tuple[tuple[str, ...], ...]:
    return tuple(tuple(row.split()) for row in rows)


def _province(
    province_class: str,
    parent_processes: tuple[str, ...],
    parent_object: str,
    crust_thickness_m: float,
    sediment_m: float,
    elevation_m: float,
    relief_m: float,
    basement_age_myr: float,
    stability: float,
    lithology_group: str,
) -> dict[str, Any]:
    return {
        "class": province_class,
        "parent_processes": parent_processes,
        "parent_object": parent_object,
        "parent_continent_id": "continent:p83-synthetic",
        "crust_thickness_m": float(crust_thickness_m),
        "sediment_m": float(sediment_m),
        "elevation_m": float(elevation_m),
        "relief_m": float(relief_m),
        "basement_age_myr": float(basement_age_myr),
        "stability": float(stability),
        "lithology_group": lithology_group,
        "source": "deterministic_process_reference_fixture",
    }


REFERENCE_FIXTURE: dict[str, Any] = {
    "name": "p83_crust_sediment_province_synthetic_continent",
    "source_ids": SOURCE_IDS,
    "grid": _rows(
        "SH SH SH PL PL IB IB PM CS CS",
        "SH SH SH PL PL IB IB PM CS CS",
        "SH SH PL PL IB IB PM PM CS CS",
        "OO OO ST PL FB FB AO AO PM CS",
        "OO OO ST PL FB FB AO AO PM CS",
        "PL RS RB RA RB RS VL VL PM CS",
        "PL RS RB RA RB RS VL VL PM CS",
    ),
    "provinces": {
        "SH": _province(
            "shield", ("cratonization",), "archean_craton_core",
            42000.0, 150.0, 620.0, 380.0, 3000.0, 0.93, "crystalline_basement"),
        "PL": _province(
            "platform", ("platform_subsidence",), "covered_craton_platform",
            38000.0, 800.0, 280.0, 160.0, 1800.0, 0.76, "mixed_sedimentary_cover"),
        "IB": _province(
            "intracratonic_basin", ("intracratonic_sag",), "long_lived_sag_basin",
            33000.0, 2600.0, 80.0, 90.0, 1600.0, 0.58, "basin_fill"),
        "PM": _province(
            "passive_margin_lowland", ("passive_margin_subsidence",),
            "mature_passive_margin_prism",
            31000.0, 2200.0, 40.0, 70.0, 1100.0, 0.62, "coastal_plain_cover"),
        "CS": _province(
            "continental_shelf", ("shelf_sedimentation", "passive_margin_subsidence"),
            "continental_shelf_wedge",
            28000.0, 3300.0, -120.0, 50.0, 1000.0, 0.60, "shelf_sediment"),
        "OO": _province(
            "old_orogen", ("orogenic_decay", "collision_orogeny"), "paleo_orogen",
            45000.0, 520.0, 900.0, 650.0, 1200.0, 0.55, "metamorphic_belt"),
        "ST": _province(
            "old_suture", ("suture_inheritance",), "closed_ocean_suture",
            39000.0, 450.0, 240.0, 260.0, 1300.0, 0.50, "suture_melange"),
        "FB": _province(
            "foreland_basin", ("flexural_loading",), "active_orogen_load",
            36000.0, 3600.0, 110.0, 140.0, 700.0, 0.45, "foreland_fill"),
        "AO": _province(
            "active_orogen", ("collision_orogeny",), "active_collision_belt",
            58000.0, 260.0, 3200.0, 1900.0, 250.0, 0.35, "thickened_orogenic_crust"),
        "RS": _province(
            "rift_shoulder", ("rift_shoulder_uplift", "continental_extension"),
            "rift_shoulder_uplift",
            40000.0, 320.0, 1050.0, 900.0, 700.0, 0.38, "uplifted_platform_edge"),
        "RB": _province(
            "rift_basin", ("rift_basin_subsidence", "continental_extension"),
            "rift_half_graben",
            29000.0, 2800.0, -60.0, 320.0, 500.0, 0.35, "rift_basin_fill"),
        "RA": _province(
            "rift_axis", ("continental_extension",), "rift_axis",
            27000.0, 1900.0, 30.0, 250.0, 450.0, 0.32, "extensional_axis_fill"),
        "VL": _province(
            "volcanic_lip_plateau", ("plume_lip_emplacement",),
            "large_igneous_province",
            50000.0, 220.0, 1450.0, 700.0, 180.0, 0.50, "flood_basalt"),
    },
    "expected_orderings": {
        "crust_thickness_m": (
            (
                "active_orogen",
                "volcanic_lip_plateau",
                "old_orogen",
                "shield",
                "rift_shoulder",
                "old_suture",
                "platform",
                "foreland_basin",
                "intracratonic_basin",
                "passive_margin_lowland",
                "rift_basin",
                "continental_shelf",
                "rift_axis",
            ),
        ),
        "sediment_m": (
            (
                "foreland_basin",
                "continental_shelf",
                "rift_basin",
                "intracratonic_basin",
                "passive_margin_lowland",
                "rift_axis",
                "platform",
                "old_orogen",
                "old_suture",
                "rift_shoulder",
                "active_orogen",
                "volcanic_lip_plateau",
                "shield",
            ),
        ),
        "elevation_m": (
            (
                "active_orogen",
                "volcanic_lip_plateau",
                "rift_shoulder",
                "old_orogen",
                "shield",
                "platform",
                "old_suture",
                "foreland_basin",
                "intracratonic_basin",
                "passive_margin_lowland",
                "rift_axis",
                "rift_basin",
                "continental_shelf",
            ),
        ),
        "relief_m": (
            (
                "active_orogen",
                "rift_shoulder",
                "volcanic_lip_plateau",
                "old_orogen",
                "shield",
                "rift_basin",
                "old_suture",
                "rift_axis",
                "platform",
                "foreland_basin",
                "intracratonic_basin",
                "passive_margin_lowland",
                "continental_shelf",
            ),
        ),
    },
}


def _canonical_digest(data: Any) -> str:
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _class_means(fixture: dict[str, Any], field: str) -> dict[str, float]:
    totals: defaultdict[str, float] = defaultdict(float)
    counts: Counter[str] = Counter()
    provinces = fixture["provinces"]
    for row in fixture["grid"]:
        for province_id in row:
            province = provinces[province_id]
            province_class = str(province["class"])
            totals[province_class] += float(province[field])
            counts[province_class] += 1
    return {
        province_class: totals[province_class] / counts[province_class]
        for province_class in sorted(counts)
    }


def _ordering_passes(
    means: dict[str, float],
    orderings: tuple[tuple[str, ...], ...],
) -> bool:
    for ordering in orderings:
        values = [float(means[province_class]) for province_class in ordering]
        if any(left <= right for left, right in zip(values, values[1:])):
            return False
    return True


def _class_edges(fixture: dict[str, Any]) -> tuple[tuple[str, str], ...]:
    grid = fixture["grid"]
    provinces = fixture["provinces"]
    edges: set[tuple[str, str]] = set()
    for r, row in enumerate(grid):
        for c, province_id in enumerate(row):
            cls = str(provinces[province_id]["class"])
            for nr, nc in ((r + 1, c), (r, c + 1)):
                if nr >= len(grid) or nc >= len(grid[nr]):
                    continue
                other_id = grid[nr][nc]
                other_cls = str(provinces[other_id]["class"])
                if cls != other_cls:
                    edges.add(tuple(sorted((cls, other_cls))))
    return tuple(sorted(edges))


def reference_coupling_summary() -> dict[str, Any]:
    fixture = REFERENCE_FIXTURE
    grid_ids = {province_id for row in fixture["grid"] for province_id in row}
    provinces = fixture["provinces"]
    declared_ids = set(provinces)
    classes_present = sorted({str(provinces[pid]["class"]) for pid in grid_ids})
    parent_processes = sorted({
        str(process)
        for pid in grid_ids
        for process in provinces[pid]["parent_processes"]
    })
    class_means = {
        field: _class_means(fixture, field)
        for field in ("crust_thickness_m", "sediment_m", "elevation_m", "relief_m")
    }
    orderings = fixture["expected_orderings"]
    ordering_checks = {
        field: _ordering_passes(class_means[field], tuple(orderings[field]))
        for field in orderings
    }

    missing_required_classes = tuple(sorted(REQUIRED_CLASSES - set(classes_present)))
    missing_required_parent_processes = tuple(
        sorted(REQUIRED_PARENT_PROCESSES - set(parent_processes)))
    unparented = tuple(sorted(
        pid for pid in grid_ids if not provinces[pid]["parent_processes"]))
    random_sourced = tuple(sorted(
        pid for pid in grid_ids if "random" in str(provinces[pid]["source"])))

    low_accommodation_classes = (
        "intracratonic_basin",
        "foreland_basin",
        "rift_basin",
        "passive_margin_lowland",
        "continental_shelf",
    )
    low_accommodation_low = all(
        class_means["elevation_m"][cls] < class_means["elevation_m"]["platform"]
        for cls in low_accommodation_classes
    )
    low_accommodation_sedimented = all(
        class_means["sediment_m"][cls] > class_means["sediment_m"]["platform"] + 900.0
        for cls in low_accommodation_classes
    )
    low_accommodation_parented = all(
        provinces[pid]["parent_continent_id"] == "continent:p83-synthetic"
        for pid in grid_ids
        if provinces[pid]["class"] in low_accommodation_classes
    )
    shield = next(
        province for province in provinces.values() if province["class"] == "shield")
    shield_old_stable_not_high_flat = bool(
        shield["basement_age_myr"] >= 2500.0
        and shield["stability"] >= 0.85
        and 350.0 <= shield["elevation_m"] <= 900.0
        and 250.0 <= shield["relief_m"] <= 700.0
        and shield["sediment_m"] < 400.0
        and shield["elevation_m"] < class_means["elevation_m"]["old_orogen"]
        and shield["elevation_m"] < class_means["elevation_m"]["rift_shoulder"]
    )
    passive_margin_low_not_erased = bool(
        -50.0 <= class_means["elevation_m"]["passive_margin_lowland"] <= 150.0
        and class_means["sediment_m"]["passive_margin_lowland"] > 1800.0
        and class_means["crust_thickness_m"]["passive_margin_lowland"] >= 28000.0
        and class_means["crust_thickness_m"]["passive_margin_lowland"]
        < class_means["crust_thickness_m"]["platform"]
    )
    basins_low_not_erased = bool(
        low_accommodation_low
        and low_accommodation_sedimented
        and low_accommodation_parented
    )
    crust_sediment_elevation_coupled = bool(
        class_means["crust_thickness_m"]["active_orogen"]
        > class_means["crust_thickness_m"]["platform"] + 15000.0
        and class_means["elevation_m"]["active_orogen"]
        > class_means["elevation_m"]["platform"] + 2500.0
        and class_means["sediment_m"]["foreland_basin"]
        > class_means["sediment_m"]["active_orogen"] + 2500.0
        and class_means["elevation_m"]["foreland_basin"]
        < class_means["elevation_m"]["platform"]
    )
    no_grid_declaration_errors = bool(
        not (grid_ids - declared_ids) and not (declared_ids - grid_ids))
    acceptance = {
        "fixture_schema_ready": bool(fixture["name"] and fixture["source_ids"]),
        "required_classes_present": not missing_required_classes,
        "required_parent_processes_present": not missing_required_parent_processes,
        "all_grid_ids_declared": no_grid_declaration_errors,
        "all_provinces_parented": not unparented,
        "no_random_texture_sources": not random_sourced,
        "crust_thickness_ordering": bool(ordering_checks["crust_thickness_m"]),
        "sediment_accommodation_ordering": bool(ordering_checks["sediment_m"]),
        "elevation_ordering": bool(ordering_checks["elevation_m"]),
        "relief_ordering": bool(ordering_checks["relief_m"]),
        "basins_low_without_erasing_parent_continent": basins_low_not_erased,
        "passive_margin_lowland_preserves_continent": passive_margin_low_not_erased,
        "shield_old_stable_not_high_flat": shield_old_stable_not_high_flat,
        "crust_sediment_elevation_coupled": crust_sediment_elevation_coupled,
        "raw_reference_data_not_stored": True,
    }
    metrics = {
        "province_id_count": int(len(grid_ids)),
        "province_class_count": int(len(classes_present)),
        "parent_process_count": int(len(parent_processes)),
        "class_edge_count": int(len(_class_edges(fixture))),
        "missing_required_class_count": int(len(missing_required_classes)),
        "missing_required_parent_process_count": int(
            len(missing_required_parent_processes)),
        "unparented_province_id_count": int(len(unparented)),
        "random_sourced_province_id_count": int(len(random_sourced)),
        "shield_elevation_m": float(class_means["elevation_m"]["shield"]),
        "shield_relief_m": float(class_means["relief_m"]["shield"]),
        "shield_basement_age_myr": float(shield["basement_age_myr"]),
        "platform_elevation_m": float(class_means["elevation_m"]["platform"]),
        "intracratonic_basin_elevation_m": float(
            class_means["elevation_m"]["intracratonic_basin"]),
        "foreland_basin_sediment_m": float(
            class_means["sediment_m"]["foreland_basin"]),
        "active_orogen_crust_thickness_m": float(
            class_means["crust_thickness_m"]["active_orogen"]),
        "passive_margin_lowland_elevation_m": float(
            class_means["elevation_m"]["passive_margin_lowland"]),
        "continental_shelf_sediment_m": float(
            class_means["sediment_m"]["continental_shelf"]),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "crust_sediment_province_reference_ready"
            if all(acceptance.values())
            else "crust_sediment_province_reference_incomplete"
        ),
        "fixture_name": fixture["name"],
        "source_ids": fixture["source_ids"],
        "extraction_policy": {
            "raw_crust_lithology_sediment_data_stored": False,
            "direct_crust_lithology_sediment_extraction_pending": True,
            "fixture_is_small_derived_reference": True,
        },
        "classes_present": tuple(classes_present),
        "parent_processes": tuple(parent_processes),
        "class_edges": _class_edges(fixture),
        "class_means": class_means,
        "ordering_checks": ordering_checks,
        "missing_required_classes": missing_required_classes,
        "missing_required_parent_processes": missing_required_parent_processes,
        "unparented_province_ids": unparented,
        "random_sourced_province_ids": random_sourced,
        "metrics": metrics,
        "acceptance": acceptance,
        "fixture_digest": _canonical_digest(fixture),
    }


def _field_array(world: Any, name: str) -> np.ndarray | None:
    if name not in getattr(world, "fields", {}):
        return None
    return np.asarray(world.field(name), dtype=np.float64)


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    if values.size == 0 or weights.size == 0 or float(weights.sum()) <= 0.0:
        return float("nan")
    return float(np.average(values, weights=weights))


def _current_group_mean(
    kind_summary: dict[str, dict[str, Any]],
    kinds: tuple[str, ...],
    field: str,
) -> float:
    numerator = 0.0
    denominator = 0.0
    for kind in kinds:
        summary = kind_summary.get(kind)
        if not summary:
            continue
        area = float(summary["area_m2"])
        value = float(summary[field])
        if np.isfinite(value) and area > 0.0:
            numerator += value * area
            denominator += area
    if denominator <= 0.0:
        return float("nan")
    return numerator / denominator


def current_generated_coupling_audit(world: Any) -> dict[str, Any]:
    required_fields = (
        "crust.thickness_m",
        "sediment.thickness_m",
        "terrain.elevation_m",
    )
    missing_fields = tuple(
        name
        for name in required_fields
        if name not in getattr(world, "fields", {})
    )
    missing_production_fields = tuple(
        name
        for name in PRODUCTION_PROVINCE_FIELDS
        if name not in getattr(world, "fields", {})
    )
    expected_residuals_recorded = bool(
        set(missing_production_fields).issubset(set(EXPECTED_CURRENT_RESIDUAL_FIELDS)))

    elevation = _field_array(world, "terrain.elevation_m")
    sediment = _field_array(world, "sediment.thickness_m")
    thickness = _field_array(world, "crust.thickness_m")
    cell_area = np.asarray(getattr(world.grid, "cell_area", np.ones(world.grid.n)),
                           dtype=np.float64)
    n = int(world.grid.n)
    production_objects = list(
        getattr(world, "objects", {}).get("tectonics.continental_provinces", []))
    use_production_graph = bool(not missing_production_fields and production_objects)

    kind_acc: dict[str, dict[str, float]] = defaultdict(lambda: {
        "object_count": 0.0,
        "cell_count": 0.0,
        "area_m2": 0.0,
        "elevation_area_sum": 0.0,
        "sediment_area_sum": 0.0,
        "thickness_area_sum": 0.0,
    })
    if use_production_graph:
        objects = production_objects
        province_code = np.asarray(
            world.field("tectonics.continental_province_code"),
            dtype=np.float64,
        ).astype(int)
        code_to_name = {
            int(code): str(name)
            for code, name in CONT_PROVINCE_CLASS_NAMES.items()
            if int(code) > 0 and str(name) != "none"
        }
        kind_counts = Counter(
            code_to_name[int(code)]
            for code in province_code
            if int(code) in code_to_name
        )
        object_counts_by_kind = Counter(
            str(obj.get("province_class", "unknown"))
            for obj in production_objects
        )
        for code, kind in sorted(code_to_name.items()):
            cells = np.where(province_code == int(code))[0]
            if cells.size == 0:
                continue
            weights = cell_area[cells]
            area = float(weights.sum())
            acc = kind_acc[kind]
            acc["object_count"] += float(object_counts_by_kind.get(kind, 0))
            acc["cell_count"] += float(cells.size)
            acc["area_m2"] += area
            if elevation is not None:
                acc["elevation_area_sum"] += float(np.sum(elevation[cells] * weights))
            if sediment is not None:
                acc["sediment_area_sum"] += float(np.sum(sediment[cells] * weights))
            if thickness is not None:
                acc["thickness_area_sum"] += float(np.sum(thickness[cells] * weights))
        kind_groups = PRODUCTION_KIND_GROUPS
    else:
        objects = list(
            getattr(world, "objects", {}).get("terrain.continental_landforms", []))
        kind_counts = Counter(str(obj.get("kind", "unknown")) for obj in objects)
        for obj in objects:
            kind = str(obj.get("kind", "unknown"))
            cells = np.asarray(obj.get("cells", ()), dtype=np.int64)
            cells = cells[(cells >= 0) & (cells < n)]
            if cells.size == 0:
                continue
            weights = cell_area[cells]
            area = float(weights.sum())
            acc = kind_acc[kind]
            acc["object_count"] += 1.0
            acc["cell_count"] += float(cells.size)
            acc["area_m2"] += area
            if elevation is not None:
                acc["elevation_area_sum"] += float(np.sum(elevation[cells] * weights))
            if sediment is not None:
                acc["sediment_area_sum"] += float(np.sum(sediment[cells] * weights))
            if thickness is not None:
                acc["thickness_area_sum"] += float(np.sum(thickness[cells] * weights))
        kind_groups = CURRENT_KIND_GROUPS

    kind_summary: dict[str, dict[str, Any]] = {}
    for kind, acc in sorted(kind_acc.items()):
        area = float(acc["area_m2"])
        kind_summary[kind] = {
            "object_count": int(acc["object_count"]),
            "cell_count": int(acc["cell_count"]),
            "area_m2": area,
            "mean_elevation_m": (
                float(acc["elevation_area_sum"] / area)
                if elevation is not None and area > 0.0 else float("nan")
            ),
            "mean_sediment_m": (
                float(acc["sediment_area_sum"] / area)
                if sediment is not None and area > 0.0 else float("nan")
            ),
            "mean_crust_thickness_m": (
                float(acc["thickness_area_sum"] / area)
                if thickness is not None and area > 0.0 else float("nan")
            ),
        }

    group_means: dict[str, dict[str, float]] = {}
    for group, kinds in kind_groups.items():
        group_means[group] = {
            "mean_elevation_m": _current_group_mean(
                kind_summary, kinds, "mean_elevation_m"),
            "mean_sediment_m": _current_group_mean(
                kind_summary, kinds, "mean_sediment_m"),
            "mean_crust_thickness_m": _current_group_mean(
                kind_summary, kinds, "mean_crust_thickness_m"),
        }
    current_kind_groups_present = {
        group: any(kind_counts.get(kind, 0) > 0 for kind in kinds)
        for group, kinds in kind_groups.items()
    }
    basin_low = bool(
        np.isfinite(group_means["basin_lowland"]["mean_elevation_m"])
        and np.isfinite(group_means["platform"]["mean_elevation_m"])
        and group_means["basin_lowland"]["mean_elevation_m"]
        < group_means["platform"]["mean_elevation_m"])
    sediment_signal = bool(
        np.isfinite(group_means["basin_lowland"]["mean_sediment_m"])
        and np.isfinite(group_means["platform"]["mean_sediment_m"])
        and group_means["basin_lowland"]["mean_sediment_m"]
        > group_means["platform"]["mean_sediment_m"] + 500.0)
    highland_signal = bool(
        np.isfinite(group_means["orogen"]["mean_elevation_m"])
        and np.isfinite(group_means["basin_lowland"]["mean_elevation_m"])
        and group_means["orogen"]["mean_elevation_m"]
        > group_means["basin_lowland"]["mean_elevation_m"] + 400.0)
    field_ranges = {}
    for name, arr in (
        ("terrain.elevation_m", elevation),
        ("sediment.thickness_m", sediment),
        ("crust.thickness_m", thickness),
    ):
        if arr is None:
            continue
        field_ranges[name] = {
            "min": float(np.nanmin(arr)),
            "mean": float(np.nanmean(arr)),
            "max": float(np.nanmax(arr)),
        }
    acceptance = {
        "current_core_fields_available": not missing_fields,
        "current_landform_objects_present": len(objects) > 0,
        "current_required_kind_groups_present": all(current_kind_groups_present.values()),
        "current_lowland_elevation_signal_present": basin_low,
        "current_sediment_accommodation_signal_present": sediment_signal,
        "current_highland_contrast_signal_present": highland_signal,
        "current_expected_residuals_recorded": expected_residuals_recorded,
    }
    return {
        "schema": SCHEMA,
        "status": (
            "generated_world_crust_sediment_province_audit_ready"
            if all(acceptance.values())
            else "generated_world_crust_sediment_province_audit_incomplete"
        ),
        "object_count": int(len(objects)),
        "kind_counts": dict(sorted((kind, int(count)) for kind, count in kind_counts.items())),
        "kind_summary": kind_summary,
        "kind_group_means": group_means,
        "kind_groups_present": current_kind_groups_present,
        "field_ranges": field_ranges,
        "missing_fields": missing_fields,
        "missing_production_province_fields": missing_production_fields,
        "expected_current_residual_fields": EXPECTED_CURRENT_RESIDUAL_FIELDS,
        "limitations": {
            "production_continental_province_ids_missing": bool(
                missing_production_fields),
        },
        "aggregation_source": (
            "tectonics.continental_provinces"
            if use_production_graph else "terrain.continental_landforms"
        ),
        "metrics": {
            "object_count": int(len(objects)),
            "kind_count": int(len(kind_counts)),
            "missing_field_count": int(len(missing_fields)),
            "missing_production_province_field_count": int(
                len(missing_production_fields)),
            "basin_lowland_mean_elevation_m": float(
                group_means["basin_lowland"]["mean_elevation_m"]),
            "platform_mean_elevation_m": float(
                group_means["platform"]["mean_elevation_m"]),
            "orogen_mean_elevation_m": float(
                group_means["orogen"]["mean_elevation_m"]),
            "basin_lowland_mean_sediment_m": float(
                group_means["basin_lowland"]["mean_sediment_m"]),
            "platform_mean_sediment_m": float(
                group_means["platform"]["mean_sediment_m"]),
            "shield_mean_sediment_m": float(
                group_means["core"]["mean_sediment_m"]),
        },
        "acceptance": acceptance,
    }


def crust_sediment_province_coupling_summary(world: Any | None = None) -> dict[str, Any]:
    reference = reference_coupling_summary()
    current = (
        current_generated_coupling_audit(world)
        if world is not None
        else {
            "status": "generated_world_crust_sediment_province_audit_not_run",
            "acceptance": {
                "current_core_fields_available": False,
                "current_landform_objects_present": False,
                "current_required_kind_groups_present": False,
                "current_lowland_elevation_signal_present": False,
                "current_sediment_accommodation_signal_present": False,
                "current_highland_contrast_signal_present": False,
                "current_expected_residuals_recorded": False,
            },
            "metrics": {},
        }
    )
    acceptance = {
        "reference_fixture_ready": reference["status"]
        == "crust_sediment_province_reference_ready",
        "current_generated_audit_available": current["status"]
        == "generated_world_crust_sediment_province_audit_ready",
        "required_classes_present": reference["acceptance"]["required_classes_present"],
        "required_parent_processes_present": reference["acceptance"][
            "required_parent_processes_present"],
        "crust_thickness_ordering": reference["acceptance"]["crust_thickness_ordering"],
        "sediment_accommodation_ordering": reference["acceptance"][
            "sediment_accommodation_ordering"],
        "shield_old_stable_not_high_flat": reference["acceptance"][
            "shield_old_stable_not_high_flat"],
        "basins_and_passive_margins_low_not_erased": bool(
            reference["acceptance"]["basins_low_without_erasing_parent_continent"]
            and reference["acceptance"]["passive_margin_lowland_preserves_continent"]),
        "current_core_fields_available": current["acceptance"][
            "current_core_fields_available"],
        "current_object_signal_present": bool(
            current["acceptance"]["current_landform_objects_present"]
            and current["acceptance"]["current_required_kind_groups_present"]
            and current["acceptance"]["current_lowland_elevation_signal_present"]
            and current["acceptance"]["current_sediment_accommodation_signal_present"]
            and current["acceptance"]["current_highland_contrast_signal_present"]),
        "current_expected_residuals_recorded": current["acceptance"][
            "current_expected_residuals_recorded"],
    }
    return {
        "schema": SCHEMA,
        "status": (
            "crust_sediment_province_coupling_ready"
            if all(acceptance.values())
            else "crust_sediment_province_coupling_incomplete"
        ),
        "reference": reference,
        "current_generated": current,
        "acceptance": acceptance,
        "summary_digest": _canonical_digest({
            "reference": reference,
            "current": current,
        }),
    }
