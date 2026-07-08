"""Reference corpus and metric scaffold for continental physiography work.

This module is intentionally offline and deterministic.  It records which
real-Earth references, theory layers, benchmark metrics, and current baseline
gaps must stay visible while the terrain and plate systems are rebuilt.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any


SCHEMA = "aevum.physiographic_reference.v1"

PHASES: tuple[dict[str, Any], ...] = (
    {
        "id": "phase_0",
        "name": "current_pipeline_inventory",
        "goal": "Freeze the current asset, metric, and code-path inventory before new changes.",
        "requires": (),
    },
    {
        "id": "phase_1",
        "name": "real_earth_reference_corpus",
        "goal": "Collect global topography, bathymetry, coastline, hydrography, and physiographic references.",
        "requires": ("phase_0",),
    },
    {
        "id": "phase_2",
        "name": "geologic_and_tectonic_province_corpus",
        "goal": "Collect crust, lithology, province, plate-boundary, fault, and reconstruction references.",
        "requires": ("phase_1",),
    },
    {
        "id": "phase_3",
        "name": "theory_model",
        "goal": "Convert Earth-science references into process-parented province and terrain rules.",
        "requires": ("phase_1", "phase_2"),
    },
    {
        "id": "phase_4",
        "name": "implementation_architecture",
        "goal": "Define province graph, drainage, ocean-margin, and terrain-coupling code boundaries.",
        "requires": ("phase_3",),
    },
    {
        "id": "phase_5",
        "name": "microbenchmarks",
        "goal": "Create fixture and generated-world gates for each province and terrain class.",
        "requires": ("phase_4",),
    },
    {
        "id": "phase_6",
        "name": "optimization_targets",
        "goal": "Set Earth-calibrated envelopes and tuning dashboards for generated worlds.",
        "requires": ("phase_5",),
    },
    {
        "id": "phase_7",
        "name": "archive_and_execution",
        "goal": "Archive decisions and keep the implementation plan executable through P70-P75.",
        "requires": ("phase_0", "phase_6"),
    },
)

REFERENCE_SOURCES: tuple[dict[str, Any], ...] = (
    {
        "id": "AEVUM_P69_GENERATED_BASELINE",
        "category": "generated_baseline",
        "phases": ("phase_0", "phase_7"),
        "url": "repo://out_bench_p69_highres_physical_ensemble_visual_audit_20260626_c",
        "uses": ("current_pipeline_inventory", "archive_execution_baseline", "p70_gap_mapping"),
    },
    {
        "id": "NOAA_ETOPO_2022",
        "category": "global_relief",
        "phases": ("phase_1", "phase_6"),
        "url": "https://www.ncei.noaa.gov/products/etopo-global-relief-model",
        "uses": ("hypsometry", "continental_relief", "ocean_bathymetry"),
    },
    {
        "id": "GEBCO_GRIDDED_BATHYMETRY",
        "category": "global_relief",
        "phases": ("phase_1", "phase_6"),
        "url": "https://www.gebco.net/data-products/gridded-bathymetry-data",
        "uses": ("shelf_slope_rise_abyss", "ridge_trench_depths"),
    },
    {
        "id": "NATURAL_EARTH_10M_PHYSICAL",
        "category": "planform",
        "phases": ("phase_1", "phase_6"),
        "url": "https://www.naturalearthdata.com/downloads/10m-physical-vectors/",
        "uses": ("coastline_complexity", "land_components", "island_arcs"),
    },
    {
        "id": "USGS_PHYSIOGRAPHIC_DIVISIONS_US",
        "category": "physiographic_provinces",
        "phases": ("phase_1", "phase_3", "phase_5"),
        "url": "https://data.usgs.gov/datacatalog/data/USGS%3Ae04ea9e9-17b6-45ae-b279-7bc35ea79539",
        "uses": ("province_graph_fixtures", "continental_internal_divisions"),
    },
    {
        "id": "NPS_PHYSIOGRAPHIC_PROVINCES",
        "category": "physiographic_provinces",
        "phases": ("phase_1", "phase_3", "phase_5"),
        "url": "https://www.nps.gov/subjects/geology/physiographic-provinces.htm",
        "uses": ("province_class_examples", "public_description_crosscheck"),
    },
    {
        "id": "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
        "category": "geologic_provinces",
        "phases": ("phase_2", "phase_3", "phase_5"),
        "url": "https://zenodo.org/records/6586972",
        "uses": ("province_parentage", "plate_context", "craton_orogen_basin_cases"),
    },
    {
        "id": "USGS_WORLD_GEOLOGIC_PROVINCES",
        "category": "geologic_provinces",
        "phases": ("phase_2", "phase_3"),
        "url": "https://data.usgs.gov/datacatalog/",
        "uses": ("global_province_names", "sedimentary_basin_cases"),
    },
    {
        "id": "GPLATES",
        "category": "plate_reconstruction",
        "phases": ("phase_2", "phase_3", "phase_4"),
        "url": "https://www.gplates.org/",
        "uses": ("plate_history_workflow", "reconstruction_validation"),
    },
    {
        "id": "PYGPLATES_EARTHBYTE",
        "category": "plate_reconstruction",
        "phases": ("phase_2", "phase_4", "phase_5"),
        "url": "https://www.earthbyte.org/category/resources/software-workflows/pygplates/",
        "uses": ("fixture_generation", "rotation_model_interop"),
    },
    {
        "id": "EARTHBYTE_RECONSTRUCTIONS",
        "category": "plate_reconstruction",
        "phases": ("phase_2", "phase_3", "phase_6"),
        "url": "https://www.earthbyte.org/category/reconstructions/",
        "uses": ("supercontinent_cycle_cases", "historical_process_coverage"),
    },
    {
        "id": "PB2002_PLATE_BOUNDARIES",
        "category": "plate_boundaries",
        "phases": ("phase_2", "phase_5", "phase_6"),
        "url": "https://peterbird.name/oldFTP/PB2002/",
        "uses": ("ridge_trench_transform_lengths", "boundary_type_fixtures"),
    },
    {
        "id": "GEM_GLOBAL_ACTIVE_FAULTS",
        "category": "plate_boundaries",
        "phases": ("phase_2", "phase_5"),
        "url": "https://github.com/GEMScienceTools/gem-global-active-faults",
        "uses": ("active_fault_density", "intraplate_reactivation_cases"),
    },
    {
        "id": "CRUST1_0",
        "category": "crust",
        "phases": ("phase_2", "phase_3", "phase_6"),
        "url": "https://ds.iris.edu/ds/products/emc-crust10/",
        "uses": ("crust_thickness", "continental_vs_oceanic_crust"),
    },
    {
        "id": "GLIM_GLOBAL_LITHOLOGY",
        "category": "lithology",
        "phases": ("phase_2", "phase_3", "phase_5"),
        "url": "https://www.geo.uni-hamburg.de/en/geologie/forschung/aquatische-geochemie/glim.html",
        "uses": ("rock_type_province_expression", "erosion_resistance_cases"),
    },
    {
        "id": "NOAA_TOTAL_SEDIMENT_THICKNESS",
        "category": "sediment",
        "phases": ("phase_2", "phase_3", "phase_6"),
        "url": "https://www.ncei.noaa.gov/products/total-sediment-thickness-oceans-seas",
        "uses": ("passive_margin_prisms", "source_to_sink_closure"),
    },
    {
        "id": "HYDROSHEDS_HYDROBASINS_HYDRORIVERS",
        "category": "drainage",
        "phases": ("phase_1", "phase_5", "phase_6"),
        "url": "https://www.hydrosheds.org/products",
        "uses": ("drainage_basin_graphs", "continental_divides", "source_to_sink_paths"),
    },
    {
        "id": "GMBA_MOUNTAIN_INVENTORY",
        "category": "mountains",
        "phases": ("phase_1", "phase_3", "phase_6"),
        "url": "https://www.earthenv.org/mountains",
        "uses": ("orogen_extent", "mountain_area_fraction", "relief_cases"),
    },
)

METRIC_SCHEMA: dict[str, tuple[dict[str, Any], ...]] = {
    "planform": (
        {"key": "land_fraction", "unit": "fraction", "target": "earthlike envelope"},
        {"key": "land_component_count", "unit": "count", "target": "few major continents plus island arcs"},
        {"key": "largest_land_component_fraction", "unit": "fraction_of_land", "target": "no single oversize megacontinent by default"},
        {"key": "land_ribbon_fraction_gt_0_5", "unit": "fraction_of_land", "target": "penalize excessive long strips"},
        {"key": "land_coastline_complexity_largest", "unit": "dimensionless", "target": "bounded fractal coastline"},
        {"key": "major_continent_count", "unit": "count", "target": "multiple major continents"},
    ),
    "hypsometry": (
        {"key": "land_elevation_mean_m", "unit": "m", "target": "Earthlike mean continental elevation"},
        {"key": "land_elevation_p95_m", "unit": "m", "target": "controlled highland tail"},
        {"key": "high_land_fraction_gt2500m", "unit": "fraction_of_land", "target": "limited plateau and active orogen area"},
        {"key": "inland_relief_p90_p10_m", "unit": "m", "target": "non-flat but not uniformly high interiors"},
        {"key": "lowland_fraction_lt500m", "unit": "fraction_of_land", "target": "broad lowland and basin expression"},
        {"key": "lowland_fraction_lt1000m", "unit": "fraction_of_land", "target": "platform and coastal plains visible"},
    ),
    "province_architecture": (
        {"key": "province_count_per_major_continent", "unit": "count", "target": "several internal provinces per continent"},
        {"key": "province_class_diversity_per_major_continent", "unit": "count", "target": "shield/platform/basin/orogen/rift mix"},
        {"key": "largest_internal_province_fraction", "unit": "fraction_of_continent", "target": "avoid one uniform interior tile"},
        {"key": "shield_share_of_continental_crust", "unit": "fraction", "target": "old cratonic cores present but bounded"},
        {"key": "platform_share_of_continental_crust", "unit": "fraction", "target": "broad stable platforms"},
        {"key": "basin_share_of_continental_crust", "unit": "fraction", "target": "visible low sedimentary basins"},
        {"key": "old_orogen_share_of_land", "unit": "fraction", "target": "eroded inherited belts"},
        {"key": "rift_basin_share_of_land", "unit": "fraction", "target": "localized extensional systems"},
        {"key": "passive_margin_lowland_share_of_land", "unit": "fraction", "target": "wide passive margin plains where appropriate"},
        {"key": "active_orogen_or_plateau_fraction_of_land", "unit": "fraction", "target": "bounded active belts and plateaus"},
    ),
    "process_parentage": (
        {"key": "parented_highland_fraction", "unit": "fraction", "target": "highlands tied to collision, rift shoulder, plume, or arc"},
        {"key": "unparented_highland_fraction", "unit": "fraction", "target": "near zero"},
        {"key": "orogen_foreland_adjacency_score", "unit": "score", "target": "orogens paired with foreland basins where expected"},
        {"key": "rift_axis_shoulder_adjacency_score", "unit": "score", "target": "rift axes paired with shoulders and basins"},
        {"key": "province_boundary_parentage_completeness", "unit": "fraction", "target": "internal boundaries carry process labels"},
    ),
    "drainage": (
        {"key": "drainage_basin_count", "unit": "count", "target": "multiple major basins per large continent"},
        {"key": "drainage_divide_province_boundary_alignment", "unit": "score", "target": "divides respond to orogens, shields, and rifts"},
        {"key": "source_to_sink_sediment_closure", "unit": "fraction", "target": "erosion and deposition approximately close"},
        {"key": "endorheic_basin_fraction", "unit": "fraction_of_land", "target": "limited but possible interior basins"},
        {"key": "lowland_fluvial_plain_fraction", "unit": "fraction_of_land", "target": "large rivers build plains and deltas"},
    ),
    "ocean_bathymetry": (
        {"key": "shelf_fraction_of_ocean", "unit": "fraction", "target": "continental shelves around passive margins"},
        {"key": "slope_rise_fraction_of_ocean", "unit": "fraction", "target": "slope/rise transition between shelf and abyss"},
        {"key": "abyss_fraction_of_ocean", "unit": "fraction", "target": "dominant but not overfilled abyssal plains"},
        {"key": "ridge_fraction_of_ocean", "unit": "fraction", "target": "connected mid-ocean ridge systems"},
        {"key": "trench_fraction_of_ocean", "unit": "fraction", "target": "localized subduction trenches"},
        {"key": "nearshore_deep_ocean_violation_fraction", "unit": "fraction", "target": "minimize ultra-deep water directly against continents"},
    ),
}

P69_EARTHLIKE_BASELINE: dict[str, Any] = {
    "source_run": "out_bench_p69_highres_physical_ensemble_visual_audit_20260626_c",
    "member": "earthlike_reference_physical_member",
    "cells": 8000,
    "metrics": {
        "land_fraction": 0.2475096054131667,
        "land_component_count": 2.0,
        "largest_land_component_fraction": 0.6141462207674615,
        "land_ribbon_fraction_gt_0_5": 0.4121528088382741,
        "land_coastline_complexity_largest": 17.435595774162692,
        "land_elevation_mean_m": 1022.1406995020552,
        "land_elevation_p95_m": 2596.835980242177,
        "high_land_fraction_gt2500m": 0.0535358015983353,
        "active_orogen_or_plateau_fraction_of_land": 0.07727272727272727,
        "shelf_fraction_of_ocean": 0.0523281786012958,
        "slope_rise_fraction_of_ocean": 0.1328795111377222,
        "abyss_fraction_of_ocean": 0.6865507959082437,
        "ridge_fraction_of_ocean": 0.07192700977953359,
        "trench_fraction_of_ocean": 0.037878888522216124,
        "shield_share_of_continental_crust": 0.1473660909843428,
        "platform_share_of_continental_crust": 0.060824787906354705,
        "basin_share_of_continental_crust": 0.04599237721835182,
    },
    "out_of_envelope": (
        "land_fraction",
        "largest_land_component_fraction",
        "land_component_count",
        "land_ribbon_fraction_gt_0_5",
        "land_coastline_complexity_largest",
        "abyss_fraction_of_ocean",
    ),
    "architectural_gaps": (
        "missing_first_class_multi_province_graph",
        "major_continents_not_required_to_have_multiple_internal_province_classes",
        "basin_lowland_expression_not_guaranteed_per_major_continent",
        "province_boundaries_not_yet_primary_terrain_drivers",
    ),
}


def _metric_to_group() -> dict[str, str]:
    mapping: dict[str, str] = {}
    for group, entries in METRIC_SCHEMA.items():
        for entry in entries:
            mapping[str(entry["key"])] = group
    return mapping


def source_inventory_summary() -> dict[str, Any]:
    phase_ids = {phase["id"] for phase in PHASES}
    categories = sorted({str(source["category"]) for source in REFERENCE_SOURCES})
    phase_counts: Counter[str] = Counter()
    category_counts: Counter[str] = Counter()
    ids: list[str] = []
    urls: list[str] = []
    invalid_phase_refs: dict[str, tuple[str, ...]] = {}
    for source in REFERENCE_SOURCES:
        ids.append(str(source["id"]))
        urls.append(str(source["url"]))
        category_counts[str(source["category"])] += 1
        source_phases = tuple(str(phase) for phase in source["phases"])
        for phase in source_phases:
            phase_counts[phase] += 1
        missing = tuple(phase for phase in source_phases if phase not in phase_ids)
        if missing:
            invalid_phase_refs[str(source["id"])] = missing

    covered_phases = {phase for phase, count in phase_counts.items() if count > 0}
    return {
        "source_count": len(REFERENCE_SOURCES),
        "category_count": len(categories),
        "categories": categories,
        "phase_coverage": dict(sorted(phase_counts.items())),
        "category_counts": dict(sorted(category_counts.items())),
        "duplicate_ids": sorted(key for key, count in Counter(ids).items() if count > 1),
        "duplicate_urls": sorted(key for key, count in Counter(urls).items() if count > 1),
        "invalid_phase_refs": invalid_phase_refs,
        "acceptance": {
            "has_minimum_source_count": len(REFERENCE_SOURCES) >= 14,
            "has_minimum_category_count": len(categories) >= 4,
            "covers_all_declared_phases": phase_ids.issubset(covered_phases),
            "has_unique_ids": len(set(ids)) == len(ids),
            "has_unique_urls": len(set(urls)) == len(urls),
            "has_no_invalid_phase_refs": not invalid_phase_refs,
        },
    }


def metric_schema_summary() -> dict[str, Any]:
    metric_to_group = _metric_to_group()
    grouped_counts = {group: len(entries) for group, entries in METRIC_SCHEMA.items()}
    duplicate_keys = sorted(
        key for key, count in Counter(metric["key"] for entries in METRIC_SCHEMA.values() for metric in entries).items()
        if count > 1)
    required_groups = {
        "planform",
        "hypsometry",
        "province_architecture",
        "process_parentage",
        "drainage",
        "ocean_bathymetry",
    }
    future_architecture_keys = {
        "province_count_per_major_continent",
        "province_class_diversity_per_major_continent",
        "largest_internal_province_fraction",
        "province_boundary_parentage_completeness",
        "drainage_basin_count",
        "drainage_divide_province_boundary_alignment",
        "source_to_sink_sediment_closure",
        "nearshore_deep_ocean_violation_fraction",
    }
    return {
        "metric_count": len(metric_to_group),
        "metric_group_count": len(METRIC_SCHEMA),
        "grouped_counts": grouped_counts,
        "duplicate_keys": duplicate_keys,
        "future_architecture_metric_count": len(future_architecture_keys & set(metric_to_group)),
        "acceptance": {
            "has_required_metric_groups": required_groups.issubset(METRIC_SCHEMA),
            "has_minimum_metric_count": len(metric_to_group) >= 30,
            "has_future_architecture_metrics": len(future_architecture_keys & set(metric_to_group)) >= 5,
            "has_unique_metric_keys": not duplicate_keys,
            "all_metrics_have_units": all(
                bool(entry.get("unit"))
                for entries in METRIC_SCHEMA.values()
                for entry in entries),
            "all_metrics_have_targets": all(
                bool(entry.get("target"))
                for entries in METRIC_SCHEMA.values()
                for entry in entries),
        },
    }


def p69_baseline_mapping() -> dict[str, Any]:
    metric_to_group = _metric_to_group()
    baseline_metrics = set(P69_EARTHLIKE_BASELINE["metrics"])
    mapped_by_group: defaultdict[str, list[str]] = defaultdict(list)
    for metric in sorted(baseline_metrics):
        if metric in metric_to_group:
            mapped_by_group[metric_to_group[metric]].append(metric)

    unmapped = sorted(metric for metric in baseline_metrics if metric not in metric_to_group)
    unmapped_out_of_envelope = sorted(
        metric for metric in P69_EARTHLIKE_BASELINE["out_of_envelope"] if metric not in metric_to_group)
    gaps = tuple(P69_EARTHLIKE_BASELINE["architectural_gaps"])
    return {
        "source_run": P69_EARTHLIKE_BASELINE["source_run"],
        "member": P69_EARTHLIKE_BASELINE["member"],
        "cells": P69_EARTHLIKE_BASELINE["cells"],
        "baseline_metric_count": len(baseline_metrics),
        "mapped_metric_count": len(baseline_metrics) - len(unmapped),
        "mapped_by_group": dict(sorted((group, metrics) for group, metrics in mapped_by_group.items())),
        "unmapped_metrics": unmapped,
        "out_of_envelope_metrics": tuple(P69_EARTHLIKE_BASELINE["out_of_envelope"]),
        "unmapped_out_of_envelope_metrics": unmapped_out_of_envelope,
        "architectural_gaps": gaps,
        "acceptance": {
            "all_baseline_metrics_mapped": not unmapped,
            "all_out_of_envelope_metrics_mapped": not unmapped_out_of_envelope,
            "explicitly_records_missing_province_graph": "missing_first_class_multi_province_graph" in gaps,
            "explicitly_records_terrain_coupling_gap": "province_boundaries_not_yet_primary_terrain_drivers" in gaps,
        },
    }


def p70_reference_scaffold_summary() -> dict[str, Any]:
    inventory = source_inventory_summary()
    metric_schema = metric_schema_summary()
    baseline = p69_baseline_mapping()
    phase_ids = [str(phase["id"]) for phase in PHASES]
    ready = (
        all(inventory["acceptance"].values())
        and all(metric_schema["acceptance"].values())
        and all(baseline["acceptance"].values())
        and len(phase_ids) == 8
    )
    return {
        "schema": SCHEMA,
        "status": "scaffold_ready" if ready else "scaffold_incomplete",
        "phases": PHASES,
        "source_inventory": inventory,
        "metric_schema": metric_schema,
        "p69_baseline_mapping": baseline,
        "next_gates": (
            "P71.province_graph_fixture_suite",
            "P72.generated_world_province_diversity_gate",
            "P73.real_earth_case_study_calibration",
            "P74.terrain_coupling_rewrite",
            "P75.release_promotion_audit",
        ),
    }
